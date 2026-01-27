#pragma once
#include <bitloop.h>
#include <vector>

#include "types.h"


SIM_BEG;

enum class MandelTransform
{
    NONE,
    FLATTEN
};

static const char* MandelFeatureNames[(int)KernelFeatures::COUNT] = {
    "NONE",
    "ITER",
    "DIST",
    "STRIPE"
};

template<typename T>
struct GammaLUT
{
    static constexpr int N = 1024;
    T table[N + 1];
    T invN; // 1.0f / N
};

template<typename T>
struct ToneFieldInfo
{
    T minv;
    T maxv;
    T mean;

    T invRange;  // 1 / (maxv - minv)
    T pivot;     // mean in normalized [0,1] space
    T span;      // max distance from pivot to either end

    void configure(T lo, T hi, T _mean)
    {
        minv = lo;
        maxv = hi;
        mean = _mean;

        T range = hi - lo;
        if (range != 0.0f)
            invRange = 1.0f / range;
        else
            invRange = 0.0f;

        // pivot in normalized [0,1] space
        pivot = (mean - lo) * invRange;
        if (pivot < 0.0f) pivot = 0.0f;
        if (pivot > 1.0f) pivot = 1.0f;

        // maximum distance to an edge
        span = std::max(pivot, 1.0f - pivot);
        if (span <= 0.0f)
            span = 1.0f; // avoid divide by zero
    }
};

struct PivotToneParams
{
    float brightness = 0.0f;
    float contrast = 1.0f;
    float gamma = 1.0f;

    bool operator==(const PivotToneParams&) const = default;
};

template<typename T>
inline GammaLUT<T> makeGammaLUT(float gamma)
{
    GammaLUT<T> lut;
    lut.invN = T{ 1 } / GammaLUT<T>::N;

    for (int i = 0; i <= GammaLUT<T>::N; ++i)
    {
        T x = i * lut.invN;
        lut.table[i] = bl::pow(x, gamma);
    }
    return lut;
}

template<typename T>
FORCE_INLINE T applyGammaLUT(const GammaLUT<T>& lut, T x)
{
    if (x < T{0}) x = T{0};
    if (x > T{1}) x = T{1};

    T f = x * GammaLUT<T>::N;
    int   i = (int)f;
    T frac = f - (T)i;

    T a = lut.table[i];
    T b = lut.table[i + 1];

    return a + (b - a) * frac;
}

template<typename T>
FORCE_INLINE T adjustFromPivot(const PivotToneParams& p,
    const ToneFieldInfo<T>& info,
    const GammaLUT<T>& lut,
    T x)
{
    // normalize input x into [0,1] based on min/max
    T t = (x - info.minv) * info.invRange;

    // normalized distance from pivot before contrast
    T d = t - info.pivot;

    // direction and magnitude
    T sign = (d >= T{ 0 }) ? T{ 1 } : T{ -1 };
    T mag = std::fabs(d);

    T span = info.span;
    T r = mag / span; // r = distance as a fraction of span

    // gamma shaping radius as a fixed fraction of span
    const T gammaFrac = T{ 1.0 };

    T r2 = r;
    if (r <= gammaFrac) {
        // Map gammaFrac -> [0,1] for LUT
        T x_norm = r / gammaFrac;
        if (x_norm < T{ 0 }) x_norm = T{ 0 };
        if (x_norm > T{ 1 }) x_norm = T{ 1 };

        T gn = applyGammaLUT(lut, x_norm);

        // Back to "fraction of span" space, but still within [0, gammaFrac]
        r2 = gn * gammaFrac;
    }
    // else: r2 == r  (no shaping outside the gamma band)

    // back to normalized distance units around pivot
    T d2 = sign * (r2 * span);

    // apply contrast after gamma, so gamma radius isn't blown up by contrast
    d2 *= p.contrast;

    // re-center around pivot and add brightness
    T t2 = info.pivot + d2;
    t2 += p.brightness;

    return t2;
}

// same as above, but returns denormalized value in original space
template<typename T>
FORCE_INLINE T adjustFromPivotDenormalized(const PivotToneParams& p,
    const ToneFieldInfo<T>& info,
    const GammaLUT<T>& lut,
    T x)
{
    T ret = adjustFromPivot(p, info, lut, x);
    ret = info.minv + ret * (info.maxv - info.minv);
    return ret;
}

template<typename T>
struct ToneManager
{
    GammaLUT<T> lut;
    ToneFieldInfo<T> info;
    PivotToneParams params;

    void setParams(PivotToneParams p) { params = p; }
    void setParams(float _brightness, float _contrast, float _gamma)
    {
        params.brightness = _brightness;
        params.contrast = _contrast;
        params.gamma = _gamma;
    }

    void configureFieldInfo(T lo, T hi, T mean)
    {
        info.configure(lo, hi, mean);
        lut = makeGammaLUT<T>(params.gamma);
    }

    T apply(T value)
    {
        return adjustFromPivot<T>(params, info, lut, value);
    }

    T applyAndDenormalize(T value)
    {
        return adjustFromPivotDenormalized<T>(params, info, lut, value);
    }
};

class GradientTexture
{
    mutable GLuint tex = 0;
    mutable bool tex_dirty = true;
    mutable int tex_width = 0;
    
    std::vector<uint32_t> colors;

    void ensureAllocated(int new_width) const
    {
        if (tex == 0)
            glGenTextures(1, &tex);

        if (tex_width == new_width)
            return;

        tex_width = new_width;

        glBindTexture(GL_TEXTURE_2D, tex);

        // linear makes the LUT behave like a smooth gradient
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

        // 2D texture of size [CACHE_SIZE x 1]
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, tex_width, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);

        glBindTexture(GL_TEXTURE_2D, 0);
    }

    void update(const uint32_t* colors_rgba_u32, int count) const
    {
        ensureAllocated(count);

        glBindTexture(GL_TEXTURE_2D, tex);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

        // upload packed RGBA8 data
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, count, 1, GL_RGBA, GL_UNSIGNED_BYTE, colors_rgba_u32);

        glBindTexture(GL_TEXTURE_2D, 0);
    }

public:

    void destroy()
    {
        if (tex)
            glDeleteTextures(1, &tex);

        tex = 0;
        tex_width = 0;
        tex_dirty = true;
    }

    int width() const
    {
        return (int)colors.size();
    }

    void set(const ImGradient& gradient)
    {
        colors.resize(gradient.cacheSize());
        std::memcpy(colors.data(), gradient.data(), gradient.cacheSize() * sizeof(uint32_t));
        tex_dirty = true;
    }

    GLuint getTexture() const
    {
        if (tex_dirty)
        {
            GLint prev_binding = 0;
            glGetIntegerv(GL_TEXTURE_BINDING_2D, &prev_binding);

            update(colors.data(), (int)colors.size());
            tex_dirty = false;

            glBindTexture(GL_TEXTURE_2D, (GLuint)prev_binding);
        }
        return tex;
    }
};

class MultiPassShader
{
    mutable std::vector<deferred_unique_ptr<ShaderSurface>> pass_surfaces;
    mutable GLuint black1x1Tex = 0;

    GLuint ensureBlack1x1Tex() const noexcept
    {
        if (black1x1Tex != 0)
            return black1x1Tex;

        const uint8_t rgba[4] = { 0, 0, 0, 255 };

        glGenTextures(1, &black1x1Tex);
        glBindTexture(GL_TEXTURE_2D, black1x1Tex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 1, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgba);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

        glBindTexture(GL_TEXTURE_2D, 0);

        return black1x1Tex;
    }

public:

    ~MultiPassShader()
    {
        GLuint blackTex = black1x1Tex;
        main_window()->threadQueue().post([blackTex]()
        {
            if (blackTex != 0)
                glDeleteTextures(1, &blackTex);
        });
    }

    void updateFragmentSource(std::string_view script)
    {
        struct PassMeta
        {
            bool enabled = true;

            bool has_scale = false;
            float scale = 1.0f;

            bool has_filter = false;
            std::string filter;

            std::string name;

            bool has_body_line = false;
            int first_body_line = 1;
        };

        struct PassBlock
        {
            PassMeta meta;
            std::string body;
        };

        static const char* base_frag_pre =
            "#version 300 es\n"
            "precision highp float;\n"
            "\n"
            "uniform vec2 u_outTexel;\n"
            "uniform sampler2D u_prev;\n"
            "uniform sampler2D u_features;\n"
            "uniform sampler2D u_gradient;\n"
            "uniform float u_time;\n"
            "\n"
            "in vec2 v_uv;\n"
            "out vec4 outColor;\n"
            "\n"
            "float clamp01(float x) { return clamp(x, 0.0, 1.0); }\n"
            "float wrap01(float x) { return fract(x); }\n"
            "vec4 samplePrev(vec2 uv) { return texture(u_prev, uv); }\n"
            "vec4 sampleGradient(float t)\n"
            "{\n"
            "    return texture(u_gradient, vec2(clamp(t, 0.0, 1.0), 0.5));\n"
            "}\n";

        static const char* base_frag_post =
            "\nvoid main() {\n"
            "    vec2 uv = v_uv;\n"
            "    vec4 f = texture(u_features, uv);\n"
            "\n"
            "    float iter = f.r;\n"
            "    float dist = f.g;\n"
            "    float stripe = f.b;\n"
            "\n"
            "    if (iter < -1.0) { outColor = vec4(1.0, 0.0, 1.0, 1.0); return; }\n"
            "    if (iter < 0.0)  { outColor = vec4(0.0, 0.0, 0.0, 1.0); return; }\n"
            "\n"
            "    vec4 prev = texture(u_prev, uv);\n"
            "    outColor = userShade(prev, iter, dist, stripe, uv);\n"
            "}\n";

        auto trim = [](std::string_view sv) -> std::string_view
        {
            while (!sv.empty() && (sv.front() == ' ' || sv.front() == '\t' || sv.front() == '\r'))
                sv.remove_prefix(1);
            while (!sv.empty() && (sv.back() == ' ' || sv.back() == '\t' || sv.back() == '\r'))
                sv.remove_suffix(1);
            return sv;
        };

        auto startsWith = [](std::string_view s, std::string_view p) -> bool
        {
            return s.size() >= p.size() && s.substr(0, p.size()) == p;
        };

        auto appendLine = [](std::string& dst, std::string_view line)
        {
            dst.append(line.data(), line.size());
            dst.push_back('\n');
        };

        auto parseDirective = [&](std::string_view line, std::string_view& out_key, std::string_view& out_arg) -> bool
        {
            line = trim(line);
            if (!startsWith(line, "//"))
                return false;

            line.remove_prefix(2);
            line = trim(line);

            if (!startsWith(line, "@"))
                return false;

            line.remove_prefix(1);
            line = trim(line);

            const size_t sp = line.find_first_of(" \t");
            if (sp == std::string_view::npos)
            {
                out_key = line;
                out_arg = {};
                return true;
            }

            out_key = trim(line.substr(0, sp));
            out_arg = trim(line.substr(sp + 1));
            return true;
        };

        auto applyPassMeta = [&](PassMeta& meta, std::string_view key, std::string_view arg)
        {
            if (key == "scale")
            {
                try {
                    meta.scale = std::stof(std::string(arg));
                    meta.has_scale = true;
                }
                catch (...) {}
            }
            else if (key == "filter")
            {
                meta.filter = std::string(arg);
                meta.has_filter = true;
            }
            else if (key == "enabled")
            {
                arg = trim(arg);
                meta.enabled = !(arg == "0" || arg == "false" || arg == "False");
            }
        };

        std::string global_code;
        std::vector<PassBlock> passes;

        PassBlock* cur_pass = nullptr;
        bool in_global_block = false;

        // ---- explicit line scan (while loop)
        std::string_view s = script;
        int line_no = 1;
        size_t parsed_prefix_len = 0;

        while (!s.empty())
        {
            parsed_prefix_len = script.size() - s.size();

            size_t eol = s.find('\n');
            std::string_view line = (eol == std::string_view::npos) ? s : s.substr(0, eol);

            std::string_view key, arg;
            const bool is_directive = parseDirective(line, key, arg);

            if (is_directive && key == "end")
            {
                // stop parsing the entire script here (do not include this line)
                // TODO: not working?
                break;
            }

            if (is_directive && key == "global")
            {
                in_global_block = true;
                cur_pass = nullptr;
            }
            else if (is_directive && key == "pass")
            {
                in_global_block = false;

                PassBlock pb;
                pb.meta.name = std::string(arg);
                passes.push_back(std::move(pb));
                cur_pass = &passes.back();
            }
            else if (is_directive && cur_pass && !cur_pass->meta.has_body_line)
            {
                // metadata lines only apply before any body content begins
                applyPassMeta(cur_pass->meta, key, arg);
            }
            else
            {
                if (cur_pass)
                {
                    if (!cur_pass->meta.has_body_line)
                    {
                        cur_pass->meta.has_body_line = true;
                        cur_pass->meta.first_body_line = line_no;
                    }
                    appendLine(cur_pass->body, line);
                }
                else
                {
                    // anything before the first pass, or inside @global, is shared
                    if (in_global_block || passes.empty())
                        appendLine(global_code, line);
                }
            }

            if (eol == std::string_view::npos)
            {
                parsed_prefix_len = script.size();
                break;
            }

            s.remove_prefix(eol + 1);
            line_no++;
        }

        // fallback: no @pass found -> treat entire script as a single pass body
        if (passes.empty())
        {
            const std::string_view parsed_script = script.substr(0, parsed_prefix_len);

            PassBlock pb;
            pb.meta.name = "pass0";
            pb.meta.has_body_line = true;
            pb.meta.first_body_line = 1;
            pb.body.assign(parsed_script.begin(), parsed_script.end());
            if (!pb.body.empty() && pb.body.back() != '\n')
                pb.body.push_back('\n');
            passes.push_back(std::move(pb));
        }

        // shrink
        while (pass_surfaces.size() > passes.size())
            pass_surfaces.pop_back();

        // grow (and provide deferred destruction)
        pass_surfaces.reserve(passes.size());
        while (pass_surfaces.size() < passes.size())
            pass_surfaces.emplace_back(make_deferred_unique<ShaderSurface>(main_window()->threadQueue()));

        for (int i = 0; i < (int)passes.size(); ++i)
        {
            const PassBlock& pb = passes[i];

            std::string frag;
            frag.reserve(4096 + global_code.size() + pb.body.size());

            frag += base_frag_pre;

            if (!global_code.empty())
            {
                frag += "\n#line 1\n";
                frag += global_code;
            }

            frag += "\nvec4 userShade(vec4 prev, float iter, float dist, float stripe, vec2 uv)\n{\n";

            frag += "#line ";
            frag += std::to_string(pb.meta.has_body_line ? pb.meta.first_body_line : 1);
            frag += "\n";

            if (!pb.meta.enabled)
                frag += "    return prev;\n";
            else
                frag += pb.body;

            frag += "}\n";
            frag += base_frag_post;

            pass_surfaces[i]->setFragmentSource(frag);

            if (pb.meta.has_scale)
                pass_surfaces[i]->setScale(pb.meta.scale);

            if (pb.meta.has_filter)
                pass_surfaces[i]->setFilter(pb.meta.filter.c_str());
        }
    }

    template<typename BindFn>
    void render(int w, int h, BindFn&& bind_inputs_and_uniforms) const
    {
        GLuint prevTex = ensureBlack1x1Tex();

        int pass_count = (int)pass_surfaces.size();
        for (int pass = 0; pass < pass_count; pass++)
        {
            auto s = pass_surfaces[pass].get();
            s->render(w, h, [&](const ShaderSurface& surf)
            {
                surf.bindTexture2D("u_prev", prevTex);
                std::forward<BindFn>(bind_inputs_and_uniforms)(surf);
            });
            prevTex = s->texture();
        }
    }

    const ShaderSurface* outputSurface() const
    {
        if (pass_surfaces.empty())
            return nullptr;
        return pass_surfaces.back().get();
    }
};

SIM_END;
