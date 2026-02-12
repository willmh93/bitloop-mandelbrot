#pragma once
#include <bitloop/nanovgx/nano_shader_surface.h>
#include "../core/types.h"

SIM_BEG;

class MultiPassShader
{
    mutable std::vector<ShaderSurfacePtr> pass_surfaces;
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

            "uniform vec2 u_outTexel;\n"
            "uniform sampler2D u_prev;\n"
            "uniform sampler2D u_features;\n"
            "uniform sampler2D u_gradient;\n"
            "uniform float min_iter;\n"
            "uniform float max_iter;\n"
            "uniform float min_dist;\n"
            "uniform float max_dist;\n"
            "uniform float min_stripe;\n"
            "uniform float max_stripe;\n"
            ///"uniform float u_time;\n" // todo: requires checking presence to determine if we should reshade every frame

            "in vec2 v_uv;\n"
            "out vec4 outColor;\n"

            "float clamp01(float x) { return clamp(x, 0.0, 1.0); }\n"
            "float wrap01(float x) { return fract(x); }\n"
            "vec4 samplePrev(vec2 uv) { return texture(u_prev, uv); }\n"
            "vec4 sampleGradient(float t)\n"
            "{\n"
            "    return texture(u_gradient, vec2(clamp(t, 0.0, 1.0), 0.5));\n"
            "}\n"

            "float _bl_dist;\n"
            "float _bl_inv_dist_range;\n"
            "\n"
            "float boundary(float thresholdA, float thresholdB)\n"
            "{\n"
            "    float x = (_bl_dist - min_dist) * _bl_inv_dist_range;\n"
            "    return 1.0 - smoothstep(thresholdA, thresholdB, x);\n"
            "}\n"
            "\n"
            "float boundary(float threshold)\n"
            "{\n"
            "    return boundary(threshold, threshold+0.1);\n"
            "}\n"
            ;

        static const char* base_frag_post =
            "\nvoid main() {\n"
            "    vec2 uv = v_uv;\n"
            "    vec4 f = texture(u_features, uv);\n"

            "    float iter = f.r;\n"
            "    float dist = f.g;\n"
            "    float stripe = f.b;\n"

            "    if (iter < -1.0) { outColor = vec4(1.0, 0.0, 1.0, 1.0); return; }\n"
            "    if (iter < 0.0)  { outColor = vec4(0.0, 0.0, 0.0, 1.0); return; }\n"

            "    vec4 prev = texture(u_prev, uv);\n"

            "    _bl_dist = dist;\n"
            "    _bl_inv_dist_range = 1.0 / (max_dist - min_dist);\n"

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

        // resize
        pass_surfaces.resize(passes.size());
        for (auto& s : pass_surfaces)
            if (!s) s = makeShaderSurface(main_window()->threadQueue());


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

    bool testCompile(
        bool& combined_vert_errors,
        bool& combined_frag_errors,
        std::string& combined_vert_log,
        std::string& combined_frag_log)
    {
        int pass_count = (int)pass_surfaces.size();

        combined_vert_errors = false;
        combined_frag_errors = false;
        combined_vert_log = "";
        combined_frag_log = "";

        for (int pass = 0; pass < pass_count; pass++)
        {
            auto s = pass_surfaces[pass].get();

            bool pass_vert_errors = false;
            bool pass_frag_errors = false;
            std::string pass_vert_log;
            std::string pass_frag_log;

            s->testCompile(
                pass_vert_errors,
                pass_frag_errors,
                pass_vert_log,
                pass_frag_log
            );

            if (pass_vert_errors) { combined_vert_errors = true; combined_vert_log += pass_vert_log + "\n"; }
            if (pass_frag_errors) { combined_frag_errors = true; combined_frag_log += pass_frag_log + "\n"; }
        }

        return !(combined_vert_errors || combined_frag_errors);
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