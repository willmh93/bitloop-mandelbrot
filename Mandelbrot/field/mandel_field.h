#pragma once
#include <bitloop.h>

#include "types.h"

SIM_BEG;

using namespace bl;

FORCE_INLINE void mul_angle(f32 c, f32 s, uint32_t n, f32& cn, f32& sn)
{
    f32 rx = 1.0f, ry = 0.0f;
    f32 bx = c, by = s;

    while (n)
    {
        if (n & 1u)
        {
            #ifdef FMA_AVAILABLE
            const f32 t0 = std::fma(rx, bx, -ry * by);
            const f32 t1 = std::fma(rx, by, ry * bx);
            #else
            const f32 t0 = rx * bx - ry * by;
            const f32 t1 = rx * by + ry * bx;
            #endif
            rx = t0;
            ry = t1;
        }

        n >>= 1;
        if (!n)
            break;

        const f32 t0 = bx * bx - by * by;
        const f32 t1 = (bx + bx) * by;
        bx = t0;
        by = t1;
    }

    cn = rx;
    sn = ry;
}

struct IterParams
{
    bool    cycle_iter_dynamic_limit = false;
    bool    cycle_iter_normalize_depth = false;
    double  cycle_iter_log1p_weight = 1.0;
    double  cycle_iter_value = 1.0f; // If dynamic, iter_lim ratio, else iter_lim

    bool operator==(const IterParams&) const = default;
};

struct DistParams
{
    bool    cycle_dist_invert = false;
    double  cycle_dist_value = 0.25;
    double  cycle_dist_sharpness = 0.9; // Used for UI (ignored during tween)

    bool operator==(const DistParams&) const = default;
};

struct StripeParams
{
    int freq = 3; // stripes per 2 pi
    float phase = 0.0; // radians

    StripeParams(int _freq = 8, float _phase = 0.0f)
        : freq(_freq), phase(_phase)
    {}

    bool operator==(const StripeParams&) const = default;
};

struct StripeAccum
{
    int  n = 1; // frequency used (needed?)
    int  W = 0; // number of samples
    f32 sum_sin = 0.0f;
    f32 sum_cos = 0.0f;
    f32 last_sin = 0.0f;
    f32 last_cos = 0.0f;

    f64 log_r2_at_escape = 0.0;
    bool escaped = false;

    StripeAccum() : n(1) {}
    StripeAccum(int freq) : n(freq) {}

    void accumulate(f32 c, f32 s)
    {
        f32 cn, sn;
        mul_angle(c, s, n, cn, sn);

        sum_cos += cn;
        sum_sin += sn;
        last_cos = cn;
        last_sin = sn;
        W++;
    }

    void escape(f64 r2)
    {
        log_r2_at_escape = std::log(r2);
        escaped = true;
    }
};

inline float stripeFromAccum(const StripeAccum& st, double log_er2, float cphi, float sphi)
{
    if (st.W <= 0) return 0.0f;

    float invW = 1.0f / (float)st.W;
    float avg = 0.5f + 0.5f * ((st.sum_sin * invW) * cphi + (st.sum_cos * invW) * sphi);

    // prev: exclude last sample to reproduce your mix
    float prev = avg;
    if (st.W > 1) {
        float invWm1 = 1.0f / (float)(st.W - 1);
        float S2 = (st.sum_sin - st.last_sin) * invWm1;
        float C2 = (st.sum_cos - st.last_cos) * invWm1;
        prev = 0.5f + 0.5f * (S2 * cphi + C2 * sphi);
    }

    // escape fraction: frac = 1 + log2( log(ER^2) / log(r2_escape) )
    float frac = 0.0f;
    if (st.escaped) {
        const double lr = std::max(st.log_r2_at_escape, 1e-300);
        double val = 1.0 + std::log2(log_er2 / lr);
        if (val < 0.0) val = 0.0; else if (val > 1.0) val = 1.0;
        frac = (float)val;
    }

    float mix = frac * avg + (1.0f - frac) * prev;

    // contrast shaping
    mix = 0.5f + 0.5f * std::tanh((mix - 0.5f));
    if (mix < 0.0f) mix = 0.0f; else if (mix > 1.0f) mix = 1.0f;
    return mix;
}

inline float stripeFromAccum_Newton(const StripeAccum& st, float cphi, float sphi)
{
    if (st.W <= 0) return 0.0f;

    float invW = 1.0f / (float)st.W;
    float avg = 0.5f + 0.5f * ((st.sum_sin * invW) * cphi +
        (st.sum_cos * invW) * sphi);

    // gentle shaping instead of strong tanh:
    float mix = avg;
    // optional mild contrast:
    // mix = 0.5f + 0.5f * std::tanh(params.contrast * (mix - 0.5f));

    return std::clamp(mix, 0.0f, 1.0f);
}


#if MANDEL_EXTENDED_FIELD_STATS
struct ExtendedFieldStats
{
    size_t slack_rebases     = 0;
    size_t dominance_rebases = 0;
    size_t iter_count        = 0;
    size_t escape_count      = 0;

    void reset()
    {
        slack_rebases        = 0;
        dominance_rebases    = 0;
        iter_count           = 0;
        escape_count         = 0;
    }

    void operator +=(const ExtendedFieldStats& rhs)
    {
        slack_rebases       += rhs.slack_rebases;
        dominance_rebases   += rhs.dominance_rebases;
        iter_count          += rhs.iter_count;
        escape_count        += rhs.escape_count;
    }
};
#endif

struct EscapeFieldPixel
{
    f64 depth;
    f64 dist;
    StripeAccum stripe;

    f32 stripe_min;
    f32 stripe_max;

    f32 final_depth;
    f32 final_dist;
    f32 final_stripe;

    #if MANDEL_EXTENDED_FIELD_STATS
    ExtendedFieldStats extended_stats;
    #endif

    void set(const EscapeFieldPixel& p)
    {
        depth = p.depth;
        dist = p.dist;
        stripe = p.stripe;
    }

    void set(f64 _depth, f64 _dist, const StripeAccum& _stripe) {
        depth  = _depth;
        dist   = _dist;
        stripe = _stripe;
    }
};

constexpr double phaseBmpScale(int phase)
{
    if (phase == 0) return 9.0;
    if (phase == 1) return 3.0;
    return 1.0;
}

struct EscapeField : public std::vector<EscapeFieldPixel>
{
    int compute_phase;

    f64 min_depth = 0;
    f64 max_depth = 0;
    f64 mean_depth = 0;

    f64 stable_min_dist = 0;
    f64 stable_max_dist = 0;

    f32 min_stripe = 0;
    f32 max_stripe = 0;
    f32 mean_stripe = 0;

    f64 log_color_cycle_iters{};
    f64 cycle_dist_value{};

    std::vector<i8> skip_flags;

    int w = 0, h = 0;

    std::vector<float> gpu_data;
    GLuint features_tex = 0;
    int tex_w = 0, tex_h = 0;

    EscapeField(int phase) : compute_phase(phase) {}

    void reset()
    {
        memset(data(), 0x80, size() * sizeof(EscapeFieldPixel)); // pattern gives each float type a negative value
        memset(skip_flags.data(), 0, skip_flags.size() * sizeof(i8));

        #if MANDEL_EXTENDED_FIELD_STATS
        for (size_t i = 0; i < size(); i++)
        {
            EscapeFieldPixel& field_pixel = std::vector<EscapeFieldPixel>::at(i);
            field_pixel.extended_stats.reset();
        }
        #endif
    }
    void setDimensions(int _w, int _h)
    {
        w = _w;
        h = _h;
        resize(w * h, { -1.0, -1.0 });
        skip_flags.resize(w * h, 0);
    }

    EscapeFieldPixel& operator ()(int x, int y)
    {
        return std::vector<EscapeFieldPixel>::at(y * w + x);
    }
    EscapeFieldPixel& at(int x, int y)
    {
        return std::vector<EscapeFieldPixel>::at(y * w + x);
    }
    EscapeFieldPixel* get(int x, int y)
    {
        int i = y * w + x;
        if (i < 0 || i >= size()) return nullptr;
        return data() + i;
    }

    void set_skip_flag(int x, int y, u8 b)
    {
        skip_flags[y * w + x] = b;
    }
    u8   get_skip_flag(int x, int y)
    {
        return skip_flags[y * w + x];
    }
    bool safe(int x, int y)
    {
        if (x < 0) return false;
        if (y < 0) return false;
        if (x >= w) return false;
        if (y >= h) return false;
        return true;
    }
    bool has_valid_depth(int x, int y)
    {
        if (safe(x, y))
        {
            EscapeFieldPixel& pixel = at(x, y);
            if (pixel.depth < INSIDE_MANDELBROT_SET_SKIPPED) return true;
        }
        return false;
    }

    void contractSkipFlags(int r = 1)
    {
        if (w <= 0 || h <= 0 || r <= 0) return;

        const int W = w;
        const int H = h;
        const size_t N = size_t(W) * size_t(H);

        // Binary mask from skip_flags
        std::vector<uint8_t> mask(N);
        for (size_t i = 0; i < N; ++i)
            mask[i] = skip_flags[i] ? 1u : 0u;

        // Summed-area table
        const int PSW = W + 1;
        const int PSH = H + 1;
        std::vector<uint32_t> ps(size_t(PSW) * PSH, 0u);

        for (int y = 0; y < H; ++y)
        {
            uint32_t rowSum = 0;
            const size_t rowOffset = size_t(y) * W;
            const size_t psRow = size_t(y + 1) * PSW;
            const size_t psPrevRow = size_t(y) * PSW;

            for (int x = 0; x < W; ++x)
            {
                rowSum += mask[rowOffset + x];
                ps[psRow + (x + 1)] = ps[psPrevRow + (x + 1)] + rowSum;
            }
        }

        const int windowSize = (2 * r + 1) * (2 * r + 1);

        for (int y = 0; y < H; ++y)
        {
            for (int x = 0; x < W; ++x)
            {
                const size_t idx = size_t(y) * W + x;

                // If original was 0, erosion must also be 0.
                if (!mask[idx])
                {
                    skip_flags[idx] = 0;
                    continue;
                }

                if (x < r || x >= W - r || y < r || y >= H - r)
                {
                    skip_flags[idx] = 0;
                    continue;
                }

                const int x0 = x - r;
                const int y0 = y - r;
                const int x1 = x + r;
                const int y1 = y + r;

                // Rectangle sum
                const uint32_t sum =
                    ps[(y1 + 1) * PSW + (x1 + 1)] -
                    ps[(y0)*PSW + (x1 + 1)] -
                    ps[(y1 + 1) * PSW + (x0)] +
                    ps[(y0)*PSW + (x0)];

                // All neighbors 1? If not, clear.
                skip_flags[idx] = (sum == uint32_t(windowSize)) ? 1 : 0;
            }
        }
    }
    void expandSkipFlags(int r = 1, bool overwrite = false)
    {
        if (w <= 0 || h <= 0) return;

        const int W = w;
        const int H = h;
        const size_t N = size_t(W) * size_t(H);

        // r <= 0 is a corner case
        if (r <= 0)
        {
            if (overwrite)
                std::fill(skip_flags.begin(), skip_flags.end(), 0);
            return;
        }

        // Read current mask from skip_flags
        std::vector<uint8_t> mask(N);
        for (size_t i = 0; i < N; ++i)
            mask[i] = skip_flags[i] ? 1u : 0u;

        // Precompute disk offsets (dx,dy) with dx*dx + dy*dy <= r*r
        struct Offset { int dx, dy; };
        static thread_local std::vector<Offset> s_offsets;
        static thread_local int s_cached_r = -1;

        if (s_cached_r != r)
        {
            s_cached_r = r;
            s_offsets.clear();
            const int r2 = r * r;
            for (int dy = -r; dy <= r; ++dy)
            {
                for (int dx = -r; dx <= r; ++dx)
                {
                    if (dx * dx + dy * dy <= r2)
                        s_offsets.push_back({ dx, dy });
                }
            }
        }

        // ring will hold the dilated mask first, then we subtract the original mask
        std::vector<uint8_t> ring(N, 0u);

        // Dilation
        for (int y = 0; y < H; ++y)
        {
            const size_t rowOffset = size_t(y) * W;
            for (int x = 0; x < W; ++x)
            {
                if (!mask[rowOffset + x]) continue;

                for (const auto& o : s_offsets)
                {
                    const int xx = x + o.dx;
                    const int yy = y + o.dy;
                    if ((unsigned)xx >= (unsigned)W || (unsigned)yy >= (unsigned)H)
                        continue;

                    ring[size_t(yy) * W + xx] = 1u;
                }
            }
        }

        // Outline = dilated AND NOT mask
        for (size_t i = 0; i < N; ++i)
            ring[i] = (ring[i] & (mask[i] ^ 1u));

        // Write back
        if (overwrite)
        {
            for (size_t i = 0; i < N; ++i)
                skip_flags[i] = ring[i] ? 1 : 0;
        }
        else
        {
            for (size_t i = 0; i < N; ++i)
            {
                if (ring[i])
                    skip_flags[i] = 1;
            }
        }
    }


    void fillFeaturesTexureData()
    {
        const size_t len = size();
        gpu_data.resize(len * 4);

        for (size_t i = 0; i < len; i++)
        {
            const EscapeFieldPixel& px = std::vector<EscapeFieldPixel>::at(i);

            const float iter = px.final_depth;
            const float dist = px.final_dist;
            const float stripe = px.final_stripe;

            gpu_data[i * 4 + 0] = iter;
            gpu_data[i * 4 + 1] = dist;
            gpu_data[i * 4 + 2] = stripe;
            gpu_data[i * 4 + 3] = 0.0f;
        }
    }


    void initFeaturesTexture()
    {
        if (features_tex == 0)
            glGenTextures(1, &features_tex);

        if (tex_w == w && tex_h == h)
            return;

        tex_w = w;
        tex_h = h;

        glBindTexture(GL_TEXTURE_2D, features_tex);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

        // reallocate storage for new resolution
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, w, h, 0, GL_RGBA, GL_FLOAT, nullptr);

        glBindTexture(GL_TEXTURE_2D, 0);
    }


    void updateFeaturesTexture()
    {
        initFeaturesTexture();

        glBindTexture(GL_TEXTURE_2D, features_tex);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

        // upload new feature data into already-allocated float texture
        glTexSubImage2D(
            GL_TEXTURE_2D,
            0,
            0, 0,
            w, h,
            GL_RGBA,
            GL_FLOAT,
            gpu_data.data()
        );

        glBindTexture(GL_TEXTURE_2D, 0);
    }

    void destroyFeaturesTexture()
    {
        if (features_tex != 0)
        {
            glDeleteTextures(1, &features_tex);
            features_tex = 0;
            tex_w = 0;
            tex_h = 0;
        }
    }

private:

    inline uint8_t dilate_at(const std::vector<uint8_t>& in, int w, int h, int x, int y, int r)
    {
        const int r2 = r * r;
        for (int dy = -r; dy <= r; ++dy)
        {
            int yy = y + dy;
            if (yy < 0 || yy >= h) continue;
            for (int dx = -r; dx <= r; ++dx)
            {
                if (dx * dx + dy * dy > r2) continue; // outside the disk
                int xx = x + dx;
                if (xx < 0 || xx >= w) continue;
                if (in[yy * w + xx]) return 1;
            }
        }
        return 0;
    }
};

struct NormalizationPixel : public EscapeFieldPixel
{
private:
    union {
        FVec2  world_pos_32;
        DVec2  world_pos_64;
        DDVec2 world_pos_128;
    };
public:

    NormalizationPixel() noexcept
        : EscapeFieldPixel{}
        , world_pos_128{}
        , stage_pos{}
        , is_final(false)
        , weight(0.0f)
    {
    }

    DVec2  stage_pos;
    bool   is_final;
    float  weight;

    template<typename T> requires is_f32<T>  [[nodiscard]] constexpr const FVec2& worldPos() { return world_pos_32; }
    template<typename T> requires is_f64<T>  [[nodiscard]] constexpr const DVec2& worldPos() { return world_pos_64; }
    template<typename T> requires is_f128<T> [[nodiscard]] constexpr const DDVec2& worldPos() { return world_pos_128; }

    void setWorldPos(FVec2 p)  { world_pos_32 = p; }
    void setWorldPos(DVec2 p)  { world_pos_64 = p; }
    void setWorldPos(DDVec2 p) { world_pos_128 = p; }
};

class NormalizationField
{
    std::vector<DVec2> local_field;
    CanvasObject128 bounds;

public:

    std::vector<NormalizationPixel> world_field;

    NormalizationField()
    {
        setShape(0.025, 2.0);
    }

    void setShape(double sample_r, double exponent)
    {
        local_field = math::delaunayMeshEllipse<f64>(0, 0, 1.0, sample_r, exponent);
        world_field.resize(local_field.size());
    }

    template<typename T>
    void updateField(const CameraInfo &camera, double scale)
    {
        T world_radius = ((T)scale / camera.relativeZoom<T>()) * T{ 2 };
        Vec2<T> cam_center_world = camera.pos<T>();

        // Get stage quad of world ellipse bounds for fast stage interpolation
        CanvasObjectBase<T> bounds;
        bounds.setCamera(camera);
        bounds.setWorldRect(cam_center_world - world_radius, world_radius * 2);
        DQuad stage_quad = bounds.stageQuad();

        for (size_t i=0; i< world_field.size(); i++)
        {
            NormalizationPixel& px = world_field[i];
            const DVec2 pt = local_field[i];

            //px.world_pos = cam_center_world + (pt * world_radius);
            px.setWorldPos(cam_center_world + (pt * world_radius));
            px.stage_pos = stage_quad.lerpPoint(0.5 + pt * 0.5);

            px.weight = 1.0f - (float)pt.mag();
        }
    }

    void clearFinalFlags()
    {
        for (NormalizationPixel& p : world_field)
            p.is_final = false;
    }

    template<typename Callback>
    void forEach(Callback&& callback, int thread_count = Thread::threadCount())
    {
        //static_assert(std::is_invocable_r_v<void, Callback, NormalizationPixel&>,
        //    "Callback must be: void(NormalizationPixel&)");

        int pixel_count = (int)world_field.size();
        std::vector<std::pair<int, int>> ranges = Thread::splitRanges<int>(pixel_count, thread_count);

        std::vector<std::future<void>> futures(thread_count);

        for (int ti = 0; ti < thread_count; ti++)
        {
            const std::pair<int, int>& range = ranges[ti];
            futures[ti] = Thread::pool().submit_task([&]()
            {
                int i0 = range.first;
                int i1 = range.second;
                if constexpr (std::is_invocable_r_v<void, Callback, NormalizationPixel&>)
                {
                    for (int i = i0; i < i1; i++)
                        callback(world_field[i]);
                }
                else
                {
                    for (int i = i0; i < i1; i++)
                        callback(world_field[i], ti);
                }
            });
        }

        for (int ti = 0; ti < thread_count; ti++)
        {
            if (futures[ti].valid())
                futures[ti].get();
        }
    }
};

SIM_END;
