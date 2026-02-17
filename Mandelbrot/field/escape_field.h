#pragma once
#include <bitloop.h>

#include "../core/types.h"

SIM_BEG;

using namespace bl;

// for accumulating stripe cos/sin values
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
    bool    iter_dynamic_limit = false;
    bool    iter_normalize_depth = false;
    double  iter_log1p_weight = 1.0;
    double  iter_cycle_value = 2.0f; // If dynamic, iter_lim ratio, else iter_lim
};

struct DistParams
{
    bool    dist_invert = false;
    double  dist_cycle_value = 0.25;
    double  cycle_dist_sharpness = 0.9; // Used for UI (ignored during tween)
};

struct StripeParams
{
    int freq = 3; // stripes per 2 pi
    float phase = 0.0; // radians

    StripeParams(int _freq = 8, float _phase = 0.0f)
        : freq(_freq), phase(_phase)
    {}
};

struct StripeAccum
{
    int  n = 1; // frequency used
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

    void escape(f32 r2)
    {
        log_r2_at_escape = std::log(r2);
        escaped = true;
    }

    void escape(f64 r2)
    {
        log_r2_at_escape = std::log(r2);
        escaped = true;
    }

    void no_escape()
    {
        log_r2_at_escape = 0.0;
    }
};

inline float stripeFromAccum(const StripeAccum& st, double log_er2, float cphi, float sphi)
{
    if (st.W <= 0) return 0.0f;

    float invW = 1.0f / (float)st.W;
    float avg = 0.5f + 0.5f * ((st.sum_sin * invW) * cphi + (st.sum_cos * invW) * sphi);

    // exclude last sample to reproduce mix
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

    f32 final_depth;
    f32 final_dist;
    f32 final_stripe;

    #if MANDEL_EXTENDED_FIELD_STATS
    ExtendedFieldStats extended_stats;
    #endif

    void set(const EscapeFieldPixel& p)
    {
        depth  = p.depth;
        dist   = p.dist;
        stripe = p.stripe;
    }

    void set(f64 _depth, f64 _dist, const StripeAccum& _stripe) {
        depth  = _depth;
        dist   = _dist;
        stripe = _stripe;
    }
};

struct EscapeField : public std::vector<EscapeFieldPixel>
{
    int w = 0, h = 0;
    int phase = -1;

    // bools indicating which pixels to skip (interior)
    std::vector<i8> skip_flags;

    // raw min/max depth (computed directly from field)
    f64 raw_min_depth   = 0;
    f64 raw_max_depth   = 0;

    // raw min/max dist is derived from formula based on zoom level
    f64 raw_min_dist    = 0;
    f64 raw_max_dist    = 0;

    // log min/max dist
    f64 log_min_dist = 0;
    f64 log_max_dist = 0;

    // raw min/max stripe computed from actual mean 
    // + magnitude derived from spline (approximate from histogram)
    f32 raw_min_stripe  = 0;
    f32 raw_max_stripe  = 0;
    f32 raw_mean_stripe = 0;
    f32 raw_mag_stripe = 0;

    // normalized iter cycle size (smaller repeats gradient at higher freq)
    f64 log_color_cycle_iters = 0;

    // final min/max values after normalization/toning/cycling
    f32 final_min_depth  = 0;
    f32 final_max_depth  = 0;
    f32 final_min_dist   = 0;
    f32 final_max_dist   = 0;
    f32 final_min_stripe = 0;
    f32 final_max_stripe = 0;

    // data for GPU upload
    std::vector<float> gpu_data;
    GLuint features_tex = 0;
    int tex_w = 0, tex_h = 0;

    EscapeField(int phase) : phase(phase) {}

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
        const size_t new_size = (size_t)_w * (size_t)_h;
        const bool changed = (_w != w) || (_h != h);

        w = _w;
        h = _h;

        resize(new_size);
        skip_flags.resize(new_size, 0);

        if (changed)
            reset(); // forces depth negative sentinel again for the new layout
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

    bool featuresTextureDirty = false;

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

        featuresTextureDirty = true;
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

        GLint prev_unpack = 4;
        glGetIntegerv(GL_UNPACK_ALIGNMENT, &prev_unpack);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, w, h, 0, GL_RGBA, GL_FLOAT, nullptr);

        glPixelStorei(GL_UNPACK_ALIGNMENT, prev_unpack);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    void updateFeaturesTexture()
    {
        initFeaturesTexture();

        glBindTexture(GL_TEXTURE_2D, features_tex);

        GLint prev_unpack = 4;
        glGetIntegerv(GL_UNPACK_ALIGNMENT, &prev_unpack);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, GL_RGBA, GL_FLOAT, gpu_data.data());

        glPixelStorei(GL_UNPACK_ALIGNMENT, prev_unpack);
        glBindTexture(GL_TEXTURE_2D, 0);

        featuresTextureDirty = false;
    }


    GLuint getTexture() const
    {
        if (featuresTextureDirty)
            const_cast<EscapeField*>(this)->updateFeaturesTexture();

        return features_tex;
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


SIM_END;
