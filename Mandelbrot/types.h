#pragma once
#include <bitloop.h>

#include <cmath>
#include <vector>

#include "build_config.h"

SIM_BEG;

using namespace bl;

constexpr double INSIDE_MANDELBROT_SET = std::numeric_limits<double>::max();
const double INSIDE_MANDELBROT_SET_SKIPPED = std::nextafter(INSIDE_MANDELBROT_SET, 0.0f);

enum MandelFlag : uint32_t
{
    // bits
    MANDEL_DYNAMIC_ITERS        = 1u << 0,
    MANDEL_SHOW_AXIS            = 1u << 1,
    MANDEL_FLATTEN              = 1u << 2,
    MANDEL_DYNAMIC_COLOR_CYCLE  = 1u << 3,
    MANDEL_NORMALIZE_DEPTH      = 1u << 4,
    MANDEL_INVERT_DIST          = 1u << 5,
    MANDEL_USE_SMOOTHING        = 1u << 6,

    // bitmasks
    //MANDEL_FLAGS_MASK   = 0x000FFFFFu, // max 24 bit-flags
    //MANDEL_SMOOTH_MASK  = 0x00F00000u, // max 16 smooth types
    MANDEL_VERSION_MASK   = 0xFF000000u, // max 255 versions

    MANDEL_VERSION_BITSHIFT = 24
};

enum struct MandelMaxDepthOptimization
{
    SLOWEST,
    SLOW,
    MEDIUM,
    FAST,
    COUNT
};

static inline const char* MandelMaxDepthOptimizationNames[(int)MandelMaxDepthOptimization::COUNT] =
{
    "Slowest (Lossless)",
    "Slow",
    "Balanced",
    "Fast (Loss at sharp angles)"
};


static inline const char* FloatingPointTypeNames[3] =
{
    "32-bit",
    "64-bit",
    "128-bit (quad-precision)"
};

struct IterParams
{
    double  cycle_iter_weight = 0.0;
    bool    cycle_iter_dynamic_limit = true;
    bool    cycle_iter_normalize_depth = true;
    double  cycle_iter_normalize_low_fact = 25;
    double  cycle_iter_normalize_high_fact = 50;
    double  cycle_iter_log1p_weight = 0.0;
    double  cycle_iter_value = 0.5f; // If dynamic, iter_lim ratio, else iter_lim

    bool operator==(const IterParams&) const = default;
};

struct DistParams
{
    double  cycle_dist_weight = 0.0;
    bool    cycle_dist_invert = false;
    double  cycle_dist_value = 0.5f;
    double  cycle_dist_sharpness = 0.9; // Used for UI (ignored during tween)

    bool operator==(const DistParams&) const = default;
};

struct StripeParams
{
    float freq = 3.0; // stripes per 2π
    float phase = 0.0; // radians
    float contrast = 3.0; // tanh shaping

    StripeParams(float _freq = 8.0f, float _phase = 0.0f, float _contrast = 3.0f)
        : freq(_freq), phase(_phase), contrast(_contrast)
    {}

    bool operator==(const StripeParams&) const = default;
};


struct EscapeFieldPixel
{
    f64 depth;

    union
    {
        f32  dist_32;
        f64  dist_64;
        f128 dist_128;
    };

    union
    {
        f32  stripe_32;
        f64  stripe_64;
        f128 stripe_128;
    };

    // todo: Reuse above vars to save space
    f32 final_depth;
    f32 final_dist;
    f32 final_stripe;

    bool flag_for_skip;

    inline void setDist(f32 d32)   { dist_32 = d32; }
    inline void setDist(f64 d64)   { dist_64 = d64; }
    inline void setDist(f128 d128) { dist_128 = d128; }

    inline void setStripe(f32 s32)   { stripe_32 = s32; }
    inline void setStripe(f64 s64)   { stripe_64 = s64; }
    inline void setStripe(f128 s128) { stripe_128 = s128; }

    template<typename T> requires std::same_as<T, f32>  constexpr T getDist() { return dist_32; }
    template<typename T> requires std::same_as<T, f64>  constexpr T getDist() { return dist_64; }
    template<typename T> requires std::same_as<T, f128> constexpr T getDist() { return dist_128; }

    template<typename T> requires std::same_as<T, f32>  constexpr T getStripe() { return stripe_32; }
    template<typename T> requires std::same_as<T, f64>  constexpr T getStripe() { return stripe_64; }
    template<typename T> requires std::same_as<T, f128> constexpr T getStripe() { return stripe_128; }
};

struct EscapeField : public std::vector<EscapeFieldPixel>
{
    int compute_phase;

    double min_depth = 0.0;
    double max_depth = 0.0;

    double assumed_iter_min = 0.0;
    double assumed_iter_max = 0.0;

    f128 stable_min_dist{};
    f128 stable_max_dist{};

    f128 min_stripe{};
    f128 max_stripe{};

    double log_color_cycle_iters;
    double cycle_dist_value;

    int w = 0, h = 0;

    EscapeField(int phase) : compute_phase(phase) {}

    void setAllDepth(double value)
    {
        for (int i = 0; i < size(); i++)
        {
            EscapeFieldPixel &p = std::vector<EscapeFieldPixel>::at(i);
            p.depth = value;
            p.dist_128 = { value, 0 };
            p.stripe_128 = { value, 0 };
            p.flag_for_skip = false;
        }
    }
    void setDimensions(int _w, int _h)
    {
        //if (size() >= (_w * _h)) return;
        w = _w;
        h = _h;
        resize(w * h, { -1.0, -1.0 });
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

    bool safe(int x, int y)
    {
        if (x < 0) return false;
        if (y < 0) return false;
        if (x >= w) return false;
        if (y >= h) return false;
        return true;
    }

    bool has_data(int x, int y)
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

        auto idx = [this](int x, int y) { return y * w + x; };

        // Snapshot current flags
        std::vector<uint8_t> in(size_t(w) * size_t(h), 0);
        for (int y = 0; y < h; ++y)
            for (int x = 0; x < w; ++x)
                in[idx(x, y)] = at(x, y).flag_for_skip ? 1 : 0;

        // Erode
        std::vector<uint8_t> out(size_t(w) * size_t(h), 0);
        for (int y = 0; y < h; ++y)
        {
            for (int x = 0; x < w; ++x)
            {
                size_t i = idx(x, y);
                if (!in[i]) { out[i] = 0; continue; }

                bool inner = true;
                for (int dy = -r; dy <= r && inner; ++dy)
                    for (int dx = -r; dx <= r; ++dx)
                    {
                        if (dx == 0 && dy == 0) continue;
                        int nx = x + dx, ny = y + dy;
                        if (!safe(nx, ny) || !in[idx(nx, ny)]) { inner = false; break; }
                    }
                
                out[i] = inner ? 1 : 0;
            }
        }

        // Write back
        for (int y = 0; y < h; ++y)
            for (int x = 0; x < w; ++x)
                at(x, y).flag_for_skip = (out[idx(x, y)] != 0);
    }

    void expandSkipFlags(int r = 1, bool overwrite = false)
    {
        if (w <= 0 || h <= 0) return;
        const int W = w, H = h;
        const size_t N = size_t(W) * size_t(H);

        // read current mask
        std::vector<uint8_t> mask(N, 0);
        for (int y = 0; y < H; ++y)
            for (int x = 0; x < W; ++x)
                mask[y * W + x] = at(x, y).flag_for_skip ? 1 : 0;

        // dilate
        std::vector<uint8_t> dil(N, 0);
        for (int y = 0; y < H; ++y)
            for (int x = 0; x < W; ++x)
                dil[y * W + x] = dilate_at(mask, W, H, x, y, r);

        // outline = dilated AND NOT mask
        std::vector<uint8_t> ring(N, 0);
        for (size_t i = 0; i < N; ++i)
            ring[i] = (dil[i] & (mask[i] ^ 1)); // 1 where new perimeter appears

        // write back
        for (int y = 0; y < H; ++y)
        {
            for (int x = 0; x < W; ++x) {
                bool v = ring[y * W + x];
                if (overwrite)
                    at(x, y).flag_for_skip = v;
                else
                    at(x, y).flag_for_skip = at(x, y).flag_for_skip || v;
            }
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
                if (dx * dx + dy * dy > r2) continue;           // outside the disk
                int xx = x + dx;
                if (xx < 0 || xx >= w) continue;
                if (in[yy * w + xx]) return 1;
            }
        }
        return 0;
    }
};

SIM_END;
