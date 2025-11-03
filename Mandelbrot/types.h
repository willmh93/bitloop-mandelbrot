#pragma once
#include <bitloop.h>

#include <cmath>
#include <vector>

enum class MandelKernelFeatures
{
    NONE = 0,
    ITER = 1,
    DIST = 2,
    STRIPES = 4,
    MIX_ALL = 7,
    COUNT
};
// todo: find way to put inside SIM_BEG ns (wasm32 error, must be in global)
bl_enable_enum_bitops(MandelKernelFeatures);

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
    bool    cycle_iter_dynamic_limit = false;
    bool    cycle_iter_normalize_depth = false;
    double  cycle_iter_log1p_weight = 1.0;
    double  cycle_iter_value = 1.0f; // If dynamic, iter_lim ratio, else iter_lim

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
    f64 dist;
    f32 stripe;

    f32 final_depth;
    f32 final_dist;
    f32 final_stripe;

    //bool flag_for_skip;
};

constexpr int phaseBmpScale(int phase)
{
    if (phase == 0) return 9;
    if (phase == 1) return 3;
    return 1;
}

struct EscapeField : public std::vector<EscapeFieldPixel>
{
    int compute_phase;
    MandelKernelFeatures mandel_features = MandelKernelFeatures::ITER;

    f64 min_depth = 0.0;
    f64 max_depth = 0.0; // unused?

    f128 stable_min_dist{};
    f128 stable_max_dist{};

    f32 min_stripe{};
    f32 max_stripe{};

    f64 log_color_cycle_iters{};
    f64 cycle_dist_value{};

    std::vector<i8> skip_flags;

    int w = 0, h = 0;

    std::vector<DVec2> plots;

    EscapeField(int phase) : compute_phase(phase) {}

    void setAllDepth(f64 value)
    {
        for (int i = 0; i < size(); i++)
        {
            EscapeFieldPixel &p = std::vector<EscapeFieldPixel>::at(i);
            p.depth = value;
            p.dist = value;
            p.stripe = (f32)value;
            //p.dist_128 = { value, 0 };
            //p.stripe_128 = { value, 0 };
            skip_flags[i] = 0;
        }
    }
    void setDimensions(int _w, int _h)
    {
        //if (size() >= (_w * _h)) return;
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

    u8 get_skip_flag(int x, int y)
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
                in[idx(x, y)] = get_skip_flag(x, y) ? 1 : 0;

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
                set_skip_flag(x, y, out[idx(x, y)] != 0);
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
                mask[y * W + x] = get_skip_flag(x, y) ? 1 : 0;

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
            for (int x = 0; x < W; ++x) 
            {
                bool v = ring[y * W + x];
                if (overwrite)
                    set_skip_flag(x, y, v);
                else
                    set_skip_flag(x, y, get_skip_flag(x, y) || v);
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

struct NormalizationPixel : public EscapeFieldPixel
{
    DVec2  stage_pos;
    DDVec2 world_pos;
    bool   is_final;
    float  weight;
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
        local_field = Math::delaunayMeshEllipse<f64>(0, 0, 1.0, sample_r, exponent);
        world_field.resize(local_field.size());
    }

    void updateField(const CameraInfo &camera, double scale)
    {
        f128   world_radius = (scale / camera.relativeZoom<f128>()) * 2.0;
        DDVec2 cam_center_world = camera.pos<f128>();

        // Get stage quad of world ellipse bounds for fast stage interpolation
        bounds.setCamera(camera);
        bounds.setWorldRect(cam_center_world - world_radius, world_radius * 2.0);
        DQuad stage_quad = bounds.stageQuad();

        for (size_t i=0; i< world_field.size(); i++)
        {
            NormalizationPixel& px = world_field[i];
            const DVec2 pt = local_field[i];

            px.world_pos = cam_center_world + (pt * world_radius);
            px.stage_pos = stage_quad.lerpPoint(0.5 + pt * 0.5);

            px.weight = 1.0f - (float)pt.magnitude();
        }
    }

    void clearFinalFlags()
    {
        for (NormalizationPixel& p : world_field)
            p.is_final = false;
    }

    template<typename Callback>
    void forEach(Callback&& callback, int thread_count = Thread::idealThreadCount())
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
