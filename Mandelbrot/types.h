#pragma once
#include <bitloop.h>

#include <cmath>
#include <vector>

enum class MandelKernelFeatures
{
    NONE    = 0b000,
    ITER    = 0b001,
    DIST    = 0b010,
    STRIPES = 0b100,
    MIX_ALL = 0b111,

    COUNT
};

enum struct MandelKernelMode : int
{
    NO_PERTURBATION            = 0,
    PERTURBATION               = 1,
    PERTURBATION_SIMD          = 2,
    PERTURBATION_SIMD_UNROLLED = 3,

    COUNT,
    AUTO = COUNT, // Same value as COUNT to avoid being included in constexpr_dispatch 'build_table'
};

// todo: find way to put inside SIM_BEG ns (wasm32 error, must be in global)
bl_enable_enum_bitops(MandelKernelFeatures);
bl_enable_enum_bitops(MandelKernelMode);

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

    void set(const EscapeFieldPixel& p)
    {
        depth = p.depth;
        dist = p.dist;
        stripe = p.stripe;
    }

    void set(f64 _depth, f64 _dist, f32 _stripe) {
        depth  = _depth;
        dist   = _dist;
        stripe = _stripe;
    }

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

    void reset()
    {
        memset(data(), 0x80, size() * sizeof(EscapeFieldPixel)); // pattern gives each float type a negative value
        memset(skip_flags.data(), 0, skip_flags.size() * sizeof(i8));
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

        // Binary mask from skip_flags (0 or 1)
        std::vector<uint8_t> mask(N);
        for (size_t i = 0; i < N; ++i)
            mask[i] = skip_flags[i] ? 1u : 0u;

        // Summed-area table (integral image) for fast rectangular window sums.
        // Size is (H+1) x (W+1)
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

                // IMPORTANT: your original code treated any pixel whose neighborhood
                // extends out of bounds as NOT "inner" (because safe(nx,ny) == false).
                // That means everything within 'r' pixels of the border is always 0.
                if (x < r || x >= W - r || y < r || y >= H - r)
                {
                    skip_flags[idx] = 0;
                    continue;
                }

                const int x0 = x - r;
                const int y0 = y - r;
                const int x1 = x + r;
                const int y1 = y + r;

                // Rectangle sum via integral image:
                // sum = ps(y1+1,x1+1) - ps(y0,x1+1) - ps(y1+1,x0) + ps(y0,x0)
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

        // r <= 0 is a corner case; preserve original behavior:
        // - The original would produce ring == 0 everywhere.
        //   * overwrite=false -> skip_flags unchanged.
        //   * overwrite=true  -> skip_flags cleared.
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

        // Precompute disk offsets (dx,dy) with dx*dx + dy*dy <= r*r.
        // Cached per thread to avoid recomputation on repeated calls.
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

        // Dilation: for every "1" pixel, turn on all pixels in its disk neighborhood.
        // This is equivalent to your previous per-output-pixel dilate_at loop,
        // but only touches neighborhoods of existing 1s.
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

        // "Outline = dilated AND NOT mask" – keep only newly added perimeter.
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

    template<typename T> requires is_f32<T>
    [[nodiscard]] constexpr const FVec2& worldPos() { return world_pos_32; }

    template<typename T> requires is_f64<T>
    [[nodiscard]] constexpr const DVec2& worldPos() { return world_pos_64; }

    template<typename T> requires is_f128<T>
    [[nodiscard]] constexpr const DDVec2& worldPos() { return world_pos_128; }

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
        local_field = Math::delaunayMeshEllipse<f64>(0, 0, 1.0, sample_r, exponent);
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
