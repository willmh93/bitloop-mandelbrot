#include "shading.h"
#include "conversions.h"

#include <cmath>
#include <algorithm>



SIM_BEG;

template<class T> struct cplx 
{ 
    T x, y; 
};

inline FloatingPointType getRequiredFloatType(MandelKernelFeatures features, f128 zoom)
{
    f128 MAX_ZOOM_FLOAT;
    f128 MAX_DOUBLE_ZOOM;

    if ((bool)(features & MandelKernelFeatures::DIST))
    {
        MAX_ZOOM_FLOAT = 10.0;
        MAX_DOUBLE_ZOOM = 2e10;
    }
    else
    {
        MAX_ZOOM_FLOAT = 10000;
        MAX_DOUBLE_ZOOM = 1e11;
    }

    if (zoom < MAX_ZOOM_FLOAT) return FloatingPointType::F32;
    if (zoom < MAX_DOUBLE_ZOOM) return FloatingPointType::F64;

    return FloatingPointType::F128;
}

template<class T>
FAST_INLINE constexpr void mandel_step(cplx<T>& z, const cplx<T>& c)
{
    T xx = z.x * z.x;
    T yy = z.y * z.y;
    T xy = (z.x * z.y);

    z.x = xx - yy + c.x;
    z.y = (xy + xy) + c.y;
}

template<class T>
FAST_INLINE constexpr void mandel_step_derivative(const cplx<T>& z, cplx<T>& dz)
{
    const T zx_dzx = z.x * dz.x;
    const T zy_dzy = z.y * dz.y;
    const T zx_dzy = z.x * dz.y;
    const T zy_dzx = z.y * dz.x;

    dz.x = ((zx_dzx - zy_dzy) + (zx_dzx - zy_dzy)) + T(1);
    dz.y = (zx_dzy + zy_dzx) + (zx_dzy + zy_dzx);
}

template<class T>
FAST_INLINE constexpr T cplx_mag2(const cplx<T>& z)
{
    return z.x * z.x + z.y * z.y;
}

template<MandelKernelFeatures Kernel_Features>
constexpr f32 escape_radius()
{
    return ((bool)(Kernel_Features & MandelKernelFeatures::DIST)) ? 512.0f : 64.0f;
}

template<MandelKernelFeatures Kernel_Features>
FAST_INLINE f64 log_escape_radius_squared()
{
    static f64 log_escape_r2 = log_as_double(escape_radius<Kernel_Features>());
    return log_escape_r2;
}


template<MandelKernelFeatures Kernel_Features>
constexpr f32 mandelbrot_smoothing_offset()
{
    constexpr f32 r2 = escape_radius<Kernel_Features>();
    return log2(log2(r2)) - 1.0f;
}

template<typename T>
FAST_INLINE bool interiorCheck(T x0, T y0)
{
    constexpr T a = T(0.25);
    constexpr T b = T(0.0625);
    constexpr T one = T(1);

    const T x_minus_a = x0 - a;
    const T q = x_minus_a * x_minus_a + y0 * y0;

    // Cardioid check (main bulb)
    if (q * (q + x_minus_a) < a * y0 * y0)
        return true;

    // Period-2 bulb check (left-side circle)
    const T x_plus_1 = x0 + one;
    if ((x_plus_1 * x_plus_1 + y0 * y0) < b)
        return true;

    return false;
}


template<class T> FAST_INLINE T clamp01(T v)
{
    return v < T(0) ? T(0) : (v > T(1) ? T(1) : v);
}

//if (interiorCheck(x0, y0)) {
//    depth = INSIDE_MANDELBROT_SET_SKIPPED;
//    if constexpr (((int)S & (int)MandelSmoothing::DIST) != 0) dist = T{ INSIDE_MANDELBROT_SET };
//    if constexpr (((int)S & (int)MandelSmoothing::STRIPES) != 0) stripes = 0.0;
//    return;
//}

template<class T>
inline constexpr bool is_flt128_v = std::is_same_v<std::remove_cvref_t<T>, bl::f128>;

// plain mandelbrot kernet for ITER + DIST + STRIPE
template<class T, MandelKernelFeatures S>
FAST_INLINE void mandel_kernel(
    const T& x0,
    const T& y0,
    int iter_lim,
    f64& depth,
    f64& dist,
    f32& stripe,
    StripeParams sp = {})
{
    constexpr bool NEED_DIST        = (bool)(S & MandelKernelFeatures::DIST);
    constexpr bool NEED_SMOOTH_ITER = (bool)(S & MandelKernelFeatures::ITER);
    constexpr bool NEED_STRIPES     = (bool)(S & MandelKernelFeatures::STRIPES);

    constexpr T escape_r2 = T(escape_radius<S>());
    constexpr T zero = T(0), one = T(1);

    int iter = 0;
    T r2 = T{ 0 };

    cplx<T> z{ zero, zero };
    cplx<T> c{ x0, y0 };
    cplx<T> dz{ one, zero };

    // stripe accumulators
    f32 sum = 0.0, last_added = 0.0;
    int sum_samples = 0;

    while (true)
    {
        if constexpr (NEED_DIST)
            mandel_step_derivative(z, dz);

        cplx<T> z0 = z;

        // step
        mandel_step(z, c);
        ++iter;

        // update radius^2
        r2 = cplx_mag2(z);

        if constexpr (NEED_STRIPES)
        {
            f32 a = (f32)Math::atan2f_fast((f32)z.y, (f32)z.x);
            last_added = 0.5f + 0.5f * std::sin(sp.freq * a + sp.phase);
            sum += last_added;
            ++sum_samples;
        }
        if (r2 > escape_r2 || iter >= iter_lim) break;
    }

    const bool escaped = (r2 > escape_r2) && (iter < iter_lim);

    if constexpr (NEED_DIST)
    {
        if (escaped) {
            const T r = sqrt(r2);
            const T dz_abs = sqrt(cplx_mag2(dz));
            dist = (dz_abs == zero) ? 0.0 : (f64)(r * log(r) / dz_abs); // log_as_double?
        }
        else {
            dist = -1.0;
        }
    }

    if (!escaped)
    {
        depth = INSIDE_MANDELBROT_SET;
        if constexpr (NEED_STRIPES) stripe = 0.0f;
        return;
    }

    if constexpr (NEED_SMOOTH_ITER)
    {
        // smooth iteration depth
        const f64 log_abs_z = 0.5 * log_as_double(r2);
        const f64 nu = (f64)iter + 1.0 - std::log2(log_abs_z);
        depth = nu - mandelbrot_smoothing_offset<S>();
    }
    else
    {
        depth = (f64)iter;
    }

    if constexpr (NEED_STRIPES)
    {
        f32 avg = (sum_samples > 0) ? (sum / (f32)sum_samples) : 0.0f;
        f32 prev = (sum_samples > 1) ? ((sum - last_added) / (f32)(sum_samples - 1)) : avg;
    
        // stripe linear interpolation
        const f64 log_r2 = std::max((f64)std::log((f64)r2), 1e-300);
        f32 frac = 1.0f + (f32)std::log2((f64)log_escape_radius_squared<S>() / log_r2);
        frac = clamp01(frac);
    
        f32 mix = frac * avg + (1.0f - frac) * prev;
        mix = 0.5f + 0.5f * std::tanh(sp.contrast * (mix - 0.5f));
        stripe = T{ clamp01(mix) };
    }
}

template<class T>
struct RefOrbit 
{
    std::vector<cplx<T>> z_pre;      // before step n+1
    std::vector<cplx<T>> z_post;     // after  step n+1
    std::vector<T>       zpost_abs2;

    T cx{}, cy{};

    int  iter_esc = 0;               // first escape iter (n+1), or == iter_lim
    bool escaped = false;
    int  max_ref = 0;
};

template<class T_lo>
struct RefOrbitLo 
{
    std::vector<T_lo> pre_x, pre_y;
    std::vector<T_lo> zpost_abs2;

    T_lo cx{}, cy{};

    int max_ref = 0;
    std::vector<T_lo> pre2x, pre2y;   // 2*Z_pre
    std::vector<T_lo> post_x, post_y; // Z_post
    std::vector<T_lo> slack2;
};


template<class T, MandelKernelFeatures S>
FAST_INLINE void build_ref_orbit(const T& cx, const T& cy, int iter_lim, RefOrbit<T>& out)
{
    constexpr T ER2 = T(escape_radius<S>());
    const T zero = T(0);

    out.cx = cx; 
    out.cy = cy;

    out.z_pre.clear(); 
    out.z_post.clear();
    out.zpost_abs2.clear();

    out.z_pre.reserve(iter_lim);
    out.z_post.reserve(iter_lim);
    out.zpost_abs2.reserve(iter_lim);

    out.iter_esc = iter_lim;
    out.escaped = false;

    cplx<T> z{ zero, zero };
    cplx<T> c{ cx, cy };

    for (int n = 0; n < iter_lim; ++n)
    {
        // pre
        out.z_pre.push_back(z);

        // step
        mandel_step(z, c);

        // post
        out.z_post.push_back(z);
        out.zpost_abs2.push_back(z.x * z.x + z.y * z.y);

        // check for escape
        if (!out.escaped && cplx_mag2(z) > ER2) 
        {
            out.iter_esc = n + 1;
            out.escaped = true;
            break;
        }
    }

    // last post-step index that is still safe to use
    out.max_ref = out.escaped ? (out.iter_esc - 1) : (int)out.z_post.size() - 1;
    if (out.max_ref < 0) out.max_ref = 0;
}

template<class Thi, class T_lo, MandelKernelFeatures S>
inline void downcast_orbit(const RefOrbit<Thi>& hi, RefOrbitLo<T_lo>& lo)
{
    const size_t N = hi.z_pre.size();
    lo.pre_x.resize(N);
    lo.pre_y.resize(N);
    lo.zpost_abs2.resize(N);

    lo.pre2x.resize(N);
    lo.pre2y.resize(N);
    lo.post_x.resize(N);
    lo.post_y.resize(N);
    lo.slack2.resize(N);

    // constants
    const T_lo ER = (T_lo)std::sqrt((f64)escape_radius<S>());

    for (size_t i = 0; i < N; ++i) {
        const T_lo zx = (T_lo)hi.z_pre[i].x;
        const T_lo zy = (T_lo)hi.z_pre[i].y;

        lo.pre_x[i] = zx;
        lo.pre_y[i] = zy;
        lo.zpost_abs2[i] = (T_lo)hi.zpost_abs2[i];

        // 2 * Z_pre
        lo.pre2x[i] = zx + zx;
        lo.pre2y[i] = zy + zy;

        // Z_post directly from high-precision orbit
        lo.post_x[i] = (T_lo)hi.z_post[i].x;
        lo.post_y[i] = (T_lo)hi.z_post[i].y;

        // avoid sqrt in hot loop
        const T_lo zabs = (T_lo)std::sqrt((f64)hi.zpost_abs2[i]);
        const T_lo s = std::max<T_lo>((T_lo)0, ER - zabs);
        lo.slack2[i] = s * s;
    }

    lo.cx = (T_lo)hi.cx;
    lo.cy = (T_lo)hi.cy;
    lo.max_ref = hi.max_ref;
}

template<class T_lo, class T_hi, MandelKernelFeatures F>
FAST_INLINE bool mandel_kernel_approx_rebase(
    const RefOrbitLo<T_lo>& ref,
    T_lo dcx, T_lo dcy,
    int iter_lim,
    f64& depth,
    f64& dist,
    f32& stripe,
    StripeParams sp = {})
{
    constexpr bool NEED_DIST = (bool)(F & MandelKernelFeatures::DIST);
    constexpr bool NEED_STRIPES = (bool)(F & MandelKernelFeatures::STRIPES);
    constexpr T_lo  ER2 = T_lo(escape_radius<F>());

    int ref_index = 0;
    const int max_ref = std::max(0, ref.max_ref);

    // --- delta state ---
    T_lo dx = 0, dy = 0;

    // --- derivative (for distance estimate) ---
    T_lo ddx = T_lo(1);
    T_lo ddy = T_lo(0);

    // --- z (previous step) ---
    T_lo zx_prev = 0, zy_prev = 0;

    // --- stripe accumulators ---
    f32 sum = 0.0, last_added = 0.0;
    int sum_samples = 0;
    f32 freq = sp.freq;
    f32 phase = sp.phase;

    for (int n = 0; n < iter_lim; ++n)
    {
        // 1) update derivative
        if constexpr (NEED_DIST) 
        {
            const T_lo tdx = T_lo(2) * zx_prev * ddx - T_lo(2) * zy_prev * ddy + T_lo(1);
            const T_lo tdy = T_lo(2) * zx_prev * ddy + T_lo(2) * zy_prev * ddx;
            ddx = tdx; 
            ddy = tdy;
        }

        // 2) read reference z_pre and advance delta
        const T_lo c2x = ref.pre2x[ref_index];
        const T_lo c2y = ref.pre2y[ref_index];

        const T_lo sx = c2x + dx;
        const T_lo sy = c2y + dy;

        const T_lo dx_sx = dx * sx;
        const T_lo dy_sy = dy * sy;
        const T_lo dx_sy = dx * sy;
        const T_lo dy_sx = dy * sx;

        dx = dx_sx - dy_sy + dcx;
        dy = dx_sy + dy_sx + dcy;

        // 3) build reference z_post (in T_lo)
        const T_lo zpostx = ref.post_x[ref_index];
        const T_lo zposty = ref.post_y[ref_index];

        const T_lo zax = zpostx + dx; // z_{n+1}.x
        const T_lo zay = zposty + dy; // z_{n+1}.y

        // 4) STRIPES sampling on z_{n+1}
        if constexpr (NEED_STRIPES) {

            f32 a = std::atan2f((f32)zay, (f32)zax);
            f32 s = 0.5f + 0.5f * std::sin(freq * a + phase);

            last_added = s;
            sum += s;
            ++sum_samples;
        }

        // 5) escape test
        const T_lo r2 = zax * zax + zay * zay;
        if (r2 > ER2)
        {
            // smooth ITER
            const f64 r2d = (f64)r2;
            const f64 log_abs_z = 0.5 * std::log(std::max(r2d, 1e-300));
            const f64 nu = f64(n + 1) + 1.0 - std::log2(std::max(log_abs_z, 1e-300));
            depth = nu - mandelbrot_smoothing_offset<F>();

            // DIST
            if constexpr (NEED_DIST) 
            {
                const T_lo r = T_lo(std::sqrt((f64)r2));
                const T_lo dzabs = T_lo(std::hypot((f64)ddx, (f64)ddy));
                dist = (dzabs == 0) ? 0 : (r * std::log(r) / dzabs);
            }
            else {
                dist = 0;
            }

            // STRIPE
            if constexpr (NEED_STRIPES) 
            {
                f32 avg = (sum_samples > 0) ? (sum / (f32)sum_samples) : 0.0f;
                f32 prev = (sum_samples > 1) ? ((sum - last_added) / (f32)(sum_samples - 1)) : avg;

                // stripe linear interpolation
                const f64 log_r2 = std::max((f64)std::log(r2), 1e-300);
                f32 frac = 1.0f + (f32)std::log2((f64)log_escape_radius_squared<F>() / log_r2);
                frac = clamp01(frac);

                f32 mix = frac * avg + (1.0f - frac) * prev;
                mix = 0.5f + 0.5f * std::tanhf(sp.contrast * (mix - 0.5f));
                stripe = clamp01(mix);
            }
            else {
                stripe = 0;
            }

            return true;
        }


        // 6) rebase
        const T_lo d2 = dx * dx + dy * dy;
        bool must_rebase = false;

        // A) after using the last safe ref step, rebase before the next iteration
        if (ref_index == max_ref)
            must_rebase = true;

        // B) boundary-aware rebase: if |z| too large compared to the remaining bailout slack at this ref step, rebase
        if (!must_rebase) 
        {
             constexpr T_lo ALPHA2 = (T_lo)(0.75 * 0.75);
             const T_lo lim2 = ALPHA2 * ref.slack2[ref_index];
             if (lim2 > (T_lo)0 && d2 > lim2)
                 must_rebase = true;
        }

        // C) pauldelbrot rebase (near critical point)
        if (!must_rebase) {
            constexpr T_lo H = (T_lo)0.95;
            if (r2 < H * d2)
                must_rebase = true;
        }

        if (must_rebase) 
        {
            // re-center on the absolute state and restart the reference at m=0
            dx = zax; 
            dy = zay; 

            ref_index = 0;

            // keep dW/dc continuous across the split change
            zx_prev = zax; 
            zy_prev = zay;

            continue;
        }

        // If no rebase, advance to the next reference step (no wrap-around)
        ++ref_index;

        // 7) carry z_{n+1} to next step as z_n
        zx_prev = zax; zy_prev = zay;
    }

    // Did not escape
    depth = INSIDE_MANDELBROT_SET;
    if constexpr (NEED_DIST)    dist   = T_lo(-1);
    if constexpr (NEED_STRIPES) stripe = T_lo(0);
    return true;
}


template<typename T, MandelKernelFeatures F, bool flatten>
bool mandelbrot_perturbation(
    CanvasImage128* bmp, 
    EscapeField* field, 
    NormalizationField& norm_field,
    int iter_lim,
    int threads, 
    int timeout,
    int& current_tile, 
    int m1, int m2, int m3,
    StripeParams stripe_params = {})
{
    typedef f64 T_lo; // f32 requires more rebasing, ends up slower

    RefOrbitLo<T_lo> orbit_lo;

    // 1) calculate high-precision reference
    Vec2<T> anchor; bmp->worldPos<T>(bmp->width() / 2, bmp->height() / 2, anchor.x, anchor.y);
    {
        RefOrbit<T> orbit_hi;
        build_ref_orbit<T, F>(anchor.x, anchor.y, iter_lim, orbit_hi);
        downcast_orbit<T, T_lo, F>(orbit_hi, orbit_lo);
    }

    // 2) quickly calclate pixels using reference
    int compute_tiles;
    switch (field->compute_phase)
    {
    case 0: compute_tiles = (threads * m1);
    case 1: compute_tiles = (threads * m2);
    case 2: compute_tiles = (threads * m3);
    default: compute_tiles = threads;
    }

    f32 tiles_sqrt = ceil(sqrt((f32)compute_tiles));
    int tile_w = (int)ceil((f32)bmp->width() / tiles_sqrt);
    int tile_h = (int)ceil((f32)bmp->height() / tiles_sqrt);

    bool frame_complete = bmp->forEachWorldTilePixel<T>(tile_w, tile_h, current_tile, [&](int x, int y, T wx, T wy)
    {
        EscapeFieldPixel& field_pixel = field->at(x, y);

        f64 depth = field_pixel.depth;
        if (depth >= 0) return;

        const T_lo dcx_lo = (T_lo)(wx - anchor.x);
        const T_lo dcy_lo = (T_lo)(wy - anchor.y);

        f64 dist{};
        f32 stripe{};
        mandel_kernel_approx_rebase<T_lo, T, F>(orbit_lo, dcx_lo, dcy_lo, iter_lim, depth, dist, stripe, stripe_params);

        field_pixel.depth = depth;
        field_pixel.dist = dist;
        field_pixel.stripe = stripe;

    }, compute_tiles, timeout);

    norm_field.forEach([&](NormalizationPixel& field_pixel)
    {
        if (field_pixel.is_final)
            return;

        const T_lo dcx_lo = (T_lo)(field_pixel.world_pos.x - anchor.x);
        const T_lo dcy_lo = (T_lo)(field_pixel.world_pos.y - anchor.y);

        // faster version of: bmp->pixelPosFromWorld(field_pixel.world_pos);
        IVec2 p = field_pixel.stage_pos / phaseBmpScale(field->compute_phase);
        EscapeFieldPixel* existing_pixel = field->get(p.x, p.y);

        if (existing_pixel)
        {
            // always update normalization pixel to nearest pixel (if available),
            // as the previous phase may lack accuracy due to grid snapping
            field_pixel.depth = existing_pixel->depth;
            field_pixel.dist = existing_pixel->dist;
            field_pixel.stripe = existing_pixel->stripe;
            return;
        }

        // coordinate lies outside of viwport rect, do calculation
        f64 depth;
        f64 dist;
        f32 stripe;
        mandel_kernel_approx_rebase<T_lo, T, F>(orbit_lo, dcx_lo, dcy_lo, iter_lim, depth, dist, stripe, stripe_params);

        field_pixel.depth = depth;
        field_pixel.dist = dist;
        field_pixel.stripe = stripe;

        // don't recalculate
        field_pixel.is_final = true;
    }, threads);

    return frame_complete;
};

template<typename T>
void calculate_normalize_info(
    EscapeField* field,
    NormalizationField& norm_field,
    CameraInfo& camera,
    const IterParams& iter_params,
    const DistParams& dist_params,
    int threads
)
{
    blPrint() << "calculate_normalize_info()";

    field->min_depth  = std::numeric_limits<f64>::max();
    field->max_depth  = std::numeric_limits<f64>::lowest();
    field->min_stripe = std::numeric_limits<f32>::max();
    field->max_stripe = std::numeric_limits<f32>::lowest();

    constexpr f32 TAIL_FRAC = 0.05;         // 5% clip per tail
    constexpr f32 MIN_RANGE = 1.0 / 256.0;  // prevent collapse
    constexpr size_t BINS = 256;// 512;
    constexpr u32 W_SCALE = 1u << 20;       // ~1e6 precision

    // --- weighted histogram over raw stripes ---
    /*std::array<std::atomic<u64>, BINS> hist_w{};
    for (auto& a : hist_w) a.store(0, std::memory_order_relaxed);
    std::atomic<u64> total_w{ 0 };*/

    struct WeightHist { 
        u64 total;
        std::array<u64, BINS> hist; 
    };

    std::vector<WeightHist> buckets(threads);

    Thread::forEachBatch(norm_field.world_field, [&](std::span<NormalizationPixel> batch, int ti)
    {
        auto& [total_w, hist_w] = buckets[ti];

        for (NormalizationPixel& field_pixel : batch)
        {
            f64 depth = field_pixel.depth;

            if (depth < INSIDE_MANDELBROT_SET_SKIPPED)
            {
                if (depth < field->min_depth) field->min_depth = depth;
                if (depth > field->max_depth) field->max_depth = depth;

                size_t bin = (size_t)std::floor(field_pixel.stripe * (BINS - 1));
                bin = std::clamp(bin, (size_t)0, BINS - 1);

                // convert weight to fixed-point
                u64 w = (u64)std::llround((f64)field_pixel.weight * (f64)W_SCALE);

                hist_w[bin] += w; // add weight to bin
                total_w     += w; // add weight to sum
            }
        }
    }, threads);

    std::array<u64, BINS> hist_w{};
    u64 total_w = 0;

    for (auto& [total, hist] : buckets)
    {
        for (int bin = 0; bin < BINS; bin++)
            hist_w[bin] += hist[bin];

        total_w += total;
    }

    /*norm_field.forEach([&](NormalizationPixel& field_pixel)
    {
        f64 depth = field_pixel.depth;

        if (depth < INSIDE_MANDELBROT_SET_SKIPPED)
        {
            if (depth < field->min_depth) field->min_depth = depth;
            if (depth > field->max_depth) field->max_depth = depth;

            size_t bin = (size_t)std::floor(field_pixel.stripe_32 * (BINS - 1));
            bin = std::clamp(bin, 0ull, BINS - 1);

            // convert weight to fixed-point
            u64 w = (u64)std::llround((f64)field_pixel.weight * (f64)W_SCALE); 

            // add weight to bin
            hist_w[(size_t)bin].fetch_add(w, std::memory_order_relaxed);

            // add weight to sum
            total_w.fetch_add(w, std::memory_order_relaxed);
        }
    }, threads);*/

    // --- weighted percentiles => min/max ---
    const u64 w_sum = total_w;// .load(std::memory_order_relaxed);
    if (w_sum == 0) {
        field->min_stripe = 0.0;
        field->max_stripe = 1.0;
    }
    else {
        const u64 rank_lo = (u64)std::floor(TAIL_FRAC * (f64)(w_sum - 1));
        const u64 rank_hi = (u64)std::floor((1.0 - TAIL_FRAC) * (f64)(w_sum - 1));

        auto wquantile = [&](u64 rank) -> f32
        {
            u64 acc = 0;
            for (int i = 0; i < BINS; ++i) 
            {
                u64 c = hist_w[(size_t)i];// .load(std::memory_order_relaxed);
                if (acc + c > rank) 
                {
                    const u64 prev = acc;
                    const u64 inside = rank - prev;
                    const f32 frac = (c > 0) ? (f32)inside / (f32)c : 0.0f;
                    return ((f32)i + frac) / (f32)(BINS - 1);
                }
                acc += c;
            }
            return 1.0f; // fallback
        };

        f32 lo = wquantile(rank_lo);
        f32 hi = wquantile(rank_hi);

        if (hi - lo < MIN_RANGE) 
        {
            const f32 mid = 0.5f * (lo + hi);
            lo = std::max(0.0f, mid - 0.5f * MIN_RANGE);
            hi = std::min(1.0f, mid + 0.5f * MIN_RANGE);
            if (hi <= lo) { lo = 0.0f; hi = 1.0f; }
        }
        field->min_stripe = lo;
        field->max_stripe = hi;
    }

    if (field->min_depth == std::numeric_limits<f64>::max()) field->min_depth = 0;

    // Calculate normalized depth/dist
    f64 sharpness_ratio = ((100.0 - dist_params.cycle_dist_sharpness) / 100.0 + 0.00001);

    f128 r1 = f128(0.0001) / camera.relativeZoom<f128>();
    f128 r2 = f128(5.0) / camera.relativeZoom<f128>();

    f128 stable_min_raw_dist = r1 * sharpness_ratio;
    f128 stable_max_raw_dist = r2;

    field->stable_min_dist = (dist_params.cycle_dist_invert ? -1 : 1) * log(stable_min_raw_dist);
    field->stable_max_dist = (dist_params.cycle_dist_invert ? -1 : 1) * log(stable_max_raw_dist); // @@ todo: make linear_log_lerp, make option

    if (iter_params.cycle_iter_dynamic_limit)
    {
        /// "color_cycle_iters" represents ratio of (assumed) iter_lim
        f64 assumed_iter_lim = mandelbrotIterLimit(camera.relativeZoom<f128>()) * 0.5;
        f64 color_cycle_iters = iter_params.cycle_iter_value * assumed_iter_lim;

        field->log_color_cycle_iters = Math::linear_log1p_lerp(color_cycle_iters, iter_params.cycle_iter_log1p_weight);
    }
    else
    {
        /// "cycle_iter_value" represents actual iter_lim
        field->log_color_cycle_iters = Math::linear_log1p_lerp(iter_params.cycle_iter_value, iter_params.cycle_iter_log1p_weight);
    }

    field->cycle_dist_value = dist_params.cycle_dist_value;
}


template<typename T, MandelKernelFeatures F, bool Normalize_Depth>
void normalize_field(
    EscapeField* field,
    CanvasImage128* bmp,
    const IterParams& iter_params,
    const DistParams& dist_params,
    int threads)
{
    f64 stable_min_dist = (f64)field->stable_min_dist;
    f64 stable_max_dist = (f64)field->stable_max_dist;
    float dist_mult = dist_params.cycle_dist_invert ? -1.0f : 1.0f;

    bmp->forEachPixel([&](int x, int y)
    {
        EscapeFieldPixel& field_pixel = field->at(x, y);
        f64 depth = field_pixel.depth;
        if (depth >= INSIDE_MANDELBROT_SET_SKIPPED || field->get_skip_flag(x, y)) return;

        f32 final_depth{}, final_dist{}, final_stripe{};

        /// ====== ITER =======
        if constexpr ((bool)(F & MandelKernelFeatures::ITER))
        {
            if constexpr (Normalize_Depth)
            {
                f64 depth_from_floor = depth - field->min_depth;
                if (depth_from_floor < 0.0) depth_from_floor = 0.0;
                final_depth = (f32)Math::linear_log1p_lerp(depth_from_floor, iter_params.cycle_iter_log1p_weight);
            }
            else
                final_depth = (f32)Math::linear_log1p_lerp(depth, iter_params.cycle_iter_log1p_weight);
        }

        /// ====== DIST =======
        if constexpr ((bool)(F & MandelKernelFeatures::DIST))
        {
            f64 dist          = field_pixel.dist;           // highest when next to mandelbrot interior
            f64 rescaled_dist = dist_mult * log(dist);
            final_dist        = dist_mult * (f32)Math::lerpFactor(rescaled_dist, stable_min_dist, stable_max_dist);

            if (!isfinite(final_dist)) final_dist = 0; // still needed? don't remove unless certain
        }

        /// ====== STRIPE =======
        if constexpr ((bool)(F & MandelKernelFeatures::STRIPES))
        {
            f32 stripe = field_pixel.stripe;
            final_stripe = (f32)Math::lerpFactor(stripe, field->min_stripe, field->max_stripe);
        }

        field_pixel.final_depth  = final_depth;
        field_pixel.final_dist   = final_dist;
        field_pixel.final_stripe = final_stripe;
    }, threads);
}


template<MandelShaderFormula F, bool maxdepth_show_optimized>
void shadeBitmap(
    EscapeField* field,
    CanvasImage128* bmp,
    ImGradient* gradient,
    f32 iter_weight,
    f32 dist_weight,
    f32 stripe_weight,
    int threads
)
{
    f64 iter_ratio, dist_ratio, stripe_ratio;
    shadingRatios(
        iter_weight, dist_weight, stripe_weight,
        iter_ratio, dist_ratio, stripe_ratio
    );

    f32 max_final_depth = (f32)field->log_color_cycle_iters;
    f32 max_final_dist  = (f32)field->cycle_dist_value;

    bmp->forEachPixel([&, field](int x, int y)
    {
        EscapeFieldPixel& field_pixel = field->at(x, y);

        if (field_pixel.depth >= INSIDE_MANDELBROT_SET_SKIPPED)
        {
            if constexpr (maxdepth_show_optimized)
            {
                if (field_pixel.depth == INSIDE_MANDELBROT_SET_SKIPPED)
                    bmp->setPixel(x, y, 0xFF7F007F);
                else
                    bmp->setPixel(x, y, 0xFF000000);
            }
            else
            {
                bmp->setPixel(x, y, 0xFF000000);
            }
            return;
        }

        uint32_t u32;

        f32 iter_v = field_pixel.final_depth / max_final_depth;// log_color_cycle_iters;
        f32 dist_v = field_pixel.final_dist / max_final_dist; // cycle_dist_value;
        f32 stripe_v = field_pixel.final_stripe;

        f32 combined_t;

        if constexpr (F == MandelShaderFormula::ITER_DIST_STRIPE)
        {
            combined_t = Math::wrap(
                ((iter_v * (f32)iter_ratio) +
                (dist_v * (f32)dist_ratio)) +
                (stripe_v * (f32)stripe_ratio),
                0.0f, 1.0f
            );
        }
        else if constexpr (F == MandelShaderFormula::ITER_DIST__MUL__STRIPE)
        {
            combined_t = Math::wrap(
                ((iter_v * (f32)iter_ratio) + (dist_v * (f32)dist_ratio)) *
                (stripe_v * (f32)stripe_ratio),
                0.0f, 1.0f
            );
        }
        else if constexpr (F == MandelShaderFormula::ITER__MUL__DIST_STRIPE)
        {
            combined_t = Math::wrap(
                (iter_v * (f32)iter_ratio) *
                ((dist_v * (f32)dist_ratio) * (stripe_v * (f32)stripe_ratio)),
                0.0f, 1.0f
            );
        }
        else if constexpr (F == MandelShaderFormula::ITER_STRIPE__MULT__DIST)
        {
            combined_t = Math::wrap(
                (stripe_v * (f32)stripe_ratio) *
                ((iter_v * (f32)iter_ratio) * (stripe_v * (f32)stripe_ratio)),
                0.0f, 1.0f
            );
        }

        gradient->unguardedRGBA(combined_t, u32);

        bmp->setPixel(x, y, u32);
    }, threads);
}

SIM_END;
