#pragma once

#include <bitloop.h>

#include "../shading/shading.h"
#include "conversions.h"

#include <cmath>
#include <algorithm>

/// official kernels
#include "../kernels/kernel_mandel.hpp"
#include "../kernels/kernel_mandel_perturb.hpp"
#include "../kernels/kernel_mandel_perturb_simd.hpp"
#include "../kernels/kernel_mandel_perturb_simd_unrolled.hpp"

/// experimental kernels
//#include "../kernels/kernel_mandel_inertial.hpp"
//#include "../kernels/kernel_newton.hpp"

SIM_BEG;

template<typename T>
FORCE_INLINE bool interiorCheck(T x0, T y0)
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

/// ----- calculate raw [iter, dist, stripe] values for each field pixel -----
template<typename T, KernelFeatures F, KernelMode K>
bool Mandelbrot_Scene::compute_mandelbrot(int timeout, double& field_compute_ms)
{
    int threads = Thread::threadCount();
    int compute_tiles = threads;

    EscapeField& pending_field = pendingField();
    WorldRasterGrid128& pending_grid = pendingRasterGrid();

    f32 tiles_sqrt = ceil(sqrt((f32)compute_tiles));
    int tile_w = (int)ceil((f32)pending_grid.rasterWidth() / tiles_sqrt);
    int tile_h = (int)ceil((f32)pending_grid.rasterHeight() / tiles_sqrt);
    bool frame_complete = false;

    if constexpr (
        K == KernelMode::PERTURBATION ||
        K == KernelMode::PERTURBATION_SIMD || 
        K == KernelMode::PERTURBATION_SIMD_UNROLLED)
    {
        // get high precision anchor world pos
        Vec2<T> anchor = camera.pos<T>();

        // calculate reference orbit in high precision, downcast to low
        //if (anchor_dd != ref_orbit_anchor)
        {
            // todo: You can only reuse last ref orbit if (iter_lim <= old iter_lim)
            //       For steady-zooms, calculate once at the target depth. For interaction,
            //       increase 'used' ceil(iter_lim, 50)

            RefOrbit<T> orbit_hi;
            build_ref_orbit<T, F>(anchor.x, anchor.y, iter_lim, orbit_hi);
            downcast_orbit<T, T_lo, F>(orbit_hi, ref_orbit_lo);

            //ref_orbit_anchor = anchor_dd;
            //blPrint() << "recalculated reference orbit";
        }

        SimpleTimer timer;

        // run perturbation kernel for remaining tiles pixels
        frame_complete = pending_grid.forEachWorldTilePixel<T>(tile_w, tile_h, P, [&](int x, int y, T wx, T wy)
        {
            EscapeFieldPixel& field_pixel = pending_field.at(x, y);

            if (interiorCheck(wx, wy)) {
                field_pixel.depth = INSIDE_MANDELBROT_SET_SKIPPED;
                return;
            }

            f64 depth = field_pixel.depth;
            if (depth >= 0) return;

            const T_lo dcx_lo = (T_lo)(wx - anchor.x);
            const T_lo dcy_lo = (T_lo)(wy - anchor.y);

            f64 dist;
            StripeAccum stripe(stripe_params.freq);

            #if MANDEL_EXTENDED_FIELD_STATS
            ExtendedFieldStats extended_stats;
            #define MANDEL_ENXTENDED_ARG ,extended_stats
            #else
            #define MANDEL_ENXTENDED_ARG
            #endif

            if constexpr (K == KernelMode::PERTURBATION_SIMD_UNROLLED)
                mandel_kernel_perturb_simd_unrolled<T_lo, F>(ref_orbit_lo, dcx_lo, dcy_lo, iter_lim, depth, dist, stripe);
            else if constexpr (K == KernelMode::PERTURBATION_SIMD)
                mandel_kernel_perturb_simd<T_lo, F>(ref_orbit_lo, dcx_lo, dcy_lo, iter_lim, depth, dist, stripe MANDEL_ENXTENDED_ARG);
            else
                mandel_kernel_perturb<T_lo, F>(ref_orbit_lo, dcx_lo, dcy_lo, iter_lim, depth, dist, stripe);
            //{
            //    f32 _depth, _dist;
            //    mandel_kernel_perturb_32<T_lo, F>(ref_orbit_lo, dcx_lo, dcy_lo, iter_lim, _depth, _dist, stripe);
            //    depth = (f64)_depth;
            //    dist = (f64)_dist;
            //}

            field_pixel.set(depth, dist, stripe);

            #if MANDEL_EXTENDED_FIELD_STATS
            field_pixel.extended_stats = extended_stats;
            #undef MANDEL_ENXTENDED_ARG
            #endif
        }, compute_tiles, timeout);

        field_compute_ms = timer.elapsed();

        // grab most up-to-date field

        // Normalization field isn't rendered. In "low accuracy" mode, use the latest complete field the moment it becomes available
        EscapeField* complete_field = &activeField();

        if (frame_complete)
            complete_field = &pending_field; // prefer high-res field
        
        // run perturbation kernel on normalization points (reusing existing "currently active" field first where acceptable)
        double bmp_scale = phaseDownscaleFactor(complete_field->phase);
        double fdp_tolerance = 9.0 * (1.0 - normalize_field_precision);
        double fdp_tolerance2 = fdp_tolerance * fdp_tolerance;

        norm_field.forEach([&](NormalizationPixel& field_pixel)
        {
            // already done 100% precise manual calculation on previous mid-compute frame? skip
            if (field_pixel.is_final) return;
            
            DVec2 fp = field_pixel.stage_pos / bmp_scale;
            if (fp.x >= 0.0 && fp.y >= 0.0 && fp.x < complete_field->w && fp.y < complete_field->h)
            {
                IVec2 p = (IVec2)fp;
                double fdpf = (fp - ((DVec2)p)).mag2(); // delta dist (pixel fraction)
                double fdp = (fdpf * bmp_scale);
            
                if (fdp < fdp_tolerance2)
                {
                    EscapeFieldPixel& existing_pixel = complete_field->at(p.x, p.y);
                    field_pixel.set(existing_pixel);
                    return;
                }
            }

            // coordinate lies outside of viwport rect, do calculation
            const Vec2<T> world_pos = field_pixel.worldPos<T>();
            const T_lo dcx_lo = (T_lo)(world_pos.x - anchor.x);
            const T_lo dcy_lo = (T_lo)(world_pos.y - anchor.y);

            f64 depth, dist; StripeAccum stripe(stripe_params.freq);


            #if MANDEL_EXTENDED_FIELD_STATS
            ExtendedFieldStats extended_stats; // placeholder
            #define MANDEL_ENXTENDED_ARG ,extended_stats
            #else
            #define MANDEL_ENXTENDED_ARG
            #endif

            if constexpr (K == KernelMode::PERTURBATION_SIMD_UNROLLED)
                mandel_kernel_perturb_simd_unrolled<T_lo, F>(ref_orbit_lo, dcx_lo, dcy_lo, iter_lim, depth, dist, stripe);
            else if constexpr (K == KernelMode::PERTURBATION_SIMD)
                mandel_kernel_perturb_simd<T_lo, F>(ref_orbit_lo, dcx_lo, dcy_lo, iter_lim, depth, dist, stripe MANDEL_ENXTENDED_ARG);
            else
                mandel_kernel_perturb<T_lo, F>(ref_orbit_lo, dcx_lo, dcy_lo, iter_lim, depth, dist, stripe);
            //{
            //    f32 _depth, _dist;
            //    mandel_kernel_perturb_32<T_lo, F>(ref_orbit_lo, dcx_lo, dcy_lo, iter_lim, _depth, _dist, stripe);
            //    depth = (f64)_depth;
            //    dist = (f64)_dist;
            //}

            field_pixel.set(depth, dist, stripe);

            // don't recalculate
            ///field_pixel.is_final = true;

            #if MANDEL_EXTENDED_FIELD_STATS
            #undef MANDEL_ENXTENDED_ARG
            #endif
        });
    }
    else // standard kernel (no perturbation)
    {
        SimpleTimer timer;

        frame_complete = pending_grid.forEachWorldTilePixel<T>(tile_w, tile_h, P, [&](int x, int y, T wx, T wy)
        {
            EscapeFieldPixel& field_pixel = pending_field.at(x, y);

            if (interiorCheck(wx, wy)) {
                field_pixel.depth = INSIDE_MANDELBROT_SET_SKIPPED;
                return;
            }
    
            f64 depth = field_pixel.depth;
            if (depth >= 0) return;
    
            f64 dist; StripeAccum stripe(stripe_params.freq);
            kernel_mandel<T, F>(wx, wy, iter_lim, depth, dist, stripe);
            //mandel_newton_alt_kernel<T, F>(wx, wy, iter_lim, depth, dist, stripe);
            //mandel_kernel_iter_only<T, F>(wx, wy, iter_lim, depth, dist, stripe);

            field_pixel.set(depth, dist, stripe);
    
        }, compute_tiles, timeout);

        field_compute_ms = timer.elapsed();
    
        norm_field.forEach([&](NormalizationPixel& field_pixel)
        {
            // already calculated on previous frame? skip
            if (field_pixel.is_final) return;
    
            IVec2 p = field_pixel.stage_pos / phaseDownscaleFactor(pending_field.phase);
            EscapeFieldPixel* existing_pixel = pending_field.get(p.x, p.y);
    
            if (existing_pixel) {
                field_pixel.set(*existing_pixel);
                return;
            }
    
            // coordinate lies outside of viwport rect, do calculation
            f64 depth, dist; StripeAccum stripe(stripe_params.freq);
            const Vec2<T> world_pos = field_pixel.worldPos<T>();
            kernel_mandel<T, F>(world_pos.x, world_pos.y, iter_lim, depth, dist, stripe);
            field_pixel.set(depth, dist, stripe);
    
            // don't recalculate
            field_pixel.is_final = true;
        });
    }

    return frame_complete;
}

template<KernelFeatures F>
FORCE_INLINE f32 resolveStripe(const StripeAccum& stripe, f32 phase)
{
    const double LOG_ER2 = log_escape_radius2<F>();
    const float cphi = std::cos(phase);
    const float sphi = std::sin(phase);
    return stripeFromAccum(stripe, LOG_ER2, cphi, sphi);
}

template<bool Normalize_Depth>
FORCE_INLINE f32 toFinalDepth(f64 depth, f64 raw_min_depth, f64 log1p_weight, f64 cycle_iters, f32 iter_ratio)
{
    if constexpr (Normalize_Depth)
    {
        f64 depth_from_floor = depth - raw_min_depth;
        if (depth_from_floor < 0.0) depth_from_floor = 0.0;
        return (f32)(math::linearLog1pLerp(depth_from_floor, log1p_weight) / cycle_iters) * iter_ratio;
    }
    else
        return (f32)(math::linearLog1pLerp(depth, log1p_weight) / cycle_iters) * iter_ratio;
}

FORCE_INLINE f64 expEase01_poly5(f64 t)
{
    // approx (exp(t)-1)/(e-1)
    constexpr f64 c1 = 0.5820155208947859;
    constexpr f64 c2 = 0.29049237125998606;
    constexpr f64 c3 = 0.09911540029626577;
    constexpr f64 c4 = 0.020296140335692475;
    constexpr f64 c5 = 0.008080567213269772;
    return t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * c5))));
}

FORCE_INLINE f64 clamp01(f64 x)
{
    return (x < 0.0) ? 0.0 : (x > 1.0 ? 1.0 : x);
}

template<bool Invert>
FORCE_INLINE f32 toFinalDist(
    f64 raw_dist, f64 min_log_dist, f64 max_log_dist,
    ToneManager<f64>& tone, f64 cycle_dist, f32 dist_ratio)
{
    constexpr f64 e_m1 = bl::exp(1.0) - 1.0;
    constexpr f64 inv_e_m1 = 1.0 / e_m1;
    constexpr f64 k_min = -inv_e_m1; // = -1/(e-1)
    constexpr f64 k_floor = -16.0; // exp(-16) ~ 1e-7

    f64 t_raw;
    if (raw_dist <= 0.0)
        t_raw = std::numeric_limits<double>::lowest();
    else
    {
        const f64 inv_range = 1.0 / (max_log_dist - min_log_dist);
        t_raw = (std::log(raw_dist) - min_log_dist) * inv_range;
    }

    f64 t;
    if (t_raw <= k_floor)  t = k_min;
    else if (t_raw >= 1.0) t = 1.0;
    else if (t_raw >= 0.0) t = expEase01_poly5(t_raw); // t = (std::exp(t) - 1.0) / e_m1;
    else                   t = (std::exp(t_raw) - 1.0) * inv_e_m1; // slow path

    t += inv_e_m1;
    t = tone.applyAndDenormalize(t);

    if constexpr (Invert)
        t = -t;

    return (f32)((t / cycle_dist) * (f64)dist_ratio);
}

template<KernelFeatures F>
FORCE_INLINE f32 toFinalStripe(float stripe, float raw_mean_stripe, ToneManager<f32>& tone, f32 stripe_ratio)
{
    //float s = resolveStripe<F>(stripe, phase);
    stripe = (stripe - raw_mean_stripe);
    stripe = tone.apply(stripe) - 0.5f;
    return stripe * stripe_ratio;
}

///template<size_t BINS>
///static void smooth_histogram(std::array<u64, BINS>& hist, int passes = 2)
///{
///    if constexpr (BINS == 0) return;
///    std::array<double, BINS> cur{}, tmp{};
///    double orig_sum = 0.0;
///    for (size_t i = 0; i < BINS; ++i) { cur[i] = (double)hist[i]; orig_sum += cur[i]; }
///    if (orig_sum == 0.0) return;
///
///    // Triangular smoothing (binomial [1,2,1] kernel), 'passes' times
///    for (int p = 0; p < passes; ++p) {
///        // edges: clamp
///        tmp[0] = (3.0 * cur[0] + cur[1]) * 0.25;
///        for (size_t i = 1; i + 1 < BINS; ++i)
///            tmp[i] = (cur[i - 1] + 2.0 * cur[i] + cur[i + 1]) * 0.25;
///        tmp[BINS - 1] = (cur[BINS - 2] + 3.0 * cur[BINS - 1]) * 0.25;
///        cur.swap(tmp);
///    }
///
///    // Renormalize to preserve total weight exactly
///    double sm_sum = 0.0;
///    for (size_t i = 0; i < BINS; ++i) sm_sum += cur[i];
///    double scale = (sm_sum > 0.0) ? (orig_sum / sm_sum) : 1.0;
///    for (size_t i = 0; i < BINS; ++i)
///        hist[i] = (u64)std::llround(cur[i] * scale);
///}
///
///template<size_t BINS>
///static float weighted_quantile_from_hist(const std::array<u64, BINS>& hist, u64 rank, u64 total_w)
///{
///    if (total_w == 0) return 0.0f;
///
///    // Build double CDF
///    std::array<double, BINS> cdf{};
///    double acc = 0.0;
///    for (size_t i = 0; i < BINS; ++i) { acc += (double)hist[i]; cdf[i] = acc; }
///
///    // Find first bin where CDF >= target
///    const double target = (double)rank + 0.5; // mid-rank for smoother behavior
///    size_t i = 0;
///    while (i + 1 < BINS && cdf[i] < target) ++i;
///
///    // Interpolate inside bin using local count
///    const double c_prev = (i == 0) ? 0.0 : cdf[i - 1];
///    const double count = std::max(1.0, cdf[i] - c_prev);
///    const double frac = std::clamp((target - c_prev) / count, 0.0, 1.0);
///
///    // Map to [0,1]
///    return (float)(((double)i + frac) / (double)(BINS - 1));
///}
///
///inline void sort_and_dedup(std::vector<std::pair<float, float>>& ideal_zf_numerator_map)
///{
///    // sort by key, then by value (so identical pairs become adjacent).
///    std::sort(ideal_zf_numerator_map.begin(), ideal_zf_numerator_map.end(),
///        [](const auto& a, const auto& b)
///    {
///        if (a.first < b.first) return true;
///        if (a.first > b.first) return false;
///        return a.second < b.second;
///    });
///
///    // remove exact duplicate pairs (both first and second equal).
///    ideal_zf_numerator_map.erase(
///        std::unique(ideal_zf_numerator_map.begin(), ideal_zf_numerator_map.end(),
///        [](const auto& a, const auto& b)
///    {
///        return a.first == b.first && a.second == b.second;
///    }),
///        ideal_zf_numerator_map.end());
///}

/// ----- determine [low, high, mean] of each feature, and normalized wrapping limits -----
template<typename T, KernelFeatures F>
void Mandelbrot_Scene::calculate_normalize_info()
{
    // todo: move min/max/mean props to normalization field? It's field-independant

    EscapeField& active_field = activeField();
    
    // Calculate normalized dist
    f64 sharpness_ratio = ((100.0 - dist_params.cycle_dist_sharpness) / 100.0 + 0.00001);
    f64 r1 = 0.1 / camera.relativeZoom<f64>();
    f64 r2 = 10.0 / camera.relativeZoom<f64>();
    f64 raw_min_dist = r1 * sharpness_ratio;
    f64 raw_max_dist = r2;
    f64 log_min_dist = std::log(raw_min_dist);
    f64 log_max_dist = std::log(raw_max_dist);

    constexpr size_t BINS = 256;
    constexpr u32 W_SCALE = 1u << 20; // ~1e6 precision

    struct WeightHist 
    { 
        // shared weights
        u64 total_w;
        std::array<u64, BINS> hist_w; // stripe histogram

        // feature * weight sums
        f64 sum_wstripe = 0.0;
    };

    std::vector<WeightHist> buckets(Thread::threadCount());

    auto raw_depth_ranges = Thread::forEachBatch(norm_field.world_field, [&](std::span<NormalizationPixel> batch)
    {
        f64 raw_min_depth = std::numeric_limits<f64>::max();
        f64 raw_max_depth = std::numeric_limits<f64>::lowest();

        for (NormalizationPixel& field_pixel : batch)
        {
            f64 depth = field_pixel.depth;
            if (depth < INSIDE_MANDELBROT_SET_SKIPPED)
            {
                if (depth < raw_min_depth) raw_min_depth = depth;
                if (depth > raw_max_depth) raw_max_depth = depth;
            }
        }

        return DRange{ raw_min_depth, raw_max_depth };
    });

    DRange raw_depth_range = DRange::mergeExpand(raw_depth_ranges);

    u64 total_w = 0;
    f64 total_wstripe = 0;

    if (!steady_zoom)
    {
        // Use normalization field to sample mean stripe value
        Thread::forEachBatch(norm_field.world_field, [&](std::span<NormalizationPixel> batch, int ti)
        {
            auto& bucket = buckets[ti];
            u64& total_w = bucket.total_w;
            auto& hist_w = bucket.hist_w;

            f64& sum_wstripe = bucket.sum_wstripe;

            for (NormalizationPixel& field_pixel : batch)
            {
                f64 depth = field_pixel.depth;

                if (depth < INSIDE_MANDELBROT_SET_SKIPPED)
                {
                    // convert weight to fixed-point
                    u64 w = (u64)std::llround((f64)field_pixel.weight * (f64)W_SCALE);
                    total_w += w;

                    f32 normalized_stripe = resolveStripe<F>(field_pixel.stripe, stripe_params.phase);

                    size_t bin = (size_t)std::floor(normalized_stripe * (BINS - 1));
                    bin = std::clamp(bin, (size_t)0, BINS - 1);
                    hist_w[bin] += w;

                    // sum up weighted features
                    sum_wstripe += w * normalized_stripe;
                }
            }
        });

        std::array<u64, BINS> hist_w{};
       
        for (auto& bucket : buckets)
        {
            std::array<u64, BINS>& bucket_hist_w = bucket.hist_w;

            for (int bin = 0; bin < BINS; bin++)
                hist_w[bin] += bucket_hist_w[bin];

            total_w += bucket.total_w;
            total_wstripe += bucket.sum_wstripe;
        }
    }

    auto stripeMagFromZoom = [](const f128& rel_zoom) -> f32
    {
        const f32 zf = (f32)bl::log10(rel_zoom);
        const f32 base_spline_mag = MandelSplines::stripe_zf_spline(zf);

        if (zf != 0.0f)
            return base_spline_mag / zf;

        // log10(1) might return exactly 0, causing division by zero and NaN stripe magnitude
        static const f32 slope0 = []() -> f32
        {
            constexpr f32 h = 1e-4f;
            const f32 yp = MandelSplines::stripe_zf_spline(h);
            const f32 ym = MandelSplines::stripe_zf_spline(-h);
            return (yp - ym) * (0.5f / h);
        }();

        return slope0;
    };


    // stripe magnitude (fairly good approximation from zoom + spline)
    // avoiding use of histogram for stripe magnitude allows for smoother 'phase' animation
    f32 stripe_mag_dynamic = stripeMagFromZoom(camera.relativeZoom<f128>());

    /// // ---- stripe magnutide (histogram version - what the below approximation tries to mimic) ----
    /// f32 stripe_mag_hist;
    /// smooth_histogram<BINS>(hist_w, 2);
    /// constexpr f32 TAIL_FRAC = 0.05f;
    /// const u64 rank_lo = (u64)std::floor(TAIL_FRAC * (f64)(total_w - 1));
    /// const u64 rank_hi = (u64)std::floor((1.0 - TAIL_FRAC) * (f64)(total_w - 1));
    /// const f32 stripe_lo = weighted_quantile_from_hist<BINS>(hist_w, rank_lo, total_w);
    /// const f32 stripe_hi = weighted_quantile_from_hist<BINS>(hist_w, rank_hi, total_w);
    /// stripe_mag_hist = stripe_hi - stripe_lo;
    /// 
    /// // ---- dev: compare normalized zoom magnitude computed by histogram, collect plots for x/y graph for manual spline approximation ----
    /// if (stripe_weight > 0.1 && zf_raw > 0.1)
    /// {
    ///     f32 ideal_numerator = (stripe_hi - stripe_lo) * zf_raw;
    ///     ideal_zf_numerator_map.push_back(std::pair(zf_raw, ideal_numerator));
    ///     sort_and_dedup(ideal_zf_numerator_map);
    /// }
    /// f32 stripe_mag = stripe_mag_from_hist ? stripe_mag_hist : stripe_mag_dynamic;

    f32 stripe_mag = stripe_mag_dynamic;
    f32 stripe_mean;

    if (!steady_zoom)
        stripe_mean = (f32)(total_wstripe / (f64)total_w);
    else
        stripe_mean = stripe_mean_locked;
    
    // set min/max/mean values
    {
        active_field.raw_min_depth = raw_depth_range.lo;//raw_min_depth;
        active_field.raw_max_depth = raw_depth_range.hi;//raw_max_depth;
        active_field.raw_min_dist = raw_min_dist;
        active_field.raw_max_dist = raw_max_dist;
        active_field.log_min_dist = log_min_dist;
        active_field.log_max_dist = log_max_dist;

        active_field.raw_min_stripe  = stripe_mean - stripe_mag;
        active_field.raw_max_stripe  = stripe_mean + stripe_mag;
        active_field.raw_mean_stripe = stripe_mean;
        active_field.raw_mag_stripe  = stripe_mag;
    }

    dist_tone   .setParams(dist_tone_params);
    stripe_tone .setParams(stripe_tone_params);
    dist_tone   .configureFieldInfo(0.0f, 1.0f, 0.0f);
    stripe_tone .configureFieldInfo(-stripe_mag, stripe_mag, 0.0f);

    if (active_field.raw_min_depth == std::numeric_limits<f64>::max()) 
        active_field.raw_min_depth = 0;

    if (iter_params.iter_dynamic_limit)
    {
        /// "color_cycle_iters" represents ratio of (assumed) iter_lim
        f64 assumed_iter_lim = mandelbrotIterLimit(camera.relativeZoom<f128>());
        f64 color_cycle_iters = iter_params.iter_cycle_value * assumed_iter_lim;

        active_field.log_color_cycle_iters = math::linearLog1pLerp(color_cycle_iters, iter_params.iter_log1p_weight);
    }
    else
    {
        /// "iter_cycle_value" represents actual iter_lim
        active_field.log_color_cycle_iters = math::linearLog1pLerp(iter_params.iter_cycle_value, iter_params.iter_log1p_weight);
    }
}


/// ----- calculates final_iter, final_dist, final_stripe & prepares data for shader -----
template<typename T, KernelFeatures F, bool Normalize_Depth, bool Invert_Dist, bool Show_Optimized>
void Mandelbrot_Scene::normalize_field()
{
    // todo: implement with GLSL shader

    EscapeField& field = activeField();
    WorldRasterGrid128& grid = activeRasterGrid();

    // raw limits / mean values
    f64 raw_min_depth   = field.raw_min_depth;
    f64 raw_max_depth   = field.raw_max_depth;

    //f64 raw_min_dist    = field->raw_min_dist;
    f64 raw_max_dist    = field.raw_max_dist;
    f64 log_min_dist    = field.log_min_dist;
    f64 log_max_dist    = field.log_max_dist;

    f32 raw_min_stripe  = field.raw_min_stripe;
    f32 raw_max_stripe  = field.raw_max_stripe;
    f32 raw_mean_stripe = field.raw_mean_stripe;

    // normalized (not final) cycle sizes for repeating patterns
    f64 cycle_iters     = field.log_color_cycle_iters;
    f32 cycle_dist      = (f32)dist_params.dist_cycle_value;

    f32 iter_ratio, dist_ratio, stripe_ratio;
    shadingRatios(
        iter_weight, dist_weight, stripe_weight,
        iter_ratio, dist_ratio, stripe_ratio
    );

    field.final_min_depth  = toFinalDepth<Normalize_Depth>(raw_min_depth, raw_min_depth, iter_params.iter_log1p_weight, cycle_iters, iter_ratio);
    field.final_max_depth  = toFinalDepth<Normalize_Depth>(raw_max_depth, raw_min_depth, iter_params.iter_log1p_weight, cycle_iters, iter_ratio);
    field.final_min_dist   = toFinalDist<Invert_Dist>(0.0, log_min_dist, log_max_dist, dist_tone, cycle_dist, dist_ratio);
    field.final_max_dist   = toFinalDist<Invert_Dist>(raw_max_dist, log_min_dist, log_max_dist, dist_tone, cycle_dist, dist_ratio);
    field.final_min_stripe = toFinalStripe<F>(raw_min_stripe, raw_mean_stripe, stripe_tone, stripe_ratio);
    field.final_max_stripe = toFinalStripe<F>(raw_max_stripe, raw_mean_stripe, stripe_tone, stripe_ratio);

    grid.forEachPixel([
        this, &field,
        raw_min_depth, raw_max_depth, log_min_dist, log_max_dist, raw_mean_stripe,
        cycle_iters, cycle_dist,
        iter_ratio, dist_ratio, stripe_ratio
    ](int x, int y)
    {
        EscapeFieldPixel& field_pixel = field.at(x, y);
        f64 depth = field_pixel.depth;

        if constexpr (Show_Optimized)
        {
            if (field.get_skip_flag(x, y))
            {
                field_pixel.final_depth = -2.0f;
                return;
            }
        }

        if (depth >= INSIDE_MANDELBROT_SET_SKIPPED)
        {
            field_pixel.final_depth = -1.0f;
            return;
        }

        f32 final_depth{};
        f32 final_dist{};
        f32 final_stripe{};

        /// ====== ITER =======
        if constexpr ((bool)(F & KernelFeatures::ITER))
            final_depth = toFinalDepth<Normalize_Depth>(depth, raw_min_depth, iter_params.iter_log1p_weight, cycle_iters, iter_ratio);

        /// ====== DIST =======
        if constexpr ((bool)(F & KernelFeatures::DIST))
            final_dist = toFinalDist<Invert_Dist>(field_pixel.dist, log_min_dist, log_max_dist, dist_tone, cycle_dist, dist_ratio);

        /// ====== STRIPE =======
        if constexpr ((bool)(F & KernelFeatures::STRIPES))
        {
            float stripe = resolveStripe<F>(field_pixel.stripe, stripe_params.phase);
            final_stripe = toFinalStripe<F>(stripe, raw_mean_stripe, stripe_tone, stripe_ratio);
        }

        field_pixel.final_depth  = final_depth;
        field_pixel.final_dist   = final_dist;
        field_pixel.final_stripe = final_stripe;
    });

    field.fillFeaturesTexureData();
}

SIM_END;
