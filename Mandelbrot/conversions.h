#pragma once
#include <bitloop.h>

SIM_BEG;

using namespace bl;

template<KernelFeatures F> constexpr f64 escape_radius() { return ((bool)(F & KernelFeatures::DIST)) ? 32.0 : 8.0; }
template<KernelFeatures F> constexpr f64 escape_radius2() { return bl::pow(escape_radius<F>(), 2.0); }
template<KernelFeatures F> constexpr f64 log_escape_radius2() { return bl::log(escape_radius2<F>()); }
template<KernelFeatures F> constexpr f32 smooth_depth_offset() { return (f32)bl::log2(bl::log2(escape_radius2<F>())) - 1.0f; }

// baselines
constexpr f128 Z0 = 1; // base zoom
constexpr f128 H0 = 1; // base height

// normalized zoom (begins at 1, increases as you zoom)
inline f128 toNormalizedZoom(f128 zoom) { return f128(1) + log10(zoom / Z0); }
inline f128 fromNormalizedZoom(f128 n)  { return Z0 * pow(f128(10), n - f128(1)); }

// height (begins at 1, decreases as you zoom - can be negative)
inline f128 toHeight(f128 zoom)  { return H0 - log10(zoom / Z0); }
inline f128 fromHeight(f128 h)   { return Z0 * pow(f128(10), H0 - h); }

// "base" iter limit
inline int mandelbrotIterLimit(f128 zoom)                   { return (int)(std::max((f128)1, log10(zoom)) * (f128)2000.0); }
inline double qualityFromIterLimit(int iter_lim, f128 zoom) { return (f64)(iter_lim) / mandelbrotIterLimit(zoom); }

// iter limit used (if dynamic, quality=ratio, otherwise raw iter limit). force dynamic during tween
inline int finalIterLimit(const CameraInfo& camera, f64 quality, bool dynamic_iter_lim, bool tweening)
{
    if (dynamic_iter_lim)
        return (int)(mandelbrotIterLimit(camera.relativeZoom<f128>()) * quality);
    else
    {
        int iters = (int)(quality);
        if (tweening)
            iters = std::min(iters, (int)(mandelbrotIterLimit(camera.relativeZoom<f128>()) * 0.25f));
        return iters;
    }
}

inline FloatingPointType getRequiredFloatType(KernelFeatures features, f128 zoom)
{
    f128 MAX_F64_ZOOM;

    // dist is more sensitive to precision, lower zoom threshold if enabled
    if ((bool)(features & KernelFeatures::DIST))
        MAX_F64_ZOOM = 2e10;
    else
        MAX_F64_ZOOM = 1e11;

    if (zoom < MAX_F64_ZOOM) return FloatingPointType::F64;
    return FloatingPointType::F128;
}

inline void shadingRatios(
    double  iter_weight, double  dist_weight, double  stripe_weight,
    double& iter_ratio,  double& dist_ratio,  double& stripe_ratio)
{
    if (iter_weight + dist_weight + stripe_weight == 0.0)
        iter_weight = 1.0;

    double sum = iter_weight + dist_weight + stripe_weight;
    iter_ratio   = iter_weight / sum;
    dist_ratio   = dist_weight / sum;
    stripe_ratio = stripe_weight / sum;
}

/*

A = lerp(iter, dist, 


*/

//inline float mixKernelFeatures(
//    float iter, 
//    float dist, 
//    float stripe,
//    float lerp_iter_dist,   // (iter * dist)   weight
//    float lerp_iter_stripe, // (iter * stripe) weight
//    float lerp_dist_stripe  // (dist * stripe) weight
//)
//{
//    float base = iter + dist + stripe;
//    float da = iter * dist   - iter - dist;
//    float db = iter * stripe - iter - stripe;
//    float dc = dist * stripe - dist - stripe;
//
//    return base
//        + lerp_iter_dist * da
//        + lerp_iter_stripe * db
//        + lerp_dist_stripe * dc;
//}

///inline float mixKernelFeatures(
///    float iter,
///    float dist,
///    float stripe,
///    float lerp_iter_dist,   // (iter, dist)   pair
///    float lerp_iter_stripe, // (iter, stripe) pair
///    float lerp_dist_stripe  // (dist, stripe) pair
///)
///{
///    float a = iter;
///    float b = dist;
///    float c = stripe;
///
///    float base = a + b + c;
///
///    // _ChatGPT_ Pairwise "multiply vs add" corrections
///    float da = a * b - a - b; // for (iter, dist)
///    float db = a * c - a - c; // for (iter, stripe)
///    float dc = b * c - b - c; // for (dist, stripe)
///
///    float wa = lerp_iter_dist;
///    float wb = lerp_iter_stripe;
///    float wc = lerp_dist_stripe;
///
///    // _ChatGPT_ Triple corrections (derived, not extra user controls)
///    // _ChatGPT_ Activates near iter*(dist+stripe)
///    float ga = std::max(0.0f, lerp_iter_dist + lerp_iter_stripe - 1.0f);
///
///    // _ChatGPT_ Activates near (iter+dist)*stripe
///    float gc = std::max(0.0f, lerp_iter_stripe + lerp_dist_stripe - 1.0f);
///
///    return base
///        + lerp_iter_dist   * da
///        + lerp_iter_stripe * db
///        + wc * dc
///        + ga * a   // fix for a*(b + c)
///        + gc * c;  // fix for (a + b)*c
///}

inline float mixKernelFeatures(
    float iter,
    float dist,
    float stripe,
    float iter_x_dist_weight,
    float dist_x_stripe_weight,
    float stripe_x_iter_weight,
    float iter_x_distStripe_weight,
    float dist_x_iterStripe_weight,
    float stripe_x_iterDist_weight)
{
    // _ChatGPT_ Common base: pure additive mix
    float base = iter + dist + stripe;

    // _ChatGPT_ Common products
    float iter_dist = iter * dist;
    float dist_stripe = dist * stripe;
    float stripe_iter = stripe * iter; // same as iter * stripe

    // _ChatGPT_ Pairwise "product vs sum" deltas (pair product + third additive)
    float d_iter_x_dist = iter_dist - iter - dist;   // (iter * dist + stripe)  - base
    float d_dist_x_stripe = dist_stripe - dist - stripe; // (dist * stripe + iter)  - base
    float d_stripe_x_iter = stripe_iter - stripe - iter; // (stripe * iter + dist)  - base

    // _ChatGPT_ Triple "sum-then-multiply" deltas
    float d_iter_x_distStripe = (iter_dist + stripe_iter) - base; // iter * (dist + stripe)  - base
    float d_dist_x_iterStripe = (iter_dist + dist_stripe) - base; // dist * (iter + stripe)  - base
    float d_stripe_x_iterDist = (stripe_iter + dist_stripe) - base; // stripe * (iter + dist)  - base

    return base
        + iter_x_dist_weight * d_iter_x_dist
        + dist_x_stripe_weight * d_dist_x_stripe
        + stripe_x_iter_weight * d_stripe_x_iter
        + iter_x_distStripe_weight * d_iter_x_distStripe
        + dist_x_iterStripe_weight * d_dist_x_iterStripe
        + stripe_x_iterDist_weight * d_stripe_x_iterDist;
}


SIM_END;
