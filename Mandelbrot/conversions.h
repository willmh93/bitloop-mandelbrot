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


SIM_END;
