#pragma once
#include <bitloop.h>

SIM_BEG;

using namespace bl;

//inline f128 toNormalizedZoom(f128 zoom)
//{
//    return log10(zoom) + 1;
//    //return log(zoom) + 1;
//}
//inline f128 fromNormalizedZoom(f128 normalized_zoom)
//{
//    return pow(normalized_zoom - 1, 10);
//    //return exp(normalized_zoom - 1);
//}
//
//inline f128 toHeight(f128 zoom)
//{
//    return f128(1.0) / toNormalizedZoom(zoom);
//}
//
//inline f128 fromHeight(f128 height)
//{
//    return fromNormalizedZoom(f128(1.0) / height);
//}

// Baselines
constexpr f128 Z0 = f128(1);   // zoom that corresponds to normalized = 1
constexpr f128 Hc = f128(1);   // choose any constant offset for height

// Normalized zoom (decades): n = 1 + log10(zoom/Z0)
inline f128 toNormalizedZoom(f128 zoom) {
    return f128(1) + log10(zoom / Z0);
}
inline f128 fromNormalizedZoom(f128 n) {
    return Z0 * pow(f128(10), n - f128(1));
}

// Height that DECREASES when you zoom in, with linear steps mapping to constant ratios:
inline f128 toHeight(f128 zoom) {
    // height = Hc - log10(zoom/Z0)
    return Hc - log10(zoom / Z0);
}
inline f128 fromHeight(f128 h) {
    // zoom = Z0 * 10^(Hc - h)
    return Z0 * pow(f128(10), Hc - h);
}

inline int mandelbrotIterLimit(f128 zoom)
{
    return (int)(std::max((f128)1, log10(zoom)) * (f128)2000.0);
    //const f128 l = log10(zoom * 400.0);
    //int iters = static_cast<int>(-19.35 * l * l + 741.0 * l - 1841.0);
    //return (100 + (std::max(0, iters))) * 3;
}

inline double qualityFromIterLimit(int iter_lim, f128 zoom)
{
    const int base = mandelbrotIterLimit(zoom);   // always >= 100
    if (base == 0) return 0.0;                      // defensive, though impossible here
    return static_cast<double>(iter_lim) / base;    // <= true quality by < 1/base
}

inline int finalIterLimit(
    CameraInfo& camera,
    double quality,
    bool dynamic_iter_lim,
    bool tweening)
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
