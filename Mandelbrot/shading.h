#pragma once
#include <bitloop.h>

#include "types.h"


SIM_BEG;



enum class MandelTransform
{
    NONE,
    FLATTEN
};

enum struct MandelShaderFormula
{
    ITER_DIST_STRIPE,          // a+b+c
    ITER_DIST__MUL__STRIPE,   // (a+b)*c
    ITER__MUL__DIST_STRIPE,  // a*(b+c)
    ITER_STRIPE__MULT__DIST,  // a*(b+c)
    COUNT
    //~ // todo: Always tween stripe weight to 0 before switching mode
};

static const char* MandelFormulaNames[(int)MandelShaderFormula::COUNT] = {
    "iter + dist + stripe",
    "iter + (dist * stripe)",
    "(iter * dist) + stripe",
    "(iter * stripe) * dist"
};

static const char* MandelSmoothingNames[(int)KernelFeatures::COUNT] = {
    "NONE",
    "ITER",
    "DIST",
    "STRIPE"
};

template<typename T>
struct GammaLUT
{
    static constexpr int N = 1024;
    T table[N + 1];
    T invN; // 1.0f / N
};

template<typename T>
struct ToneFieldInfo
{
    T minv;
    T maxv;
    T mean;

    T invRange;  // 1 / (maxv - minv)
    T pivot;     // mean in normalized [0,1] space
    T span;      // max distance from pivot to either end

    void configure(T lo, T hi, T _mean)
    {
        minv = lo;
        maxv = hi;
        mean = _mean;

        T range = hi - lo;
        if (range != 0.0f)
            invRange = 1.0f / range;
        else
            invRange = 0.0f;

        // pivot in normalized [0,1] space
        pivot = (mean - lo) * invRange;
        if (pivot < 0.0f) pivot = 0.0f;
        if (pivot > 1.0f) pivot = 1.0f;

        // maximum distance to an edge
        span = std::max(pivot, 1.0f - pivot);
        if (span <= 0.0f)
            span = 1.0f; // avoid divide by zero
    }
};

struct PivotToneParams
{
    float brightness = 0.0f;
    float contrast = 1.0f;
    float gamma = 1.0f;

    bool operator==(const PivotToneParams&) const = default;
};

template<typename T>
inline GammaLUT<T> makeGammaLUT(float gamma)
{
    GammaLUT<T> lut;
    lut.invN = T{ 1 } / GammaLUT<T>::N;

    for (int i = 0; i <= GammaLUT<T>::N; ++i)
    {
        T x = i * lut.invN;
        lut.table[i] = bl::pow(x, gamma);
    }
    return lut;
}

template<typename T>
FORCE_INLINE T applyGammaLUT(const GammaLUT<T>& lut, T x)
{
    if (x < T{0}) x = T{0};
    if (x > T{1}) x = T{1};

    T f = x * GammaLUT<T>::N;
    int   i = (int)f;
    T frac = f - (T)i;

    T a = lut.table[i];
    T b = lut.table[i + 1];

    return a + (b - a) * frac;
}

template<typename T>
FORCE_INLINE T adjustFromPivot(const PivotToneParams& p,
    const ToneFieldInfo<T>& info,
    const GammaLUT<T>& lut,
    T x)
{
    // normalize input x into [0,1] based on min/max
    T t = (x - info.minv) * info.invRange;

    // normalized distance from pivot before contrast
    T d = t - info.pivot;

    // direction and magnitude
    T sign = (d >= T{ 0 }) ? T{ 1 } : T{ -1 };
    T mag = std::fabs(d);

    T span = info.span;
    T r = mag / span; // r = distance as a fraction of span

    // gamma shaping radius as a fixed fraction of span
    const T gammaFrac = T{ 1.0 };

    T r2 = r;
    if (r <= gammaFrac) {
        // Map gammaFrac -> [0,1] for LUT
        T x_norm = r / gammaFrac;
        if (x_norm < T{ 0 }) x_norm = T{ 0 };
        if (x_norm > T{ 1 }) x_norm = T{ 1 };

        T gn = applyGammaLUT(lut, x_norm);

        // Back to "fraction of span" space, but still within [0, gammaFrac]
        r2 = gn * gammaFrac;
    }
    // else: r2 == r  (no shaping outside the gamma band)

    // back to normalized distance units around pivot
    T d2 = sign * (r2 * span);

    // apply contrast after gamma, so gamma radius isn't blown up by contrast
    d2 *= p.contrast;

    // re-center around pivot and add brightness
    T t2 = info.pivot + d2;
    t2 += p.brightness;

    return t2;
}

// same as above, but returns denormalized value in original space
template<typename T>
FORCE_INLINE T adjustFromPivotDenormalized(const PivotToneParams& p,
    const ToneFieldInfo<T>& info,
    const GammaLUT<T>& lut,
    T x)
{
    T ret = adjustFromPivot(p, info, lut, x);
    ret = info.minv + ret * (info.maxv - info.minv);
    return ret;
}

template<typename T>
struct ToneManager
{
    GammaLUT<T> lut;
    ToneFieldInfo<T> info;
    PivotToneParams params;

    void setParams(PivotToneParams p) { params = p; }
    void setParams(float _brightness, float _contrast, float _gamma)
    {
        params.brightness = _brightness;
        params.contrast = _contrast;
        params.gamma = _gamma;
    }

    void configureFieldInfo(T lo, T hi, T mean)
    {
        info.configure(lo, hi, mean);
        lut = makeGammaLUT<T>(params.gamma);
    }

    T apply(T value)
    {
        return adjustFromPivot<T>(params, info, lut, value);
    }

    T applyAndDenormalize(T value)
    {
        return adjustFromPivotDenormalized<T>(params, info, lut, value);
    }
};



SIM_END;
