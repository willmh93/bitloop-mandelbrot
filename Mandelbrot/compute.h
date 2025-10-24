#ifdef _MSC_VER
#pragma float_control(precise, off)
#endif

#include "build_config.h"

#include "shading.h"
#include "conversions.h"

#include <cmath>
#include <complex>
#include <algorithm>

SIM_BEG;

inline FloatingPointType getRequiredFloatType(MandelKernelFeatures smoothing, f128 zoom)
{
    f128 MAX_ZOOM_FLOAT;
    f128 MAX_DOUBLE_ZOOM;

    if (((int)smoothing & (int)MandelKernelFeatures::DIST))
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



inline int mandelbrot_depth(double x0, double y0, int iter_lim)
{
    double x = 0.0, y = 0.0, xx = 0.0, yy = 0.0;
    int iter = 0;

    while (xx + yy <= 4.0 && iter < iter_lim)
    {
        y = (2.0 * x * y + y0);
        x = (xx - yy + x0);
        xx = x * x;
        yy = y * y;
        iter++;
    }

    return iter;
}

inline double mandelbrot_dist(double x0, double y0, int iter_lim)
{
    double x = 0.0, y = 0.0;        // z
    double dx = 1.0, dy = 0.0;        // dz/dc

    for (int i = 0; i < iter_lim; ++i)
    {
        double r2 = x * x + y * y;
        if (r2 > 512.0)
        {
            double r = std::sqrt(r2);                     // |z|
            double dz = std::sqrt(dx * dx + dy * dy);      // |dz|
            if (dz == 0.0) return 0.0;
            return r * std::log(r) / dz;                   // distance
        }

        // save current z for the derivative update
        double xold = x;
        double yold = y;

        // z_{n+1} = z_n^2 + c
        x = xold * xold - yold * yold + x0;
        y = 2.0 * xold * yold + y0;

        // dz_{n+1} = 2 z_n dz_n + 1   (uses xold,yold)
        double dx_new = 2.0 * (xold * dx - yold * dy) + 1.0;
        double dy_new = 2.0 * (xold * dy + yold * dx);

        dx = dx_new;
        dy = dy_new;
    }
    return INSIDE_MANDELBROT_SET;                         // inside set
}

namespace detail
{
    // Complex helpers

    template<class T> struct cplx { T x, y; };

    template<class T>
    FAST_INLINE constexpr void step(cplx<T>& z, const cplx<T>& c)
    {
        T xx = z.x * z.x;
        T yy = z.y * z.y;
        T xy = (z.x * z.y);

        z.x = xx - yy + c.x;
        z.y = (xy + xy) + c.y;
    }

    template<class T>
    FAST_INLINE constexpr void step_d(const cplx<T>& z, cplx<T>& dz)
    {
        const T zx_dzx = z.x * dz.x;
        const T zy_dzy = z.y * dz.y;
        const T zx_dzy = z.x * dz.y;
        const T zy_dzx = z.y * dz.x;

        dz.x = ((zx_dzx - zy_dzy) + (zx_dzx - zy_dzy)) + T(1);
        dz.y = (zx_dzy + zy_dzx) + (zx_dzy + zy_dzx);
    }

    template<class T>
    FAST_INLINE constexpr T mag2(const cplx<T>& z)
    {
        return z.x * z.x + z.y * z.y;
    }


} // namespace detail


template<MandelKernelFeatures Kernel_Features>
constexpr float escape_radius()
{
    return (((int)Kernel_Features & (int)MandelKernelFeatures::DIST) ? 512.0 : 64.0);
}

template<MandelKernelFeatures Kernel_Features>
FAST_INLINE double log_escape_radius_squared()
{
    static double log_escape_r2 = log_as_double(escape_radius<Kernel_Features>());
    return log_escape_r2;
}


template<MandelKernelFeatures Kernel_Features>
constexpr float mandelbrot_smoothing_offset()
{
    constexpr float r2 = escape_radius<Kernel_Features>();
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

template<class T> inline constexpr bool is_flt128_v = false;
template<> inline constexpr bool is_flt128_v<f128> = true;

template<class T, MandelKernelFeatures S>
FAST_INLINE void mandel_kernel(
    const T& x0,
    const T& y0,
    int iter_lim,
    double& depth,
    T& dist,
    T& stripe,
    StripeParams sp = {})
{
    //if (interiorCheck(x0, y0)) {
    //    depth = INSIDE_MANDELBROT_SET_SKIPPED;
    //    if constexpr (((int)S & (int)MandelSmoothing::DIST) != 0) dist = T{ INSIDE_MANDELBROT_SET };
    //    if constexpr (((int)S & (int)MandelSmoothing::STRIPES) != 0) stripes = 0.0;
    //    return;
    //}

    /// @@@ TODO @@@
    /// 
    /// DIST calculation only needs 'float' up to zoom 10000, and 'double' up to zoom 1e26
    /// Currently, you switch to f128 at 2e12 zoom for DIST mode.
    ///
    /// Switch from using single 'T' to 'IterT' and 'DistT'
    /// 
    /// which is actually MUCH higher than where you currently switch to f128,
    /// big gains to be made inbetween 2e12 and 1e26 for DIST mode

    using detail::cplx;
    constexpr bool NEED_DIST        = (bool)((int)S & (int)MandelKernelFeatures::DIST);
    constexpr bool NEED_SMOOTH_ITER = (bool)((int)S & (int)MandelKernelFeatures::ITER);
    constexpr bool NEED_STRIPES     = (bool)((int)S & (int)MandelKernelFeatures::STRIPES);

    constexpr T escape_r2 = T(escape_radius<S>());
    constexpr T zero = T(0), one = T(1);

    int iter = 0;
    T r2 = T{ 0 };

    cplx<T> z{ zero, zero };
    cplx<T> c{ x0, y0 };
    cplx<T> dz{ one, zero };

    // stripe accumulators
    float sum = 0.0, last_added = 0.0;
    int sum_samples = 0;

    while (true)
    {
        if constexpr (NEED_DIST)
            detail::step_d(z, dz);

        cplx<T> z0 = z;

        // step
        detail::step(z, c);
        ++iter;


        // update radius^2
        r2 = detail::mag2(z);

        if constexpr (NEED_STRIPES)
        {
            float a = (float)Math::atan2f_fast((float)z.y, (float)z.x);
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
            const T dz_abs = sqrt(detail::mag2(dz));
            dist = (dz_abs == zero) ? T{ 0 } : (r * log(r) / dz_abs); // log_as_double?
        }
        else {
            dist = T{ -1 };// T{ INSIDE_MANDELBROT_SET };
        }
    }

    if (!escaped)
    {
        depth = INSIDE_MANDELBROT_SET;
        if constexpr (NEED_STRIPES) stripe = T{ 0 };
        return;
    }

    // smooth iteration depth
    if constexpr (NEED_SMOOTH_ITER)
    {
        const double log_abs_z = 0.5 * log_as_double(r2);// (double)log(r2);
        //const double nu = (double)iter + 1.0 - std::log2(std::max(log_abs_z, 1e-30));
        const double nu = (double)iter + 1.0 - std::log2(log_abs_z);
        depth = nu - mandelbrot_smoothing_offset<S>();
    }
    else
    {
        depth = (double)iter;
    }

    if constexpr (NEED_STRIPES)
    {
        float avg = (sum_samples > 0) ? (sum / (float)sum_samples) : 0.0f;
        float prev = (sum_samples > 1) ? ((sum - last_added) / (float)(sum_samples - 1)) : avg;
    
        // stripeAC interpolation weight (fraction inside the last band)
        // frac = 1 + log2( log(ER^2) / log(|z|^2) ), clamped to [0,1]
        float frac = 1.0f + (float)std::log2(log_escape_radius_squared<S>() / std::max(log_as_double(r2), 1e-300));
        frac = clamp01(frac);
    
        float mix = frac * avg + (1.0f - frac) * prev; // linear interpolation
        mix = 0.5f + 0.5f * std::tanh(sp.contrast * (mix - 0.5f)); // optional shaping
        stripe = T{ clamp01(mix) };
    }
}


/*template<class T, MandelSmoothing S>
FAST_INLINE void mandel_kernel(
    const T& x0,
    const T& y0,
    int iter_lim,
    double& depth,
    T& dist,
    float& stripes,
    StripeParams sp = {})
{
    //if (interiorCheck(x0, y0)) {
    //    depth = INSIDE_MANDELBROT_SET_SKIPPED;
    //    if constexpr (((int)S & (int)MandelSmoothing::DIST) != 0) dist = T{ INSIDE_MANDELBROT_SET };
    //    if constexpr (((int)S & (int)MandelSmoothing::STRIPES) != 0) stripes = 0.0;
    //    return;
    //}

    using detail::cplx;
    constexpr bool NEED_DIST = (bool)((int)S & (int)MandelSmoothing::DIST);
    constexpr bool NEED_SMOOTH_ITER = (bool)((int)S & (int)MandelSmoothing::ITER);
    constexpr bool NEED_STRIPES = (bool)((int)S & (int)MandelSmoothing::STRIPES);

    constexpr T escape_r2 = T(escape_radius<S>());
    constexpr T zero = T(0), one = T(1);

    int iter = 0;
    T r2 = T{ 0 };

    cplx<T> z{ zero, zero };
    cplx<T> c{ x0, y0 };
    cplx<T> dz{ one, zero };

    // stripe accumulators
    //float sum = 0.0, last_added = 0.0;
    //int sum_samples = 0;

    float mean_s = 0.5f, m2_s = 0.0f, sum_w = 0.0f;

    while (true)
    {
        if constexpr (NEED_DIST)
            detail::step_d(z, dz);

        cplx<T> z0 = z;

        // step
        detail::step(z, c);
        ++iter;


        // update radius^2
        r2 = detail::mag2(z);

        //if constexpr (NEED_STRIPES)
        //{
        //    //float a = (float)atan2(z.y, z.x);
        //    float a = (float)Math::atan2f_fast((float)z.y, (float)z.x);
        //    last_added = 0.5f + 0.5f * std::sin(sp.freq * a + sp.phase);
        //    sum += last_added;
        //    ++sum_samples;
        //}
        if constexpr (NEED_STRIPES)
        {
            float a = (float)Math::atan2f_fast((float)z.y, (float)z.x);
            float s = 0.5f + 0.5f * std::sin(sp.freq * a + sp.phase);

            // Tail-emphasis EMA (effective window ~ L)
            // Choose L ~ 16..32. alpha = 2/(L+1). You can expose this via sp.
            float L = sp.L_short;
            float alpha = 2.0f / (L + 1.0f);

            // Optional magnitude weighting to suppress far-field dilution
            // Set p in [0.5, 2.0]; p=1 is a good start
            float p = sp.p;
            float w_mag = 1.0f;
            if constexpr (true) { // set to false to disable |z|-weighting
                // w_mag ≈ (ER^2 / r2)^p, clamped to avoid INF/NaN
                double ratio = std::max(1e-30, (double)escape_r2 / std::max((double)r2, 1e-30));
                w_mag = (float)std::pow(ratio, p);
            }

            float w = alpha * w_mag;

            // Keep a few accumulators outside the loop (init to 0 before the loop)
            // sum_w, mean_s, m2_s (for variance)
            // Update running weighted mean/variance (Welford with weights)
            //float prev_sum_w = sum_w;
            sum_w = sum_w + w;
            float delta = s - mean_s;
            float R = (sum_w > 0.0f) ? (w / sum_w) : 0.0f;
            mean_s = mean_s + R * delta;
            float delta2 = s - mean_s;
            m2_s = m2_s + w * delta * delta2; // weighted second moment
        }
        if (r2 > escape_r2 || iter >= iter_lim) break;
    }

    const bool escaped = (r2 > escape_r2) && (iter < iter_lim);

    if constexpr (NEED_DIST)
    {
        if (escaped) {
            const T r = sqrt(r2);
            const T dz_abs = sqrt(detail::mag2(dz));
            dist = (dz_abs == zero) ? T{ 0 } : (r * log(r) / dz_abs); // log_as_double?
        }
        else {
            dist = T{ -1 };// T{ INSIDE_MANDELBROT_SET };
        }
    }

    if (!escaped)
    {
        depth = INSIDE_MANDELBROT_SET;
        if constexpr (NEED_STRIPES) stripes = 0.0;
        return;
    }

    // smooth iteration depth
    if constexpr (NEED_SMOOTH_ITER)
    {
        const double log_abs_z = 0.5 * log_as_double(r2);// (double)log(r2);
        //const double nu = (double)iter + 1.0 - std::log2(std::max(log_abs_z, 1e-30));
        const double nu = (double)iter + 1.0 - std::log2(log_abs_z);
        depth = nu - mandelbrot_smoothing_offset<S>();
    }
    else
    {
        depth = (double)iter;
    }

    if constexpr (NEED_STRIPES)
    {
        // mean
        float mu = (sum_w > 0.0f) ? mean_s : 0.5f;

        // Optional variance normalization for consistent contrast
        // Comment out if you prefer just the mean.
        //float var = (sum_w > 0.0f) ? (m2_s / sum_w) : 0.0f;
        //float sigma = std::sqrt(std::max(var, 1e-8f));

        // Use smooth-iteration "frac" as a final temporal blend if you like
        float frac = 1.0f + (float)std::log2(log_as_double(escape_r2) / std::max(log_as_double(r2), 1e-300));
        frac = clamp01(frac);

        float value = mu; // or (0.5f + 0.5f * (mu - 0.5f) / (k*sigma + eps)) for variance normalization
        // Simple variance normalization with gain k:
        // constexpr float k = 1.2f; value = 0.5f + 0.5f * clamp((mu - 0.5f) / (k*sigma + 1e-6f), -1.0f, 1.0f);

        value = 0.5f + 0.5f * std::tanh(sp.contrast * (value - 0.5f));
        stripes = clamp01(value);
    }

    //if constexpr (NEED_STRIPES)
    //{
    //    float avg = (sum_samples > 0) ? (sum / (float)sum_samples) : 0.0f;
    //    float prev = (sum_samples > 1) ? ((sum - last_added) / (float)(sum_samples - 1)) : avg;
    //
    //    // stripeAC interpolation weight (fraction inside the last band)
    //    // frac = 1 + log2( log(ER^2) / log(|z|^2) ), clamped to [0,1]
    //    float frac = 1.0f + (float)std::log2(log_as_double(escape_r2) / std::max(log_as_double(r2), 1e-300));
    //    frac = clamp01(frac);
    //
    //    float mix = frac * avg + (1.0f - frac) * prev; // linear interpolation
    //    mix = 0.5f + 0.5f * std::tanh(sp.contrast * (mix - 0.5f)); // optional shaping
    //    stripes = clamp01(mix);
    //}
}*/


/*
template<class T, MandelSmoothing S>
FAST_INLINE void mandel_kernel(
    const T& x0,
    const T& y0,
    int iter_lim,
    double& depth,
    T& dist,
    float& stripes,
    StripeParams sp = {})
{
    using detail::cplx;
    constexpr bool NEED_DIST = (bool)((int)S & (int)MandelSmoothing::DIST);
    constexpr bool NEED_SMOOTH_ITER = (bool)((int)S & (int)MandelSmoothing::ITER);
    constexpr bool NEED_STRIPES = (bool)((int)S & (int)MandelSmoothing::STRIPES);

    constexpr T escape_r2 = T(escape_radius<S>());
    constexpr T zero = T(0), one = T(1);

    int iter = 0;
    T r2 = T{ 0 };

    cplx<T> z{ zero, zero };
    cplx<T> c{ x0, y0 };
    cplx<T> dz{ one, zero };      // for distance
    cplx<T> z0;                 // previous z (optional)

    // --- helpers ---
    auto unwrap = [](float prev, float a)->float {
        float d = a - prev;
        while (d > 3.14159265f) d -= 6.28318531f;
        while (d < -3.14159265f) d += 6.28318531f;
        return prev + d;
    };

    // --- stripe angle tracking (unwrapped) ---
    float theta_prev = 0.0f;    // θ at step n-1
    float theta_curr = 0.0f;    // θ at step n
    bool  have_theta = false;

    float sum = 0.0, last_added = 0.0;
    int sum_samples = 0;

    float theta_u = 0.0f, theta_p = 0.0f;


    // --- iterate ---
    while (true)
    {
        if constexpr (NEED_DIST)
            detail::step_d(z, dz);

        z0 = z;

        detail::step(z, c);
        ++iter;

        r2 = detail::mag2(z);

        if constexpr (NEED_STRIPES) {
            float a_raw = (float)Math::atan2f_fast((float)z.y, (float)z.x); // use same fn as in-loop
            if (!have_theta) {
                theta_p = theta_u = a_raw;
                have_theta = true;
            }
            else
            {
                theta_p = theta_u;
                theta_u = unwrap(theta_u, a_raw);
            }
        }

        if (r2 > escape_r2 || iter >= iter_lim)
            break;
    }

    const bool escaped = (r2 > escape_r2) && (iter < iter_lim);

    // --- distance estimate (unchanged) ---
    if constexpr (NEED_DIST)
    {
        if (escaped) {
            const T r = sqrt(r2);
            const T dz_abs = sqrt(detail::mag2(dz));
            dist = (dz_abs == zero) ? T{ 0 } : (r * log(r) / dz_abs);
        }
        else {
            dist = T{ -1 };
        }
    }

    // inside? (no smooth iter / stripes)
    if (!escaped)
    {
        depth = INSIDE_MANDELBROT_SET;
        if constexpr (NEED_STRIPES) stripes = 0.0f;
        return;
    }

    // --- smooth iteration depth (unchanged) ---
    if constexpr (NEED_SMOOTH_ITER)
    {
        const double log_abs_z = 0.5 * (double)log(r2);
        const double nu = (double)iter + 1.0 - std::log2(std::max(log_abs_z, 1e-30));
        depth = nu - mandelbrot_smoothing_offset<S>();
    }
    else
    {
        depth = (double)iter;
    }

    if constexpr (NEED_STRIPES)
    {
        auto arg_fn = [](float y, float x)->float { return (float)Math::atan2f_fast(y, x); };

        // Smooth fractional iterate t. If we hit the cap without escape, treat as full step.
        const double log_abs_z = 0.5 * (double)std::log(std::max((double)r2, 1e-300));
        const double nu = (double)iter + 1.0 - std::log2(std::max(log_abs_z, 1e-300));
        float t = (float)(nu - (double)iter);
        if (iter >= iter_lim && !(r2 > escape_r2)) t = 1.0f;
        if (t < 0.0f) t = 0.0f; else if (t > 1.0f) t = 1.0f;

        // Angles at z_{n-1}, z_n using the SAME arg as in-loop
        const float a0 = arg_fn((float)z0.y, (float)z0.x);
        const float a1 = arg_fn((float)z.y, (float)z.x);

        // Stable signed delta via cross/dot (consistent across neighbors)
        const float u0x = (float)z0.x, u0y = (float)z0.y;
        const float u1x = (float)z.x, u1y = (float)z.y;
        const float inv0 = 1.0f / std::max(std::sqrt(u0x * u0x + u0y * u0y), 1e-30f);
        const float inv1 = 1.0f / std::max(std::sqrt(u1x * u1x + u1y * u1y), 1e-30f);
        const float x0 = u0x * inv0, y0 = u0y * inv0;
        const float x1 = u1x * inv1, y1 = u1y * inv1;
        const float dot = x0 * x1 + y0 * y1;
        const float crs = x0 * y1 - y0 * x1;
        const float dA = std::atan2(crs, dot);   // (-pi, pi]

        // Tail-midpoint angle (anchored to the endpoint you sampled in-loop)
        const float a_tail_mid = a1 - 0.5f * t * dA;

        // Fractional micro-average (use midpoint; you can swap for 2–3 samples if you like)
        const float frac_avg = 0.5f + 0.5f * std::sin(sp.freq * a_tail_mid + sp.phase);

        // Replace the last FULL sample with a FRACTIONAL contribution and a FRACTIONAL denominator
        const float S_prev = sum - last_added;           // remove the last full-sample you added
        const float N_prev = (float)(sum_samples - 1);   // its count
        const float mean = (sum_samples <= 1)
            ? frac_avg
            : (S_prev + t * frac_avg) / (N_prev + t);    // <-- critical: + t, not + 1

        int n_int = (int)std::floor(nu);

        // if you want this to respect your "Frequency" control for integer values:
        int flip = (int)std::floor(sp.freq) * n_int;   // number of π flips
        bool odd = (flip & 1) != 0;

        // evaluate your stripe at the phase you already computed (a_tail_mid)
        float s = std::sin(sp.freq * a_tail_mid + sp.phase);

        // flip per layer if needed (removes the checkerboard fill)
        if (odd) s = -s;

        float val = 0.5f + 0.5f * s;
        val = 0.5f + 0.5f * std::tanh(sp.contrast * (val - 0.5f)); // optional shaping
        stripes = clamp01(val);
    }
}
*/

template<bool smooth>
inline double mandelbrot_spline_iter(double x0, double y0, int iter_lim, ImSpline::Spline& x_spline, ImSpline::Spline& y_spline)
{
    double x = 0.0, y = 0.0, xx = 0.0, yy = 0.0;
    int iter = 0;
    while (xx + yy <= 4.0 && iter < iter_lim)
    {
        y = (2.0 * x * y + y0);
        x = (xx - yy + x0);

        float spline_xx = x_spline(static_cast<float>(x * x));
        float spline_yy = y_spline(static_cast<float>(y * y));

        xx = spline_xx;
        yy = spline_yy;

        iter++;
    }

    // Ensures black for deep-set points
    if (iter == iter_lim)
        return iter_lim;

    if constexpr (smooth)
        return iter + (1.0 - log2(log2(xx + yy) / 2.0));
    else
        return iter;
}

/*template<
    bool Smooth,
    bool Show_Period2_Bulb
>
bool radialMandelbrot()
{
    //double f_max_iter = static_cast<double>(iter_lim);
    return pending_bmp->forEachWorldPixel(camera, current_row, [&](int x, int y, double angle, double point_dist)
    {
        DVec2 polard_coord = cardioid_lerper.originalPolarCoordinate(angle, point_dist, cardioid_lerp_amount);

        // Below x-axis
        if (polard_coord.y < 0)
        {
            pending_bmp->setPixel(x, y, 0, 0, 0, 255);
            return;
        }

        DVec2 mandel_pt = Cardioid::fromPolarCoordinate(polard_coord.x, polard_coord.y);
        double recalculated_orig_angle = cardioid_lerper.originalPolarCoordinate(mandel_pt.x, mandel_pt.y, 1.0).x;

        // Avoid picking pixels from opposite side of the cardioid which would be sampled twice
        bool hide =
            (polard_coord.x < Math::PI && recalculated_orig_angle > Math::PI * 1.1) ||
            (polard_coord.x > Math::PI && recalculated_orig_angle < Math::PI * 0.9);

        if (hide)
        {
            pending_bmp->setPixel(x, y, 0, 0, 0, 255);
            return;
        }

        // Optionally hide everything left of main cardioid
        if constexpr (!Show_Period2_Bulb)
        {
            if (mandel_pt.x < -0.75)
            {
                pending_bmp->setPixel(x, y, 0, 0, 0, 255);
                return;
            }
        }

        double smooth_iter = mandelbrot_iter<Smooth>(mandel_pt.x, mandel_pt.y, iter_lim);

        uint32_t u32;
        iter_gradient_color(smooth_iter, u32);

        //double ratio = smooth_iter / f_max_iter;
        //iter_ratio_color(ratio, r, g, b);

        pending_bmp->setPixel(x, y, u32);
    });
};*/


/// --------------------------------------------------------------
/// Uses the bmp's CanvasObject::worldQuad to to loop over pixels,
/// and calculates mandelbrot set for all pixels which don't already
/// have a depth > 0 (i.e. Filling in the blanks)
/// --------------------------------------------------------------

template<typename T, MandelKernelFeatures Smoothing, bool flatten>
bool mandelbrot(CanvasImage128* bmp, EscapeField* field, int iter_lim, int threads, int timeout, int& current_row, StripeParams stripe_params={})
{
    bool frame_complete = bmp->forEachWorldPixel<T>(current_row, [&](int x, int y, T wx, T wy)
    {
        EscapeFieldPixel& field_pixel = field->at(x, y);
        if (field_pixel.flag_for_skip)
        {
            field_pixel.depth = INSIDE_MANDELBROT_SET_SKIPPED;
            return;
        }

        double depth = field_pixel.depth;
        if (depth >= 0) return;

        T dist{};
        T stripe{};
        mandel_kernel<T, Smoothing>(wx, wy, iter_lim, depth, dist, stripe, stripe_params);

        // Guard against infinite distance, cap at finite min/max
        if (!isfinite(dist))
        {
            if (dist < T{ 0 })
                dist = -std::numeric_limits<T>::min();
            else
                dist = -std::numeric_limits<T>::max();
        }

        field_pixel.depth = depth;    // not likely not need more than float precision
        field_pixel.setStripe(stripe); // store stripe with appropriate precision
        field_pixel.setDist(dist);  // store dist with appropriate precision

    }, threads, timeout);

    return frame_complete;
};

template<typename T>
void normalize_shading_limits(
    EscapeField* field,
    CanvasImage128* bmp,
    CameraInfo& camera,
    const IterParams& iter_params,
    const DistParams& dist_params
)
{
    field->min_depth = std::numeric_limits<f64>::max();
    field->max_depth = std::numeric_limits<f64>::lowest();
    field->min_stripe = std::numeric_limits<f128>::max();
    field->max_stripe = std::numeric_limits<f128>::lowest();

    // Redetermine minimum depth for entire visible field
    bmp->forEachPixel([&](int x, int y)
    {
        EscapeFieldPixel& field_pixel = field->at(x, y);
        double depth = field_pixel.depth;

        if (depth >= INSIDE_MANDELBROT_SET_SKIPPED || field_pixel.flag_for_skip) return;

        if (depth < field->min_depth) field->min_depth = depth;
        if (depth > field->max_depth) field->max_depth = depth;

        T stripe = field_pixel.getStripe<T>();
        if (stripe < field->min_stripe) field->min_stripe = stripe;
        if (stripe > field->max_stripe) field->max_stripe = stripe;
    }, 0);


    // rather than hard min_depth, get average of lowest 25%, and divide by 2?
    // The issue is the noisiness AS you zoom, even with a mean, the diff with that minimum 

    if (field->min_depth == std::numeric_limits<double>::max()) field->min_depth = 0;

    //double avg_depth = sum_depth / (double)(sum_depth_samples);

    double low_ratio = 0;/// 0.25;  // iter_params.cycle_iter_normalize_low_fact / 100.0;
    double high_ratio = 0.5;/// 0.5;  // iter_params.cycle_iter_normalize_high_fact / 100.0;
    /////field->assumed_iter_min = mandelbrotIterLimit(cam_view.relativeZoom<f128>()) * low_ratio;
    field->assumed_iter_min = mandelbrotIterLimit(camera.relativeZoom<f128>()) * low_ratio;
    field->assumed_iter_max = mandelbrotIterLimit(camera.relativeZoom<f128>()) * high_ratio;


    // Calculate normalized depth/dist
    double dist_min_pixel_ratio = ((100.0 - dist_params.cycle_dist_sharpness) / 100.0 + 0.00001);
    f128 stable_min_raw_dist = camera.getTransform().toWorldOffset<f128>(dist_min_pixel_ratio, 0.0).magnitude(); // fraction of a pixel
    f128 stable_max_raw_dist = f128{ (bmp->worldSize().magnitude()) } / f128{ 2.0 }; // half diagonal world viewport size

    // @@ todo: Use downgraded type from: stable_min_raw_dist?
    field->stable_min_dist = (dist_params.cycle_dist_invert ? -1 : 1) * log(stable_min_raw_dist);
    field->stable_max_dist = (dist_params.cycle_dist_invert ? -1 : 1) * log(stable_max_raw_dist);

    ///blPrint() << "stable_min_dist: " << field->stable_min_dist << "   stable_max_dist: " << field->stable_max_dist;

    //double color_cycle_iters;
    //if (iter_params.cycle_iter_dynamic_limit)
    //    color_cycle_iters = (iter_params.cycle_iter_value * (field->assumed_iter_max - (iter_params.cycle_iter_normalize_depth ? field->assumed_iter_min : 0)));
    //else
    //    color_cycle_iters = iter_params.cycle_iter_value;

    if (iter_params.cycle_iter_dynamic_limit)
    {
        double assumed_iter_lim = mandelbrotIterLimit(camera.relativeZoom<f128>()) * 0.5;

        /// "cycle_iter_value" represents ratio of iter_lim
        double color_cycle_iters;
        if (iter_params.cycle_iter_normalize_depth)
            color_cycle_iters = (iter_params.cycle_iter_value * (field->assumed_iter_max - field->assumed_iter_min));
        else
            color_cycle_iters = iter_params.cycle_iter_value * assumed_iter_lim;

        field->log_color_cycle_iters = Math::linear_log1p_lerp(color_cycle_iters, iter_params.cycle_iter_log1p_weight);
    }
    else
    {
        /// "cycle_iter_value" represents actual iter_lim
        field->log_color_cycle_iters = Math::linear_log1p_lerp(iter_params.cycle_iter_value, iter_params.cycle_iter_log1p_weight);
    }

    field->cycle_dist_value = dist_params.cycle_dist_value;

    ///#ifdef BL_DEBUG // todo: Remove when certain bug fixed
    ///if (isnan(field->log_color_cycle_iters))
    ///{
    ///    DebugBreak();
    ///    if (cycle_iter_dynamic_limit)
    ///    {
    ///        double assumed_iter_lim = mandelbrotIterLimit(cam_view.zoom) * 0.5;
    ///
    ///        /// "cycle_iter_value" represents ratio of iter_lim
    ///        double color_cycle_iters = (cycle_iter_value * (assumed_iter_lim - (cycle_iter_normalize_depth ? field->min_depth : 0)));
    ///        field->log_color_cycle_iters = Math::linear_log1p_lerp(color_cycle_iters, cycle_iter_log1p_weight);
    ///    }
    ///    else
    ///    {
    ///        /// "cycle_iter_value" represents actual iter_lim
    ///        field->log_color_cycle_iters = Math::linear_log1p_lerp(cycle_iter_value, cycle_iter_log1p_weight);
    ///    }
    ///}
    ///#endif
}

inline double wrap_range(double x, double low, double high) {
    const double width = high - low;
    //if (!(width > 0)) return low;            // or handle error
    double y = std::fmod(x - low, width);    // wrap around 0..width
    if (y < 0) y += width;                   // fix negatives
    return low + y;                          // back to low..high)
}

template<typename T>
void refreshFieldDepthNormalized(
    EscapeField* pending_field,
    CanvasImage128* pending_bmp,
    MandelKernelFeatures smoothing_type,
    const IterParams& iter_params,
    const DistParams& dist_params,
    ///[[maybe_unused]] bool   cycle_iter_normalize_depth,
    ///double cycle_iter_normalize_low_fact,
    ///[[maybe_unused]] double cycle_iter_normalize_high_fact,
    ///double cycle_iter_log1p_weight,
    ///bool   cycle_dist_invert,
    //double cycle_dist_sharpness,
    int threads)
{
    ///double assumed_iter_min = pending_field->assumed_iter_min;
    ///double assumed_iter_max = pending_field->assumed_iter_max;

    T stable_min_dist = (T)pending_field->stable_min_dist;
    T stable_max_dist = (T)pending_field->stable_max_dist;

    pending_bmp->forEachPixel([&](int x, int y)
    {
        EscapeFieldPixel& field_pixel = pending_field->at(x, y);
        double depth = field_pixel.depth;
        if (depth >= INSIDE_MANDELBROT_SET_SKIPPED || field_pixel.flag_for_skip) return;

        /// ====== ITER =======
        //double low_ratio = 0.25;
        //double high_ratio = 0.5;
        //double floor_depth = iter_params.cycle_iter_normalize_depth ? (pending_field->min_depth * low_ratio) : 0;
        //double ceil_depth = iter_params.cycle_iter_normalize_depth ? (pending_field->max_depth * high_ratio) : 0;

        f64 final_depth;
        
        if (iter_params.cycle_iter_normalize_depth)
        {
            //f64 assumed_log_depth_min = Math::linear_log1p_lerp(pending_field->assumed_iter_min,  iter_params.cycle_iter_log1p_weight);
            //f64 assumed_log_depth_max = Math::linear_log1p_lerp(pending_field->assumed_iter_max,  iter_params.cycle_iter_log1p_weight);

            final_depth = Math::linear_log1p_lerp(depth - pending_field->min_depth, iter_params.cycle_iter_log1p_weight);
        }
        else
            final_depth = Math::linear_log1p_lerp(depth, iter_params.cycle_iter_log1p_weight);

        //double floor_depth = iter_params.cycle_iter_normalize_depth ? pending_field->assumed_iter_min : 0;

        //double depth_window = pending_field->assumed_iter_max - pending_field->assumed_iter_min;
        ///double normalized_log_depth = iter_params.cycle_iter_normalize_depth ?
        ///    wrap_range(log_depth, assumed_log_depth_min, assumed_log_depth_max)
        ///    : log_depth;


        //float final_depth = (float)Math::linear_log1p_lerp(std::max(-0.99999, depth - floor_depth), iter_params.cycle_iter_log1p_weight);
        //float final_depth = (float)normalized_log_depth;// (float)Math::linear_log1p_lerp(std::max(-0.99999, normalized_depth), iter_params.cycle_iter_log1p_weight);

        /// ====== DIST =======
        T raw_dist = field_pixel.getDist<T>();
        if (raw_dist < std::numeric_limits<T>::epsilon())
            raw_dist = std::numeric_limits<T>::epsilon();

        // raw_dist is highest next to mandelbrot set
        T dist = ((int)smoothing_type & (int)MandelKernelFeatures::DIST) ?
            ((dist_params.cycle_dist_invert ? T{-1.0} : T{1.0}) * log(raw_dist))
            : T{0};

        float dist_factor = (float)(dist_params.cycle_dist_invert ?
             -Math::lerpFactor(dist, stable_min_dist, stable_max_dist) :
              Math::lerpFactor(dist, stable_min_dist, stable_max_dist)
        );

        float final_dist = dist_factor;

        /// ====== STRIPE =======
        T stripe = field_pixel.getStripe<T>();

        float final_stripe;
        if (pending_field->min_stripe == pending_field->max_stripe)
            final_stripe = 0.0f;
        else
            final_stripe = (float)Math::lerpFactor<T>(stripe, (T)pending_field->min_stripe, (T)pending_field->max_stripe);

        field_pixel.final_depth = (f32)final_depth;
        field_pixel.final_dist = final_dist;
        field_pixel.final_stripe = final_stripe;
    }, threads);
}


template<MandelShaderFormula F, bool maxdepth_show_optimized>
void shadeBitmap(
    EscapeField* field,
    CanvasImage128* bmp,
    ImGradient* gradient,
    //float max_final_depth,
    //float max_final_dist,
    float iter_weight,
    float dist_weight,
    float stripe_weight,
    int threads
)
{
    double iter_ratio, dist_ratio, stripe_ratio;
    shadingRatios(
        iter_weight, dist_weight, stripe_weight,
        iter_ratio, dist_ratio, stripe_ratio
    );

    //blPrint() << "max_final_depth: " << max_final_depth;
    //blPrint() << "max_final_dist: " << max_final_dist;

    float max_final_depth = (float)field->log_color_cycle_iters;
    float max_final_dist  = (float)field->cycle_dist_value;

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

        float iter_v = field_pixel.final_depth / max_final_depth;// log_color_cycle_iters;
        float dist_v = field_pixel.final_dist / max_final_dist; // cycle_dist_value;
        float stripe_v = field_pixel.final_stripe;

        float combined_t;

        if constexpr (F == MandelShaderFormula::ITER_DIST_STRIPE)
        {
            combined_t = Math::wrap(
                ((iter_v * (float)iter_ratio) +
                (dist_v * (float)dist_ratio)) +
                (stripe_v * (float)stripe_ratio),
                0.0f, 1.0f
            );
        }
        else if constexpr (F == MandelShaderFormula::ITER_DIST__MUL__STRIPE)
        {
            combined_t = Math::wrap(
                ((iter_v * (float)iter_ratio) + (dist_v * (float)dist_ratio)) *
                (stripe_v * (float)stripe_ratio),
                0.0f, 1.0f
            );
        }
        else if constexpr (F == MandelShaderFormula::ITER__MUL__DIST_STRIPE)
        {
            combined_t = Math::wrap(
                (iter_v * (float)iter_ratio) *
                ((dist_v * (float)dist_ratio) * (stripe_v * (float)stripe_ratio)),
                0.0f, 1.0f
            );
        }
        else if constexpr (F == MandelShaderFormula::ITER_STRIPE__MULT__DIST)
        {
            combined_t = Math::wrap(
                (stripe_v * (float)stripe_ratio) *
                ((iter_v * (float)iter_ratio) * (stripe_v * (float)stripe_ratio)),
                0.0f, 1.0f
            );
        }

        gradient->unguardedRGBA(combined_t, u32);

        bmp->setPixel(x, y, u32);
    }, threads);
}

SIM_END;
