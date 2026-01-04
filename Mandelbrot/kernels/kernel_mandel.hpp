#pragma once

SIM_BEG;

// --- standard mandelbrot kernel (with ITER/DIST/STRIPE) ---

template<class T>
FORCE_INLINE constexpr void mandel_step(Complex<T>& z, const Complex<T>& c)
{
    // optimized: z = z*z + c;
    const T xx = z.x * z.x;
    const T yy = z.y * z.y;
    const T xy = z.x * z.y;
    
    z.x = xx - yy + c.x;
    z.y = xy + xy + c.y;
}

template<class T>
FORCE_INLINE constexpr void mandel_step_and_derivative(Complex<T>& z, Complex<T>& dz, const Complex<T>& c)
{
    // optimized: dz = (z * dz) * 2 + 1;
    const T zx = z.x;
    const T zy = z.y;

    const T zx_dzx = zx * dz.x;
    const T zy_dzy = zy * dz.y;
    const T zx_dzy = zx * dz.y;
    const T zy_dzx = zy * dz.x;

    const T re = zx_dzx - zy_dzy;
    const T im = zx_dzy + zy_dzx;

    dz.x = re + re + 1;
    dz.y = im + im;

    // optimized: z = z*z + c;
    const T xx = zx * zx;
    const T yy = zy * zy;
    const T xy = zx * zy;

    z.x = xx - yy + c.x;
    z.y = xy + xy + c.y;
}

template<class T, KernelFeatures S>
FORCE_INLINE void kernel_mandel(T x0, T y0, int iter_lim, f64& depth, f64& dist, StripeAccum& stripe)
{
    constexpr bool NEED_DIST        = (bool)(S & KernelFeatures::DIST);
    constexpr bool NEED_SMOOTH_ITER = (bool)(S & KernelFeatures::ITER);
    constexpr bool NEED_STRIPES     = (bool)(S & KernelFeatures::STRIPES);

    constexpr T escape_r2 = T(escape_radius2<S>());
    constexpr T zero = T(0), one = T(1);

    Complex<T> z{ zero, zero };
    Complex<T> c{ x0, y0 };
    Complex<T> dz{ one, zero };
    
    int iter = 0;
    T r2 = T{ 0 };

    while (true)
    {
        if constexpr (NEED_DIST)
            mandel_step_and_derivative(z, dz, c);
        else
            mandel_step(z, c);

        ++iter;

        // update r^2
        r2 = z.mag2();

        if constexpr (NEED_STRIPES)
        {
            float invr = 1.0f / std::sqrt((f32)r2);
            float c = (float)z.x * invr; // cos
            float s = (float)z.y * invr; // sin
            stripe.accumulate(c, s);
        }
        if (r2 > escape_r2|| iter >= iter_lim) break;
    }

    const bool escaped = (r2 > escape_r2) && (iter < iter_lim);

    if constexpr (NEED_DIST)
    {
        if (escaped) {
            const T r = sqrt(r2);
            const T dz_abs = sqrt(dz.mag2());
            dist = (dz_abs == zero) ? 0.0 : (f64)(r * log(r) / dz_abs);
        }
        else
            dist = -1.0;
    }

    if (!escaped)
    {
        depth = INSIDE_MANDELBROT_SET;
        if constexpr (NEED_STRIPES) stripe.escape(0.0);
        return;
    }

    if constexpr (NEED_SMOOTH_ITER)
    {
        // smooth iteration depth
        const f64 log_abs_z = 0.5 * log_as_double(r2);
        const f64 nu = (f64)iter + 1.0 - std::log2(log_abs_z);
        depth = nu - smooth_depth_offset<S>();
        if (depth < 0) depth = 0;
    }
    else
    {
        depth = (f64)iter;
    }

    if constexpr (NEED_STRIPES)
        stripe.escape((f64)r2);
}

SIM_END;