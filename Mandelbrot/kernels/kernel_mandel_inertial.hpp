#pragma once

SIM_BEG;

template<class T>
FORCE_INLINE constexpr void mandel_step_inertial(Complex<T>& z, Complex<T>& z_vel, const Complex<T>& c)
{
    z = z * z + c + z_vel;
    z_vel += z * T(0.02);
    z_vel *= T(0.999);
}

template<class T, KernelFeatures F>
FORCE_INLINE void mandel_kernel_inertial(
    T x0, T y0,
    int iter_lim,
    f64& depth,
    f64& dist,
    StripeAccum& stripe)
{
    constexpr bool NEED_SMOOTH_ITER = (bool)(F & KernelFeatures::ITER);
    constexpr T escape_r2 = T(escape_radius2<F>());
    constexpr T zero = T(0);

    Complex<T> c{ x0, y0 };
    Complex<T> z{ zero, zero };
    Complex<T> z_vel{ 0.0, 0.0 };

    int iter = 0;
    T r2 = T(0);

    while (true)
    {
        mandel_step_inertial(z, z_vel, c);
        ++iter;

        r2 = z.mag2();
        if (r2 > escape_r2 || iter >= iter_lim)
            break;
    }

    const bool escaped = (r2 > escape_r2) && (iter < iter_lim);

    // unused in this kernel
    dist = -1.0;
    stripe.escape(0.0);

    if (!escaped)
    {
        depth = INSIDE_MANDELBROT_SET;
        return;
    }

    if constexpr (NEED_SMOOTH_ITER)
    {
        const f64 log_abs_z = 0.5 * log_as_double(r2);
        const f64 nu = (f64)iter + 1.0 - std::log2(log_abs_z);
        depth = nu - smooth_depth_offset<F>();
        if (depth < 0.0) depth = 0.0;
    }
    else
    {
        depth = (f64)iter;
    }
}

SIM_END;