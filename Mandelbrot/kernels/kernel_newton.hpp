#pragma once

SIM_BEG;

template<class T, KernelFeatures S>
FORCE_INLINE void newton_kernel(
    T x0, T y0,
    int iter_lim,
    f64& depth,
    f64& dist,
    StripeAccum& stripe)
{
    constexpr bool NEED_DIST = (bool)(S & KernelFeatures::DIST);
    constexpr bool NEED_SMOOTH_ITER = (bool)(S & KernelFeatures::ITER);
    constexpr bool NEED_STRIPES = (bool)(S & KernelFeatures::STRIPES);

    constexpr T zero = T(0);
    constexpr T one = T(1);

    const T   eps2 = T(1e-12);
    const f64 eps = std::sqrt((f64)eps2);

    int iter = 0;

    Complex<T> z{ x0, y0 };
    Complex<T> fz{ zero, zero };
    Complex<T> dfz{ zero, zero };
    Complex<T> dfz_prev{ zero, zero };

    T   f2 = T(0);
    T   f2_prev = T(0);

    bool have_cross = false;
    f64 smooth_iter = 0.0;
    f64 t_cross = 1.0;

    bool have_dist_cross = false;
    f64 dist_smooth = 0.0;

    while (true)
    {
        // keep previous residuals / derivative
        f2_prev = f2;
        dfz_prev = dfz;
        const bool prev_below = (f2_prev <= eps2);

        // f(z) = z^3 - 1,  f'(z) = 3 z^2
        Complex<T> z2 = z * z;
        Complex<T> z3 = z2 * z;

        fz = z3 - Complex<T>(one, zero);
        dfz = T(3) * z2;

        f2 = cplx_mag2(fz);
        const bool converged = (f2 <= eps2);
        const T    df2 = cplx_mag2(dfz);
        const bool bad_derivative = (df2 == zero);

        // --------- smooth ITER (fractional iteration) ----------
        if constexpr (NEED_SMOOTH_ITER)
        {
            if (!prev_below && converged && iter > 0 && !have_cross)
            {
                const f64 fp = std::sqrt((f64)f2_prev);
                const f64 fc = std::sqrt((f64)f2);
                const f64 log_fp = std::log(fp);
                const f64 log_fc = std::log((fc > 0.0) ? fc : eps * 1e-6);
                const f64 log_eps = std::log(eps);

                f64 t = 0.0;
                if (log_fc != log_fp)
                    t = (log_eps - log_fp) / (log_fc - log_fp);

                if (t < 0.0) t = 0.0;
                else if (t > 1.0) t = 1.0;

                smooth_iter = (f64)(iter - 1) + t;
                t_cross = t;
                have_cross = true;
            }
        }

        // --------- smooth DIST across the same crossing ----------
        if constexpr (NEED_DIST)
        {
            if (!prev_below && converged && iter > 0 && !have_dist_cross)
            {
                const f64 fp = std::sqrt((f64)f2_prev);
                const f64 fc = std::sqrt((f64)f2);
                const f64 dfp_abs = std::sqrt((f64)cplx_mag2(dfz_prev));
                const f64 dfc_abs = std::sqrt((f64)cplx_mag2(dfz));

                const f64 d_prev = (dfp_abs == 0.0) ? 0.0 : fp / dfp_abs;
                const f64 d_curr = (dfc_abs == 0.0) ? 0.0 : fc / dfc_abs;

                f64 t = 0.0;
                if (have_cross)
                {
                    t = t_cross;
                }
                else
                {
                    const f64 log_fp = std::log(fp);
                    const f64 log_fc = std::log((fc > 0.0) ? fc : eps * 1e-6);
                    const f64 log_eps = std::log(eps);

                    if (log_fc != log_fp)
                        t = (log_eps - log_fp) / (log_fc - log_fp);

                    if (t < 0.0) t = 0.0;
                    else if (t > 1.0) t = 1.0;
                }

                dist_smooth = (1.0 - t) * d_prev + t * d_curr;
                have_dist_cross = true;
            }
        }

        // --------- STRIPE accumulation ----------
        if constexpr (NEED_STRIPES)
        {
            T r2 = cplx_mag2(z);
            if (r2 > zero)
            {
                float invr = 1.0f / std::sqrt((f32)r2);
                float c = (float)z.x * invr;
                float s = (float)z.y * invr;
                stripe.accumulate(c, s);
            }
        }

        if (converged || bad_derivative || iter >= iter_lim)
            break;

        // Newton step
        Complex<T> step = fz / dfz;
        z -= step;

        ++iter;
    }

    const bool converged = (f2 <= eps2) && (iter < iter_lim);

    // --------- final DIST ----------
    if constexpr (NEED_DIST)
    {
        if (converged)
        {
            if (have_dist_cross)
            {
                dist = dist_smooth;
            }
            else
            {
                const T f_abs = sqrt(f2);
                const T df_abs = sqrt(cplx_mag2(dfz));
                dist = (df_abs == zero) ? 0.0 : (f64)(f_abs / df_abs);
            }
        }
        else
        {
            dist = -1.0;
        }
    }

    // bail out early for non-converged points
    if (!converged)
    {
        depth = INSIDE_MANDELBROT_SET;
        if constexpr (NEED_STRIPES)
            stripe.escape(0.0);      // keep your existing "no convergence" behaviour
        return;
    }

    // --------- final ITER depth ----------
    if constexpr (NEED_SMOOTH_ITER)
    {
        f64 it = have_cross ? smooth_iter : (f64)iter;
        if (it < 0.0) it = 0.0;
        depth = it;
    }
    else
    {
        depth = (f64)iter;
    }

    // --------- final STRIPE encoding (NO banding) ----------
    if constexpr (NEED_STRIPES)
    {
        const f64 log_er2 = log_escape_radius2<S>();

        if (converged)
        {
            const f64 t = have_cross ? t_cross : 1.0;

            const f64 scale = std::pow(2.0, t - 1.0);
            const f64 lr = log_er2 / scale;

            stripe.log_r2_at_escape = lr;
            stripe.escaped = true;
        }
        else
        {
            stripe.escape(0.0);
        }
    }
}


/*template<class T, KernelFeatures S>
FORCE_INLINE void mandel_newton_alt_kernel(
    T x0, T y0,
    int iter_lim,
    f64& depth,
    f64& dist,
    StripeAccum& stripe)
{
    constexpr bool NEED_DIST = (bool)(S & KernelFeatures::DIST);
    constexpr bool NEED_SMOOTH_ITER = (bool)(S & KernelFeatures::ITER);
    constexpr bool NEED_STRIPES = (bool)(S & KernelFeatures::STRIPES);

    constexpr T escape_r2 = T(escape_radius2<S>());
    constexpr T zero = T(0), one = T(1);

    // Newton convergence threshold (same spirit as your newton_kernel)
    const T   eps2 = T(1e-12);
    const f64 eps = std::sqrt((f64)eps2);

    int step_i = 0;          // total alternating steps executed
    int newt_i = 0;          // number of completed Newton updates (for smoothing bookkeeping)

    Complex<T> z{ zero, zero }; // Mandel-style start (parameter-plane hybrid)
    Complex<T> c{ x0, y0 };

    // derivative accumulator (same convention as your kernel_mandel)
    Complex<T> dz{ one, zero };

    // book-keeping
    T r2 = T(0);

    // Newton residual / derivative state for smoothing
    Complex<T> fz{ zero, zero };
    Complex<T> dfz{ zero, zero };
    Complex<T> dfz_prev{ zero, zero };
    T f2 = T(0);
    T f2_prev = T(0);

    bool have_cross = false;
    f64  t_cross = 1.0;
    f64  smooth_newt = 0.0;      // in "Newton-evaluation index" units (like your newton_kernel)

    bool have_dist_cross = false;
    f64  dist_smooth = 0.0;

    bool escaped = false;
    bool converged = false;
    bool bad_newton = false;

    while (true)
    {
        if (step_i >= iter_lim)
            break;

        const bool do_mandel = (step_i % 3) == 0;// ((step_i & 1) == 0);

        if (do_mandel)
        {
            if constexpr (NEED_DIST)
                mandel_step_derivative(z, dz);

            mandel_step(z, c);
            ++step_i;

            r2 = cplx_mag2(z);

            if constexpr (NEED_STRIPES)
            {
                float invr = 1.0f / std::sqrt((f32)r2);
                float cc = (float)z.x * invr;
                float ss = (float)z.y * invr;
                stripe.accumulate(cc, ss);
            }

            if (r2 > escape_r2)
            {
                escaped = true;
                break;
            }
        }
        else
        {
            // --- evaluate Newton state at current z (this is the "odd" phase) ---
            f2_prev = f2;
            dfz_prev = dfz;
            const bool prev_below = (f2_prev <= eps2);

            // f(z) = z^3 - 1,  f'(z) = 3 z^2,  f''(z) = 6 z
            Complex<T> z2 = z * z;
            Complex<T> z3 = z2 * z;

            fz = z3 - Complex<T>(one, zero);
            dfz = T(3) * z2;

            f2 = cplx_mag2(fz);
            const bool now_conv = (f2 <= eps2);

            const T df2 = cplx_mag2(dfz);
            const bool bad_derivative = (df2 == zero);

            // --- smooth ITER (fractional) across residual crossing (between Newton evaluations) ---
            if constexpr (NEED_SMOOTH_ITER)
            {
                if (!prev_below && now_conv && newt_i > 0 && !have_cross)
                {
                    const f64 fp = std::sqrt((f64)f2_prev);
                    const f64 fc = std::sqrt((f64)f2);

                    const f64 log_fp = std::log(fp);
                    const f64 log_fc = std::log((fc > 0.0) ? fc : eps * 1e-6);
                    const f64 log_eps = std::log(eps);

                    f64 t = 0.0;
                    if (log_fc != log_fp)
                        t = (log_eps - log_fp) / (log_fc - log_fp);

                    if (t < 0.0) t = 0.0;
                    else if (t > 1.0) t = 1.0;

                    smooth_newt = (f64)(newt_i - 1) + t;
                    t_cross = t;
                    have_cross = true;
                }
            }

            // --- smooth DIST across the same crossing ---
            if constexpr (NEED_DIST)
            {
                if (!prev_below && now_conv && newt_i > 0 && !have_dist_cross)
                {
                    const f64 fp = std::sqrt((f64)f2_prev);
                    const f64 fc = std::sqrt((f64)f2);

                    const f64 dfp_abs = std::sqrt((f64)cplx_mag2(dfz_prev));
                    const f64 dfc_abs = std::sqrt((f64)cplx_mag2(dfz));

                    const f64 d_prev = (dfp_abs == 0.0) ? 0.0 : fp / dfp_abs;
                    const f64 d_curr = (dfc_abs == 0.0) ? 0.0 : fc / dfc_abs;

                    f64 t = 0.0;
                    if (have_cross)
                    {
                        t = t_cross;
                    }
                    else
                    {
                        const f64 log_fp = std::log(fp);
                        const f64 log_fc = std::log((fc > 0.0) ? fc : eps * 1e-6);
                        const f64 log_eps = std::log(eps);

                        if (log_fc != log_fp)
                            t = (log_eps - log_fp) / (log_fc - log_fp);

                        if (t < 0.0) t = 0.0;
                        else if (t > 1.0) t = 1.0;
                    }

                    dist_smooth = (1.0 - t) * d_prev + t * d_curr;
                    have_dist_cross = true;
                }
            }

            // --- STRIPE accumulation at the evaluation point (matches your Newton kernel behaviour) ---
            if constexpr (NEED_STRIPES)
            {
                T zr2 = cplx_mag2(z);
                if (zr2 > zero)
                {
                    float invr = 1.0f / std::sqrt((f32)zr2);
                    float cc = (float)z.x * invr;
                    float ss = (float)z.y * invr;
                    stripe.accumulate(cc, ss);
                }
            }

            if (now_conv)
            {
                converged = true;
                break;
            }

            if (bad_derivative)
            {
                bad_newton = true;
                break;
            }

            // --- propagate derivative through Newton map (for escape DE) ---
            if constexpr (NEED_DIST)
            {
                // N'(z) = f(z) f''(z) / (f'(z))^2, with f''(z) = 6 z
                Complex<T> den = dfz * dfz;                 // (f')^2
                const T den2 = cplx_mag2(den);
                if (den2 == zero)
                {
                    bad_newton = true;
                    break;
                }

                Complex<T> num = fz * (T(6) * z);          // f * f''
                Complex<T> nd = num / den;                // N'(z)
                dz *= nd;
            }

            // Newton update
            Complex<T> step = fz / dfz;
            z -= step;

            ++newt_i;
            ++step_i;

            r2 = cplx_mag2(z);
            if (r2 > escape_r2)
            {
                escaped = true;
                break;
            }
        }
    }

    // ---------- final DIST ----------
    if constexpr (NEED_DIST)
    {
        if (escaped)
        {
            const T r = sqrt(r2);
            const T dz_abs = sqrt(cplx_mag2(dz));
            dist = (dz_abs == zero) ? 0.0 : (f64)(r * log(r) / dz_abs);
        }
        else if (converged)
        {
            if (have_dist_cross)
            {
                dist = dist_smooth;
            }
            else
            {
                const T f_abs = sqrt(f2);
                const T df_abs = sqrt(cplx_mag2(dfz));
                dist = (df_abs == zero) ? 0.0 : (f64)(f_abs / df_abs);
            }
        }
        else
        {
            dist = -1.0;
        }
    }

    // ---------- resolve "no result" ----------
    if (!escaped && !converged)
    {
        depth = INSIDE_MANDELBROT_SET;
        if constexpr (NEED_STRIPES) stripe.escape(0.0);
        return;
    }

    // ---------- final ITER depth ----------
    if (escaped)
    {
        if constexpr (NEED_SMOOTH_ITER)
        {
            const f64 log_abs_z = 0.5 * log_as_double(r2);
            const f64 nu = (f64)step_i + 1.0 - std::log2(log_abs_z);
            depth = nu - smooth_depth_offset<S>();
            if (depth < 0) depth = 0;
        }
        else
        {
            depth = (f64)step_i;
        }
    }
    else // converged
    {
        if constexpr (NEED_SMOOTH_ITER)
        {
            // Newton evaluations happen at odd step indices (1,3,5,...) in this schedule.
            // Convert smooth_newt (Newton-eval index space) -> step-space: step = 2*smooth_newt + 1
            f64 it_steps = have_cross ? (2.0 * smooth_newt + 1.0) : (f64)step_i;
            if (it_steps < 0.0) it_steps = 0.0;
            depth = it_steps;
        }
        else
        {
            depth = (f64)step_i;
        }
    }

    // ---------- final STRIPE encoding ----------
    if constexpr (NEED_STRIPES)
    {
        if (escaped)
        {
            stripe.escape((f64)r2);
        }
        else // converged
        {
            const f64 log_er2 = log_escape_radius2<S>();
            const f64 t = have_cross ? t_cross : 1.0;

            const f64 scale = std::pow(2.0, t - 1.0);
            const f64 lr = log_er2 / scale;

            stripe.log_r2_at_escape = lr;
            stripe.escaped = true;
        }
    }
}*/

template<class T, KernelFeatures S>
FORCE_INLINE void mandel_newton_alt_kernel(
    T x0, T y0,
    int iter_lim,
    f64& depth,
    f64& dist,
    StripeAccum& stripe,
    StripeParams sp = {})
{
    constexpr bool NEED_DIST = (bool)(S & KernelFeatures::DIST);
    constexpr bool NEED_SMOOTH_ITER = (bool)(S & KernelFeatures::ITER);
    constexpr bool NEED_STRIPES = (bool)(S & KernelFeatures::STRIPES);

    constexpr T escape_r2 = T(escape_radius2<S>());
    constexpr T zero = T(0), one = T(1);

    const T   eps2 = T(1e-12);
    const f64 eps = std::sqrt((f64)eps2);

    int step_i = 0;   // total alternating steps executed
    int newt_i = 0;   // number of completed Newton updates

    Complex<T> z{ zero, zero }; // parameter-plane: z0 = 0
    Complex<T> c{ x0, y0 };

    Complex<T> dz{ one, zero }; // derivative accumulator (for escape DE)
    T r2 = T(0);

    // Newton state for smoothing / distance
    Complex<T> fz{ zero, zero };
    Complex<T> dfz{ zero, zero };
    Complex<T> dfz_prev{ zero, zero };
    T f2 = T(0);
    T f2_prev = T(0);

    bool have_cross = false;
    f64  t_cross = 1.0;
    f64  smooth_newt = 0.0;

    bool have_dist_cross = false;
    f64  dist_smooth = 0.0;

    bool escaped = false;
    bool escaped_on_mandel = false;
    bool converged = false;
    bool bad_newton = false;

    //if constexpr (NEED_STRIPES)
    //    stripe.n = (int)sp.freq;

    while (true)
    {
        if (step_i >= iter_lim)
            break;

        //const bool do_mandel = ((step_i & 1) == 0);
        const bool do_mandel = (step_i % 3) == 0;// ((step_i & 1) == 0);

        if (do_mandel)
        {
            if constexpr (NEED_DIST)
                mandel_step_derivative(z, dz);

            mandel_step(z, c);
            ++step_i;

            r2 = cplx_mag2(z);

            if constexpr (NEED_STRIPES)
            {
                float invr = 1.0f / std::sqrt((f32)r2);
                float cc = (float)z.x * invr;
                float ss = (float)z.y * invr;
                stripe.accumulate(cc, ss);
            }

            if (r2 > escape_r2)
            {
                escaped = true;
                escaped_on_mandel = true;
                break;
            }
        }
        else
        {
            f2_prev = f2;
            dfz_prev = dfz;
            const bool prev_below = (f2_prev <= eps2);

            // f(z) = z^3 - 1,  f'(z) = 3 z^2,  f''(z) = 6 z
            Complex<T> z2 = z * z;
            Complex<T> z3 = z2 * z;

            fz = z3 - Complex<T>(one, zero);
            dfz = T(3) * z2;

            f2 = cplx_mag2(fz);
            const bool now_conv = (f2 <= eps2);

            const T df2 = cplx_mag2(dfz);
            const bool bad_derivative = (df2 == zero);

            if constexpr (NEED_SMOOTH_ITER)
            {
                if (!prev_below && now_conv && newt_i > 0 && !have_cross)
                {
                    const f64 fp = std::sqrt((f64)f2_prev);
                    const f64 fc = std::sqrt((f64)f2);

                    const f64 log_fp = std::log(fp);
                    const f64 log_fc = std::log((fc > 0.0) ? fc : eps * 1e-6);
                    const f64 log_eps = std::log(eps);

                    f64 t = 0.0;
                    if (log_fc != log_fp)
                        t = (log_eps - log_fp) / (log_fc - log_fp);

                    if (t < 0.0) t = 0.0;
                    else if (t > 1.0) t = 1.0;

                    smooth_newt = (f64)(newt_i - 1) + t;
                    t_cross = t;
                    have_cross = true;
                }
            }

            if constexpr (NEED_DIST)
            {
                if (!prev_below && now_conv && newt_i > 0 && !have_dist_cross)
                {
                    const f64 fp = std::sqrt((f64)f2_prev);
                    const f64 fc = std::sqrt((f64)f2);

                    const f64 dfp_abs = std::sqrt((f64)cplx_mag2(dfz_prev));
                    const f64 dfc_abs = std::sqrt((f64)cplx_mag2(dfz));

                    const f64 d_prev = (dfp_abs == 0.0) ? 0.0 : fp / dfp_abs;
                    const f64 d_curr = (dfc_abs == 0.0) ? 0.0 : fc / dfc_abs;

                    f64 t = 0.0;
                    if (have_cross)
                    {
                        t = t_cross;
                    }
                    else
                    {
                        const f64 log_fp = std::log(fp);
                        const f64 log_fc = std::log((fc > 0.0) ? fc : eps * 1e-6);
                        const f64 log_eps = std::log(eps);

                        if (log_fc != log_fp)
                            t = (log_eps - log_fp) / (log_fc - log_fp);

                        if (t < 0.0) t = 0.0;
                        else if (t > 1.0) t = 1.0;
                    }

                    dist_smooth = (1.0 - t) * d_prev + t * d_curr;
                    have_dist_cross = true;
                }
            }

            if constexpr (NEED_STRIPES)
            {
                T zr2 = cplx_mag2(z);
                if (zr2 > zero)
                {
                    float invr = 1.0f / std::sqrt((f32)zr2);
                    float cc = (float)z.x * invr;
                    float ss = (float)z.y * invr;
                    stripe.accumulate(cc, ss);
                }
            }

            if (now_conv)
            {
                converged = true;
                break;
            }

            if (bad_derivative)
            {
                bad_newton = true;
                break;
            }

            if constexpr (NEED_DIST)
            {
                // N'(z) = f(z) f''(z) / (f'(z))^2
                Complex<T> den = dfz * dfz;
                const T den2 = cplx_mag2(den);
                if (den2 == zero)
                {
                    bad_newton = true;
                    break;
                }

                Complex<T> num = fz * (T(6) * z);
                Complex<T> nd = num / den;
                dz *= nd;
            }

            Complex<T> step = fz / dfz;
            z -= step;

            ++newt_i;
            ++step_i;

            r2 = cplx_mag2(z);
            if (r2 > escape_r2)
            {
                escaped = true;
                escaped_on_mandel = false;
                break;
            }
        }
    }

    if constexpr (NEED_DIST)
    {
        if (escaped)
        {
            const T r = sqrt(r2);
            const T dz_abs = sqrt(cplx_mag2(dz));
            dist = (dz_abs == zero) ? 0.0 : (f64)(r * log(r) / dz_abs);
        }
        else if (converged)
        {
            if (have_dist_cross)
            {
                dist = dist_smooth;
            }
            else
            {
                const T f_abs = sqrt(f2);
                const T df_abs = sqrt(cplx_mag2(dfz));
                dist = (df_abs == zero) ? 0.0 : (f64)(f_abs / df_abs);
            }
        }
        else
        {
            dist = -1.0;
        }
    }

    if (!escaped && !converged)
    {
        depth = INSIDE_MANDELBROT_SET;
        if constexpr (NEED_STRIPES) stripe.escape(0.0);
        return;
    }

    if (escaped)
    {
        if constexpr (NEED_SMOOTH_ITER)
        {
            const f64 log_abs_z = 0.5 * log_as_double(r2);

            const int mandel_steps =
                escaped_on_mandel ? ((step_i + 1) >> 1) : (step_i >> 1);

            f64 nu_m = (f64)mandel_steps + 1.0 - std::log2(log_abs_z);
            nu_m -= smooth_depth_offset<S>();
            if (nu_m < 0.0) nu_m = 0.0;

            depth = escaped_on_mandel ? (2.0 * nu_m - 1.0) : (2.0 * nu_m);
            if (depth < 0.0) depth = 0.0;
        }
        else
        {
            depth = (f64)step_i;
        }
    }
    else
    {
        if constexpr (NEED_SMOOTH_ITER)
        {
            f64 it_steps = have_cross ? (2.0 * smooth_newt + 1.0) : (f64)step_i;
            if (it_steps < 0.0) it_steps = 0.0;
            depth = it_steps;
        }
        else
        {
            depth = (f64)step_i;
        }
    }

    if constexpr (NEED_STRIPES)
    {
        if (escaped)
        {
            stripe.escape((f64)r2);
        }
        else
        {
            const f64 log_er2 = log_escape_radius2<S>();
            const f64 t = have_cross ? t_cross : 1.0;

            const f64 scale = std::pow(2.0, t - 1.0);
            const f64 lr = log_er2 / scale;

            stripe.log_r2_at_escape = lr;
            stripe.escaped = true;
        }
    }
}

SIM_END;