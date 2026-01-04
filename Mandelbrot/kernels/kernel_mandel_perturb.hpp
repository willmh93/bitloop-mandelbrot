#pragma once

#include "ref_orbit.hpp"

SIM_BEG;

// Perturbation (no SIMD)
template<class T_lo, class T_hi, KernelFeatures F>
FORCE_INLINE bool mandel_kernel_perturb(
    const RefOrbitLo<T_lo>& ref,
    T_lo dcx, T_lo dcy,
    int iter_lim,
    f64& depth,
    f64& dist,
    StripeAccum& stripe)
{
    using V2 = simd2::v2<T_lo>;

    constexpr bool NEED_DIST    = (bool)(F & KernelFeatures::DIST);
    constexpr bool NEED_STRIPES = (bool)(F & KernelFeatures::STRIPES);
    constexpr T_lo  ER2 = T_lo(escape_radius2<F>());

    int ref_i = 0;
    const int max_ref = std::max(0, ref.max_ref);

    // --- delta state ---
    T_lo dx = 0, dy = 0;

    // --- derivative (for distance estimate) ---
    T_lo ddx = T_lo(1);
    T_lo ddy = T_lo(0);

    // --- z (previous step) ---
    T_lo zx_prev = 0, zy_prev = 0;

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

        // 2) advance delta
        const RefStep<T_lo>& s = ref.step[ref_i];

        const T_lo c2x = s.pre2x;
        const T_lo c2y = s.pre2y;

        const T_lo sx = c2x + dx;
        const T_lo sy = c2y + dy;

        const T_lo dx_sx = dx * sx;
        const T_lo dy_sy = dy * sy;
        const T_lo dx_sy = dx * sy;
        const T_lo dy_sx = dy * sx;

        dx = dx_sx - dy_sy + dcx;
        dy = dx_sy + dy_sx + dcy;

        // 3) build reference z_post
        const T_lo zpostx = s.postx;
        const T_lo zposty = s.posty;

        const T_lo zax = zpostx + dx;
        const T_lo zay = zposty + dy;

        const T_lo r2 = zax * zax + zay * zay;

        // 4) STRIPES sampling on z_{n+1}
        if constexpr (NEED_STRIPES) {

            float invr = 1.0f / std::sqrt((f32)r2);
            float c = (float)zax * invr; // cos
            float s = (float)zay * invr; // sin
            stripe.accumulate(c, s);
        }

        // 5) escape test
        if (r2 > ER2)
        {
            // smooth ITER
            const f64 r2d = (f64)r2;
            const f64 log_abs_z = 0.5 * std::log(std::max(r2d, 1e-300));
            const f64 nu = f64(n + 1) + 1.0 - std::log2(std::max(log_abs_z, 1e-300));
            depth = nu - smooth_depth_offset<F>();
            if (depth < 0) depth = 0;

            // DIST
            if constexpr (NEED_DIST)
            {
                const T_lo r = T_lo(std::sqrt((f64)r2));
                const T_lo dzabs = T_lo(std::hypot((f64)ddx, (f64)ddy));
                dist = (dzabs == 0) ? 0 : (r * std::log(r) / dzabs);
            }
            else 
            {
                dist = 0;
            }

            // STRIPE
            if constexpr (NEED_STRIPES)
                stripe.escape(r2d);

            return true;
        }


        // 6) rebase
        const T_lo d2 = dx * dx + dy * dy;
        bool must_rebase = false;

        // A) after using the last safe ref step, rebase before the next iteration
        if (ref_i == max_ref)
            must_rebase = true;

        // B) boundary-aware rebase
        if (!must_rebase)
        {
            constexpr T_lo ALPHA2 = (T_lo)(0.75 * 0.75);
            const T_lo lim2 = ALPHA2 * s.slack2;
            if (d2 > lim2)
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
            // re-center
            dx = zax;
            dy = zay;

            ref_i = 0;

            zx_prev = zax;
            zy_prev = zay;

            continue;
        }

        ++ref_i;
        zx_prev = zax;
        zy_prev = zay;
    }

    // Did not escape
    depth = INSIDE_MANDELBROT_SET;
    if constexpr (NEED_DIST)    dist = T_lo(-1);
    if constexpr (NEED_STRIPES) stripe.escape(0.0);
    return true;
}

SIM_END;