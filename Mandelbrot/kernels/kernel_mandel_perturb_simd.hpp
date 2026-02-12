#pragma once

#include "ref_orbit.hpp"

SIM_BEG;

// Perturbation (SIMD) - Not actively used, but kept for debugging/experimenting/timing
// since it's easier to work with than unrolled version (which is actively used)

template<class T_lo, KernelFeatures F>
FORCE_INLINE bool mandel_kernel_perturb_simd(
    const RefOrbitLo<T_lo>& ref,
    T_lo dcx, T_lo dcy,
    int iter_lim,
    f64& depth,
    f64& dist,
    StripeAccum& stripe
    #if MANDEL_EXTENDED_FIELD_STATS
    , ExtendedFieldStats& extended_stats
    #endif
)
{
    using V2 = simd2::v2<T_lo>;

    constexpr bool NEED_DIST = (bool)(F & KernelFeatures::DIST);
    constexpr bool NEED_STRIPES = (bool)(F & KernelFeatures::STRIPES);
    constexpr T_lo ER2 = T_lo(escape_radius2<F>());

    int ref_i = 0;
    const int max_ref = std::max(0, ref.max_ref);

    // --- delta state ---
    V2 v_d = V2::set(0, 0);

    // --- derivative (for distance estimate) ---
    T_lo ddx = T_lo(1);
    T_lo ddy = T_lo(0);

    // --- z (previous step) ---
    T_lo zx_prev = 0, zy_prev = 0;

    #if defined(__wasm_simd128__)
    #pragma clang loop unroll_count(2)
    #endif
    for (int n = 0; n < iter_lim; ++n)
    {
        //timer_sample(ITER);
        //timer_sample(DIST);
        //timer_sample(STRIPE);
        //timer_sample(ESCAPE);
        //timer_sample(REBASE);

        ///timer_t0(NORM_FIELD);
        ///timer_t1(NORM_FIELD);

        #if MANDEL_EXTENDED_FIELD_STATS
        extended_stats.iter_count++;
        #endif

        // 1) update derivative
        {
            ///timer_t0(DIST);
            if constexpr (NEED_DIST)
            {
                auto v_z = V2::set(zx_prev, zy_prev);
                auto v_dd = V2::set(ddx, ddy);
                auto v_twoz = V2::mul(V2::set(T_lo(2), T_lo(2)), v_z);
                auto v_t = V2::add(V2::cmul(v_twoz, v_dd), V2::set(T_lo(1), T_lo(0)));
                ddx = v_t.x();
                ddy = v_t.y();
            }
            ///timer_t1(DIST);
        }

        ///timer_t0(ITER);
        const RefStep<T_lo>& s = ref.step[ref_i];

        auto v_pre2 = V2::load(&s.pre2x);
        auto v_post = V2::load(&s.postx);

        // 2) advance delta
        auto v_s = V2::add(v_pre2, v_d);
        v_d = V2::add(V2::cmul(v_d, v_s), V2::set(dcx, dcy));

        // 3) build ref
        auto v_zn1 = V2::add(v_post, v_d);
        const T_lo zax = v_zn1.x();
        const T_lo zay = v_zn1.y();
        const T_lo r2 = V2::dot(v_zn1);
        ///timer_t1(ITER);

        // 4) STRIPES: accumulate unit direction for circular mean
        {
            timer_t0(STRIPE);
            if constexpr (NEED_STRIPES)
            {
                float invr = 1.0f / std::sqrt((f32)r2);
                float c = (float)zax * invr; // cos
                float s = (float)zay * invr; // sin
                stripe.accumulate(c, s);
            }
            timer_t1(STRIPE);
        }

        // 5) escape test
        if (r2 > ER2)
        {
            ///timer_t0(ESCAPE);

            #if MANDEL_EXTENDED_FIELD_STATS
            extended_stats.escape_count++;
            #endif

            // smooth ITER
            const f64 r2d = (f64)r2;
            const f64 log_abs_z = 0.5 * std::log(std::max(r2d, 1e-300));
            const f64 nu = f64(n + 1) + 1.0 - std::log2(std::max(log_abs_z, 1e-300));
            depth = nu - smooth_depth_offset<F>();
            if (depth < 0) depth = 0;

            // DIST
            {
                if constexpr (NEED_DIST)
                {
                    const T_lo r = T_lo(std::sqrt((f64)r2));
                    const T_lo dzabs = T_lo(std::hypot((f64)ddx, (f64)ddy));
                    dist = (dzabs == 0) ? 0 : (r * std::log(r) / dzabs);
                }

                if constexpr (!NEED_DIST)
                    dist = 0;
            }

            // STRIPE
            if constexpr (NEED_STRIPES)
                stripe.escape((f32)r2);
            ///timer_t1(ESCAPE);

            //timer_t1(NORM_FIELD)
            return true;
        }

        ///timer_t0(REBASE);
        // 6) rebase checks
        const T_lo d2 = V2::dot(v_d);
        bool must_rebase = (ref_i == max_ref);

        constexpr T_lo ALPHA2 = (T_lo)(0.75 * 0.75);
        constexpr T_lo H = (T_lo)0.95;

        if (!must_rebase)
        {
            const T_lo lim2 = ALPHA2 * s.slack2;
            if (lim2 > (T_lo)0 && d2 > lim2)
            {
                must_rebase = true;

                #if MANDEL_EXTENDED_FIELD_STATS
                extended_stats.slack_rebases++;
                #endif
            }
        }

        if (!must_rebase)
        {
            if (r2 < H * d2)
            {
                must_rebase = true;

                #if MANDEL_EXTENDED_FIELD_STATS
                extended_stats.dominance_rebases++;
                #endif
            }
        }

        if (must_rebase)
        {
            // re-center
            v_d = V2::set(zax, zay);
            ref_i = 0;

            zx_prev = zax;
            zy_prev = zay;

            //timer_t1(NORM_FIELD);

            ///timer_t1(REBASE);
            continue;
        }
        ///timer_t1(REBASE);

        ++ref_i;
        zx_prev = zax;
        zy_prev = zay;

        //timer_t1(NORM_FIELD);
    }

    // Did not escape
    depth = INSIDE_MANDELBROT_SET;

    if constexpr (NEED_DIST) dist = T_lo(-1);
    if constexpr (NEED_STRIPES) stripe.escape(0.0f);

    return true;
}

SIM_END;