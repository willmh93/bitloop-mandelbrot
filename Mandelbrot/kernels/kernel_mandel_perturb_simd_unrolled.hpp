#pragma once

#include "ref_orbit.hpp"

SIM_BEG;

// Perturbation (SIMD, unrolled)
template<class T_lo, KernelFeatures F>
FORCE_INLINE bool mandel_kernel_perturb_simd_unrolled(
    const RefOrbitLo<T_lo>& ref,
    T_lo dcx, T_lo dcy, int iter_lim,
    f64& depth, f64& dist, StripeAccum& stripe)
{
    using V2 = simd2::v2<T_lo>;

    constexpr bool NEED_DIST    = (bool)(F & KernelFeatures::DIST);
    constexpr bool NEED_STRIPES = (bool)(F & KernelFeatures::STRIPES);
    constexpr T_lo ER2 = T_lo(escape_radius2<F>());
    constexpr T_lo MIN_FLT = (T_lo)1.17549e-38;
    const int max_ref = std::max(0, ref.max_ref);

    int ref_i = 0, n = 0;
    T_lo dx = 0, dy = 0;
    T_lo ddx = T_lo(1), ddy = T_lo(0);
    T_lo zx_prev = 0, zy_prev = 0;

    // unroll manually
    #if defined(__clang__)
    #pragma clang loop unroll(disable)
    #endif

    while (n < iter_lim)
    {
        #define MANDEL_KERNEL_UNROLLED_BODY                                                             \
        {                                                                                               \
            if constexpr (NEED_DIST) {                                                                  \
                auto v_z = V2::set(zx_prev, zy_prev);                                                   \
                auto v_dd = V2::set(ddx, ddy);                                                          \
                auto v_twoz = V2::mul(V2::set(T_lo(2), T_lo(2)), v_z);                                  \
                auto v_t = V2::add(V2::cmul(v_twoz, v_dd), V2::set(T_lo(1), T_lo(0)));                  \
                ddx = v_t.x(); ddy = v_t.y();                                                           \
            }                                                                                           \
                                                                                                        \
            const RefStep<T_lo>& s = ref.step[ref_i];                                                   \
            auto v_pre2 = V2::load(&s.pre2x);                                                           \
            auto v_post = V2::load(&s.postx);                                                           \
            auto v_dx = V2::set(dx, dy);                                                                \
            auto v_s = V2::add(v_pre2, v_dx);                                                           \
            v_dx = V2::add(V2::cmul(v_dx, v_s), V2::set(dcx, dcy));                                     \
            dx = v_dx.x(); dy = v_dx.y();                                                               \
                                                                                                        \
            auto v_zn1 = V2::add(v_post, V2::set(dx, dy));                                              \
            const T_lo r2 = V2::dot(v_zn1);                                                             \
            const T_lo zax = v_zn1.x(), zay = v_zn1.y();                                                \
                                                                                                        \
            if constexpr (NEED_STRIPES) {                                                               \
                f32 invr = 1.0f / std::sqrt((f32)r2); /* sqrt f32 faster? */                            \
                f32 c = (f32)zax * invr;                                                                \
                f32 s = (f32)zay * invr;                                                                \
                stripe.accumulate(c, s);                                                                \
            }                                                                                           \
            if (r2 > ER2) {                                                                             \
                const T_lo log_abs_z = T_lo(0.5) * std::log(std::max(r2, MIN_FLT));                     \
                const T_lo nu = T_lo(n + 2) - std::log2(std::max(log_abs_z, MIN_FLT));                  \
                depth = (f64)(nu - smooth_depth_offset_d<F>());                                         \
                if (depth < 0) depth = 0;                                                               \
                                                                                                        \
                if constexpr (NEED_DIST) {                                                              \
                    const T_lo r = std::sqrt(r2);                                                       \
                    const T_lo dzabs = T_lo(std::hypot(ddx, ddy));                                      \
                    dist = (f64)((dzabs == 0) ? 0 : (r * std::log(r) / dzabs));                         \
                }                                                                                       \
                else dist = 0;                                                                          \
                                                                                                        \
                if constexpr (NEED_STRIPES)                                                             \
                    stripe.escape(r2);                                                                  \
                                                                                                        \
                return true;                                                                            \
            }                                                                                           \
                                                                                                        \
            const T_lo d2 = V2::dot(V2::set(dx, dy));                                                   \
            bool must_rebase = (ref_i == max_ref);                                                      \
                                                                                                        \
            if (!must_rebase) {                                                                         \
                constexpr T_lo ALPHA2 = (T_lo)(0.75 * 0.75);                                            \
                const T_lo lim2 = ALPHA2 * s.slack2;                                                    \
                if (lim2 > (T_lo)0 && d2 > lim2) must_rebase = true;                                    \
            }                                                                                           \
            if (!must_rebase) {                                                                         \
                constexpr T_lo H = (T_lo)0.95;                                                          \
                if (r2 < H * d2) must_rebase = true;                                                    \
            }                                                                                           \
            if (must_rebase) {                                                                          \
                dx = zax; dy = zay;                                                                     \
                ++n; ref_i = 0;                                                                         \
                zx_prev = zax; zy_prev = zay;                                                           \
                continue;                                                                               \
            }                                                                                           \
                                                                                                        \
            zx_prev = zax; zy_prev = zay;                                                               \
            ++n; ++ref_i;                                                                               \
        }

        MANDEL_KERNEL_UNROLLED_BODY; if (n >= iter_lim) break;
        MANDEL_KERNEL_UNROLLED_BODY; if (n >= iter_lim) break;
        MANDEL_KERNEL_UNROLLED_BODY; if (n >= iter_lim) break;
        MANDEL_KERNEL_UNROLLED_BODY; if (n >= iter_lim) break;
        MANDEL_KERNEL_UNROLLED_BODY; if (n >= iter_lim) break;
        MANDEL_KERNEL_UNROLLED_BODY; if (n >= iter_lim) break;
        MANDEL_KERNEL_UNROLLED_BODY; if (n >= iter_lim) break;
        MANDEL_KERNEL_UNROLLED_BODY;
    }

    // Did not escape
    depth = INSIDE_MANDELBROT_SET;
    if constexpr (NEED_DIST)    dist = T_lo(-1);
    if constexpr (NEED_STRIPES) stripe.no_escape();
    return true;
}

SIM_END;