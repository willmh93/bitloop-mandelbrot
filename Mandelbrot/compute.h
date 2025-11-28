#pragma once

#include "shading.h"
#include "conversions.h"
#include <cmath>
#include <algorithm>

#include <bitloop.h>

SIM_BEG;

template<class T>
using cplx = Vec2<T>;

template<class T>
FORCE_INLINE constexpr void mandel_step(cplx<T>& z, const cplx<T>& c)
{
    T xx = z.x * z.x;
    T yy = z.y * z.y;
    T xy = (z.x * z.y);

    z.x = xx - yy + c.x;
    z.y = (xy + xy) + c.y;
}

template<class T>
FORCE_INLINE constexpr void mandel_step_derivative(const cplx<T>& z, cplx<T>& dz)
{
    const T zx_dzx = z.x * dz.x;
    const T zy_dzy = z.y * dz.y;
    const T zx_dzy = z.x * dz.y;
    const T zy_dzx = z.y * dz.x;

    dz.x = ((zx_dzx - zy_dzy) + (zx_dzx - zy_dzy)) + T(1);
    dz.y = (zx_dzy + zy_dzx) + (zx_dzy + zy_dzx);
}

template<class T>
FORCE_INLINE constexpr T cplx_mag2(const cplx<T>& z)
{
    return z.x * z.x + z.y * z.y;
}

template<typename T>
FORCE_INLINE bool interiorCheck(T x0, T y0)
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

template<class T>
FORCE_INLINE T clamp01(T v)
{
    return v < T(0) ? T(0) : (v > T(1) ? T(1) : v);
}

// helper for looking ahead multiple steps (experimental)
template<typename T>
inline cplx<T> process_z(int steps, cplx<T> z, const cplx<T>& c)
{
    for (int i = 0; i < steps; i++)
        mandel_step(z, c);

    return z;
}

// plain mandelbrot kernet for ITER + DIST + STRIPE
template<class T, KernelFeatures S>
FORCE_INLINE void mandel_kernel(
    T x0, T y0,
    int iter_lim,
    f64& depth,
    f64& dist,
    StripeAccum& stripe,
    StripeParams sp = {})
{
    constexpr bool NEED_DIST        = (bool)(S & KernelFeatures::DIST);
    constexpr bool NEED_SMOOTH_ITER = (bool)(S & KernelFeatures::ITER);
    constexpr bool NEED_STRIPES     = (bool)(S & KernelFeatures::STRIPES);

    constexpr T escape_r2 = T(escape_radius2<S>());
    constexpr T zero = T(0), one = T(1);

    int iter = 0;
    T r2 = T{ 0 };

    cplx<T> z{ zero, zero };
    cplx<T> c{ x0, y0 };
    cplx<T> dz{ one, zero };

    // stripe accumulators
    if constexpr (NEED_STRIPES)
        stripe.n = (int)sp.freq;  // fix n for this compute

    while (true)
    {
        if constexpr (NEED_DIST)
            mandel_step_derivative(z, dz);

        cplx<T> z0 = z;

        // step
        mandel_step(z, c);
        ++iter;

        // update radius^2
        r2 = cplx_mag2(z);

        if constexpr (NEED_STRIPES)
        {
            float invr = 1.0f / std::sqrt((f32)r2);
            float c = (float)z.x * invr; // cos
            float s = (float)z.y * invr; // sin
            stripe.accumulate(c, s);
        }
        if (r2 > escape_r2 || iter >= iter_lim) break;
    }

    const bool escaped = (r2 > escape_r2) && (iter < iter_lim);

    if constexpr (NEED_DIST)
    {
        if (escaped) {
            const T r = sqrt(r2);
            const T dz_abs = sqrt(cplx_mag2(dz));
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

template<class T>
struct alignas(16) RefStep {
    T pre2x, pre2y;
    T postx, posty;
    T slack2;      
    T _pad;
};

template<class T>
struct RefOrbit 
{
    std::vector<cplx<T>> z_pre;  // before step n+1
    std::vector<cplx<T>> z_post; // after  step n+1
    std::vector<T>       zpost_abs2;

    T cx{}, cy{};

    int  iter_esc = 0;
    bool escaped = false;
    int  max_ref = 0;
};

template<class T_lo>
struct RefOrbitLo 
{
    T_lo cx{}, cy{};
    int max_ref = 0;

    std::vector<RefStep<T_lo>> steps;
    const RefStep<T_lo>* step = nullptr;
};

template<class T, KernelFeatures S>
FORCE_INLINE void build_ref_orbit(const T& cx, const T& cy, int iter_lim, RefOrbit<T>& out)
{
    constexpr T ER2 = T(escape_radius2<S>());
    const T zero = T(0);

    out.cx = cx; 
    out.cy = cy;

    out.z_pre.clear(); 
    out.z_post.clear();
    out.zpost_abs2.clear();

    out.z_pre.reserve(iter_lim);
    out.z_post.reserve(iter_lim);
    out.zpost_abs2.reserve(iter_lim);

    out.iter_esc = iter_lim;
    out.escaped = false;

    cplx<T> z{ zero, zero };
    cplx<T> c{ cx, cy };

    for (int n = 0; n < iter_lim; ++n)
    {
        // pre
        out.z_pre.push_back(z);

        // step
        mandel_step(z, c);

        // post
        out.z_post.push_back(z);
        out.zpost_abs2.push_back(z.x * z.x + z.y * z.y);

        // check for escape
        if (!out.escaped && cplx_mag2(z) > ER2) 
        {
            out.iter_esc = n + 1;
            out.escaped = true;
            break;
        }
    }

    // last post-step index that is still safe to use
    out.max_ref = out.escaped ? (out.iter_esc - 1) : (int)out.z_post.size() - 1;
    if (out.max_ref < 0) out.max_ref = 0;
}

template<class Thi, class T_lo, KernelFeatures F>
inline void downcast_orbit(const RefOrbit<Thi>& hi, RefOrbitLo<T_lo>& lo)
{
    const size_t N = hi.z_pre.size();
    lo.steps.resize(N);

    const T_lo ER = (T_lo)escape_radius<F>();

    for (size_t i = 0; i < N; ++i) 
    {
        const T_lo zx = (T_lo)hi.z_pre[i].x;
        const T_lo zy = (T_lo)hi.z_pre[i].y;

        const T_lo postx = (T_lo)hi.z_post[i].x;
        const T_lo posty = (T_lo)hi.z_post[i].y;

        // precompute slack^2
        const T_lo zabs = (T_lo)std::sqrt((f64)hi.zpost_abs2[i]);
        const T_lo s = std::max<T_lo>(T_lo(0), ER - zabs);
        const T_lo slack2 = s * s;

        // fill packed step
        RefStep<T_lo>& st = lo.steps[i];
        st.pre2x = zx + zx;
        st.pre2y = zy + zy;
        st.postx = postx;
        st.posty = posty;
        st.slack2 = slack2;
        st._pad = T_lo(0);
    }

    lo.cx = (T_lo)hi.cx;
    lo.cy = (T_lo)hi.cy;
    lo.max_ref = hi.max_ref;

    // expose raw pointer for kernel
    lo.step = lo.steps.data();
}

// Perturbation (no SIMD)
template<class T_lo, class T_hi, KernelFeatures F>
FORCE_INLINE bool mandel_kernel_perturb(
    const RefOrbitLo<T_lo>& ref,
    T_lo dcx, T_lo dcy,
    int iter_lim,
    f64& depth,
    f64& dist,
    StripeAccum& stripe,
    StripeParams sp = {})
{
    using V2 = simd2::v2<T_lo>;

    constexpr bool NEED_DIST = (bool)(F & KernelFeatures::DIST);
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

    // --- stripe accumulators ---
    if constexpr (NEED_STRIPES)
        stripe.n = (int)sp.freq;  // fix n for this compute

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
            else {
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


// Perturbation (SIMD) - Not actively used, but kept for debugging/experimenting with new features
template<class T_lo, class T_hi, KernelFeatures F>
FORCE_INLINE bool mandel_kernel_perturb_simd(
    const RefOrbitLo<T_lo>& ref,
    T_lo dcx, T_lo dcy,
    int iter_lim,
    f64& depth,
    f64& dist,
    StripeAccum& stripe,
    StripeParams sp
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

    // --- stripe ---
    if constexpr (NEED_STRIPES) 
        stripe.n = (int)sp.freq;  // fix n for this compute

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
                stripe.escape((f64)r2);
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
    if constexpr (NEED_STRIPES) stripe.escape(0.0);

    return true;
}

// Perturbation (SIMD, unrolled)
template<class T_lo, class T_hi, KernelFeatures F>
FORCE_INLINE bool mandel_kernel_perturb_simd_unrolled(
    const RefOrbitLo<T_lo>& ref,
    T_lo dcx, T_lo dcy, int iter_lim,
    f64& depth, f64& dist, StripeAccum& stripe,
    StripeParams sp = {})
{
    using V2 = simd2::v2<T_lo>;

    constexpr bool NEED_DIST = (bool)(F & KernelFeatures::DIST);
    constexpr bool NEED_STRIPES = (bool)(F & KernelFeatures::STRIPES);
    constexpr T_lo ER2 = T_lo(escape_radius2<F>());
    const int max_ref = std::max(0, ref.max_ref);

    int ref_i = 0, n = 0;
    T_lo dx = 0, dy = 0;
    T_lo ddx = T_lo(1), ddy = T_lo(0);
    T_lo zx_prev = 0, zy_prev = 0;

    // --- stripes ---
    if constexpr (NEED_STRIPES)
        stripe.n = (int)sp.freq;

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
                float invr = 1.0f / std::sqrt((f32)r2);                                                 \
                float c = (float)zax * invr;                                                            \
                float s = (float)zay * invr;                                                            \
                stripe.accumulate(c, s);                                                                \
            }                                                                                           \
            if (r2 > ER2) {                                                                             \
                const f64 r2d = (f64)r2;                                                                \
                const f64 log_abs_z = 0.5 * std::log(std::max(r2d, 1e-300));                            \
                const f64 nu = f64(n + 1) + 1.0 - std::log2(std::max(log_abs_z, 1e-300));               \
                depth = nu - smooth_depth_offset<F>();                                                  \
                if (depth < 0) depth = 0;                                                               \
                                                                                                        \
                if constexpr (NEED_DIST) {                                                              \
                    const T_lo r = T_lo(std::sqrt((f64)r2));                                            \
                    const T_lo dzabs = T_lo(std::hypot((f64)ddx, (f64)ddy));                            \
                    dist = (dzabs == 0) ? 0 : (r * std::log(r) / dzabs);                                \
                }                                                                                       \
                else dist = 0;                                                                          \
                                                                                                        \
                if constexpr (NEED_STRIPES)                                                             \
                    stripe.escape(r2d);                                                                 \
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
    if constexpr (NEED_STRIPES) stripe.escape(0.0);
    return true;
}

/// ----- calculate raw [iter, dist, stripe] values for each field pixel -----
template<typename T, KernelFeatures F, KernelMode K>
bool Mandelbrot_Scene::compute_mandelbrot(int timeout)
{
    typedef f64 T_lo;

    int threads = Thread::threadCount();
    int compute_tiles = threads;

    f32 tiles_sqrt = ceil(sqrt((f32)compute_tiles));
    int tile_w = (int)ceil((f32)pending_bmp->width() / tiles_sqrt);
    int tile_h = (int)ceil((f32)pending_bmp->height() / tiles_sqrt);
    bool frame_complete = false;

    if constexpr ((bool)(K & KernelMode::PERTURBATION_MASK))
    {
        // get high precision anchor world pos
        RefOrbitLo<T_lo> orbit_lo;
        Vec2<T> anchor = camera.pos<T>();

        // calculate reference orbit in high precision, downcast to low
        {
            RefOrbit<T> orbit_hi;
            build_ref_orbit<T, F>(anchor.x, anchor.y, iter_lim, orbit_hi);
            downcast_orbit<T, T_lo, F>(orbit_hi, orbit_lo);
        }

        // run perturbation kernel for remaining tiles pixels
        frame_complete = pending_bmp->forEachWorldTilePixel<T>(tile_w, tile_h, P, [&](int x, int y, T wx, T wy)
        {
            EscapeFieldPixel& field_pixel = pending_field->at(x, y);

            if (interiorCheck(wx, wy)) {
                field_pixel.depth = INSIDE_MANDELBROT_SET_SKIPPED;
                return;
            }

            f64 depth = field_pixel.depth;
            if (depth >= 0) return;

            const T_lo dcx_lo = (T_lo)(wx - anchor.x);
            const T_lo dcy_lo = (T_lo)(wy - anchor.y);

            f64 dist;
            StripeAccum stripe;

            #if MANDEL_EXTENDED_FIELD_STATS
            ExtendedFieldStats extended_stats;
            #define MANDEL_ENXTENDED_ARG ,extended_stats
            #else
            #define MANDEL_ENXTENDED_ARG
            #endif

            if constexpr (K == KernelMode::PERTURBATION_SIMD_UNROLLED)
                mandel_kernel_perturb_simd_unrolled<T_lo, T, F>(orbit_lo, dcx_lo, dcy_lo, iter_lim, depth, dist, stripe, stripe_params);
            else if constexpr (K == KernelMode::PERTURBATION_SIMD)
                mandel_kernel_perturb_simd<T_lo, T, F>(orbit_lo, dcx_lo, dcy_lo, iter_lim, depth, dist, stripe, stripe_params MANDEL_ENXTENDED_ARG);
            else
                mandel_kernel_perturb<T_lo, T, F>(orbit_lo, dcx_lo, dcy_lo, iter_lim, depth, dist, stripe, stripe_params);


            field_pixel.set(depth, dist, stripe);

            #if MANDEL_EXTENDED_FIELD_STATS
            field_pixel.extended_stats = extended_stats;
            #undef MANDEL_ENXTENDED_ARG
            #endif
        }, compute_tiles, timeout);

        // grab most up-to-date field
        EscapeField* complete_field = active_field;
        if (frame_complete)
            complete_field = pending_field; // prefer high-res field
        else if (!complete_field)
            blBreak();                      // shouldn't reach here - phase 0 is always one-shot

        // run perturbation kernel on normalization points (reusing existing "currently active" field first where acceptable)
        double bmp_scale = phaseBmpScale(complete_field->compute_phase);
        double fdp_tolerance = 9.0 * (1.0 - normalize_field_precision);
        double fdp_tolerance2 = fdp_tolerance * fdp_tolerance;

        std::atomic<u64> norm_full_computes(0);
        norm_field.forEach([&](NormalizationPixel& field_pixel)
        {
            // already done 100% precise manual calculation on previous frame? skip
            if (field_pixel.is_final) return;

            DVec2 fp = field_pixel.stage_pos / bmp_scale;
            if (fp.x >= 0.0 && fp.y >= 0.0 && fp.x < complete_field->w && fp.y < complete_field->h)
            {
                IVec2 p = (IVec2)fp;
                double fdpf = (fp - ((DVec2)p)).mag2(); // delta dist (pixel fraction)
                double fdp = (fdpf * bmp_scale);

                if (fdp < fdp_tolerance2)
                {
                    EscapeFieldPixel& existing_pixel = complete_field->at(p.x, p.y);
                    field_pixel.set(existing_pixel);
                    return;
                }
            }

            // coordinate lies outside of viwport rect, do calculation
            const Vec2<T> world_pos = field_pixel.worldPos<T>();
            const T_lo dcx_lo = (T_lo)(world_pos.x - anchor.x);
            const T_lo dcy_lo = (T_lo)(world_pos.y - anchor.y);

            f64 depth, dist; StripeAccum stripe;


            #if MANDEL_EXTENDED_FIELD_STATS
            ExtendedFieldStats extended_stats; // placeholder
            #define MANDEL_ENXTENDED_ARG ,extended_stats
            #else
            #define MANDEL_ENXTENDED_ARG
            #endif

            norm_full_computes.fetch_add(1);

            if constexpr (K == KernelMode::PERTURBATION_SIMD_UNROLLED)
                mandel_kernel_perturb_simd_unrolled<T_lo, T, F>(orbit_lo, dcx_lo, dcy_lo, iter_lim, depth, dist, stripe, stripe_params);
            else if constexpr (K == KernelMode::PERTURBATION_SIMD)
                mandel_kernel_perturb_simd<T_lo, T, F>(orbit_lo, dcx_lo, dcy_lo, iter_lim, depth, dist, stripe, stripe_params MANDEL_ENXTENDED_ARG);
            else
                mandel_kernel_perturb<T_lo, T, F>(orbit_lo, dcx_lo, dcy_lo, iter_lim, depth, dist, stripe, stripe_params);

            field_pixel.set(depth, dist, stripe);

            // don't recalculate
            field_pixel.is_final = true;

            #if MANDEL_EXTENDED_FIELD_STATS
            #undef MANDEL_ENXTENDED_ARG
            #endif
        });
    }
    else // standard kernel (no perturbation)
    {
        frame_complete = pending_bmp->forEachWorldTilePixel<T>(tile_w, tile_h, P, [&](int x, int y, T wx, T wy)
        {
            EscapeFieldPixel& field_pixel = pending_field->at(x, y);

            if (interiorCheck(wx, wy)) {
                field_pixel.depth = INSIDE_MANDELBROT_SET_SKIPPED;
                return;
            }
    
            f64 depth = field_pixel.depth;
            if (depth >= 0) return;
    
            f64 dist; StripeAccum stripe;
            mandel_kernel<T, F>(wx, wy, iter_lim, depth, dist, stripe, stripe_params);
            field_pixel.set(depth, dist, stripe);
    
        }, compute_tiles, timeout);
    
        norm_field.forEach([&](NormalizationPixel& field_pixel)
        {
            // already calculated on previous frame? skip
            if (field_pixel.is_final) return;
    
            IVec2 p = field_pixel.stage_pos / phaseBmpScale(pending_field->compute_phase);
            EscapeFieldPixel* existing_pixel = pending_field->get(p.x, p.y);
    
            if (existing_pixel) {
             field_pixel.set(*existing_pixel);
             return;
         }
    
            // coordinate lies outside of viwport rect, do calculation
            f64 depth, dist; StripeAccum stripe;
            const Vec2<T> world_pos = field_pixel.worldPos<T>();
            mandel_kernel<T, F>(world_pos.x, world_pos.y, iter_lim, depth, dist, stripe, stripe_params);
            field_pixel.set(depth, dist, stripe);
    
            // don't recalculate
            field_pixel.is_final = true;
        });
    }

    return frame_complete;
}


template<bool Normalize_Depth>
FORCE_INLINE f64 toNormalizedDepth(f64 depth, f64 min_depth, f64 log1p_weight)
{
    if constexpr (Normalize_Depth)
    {
        f64 depth_from_floor = depth - min_depth;
        if (depth_from_floor < 0.0) depth_from_floor = 0.0;
        return Math::linear_log1p_lerp(depth_from_floor, log1p_weight);
    }
    else
        return Math::linear_log1p_lerp(depth, log1p_weight);
}

// returns signed log-normalized distance (depending on inversion)
template<bool Invert>
FORCE_INLINE f64 toNormalizedDist(f64 dist, f64 stable_min_dist, f64 stable_max_dist)
{
    //assert(dist >= 0);
    f64 rescaled_dist = log(dist);
    return Math::lerpFactor(rescaled_dist, stable_min_dist, stable_max_dist);
}

template<KernelFeatures F>
FORCE_INLINE f32 toNormalizedStripe(const StripeAccum& stripe, f32 phase)
{
    const double LOG_ER2 = log_escape_radius2<F>();
    const float cphi = std::cos(phase);
    const float sphi = std::sin(phase);
    return stripeFromAccum(stripe, LOG_ER2, cphi, sphi);
}

/// ----- determine [low, high, mean] of each feature, and normalized wrapping limits -----
template<typename T, KernelFeatures F, bool Normalize_Depth, bool Invert_Dist>
void Mandelbrot_Scene::calculate_normalize_info()
{
    // todo: move min/max/mean props to normalization field? It's field-independant
    // Calculate normalized depth/dist
    f32 zf = std::max(1.0f, (f32)bl::log10(camera.relativeZoom<f128>()));
    f64 sharpness_ratio = ((100.0 - dist_params.cycle_dist_sharpness) / 100.0 + 0.00001);
    f128 r1 = f128(0.1) / camera.relativeZoom<f128>();
    f128 r2 = f128(10.0) / camera.relativeZoom<f128>();
    f128 stable_min_raw_dist = r1 * sharpness_ratio;
    f128 stable_max_raw_dist = r2;
    f64  stable_min_dist = (f64)log(stable_min_raw_dist);
    f64  stable_max_dist = (f64)log(stable_max_raw_dist);

    constexpr size_t BINS = 256;
    constexpr u32 W_SCALE = 1u << 20; // ~1e6 precision

    struct WeightHist 
    { 
        // shared weights
        u64 total_w;
        std::array<u64, BINS> hist_w; // stripe histogram

        // feature * weight sums
        f64 sum_wstripe = 0.0;
    };

    std::vector<WeightHist> buckets(Thread::threadCount());

    // first-pass: establish ITER limits (todo: use histogram here instead of hard limits)
    active_field->min_depth = std::numeric_limits<f64>::max();
    active_field->max_depth = std::numeric_limits<f64>::lowest();

    norm_field.forEach([&](NormalizationPixel& field_pixel) 
    {
        f64 depth = field_pixel.depth;
        if (depth < INSIDE_MANDELBROT_SET_SKIPPED)
        {
            if (depth < active_field->min_depth) active_field->min_depth = depth;
            if (depth > active_field->max_depth) active_field->max_depth = depth;
        }
    });


    Thread::forEachBatch(norm_field.world_field, [&](std::span<NormalizationPixel> batch, int ti)
    {
        auto& bucket = buckets[ti];
        u64&  total_w = bucket.total_w;
        auto& hist_w  = bucket.hist_w;

        f64& sum_wstripe = bucket.sum_wstripe;

        for (NormalizationPixel& field_pixel : batch)
        {
            f64 depth = field_pixel.depth;

            if (depth < INSIDE_MANDELBROT_SET_SKIPPED)
            {
                // convert weight to fixed-point
                u64 w = (u64)std::llround((f64)field_pixel.weight * (f64)W_SCALE);
                total_w += w;

                f32 normalized_stripe = toNormalizedStripe<F>(field_pixel.stripe, stripe_params.phase);

                size_t bin = (size_t)std::floor(normalized_stripe * (BINS - 1));
                bin = std::clamp(bin, (size_t)0, BINS - 1);
                hist_w[bin] += w;

                // sum up weighted features
                sum_wstripe += w * normalized_stripe;
            }
        }
    });

    std::array<u64, BINS> hist_w{};
    u64 total_w = 0;
    f64 total_wstripe = 0;

    for (auto& bucket : buckets)
    {
        std::array<u64, BINS>& bucket_hist_w = bucket.hist_w;

        for (int bin = 0; bin < BINS; bin++)
            hist_w[bin] += bucket_hist_w[bin];

        total_w       += bucket.total_w;
        total_wstripe += bucket.sum_wstripe;
    }

    f32 stripe_mean = (f32)(total_wstripe / (f64)total_w);
    //f32 stripe_mag = 0.5f / zf;
    f32 stripe_mag = 0.1f / zf;

    active_field->stable_min_dist = stable_min_dist;
    active_field->stable_max_dist = stable_max_dist;
    active_field->mean_stripe = stripe_mean;

    dist_tone   .setParams(dist_tone_params);
    stripe_tone .setParams(stripe_tone_params);
    dist_tone   .configureFieldInfo(0.0f, 1.0f, 0.0f);
    stripe_tone .configureFieldInfo(-stripe_mag, stripe_mag, 0.0f);

    if (active_field->min_depth == std::numeric_limits<f64>::max()) 
        active_field->min_depth = 0;

    if (iter_params.cycle_iter_dynamic_limit)
    {
        /// "color_cycle_iters" represents ratio of (assumed) iter_lim
        f64 assumed_iter_lim = mandelbrotIterLimit(camera.relativeZoom<f128>()) * 0.5;
        f64 color_cycle_iters = iter_params.cycle_iter_value * assumed_iter_lim;

        active_field->log_color_cycle_iters = Math::linear_log1p_lerp(color_cycle_iters, iter_params.cycle_iter_log1p_weight);
    }
    else
    {
        /// "cycle_iter_value" represents actual iter_lim
        active_field->log_color_cycle_iters = Math::linear_log1p_lerp(iter_params.cycle_iter_value, iter_params.cycle_iter_log1p_weight);
    }

    active_field->cycle_dist_value = dist_params.cycle_dist_value;
}


/// ----- calculates final_iter, final_dist, final_stripe -----
template<typename T, KernelFeatures F, bool Normalize_Depth, bool Invert_Dist>
void Mandelbrot_Scene::normalize_field()
{
    f64 stable_min_dist = active_field->stable_min_dist;
    f64 stable_max_dist = active_field->stable_max_dist;

    f64 cycle_iters = active_field->log_color_cycle_iters;
    f32 cycle_dist = (f32)active_field->cycle_dist_value;

    //f32 zf = std::max(1.0f, (f32)bl::log10(camera.relativeZoom<f128>()));

    active_bmp->forEachPixel([&](int x, int y)
    {
        EscapeFieldPixel& field_pixel = active_field->at(x, y);
        f64 depth = field_pixel.depth;
        if (depth >= INSIDE_MANDELBROT_SET_SKIPPED || active_field->get_skip_flag(x, y)) return;

        f64 final_depth{};
        f64 final_dist{};
        f32 final_stripe{};

        /// ====== ITER =======
        if constexpr ((bool)(F & KernelFeatures::ITER))
        {
            final_depth = toNormalizedDepth<Normalize_Depth>(depth, active_field->min_depth, iter_params.cycle_iter_log1p_weight);
            final_depth /= cycle_iters;
        }

        /// ====== DIST =======
        if constexpr ((bool)(F & KernelFeatures::DIST))
        {
            final_dist = toNormalizedDist<Invert_Dist>(field_pixel.dist, stable_min_dist, stable_max_dist);
            final_dist = dist_tone.applyAndDenormalize(final_dist);
            final_dist *= Invert_Dist ? -1.0 : 1.0;
            final_dist /= cycle_dist;
        }

        /// ====== STRIPE =======
        if constexpr ((bool)(F & KernelFeatures::STRIPES))
        {
            float s = toNormalizedStripe<F>(field_pixel.stripe, stripe_params.phase);
            s = (s - active_field->mean_stripe);// *zf;
            final_stripe = stripe_tone.apply(s) - 0.5f;
        }

        field_pixel.final_depth  = (f32)final_depth;
        field_pixel.final_dist   = (f32)final_dist;
        field_pixel.final_stripe = final_stripe;
    });
}

/// ----- shades pixels based on final_iter, final_dist, final_stripe -----
template<MandelShaderFormula F, bool Highlight_Optimized_Interior>
void Mandelbrot_Scene::shadeBitmap()
{
    f64 iter_ratio, dist_ratio, stripe_ratio;
    shadingRatios(
        iter_weight, dist_weight, stripe_weight,
        iter_ratio, dist_ratio, stripe_ratio
    );

    active_bmp->forEachPixel([&](int x, int y)
    {
        EscapeFieldPixel& field_pixel = active_field->at(x, y);

        if (field_pixel.depth >= INSIDE_MANDELBROT_SET_SKIPPED)
        {
            if constexpr (Highlight_Optimized_Interior)
            {
                if (field_pixel.depth == INSIDE_MANDELBROT_SET_SKIPPED)
                    active_bmp->setPixel(x, y, 0xFF7F007F);
                else
                    active_bmp->setPixel(x, y, 0xFF000000);
            }
            else
            {
                active_bmp->setPixel(x, y, 0xFF000000);
            }
            return;
        }

        uint32_t u32;

        f32 iter_v = field_pixel.final_depth;
        f32 dist_v = field_pixel.final_dist;
        f32 stripe_v = field_pixel.final_stripe;

        f32 w_iter_v   = (iter_v   * (f32)iter_ratio);
        f32 w_dist_v   = (dist_v   * (f32)dist_ratio);
        f32 w_stripe_v = (stripe_v * (f32)stripe_ratio);

        f32 combined_t;

        if constexpr (F == MandelShaderFormula::ITER_DIST_STRIPE)
        {
            //combined_t = std::clamp(w_iter_v + w_dist_v + w_stripe_v, 0.0f, 1.0f);
            combined_t = Math::wrap01(w_iter_v + w_dist_v + w_stripe_v);
        }
        else if constexpr (F == MandelShaderFormula::ITER_DIST__MUL__STRIPE)
        {
            combined_t = Math::wrap01((w_iter_v + w_dist_v) * w_stripe_v);
        }
        else if constexpr (F == MandelShaderFormula::ITER__MUL__DIST_STRIPE)
        {
            combined_t = Math::wrap01(w_iter_v * (w_dist_v + w_stripe_v));
        }
        else if constexpr (F == MandelShaderFormula::ITER_STRIPE__MULT__DIST)
        {
            combined_t = Math::wrap01((w_iter_v + w_stripe_v) * w_dist_v);
        }

        if (isfinite(combined_t))
        {
            gradient_shifted.unguardedRGBA(combined_t, u32);
            active_bmp->setPixel(x, y, u32);
        }
    });
}

SIM_END;
