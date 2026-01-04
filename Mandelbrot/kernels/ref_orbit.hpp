#pragma once

SIM_BEG;

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
    std::vector<Complex<T>> z_pre;  // before step n+1
    std::vector<Complex<T>> z_post; // after  step n+1
    std::vector<T>          zpost_abs2;

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

template<class T, KernelFeatures F>
FORCE_INLINE void build_ref_orbit(const T& cx, const T& cy, int iter_lim, RefOrbit<T>& out)
{
    constexpr T ER2 = T(escape_radius2<F>());
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

    Complex<T> z{ zero, zero };
    Complex<T> c{ cx, cy };

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
        if (!out.escaped && z.mag2() > ER2)
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

SIM_END;