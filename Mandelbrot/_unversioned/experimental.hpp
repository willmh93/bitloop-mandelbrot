#if MANDEL_UNVERSIONED_EXPERIMENTAL
#include "../build_config.h"
#include <bitloop.h>

SIM_BEG;
using namespace bl;

struct Cx { double re, im; };
static inline Cx  cx_add(Cx a, Cx b) { return { a.re + b.re, a.im + b.im }; }
static inline Cx  cx_sub(Cx a, Cx b) { return { a.re - b.re, a.im - b.im }; }
static inline Cx  cx_mul(Cx a, Cx b) { return { a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re }; }
static inline Cx  cx_div(Cx a, Cx b) {
    double d = b.re * b.re + b.im * b.im;
    return { (a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d };
}
static inline Cx  cx_scale(Cx a, double s) { return { a.re * s, a.im * s }; }
static inline double cx_abs2(Cx a) { return a.re * a.re + a.im * a.im; }
static inline double cx_abs(Cx a) { return std::sqrt(cx_abs2(a)); }
static inline bool   finite(Cx a) { return std::isfinite(a.re) && std::isfinite(a.im); }

static inline Cx cx_log_principal(Cx z) {
    double r = std::sqrt(cx_abs2(z));
    return { std::log(r), std::atan2(z.im, z.re) };
}

// Compute logPhi(c) and dlogPhi(c)=d/dc logPhi(c) for parameter c outside M.
// Returns false if c does not escape within max_iter (treat as inside/too close).
//
// Core idea:
//   After escape, z grows by squaring, so logPhi ≈ log(z_n) / 2^(n-1)
//   and dlogPhi ≈ (dz_n/z_n) / 2^(n-1).
static bool compute_logphi_dlogphi(Cx c,
    Cx& out_logPhi,
    Cx& out_dlogPhi,
    int& out_escape_iter,
    int max_iter = 200000)
{
    constexpr double TAU = 6.283185307179586476925286766559;
    constexpr double ESC_R2 = 4.0;     // |z| > 2 => escaped
    constexpr double Z_BIG2 = 1e24;    // once |z|~1e12, asymptotic approx is excellent
    constexpr double Z_STOP2 = 1e300;   // stop before overflow in squaring

    Cx z{ 0.0, 0.0 };
    Cx dz{ 0.0, 0.0 };  // dz/dc
    int n = 0;

    auto step = [&]() {
        // dz <- 2*z*dz + 1
        Cx two_z{ 2.0 * z.re, 2.0 * z.im };
        dz = cx_add(cx_mul(two_z, dz), Cx{ 1.0, 0.0 });
        // z <- z^2 + c
        z = cx_add(cx_mul(z, z), c);
        ++n;
    };

    // Iterate until escape
    while (n < max_iter && cx_abs2(z) <= ESC_R2) {
        step();
        if (!finite(z) || !finite(dz)) return false;
    }
    if (n >= max_iter) return false;

    out_escape_iter = n;

    // Unwrapped log(z) tracking
    Cx L = cx_log_principal(z);
    if (!finite(L)) return false;

    // Keep iterating a few steps until |z| is "big enough" (but not so big it overflows).
    // This keeps logPhi/dlogPhi accurate without fixed "refine_iters" that can explode.
    while (n < max_iter && cx_abs2(z) < Z_BIG2 && cx_abs2(z) < Z_STOP2) {
        Cx z_prev = z;
        double Lprev_im = L.im;

        step();
        if (!finite(z) || !finite(dz)) return false;

        Cx Ln = cx_log_principal(z);
        // unwrap imag(Ln) to be close to 2*imag(L)
        double target = 2.0 * Lprev_im;
        double k = std::nearbyint((target - Ln.im) / TAU);
        Ln.im += k * TAU;
        L = Ln;

        if (cx_abs2(z_prev) > Z_STOP2) break;
    }

    // scale = 2^{-(n-1)} (exact for powers of two)
    const int shift = -(n - 1);
    const double scale = std::scalbn(1.0, shift);
    if (scale == 0.0) return false;

    out_logPhi = { L.re * scale,  L.im * scale };

    // dlogPhi ≈ (dz/z) * scale
    if (cx_abs2(z) < 1e-300) return false;
    Cx dz_over_z = cx_div(dz, z);
    out_dlogPhi = { dz_over_z.re * scale, dz_over_z.im * scale };

    return finite(out_logPhi) && finite(out_dlogPhi);
}

// Follow a “ray” in parameter space by descending the potential g = Re(logPhi).
// We use k as the number of steps, each aiming to reduce g by DG (adaptive near boundary).
//
// This gives you a stable, monotone “inward” path from |c|=8 that matches the flow-line intuition
// in your screenshot without any inversion solves.
/*bool mandel_angle_iter_to_xy(double angle_rad, int k, double& out_x, double& out_y)
{
    // Your shell radius (iter=1 with bailout 8): start point c0 = 8 * e^{i angle}
    constexpr double START_R = 4.0;

    // Potential drop per step (tune to taste):
    // Larger => moves inward faster but less safe near boundary.
    constexpr double DG = 0.05;

    // Stop before we are “too close” to the boundary for double precision.
    constexpr double G_MIN = 1e-12;

    angle_rad = std::fmod(angle_rad, 6.283185307179586476925286766559);
    if (angle_rad < 0) angle_rad += 6.283185307179586476925286766559;

    Cx c{ START_R * std::cos(angle_rad), START_R * std::sin(angle_rad) };

    for (int i = 0; i < k; ++i)
    {
        Cx logPhi{}, dlogPhi{};
        int escIter = 0;
        if (!compute_logphi_dlogphi(c, logPhi, dlogPhi, escIter))
            return false; // should not happen for exterior points unless extremely near boundary

        double g = logPhi.re; // potential
        if (!(g > G_MIN)) {
            // Essentially at the boundary; further steps are not meaningful in doubles.
            break;
        }

        // Gradient of g from analytic derivative f'(c)=d/dc logPhi:
        // For f=u+iv analytic: grad(u) = (Re f', -Im f')
        const double gx = dlogPhi.re;
        const double gy = -dlogPhi.im;

        const double grad2 = gx * gx + gy * gy;
        if (!(grad2 > 0.0) || !std::isfinite(grad2))
            return false;

        // Desired potential drop this step; don’t overshoot past g=0.
        double dg = std::min(DG, 0.5 * g);

        // Backtracking to keep us safely outside and monotone in g
        bool accepted = false;
        for (int bt = 0; bt < 12; ++bt)
        {
            // Choose dc so that g decreases by approximately dg:
            // dg ≈ grad · dc => with dc = -(dg/|grad|^2) grad => dg ≈ -dg
            const double s = dg / grad2;
            Cx c_try{ c.re - gx * s, c.im - gy * s };

            Cx logPhi2{}, dlogPhi2{};
            int esc2 = 0;
            if (compute_logphi_dlogphi(c_try, logPhi2, dlogPhi2, esc2))
            {
                // enforce monotone decrease of potential
                if (logPhi2.re < g) {
                    c = c_try;
                    accepted = true;
                    break;
                }
            }

            dg *= 0.5; // shrink step and retry
            if (dg < 1e-18) break;
        }

        if (!accepted) return false;
    }

    out_x = c.re;
    out_y = c.im;
    return true;
}
*/

static inline double unwrap_to(double value, double target)
{
    constexpr double TAU = 6.283185307179586476925286766559;
    double k = std::nearbyint((target - value) / TAU);
    return value + k * TAU;
}

static inline double escapeDepth_to_g(int N, double bailout_R = 2.0)
{
    // g_N = ln(R) * 2^{-(N-1)}
    if (N <= 1) return std::log(bailout_R);
    return std::scalbn(std::log(bailout_R), -(N - 1));
}

static inline bool eval_logphi_locked_theta(
    Cx c,
    double theta_locked,
    Cx& out_logPhi,
    Cx& out_dlogPhi,
    int& out_escape_iter,
    int max_iter)
{
    if (!compute_logphi_dlogphi(c, out_logPhi, out_dlogPhi, out_escape_iter, max_iter))
        return false;

    // Lock imag(logPhi) to the same branch as theta_locked to prevent 2π jumps.
    out_logPhi.im = unwrap_to(out_logPhi.im, theta_locked);

    return finite(out_logPhi) && finite(out_dlogPhi);
}

bool mandel_angle_iter_to_xy(double angle_rad, int k, double& out_x, double& out_y)
{
    constexpr double TAU = 6.283185307179586476925286766559;
    constexpr double BAILOUT_R = 2.0;

    // Use a larger start radius: makes theta determination and initial guess more well-conditioned.
    constexpr double START_R = 16.0;

    // Practical double floor
    constexpr double G_MIN = 1e-15;

    // Residual tolerances in logPhi space
    constexpr double G_TOL = 2e-14;
    constexpr double TH_TOL = 2e-12;

    constexpr int MAX_NEWTON_ITERS = 20;
    constexpr int MAX_BACKTRACK = 14;

    angle_rad = std::fmod(angle_rad, TAU);
    if (angle_rad < 0) angle_rad += TAU;

    int N = std::max(1, k);
    double g_target = escapeDepth_to_g(N, BAILOUT_R);
    if (g_target < G_MIN) g_target = G_MIN;

    int max_iter = std::min(200000, std::max(256, 6 * N + 128));

    Cx c0{ START_R * std::cos(angle_rad), START_R * std::sin(angle_rad) };

    // Determine locked theta from the starting point
    Cx logPhi0{}, dlogPhi0{};
    int esc0 = 0;
    if (!compute_logphi_dlogphi(c0, logPhi0, dlogPhi0, esc0, max_iter))
        return false;

    if (!(logPhi0.re > 0.0) || !finite(logPhi0))
        return false;

    const double theta = logPhi0.im;

    // Initial guess: radial scaling in c-plane (good enough when START_R is large)
    double scale = std::exp(std::clamp(g_target - logPhi0.re, -8.0, 8.0));
    Cx c = cx_scale(c0, scale);

    Cx target{ g_target, theta };

    // Trust region cap: prevents catastrophic jumps. (Tune if you like.)
    auto max_step_for = [](Cx ccur) -> double {
        double r = cx_abs(ccur);
        return std::max(0.05, 0.25 * r); // cap to 25% of radius, min 0.05
    };

    for (int it = 0; it < MAX_NEWTON_ITERS; ++it)
    {
        Cx logPhi{}, dlogPhi{};
        int esc = 0;
        if (!eval_logphi_locked_theta(c, theta, logPhi, dlogPhi, esc, max_iter))
            return false;

        if (!(logPhi.re > G_MIN))
            return false;

        Cx F = cx_sub(logPhi, target);

        if (std::abs(F.re) <= G_TOL && std::abs(F.im) <= TH_TOL) {
            out_x = c.re;
            out_y = c.im;
            return true;
        }

        double d2 = cx_abs2(dlogPhi);
        if (!(d2 > 0.0) || !std::isfinite(d2))
            return false;

        // Newton step in c: delta = F / dlogPhi
        Cx delta = cx_div(F, dlogPhi);

        // Trust region cap
        double step_len = cx_abs(delta);
        double step_cap = max_step_for(c);
        if (step_len > step_cap) {
            delta = cx_scale(delta, step_cap / step_len);
        }

        // Backtracking: accept only if residual decreases AND stays exterior
        bool accepted = false;
        double t = 1.0;

        double err = std::max(std::abs(F.re), std::abs(F.im));

        for (int bt = 0; bt < MAX_BACKTRACK; ++bt)
        {
            Cx c_try = cx_sub(c, cx_scale(delta, t));

            Cx logPhi2{}, dlogPhi2{};
            int esc2 = 0;
            if (!eval_logphi_locked_theta(c_try, theta, logPhi2, dlogPhi2, esc2, max_iter)) {
                t *= 0.5;
                continue;
            }

            if (!(logPhi2.re > G_MIN)) {
                t *= 0.5;
                continue;
            }

            Cx F2 = cx_sub(logPhi2, target);
            double err2 = std::max(std::abs(F2.re), std::abs(F2.im));

            if (err2 < err) {
                c = c_try;
                accepted = true;
                break;
            }

            t *= 0.5;
        }

        if (!accepted)
            return false;
    }

    return false;
}

//static inline DVec2 operator+(DVec2 a, DVec2 b) { return { a.x + b.x, a.y + b.y }; }
//static inline DVec2 operator-(DVec2 a, DVec2 b) { return { a.x - b.x, a.y - b.y }; }
//static inline DVec2 operator*(DVec2 a, double s) { return { a.x * s, a.y * s }; }
//
//static inline double dot(DVec2 a, DVec2 b) { return a.x * b.x + a.y * b.y; }
//static inline double len2(DVec2 a) { return dot(a, a); }
//static inline double len(DVec2 a) { return std::sqrt(len2(a)); }

static double point_to_segment_distance(DVec2 p, DVec2 a, DVec2 b)
{
    DVec2 ab = b - a;
    double ab2 = ab.mag2();// len2(ab);
    //if (ab2 <= 0.0) return len(p - a);
    if (ab2 <= 0.0) return (p - a).mag();
    double t = (p - a).dot(ab) / ab2;
    t = std::clamp(t, 0.0, 1.0);
    DVec2 q = a + ab * t;
    //return len(p - q);
    return (p - q).mag();
}

static bool sample_point(double angle, int input_iters, DVec2& out)
{
    double x, y;
    if (!mandel_angle_iter_to_xy(angle, input_iters, x, y))
        return false;
    out = { x, y };
    return true;
}

// Adaptive subdivision from (a0,p0) to (a1,p1); appends points to 'out' in order.
// 'out' must already contain p0.
static bool emit_adaptive_segment(
    std::vector<DVec2>& out,
    double a0, DVec2 p0,
    double a1, DVec2 p1,
    int input_iters,
    double max_chord,
    double max_sag,
    int max_depth)
{
    struct Seg { double a0, a1; DVec2 p0, p1; int d; };

    std::vector<Seg> stack;
    stack.reserve(256);
    stack.push_back({ a0, a1, p0, p1, 0 });

    while (!stack.empty())
    {
        Seg s = stack.back();
        stack.pop_back();

        double chord = (s.p1 - s.p0).mag();

        if (s.d >= max_depth)
        {
            out.push_back(s.p1);
            continue;
        }

        // Midpoint
        double am = 0.5 * (s.a0 + s.a1);
        DVec2 pm;
        if (!sample_point(am, input_iters, pm))
            return false;

        double sag = point_to_segment_distance(pm, s.p0, s.p1);

        if (chord > max_chord || sag > max_sag)
        {
            // Subdivide: push right then left so left is processed first.
            stack.push_back({ am, s.a1, pm, s.p1, s.d + 1 });
            stack.push_back({ s.a0, am, s.p0, pm, s.d + 1 });
        }
        else
        {
            out.push_back(s.p1);
        }
    }

    return true;
}

// Builds a closed polyline around [0, 2pi], adaptive in world-space.
// Returns empty vector on failure.
std::vector<DVec2> build_mandel_outline_path(int input_iters, int quality, int thread_count)
{
    //constexpr double PI = 3.1415926535897932384626433832795;
    constexpr double TAU = 6.283185307179586476925286766559;

    quality = std::clamp(quality, 1, 50);
    thread_count = std::max(1, thread_count);

    // Heuristics: higher quality => smaller chord/sag thresholds and deeper subdivision.
    // baseSeg governs the overall density; adaptive subdivision adds detail where needed.
    const int baseSeg = 256 * quality;

    const double max_chord = 6.0 / double(baseSeg);       // world units
    const double max_sag = 0.20 * max_chord;            // sagitta tolerance
    const int    max_depth = std::clamp(8 + 2 * quality, 12, 40);

    // Split angle domain into contiguous chunks (deterministic merge).
    const int tc = std::min(thread_count, baseSeg); // avoid more threads than work
    std::vector<double> bounds(tc + 1);
    for (int i = 0; i <= tc; ++i)
        bounds[i] = TAU * (double(i) / double(tc));

    // Precompute boundary points sequentially (stabilizes joins).
    std::vector<DVec2> bound_pts(tc + 1);
    for (int i = 0; i <= tc; ++i)
    {
        if (!sample_point(bounds[i], input_iters, bound_pts[i]))
            return {}; // cannot build a closed path for these parameters
    }

    std::atomic<bool> ok{ true };
    std::vector<std::vector<DVec2>> per(tc);
    std::vector<std::future<void>> futures(tc);

    for (int ti = 0; ti < tc; ++ti)
    {
        futures[ti] = Thread::pool().submit_task([&, ti]()
        {
            if (!ok.load(std::memory_order_relaxed))
                return;

            const double a0 = bounds[ti];
            const double a1 = bounds[ti + 1];
            const DVec2   p0 = bound_pts[ti];
            const DVec2   p1 = bound_pts[ti + 1];

            // Start with p0; adaptive emitter appends to p1
            std::vector<DVec2> local;
            local.reserve((baseSeg / tc) * 4);
            local.push_back(p0);

            // Further split each chunk into a few coarse segments before adapting (reduces recursion pressure).
            const int coarse = std::max(4, (baseSeg / tc) / 4);
            double da = (a1 - a0) / double(coarse);

            double ac = a0;
            DVec2 pc = p0;

            for (int i = 1; i <= coarse; ++i)
            {
                double an = (i == coarse) ? a1 : (a0 + da * double(i));
                DVec2 pn;
                if (!sample_point(an, input_iters, pn)) { ok.store(false); return; }

                if (!emit_adaptive_segment(local, ac, pc, an, pn, input_iters,
                    max_chord, max_sag, max_depth))
                {
                    ok.store(false); return;
                }

                ac = an;
                pc = pn;
            }

            per[ti] = std::move(local);
        });
    }

    for (int ti = 0; ti < tc; ++ti)
    {
        if (futures[ti].valid())
            futures[ti].get();
    }

    if (!ok.load())
        return {};

    // Merge in order; drop duplicate first point of each subsequent chunk.
    std::vector<DVec2> path;
    path.reserve(baseSeg * 2);

    for (int ti = 0; ti < tc; ++ti)
    {
        const auto& seg = per[ti];
        if (seg.empty())
            return {};

        if (ti == 0) path.insert(path.end(), seg.begin(), seg.end());
        else         path.insert(path.end(), seg.begin() + 1, seg.end());
    }

    // Ensure closed (last point equals first). If not, append first.
    if (!path.empty())
    {
        DVec2 a = path.front(), b = path.back();
        if (std::abs(a.x - b.x) > 0.0 || std::abs(a.y - b.y) > 0.0)
            path.push_back(a);
    }

    return path;
}

void Mandelbrot_Scene::processExperimental()
{
    if (first_frame || Changed(input_iters, input_quality))
    {
        double inc = math::pi / 100.0;
        iter_paths.clear();
        for (double angle = 0.0; angle < math::pi * 2; angle += inc)
        {
            double px, py;
            int i = 0;

            iter_paths.push_back(std::vector<DVec2>());
            std::vector<DVec2>& path = iter_paths.back();

            mandel_angle_iter_to_xy(angle, i++, px, py);
            path.push_back({ px, py });

            for (; i < input_iters; i++)
            {
                mandel_angle_iter_to_xy(angle, i, px, py);
                path.push_back({ px, py });
            }
        }

        ring_paths.clear();
        for (int i = 1; i < input_iters; i++)
            ring_paths.push_back(build_mandel_outline_path(i, input_quality, Thread::threadCount()));
    }
}

void Mandelbrot_Scene::drawExperimental(Viewport* ctx) const
{
    //double inc = math::pi / 50000.0;
    //double px, py;
    //mandel_angle_iter_to_xy(0.0, input_iters, px, py);
    //
    //ctx->worldHudMode();
    //ctx->setLineWidth(1);
    //ctx->setStrokeStyle(Color::red);
    //ctx->beginPath();
    //ctx->moveTo(px, py);
    //for (double angle = 0.0; angle < math::pi * 2; angle += inc)
    //{
    //    mandel_angle_iter_to_xy(angle, input_iters, px, py);
    //    ctx->lineTo(px, py);
    //}
    //ctx->stroke();
    //ctx->worldMode();

    ctx->worldHudMode();
    ctx->setLineWidth(1);
    ctx->setStrokeStyle(Color::red);

    for (auto& path : iter_paths)
        ctx->strokePath(path);

    for (auto& path : ring_paths)
        ctx->strokePath(path);

    ctx->worldMode();
}

void Mandelbrot_Scene::UI::populateExperimental()
{
    if (ImGui::Section("Experimental")) {
    {
        bl_scoped(input_angle, input_iters, input_quality);
        ImGui::SliderDouble("input_angle", &input_angle, 0.0, math::pi * 2.0);
        ImGui::SliderInt("input_iters", &input_iters, 1, 200);
        ImGui::SliderInt("input_quality", &input_quality, 1, 200);
    }
}

void Mandelbrot_Scene::UI::populateSplinesDev()
{
    bl_scoped(tween_pos_spline);
    bl_scoped(tween_zoom_lift_spline);
    bl_scoped(tween_base_zoom_spline);

    static ImRect vr = { 0.0f, 1.0f, 1.0f, 0.0f };

    ImGui::SeparatorText("Position Tween");
    ImSpline::SplineEditor("tween_pos", &tween_pos_spline, &vr);

    ImGui::SeparatorText("Lift Tween");
    ImSpline::SplineEditor("tween_zoom_lift", &tween_zoom_lift_spline, &vr);

    ImGui::SeparatorText("Base Zoom Tween");
    ImSpline::SplineEditor("tween_base_zoom", &tween_base_zoom_spline, &vr);

    //ImGui::InputTextMultiline("###pos_buf", pos_tween_buf, 1024, ImVec2(0, 0), ImGuiInputTextFlags_AllowTabInput));
    if (ImGui::Button("Copy position spline")) ImGui::SetClipboardText(tween_pos_spline.serialize(SplineSerializationMode::CPP_ARRAY, 3).c_str());

    //ImGui::InputTextMultiline("###pos_buf", zoom_tween_buf, 1024, ImVec2(0, 0), ImGuiInputTextFlags_AllowTabInput));
    if (ImGui::Button("Copy lift spline")) ImGui::SetClipboardText(tween_zoom_lift_spline.serialize(SplineSerializationMode::CPP_ARRAY, 3).c_str());

    if (ImGui::Button("Copy base zoom spline")) ImGui::SetClipboardText(tween_base_zoom_spline.serialize(SplineSerializationMode::CPP_ARRAY, 3).c_str());

    ///static std::string results;
    ///static std::string spline_results;
    ///
    ///bl_scoped(stripe_mag_numerator, stripe_mag_from_hist);
    ///ImGui::Checkbox("Stripe Use Histogram", &stripe_mag_from_hist);
    ///ImGui::SliderFloat("stripe_mag_numerator", &stripe_mag_numerator, 0.01f, 0.3f);
    ///
    ///bl_scoped(stripe_zf_spline, stripe_zf_spline_rect);
    ///bl_pull(ideal_zf_numerator_map);
    ///
    ///if (ImSpline::BeginSplineEditor("ZF Spline", &stripe_zf_spline, &stripe_zf_spline_rect, 300.0f, ImSplineFlags_InvertY))
    ///{
    ///    if (ideal_zf_numerator_map.size()) {
    ///        for (const auto& [key, value] : ideal_zf_numerator_map)
    ///            ImSpline::PlotPoint(key, value);
    ///    }
    ///    ImSpline::EndSplineEditor();
    ///}
    ///
    ///if (ImGui::Button("Update Results"))
    ///{
    ///    std::stringstream ss;
    ///    for (const auto& [key, value] : ideal_zf_numerator_map)
    ///        ss << key << " = " << value << "\n";
    ///
    ///    results = ss.str();
    ///    spline_results = stripe_zf_spline.serialize(SplineSerializationMode::CPP_ARRAY, 3);
    ///}
    ///ImGui::InputTextMultiline("Results", &results);
    ///ImGui::InputTextMultiline("Spline", &spline_results);
}


SIM_END;
#endif