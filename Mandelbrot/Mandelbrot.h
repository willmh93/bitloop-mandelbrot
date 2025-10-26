#pragma once
#include <bitloop.h>
#include <cmath>

#include "Cardioid/Cardioid.h"

// Mandelbrot includes
#include "build_config.h"

#include "types.h"
#include "mandel_state.h"
#include "compute.h"
#include "shading.h"
#include "gradient.h"
#include "tween_splines.h"
#include "examples.h"
#include "mandel_stats.h"

SIM_BEG;

using namespace bl;



//  ──────  MandelState inheritance  ────── 
//  Lerp between state_a ─> state_b saves result in "this" inherited MandelState

struct Mandelbrot_Scene : public MandelState, public Scene<Mandelbrot_Scene>
{
    // ────── config ──────
    struct Config { std::string load_example_name; };
    Mandelbrot_Scene(Config& config) : load_example_name(config.load_example_name) {}

    static inline MandelExampleMap mandel_presets = generateMandelPresets();

    std::string load_example_name;

    // ────── threads ──────
    static constexpr int MAX_THREADS = 0; // 0 = Use max threads
    inline int numThreads() const {
        if constexpr (MAX_THREADS > 0) return MAX_THREADS;
        return Thread::idealThreadCount();
    }

    // ────── resources ──────
    NanoFont font;

    // ────── tweening ──────
    bool tweening = false;
    f64  tween_progress = 0; // 0..1
    f64  tween_duration = 0;
    f128 tween_lift = 0;

    MandelState  state_a;
    MandelState  state_b;
    DDAngledRect tween_r1;
    DDAngledRect tween_r2;

    // For any given MandelState, we need a way of predicting what the destination world-quad will be for tweening.
    // That means we need to use the active viewport size + the target MandelState's camera info.
    DDAngledRect getAngledRect(const MandelState& s) const {
        return DDAngledRect(s.camera.pos<f128>(), camera.viewportWorldSize<f128>(), (f128)s.camera.rotation());
    }

    // ────── gradients ──────
    ImGradient gradient_shifted;


    // ────── Splines ──────
    ImSpline::Spline x_spline               = MandelSplines::x_spline;
    ImSpline::Spline y_spline               = MandelSplines::y_spline;
    ImSpline::Spline tween_pos_spline       = MandelSplines::tween_pos_spline;
    ImSpline::Spline tween_zoom_lift_spline = MandelSplines::tween_zoom_lift_spline;
    ImSpline::Spline tween_base_zoom_spline = MandelSplines::tween_base_zoom_spline;
    ImSpline::Spline tween_color_cycle      = MandelSplines::tween_color_cycle;

    // ────── saving & loading / URL ──────
    //std::string data_buf = ""; // data for text box (save data / url)
    //void updateConfigBuffer() { data_buf = serialize(); }
    //void loadConfigBuffer()   { deserialize(data_buf);  }

    std::string getStateData() const { return serialize(); }
    void loadState(std::string data) { deserialize(data); }

    //void onSavefileChanged();
    std::string getURL() const;

    // ────── fields & bitmaps ──────
    EscapeField field_9x9 = EscapeField(0); // processed in a single frame (fast)
    EscapeField field_3x3 = EscapeField(1); // processed over multiple frames
    EscapeField field_1x1 = EscapeField(2); // processed over multiple frames

    CanvasImage128 bmp_9x9;
    CanvasImage128 bmp_3x3;
    CanvasImage128 bmp_1x1;

    CanvasImage128* pending_bmp = nullptr;
    CanvasImage128* active_bmp = nullptr;

    EscapeField* pending_field = nullptr;
    EscapeField* active_field = nullptr;

    MandelStats stats;

    // ────── dynamicly set at runtime ──────
    bool    display_intro = true;
    double  log_color_cycle_iters = 0.0;
    int     iter_lim = 0; // Actual iter limit
    int     current_row = 0;

    // phase 0 = 9x smaller
    // phase 1 = 3x smaller
    // phase 2 = full resolution
    int  computing_phase = 0;
    bool first_frame = true;
    bool frame_complete = false;  // Similar to finished_compute, but not cleared until next compute starts
    bool final_frame_complete = true;
    DDQuad world_quad{};

    // ────── expensive interior "forwarding" optimization ──────
    struct {
        int c1 = 5, e1 = 3, c2 = 7, e2 = 3;
    }  interior_phases_contract_expand;

    bool maxdepth_show_optimized = false;

    // ────── timers ──────
    std::chrono::steady_clock::time_point compute_t0;
    Math::MovingAverage::MA<double> timer_ma = Math::MovingAverage::MA<double>(1);
    double dt_avg = 0;

    // ────── tweening ──────
    void startTween(const MandelState& target);
    void lerpState(const MandelState& a, const MandelState& b, double f, bool complete);

    // ────── camera navigation easing ──────
    CameraNavigator navigator;

    Math::MovingAverage::MA<DDVec2> avg_vel_pos = Math::MovingAverage::MA<DDVec2>(8);
    Math::MovingAverage::MA<double> avg_vel_zoom = Math::MovingAverage::MA<double>(8);

    DDVec2 camera_vel_pos{};
    double camera_vel_zoom{1};

    // ────── deep-zoom animation (temporary until timeline feature added) ──────
    bool steady_zoom = false;
    f128 steady_zoom_mult_speed{ 0.01 };
    int tween_frames_elapsed = 0;
    int tween_expected_frames = 0;
    Math::MovingAverage::MA<double> expected_time_left_ma = Math::MovingAverage::MA<double>(5);

    // ────── cardioid ──────
    bool show_period2_bulb = true;

    #if MANDEL_FEATURE_INTERACTIVE_CARDIOID
    bool show_interactive_cardioid = false;
    bool   animate_cardioid_angle = true;
    double ani_angle = 0.0;
    double ani_inc = Math::toRadians(2.0);
    #endif

    double cardioid_lerp_amount = 1.0; // 1 - flatten
    Cardioid::CardioidLerper cardioid_lerper;

    // ────── user interface ──────
    struct UI : Interface
    {
        using Interface::Interface;
        void sidebar();


        // dialog flags
        bool show_save_dialog = false;
        bool show_load_dialog = false;
        bool show_share_dialog = false;
        #ifdef __EMSCRIPTEN__
        bool opening_load_popup = false;
        #endif

        // textbox buffers
        std::string url_buf;
        std::string data_buf;

        // stats
        std::vector<int> depth_xs, depth_ys;
        //std::vector<double> max_depth_xs, max_depth_ys;

        // UI sections
        void populateSavingLoading();
        void populateExamples();
        void populateCameraView();
        void populateQualityOptions();
        void populateColorCycleOptions();
        void populateGradientShiftOptions();
        void populateGradientPicker();
        void populateStats();

        void populateExperimental();
        void populateSplinesDev();
    };

    // ────── simulation processing ──────
    void sceneStart() override;
    void sceneMounted(Viewport* viewport) override;

    // ────── viewport handling ──────
    bool mandelChanged();
    bool shadingFormulaChanged();

    void updateAnimation();
    void updateTweening(double dt);
    bool updateGradient();
    void updateCameraView();
    void updateFieldSizes(Viewport* ctx);
    void updateEnabledKernelFeatures();
    void updateActivePhaseAndField();
    bool processCompute();

    void viewportProcess(Viewport* ctx, double dt) override;
    void viewportDraw(Viewport* ctx) const override;

    void collectStats();

    // ────── input ──────
    void onEvent(Event e) override;
};

struct Mandelbrot_Project : public BasicProject
{
    static ProjectInfo info() {
        return ProjectInfo({ "Fractal", "Mandelbrot", "Mandelbrot Viewer" });
    }

    void projectPrepare(Layout& layout) override;
};

#undef DEV_MODE

SIM_END;
