#pragma once

// Include build_config.h before bitloop.h
#include "build_config.h"

#include <bitloop.h>
#include <cmath>

#include "Cardioid/Cardioid.h"

// Mandelbrot includes

#include "types.h"
#include "conversions.h"
#include "mandel_state.h"
#include "shading.h"
#include "gradient.h"
#include "tween_splines.h"
#include "examples.h"
#include "mandel_stats.h"

SIM_BEG;

using namespace bl;



//  ──────  MandelState inheritance  ────── 
//  Tween lerps from state_a ─> state_b, saves result in inherited MandelState

struct Mandelbrot_Scene : public Scene<Mandelbrot_Scene>, public MandelState
{
    // ────── config ──────
    struct Config { std::string load_example_name; };
    Mandelbrot_Scene(Config& config) : load_example_name(config.load_example_name) {}
    std::string load_example_name;

    // ────── examples ──────
    static inline const MandelExampleMap mandel_presets = generateMandelPresets();
    bool rendering_examples = false;
    int rendering_example_i = 0;
    std::string render_batch_name;


    // ────── threads ──────
    static constexpr int MAX_THREADS = 0; // 0 = max threads

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
    // This means we need to use the active viewport size + the target MandelState's camera info.
    DDAngledRect getAngledRect(const MandelState& s) const {
        return DDAngledRect(s.camera.pos<f128>(), camera.viewportWorldSize<f128>(), (f128)s.camera.rotation());
    }

    // ────── gradients ──────
    ImGradient gradient_shifted;

    // ────── shading ──────
    ToneManager<f64> dist_tone;
    ToneManager<f32> stripe_tone;


    // ────── Splines ──────
    ImSpline::Spline x_spline               = MandelSplines::x_spline;
    ImSpline::Spline y_spline               = MandelSplines::y_spline;
    ImSpline::Spline tween_pos_spline       = MandelSplines::tween_pos_spline;
    ImSpline::Spline tween_zoom_lift_spline = MandelSplines::tween_zoom_lift_spline;
    ImSpline::Spline tween_base_zoom_spline = MandelSplines::tween_base_zoom_spline;
    ImSpline::Spline tween_color_cycle      = MandelSplines::tween_color_cycle;

    // ────── saving & loading / URL ──────
    std::string serializeState() const { return serialize(); }
    void loadState(std::string_view data) { deserialize(data); }

    std::string getURL() const;


    // ────── fields for each phase (raw float data) ──────
    EscapeField field_9x9 = EscapeField(0); // processed in a single frame (usually fast enough for smooth interaction)
    EscapeField field_3x3 = EscapeField(1); // processed over multiple frames
    EscapeField field_1x1 = EscapeField(2); // processed over multiple frames (final compute phase)

    EscapeField* pending_field = nullptr; // ptr to unfinished field
    EscapeField* active_field = nullptr;  // ptr to finished normalized field (shades to active_bmp)


    // ────── bitmaps for each phase (RGBA data, with multithreaded forEachWorldTilePixel helpers) ──────
    CanvasImage128 bmp_9x9;
    CanvasImage128 bmp_3x3;
    CanvasImage128 bmp_1x1;

    CanvasImage128* pending_bmp = nullptr; // ptr to the unfinished bitmap
    CanvasImage128* active_bmp = nullptr;  // ptr to the finished displayed bitmap


    bool iter_hist_visible = false;
    bool dist_hist_visible = false;
    bool stripe_hist_visible = false;

    // ────── normalization ──────
    NormalizationField norm_field;

    /// todo: lerped means/ranges for smoother zooms/pans
    //NormalizationField norm_field_zoom_out;
    //NormalizationField norm_field_zoom_in;
    //NormalizationField norm_field_pan_tl;
    //NormalizationField norm_field_pan_tr;
    //NormalizationField norm_field_pan_bl;
    //NormalizationField norm_field_pan_br;

    bool preview_normalization_field = false;

    /// dev: stripe normalization debug info
    //float stripe_mag_numerator = 0.1f;
    //bool  stripe_mag_from_hist = false;
    //std::map<float, float> ideal_zf_numerator_map;
    //std::vector<std::pair<float, float>> ideal_zf_numerator_map;
    //ImRect stripe_zf_spline_rect = ImRect(0.0f, 0.0f, 30.0f, 1.0f);
    
    // ────── stats ──────
    MandelStats stats;
    

    // ────── dynamic ──────
    bool           display_intro = true;
    double         log_color_cycle_iters = 0.0;
    int            iter_lim = 0; // real iteration limit fed to mandelbrot kernel
    KernelFeatures mandel_features = KernelFeatures::ITER;
    KernelMode     kernel_mode = KernelMode::AUTO;

    // multi-frame tile progress tracker for the pending phase
    TileBlockProgress P;

    // phase 0 = 9x smaller field
    // phase 1 = 3x smaller field
    // phase 2 = full resolution field
    int  computing_phase = 0;
    bool first_frame = true;
    bool frame_complete = false;  // Similar to finished_compute, but not cleared until next compute starts
    bool final_frame_complete = true;
    DDQuad world_quad{};

    // ────── interior-forwarding optimization ──────
    struct {
        // c1, e1: [contract, expand] phase 0 interior => forward to phase 1
        // c2, e2: [contract, expand] phase 1 interior => forward to phase 2
        // This prevents forwarding results of interior edges (but can fail at sharp angles if too aggressive)
        int c1 = 5, e1 = 3, c2 = 7, e2 = 3;
    }  interior_phases_contract_expand;

    bool maxdepth_show_optimized = false;

    // ────── timers ──────
    std::chrono::steady_clock::time_point compute_t0;
    math::SMA<double> timer_ma = math::SMA<double>(1);
    double dt_avg = 0;

    // ────── tweening ──────
    void startTween(const MandelState& target);
    double lerpState(const MandelState& a, const MandelState& b, double f, bool complete);

    // ────── camera navigation easing ──────
    CameraNavigator navigator;

    math::SMA<DDVec2> avg_vel_pos = math::SMA<DDVec2>(8);
    math::SMA<double> avg_vel_zoom = math::SMA<double>(8);

    DDVec2 camera_vel_pos{};
    double camera_vel_zoom{1};

    // ────── deep-zoom animation (temporary until timeline feature added) ──────
    bool steady_zoom = false;
    f128 steady_zoom_mult_speed{ 0.01 };
    //f128 steady_zoom_mult_speed{ 0.05 }; // faster test zoom
    float steady_zoom_pct = 0;
    int tween_frames_elapsed = 0;
    int tween_expected_frames = 0;
    math::SMA<f64> expected_time_left_ma = math::SMA<f64>(5); // TODO: Use or remove

    // ────── undo/redo ──────
    const int history_count = 100;
    bool pending_checkpoint_flag = false;
    bool loading_history_state = false;
    std::vector<std::pair<int, std::string>> history;
    int current_history_i = 0;

    // ────── cardioid ──────
    bool show_period2_bulb = true;

    #if MANDEL_FEATURE_INTERACTIVE_CARDIOID
    bool show_interactive_cardioid = false;
    bool   animate_cardioid_angle = true;
    double ani_angle = 0.0;
    double ani_inc = math::toRadians(2.0);
    #endif

    // ────── extra features ──────
    double cardioid_lerp_amount = 1.0; // 0 = flatten
    Cardioid::CardioidLerper cardioid_lerper;

    // ────── user interface ──────
    struct UI : ViewModel
    {
        using ViewModel::ViewModel;


        MandelBookmarkManager bookmark_manager;
        MandelBookmarkList examples_list;

        // dialog flags
        bool show_save_dialog = false;
        bool show_load_dialog = false;
        bool show_share_dialog = false;
        #ifdef __EMSCRIPTEN__
        bool opening_load_popup = false;
        #endif
        
        bool show_iter_hist = false;
        bool show_dist_hist = false;
        bool show_stripe_hist = false;


        // textbox buffers
        std::string url_buf;
        std::string data_buf;

        double orbit_padding = 0.3;

        int selected_preset_i = 0;

        void init() override;
        void overlay() override;
        void sidebar() override;

        void launchBookmark(std::string_view data);


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
        void populateMouseOrbit();
        void populateHistory();

        void populateSplinesDev();


        void populateCaptureOptions();
    };

    // ────── simulation processing ──────
    void sceneStart() override;
    void sceneDestroy() override;
    void sceneMounted(Viewport* viewport) override;

    // ────── compute ──────
    template<typename T, KernelFeatures F, KernelMode K>  
    bool compute_mandelbrot(int timeout);

    template<typename T, KernelFeatures F>
    void calculate_normalize_info();

    template<typename T, KernelFeatures F, bool Normalize_Depth, bool Invert_Dist>
    void normalize_field();

    template<bool Highlight_Optimized_Interior>
    void shadeBitmap();

    // ────── viewport handling ──────
    bool mandelChanged();
    bool normalizationOptionsChanged();

    void processBatchSnapshot();
    void updateAnimation();
    void updateTweening(double dt);
    bool updateGradient();
    void updateCameraView();
    void updateFieldSizes(Viewport* ctx);
    void updateEnabledKernelFeatures();
    void updateActivePhaseAndField();
    bool processCompute();
    void processCapturing(bool finished_compute, bool reshade);
    void processUndoRedo(bool normalization_opts_changed);

    #if MANDEL_EXPERIMENTAL_TESTS
    double input_angle = 0.0;
    int    input_iters = 20;
    int    input_quality = 2;

    double plot_x, plot_y;
    std::vector<std::vector<DVec2>> iter_paths;
    std::vector<std::vector<DVec2>> ring_paths;

    void processExperimental();
    void drawExperimental(Viewport* ctx) const;
    #endif

    void viewportProcess(Viewport* ctx, double dt) override;
    void viewportDraw(Viewport* ctx) const override;

    void collectStats(bool renormalized);

    //void onEncodeFrame(bytebuf& data, int request_id, const SnapshotPreset& preset) override
    //{}

    // ────── input ──────
    void onEvent(Event e) override;
    void onKeyDown(KeyEvent e) override;

    void historyKeyEvent(KeyEvent e);
};

struct Mandelbrot_Project : public BasicProject
{
    static ProjectInfo info() {
        return ProjectInfo({ "Fractal", "Mandelbrot", "Mandelbrot Viewer" });
    }

    void projectPrepare(Layout& layout) override;
};

SIM_END;

