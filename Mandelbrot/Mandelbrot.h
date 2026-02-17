#pragma once

// Include build_config.h before bitloop.h
#include "config/build_config.h"

#include <bitloop.h>
#include <cmath>

// subprojects
#include "Cardioid/Cardioid.h"

// external
#include "external/TextEditor.h"

// Mandelbrot includes
#include "core/types.h"
#include "core/conversions.h"
#include "state/mandel_state.h"
#include "shading/shading.h"
#include "shading/gradient.h"
#include "config/tween_splines.h"
#include "state/examples.h"
#include "field/escape_field.h"
#include "field/normalization_field.h"
#include "stats/mandel_stats.h"

#include "kernels/kernel_mandel.hpp"
#include "kernels/ref_orbit.hpp"

#if MANDEL_UNVERSIONED_EXPERIMENTAL
#include "_unversioned/experimental.h"
#endif

SIM_BEG;

using namespace bl;

// each phase given it's own field / raster grid mapping / shader
struct MandelPhaseData
{
    int phase = -1;

    EscapeField         field;
    WorldRasterGrid128  grid;
    MultiPassShader     shader;

    MandelPhaseData(int _phase) : 
        phase(_phase), 
        field(_phase),
        grid() 
    {}
};

struct Mandelbrot_Scene : public Scene<Mandelbrot_Scene>, public MandelState
{
    // ────── config ──────
    struct Config { std::string load_example_name; };
    Mandelbrot_Scene(Config& config) : load_example_name(config.load_example_name) {}
    std::string load_example_name;

    // ────── examples ──────
    MandelBookmarkManager bookmark_manager;
    bool rendering_examples = false;
    int rendering_example_i = 0;
    bool ignore_preset_filters = false;
    std::string render_batch_name;

    // ────── threads ──────
    static constexpr int MAX_THREADS = 0; // 0 = max threads

    // ────── minimum f32 precision ──────
    typedef f64 T_lo; // minimum f32 precision (f32 supported but slower, may be useful for GPU shading)

    // ────── cached compute data ──────
    RefOrbitLo<T_lo> ref_orbit_lo;
    //DDVec2 ref_orbit_anchor = DDVec2::highest();

    // ────── resources ──────
    NanoFont font;

    // ────── tweening / steady zoom ──────
    bool tweening = false;
    MandelState  state_a;
    MandelState  state_b;

    // ────── gradients ──────
    ImGradient gradient_shifted;

    // ────── shading ──────
    ToneManager<f64> dist_tone;
    ToneManager<f32> stripe_tone;

    // ────── saving & loading ──────
    std::string serializeState() const;
    void loadState(std::string_view data, bool is_big_jump = false);

    // ────── phases ──────
    int surface_w; /// = (ceil(w / 9)) * 9;
    int surface_h; /// = (ceil(h / 9)) * 9;

    MandelPhaseData phases[PHASE_COUNT] = {
        MandelPhaseData(PHASE_27X),
        MandelPhaseData(PHASE_9X),
        MandelPhaseData(PHASE_3X),
        MandelPhaseData(PHASE_1X)
    };

    // ────── fields for each phase (raw f32 data) ──────
    EscapeField& activeField() { return phases[active_phase].field; }
    EscapeField& pendingField() { return phases[computing_phase].field; }
    const EscapeField& activeField()  const { return phases[active_phase].field; }
    const EscapeField& pendingField() const { return phases[computing_phase].field; }


    // ────── raster grids for each phase (no RGBA data, just forEachPixel helpers) ──────
    WorldRasterGrid128& activeRasterGrid()  { return phases[active_phase].grid; }
    WorldRasterGrid128& pendingRasterGrid() { return phases[computing_phase].grid; }
    const WorldRasterGrid128& activeRasterGrid()  const { return phases[active_phase].grid; }
    const WorldRasterGrid128& pendingRasterGrid() const { return phases[computing_phase].grid; }


    // ────── shading ──────
    static inline const char* init_shader_source =
        "// @pass base\n"
        "return sampleGradient(wrap01(iter + dist + stripe));\n";

    bool update_editor_shader_source = false;

    GradientTexture gradient_tex;

    // don't sync histogram info unless visible
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
    //f32 stripe_mag_numerator = 0.1f;
    //bool  stripe_mag_from_hist = false;
    //std::map<f32, f32> ideal_zf_numerator_map;
    //std::vector<std::pair<f32, f32>> ideal_zf_numerator_map;
    //ImRect stripe_zf_spline_rect = ImRect(0.0f, 0.0f, 30.0f, 1.0f);
    
    // ────── stats ──────
    MandelStats stats;
    
    // ────── compute ──────
    int                iter_lim            = 0; // iteration limit provided to mandelbrot kernel
    KernelFeatures     mandel_features     = KernelFeatures::ITER;
    KernelMode         kernel_mode         = KernelMode::AUTO;
    KernelMode         deduced_kernel_mode = KernelMode::NO_PERTURBATION;
    FloatingPointType  float_type          = FloatingPointType::F64;

    // multi-frame tile progress tracker for the pending phase
    TileBlockProgress P;

    // phase 0 = 9x smaller field
    // phase 1 = 3x smaller field
    // phase 2 = full resolution field
    int  starting_phase = 0; // may start at later phase if performance is good
    int  computing_phase = 0;
    int  active_phase = 0;
    int  next_starting_phase = 0; // signal to start on this phase next full compute if performance is good

    bool display_intro = true;
    bool display_alignment_overlay = false;
    bool first_frame = true;
    bool frame_complete = false;  // Similar to finished_compute, but not cleared until next compute starts
    bool final_frame_complete = true;

    bool renormalize = false;
    bool reshade = false;

    DDQuad world_quad{};
    
    // ────── interior-forwarding optimization ──────
    struct ContractExpandPhase { int contract, expand; };
    ContractExpandPhase contract_expand_phases[PHASE_COUNT - 1] = {
        // phase [n-1] => phase [n]  // todo: Move each contract/expand pair to MandelPhase
        { 5, 0 }, 
        { 5, 3 },
        { 7, 3 }
    };

    bool maxdepth_show_optimized = false;

    // ────── timers ──────
    f64 dt_last = 0;

    SimpleTimer compute_timer;

    // compute time estimation
    std::array<f64, PHASE_COUNT>            phase_timers;
    std::array<math::SMA<f64>, PHASE_COUNT> phase_elapsed_mult_sma_list;
    std::array<f64, PHASE_COUNT>            phase_elapsed_mult_results;

    // todo: Move to MandelPhase
    math::SMA<f64>                          phase_elapsed_mult_sma = math::SMA<f64>(2);
    f64                                     phase_elapsed_mult_sma_result = 6.5;
    f64                                     phase_elapsted_estimated_final = 0;

    int compute_timeout = 0;

    SimpleTimer cam_change_timer;
    bool adjustingCamera() const { return (cam_change_timer.elapsed() < 500.0); }
    int  updateComputeTimeout();

    // ────── camera navigation easing ──────
    CameraNavigator navigator;

    math::SMA<DDVec2> avg_vel_pos = math::SMA<DDVec2>(8);
    math::SMA<f64> avg_vel_zoom = math::SMA<f64>(8);

    DDVec2 cam_vel_pos{};
    f64 cam_vel_zoom{1};

    // ────── deep-zoom animation (temporary until timeline feature added) ──────
    bool   steady_zoom = false;
    bool   record_steady_zoom = false;
    f64    steady_zoom_mult_speed = 0.01;
    f64    steady_zoom_rotate_speed = 0.0;
    f32    steady_zoom_pct = 0;
    int    steady_zoom_frames_elapsed = 0;
    int    steady_zoom_expected_frames = 0;
    math::SMA<f64> expected_time_left_ma = math::SMA<f64>(5); // todo: Use or remove

    // use same mean stream during deep zoom. todo: Issue if phase is animated? Or does mean change
    f32  stripe_mean_locked = 0.0;

    // ────── undo/redo ──────
    const int history_count = 100;
    bool pending_checkpoint_flag = false;
    bool loading_history_state = false;
    std::vector<std::pair<int, std::string>> history;
    int current_history_i = 0;

    // ────── cardioid ──────
    #if MANDEL_FEATURE_INTERACTIVE_CARDIOID
    bool   show_interactive_cardioid = false;
    bool   animate_cardioid_angle = true;
    f64 ani_angle = 0.0;
    f64 ani_inc = math::toRadians(2.0);
    #endif

    // ────── user interface ──────
    struct UI : BufferedInterfaceModel
    {
        using BufferedInterfaceModel::BufferedInterfaceModel;

        // dialog flags
        bool show_save_dialog  = false;
        bool show_load_dialog  = false;
        bool show_share_dialog = false;
        
        // histogram visibilities
        bool show_iter_hist    = false;
        bool show_dist_hist    = false;
        bool show_stripe_hist  = false;

        // textbox buffers
        std::string url_buf;
        std::string data_buf;

        // preset filters
        int selected_preset_i = 0;

        // shader editor
        TextEditor  editor;
        f32       editor_height;
        bool        combined_vert_errors = false;
        bool        combined_frag_errors = false;
        std::string combined_vert_log = "";
        std::string combined_frag_log = "";
        bool        auto_apply_shader = true;

        // mouse orbit
        f64 orbit_padding = 0.3;

        // UI init
        void init() override;
        void overlay() override;
        void sidebar() override;

        void launchBookmark(std::string_view data);

        // UI sections
        void populateSavingLoading();
        void populateExamples();
        void populateCameraView();
        void populateQualityOptions();
        void populateInputOptions();
        void populateShaderEditor();
        void populateGradientOptions();
        void populateAnimation();
        void populateStats();
        void populateCaptureOptions();

        void populateMouseOrbit();
        void populateHistory();

        #if MANDEL_UNVERSIONED_EXPERIMENTAL
        void populateExperimental();
        void populateSplinesDev();
        #endif
    };

    // ────── simulation processing ──────
    void sceneStart() override;
    void sceneDestroy() override;
    void sceneMounted(Viewport* viewport) override;

    // ────── compute ──────
    template<typename T, KernelFeatures F, KernelMode K>
    bool compute_mandelbrot(int timeout, f64& field_compute_ms);

    template<typename T, KernelFeatures F>
    void calculate_normalize_info();

    template<typename T, KernelFeatures F, bool Normalize_Depth, bool Invert_Dist, bool Show_Optimized>
    void normalize_field();


    // ────── viewport handling ──────
    bool mandelChanged();
    bool normalizationOptionsChanged();

    void beginSteadyZoom();
    void updateSteadyZoom(f64 dt);
    void endSteadyZoom();

    // process pipeline
    void loadNextBatchSnapshotState();
    void updateAnimation();
    void updateCameraView(Viewport* ctx);
    void updateQuality();
    void resetCompute();
    void updateFieldSizes();
    void updateEnabledKernelFeatures();
    bool processCompute();
    bool updateGradient();
    void processCapturing(bool finished_compute, bool reshade);
    void processUndoRedo(bool normalization_opts_changed, bool gradient_changed);

    #if MANDEL_UNVERSIONED_EXPERIMENTAL
    Mandel_Experimental experimental;
    #endif

    // processing (worker)
    void viewportProcess(Viewport* ctx, f64 dt) override;
    void collectStats(bool renormalized);

    // drawing (GUI)
    void renderShaderChain(Viewport* ctx) const;
    void viewportDraw(Viewport* ctx) const override;
    void onEncodeFrame(EncodeFrame& frame, const CapturePreset& preset) override;

    // ────── input ──────
    void onEvent(Event e) override;
    void onKeyDown(KeyEvent e) override;

    void historyOnKeyDown(KeyEvent e);
};

struct Mandelbrot_Project : public BasicProject
{
    static ProjectInfo info() {
        return ProjectInfo({ "Fractal", "Mandelbrot", "Mandelbrot Viewer" });
    }

    void projectPrepare(Layout& layout) override;
    bool sidebarVisible() const override { return false; }
};

SIM_END;

