
#include "Mandelbrot.h"

#include "core/mandel_process.hpp"
#include "ui/history.hpp"

#if MANDEL_UNVERSIONED_EXPERIMENTAL
#include "_unversioned/experimental.hpp"
#endif

SIM_BEG;

/// ─────────────────────── Project ───────────────────────

void Mandelbrot_Project::projectPrepare(Layout& layout)
{
    Mandelbrot_Scene::Config config;
    create<Mandelbrot_Scene>(config)->mountTo(layout);
}

/// ─────────────────────── Scene ───────────────────────

void Mandelbrot_Scene::sceneStart()
{
    Thread::setMaxThreads(MAX_THREADS);

    SnapshotPresetList& all_presets = main_window()->getSnapshotPresetManager()->allPresets();
    valid_presets.clear();
    for (const auto& preset : all_presets)
        valid_presets[preset.hashedAlias()] = true;

    generateGradientFromPreset(gradient, GradientPreset::CLASSIC);

    font = NanoFont::create("/data/fonts/DroidSans.ttf");

    for (int phase = 0; phase < PHASE_COUNT; phase++)
        phases[phase].shader.updateFragmentSource(init_shader_source);

    phase_elapsed_mult_sma.setLength(2);
    phase_elapsed_mult_sma_result = 6.5;

    for (int phase = 0; phase < PHASE_COUNT; phase++)
    {
        auto& timer = phase_elapsed_mult_sma_list[phase];
        timer.setLength(3);

        phase_elapsed_mult_results[phase] = 6.5;
    }
}

void Mandelbrot_Scene::sceneDestroy()
{
    main_window()->threadQueue().post([&]()
    {
        for (MandelPhaseData& phase : phases)
            phase.field.destroyFeaturesTexture();

        gradient_tex.destroy();
        bookmark_manager.destroyAllTextures();
    });
}

void Mandelbrot_Scene::sceneMounted(Viewport* ctx)
{
    //camera->setCameraStageSnappingSize(1);
    camera.setSurface(ctx);
    camera.setOriginViewportAnchor(Anchor::CENTER);
    camera.focusWorldRect(-2.25, -1.25, 0.75, 1.25);
    camera.uiSetCurrentAsDefault();

    for (int phase = 0; phase < PHASE_COUNT; phase++)
        phases[phase].grid.setCamera(camera);

    navigator.setTarget(camera);
    navigator.setDirectCameraPanning(true);
    navigator.restrictRelativeZoomRange(0.0001, 5.0e28);

    #ifdef __EMSCRIPTEN__
    // If URL has encoded state data, load on startup
    if (platform()->url_has("data"))
        loadState(platform()->url_get_string("data"));
    #endif

    history.push_back({ 0, serialize() });
}

std::string Mandelbrot_Scene::serializeState() const 
{
    return MandelState::serialize(); 
}

void Mandelbrot_Scene::loadState(std::string_view data, bool is_big_jump)
{
    deserialize(data);
    update_editor_shader_source = true;

    // on a drastic change, assume the worst about the difficulty of the target state,
    // and restart on lowest resolution phase
    if (is_big_jump)
    {
        next_starting_phase = 0;

        for (int phase = 0; phase < PHASE_COUNT; phase++)
        {
            auto& timer = phase_elapsed_mult_sma_list[phase];
            timer.clear();
            phase_elapsed_mult_results[phase] = 6.5;
        }
    }
}

int Mandelbrot_Scene::updateComputeTimeout()
{
    if (steady_zoom || computing_phase == starting_phase)
        // always one-shot the first phase, or the entire compute if steady-zoom (fixed resolution + interaction disabled)
        compute_timeout = 0;

    else if (adjustingCamera() || mouse->buttonClicked(MouseButton::LEFT))
        // smooth panning / zooming
        compute_timeout = 128;

    else
    {
        // potential multi-frame compute, optimize for balance between speed / responsiveness
        f64 expected_next_phase_dt = phase_timers[std::max(0, computing_phase - 1)] * 8.0;
        if (expected_next_phase_dt > 256.0 || float_type == FloatingPointType::F128)
            compute_timeout = 256; // slightly larger timeout for difficult frames = faster compute, more interaction lag
        else
            compute_timeout = 512; // enough to one-shot easy phases < 256ms expected (timeout as a last resort)
    }
    return compute_timeout;
}

void Mandelbrot_Scene::viewportProcess(Viewport* ctx, double dt)
{
    // ────── disable capture this frame, re-enable depending on mode ──────
    permitCaptureFrame(false);

    if (rendering_examples)
        loadNextBatchSnapshotState();
    
    // ────── automatic updates (animation / tweening) ──────
    updateAnimation();
    updateSteadyZoom(dt);

    // ────── update states which can determine whether we need to recompute ──────
    updateCameraView(ctx);
    updateQuality();
    updateEnabledKernelFeatures();

    if (mandelChanged())
        resetCompute();

    updateFieldSizes();

    // ────── compute unnormalized field values (ITER, DIST, STRIPE) ──────
    bool finished_compute = false;
    if (!final_frame_complete)
        finished_compute = processCompute();
    
    // ────── check if renormalizing / reshading ──────
    bool gradient_changed = updateGradient();
    bool normalization_opts_changed = normalizationOptionsChanged();

    renormalize = finished_compute || normalization_opts_changed;
    reshade = renormalize || (gradient_changed && frame_complete);

    // ────── renormalize if flagged  ──────
    if (renormalize)
    {
        bool normalize_depth = iter_params.iter_normalize_depth;
        bool invert_dist = dist_params.dist_invert;

        table_invoke(dispatch_table(calculate_normalize_info), float_type, mandel_features);
        table_invoke(dispatch_table(normalize_field),          float_type, mandel_features, normalize_depth, invert_dist, maxdepth_show_optimized);
    }

    // ────── reshade on shader source change ──────
    // todo: fine for short scripts, but check flag for larger script
    if (Changed(shader_source_txt))
    {
        for (int phase = 0; phase < PHASE_COUNT; phase++)
            phases[phase].shader.updateFragmentSource(shader_source_txt);
        
        reshade = true;
    }

    // ────── check if we should permit frame capture on next draw (snapshot/record) ──────
    processCapturing(finished_compute, reshade);

    // undo/redo logic
    if (!platform()->is_mobile())
        processUndoRedo(normalization_opts_changed, gradient_changed);

    // Gather stats / realtime info
    collectStats(renormalize);

    #if MANDEL_UNVERSIONED_EXPERIMENTAL
    processExperimental();
    #endif

    first_frame = false;

    bool toggled_overelays = Changed(show_axis, display_alignment_overlay, show_interactive_cardioid, preview_normalization_field);
    if (reshade || toggled_overelays)
        requestRedraw(true);

    if (mouse->buttonDown(MouseButton::LEFT) || cam_vel_zoom != 1 || !cam_vel_pos.isZero())
        requestImmediateUpdate();
}

void Mandelbrot_Scene::renderShaderChain(Viewport* ctx) const
{
    // todo: Not right when you do high SSAA renders?
    FVec2 inv_size = FVec2(1.0f, 1.0f) / (FVec2)ctx->outputSize();
    
    const EscapeField& field = activeField();
    const WorldRasterGrid128& grid = activeRasterGrid();

    phases[active_phase].shader.render(grid.rasterWidth(), grid.rasterHeight(), [&](const ShaderSurface& s)
    {
        s.setUniform1f("min_iter",   field.final_min_depth);
        s.setUniform1f("max_iter",   field.final_max_depth);
        s.setUniform1f("min_dist",   field.final_min_dist);
        s.setUniform1f("max_dist",   field.final_max_dist);
        s.setUniform1f("min_stripe", field.final_min_stripe);
        s.setUniform1f("max_stripe", field.final_max_stripe);

        //s.setUniform1f("u_time", (f32)running_dt());
        s.setUniform2f("u_outTexel", inv_size.x, inv_size.y);
        s.bindTexture2D("u_features", field.getTexture());
        s.bindTexture2D("u_gradient", gradient_tex.getTexture());
    });
}

void Mandelbrot_Scene::viewportDraw(Viewport* ctx) const
{
    if (reshade)
        renderShaderChain(ctx);

    // apply 128-bit camera world transformation
    ctx->setTransform(camera.getTransform());

    const ShaderSurface* composite_surface = phases[active_phase].shader.outputSurface();

    // Draw active phase bitmap
    ///ctx->setGlobalAlpha(0.5);
    ctx->drawShaderSurface(*composite_surface, 0.0f, 0.0f, (f32)surface_w, (f32)surface_h);
    ///ctx->setGlobalAlpha(1.0);
    
    double zoom_mag = camera.relativeZoom<f64>();

    if (show_axis && zoom_mag < 1.0e7)
        ctx->drawWorldAxis(0.5, 0, 0.5);

    if (display_alignment_overlay)
    {
        ctx->setStrokeStyle(Color::white);
        ctx->setLineWidth(1);
        ctx->stageMode();
        ctx->beginPath();
        ctx->moveTo(ctx->width() / 2.0, 0.0);
        ctx->lineTo(ctx->width() / 2.0, ctx->height());
        ctx->moveTo(0.0, ctx->height() / 2.0);
        ctx->lineTo(ctx->width(), ctx->height() / 2.0);
        ctx->stroke();
    }

    // if "Normalization Sampling" -> "Preview" ticked
    if (preview_normalization_field)
    {
        ctx->stageMode();
        ctx->setFillStyle(255, 255, 255, 70);
        for (auto& p : norm_field.world_field)
        {
            float rad = std::max(1.0f, p.weight * 4.0f);
            ctx->fillEllipse<f128>(p.stage_pos, rad);
        }
		ctx->worldMode();
    }

    #if MANDEL_UNVERSIONED_EXPERIMENTAL
    drawExperimental(ctx);
    #endif
    
    #if MANDEL_FEATURE_INTERACTIVE_CARDIOID
    if (show_interactive_cardioid && zoom_mag < 1000000.0)
    {
        // interactive cardioid
        Cardioid::plot(this, ctx, false);
        Cardioid::animatePlot(ctx, 1.0, 0.0, ani_angle);
    }
    #endif

    ctx->worldHudMode();

    ctx->stageMode();

    #ifndef BL_DEBUG
    if (display_intro)
    {
        // intro message on startup
        ctx->setFont(font);
        ctx->setFontSize(20);
        ctx->setFillStyle(Color::white);

        if (platform()->is_mobile())
        {
            // touchscreen instructions
            ctx->fillText("Controls:", scale_size(10.0), scale_size(10.0));
            ctx->fillText("  - Touch & drag to move", scale_size(10.0), scale_size(35.0));
            ctx->fillText("  - Pinch to zoom / rotate", scale_size(10.0), scale_size(60.0));
        }
        else
        {
            // desktop instructions
            ctx->fillText("Controls:", scale_size(10.0), scale_size(10.0));
            ctx->fillText("  - Click & drag to move", scale_size(10.0), scale_size(35.0));
            ctx->fillText("  - Mouse wheel to zoom", scale_size(10.0), scale_size(60.0));
        }

        ctx->setTextAlign(TextAlign::ALIGN_LEFT);
        ctx->setTextBaseline(TextBaseline::BASELINE_BOTTOM);
        ctx->fillText("Contact:  will.hemsworth@bitloop.dev", scale_size(10.0), ctx->height() - scale_size(10.0));
    }
    #endif
}

void Mandelbrot_Scene::onEncodeFrame(EncodeFrame& frame, const CapturePreset& preset)
{
    // attach XMP metadata to recover MandelState from any snapshot
    if (preset.isImage())
        frame.payload = serializeState();
}

void Mandelbrot_Scene::onEvent(Event e)
{
    if (!e.ownedBy(this))
        return;
    
    // don't permit world navigation while tweening
    if (tweening || isRecording())
        return;

    //if (e.sdl()->type == SDL_EVENT_DROP_FILE)
    //{
    //    std::string path = e.sdl()->drop.data; // file path (UTF-8)
    //    blPrint() << "Mandelbrot: " << path;
    //}

    DDVec2 old_pos = camera.pos<f128>();
    f128 old_zoom = camera.relativeZoom<f128>();
    
    if (navigator.handleWorldNavigation(e, true, true))
    {
        // camera altered - calculate velocities by comparing with previous frame, smooth with SMA
        f128 zoom = camera.relativeZoom<f128>();
        avg_vel_pos.push(camera.pos<f128>() - old_pos);
        avg_vel_zoom.push( (double)(zoom / old_zoom) );

        cam_vel_pos = avg_vel_pos.average();
        cam_vel_zoom = (avg_vel_zoom.average() - 1) * 0.6 + 1;

        avg_vel_zoom.clear();

        cam_change_timer.begin();
    }

    if (e.type() == SDL_EVENT_MOUSE_BUTTON_UP ||
        e.type() == SDL_EVENT_FINGER_UP)
    {
        avg_vel_zoom.clear();
        avg_vel_pos.clear();
    }
}

void Mandelbrot_Scene::onKeyDown(KeyEvent e)
{
    historyOnKeyDown(e);
}

SIM_END;
