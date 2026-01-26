
#include "Mandelbrot.h"
#include "mandel_process.hpp"
#include "history.hpp"
#include "experimental.h"

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

    // todo: Stop this getting called twice on startup
    generateGradientFromPreset(gradient, GradientPreset::CLASSIC);

    cardioid_lerper.create(math::tau / 5760.0, 0.005);

    font = NanoFont::create("/data/fonts/DroidSans.ttf");

    createShaderPasses(init_shader_source, pass_surfaces);
}

void Mandelbrot_Scene::sceneDestroy()
{
    blPrint() << "Mandelbrot_Scene::sceneDestroy()";

    main_window()->deferredGuiDestructionQueue().enqueue([&]() noexcept
    {
        if (black1x1Tex != 0) 
            glDeleteTextures(1, &black1x1Tex);

        field_1x1.destroyFeaturesTexture();
        field_3x3.destroyFeaturesTexture();
        field_9x9.destroyFeaturesTexture();

        gradient_tex.destroy();
    });
}

void Mandelbrot_Scene::sceneMounted(Viewport* ctx)
{
    blPrint() << "Mandelbrot_Scene::sceneMounted()";

    //camera->setCameraStageSnappingSize(1);
    camera.setSurface(ctx);
    camera.setOriginViewportAnchor(Anchor::CENTER);
    camera.focusWorldRect(-2.25, -1.25, 0.75, 1.25);
    camera.uiSetCurrentAsDefault();

    bmp_9x9.setCamera(camera);
    bmp_3x3.setCamera(camera);
    bmp_1x1.setCamera(camera);

    navigator.setTarget(camera);
    navigator.setDirectCameraPanning(true);
    navigator.restrictRelativeZoomRange(0.0001, 5.0e28);

    #ifdef __EMSCRIPTEN__
    // If URL has encoded state data, load on startup
    if (platform()->url_has("data"))
        loadState(platform()->url_get_string("data"));
    #endif

    /*if (!load_example_name.empty())
    {
        data_buf = mandel_presets[load_example_name];
        if (!data_buf.empty())
            loadConfigBuffer();
    }*/

    history.push_back({ 0, serialize() });
}

void Mandelbrot_Scene::viewportProcess(Viewport* ctx, double dt)
{
    // ────── disable capture this frame, re-enable depending on mode ──────
    permitCaptureFrame(false);

    // if "Render All" was clicked, call beginSnapshot() for each preset until done
    processBatchSnapshot();
    
    // ────── automatic updates (animation / tweening) ──────
    updateAnimation();
    updateTweening(dt);

    // ────── update states needed for compute ──────
    updateCameraView();
    updateFieldSizes(ctx);
    updateEnabledKernelFeatures();
    updateActivePhaseAndField();

    // ────── compute unnormalized field values (ITER, DIST, STRIPE) ──────
    bool finished_compute = false;
    if (!final_frame_complete)
    {
        finished_compute = processCompute();
    }
    
    // ────── check if renormalizing / reshading ──────
    bool gradient_changed = updateGradient();
    bool normalization_opts_changed = normalizationOptionsChanged();

    renormalize = finished_compute || normalization_opts_changed;
    reshade = renormalize || (gradient_changed && frame_complete);

    // ────── renormalize if flagged  ──────
    if (renormalize)
    {
        FloatingPointType float_type = getRequiredFloatType(mandel_features, camera.relativeZoom<f128>());
        bool normalize_depth = iter_params.cycle_iter_normalize_depth;
        bool invert_dist = dist_params.cycle_dist_invert;

        table_invoke(dispatch_table(calculate_normalize_info), float_type, mandel_features);
        table_invoke(dispatch_table(normalize_field),          float_type, mandel_features, normalize_depth, invert_dist, maxdepth_show_optimized);
    }

    // ────── reshade if flagged ──────
    if (Changed(shader_source_txt))
    {
        createShaderPasses(shader_source_txt, pass_surfaces);
        reshade = true;
    }

    // check if we should permit frame capture on next draw (snapshot/record) - condition varies depending on mode
    processCapturing(finished_compute, reshade);

    // undo/redo logic
    if (!platform()->is_mobile())
        processUndoRedo(normalization_opts_changed, gradient_changed);

    // Gather stats / realtime info
    collectStats(renormalize);

    #if MANDEL_EXPERIMENTAL_TESTS
    processExperimental();
    #endif

    first_frame = false;
}

void Mandelbrot_Scene::renderShaderChain(Viewport* ctx) const
{
    GLuint prevTex = ensureBlack1x1Tex();
    FVec2 inv_size = FVec2(1.0f,1.0f) / (FVec2)ctx->outputSize();

    for (int pass = 0; pass < pass_surfaces.size(); pass++)
    {
        const ShaderSurface* surface = getPassSurfaces(pass)->getSurface(active_phase);
        surface->render(active_bmp->width(), active_bmp->height(), [&](const ShaderSurface& s)
        {
            s.setUniform2f("u_outTexel", inv_size.x, inv_size.y);
            s.bindTexture2D("u_prev", prevTex, 0);
            s.bindTexture2D("u_features", active_field->features_tex, 1);
            s.bindTexture2D("u_gradient", gradient_tex.getTexture(), 2);
        });
        prevTex = surface->texture();
    }
}

void Mandelbrot_Scene::viewportDraw(Viewport* ctx) const
{
    if (renormalize)
        active_field->updateFeaturesTexture();

    if (reshade)
        renderShaderChain(ctx);

    // apply 128-bit camera world transformation
    ctx->setTransform(camera.getTransform());

    const ShaderSurface* composite_surface = getPassSurfaces((int)pass_surfaces.size() - 1)->getSurface(active_phase);

    // Draw active phase bitmap
    if (composite_surface)
    {
        ///ctx->setGlobalAlpha(0.5);
        ctx->drawShaderSurface(*composite_surface, 0.0f, 0.0f, (f32)ctx->width(), (f32)ctx->height());
        ///ctx->setGlobalAlpha(1.0);
    }

    double zoom_mag = camera.relativeZoom<double>();

    if (show_axis && zoom_mag < 1.0e7)
        ctx->drawWorldAxis(0.5, 0, 0.5);

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

    #if MANDEL_EXPERIMENTAL_TESTS
    drawExperimental(ctx);
    #endif
    
    #if MANDEL_FEATURE_INTERACTIVE_CARDIOID
    if (show_interactive_cardioid && zoom_mag < 1000000.0)
    {
        if (!flatten)
        {
            // interactive cardioid
            Cardioid::plot(this, ctx, false);
            Cardioid::animatePlot(ctx, 1.0, 0.0, ani_angle);

            /// Interesting...
            /// for (double m=0; m<=1.0; m+=0.01)
            /// Cardioid::animatePlot(this, ctx, 1.0, 0.0, math::pi*m);
        }
        else
        {
            // Lerped/Flattened Mandelbrot
            ctx->scalingLines(false);
            ctx->setLineWidth(1);
            ctx->beginPath();
            ctx->drawPath(cardioid_lerper.lerped(cardioid_lerp_amount));
            ctx->stroke();
        }
    }
    #endif

    ctx->stageMode();

    #ifndef BL_DEBUG
    if (display_intro)
    {
        // intro message on startup
        ctx->setFont(font);
        ctx->setFontSize(20);

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

void Mandelbrot_Scene::onEvent(Event e)
{
    if (!this->ownsEvent(e))
        return;
    
    #ifdef BL_RELEASE
    // don't permit world navigation while tweening
    if (tweening || isRecording())
        return;
    #endif

    // enter/exit fullscreen hotkeys (F11 / Escape)
    if (e.type() == SDL_EVENT_KEY_DOWN)
    {
        if (e.sdl()->key.key == SDLK_F11)
        {
            SDL_WindowFlags flags = SDL_GetWindowFlags(platform()->sdl_window());
            bool fullscreen = (flags & SDL_WINDOW_FULLSCREEN);

            SDL_SetWindowFullscreen(platform()->sdl_window(), !fullscreen);
        }
        else if (e.sdl()->key.key == SDLK_ESCAPE)
        {
            SDL_SetWindowFullscreen(platform()->sdl_window(), false);
        }
    }

    
    DDVec2 old_pos = camera.pos<f128>();
    f128 old_zoom = camera.relativeZoom<f128>();
    
    if (navigator.handleWorldNavigation(e, true, true))
    {
        f128 zoom = camera.relativeZoom<f128>();
        avg_vel_pos.push(camera.pos<f128>() - old_pos);
        avg_vel_zoom.push( (double)(zoom / old_zoom) );

        camera_vel_pos = avg_vel_pos.average();
        camera_vel_zoom = (avg_vel_zoom.average() - 1) * 0.6 + 1;

        //if (!platform()->is_mobile()) // Stopped wheeling mouse
            avg_vel_zoom.clear();
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
    historyKeyEvent(e);
}

SIM_END;
