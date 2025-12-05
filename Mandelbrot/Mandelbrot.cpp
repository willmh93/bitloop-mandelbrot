
#include "Mandelbrot.h"
#include "mandel_process.hpp"

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

    // todo: Stop this getting called twice on startup
    generateGradientFromPreset(gradient, GradientPreset::CLASSIC);

    cardioid_lerper.create(Math::TWO_PI / 5760.0, 0.005);

    font = NanoFont::create("/data/fonts/DroidSans.ttf");
}

void Mandelbrot_Scene::sceneDestroy()
{
    blPrint() << "Mandelbrot_Scene::sceneDestroy()";
}

void Mandelbrot_Scene::sceneMounted(Viewport* ctx)
{
    blPrint() << "Mandelbrot_Scene::sceneMounted()";

    //camera->setCameraStageSnappingSize(1);
    camera.setSurface(ctx);
    camera.setOriginViewportAnchor(Anchor::CENTER);
    camera.focusWorldRect(-2, -1.25, 1, 1.25);
    camera.uiSetCurrentAsDefault();

    bmp_9x9.setCamera(camera);
    bmp_3x3.setCamera(camera);
    bmp_1x1.setCamera(camera);

    navigator.setTarget(camera);
    navigator.setDirectCameraPanning(true);
    //camera->restrictRelativeZoomRange(0.5, 1e+300);

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
}

void Mandelbrot_Scene::viewportProcess(Viewport* ctx, double dt)
{
    /// Never record a frame unless we finish computing a full frame
    captureFrame(false);

    // ────── automatic updates (animation / tweening) ──────
    updateAnimation();
    updateTweening(dt);

    // ────── update states needed for compute ──────
    updateCameraView();
    updateFieldSizes(ctx);
    updateEnabledKernelFeatures();
    updateActivePhaseAndField();

    // ────── compute ──────
    bool finished_compute = false;
    if (!final_frame_complete)
    {
        // compute the unnormalized base data needed for the normalization/reshading stages
        finished_compute = processCompute();
    }

    // ────── color cycle changed? ──────
    bool gradient_changed = updateGradient();
    bool shading_formula_changed = shadingFormulaChanged();

    bool renormalize = finished_compute || shading_formula_changed;
    bool reshade = renormalize || (gradient_changed && frame_complete);


    // ────── renormalize cached data if necessary  ──────
    if (renormalize)
    {
        FloatingPointType float_type = getRequiredFloatType(mandel_features, camera.relativeZoom<f128>());
        bool normalize_depth = iter_params.cycle_iter_normalize_depth;
        bool invert_dist = dist_params.cycle_dist_invert;

        table_invoke(dispatch_table(calculate_normalize_info), float_type, mandel_features, normalize_depth, invert_dist);
        table_invoke(dispatch_table(normalize_field),          float_type, mandel_features, normalize_depth, invert_dist);
    }

    // ────── reshade if data renormalized, or gradient changed ──────
    if (reshade)
    {
        table_invoke(dispatch_table(shadeBitmap), (MandelShaderFormula)shade_formula, maxdepth_show_optimized);
    }

    // unless we're doing a steady-zoom animation, capture basic animated shading every frame
    if (!steady_zoom && reshade && final_frame_complete)
        captureFrame(true);

    first_frame = false;

    if (!steady_zoom && show_interactive_cardioid && animate_cardioid_angle)
    {
        // Cardioid animation
        ani_angle += ani_inc;
        ani_angle = Math::wrapRadians2PI(ani_angle);

        // If not tweening and we're animating the cardioid -> record frame
        if (isRecording())
            captureFrame(true);
    }

    // Gather stats / realtime info
    collectStats(renormalize);

}

void Mandelbrot_Scene::viewportDraw(Viewport* ctx) const
{
    ctx->setTransform(camera.getTransform());

    // Draw active phase bitmap
    if (active_bmp)
    {
        ///ctx->setGlobalAlpha(0.5);
        ctx->drawImage(*active_bmp, active_bmp->worldQuad());
        ///ctx->setGlobalAlpha(1.0);
    }

    double zoom_mag = camera.relativeZoom<double>();

    if (show_axis && zoom_mag < 1.0e7)
        ctx->drawWorldAxis(0.5, 0, 0.5);

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
            /// Cardioid::animatePlot(this, ctx, 1.0, 0.0, Math::PI*m);
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

SIM_END;
