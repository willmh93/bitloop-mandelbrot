#include "Mandelbrot.h"

SIM_BEG;

void Mandelbrot_Scene::updateAnimation()
{
    //if (Changed(show_color_animation_options, show_axis, gradient_shift_step, hue_shift_step))
    //    savefile_changed = true;

    double ani_mult = fpsFactor();

    // ────── Progressing animation ──────
    if (!steady_zoom || capturedLastFrame())
    {
        if (show_color_animation_options)
        {
            // Animation
            if (fabs(gradient_shift_step) > 1.0e-4)
                gradient_shift = Math::wrap(gradient_shift + gradient_shift_step * ani_mult, 0.0, 1.0);

            if (fabs(hue_shift_step) > 1.0e-4)
                hue_shift = Math::wrap(hue_shift + hue_shift_step * ani_mult, 0.0, 360.0);
        }
        else
        {
            //if (Changed(gradient_shift, hue_shift, gradient_shift_step, hue_shift_step))
            //    savefile_changed = true;
        }
    }
}

bool Mandelbrot_Scene::updateGradient()
{
    if (first_frame || Changed(gradient, gradient_shift, hue_shift))
    {
        transformGradient(gradient_shifted, gradient, (float)gradient_shift, (float)hue_shift);
        //savefile_changed = true;
        ///colors_updated = true;
        return true;
    }
    return false;
}

void Mandelbrot_Scene::updateCameraView()
{
    #if MANDEL_FEATURE_FLATTEN_MODE
    if (Changed(flatten))
    {
        if (flatten) flatten_amount = 0.0;
        camera->focusWorldRect(-2, -1, 2, 1);
        camera.zoom = camera->relativeZoom<f128>().x;
    }

    if (Changed(flatten_amount))
    {
        using namespace Math;
        double t = flatten_amount;
        DRect r;

        if (t < 0.5)      r = lerp(DRect(-2.0, -1.5, 0.5, 1.5), DRect(-2.0, -0.2, 1.5, 3.5), lerpFactor(t, 0.0, 0.5));
        else if (t < 0.7) r = lerp(DRect(-2.0, -0.2, 1.5, 3.5), DRect(-1.5, -0.2, 4.0, 3.5), lerpFactor(t, 0.5, 0.7));
        else              r = lerp(DRect(-1.5, -0.2, 4.0, 3.5), DRect(0.0, -1.5, 4.5, 0.5), lerpFactor(t, 0.7, 1.0));

        camera->focusWorldRect(r, false);
        camera.zoom = camera->relativeZoom<f128>().x;
    }
    #endif

    if (!tweening) {

        #if MANDEL_FEATURE_CAMERA_EASING
        // Process camera velocity
        if (mouse->pressed)
        {
            camera_vel_pos = DDVec2{};

            // Stop zoom velocity on touch
            camera_vel_zoom = 1.0;
        }
        else
        {
            f128 threshold = camera.relativeZoom<f128>() / 1000.0;
            if (camera_vel_pos.magnitude() > threshold)
                camera.setPos(camera.pos<f128>() + camera_vel_pos);
            else
                camera_vel_pos = DDVec2{};

            if (fabs(camera_vel_zoom - 1.0) > 0.01)
                camera.setRelativeZoom(camera.relativeZoom<f128>() * camera_vel_zoom);
            else
                camera_vel_zoom = 1.0;

            camera_vel_pos *= 0.8;
            camera_vel_zoom += (1.0 - camera_vel_zoom) * 0.2; // Ease back to 1x
        }
        #endif
    }
}

void Mandelbrot_Scene::updateFieldSizes(Viewport* ctx)
{
    int iw = (int)(ceil(ctx->width() / 9)) * 9;
    int ih = (int)(ceil(ctx->height() / 9)) * 9;

    // set field sizes
    field_9x9.setDimensions(iw / 9, ih / 9);
    field_3x3.setDimensions(iw / 3, ih / 3);
    field_1x1.setDimensions(iw, ih);

    // set bitmap sizes
    bmp_9x9.setBitmapSize(iw / 9, ih / 9);
    bmp_3x3.setBitmapSize(iw / 3, ih / 3);
    bmp_1x1.setBitmapSize(iw, ih);

    // determine world quad from viewport rect
    bmp_9x9.setStageRect(0, 0, iw, ih);
    bmp_3x3.setStageRect(0, 0, iw, ih);
    bmp_1x1.setStageRect(0, 0, iw, ih);

    world_quad = bmp_1x1.worldQuad();
}

void Mandelbrot_Scene::updateEnabledKernelFeatures()
{
    bool mandel_changed = mandelChanged();

    constexpr double eps = std::numeric_limits<double>::epsilon();
    MandelKernelFeatures old_smoothing = mandel_features;
    MandelKernelFeatures new_smoothing = MandelKernelFeatures::NONE;

    if (iter_weight > eps)   new_smoothing |= MandelKernelFeatures::ITER;
    if (dist_weight > eps)   new_smoothing |= MandelKernelFeatures::DIST;
    if (stripe_weight > eps) new_smoothing |= MandelKernelFeatures::STRIPES;

    bool force_upgrade       = (bool)( new_smoothing & ~old_smoothing);
    bool downgrade_on_change = (bool)(~new_smoothing &  old_smoothing);

    if (force_upgrade || (mandel_changed && downgrade_on_change))
    {
        mandel_changed = true;
        mandel_features = new_smoothing;
    }
}

void Mandelbrot_Scene::updateActivePhaseAndField()
{
    bool mandel_changed = mandelChanged();

    // ────── presented mandelbrot changed? restart at phase 0 ──────
    if (mandel_changed)
    {
        bmp_9x9.clear(0, 255, 0, 255);
        bmp_3x3.clear(0, 255, 0, 255);
        bmp_1x1.clear(0, 255, 0, 255);

        // hide intro on first change
        if (!first_frame) display_intro = false;

        computing_phase = 0;
        current_row = 0;
        final_frame_complete = false;
        field_9x9.setAllDepth(-1.0);

        compute_t0 = std::chrono::steady_clock::now();
    }

    // ────── set pending bitmap & pending field ──────
    {
        switch (computing_phase)
        {
        case 0:  pending_bmp = &bmp_9x9;  pending_field = &field_9x9;  break;
        case 1:  pending_bmp = &bmp_3x3;  pending_field = &field_3x3;  break;
        case 2:  pending_bmp = &bmp_1x1;  pending_field = &field_1x1;  break;
        }
    }
}

bool Mandelbrot_Scene::processCompute()
{
    bool finished_compute = false;

    iter_lim = finalIterLimit(camera, quality, dynamic_iter_lim, tweening);

    FloatingPointType float_type = getRequiredFloatType(mandel_features, camera.relativeZoom<f128>());

    // Calculate first low-res phase in one-shot (no timeout)
    // For high-res phases, break up work across multiple frames if necessary (kept track of with current_row)
    int timeout = computing_phase == 0 ? 0 : 16;

    //if (!flatten)
    {
        //DQuad quad = ctx->worldQuad();
        //bool x_axis_visible = quad.intersects({ {quad.minX(), 0}, {quad.maxX(), 0} });

        // Run appropriate kernel for given settings
        finished_compute = frame_complete = table_invoke(
            build_table(mandelbrot, [&], pending_bmp, pending_field, iter_lim, numThreads(), timeout, current_row, stripe_params),
            float_type, mandel_features, flatten
        );
    }
    //else
    //{
    //    // Flat lerp
    //    finished_compute = dispatchBooleans(
    //        boolsTemplate(radialMandelbrot, [&]),
    //        mandel_features != MandelSmoothing::NONE,
    //        show_period2_bulb
    //    );
    //}

    if (finished_compute)
    {
        // Finished computing pending_field, set as active field and use for future color updates
        active_bmp = pending_bmp;
        active_field = pending_field;

        // ======== Result forwarding ========

        switch (computing_phase) {
            case 0:
                /// finished first 9x9 low-res phase, forward computed pixels to 3x3 phase
                field_3x3.setAllDepth(-1.0);
                bmp_9x9.forEachPixel([this](int x, int y)
                {
                    field_3x3(x*3+1, y*3+1) = field_9x9(x, y);

                    // Interior skipping optimization
                    if (!field_9x9.has_data(x, y))
                        field_9x9(x, y).flag_for_skip = true;
                });

                field_9x9.contractSkipFlags(interior_phases_contract_expand.c1);
                field_9x9.expandSkipFlags(interior_phases_contract_expand.e1);
                bmp_9x9.forEachPixel([this](int x, int y)
                {
                    if (field_9x9(x, y).flag_for_skip)
                    {
                        int x0 = x * 3, y0 = y * 3;
                        for (int py = y0; py <= y0 + 3; py++)
                            for (int px = x0; px <= x0 + 3; px++)
                                field_3x3(px, py).flag_for_skip = true;
                    }
                });

                break;

            case 1:
                /// finished second 3x3 low-res phase, forward computed pixels to final 1x1 phase
                field_1x1.setAllDepth(-1.0);
                bmp_3x3.forEachPixel([this](int x, int y)
                {
                    field_1x1(x*3+1, y*3+1) = field_3x3(x, y);

                    // Interior skipping optimization
                    if (!field_3x3.has_data(x, y))
                        field_3x3(x, y).flag_for_skip = true;
                });

                field_3x3.contractSkipFlags(interior_phases_contract_expand.c2);
                field_3x3.expandSkipFlags(interior_phases_contract_expand.e2);
                bmp_3x3.forEachPixel([this](int x, int y)
                {
                    if (field_3x3(x, y).flag_for_skip)
                    {
                        int x0 = x * 3, y0 = y * 3;
                        for (int py = y0; py <= y0 + 3; py++)
                        {
                            for (int px = x0; px <= x0 + 3; px++)
                            {
                                EscapeFieldPixel& pixel = field_1x1(px, py);
                                pixel.depth = INSIDE_MANDELBROT_SET_SKIPPED;
                                pixel.flag_for_skip = true;
                            }
                        }
                    }
                });

                break;

            case 2:
                /// finished final 1x1 phase, permit capture of this frame
                captureFrame(true);
                final_frame_complete = true;

                auto elapsed = std::chrono::steady_clock::now() - compute_t0;
                double dt = std::chrono::duration<double, std::milli>(elapsed).count();
                dt_avg = timer_ma.push(dt);
                break;
        }

        if (computing_phase < 2)
            computing_phase++;
    }

    return finished_compute;
}

bool Mandelbrot_Scene::mandelChanged()
{
    bool view_changed         = Changed(world_quad);
    bool quality_opts_changed = Changed(quality, dynamic_iter_lim, mandel_features, maxdepth_optimize, maxdepth_show_optimized);
    bool compute_opts_changed = Changed(stripe_params);
    bool splines_changed      = Changed(x_spline.hash(), y_spline.hash());
    bool flatten_changed      = Changed(flatten, show_period2_bulb, cardioid_lerp_amount);
    bool features_changed     = Changed(mandel_features);
    return (first_frame || view_changed || quality_opts_changed || compute_opts_changed || splines_changed || flatten_changed || features_changed);
}

bool Mandelbrot_Scene::shadingFormulaChanged()
{
    bool shade_formula_changed   = Changed(shade_formula);
    bool weights_changed         = Changed(iter_weight, dist_weight, stripe_weight);
    bool cycle_iter_opts_changed = Changed(iter_params);
    bool cycle_dist_opts_changed = Changed(dist_params);
    return (shade_formula_changed || weights_changed || cycle_iter_opts_changed || cycle_dist_opts_changed);
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

    // ────── resume unfinished compute ──────
    bool finished_compute = false;
    if (!final_frame_complete)
        finished_compute = processCompute();

    // ────── color cycle changed? ──────
    bool gradient_changed        = updateGradient();
    bool shading_formula_changed = shadingFormulaChanged();

    bool renormalize = finished_compute || shading_formula_changed;
    bool reshade     = renormalize || (gradient_changed && frame_complete);


    // ────── renormalize cached data if necessary  ──────
    if (renormalize)
    {
        FloatingPointType float_type = getRequiredFloatType(mandel_features, camera.relativeZoom<f128>());

        table_invoke( build_table(normalize_shading_limits, [&], active_field, active_bmp, camera, iter_params, dist_params, numThreads()),
            float_type);

        table_invoke( build_table(refreshFieldDepthNormalized, [&], active_field, active_bmp, mandel_features, iter_params, dist_params, numThreads()),
            float_type);

        if (final_frame_complete)
            collectStats();
    }

    // ────── reshade if data renormalized, or gradient changed ──────
    if (reshade)
    {
        table_invoke(
            build_table(shadeBitmap, [&], active_field, active_bmp, &gradient_shifted, (float)iter_weight, (float)dist_weight, (float)stripe_weight, numThreads()),
            (MandelShaderFormula)shade_formula, maxdepth_show_optimized
        );
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

    //if (savefile_changed)
    //    onSavefileChanged();

    // stats for mouse position
    if (active_field && active_bmp)
    {
        if (mouse->stage_x >= 0 && mouse->stage_y >= 0 && mouse->stage_x < active_bmp->width() && mouse->stage_y < active_bmp->height())
        {
            IVec2 pos = active_bmp->pixelPosFromWorld(DDVec2(mouse->world_x, mouse->world_y));
            EscapeFieldPixel* p = active_field->get(pos.x, pos.y);
            if (p) stats.hovered_field_pixel = *p;
        }
    }
}

SIM_END;