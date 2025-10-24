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

void Mandelbrot_Scene::updateTweening(double dt)
{
    if (tweening)
    {
        double ani_mult = fpsFactor();

        if (steady_zoom)
        {
            bool finished_frame = capturedLastFrame();
            if (finished_frame)
            {
                if (camera.relativeZoom<f128>() < state_b.camera.relativeZoom<f128>())
                {
                    ///blPrint() << "FINISHED FRAME. PROGRESSING";

                    auto stepsToReach = [](f128 A, f128 B, f128 factor) {
                        double n = (double)(log(B / A) / log(factor));
                        return (int)ceil(n);
                    };

                    tween_frames_elapsed++;

                    tween_expected_frames = stepsToReach(
                        state_a.camera.relativeZoom<f128>(),
                        state_b.camera.relativeZoom<f128>(),
                        1.0 + steady_zoom_mult_speed); // seconds


                    //static double ease_duration = 1.0; // seconds
                    //double ease_pct_of_total = 

                    // double ease_in_mult = std::min(tween_frames_elapsed / ease_duration, 1.0);
                    // double ease_out_mult = std::min((tween_expected_frames - tween_frames_elapsed) / ease_duration, 1.0);
                    // double ease_mult = std::min(ease_in_mult, ease_out_mult);

                    //camera.zoom *= 1 + (steady_zoom_mult_speed * f128{ ease_mult });

                    ///camera.zoom *= 1.0 + steady_zoom_mult_speed;
                    camera.setRelativeZoom(camera.relativeZoom<f128>() * (1.0 + steady_zoom_mult_speed));
                    camera.setRotation(camera.rotation() + Math::toRadians(0.075));

                    // Update estimated time remaining
                    double expected_time_left = dt * (tween_expected_frames - tween_frames_elapsed);
                    expected_time_left_ma.push(expected_time_left);
                }
                else
                {
                    camera.setRelativeZoom(state_b.camera.relativeZoom<f128>());
                    tweening = false;
                    queueEndRecording();
                }
            }
        }
        else
        {
            double speed = 0.01 / tween_duration;
            tween_progress += speed * ani_mult;

            if (tween_progress < 1.0)
                lerpState(state_a, state_b, tween_progress, false);
            else
            {
                lerpState(state_a, state_b, 1.0, true);

                tween_progress = 0.0;
                tweening = false;
            }
        }
    }
}

void Mandelbrot_Scene::updateGradient()
{
    if (first_frame || Changed(gradient, gradient_shift, hue_shift))
    {
        transformGradient(gradient_shifted, gradient, (float)gradient_shift, (float)hue_shift);
        //savefile_changed = true;
        colors_updated = true;
    }
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
    double rw = round(ctx->width());
    double rh = round(ctx->height());
    int iw = (int)(ceil(rw / 9)) * 9;
    int ih = (int)(ceil(rh / 9)) * 9;
    world_quad = camera.getTransform().toWorldQuad<f128>(0, 0, iw, ih);

    bmp_9x9.setCamera(camera);
    bmp_3x3.setCamera(camera);
    bmp_1x1.setCamera(camera);

    bmp_9x9.setStageRect(0, 0, iw, ih);
    bmp_3x3.setStageRect(0, 0, iw, ih);
    bmp_1x1.setStageRect(0, 0, iw, ih);

    bmp_9x9.setBitmapSize(iw / 9, ih / 9);
    bmp_3x3.setBitmapSize(iw / 3, ih / 3);
    bmp_1x1.setBitmapSize(iw, ih);

    field_9x9.setDimensions(iw / 9, ih / 9);
    field_3x3.setDimensions(iw / 3, ih / 3);
    field_1x1.setDimensions(iw, ih);
}

void Mandelbrot_Scene::updateKernelMode(bool mandel_changed)
{
    constexpr double eps = std::numeric_limits<double>::epsilon();
    int old_smoothing = smoothing_type;
    int new_smoothing = 0;

    if (iter_weight > eps)   new_smoothing |= (int)MandelSmoothing::ITER;
    if (dist_weight > eps)   new_smoothing |= (int)MandelSmoothing::DIST;
    if (stripe_weight > eps) new_smoothing |= (int)MandelSmoothing::STRIPES;

    bool force_upgrade = new_smoothing & ~old_smoothing;
    bool downgrade_on_change = ~new_smoothing & old_smoothing;

    if (force_upgrade || (mandel_changed && downgrade_on_change))
    {
        mandel_changed = true;
        smoothing_type = new_smoothing;
    }
}

void Mandelbrot_Scene::updateActiveField(bool mandel_changed)
{
    // ────── Presented Mandelbrot *actually* changed? Restart on 9x9 bmp (phase 0) ──────
    if (mandel_changed)
    {
        bmp_9x9.clear(0, 255, 0, 255);
        bmp_3x3.clear(0, 255, 0, 255);
        bmp_1x1.clear(0, 255, 0, 255);

        // Hide intro
        if (!first_frame)
            display_intro = false;

        computing_phase = 0;
        current_row = 0;
        final_frame_complete = false;
        field_9x9.setAllDepth(-1.0);

        compute_t0 = std::chrono::steady_clock::now();
    }

    // ────── Prepare bmp / depth-field dimensions and view rect ──────
    {
        // Set pending bitmap & pending field
        switch (computing_phase)
        {
        case 0:  pending_bmp = &bmp_9x9;  pending_field = &field_9x9;  break;
        case 1:  pending_bmp = &bmp_3x3;  pending_field = &field_3x3;  break;
        case 2:  pending_bmp = &bmp_1x1;  pending_field = &field_1x1;  break;
        }
    }
}

bool Mandelbrot_Scene::shouldCompute(bool mandel_changed)
{
    bool phase_changed = Changed(computing_phase);

    if (mandel_changed ||  // Has ANY option would would alter the final mandelbrot changed?
        phase_changed ||   // Just computed last phase, begin computing next phase
        current_row != 0)  // Not finished computing current phase, resume computing current phase
    {
        return true;
    }
    return false;
}

bool Mandelbrot_Scene::processCompute()
{
    bool finished_compute = false;

    iter_lim = finalIterLimit(camera, quality, dynamic_iter_lim, tweening);

    //if (!flatten)
    {
        //DQuad quad = ctx->worldQuad();
        //bool x_axis_visible = quad.intersects({ {quad.minX(), 0}, {quad.maxX(), 0} });

        // Run appropriate kernel for given settings
        finished_compute = frame_complete = compute_mandelbrot(pending_field, pending_bmp);
    }
    //else
    //{
    //    // Flat lerp
    //    finished_compute = dispatchBooleans(
    //        boolsTemplate(radialMandelbrot, [&]),
    //        smoothing_type != MandelSmoothing::NONE,
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
                // ======== Finished first 9x9 low-res phase, forward computed pixels to 3x3 phase ========
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
                // ======== Finished second 3x3 low-res phase, forward computed pixels to final 1x1 phase ========
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
                // ======== Finished final 1x1 phase, permit image/video capture of this frame ========
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

void Mandelbrot_Scene::viewportProcess(Viewport* ctx, double dt)
{
    /// Process Viewports running this Scene
    //bool savefile_changed = false;

    // Never record a frame unless we finish computing a full frame
    captureFrame(false);

    updateAnimation();
    updateTweening(dt);
    updateGradient();
    updateCameraView();
    updateFieldSizes(ctx);

    // ────── any properties changed that require restarting compute? ──────
    bool started_tween = Changed(tweening);
    bool view_changed = Changed(world_quad);
    bool quality_opts_changed = Changed(quality, dynamic_iter_lim, smoothing_type, maxdepth_optimize, maxdepth_show_optimized);
    bool splines_changed = Changed(x_spline.hash(), y_spline.hash());
    bool flatten_changed = Changed(flatten, show_period2_bulb, cardioid_lerp_amount);
    bool compute_opts_changed = Changed(stripe_params);

    bool mandel_changed = (first_frame || view_changed || started_tween || quality_opts_changed || compute_opts_changed || splines_changed || flatten_changed);

    //if (mandel_changed)
    //    savefile_changed = true;

    updateKernelMode(mandel_changed);
    updateActiveField(mandel_changed);

    // Do compute if mandel changed
    // ────── do compute (if mandel changed) ──────
    bool finished_compute = false;
    if (shouldCompute(mandel_changed))
        finished_compute = processCompute();

    // ────── Color cycle changed? ──────
    bool shade_formula_changed = Changed(shade_formula);
    bool weights_changed = Changed(iter_weight, dist_weight, stripe_weight);
    bool cycle_iter_opts_changed = Changed(iter_params);
    bool cycle_dist_opts_changed = Changed(dist_params);

    if (shade_formula_changed || weights_changed || cycle_iter_opts_changed || cycle_dist_opts_changed)
    {
        // Reshade active_bmp
        colors_updated = true;
        //savefile_changed = true;
    }

    bool renormalize = finished_compute || cycle_iter_opts_changed || cycle_dist_opts_changed || shade_formula_changed || weights_changed;
    bool reshade_bmp = renormalize || (colors_updated && frame_complete);

    // ────── Renormalize cached data if necessary  ──────
    if (renormalize)
    {
        normalize_field(active_field, active_bmp);

        if (final_frame_complete)
            collectStats();
    }

    // ────── Refresh normalized values if styles change, or if we do full compute ──────
    if (reshade_bmp)
    {
        table_invoke(
            build_table(shadeBitmap, [&], active_field, active_bmp, &gradient_shifted, (float)iter_weight, (float)dist_weight, (float)stripe_weight, numThreads()),
            (MandelShaderFormula)shade_formula, maxdepth_show_optimized
        );

        colors_updated = false;
        //savefile_changed = true;
    }

    // Unless we're doing a steady-zoom animation, just "record what we see" and capture every frame
    /// todo: Disable for high quality snapshots
    if (!steady_zoom && reshade_bmp && final_frame_complete)
        captureFrame(true);

    first_frame = false;

    if (!steady_zoom && show_interactive_cardioid && animate_cardioid_angle)
    {
        // Cardioid animation
        ani_angle += ani_inc;
        ani_angle = Math::wrapRadians2PI(ani_angle);

        if (isRecording())
            captureFrame(true);
    }

    //if (savefile_changed)
    //    onSavefileChanged();

    // stats for mouse position
    if (active_field && active_bmp)
    {
        //ctx->print() << "\nactive_field->max_depth: " << active_field->max_depth;
        ///ctx->print() << "\nactive_field->min_dist: " << active_field->min_dist;
        ///ctx->print() << "\nactive_field->max_dist: " << active_field->max_dist;

        int px = (int)mouse->stage_x;
        int py = (int)mouse->stage_y;
        if (px >= 0 && py >= 0 && px < active_bmp->width() && py < active_bmp->height())
        {
            IVec2 pos = active_bmp->pixelPosFromWorld(DDVec2(mouse->world_x, mouse->world_y));
            EscapeFieldPixel* p = active_field->get(pos.x, pos.y);

            if (p)
            {
                stats.hovered_field_pixel = *p;
                //double raw_depth = p->depth;
                //double raw_dist = p->getDist

                //double dist = log(raw_dist);

                ///double dist_factor = Math::lerpFactor(dist, active_field->min_dist, active_field->max_dist);

                ///ctx->print() << "\nraw_depth: " << raw_depth << "\n";
                ///ctx->print() << "\nraw_dist: " << raw_dist << "\n";
                /// 
                //ctx->print() << "log_dist: " << dist << "\n";
                ///ctx->print() << "dist_factor: " << dist_factor << "\n\n";

                /*f128 stable_min_raw_dist = camera.getTransform().toWorldOffset<f128>(DDVec2{ 0.5, 0 }).magnitude();
                f128 stable_max_raw_dist = active_bmp->worldSize().magnitude() / 2.0;

                f128 stable_min_dist = log(stable_min_raw_dist);
                f128 stable_max_dist = log(stable_max_raw_dist);

                ctx->print() << "stable_min_raw_dist: " << stable_min_raw_dist << "\n";
                ctx->print() << "stable_max_raw_dist: " << stable_max_raw_dist << "\n\n";

                ctx->print() << "min_possible_dist: " << stable_min_dist << "\n";
                ctx->print() << "max_possible_dist: " << stable_max_dist << "\n\n";

                ctx->print() << "Stabilized factor: " << Math::lerpFactor(dist, stable_min_dist, stable_max_dist);


                double lower_depth_bound = cycle_iter_normalize_depth ? pending_field->min_depth : 0;

                double normalized_depth = Math::lerpFactor(raw_depth, lower_depth_bound, (double)iter_lim);
                double normalized_dist = Math::lerpFactor(raw_dist, pending_field->min_dist, pending_field->max_dist);

                ctx->print() << std::setprecision(18);
                ctx->print() << "\n\nnormalized_depth: " << normalized_depth;
                ctx->print() << "\nnormalized_dist: " << normalized_dist;

                ctx->print() << "\n\nraw_iters: " << raw_depth;
                ctx->print() << "\nraw_dist: " << raw_dist << "\n";

                double single_pixel_raw_dist = camera->toWorldOffset(DVec2{ 0.001, 0 }).magnitude();
                double min_possible_dist = log(single_pixel_raw_dist);

                double max_possible_raw_dist = ctx->worldSize().magnitude();
                double max_possible_dist = log(max_possible_raw_dist);
                //double max_possible_dist = de_cap_from_view(
                //    camera->x(),
                //    camera->y(),
                //    ctx->worldSize().x / 2,
                //    ctx->worldSize().y,
                //    sqrt(escape_radius<MandelSmoothing::DIST>())
                //);

                ctx->print() << "\nlog_dist: " << log(raw_dist) << "\n";

                //ctx->print() << "max_possible_raw_dist: " << max_possible_raw_dist << "\n";
                ctx->print() << "single_pixel_raw_dist: " << single_pixel_raw_dist << "\n";
                ctx->print() << "max_possible_raw_dist: " << max_possible_raw_dist << "\n";

                ctx->print() << "min_possible_dist: " << min_possible_dist << "\n\n";
                ctx->print() << "max_possible_dist: " << max_possible_dist << "\n\n";

                double stable_normalized_dist = Math::lerpFactor(log(raw_dist), pending_field->min_dist, max_possible_dist);
                ctx->print() << "stable_normalized_dist: " << stable_normalized_dist << "\n";
                */
            }
        }
    }
}

SIM_END;