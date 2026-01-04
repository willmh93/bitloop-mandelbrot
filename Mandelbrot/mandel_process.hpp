#include "Mandelbrot.h"

#include "compute.h"

SIM_BEG;

void Mandelbrot_Scene::processBatchSnapshot()
{
    if (rendering_examples)
    {
        if (!isSnapshotting())
        {
            // no active snapshot processing, check if there is one pending...

            if (rendering_example_i >= mandel_presets.size())
            {
                // all examples rendered => end
                rendering_examples = false;
                rendering_example_i = 0;
            }
            else
            {
                // begin new snapshot...

                // grab example name / state data
                auto [key, data] = *std::next(mandel_presets.begin(), rendering_example_i);
                std::string name = key.name;

                // convert example name to lowercase for savefile, replace spaces with underscores
                std::transform(name.begin(), name.end(), name.begin(), [](char ch) {
                    return (ch == ' ') ? '_' : std::tolower(ch);
                });

                std::filesystem::path filepath = render_batch_name;
                filepath /= name;

                loadState(data);

                SnapshotPresetList all_presets = main_window()->getSnapshotPresetManager()->allPresets();
                SnapshotPresetList filtered    = all_presets.filtered(valid_presets);

                beginSnapshotList(filtered, filepath.string().c_str());
            }

            rendering_example_i++;
        }
    }
}

void Mandelbrot_Scene::updateAnimation()
{
    //if (Changed(show_color_animation_options, show_axis, gradient_shift_step, hue_shift_step))
    //    savefile_changed = true;

    double ani_mult = fpsFactor();

    // ────── color cycle animation ──────
    if (!steady_zoom || capturedLastFrame())
    {
        if (final_frame_complete)
        {
            if (show_color_animation_options)
            {
                // Animation
                if (fabs(gradient_shift_step) > 1.0e-4)
                    gradient_shift = math::wrap(gradient_shift + gradient_shift_step * ani_mult, 0.0, 1.0);

                if (fabs(hue_shift_step) > 1.0e-4)
                    hue_shift = math::wrap(hue_shift + hue_shift_step * ani_mult, 0.0, 360.0);
            }
            else
            {
                //if (Changed(gradient_shift, hue_shift, gradient_shift_step, hue_shift_step))
                //    savefile_changed = true;
            }
        }
    }

    // cardioid animation
    if (!steady_zoom && show_interactive_cardioid && animate_cardioid_angle)
    {
        ani_angle += ani_inc;
        ani_angle = math::wrapRadians2PI(ani_angle);
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
            f128 threshold = 0.001 / camera.relativeZoom<f128>();
            if (camera_vel_pos.mag() > threshold)
                camera.setPos(camera.pos<f128>() + camera_vel_pos);
            else
                camera_vel_pos = DDVec2{};

            if (fabs(camera_vel_zoom - 1.0) > 0.01)
                camera.setRelativeZoom(camera.relativeZoom<f128>() * camera_vel_zoom);
            else
                camera_vel_zoom = 1.0;

            camera_vel_pos *= 0.8 / fpsFactor();
            camera_vel_zoom += (1.0 - camera_vel_zoom) * std::min(1.0, 0.2 * fpsFactor()); // Ease back to 1x
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

    bool changed_normalization_field_shape = false;
    if (first_frame || Changed(normalize_field_quality) || Changed(normalize_field_exponent))
    {
        double spacing = (1.0 / normalize_field_quality) / 500.0;
        norm_field.setShape(spacing, normalize_field_exponent);
        changed_normalization_field_shape = true;
    }

    if (first_frame ||
        changed_normalization_field_shape ||
        Changed(world_quad) ||
        Changed(normalize_field_scale))
    {
        FloatingPointType float_type = getRequiredFloatType(mandel_features, camera.relativeZoom<f128>());
        table_invoke(dispatch_table_targ(norm_field, NormalizationField::updateField, camera, normalize_field_scale), float_type);
    }

    //if (Changed(world_quad))
    //{
    //    FloatingPointType float_type = getRequiredFloatType(mandel_features, camera.relativeZoom<f128>());
    //    table_invoke(dispatch_table_targ(norm_field_zoom_out, NormalizationField::updateField, camera, normalize_field_scale), float_type);
    //}
}

void Mandelbrot_Scene::updateEnabledKernelFeatures()
{
    bool mandel_changed = mandelChanged();

    constexpr double eps = std::numeric_limits<double>::epsilon();
    KernelFeatures old_smoothing = mandel_features;
    KernelFeatures new_smoothing = KernelFeatures::NONE;

    if (iter_weight > eps)   new_smoothing |= KernelFeatures::ITER;
    if (dist_weight > eps)   new_smoothing |= KernelFeatures::DIST;
    if (stripe_weight > eps) new_smoothing |= KernelFeatures::STRIPES;

    bool force_upgrade       = (bool)( new_smoothing & ~old_smoothing);
    bool downgrade_on_change = (bool)(~new_smoothing &  old_smoothing);

    if (force_upgrade || (mandel_changed && downgrade_on_change))
        mandel_features = new_smoothing;
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
        P.reset_progress_only();
        final_frame_complete = false;

        field_9x9.reset();
        field_3x3.reset();
        field_1x1.reset();

        norm_field.clearFinalFlags();

        compute_t0 = std::chrono::steady_clock::now();

        timer_begin_group(ITER);
        timer_begin_group(DIST);
        timer_begin_group(STRIPE);
        timer_begin_group(ESCAPE);
        timer_begin_group(REBASE);
        timer_begin_group(NORM_FIELD);
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
    // For high-res phases, break up work across multiple frames if necessary to not lock panning
    int timeout = 0;
    if (!steady_zoom && computing_phase > 0) 
    {
        switch (float_type) {
            case FloatingPointType::F64:  timeout = 128; break;
            case FloatingPointType::F128: timeout = 256; break;
            default: timeout = 0;
        }
    }

    //if (!flatten)
    {
        //DQuad quad = ctx->worldQuad();
        //bool x_axis_visible = quad.intersects({ {quad.minX(), 0}, {quad.maxX(), 0} });

        // Run appropriate kernel for given settings
        KernelMode deduced_kernel_mode = kernel_mode;

        // if "AUTO" => use perturbation at f128 depths
        // otherwise no performance gain => use simple kernel
        if (deduced_kernel_mode == KernelMode::AUTO)
            deduced_kernel_mode = (float_type >= FloatingPointType::F128) ? KernelMode::PERTURBATION_SIMD_UNROLLED : KernelMode::NO_PERTURBATION;

        finished_compute = frame_complete = 
            table_invoke(dispatch_table(compute_mandelbrot, timeout), float_type, mandel_features, deduced_kernel_mode);


        //finished_compute = frame_complete = table_invoke(
        //    build_table(mandelbrot_perturbation, [&], pending_bmp, pending_field, active_field, norm_field, normalize_field_precision, iter_lim, timeout, P, stripe_params),
        //    float_type, mandel_features, deduced_kernel_mode
        //);
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
            {
                /// finished first 9x9 low-res phase, forward computed pixels to 3x3 phase
                int interior_c = 0;
                bmp_9x9.forEachPixel([&](int x, int y)
                {
                    field_3x3(x*3+1, y*3+1) = field_9x9(x, y);
                
                    // Interior skipping optimization
                    if (!field_9x9.has_valid_depth(x, y))
                    {
                        field_9x9.set_skip_flag(x, y, 1);
                        interior_c++;
                    }
                });
                int min_c = interior_phases_contract_expand.c1 * interior_phases_contract_expand.c1;
                if (interior_c > min_c)
                {
                    field_9x9.contractSkipFlags(interior_phases_contract_expand.c1);
                    field_9x9.expandSkipFlags(interior_phases_contract_expand.e1);
                    bmp_9x9.forEachPixel([this](int x, int y)
                    {
                        if (field_9x9.get_skip_flag(x, y))
                        {
                            int x0 = x * 3, y0 = y * 3;
                            for (int py = y0; py <= y0 + 3; py++)
                                for (int px = x0; px <= x0 + 3; px++)
                                    field_3x3.set_skip_flag(px, py, 1);
                        }
                    });
                }
            }
            break;

            case 1:
            {
                /// finished second 3x3 low-res phase, forward computed pixels to final 1x1 phase
                int interior_c = 0;

                bmp_3x3.forEachPixel([&](int x, int y)
                {
                    field_1x1(x*3+1, y*3+1) = field_3x3(x, y);
                
                    // Interior skipping optimization
                    if (!field_3x3.has_valid_depth(x, y))
                    {
                        field_3x3.set_skip_flag(x, y, 1);
                        interior_c++;
                    }
                });
                
                int min_c = interior_phases_contract_expand.c2 * interior_phases_contract_expand.c1;
                if (interior_c > min_c)
                {
                    field_3x3.contractSkipFlags(interior_phases_contract_expand.c2);
                    field_3x3.expandSkipFlags(interior_phases_contract_expand.e2);
                    bmp_3x3.forEachPixel([this](int x, int y)
                    {
                        if (field_3x3.get_skip_flag(x, y))
                        {
                            int x0 = x * 3, y0 = y * 3;
                            for (int py = y0; py <= y0 + 3; py++)
                            {
                                for (int px = x0; px <= x0 + 3; px++)
                                {
                                    EscapeFieldPixel& pixel = field_1x1(px, py);
                                    pixel.depth = INSIDE_MANDELBROT_SET_SKIPPED;
                                    field_1x1.set_skip_flag(px, py, 1);
                                }
                            }
                        }
                    });
                }
            }
            break;

            case 2:
                /// finished final 1x1 phase
                final_frame_complete = true;

                //timer_end_group(ITER);
                //timer_end_group(DIST);
                //timer_end_group(STRIPE);
                //timer_end_group(ESCAPE);
                //timer_end_group(REBASE);
                //timer_end_group(NORM_FIELD);
                //blPrint("\n");

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

void Mandelbrot_Scene::processCapturing(bool finished_compute, bool reshade)
{
    if (final_frame_complete)
    {
        // just finished computing final 1x1 phase? Always permit capture of this frame
        if (finished_compute)
            permitCaptureFrame(true);

        if (!steady_zoom)
        {
            // unless we're doing a deep-zoom, capture basic basic color-cycle animation every frame
            if (reshade)
                permitCaptureFrame(true);

            // If animating the cardioid (and not doing a deep zoom) -> record frame
            if (show_interactive_cardioid && animate_cardioid_angle)
            {
                if (isRecording())
                    permitCaptureFrame(true);
            }
        }
    }
}



bool Mandelbrot_Scene::mandelChanged()
{
    bool view_changed            = Changed(world_quad);
    bool quality_opts_changed    = Changed(quality, dynamic_iter_lim, interior_forwarding, maxdepth_show_optimized);
    bool compute_opts_changed    = Changed(stripe_params.freq);
    bool features_changed        = Changed(mandel_features);
    bool normalize_field_changed = Changed(normalize_field_quality, normalize_field_exponent, normalize_field_scale, normalize_field_precision);
    bool kernel_mode_changed     = Changed(kernel_mode);

    //bool debug_opts_changed = Changed(stripe_mag_numerator, stripe_mag_from_hist);

    return (
        first_frame ||
        view_changed || 
        quality_opts_changed || 
        compute_opts_changed ||
        features_changed || 
        normalize_field_changed ||
        kernel_mode_changed 
        //|| debug_opts_changed
    );
}

bool Mandelbrot_Scene::normalizationOptionsChanged()
{
    bool weights_changed            = Changed(iter_weight, dist_weight, stripe_weight);
    bool shade_formula_changed      = Changed(iter_x_dist_weight,       dist_x_stripe_weight,     stripe_x_iter_weight,
                                              iter_x_distStripe_weight, dist_x_iterStripe_weight, stripe_x_iterDist_weight);

    bool cycle_iter_opts_changed    = Changed(iter_params);
    bool cycle_dist_opts_changed    = Changed(dist_params, dist_tone_params);
    bool cycle_stripe_opts_changed  = Changed(stripe_params.phase, stripe_tone_params);

    return (shade_formula_changed || weights_changed || cycle_iter_opts_changed || cycle_dist_opts_changed || cycle_stripe_opts_changed);
}


SIM_END;