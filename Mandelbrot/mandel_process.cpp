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
            if (camera_vel_pos.mag() > threshold)
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

template<typename T>
void updateNormalizeFieldShape(NormalizationField& norm_field, const CameraInfo &camera, double normalize_field_scale)
{
    norm_field.updateField<T>(camera, normalize_field_scale);
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
        double spacing = (1.0 / normalize_field_quality) / 200.0;
        norm_field.setShape(spacing, normalize_field_exponent);
        changed_normalization_field_shape = true;
    }

    if (first_frame ||
        changed_normalization_field_shape ||
        Changed(world_quad) ||
        Changed(normalize_field_scale))
    {
        FloatingPointType float_type = getRequiredFloatType(mandel_features, camera.zoom<f128>());
        table_invoke(build_table(updateNormalizeFieldShape, [&], norm_field, camera, normalize_field_scale), float_type);
    }

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
        bool use_perturbation = true;

        if (kernel_mode == MandelKernelMode::AUTO)
            use_perturbation = (float_type >= FloatingPointType::F128);
        else if (kernel_mode == MandelKernelMode::FULL)
            use_perturbation = false;

        finished_compute = frame_complete = table_invoke(
            build_table(mandelbrot_perturbation, [&], pending_bmp, pending_field, norm_field, iter_lim, timeout, P, m1, m2, m3, stripe_params),
            float_type, mandel_features, use_perturbation, flatten
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
            {
                /// finished first 9x9 low-res phase, forward computed pixels to 3x3 phase
                //field_3x3.setAllDepth(-1.0);
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
                //field_1x1.setAllDepth(-1.0);
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
    //bool stripe_changed = Changed(stripe_weight);

    bool view_changed            = Changed(world_quad);
    bool quality_opts_changed    = Changed(quality, dynamic_iter_lim, maxdepth_optimize, maxdepth_show_optimized);
    bool compute_opts_changed    = Changed(stripe_params);
    //bool splines_changed         = Changed(x_spline.hash(), y_spline.hash());
    //bool flatten_changed         = Changed(flatten, show_period2_bulb, cardioid_lerp_amount);
    bool features_changed        = Changed(mandel_features);
    bool normalize_field_changed = Changed(normalize_field_quality) || Changed(normalize_field_exponent) || Changed(normalize_field_scale);
    bool kernel_mode_changed     = Changed(kernel_mode);
    bool m_changed = Changed(m1, m2, m3);
    return (
        first_frame ||
        view_changed || 
        quality_opts_changed || 
        compute_opts_changed || 
        //splines_changed || 
        //flatten_changed || 
        features_changed || 
        normalize_field_changed ||
        kernel_mode_changed ||
        m_changed
    );
}

bool Mandelbrot_Scene::shadingFormulaChanged()
{
    bool shade_formula_changed   = Changed(shade_formula);
    bool weights_changed         = Changed(iter_weight, dist_weight, stripe_weight);
    bool cycle_iter_opts_changed = Changed(iter_params);
    bool cycle_dist_opts_changed = Changed(dist_params);
    return (shade_formula_changed || weights_changed || cycle_iter_opts_changed || cycle_dist_opts_changed);
}



SIM_END;