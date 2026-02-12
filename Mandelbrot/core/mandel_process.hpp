#include "../Mandelbrot.h"

#include "compute.h"

SIM_BEG;

void Mandelbrot_Scene::loadNextBatchSnapshotState()
{
    // if "Render All" was clicked, call beginSnapshot() for each preset until done
    assert(rendering_examples);

    if (!isSnapshotting())
    {
        // no active snapshot processing, check if there is one pending...
        const auto& examples = bookmark_manager.find("Examples").getItems();

        if (rendering_example_i >= examples.size())
        {
            // all examples rendered => end
            rendering_examples = false;
            rendering_example_i = 0;
        }
        else
        {
            // begin new snapshot...
            
            // grab example name / state data
            const auto& example = examples[rendering_example_i];

            const std::string& data = example.data;
            const std::string name = example.thumbName();
           
            std::filesystem::path filepath = render_batch_name;
            filepath /= name + "_%s";
            
            loadState(data);
            
            SnapshotPresetList all_presets = main_window()->getSettingsConfig()->enabledImagePresets();
            SnapshotPresetList filtered    = ignore_preset_filters ?
                all_presets : 
                all_presets.filtered(valid_presets);
            
            beginSnapshotList(filtered, filepath.string().c_str());
        }
            
        rendering_example_i++;
    }
}

void Mandelbrot_Scene::updateAnimation()
{
    double ani_mult = fpsFactor();

    // ────── color cycle animation ──────
    if (!steady_zoom || capturedLastFrame())
    {
        if (final_frame_complete)
        {
            bool stepped = false;

            // Animation
            if (animate_gradient_shift)
            {
                if (fabs(gradient_shift_step) > 1.0e-4)
                {
                    gradient_shift = math::wrap(gradient_shift + gradient_shift_step * ani_mult, 0.0, 1.0);
                    stepped = true;
                }
            }

            if (animate_gradient_hue)
            {
                if (fabs(hue_shift_step) > 1.0e-4)
                {
                    hue_shift = math::wrap(hue_shift + hue_shift_step * ani_mult, 0.0, 360.0);
                    stepped = true;
                }
            }

            if (animate_stripe_phase)
            {
                if (fabs(phase_step) > 1.0e-4)
                {
                    stripe_params.phase = math::wrap(stripe_params.phase + phase_step * (f32)ani_mult, 0.0f, math::tau_f);
                    stepped = true;
                }
            }

            if (stepped)
                requestRedraw(true);
        }
    }

    // cardioid animation
    if (!steady_zoom && show_interactive_cardioid && animate_cardioid_angle)
    {
        ani_angle += ani_inc * ani_mult;
        ani_angle = math::wrapRadians2PI(ani_angle);
    }

    if (Changed(ani_angle))
        requestRedraw(true);
}

bool Mandelbrot_Scene::updateGradient()
{
    if (first_frame || Changed(gradient, gradient_shift, hue_shift))
    {
        transformGradient(gradient_shifted, gradient, (float)gradient_shift, (float)hue_shift);
        gradient_tex.set(gradient_shifted);

        return true;
    }
    return false;
}

void Mandelbrot_Scene::updateCameraView(Viewport* ctx)
{
    #if MANDEL_FEATURE_CAMERA_EASING
    if (!tweening)
    {
        // Process camera velocity
        if (mouse->buttonDown(MouseButton::LEFT))
        {
            cam_vel_pos = DDVec2{};

            // Stop zoom velocity on touch
            cam_vel_zoom = 1.0;
        }
        else
        {
            f128 threshold = 0.001 / camera.relativeZoom<f128>();
            if (cam_vel_pos.mag() > threshold)
                camera.setPos(camera.pos<f128>() + cam_vel_pos);
            else
                cam_vel_pos = DDVec2{};

            if (fabs(cam_vel_zoom - 1.0) > 0.01)
                camera.setRelativeZoom(camera.relativeZoom<f128>() * cam_vel_zoom);
            else
                cam_vel_zoom = 1.0;

            const f64 alpha = 0.85;
            const f64 f = fpsFactor();

            cam_vel_pos  *= std::pow(alpha, f);
            cam_vel_zoom = 1.0 - (1.0 - cam_vel_zoom) * std::pow(alpha, f); // Ease back to 1x
        }
    }
    #endif

    // update surface size (snapped to lowest resolution) + world quad
    int max_division = (int)phaseDownscaleFactor(0);
    surface_w = (int)(ceil(ctx->width() / max_division)) * max_division;
    surface_h = (int)(ceil(ctx->height() / max_division)) * max_division;
    world_quad = camera.getTransform().toWorldQuad<f128>(DRect{ surface_w, surface_h });
}

void Mandelbrot_Scene::updateFieldSizes()
{
    // set field/raster grid sizes + target render rect
    int cur_division = (int)phaseDownscaleFactor(0);
    for (int i = 0; i < PHASE_COUNT; i++)
    {
        int phase_w = surface_w / cur_division;
        int phase_h = surface_h / cur_division;

        phases[i].field.setDimensions(phase_w, phase_h);
        phases[i].grid.setRasterSize(phase_w, phase_h);
        phases[i].grid.setStageRect(0, 0, surface_w, surface_h);

        cur_division /= 3;
    }

    // update normalization field shape
    bool changed_normalization_field_shape = false;
    if (first_frame || Changed(normalize_field_quality) || Changed(normalize_field_exponent))
    {
        double spacing = (1.0 / normalize_field_quality) / 500.0;
        norm_field.setShape(spacing, normalize_field_exponent);
        changed_normalization_field_shape = true;
    }

    // update normalization field sample points
    if (first_frame ||
        Changed(world_quad) ||
        changed_normalization_field_shape ||
        Changed(normalize_field_scale/*, normalize_field_precision*/)
        )
    {
        table_invoke(dispatch_table_targ(norm_field, NormalizationField::updateField, camera, normalize_field_scale), float_type);
    }
}

void Mandelbrot_Scene::updateEnabledKernelFeatures()
{
    bool mandel_changed = mandelChanged();

    constexpr double eps = std::numeric_limits<double>::epsilon();
    KernelFeatures old_features = mandel_features;
    KernelFeatures new_features = KernelFeatures::NONE;

    if (iter_weight > eps)   new_features |= KernelFeatures::ITER;
    if (dist_weight > eps)   new_features |= KernelFeatures::DIST;
    if (stripe_weight > eps) new_features |= KernelFeatures::STRIPES;

    // upgrade immediately and recompute when a shader requires a new features (even if view doesn't change),
    // but don't recompute when downgrading available features until it becomes necessary (e.g. a view change)
    bool force_upgrade       = (bool)( new_features & ~old_features);
    bool downgrade_on_change = (bool)(~new_features &  old_features);

    if (force_upgrade || (mandel_changed && downgrade_on_change))
        mandel_features = new_features;
}

void Mandelbrot_Scene::updateQuality()
{
    bool view_changed = Changed(camera);
    bool quality_opts_changed = Changed(quality, dynamic_iter_lim, tweening);

    if (first_frame || view_changed || quality_opts_changed)
    {
        float_type = getRequiredFloatType(mandel_features, camera.relativeZoom<f128>());
        iter_lim = finalIterLimit(camera, quality, dynamic_iter_lim, tweening);
    }
}

void Mandelbrot_Scene::resetCompute()
{
    // hide intro on first change
    if (!first_frame) 
        display_intro = false;

    // reset phases
    starting_phase = next_starting_phase;
    computing_phase = starting_phase;

    // reset forEachPixel progress tracker
    P.reset_progress_only();
    final_frame_complete = false;

    //// reset fields
    norm_field.clearFinalFlags();
    for (MandelPhaseData& phase : phases)
        phase.field.reset();

    // reset timers
    compute_timer.begin();
    for (int phase = 0; phase < PHASE_COUNT; phase++)
        phase_timers[phase] = 0;
}

bool Mandelbrot_Scene::processCompute()
{
    bool finished_compute = false; // todo: rename to phase_complete?

    // Calculate first low-res phase in one-shot (no timeout)
    // For high-res phases, break up work across multiple frames if necessary to not lock panning
    int timeout = updateComputeTimeout();

    // Run appropriate kernel for given settings
    deduced_kernel_mode = kernel_mode;

    // if "AUTO" => use perturbation at f128 depths
    // otherwise no performance gain => use simple kernel
    if (deduced_kernel_mode == KernelMode::AUTO)
    {
        // todo: Make sure we don't switch mode inbetween phases? Some kernels might rely on cached info
        //       from the same kernel 
        if (float_type >= FloatingPointType::F128)
            deduced_kernel_mode = KernelMode::PERTURBATION_SIMD_UNROLLED;
        else
            deduced_kernel_mode = KernelMode::NO_PERTURBATION;
    }

    // resume computing this phase
    {
        double elapsed = 0.0;

        finished_compute = frame_complete = 
            table_invoke(dispatch_table(compute_mandelbrot, timeout, elapsed), float_type, mandel_features, deduced_kernel_mode);

        phase_timers[computing_phase] += elapsed;
    }

    if (finished_compute)
    {
        if (computing_phase > starting_phase)
        {
            auto& src_sma = phase_elapsed_mult_sma_list[computing_phase];
            f64 mult = phase_timers[computing_phase] / phase_timers[computing_phase - 1];

            phase_elapsed_mult_results[computing_phase] = src_sma.push(mult);
            phase_elapsed_mult_sma_result = phase_elapsed_mult_sma.push(mult);
        }

        if (computing_phase == starting_phase)
        {
            phase_elapsted_estimated_final = predictFinalPhaseDuration(computing_phase, phase_timers[computing_phase], phase_elapsed_mult_results);
            if (starting_phase == 0)
            {
                if (phase_elapsted_estimated_final < 1500)
                    next_starting_phase = 1;
            }
            else
            {
                if (phase_elapsted_estimated_final > 1750)
                    next_starting_phase = 0;
            }
        }

        // Finished computing pending_field, set as active field and use for future color updates
        active_phase = computing_phase;

        // ======== Result forwarding ========
        if (interior_forwarding != (int)MandelInteriorForwarding::SLOWEST)
        {
            int next_phase = computing_phase + 1;
            WorldRasterGrid128& src_grid = phases[computing_phase].grid;
            EscapeField& src_field = phases[computing_phase].field;
            EscapeField& dst_field = phases[next_phase].field;

            if (computing_phase < PHASE_COUNT - 1)
            {
                ContractExpandPhase& contact_expand_phase = contract_expand_phases[computing_phase];
                int interior_c = 0;

                src_grid.forEachPixel([&](int x, int y)
                {
                    dst_field(x * 3 + 1, y * 3 + 1) = src_field(x, y);

                    // Interior skipping optimization
                    if (!src_field.has_valid_depth(x, y)) {
                        src_field.set_skip_flag(x, y, 1);
                        interior_c++;
                    }
                });

                // Only apply optimization if there's a significant amount of interior that can be forwarded
                int min_c = contact_expand_phase.contract * contact_expand_phase.contract;
                if (interior_c > min_c)
                {
                    src_field.contractSkipFlags(contact_expand_phase.contract);
                    src_field.expandSkipFlags(contact_expand_phase.expand);

                    bool final_forward = (computing_phase == PHASE_COUNT - 2);
                    table_invoke(lambda_table([&]<bool Final_Step>()
                    {
                        src_grid.forEachPixel([&](int x, int y)
                        {
                            if (!src_field.get_skip_flag(x, y)) return;

                            int x0 = x * 3;
                            int y0 = y * 3;
                            for (int py = y0; py <= y0 + 3; py++) {
                                for (int px = x0; px <= x0 + 3; px++) {
                                    EscapeFieldPixel& pixel = dst_field(px, py);
                                    if constexpr (Final_Step)
                                        pixel.depth = INSIDE_MANDELBROT_SET_SKIPPED;
                                    dst_field.set_skip_flag(px, py, 1);
                                }
                            }
                        });
                    }), final_forward);
                };
            }
        }

        if (computing_phase >= PHASE_COUNT - 1)
        {
            // finished final phase
            final_frame_complete = true;
            dt_last = compute_timer.elapsed();
        }
        else
        {
            // more phases to go
            computing_phase++;
            requestImmediateUpdate();
        }
    }
    else
    {
        if (!P.finished())
            requestImmediateUpdate();
    }

    return finished_compute;
}

void Mandelbrot_Scene::processCapturing(bool finished_compute, bool reshade)
{
    if (final_frame_complete)
    {
        // just finished computing final phase? Always permit capture of this frame
        if (finished_compute)
            permitCaptureFrame(true);

        // did we *begin* capture on final phase without needing to recompute? 
        // Previous finished_compute condition won't trigger
        if (capturedFrameCount() == 0)
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
    bool quality_opts_changed    = Changed(iter_lim, dynamic_iter_lim, interior_forwarding, maxdepth_show_optimized);
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
    bool cycle_iter_opts_changed    = Changed(iter_params);
    bool cycle_dist_opts_changed    = Changed(dist_params, dist_tone_params);
    bool cycle_stripe_opts_changed  = Changed(stripe_params.phase, stripe_tone_params);
    bool norm_field_changed         = Changed(normalize_field_scale, normalize_field_quality, normalize_field_exponent, normalize_field_precision);

    return (weights_changed || cycle_iter_opts_changed || cycle_dist_opts_changed || cycle_stripe_opts_changed || norm_field_changed);
}


SIM_END;