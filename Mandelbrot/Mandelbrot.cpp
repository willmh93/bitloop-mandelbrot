#include <bitloop.h>

#include "Mandelbrot.h"


#include "tween.h"
#include "conversions.h"

#ifdef __EMSCRIPTEN__
#include <bitloop/platform/emscripten_browser_clipboard.h>
#endif


SIM_BEG;


/// ─────────────────────── Project ───────────────────────

void Mandelbrot_Project::projectPrepare(Layout& layout)
{
    Mandelbrot_Scene::Config config;
    create<Mandelbrot_Scene>(config)->mountTo(layout);
}

/// ─────────────────────── Scene ───────────────────────

// ────── Compute/Normalize dispatch helpers ──────

bool Mandelbrot_Scene::compute_mandelbrot(EscapeField* field, CanvasImage128* bmp)
{
    MandelSmoothing   smoothing  = static_cast<MandelSmoothing>(smoothing_type);
    FloatingPointType float_type = getRequiredFloatType(smoothing, camera.getRelativeZoom<f128>());

    // Calculate first low-res phase in one-shot (no timeout)
    // For high-res phases, break up work across multiple frames if necessary (kept track of with current_row)
    int timeout = computing_phase == 0 ? 0 : 16;

    return floatInvoke(float_type, [&]<typename T>()
    {
        // compute mandelbrot with first template arg T = [float / double / f128] as determined by float_type...

        return table_invoke<T>(
            /* function arguments    */ build_table(mandelbrot, [&], bmp, field, iter_lim, numThreads(), timeout, current_row, stripe_params),
            /* rest of template args */ smoothing, flatten
        );
    });
}

void Mandelbrot_Scene::normalize_field(EscapeField* field, CanvasImage128* bmp)
{
    MandelSmoothing    smoothing = static_cast<MandelSmoothing>(smoothing_type);
    FloatingPointType float_type = getRequiredFloatType(smoothing, camera.getRelativeZoom<f128>());

    floatInvoke(float_type, [&]<typename T>() { normalize_shading_limits<T>(field, bmp, camera, iter_params, dist_params); });
    floatInvoke(float_type, [&]<typename T>() { refreshFieldDepthNormalized<T>(field, bmp, smoothing, iter_params, dist_params, numThreads()); });
}

// ────── UI ──────

void Mandelbrot_Scene::UI::sidebar()
{
    bl_scoped(camera);
    bl_scoped(tweening);
    bl_scoped(dynamic_iter_lim);
    bl_scoped(quality);
    bl_scoped(iter_lim);
    bl_scoped(colors_updated);

    bl_pull(gradient);

    // experimental
    bl_scoped(flatten);

    {
        bl_scoped(stats);

        
        xs.reserve(stats.depth_histogram.size());
        ys.reserve(stats.depth_histogram.size());
        xs.clear();
        ys.clear();

        for (const auto& [k, v] : stats.depth_histogram) {
            //if (v <= 1) continue; // Hide very small counts
            xs.push_back(k);
            ys.push_back(v);
        }

        if (ImPlot::BeginPlot("Depth Histogram"))
        {
            ImPlot::SetupAxis(ImAxis_X1, "depth", ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxis(ImAxis_Y1, "pixel count", ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            if (stats.depth_histogram.size())
            {
                ImPlot::PlotLine("##depth_hist", xs.data(), ys.data(), (int)xs.size());
            }
            ImPlot::EndPlot();
        }
    }
    
    bool is_tweening = tweening;
    if (is_tweening)
        ImGui::BeginDisabled();

    if (ImGui::Section("Saving & Loading", true, 0))
    {
        bl_scoped(config_buf);

        if (ImGui::Button("Save"))
        {
            show_save_dialog = true;
            //strcpy(config_buf_name, "");
            bl_schedule([](Mandelbrot_Scene& scene)
            {
                scene.updateConfigBuffer();
            });

            ImGui::OpenPopup("Save Data"); // open on this frame
        }

        ImGui::SameLine();
        if (ImGui::Button("Load"))
        {
            //show_load_dialog = true;

            // Attempt to load immediately from the clipboard (if valid save data)
            //emscripten_browser_clipboard::paste_now(&on_paste, &config_buf);

            ///#ifdef __EMSCRIPTEN__
            ///emscripten_browser_clipboard::paste_now([&](std::string&& buf)
            ///{
            ///    config_buf = buf;
            ///    blPrint() << "config_buf: " << config_buf;
            ///
            ///    opening_load_popup = true;
            ///});
            ///#else
            config_buf = "";
            show_load_dialog = true;
            ImGui::OpenPopup("Load Data");
            ///#endif
        }

        #ifdef __EMSCRIPTEN__
        if (opening_load_popup)
        {
            opening_load_popup = false;
            show_load_dialog = true;
            ImGui::OpenPopup("Load Data");
        }
        #endif

        ImGui::SameLine();
        if (ImGui::Button("Share"))
        {
            show_share_dialog = true;
            url = dataToURL(config_buf);
            ImGui::OpenPopup("Share URL");

        }

        // Save Dialog
        ImGui::SetNextWindowSize(scale_size(350, 300), ImGuiCond_FirstUseEver);
        if (ImGui::BeginPopupModal("Save Data", &show_save_dialog))
        {
            ImVec2 avail = ImGui::GetContentRegionAvail();
            avail.y -= ImGui::GetFrameHeightWithSpacing(); // leave room for buttons/input

            //if (!platform()->is_mobile())
            //{
            //    avail.y -= ImGui::GetFrameHeightWithSpacing();
            //    ImGui::AlignTextToFramePadding();
            //    ImGui::Text("Name:");
            //    ImGui::SameLine();
            //    if (ImGui::InputText("###mandel_name", config_buf_name, 28))
            //        updateConfigBuffer();
            //}

            ImGui::PushFont(main_window()->monoFont());
            ImGui::InputTextMultiline("###Config", &config_buf, avail, ImGuiInputTextFlags_ReadOnly);
            ImGui::PopFont();

            if (ImGui::Button("Copy to Clipboard"))
                ImGui::SetClipboardText(config_buf.c_str());

            ImGui::SameLine();
            if (ImGui::Button("Close")) ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
        }

        // Load Dialog
        ImGui::SetNextWindowSize(scale_size(350, 300), ImGuiCond_FirstUseEver);
        if (ImGui::BeginPopupModal("Load Data", &show_load_dialog))
        {
            ImVec2 avail = ImGui::GetContentRegionAvail();
            avail.y -= ImGui::GetFrameHeightWithSpacing(); // leave room for buttons

            ImGui::PushFont(main_window()->monoFont());
            //blPrint() << "InputTextMultiline @ config_buf: " << config_buf;
            ImGui::InputTextMultiline("###Config", &config_buf, avail, ImGuiInputTextFlags_AlwaysOverwrite);
            ImGui::PopFont();

            if (ImGui::Button("Paste"))
            {
                #ifdef __EMSCRIPTEN__
                emscripten_browser_clipboard::paste_now([&](std::string&& buf) {
                    config_buf = buf;
                });
                #else
                size_t s;
                config_buf = (char*)SDL_GetClipboardData("text/plain", &s);
                #endif
            }

            ImGui::SameLine();
            ImGui::Dummy(ImVec2(10.0f, 0.0f));
            ImGui::SameLine();

            if (ImGui::Button("Load"))
            {
                bl_schedule([](Mandelbrot_Scene& scene)
                {
                    scene.loadConfigBuffer();
                });
                ImGui::CloseCurrentPopup();
            }

            ImGui::SameLine();
            if (ImGui::Button("Cancel")) ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
        }

        // Share Dialog
        ImGui::SetNextWindowSize(scale_size(350, 120), ImGuiCond_FirstUseEver);
        if (ImGui::BeginPopupModal("Share URL", &show_share_dialog))
        {
            ImVec2 avail = ImGui::GetContentRegionAvail();
            avail.y -= ImGui::GetFrameHeightWithSpacing() * 1; // leave room for buttons/input
            ImGui::PushFont(main_window()->monoFont());
            ImGui::InputTextMultiline("###url", &url, avail, ImGuiInputTextFlags_ReadOnly);
            ImGui::PopFont();
            if (ImGui::Button("Copy to Clipboard"))
                ImGui::SetClipboardText(url.c_str());
            ImGui::SameLine();
            if (ImGui::Button("Close")) ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
        }
    }

    if (ImGui::Section("Examples", true, 0.0f, 2.0f))
    {
        bl_scoped(steady_zoom);

        if (platform()->is_desktop_native())
            ImGui::Checkbox("Steady zoom", &steady_zoom); // For desktop deep zoom recordings


        ImGuiStyle& style = ImGui::GetStyle();
        float avail_full = ImGui::GetContentRegionAvail().x;
        float min_btn_w = scale_size(100.0f);
        int   cols = (int)((avail_full + style.ItemSpacing.x) / (min_btn_w + style.ItemSpacing.x));
        cols = cols < 1 ? 1 : cols;
    
        if (ImGui::BeginTable("preset_grid", cols, ImGuiTableFlags_SizingStretchProp))
        {
            int i = 0;
            for (const auto& [key, data] : mandel_presets)
            {
                const auto& [category, name] = key;

                ImGui::TableNextColumn();
                ImGui::PushID(i);
                if (ImGui::Button(category))
                {
                    //MandelState dest;
                    //dest.deserialize(preset.data);

                    
                    if (steady_zoom)
                    {
                        bl_schedule([data](Mandelbrot_Scene& scene)
                        {
                            scene.tween_frames_elapsed = 0;

                            // Give destination same reference zoom level
                            scene.state_b.camera.setReferenceZoom(scene.camera.getReferenceZoom<f128>());

                            // Set destination
                            scene.state_b.deserialize(data);

                            // Set current state to match, but reset back to current zoom
                            f128 current_zoom = scene.camera.getRelativeZoom<f128>();
                            static_cast<MandelState&>(scene) = scene.state_b;
                            scene.camera.setRelativeZoom(current_zoom);

                            // Mark starting state for checking lerp progress
                            scene.state_a = static_cast<MandelState&>(scene);

                            // Begin tweening
                            scene.tweening = true;
                        });
                    }
                    else
                    {
                        bl_schedule([data](Mandelbrot_Scene& scene)
                        {
                            // Give destination same reference zoom level
                            scene.state_b.camera.setReferenceZoom(scene.camera.getReferenceZoom<f128>());

                            scene.state_b.deserialize(data);
                            startTween(scene);
                        });
                    }
                }
                ImGui::PopID();

                i++;
            }
            ImGui::EndTable();
        }
    }

    /// --------------------------------------------------------------
    //static bool show_view_by_default = !platform()->is_mobile(); // Expect navigate by touch for mobile
    if (ImGui::Section("View", true, 5.0f, 2.0f))
    {
        bl_scoped(show_axis);

        if (ImGui::Button("Fullscreen"))
        {
            SDL_SetWindowFullscreen(platform()->sdl_window(), true);
        }
        ImGui::SameLine();
        ImGui::Checkbox("Show Axis", &show_axis);

        #if MANDEL_FEATURE_INTERACTIVE_CARDIOID
        bl_scoped(show_interactive_cardioid);
        bl_scoped(ani_angle);
        bl_scoped(ani_inc);
        bl_scoped(animate_cardioid_angle);

        if (!flatten && !platform()->is_mobile())
        {
            ImGui::SameLine();
            ImGui::Dummy(ImVec2(10.0f, 0.0f));
            ImGui::SameLine();
            ImGui::Checkbox("Show Interactive Cardioid", &show_interactive_cardioid);
        }
        #endif
        
        camera.populateUI({-5.0, -5.0, 5.0, 5.0});

        #if MANDEL_FEATURE_INTERACTIVE_CARDIOID
        if (show_interactive_cardioid)
        {
            ImGui::GroupBox box("cardioid_group", "Interactive Cardioid", scale_size(13.0f), scale_size(20.0f));
            ImGui::Checkbox("Apply Animation", &animate_cardioid_angle);
            ImGui::SliderAngle("Angle", &ani_angle);
            if (animate_cardioid_angle)
                ImGui::SliderAngle("Angle Inc", &ani_inc, -Math::HALF_PI, Math::HALF_PI, 1);
        }
        #endif
    }

    /// --------------------------------------------------------------
    if (ImGui::Section("Quality", false, 5.0f, 2.0f))
    {
        if (ImGui::Checkbox("Dynamic Iteration Limit", &dynamic_iter_lim))
        {
            if (!dynamic_iter_lim)
            {
                quality = iter_lim;
            }
            else
            {
                quality = qualityFromIterLimit(iter_lim, camera.getRelativeZoom<f128>());
            }
        }

        if (dynamic_iter_lim)
        {
            double quality_pct = quality * 100.0;
            static double initial_quality = quality_pct;
            ImGui::PushID("QualitySlider");
            ImGui::RevertableDragDouble("###Quality", &quality_pct, &initial_quality, 1, 1.0, 500.0, "%.0f%%");
            ImGui::PopID();
            quality = quality_pct / 100.0;

            ImGui::SameLine();
            ImGui::Text("= %d Iters", finalIterLimit(camera, quality, dynamic_iter_lim, tweening));
        }
        else
            ImGui::DragDouble("Max Iterations", &quality, 1000.0, 1.0, 1000000.0, "%.0f", ImGuiSliderFlags_Logarithmic);

        bl_scoped(maxdepth_optimize);
        bl_scoped(interior_phases_contract_expand);
        bl_scoped(maxdepth_show_optimized);

        ImGui::Spacing(); ImGui::Spacing();
        ImGui::Text("Max-depth result forwarding");
        if (ImGui::Combo("###MandelMaxDepthOptimization", &maxdepth_optimize, MandelMaxDepthOptimizationNames, (int)MandelMaxDepthOptimization::COUNT))
        {
            switch ((MandelMaxDepthOptimization)maxdepth_optimize)
            {
            case MandelMaxDepthOptimization::SLOWEST: interior_phases_contract_expand = { 10, 0, 10, 0 }; break;
            case MandelMaxDepthOptimization::SLOW:    interior_phases_contract_expand = { 6, 2, 8, 2 };   break;
            case MandelMaxDepthOptimization::MEDIUM:  interior_phases_contract_expand = { 5, 3, 7, 3 };   break;
            case MandelMaxDepthOptimization::FAST:    interior_phases_contract_expand = { 1, 0, 1, 0 };   break;
            default: break;
            }
        }

        ImGui::Checkbox("Highlight optimized regions", &maxdepth_show_optimized);
    } // End Header

    //if (ImGui::Section("Instant Styles", false, 5.0f, 2.0f))
    //{
    //    if (ImGui::Button("High Contrast"))
    //    {
    //        iter_dist_mix = 1;
    //        //cycle_iter_dynamic_limit = true;
    //        //cycle_iter_normalize_depth = true;
    //        //cycle_iter_value = 4;
    //        cycle_dist_invert = false;
    //        cycle_dist_value = 1;
    //        cycle_dist_sharpness = 100;
    //        show_color_animation_options = false;
    //        gradient_shift = 0.0;
    //        hue_shift = 0.0;
    //        loadGradientPreset(GradientPreset::CLASSIC);
    //    }
    //}

    /// --------------------------------------------------------------
    if (ImGui::Section("Colour Cycle", true, 5.0f, 2.0f))
    {
        // Weights / Mix formula
        bl_scoped(iter_weight);
        bl_scoped(dist_weight);
        bl_scoped(stripe_weight);
        bl_scoped(shade_formula);

        // Iter
        bl_scoped(use_smoothing); /// todo

        // Stripe
        bl_scoped(iter_params);
        bl_scoped(dist_params);
        bl_scoped(stripe_params);


        double iter_ratio, dist_ratio, stripe_ratio;
        shadingRatios(
            iter_weight, dist_weight, stripe_weight,
            iter_ratio, dist_ratio, stripe_ratio
        );

        char iter_header[64], dist_header[64], stripe_header[64];
        sprintf(iter_header, "Iteration  -  %d%% Weight", (int)(iter_ratio * 100.0));
        sprintf(dist_header, "Distance  -  %d%% Weight",  (int)(dist_ratio * 100.0));
        sprintf(stripe_header, "Stripe  -  %d%% Weight",  (int)(stripe_ratio * 100.0));

        // ────── Color Cycle: ITER ──────
        {
            ImGui::GroupBox box("iter_group", iter_header, scale_size(13.0f), scale_size(20.0f));
            ImGui::SliderDouble("Iter Ratio", &iter_weight, 0.0, 1.0, "%.2f");
            ImGui::Spacing(); ImGui::Spacing();

            if (ImGui::Checkbox("% of Max Iters", &iter_params.cycle_iter_dynamic_limit))
            {
                if (iter_params.cycle_iter_dynamic_limit)
                    iter_params.cycle_iter_value /= iter_lim;
                else
                    iter_params.cycle_iter_value *= iter_lim;
            }
            ImGui::SameLine(); ImGui::Dummy(ImVec2(10.0f, 0.0f)); ImGui::SameLine();

            ImGui::Checkbox("Smooth", &use_smoothing);
            ImGui::Checkbox("Normalize to Zoom", &iter_params.cycle_iter_normalize_depth);

            if (iter_params.cycle_iter_normalize_depth)
            {
                if (iter_params.cycle_iter_normalize_depth)
                {
                    ImGui::SliderDouble("Low %", &iter_params.cycle_iter_normalize_low_fact, 0.001, 100.0, "%.4f%%", ImGuiSliderFlags_AlwaysClamp);
                    ImGui::SliderDouble("High %", &iter_params.cycle_iter_normalize_high_fact, 0.001, 100.0, "%.4f%%", ImGuiSliderFlags_AlwaysClamp);
                }
                else
                {

                }
            }


            float required_space = 0.0f;
            box.IncreaseRequiredSpaceForLabel(required_space, iter_params.cycle_iter_dynamic_limit ? "% Iterations" : "Iterations");
            box.IncreaseRequiredSpaceForLabel(required_space, "Logarithmic");

            double raw_cycle_iters;
            if (iter_params.cycle_iter_dynamic_limit)
            {
                double cycle_pct = iter_params.cycle_iter_value * 100.0;

                box.SetNextItemWidthForSpace(required_space);
                ImGui::SliderDouble("% Iterations", &cycle_pct, 0.001, 100.0, "%.4f%%",
                    (/*color_cycle_use_log1p ?*/ ImGuiSliderFlags_Logarithmic /*: 0*/) |
                    ImGuiSliderFlags_AlwaysClamp);

                iter_params.cycle_iter_value = cycle_pct / 100.0;

                raw_cycle_iters = finalIterLimit(camera, quality, dynamic_iter_lim, tweening) * iter_params.cycle_iter_value;

            }
            else
            {
                box.SetNextItemWidthForSpace(required_space);
                ImGui::SliderDouble("Iterations", &iter_params.cycle_iter_value, 1.0, (float)iter_lim, "%.3f",
                    (/*color_cycle_use_log1p ?*/ ImGuiSliderFlags_Logarithmic /*: 0*/) |
                    ImGuiSliderFlags_AlwaysClamp);

                raw_cycle_iters = iter_params.cycle_iter_value;
            }

            double log_pct = iter_params.cycle_iter_log1p_weight * 100.0;
            box.SetNextItemWidthForSpace(required_space);
            if (ImGui::SliderDouble_InvLog("Logarithmic", &log_pct, 0.0, 100.0, "%.2f%%", ImGuiSliderFlags_NoRoundToFormat | ImGuiSliderFlags_AlwaysClamp))
            {
                iter_params.cycle_iter_log1p_weight = log_pct / 100.0;
            }


            ///// Print cycle iter values
            ///{
            ///    if (cycle_iter_dynamic_limit)
            ///        ImGui::Text("raw_cycle_iters = %.1f", raw_cycle_iters);
            ///
            ///    double log_cycle_iters = Math::linear_log1p_lerp(raw_cycle_iters, cycle_iter_log1p_weight);
            ///    ImGui::Text("log_cycle_iters = %.1f", log_cycle_iters);
            ///}
        }

        // ────── Color Cycle: DIST ──────
        {
            ImGui::GroupBox box("dist_group", dist_header);
            ImGui::SliderDouble("Dist Weight", &dist_weight, 0.0, 1.0, "%.2f");
            ImGui::Spacing(); ImGui::Spacing();

            float required_width = 0.0f;
            box.IncreaseRequiredSpaceForLabel(required_width, "Dist");
            box.IncreaseRequiredSpaceForLabel(required_width, "Sharpness");

            ImGui::Checkbox("Invert", &dist_params.cycle_dist_invert);

            box.SetNextItemWidthForSpace(required_width);
            ImGui::SliderDouble("Dist", &dist_params.cycle_dist_value, 0.001, 1.0, "%.5f", ImGuiSliderFlags_Logarithmic | ImGuiSliderFlags_AlwaysClamp);
            ImGui::PushID("DistLog");
            ImGui::PopID();

            box.SetNextItemWidthForSpace(required_width);
            ImGui::SliderDouble_InvLog("Sharpness", &dist_params.cycle_dist_sharpness, 0.0, 100.0, "%.4f%%", ImGuiSliderFlags_NoRoundToFormat | ImGuiSliderFlags_AlwaysClamp);
        }

        // ────── Color Cycle: STRIPE ──────
        {
            ImGui::GroupBox box("stripe_group", stripe_header);
            ImGui::SliderDouble("Stripe Weight", &stripe_weight, 0.0, 1.0, "%.2f");
            ImGui::Spacing();
            ImGui::Spacing();

            ImGui::SliderFloat("Frequency", &stripe_params.freq, 1.0f, 100.0f, "%.0f", ImGuiSliderFlags_AlwaysClamp);
            ImGui::SliderFloat("Phase",     &stripe_params.phase, 0.0f, Math::PI*2.0f, "%.5f", ImGuiSliderFlags_AlwaysClamp);
            ImGui::SliderFloat("Contrast",  &stripe_params.contrast, 1.0f, 100.0f, "%.3f", ImGuiSliderFlags_AlwaysClamp);
        }

        ImGui::Combo("###MandelFormula", &shade_formula, MandelFormulaNames, (int)MandelShaderFormula::COUNT);
    }

    if (ImGui::Section("Color Offset + Animation", true, 5.0f, 2.0f))
    {
        // Shift
        bl_pull(gradient_shifted); // Recalculated in viewportProcess, no need to push
        bl_scoped(hue_shift);
        bl_scoped(gradient_shift);

        // Animation
        bl_scoped(show_color_animation_options);
        bl_scoped(gradient_shift_step);
        bl_scoped(hue_shift_step);

        ImGui::Checkbox("Animate", &show_color_animation_options);

        float required_space = 0.0f;
        ImGui::IncreaseRequiredSpaceForLabel(required_space, "Gradient shift");

        ImGui::SetNextItemWidthForSpace(required_space);
        if (ImGui::DragDouble("Gradient shift", &gradient_shift, 0.01, -100.0, 100.0, " %.3f", ImGuiSliderFlags_AlwaysClamp))
        {
            gradient_shift = Math::wrap(gradient_shift, 0.0, 1.0);
            colors_updated = true;

            transformGradient(gradient_shifted, gradient, (float)gradient_shift, (float)hue_shift);
        }

        if (show_color_animation_options)
        {
            ImGui::Indent();
            ImGui::PushID("gradient_increment");
            ImGui::SetNextItemWidthForSpace(required_space);
            ImGui::SliderDouble("Increment", &gradient_shift_step, -0.02, 0.02, "%.4f");
            ImGui::PopID();
            ImGui::Unindent();
        }

        ImGui::SetNextItemWidthForSpace(required_space);
        if (ImGui::SliderDouble("Hue shift", &hue_shift, 0.0, 360, "%.3f"))
        {
            colors_updated = true;
            transformGradient(gradient_shifted, gradient, (float)gradient_shift, (float)hue_shift);
        }

        if (show_color_animation_options)
        {
            ImGui::Indent();
            ImGui::PushID("hue_increment");
            ImGui::SetNextItemWidthForSpace(required_space);
            ImGui::SliderDouble("Increment", &hue_shift_step, -5.0, 5.0, "%.3f", ImGuiSliderFlags_AlwaysClamp);
            ImGui::PopID();
            ImGui::Unindent();
        }


        ImGui::Spacing();
        ImGui::Text("Live preview");
        ImGui::GradientButton(&gradient_shifted, platform()->dpr());
    }

    if (ImGui::Section("Base Color Gradient", true, 5.0f, 2.0f))
    {
        ImGui::Text("Load Preset");
        static int selecting_template = -1;
        if (ImGui::Combo("###ColorTemplate", &selecting_template, ColorGradientNames, (int)GradientPreset::COUNT))
        {
            generateGradientFromPreset(gradient, (GradientPreset)selecting_template);

            selecting_template = -1;
            colors_updated = true;
        }

        if (ImGui::Button("Copy gradient C++ code"))
        {
            ImGui::SetClipboardText(gradient.to_cpp_marks().c_str());
        }

        ImGui::Dummy(scale_size(0, 8));

        if (ImGui::GradientEditor(&gradient,
            platform()->ui_scale_factor(), 
            platform()->ui_scale_factor(2.0f)))
        {
            colors_updated = true;

            // Shift
            bl_pull(gradient_shifted);
            bl_pull(hue_shift);
            bl_pull(gradient_shift);

            transformGradient(gradient_shifted, gradient, (float)gradient_shift, (float)hue_shift);
        }

    } // End Header

    /*if (ImGui::Section("Experimental")) {

        bl_scoped(show_period2_bulb);
        bl_scoped(flatten_amount);

        ImGui::Checkbox("Flatten", &flatten);

        if (flatten)
        {
            ImGui::Indent();
            if (ImGui::SliderDouble("Flatness", &flatten_amount, 0.0, 1.0, "%.2f"))
                cardioid_lerp_amount = 1.0 - flatten_amount;

            ImGui::Checkbox("Show period-2 bulb", &show_period2_bulb);
            ImGui::Unindent();
            ImGui::Dummy(scale_size(0, 10));
        }

        //static ImRect vr = { 0.0f, 0.8f, 0.8f, 0.0f };
        //ImGui::SeparatorText("Iteration Spline Mapping");
        //if (ImSpline::SplineEditor("iter_spline", &iter_gradient_spline, &vr))
        //{
        //    colors_updated = true;
        //}

        if (!flatten)
        {
            /// --------------------------------------------------------------
            ImGui::SeparatorText("XX, YY Spline Relationship");
            /// --------------------------------------------------------------

            bl_scoped(x_spline);
            bl_scoped(y_spline);

            static ImRect vr = { 0.0f, 0.8f, 0.8f, 0.0f };
            ImSpline::SplineEditorPair("X/Y Spline", &x_spline, &y_spline, &vr, 900.0f);
        }

    } // End Header
    */
    
    if (is_tweening)
        ImGui::EndDisabled();

    // ======== Developer ========
    #if MANDEL_DEV_EDIT_TWEEN_SPLINES
    {
        bl_scoped(tween_pos_spline);
        bl_scoped(tween_zoom_lift_spline);
        bl_scoped(tween_base_zoom_spline);

        static ImRect vr = { 0.0f, 1.0f, 1.0f, 0.0f };

        ImGui::SeparatorText("Position Tween");
        ImSpline::SplineEditor("tween_pos", &tween_pos_spline, &vr);
        ImSpline::SplineEditor("tween_zoom_lift", &tween_zoom_lift_spline, &vr);
        ImSpline::SplineEditor("tween_base_zoom", &tween_base_zoom_spline, &vr);

        //ImGui::InputTextMultiline("###pos_buf", pos_tween_buf, 1024, ImVec2(0, 0), ImGuiInputTextFlags_AllowTabInput));
        if (ImGui::Button("Copy position spline")) ImGui::SetClipboardText(tween_pos_spline.serialize(SplineSerializationMode::CPP_ARRAY, 3).c_str());

        //ImGui::InputTextMultiline("###pos_buf", zoom_tween_buf, 1024, ImVec2(0, 0), ImGuiInputTextFlags_AllowTabInput));
        if (ImGui::Button("Copy lift spline")) ImGui::SetClipboardText(tween_zoom_lift_spline.serialize(SplineSerializationMode::CPP_ARRAY, 3).c_str());

        if (ImGui::Button("Copy base zoom spline")) ImGui::SetClipboardText(tween_base_zoom_spline.serialize(SplineSerializationMode::CPP_ARRAY, 3).c_str());
    }
    #endif


    if (colors_updated)
        bl_push(gradient);
}

void Mandelbrot_Scene::updateConfigBuffer()
{
    config_buf = serialize();
}

void Mandelbrot_Scene::loadConfigBuffer()
{
    deserialize(config_buf);
}

void Mandelbrot_Scene::onSavefileChanged()
{
    config_buf = serialize();
    #ifdef __EMSCRIPTEN__
    //platform()->url_set_string("data", config_buf.c_str());
    #endif
}


// ────── Scene overrides ──────

void Mandelbrot_Scene::sceneStart()
{
    // todo: Stop this getting called twice on startup
    //loadGradientPreset(GradientPreset::CLASSIC);
    generateGradientFromPreset(gradient, GradientPreset::CLASSIC);

    cardioid_lerper.create(Math::TWO_PI / 5760.0, 0.005);

    font = NanoFont::create("/data/fonts/DroidSans.ttf");
}

void Mandelbrot_Scene::sceneMounted(Viewport* ctx)
{
    //camera->setCameraStageSnappingSize(1);
    camera.setSurface(ctx);
    camera.setOriginViewportAnchor(Anchor::CENTER);
    camera.focusWorldRect(-2, -1.25, 1, 1.25);
    camera.uiSetCurrentAsDefault();

    navigator.setTarget(camera);
    navigator.setDirectCameraPanning(true);
    //camera->restrictRelativeZoomRange(0.5, 1e+300);


    reference_zoom = camera.getReferenceZoom<f128>();
    ctx_stage_size = ctx->size();

    #ifdef __EMSCRIPTEN__
    if (platform()->url_has("data"))
    {
        config_buf = platform()->url_get_string("data");
        loadConfigBuffer();
    }
    #endif

    /*if (!load_example_name.empty())
    {
        config_buf = mandel_presets[load_example_name];
        if (!config_buf.empty())
            loadConfigBuffer();
    }*/
}

void Mandelbrot_Scene::viewportProcess(Viewport* ctx, double dt)
{
    //blPrint() << "------------------- frame ------------------------";

    //double ani_dt = dt;
    //if (steady_zoom) ani_dt = 1.0 / 60.0;
    double ani_mult = fpsFactor();// ani_dt / 0.016;

    /// Process Viewports running this Scene
    ctx_stage_size = ctx->size();

    
    
    bool savefile_changed = false;

    if (Changed(show_color_animation_options, show_axis, gradient_shift_step, hue_shift_step))
        savefile_changed = true;

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
            if (Changed(gradient_shift, hue_shift, gradient_shift_step, hue_shift_step))
                savefile_changed = true;
        }
    }

    //if (input.touch.now().fingers.size())
    //    blPrint() << input.touch.now().fingers[0].stage_x;

    // ────── Tweening ──────
    if (tweening)
    {
        if (steady_zoom)
        {
            bool finished_frame = capturedLastFrame();
            if (finished_frame)
            {
                if (camera.getRelativeZoom<f128>() < state_b.camera.getRelativeZoom<f128>())
                {
                    ///blPrint() << "FINISHED FRAME. PROGRESSING";

                    auto stepsToReach = [](f128 A, f128 B, f128 factor) {
                        double n = (double)(log(B / A) / log(factor));
                        return (int)ceil(n);
                    };

                    tween_frames_elapsed++;

                    tween_expected_frames = stepsToReach(
                        state_a.relative_zoom(),
                        state_b.relative_zoom(),
                        1.0 + steady_zoom_mult_speed); // seconds


                    //static double ease_duration = 1.0; // seconds
                    //double ease_pct_of_total = 

                    // double ease_in_mult = std::min(tween_frames_elapsed / ease_duration, 1.0);
                    // double ease_out_mult = std::min((tween_expected_frames - tween_frames_elapsed) / ease_duration, 1.0);
                    // double ease_mult = std::min(ease_in_mult, ease_out_mult);

                    //camera.zoom *= 1 + (steady_zoom_mult_speed * f128{ ease_mult });

                    ///camera.zoom *= 1.0 + steady_zoom_mult_speed;
                    camera.setRelativeZoom(camera.getRelativeZoom<f128>() * (1.0 + steady_zoom_mult_speed));
                    camera.setRotation(camera.rotation() + Math::toRadians(0.075));

                    // Update estimated time remaining
                    double expected_time_left = dt * (tween_expected_frames - tween_frames_elapsed);
                    expected_time_left_ma.push(expected_time_left);
                }
                else
                {
                    camera.setRelativeZoom(state_b.camera.getRelativeZoom<f128>());
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
                lerpState(*this, state_a, state_b, tween_progress, false);
            else
            {
                lerpState(*this, state_a, state_b, 1.0, true);

                tween_progress = 0.0;
                tweening = false;
            }
        }
    }

    // ────── Updating Shifted Gradient ──────
    if (first_frame || Changed(gradient, gradient_shift, hue_shift))
    {
        transformGradient(gradient_shifted, gradient, (float)gradient_shift, (float)hue_shift);
        savefile_changed = true;
        colors_updated = true;
    }

    // ────── Calculate Depth Limit ──────
    iter_lim = finalIterLimit(camera, quality, dynamic_iter_lim, tweening);

    // ────── Flattening ──────
    #if MANDEL_FEATURE_FLATTEN_MODE
    if (Changed(flatten))
    {
        if (flatten) flatten_amount = 0.0;
        camera->focusWorldRect(-2, -1, 2, 1);
        camera.zoom = camera->getRelativeZoom<f128>().x;
    }

    if (Changed(flatten_amount))
    {
        using namespace Math;
        double t = flatten_amount;
        DRect r;

        if (t < 0.5)      r = lerp(DRect(-2.0,-1.5,0.5,1.5), DRect(-2.0,-0.2,1.5,3.5), lerpFactor(t,0.0,0.5));
        else if (t < 0.7) r = lerp(DRect(-2.0,-0.2,1.5,3.5), DRect(-1.5,-0.2,4.0,3.5), lerpFactor(t,0.5,0.7));
        else              r = lerp(DRect(-1.5,-0.2,4.0,3.5), DRect( 0.0,-1.5,4.5,0.5), lerpFactor(t,0.7,1.0));

        camera->focusWorldRect(r, false);
        camera.zoom = camera->getRelativeZoom<f128>().x;
    }
    #endif

    // ────── Camera View ──────
    if (!tweening) {

        #if MANDEL_FEATURE_CAMERA_EASING
        // Process camera velocity
        if (mouse->pressed)
        {
            camera_vel_pos = DDVec2{};

            // Stop zoom velocity on touch
            //if (platform()->is_mobile())
            camera_vel_zoom = 1.0;
        }
        else
        {
            f128 threshold = camera.getRelativeZoom<f128>() / 1000.0;
            if (camera_vel_pos.magnitude() > threshold)
                camera.setPos(camera.pos<f128>() + camera_vel_pos);
            else
                camera_vel_pos = DDVec2{};

            if (fabs(camera_vel_zoom - 1.0) > 0.01)
                camera.setRelativeZoom(camera.getRelativeZoom<f128>() * camera_vel_zoom);
            else
                camera_vel_zoom = 1.0;

            camera_vel_pos *= 0.8;
            camera_vel_zoom += (1.0 - camera_vel_zoom) * 0.2; // Ease back to 1x
        }
        #endif
    }

    // ────── Ensure size divisble by 9 for result forwarding from: [9x9] to [3x3] to [1x1] ──────
    double upscale = 1.0;
    double rw = ceil(ctx->width() / upscale);
    double rh = ceil(ctx->height() / upscale);
    int iw = (int)(ceil(rw / 9)) * 9;
    int ih = (int)(ceil(rh / 9)) * 9;

    ///blPrint() << std::setprecision(30);
    ///DDQuad old_quad = world_quad;
    world_quad = camera.getTransform().toWorldQuad<f128>(0, 0, iw * upscale, ih * upscale);
    ///
    ///blPrint() << "\n\nOld world_quad:" << old_quad <<
    ///             "\nNew world_quad:" << world_quad <<
    ///             (world_quad != old_quad ? " (DIFFERENT)" : "") << "\n\n";


    // ────── Any properties changed that require restarting compute? ──────
    bool started_tween        = Changed(tweening);
    bool view_changed         = Changed(world_quad);
    bool quality_opts_changed = Changed(quality, dynamic_iter_lim, smoothing_type, maxdepth_optimize, maxdepth_show_optimized);
    bool splines_changed      = Changed(x_spline.hash(), y_spline.hash());
    bool flatten_changed      = Changed(flatten, show_period2_bulb, cardioid_lerp_amount);
    bool compute_opts_changed = Changed(stripe_params);

    bool mandel_changed = (first_frame || view_changed || started_tween || quality_opts_changed || compute_opts_changed || splines_changed || flatten_changed);

    if (mandel_changed)
        savefile_changed = true;

    // ────── Upgrade mode if new feature, downgrade if lost feature & mandel changed  ──────
    {
        constexpr double eps = std::numeric_limits<double>::epsilon();
        int old_smoothing = smoothing_type;
        int new_smoothing = 0;

        if (iter_weight > eps)   new_smoothing |= (int)MandelSmoothing::ITER;
        if (dist_weight > eps)   new_smoothing |= (int)MandelSmoothing::DIST;
        if (stripe_weight > eps) new_smoothing |= (int)MandelSmoothing::STRIPES;

        bool force_upgrade       = new_smoothing & ~old_smoothing;
        bool downgrade_on_change = ~new_smoothing & old_smoothing;

        if (force_upgrade || (mandel_changed && downgrade_on_change))
        {
            mandel_changed = true;
            smoothing_type = new_smoothing;
        }
    }

    // ────── Presented Mandelbrot *actually* changed? Restart on 9x9 bmp (phase 0) ──────
    if (mandel_changed)
    {
        bmp_9x9.clear(0, 255, 0, 255);
        bmp_3x3.clear(0, 255, 0, 255);
        bmp_1x1.clear(0, 255, 0, 255);

        // Hide intro
        if (!first_frame)
            display_intro = false;

        blPrint() << "RESETTING PHASE to 0";

        computing_phase = 0;
        current_row = 0;
        final_frame_complete = false;
        field_9x9.setAllDepth(-1.0);

        #if MANDEL_DEV_PERFORMANCE_TIMERS
        compute_t0 = std::chrono::steady_clock::now();
        #endif
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

        pending_bmp->setCamera(camera);
        pending_bmp->setStageRect(0, 0, iw * upscale, ih * upscale);

        bmp_9x9.setBitmapSize(iw / 9, ih / 9);
        bmp_3x3.setBitmapSize(iw / 3, ih / 3);
        bmp_1x1.setBitmapSize(iw, ih);

        field_9x9.setDimensions(iw / 9, ih / 9);
        field_3x3.setDimensions(iw / 3, ih / 3);
        field_1x1.setDimensions(iw, ih);
    }


    bool do_compute = false;
    bool finished_compute = false;

    bool phase_changed = Changed(computing_phase);

    // ────── Start new (or resume an ongoing) compute? ──────
    if (mandel_changed   ||  // Has ANY option would would alter the final mandelbrot changed?
        phase_changed    ||  // Finished computing last phase, begin computing next phase
        current_row != 0)    // Not finished computing current phase, resume computing current phase
    {
        do_compute = true;
    }

    // Never record a frame unless we finish computing a full frame
    captureFrame(false);

    // ────── Compute Mandelbrot & Normalize ──────
    if (do_compute)
    {
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

                    #if MANDEL_DEV_PERFORMANCE_TIMERS
                    auto elapsed = std::chrono::steady_clock::now() - compute_t0;
                    double dt = std::chrono::duration<double, std::milli>(elapsed).count();
                    dt_avg = timer_ma.push(dt);
                    #endif
                    break;
            }

            if (computing_phase < 2)
                computing_phase++;
        }
    }

    // ────── Color cycle changed? ──────
    bool shade_formula_changed    = Changed(shade_formula);
    bool weights_changed          = Changed(iter_weight, dist_weight, stripe_weight);
    bool cycle_iter_opts_changed  = Changed(iter_params);
    bool cycle_dist_opts_changed  = Changed(dist_params);

    if (shade_formula_changed || weights_changed || cycle_iter_opts_changed || cycle_dist_opts_changed)
    {
        // Reshade active_bmp
        colors_updated = true;
        savefile_changed = true;
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
        savefile_changed = true;
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

    if (savefile_changed)
        onSavefileChanged();
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

    double zoom_mag = camera.getRelativeZoom<double>();

    if (show_axis && zoom_mag < 1.0e7)
        ctx->drawWorldAxis(0.5, 0, 0.5);

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

    ///ctx->print() << "Shift: " << gradient_shift;
    ///ctx->print() << "\nFPS Factor: " << fpsFactor();

    #ifndef BL_DEBUG
    if (display_intro)
    {
        ctx->setFont(font);  
        ctx->setFontSize(20);

        //if (!platform()->is_desktop_native())
        {
            ctx->fillText("Controls:", scale_size(10.0), scale_size(10.0));
            ctx->fillText("  - Touch & drag to move", scale_size(10.0), scale_size(35.0));
            ctx->fillText("  - Pinch to zoom / rotate", scale_size(10.0), scale_size(60.0));
        }

        ctx->setTextAlign(TextAlign::ALIGN_LEFT);
        ctx->setTextBaseline(TextBaseline::BASELINE_BOTTOM);
        //ctx->fillText("Developer:    Will Hemsworth", scale_size(10), ctx->height() - scale_size(32));
        ctx->fillText("Contact:  will.hemsworth@bitloop.dev", scale_size(10.0), ctx->height() - scale_size(10.0));
    }
    #endif

    #if MANDEL_DEV_PRINT_ACTIVE_FLOAT_TYPE
    {
        ctx->print() << "Quality:   " << FloatingPointTypeNames[(int)getRequiredFloatType((MandelSmoothing)smoothing_type, camera.zoom)] << "\n";
    }
    #endif


    #if MANDEL_DEV_PRINT_ACTIVE_COMPUTE_FEATURES
    {
        ctx->print() << "Computing: ";
        for (int i = 0, j = 0; i <= 4; i++) {
            if (smoothing_type & (1 << i)) {
                if (j++ != 0) ctx->print() << " | ";
                ctx->print() << MandelSmoothingNames[i + 1];
            }
        }
        ctx->print() << "\n";
    }
    #endif

    #if MANDEL_DEV_PERFORMANCE_TIMERS
    {
        ctx->print() << "Compute timer: " << dt_avg << "\n";
    }
    #endif

    /*ctx->print() << "\ntween_frames_elapsed: " << tween_frames_elapsed << "\n";
    ctx->print() << "tween_expected_frames: " << tween_expected_frames << "\n";
    ctx->print() << "time remaining: " << expected_time_left_ma.average() << "\n";*/

    //for (auto pair : ui_stage)
    //{
    //    auto entry = pair.second;
    //    ctx->print() << entry.name << " live value: " << entry.to_string_value() << "\n";
    //    //ctx->print() << entry.name << " live mrked: " << entry.to_string_marked_live() << "\n";
        //ctx->print() << entry.name << " shdw mrked: " << entry.to_string_marked_shadow() << "\n";
    //}

    /*if (active_field && active_bmp)
    {
        //ctx->print() << "\nactive_field->max_depth: " << active_field->max_depth;
        ///ctx->print() << "\nactive_field->min_dist: " << active_field->min_dist;
        ///ctx->print() << "\nactive_field->max_dist: " << active_field->max_dist;

        int px = (int)mouse->stage_x;
        int py = (int)mouse->stage_y;
        if (px >= 0 && py >= 0 && px < active_bmp->width() && py < active_bmp->height())
        {
            IVec2 pos = active_bmp->pixelPosFromWorld(DVec2(mouse->world_x, mouse->world_y));
            EscapeFieldPixel* p = active_field->get(pos.x, pos.y);

            if (p)
            {
                double raw_depth = p->depth;
                double raw_dist = p->dist;

                double dist = log(raw_dist);

                ///double dist_factor = Math::lerpFactor(dist, active_field->min_dist, active_field->max_dist);

                ctx->print() << "\nraw_depth: " << raw_depth << "\n";
                ctx->print() << "\nraw_dist: " << raw_dist << "\n";
                //ctx->print() << "log_dist: " << dist << "\n";
                ///ctx->print() << "dist_factor: " << dist_factor << "\n\n";

                double stable_min_raw_dist = camera->toWorldOffset(DVec2{ 0.5, 0 }).magnitude();
                double stable_max_raw_dist = active_bmp->worldSize().magnitude() / 2.0;

                double stable_min_dist = log(stable_min_raw_dist);
                double stable_max_dist = log(stable_max_raw_dist);

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

            }
        }
    }
    */

    //DQuad quad = ctx->worldQuad();
    //bool x_axis_visible = quad.intersects({ {quad.minX(), 0}, {quad.maxX(), 0}});
    //ctx->print() << "\nx_axis_visible: " << (x_axis_visible ? "true" : "false");

    //ctx->printTouchInfo();



    //ctx->print() << "camera_vel: " << camera_vel_pos << "\n";
   
    //if (input.touch.finger(0).pressed)
    //{
    //  minimizeMaximizeButton(ctx, false);
    //}

    //if (Input::Touch().Fingers(0).Released(minimizeMaximizeButton(ctx, false)))
    //{
    //
    //}
    //
    //if (input.touch.now().fingers(0).clicked(minimizeMaximizeButton(ctx, false)))
    //{
    //}
}

void Mandelbrot_Scene::collectStats()
{
    if (!active_field) return;
    EscapeField& field = *active_field;

    memset(&stats.dirty, 0, sizeof(stats.dirty));

    // Depth Histogram
    {
        int bucket_size = std::max(1, (int)((field.max_depth - field.min_depth) / 1000.0));
        if (bucket_size > 0)
        {
            auto& hist = stats.depth_histogram;
            hist.clear();
            for (size_t i = 0; i < field.size(); i++)
            {
                EscapeFieldPixel& p = field[i];
                if (p.depth > 1e10 || isnan(p.depth)) continue;

                int bucket_depth = (int)(p.depth / bucket_size) * bucket_size;
                hist[bucket_depth]++;
            }

            // Trim small isolated entries from back of entries
            if (hist.size() >= 2)
            {
                int min_depth = hist.begin()->first;
                int max_depth = (--hist.end())->first;

                i64 sum_count = 0;
                for (auto [k, c] : hist) sum_count += c;
                int mean_px_count = (int)((double)sum_count / (double)hist.size());
                int min_px_count = mean_px_count / 50;
                auto valid_entry = [&](int d) { return hist.count(d) && hist[d] >= min_px_count; };

                // Find 'x' consecutive valid entries from the back and stop
                const int consecutive_entries = 5;
                int trim_from_depth;
                for (trim_from_depth = max_depth; trim_from_depth >= min_depth; trim_from_depth -= bucket_size)
                {
                    int d = (int)(trim_from_depth / bucket_size) * bucket_size;
                    bool gaps = false;
                    for (int j = 0; j < consecutive_entries; j++)
                    {
                        if (!valid_entry(d - bucket_size * j)) {
                            gaps = true;
                            break;
                        }
                    }
                    if (!gaps) break;
                }

                auto it = hist.upper_bound(trim_from_depth);
                hist.erase(it, hist.end());
                max_depth = trim_from_depth;

                for (int d = min_depth + bucket_size; d < max_depth; d += bucket_size)
                {
                    int bucket_depth = (int)(d / bucket_size) * bucket_size;
                    if (hist.count(bucket_depth) == 0)
                        hist[bucket_depth] = 0;
                }
            }
            stats.dirty.depth_histogram = true;
        }
    }
}

void Mandelbrot_Scene::onEvent(Event e)
{
    #ifdef BL_RELEASE
    if (tweening || isRecording())
        return;
    #endif

    // enter/exit fullscreen hotkeys (F11/Escape)
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
    f128 old_zoom = camera.getRelativeZoom<f128>();
    
    if (navigator.handleWorldNavigation(e, true, true))
    {
        f128 zoom = camera.getRelativeZoom<f128>();
        avg_vel_pos.push(camera.pos<f128>() - old_pos);
        avg_vel_zoom.push( (double)(zoom / old_zoom) );

        camera_vel_pos = avg_vel_pos.average();
        camera_vel_zoom = (avg_vel_zoom.average() - 1) * 0.6 + 1;

        if (!platform()->is_mobile())
        {
            // Stopped wheeling mouse
            avg_vel_zoom.clear();
        }
    }

    if (e.type() == SDL_EVENT_MOUSE_BUTTON_UP)
    {
        avg_vel_zoom.clear();
        avg_vel_pos.clear();
    }
}

SIM_END;
