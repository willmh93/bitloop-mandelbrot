﻿#include "Mandelbrot.h"
#ifdef __EMSCRIPTEN__
#include <bitloop/platform/emscripten_browser_clipboard.h>
#endif

SIM_BEG;

std::string dataToURL(std::string data_buf)
{
    #ifdef __EMSCRIPTEN__
    return platform()->url_get_base() + "?data=" + data_buf;
    #else
    return "https://bitloop.dev/Mandelbrot?data=" + data_buf;
    #endif
}

void Mandelbrot_Scene::UI::populateSavingLoading()
{
    if (ImGui::Section("Saving & Loading", true, 0))
    {
        if (ImGui::Button("Save"))
        {
            show_save_dialog = true;
            bl_schedule([&](Mandelbrot_Scene& scene)
            {
                data_buf = scene.getStateData();
            });

            ImGui::OpenPopup("Save Data"); // open on this frame
        }

        ImGui::SameLine();
        if (ImGui::Button("Load"))
        {
            //show_load_dialog = true;

            // Attempt to load immediately from the clipboard (if valid save data)
            //emscripten_browser_clipboard::paste_now(&on_paste, &data_buf);

            // paste immediately on opening dialog
            ///#ifdef __EMSCRIPTEN__
            ///emscripten_browser_clipboard::paste_now([&](std::string&& buf)
            ///{
            ///    data_buf = buf;
            ///    blPrint() << "data_buf: " << data_buf;
            ///
            ///    opening_load_popup = true;
            ///});
            ///#else
            /// 
            data_buf = "";
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
            bl_schedule([&](Mandelbrot_Scene& scene)
            {
                url_buf = dataToURL(scene.getStateData());
            });
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
            ImGui::InputTextMultiline("###Config", &data_buf, avail, ImGuiInputTextFlags_ReadOnly);
            ImGui::PopFont();

            if (ImGui::Button("Copy to Clipboard"))
                ImGui::SetClipboardText(data_buf.c_str());

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
            //blPrint() << "InputTextMultiline @ data_buf: " << data_buf;
            ImGui::InputTextMultiline("###Config", &data_buf, avail, ImGuiInputTextFlags_AlwaysOverwrite);
            ImGui::PopFont();

            if (ImGui::Button("Paste"))
            {
                #ifdef __EMSCRIPTEN__
                emscripten_browser_clipboard::paste_now([&](std::string&& buf) {
                    data_buf = buf;
                });
                #else
                size_t s;
                data_buf = (char*)SDL_GetClipboardData("text/plain", &s);
                #endif
            }

            ImGui::SameLine();
            ImGui::Dummy(ImVec2(10.0f, 0.0f));
            ImGui::SameLine();

            if (ImGui::Button("Load"))
            {
                bl_schedule([&](Mandelbrot_Scene& scene) {
                    scene.loadState(data_buf);
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
            ImGui::InputTextMultiline("###url", &url_buf, avail, ImGuiInputTextFlags_ReadOnly);
            ImGui::PopFont();
            if (ImGui::Button("Copy to Clipboard"))
                ImGui::SetClipboardText(url_buf.c_str());
            ImGui::SameLine();
            if (ImGui::Button("Close")) ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
        }
    }
}
void Mandelbrot_Scene::UI::populateExamples()
{
    if (ImGui::Section("Examples", true))
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
                            f128 current_zoom = scene.camera.relativeZoom<f128>();
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
                            MandelState target;
                            target.deserialize(data);
                            scene.startTween(target);
                        });
                    }
                }
                ImGui::PopID();

                i++;
            }
            ImGui::EndTable();
        }
    }
}
void Mandelbrot_Scene::UI::populateCameraView()
{
    bl_scoped(flatten);
    bl_scoped(camera);

    //static bool show_view_by_default = !platform()->is_mobile(); // Expect navigate by touch for mobile
    if (ImGui::Section("View", true))
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

        camera.populateUI({ -5.0, -5.0, 5.0, 5.0 });

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
}
void Mandelbrot_Scene::UI::populateQualityOptions()
{
    bl_scoped(dynamic_iter_lim);
    bl_scoped(quality);
    bl_scoped(iter_lim);
    bl_scoped(camera);
    bl_scoped(tweening);

    if (ImGui::Section("Quality", false))
    {
        if (ImGui::Checkbox("Dynamic Iteration Limit", &dynamic_iter_lim))
        {
            if (!dynamic_iter_lim)
            {
                quality = iter_lim;
            }
            else
            {
                quality = qualityFromIterLimit(iter_lim, camera.relativeZoom<f128>());
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
    }
}
void Mandelbrot_Scene::UI::populateColorCycleOptions()
{
    if (ImGui::Section("Colour Cycle", true))
    {
        bl_scoped(tweening);
        bl_scoped(iter_lim);
        bl_scoped(camera);
        bl_scoped(quality);
        bl_scoped(dynamic_iter_lim);

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
        sprintf(dist_header, "Distance  -  %d%% Weight", (int)(dist_ratio * 100.0));
        sprintf(stripe_header, "Stripe  -  %d%% Weight", (int)(stripe_ratio * 100.0));

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
            ImGui::SliderFloat("Phase", &stripe_params.phase, 0.0f, Math::PI * 2.0f, "%.5f", ImGuiSliderFlags_AlwaysClamp);
            ImGui::SliderFloat("Contrast", &stripe_params.contrast, 1.0f, 100.0f, "%.3f", ImGuiSliderFlags_AlwaysClamp);
        }

        ImGui::Combo("###MandelFormula", &shade_formula, MandelFormulaNames, (int)MandelShaderFormula::COUNT);
    }
}
void Mandelbrot_Scene::UI::populateGradientShiftOptions()
{
    if (ImGui::Section("Gradient Offset + Animation", true))
    {
        // Shift
        bl_pull(gradient_shifted); // Recalculated in viewportProcess, no need to push
        bl_pull(gradient);

        bl_scoped(hue_shift);
        bl_scoped(gradient_shift);

        // Animation
        bl_scoped(show_color_animation_options);
        bl_scoped(gradient_shift_step);
        bl_scoped(hue_shift_step);

        bl_scoped(colors_updated);

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
}
void Mandelbrot_Scene::UI::populateGradientPicker()
{
    if (ImGui::Section("Base Color Gradient", true, 2.0f))
    {
        bl_scoped(gradient);
        bl_scoped(colors_updated);

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

        ImGui::Spacing();

    }
}
void Mandelbrot_Scene::UI::populateStats()
{
    if (ImGui::Section("Statistics", true))
    {
        bl_scoped(stats);

        // --- compute info ---
        {
            ImGui::GroupBox box("compute_info", "Compute info", scale_size(13.0f), scale_size(20.0f));

            bl_pull(smoothing_type);
            bl_pull(camera);
            bl_pull(dt_avg);
            bl_pull(computing_phase);

            // float precision
            {

                const char* float_precision = FloatingPointTypeNames[(int)getRequiredFloatType((MandelSmoothing)smoothing_type, camera.relativeZoom<f128>())];
                ImGui::Text("Active float precision:  %s", float_precision);
            }

            // current phase
            {
                
                ImGui::Text("Current phase:  %d/2", computing_phase);

                std::stringstream ss;
                for (int i = 0, j = 0; i <= 4; i++) {
                    if (smoothing_type & (1 << i)) {
                        if (j++ != 0) ss << " | ";
                        ss << MandelSmoothingNames[i + 1];
                    }
                }
            }

            // compute timer
            {
                ImGui::Text("Full compute time:  %.0f ms", dt_avg);
            }
            
            ImGui::Spacing();
            ImGui::Spacing();
            
            // mouse pixel info
            {
                ImGui::Text("Mouse info:");

                if (stats.hovered_field_pixel.depth < INSIDE_MANDELBROT_SET_SKIPPED)
                {
                    ImGui::Text("Mouse depth: %.0f iterations", stats.hovered_field_pixel.depth);
                    ImGui::Text("Log depth: %.2f iterations", stats.hovered_field_pixel.final_depth);
                }
                else
                {
                    ImGui::Text("Mouse depth: <escaped> iterations");
                    ImGui::Text("Log depth: <escaped> iterations");
                }
            }
        }

        // --- plots ---
        ImGui::Spacing();
        ImGui::Spacing();

        // depth histogram
        {
            depth_xs.reserve(stats.depth_histogram.size());
            depth_ys.reserve(stats.depth_histogram.size());
            depth_xs.clear();
            depth_ys.clear();

            for (const auto& [k, v] : stats.depth_histogram)
            {
                //if (v <= 1) continue; // Hide very small counts
                depth_xs.push_back(k);
                depth_ys.push_back(v);
            }

            if (ImPlot::BeginPlot("Depth Histogram"))
            {
                ImPlot::SetupAxis(ImAxis_X1, "depth", ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxis(ImAxis_Y1, "pixel count", ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                if (stats.depth_histogram.size())
                {
                    ImPlot::PlotLine("##depth_hist", depth_xs.data(), depth_ys.data(), (int)depth_xs.size());

                    double max_depth = (--stats.depth_histogram.end())->first;
                    if (stats.hovered_field_pixel.depth < max_depth)
                    {

                        auto nearest_it = stats.depth_histogram.lower_bound((int)stats.hovered_field_pixel.depth);
                        if (nearest_it != stats.depth_histogram.end())
                        {
                            double px = nearest_it->first;
                            double py = nearest_it->second;
                            ImPlot::PlotInfLines("##depth_mark", &px, 1); // vertical line at X = mark
                            ImPlot::PlotScatter("##depth_point", &px, &py, 1);
                        }
                    }
                }
                ImPlot::EndPlot();
            }
        }
    }
}

void Mandelbrot_Scene::UI::populateExperimental()
{
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

    if (ImGui::Section("Experimental")) {

        bl_scoped(show_period2_bulb);
        bl_scoped(flatten_amount);
        bl_scoped(flatten);
        bl_scoped(cardioid_lerp_amount);

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
    
}
void Mandelbrot_Scene::UI::populateSplinesDev()
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

void Mandelbrot_Scene::UI::sidebar()
{
    bl_scoped(tweening);
    bl_scoped(colors_updated);
    bl_pull(gradient);

    bool is_tweening = tweening;
    if (is_tweening)
        ImGui::BeginDisabled();
        
    populateSavingLoading();
    populateExamples();
    populateCameraView();
    populateQualityOptions();
    populateColorCycleOptions();
    populateGradientShiftOptions();
    populateGradientPicker();
    populateStats();
    //populateExperimental();
        
    if (is_tweening)
        ImGui::EndDisabled();

    // ======== Developer ========
    #if MANDEL_DEV_EDIT_TWEEN_SPLINES
    populateSplinesDev();
    #endif

    if (colors_updated)
        bl_push(gradient);
}



SIM_END;