#include "Mandelbrot.h"
#include "compute.h"

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

// Call once per frame where you want the button to appear.
// Returns true if it toggled.
// Example:
//   static bool is_fullscreen = false;
//   DrawFullscreenOverlayButton(&is_fullscreen, ImVec2(24, 24), 30.0f, 0.66f);

void Mandelbrot_Scene::UI::init()
{
    ImGuiStyle& style = ImGui::GetStyle();
    style.HoverStationaryDelay = 0.50f;
}

void Mandelbrot_Scene::UI::sidebar()
{
    bl_scoped(tweening);

    if (tweening)
        ImGui::BeginDisabled();

    populateSavingLoading();
    populateCameraView();
    populateExamples();
    populateQualityOptions();
    populateColorCycleOptions();
    populateGradientShiftOptions();
    populateGradientPicker();
    populateStats();
    //populateExperimental();
    if (!platform()->is_mobile())
        populateMouseOrbit();

    if (tweening)
        ImGui::EndDisabled();

    // ======== Developer ========
    #if MANDEL_DEV_EDIT_TWEEN_SPLINES
    populateSplinesDev();
    #endif
}


bool DrawFullscreenOverlayButton(bool* fullscreen, bool *was_held, ImVec2 screen_pos, float size = 30.0f, float alpha = 0.66f)
{
    ImDrawList* dl = ImGui::GetForegroundDrawList();
    ImGuiIO& io = ImGui::GetIO();

    const ImRect r(screen_pos, screen_pos + ImVec2(size, size));

    // Colors
    const ImU32 col_bg = ImGui::GetColorU32(ImVec4(0, 0, 0, alpha));
    const ImU32 col_bg_h = ImGui::GetColorU32(ImVec4(0, 0, 0, alpha + 0.15f));
    const ImU32 col_bg_p = ImGui::GetColorU32(ImVec4(0, 0, 0, alpha + 0.25f));
    const ImU32 col_bd = ImGui::GetColorU32(ImVec4(1, 1, 1, 0.25f));
    const ImU32 col_fg = ImGui::GetColorU32(ImGui::GetStyle().Colors[ImGuiCol_Text]);

    // Hit test
    const bool hovered = r.Contains(io.MousePos);
    const bool held = hovered && ImGui::IsMouseDown(ImGuiMouseButton_Left);
    const bool clicked = hovered && *was_held && ImGui::IsMouseReleased(ImGuiMouseButton_Left);
    *was_held = held;

    // Background
    const float rounding = size / 6.0f;
    dl->AddRectFilled(r.Min, r.Max, held ? col_bg_p : (hovered ? col_bg_h : col_bg), rounding);
    dl->AddRect(r.Min, r.Max, col_bd, rounding);

    // Icon geometry
    const float t = scale_size(2.0f);     // stroke thickness
    const float inset = size * 0.08f;         // distance from edges
    const float len = size * 0.15f;         // target leg length

    const float max_len = ImMax(0.0f, (size - 2.0f * inset) * 0.5f);
    const float L = ImMin(len, max_len);

    auto lineL = [&](const ImVec2& c, const ImVec2& dx, const ImVec2& dy)
    {
        dl->AddLine(c, c + dx * L, col_fg, t);
        dl->AddLine(c, c + dy * L, col_fg, t);
    };

    // Outer inset corners
    const ImVec2 tl = r.Min + ImVec2(inset, inset);
    const ImVec2 tr = ImVec2(r.Max.x - inset, r.Min.y + inset);
    const ImVec2 bl = ImVec2(r.Min.x + inset, r.Max.y - inset);
    const ImVec2 br = r.Max - ImVec2(inset, inset);

    if (!*fullscreen)
    {
        // Expand: inner box corners (open outward)
        const ImVec2 ia = r.Min + ImVec2(inset + L, inset + L);
        const ImVec2 ib = r.Max - ImVec2(inset + L, inset + L);

        lineL(ImVec2(ia.x, ia.y), ImVec2(+1, 0), ImVec2(0, +1)); // top-left
        lineL(ImVec2(ib.x, ia.y), ImVec2(-1, 0), ImVec2(0, +1)); // top-right
        lineL(ImVec2(ia.x, ib.y), ImVec2(+1, 0), ImVec2(0, -1)); // bottom-left
        lineL(ImVec2(ib.x, ib.y), ImVec2(-1, 0), ImVec2(0, -1)); // bottom-right
    }
    else
    {
        // Collapse: from outer corners inward
        const float inset2 = inset * 2.0f;
        const ImVec2 ia = r.Min + ImVec2(inset2 + L, inset2 + L);
        const ImVec2 ib = r.Max - ImVec2(inset2 + L, inset2 + L);

        lineL(ImVec2(ia.x, ia.y), ImVec2(-1, 0), ImVec2(0, -1)); // top-left inward
        lineL(ImVec2(ib.x, ia.y), ImVec2(+1, 0), ImVec2(0, -1)); // top-right inward
        lineL(ImVec2(ia.x, ib.y), ImVec2(-1, 0), ImVec2(0, +1)); // bottom-left inward
        lineL(ImVec2(ib.x, ib.y), ImVec2(+1, 0), ImVec2(0, +1)); // bottom-right inward
    }

    if (clicked)
    {
        *fullscreen = !*fullscreen;
        return true;
    }
    return false;
}

void Mandelbrot_Scene::UI::overlay()
{
    IVec2 ctx_size = main_window()->viewportSize();
    float btn_size = 60.0f;
    float btn_space = 10.0f;

    SDL_WindowFlags flags = SDL_GetWindowFlags(platform()->sdl_window());
    bool fullscreen = (flags & SDL_WINDOW_FULLSCREEN);
    static bool fullscreen_btn_held = false;
    if (DrawFullscreenOverlayButton(&fullscreen, &fullscreen_btn_held, scale_size(ctx_size.x - btn_size - btn_space, btn_space), scale_size(btn_size)))
    {
        SDL_SetWindowFullscreen(platform()->sdl_window(), fullscreen);
        main_window()->setSidebarVisible(!fullscreen);
    }
}

void Mandelbrot_Scene::UI::populateSavingLoading()
{
    if (ImGui::CollapsingHeaderBox("Saving & Loading", true))
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

        ImGui::EndCollapsingHeaderBox();
    }
}
void Mandelbrot_Scene::UI::populateCameraView()
{
    bl_scoped(flatten);
    bl_scoped(camera);

    //static bool show_view_by_default = !platform()->is_mobile(); // Expect navigate by touch for mobile
    if (ImGui::CollapsingHeaderBox("View", true))
    {
        bl_scoped(show_axis);

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
            //ImGui::GroupBox box("cardioid_group", "Interactive Cardioid", scale_size(13.0f), scale_size(20.0f));
            ImGui::BeginLabelledBox("Interactive Cardioid");

            ImGui::Checkbox("Apply Animation", &animate_cardioid_angle);
            ImGui::SliderAngle("Angle", &ani_angle);
            if (animate_cardioid_angle)
                ImGui::SliderAngle("Angle Step", &ani_inc, -Math::HALF_PI, Math::HALF_PI, 1);

            ImGui::EndLabelledBox();
        }
        #endif

        ImGui::EndCollapsingHeaderBox();
    }
}
void Mandelbrot_Scene::UI::populateExamples()
{
    if (ImGui::CollapsingHeaderBox("Examples", true))
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
                            main_window()->setFixedFrameTimeDelta(true);

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
                            scene.loadState(data);

                            // Give destination same reference zoom level
                            ///MandelState target;
                            ///target.camera.setReferenceZoom(scene.camera.getReferenceZoom<f128>());
                            ///target.deserialize(data);
                            ///scene.startTween(target);
                        });
                    }
                }
                ImGui::PopID();

                i++;
            }
            ImGui::EndTable();
        }

        ImGui::EndCollapsingHeaderBox();
    }
}
void Mandelbrot_Scene::UI::populateQualityOptions()
{
    if (ImGui::CollapsingHeaderBox("Quality", false))
    {
        bl_view(camera, iter_lim, tweening);
        bl_scoped(dynamic_iter_lim, quality);

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

        bl_scoped(interior_forwarding);
        bl_scoped(interior_phases_contract_expand);
        bl_scoped(maxdepth_show_optimized);

        ImGui::Spacing(); ImGui::Spacing();
        ImGui::Text("Interior forwarding");
        if (ImGui::Combo("###MandelInteriorForwarding", &interior_forwarding, MandelMaxDepthOptimizationNames, (int)MandelInteriorForwarding::COUNT))
        {
            switch ((MandelInteriorForwarding)interior_forwarding)
            {
            case MandelInteriorForwarding::SLOWEST: interior_phases_contract_expand = { 10, 0, 10, 0 }; break;
            case MandelInteriorForwarding::SLOW:    interior_phases_contract_expand = { 6, 2, 8, 2 };   break;
            case MandelInteriorForwarding::MEDIUM:  interior_phases_contract_expand = { 5, 3, 7, 3 };   break;
            case MandelInteriorForwarding::FAST:    interior_phases_contract_expand = { 1, 0, 1, 0 };   break;
            default: break;
            }
        }

        ImGui::Checkbox("Highlight optimized regions", &maxdepth_show_optimized);
        ImGui::Spacing(); ImGui::Spacing();

        bl_scoped(kernel_mode);
        ImGui::Text("Kernel Mode");
        ///--------------------------------------------------------------------------------------------------------------------------
        ImGui::RadioButton("Auto",                            &(int&)kernel_mode, (int)KernelMode::AUTO);
        ImGui::RadioButton("No Perturbation",                 &(int&)kernel_mode, (int)KernelMode::NO_PERTURBATION);
        ImGui::RadioButton("Perturbation (no SIMD)",          &(int&)kernel_mode, (int)KernelMode::PERTURBATION);
        ImGui::RadioButton("Perturbation (SIMD)",             &(int&)kernel_mode, (int)KernelMode::PERTURBATION_SIMD);
        ImGui::RadioButton("Perturbation (SIMD, unrolled)",   &(int&)kernel_mode, (int)KernelMode::PERTURBATION_SIMD_UNROLLED);
        ///--------------------------------------------------------------------------------------------------------------------------

        ImGui::EndCollapsingHeaderBox();
    }

    if (ImGui::CollapsingHeaderBox("Normalization Sampling", false))
    {
        bl_scoped(normalize_field_scale);
        bl_scoped(normalize_field_quality);
        bl_scoped(normalize_field_exponent);
        bl_scoped(normalize_field_precision);
        bl_scoped(preview_normalization_field);
        static double init_normalize_field_scale = normalize_field_scale;
        static double init_normalize_field_quality = normalize_field_quality;
        static double init_normalize_field_exponent = normalize_field_exponent;
        static double init_normalize_field_precision = normalize_field_precision;
        ImGui::Checkbox("Preview", &preview_normalization_field);

        ImGui::RevertableSliderDouble("Scale", &normalize_field_scale, &init_normalize_field_scale, 0.5, 4.0, "%.2f");
        ImGui::SetItemTooltip("Controls size of normalization field.\n(samples outside of viewport must be calculated)");

        ImGui::RevertableSliderDouble("Quality", &normalize_field_quality, &init_normalize_field_quality, 0.1, 1.0, "%.2f");
        ImGui::SetItemTooltip("Controls the number of samples taken.\n(Keep low for better performance, especially if Scale is high)");

        ImGui::RevertableSliderDouble("Exponent", &normalize_field_exponent, &init_normalize_field_exponent, 1.0, 4.0, "%.2f");
        ImGui::SetItemTooltip("Higher = Sample more heavily near camera center.\n(useful for zoom animations, otherwise keep at 1.0)");
        
        ImGui::RevertableSliderDouble("Precision", &normalize_field_precision, &normalize_field_precision, 0.0, 1.0, "%.2f", ImGuiSliderFlags_AlwaysClamp);
        ImGui::SetItemTooltip("0 = Always reuse nearest cached pixel where possible (faster)\n1 = Always recalculate at exact sample coordinate (reduces flicker)");

        if (normalize_field_precision < 0.2)
        {
            if (normalize_field_scale > 1.5 && normalize_field_quality > 0.4)
            {
                ImGui::Spacing();
                ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0, 0.2, 0.2, 1.0));
                ImGui::TextUnformatted("Warning:\n   Sampling heavily outside of viewport,\n   performance may start to suffer.");
                ImGui::PopStyleColor();
            }
        }
        else
        {
            if (normalize_field_quality > 0.25)
            {
                ImGui::Spacing(); 
                ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0, 0.2, 0.2, 1.0));
                ImGui::TextUnformatted("Warning:\n   Large sample count,\n   performance may start to suffer.");
                ImGui::PopStyleColor();
            }
        }

        ImGui::EndCollapsingHeaderBox();
    }
    //if (ImGui::BeginSectionBox("AAA", false))
    //{
    //    ImGui::EndSectionBox();
    //}
    //if (ImGui::BeginSectionBox("BBB", false))
    //{
    //    ImGui::EndSectionBox();
    //}
    //if (ImGui::BeginSectionBox("CCC", false))
    //{
    //    ImGui::EndSectionBox();
    //}
}
void Mandelbrot_Scene::UI::populateColorCycleOptions()
{
    if (ImGui::CollapsingHeaderBox("Mandelbrot Parameters", false))
    {
        // Weights / Mix formula
        bl_scoped(iter_weight, dist_weight, stripe_weight);
        bl_scoped(shade_formula);


        double iter_ratio, dist_ratio, stripe_ratio;
        shadingRatios(
            iter_weight, dist_weight, stripe_weight,
            iter_ratio, dist_ratio, stripe_ratio
        );

        char iter_header[64], dist_header[64], stripe_header[64];
        sprintf(iter_header, "Iteration  -  %d%% Weight", (int)(iter_ratio * 100.0));
        sprintf(dist_header, "Distance  -  %d%% Weight", (int)(dist_ratio * 100.0));
        sprintf(stripe_header, "Stripe  -  %d%% Weight", (int)(stripe_ratio * 100.0));


        ImGui::SliderDouble("ITER Weight", &iter_weight, 0.0, 1.0, "%.3f", ImGuiSliderFlags_AlwaysClamp);
        ImGui::SliderDouble("DIST Weight", &dist_weight, 0.0, 1.0, "%.3f", ImGuiSliderFlags_AlwaysClamp);
        ImGui::SliderDouble("STRIPE Weight", &stripe_weight, 0.0, 1.0, "%.3f", ImGuiSliderFlags_AlwaysClamp);

        ImGui::Combo("###MandelFormula", &shade_formula, MandelFormulaNames, (int)MandelShaderFormula::COUNT);

        ImGui::EndCollapsingHeaderBox();

        if (ImGui::BeginTabBar("params_tabs"))
        {
            bl_scoped(iter_weight, dist_weight, stripe_weight);
            bl_view(stats, camera, tweening);

            if (iter_weight > 0.0001f && ImGui::TabBox("Iteration"))
            {
                //ImGui::Spacing();
                //ImGui::GroupBox box("iter_group", nullptr);


                float required_space = 0.0f;
                ImGui::IncreaseRequiredSpaceForLabel(required_space, "% Iterations");

                // ────── Color Cycle: ITER ──────
                {
                    bl_scoped(iter_params);
                    bl_scoped(iter_lim);
                    bl_scoped(quality);
                    bl_scoped(dynamic_iter_lim);

                    if (ImGui::Checkbox("% of Max Iters", &iter_params.cycle_iter_dynamic_limit))
                    {
                        if (iter_params.cycle_iter_dynamic_limit)
                            iter_params.cycle_iter_value /= iter_lim;
                        else
                            iter_params.cycle_iter_value *= iter_lim;
                    }
                    ImGui::SameLine(); ImGui::Dummy(ImVec2(10.0f, 0.0f)); ImGui::SameLine();

                    /// todo
                    //bl_scoped(use_smoothing); 
                    //ImGui::Checkbox("Smooth", &use_smoothing);

                    ImGui::Checkbox("Normalize to Zoom", &iter_params.cycle_iter_normalize_depth);

                    double raw_cycle_iters;
                    if (iter_params.cycle_iter_dynamic_limit)
                    {
                        double cycle_pct = iter_params.cycle_iter_value * 100.0;

                        ImGui::SetNextItemWidthForSpace(required_space);
                        ImGui::SliderDouble("% Iterations", &cycle_pct, 0.001, 100.0, "%.4f%%",
                            (/*color_cycle_use_log1p ?*/ ImGuiSliderFlags_Logarithmic /*: 0*/) |
                            ImGuiSliderFlags_AlwaysClamp);

                        iter_params.cycle_iter_value = cycle_pct / 100.0;

                        raw_cycle_iters = finalIterLimit(camera, quality, dynamic_iter_lim, tweening) * iter_params.cycle_iter_value;

                    }
                    else
                    {
                        ImGui::SetNextItemWidthForSpace(required_space);
                        ImGui::SliderDouble("Iterations", &iter_params.cycle_iter_value, 1.0, (float)iter_lim, "%.3f",
                            (/*color_cycle_use_log1p ?*/ ImGuiSliderFlags_Logarithmic /*: 0*/) |
                            ImGuiSliderFlags_AlwaysClamp);

                        raw_cycle_iters = iter_params.cycle_iter_value;
                    }

                    double log_pct = iter_params.cycle_iter_log1p_weight * 100.0;
                    ImGui::SetNextItemWidthForSpace(required_space);
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

                if (!platform()->is_mobile())
                {
                    ImGui::Spacing();
                    ImGui::SeparatorText("Histogram");
                    {
                        float w = ImGui::GetContentRegionAvail().x;
                        if (ImPlot::BeginPlot("Depth", ImVec2(w, 0)))
                        {
                            ImPlot::SetupAxis(ImAxis_X1, "Value", ImPlotAxisFlags_AutoFit);
                            ImPlot::SetupAxis(ImAxis_Y1, "Pixels", ImPlotAxisFlags_AutoFit);
                            ImPlot::PlotHistogram("##iter_hist",
                                stats.iter_histogram.data(),
                                (int)stats.iter_histogram.size(),
                                100, 1.0);
                            ImPlot::EndPlot();
                        }
                    }
                }
                ImGui::EndTabBox();
            }

            if (dist_weight > 0.0001f && ImGui::TabBox("Distance"))
            {
                ImGui::Spacing();

                float required_width = 0.0f;
                ImGui::IncreaseRequiredSpaceForLabel(required_width, "Sharpness");


                // ────── Color Cycle: DIST ──────
                {
                    bl_scoped(dist_params);

                    ImGui::Checkbox("Invert", &dist_params.cycle_dist_invert);

                    ImGui::SetNextItemWidthForSpace(required_width);
                    ImGui::SliderDouble("Dist", &dist_params.cycle_dist_value, 0.001, 1.0, "%.5f", ImGuiSliderFlags_Logarithmic);

                    ImGui::SetNextItemWidthForSpace(required_width);
                    ImGui::SliderDouble_InvLog("Sharpness", &dist_params.cycle_dist_sharpness, 0.0, 100.0, "%.4f%%",
                        ImGuiSliderFlags_NoRoundToFormat | ImGuiSliderFlags_AlwaysClamp);
                }

                // ────── Tone: DIST ──────
                {
                    bl_scoped(dist_tone_params);

                    ImGui::Spacing();
                    ImGui::SeparatorText("Tone Adjustments");

                    // Contrast (disabled as it does essentially the same thing as "Dist" slider)
                    ///ImGui::SetNextItemWidthForSpace(required_width);
                    ///ImGui::SliderFloat("Contrast", &dist_tone_params.contrast, 0.01f, 10.0f, "%.3f", ImGuiSliderFlags_Logarithmic);

                    // Gamma
                    float gamma_exp = 1.0f / dist_tone_params.gamma;
                    ImGui::SetNextItemWidthForSpace(required_width);
                    ImGui::SliderFloat("Gamma", &gamma_exp, 0.0f, 2.0f, "%.3f");
                    dist_tone_params.gamma = 1.0f / gamma_exp;

                    // Offset
                    ImGui::SetNextItemWidthForSpace(required_width);
                    ImGui::SliderFloat("Offset", &dist_tone_params.brightness, -1.0f, 1.0f, "%.3f");

                    if (!platform()->is_mobile())
                    {
                        // Histogram
                        ImGui::Spacing();
                        ImGui::SeparatorText("Histogram");
                        {
                            float w = ImGui::GetContentRegionAvail().x;
                            if (ImPlot::BeginPlot("Distance", ImVec2(w, 0)))
                            {
                                ImPlot::SetupAxis(ImAxis_X1, "Value", ImPlotAxisFlags_AutoFit);
                                ImPlot::SetupAxis(ImAxis_Y1, "Pixels", ImPlotAxisFlags_AutoFit);
                                ImPlot::PlotHistogram("##dist_hist",
                                    stats.dist_histogram.data(),
                                    (int)stats.dist_histogram.size(),
                                    100, 1.0
                                );
                                ImPlot::EndPlot();
                            }
                        }
                    }
                }

                ImGui::EndTabBox();
            }

            if (stripe_weight > 0.0001f && ImGui::TabBox("Stripe"))
            {
                ImGui::Spacing();

                float required_width = 0.0f;
                ImGui::IncreaseRequiredSpaceForLabel(required_width, "Frequency");

                // ────── Color Cycle: STRIPE ──────
                {
                    bl_scoped(stripe_params);

                    ImGui::SetNextItemWidthForSpace(required_width);
                    ImGui::SliderFloat("Frequency", &stripe_params.freq, 1.0f, 100.0f, "%.0f", ImGuiSliderFlags_AlwaysClamp);

                    ImGui::SetNextItemWidthForSpace(required_width);
                    ImGui::SliderAngle("Phase", &stripe_params.phase, 0.0f, Math::FPI * 2.0f, 0, ImGuiSliderFlags_AlwaysClamp);
                }

                // ────── Tone: STRIPE ──────
                {
                    bl_scoped(stripe_tone_params);

                    ImGui::Spacing();
                    ImGui::SeparatorText("Tone Adjustments");

                    // Contrast
                    ImGui::SetNextItemWidthForSpace(required_width);
                    ImGui::SliderFloat("Contrast", &stripe_tone_params.contrast, 0.01f, 10.0f, "%.3f", ImGuiSliderFlags_Logarithmic);

                    // Gamma
                    float gamma_exp = 1.0f / stripe_tone_params.gamma;
                    ImGui::SetNextItemWidthForSpace(required_width);
                    ImGui::SliderFloat("Gamma", &gamma_exp, 0.0f, 2.0f, "%.3f");
                    stripe_tone_params.gamma = 1.0f / gamma_exp;

                    // Offset
                    ImGui::SetNextItemWidthForSpace(required_width);
                    ImGui::SliderFloat("Offset", &stripe_tone_params.brightness, -1.0f, 1.0f, "%.3f");

                    if (!platform()->is_mobile())
                    {
                        // Histogram
                        ImGui::Spacing();
                        ImGui::SeparatorText("Histogram");
                        {
                            ///ImGui::Text("Mean STRIPE: %.6f", stats.field_info.mean_stripe);
                            float w = ImGui::GetContentRegionAvail().x;
                            if (ImPlot::BeginPlot("Stripe", ImVec2(w, 0)))
                            {
                                ImPlot::SetupAxis(ImAxis_X1, "Value", ImPlotAxisFlags_AutoFit);
                                ImPlot::SetupAxis(ImAxis_Y1, "Pixels", ImPlotAxisFlags_AutoFit);
                                ImPlot::PlotHistogram("##stripe_hist",
                                    stats.stripe_histogram.data(),
                                    (int)stats.stripe_histogram.size(),
                                    100, 1.0
                                );
                                ImPlot::EndPlot();
                            }
                        }
                    }
                }

                ImGui::EndTabBox();
            }

            ImGui::EndTabBar();
        }
    }
}
void Mandelbrot_Scene::UI::populateGradientShiftOptions()
{
    if (ImGui::CollapsingHeaderBox("Gradient Offset + Animation", false))
    {
        // Shift
        bl_view(gradient);
        bl_pull(gradient_shifted); // no need to push (calculated already in viewportProcess)
        bl_scoped(hue_shift, gradient_shift);

        // Animation
        bl_scoped(show_color_animation_options);
        bl_scoped(gradient_shift_step, hue_shift_step);

        ImGui::Checkbox("Animate", &show_color_animation_options);

        float required_space = 0.0f;
        ImGui::IncreaseRequiredSpaceForLabel(required_space, "Gradient shift");

        ImGui::SetNextItemWidthForSpace(required_space);
        if (ImGui::DragDouble("Gradient shift", &gradient_shift, 0.01, -100.0, 100.0, " %.3f", ImGuiSliderFlags_AlwaysClamp))
        {
            gradient_shift = Math::wrap(gradient_shift, 0.0, 1.0);
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
        double hue_shift_rad = Math::toRadians(hue_shift);
        static double initial_hue_shift_rad = hue_shift_rad;
        if (ImGui::RevertableSliderAngle("Hue shift", &hue_shift_rad, &initial_hue_shift_rad, 0.0, Math::TWO_PI, 3))
        {
            hue_shift = Math::toDegrees(hue_shift_rad);
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

        ImGui::EndCollapsingHeaderBox();
    }
}
void Mandelbrot_Scene::UI::populateGradientPicker()
{
    if (ImGui::CollapsingHeaderBox("Base Color Gradient", false))
    {
        bl_scoped(gradient);

        ImGui::Text("Load Preset");
        static int selecting_template = -1;
        if (ImGui::Combo("###ColorTemplate", &selecting_template, ColorGradientNames, (int)GradientPreset::COUNT))
        {
            generateGradientFromPreset(gradient, (GradientPreset)selecting_template);
            selecting_template = -1;
        }

        #if MANDEL_DEV_MODE
        if (ImGui::Button("Copy gradient C++ code"))
        {
            ImGui::SetClipboardText(gradient.to_cpp_marks().c_str());
        }
        #endif

        ImGui::Dummy(scale_size(0, 8));

        if (ImGui::GradientEditor(&gradient,
            platform()->dpr(),
            platform()->dpr() * platform()->thumbScale()))
        {
            // Shift
            ///bl_pull(gradient_shifted, hue_shift, gradient_shift);
            ///transformGradient(gradient_shifted, gradient, (float)gradient_shift, (float)hue_shift);
        }

        ImGui::Spacing();
        ImGui::EndCollapsingHeaderBox();
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
        bl_scoped(flatten, flatten_amount);
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

            bl_scoped(x_spline, y_spline);

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

    ImGui::SeparatorText("Lift Tween");
    ImSpline::SplineEditor("tween_zoom_lift", &tween_zoom_lift_spline, &vr);

    ImGui::SeparatorText("Base Zoom Tween");
    ImSpline::SplineEditor("tween_base_zoom", &tween_base_zoom_spline, &vr);

    //ImGui::InputTextMultiline("###pos_buf", pos_tween_buf, 1024, ImVec2(0, 0), ImGuiInputTextFlags_AllowTabInput));
    if (ImGui::Button("Copy position spline")) ImGui::SetClipboardText(tween_pos_spline.serialize(SplineSerializationMode::CPP_ARRAY, 3).c_str());

    //ImGui::InputTextMultiline("###pos_buf", zoom_tween_buf, 1024, ImVec2(0, 0), ImGuiInputTextFlags_AllowTabInput));
    if (ImGui::Button("Copy lift spline")) ImGui::SetClipboardText(tween_zoom_lift_spline.serialize(SplineSerializationMode::CPP_ARRAY, 3).c_str());

    if (ImGui::Button("Copy base zoom spline")) ImGui::SetClipboardText(tween_base_zoom_spline.serialize(SplineSerializationMode::CPP_ARRAY, 3).c_str());
}

struct Limits2D
{
    double xmin, xmax;
    double ymin, ymax;
};

Limits2D computeTaperedFitLimits(
    const std::vector<double>& xs,
    const std::vector<double>& ys,
    double hard_frac = 0.70,
    double pad_frac = 0.3)
{
    const size_t n = xs.size();
    if (n == 0) return { -1, 1, -1, 1 };

    size_t k = (size_t)std::floor(hard_frac * (double)n);
    if (k < 1) k = 1;
    if (k > n) k = n;

    // hard include up to hard_frac
    double base_xmin = xs[0], base_xmax = xs[0];
    double base_ymin = ys[0], base_ymax = ys[0];
    for (size_t i = 0; i < k; ++i) {
        base_xmin = std::min(base_xmin, xs[i]);
        base_xmax = std::max(base_xmax, xs[i]);
        base_ymin = std::min(base_ymin, ys[i]);
        base_ymax = std::max(base_ymax, ys[i]);
    }

    double xmin = base_xmin, xmax = base_xmax;
    double ymin = base_ymin, ymax = base_ymax;

    // tapered influence
    if (k < n) 
    {
        const double denom = (double)(n - 1 - k);
        for (size_t i = k; i < n; ++i) 
        {
            const double t = (denom > 0.0) ? (double)(i - k) / denom : 1.0;
            const double w = 1.0 - t;

            const double x = xs[i];
            const double y = ys[i];

            if (x > base_xmax) xmax = std::max(xmax, base_xmax + (x - base_xmax) * w);
            if (x < base_xmin) xmin = std::min(xmin, base_xmin + (x - base_xmin) * w);

            if (y > base_ymax) ymax = std::max(ymax, base_ymax + (y - base_ymax) * w);
            if (y < base_ymin) ymin = std::min(ymin, base_ymin + (y - base_ymin) * w);
        }
    }

    // padding
    const double dx = (xmax - xmin);
    const double dy = (ymax - ymin);
    const double padx = (dx > 0.0) ? dx * pad_frac : 1.0;
    const double pady = (dy > 0.0) ? dy * pad_frac : 1.0;

    xmin -= padx; xmax += padx;
    ymin -= pady; ymax += pady;

    return { xmin, xmax, ymin, ymax };
}

template<typename T>
void computeOrbit(f128 x0, f128 y0, int iter_lim, std::vector<f64>& xs, std::vector<f64>& ys)
{
    using cplx = cplx<T>;
    constexpr T escape_r = T(16.0);

    cplx z{ 0.0, 0.0 };
    cplx c{ T(x0), T(y0) };

    for (int i = 0; i < iter_lim; i++)
    {
        mandel_step(z, c);
        xs.push_back((f64)z.x);
        ys.push_back((f64)z.y);

        T r2 = z.x * z.x + z.y * z.y;
        if (r2 > escape_r)
            break;
    }
}

void Mandelbrot_Scene::UI::populateMouseOrbit()
{
    if (ImGui::CollapsingHeaderBox("Mouse Orbit", false))
    {
        bl_view(camera);
        bl_view(stats);
        bl_view(iter_lim);
        bl_view(mandel_features);

        ImGui::SliderDouble("Padding", &orbit_padding, 0.0, 1.0, "%.2f", ImGuiSliderFlags_AlwaysClamp);

        std::vector<f64> xs, ys;
        xs.push_back(0.0);
        ys.push_back(0.0);

        auto [x0, y0] = stats.hovered_field_world_pos;

        FloatingPointType float_type = getRequiredFloatType(mandel_features, camera.relativeZoom<f128>());
        table_invoke(dispatch_table(computeOrbit, x0, y0, iter_lim, xs, ys), float_type);

        if (ImPlot::BeginPlot("Mouse orbit"))
        {
            auto lims = computeTaperedFitLimits(xs, ys, 0.8, orbit_padding);
            ImPlot::SetupAxisLimits(ImAxis_X1, lims.xmin, lims.xmax, ImGuiCond_Always);
            ImPlot::SetupAxisLimits(ImAxis_Y1, lims.ymin, lims.ymax, ImGuiCond_Always);
            if (xs.size())
            {
                ImPlot::PlotLine("##mouse_plot", xs.data(), ys.data(), (int)xs.size());
            }
            ImPlot::EndPlot();
        }

        ImGui::EndCollapsingHeaderBox();
    }
}

SIM_END;