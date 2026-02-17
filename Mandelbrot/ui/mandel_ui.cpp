#include "../Mandelbrot.h"

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

void Mandelbrot_Scene::UI::init()
{
    ImGuiStyle& style = ImGui::GetStyle();
    style.HoverStationaryDelay = 0.50f;

    bl_scoped(bookmark_manager);
    bookmark_manager.loadCategoryDirs("/data/bookmarks");
    bookmark_manager.ensureCategory("Examples");

    #if MANDEL_UPDATE_ALL_BOOKMARKS
    {
        for (int i = 0; i < bookmark_manager.size(); i++)
        {
            auto& list = bookmark_manager.at(i).second;
            std::vector<MandelBookmark>& items = list.getItems();
            for (auto& bookmark : items)
            {
                MandelState state;
                state.deserialize(bookmark.data);
                bookmark.data = state.serialize();
            }
        }
        bookmark_manager.saveAll("Updated");
    }
    #endif

    #if MANDEL_RUN_DICTIONARY_TUNINGS
    {
        int current_dictionary = 2;

        DictTuneConfig cfg;
        cfg.quality = 11;
        cfg.window = 22;
        cfg.score_mode = DictScoreMode::Sum;
        cfg.start_with_all_tokens = true;
        cfg.enable_pair_removal = true;
    
        auto tuned = tune_brotli_dictionary(bookmark_manager.shaders, MandelState::getDictionaryTokens(current_dictionary), cfg);
        std::string var_name = "kDictTokensV";
        var_name += std::to_string(current_dictionary);
        blPrint() << tuned.to_initializer(var_name) << "\n";
    }
    #endif

    auto pal = TextEditor::GetDarkPalette();
    pal[(int)TextEditor::PaletteIndex::PreprocIdentifier] = IM_COL32(255, 80, 80, 255);
    editor.SetPalette(pal);

    auto lang = TextEditor::LanguageDefinition::GLSL();

    // shader-pass tags
    lang.mTokenRegexStrings.insert(lang.mTokenRegexStrings.begin(),
        { "@[^\\r\\n]*", TextEditor::PaletteIndex::PreprocIdentifier });

    static constexpr const char* kGlslTypes[] = {
        // scalars
        "void", "bool", "int", "uint", "float", "double",

        // vectors
        "vec2", "vec3", "vec4",
        "bvec2", "bvec3", "bvec4",
        "ivec2", "ivec3", "ivec4",
        "uvec2", "uvec3", "uvec4",
        "dvec2", "dvec3", "dvec4",

        // matrices
        "mat2", "mat3", "mat4",
        "mat2x2", "mat2x3", "mat2x4",
        "mat3x2", "mat3x3", "mat3x4",
        "mat4x2", "mat4x3", "mat4x4",
        "dmat2", "dmat3", "dmat4",
        "dmat2x2", "dmat2x3", "dmat2x4",
        "dmat3x2", "dmat3x3", "dmat3x4",
        "dmat4x2", "dmat4x3", "dmat4x4",

        // samplers (add/remove as needed)
        "sampler1D", "sampler2D", "sampler3D", "samplerCube",
        "sampler2DShadow", "samplerCubeShadow",
        "sampler2DArray", "sampler2DArrayShadow",
        "isampler2D", "usampler2D",
    };

    for (const char* t : kGlslTypes)
        lang.mKeywords.insert(t);

    auto& rx = lang.mTokenRegexStrings;

    auto itIdentifierRule = std::find_if(rx.begin(), rx.end(),
        [](const TextEditor::LanguageDefinition::TokenRegexString& r)
    {
        return r.first == "[a-zA-Z_][a-zA-Z0-9_]*";
    });

    // 2) iter/dist/stripe
    auto itIdent = std::find_if(rx.begin(), rx.end(),
        [](const TextEditor::LanguageDefinition::TokenRegexString& r)
    {
        return r.first == "[a-zA-Z_][a-zA-Z0-9_]*";
    });

    rx.insert(itIdent, { R"(\b(iter|dist|stripe)\b)", TextEditor::PaletteIndex::KnownIdentifier });

    // functions
    {
        TextEditor::Identifier id; id.mDeclaration = "boundary(falloff_begin, falloff_end [optional])";
        lang.mIdentifiers.insert(std::make_pair(std::string("boundary"), id));
        editor.SetPalette(pal);
    }

    editor_height = scale_size(250.0f);
    editor.SetShowWhitespaces(false);
    editor.SetLanguageDefinition(lang);
    editor.SetText(Mandelbrot_Scene::init_shader_source);
}

void Mandelbrot_Scene::UI::sidebar()
{
    bl_scoped(tweening);

    if (tweening)
        ImGui::BeginCollapsingHeaderContentsDisabled();

    ImGui::SeparatorText("Controls");
    populateSavingLoading();
    populateExamples();

    ImGui::SeparatorText("Active State");
    populateCameraView();
    populateInputOptions();
    populateShaderEditor();
    populateGradientOptions();
    populateAnimation();
    populateQualityOptions();
    //populateGradientPicker();
    populateCaptureOptions();
    //populateExperimental();

    ImGui::SeparatorText("Other");
    populateStats();
    if (!platform()->is_mobile())
    {
        populateMouseOrbit();

        #ifdef MANDEL_DEV_MODE
        populateHistory();
        #endif
    }

    if (tweening)
        ImGui::EndCollapsingHeaderContentsDisabled();

    // ======== Developer ========
    #if MANDEL_DEV_EDIT_TWEEN_SPLINES
    populateSplinesDev();
    #endif
}

void Mandelbrot_Scene::UI::overlay()
{
}

void Mandelbrot_Scene::UI::launchBookmark(std::string_view data)
{
    bl_schedule([data](Mandelbrot_Scene& scene)
    {
        scene.loadState(data, true);
    });
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
                data_buf = scene.serializeState();
            });

            ImGui::OpenPopup("Save Data"); // open on this frame
        }

        ImGui::SameLine();
        if (ImGui::Button("Load"))
        {
            data_buf = "";
            show_load_dialog = true;
            ImGui::OpenPopup("Load Data");
        }

        ImGui::SameLine();
        if (ImGui::Button("Share"))
        {
            show_share_dialog = true;
            bl_schedule([&](Mandelbrot_Scene& scene)
            {
                url_buf = dataToURL(scene.serializeState());
            });
            ImGui::OpenPopup("Share URL");
        }

        // Save Dialog
        ImGui::SetNextWindowSize(scale_size(350, 300), ImGuiCond_FirstUseEver);
        if (ImGui::BeginPopupModal("Save Data", &show_save_dialog))
        {
            ImVec2 avail = ImGui::GetContentRegionAvail();
            avail.y -= ImGui::GetFrameHeightWithSpacing(); // leave room for buttons/input

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
                    scene.loadState(data_buf, true);
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
    bl_scoped(camera);

    //static bool show_view_by_default = !platform()->is_mobile(); // Expect navigate by touch for mobile
    if (ImGui::CollapsingHeaderBox("View", true))
    {
        bl_scoped(show_axis, display_alignment_overlay);
        ImGui::Checkbox("Show Axis", &show_axis);
        ImGui::SameLine();
        ImGui::Dummy(scale_size(4, 0));
        ImGui::SameLine();
        ImGui::Checkbox("Alignment Lines", &display_alignment_overlay);

        #if MANDEL_FEATURE_INTERACTIVE_CARDIOID
        bl_scoped(show_interactive_cardioid);
        ImGui::Checkbox("Interactive Cardioid", &show_interactive_cardioid);
        #endif

        camera.populateUI({ -5.0, -5.0, 5.0, 5.0 });

        #if MANDEL_FEATURE_INTERACTIVE_CARDIOID
        if (show_interactive_cardioid)
        {
            ImGui::BeginLabelledBox("Interactive Cardioid");

            bl_scoped(ani_angle);
            bl_scoped(ani_inc);
            bl_scoped(animate_cardioid_angle);

            ImGui::Checkbox("Apply Animation", &animate_cardioid_angle);
            ImGui::SliderAngle("Angle", &ani_angle);
            if (animate_cardioid_angle)
                ImGui::SliderAngle("Step", &ani_inc, -math::half_pi, math::half_pi, 1);

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
        bl_scoped(bookmark_manager);

        const float sx = ImGui::GetStyle().ItemSpacing.x;
        const float sy = ImGui::GetStyle().ItemSpacing.y;

        const float btn_w = 128.0f, btn_h = 72.0f;
        const float btn_sw = scale_size(btn_w), btn_sh = scale_size(btn_h);

        // red headers
        ImGui::PushStyleColor(ImGuiCol_Header, ImVec4(0.5f, 0.2f, 0.2f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_HeaderHovered, ImVec4(0.6f, 0.3f, 0.3f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_HeaderActive, ImVec4(1.0f, 0.4f, 0.4f, 1.0f));

        {
            for (int category_i = 0; category_i < bookmark_manager.size(); category_i++)
            {
                auto [category_name, list] = bookmark_manager.at(category_i);
                auto& bookmarks = list.getItems();

                ImGui::PushID(category_i);
                if (ImGui::CollapsingHeaderBox(category_name.c_str()))
                {
                    const int count = (int)bookmarks.size();
                    const int rows_visible = 3;
                    const float strip_h = rows_visible * btn_sh + (rows_visible - 1) * sy;

                    if (ImGui::BeginChild("##examples_region", ImVec2(0.0f, strip_h), 0))
                    {
                        const float window_visible_x2 = ImGui::GetWindowPos().x + ImGui::GetWindowContentRegionMax().x;

                        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
                        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0, 0, 0, 0));
                        ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0, 0, 0, 0));
                        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0, 0));

                        for (int idx = 0; idx < count; ++idx)
                        {
                            MandelBookmark& bookmark = bookmarks[idx];

                            ImGui::PushID(idx);
                            if (ImGui::ImageButton("##bookmark", bookmark.thumbTexture(), ImVec2(btn_sw, btn_sh)))
                                launchBookmark(bookmark.data);
                            ImGui::PopID();

                            const float last_button_x2 = ImGui::GetItemRectMax().x;
                            const float next_button_x2 = last_button_x2 + sx + btn_sw;

                            if (idx + 1 < count && next_button_x2 < window_visible_x2)
                                ImGui::SameLine(0.0f, sx);
                        }

                        ImGui::PopStyleVar();
                        ImGui::PopStyleColor(3);

                    }
                    ImGui::EndChild();

                    bool allow_add_new = false;

                    #if MANDEL_DEV_MODE
                    allow_add_new = true;
                    #endif

                    #ifndef BITLOOP_DEV_MODE
                    // if RELEASE build, don't allow modifying bundled example bookmarks, even if in mandel dev mode
                    if (category_name == "Examples")
                        allow_add_new = false;
                    #endif

                    if (allow_add_new)
                    {
                        if (ImGui::Button("Bookmark Active"))
                        {
                            // Generate bookmark on worker thread (since we also need to generate a thumbnail)
                            bl_schedule([&, category_i, category_name](Mandelbrot_Scene& scene)
                            {
                                // Serialize Mandelbrot state
                                std::string state_data = scene.serialize();

                                // Grab a suitable thumbnail preset
                                SnapshotPresetList all_presets = main_window()->getSnapshotPresetManager()->allPresets();
                                CapturePreset* preset = all_presets.findByAlias("thumb128x72_hd");
                                assert(preset != nullptr);

                                std::string thumb_path = ProjectBase::activeProject()->root_path(
                                    "data/bookmarks/" + category_name + "/" + MandelBookmark(state_data).thumbFilename()
                                );

                                // Create bookmark, generate thumbnail image, load direct from memory (also saves to data/thumbnails/ for embedding)
                                scene.beginSnapshot(*preset, thumb_path, [&, state_data, category_i](bytebuf& thumb_data, const CapturePreset& preset)
                                {
                                    MandelBookmark bookmark(state_data);

                                    // Load thumbnail from memory
                                    bookmark.loadThumbnail(thumb_data, preset.width(), preset.height());

                                    // Finally, add bookmark to the target list
                                    MandelBookmarkList& list = bookmark_manager.at(category_i).second; /// todo: reaccessing UI from worker, maybe not thread-safe
                                    list.addItem(bookmark);
                                });
                            });
                        }
                    }

                    ImGui::EndCollapsingHeaderBox();
                }
                ImGui::PopID();
            }
        }

        ImGui::PopStyleColor(3);

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

            ImGui::Text("  = %d Iters", finalIterLimit(camera, quality, dynamic_iter_lim, tweening));
        }
        else
        {
            ImGui::DragDouble("###Quality", &quality, 1000.0, 1.0, 1000000.0, "%.0f", ImGuiSliderFlags_Logarithmic);
            ImGui::Text("  = %d Iters", finalIterLimit(camera, quality, dynamic_iter_lim, tweening));
        }

        if (!platform()->is_mobile())
        {
            bl_scoped(interior_forwarding);
            bl_scoped(contract_expand_phases);
            bl_scoped(maxdepth_show_optimized);

            ImGui::Spacing();
            ImGui::Spacing();
            ImGui::Spacing();
            ImGui::Text("Interior forwarding");
            if (ImGui::Combo("###MandelInteriorForwarding", &interior_forwarding, MandelMaxDepthOptimizationNames, (int)MandelInteriorForwarding::COUNT))
            {
                using CE = ContractExpandPhase;
                switch ((MandelInteriorForwarding)interior_forwarding)
                {
                //case MandelInteriorForwarding::SLOWEST: 
                //    contract_expand_phases = { CE{16, 0}, CE{16, 0}, CE{16, 0} }; break;
                
                case MandelInteriorForwarding::SLOW:
                    contract_expand_phases = { CE{7, 0}, CE{6, 2}, CE{8, 2} };   break;
                
                case MandelInteriorForwarding::MEDIUM:
                    contract_expand_phases = { CE{5, 0}, CE{5, 3}, CE{7, 3} };   break;
                
                case MandelInteriorForwarding::FAST:
                    contract_expand_phases = { CE{2, 0}, CE{1, 0}, CE{1, 0} };   break;
                
                default: break;
                }
            }

            ImGui::Checkbox("Highlight optimized regions", &maxdepth_show_optimized);
            ImGui::Spacing(); ImGui::Spacing();

            #if MANDEL_DEV_MODE
            bl_scoped(kernel_mode);
            ImGui::Text("Kernel Mode");
            ///--------------------------------------------------------------------------------------------------------------------------
            ImGui::RadioButton("Auto", &(int&)kernel_mode, (int)KernelMode::AUTO);
            ImGui::RadioButton("No Perturbation", &(int&)kernel_mode, (int)KernelMode::NO_PERTURBATION);
            ImGui::RadioButton("Perturbation (no SIMD)", &(int&)kernel_mode, (int)KernelMode::PERTURBATION);
            ImGui::RadioButton("Perturbation (SIMD)", &(int&)kernel_mode, (int)KernelMode::PERTURBATION_SIMD);
            ImGui::RadioButton("Perturbation (SIMD, unrolled)", &(int&)kernel_mode, (int)KernelMode::PERTURBATION_SIMD_UNROLLED);
            ///--------------------------------------------------------------------------------------------------------------------------
            #endif
        }

        ImGui::EndCollapsingHeaderBox();
    }

    if (ImGui::CollapsingHeaderBox("Normalization Field", false))
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
        ImGui::SetItemTooltip("Radius of normalization field");

        ImGui::RevertableSliderDouble("Quality", &normalize_field_quality, &init_normalize_field_quality, 0.1, 1.0, "%.2f");
        ImGui::SetItemTooltip("Controls the number of samples taken");

        ImGui::RevertableSliderDouble("Exponent", &normalize_field_exponent, &init_normalize_field_exponent, 1.0, 4.0, "%.2f");
        ImGui::SetItemTooltip("Higher = Sample more heavily near camera center.\n(useful for zoom animations, otherwise keep at 1.0)");
        
        ImGui::RevertableSliderDouble("Accuracy", &normalize_field_precision, &normalize_field_precision, 0.0, 1.0, "%.2f", ImGuiSliderFlags_AlwaysClamp);
        ImGui::SetItemTooltip("0 = Reuse nearest cached pixel where possible (fast)\n1 = Force calculate at exact sample coordinate (less flicker)");

        if (normalize_field_precision < 0.2)
        {
            if (normalize_field_scale > 1.5 && normalize_field_quality > 0.4)
            {
                ImGui::Spacing();
                ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.2f, 0.2f, 1.0f));
                ImGui::TextUnformatted("Warning:\n   Sampling heavily outside of viewport,\n   performance may start to suffer.");
                ImGui::PopStyleColor();
            }
        }
        else
        {
            if (normalize_field_quality > 0.25)
            {
                ImGui::Spacing(); 
                ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.2f, 0.2f, 1.0f));
                ImGui::TextUnformatted("Warning:\n   Large sample count,\n   performance may start to suffer.");
                ImGui::PopStyleColor();
            }
        }

        ImGui::EndCollapsingHeaderBox();
    }
}
void Mandelbrot_Scene::UI::populateInputOptions()
{
    if (ImGui::CollapsingHeaderBox("Inputs", false))
    {
        // Weights / Mix formula
        bl_scoped(iter_weight, dist_weight, stripe_weight);

        f32 iter_ratio, dist_ratio, stripe_ratio;
        shadingRatios(
            iter_weight, dist_weight, stripe_weight,
            iter_ratio, dist_ratio, stripe_ratio
        );

        char iter_header[32], dist_header[32], stripe_header[32];
        sprintf(iter_header, "iter - %d%%###iter_tab", (int)(iter_ratio * 100.0f));
        sprintf(dist_header, "dist - %d%%###dist_tab", (int)(dist_ratio * 100.0f));
        sprintf(stripe_header, "stripe - %d%%###stripe_tab", (int)(stripe_ratio * 100.0f));

        ImGui::SeparatorText("Features");
        ImGui::SliderDouble("iter", &iter_weight, 0.0, 1.0, "%.3f", ImGuiSliderFlags_AlwaysClamp);
        ImGui::SliderDouble("dist", &dist_weight, 0.0, 1.0, "%.3f", ImGuiSliderFlags_AlwaysClamp);
        ImGui::SliderDouble("stripe", &stripe_weight, 0.0, 1.0, "%.3f", ImGuiSliderFlags_AlwaysClamp);

        ImGui::Spacing();
        ImGui::Spacing();
        if (ImGui::BeginTabBar("params_tabs"))
        {
            bl_scoped(iter_weight, dist_weight, stripe_weight);
            bl_view(stats, camera, tweening);
            bl_scoped(iter_hist_visible, dist_hist_visible, stripe_hist_visible);

            ImGui::PushID("iter_tab");
            if (iter_weight > 0.0001f && ImGui::TabBox(iter_header))
            {
                float required_space = 0.0f;
                ImGui::IncreaseRequiredSpaceForLabel(required_space, "% Iterations");

                iter_hist_visible = false;
                dist_hist_visible = false;
                stripe_hist_visible = false;

                // ────── Color Cycle: ITER ──────
                {
                    bl_scoped(iter_params);
                    bl_scoped(iter_lim);
                    bl_scoped(quality);
                    bl_scoped(dynamic_iter_lim);

                    if (ImGui::Checkbox("% of Max Iters", &iter_params.iter_dynamic_limit))
                    {
                        if (iter_params.iter_dynamic_limit)
                            iter_params.iter_cycle_value /= iter_lim;
                        else
                            iter_params.iter_cycle_value *= iter_lim;
                    }

                    /// todo
                    //bl_scoped(use_smoothing); 
                    //ImGui::Checkbox("Smooth", &use_smoothing);

                    ImGui::Checkbox("Use Normalization Field", &iter_params.iter_normalize_depth);

                    double raw_cycle_iters;
                    if (iter_params.iter_dynamic_limit)
                    {
                        double cycle_pct = iter_params.iter_cycle_value * 100.0;

                        ImGui::SetNextItemWidthForSpace(required_space);
                        ImGui::SliderDouble("% Iterations", &cycle_pct, 0.001, 100.0, "%.4f%%",
                            (/*color_cycle_use_log1p ?*/ ImGuiSliderFlags_Logarithmic /*: 0*/) |
                            ImGuiSliderFlags_AlwaysClamp);

                        iter_params.iter_cycle_value = cycle_pct / 100.0;

                        raw_cycle_iters = finalIterLimit(camera, quality, dynamic_iter_lim, tweening) * iter_params.iter_cycle_value;

                    }
                    else
                    {
                        ImGui::SetNextItemWidthForSpace(required_space);
                        ImGui::SliderDouble("Iterations", &iter_params.iter_cycle_value, 1.0, (float)iter_lim, "%.3f",
                            (/*color_cycle_use_log1p ?*/ ImGuiSliderFlags_Logarithmic /*: 0*/) |
                            ImGuiSliderFlags_AlwaysClamp);

                        raw_cycle_iters = iter_params.iter_cycle_value;
                    }

                    double log_pct = iter_params.iter_log1p_weight * 100.0;
                    ImGui::SetNextItemWidthForSpace(required_space);
                    if (ImGui::SliderDouble_InvLog("Logarithmic", &log_pct, 0.0, 100.0, "%.2f%%", ImGuiSliderFlags_NoRoundToFormat | ImGuiSliderFlags_AlwaysClamp))
                    {
                        iter_params.iter_log1p_weight = log_pct / 100.0;
                    }

                    ///// Print cycle iter values
                    ///{
                    ///    if (iter_dynamic_limit)
                    ///        ImGui::Text("raw_cycle_iters = %.1f", raw_cycle_iters);
                    ///
                    ///    double log_cycle_iters = math::linearLog1pLerp(raw_cycle_iters, iter_log1p_weight);
                    ///    ImGui::Text("log_cycle_iters = %.1f", log_cycle_iters);
                    ///}
                }

                if (!platform()->is_mobile())
                {
                    ImGui::Spacing();

                    if (!show_iter_hist) {
                        if (ImGui::Button("Show Histogram"))
                            show_iter_hist = true;
                    }
                    else if (ImGui::Button("Hide Histogram"))
                        show_iter_hist = false;
                    
                    if (show_iter_hist)
                    {
                        float w = ImGui::GetContentRegionAvail().x;
                        if (ImPlot::BeginPlot("Depth", ImVec2(w, 0)))
                        {
                            iter_hist_visible = true;

                            ImPlot::SetupAxis(ImAxis_X1, "Value", ImPlotAxisFlags_AutoFit);
                            ImPlot::SetupAxis(ImAxis_Y1, "Pixels", ImPlotAxisFlags_AutoFit);
                            ImPlot::PlotHistogram("##iter_hist",
                                stats.iter_histogram.data(),
                                (int)stats.iter_histogram.size(),
                                100, 1.0);

                            //ImPlot::SetNextLineStyle(ImVec4(1, 0, 0, 1), 2.0f);
                            //ImPlot::PlotInfLines("#min_depth_v", &stats.field_info.final_min_depth, 1);
                            //ImPlot::SetNextLineStyle(ImVec4(1, 0, 0, 1), 2.0f);
                            //ImPlot::PlotInfLines("#max_depth_v", &stats.field_info.final_max_depth, 1);
                            ImPlot::EndPlot();
                        } 
                    }
                }
                ImGui::EndTabBox();
            }
            ImGui::PopID();

            if (dist_weight > 0.0001f && ImGui::TabBox(dist_header))
            {
                ImGui::Spacing();

                float required_width = 0.0f;
                ImGui::IncreaseRequiredSpaceForLabel(required_width, "Sharpness");


                // ────── Color Cycle: DIST ──────
                {
                    bl_scoped(dist_params);

                    ImGui::Checkbox("Invert", &dist_params.dist_invert);

                    ImGui::SetNextItemWidthForSpace(required_width);
                    ImGui::SliderDouble("Distance", &dist_params.dist_cycle_value, 0.001, 1.0, "%.5f", ImGuiSliderFlags_Logarithmic);

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
                    ImGui::SliderFloat("Gamma", &gamma_exp, 0.001f, 2.0f, "%.3f");
                    dist_tone_params.gamma = 1.0f / gamma_exp;

                    // Offset
                    ImGui::SetNextItemWidthForSpace(required_width);
                    ImGui::SliderFloat("Offset", &dist_tone_params.brightness, -1.0f, 1.0f, "%.3f");
                }

                if (!platform()->is_mobile())
                {
                    // Histogram
                    ImGui::Spacing();

                    if (!show_dist_hist) {
                        if (ImGui::Button("Show Histogram"))
                            show_dist_hist = true;
                    }
                    else if (ImGui::Button("Hide Histogram"))
                        show_dist_hist = false;

                    if (show_dist_hist)
                    {
                        float w = ImGui::GetContentRegionAvail().x;
                        if (ImPlot::BeginPlot("Dist", ImVec2(w, 0)))
                        {
                            dist_hist_visible = true;

                            ImPlot::SetupAxis(ImAxis_X1, "Value", ImPlotAxisFlags_AutoFit);
                            ImPlot::SetupAxis(ImAxis_Y1, "Pixels", ImPlotAxisFlags_AutoFit);
                            ImPlot::PlotHistogram("##dist_hist",
                                stats.dist_histogram.data(),
                                (int)stats.dist_histogram.size(),
                                100, 1.0
                            );
                            //ImPlot::SetNextLineStyle(ImVec4(1, 0, 0, 1), 2.0f);
                            //ImPlot::PlotInfLines("#min_dist_v", &stats.field_info.final_min_dist, 1);
                            //ImPlot::SetNextLineStyle(ImVec4(1, 0, 0, 1), 2.0f);
                            //ImPlot::PlotInfLines("#max_dist_v", &stats.field_info.final_max_dist, 1);
                            ImPlot::EndPlot();
                        }
                    }
                }

                ImGui::EndTabBox();
            }

            if (stripe_weight > 0.0001f && ImGui::TabBox(stripe_header))
            {
                ImGui::Spacing();

                float required_width = 0.0f;
                ImGui::IncreaseRequiredSpaceForLabel(required_width, "Frequency");

                //bl_scoped(stripe_mag_from_hist);
                //ImGui::Checkbox("Histogram normalize", &stripe_mag_from_hist);

                // ────── Color Cycle: STRIPE ──────
                {
                    bl_scoped(stripe_params);

                    ImGui::SetNextItemWidthForSpace(required_width);
                    ImGui::SliderInt("Frequency", &stripe_params.freq, 1, 20, "%d", ImGuiSliderFlags_AlwaysClamp);

                    ImGui::SetNextItemWidthForSpace(required_width);
                    ImGui::SliderAngle("Phase", &stripe_params.phase, 0.0f, math::pi_f * 2.0f, 0, ImGuiSliderFlags_AlwaysClamp);
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
                    ImGui::SliderFloat("Gamma", &gamma_exp, 0.001f, 2.0f, "%.3f");
                    stripe_tone_params.gamma = 1.0f / gamma_exp;

                    // Offset
                    ImGui::SetNextItemWidthForSpace(required_width);
                    ImGui::SliderFloat("Offset", &stripe_tone_params.brightness, -1.0f, 1.0f, "%.3f");
                }

                if (!platform()->is_mobile())
                {
                    // Histogram
                    ImGui::Spacing();
                    if (!show_stripe_hist) {
                        if (ImGui::Button("Show Histogram"))
                            show_stripe_hist = true;
                    }
                    else if (ImGui::Button("Hide Histogram"))
                        show_stripe_hist = false;

                    if (show_stripe_hist)
                    {
                        ///ImGui::Text("Mean STRIPE: %.6f", stats.field_info.raw_mean_stripe);
                        float w = ImGui::GetContentRegionAvail().x;
                        if (ImPlot::BeginPlot("Stripe", ImVec2(w, 0)))
                        {
                            stripe_hist_visible = true;

                            ImPlot::SetupAxis(ImAxis_X1, "Value", ImPlotAxisFlags_AutoFit);
                            ImPlot::SetupAxis(ImAxis_Y1, "Pixels", ImPlotAxisFlags_AutoFit);
                            ImPlot::PlotHistogram("##stripe_hist",
                                stats.stripe_histogram.data(),
                                (int)stats.stripe_histogram.size(),
                                100, 1.0
                            );

                            //ImPlot::SetNextLineStyle(ImVec4(1, 0, 0, 1), 2.0f);
                            //ImPlot::PlotInfLines("#min_stripe_v", &stats.field_info.final_min_stripe, 1);
                            //ImPlot::SetNextLineStyle(ImVec4(1, 0, 0, 1), 2.0f);
                            //ImPlot::PlotInfLines("#max_stripe_v", &stats.field_info.final_max_stripe, 1);
                            ImPlot::EndPlot();
                        }
                    }
                }

                ImGui::EndTabBox();
            }

            ImGui::EndTabBar();
        }

        ImGui::EndCollapsingHeaderBox();
    }
}
void Mandelbrot_Scene::UI::populateShaderEditor()
{
    if (platform()->is_mobile())
        return; // shader editing on mobile is not really feasible. 
                // todo: Provide template shaders and hook up shader inputs to UI controls (previous section)

    if (ImGui::CollapsingHeaderBox("Shader", false))
    {
        ImGuiContext& g = *ImGui::GetCurrentContext();
        auto cpos = editor.GetCursorPosition();
        g.PlatformImeData.WantVisible = true;
        g.PlatformImeData.InputPos = ImVec2(0, 0);// ImVec2(cpos.x - 1.0f, cpos.y - g.FontSize);
        g.PlatformImeData.InputLineHeight = g.FontSize;
        g.PlatformImeViewport = ImGui::GetCurrentWindow()->Viewport->ID;

        // update editor text on request from worker
        {
            bl_scoped(update_editor_shader_source);
            if (update_editor_shader_source)
            {
                bl_scoped(shader_source_txt);
                editor.SetText(shader_source_txt);

                combined_frag_log = "";

                update_editor_shader_source = false;
            }
        }

        ImGui::TextUnformatted("GLSL fragment shader:");

        ImVec2 avail = ImGui::GetContentRegionAvail();

        const float splitterThickness = 6.0f;
        const float minEditorHeight = 120.0f;

        float maxEditorHeight = 800.0f;
        if (maxEditorHeight < minEditorHeight)
            maxEditorHeight = minEditorHeight;

        // clamp
        if (editor_height < minEditorHeight) editor_height = minEditorHeight;
        if (editor_height > maxEditorHeight) editor_height = maxEditorHeight;

        // editor takes a controlled height
        ImGui::PushFont(main_window()->monoFont());
        editor.Render("##editor", ImVec2(avail.x, editor_height), true);
        ImGui::PopFont();

        if (editor.IsTextChanged())
        {
            MultiPassShader test_shader;

            if (auto_apply_shader)
            {
                bl_scoped(shader_source_txt);
                shader_source_txt = editor.GetText();
                test_shader.updateFragmentSource(shader_source_txt);
            }
            else
            {
                test_shader.updateFragmentSource(editor.GetText());
            }

            test_shader.testCompile(
                combined_vert_errors,
                combined_frag_errors,
                combined_vert_log,
                combined_frag_log
            );
        }

        // draggable splitter
        ImGui::InvisibleButton("##editor_splitter", ImVec2(avail.x, splitterThickness));
        if (ImGui::IsItemHovered() || ImGui::IsItemActive())
            ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeNS);

        if (ImGui::IsItemActive())
            editor_height += ImGui::GetIO().MouseDelta.y;

        // visible separator line
        {
            ImDrawList* dl = ImGui::GetWindowDrawList();
            ImVec2 a = ImGui::GetItemRectMin();
            ImVec2 b = ImGui::GetItemRectMax();
            float y = (a.y + b.y) * 0.5f;
            dl->AddLine(ImVec2(a.x, y), ImVec2(b.x, y), ImGui::GetColorU32(ImGuiCol_Separator), 1.0f);
        }

        if (!auto_apply_shader)
        {
            if (ImGui::Button("Apply"))
            {
                bl_scoped(shader_source_txt);
                shader_source_txt = editor.GetText();
            }
            ImGui::SameLine();
            if (ImGui::Button("Revert"))
            {
                bl_scoped(shader_source_txt);
                editor.SetText(shader_source_txt);
            }
            ImGui::SameLine();
        }
        if (ImGui::Checkbox("Auto Apply", &auto_apply_shader))
        {
            if (auto_apply_shader)
            {
                bl_scoped(shader_source_txt);
                shader_source_txt = editor.GetText();
            }
        }

        // bottom panel uses the remaining space
        ImGui::BeginChild("##bottom_panel", scale_size(0, 100), true);
        ImGui::InputTextMultiline("##error_log", &combined_frag_log, ImVec2(-FLT_MIN, 0), ImGuiInputTextFlags_ReadOnly);
        ImGui::EndChild();
        //ImGui::TextUnformatted(combined_frag_log.c_str());

        ImGui::EndCollapsingHeaderBox();
    }
}
void Mandelbrot_Scene::UI::populateGradientOptions()
{
    if (ImGui::CollapsingHeaderBox("Gradient", false))
    {
        ImGui::BeginLabelledBox("Base Gradient");
        {
            bl_scoped(gradient);

            ImGui::Text("Load Preset");
            static int selecting_template = -1;
            if (ImGui::Combo("###ColorTemplate", &selecting_template, ColorGradientNames, (int)GradientPreset::COUNT))
            {
                generateGradientFromPreset(gradient, (GradientPreset)selecting_template);
                selecting_template = -1;
            }


            ImGui::Dummy(scale_size(0, 8));

            if (ImGui::GradientEditor(&gradient,
                platform()->dpr(),
                platform()->dpr() * platform()->thumbScale(),
                scale_size(250.0f)))
            {
                // Shift (not needed?)
                ///bl_pull(gradient_shifted, hue_shift, gradient_shift);
                ///transformGradient(gradient_shifted, gradient, (float)gradient_shift, (float)hue_shift);
            }

            

            #if MANDEL_DEV_MODE
            if (ImGui::Button("Copy gradient C++ code"))
            {
                ImGui::SetClipboardText(gradient.to_cpp_marks().c_str());
            }
            #endif

            ImGui::Spacing();
        }
        ImGui::EndLabelledBox();

        ImGui::BeginLabelledBox("Transform");
        {
            bl_scoped(hue_shift, gradient_shift);

            // Shift
            {
                bl_view(gradient);
                bl_pull(gradient_shifted); // no need to push (calculated already in viewportProcess)

                ImGui::TextUnformatted("Shift:");
                float required_space = 0.0f;
                ImGui::IncreaseRequiredSpaceForLabel(required_space, "Gradient    ");

                ImGui::SetNextItemWidthForSpace(required_space);
                static double initial_gradient_shift = 0.0;
                if (ImGui::RevertableDragDouble("Gradient", &gradient_shift, &initial_gradient_shift, 0.01, -100.0, 100.0, " %.3f", ImGuiSliderFlags_AlwaysClamp))
                {
                    gradient_shift = math::wrap(gradient_shift, 0.0, 1.0);

                    // immediately preview
                    transformGradient(gradient_shifted, gradient, (float)gradient_shift, (float)hue_shift);
                }

                
                double hue_shift_rad = math::toRadians(hue_shift);
                static double initial_hue_shift_rad = 0.0;
                ImGui::SetNextItemWidthForSpace(required_space);
                if (ImGui::RevertableSliderAngle("Hue", &hue_shift_rad, &initial_hue_shift_rad, 0.0, math::tau, 3))
                {
                    hue_shift = math::toDegrees(hue_shift_rad);
                    transformGradient(gradient_shifted, gradient, (float)gradient_shift, (float)hue_shift);
                }

                

                ImGui::Spacing();
                ImGui::Text("Live preview");
                ImGui::GradientButton(&gradient_shifted, platform()->dpr());
            }

            if (ImGui::Button("Set as base gradient"))
            {
                bl_scoped(gradient);
                bl_pull_temp(gradient_shifted);
                gradient = gradient_shifted;

                gradient_shift = 0;
                hue_shift = 0;
            }

            ImGui::SameLine();
            if (ImGui::Button("Flatten hue"))
            {
                bl_scoped(gradient);
                auto& marks = gradient.getMarks();
                std::vector<Color> colors;

                for (float x = 0.0f; x <= 1.0f; x += 0.01f) {
                    float color[4];
                    gradient.getColorAt(x, color);
                    colors.push_back(Color(color));
                }
                float avg_hue = Color::avgHueEstimate(colors);
                for (auto& mark : marks)
                {
                    auto adjusted = Color(mark.color).setHue(avg_hue).vec4();
                    memcpy(mark.color, adjusted.data(), sizeof(float) * 4);
                }
                gradient.refreshCache();
            }

            ImGui::SameLine();
            if (ImGui::Button("Invert"))
            {
                bl_scoped(gradient);
                auto& marks = gradient.getMarks();

                for (auto& mark : marks)
                {
                    auto c = Color(mark.color);
                    c.r = 255 - c.r;
                    c.g = 255 - c.g;
                    c.b = 255 - c.b;
                    memcpy(mark.color, c.vec4().data(), sizeof(float) * 4);
                }

                gradient.refreshCache();
            }

            ImGui::SameLine();
            if (ImGui::Button("Reverse"))
            {
                bl_scoped(gradient);
                auto& marks = gradient.getMarks();
 
                for (auto& mark : marks)
                    mark.position = 1.0f - mark.position;

                gradient.refreshCache();
            }

            
        }
        ImGui::EndLabelledBox();

        ImGui::EndCollapsingHeaderBox();
    }
}
void Mandelbrot_Scene::UI::populateAnimation()
{
    if (ImGui::CollapsingHeaderBox("Animate", false))
    {
        bl_scoped(gradient_shift_step, hue_shift_step);
        bl_scoped(animate_gradient_shift, animate_gradient_hue);

        ImGui::TextUnformatted("Gradient:");

        float required_space = 0.0f;
        ImGui::IncreaseRequiredSpaceForLabel(required_space, "Gradient    ");

        // Gradient Shift
        ImGui::Checkbox("##gradient_shift_ani", &animate_gradient_shift);
        ImGui::SameLine();

        if (!animate_gradient_shift) ImGui::BeginDisabled();
        ImGui::PushID("gradient_increment");
        ImGui::SetNextItemWidthForSpace(required_space);
        ImGui::SliderDouble("Shift", &gradient_shift_step, -0.02, 0.02, "%.4f");
        ImGui::PopID();
        if (!animate_gradient_shift) ImGui::EndDisabled();

        // Gradient Hue
        ImGui::Checkbox("##gradient_hue_ani", &animate_gradient_hue);
        ImGui::SameLine();

        if (!animate_gradient_hue) ImGui::BeginDisabled();
        ImGui::PushID("hue_increment");
        ImGui::SetNextItemWidthForSpace(required_space);
        ImGui::SliderDouble("Hue", &hue_shift_step, -5.0, 5.0, "%.3f", ImGuiSliderFlags_AlwaysClamp);
        ImGui::PopID();
        if (!animate_gradient_hue) ImGui::EndDisabled();

        // Stripe
        ImGui::TextUnformatted("Stripe:");

        bl_scoped(animate_stripe_phase);
        bl_scoped(phase_step);

        ImGui::Checkbox("##phase_ani", &animate_stripe_phase);
        ImGui::SameLine();

        if (!animate_stripe_phase) ImGui::BeginDisabled();
        ImGui::PushID("phase_increment");
        ImGui::SetNextItemWidthForSpace(required_space);
        ImGui::SliderAngle("Phase", &phase_step, -math::pi_f / 100.0f, math::pi_f / 100.0f, 1, ImGuiSliderFlags_AlwaysClamp);
        ImGui::PopID();
        if (!animate_stripe_phase) ImGui::EndDisabled();
  
        ImGui::EndCollapsingHeaderBox();
    }
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
    using Complex = Complex<T>;
    constexpr T escape_r = T(16.0);

    Complex z{ 0.0, 0.0 };
    Complex c{ T(x0), T(y0) };

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

template<typename T>
void computeOrbitVel(f128 x0, f128 y0, int iter_lim, std::vector<f64>& xs, std::vector<f64>& ys)
{
    using Complex = Complex<T>;
    constexpr T escape_r = T(16.0);

    Complex z{ 0.0, 0.0 };
    Complex c{ T(x0), T(y0) };
    //Complex z_vel{ 0.0, 0.0 };

    for (int i = 0; i < iter_lim; i++)
    {
        //mandel_step_inertial(z, z_vel, c);
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
        bl_view(stats, iter_lim, float_type);

        ImGui::SliderDouble("Padding", &orbit_padding, 0.0, 1.0, "%.2f", ImGuiSliderFlags_AlwaysClamp);

        std::vector<f64> xs, ys;
        xs.push_back(0.0);
        ys.push_back(0.0);

        auto [x0, y0] = stats.hovered_field_world_pos;

        table_invoke(dispatch_table(computeOrbitVel, x0, y0, iter_lim, xs, ys), float_type);

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

void Mandelbrot_Scene::UI::populateCaptureOptions()
{
    if (platform()->is_mobile())
        return;

    if (ImGui::CollapsingHeaderBox("Capture Options", false))
    {
        // ---------------- Permitted capture presets during batch snapshot ----------------
        ImGui::BeginLabelledBox("Allowed presets (for current state)");
        auto& standard_presets = main_window()->getSnapshotPresetManager()->allPresets();

        bl_scoped(valid_presets);
        populateCapturePresetsList<CapturePresetsSelectMode::MULTI>([&](int i) -> CapturePreset& {
            return standard_presets[i];
        }, (int)standard_presets.size(), &valid_presets, selected_preset_i);

        ImGui::EndLabelledBox();

        // ---------------- Batch snapshot examples ----------------
        ImGui::Spacing();
        ImGui::BeginLabelledBox("Batch snapshot");
        ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + std::max(scale_size(150.0f), ImGui::GetContentRegionAvail().x));
        ImGui::TextWrapped(
            "Renders all examples with the global settings "
            "filtered by the allowed per-state preset list above."
        );

        /// TODO: Have an option: "skip if folder already contains snapshot with matching XMP state"
        /// - works with non-deterministic filenames (snap1, snap2, etc)
        /// - but harder to do than matching hash-based filenames

        bl_scoped(rendering_examples, rendering_example_i, render_batch_name, ignore_preset_filters);
        if (render_batch_name.empty()) ImGui::BeginDisabled();
        if (ImGui::Button("Render ALL bookmarks"))
        {
            rendering_examples = true;
            rendering_example_i = 0;
        }
        if (render_batch_name.empty()) ImGui::EndDisabled();

        float required_space = 0.0f;
        ImGui::IncreaseRequiredSpaceForLabel(required_space, "Folder name");
        ImGui::SetNextItemWidthForSpace(required_space);

        ImGui::SetNextItemWidthForSpace(required_space);
        ImGui::InputText("Folder name", &render_batch_name);
        ImGui::Checkbox("Ignore filters", &ignore_preset_filters);

        ImGui::EndLabelledBox();

        // ---------------- Steady Zoom Animation ----------------
        ImGui::Spacing();
        ImGui::BeginLabelledBox("Steady Zoom Animation");
        ImGui::TextWrapped("Begins a steady zoom animation from zoom=1 to the current view");
        
        bl_scoped(record_steady_zoom, steady_zoom_mult_speed, steady_zoom_rotate_speed);
        ImGui::Checkbox("Record", &record_steady_zoom);

        // begin animation button
        if (ImGui::Button("Begin"))
            bl_schedule([](Mandelbrot_Scene& s) { s.beginSteadyZoom(); });
        
        // end animation button
        bool contents_disabled = ImGui::CollapsingHeaderContentsDisabled();
        if (contents_disabled) ImGui::EndDisabled(); else ImGui::BeginDisabled();
        ImGui::SameLine();
        if (ImGui::Button("Stop"))
            bl_schedule([](Mandelbrot_Scene& s) { s.endSteadyZoom(); });
        
        if (contents_disabled) ImGui::BeginDisabled(); else ImGui::EndDisabled();



        ImGui::SetNextItemWidthForSpace(required_space);
        ImGui::SliderDouble("Zoom rate", &steady_zoom_mult_speed, 0.01, 0.1,
            "%.2f", ImGuiSliderFlags_AlwaysClamp);

        constexpr double spin_mag = math::toRadians(0.1);
        ImGui::SetNextItemWidthForSpace(required_space);
        ImGui::SliderAngle("Spin rate", &steady_zoom_rotate_speed, -spin_mag, spin_mag,
            2, ImGuiSliderFlags_AlwaysClamp);

        ImGui::EndLabelledBox();

        ImGui::EndCollapsingHeaderBox();
    }
}

SIM_END;