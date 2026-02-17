#include "../Mandelbrot.h"
#include <format>

SIM_BEG;

void SetupTable2(float col1_w)
{
    ImGui::TableSetupColumn("Label", ImGuiTableColumnFlags_WidthFixed, col1_w);
    ImGui::TableSetupColumn("Value", ImGuiTableColumnFlags_WidthStretch);
}

static void TableRowText2(const char* label, const char* fmt, ...)
{
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::TextUnformatted(label);

    ImGui::TableSetColumnIndex(1);
    va_list args;
    va_start(args, fmt);
    ImGui::TextV(fmt, args);
    va_end(args);
}

static void TableRowSpacing2()
{
    ImGui::TableNextRow();
    ImGui::TableNextColumn(); ImGui::Dummy(ImVec2(0, ImGui::GetStyle().ItemSpacing.y));
    ImGui::TableNextColumn(); ImGui::Dummy(ImVec2(0, ImGui::GetStyle().ItemSpacing.y));
}

void Mandelbrot_Scene::UI::populateStats()
{
    if (ImGui::CollapsingHeaderBox("Statistics", false))
    {
        ImVec2 cellPad = ImGui::GetStyle().CellPadding;
        int flags = 
            ImGuiTableFlags_SizingStretchProp |
            ImGuiTableFlags_PadOuterX |
            ImGuiTableFlags_NoSavedSettings;

        ImGui::PushStyleVar(ImGuiStyleVar_CellPadding, ImVec2(cellPad.x + 8.0f, cellPad.y));

        #if MANDEL_EXTENDED_FIELD_STATS
        const char* col1_longest = "Total dominance rebases";
        #else
        const char* col1_longest = "Kernel Features";
        #endif

        const float col1_w = ImGui::CalcTextSize(col1_longest).x + ImGui::GetStyle().CellPadding.x;

        bl_pull(stats, camera, deduced_kernel_mode, float_type, mandel_features, steady_zoom_pct, dt_last, computing_phase);

        int decimals = camera.getPositionDecimals();

        // --- compute info ---
        {
            ImGui::BeginLabelledBox("Compute");
            
            if (ImGui::BeginTable("##compute_stats", 2, flags, ImVec2(-FLT_MIN, 0.0f)))
            {
                SetupTable2(col1_w);

                // float precision
                const char* float_precision = FloatingPointTypeNames[(int)float_type];
                TableRowText2("Float precision", float_precision);
                TableRowText2("Current phase", "%d/%d", computing_phase, PHASE_COUNT-1);
                TableRowText2("Kernel", KernelModeNames[(int)deduced_kernel_mode]);

                #if MANDEL_EXTENDED_FIELD_STATS
                ImGui::Spacing();
                TableRowText2("Total slack rebases", "%s", text::format_human_u64(stats.field_info.extended.slack_rebases).c_str());
                TableRowText2("Total dominance rebases",  "%s", text::format_human_u64(stats.field_info.extended.dominance_rebases).c_str());
                TableRowText2("Total iterations", "%s", text::format_human_u64(stats.field_info.extended.iter_count).c_str());
                TableRowText2("Total escapes", "%s", text::format_human_u64(stats.field_info.extended.escape_count).c_str());
                #endif

                // active features
                {
                    std::stringstream ss;
                    for (int i = 0, j = 0; i <= 4; i++) {
                        if ((bool)(mandel_features & (1 << i))) {
                            if (j++ != 0) ss << " | ";
                            ss << MandelFeatureNames[i + 1];
                        }
                    }
                    TableRowText2("Kernel Features", "%s", ss.str().c_str());
                }

                #if MANDEL_DEV_MODE
                bl_pull(compute_timeout);
                TableRowText2("Compute timeout", "%d ms", compute_timeout);
                #endif

                // compute timer
                TableRowText2("Compute timer", "%.0f ms", dt_last);

                #if MANDEL_DEV_MODE
                TableRowSpacing2();
                TableRowText2("Phase timers:", "");

                bl_pull(phase_timers, starting_phase, final_frame_complete);
                for (int phase = 0; phase < PHASE_COUNT; phase++)
                {
                    std::string phase_txt = "Phase: ";
                    phase_txt += std::to_string(phase);

                    if (phase > starting_phase && (phase < computing_phase || final_frame_complete))
                    {
                        double mult = (phase_timers[phase] / phase_timers[phase-1]);
                        TableRowText2(phase_txt.c_str(), "%.2f ms (%.1fx)", phase_timers[phase], mult);
                    }
                    else
                    {
                        if (phase >= starting_phase)
                            TableRowText2(phase_txt.c_str(), "%.2f ms", phase_timers[phase]);
                        else
                            TableRowText2(phase_txt.c_str(), "-");
                    }
                }

                bl_pull(phase_elapsted_estimated_final, phase_elapsed_mult_results, phase_elapsed_mult_sma_result);

                TableRowSpacing2();
                TableRowText2("Average mults:", "");
                for (int phase = 0; phase < PHASE_COUNT; phase++)
                {
                    std::string phase_txt = "Phase: ";
                    phase_txt += std::to_string(phase);

                    TableRowText2(phase_txt.c_str(), "%.1fx", phase_elapsed_mult_results[phase]);
                }

                TableRowSpacing2();
                TableRowText2("Avg mult", "%.1fx", phase_elapsed_mult_sma_result);
                TableRowText2("Expected final phase", "%.4f", phase_elapsted_estimated_final);
                #endif
                
                ImGui::EndTable();
            }

            ImGui::EndLabelledBox();
        }

        if (!platform()->is_mobile())
        {
            // --- normalize info ---
            ImGui::BeginLabelledBox("Normalization Field");
            if (ImGui::BeginTable("##field_stats", 2, flags, ImVec2(-FLT_MIN, 0.0f)))
            {
                SetupTable2(col1_w);

                TableRowText2("Raw:", "");
                TableRowText2("    ITER min", std::to_string(stats.field_info.raw_min_depth).c_str());
                TableRowText2("    ITER max", std::to_string(stats.field_info.raw_max_depth).c_str());

                TableRowSpacing2();
                TableRowText2("    DIST min", std::to_string(stats.field_info.raw_min_dist).c_str());
                TableRowText2("    DIST max", std::to_string(stats.field_info.raw_max_dist).c_str());

                TableRowSpacing2();
                TableRowText2("    STRIPE min", std::to_string(stats.field_info.raw_min_stripe).c_str());
                TableRowText2("    STRIPE max", std::to_string(stats.field_info.raw_max_stripe).c_str());

                TableRowSpacing2();
                TableRowText2("    STRIPE mean",       std::to_string(stats.field_info.raw_mean_stripe).c_str());
                TableRowText2("    STRIPE magnitude",  std::to_string(stats.field_info.raw_mag_stripe).c_str());

                //

                TableRowSpacing2();
                TableRowText2("Final:", "");
                TableRowText2("    ITER min", std::to_string(stats.field_info.final_min_depth).c_str());
                TableRowText2("    ITER max", std::to_string(stats.field_info.final_max_depth).c_str());

                TableRowSpacing2();
                TableRowText2("    DIST min", std::to_string(stats.field_info.final_min_dist).c_str());
                TableRowText2("    DIST max", std::to_string(stats.field_info.final_max_dist).c_str());

                TableRowSpacing2();
                TableRowText2("    STRIPE min", std::to_string(stats.field_info.final_min_stripe).c_str());
                TableRowText2("    STRIPE max", std::to_string(stats.field_info.final_max_stripe).c_str());

                ImGui::EndTable();
            }
            ImGui::EndLabelledBox();

            // --- mouse info ---
            ImGui::BeginLabelledBox("Mouse");
            if (ImGui::BeginTable("##hover_stats", 2, flags, ImVec2(-FLT_MIN, 0.0f)))
            {
                bl_pull(stripe_params);

                auto px = stats.hovered_field_pixel;

                float phase = stripe_params.phase;
                const f64 LOG_ER2 = table_invoke(dispatch_table(log_escape_radius2), mandel_features);
                const f32 cphi    = std::cos(phase);
                const f32 sphi    = std::sin(phase);
                const f32 stripe  = stripeFromAccum(px.stripe, LOG_ER2, cphi, sphi);

                bool escaped = (px.depth < INSIDE_MANDELBROT_SET_SKIPPED);
                std::string x_str = to_string(stats.hovered_field_world_pos.x, decimals, true);
                std::string y_str = to_string(stats.hovered_field_world_pos.y, decimals, true);
                std::string dist_str = to_string(px.dist, decimals, true);

                SetupTable2(col1_w);

                TableRowText2("Coordinate:", "");
                TableRowText2("    Real", x_str.c_str());
                TableRowText2("    Imaginary", y_str.c_str());

                TableRowSpacing2();
                TableRowText2("Raw:", "");
                TableRowText2("    DEPTH",  (escaped ? std::format("{:.2f} iterations", px.depth) : "<interior>").c_str());
                TableRowText2("    DIST",   (escaped ? dist_str : "<interior>").c_str());
                TableRowText2("    STRIPE", (escaped ? std::format("{:.4f}", stripe) : "<interior>").c_str());

                TableRowSpacing2();
                TableRowText2("Final:", "");
                TableRowText2("    DEPTH",  (escaped ? std::format("{:.3f}", px.final_depth) : "<interior>").c_str());
                TableRowText2("    DIST",   (escaped ? std::format("{:.3f}", px.final_dist) : "<interior>").c_str());
                TableRowText2("    STRIPE", (escaped ? std::format("{:.3f}", px.final_stripe) : "<interior>").c_str());

                ImGui::EndTable();
            }
            ImGui::EndLabelledBox();
        }

        // --- camera info ---
        {
            ImGui::BeginLabelledBox("Camera");
            if (ImGui::BeginTable("##camera_stats", 2, flags, ImVec2(-FLT_MIN, 0.0f)))
            {
                SetupTable2(col1_w);
                
                auto format_zoom = [](f128 f) -> std::string { return to_string(f, 5, false, true, true).c_str(); };
                TableRowText2("Zoom",                "%s",   format_zoom(camera.relativeZoom<f128>()).c_str());
                #if MANDEL_DEV_MODE                  
                TableRowText2("Height",              "%.4f", (f64)toHeight(camera.relativeZoom<f128>()));
                TableRowText2("Normalized Zoom",     "%.4f", (f64)toNormalizedZoom(camera.relativeZoom<f128>()));
                #endif

                ImGui::EndTable();
            }
            ImGui::EndLabelledBox();
        }

        // --- capture info ---
        if (!platform()->is_mobile())
        {
            ImGui::BeginLabelledBox("Capture");
            if (ImGui::BeginTable("##capture_stats", 2, flags, ImVec2(-FLT_MIN, 0.0f)))
            {
                SetupTable2(col1_w);

                TableRowText2("Capture progress", "%.1f%%", steady_zoom_pct);

                ImGui::EndTable();
            }
            ImGui::EndLabelledBox();
        }

        ImGui::PopStyleVar();

        ImGui::EndCollapsingHeaderBox();
    }
}

template<bool I, bool D, bool S> void syncHist(EscapeField& field, MandelStats& stats)
{
    if constexpr (!I && !D && !S) return;
    for (size_t i = 0; i < field.size(); i++)
    {
        EscapeFieldPixel& p = field[i];
        if (p.depth > 1e10 || isnan(p.depth)) continue;
        if constexpr (I) stats.iter_histogram.push_back(p.final_depth);
        if constexpr (D) stats.dist_histogram.push_back(p.final_dist);
        if constexpr (S) stats.stripe_histogram.push_back(p.final_stripe);
    }
};

void Mandelbrot_Scene::collectStats(bool renormalized)
{
    // reset stats we need syncing
    stats.dirty = 0;

    // always sync field info with UI
    stats.sync(MandelStats::Dirty::field_info);

    bool finalized_compute = (renormalized && final_frame_complete);
    (void)finalized_compute; // hide warning when no sections use bool

    EscapeField& active_field = activeField();
    WorldRasterGrid128& active_grid = activeRasterGrid();

    // stats for mouse position
    {
        // raw min/max values
        stats.field_info.raw_min_depth     = active_field.raw_min_depth;
        stats.field_info.raw_max_depth     = active_field.raw_max_depth;
        stats.field_info.raw_min_dist      = active_field.raw_min_dist;
        stats.field_info.raw_max_dist      = active_field.raw_max_dist;
        stats.field_info.raw_min_stripe    = active_field.raw_min_stripe;
        stats.field_info.raw_max_stripe    = active_field.raw_max_stripe;
        stats.field_info.raw_mean_stripe   = active_field.raw_mean_stripe;
        stats.field_info.raw_mag_stripe    = active_field.raw_mag_stripe;

        // final min/max values after normalization/toning/cycling
        stats.field_info.final_min_depth   = active_field.final_min_depth;
        stats.field_info.final_max_depth   = active_field.final_max_depth;
        stats.field_info.final_min_dist    = active_field.final_min_dist;
        stats.field_info.final_max_dist    = active_field.final_max_dist;
        stats.field_info.final_min_stripe  = active_field.final_min_stripe;
        stats.field_info.final_max_stripe  = active_field.final_max_stripe;

        if (mouse->stage_x >= 0 &&
            mouse->stage_y >= 0 &&
            mouse->stage_x < active_grid.rasterWidth() &&
            mouse->stage_y < active_grid.rasterHeight())
        {
            stats.sync(MandelStats::Dirty::hovered_field_pixel);

            IVec2 pos = active_grid.pixelPosFromWorld(DDVec2(mouse->world_x, mouse->world_y));
            EscapeFieldPixel* p = active_field.get(pos.x, pos.y);
            if (p) stats.hovered_field_pixel = *p;

            active_grid.pixelWorldPos<f128>((int)mouse->client_x, (int)mouse->client_y, stats.hovered_field_world_pos.x, stats.hovered_field_world_pos.y);
        }
    }

    #if MANDEL_EXTENDED_FIELD_STATS
    if (finalized_compute)
    {
        /// Grab extended info in batches, and merge
        std::vector<ExtendedFieldStats> batches = Thread::forEachBatch(*active_field, [&](std::span<EscapeFieldPixel> batch)
        {
            ExtendedFieldStats s;
            for (EscapeFieldPixel& field_pixel : batch) 
                s += field_pixel.extended_stats;
            return s;
        });

        ExtendedFieldStats merged;
        for (const ExtendedFieldStats& s : batches) 
            merged += s;

        stats.field_info.extended = merged;
    }
    #endif

    if (finalized_compute || 
        Changed(iter_hist_visible) ||
        Changed(dist_hist_visible) ||
        Changed(stripe_hist_visible))
    { 
        if (iter_hist_visible)
        {
            stats.sync(MandelStats::Dirty::iter_histogram);
            stats.iter_histogram.clear();
        }
        
        if (dist_hist_visible)
        {
            stats.sync(MandelStats::Dirty::dist_histogram);
            stats.dist_histogram.clear();
        }

        if (stripe_hist_visible)
        {
            stats.sync(MandelStats::Dirty::stripe_histogram);
            stats.stripe_histogram.clear();
        }
  
        table_invoke(dispatch_table(syncHist, active_field, stats),
            iter_hist_visible,
            dist_hist_visible,
            stripe_hist_visible
        );
    }
}

SIM_END;