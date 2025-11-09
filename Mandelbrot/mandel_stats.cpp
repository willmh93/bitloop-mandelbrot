#include "Mandelbrot.h"

SIM_BEG;

void Mandelbrot_Scene::UI::populateStats()
{
    //return;

    if (ImGui::Section("Statistics", true))
    {
        bl_scoped(stats);

        // --- compute info ---
        {
            ImGui::GroupBox box("compute_info", "Compute info", scale_size(13.0f), scale_size(20.0f));

            bl_pull(mandel_features);
            bl_pull(camera);
            bl_pull(dt_avg);
            bl_pull(computing_phase);

            // float precision
            const char* float_precision = FloatingPointTypeNames[(int)getRequiredFloatType((MandelKernelFeatures)mandel_features, camera.relativeZoom<f128>())];
            ImGui::Text("Active float precision:  %s", float_precision);

            // current phase/row
            ImGui::Text("Current phase:  %d/2", computing_phase);

            // active features
            {
                std::stringstream ss;
                for (int i = 0, j = 0; i <= 4; i++) {
                    if ((bool)(mandel_features & (1 << i))) {
                        if (j++ != 0) ss << " | ";
                        ss << MandelSmoothingNames[i + 1];
                    }
                }
                ImGui::Text("Computing Features:  %s", ss.str().c_str());
            }

            // compute timer
            ImGui::Text("Full compute time:  %.0f ms", dt_avg);

            ImGui::Spacing();
            ImGui::Spacing();

            // mouse pixel info
            {
                ImGui::Text("Mouse world pos:   [%.2f, %.2f]", stats.hovered_field_world_pos.x, stats.hovered_field_world_pos.y);

                if (stats.hovered_field_pixel.depth < INSIDE_MANDELBROT_SET_SKIPPED)
                {
                    ImGui::Text("Mouse depth:        %.1f iterations", stats.hovered_field_pixel.depth);
                    ImGui::Text("Mouse depth (norm): %.4f iterations", stats.hovered_field_pixel.final_depth);
                }
                else
                {
                    ImGui::Text("Mouse depth:        <interior>");
                    ImGui::Text("Mouse depth (norm): <interior>");
                }
            }
        }

        {
            bl_pull(camera);

            ImGui::GroupBox box("cam_info", "Camera info", scale_size(13.0f), scale_size(20.0f));
            ImGui::Text("Height: %.4f", (double)toHeight(camera.relativeZoom<f128>()));
            ImGui::Text("Normalized Zoom: %.4f", (double)toNormalizedZoom(camera.relativeZoom<f128>()));
        }

        // --- capture info ---
        {
            bl_pull(steady_zoom_pct);

            ImGui::GroupBox box("capture_info", "Capture info", scale_size(13.0f), scale_size(20.0f));
            ImGui::Text("Capture progress: %.1f%%", steady_zoom_pct);
            
        }

        // --- plots ---
        ImGui::Spacing();
        ImGui::Spacing();

        // depth histogram
        ///{
        ///
        ///    if (ImPlot::BeginPlot("Depth Histogram"))
        ///    {
        ///        depth_xs.reserve(stats.depth_histogram.size());
        ///        depth_ys.reserve(stats.depth_histogram.size());
        ///        depth_xs.clear();
        ///        depth_ys.clear();
        ///
        ///        double bucket_size = stats.depth_histogram_bucket_size;
        ///        double min_bucket_depth = (int)(stats.field_info.min_depth / bucket_size) * bucket_size;
        ///        for (int i=0; i < stats.depth_histogram.size(); i++)
        ///        {
        ///            //if (v <= 1) continue; // Hide very small counts
        ///            int d = (int)(min_bucket_depth + (i * stats.depth_histogram_bucket_size));
        ///            depth_xs.push_back(d);
        ///            depth_ys.push_back(stats.depth_histogram[i]);
        ///        }
        ///
        ///        ImPlot::SetupAxis(ImAxis_X1, "depth", ImPlotAxisFlags_AutoFit);
        ///        ImPlot::SetupAxis(ImAxis_Y1, "pixel count", ImPlotAxisFlags_AutoFit);
        ///        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
        ///        if (stats.depth_histogram.size())
        ///        {
        ///            ImPlot::PlotLine("##depth_hist", depth_xs.data(), depth_ys.data(), (int)depth_xs.size());
        ///
        ///            /*double max_depth = (--stats.depth_histogram.end())->first;
        ///            if (stats.hovered_field_pixel.depth < max_depth)
        ///            {
        ///
        ///                auto nearest_it = stats.depth_histogram.lower_bound((int)stats.hovered_field_pixel.depth);
        ///                if (nearest_it != stats.depth_histogram.end())
        ///                {
        ///                    double px = nearest_it->first;
        ///                    double py = nearest_it->second;
        ///                    ImPlot::PlotInfLines("##depth_mark", &px, 1); // vertical line at X = mark
        ///                    ImPlot::PlotScatter("##depth_point", &px, &py, 1);
        ///                }
        ///            }*/
        ///        }
        ///        ImPlot::EndPlot();
        ///    }
        ///}

        // zoom - max_depth graph
        /*{
            if (ImPlot::BeginPlot("Max depth"))
            {
                max_depth_xs.clear();
                max_depth_ys.clear();
                int i = 1;
                for (double z = 1; z < 1e30; z *= 10)
                {
                    double max_depth = mandelbrotIterLimit(z);
                    max_depth_xs.push_back(i);
                    max_depth_ys.push_back(max_depth);
                    i++;
                }
                ImPlot::SetupAxis(ImAxis_X1, "zoom e10", ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxis(ImAxis_Y1, "max iter", ImPlotAxisFlags_AutoFit);
                //ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
                //ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                ImPlot::PlotLine("##maxdepth_graph", max_depth_xs.data(), max_depth_ys.data(), (int)max_depth_xs.size());
                ImPlot::EndPlot();
            }
        }*/

        /* // todo: crashes since perturbation
		{
            if (ImPlot::BeginPlot("Stripe Angle"))
            {
                ImPlot::PlotHistogram(
                    "##stripe_hist",
                    stats.stripe_histogram.data(),
                    (int)stats.stripe_histogram.size(),
                    1000, 1.0, ImPlotRange(0.0, 1.0),
                    //ImPlotHistogramFlags_NoOutliers
                );
                ImPlot::EndPlot();

                //stripe_xs.clear();
                ////max_depth_ys.clear();
                //int i = 0;
                //for (int i = 0; i<360; i++)
                //{
                //    stripe_xs.push_back(i);
                //}
                //
                //ImPlot::SetupAxis(ImAxis_X1, "angle", ImPlotAxisFlags_AutoFit);
                //ImPlot::SetupAxis(ImAxis_Y1, "px count", ImPlotAxisFlags_AutoFit);
                ////ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
                ////ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                //ImPlot::PlotLine("##stripe_graph", stripe_xs.data(), stats.stripe_histogram, 360);
                //ImPlot::EndPlot();
            }
        }*/
    }
}

void Mandelbrot_Scene::collectStats(bool renormalized)
{
    // stats for mouse position
    if (active_field && active_bmp)
    {
        if (mouse->stage_x >= 0 && mouse->stage_y >= 0 && mouse->stage_x < active_bmp->width() && mouse->stage_y < active_bmp->height())
        {
            IVec2 pos = active_bmp->pixelPosFromWorld(DDVec2(mouse->world_x, mouse->world_y));
            EscapeFieldPixel* p = active_field->get(pos.x, pos.y);
            if (p) stats.hovered_field_pixel = *p;

            active_bmp->worldPos((int)mouse->client_x, (int)mouse->client_y, stats.hovered_field_world_pos.x, stats.hovered_field_world_pos.y);

        }
    }

    UNUSED(renormalized);

    /*if (renormalized && final_frame_complete)
    { 
        if (!this->active_field) return;
        EscapeField& field = *this->active_field;

        memset(&stats.dirty, 0, sizeof(stats.dirty));

        stats.field_info.min_depth = field.min_depth;
        stats.field_info.min_depth = field.max_depth;

        double trimmed_min_depth = field.min_depth, trimmed_max_depth = field.max_depth;

        if (trimmed_max_depth > trimmed_min_depth)
        {
            // Depth Histogram
            for (int k = 0; k < 2; k++)
            {
                double bucket_size = (trimmed_max_depth - trimmed_min_depth) / 500.0;
                stats.depth_histogram_bucket_size = bucket_size;

                if (bucket_size > std::numeric_limits<double>::epsilon())
                {
                    auto& hist = stats.depth_histogram;
                    hist.clear();
                    hist.resize(501);
                    for (size_t i = 0; i < field.size(); i++)
                    {
                        EscapeFieldPixel& p = field[i];
                        if (p.depth > 1e10 || isnan(p.depth)) continue;

                        int bucket_depth = (int)((p.depth - trimmed_min_depth) / bucket_size);// *bucket_size;
                        hist[bucket_depth]++;
                    }

                    // Trim small isolated entries from back of entries
                    //if (hist.size() >= 2)
                    //{
                    //    int min_depth = hist.front();
                    //    int max_depth = hist.back();
                    //
                    //    i64 sum_count = 0;
                    //    for (auto c : hist) sum_count += c;
                    //    int mean_px_count = (int)((double)sum_count / (double)hist.size());
                    //    int min_px_count = mean_px_count / 50;
                    //    auto valid_entry = [&](int d) { return hist[d] >= min_px_count; };
                    //    
                    //    // Find 'x' consecutive valid entries from the back and stop
                    //    const int consecutive_entries = 5;
                    //    //int trim_from_depth;
                    //    int i;
                    //    //for (trim_from_depth = max_depth; trim_from_depth >= min_depth; trim_from_depth -= bucket_size)
                    //    for (i = 499; i >= 0; --i)
                    //    {
                    //        //int d = (int)(trim_from_depth / bucket_size);// *bucket_size;
                    //        bool gaps = false;
                    //        for (int j = 0; j < consecutive_entries; j++)
                    //        {
                    //            if (!valid_entry(d -  j)) {
                    //                gaps = true;
                    //                break;
                    //            }
                    //        }
                    //        if (!gaps) break;
                    //    }
                    //
                    //    auto it = hist.upper_bound(trim_from_depth);
                    //    hist.erase(it, hist.end());
                    //    max_depth = trim_from_depth;
                    //
                    //    for (int d = min_depth + bucket_size; d < max_depth; d += bucket_size)
                    //    {
                    //        int bucket_depth = (int)(d / bucket_size) * bucket_size;
                    //        if (hist.count(bucket_depth) == 0)
                    //            hist[bucket_depth] = 0;
                    //    }
                    //
                    //    trimmed_max_depth = max_depth;
                    //}
                    stats.dirty.depth_histogram = true;
                }
            }
        }

        // Stripe histogram
        stats.stripe_histogram.clear();
        //memset(&stats.stripe_histogram, 0, sizeof(stats.stripe_histogram));
        for (size_t i = 0; i < field.size(); i++)
        {
            EscapeFieldPixel& p = field[i];
            if (p.depth > 1e10 || isnan(p.depth)) continue;
            stats.stripe_histogram.push_back(p.final_stripe);
            //int angle_index = (int)Math::toDegrees(p.stripe_32);
            //stats.stripe_histogram[angle_index]++;
        }
    }*/
}

SIM_END;