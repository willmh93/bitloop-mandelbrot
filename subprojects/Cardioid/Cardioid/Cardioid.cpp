#include "Cardioid.h"
//#include "imgui_spline.h"

SIM_BEG;

// Exported
void plot(const SceneBase *scene, Viewport* ctx, bool interactive, int segments, double ox)
{
    ctx->saveCameraTransform();
    ctx->worldHudMode();

    // Full plot
    ctx->setLineWidth(2);
    ctx->setStrokeStyle(255, 0, 0);
    ctx->beginPath();

    double angle_step = Math::TWO_PI / static_cast<double>(segments);
    bool first = true;

    for (double angle = 0.0; angle < Math::TWO_PI; angle += angle_step)
    {
        double plot_x = 0.5 * cos(angle) - 0.25 * cos(angle * 2.0);
        double plot_y = 0.5 * sin(angle) - 0.25 * sin(angle * 2.0);

        if (first)
        {
            ctx->moveTo(plot_x, plot_y);
            first = false;
        }
        else
            ctx->lineTo(plot_x, plot_y);
    }

    ctx->stroke();

    if (interactive)
    {
        double tx1, ty1, tx2, ty2, ta, oa, d;
        cardioidPolarCoord((f64)scene->mouse->world_x, (f64)scene->mouse->world_y, tx1, ty1, ta, d, oa);

        ctx->scalingSizes(true);

        // Big blue circle
        double r1 = 0.5;
        ctx->setStrokeStyle(255, 0, 255);
        ctx->strokeEllipse(ox, 0.0, r1);


        // p1 dot
        double p1_x = r1 * cos(oa) + ox;
        double p1_y = r1 * sin(oa);

        // Med green circle
        double r2 = r1 / 4.0;
        double angle2 = 2.0 * oa;
        double c2_x = p1_x - r2 * cos(angle2);
        double c2_y = p1_y - r2 * sin(angle2);
        ctx->setStrokeStyle(255, 0, 255);
        ctx->strokeEllipse(c2_x, c2_y, r2);

        ctx->scalingSizes(false);

        DVec2 off = ctx->Offset(5, 5);

        // draw p1 dot
        ctx->setFillStyle(255, 0, 255);
        ctx->fillEllipse(p1_x, p1_y, 8.0);
        ctx->fillText("p1", p1_x + off.x, p1_y + off.y);

        // draw p2 dot
        double p2_x = c2_x - r2 * cos(angle2);
        double p2_y = c2_y - r2 * sin(angle2);
        ctx->setFillStyle(255, 0, 255);
        ctx->fillEllipse(p2_x, p2_y, 8.0);
        ctx->fillText("p2", p2_x + off.x, p2_y + off.y);

        // draw tangent arrow
        double tangent_head_x = p2_x + cos(ta) * 0.2;
        double tangent_head_y = p2_y + sin(ta) * 0.2;
        ctx->drawArrow(DVec2{ p2_x, p2_y }, DVec2{ tangent_head_x, tangent_head_y }, Color(255, 255, 255));
        //ctx->beginPath();
        //ctx->arrowMoveTo(p2_x, p2_y);
        //ctx->arrowDrawTo(tangent_head_x, tangent_head_y);
        //ctx->stroke();

        // Draw perpendicular arrow (to mouse)
        double perp_angle = Math::wrapRadians2PI(ta - Math::HALF_PI);
        tx2 = tx1 + cos(perp_angle) * d;
        ty2 = ty1 + sin(perp_angle) * d;

        //ctx->beginPath();
        //ctx->arrowMoveTo(tx1, ty1);
        //ctx->arrowDrawTo(tx2, ty2);
        //ctx->stroke();
        ctx->drawArrow(DVec2{ tx1, ty1 }, DVec2{ tx2, ty2 }, Color(255, 255, 255));
    }
    ctx->restoreCameraTransform();
}

void animatePlot(Viewport* ctx, double scale, double ox, double orig_angle, double dist)
{
    ctx->scalingSizes(true);
    ctx->scalingLines(false);
    ctx->setLineWidth(2);

    // Big blue circle
    double r1 = scale * 0.5;
    ctx->setStrokeStyle(80, 80, 255);
    ctx->strokeEllipse(ox * scale, 0.0, r1);

    // p1 dot
    double p1_x = r1 * cos(orig_angle) + ox * scale;
    double p1_y = r1 * sin(orig_angle);

    // Med green circle
    double r2 = r1 / 4.0;
    double angle2 = 2.0 * orig_angle;
    double c2_x = p1_x - r2 * cos(angle2);
    double c2_y = p1_y - r2 * sin(angle2);
    ctx->setStrokeStyle(0, 255, 0);
    ctx->strokeEllipse(c2_x, c2_y, r2);

    DVec2 off = ctx->Offset(5, 5);

    //camera->scalingLines(false);

    ctx->scalingSizes(false);

    // draw p1 dot
    ctx->setFillStyle(0, 255, 255);
    ctx->fillEllipse(p1_x, p1_y, 5.0);
    ctx->fillText("p1", p1_x + off.x, p1_y + off.y);

    // draw p2 dot
    double p2_x = c2_x - r2 * cos(angle2);
    double p2_y = c2_y - r2 * sin(angle2);
    ctx->setFillStyle(0, 255, 255);
    ctx->fillEllipse(p2_x, p2_y, 5.0);
    ctx->fillText("p2", p2_x + off.x, p2_y + off.y);

    double tangent_angle = 1.5 * orig_angle;

    // draw tangent arrow
    double tangent_head_x = p2_x + cos(tangent_angle) * 0.2;
    double tangent_head_y = p2_y + sin(tangent_angle) * 0.2;
    ctx->drawArrow(DVec2{ p2_x, p2_y }, DVec2{ tangent_head_x, tangent_head_y }, Color(255, 255, 255));
    //ctx->beginPath();
    //ctx->arrowMoveTo(p2_x, p2_y);
    //ctx->arrowDrawTo(tangent_head_x, tangent_head_y);
    //ctx->stroke();


    // Draw perpendicular arrow (to mouse)
    //double perp_angle = Math::wrapRadians2PI(tangent_angle - M_PI / 2.0);
    ///tx2 = tx1 + cos(perp_angle) * d;
    ///ty2 = ty1 + sin(perp_angle) * d;
    DVec2 tp = fromPolarCoordinate(orig_angle, dist);

    double tx1 = 0.5 * cos(orig_angle) - 0.25 * cos(orig_angle * 2.0);
    double ty1 = 0.5 * sin(orig_angle) - 0.25 * sin(orig_angle * 2.0);
    double tx2 = tp.x;
    double ty2 = tp.y;

    //ctx->beginPath();
    //ctx->arrowMoveTo(tx1, ty1);
    //ctx->arrowDrawTo(tx2, ty2);
    //ctx->stroke();

    ctx->drawArrow(DVec2{ tx1, ty1 }, DVec2{ tx2, ty2 }, Color(255, 255, 255));

    ///ctx->print() << "p1 angle: " << QString::asprintf("%.0f deg", oa * 180.0 / M_PI);
    ///ctx->print() << "\np2 tangent angle: " << QString::asprintf("%.0f deg", ta * 180.0 / M_PI);
    ///ctx->print() << "\np2 perp angle: " << QString::asprintf("%.0f deg", perp_angle * 180.0 / M_PI);


    //double ta_perp = wrapRadians2PI(ta - M_PI / 2.0);
    //double ta_perp = originalAngleFromPerpAngle(ta);

    /*QString coord_txt = QString::asprintf("(%.1fd, %.2f)", perp_angle *180.0/M_PI, d);
    ctx->fillText(coord_txt, tx2 + off.x, ty2 + off.y);

    DVec2 pt = fromPolarCoordinate(perp_angle, 1);
    ctx->setFillStyle(255, 255, 255);
    ctx->fillEllipse(pt.x, pt.y, 0.015);*/
}

/// BasicProject ///

void Cardioid_Project::projectPrepare(Layout& layout)
{
    layout << create<Cardioid_Scene>();
}

/// Scene ///

void Cardioid_Scene::UI::sidebar()
{
    bl_scoped(flatten);
    bl_scoped(interactive);
    bl_scoped(interact_angle_step);
    bl_scoped(interact_spin_mult);
    bl_scoped(interact_angle);
    bl_scoped(interact_dist);

    bl_scoped(ani_angle);
    bl_scoped(ani_inc);
    bl_scoped(animate);

    //ImGui::SliderDouble("Angle", &interact_angle, 0, (2 * M_PI));
    ImGui::Checkbox("Flatten", &flatten);
    ImGui::Checkbox("Interactive", &interactive);
    ImGui::SliderDouble("Rotate Speed", &interact_angle_step, 0, Math::TWO_PI / 100.0);


    ImGui::SliderDouble("Spin multiplier", &interact_spin_mult, 0.0, 1.0);
    ImGui::SliderDouble("Distance", &interact_dist, 0.0, 1.0);


    ImGui::Checkbox("Animate Rotation", &animate);
    ImGui::SliderAngle("Angle", &ani_angle);
    ImGui::SliderAngle("Angle Inc", &ani_inc, -Math::HALF_PI, Math::HALF_PI, 1);

    //ImGui::Checkbox("show offset", &show_offset);
    //ImGui::Checkbox("show original", &show_original);
    //ImGui::Checkbox("show alternative", &show_alternative);

    ///static float v[5] = { 0.390f, 0.575f, 0.565f, 1.000f };
    ///ImGui::Bezier( "easeOutSine", v );       // draw
    ///float y = ImGui::BezierValue( 0.5f, v ); // x delta in [0..1] range
}

void Cardioid_Scene::sceneStart()
{
    //cumulative_cardioid = computeCumulativeCardioid((2 * M_PI) / 72000.0);
    cumulative_cardioid_lookup.create(Math::TWO_PI / 3600.0, 0.01);

}

void Cardioid_Scene::sceneMounted(Viewport* ctx)
{
    /// Initialize viewport (after sceneStart)
    camera.setSurface(ctx);
    camera.setOriginViewportAnchor(Anchor::CENTER);
    camera.focusWorldRect(-0.9, -1, 0.6, 1);

    navigator.setTarget(camera);
    navigator.setDirectCameraPanning(true);
    //camera.restrictRelativeZoomRange(0.001, 100);
}

void Cardioid_Scene::sceneProcess()
{
    if (isRecording())
        captureFrame(false);

    if (animate)
    {
        if (interactive)
        {

        }
        else
        {
            ani_angle += ani_inc;
            ani_angle = Math::wrapRadians2PI(ani_angle);

            if (isRecording())
                captureFrame(true);
        }
    }
}

double originalAngleFromPerpAngle(double perp_angle)
{
    return Math::wrapRadians2PI((perp_angle + Math::HALF_PI) / 1.5);
}

void Cardioid_Scene::viewportProcess(Viewport*, double)
{
    /// Process Viewports running this Scene
    //interact_angle = originalAngleFromPoint(mouse->world_x, mouse->world_y);

}

void Cardioid_Scene::viewportDraw(Viewport* ctx) const
{
    ctx->setTransform(camera.getTransform());

    /// Draw Scene to Viewport
    ctx->drawWorldAxis();
    ctx->worldHudMode();

    double ox = show_offset ? -0.25 : 0;

    
    //plotCumulativeCardioid(ctx, cumulative_cardioid, interact_spin_mult);


    ctx->scalingSizes(false);
    ctx->scalingLines(false);

    if (animate)
    {
        ///if (show_alternative)
        ///    fullPlotAlternative(ctx, 1, ox);
         /// 
        if (interactive)
        {
            double tx1, ty1, ta, oa, d;
            cardioidPolarCoord((f64)mouse->world_x, (f64)mouse->world_y, tx1, ty1, ta, d, oa);

            fullPlot(ctx, 1, ox);
            animatePlot(ctx, 1, ox, oa, d);
        }
        else
        {
            fullPlot(ctx, 1, ox);
            animatePlot(ctx, 1, ox, ani_angle, interact_dist);
        }
    }
    else
    {
        if (flatten)
        {
            if (interactive)
            {
                ctx->beginPath();
                ctx->drawPath(cumulative_cardioid_lookup.lerped(interact_spin_mult));
                ctx->stroke();

                // Draw red projected dot
                DVec2 p = cumulative_cardioid_lookup.project(interact_angle, interact_dist, interact_spin_mult);
                ctx->setFillStyle(255, 0, 0);
                ctx->fillEllipse(p.x, p.y, 5.0);

                // Mouse project test (should follow mouse)
                DVec2 p2 = cumulative_cardioid_lookup.originalPolarCoordinate((f64)mouse->world_x, (f64)mouse->world_y, interact_spin_mult);
                DVec2 p3 = cumulative_cardioid_lookup.project(p2.x, p2.y, interact_spin_mult);

                ctx->setFillStyle(255, 0, 0);
                ctx->fillEllipse(p3.x, p3.y, 5.0);
            }
        }
        else
        {
            fullPlot(ctx, 1, ox);
            animatePlot(ctx, 1, ox, ani_angle, interact_dist);
        }
    }

    //ctx->print() << "\nfps: " << fps(60) << " ms";
}

void Cardioid_Scene::onEvent(Event e)
{
    navigator.handleWorldNavigation(e, true);
}

double originalAngle(double p2_x, double p2_y)
{
    double c2 = (1.0 - std::sqrt(3.0 - 8.0 * p2_x)) / 2.0;
    return std::atan2((2.0 * p2_y) / (1.0 - c2), c2);
}

void Cardioid_Scene::animatePlot(Viewport* ctx, double scale, double ox, double orig_angle, double dist) const
{
    ctx->scalingSizes(true);
    ctx->scalingLines(false);

    ctx->setLineWidth(2);
    
    // Big blue circle
    double r1 = scale * 0.5;
    ctx->setStrokeStyle(80, 80, 255);
    ctx->strokeEllipse(ox * scale, 0.0, r1);

    // p1 dot
    double p1_x = r1 * cos(orig_angle) + ox * scale;
    double p1_y = r1 * sin(orig_angle);

    // Med green circle
    double r2 = r1 / 4.0;
    double angle2 = 2.0 * orig_angle;
    double c2_x = p1_x - r2 * cos(angle2);
    double c2_y = p1_y - r2 * sin(angle2);

    ctx->setStrokeStyle(0, 255, 0);
    ctx->strokeEllipse(c2_x, c2_y, r2);

    DVec2 off = ctx->Offset(5, 5);

    ctx->scalingSizes(false);

    // draw p1 dot
    ctx->setFillStyle(0, 255, 255);
    ctx->fillEllipse(p1_x, p1_y, 5.0);
    ctx->fillText("p1", p1_x + off.x, p1_y + off.y);

    // draw p2 dot
    double p2_x = c2_x - r2 * cos(angle2);
    double p2_y = c2_y - r2 * sin(angle2);
    ctx->setFillStyle(0, 255, 255);
    ctx->fillEllipse(p2_x, p2_y, 5.0);
    ctx->fillText("p2", p2_x + off.x, p2_y + off.y);

    double tangent_angle = 1.5 * orig_angle;

    // draw tangent arrow
    double tangent_head_x = p2_x + cos(tangent_angle) * 0.2;
    double tangent_head_y = p2_y + sin(tangent_angle) * 0.2;
    ctx->drawArrow(DVec2{ p2_x, p2_y }, DVec2{ tangent_head_x, tangent_head_y }, Color(255, 255, 255));
    //ctx->beginPath();
    //ctx->arrowMoveTo(p2_x, p2_y);
    //ctx->arrowDrawTo(tangent_head_x, tangent_head_y);
    //ctx->stroke();
    
    
    // Draw perpendicular arrow (to mouse)
    //double perp_angle = Math::wrapRadians2PI(tangent_angle - M_PI / 2.0);
    ///tx2 = tx1 + cos(perp_angle) * d;
    ///ty2 = ty1 + sin(perp_angle) * d;
    DVec2 tp = fromPolarCoordinate(orig_angle, dist);

    double tx1 = 0.5 * cos(orig_angle) - 0.25 * cos(orig_angle * 2.0);
    double ty1 = 0.5 * sin(orig_angle) - 0.25 * sin(orig_angle * 2.0);
    double tx2 = tp.x;
    double ty2 = tp.y;

    //ctx->beginPath();
    //ctx->arrowMoveTo(tx1, ty1);
    //ctx->arrowDrawTo(tx2, ty2);
    //ctx->stroke();

    ctx->drawArrow(DVec2{ tx1, ty1 }, DVec2{ tx2, ty2 }, Color(255, 255, 255));

    ///ctx->print() << "p1 angle: " << QString::asprintf("%.0f deg", oa * 180.0 / M_PI);
    ///ctx->print() << "\np2 tangent angle: " << QString::asprintf("%.0f deg", ta * 180.0 / M_PI);
    ///ctx->print() << "\np2 perp angle: " << QString::asprintf("%.0f deg", perp_angle * 180.0 / M_PI);


    //double ta_perp = wrapRadians2PI(ta - M_PI / 2.0);
    //double ta_perp = originalAngleFromPerpAngle(ta);

    /*QString coord_txt = QString::asprintf("(%.1fd, %.2f)", perp_angle *180.0/M_PI, d);
    ctx->fillText(coord_txt, tx2 + off.x, ty2 + off.y);

    DVec2 pt = fromPolarCoordinate(perp_angle, 1);
    ctx->setFillStyle(255, 255, 255);
    ctx->fillEllipse(pt.x, pt.y, 0.015);*/
}

void Cardioid_Scene::plotCumulativeCardioid(
    Viewport* ctx, 
    const CumulativeCardioid& cardioid,
    double angle_mult)
{
    double plot_x = 0.25;
    double plot_y = 0.0;
    double plot_direction = 0.0;

    ctx->beginPath();
    ctx->moveTo(plot_x, plot_y);

    for (int i = 0; i < cardioid.size(); i++)
    {
        double step_angle = cardioid[i].step_angle;
        double step_dist = cardioid[i].step_dist;

        double scaled_step_angle = step_angle * angle_mult;
        plot_direction += scaled_step_angle;

        plot_x += cos(plot_direction) * step_dist;
        plot_y += sin(plot_direction) * step_dist;

        ctx->lineTo(plot_x, plot_y);
    }

    ctx->stroke();
}

void Cardioid_Scene::fullPlot(Viewport* ctx, double scale, double ox) const
{
    //camera->stageMode();

    // Full plot
    ctx->setStrokeStyle(255, 0, 0);
    ctx->beginPath();
    ctx->setLineWidth(3);

    double r1 = scale * 0.5;
    double r2 = r1 / 4.0;

    double angle_step = Math::TWO_PI / 100.0;

    bool first = true;
    for (double angle = 0.0; angle < Math::TWO_PI; angle += angle_step)
    {
        //double straight_len = (angle / (Math::TWO_PI)) * 3.837;

        double angle2 = 2.0 * angle;

        double c1_x = r1 * cos(angle) + ox * scale;
        double c1_y = r1 * sin(angle);
        double c2_x = c1_x - r2 * cos(angle2);
        double c2_y = c1_y - r2 * sin(angle2);

        double plot_x = c2_x - r2 * cos(angle2);
        double plot_y = c2_y - r2 * sin(angle2);

        //double plot_x = 0.5 * cos(angle) - 0.25 * cos(angle * 2.0);
        //double plot_y = 0.5 * sin(angle) - 0.25 * sin(angle * 2.0);

        if (first)
        {
            ctx->moveTo(plot_x, plot_y);
            first = false;
        }
        else
        {
            ctx->lineTo(plot_x, plot_y);
        }
    }

    ctx->stroke();
}

void Cardioid_Scene::fullPlotAlternative(Viewport* ctx, double scale, double ox) const
{
    // Full plot
    ctx->setStrokeStyle(255, 0, 255);
    ctx->beginPath();

    ctx->setLineWidth(4);
    ctx->setStrokeStyle(255, 0, 255);

    double angle_step = Math::TWO_PI / 100.0;

    bool first = true;
    for (double angle = 0.0; angle < Math::TWO_PI; angle += angle_step)
    {
        //double a = sin(angle * 0.5);
        double plot_x = cos(angle) * pow(sin(angle * 0.5), 2) + (ox + 0.25) * scale;
        double plot_y = sin(angle) * pow(sin(angle * 0.5), 2);

        if (first)
        {
            ctx->moveTo(plot_x, plot_y);
            first = false;
        }
        else
        {
            ctx->lineTo(plot_x, plot_y);
        }
    }
    ctx->stroke();
}

void Cardioid_Graph_Scene::sceneStart()
{
}

void Cardioid_Graph_Scene::sceneMounted(Viewport* ctx)
{
    camera.setSurface(ctx);
    camera.setOriginViewportAnchor(Anchor::CENTER);
    camera.focusWorldRect(-0.9, -1, 0.6, 1);
}

void Cardioid_Graph_Scene::viewportProcess(Viewport* ctx, double)
{

    int iw = static_cast<int>(ctx->width() / 2);
    int ih = static_cast<int>(ctx->height() / 2);

    bmp.setBitmapSize(iw, ih);
    bmp.setStageRect(0, 0, ctx->width(), ctx->height());

    if (bmp.needsReshading())
    {
        // Tangent angle heatmap
        int current_row = 0;
        bmp.forEachWorldPixel<double>(current_row, [this](int x, int y, double wx, double wy)
        {
            double tx1, ty1, ta, oa, d;
            cardioidPolarCoord(wx, wy, tx1, ty1, ta, d, oa);

            if (d >= 0)
            {
                double perp_angle = Math::wrapRadians(ta - Math::HALF_PI);

                double neg_angle = std::max(0.0, -perp_angle);
                double pos_angle = std::max(0.0, perp_angle);

                int neg_col = (int)((std::max(0.0, std::min(Math::PI, neg_angle)) / Math::PI) * 255.0);
                int pos_col = (int)((std::max(0.0, std::min(Math::PI, pos_angle)) / Math::PI) * 255.0);

                bmp.setPixel(x, y, neg_col, pos_col, 0, 255);
            }
            else
            {
                bmp.setPixel(x, y, 0, 0, 0, 255);
            }
        });

        // Tangent angle heatmap
        /*std::vector<QFuture<void>> futures(thread_count);

        int pixel_count = iw * ih;
        auto pixel_ranges = splitRanges(pixel_count, thread_count);

        for (int ti = 0; ti < thread_count; ti++)
        {
            auto& pixel_range = pixel_ranges[ti];
            futures[ti] = QtConcurrent::run([&]()
            {
                for (int i = pixel_range.first; i < pixel_range.second; i++)
                {
                    int x = (i % iw), y = (i / iw);
                    double wx = bmp.wx(x);
                    double wy = bmp.wy(y);

                    double tx1, ty1, tx2, ty2, ta, oa, d;
                    cardioidPolarCoord(wx, wy, tx1, ty1, ta, d, oa);

                    if (d >= 0)
                    {
                        //double orig_angle = originalAngleFromPoint(wx, wy);
                        //double perp_angle = wrapRadians2PI(1.5 * orig_angle - M_PI / 2.0) - M_PI;
                        //double perp_angle = ta - M_PI;// wrapRadians2PI(ta - M_PI / 2.0) - M_PI;
                        double perp_angle = wrapRadians(ta - M_PI/2.0);// wrapRadians2PI(ta - M_PI / 2.0) - M_PI;

                        double neg_angle = max(0.0, -perp_angle);
                        double pos_angle = max(0.0, perp_angle);

                        int neg_col = (int)((max(0.0, min(M_PI, neg_angle)) / M_PI) * 255.0);
                        int pos_col = (int)((max(0.0, min(M_PI, pos_angle)) / M_PI) * 255.0);

                        bmp.setPixel(x, y, neg_col, pos_col, 0, 255);
                    }
                    else
                    {
                        bmp.setPixel(x, y, 0, 0, 0, 255);
                    }
                }
            });
        }
 
        for (auto& future : futures)
            future.waitForFinished();*/
    }
}

void Cardioid_Graph_Scene::viewportDraw(Viewport* ctx) const
{
    /// Bmp plotting
    //ctx->drawImage(bmp);
    ctx->drawWorldAxis();
    ctx->worldHudMode();
    
    double angle_step = Math::TWO_PI / 720.0;
    bool first;

    /// Show unstable region bound
    /*ctx->setLineWidth(2);
    ctx->setStrokeStyle(0, 0, 0);
    ctx->beginPath();
    ctx->moveTo(0.25, 0);
    for (double x = 0.25; x < 50.0; x += 0.01)
    {
        double y_bound = 0.001 + 4 * pow(x - 0.24, 2);
        ctx->lineTo(x, y_bound);
    }
    ctx->stroke();*/

    // Full plot
    ctx->setLineWidth(2);
    ctx->setStrokeStyle(255, 255, 255);
    ctx->beginPath();

    first = true;
    for (double angle = 0.0; angle < Math::TWO_PI; angle += angle_step)
    {
        double plot_x = 0.5 * cos(angle) - 0.25 * cos(angle * 2.0);
        double plot_y = 0.5 * sin(angle) - 0.25 * sin(angle * 2.0);

        if (first)
        {
            ctx->moveTo(plot_x, plot_y);
            first = false;
        }
        else
        {
            ctx->lineTo(plot_x, plot_y);
        }
    }

    ctx->stroke();
}


SIM_END;
