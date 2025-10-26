
#include "Mandelbrot.h"

SIM_BEG;

/// ─────────────────────── Project ───────────────────────

void Mandelbrot_Project::projectPrepare(Layout& layout)
{
    Mandelbrot_Scene::Config config;
    create<Mandelbrot_Scene>(config)->mountTo(layout);
}

/// ─────────────────────── Scene ───────────────────────

void Mandelbrot_Scene::sceneStart()
{
    // todo: Stop this getting called twice on startup
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

    bmp_9x9.setCamera(camera);
    bmp_3x3.setCamera(camera);
    bmp_1x1.setCamera(camera);

    navigator.setTarget(camera);
    navigator.setDirectCameraPanning(true);
    //camera->restrictRelativeZoomRange(0.5, 1e+300);


    #ifdef __EMSCRIPTEN__
    // If URL has encoded state data, load on startup
    if (platform()->url_has("data"))
        loadState(platform()->url_get_string("data"));
    #endif

    /*if (!load_example_name.empty())
    {
        data_buf = mandel_presets[load_example_name];
        if (!data_buf.empty())
            loadConfigBuffer();
    }*/
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

    double zoom_mag = camera.relativeZoom<double>();

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
    // Show intro message on startup
    if (display_intro)
    {
        ctx->setFont(font);
        ctx->setFontSize(20);

        if (platform()->is_mobile())
        {
            // touchscreen instructions
            ctx->fillText("Controls:", scale_size(10.0), scale_size(10.0));
            ctx->fillText("  - Touch & drag to move", scale_size(10.0), scale_size(35.0));
            ctx->fillText("  - Pinch to zoom / rotate", scale_size(10.0), scale_size(60.0));
        }
        else
        {
            // desktop instructions
            ctx->fillText("Controls:", scale_size(10.0), scale_size(10.0));
            ctx->fillText("  - Click & drag to move", scale_size(10.0), scale_size(35.0));
            ctx->fillText("  - Mouse wheel to zoom", scale_size(10.0), scale_size(60.0));
        }

        ctx->setTextAlign(TextAlign::ALIGN_LEFT);
        ctx->setTextBaseline(TextBaseline::BASELINE_BOTTOM);
        ctx->fillText("Contact:  will.hemsworth@bitloop.dev", scale_size(10.0), ctx->height() - scale_size(10.0));
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


//void Mandelbrot_Scene::onSavefileChanged()
//{
//    data_buf = serialize();
//    #ifdef __EMSCRIPTEN__
//    // Appears to slow down browser during animation - call less often or rely on save/load/share buttons only
//    //platform()->url_set_string("data", data_buf.c_str());
//    #endif
//}


void Mandelbrot_Scene::onEvent(Event e)
{
    #ifdef BL_RELEASE
    // don't permit world navigation while tweening
    if (tweening || isRecording())
        return;
    #endif

    // enter/exit fullscreen hotkeys (F11 / Escape)
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
    f128 old_zoom = camera.relativeZoom<f128>();
    
    if (navigator.handleWorldNavigation(e, true, true))
    {
        f128 zoom = camera.relativeZoom<f128>();
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
