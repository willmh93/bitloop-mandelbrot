#pragma once
#include <bitloop.h>

SIM_BEG;

using namespace bl;

namespace MandelSplines
{
    static ImSpline::Spline x_spline = ImSpline::Spline(100, {
        {0.0f, 0.0f}, {0.1f, 0.1f}, {0.2f, 0.2f},
        {0.3f, 0.3f}, {0.4f, 0.4f}, {0.5f, 0.5f},
        {0.6f, 0.6f}, {0.7f, 0.7f}, {0.8f, 0.8f}
    });
    static ImSpline::Spline y_spline = ImSpline::Spline(100, {
        {0.0f, 0.0f}, {0.1f, 0.1f}, {0.2f, 0.2f},
        {0.3f, 0.3f}, {0.4f, 0.4f}, {0.5f, 0.5f},
        {0.6f, 0.6f}, {0.7f, 0.7f}, {0.8f, 0.8f}
    });
    static ImSpline::Spline tween_pos_spline = ImSpline::Spline(100, {
        {-0.147f, 0.000f}, {0.253f,0.000f}, {0.553f,0.000f},
        {0.3710f, 1.000f}, {0.675f,1.000f}, {1.070f,1.000f}
    });
    static ImSpline::Spline tween_zoom_lift_spline = ImSpline::Spline(100, {
        {-0.100f, 0.000f}, {0.000f, 0.000f}, {0.100f, 0.000f},
        {0.2010f, 0.802f}, {0.251f, 1.002f}, {0.351f, 1.402f},
        {0.5720f, 1.419f}, {0.672f, 0.799f}, {0.800f, 0.000f},
        {0.8000f, 0.000f}, {1.000f, 0.000f}, {1.104f, 0.000f}
    });
    //static ImSpline::Spline tween_base_zoom_spline = ImSpline::Spline(100, {
    //    {-0.3f, 0.0f}, {0.251f, 0.0f}, {0.5f, 0.0f},
    //    {0.50f, 1.0f}, {1.0f, 1.0f}, {1.5f, 1.0f}
    //});
    static ImSpline::Spline tween_base_zoom_spline = ImSpline::Spline(100,
        { { -0.5, 0 }, { 0,0 }, { 0.5,0 }, { 0.5,1 }, { 1,1 }, { 1.5,1 } }
    );
    static ImSpline::Spline tween_color_cycle = ImSpline::Spline(100, {
        {0.072f, 0.000f}, {0.500f, 0.000f}, {0.845f, 0.000f},
        {0.750f, 1.000f}, {1.000f, 1.000f}, {1.250f, 1.000f}
    });
    static ImSpline::Spline stripe_zf_spline = ImSpline::Spline(100, {
        {-0.48f,-0.13f}, {0.0f,0.0f},   {1.49f,0.33f},
        {2.21f,0.07f},   {11.99f,0.2f}, {29.42f,0.43f},
        {30.94f,0.45f},  {38.8f,0.56f}, {48.85f,0.7f}
    });
}

SIM_END;
