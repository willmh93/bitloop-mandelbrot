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
        {-0.147, 0.000}, {0.253,0.000f}, {0.553,0.000f},
        {0.3710, 1.000}, {0.675,1.000f}, {1.070,1.000f}
    });
    static ImSpline::Spline tween_zoom_lift_spline = ImSpline::Spline(100, {
        {-0.100,0.000}, {0.000, 0.000}, {0.100, 0.000},
        {0.2010,0.802}, {0.251, 1.002}, {0.351, 1.402},
        {0.5720,1.419}, {0.672, 0.799}, {0.800, 0.000},
        {0.8000,0.000}, {1.000, 0.000}, {1.104, 0.000}
    });
    static ImSpline::Spline tween_base_zoom_spline = ImSpline::Spline(100, {
        {-0.3f, 0.0f}, {0.0f, 0.0f}, {0.5f, 0.0f},
        {0.50f, 1.0f}, {1.0f, 1.0f}, {1.5f, 1.0f}
    });
    static ImSpline::Spline tween_color_cycle = ImSpline::Spline(100, {
        {0.072f, 0.000f}, {0.500f, 0.000f}, {0.845f, 0.000f},
        {0.750f, 1.000f}, {1.000f, 1.000f}, {1.250f, 1.000f}
    });
}

SIM_END;
