#pragma once
#include <bitloop.h>

SIM_BEG;

using namespace bl;

namespace MandelSplines
{
    static ImSpline::Spline stripe_zf_spline = ImSpline::Spline(100, {
        {-0.48f,-0.13f}, {0.0f,0.0f},   {1.49f,0.33f},
        {2.21f,0.07f},   {11.99f,0.2f}, {29.42f,0.43f},
        {30.94f,0.45f},  {38.8f,0.56f}, {48.85f,0.7f}
    });
}

SIM_END;
