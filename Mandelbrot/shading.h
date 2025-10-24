#pragma once
#include <bitloop.h>

#include "types.h"
#include "mandel_state.h"

SIM_BEG;

enum class MandelKernelFeatures
{
    NONE = 0,
    ITER = 1,
    DIST = 2,
    STRIPES = 4,
    MIX_ALL = 7,
    COUNT
};

enum class MandelTransform
{
    NONE,
    FLATTEN
};

enum struct MandelShaderFormula
{
    ITER_DIST_STRIPE,          // a+b+c
    ITER_DIST__MUL__STRIPE,   // (a+b)*c
    ITER__MUL__DIST_STRIPE,  // a*(b+c)
    ITER_STRIPE__MULT__DIST,  // a*(b+c)
    COUNT
    //~ // todo: Always tween stripe weight to 0 before switching mode
};

static const char* MandelFormulaNames[(int)MandelShaderFormula::COUNT] = {
    "iter + dist + stripe",
    "iter + (dist * stripe)",
    "(iter * dist) + stripe",
    "(iter * stripe) * dist"
};

static const char* MandelSmoothingNames[(int)MandelKernelFeatures::COUNT] = {
    "NONE",
    "ITER",
    "DIST",
    "STRIPE"
};

SIM_END;
