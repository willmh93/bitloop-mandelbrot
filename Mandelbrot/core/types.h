#pragma once
#include <bitloop.h>

#include <cmath>
#include <vector>

enum class KernelFeatures
{
    NONE    = 0b000,
    ITER    = 0b001,
    DIST    = 0b010,
    STRIPES = 0b100,
    MIX_ALL = 0b111,

    COUNT
};

enum struct KernelMode : int
{
    NO_PERTURBATION              = 0,
    PERTURBATION                 = 1,
    PERTURBATION_SIMD            = 2, // kept for development / easier debugging
    PERTURBATION_SIMD_UNROLLED   = 3,

    COUNT,
    AUTO = COUNT, // Same value as COUNT to avoid being included in constexpr_dispatch 'build_table'
    PERTURBATION_MASK = PERTURBATION | PERTURBATION_SIMD | PERTURBATION_SIMD_UNROLLED
};

// todo: find way to put inside SIM_BEG ns (wasm32 error, must be in global)
bl_enable_enum_bitops(KernelFeatures);
bl_enable_enum_bitops(KernelMode);

SIM_BEG;

using namespace bl;

constexpr double INSIDE_MANDELBROT_SET = std::numeric_limits<double>::max();
const double INSIDE_MANDELBROT_SET_SKIPPED = std::nextafter(INSIDE_MANDELBROT_SET, 0.0f);

enum MandelFlag : uint32_t
{
    // bits
    MANDEL_DYNAMIC_ITERS        = 1u << 0,
    MANDEL_SHOW_AXIS            = 1u << 1,
    MANDEL_FLATTEN              = 1u << 2,
    MANDEL_DYNAMIC_COLOR_CYCLE  = 1u << 3,
    MANDEL_NORMALIZE_DEPTH      = 1u << 4,
    MANDEL_INVERT_DIST          = 1u << 5,
    MANDEL_USE_SMOOTHING        = 1u << 6,

    // bitmasks
    //MANDEL_FLAGS_MASK   = 0x000FFFFFu, // max 24 bit-flags
    //MANDEL_SMOOTH_MASK  = 0x00F00000u, // max 16 smooth types
    MANDEL_VERSION_MASK   = 0xFF000000u, // max 255 versions

    MANDEL_VERSION_BITSHIFT = 24
};

enum struct MandelInteriorForwarding
{
    SLOWEST,
    SLOW,
    MEDIUM,
    FAST,
    COUNT
};

static inline const char* MandelMaxDepthOptimizationNames[(int)MandelInteriorForwarding::COUNT] =
{
    "Slowest (Lossless)",
    "Slow",
    "Balanced",
    "Fast (Loss at sharp angles)"
};

static inline const char* FloatingPointTypeNames[3] =
{
    "32-bit",
    "64-bit",
    "128-bit (quad-precision)"
};


SIM_END;
