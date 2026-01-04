#pragma once

#define BL_TIMERS_ENABLED                           0
#define BL_SIMD_FORCE_SCALAR                        1 // faster as scalar for current implementation

/// ──────  Enable/Disable all Dev options ──────
#define MANDEL_DEV_MODE                             0


#if MANDEL_DEV_MODE // && defined(BL_DEBUG)
    /// ──────  Enabled Dev Options ──────
    #define MANDEL_DEV_EDIT_TWEEN_SPLINES           0
    #define MANDEL_EXTENDED_FIELD_STATS             0
    #define MANDEL_EXPERIMENTAL_TESTS               0
#endif


/// ──────  Feature toggles ──────
#define MANDEL_FEATURE_CAMERA_EASING                1
#define MANDEL_FEATURE_INTERACTIVE_CARDIOID         1
#define MANDEL_FEATURE_SPLINE_MODE                  0 // todo
#define MANDEL_FEATURE_FLATTEN_MODE                 0 // todo

