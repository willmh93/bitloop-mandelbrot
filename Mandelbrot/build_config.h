#pragma once

/// ──────  Dev Mode ON/OFF ──────
#define MANDEL_DEV_MODE 1

#if MANDEL_DEV_MODE //&& defined(BL_DEBUG)

    /// ──────  Enabled Dev Options ──────
#define MANDEL_DEV_PRINT_ACTIVE_FLOAT_TYPE          0
#define MANDEL_DEV_PRINT_ACTIVE_COMPUTE_FEATURES    0
#define MANDEL_DEV_PERFORMANCE_TIMERS               0
#define MANDEL_DEV_EDIT_TWEEN_SPLINES               0
#define MANDEL_DEV_HIGHLIGHT_MAXDEPTH_FORWARDING    1

#endif


/// ──────  Feature toggles ──────
#define MANDEL_FEATURE_CAMERA_EASING            1
#define MANDEL_FEATURE_INTERACTIVE_CARDIOID     1
#define MANDEL_FEATURE_SPLINE_MODE              0 // todo
#define MANDEL_FEATURE_FLATTEN_MODE             0 // todo

