#pragma once

/// ──────  Dev Mode ON/OFF ──────
#define MANDEL_DEV_MODE 1

#if MANDEL_DEV_MODE //&& defined(BL_DEBUG)

    /// ──────  Enabled Dev Options ──────
#define MANDEL_DEV_EDIT_TWEEN_SPLINES           0

#endif


/// ──────  Feature toggles ──────
#define MANDEL_FEATURE_CAMERA_EASING            1
#define MANDEL_FEATURE_INTERACTIVE_CARDIOID     1
#define MANDEL_FEATURE_SPLINE_MODE              0 // todo
#define MANDEL_FEATURE_FLATTEN_MODE             0 // todo

