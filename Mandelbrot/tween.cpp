#include "Mandelbrot.h"
#include "conversions.h"

SIM_BEG;

f128 tweenDistance(
    MandelState& a,
    MandelState& b)
{
    f128 dh = toHeight(b.camera.getRelativeZoom<f128>()) - toHeight(a.camera.getRelativeZoom<f128>());
    f128 dx = b.camera.x<f128>() - a.camera.x<f128>();
    f128 dy = b.camera.y<f128>() - a.camera.y<f128>();
    f128 d = sqrt(dx * dx + dy * dy + dh * dh);
    return d;
}

void startTween(Mandelbrot_Scene &scene)
{
    // Now, set start/end tween states
    scene.state_a = scene;
    //scene.state_b = target;

    // Switch to raw iter_lim for tween (and switch to quality mode on finish)
    scene.dynamic_iter_lim = false;
    scene.quality = scene.iter_lim;


    scene.tween_r1 = scene.getAngledRect(scene.state_a);
    scene.tween_r2 = scene.getAngledRect(scene.state_b);

    DDAngledRect encompassing;
    encompassing.fitTo(scene.tween_r1, scene.tween_r2, scene.tween_r1.aspectRatio());

    f128 encompassing_zoom = (scene.stageWorldSize() / encompassing.size).average();
    f128 encompassing_height = std::min(f128(1.0), toHeight(encompassing_zoom)); // Cap at max height
    scene.tween_lift = encompassing_height - std::max(
        toHeight(scene.state_a.camera.getRelativeZoom<f128>()),
        toHeight(scene.state_b.camera.getRelativeZoom<f128>())
    );

    f128 max_lift = 1.0 - toHeight(scene.state_b.camera.getRelativeZoom<f128>());
    if (max_lift < 0.0) max_lift = 0;
    if (scene.tween_lift < 0) scene.tween_lift = 0;
    if (scene.tween_lift > max_lift) scene.tween_lift = max_lift;


    // Begin tween
    scene.tween_progress = 0.0;
    scene.tweening = true;
    scene.tween_duration = 2;// tweenDistance(scene_data.state_a, scene_data.state_b);
    scene.dist_params.cycle_dist_invert = scene.state_b.dist_params.cycle_dist_invert;
}

void lerpState(
    Mandelbrot_Scene& scene,
    MandelState& a,
    MandelState& b,
    double f,
    bool complete)
{
    MandelState& dst = scene;

    float pos_f = scene.tween_pos_spline((float)f);

    // === Get true src/dst iter_lim for consistent tweening (ignore dynamic calculation until tween complete) ===
    double src_iter_lim = a.dynamic_iter_lim ?
        (mandelbrotIterLimit(a.camera.getRelativeZoom<f128>()) * a.quality) :
        a.quality;

    double dst_iter_lim = b.dynamic_iter_lim ?
        (mandelbrotIterLimit(b.camera.getRelativeZoom<f128>()) * b.quality) :
        b.quality;

    // === Lerp Camera View and normalized zoom from "height" ===
    double lift_weight = (double)scene.tween_zoom_lift_spline((float)f);
    f128 lift_height = scene.tween_lift * lift_weight;
    f128 a_height = toHeight(a.camera.getRelativeZoom<f128>());
    f128 b_height = toHeight(b.camera.getRelativeZoom<f128>());

    double base_zoom_f = scene.tween_base_zoom_spline((float)f);
    f128 dst_height = Math::lerp(a_height, b_height, base_zoom_f) + lift_height;
    CameraInfo::lerp(dst.camera, a.camera, b.camera, pos_f);
    dst.camera.setRelativeZoom(fromHeight(dst_height)); // override zoom from computed "height"

    // === Quality ===
    dst.quality = Math::lerp(src_iter_lim, dst_iter_lim, pos_f);

    // === Color Cycle ===
    float color_cycle_f = scene.tween_color_cycle((float)f);

    // === Shader weights ===
    dst.iter_weight = Math::lerp(a.iter_weight, b.iter_weight, color_cycle_f);
    dst.dist_weight = Math::lerp(a.dist_weight, b.dist_weight, color_cycle_f);
    dst.stripe_weight = Math::lerp(a.stripe_weight, b.stripe_weight, color_cycle_f);

    // === Iter shader === 
    dst.iter_params.cycle_iter_value = Math::lerp(a.iter_params.cycle_iter_value, b.iter_params.cycle_iter_value, color_cycle_f);
    dst.iter_params.cycle_iter_log1p_weight = Math::lerp(a.iter_params.cycle_iter_log1p_weight, b.iter_params.cycle_iter_log1p_weight, color_cycle_f);

    // === Dist Shader
    dst.dist_params.cycle_dist_value = Math::lerp(a.dist_params.cycle_dist_value, b.dist_params.cycle_dist_value, color_cycle_f);
    dst.dist_params.cycle_dist_sharpness = Math::lerp(a.dist_params.cycle_dist_sharpness, b.dist_params.cycle_dist_sharpness, color_cycle_f);

    // === Stripe Shader
    dst.stripe_params.freq = Math::lerp(a.stripe_params.freq, b.stripe_params.freq, color_cycle_f);
    dst.stripe_params.phase = Math::lerp(a.stripe_params.phase, b.stripe_params.phase, color_cycle_f);
    dst.stripe_params.contrast = Math::lerp(a.stripe_params.contrast, b.stripe_params.contrast, color_cycle_f);

    // === Lerp Gradient/Hue Shift===
    dst.gradient_shift = Math::lerp(a.gradient_shift, b.gradient_shift, pos_f);
    dst.hue_shift = Math::lerp(a.hue_shift, b.hue_shift, pos_f);

    // === Lerp animation speed for Gradient/Hue Shift ===
    dst.gradient_shift_step = Math::lerp(a.gradient_shift_step, b.gradient_shift_step, pos_f);
    dst.hue_shift_step = Math::lerp(a.hue_shift_step, b.hue_shift_step, pos_f);

    // === Lerp Color Gradient ===
    ImGradient::lerp(scene.gradient, a.gradient, b.gradient, (float)f);

    if (complete)
    {
        dst.dynamic_iter_lim = b.dynamic_iter_lim;
        dst.quality = b.quality;

        dst.iter_params.cycle_iter_dynamic_limit = b.iter_params.cycle_iter_dynamic_limit;

        dst.show_axis = b.show_axis;
        dst.show_color_animation_options = b.show_color_animation_options;
    }

    // === Lerp quality ===
    //dst->iter_lim = Math::lerp(state_a.iter_lim, state_b.iter_lim, f);

    // === Lerp Cardioid Flattening Factor ===
    //dst->cardioid_lerp_amount = Math::lerp(state_a.cardioid_lerp_amount, state_b.cardioid_lerp_amount, f);

    // Spline Data
    //memcpy(dst->x_spline_point, Math::lerp(state_a.x_spline_points, state_b.x_spline_points, f));
}

SIM_END;
