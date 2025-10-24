#include "Mandelbrot.h"
#include "conversions.h"

SIM_BEG;

/*f128 tweenDistance(
    MandelState& a,
    MandelState& b)
{
    f128 dh = toHeight(b.camera.relativeZoom<f128>()) - toHeight(a.camera.relativeZoom<f128>());
    f128 dx = b.camera.x<f128>() - a.camera.x<f128>();
    f128 dy = b.camera.y<f128>() - a.camera.y<f128>();
    f128 d = sqrt(dx * dx + dy * dy + dh * dh);
    return d;
}*/

void Mandelbrot_Scene::startTween(const MandelState& target)
{
    // Now, set start/end tween states
    this->state_a = *this;
    this->state_b = target;

    this->state_b.camera.setReferenceZoom(camera.getReferenceZoom<f128>());

    // Switch to raw iter_lim mode for tween (and switch to target mode on finish)
    this->dynamic_iter_lim = false;
    this->quality = this->iter_lim;

    this->tween_r1 = this->getAngledRect(this->state_a);
    this->tween_r2 = this->getAngledRect(this->state_b);

    // Find the "quad" (angled rect) that fits neatly around the start-quad and end-quad
    DDAngledRect encompassing;
    encompassing.fitTo(this->tween_r1, this->tween_r2, this->tween_r1.aspectRatio());

    // Use that encompassing quad to work out the zoom level when the camera zooms-out
    f128 encompassing_zoom = (camera.viewportWorldSize<f128>() / encompassing.size).average();

    // Convert to "height" (zoom, but on a log scale)
    f128 encompassing_height = std::min(f128(1.0), toHeight(encompassing_zoom));

    // From the current height, how much do we need to "lift" the camera up? (again, on linear scale)
    this->tween_lift = encompassing_height - std::max(
        toHeight(this->state_a.camera.relativeZoom<f128>()),
        toHeight(this->state_b.camera.relativeZoom<f128>())
    );

    // Make sure we don't lift up too high
    f128 max_lift = 1.0 - toHeight(this->state_b.camera.relativeZoom<f128>());
    if (max_lift < 0.0) max_lift = 0;
    if (this->tween_lift < 0) this->tween_lift = 0;
    if (this->tween_lift > max_lift) this->tween_lift = max_lift;

    // Begin tween
    this->tween_progress = 0.0;
    this->tweening = true;
    this->tween_duration = 2;// tweenDistance(scene_data.state_a, scene_data.state_b);
    this->dist_params.cycle_dist_invert = this->state_b.dist_params.cycle_dist_invert;
}

void Mandelbrot_Scene::lerpState(
    const MandelState& a,
    const MandelState& b,
    double f,
    bool complete)
{
    float pos_f = this->tween_pos_spline((float)f);

    // get true src/dst iter_lim (ignore dynamic calculation until tween complete)
    double src_iter_lim = a.dynamic_iter_lim ? (mandelbrotIterLimit(a.camera.relativeZoom<f128>()) * a.quality) : a.quality;
    double dst_iter_lim = b.dynamic_iter_lim ? (mandelbrotIterLimit(b.camera.relativeZoom<f128>()) * b.quality) : b.quality;

    // lerp camera
    double lift_weight = (double)this->tween_zoom_lift_spline((float)f);
    f128 lift_height = this->tween_lift * lift_weight;
    f128 a_height = toHeight(a.camera.relativeZoom<f128>());
    f128 b_height = toHeight(b.camera.relativeZoom<f128>());

    double base_zoom_f = this->tween_base_zoom_spline((float)f);
    f128 dst_height = Math::lerp(a_height, b_height, base_zoom_f) + lift_height;
    CameraInfo::lerp(this->camera, a.camera, b.camera, pos_f);
    this->camera.setRelativeZoom(fromHeight(dst_height)); // override zoom from computed "height"

    // quality
    this->quality = Math::lerp(src_iter_lim, dst_iter_lim, pos_f);

    // color cycle
    float color_cycle_f = this->tween_color_cycle((float)f);

    // shader weights
    this->iter_weight = Math::lerp(a.iter_weight, b.iter_weight, color_cycle_f);
    this->dist_weight = Math::lerp(a.dist_weight, b.dist_weight, color_cycle_f);
    this->stripe_weight = Math::lerp(a.stripe_weight, b.stripe_weight, color_cycle_f);

    // iter shader
    this->iter_params.cycle_iter_value = Math::lerp(a.iter_params.cycle_iter_value, b.iter_params.cycle_iter_value, color_cycle_f);
    this->iter_params.cycle_iter_log1p_weight = Math::lerp(a.iter_params.cycle_iter_log1p_weight, b.iter_params.cycle_iter_log1p_weight, color_cycle_f);

    // dist shader
    this->dist_params.cycle_dist_value = Math::lerp(a.dist_params.cycle_dist_value, b.dist_params.cycle_dist_value, color_cycle_f);
    this->dist_params.cycle_dist_sharpness = Math::lerp(a.dist_params.cycle_dist_sharpness, b.dist_params.cycle_dist_sharpness, color_cycle_f);

    // stripe shader
    this->stripe_params.freq = Math::lerp(a.stripe_params.freq, b.stripe_params.freq, color_cycle_f);
    this->stripe_params.phase = Math::lerp(a.stripe_params.phase, b.stripe_params.phase, color_cycle_f);
    this->stripe_params.contrast = Math::lerp(a.stripe_params.contrast, b.stripe_params.contrast, color_cycle_f);

    // lerp gradient/hue shift
    this->gradient_shift = Math::lerp(a.gradient_shift, b.gradient_shift, pos_f);
    this->hue_shift = Math::lerp(a.hue_shift, b.hue_shift, pos_f);

    // lerp animation speed for gradient/hue shift
    this->gradient_shift_step = Math::lerp(a.gradient_shift_step, b.gradient_shift_step, pos_f);
    this->hue_shift_step = Math::lerp(a.hue_shift_step, b.hue_shift_step, pos_f);

    // lerp color gradient
    ImGradient::lerp(this->gradient, a.gradient, b.gradient, (float)f);

    if (complete)
    {
        this->dynamic_iter_lim = b.dynamic_iter_lim;
        this->quality = b.quality;

        this->iter_params.cycle_iter_dynamic_limit = b.iter_params.cycle_iter_dynamic_limit;

        this->show_axis = b.show_axis;
        this->show_color_animation_options = b.show_color_animation_options;
    }

    // === Lerp quality ===
    //this->iter_lim = Math::lerp(state_a.iter_lim, state_b.iter_lim, f);

    // === Lerp Cardioid Flattening Factor ===
    //this->cardioid_lerp_amount = Math::lerp(state_a.cardioid_lerp_amount, state_b.cardioid_lerp_amount, f);

    // Spline Data
    //memcpy(this->x_spline_point, Math::lerp(state_a.x_spline_points, state_b.x_spline_points, f));
}

SIM_END;
