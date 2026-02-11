#pragma once
#include <bitloop.h>
#include <array>
#include <string_view>
#include "shading.h"
#include "gradient.h"
#include "mandel_field.h"

SIM_BEG;
using namespace bl;


struct MandelState
{
    // brotli dictionaries (versioned)
    static std::span<const std::string_view> getDictionaryTokens(int v);
    static const compression::BrotliDict*    getDictionary(int v);

    /// ─────────────────────── Saveable Info ───────────────────────

    CameraInfo       camera;
                     
    bool             show_axis                      = true;
    bool             dynamic_iter_lim               = true;
    double           quality                        = 0.5; // Used for UI (ignored during tween, represents iter_lim OR % of iter_lim)
    int              interior_forwarding            = (int)MandelInteriorForwarding::MEDIUM;
                                                    
    bool             use_smoothing                  = true;
                     
    double           iter_weight                    = 1.0;
    double           dist_weight                    = 0.0;
    double           stripe_weight                  = 0.0;
                     
    float            iter_x_dist_weight             = 0.0f; // (iter * dist) + stripe
    float            dist_x_stripe_weight           = 0.0f; // (dist * stripe) + iter
    float            stripe_x_iter_weight           = 0.0f; // (stripe * iter) + dist
    float            iter_x_distStripe_weight       = 0.0f; // iter * (dist + stripe)
    float            dist_x_iterStripe_weight       = 0.0f; // dist * (iter + stripe)
    float            stripe_x_iterDist_weight       = 0.0f; // stripe * (iter + dist)

    std::string shader_source_txt =
        "vec4 userShade(float iter, float dist, float stripe, vec2 uv)\n{\n"
        "    return sampleGradient(wrap01(iter + dist + stripe));\n"
        "}";
                     
    IterParams       iter_params;                              
    DistParams       dist_params;
    StripeParams     stripe_params;
    
    PivotToneParams  dist_tone_params;
    PivotToneParams  stripe_tone_params;

    double           gradient_shift                 = 0.0;
    double           hue_shift                      = 0.0;
                                                    
    double           gradient_shift_step            = 0.0078;
    double           hue_shift_step                 = 0.136;
    float            phase_step                     = 0.0f;
               
    //int            shade_formula = (int)MandelShaderFormula::ITER_DIST_STRIPE;

    ImGradient       gradient;
    
    // animate   
    //bool           show_color_animation_options = false;
    bool             animate_gradient_shift = false;
    bool             animate_gradient_hue = false;
    bool             animate_stripe_phase = false;
    
    // Flatten   
    bool         flatten = false;
    double       flatten_amount = 0.0;

    /// ─────────────────────── Normalization Info ───────────────────────

    double normalize_field_scale = 2.0;
    double normalize_field_quality = 0.1;
    double normalize_field_exponent = 1.0;
    double normalize_field_precision = 0.0; // 0=never recalculate when overlapping nearest existing pixel, 1=always use EXACT sample coordinate

    SnapshotPresetHashMap valid_presets;

    ///////////////////////////////////////////////////////////////////////////

    bool operator==(const MandelState&) const = default;

    std::string serialize() const;
    bool deserialize(std::string_view sv)
    {
        // first, try decoding compressed
        if (_deserialize(sv, true)) return true;
        if (_deserialize(sv, false)) return true;
        return true;
    }

private:

    bool _deserialize(std::string_view sv, bool COMPRESS_CONFIG);

    static MandelState getDefaults(int version)
    {
        MandelState s;

        if (version == 0)
        {
            s.interior_forwarding = (int)MandelInteriorForwarding::MEDIUM;

            s.iter_weight = 1.0;
            s.dist_weight = 0.0;
            s.stripe_weight = 0.0;

            s.iter_x_dist_weight = 0.0f;
            s.dist_x_stripe_weight = 0.0f;
            s.stripe_x_iter_weight = 0.0f;
            s.iter_x_distStripe_weight = 0.0f;
            s.dist_x_iterStripe_weight = 0.0f;
            s.stripe_x_iterDist_weight = 0.0f;

            s.shader_source_txt = 
                "// @pass base\n"
                "return sampleGradient(wrap01(iter + dist + stripe));\n";

            s.iter_params.cycle_iter_dynamic_limit = false;
            s.iter_params.cycle_iter_normalize_depth = false;
            s.iter_params.cycle_iter_log1p_weight = false;
            s.iter_params.cycle_iter_value = 1.0f;

            s.dist_params.cycle_dist_invert = false;
            s.dist_params.cycle_dist_value = 0.25;
            s.dist_params.cycle_dist_sharpness = 0.9;

            s.stripe_params.freq = 3;
            s.stripe_params.phase = 0.0;

            s.dist_tone_params.brightness = 0.0f;
            //s.dist_tone_params.contrast = 1.0f; // dist has no contrast
            s.dist_tone_params.gamma = 1.0f;

            s.stripe_tone_params.brightness = 0.0f;
            s.stripe_tone_params.contrast = 1.0f;
            s.stripe_tone_params.gamma = 1.0f;

            s.gradient_shift = 0.0;
            s.hue_shift = 0.0;

            s.gradient_shift_step = 0.0078;
            s.hue_shift_step = 0.136;
            s.phase_step = math::toRadians(1.0f);

            generateGradientFromPreset(s.gradient, GradientPreset::CLASSIC);
            //s.show_color_animation_options = false;
            s.animate_gradient_shift = false;
            s.animate_gradient_hue = false;
            s.animate_stripe_phase = false;

            s.normalize_field_scale = 2.0;
            s.normalize_field_quality = 0.1;
            s.normalize_field_exponent = 1.0;
            s.normalize_field_precision = 0.0;
        }

        return s;
    }
};

SIM_END;
