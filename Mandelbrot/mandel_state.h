#pragma once
#include <bitloop.h>
#include "shading.h"
#include "mandel_field.h"

SIM_BEG;
using namespace bl;

struct MandelState
{
    /// ─────────────────────── Saveable Info ───────────────────────

    bool show_axis = true;
    CameraInfo camera;

    /// ─────────────────────── Tweenable Info ───────────────────────

    bool         dynamic_iter_lim               = true;
    double       quality                        = 0.5; // Used for UI (ignored during tween, represents iter_lim OR % of iter_lim)
    int          interior_forwarding              = (int)MandelInteriorForwarding::MEDIUM;
                                                
    bool         use_smoothing                  = true;
                 
    double       iter_weight                    = 1.0;
    double       dist_weight                    = 0.0;
    double       stripe_weight                  = 0.0;
    
    IterParams   iter_params;                              
    DistParams   dist_params;
    StripeParams stripe_params;
    
    PivotToneParams dist_tone_params;
    PivotToneParams stripe_tone_params;

    double       gradient_shift                 = 0.0;
    double       hue_shift                      = 0.0;
                                                
    double       gradient_shift_step            = 0.0078;
    double       hue_shift_step                 = 0.136;
               
    int          shade_formula = (int)MandelShaderFormula::ITER_DIST_STRIPE;

    ImGradient   gradient;
                 
    // animate   
    bool         show_color_animation_options   = false;
                 
    // Flatten   
    bool         flatten = false;
    double       flatten_amount = 0.0;

    /// ─────────────────────── Normalization Info ───────────────────────

    double normalize_field_scale = 2.0;
    double normalize_field_quality = 0.1;
    double normalize_field_exponent = 1.0;
    double normalize_field_precision = 0.0; // 0=never recalculate when overlapping nearest existing pixel, 1=always use EXACT sample coordinate

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
};

SIM_END;
