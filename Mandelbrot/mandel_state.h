#pragma once
#include <bitloop.h>
#include "shading.h"

SIM_BEG;
using namespace bl;

struct MandelState
{
    /// ─────────────────────── Saveable Info ───────────────────────

    bool show_axis = true;
    CameraInfo camera;

    /// ─────────────────────── Tweenable Info ───────────────────────

    bool        dynamic_iter_lim               = true;
    double      quality                        = 0.5; // Used for UI (ignored during tween, represents iter_lim OR % of iter_lim)
    int         maxdepth_optimize              = (int)MandelMaxDepthOptimization::MEDIUM;
                                               
    bool        use_smoothing                  = true;

    double      iter_weight                    = 1.0;
    double      dist_weight                    = 0.0;
    double      stripe_weight                  = 0.0;
    
    IterParams   iter_params;                              
    DistParams   dist_params;
    StripeParams stripe_params; // todo: Separate StripeComputeParams/StripeShadeParams
                
    double      gradient_shift                 = 0.0;
    double      hue_shift                      = 0.0;
                                               
    double      gradient_shift_step            = 0.0078;
    double      hue_shift_step                 = 0.136;
                
    int         smoothing_type = (int)MandelKernelFeatures::ITER; // this is being dynamically set, no need to save
    int         shade_formula = (int)MandelShaderFormula::ITER_DIST_STRIPE;

    ImGradient  gradient;

    // animate
    bool        show_color_animation_options   = false;

    // Flatten
    bool        flatten = false;
    double      flatten_amount = 0.0;

    ///////////////////////////////////////////////////////////////////////////

    bool operator==(const MandelState&) const = default;

    std::string serialize() const;
    bool deserialize(std::string_view sv)
    {
        // First, try decoding compressed
        if (_deserialize(sv, true)) return true;
        if (_deserialize(sv, false)) return true;
        return true;
    }

private:

    bool _deserialize(std::string_view sv, bool COMPRESS_CONFIG);
};

SIM_END;
