#pragma once
#include <bitloop.h>
#include <array>
#include <string_view>

#include "../shading/shading.h"
#include "../shading/gradient.h"
#include "../field/escape_field.h"

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

    ImGradient       gradient;
    
    // animate   
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

    static MandelState getDefaults(int version);
};

SIM_END;
