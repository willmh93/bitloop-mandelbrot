#include "mandel_state.h"

SIM_BEG;

std::string MandelState::serialize() const
{
    constexpr bool COMPRESS_CONFIG = true;
    constexpr int defaults_version = 0;

    MandelState defaults = getDefaults(defaults_version);

    // Increment version each time the format changes
    uint32_t version = 1;
    uint32_t flags = 0;

    if (dynamic_iter_lim)                       flags |= MANDEL_DYNAMIC_ITERS;
    if (show_axis)                              flags |= MANDEL_SHOW_AXIS;
    if (flatten)                                flags |= MANDEL_FLATTEN;
    if (iter_params.cycle_iter_dynamic_limit)   flags |= MANDEL_DYNAMIC_COLOR_CYCLE;
    if (iter_params.cycle_iter_normalize_depth) flags |= MANDEL_NORMALIZE_DEPTH;
    if (dist_params.cycle_dist_invert)          flags |= MANDEL_INVERT_DIST;
    if (use_smoothing)                          flags |= MANDEL_USE_SMOOTHING;
    flags |= (version << MANDEL_VERSION_BITSHIFT);

    JSON::json info;
    int decimals = camera.getPositionDecimals();

    // meta
    info["_"] = defaults_version;
    info["f"] = flags;

    if (version >= 0)
    {
        /// Mark JSON floats with "CLEANFLOAT()" to remove trailing zeros/rounding errors when serializing

        // View
        info["x"] = to_string(camera.x<f128>(), decimals);
        info["y"] = to_string(camera.y<f128>(), decimals);
        info["z"] = to_string(camera.relativeZoom<f128>(), 5, false, true, true);
        if (camera.stretchX() != 1.0) info["a"] = JSON::markCleanFloat(camera.stretchX(), 3);
        if (camera.stretchY() != 1.0) info["b"] = JSON::markCleanFloat(camera.stretchY(), 3);
        if (camera.rotation() != 0.0) info["r"] = JSON::markCleanFloat(math::toDegrees(camera.rotation()), 0);

        // Quality
        info["q"] = JSON::markCleanFloat(quality, 3);
        if (interior_forwarding != defaults.interior_forwarding) info["Q"] = interior_forwarding;

        // Color cycle
        {
            /// --- weights ---
            if (iter_weight != defaults.iter_weight)     info["m"] = JSON::markCleanFloat(iter_weight, 2);
            if (dist_weight != defaults.dist_weight)     info["n"] = JSON::markCleanFloat(dist_weight, 2);
            if (stripe_weight != defaults.stripe_weight) info["o"] = JSON::markCleanFloat(stripe_weight, 2);

            if (iter_x_dist_weight       != defaults.iter_x_dist_weight)       info["C"] = JSON::markCleanFloat(iter_x_dist_weight, 3);
            if (dist_x_stripe_weight     != defaults.dist_x_stripe_weight)     info["D"] = JSON::markCleanFloat(dist_x_stripe_weight, 3);
            if (stripe_x_iter_weight     != defaults.stripe_x_iter_weight)     info["F"] = JSON::markCleanFloat(stripe_x_iter_weight, 3);
            if (iter_x_distStripe_weight != defaults.iter_x_distStripe_weight) info["I"] = JSON::markCleanFloat(iter_x_distStripe_weight, 3);
            if (dist_x_iterStripe_weight != defaults.dist_x_iterStripe_weight) info["J"] = JSON::markCleanFloat(dist_x_iterStripe_weight, 3);
            if (stripe_x_iterDist_weight != defaults.stripe_x_iterDist_weight) info["L"] = JSON::markCleanFloat(stripe_x_iterDist_weight, 3);

            /// --- features ---

            // ITER
            if (iter_params.cycle_iter_value        != defaults.iter_params.cycle_iter_value)        info["i"] = JSON::markCleanFloat(iter_params.cycle_iter_value);
            if (iter_params.cycle_iter_log1p_weight != defaults.iter_params.cycle_iter_log1p_weight) info["l"] = JSON::markCleanFloat(iter_params.cycle_iter_log1p_weight);

            // DIST
            if (dist_params.cycle_dist_value        != defaults.dist_params.cycle_dist_value)        info["d"] = JSON::markCleanFloat(dist_params.cycle_dist_value, 5);
            if (dist_params.cycle_dist_sharpness    != defaults.dist_params.cycle_dist_sharpness)    info["s"] = JSON::markCleanFloat(dist_params.cycle_dist_sharpness, 5);

            // STRIPE
            if (stripe_params.freq                  != defaults.stripe_params.freq)                  info["v"] = stripe_params.freq;
            if (stripe_params.phase                 != defaults.stripe_params.phase)                 info["j"] = JSON::markCleanFloat(stripe_params.phase, 4);

            /// --- tone ---

            // DIST tone
            if (dist_tone_params.brightness != defaults.dist_tone_params.brightness)     info["E"] = JSON::markCleanFloat(dist_tone_params.brightness, 3);
            if (dist_tone_params.gamma      != defaults.dist_tone_params.gamma)          info["K"] = JSON::markCleanFloat(dist_tone_params.gamma, 3);

            // STRIPE tone
            if (stripe_tone_params.contrast   != defaults.stripe_tone_params.contrast)   info["c"] = JSON::markCleanFloat(stripe_tone_params.contrast, 3);
            if (stripe_tone_params.brightness != defaults.stripe_tone_params.brightness) info["e"] = JSON::markCleanFloat(stripe_tone_params.brightness, 3);
            if (stripe_tone_params.gamma      != defaults.stripe_tone_params.gamma)      info["k"] = JSON::markCleanFloat(stripe_tone_params.gamma, 3);
        }

        // Shift
        if (gradient_shift != defaults.gradient_shift) info["g"] = JSON::markCleanFloat(gradient_shift);
        if (hue_shift != defaults.hue_shift)           info["h"] = JSON::markCleanFloat(hue_shift);

        // Shift increment
        if (show_color_animation_options != defaults.show_color_animation_options) info["A"] = show_color_animation_options ? 1 : 0;
        if (gradient_shift_step          != defaults.gradient_shift_step)          info["G"] = JSON::markCleanFloat(gradient_shift_step);
        if (hue_shift_step               != defaults.hue_shift_step)               info["H"] = JSON::markCleanFloat(hue_shift_step);

        // Gradient
        if (gradient != defaults.gradient)  info["p"] = gradient.serialize();

        if (normalize_field_scale     != defaults.normalize_field_scale)     info["M"] = JSON::markCleanFloat(normalize_field_scale);
        if (normalize_field_quality   != defaults.normalize_field_quality)   info["N"] = JSON::markCleanFloat(normalize_field_quality);
        if (normalize_field_exponent  != defaults.normalize_field_exponent)  info["O"] = JSON::markCleanFloat(normalize_field_exponent);
        if (normalize_field_precision != defaults.normalize_field_precision) info["P"] = JSON::markCleanFloat(normalize_field_precision);

        // Capture filters
        std::string concatenated_filters_b64;
        auto& all_presets = main_window()->getSnapshotPresetManager()->allPresets();
        for (const auto [hash, enabled] : valid_presets)
        {
            auto *p = all_presets.find(hash);
            concatenated_filters_b64 += p->getAlias();
            concatenated_filters_b64 += ";";
        }
        if (concatenated_filters_b64 != "")
            info["R"] = concatenated_filters_b64;
    }
    else if (version >= 1)
    {
        //info["A"] = x_spline.serialize(SplineSerializationMode::COMPRESS_SHORTEST);
        //info["B"] = y_spline.serialize(SplineSerializationMode::COMPRESS_SHORTEST);
    }

    /// Generate JSON with entries containing "CLEANFLOAT(...)" strings, then strip it away leaving just the raw floats
    /// Guaranteed to contain no rounding errors
    std::string raw_json = JSON::unmarkCleanFloats(info.dump());

    if constexpr (COMPRESS_CONFIG)
    {
        // Post-process: Remove CLEANFLOAT() markers, compress & wrap lines
        std::string json_unquoted = JSON::json_remove_key_quotes(raw_json);

        blPrint() << "Saving: " << json_unquoted;
        std::string compressed_unquoted = compression::brotli_ascii_compress(json_unquoted);

        return compressed_unquoted;
    }
    else
    {
        blPrint() << "Saving: " << raw_json;
        return raw_json;
    }
}

bool MandelState::_deserialize(std::string_view sv, bool COMPRESS_CONFIG)
{
    sv = text::trimStringView(sv);
    std::string txt = sv.data();

    std::string uncompressed;
    if (COMPRESS_CONFIG)
    {
        if (!compression::valid_b62(txt)) return false;
        uncompressed = compression::brotli_ascii_decompress(text::unwrapString(txt));
    }
    else
    {
        uncompressed = txt;
    }

    blPrint() << "uncleaned: " << uncompressed;
    uncompressed = JSON::json_add_key_quotes(uncompressed);
    uncompressed = JSON::json_add_leading_zeros(uncompressed);
    blPrint() << "decoded: " << uncompressed;

    nlohmann::json info = nlohmann::json::parse(uncompressed, nullptr, false);
    if (info.is_discarded())
        return false; // Failed to parse

    uint32_t flags = info.value("f", 0ul);
    uint32_t version = (flags & MANDEL_VERSION_MASK) >> MANDEL_VERSION_BITSHIFT;
    int defaults_version = info.value("_", 0);

    MandelState defaults = getDefaults(defaults_version);

    if (version >= 0)
    {
        //mandel_features = ((flags & MANDEL_SMOOTH_MASK) >> MANDEL_SMOOTH_BITSHIFT);

        dynamic_iter_lim                       = flags & MANDEL_DYNAMIC_ITERS;
        show_axis                              = flags & MANDEL_SHOW_AXIS;
        flatten                                = flags & MANDEL_FLATTEN;
        iter_params.cycle_iter_dynamic_limit   = flags & MANDEL_DYNAMIC_COLOR_CYCLE;
        iter_params.cycle_iter_normalize_depth = flags & MANDEL_NORMALIZE_DEPTH;
        dist_params.cycle_dist_invert          = flags & MANDEL_INVERT_DIST;
        use_smoothing                          = flags & MANDEL_USE_SMOOTHING;

        // View
        if (info.contains("x") && info.find("x").value().is_string())
            camera.setX(from_string(info.value("x", "0").c_str())); // f128 path (deserialize from string)
        else
            camera.setX(info.value("x", 0.0)); // old double path (raw)

        if (info.contains("y") && info.find("y").value().is_string())
            camera.setY(from_string(info.value("y", "0").c_str())); // f128 path (deserialize from string)
        else
            camera.setY(info.value("y", 0.0)); // old double path (raw)

        if (info.contains("z") && info.find("z").value().is_string())
            camera.setRelativeZoom(from_string(info.value("z", "0").c_str())); // f128 path (deserialize from string)
        else
            camera.setRelativeZoom(info.value("z", 0.0)); // old double path (raw)

        camera.setStretchX(info.value("a", 1.0));
        camera.setStretchY(info.value("b", 1.0));
        camera.setRotation(math::toRadians(info.value("r", 0.0)));

        // Quality
        quality = info.value("q", quality);

        // Color cycle
        {
            // Mix ratio (v0 - recover new weights from old format)
            if (version == 0)
            {
                /// TODO: deprecate
                iter_weight   = 1.0 - info.value("m", defaults.iter_weight); 
                dist_weight   = info.value("n", 1.0 - defaults.iter_weight);
                stripe_weight = info.value("o", 0.0);
            }
            else if (version >= 1)
            {
                iter_weight   = info.value("m", defaults.iter_weight);
                dist_weight   = info.value("n", defaults.dist_weight);
                stripe_weight = info.value("o", defaults.stripe_weight);
            }

            iter_x_dist_weight       = info.value("C", defaults.iter_x_dist_weight);
            dist_x_stripe_weight     = info.value("D", defaults.dist_x_stripe_weight);
            stripe_x_iter_weight     = info.value("F", defaults.stripe_x_iter_weight);
            iter_x_distStripe_weight = info.value("I", defaults.iter_x_distStripe_weight);
            dist_x_iterStripe_weight = info.value("J", defaults.dist_x_iterStripe_weight);
            stripe_x_iterDist_weight = info.value("L", defaults.stripe_x_iterDist_weight);

            // ITER
            iter_params.cycle_iter_value        = info.value("i", defaults.iter_params.cycle_iter_value);
            iter_params.cycle_iter_log1p_weight = info.value("l", defaults.iter_params.cycle_iter_log1p_weight);

            // DIST
            dist_params.cycle_dist_value     = info.value("d", defaults.dist_params.cycle_dist_value);
            dist_params.cycle_dist_sharpness = info.value("s", defaults.dist_params.cycle_dist_sharpness);

            // STRIPE
            stripe_params.freq  = info.value("v", defaults.stripe_params.freq);
            stripe_params.phase = info.value("j", defaults.stripe_params.phase);

            // DIST tone
            dist_tone_params.brightness = info.value("E", defaults.dist_tone_params.brightness);
            dist_tone_params.gamma      = info.value("K", defaults.dist_tone_params.gamma);

            // STRIPE tone
            stripe_tone_params.contrast   = info.value("c", defaults.stripe_tone_params.contrast);
            stripe_tone_params.brightness = info.value("e", defaults.stripe_tone_params.brightness);
            stripe_tone_params.gamma      = info.value("k", defaults.stripe_tone_params.gamma);
        }

        // Shift
        show_color_animation_options = info.value("A", defaults.show_color_animation_options ? 1 : 0) != 0;
        gradient_shift      = info.value("g", defaults.gradient_shift);
        hue_shift           = info.value("h", defaults.hue_shift);

        // Shift increment
        gradient_shift_step = info.value("G", defaults.gradient_shift_step);
        hue_shift_step      = info.value("H", defaults.hue_shift_step);

        // Gradient
        if (info.contains("p"))
            gradient.deserialize(info.value("p", ""));
        else
            gradient = defaults.gradient;

        // Normalization
        normalize_field_scale     = info.value("M", defaults.normalize_field_scale);
        normalize_field_quality   = info.value("N", defaults.normalize_field_quality);
        normalize_field_exponent  = info.value("O", defaults.normalize_field_exponent);
        normalize_field_precision = info.value("P", defaults.normalize_field_precision);

        // Capture filters
        std::string          concat_aliases = info.value("R", std::string(""));
        SnapshotPresetList&  all_presets    = main_window()->getSnapshotPresetManager()->allPresets();
        string_view_list     aliases        = text::split(concat_aliases, ';', true);

        valid_presets.clear();
        for (const auto& alias : aliases)
        {
            SnapshotPreset* preset = all_presets.findByAlias(alias);
            if (preset) valid_presets[preset->hashedAlias()] = true;
        }

        //std::string concatenated_filters_b64 = info.value("R", std::string(""));
        //valid_presets.clear();
        //for (int i = 0; i < concatenated_filters_b64.size(); i += 11)
        //{
        //    bl::hash_t hash = compression::b64_decode_u64(concatenated_filters_b64.substr(i, 11));
        //    valid_presets[hash] = true;
        //}
    }

    if (version >= 1)
    {
        //if (info.contains("A")) x_spline.deserialize(info["A"].get<std::string>());
        //if (info.contains("B")) y_spline.deserialize(info["B"].get<std::string>());
    }

    return true;
}


SIM_END;