#include "mandel_state.h"

SIM_BEG;

std::string MandelState::serialize() const
{
    constexpr bool COMPRESS_CONFIG = true;

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

    //flags |= ((uint32_t)mandel_features << MANDEL_SMOOTH_BITSHIFT);
    flags |= (version << MANDEL_VERSION_BITSHIFT);

    JSON::json info;
    int decimals = camera.getPositionDecimals();

    if (version >= 0)
    {
        info["f"] = flags;

        /// Mark JSON floats with "CLEANFLOAT()" to remove trailing zeros/rounding errors when serializing

        // View
        info["x"] = to_string(camera.x<f128>(), decimals);
        info["y"] = to_string(camera.y<f128>(), decimals);
        info["z"] = to_string(camera.relativeZoom<f128>(), 32);
        info["a"] = JSON::markCleanFloat(camera.stretchX(), 3);
        info["b"] = JSON::markCleanFloat(camera.stretchY(), 3);
        info["r"] = JSON::markCleanFloat(Math::toDegrees(camera.rotation()), 0);

        // Quality
        info["q"] = JSON::markCleanFloat(quality, 3);

        // Color cycle
        {
            // Mix ratio
            info["m"] = JSON::markCleanFloat(iter_weight, 2);
            info["n"] = JSON::markCleanFloat(dist_weight, 2);
            info["o"] = JSON::markCleanFloat(stripe_weight, 2);

            // ITER
            info["i"] = JSON::markCleanFloat(iter_params.cycle_iter_value);
            info["l"] = JSON::markCleanFloat(iter_params.cycle_iter_log1p_weight);

            // DIST
            info["d"] = JSON::markCleanFloat(dist_params.cycle_dist_value, 5);
            info["s"] = JSON::markCleanFloat(dist_params.cycle_dist_sharpness, 5);

            // STRIPE
            info["v"] = (int)stripe_params.freq;
            info["j"] = JSON::markCleanFloat(stripe_params.phase, 4);
            info["c"] = JSON::markCleanFloat(stripe_params.contrast, 3);
        }

        // Shift
        info["g"] = JSON::markCleanFloat(gradient_shift);
        info["h"] = JSON::markCleanFloat(hue_shift);

        // Shift increment
        info["A"] = show_color_animation_options ? 1 : 0;
        info["G"] = JSON::markCleanFloat(gradient_shift_step);
        info["H"] = JSON::markCleanFloat(hue_shift_step);

        // Gradient
        info["p"] = gradient.serialize();
    }
    else if (version >= 1)
    {
        //info["A"] = x_spline.serialize(SplineSerializationMode::COMPRESS_SHORTEST);
        //info["B"] = y_spline.serialize(SplineSerializationMode::COMPRESS_SHORTEST);
    }

    std::string raw_json = JSON::unmarkCleanFloats(info.dump());

    if constexpr (COMPRESS_CONFIG)
    {
        // Post-process: Remove CLEANFLOAT() markers, compress & wrap lines
        std::string json_unquoted = JSON::json_remove_key_quotes(raw_json);
        std::string compressed_unquoted = Compression::brotli_ascii_compress(json_unquoted);
        return compressed_unquoted;
    }
    else
    {
        return raw_json;
    }
}

bool MandelState::_deserialize(std::string_view sv, bool COMPRESS_CONFIG)
{
    sv = TextUtil::trim_view(sv);
    std::string txt = sv.data();

    std::string uncompressed;
    if (COMPRESS_CONFIG)
    {
        if (!Compression::valid_b62(txt)) return false;
        uncompressed = Compression::brotli_ascii_decompress(TextUtil::unwrapString(txt));
    }
    else
    {
        uncompressed = txt;
    }

    uncompressed = JSON::json_add_key_quotes(uncompressed);
    uncompressed = JSON::json_add_leading_zeros(uncompressed);

    blPrint() << "decoded: " << uncompressed;

    nlohmann::json info = nlohmann::json::parse(uncompressed, nullptr, false);
    if (info.is_discarded())
        return false; // Failed to parse

    uint32_t flags = info.value("f", 0ul);
    uint32_t version = (flags & MANDEL_VERSION_MASK) >> MANDEL_VERSION_BITSHIFT;

    if (version >= 0)
    {
        //mandel_features = ((flags & MANDEL_SMOOTH_MASK) >> MANDEL_SMOOTH_BITSHIFT);

        dynamic_iter_lim = flags & MANDEL_DYNAMIC_ITERS;
        show_axis = flags & MANDEL_SHOW_AXIS;
        flatten = flags & MANDEL_FLATTEN;
        iter_params.cycle_iter_dynamic_limit = flags & MANDEL_DYNAMIC_COLOR_CYCLE;
        iter_params.cycle_iter_normalize_depth = flags & MANDEL_NORMALIZE_DEPTH;
        dist_params.cycle_dist_invert = flags & MANDEL_INVERT_DIST;
        use_smoothing = flags & MANDEL_USE_SMOOTHING;

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
        camera.setRotation(Math::toRadians(info.value("r", 0.0)));

        // Quality
        quality = info.value("q", quality);

        // Color cycle
        {
            // Mix ratio
            //iter_dist_mix = info.value("m", iter_dist_mix);

            // old iter_dist_mix implied 0=full iter, 1=full dist


            iter_weight = 1.0 - info.value("m", iter_weight); // recover new weights from old format
            dist_weight = info.value("n", 1.0 - iter_weight);
            stripe_weight = info.value("o", 0.0);

            // ITER
            iter_params.cycle_iter_value = info.value("i", iter_params.cycle_iter_value);
            iter_params.cycle_iter_log1p_weight = info.value("l", iter_params.cycle_iter_log1p_weight);

            // DIST
            dist_params.cycle_dist_value = info.value("d", dist_params.cycle_dist_value);
            dist_params.cycle_dist_sharpness = info.value("s", dist_params.cycle_dist_sharpness);

            // STRIPE
            stripe_params.freq = info.value("v", 0.0f);
            stripe_params.phase = info.value("j", 0.0f);
            stripe_params.contrast = info.value("c", 3.0f);
        }

        // Shift
        show_color_animation_options = info.value("A", show_color_animation_options ? 1 : 0) != 0;
        gradient_shift = info.value("g", 0.0);
        hue_shift = info.value("h", 0.0);

        // Shift increment
        gradient_shift_step = info.value("G", gradient_shift_step);
        hue_shift_step = info.value("H", hue_shift_step);

        // Gradient
        if (info.contains("p"))
            gradient.deserialize(info.value("p", ""));
    }

    if (version >= 1)
    {
        iter_weight = info.value("m", iter_weight);
        dist_weight = info.value("n", dist_weight);
        stripe_weight = info.value("o", stripe_weight);

        //if (info.contains("A")) x_spline.deserialize(info["A"].get<std::string>());
        //if (info.contains("B")) y_spline.deserialize(info["B"].get<std::string>());
    }

    return true;
}


SIM_END;