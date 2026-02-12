#pragma once
#include <bitloop.h>

#include "../config/build_config.h"
#include "../field/escape_field.h"

SIM_BEG;

using namespace bl;

struct MandelStats
{
    // dirty flags for different stats
    enum struct Dirty : u32 {
        iter_histogram      = 1u << 0,
        dist_histogram      = 1u << 1,
        stripe_histogram    = 1u << 2,
        field_info          = 1u << 3,
        hovered_field_pixel = 1u << 4,
    };

    u32 dirty;
    void sync(Dirty type) { dirty |= (u32)type; }

    struct FieldInfo
    {
        // raw min/max depth (computed directly from field)
        f64 raw_min_depth = 0;
        f64 raw_max_depth = 0;
        f64 raw_min_dist = 0;
        f64 raw_max_dist = 0;
        f32 raw_min_stripe = 0;
        f32 raw_max_stripe = 0;

        // final min/max values after normalization/toning/cycling
        f32 final_min_depth = 0;
        f32 final_max_depth = 0;
        f32 final_min_dist = 0;
        f32 final_max_dist = 0;
        f32 final_min_stripe = 0;
        f32 final_max_stripe = 0;
        
        #if MANDEL_EXTENDED_FIELD_STATS
        ExtendedFieldStats extended;
        #endif
    };

    // feature ==> pixel count histograms
    std::vector<f32> iter_histogram;
    std::vector<f32> dist_histogram;
    std::vector<f32> stripe_histogram;

    // field_info
    FieldInfo field_info;

    // mandelbrot stats for the pixel the mouse is hovered over
    EscapeFieldPixel hovered_field_pixel;
    DDVec2 hovered_field_world_pos;

    // always run sync (but only sync dirty states by flag)
    bool operator==(const MandelStats&) const { return false; }

    MandelStats& operator=(const MandelStats& rhs)
    {
        if (rhs.dirty & (u32)Dirty::iter_histogram)       iter_histogram = rhs.iter_histogram;
        if (rhs.dirty & (u32)Dirty::dist_histogram)       dist_histogram = rhs.dist_histogram;
        if (rhs.dirty & (u32)Dirty::stripe_histogram)     stripe_histogram = rhs.stripe_histogram;

        if (rhs.dirty & (u32)Dirty::field_info)           field_info = rhs.field_info;
        if (rhs.dirty & (u32)Dirty::hovered_field_pixel) { 
            hovered_field_pixel = rhs.hovered_field_pixel;
            hovered_field_world_pos = rhs.hovered_field_world_pos;
        }
        dirty = 0;

        return *this;
    }
};

SIM_END;