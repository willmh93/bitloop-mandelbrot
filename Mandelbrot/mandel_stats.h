#pragma once
#include <bitloop.h>

SIM_BEG;

struct MandelStats
{
    // dirty flags for different stats
    struct {
        bool depth_histogram = false;      // {1}
        bool stripe_histogram = false;     // {2}
        bool field_info = false;           // {3}
        bool hovered_field_pixel = false;  // {4}
    } dirty;

    struct FieldInfo
    {
        double min_depth,     max_depth;                      // true min/max field depth
        double min_log_depth, max_log_depth;                  // lerped log of true min/max depth
        double assumed_min_depth, assumed_max_depth;          // for current zoom, rough "estimated" min/max depth
        double assumed_log_min_depth, assumed_log_max_depth;  // for current zoom, lerped log of rough "estimated" min/max depth
    };

    // {1} bucket_depth ==> pixel count
    std::vector<int> depth_histogram;
    double depth_histogram_bucket_size;

    // {2} pixel stripe angle ==> pixel count
    //int stripe_histogram[360];
    std::vector<double> stripe_histogram;

    // {3} field_info
    FieldInfo field_info;

    // {4} mandelbrot stats for the pixel the mouse is hovered over
    EscapeFieldPixel hovered_field_pixel;
    DVec2 hovered_field_world_pos;

    // Check if any flags are dirty
    bool operator==(const MandelStats& rhs) const {
        return memcmp(&dirty, &rhs.dirty, sizeof(dirty)) == 0;
    }

    MandelStats& operator=(const MandelStats& rhs)
    {
        if (rhs.dirty.depth_histogram)      depth_histogram = rhs.depth_histogram;
        //if (rhs.dirty.stripe_histogram)     memcpy(stripe_histogram, rhs.stripe_histogram, sizeof(stripe_histogram));
        if (rhs.dirty.stripe_histogram)     stripe_histogram = rhs.stripe_histogram;
        if (rhs.dirty.field_info)           field_info = rhs.field_info;
        if (rhs.dirty.hovered_field_pixel)  hovered_field_pixel = rhs.hovered_field_pixel;

        return *this;
    }
};


// todo: operator== memcmp really be needed, but returning 'false' every time causes a freeze which shouldn't happen.
//       This must be because syncing "everything else" is paused as long as "something is changing" which
//       really points to an issue with your VarBuffer solution. One variable should be able to sync
//       indefinitely both ways without freezing the whole app. You did do it that way for a reason though...
//       predictable cause and affect. e.g. Syncing isn't back-forth-back-forth, UI runs more often than Live,
//       so when the live buffer gets changed, the UI can just override it before it has a chance to sync.
//
//        if (!shadow_changed)
//            pushDataToShadow(); // If shadow is *always* changing, no live data EVER gets send to it. Maybe
//                                   experiment with doing this per-variable again (which wouldn't fix this issue,
//                                   but would prevent the bug from leaking to other variables). A real fix would
//                                   be to somehow fix the order of syncs when a render or worker frame is about
//                                   to occur
//

SIM_END;