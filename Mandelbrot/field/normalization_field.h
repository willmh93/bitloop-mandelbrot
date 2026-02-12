#pragma once
#include <bitloop.h>

#include "escape_field.h"

SIM_BEG;

struct NormalizationPixel : public EscapeFieldPixel
{
private:
    union {
        FVec2  world_pos_32;
        DVec2  world_pos_64;
        DDVec2 world_pos_128;
    };
public:

    NormalizationPixel() noexcept
        : EscapeFieldPixel{}
        , world_pos_128{}
        , stage_pos{}
        , is_final(false)
        , weight(0.0f)
    {
    }

    DVec2  stage_pos;
    bool   is_final;
    float  weight;

    template<typename T> requires is_f32<T> [[nodiscard]] constexpr const FVec2& worldPos() { return world_pos_32; }
    template<typename T> requires is_f64<T> [[nodiscard]] constexpr const DVec2& worldPos() { return world_pos_64; }
    template<typename T> requires is_f128<T> [[nodiscard]] constexpr const DDVec2& worldPos() { return world_pos_128; }

    void setWorldPos(FVec2 p) { world_pos_32 = p; }
    void setWorldPos(DVec2 p) { world_pos_64 = p; }
    void setWorldPos(DDVec2 p) { world_pos_128 = p; }
};

class NormalizationField
{
    std::vector<DVec2> local_field;
    WorldObject128 bounds;

public:

    std::vector<NormalizationPixel> world_field;

    NormalizationField()
    {
        setShape(0.025, 2.0);
    }

    void setShape(double sample_r, double exponent)
    {
        local_field = math::delaunayMeshEllipse<f64>(0, 0, 1.0, sample_r, exponent);
        world_field.resize(local_field.size());
    }

    template<typename T>
    void updateField(const CameraInfo& camera, double scale)
    {
        T world_radius = ((T)scale / camera.relativeZoom<T>()) * T{ 2 };
        Vec2<T> cam_center_world = camera.pos<T>();

        // Get stage quad of world ellipse bounds for fast stage interpolation
        WorldObjectT<T> bounds;
        bounds.setCamera(camera);
        bounds.setWorldRect(cam_center_world - world_radius, world_radius * 2);
        DQuad stage_quad = bounds.stageQuad();

        for (size_t i = 0; i < world_field.size(); i++)
        {
            NormalizationPixel& px = world_field[i];
            const DVec2 pt = local_field[i];

            //px.world_pos = cam_center_world + (pt * world_radius);
            px.setWorldPos(cam_center_world + (pt * world_radius));
            px.stage_pos = stage_quad.lerpPoint(0.5 + pt * 0.5);

            px.weight = 1.0f - (float)pt.mag();
        }
    }

    void clearFinalFlags()
    {
        for (NormalizationPixel& p : world_field)
            p.is_final = false;
    }

    template<typename Callback>
    void forEach(Callback&& callback, int thread_count = Thread::threadCount())
    {
        //static_assert(std::is_invocable_r_v<void, Callback, NormalizationPixel&>,
        //    "Callback must be: void(NormalizationPixel&)");

        int pixel_count = (int)world_field.size();
        std::vector<std::pair<int, int>> ranges = Thread::splitRanges<int>(pixel_count, thread_count);

        std::vector<std::future<void>> futures(thread_count);

        for (int ti = 0; ti < thread_count; ti++)
        {
            const std::pair<int, int>& range = ranges[ti];
            futures[ti] = Thread::pool().submit_task([&]()
            {
                int i0 = range.first;
                int i1 = range.second;
                if constexpr (std::is_invocable_r_v<void, Callback, NormalizationPixel&>)
                {
                    for (int i = i0; i < i1; i++)
                        callback(world_field[i]);
                }
                else
                {
                    for (int i = i0; i < i1; i++)
                        callback(world_field[i], ti);
                }
            });
        }

        for (int ti = 0; ti < thread_count; ti++)
        {
            if (futures[ti].valid())
                futures[ti].get();
        }
    }
};

SIM_END;
