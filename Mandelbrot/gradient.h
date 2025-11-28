#pragma once
#include <bitloop.h>

SIM_BEG;

enum struct GradientPreset
{
    CLASSIC,
    RAINBOW,
    WAVES,
    PHOTON,
    CLOUDS,
    RGB_BANDS,
    ///////////
    COUNT
};

static const char* ColorGradientNames[(int)GradientPreset::COUNT] = {
    "CLASSIC",
    "RAINBOW",
    "WAVES",
    "REFRACTION",
    "CLOUDS",
    "RGB_BANDS"
};


template<GradientPreset type>
inline void colorGradientTemplate(double t, uint8_t& r, uint8_t& g, uint8_t& b)
{
    if constexpr (type == GradientPreset::CLASSIC)
    {
        r = (uint8_t)(9 * (1 - t) * t * t * t * 255);
        g = (uint8_t)(15 * (1 - t) * (1 - t) * t * t * 255);
        b = (uint8_t)(8.5 * (1 - t) * (1 - t) * (1 - t) * t * 255);
    }
    else if constexpr (type == GradientPreset::RAINBOW)
    {
        float a = (float)t * 3.14159265f;
        r = (uint8_t)(sin(a) * sin(a) * 255);
        g = (uint8_t)(sin(a + 2.0944f) * sin(a + 2.0944f) * 255);
        b = (uint8_t)(sin(a + 4.1888f) * sin(a + 4.1888f) * 255);
    }
}

inline void transformGradient(ImGradient& dest, const ImGradient& src, float gradient_shift, float hue_shift)
{
    auto& marks = src.getMarks();
    auto& shifted_marks = dest.getMarks();
    shifted_marks.resize(marks.size());
    for (size_t i = 0; i < marks.size(); i++)
    {
        auto adjusted = Color(marks[i].color).adjustHue((float)hue_shift).vec4();
        memcpy(shifted_marks[i].color, adjusted.data(), sizeof(float) * 4);
        shifted_marks[i].position = Math::wrap(marks[i].position + (float)gradient_shift, 0.0f, 1.0f);
    }
    dest.refreshCache();
}

inline void colorGradientTemplate(GradientPreset type, double t, uint8_t& r, uint8_t& g, uint8_t& b)
{
    switch (type)
    {
    case GradientPreset::CLASSIC: colorGradientTemplate<GradientPreset::CLASSIC>(t, r, g, b); break;
    case GradientPreset::RAINBOW: colorGradientTemplate<GradientPreset::RAINBOW>(t, r, g, b); break;
    default: break;
    }
}

inline void generateGradientFromPreset(
    ImGradient& grad,
    GradientPreset type,
    float hue_threshold = 0.3f,
    float sat_threshold = 0.3f,
    float val_threshold = 0.3f)
{
    grad.getMarks().clear();

    switch (type)
    {
    case GradientPreset::CLASSIC:
    {
        grad.addMark(0.0f, ImColor(0, 0, 0));
        grad.addMark(0.2f, ImColor(39, 39, 214));
        grad.addMark(0.4f, ImColor(0, 143, 255));
        grad.addMark(0.6f, ImColor(255, 255, 68));
        grad.addMark(0.8f, ImColor(255, 30, 0));
    }
    break;

    case GradientPreset::WAVES:
    {
        grad.addMark(0.0f, ImColor(0, 0, 0));
        grad.addMark(0.3f, ImColor(73, 54, 254));
        grad.addMark(0.47f, ImColor(242, 22, 116));
        grad.addMark(0.53f, ImColor(255, 56, 41));
        grad.addMark(0.62f, ImColor(208, 171, 1));
        grad.addMark(0.62001f, ImColor(0, 0, 0));
    }
    break;

    case GradientPreset::PHOTON:
    {
        grad.addMark(0.0f, ImColor(0, 0, 0));
        grad.addMark(0.45f, ImColor(70, 70, 70));
        grad.addMark(0.5f, ImColor(255, 255, 255));
        grad.addMark(0.55f, ImColor(50, 50, 50));
    }
    break;

    case GradientPreset::CLOUDS:
    {
        grad.addMark(0.00f, ImColor(255, 255, 255));
        grad.addMark(0.27f, ImColor(39, 39, 214));
        grad.addMark(0.40f, ImColor(0, 176, 255));
        grad.addMark(0.60f, ImColor(255, 255, 255));
        grad.addMark(0.80f, ImColor(255, 255, 255));
    }
    break;

    case GradientPreset::RGB_BANDS:
    {
        grad.addMark(0.06f, ImColor(0, 0, 0));
        grad.addMark(0.16f, ImColor(0, 255, 252));
        grad.addMark(0.27f, ImColor(0, 0, 0));
		
        grad.addMark(0.40f, ImColor(0, 0, 0));
        grad.addMark(0.50f, ImColor(255, 0, 0));
        grad.addMark(0.61f, ImColor(0, 0, 0));
		
        grad.addMark(0.73f, ImColor(0, 0, 0));
        grad.addMark(0.84f, ImColor(3, 255, 0));
        grad.addMark(0.94f, ImColor(0, 0, 0));
    }
    break;

    default: // Rely on procedurally generated gradient (less ideal as worse marks placement)
    {
        uint8_t last_r, last_g, last_b;
        colorGradientTemplate(type, 0.0f, last_r, last_g, last_b);
        grad.addMark(0.0f, ImColor(last_r, last_g, last_b));

        float last_h, last_s, last_v;
        Color::RGBtoHSV(last_r, last_g, last_b, last_h, last_s, last_v);

        for (float x = 0.0f; x < 1.0f; x += 0.01f)
        {
            uint8_t r, g, b;
            float h, s, v;

            colorGradientTemplate(type, x, r, g, b);
            Color::RGBtoHSV(r, g, b, h, s, v);

            float h_ratio = Math::absAvgRatio(last_h, h);
            float s_ratio = Math::absAvgRatio(last_s, s);
            float v_ratio = Math::absAvgRatio(last_v, v);

            //float avg_ratio = Math::avg(h_ratio, s_ratio, v_ratio);
            if (h_ratio > hue_threshold ||
                s_ratio > sat_threshold ||
                v_ratio > val_threshold)
            {
                grad.addMark(x, ImColor(r, g, b));
                last_h = h;
                last_s = s;
                last_v = v;
            }
        }
    }
    break;
    }
}

SIM_END;
