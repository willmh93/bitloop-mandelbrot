#pragma once
#include <bitloop.h>
#include <complex>

SIM_BEG;

using namespace bl;


void plot(const SceneBase* scene, Viewport* ctx, bool interactive, int segments = 720, double ox=0);
void animatePlot(Viewport* ctx, double scale, double ox, double orig_angle, double dist=0.2);

inline double cardioidSquaredDistance(double angle, double mx, double my)
{
    double cos_a = cos(angle);
    double sin_a = sin(angle);
    double cos_2a = 2.0 * cos_a * cos_a - 1.0;
    double sin_2a = 2.0 * sin_a * cos_a;

    double px = 0.5 * cos_a - 0.25 * cos_2a;
    double py = 0.5 * sin_a - 0.25 * sin_2a;
    double dx = mx - px;
    double dy = my - py;
    return dx * dx + dy * dy;
}

inline double originalAngleBinarySearch(double mx, double my)
{
    double left = 0.0;
    double right = 2.0 * math::pi;
    const double eps = 0.02;// 1e-13;

    while (right - left > eps)
    {
        double m1 = left + (right - left) / 3.0;
        double m2 = right - (right - left) / 3.0;

        if (cardioidSquaredDistance(m1, mx, my) < cardioidSquaredDistance(m2, mx, my))
            right = m2;
        else
            left = m1;
    }

    return (left + right) / 2.0;
}

// Newton's Method to find theta that minimizes distance^2 from (px,py)
inline double originalAngleFromPoint(
    double px, 
    double py, 
    double tolerance = 1e-10, 
    int maxIter = 50)
{
    // Check if inside unstable region for Newton's method
    double y_bound = 0.001 + 4 * pow(px - 0.23, 2.0);
    if (px >= 0.25 && py < y_bound && py > -y_bound)
        return originalAngleBinarySearch(px, py);

    double initialGuess = atan2(py, px - 0.25);
    double theta = initialGuess;

    for (int i = 0; i < maxIter; ++i)
    {
        double _cos = cos(theta);
        double _sin = sin(theta);
        double _cos2 = cos(2.0 * theta);
        double _sin2 = sin(2.0 * theta);

        // Compute Delta x and Delta y
        double X = (0.5 * _cos - 0.25 * _cos2) - px;
        double Y = (0.5 * _sin - 0.25 * _sin2) - py;

        // First derivative f'(theta)
        double d1 = 2.0 * (
            X * (-0.5 * _sin + 0.5 * _sin2)
            + Y * (0.5 * _cos - 0.5 * _cos2)
            );

        // Second derivative f''(theta)
        double d2 = 2.0 * (((-0.5 * _sin + 0.5 * _sin2) * (-0.5 * _sin + 0.5 * _sin2)) + X * (-0.5 * _cos + _cos2))
            + 2.0 * (((0.5 * _cos - 0.5 * _cos2) * (0.5 * _cos - 0.5 * _cos2)) + Y * (-0.5 * _sin + _sin2));

        // If denominator is near zero, break (avoid divide by zero)
        if (fabs(d2) < 1e-14) break;

        // Newton update
        double step = d1 / d2;
        theta -= step;

        // Wrap angle
        if (theta < 0)
            theta += 2.0 * math::pi;
        else if (theta >= 2.0 * math::pi)
            theta -= 2.0 * math::pi;

        // If the step is sufficiently small, we consider we've converged
        if (fabs(step) < tolerance)
            break;
    }

    return theta;
}

inline void cardioidPolarCoord(
    double px,
    double py,
    double& tangent_x,
    double& tangent_y,
    double& tangent_angle,
    double& dist,
    double& orig_angle)
{
    if (px * px + py * py < 0.0625)
    {
        tangent_x = 0.25;
        tangent_y = 0;
        orig_angle = 0;
        tangent_angle = 0;
        dist = -1e-13;
        return;
    }

    orig_angle = math::wrapRadians2PI(originalAngleFromPoint(px, py));
    tangent_angle = math::wrapRadians2PI(1.5 * orig_angle);
    tangent_x = 0.5 * cos(orig_angle) - 0.25 * cos(orig_angle * 2.0);
    tangent_y = 0.5 * sin(orig_angle) - 0.25 * sin(orig_angle * 2.0);

    double dx = px - tangent_x;
    double dy = py - tangent_y;

    //double x1 = dx * cos(tangent_angle) + dy * sin(tangent_angle);
    double y1 = dy * cos(tangent_angle) - dx * sin(tangent_angle);

    dist = -y1;

    if (dist < 0)
    {
        tangent_x = 0.25;
        tangent_y = 0;
        orig_angle = 0;
        tangent_angle = 0;
        dist = -1e-13;
    }
}


inline DVec2 fromPolarCoordinate(double angle, double dist)
{
    //double angle = (perp_angle + M_PI / 2.0) / 1.5;
    double tangent_angle = 1.5 * angle;
    double perp_angle = tangent_angle - math::pi / 2.0;
    double tx = 0.5 * cos(angle) - 0.25 * cos(angle * 2.0);
    double ty = 0.5 * sin(angle) - 0.25 * sin(angle * 2.0);
    double px = tx + cos(perp_angle) * dist;
    double py = ty + sin(perp_angle) * dist;
    return { px, py };
}


double originalAngleFromPerpAngle(
    double perp_angle);

struct CardioidStepInfo
{
    double step_angle;
    double step_dist;
};

struct CardioidSegment : public DVec2
{
    double angle = 0.0;

    static CardioidSegment lerp(const CardioidSegment& a, const CardioidSegment& b, double ratio)
    {
        CardioidSegment ret;
        ret.x = (a.x + (b.x - a.x) * ratio);
        ret.y = (a.y + (b.y - a.y) * ratio);
        ret.angle = a.angle + math::closestAngleDifference(a.angle, b.angle) * ratio;
        return ret;
    }
};

typedef std::vector<CardioidSegment> LerpedCardioid;

struct CumulativeCardioid : public std::vector<CardioidStepInfo>
{
    LerpedCardioid computeSegments(double angle_mult)
    {
        LerpedCardioid ret;

        double plot_x = 0.25;
        double plot_y = 0.0;
        double plot_direction = 0.0;

        for (int i = 0; i < size(); i++)
        {
            double step_angle = at(i).step_angle;
            double step_dist = at(i).step_dist;

            double scaled_step_angle = step_angle * angle_mult;
            plot_direction += scaled_step_angle;

            plot_x += cos(plot_direction) * step_dist;
            plot_y += sin(plot_direction) * step_dist;


            CardioidSegment segment;
            segment.x = plot_x;
            segment.y = plot_y;
            segment.angle = plot_direction;
            ret.push_back(segment);
        }

        return ret;
    }
};


struct CardioidLerper : public std::vector<LerpedCardioid>
{
    CumulativeCardioid createCumulativeCardioid(double angle_step)
    {
        CumulativeCardioid ret;
        double old_plot_x = 0.25;
        double old_plot_y = 0.0;
        double old_segment_angle = 0.0;
        int steps = static_cast<int>(round(math::tau / angle_step));
        double delta_angle;
        double delta_dist;
        double angle = 0.0;
        for (int i = 0; i <= steps; i++)
        {
            double plot_x = 0.5 * cos(angle) - 0.25 * cos(angle * 2.0);
            double plot_y = 0.5 * sin(angle) - 0.25 * sin(angle * 2.0);

            double dx = plot_x - old_plot_x;
            double dy = plot_y - old_plot_y;
            double segment_angle = atan2(dy, dx);
            delta_angle = math::closestAngleDifference(old_segment_angle, segment_angle);
            delta_dist = sqrt(dx * dx + dy * dy);

            ret.push_back({
                delta_angle,
                delta_dist,
                });

            old_segment_angle = segment_angle;
            old_plot_x = plot_x;
            old_plot_y = plot_y;
            angle += angle_step;
        }
        return ret;
    }

    void create(double angle_step, double scale_step)
    {
        CumulativeCardioid cumulative_cardioid = createCumulativeCardioid(angle_step);
        int steps = static_cast<int>(round(1.0 / scale_step));
        double scale = 0.0;
        clear();
        for (int i = 0; i <= steps; i++)
        {
            push_back(cumulative_cardioid.computeSegments(scale));
            scale += scale_step;
        }
    }

    LerpedCardioid lerped(double weight) const
    {
        double f_index = weight * static_cast<double>(size()-1);
        size_t i0 = static_cast<int>(f_index);
        size_t i1 = i0 + 1;
        if (i1 >= size())
            i1 = size() - 1;
        double lerp_ratio = f_index - static_cast<double>(i0);
        const LerpedCardioid& a = at(i0);
        const LerpedCardioid& b = at(i1);
        size_t len = a.size();
        LerpedCardioid ret;
        for (int i = 0; i < len; i++)
        {
            ret.push_back(CardioidSegment::lerp(a[i], b[i], lerp_ratio));
        }
        return ret;
    }


    int nearestPoint(double x, double y, double weight) const
    {
        double f_cardioid_index = weight * static_cast<double>(size() - 1);
        int cardioid_i0 = static_cast<int>(f_cardioid_index);
        const LerpedCardioid& cardioid = at(cardioid_i0);

        int lower = 0;
        int upper = static_cast<int>(cardioid.size()) - 1;
        int best_index = -1;
        double best_dist = std::numeric_limits<double>::max();
        double d1, d2, dx, dy;

        // Ternary search-like loop
        while (upper - lower > 2)
        {
            int mid1 = lower + (upper - lower) / 3;
            int mid2 = upper - (upper - lower) / 3;

            const auto& p1 = cardioid[mid1];
            const auto& p2 = cardioid[mid2];

            dx = (x - p1.x);
            dy = (y - p1.y);
            d1 = dx * dx + dy * dy;

            dx = (x - p2.x);
            dy = (y - p2.y);
            d2 = dx * dx + dy * dy;

            if (d1 < d2)
                upper = mid2;
            else
                lower = mid1;
        }

        // Final linear pass over narrowed range
        for (int j = lower; j <= upper; ++j)
        {
            const auto& p = cardioid[j];
            dx = x - p.x;
            dy = y - p.y;
            double d = dx * dx + dy * dy;

            if (d < best_dist)
            {
                best_dist = d;
                best_index = j;
            }
        }
        return best_index;
    }


    DVec2 originalPolarCoordinate(double x, double y, double weight) const
    {
        double f_cardioid_index = weight * static_cast<double>(size() - 1);
        int cardioid_i0 = static_cast<int>(f_cardioid_index);
        const LerpedCardioid& cardioid = at(cardioid_i0);

        int i = nearestPoint(x, y, weight);

        const CardioidSegment& p = cardioid[i];
        double angle = (i / (double)cardioid.size()) * (math::tau);
        double tangent_angle = p.angle;

        double dx = x - p.x;
        double dy = y - p.y;
        double y1 = dy * cos(tangent_angle) - dx * sin(tangent_angle);

        return {
            angle,
            -y1
        };
    }

    DVec2 project(double angle, double dist, double weight) const
    {
        double f_cardioid_index = weight * static_cast<double>(size()-1);
        size_t cardioid_i0 = static_cast<int>(f_cardioid_index);
        size_t cardioid_i1 = cardioid_i0 + 1;
        if (cardioid_i1 >= size())
            cardioid_i1 = size() - 1;

        double cardioid_lerp_ratio = f_cardioid_index - static_cast<double>(cardioid_i0);

        const LerpedCardioid& a = at(cardioid_i0);
        const LerpedCardioid& b = at(cardioid_i1);
        double f_len = static_cast<double>(a.size());

        double segment_i_ratio = angle / math::tau;
        int segment_i = static_cast<int>(segment_i_ratio * f_len);

        CardioidSegment ret = CardioidSegment::lerp(a[segment_i], b[segment_i], cardioid_lerp_ratio);
        double perp_angle = ret.angle - math::half_pi;
        
        return {
            ret.x + cos(perp_angle) * dist,
            ret.y + sin(perp_angle) * dist
        };
    }
};

struct Cardioid_Scene : public Scene<Cardioid_Scene>
{
    // --- Scene management ---
    CameraInfo          camera;
    CameraNavigator navigator;

    bool   show_offset = false;
    bool   flatten = false;
    bool   interactive = false;
    bool   animate = true;

    double ani_angle = 0.0;
    double ani_inc = math::toRadians(2.0);

    double interact_angle_step = math::tau / 720.0;
    double interact_spin_mult = 1.0;
    double interact_angle = 0.0;
    double interact_dist = 0.0;

    struct UI : ViewModel
    {
        using ViewModel::ViewModel;
        void sidebar();
    };
    
    void sceneStart() override;
    void sceneMounted(Viewport* viewport) override;
    void sceneProcess() override;

    // --- Update methods ---

    void plotCumulativeCardioid(Viewport* ctx, const CumulativeCardioid &items, double angle_mult=1.0);

    CumulativeCardioid cumulative_cardioid;
    CardioidLerper cumulative_cardioid_lookup;

    void fullPlot(Viewport* ctx, double scale, double ox) const;
    void fullPlotAlternative(Viewport* ctx, double scale, double ox) const;
    //void fullPlotAlternative2(Viewport* ctx, double scale, double ox);
    void animatePlot(Viewport* ctx, double scale, double ox, double orig_angle, double dist) const;

    // --- Viewport ---
    void viewportProcess(Viewport* ctx, double dt) override;
    void viewportDraw(Viewport* ctx) const override;

    void onEvent(Event e) override;
};

struct Cardioid_Graph_Scene : public BasicScene
{
    struct Config {
        Cardioid_Scene* main_scene;
    };

    Cardioid_Graph_Scene(Config& info) :
        main_scene(info.main_scene)
    {}

    CameraInfo          camera;
    CameraNavigator navigator;

    Cardioid_Scene* main_scene;
    CanvasImage bmp;

    void sceneStart() override;

    void sceneMounted(Viewport* viewport) override;
    void viewportProcess(Viewport* ctx, double dt) override;
    void viewportDraw(Viewport* ctx) const override;
};

struct Cardioid_Project : public BasicProject
{
    //static std::vector<std::string> categorize() {
    //    return { "Fractal", "Mandelbrot", "Main Cardioid" };
    //}

    static ProjectInfo info()
    {
        return ProjectInfo({ "Fractal", "Mandelbrot", "Main Cardioid" });
    }

    void projectPrepare(Layout& layout) override;
};

SIM_END;
