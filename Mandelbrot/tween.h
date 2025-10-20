#include <bitloop.h>

SIM_BEG;

struct Mandelbrot_Scene_Data;
struct MandelState;


void startTween(Mandelbrot_Scene& scene_data);

void lerpState(
    Mandelbrot_Scene& scene_data,
    MandelState& a,
    MandelState& b,
    double f, bool complete);

SIM_END;
