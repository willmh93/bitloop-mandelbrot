#include "Mandelbrot.h"

SIM_BEG;

void Mandelbrot_Scene::UI::populateHistory()
{
    if (ImGui::CollapsingHeaderBox("History"))
    {
        bl_pull(history, current_history_i);
        if (ImGui::BeginListBox("History"))
        {
            for (size_t i = 0; i < history.size(); i++)
            {
                std::string txt = "Item " + std::to_string(i);
                if (ImGui::Selectable(txt.c_str(), i==current_history_i))
                {
                    bl_schedule([i](Mandelbrot_Scene& scene)
                    {
                        scene.current_history_i = (int)i;
                        scene.loadState(scene.history[scene.current_history_i].second);
                        scene.loading_history_state = true;
                    });
                }
            }
            ImGui::EndListBox();
        }
        ImGui::EndCollapsingHeaderBox();
    }
}

void Mandelbrot_Scene::processUndoRedo(bool normalization_opts_changed [[maybe_unused]], bool gradient_changed [[maybe_unused]] )
{
    // make sure we don't "record" undo-events
    if (!loading_history_state)
    {
        bool recomputed = mandelChanged();

        // check if any normalization opts changed due to user input (ignore changes where animation is in effect)
        bool weights_changed = Changed(iter_weight, dist_weight, stripe_weight);
        bool shade_formula_changed = Changed(
            iter_x_dist_weight, dist_x_stripe_weight, stripe_x_iter_weight,
            iter_x_distStripe_weight, dist_x_iterStripe_weight, stripe_x_iterDist_weight);

        bool cycle_iter_opts_changed = Changed(iter_params);
        bool cycle_dist_opts_changed = Changed(dist_params, dist_tone_params);
        bool cycle_stripe_opts_changed = Changed(stripe_tone_params) || (!animate_stripe_phase && Changed(stripe_params.phase));

        bool normalization_opts_changed = 
            (shade_formula_changed || weights_changed || cycle_iter_opts_changed || cycle_dist_opts_changed || cycle_stripe_opts_changed);

        // check if base gradient changed, or gradient shift changed (only if no animation present)
        bool gradient_opts_changed = Changed(gradient) || (!animate_gradient_shift && Changed(gradient_shift));
        bool animate_opts_changed = Changed(animate_gradient_shift, animate_gradient_hue, animate_stripe_phase);

        // track how long since we last changed
        if (recomputed || normalization_opts_changed || gradient_opts_changed || animate_opts_changed)
            pending_checkpoint_flag = true;

        if (pending_checkpoint_flag &&       // check a change significant *could* have occured
            camera_vel_pos.isZero() &&       // check no pan inertia
            (camera_vel_zoom == 1) &&        // check no zoom inertia
            !mouse->pressed &&               // check not currently panning
            !main_window()->isEditingUI())   // check not presently dragging imgui value slider
        {
            // serialize mandel state
            std::string serialized_state = serialize();

            // check if the "new" state is different to the "current" state
            if (history.empty() || history[current_history_i].second != serialized_state)
            {
                // discard "future" changes
                if (current_history_i < history.size() - 1)
                    history.erase(history.begin() + current_history_i + 1, history.end());

                // start new branch, add current state as "now"
                history.push_back(std::make_pair(current_history_i, serialized_state));
                current_history_i++;

                // limit history size
                if (history.size() > history_count)
                {
                    history.erase(history.begin());
                    current_history_i--;
                }
            }

            pending_checkpoint_flag = false;
        }
    }

    // start tracking history changes again (now that the previous load's knock-on changes have applied)
    loading_history_state = false;
}

void Mandelbrot_Scene::historyKeyEvent(KeyEvent e)
{
    if (!platform()->is_mobile())
    {
        if (e.keyMod() & SDL_KMOD_LCTRL)
        {
            if (e.scanCode() == SDL_SCANCODE_Y ||
                (e.scanCode() == SDL_SCANCODE_Z && e.keyMod() & SDL_KMOD_LSHIFT))
            {
                if (current_history_i < history.size() - 1)
                {
                    current_history_i++;
                    loadState(history[current_history_i].second);
                    loading_history_state = true;
                }
            }
            else if (e.scanCode() == SDL_SCANCODE_Z)
            {
                if (current_history_i > 0)
                {
                    current_history_i--;
                    loading_history_state = true;
                    loadState(history[current_history_i].second);
                }
            }
        }
    }
}

SIM_END;