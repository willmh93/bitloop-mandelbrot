#include "Mandelbrot.h"

SIM_BEG;

void Mandelbrot_Scene::collectStats()
{
    if (!this->active_field) return;
    EscapeField& field = *this->active_field;

    memset(&stats.dirty, 0, sizeof(stats.dirty));

    // Depth Histogram
    {
        int bucket_size = std::max(1, (int)((field.max_depth - field.min_depth) / 1000.0));
        if (bucket_size > 0)
        {
            auto& hist = stats.depth_histogram;
            hist.clear();
            for (size_t i = 0; i < field.size(); i++)
            {
                EscapeFieldPixel& p = field[i];
                if (p.depth > 1e10 || isnan(p.depth)) continue;

                int bucket_depth = (int)(p.depth / bucket_size) * bucket_size;
                hist[bucket_depth]++;
            }

            // Trim small isolated entries from back of entries
            if (hist.size() >= 2)
            {
                int min_depth = hist.begin()->first;
                int max_depth = (--hist.end())->first;

                i64 sum_count = 0;
                for (auto [k, c] : hist) sum_count += c;
                int mean_px_count = (int)((double)sum_count / (double)hist.size());
                int min_px_count = mean_px_count / 50;
                auto valid_entry = [&](int d) { return hist.count(d) && hist[d] >= min_px_count; };

                // Find 'x' consecutive valid entries from the back and stop
                const int consecutive_entries = 5;
                int trim_from_depth;
                for (trim_from_depth = max_depth; trim_from_depth >= min_depth; trim_from_depth -= bucket_size)
                {
                    int d = (int)(trim_from_depth / bucket_size) * bucket_size;
                    bool gaps = false;
                    for (int j = 0; j < consecutive_entries; j++)
                    {
                        if (!valid_entry(d - bucket_size * j)) {
                            gaps = true;
                            break;
                        }
                    }
                    if (!gaps) break;
                }

                auto it = hist.upper_bound(trim_from_depth);
                hist.erase(it, hist.end());
                max_depth = trim_from_depth;

                for (int d = min_depth + bucket_size; d < max_depth; d += bucket_size)
                {
                    int bucket_depth = (int)(d / bucket_size) * bucket_size;
                    if (hist.count(bucket_depth) == 0)
                        hist[bucket_depth] = 0;
                }
            }
            stats.dirty.depth_histogram = true;
        }
    }
}

SIM_END;