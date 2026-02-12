#pragma once
#include <bitloop.h>
#include <fstream>
#include <tsl/ordered_map.h>

#include "../config/build_config.h"

SIM_BEG;

using namespace bl;

struct MandelBookmark : public Hashable
{
    std::string data;
    std::string loaded_path;

    GLuint thumb_tex = 0; // must be loaded on UI-side

    MandelBookmark() {}
    MandelBookmark(const MandelBookmark& rhs)
    { 
        data = rhs.data;
        thumb_tex = rhs.thumb_tex;
        loaded_path = rhs.loaded_path;
    }
    MandelBookmark(std::string_view _data)
    {
        data = _data;
        invalidate_hash();
    }

    void destroyTexture()
    {
        if (thumb_tex != 0)
        {
            glDeleteTextures(1, (GLuint*)&thumb_tex);
            thumb_tex = 0;
        }
    }
    
    hash_t compute_hash() const noexcept override
    {
        return Hasher::hash_string(data);
    }

    bool loadFromThumbnail(std::string_view webp_path)
    {
        std::ifstream in(webp_path.data(), std::ios::binary);
        if (!in)
            return false;

        std::string txt;
        in.seekg(0, std::ios::end);
        txt.resize((std::size_t)in.tellg());
        in.seekg(0, std::ios::beg);
        in.read(txt.data(), (std::streamsize)txt.size());

        bytebuf webp_data = bytebuf(txt.data(), txt.data() + txt.size());
        webp_extract_save_from_xmp(webp_data, data);

        invalidate_hash();

        loaded_path = webp_path;

        return true;
    }

    // saves existing thumbnail with updated XMP data
    void save(std::string_view webp_path)
    {
        std::ifstream in(loaded_path, std::ios::binary);
        if (!in) throw std::runtime_error("failed to open webp for read");

        in.seekg(0, std::ios::end);
        std::string bytes;
        bytes.resize(static_cast<std::size_t>(in.tellg()));
        in.seekg(0, std::ios::beg);
        in.read(bytes.data(), static_cast<std::streamsize>(bytes.size()));

        bytebuf webp(bytes.begin(), bytes.end());

        if (!webp_set_save_string_as_xmp_inplace(webp, data))
            throw std::runtime_error("failed to update XMP chunk");

        std::error_code ec;
        std::filesystem::create_directories(std::filesystem::path(webp_path).parent_path(), ec);

        std::ofstream out(std::string(webp_path), std::ios::binary | std::ios::trunc);
        if (!out) throw std::runtime_error("failed to open webp for write");
        out.write(reinterpret_cast<const char*>(webp.data()), static_cast<std::streamsize>(webp.size()));
    }

    std::string thumbName() const  // stem only
    { return compression::b64_encode_u64(stable_hash()); }          

    std::string thumbFilename() const  // + extension
    { return thumbName() + ".webp"; }                            

    std::string thumbPath() const  // + full parent path
    { return ProjectBase::activeProject()->root_path("data/bookmarks/" + thumbFilename()); } 

    void loadThumbnail()
    {
        if (thumb_tex == 0)
            thumb_tex = loadGLTextureRGBA8(thumbPath().c_str());
    }
    void loadThumbnail(std::string_view path)
    {
        // used when loading from directory
        if (thumb_tex == 0)
            thumb_tex = loadGLTextureRGBA8(path.data());
    }
    void loadThumbnail(const bytebuf& data, int w, int h)
    {
        // used when inserting bookmark from capture
        thumb_tex = createGLTextureRGBA8(data.data(), w, h);
    }

    int thumbTexture()
    {
        return thumb_tex;
    }
};

class MandelBookmarkList : public Hashable
{
    std::vector<MandelBookmark> items;

public:

    MandelBookmarkList() {}
    MandelBookmarkList(const MandelBookmarkList& rhs) 
    {
        items = rhs.items; 
    }
    MandelBookmarkList& operator =(const MandelBookmarkList& rhs) { 
        items = rhs.items;
        return *this; 
    }

    hash_t compute_hash() const noexcept override
    {
        StableHasher h;
        for (size_t i = 0; i < items.size(); i++) 
            h.add(items[i].hash());
        return h.finish();
    }

    std::vector<MandelBookmark>&       getItems() { return items; }
    const std::vector<MandelBookmark>& getItems() const { return items; }

    void addItem(std::string_view data)
    {
        items.push_back(MandelBookmark(data));
        invalidate_hash();
    }

    void addItem(const MandelBookmark& bookmark)
    {
        items.push_back(bookmark);
        invalidate_hash();
    }

    bool loadDirectoryBookmarks(const char* dir)
    {
        namespace fs = std::filesystem;

        std::error_code ec;

        if (!fs::exists(dir, ec) || !fs::is_directory(dir, ec))
            return false;

        items.clear();
        for (const fs::directory_entry& e : fs::directory_iterator(dir, fs::directory_options::skip_permission_denied, ec))
        {
            if (ec) break;
            if (!e.is_regular_file(ec) || ec) { ec.clear(); continue; }

            const fs::path& p = e.path();
            //blPrint() << "loadDirectoryBookmarks p: " << p.string();

            // Only *.webp
            auto ext = p.extension().string();
            for (char& c : ext) c = (char)std::tolower((unsigned char)c);
            if (ext != ".webp")
                continue;

            std::string webp_path = p.lexically_normal().string().c_str();

            MandelBookmark bookmark;
            bookmark.loadThumbnail(webp_path);
            bookmark.loadFromThumbnail(webp_path);
            items.push_back(bookmark);
        }

        invalidate_hash();
        return true;
    }
};

class MandelBookmarkManager : public Hashable
{
    tsl::ordered_map<std::string, MandelBookmarkList> lists;

public:

    #if MANDEL_RUN_DICTIONARY_TUNINGS
    std::vector<std::string> shaders;
    #endif

    MandelBookmarkManager() {}
    MandelBookmarkManager(const MandelBookmarkManager& rhs)
    { 
        lists = rhs.lists; 
    }
    MandelBookmarkManager& operator =(const MandelBookmarkManager& rhs) { 
        lists = rhs.lists; 
        return *this;
    }
    
    hash_t compute_hash() const noexcept override
    {
        StableHasher h;
        for (auto it = lists.begin(); it != lists.end(); it++)
        {
            h.add_string(it->first);
            h.add(it->second.hash());
        }
        return h.finish();
    }

    MandelBookmarkList& ensureCategory(std::string name)
    {
        if (!lists.count(name)) lists[name] = MandelBookmarkList();
        invalidate_hash();
        return lists[name];
    }

    MandelBookmarkList& find(std::string name) 
    { 
        MandelBookmarkList& ret = lists[name];
        invalidate_hash();
        return ret;
    }

    std::pair<const std::string&, MandelBookmarkList&> at(int i)
    {
        auto it = lists.begin() + static_cast<std::ptrdiff_t>(i);
        return { it->first, it.value() };
    }

    int size() { return (int)lists.size(); }

    void loadCategoryDirs(std::string_view parent_dir)
    {
        namespace fs = std::filesystem;

        std::error_code ec;
        std::string root_dir = ProjectBase::activeProject()->root_path(parent_dir);

        blPrint() << "loadCategoryDirs root dir: " << root_dir;

        if (!fs::exists(root_dir, ec) || !fs::is_directory(root_dir, ec))
            return;

        for (const fs::directory_entry& e : fs::directory_iterator(root_dir, fs::directory_options::skip_permission_denied, ec))
        {
            if (ec) break;
            if (!e.is_directory(ec) || ec) { ec.clear(); continue; }

            // Load bookmarks from directory, using directory name as list name
            const fs::path& p = e.path();

            blPrint() << "loadCategoryDirs p: " << p.string();

            const std::string list_name = p.stem().string();
            lists[list_name].loadDirectoryBookmarks(p.lexically_normal().string().c_str());
        }

        invalidate_hash();

        // debug compression (breaks emscripten builds, run on desktop only)
        #if MANDEL_RUN_DICTIONARY_TUNINGS
        for (int i = 0; i < lists.size(); i++)
        {
            auto [list_name, list] = at(i);
            std::vector<MandelBookmark>& bookmarks = list.getItems();
            for (auto& bookmark : bookmarks)
            {
                MandelState test_compression;
                test_compression.deserialize(bookmark.data);
                std::string& shader_source_txt = test_compression.shader_source_txt;
                bool found = false;
                for (std::string& shader : shaders)
                {
                    if (shader == shader_source_txt)
                    {
                        found = true;
                        break;
                    }
                }
                if (!found)
                    shaders.push_back(shader_source_txt);
            }
        }
        #endif
    }

    void saveAll(std::string_view dir)
    {
        for (int i = 0; i < lists.size(); i++)
        {
            auto [list_name, list] = at(i);
            std::vector<MandelBookmark>& bookmarks = list.getItems();
            for (auto& bookmark : bookmarks)
            {
                std::filesystem::path filepath = ProjectBase::activeProject()->root_path("data/bookmarks/");
                filepath /= dir;
                filepath /= list_name;
                filepath /= bookmark.thumbFilename();
                bookmark.save(filepath.lexically_normal().string());
            }
        }
    }

    void loadMissingThumbnails()
    {
        for (int i=0; i<lists.size(); i++)
        {
            auto [list_name, list] = at(i);
            std::vector<MandelBookmark>& bookmarks = list.getItems();
            for (auto& bookmark : bookmarks)
            {
                bookmark.loadThumbnail();
            }
        }
    }

    void destroyAllTextures()
    {
        for (int i = 0; i < lists.size(); i++)
        {
            auto [list_name, list] = at(i);
            std::vector<MandelBookmark>& bookmarks = list.getItems();
            for (auto& bookmark : bookmarks)
            {
                bookmark.destroyTexture();
            }
        }
    }
};

enum class DictScoreMode
{
    Sum,      // sum of compressed lengths across corpus
    Max,      // worst-case compressed length across corpus
};

struct DictTuneConfig
{
    int quality = 11;
    int window = 22;

    DictScoreMode score_mode = DictScoreMode::Sum;

    bool start_with_all_tokens = true;

    // enables expensive but sometimes helpful moves that remove two tokens at once
    bool enable_pair_removal = false;
    std::size_t pair_removal_budget_per_pass = 4096;

    std::size_t max_passes = 128;

    // if true, prints progress to stdout.
    bool verbose = true;
};

struct DictTuneResult
{
    std::vector<std::string_view> kept;
    std::vector<std::string_view> removed;

    std::size_t baseline_score = 0;
    std::size_t final_score = 0;

    std::vector<std::string> log;

    std::string to_initializer(std::string_view var_name) const
    {
        constexpr std::size_t kMaxCols = 150;

        auto escape_sv = [](std::string_view s) -> std::string
        {
            std::string out;
            out.reserve(s.size());
            for (char c : s)
            {
                if (c == '\\') out += "\\\\";
                else if (c == '\"') out += "\\\"";
                else if (c == '\n') out += "\\n";
                else if (c == '\t') out += "\\t";
                else out += c;
            }
            return out;
        };

        std::string out;
        out += "static constexpr auto ";
        out += var_name;
        out += " = std::to_array<std::string_view>({\n";

        const std::string indent = "    ";
        std::string line = indent;

        for (std::size_t i = 0; i < kept.size(); ++i)
        {
            std::string token;
            token.reserve(kept[i].size() + 8);

            token += '"';
            token += escape_sv(kept[i]);
            token += '"';

            // omit trailing comma on last element for cleaner diffs
            if (i + 1 != kept.size())
                token += ",";

            const std::string piece = (line.size() == indent.size()) ? token : (" " + token);

            if (line.size() + piece.size() > kMaxCols && line.size() > indent.size())
            {
                out += line;
                out += "\n";
                line = indent;
                line += token;
            }
            else
            {
                line += piece;
            }
        }

        if (line.size() > indent.size())
        {
            out += line;
            out += "\n";
        }

        out += "});\n";
        return out;
    }

};

static inline std::uint64_t fnv1a64(std::span<const std::uint64_t> data)
{
    std::uint64_t h = 14695981039346656037ull;
    for (std::uint64_t x : data)
    {
        h ^= x;
        h *= 1099511628211ull;
    }
    return h;
}

static inline std::uint64_t bitmask_hash(const std::vector<std::uint8_t>& enabled)
{
    // hashes enabled mask in 64-bit chunks for caching
    std::vector<std::uint64_t> packed;
    packed.reserve((enabled.size() + 7) / 8);
    std::uint64_t acc = 0;
    int bit = 0;

    for (std::size_t i = 0; i < enabled.size(); ++i)
    {
        acc |= (std::uint64_t(enabled[i] ? 1 : 0) << bit);
        ++bit;
        if (bit == 64)
        {
            packed.push_back(acc);
            acc = 0;
            bit = 0;
        }
    }
    if (bit != 0) packed.push_back(acc);

    return fnv1a64(std::span<const std::uint64_t>(packed.data(), packed.size()));
}

static inline std::vector<std::string_view> collect_tokens(
    std::span<const std::string_view> pool,
    const std::vector<std::uint8_t>& enabled)
{
    std::vector<std::string_view> out;
    out.reserve(pool.size());
    for (std::size_t i = 0; i < pool.size(); ++i)
        if (enabled[i]) out.push_back(pool[i]);
    return out;
}

static inline std::size_t score_dictionary(
    const std::vector<std::string>& corpus,
    std::span<const std::string_view> dict_tokens,
    const DictTuneConfig& cfg)
{
    compression::BrotliDict dict(dict_tokens);

    std::size_t sum = 0;
    std::size_t mx = 0;

    for (const std::string& s : corpus)
    {
        const std::size_t n = compression::brotli_ascii_compress_with_dict(s, cfg.quality, cfg.window, &dict).size();
        sum += n;
        mx = std::max(mx, n);
    }

    return (cfg.score_mode == DictScoreMode::Sum) ? sum : mx;
}

static DictTuneResult tune_brotli_dictionary(
    const std::vector<std::string>& corpus,
    std::span<const std::string_view> token_pool,
    const DictTuneConfig& cfg)
{
    if (corpus.empty())
        throw std::runtime_error("tune_brotli_dictionary: corpus is empty");
    if (token_pool.empty())
        throw std::runtime_error("tune_brotli_dictionary: token_pool is empty");

    std::vector<std::uint8_t> enabled(token_pool.size(), cfg.start_with_all_tokens ? 1 : 0);

    std::unordered_map<std::uint64_t, std::size_t> cache;
    cache.reserve(8192);

    auto eval = [&](const std::vector<std::uint8_t>& mask) -> std::size_t
    {
        const std::uint64_t key = bitmask_hash(mask);
        if (auto it = cache.find(key); it != cache.end())
            return it->second;

        std::vector<std::string_view> dict_tokens = collect_tokens(token_pool, mask);
        const std::size_t score = score_dictionary(corpus, dict_tokens, cfg);
        cache.emplace(key, score);
        return score;
    };

    DictTuneResult result;

    const std::size_t baseline = eval(enabled);
    result.baseline_score = baseline;
    std::size_t best_score = baseline;

    if (cfg.verbose)
        blPrint() << "baseline score: " << baseline << "\n";

    auto log_step = [&](std::string msg)
    {
        if (cfg.verbose) blPrint() << msg << "\n";
        result.log.push_back(std::move(msg));
    };

    for (std::size_t pass = 0; pass < cfg.max_passes; ++pass)
    {
        bool improved_any = false;

        // ---- best single removal ----
        {
            std::size_t best_i = std::numeric_limits<std::size_t>::max();
            std::size_t best_candidate = best_score;

            for (std::size_t i = 0; i < enabled.size(); ++i)
            {
                if (!enabled[i]) continue;
                enabled[i] = 0;
                const std::size_t s = eval(enabled);
                enabled[i] = 1;

                if (s < best_candidate)
                {
                    best_candidate = s;
                    best_i = i;
                }
            }

            if (best_i != std::numeric_limits<std::size_t>::max())
            {
                enabled[best_i] = 0;
                best_score = best_candidate;
                improved_any = true;
                log_step("pass " + std::to_string(pass) + ": remove [" + std::to_string(best_i) +
                    "] -> score " + std::to_string(best_score));
            }
        }

        // ---- pair removal ----
        if (cfg.enable_pair_removal)
        {
            std::size_t budget = cfg.pair_removal_budget_per_pass;

            std::size_t best_i = std::numeric_limits<std::size_t>::max();
            std::size_t best_j = std::numeric_limits<std::size_t>::max();
            std::size_t best_candidate = best_score;

            for (std::size_t i = 0; i < enabled.size(); ++i)
            {
                if (!enabled[i]) continue;
                for (std::size_t j = i + 1; j < enabled.size(); ++j)
                {
                    if (!enabled[j]) continue;

                    enabled[i] = 0;
                    enabled[j] = 0;
                    const std::size_t s = eval(enabled);
                    enabled[i] = 1;
                    enabled[j] = 1;

                    if (s < best_candidate)
                    {
                        best_candidate = s;
                        best_i = i;
                        best_j = j;
                    }

                    if (--budget == 0) break;
                }
                if (budget == 0) break;
            }

            if (best_i != std::numeric_limits<std::size_t>::max())
            {
                enabled[best_i] = 0;
                enabled[best_j] = 0;
                best_score = best_candidate;
                improved_any = true;
                log_step("pass " + std::to_string(pass) + ": remove pair [" + std::to_string(best_i) +
                    "," + std::to_string(best_j) + "] -> score " + std::to_string(best_score));
            }
        }

        // ---- best single addition ----
        {
            std::size_t best_i = std::numeric_limits<std::size_t>::max();
            std::size_t best_candidate = best_score;

            for (std::size_t i = 0; i < enabled.size(); ++i)
            {
                if (enabled[i]) continue;
                enabled[i] = 1;
                const std::size_t s = eval(enabled);
                enabled[i] = 0;

                if (s < best_candidate)
                {
                    best_candidate = s;
                    best_i = i;
                }
            }

            if (best_i != std::numeric_limits<std::size_t>::max())
            {
                enabled[best_i] = 1;
                best_score = best_candidate;
                improved_any = true;
                log_step("pass " + std::to_string(pass) + ": add [" + std::to_string(best_i) +
                    "] -> score " + std::to_string(best_score));
            }
        }

        if (!improved_any)
            break;
    }

    result.final_score = best_score;

    // finalize kept/removed lists
    result.kept.reserve(token_pool.size());
    result.removed.reserve(token_pool.size());
    for (std::size_t i = 0; i < token_pool.size(); ++i)
        (enabled[i] ? result.kept : result.removed).push_back(token_pool[i]);

    if (cfg.verbose)
        blPrint() << "final score: " << result.final_score
        << " (delta " << (int64_t(result.final_score) - int64_t(result.baseline_score)) << ")\n";

    return result;
}

SIM_END;
