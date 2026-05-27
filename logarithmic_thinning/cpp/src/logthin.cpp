#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cfenv>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <cstdlib>

namespace {

struct CliOptions {
    std::filesystem::path input_path;
    std::filesystem::path output_path;
    double thinning_factor{1.0003};
};

void print_usage(const char *program) {
    std::cerr << "Usage: " << program
              << " --input <sorted_csv> --output <thinned_csv> --factor <float>\n";
}

CliOptions parse_arguments(int argc, char **argv) {
    CliOptions options;
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--input" && i + 1 < argc) {
            options.input_path = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            options.output_path = argv[++i];
        } else if (arg == "--factor" && i + 1 < argc) {
            options.thinning_factor = std::stod(argv[++i]);
        } else if (arg == "--help") {
            print_usage(argv[0]);
            std::exit(EXIT_SUCCESS);
        } else {
            print_usage(argv[0]);
            throw std::runtime_error("Unknown or incomplete argument: " + arg);
        }
    }

    if (options.input_path.empty()) {
        throw std::runtime_error("--input is required.");
    }
    if (options.output_path.empty()) {
        throw std::runtime_error("--output is required.");
    }
    if (options.thinning_factor <= 1.0) {
        throw std::runtime_error("--factor must be greater than 1.0.");
    }
    if (!options.output_path.has_parent_path()) {
        // remain relative to current directory
    } else {
        std::filesystem::create_directories(options.output_path.parent_path());
    }
    return options;
}

std::size_t determine_threshold(double factor) {
    const double value = factor / (factor - 1.0);
    return static_cast<std::size_t>(std::llround(value));
}

std::vector<std::size_t> logarithmic_thinning(std::size_t total_rows, double factor) {
    std::vector<std::size_t> indices;
    if (total_rows == 0) {
        return indices;
    }
    const std::size_t threshold = determine_threshold(factor);

    std::size_t index = total_rows;
    indices.push_back(index);
    while (index > threshold) {
        const double divided = std::floor(static_cast<double>(index) / factor);
        if (divided < 1.0) {
            break;
        }
        index = static_cast<std::size_t>(divided);
        if (index == 0) {
            break;
        }
        indices.push_back(index);
    }

    for (std::size_t i = threshold; i-- > 1;) {
        indices.push_back(i);
    }

    std::sort(indices.begin(), indices.end());
    indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
    return indices;
}

void apply_thinning(const CliOptions &options) {
    std::ifstream input(options.input_path);
    if (!input) {
        throw std::runtime_error("Failed to reopen input CSV: " + options.input_path.string());
    }

    std::ofstream output(options.output_path);
    if (!output) {
        throw std::runtime_error("Failed to open output CSV: " + options.output_path.string());
    }

    std::string header;
    if (!std::getline(input, header)) {
        output << "chr,pos,pp,original_index,row_number\n";
        return;
    }

    struct Row {
        int chr{};
        int pos{};
        double pp{};
        std::int64_t pp_key{};
        int original_index{};
        std::size_t original_row{};
        std::string line;
    };

    std::vector<Row> rows;
    rows.reserve(1024);

    std::fesetround(FE_TONEAREST);

    std::string line;
    std::size_t original_row_counter = 0;
    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        ++original_row_counter;

        Row row{};
        row.line = line;
        row.original_row = original_row_counter;

        std::stringstream ss(line);
        std::string field;

        if (!std::getline(ss, field, ',')) {
            throw std::runtime_error("Malformed CSV line: missing chr field.");
        }
        row.chr = std::stoi(field);

        if (!std::getline(ss, field, ',')) {
            throw std::runtime_error("Malformed CSV line: missing pos field.");
        }
        row.pos = std::stoi(field);

        if (!std::getline(ss, field, ',')) {
            throw std::runtime_error("Malformed CSV line: missing pp field.");
        }
        row.pp = std::stod(field);

        if (!std::getline(ss, field, ',')) {
            throw std::runtime_error("Malformed CSV line: missing original_index field.");
        }
        row.original_index = std::stoi(field);

        row.pp_key = static_cast<std::int64_t>(std::nearbyint(row.pp * 1'000'000.0));

        rows.push_back(std::move(row));
    }

    const std::size_t total_rows = rows.size();
    std::vector<std::size_t> keep_rows = logarithmic_thinning(total_rows, options.thinning_factor);
    const std::size_t top_n = std::min<std::size_t>(total_rows, 2000);
    for (std::size_t rank = 1; rank <= top_n; ++rank) {
        keep_rows.push_back(rank);
    }
    std::sort(keep_rows.begin(), keep_rows.end());
    keep_rows.erase(std::unique(keep_rows.begin(), keep_rows.end()), keep_rows.end());

    std::vector<std::size_t> order(total_rows);
    std::iota(order.begin(), order.end(), 0);

    std::sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
        const Row &a = rows[lhs];
        const Row &b = rows[rhs];
        if (a.pp_key != b.pp_key) return a.pp_key > b.pp_key;
        if (a.chr != b.chr) return a.chr < b.chr;
        if (a.pos != b.pos) return a.pos < b.pos;
        if (a.original_index != b.original_index) return a.original_index < b.original_index;
        return a.original_row < b.original_row;
    });

    std::vector<std::pair<std::size_t, const Row *>> kept;
    kept.reserve(keep_rows.size());
    for (std::size_t rank : keep_rows) {
        if (rank == 0 || rank > order.size()) {
            continue;
        }
        const std::size_t row_index = order[rank - 1];
        kept.emplace_back(rank, &rows[row_index]);
    }

    std::sort(kept.begin(), kept.end(), [](const auto &lhs, const auto &rhs) {
        return lhs.first < rhs.first;
    });

    output << header << ",row_number\n";
    for (const auto &entry : kept) {
        output << entry.second->line << ',' << entry.first << '\n';
    }
}

}  // namespace

int main(int argc, char **argv) {
    try {
        const CliOptions options = parse_arguments(argc, argv);
        apply_thinning(options);
        return EXIT_SUCCESS;
    } catch (const std::exception &ex) {
        std::cerr << "logthin error: " << ex.what() << '\n';
        return EXIT_FAILURE;
    }
}
