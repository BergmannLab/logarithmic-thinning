#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <memory>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <cstdlib>

namespace {

struct Record {
    double pp;
    std::int32_t original_index;
    std::int32_t chr;
    std::int32_t pos;
};

struct ChunkMetadata {
    std::filesystem::path path;
    std::uint64_t record_count;
};

struct HeapItem {
    Record record{};
    std::size_t chunk_index{};
};

class ChunkReader {
public:
    explicit ChunkReader(const std::filesystem::path &path) : stream_(path, std::ios::binary) {
        if (!stream_) {
            throw std::runtime_error("Failed to open chunk file: " + path.string());
        }
        stream_.read(reinterpret_cast<char *>(&remaining_), sizeof(remaining_));
        if (!stream_) {
            throw std::runtime_error("Failed to read chunk metadata from: " + path.string());
        }
        fetch_next();
    }

    bool has_value() const { return has_value_; }

    const Record &current() const { return current_; }

    void advance() { fetch_next(); }

private:
    void fetch_next() {
        if (remaining_ == 0) {
            has_value_ = false;
            return;
        }
        stream_.read(reinterpret_cast<char *>(&current_), sizeof(Record));
        if (!stream_) {
            throw std::runtime_error("Failed to read record from chunk.");
        }
        --remaining_;
        has_value_ = true;
    }

    std::ifstream stream_;
    std::uint64_t remaining_{0};
    Record current_{};
    bool has_value_{false};
};

std::vector<std::string> split_whitespace(const std::string &line) {
    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::string token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

double parse_double(const std::string &token) {
    if (token.empty()) {
        return -std::numeric_limits<double>::infinity();
    }
    if (token == "NA" || token == "NaN") {
        return -std::numeric_limits<double>::infinity();
    }
    char *end = nullptr;
    const double value = std::strtod(token.c_str(), &end);
    if (end == token.c_str()) {
        throw std::runtime_error("Failed to parse floating point value from token: " + token);
    }
    return value;
}

std::int32_t parse_int32(const std::string &token, const char *context) {
    if (token.empty()) {
        throw std::runtime_error(std::string("Unexpected empty token while parsing ") + context);
    }
    char *end = nullptr;
    const long value = std::strtol(token.c_str(), &end, 10);
    if (end == token.c_str()) {
        throw std::runtime_error("Failed to parse integer value from token: " + token);
    }
    if (value < std::numeric_limits<std::int32_t>::min() || value > std::numeric_limits<std::int32_t>::max()) {
        throw std::runtime_error("Value out of 32-bit range while parsing " + std::string(context));
    }
    return static_cast<std::int32_t>(value);
}

struct CliOptions {
    std::filesystem::path input_path;
    std::filesystem::path output_path;
    std::filesystem::path temporary_directory;
    double threshold{0.0};
    std::size_t rows_per_chunk{1'000'000};
    bool keep_temporary{false};
    bool disable_binary_output{false};
};

void print_usage(const char *program) {
    std::cerr
        << "Usage: " << program
        << " --input <bgenie_file> --output <sorted_csv>"
        << " [--threshold <float>] [--chunk-rows <rows>] [--tmpdir <dir>]"
        << " [--keep-temp] [--disable-binary]\n";
}

CliOptions parse_arguments(int argc, char **argv) {
    CliOptions options;
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--input" && i + 1 < argc) {
            options.input_path = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            options.output_path = argv[++i];
        } else if (arg == "--threshold" && i + 1 < argc) {
            options.threshold = std::stod(argv[++i]);
        } else if (arg == "--chunk-rows" && i + 1 < argc) {
            long long rows = std::stoll(argv[++i]);
            if (rows <= 0) {
                throw std::runtime_error("--chunk-rows must be positive.");
            }
            options.rows_per_chunk = static_cast<std::size_t>(rows);
        } else if (arg == "--tmpdir" && i + 1 < argc) {
            options.temporary_directory = argv[++i];
        } else if (arg == "--keep-temp") {
            options.keep_temporary = true;
        } else if (arg == "--disable-binary") {
            options.disable_binary_output = true;
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
    if (options.temporary_directory.empty()) {
        options.temporary_directory = options.output_path.parent_path();
    }
    if (options.temporary_directory.empty()) {
        options.temporary_directory = std::filesystem::current_path();
    }

    return options;
}

std::vector<std::size_t> detect_log_columns(const std::vector<std::string> &header_tokens) {
    std::vector<std::size_t> indices;
    for (std::size_t i = 0; i < header_tokens.size(); ++i) {
        const std::string &column = header_tokens[i];
        const std::string suffix = "true-log10p";
        if (column.size() >= suffix.size() &&
            column.compare(column.size() - suffix.size(), suffix.size(), suffix) == 0) {
            indices.push_back(i);
        }
    }
    return indices;
}

std::string timestamped_chunk_name(std::size_t index) {
    const auto now = std::chrono::steady_clock::now().time_since_epoch();
    const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(now).count();
    return "chunk_" + std::to_string(micros) + "_" + std::to_string(index) + ".bin";
}

void write_chunk(const std::filesystem::path &path, std::vector<Record> &records,
                 std::vector<ChunkMetadata> &chunks) {
    if (records.empty()) {
        return;
    }
    std::sort(records.begin(), records.end(), [](const Record &lhs, const Record &rhs) {
        return lhs.pp > rhs.pp;
    });

    std::ofstream stream(path, std::ios::binary);
    if (!stream) {
        throw std::runtime_error("Failed to open chunk file for writing: " + path.string());
    }
    const std::uint64_t count = static_cast<std::uint64_t>(records.size());
    stream.write(reinterpret_cast<const char *>(&count), sizeof(count));
    stream.write(reinterpret_cast<const char *>(records.data()), sizeof(Record) * records.size());
    if (!stream) {
        throw std::runtime_error("Failed to write chunk file: " + path.string());
    }
    chunks.push_back({path, count});
    records.clear();
}

void merge_chunks(const std::vector<ChunkMetadata> &chunks, const std::filesystem::path &csv_output,
                  bool disable_binary_output) {
    if (chunks.empty()) {
        std::ofstream csv(csv_output);
        if (!csv) {
            throw std::runtime_error("Failed to open output CSV: " + csv_output.string());
        }
        csv << "chr,pos,pp,original_index\n";
        return;
    }

    const auto binary_output = csv_output.parent_path() / "final_sorted_data.bin";

    std::priority_queue<HeapItem, std::vector<HeapItem>,
                        bool (*)(const HeapItem &, const HeapItem &)>
        heap([](const HeapItem &lhs, const HeapItem &rhs) { return lhs.record.pp < rhs.record.pp; });

    std::vector<std::unique_ptr<ChunkReader>> readers;
    readers.reserve(chunks.size());
    for (const auto &chunk : chunks) {
        auto reader = std::make_unique<ChunkReader>(chunk.path);
        if (reader->has_value()) {
            heap.push({reader->current(), readers.size()});
        }
        readers.push_back(std::move(reader));
    }

    std::ofstream csv(csv_output);
    if (!csv) {
        throw std::runtime_error("Failed to open output CSV: " + csv_output.string());
    }
    csv << "chr,pos,pp,original_index\n";
    csv.setf(std::ios::fixed);
    csv.precision(6);

    std::ofstream bin;
    if (!disable_binary_output) {
        bin.open(binary_output, std::ios::binary);
        if (!bin) {
            throw std::runtime_error("Failed to open binary output: " + binary_output.string());
        }
        const char magic[8] = {'L', 'O', 'G', 'S', 'O', 'R', 'T', '\0'};
        bin.write(magic, sizeof(magic));
    }

    while (!heap.empty()) {
        HeapItem item = heap.top();
        heap.pop();

        csv << item.record.chr << ',' << item.record.pos << ','
            << item.record.pp << ',' << item.record.original_index << '\n';
        if (bin) {
            bin.write(reinterpret_cast<const char *>(&item.record), sizeof(Record));
        }

        auto &reader = readers[item.chunk_index];
        reader->advance();
        if (reader->has_value()) {
            heap.push({reader->current(), item.chunk_index});
        }
    }
}

}  // namespace

int main(int argc, char **argv) {
    try {
        const CliOptions options = parse_arguments(argc, argv);

        if (!std::filesystem::exists(options.input_path)) {
            throw std::runtime_error("Input file does not exist: " + options.input_path.string());
        }

        std::ifstream input(options.input_path);
        if (!input) {
            throw std::runtime_error("Failed to open input file: " + options.input_path.string());
        }

        std::filesystem::create_directories(options.temporary_directory);
        std::filesystem::create_directories(options.output_path.parent_path());

        std::string header_line;
        if (!std::getline(input, header_line)) {
            throw std::runtime_error("Input file is empty.");
        }
        std::vector<std::string> header_tokens = split_whitespace(header_line);
        if (header_tokens.size() < 3) {
            throw std::runtime_error("Unexpected header format; expected at least chr, pos and one trait.");
        }

        const std::vector<std::size_t> feature_indices = detect_log_columns(header_tokens);
        if (feature_indices.empty()) {
            throw std::runtime_error("No columns ending with 'true-log10p' found in header.");
        }

        if (feature_indices.front() <= 1) {
            // indices[0] == 0 would imply chr column flagged, so we sanity check
            for (std::size_t idx : feature_indices) {
                if (idx <= 1) {
                    throw std::runtime_error("Detected trait column index overlapping structural columns.");
                }
            }
        }

        std::vector<Record> buffer;
        buffer.reserve(options.rows_per_chunk * feature_indices.size());
        std::vector<ChunkMetadata> chunks;
        chunks.reserve(64);

        std::size_t current_chunk_rows = 0;
        std::size_t chunk_counter = 0;
        std::int64_t global_row_index = 0;
        std::string line;

        while (std::getline(input, line)) {
            if (line.empty()) {
                continue;
            }
            std::vector<std::string> tokens = split_whitespace(line);
            if (tokens.size() < header_tokens.size()) {
                throw std::runtime_error("Encountered row with fewer columns than header.");
            }

            const std::int32_t chr = parse_int32(tokens[0], "chr");
            const std::int32_t pos = parse_int32(tokens[1], "pos");

            for (std::size_t column_index : feature_indices) {
                const double value = parse_double(tokens[column_index]);
                if (value > options.threshold) {
                    Record record{};
                    record.pp = value;
                    record.original_index = static_cast<std::int32_t>(global_row_index);
                    record.chr = chr;
                    record.pos = pos;
                    buffer.push_back(record);
                }
            }

            ++global_row_index;
            ++current_chunk_rows;

            if (current_chunk_rows >= options.rows_per_chunk) {
                const std::filesystem::path chunk_path =
                    options.temporary_directory / timestamped_chunk_name(chunk_counter++);
                write_chunk(chunk_path, buffer, chunks);
                current_chunk_rows = 0;
            }
        }

        if (!buffer.empty()) {
            const std::filesystem::path chunk_path =
                options.temporary_directory / timestamped_chunk_name(chunk_counter++);
            write_chunk(chunk_path, buffer, chunks);
        }

        merge_chunks(chunks, options.output_path, options.disable_binary_output);

        if (!options.keep_temporary) {
            for (const auto &chunk : chunks) {
                std::error_code ec;
                std::filesystem::remove(chunk.path, ec);
            }
        }

        return EXIT_SUCCESS;
    } catch (const std::exception &ex) {
        std::cerr << "logsort error: " << ex.what() << '\n';
        return EXIT_FAILURE;
    }
}
