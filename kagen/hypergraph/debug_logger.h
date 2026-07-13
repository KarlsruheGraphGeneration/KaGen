#pragma once

#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>

namespace kagen {

class CsvDebugLogger {
public:
    CsvDebugLogger(const std::string& filename, const bool append, std::string_view header)
        : out_(filename, append ? (std::ios::out | std::ios::app) : (std::ios::out | std::ios::trunc)) {
        const bool should_write_header =
            !append || !std::filesystem::exists(filename) || std::filesystem::file_size(filename) == 0;

        if (should_write_header) {
            out_ << header << '\n';
        }
    }

protected:
    std::ofstream& out() {
        return out_;
    }

private:
    std::ofstream out_;
};

} // namespace kagen