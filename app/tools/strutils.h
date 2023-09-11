#pragma once

#include <string>

namespace kagen {
inline std::string ExtractFilename(const std::string& filename) {
    const auto pos = filename.find_last_of('/');
    if (pos == std::string::npos) {
        return filename;
    }
    return filename.substr(pos + 1);
}

inline std::string StripExtension(const std::string& filename) {
    const auto pos = filename.find_last_of('.');
    if (pos == std::string::npos) {
        return filename;
    }
    return filename.substr(0, pos);
}
} // namespace kagen
