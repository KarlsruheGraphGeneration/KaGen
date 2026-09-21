#pragma once

#include "kagen/kagen.h"

#include <exception>
#include <string>
#include <utility>

namespace kagen {
constexpr PEID        ROOT             = 0;
constexpr std::size_t NEWTON_MAX_ITERS = 10000;
constexpr double      NEWTON_EPS       = 0.001;

enum Direction {
    Up,
    Down,
    Left,
    Right,
    Front,
    Back,
};

class ConfigurationError : public std::exception {
public:
    ConfigurationError(std::string what) : _what(std::move(what)) {}

    const char* what() const noexcept override {
        return _what.c_str();
    }

private:
    std::string _what;
};
} // namespace kagen
