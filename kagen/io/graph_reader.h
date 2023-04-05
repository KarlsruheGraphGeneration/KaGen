#pragma once

#include <stdexcept>
#include <string>
#include <utility>

#include "kagen/definitions.h"

namespace kagen {
class IOError : public std::exception {
public:
    IOError(std::string what) : _what(std::move(what)) {}

    const char* what() const noexcept override {
        return _what.c_str();
    }

private:
    std::string _what;
};

using GraphSize = std::pair<SInt, SInt>;

class GraphReader {
public:
    virtual ~GraphReader() = default;

    virtual GraphSize ReadSize() = 0;

    virtual Graph Read(SInt from_from, SInt to_node, SInt to_edge, GraphRepresentation representation) = 0;

    virtual SInt FindNodeByEdge(SInt edge) = 0;
};
} // namespace kagen
