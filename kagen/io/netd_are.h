#pragma once

#include "kagen/io/graph_format.h"
#include "kagen/io/mmap_toker.h"
#include "kagen/kagen.h"

#include <mpi.h>

#include <optional>

namespace kagen {
class NetDAreFactory : public FileFormatFactory {
public:
    std::vector<std::string> DefaultExtensions() const final {
        return {"netD"};
    }

    std::unique_ptr<GraphReader> CreateReader(const InputGraphConfig& config, PEID rank, PEID size) const final;
};

class NetDAreReader : public GraphReader {
public:
    NetDAreReader(const std::string& filename, bool noop);

    GraphSize ReadSize() final;

    Graph Read(SInt from_vertex, SInt to_vertex, SInt to_edge, GraphRepresentation representation) final;

    SInt FindNodeByEdge(SInt edge) final;

    int Deficits() const final;

private:
    SInt Remap(char type, SInt id) const;
    SInt RemapPad(SInt pad) const;
    SInt RemapCell(SInt cell) const;

    bool                           noop_;
    MappedFileToker                toker_;
    std::optional<MappedFileToker> are_toker_;

    std::uint64_t num_pins_;
    std::uint64_t num_nets_;
    std::uint64_t num_modules_;
    std::uint64_t pad_offset_;
};
} // namespace kagen

