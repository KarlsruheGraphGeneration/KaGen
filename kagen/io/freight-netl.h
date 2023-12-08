#pragma once

#include "kagen/io/graph_format.h"

#include <mpi.h>

#include <string>

namespace kagen {
class FreightNetlWriter : public StandardGraphWriter {
public:
    FreightNetlWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size);

protected:
    void WriteHeader(const std::string& filename, SInt n, SInt m) final;

    bool WriteBody(const std::string& filename) final;
};

class FreightNetlFactory : public FileFormatFactory {
public:
    std::vector<std::string> DefaultExtensions() const final {
        return {"netl"};
    }

    std::unique_ptr<GraphWriter> CreateWriter(
        const OutputGraphConfig& config, Graph& graph, GraphInfo info, const PEID rank, const PEID size) const final;
};

class FreightNetlEpWriter : public StandardGraphWriter {
public:
    FreightNetlEpWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size);

protected:
    void WriteHeader(const std::string& filename, SInt n, SInt m) final;

    bool WriteBody(const std::string& filename) final;
};

class FreightNetlEpFactory : public FileFormatFactory {
public:
    std::vector<std::string> DefaultExtensions() const final {
        return {"ep.netl"};
    }

    std::unique_ptr<GraphWriter> CreateWriter(
        const OutputGraphConfig& config, Graph& graph, GraphInfo info, const PEID rank, const PEID size) const final;
};
} // namespace kagen
