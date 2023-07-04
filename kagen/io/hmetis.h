#pragma once

#include <string>

#include <mpi.h>

#include "kagen/io/graph_format.h"

namespace kagen {
class HmetisWriter : public StandardGraphWriter {
public:
    HmetisWriter(bool directed, const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size);

protected:
    void WriteHeader(const std::string& filename, SInt n, SInt m) final;

    bool WriteBody(const std::string& filename) final;

    void WriteFooter(const std::string& filename) final;

private:
    bool directed_;
};

class HmetisFactory : public FileFormatFactory {
public:
    std::string DefaultExtension() const final {
        return "hmetis";
    }

    std::unique_ptr<GraphWriter> CreateWriter(
        const OutputGraphConfig& config, Graph& graph, GraphInfo info, const PEID rank, const PEID size) const final;
};

class DirectedHmetisFactory : public FileFormatFactory {
public:
    std::string DefaultExtension() const final {
        return "hmetis";
    }

    std::unique_ptr<GraphWriter> CreateWriter(
        const OutputGraphConfig& config, Graph& graph, GraphInfo info, const PEID rank, const PEID size) const final;
};
} // namespace kagen
