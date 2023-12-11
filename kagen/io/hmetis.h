#pragma once

#include "kagen/io/graph_format.h"

#include <mpi.h>

#include <string>

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
    std::vector<std::string> DefaultExtensions() const final {
        return {"hmetis"};
    }

    std::unique_ptr<GraphWriter> CreateWriter(
        const OutputGraphConfig& config, Graph& graph, GraphInfo info, const PEID rank, const PEID size) const final;
};

class DirectedHmetisFactory : public FileFormatFactory {
public:
    std::vector<std::string> DefaultExtensions() const final {
        return {"hmetis"};
    }

    std::unique_ptr<GraphWriter> CreateWriter(
        const OutputGraphConfig& config, Graph& graph, GraphInfo info, const PEID rank, const PEID size) const final;
};

class HmetisEpWriter : public StandardGraphWriter {
public:
    HmetisEpWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size);

protected:
    void WriteHeader(const std::string& filename, SInt n, SInt m) final;

    bool WriteBody(const std::string& filename) final;
};

class HmetisEpFactory : public FileFormatFactory {
public:
    std::vector<std::string> DefaultExtensions() const final {
        return {"ep.hmetis"};
    }

    std::unique_ptr<GraphWriter> CreateWriter(
        const OutputGraphConfig& config, Graph& graph, GraphInfo info, const PEID rank, const PEID size) const final;
};
} // namespace kagen
