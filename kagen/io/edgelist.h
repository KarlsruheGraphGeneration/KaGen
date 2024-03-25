#pragma once

#include "kagen/io/graph_format.h"
#include "kagen/io/mmap_toker.h"
#include "kagen/kagen.h"

#include <mpi.h>

#include <fstream>

namespace kagen {
//
// Text-based edge list 
//

class EdgelistWriter : public StandardGraphWriter {
public:
    EdgelistWriter(
        bool header, bool directed, const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank,
        PEID size);

protected:
    void WriteHeader(const std::string& filename, SInt n, SInt m) final;

    bool WriteBody(const std::string& filename) final;

private:
    bool header_;
    bool directed_;
};

class EdgelistFactory : public FileFormatFactory {
public:
    std::vector<std::string> DefaultExtensions() const final {
        return {"edgelist"};
    }

    std::unique_ptr<GraphWriter>
    CreateWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size) const final;
};

class UndirectedEdgelistFactory : public FileFormatFactory {
public:
    std::vector<std::string> DefaultExtensions() const final {
        return {"edgelist"};
    }

    std::unique_ptr<GraphWriter>
    CreateWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size) const final;
};

//
// Binary edge list
//

class BinaryEdgelistWriter : public StandardGraphWriter {
public:
    BinaryEdgelistWriter(
        bool header, bool directed, int datatype_size, const OutputGraphConfig& config, Graph& graph, GraphInfo info,
        PEID rank, PEID size);

protected:
    void WriteHeader(const std::string& filename, SInt n, SInt m) final;

    bool WriteBody(const std::string& filename) final;

private:
    bool header_;
    bool directed_;
    int  datatype_size_;
};

class BinaryEdgelistFactory : public FileFormatFactory {
public:
    std::vector<std::string> DefaultExtensions() const final {
        return {"binary-edgelist"};
    }

    std::unique_ptr<GraphWriter>
    CreateWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size) const final;
};

class UndirectedBinaryEdgelistFactory : public FileFormatFactory {
public:
    std::vector<std::string> DefaultExtensions() const final {
        return {"undirected-binary-edgelist"};
    }

    std::unique_ptr<GraphWriter>
    CreateWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size) const final;
};

//
// XtraPuLP edge list
//

class XtrapulpFactory : public FileFormatFactory {
public:
    std::vector<std::string> DefaultExtensions() const final {
        return {"xtrapulp"};
    }

    std::unique_ptr<GraphWriter>
    CreateWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size) const final;
};

// 
// Text edge list without any headers / line prefixes etc.
// Just lines 
// <source> <target>
//

class PlainEdgelistReader : public GraphReader {
public:
    PlainEdgelistReader(const std::string& filename);

    GraphSize ReadSize() final;

    Graph Read(SInt from_vertex, SInt to_vertex, SInt to_edge, GraphRepresentation representation) final;

    SInt FindNodeByEdge(SInt edge) final;

    int Deficits() const final;

private:
    MappedFileToker toker_;
};

class PlainEdgelistWriter : public StandardGraphWriter {
public:
    PlainEdgelistWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size);

protected:
    void WriteHeader(const std::string& filename, SInt n, SInt m) final;

    bool WriteBody(const std::string& filename) final;
};

class PlainEdgelistFactory : public FileFormatFactory {
public:
    std::vector<std::string> DefaultExtensions() const final {
        return {"plain-edgelist"};
    }

    std::unique_ptr<GraphReader> CreateReader(const InputGraphConfig& config, PEID rank, PEID size) const final;

    std::unique_ptr<GraphWriter>
    CreateWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size) const final;
};

//
// Binary edge list without header etc., but with edge weights:
// <source><target><weight>
//
// Ctor takes number of bits for vertex IDs / edge weights
//

class WeightedBinaryEdgelistWriter : public StandardGraphWriter {
public:
    WeightedBinaryEdgelistWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size);

protected:
    void WriteHeader(const std::string& filename, SInt n, SInt m) final;

    bool WriteBody(const std::string& filename) final;
};

class WeightedBinaryEdgelistReader : public GraphReader {
public:
    WeightedBinaryEdgelistReader(const std::string& filename, SInt vtx_width = 32, SInt adjwgt_width = 8);

    GraphSize ReadSize() final;

    Graph Read(SInt from_vertex, SInt to_vertex, SInt to_edge, GraphRepresentation representation) final;

    SInt FindNodeByEdge(SInt edge) final;

    int Deficits() const final;

private:
    SInt StepSize() const;

    std::ifstream in_;

    SInt vtx_width_;
    SInt adjwgt_width_;
};

class WeightedBinaryEdgelistFactory : public FileFormatFactory {
public:
    std::vector<std::string> DefaultExtensions() const final {
        return {"weighted-binary-edgelist"};
    }

    std::unique_ptr<GraphReader> CreateReader(const InputGraphConfig& config, PEID rank, PEID size) const final;

    std::unique_ptr<GraphWriter>
    CreateWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size) const final;
};
} // namespace kagen
