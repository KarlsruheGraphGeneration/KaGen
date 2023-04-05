#pragma once

#include <string>

#include <mpi.h>

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/io/graph_format.h"

namespace kagen {
class SequentialGraphWriter : public GraphWriter {
protected:
    enum Requirement {
        NONE              = 0,
        SORTED_EDGES      = 1 << 1,
        COORDINATES       = 1 << 2,
        COORDINATES_2D    = 1 << 3,
        COORDINATES_3D    = 1 << 4,
        NO_VERTEX_WEIGHTS = 1 << 5,
        NO_EDGE_WEIGHTS   = 1 << 6,
    };

public:
    SequentialGraphWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm);

    virtual void Write(bool report_progress) override;

protected:
    virtual void AppendHeaderTo(const std::string& filename, SInt n, SInt m) = 0;

    virtual void AppendTo(const std::string& filename) = 0;

    virtual void AppendFooterTo(const std::string& filename);

    virtual int Requirements() const;

private:
    static void CreateFile(const std::string& filename);
};
} // namespace kagen
