#include "kagen/io/graph_format.h"


namespace kagen{
class HmetisHypergraphWriter : public StandardGraphWriter {
public:
    HmetisHypergraphWriter(const OutputGraphConfig& config,
        Graph& graph,
        GraphInfo info,
        PEID rank,
        PEID size
    ) : StandardGraphWriter(config, graph, info, rank, size) {}

protected:
    void WriteHeader(const std::string& filename, SInt n, SInt m) override;
    bool WriteBody(const std::string& filename) override;
    void WriteFooter(const std::string& filename) override;
};

class HmetisHypergraphFactory : public FileFormatFactory {
public:
    std::vector<std::string> DefaultExtensions() const final {
        return {"hmetisHyper"};
    }

    std::unique_ptr<GraphWriter> CreateWriter(const OutputGraphConfig& config, Graph& graph, GraphInfo info, PEID rank, PEID size) const final;
};
}