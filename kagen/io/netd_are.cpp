#include "kagen/io/netd_are.h"

#include "kagen/io.h"
#include "kagen/kagen.h"

#include <fstream>
#include <string>

namespace kagen {
namespace {
constexpr static bool kDebug = false;

inline std::string StripExtension(const std::string& filename) {
    const auto pos = filename.find_last_of('.');
    if (pos == std::string::npos) {
        return filename;
    }
    return filename.substr(0, pos);
}
} // namespace

std::unique_ptr<GraphReader> NetDAreFactory::CreateReader(const InputGraphConfig& config, const PEID rank, PEID) const {
    return std::make_unique<NetDAreReader>(config.filename, rank != 0);
}

NetDAreReader::NetDAreReader(const std::string& filename, const bool noop) : noop_(noop), toker_(filename) {
    const std::string basename     = StripExtension(filename);
    const std::string are_filename = basename + ".are";
    if (std::ifstream(are_filename).good()) {
        are_toker_.emplace(are_filename);
    }
}

GraphReader::GraphSize NetDAreReader::ReadSize() {
    toker_.SkipLine(); // first line is ignored
    num_pins_ = toker_.ScanUnsigned();
    toker_.ConsumeChar('\n');
    num_nets_ = toker_.ScanUnsigned();
    toker_.ConsumeChar('\n');
    num_modules_ = toker_.ScanUnsigned();
    toker_.ConsumeChar('\n');
    pad_offset_ = toker_.ScanUnsigned();
    toker_.ConsumeChar('\n');

    if constexpr (kDebug) {
        std::cout << "[Debug] netD header: #pins: " << num_pins_ << ", #nets: " << num_nets_
                  << ", #modules: " << num_modules_ << ", pad offset: " << pad_offset_ << std::endl;
    }

    return {1, 1};
}

Graph NetDAreReader::Read(SInt, SInt, SInt, GraphRepresentation) {
    if (noop_) {
        return {};
    }

    Graph graph;

    if (are_toker_) {
        graph.vertex_weights.resize(num_modules_);
        if constexpr (kDebug) {
            std::cout << "[Debug] Found associated *.are file: reading " << num_modules_ << " weights" << std::endl;
        }
        for (std::uint64_t mod = 0; mod < num_modules_; ++mod) {
            const char type   = are_toker_->ScanChar();
            const SInt id     = are_toker_->ScanUnsigned();
            const SInt weight = are_toker_->ScanUnsigned();
            are_toker_->ConsumeChar('\n');

            graph.vertex_weights[Remap(type, id)] = weight;
        }
        if (!are_toker_->ValidPosition() && kDebug) {
            std::cout << "[Debug] Finished reading *.are file, but have not reached EOF (next integer: "
                      << are_toker_->ScanUnsigned() << "; EOF': " << !are_toker_->ValidPosition() << ")" << std::endl;
        }
    }

    SInt                               flushed_nets = 0;
    std::vector<std::pair<char, SInt>> current_net;
    auto                               flush_current_net = [&] {
        for (const auto& [sdir, sid]: current_net) {
            for (const auto& [tdir, tid]: current_net) {
                if (sid == tid) {
                    continue;
                }
                if (sdir != tdir || (sdir == 'B' && tdir == 'B')) {
                    graph.edges.emplace_back(sid, tid);
                }
            }
        }
        current_net.clear();
        ++flushed_nets;
    };

    while (true) {
        const char type = toker_.ScanChar();
        if (type == EOF) {
            break;
        }

        const SInt id        = toker_.ScanUnsigned();
        const char marker    = toker_.ScanChar();
        const char direction = toker_.ScanChar();
        toker_.ConsumeChar('\n');

        if (marker != 's' && marker != 'l') {
            throw IOError(std::string("unexpected marker ") + marker + " in *.netD file");
        }
        if (direction != 'O' && direction != 'B' && direction != 'I') {
            throw IOError(std::string("unexpected direction ") + direction + " in *.netD file");
        }

        if (marker == 's') {
            flush_current_net();
        }
        current_net.emplace_back(direction, Remap(type, id));
    }
    flush_current_net();

    if constexpr (kDebug) {
        std::cout << "[Debug] Read " << graph.edges.size() << " edges and " << graph.vertex_weights.size()
                  << " vertex weights" << std::endl;
        std::cout << "[Debug] Number of nets processed: " << flushed_nets << std::endl;
    }

    graph.vertex_range   = {0, num_modules_};
    graph.representation = GraphRepresentation::EDGE_LIST;

    return graph;
}

SInt NetDAreReader::Remap(const char type, const SInt id) const {
    if (type != 'p' && type != 'a') {
        throw IOError(std::string("unexpected module type ") + type + " in *.netD file");
    }
    return type == 'p' ? RemapPad(id) : RemapCell(id);
}

SInt NetDAreReader::RemapPad(const SInt pad) const {
    if (pad < 1 || pad > num_modules_ - pad_offset_ - 1) {
        throw IOError(
            std::string("out-of-bounds pad ID ") + std::to_string(pad) + " (expected range: [1, "
            + std::to_string(num_modules_ - pad_offset_ - 1) + "])");
    }
    return pad + pad_offset_;
}

SInt NetDAreReader::RemapCell(const SInt cell) const {
    if (cell > pad_offset_) {
        throw IOError(
            std::string("out-of-bounds cell ID ") + std::to_string(cell) + " (expected range: [0, "
            + std::to_string(pad_offset_) + "])");
    }
    return cell;
}

SInt NetDAreReader::FindNodeByEdge(SInt) {
    return 0;
}

int NetDAreReader::Deficits() const {
    return ReaderDeficits::REQUIRES_REDISTRIBUTION | ReaderDeficits::EDGE_LIST_ONLY
           | ReaderDeficits::UNKNOWN_NUM_VERTICES | ReaderDeficits::UNKNOWN_NUM_EDGES;
}
} // namespace kagen
