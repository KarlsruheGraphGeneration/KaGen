/*******************************************************************************
 * include/io/generator_io.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Daniel Funke <funke@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <vector>

#include <fcntl.h>
#include <mpi.h>
#include <sys/stat.h>
#include <unistd.h>

#include "kagen/generator_config.h"

namespace kagen {
namespace internal {
template <std::size_t kBufferSize = 1024 * 1024, std::size_t kBufferSizeLimit = kBufferSize - 1024>
class BufferedTextOutput {
public:
    BufferedTextOutput(const std::string& filename)
        : _fd{open(filename.c_str(), O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR)} {
        if (_fd < 0) {
            std::cout << "cannot write to " << filename << std::endl;
            std::exit(0);
        }
    }

    ~BufferedTextOutput() {
        force_flush();
        close(_fd);
    }

    BufferedTextOutput& write_char(const char ch) {
        *(_buffer_pos)++ = ch;
        return *this;
    }

    template <typename Int>
    BufferedTextOutput& write_int(Int value) {
        static char rev_buffer[80];

        int pos = 0;
        do {
            rev_buffer[pos++] = value % 10;
            value /= 10;
        } while (value > 0);

        while (pos > 0) {
            *(_buffer_pos++) = '0' + rev_buffer[--pos];
        }
        return *this;
    }

    BufferedTextOutput& flush() {
        if (static_cast<std::size_t>(_buffer_pos - _buffer) >= kBufferSizeLimit) {
            force_flush();
        }
        return *this;
    }

private:
    void force_flush() {
        write(_fd, _buffer, _buffer_pos - _buffer);
        _buffer_pos = _buffer;
    }

    int   _fd;
    char  _buffer[kBufferSize]{0};
    char* _buffer_pos{_buffer};
};
} // namespace internal

class GeneratorIO {
public:
    using Edge = std::tuple<SInt, SInt>;

    GeneratorIO(PGeneratorConfig& config) : config_(config), local_num_edges_(0) {
        dist_.resize(config_.dist_size);
    }

    inline void UpdateDist(SInt node_id) {
        if (node_id < dist_.size()) {
            dist_[node_id]++;
        }
        local_num_edges_++;
    }

    void ReserveEdges(SInt num_edges) {
        edges_.reserve(num_edges);
    }

    template <typename... Args>
    inline void PushEdge(Args... args) {
        edges_.emplace_back(std::make_tuple(args...));
        local_num_edges_++;
    }

    auto& GetEdges() {
        return edges_;
    }

    SInt NumEdges() const {
        return edges_.size() > 0 ? edges_.size() : local_num_edges_ / 2;
    }

    void OutputDist() const {
        // Exchange local dist
        PEID rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        std::vector<SInt> global_dist(dist_.size(), 0);
        MPI_Reduce(&dist_[0], &global_dist[0], dist_.size(), MPI_LONG, MPI_SUM, ROOT, MPI_COMM_WORLD);
        if (rank == ROOT) {
            FILE* fout = fopen(config_.output_file.c_str(), "w+");
            for (SInt i = 0; i < global_dist.size(); ++i) {
                fprintf(fout, "%llu\n", global_dist[i]);
            }
            fclose(fout);
        }
    }

    void OutputEdges() const {
        if (config_.output_single_file) {
            GatherPrint();
        } else {
            Print();
        }
    }

private:
    PGeneratorConfig& config_;

    std::vector<SInt> dist_;
    std::vector<Edge> edges_;

    SInt local_num_edges_;

    void GatherPrint() const {
        // Exchange local dist
        PEID rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        // Gather number of edges for each PE
        std::vector<int> displ(size);
        std::vector<int> num_edges(size);
        int              lSize = NumEdges();
        MPI_Gather(&lSize, 1, MPI_INT, num_edges.data(), 1, MPI_INT, ROOT, MPI_COMM_WORLD);
        int current_displ   = 0;
        int total_num_edges = 0;
        if (rank == ROOT) {
            for (SInt i = 0; i < num_edges.size(); ++i) {
                displ[i] = current_displ;
                total_num_edges += num_edges[i];
                current_displ = total_num_edges;
            }
        }

        // Gather actual edges
        MPI_Datatype MPI_EDGE;
        MPI_Type_vector(1, 2, 0, MPI_LONG, &MPI_EDGE);
        MPI_Type_commit(&MPI_EDGE);
        std::vector<Edge> edges(total_num_edges);
        MPI_Gatherv(
            edges_.data(), lSize, MPI_EDGE, edges.data(), num_edges.data(), displ.data(), MPI_EDGE, ROOT,
            MPI_COMM_WORLD);

        if (rank == ROOT) {
            // Sort edges and remove duplicates
            std::sort(std::begin(edges), std::end(edges));
            // SInt total_edges = edges.size();
            edges.erase(unique(edges.begin(), edges.end()), edges.end());

            switch (config_.output_format) {
                case OutputFormat::BINARY_EDGE_LIST:
                    WriteBinaryEdgeList(config_.output_file, edges.size(), edges);
                    break;

                case OutputFormat::EDGE_LIST:
                    WriteEdgeList(config_.output_file, edges.size(), edges);
                    break;
            }
        }
    }

    // node id output
    void Print() const {
        PEID rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        SInt num_edges       = edges_.size();
        SInt total_num_edges = 0;
        MPI_Allreduce(&num_edges, &total_num_edges, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

        const std::string filename = config_.output_file + "_" + std::to_string(rank);
        switch (config_.output_format) {
            case OutputFormat::BINARY_EDGE_LIST:
                WriteBinaryEdgeList(filename, total_num_edges, edges_);
                break;

            case OutputFormat::EDGE_LIST:
                WriteEdgeList(filename, total_num_edges, edges_);
                break;
        }
    }

    void WriteBinaryEdgeList(const std::string& filename, const SInt total_m, const std::vector<Edge>& edges) const {
        FILE* fout = fopen(filename.c_str(), "wb+");
        if (!fout) {
            std::cerr << "Error: cannot write to " << filename << "\n";
            std::exit(1);
        }

        if (config_.output_header) {
            fwrite(&config_.n, sizeof(SInt), 1, fout);
            fwrite(&total_m, sizeof(SInt), 1, fout);
        }
        for (auto edge: edges) {
            const SInt source = std::get<0>(edge) + 1;
            const SInt target = std::get<1>(edge) + 1;
            fwrite(&source, sizeof(SInt), 1, fout);
            fwrite(&target, sizeof(SInt), 1, fout);
        }

        fclose(fout);
    }

    void WriteEdgeList(const std::string& filename, const SInt total_m, const std::vector<Edge>& edges) const {
        internal::BufferedTextOutput<> out(filename);
        if (config_.output_header) {
            out.write_char('p')
                .write_char(' ')
                .write_int(config_.n)
                .write_char(' ')
                .write_int(total_m)
                .write_char('\n')
                .flush();
        }
        for (const auto& edge: edges) {
            out.write_char('e')
                .write_char(' ')
                .write_int(std::get<0>(edge) + 1)
                .write_char(' ')
                .write_int(std::get<1>(edge) + 1)
                .write_char('\n')
                .flush();
        }
    }
};
} // namespace kagen
