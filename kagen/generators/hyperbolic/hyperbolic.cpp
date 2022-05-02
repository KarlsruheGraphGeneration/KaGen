#include "kagen/generators/hyperbolic/hyperbolic.h"

#include "kagen/generators/generator.h"

namespace kagen {
Hyperbolic::Hyperbolic(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size),
      rng_(config) {
    // Globals
    alpha_         = (config_.plexp - 1) / 2;
    target_r_      = PGGeometry::GetTargetRadius(config_.n, config_.n * config_.avg_degree / 2, alpha_);
    cosh_target_r_ = cosh(target_r_);
    pdm_target_r_  = (cosh_target_r_ - 1) / 2;
    // clique_thres_ = target_r_ / 2.0;
    clique_thres_ = 0;

    // PE-specific
    total_annuli_        = floor(alpha_ * target_r_ / log(2));
    SInt chunks_per_pe   = config_.k / size_;
    SInt leftover_chunks = config_.k % size_;
    local_chunks_        = chunks_per_pe + ((SInt)rank_ < leftover_chunks);

    // Compute local chunk range
    local_chunk_start_ = local_chunks_ * rank_ + ((SInt)rank_ >= leftover_chunks ? leftover_chunks : 0);
    local_chunk_end_   = local_chunk_start_ + local_chunks_;

    LPFloat phi_per_chunk = 2 * M_PI / config_.k;
    pe_min_phi_           = local_chunk_start_ * phi_per_chunk;
    pe_max_phi_           = local_chunk_end_ * phi_per_chunk;

    // Init data structures
    annuli_.set_empty_key(total_annuli_ * config_.k);
    chunks_.set_empty_key(config_.k);
    boundaries_.resize(total_annuli_);

    // Compute number of cells_
    SInt total_cells = 0;
    cells_per_annulus_.resize(total_annuli_, std::numeric_limits<SInt>::max());
    for (SInt i = 0; i < total_annuli_; ++i) {
        global_cell_ids_.emplace_back(total_cells);
        total_cells += GridSizeForAnnulus(i) * size_;
    }
    cells_.set_empty_key(total_cells + 1);
    vertices_.set_empty_key(total_cells + 1);

    // Epsilon comparison
    chunk_eps_ = phi_per_chunk / 1000;
    cell_eps_  = (2 * M_PI / GridSizeForAnnulus(total_annuli_ - 1)) / 1000;
    point_eps_ = std::numeric_limits<LPFloat>::epsilon();

    // Vertex range
    start_node_ = std::numeric_limits<SInt>::max();
    num_nodes_  = 0;
}

int Hyperbolic::Requirements() const {
    return GeneratorRequirement::POWER_OF_TWO_COMMUNICATOR_SIZE;
}

bool Hyperbolic::AlmostUndirected() const {
    return true;
}

void Hyperbolic::GenerateImpl() {
    // Compute local chunks
    for (SInt i = local_chunk_start_; i < local_chunk_end_; ++i) {
        ComputeChunk(i);
        ComputeAnnuli(i);
    }

    // if (rank_ == ROOT)
    //   std::cout << "computed chunks" << std::endl;

    // Local points
    for (SInt i = local_chunk_start_; i < local_chunk_end_; ++i) {
        for (SInt j = 0; j < total_annuli_; ++j) {
            // if (rank_ == ROOT)
            //   std::cout << "gen cells " << j << " " << i << std::endl;
            GenerateCells(j, i);
            for (SInt k = 0; k < GridSizeForAnnulus(j); ++k)
                GenerateVertices(j, i, k);
        }
    }

    // if (rank_ == ROOT)
    //   std::cout << "generated vertices" << std::endl;

    // Local edges
    for (SInt i = local_chunk_start_; i < local_chunk_end_; ++i) {
        for (SInt j = 0; j < total_annuli_; ++j) {
            GenerateEdges(j, i);
        }
    }

    SetVertexRange(start_node_, start_node_ + num_nodes_);
}

void Hyperbolic::ComputeAnnuli(const SInt chunk_id) {
    SInt n      = std::get<0>(chunks_[chunk_id]);
    SInt offset = std::get<3>(chunks_[chunk_id]);

    LPFloat min_r      = 0;
    LPFloat total_area = PGGeometry::RadiusToHyperbolicArea(alpha_ * target_r_);

    for (SInt i = 1; i < total_annuli_ + 1; i++) {
        // Distribute points
        LPFloat max_r = i * target_r_ / total_annuli_;
        LPFloat ring_area =
            PGGeometry::RadiusToHyperbolicArea(alpha_ * max_r) - PGGeometry::RadiusToHyperbolicArea(alpha_ * min_r);

        // Variate
        SInt h = sampling::Spooky::hash(config_.seed + total_annuli_ * config_.k + chunk_id * total_annuli_ + i);
        SInt n_annulus = rng_.GenerateBinomial(h, n, ring_area / total_area);

        // Push annuli_
        annuli_[ComputeGlobalChunkId(i - 1, chunk_id)] = std::make_tuple(n_annulus, min_r, max_r, false, offset);
        boundaries_[i - 1]                             = std::make_pair(cosh(min_r), sinh(min_r));
        // if (rank_ == 2)
        //   printf("a %llu %f %f %f %f %llu\n", n_annulus, min_r, max_r, std::get<1>(chunks_[chunk_id]),
        //   std::get<2>(chunks_[chunk_id]), offset);
        min_r = max_r;
        n -= n_annulus;
        offset += n_annulus;
        total_area -= ring_area;
    }
    // if (config_.thres > 0)
    //   clique_thres_ = std::get<2>(annuli_[config_.thres]);
}

void Hyperbolic::ComputeChunk(const SInt chunk_id) {
    ComputeChunk(chunk_id, config_.n, config_.k, 0, 2 * M_PI, 0, 1, 0);
}

void Hyperbolic::ComputeChunk(
    const SInt chunk_id, const SInt n, const SInt k, const LPFloat min_phi, const LPFloat max_phi,
    const SInt chunk_start, const SInt level, const SInt offset) {
    // Base case
    if (k == 1) {
        chunks_[chunk_id] = std::make_tuple(n, min_phi, max_phi, offset);
        // if (rank_ == 2)
        //   printf("c %llu %f %f %llu\n", n, min_phi, max_phi, offset);
        return;
    }

    // Compute point distribution
    SInt midk = (k + 1) / 2;

    // Generate variate
    SInt h                = sampling::Spooky::hash(config_.seed + level * config_.k + chunk_start);
    SInt splitter_variate = rng_.GenerateBinomial(h, n, (LPFloat)midk / k);

    // Compute splitter
    LPFloat middlePhi = (max_phi - min_phi) * ((LPFloat)midk / k) + min_phi;
    // Manuel fix
    if (-1e-8 < middlePhi && middlePhi <= 0.0)
        middlePhi = 0;

    // Recurse
    if (chunk_id < chunk_start + midk)
        ComputeChunk(chunk_id, splitter_variate, midk, min_phi, middlePhi, chunk_start, level + 1, offset);
    else
        ComputeChunk(
            chunk_id, n - splitter_variate, k - midk, middlePhi, max_phi, chunk_start + midk, level + 1,
            offset + splitter_variate);
}

void Hyperbolic::GenerateCells(const SInt annulus_id, SInt chunk_id) {
    bool clique = false;
    // auto &annulus = annuli_[ComputeGlobalChunkId(annulus_id, chunk_id)];
    // if (std::get<1>(annulus) < clique_thres_) {
    //   chunk_id = local_chunk_start_;
    //   clique = true;
    // }

    // Lazily compute chunk
    if (chunks_.find(chunk_id) == end(chunks_)) {
        ComputeChunk(chunk_id);
        ComputeAnnuli(chunk_id);
    }
    auto& chunk   = chunks_[chunk_id];
    auto& annulus = annuli_[ComputeGlobalChunkId(annulus_id, chunk_id)];

    // Stop if cell distribution already generated
    if (std::get<3>(annulus))
        return;

    SInt    n, offset, seed = 0;
    LPFloat min_phi, max_phi;

    // Retrieve parameters
    if (clique) {
        n       = std::get<0>(annulus);
        offset  = std::get<3>(annulus);
        min_phi = 0.0;
        max_phi = 2 * M_PI;
        seed    = config_.seed + annulus_id + config_.n;
    } else {
        n       = std::get<0>(annulus);
        offset  = std::get<4>(annulus);
        min_phi = std::get<1>(chunk);
        max_phi = std::get<2>(chunk);
    }
    // Manuel fix
    if (-1e-8 < min_phi && min_phi <= 0.0)
        min_phi = 0;

    LPFloat total_phi = max_phi - min_phi;
    LPFloat grid_phi  = total_phi / GridSizeForAnnulus(annulus_id);
    // if (rank_ == ROOT)
    //   std::cout << "grid size " << GridSizeForAnnulus(annulus_id) << std::endl;
    for (SInt i = 0; i < GridSizeForAnnulus(annulus_id); ++i) {
        // Variate
        if (!clique)
            seed = config_.seed + annulus_id * config_.k + chunk_id + i + config_.n;
        SInt h      = sampling::Spooky::hash(seed);
        SInt n_cell = rng_.GenerateBinomial(h, n, grid_phi / total_phi);

        SInt global_cell_id = ComputeGlobalCellId(annulus_id, chunk_id, i);
        cells_[global_cell_id] =
            std::make_tuple(n_cell, min_phi + (grid_phi * i), min_phi + (grid_phi * (i + 1)), false, offset);
        // if (rank_ == 2)
        //   printf("g %llu %f %f %f %f %llu\n", n_cell, std::get<1>(annulus), std::get<2>(annulus), min_phi +
        //   (grid_phi * i), min_phi + (grid_phi * (i + 1)), offset);
        n -= n_cell;
        offset += n_cell;
        total_phi -= grid_phi;
    }
    std::get<3>(annulus) = true;
}

void Hyperbolic::GenerateVertices(const SInt annulus_id, SInt chunk_id, const SInt cell_id) {
    bool clique = false;
    // auto &annulus = annuli_[annulus_id];
    // if (std::get<1>(annulus) < clique_thres_) {
    //   chunk_id = local_chunk_start_;
    //   clique = true;
    // }

    // Lazily compute chunk
    if (chunks_.find(chunk_id) == end(chunks_)) {
        ComputeChunk(chunk_id);
        ComputeAnnuli(chunk_id);
    }
    // auto &chunk = chunks_[chunk_id];
    auto& annulus = annuli_[ComputeGlobalChunkId(annulus_id, chunk_id)];

    // Lazily compute cells distribution
    if (!std::get<3>(annulus))
        GenerateCells(annulus_id, chunk_id);

    // Check if cell was generated
    SInt  global_cell_id = ComputeGlobalCellId(annulus_id, chunk_id, cell_id);
    auto& cell           = cells_[global_cell_id];
    if (std::get<3>(cell))
        return;

    // Compute vertex distribution
    SInt    n       = std::get<0>(cell);
    SInt    offset  = std::get<4>(cell);
    LPFloat min_phi = std::get<1>(cell);
    LPFloat max_phi = std::get<2>(cell);
    LPFloat min_r   = std::get<1>(annulus);
    LPFloat max_r   = std::get<2>(annulus);

    SInt seed = 0;
    if (clique)
        seed = config_.seed + annulus_id * config_.k * GridSizeForAnnulus(annulus_id);
    else
        seed = config_.seed + annulus_id * config_.k * GridSizeForAnnulus(annulus_id)
               + chunk_id * GridSizeForAnnulus(annulus_id) + cell_id + config_.n;

    SInt h = sampling::Spooky::hash(seed);
    mersenne.RandomInit(h);
    sorted_mersenne.RandomInit(h, n);
    LPFloat              mincdf        = cosh(alpha_ * min_r);
    LPFloat              maxcdf        = cosh(alpha_ * max_r);
    std::vector<Vertex>& cell_vertices = vertices_[global_cell_id];
    cell_vertices.reserve(n);
    for (SInt i = 0; i < n; i++) {
        // Compute coordinates
        LPFloat angle  = sorted_mersenne.Random() * (max_phi - min_phi) + min_phi;
        LPFloat radius = acosh(mersenne.Random() * (maxcdf - mincdf) + mincdf) / alpha_;

        // Perform pdm transformation
        LPFloat inv_len    = (cosh(radius) + 1.0) / 2.0;
        LPFloat pdm_radius = sqrt(1.0 - 1.0 / inv_len);
        LPFloat x          = pdm_radius * sin(angle);
        LPFloat y          = pdm_radius * cos(angle);
        LPFloat gamma      = 1.0 / (1.0 - pdm_radius * pdm_radius);
        cell_vertices.emplace_back(angle, radius, x, y, gamma, offset + i);
        // if (rank_ == 2)
        //   printf("p %lld %f %f %d\n", offset + i, radius, angle, rank_);
        if ((start_node_ > offset + i) && (pe_min_phi_ <= angle && pe_max_phi_ > angle))
            start_node_ = offset;
        if (pe_min_phi_ <= angle && pe_max_phi_ > angle)
            num_nodes_++;

        if (config_.coordinates) {
            PushCoordinate(x, y);
        }
    }
    std::get<3>(cell) = true;
}

void Hyperbolic::GenerateEdges(const SInt annulus_id, const SInt chunk_id) {
    current_annulus_ = annulus_id;
    current_chunk_   = chunk_id;
    for (SInt cell_id = 0; cell_id < GridSizeForAnnulus(annulus_id); ++cell_id) {
        SInt global_cell_id = ComputeGlobalCellId(annulus_id, chunk_id, cell_id);
        if (std::get<0>(cells_[global_cell_id]) == 0)
            continue;
        current_cell_ = cell_id;
        for (SInt i = 0; i < vertices_[global_cell_id].size(); ++i) {
            // const Vertex &v = cell_vertices[i];
            // Need copy because of hash movement
            const Vertex v = vertices_[global_cell_id][i];
            if (pe_min_phi_ > std::get<0>(v) || pe_max_phi_ < std::get<0>(v))
                continue;
            // if (rank_ == 2)
            //   printf("qp %lld %f %f %f\n", std::get<5>(v), std::get<1>(v), std::get<0>(v), target_r_);
            QueryBoth(annulus_id, chunk_id, cell_id, v);
        }
    }
}

void Hyperbolic::QueryBoth(const SInt annulus_id, const SInt chunk_id, const SInt cell_id, const Vertex& q) {
    Query(annulus_id, chunk_id, cell_id, q);
    if (config_.query_both && annulus_id > 0) {
        auto&   chunk         = chunks_[chunk_id];
        LPFloat min_chunk_phi = std::get<1>(chunk);
        LPFloat max_chunk_phi = std::get<2>(chunk);
        LPFloat grid_phi      = (max_chunk_phi - min_chunk_phi) / GridSizeForAnnulus(annulus_id - 1);
        SInt    next_cell_id  = floor((std::get<0>(q) - min_chunk_phi) / grid_phi);
        Query(annulus_id - 1, chunk_id, next_cell_id, q, false);
    }
}

void Hyperbolic::Query(
    const SInt annulus_id, const SInt chunk_id, const SInt cell_id, const Vertex& q, bool search_down) {
    // Boundaries
    auto& annulus        = annuli_[ComputeGlobalChunkId(annulus_id, chunk_id)];
    auto  current_bounds = GetBoundaryPhis(std::get<0>(q), std::get<1>(q), annulus_id);
    current_min_phi_     = std::get<0>(current_bounds);
    current_max_phi_     = std::get<1>(current_bounds);

    LPFloat min_cell_phi = std::get<1>(cells_[ComputeGlobalCellId(annulus_id, chunk_id, cell_id)]);
    LPFloat max_cell_phi = std::get<2>(cells_[ComputeGlobalCellId(annulus_id, chunk_id, cell_id)]);

    // if (rank_ == 2) {
    //   std::cout << "go down " << chunk_id << " " << annulus_id << " " << cell_id << std::endl;
    //   std::cout << "min phi " << current_min_phi_ << " max phi " << current_max_phi_ << std::endl;
    //   std::cout << "min cell phi " << min_cell_phi << " max cell phi " << max_cell_phi << std::endl;
    // }

    right_processed_chunk_ = chunk_id;
    right_processed_cell_  = cell_id;

    // Iterate over cell
    if (search_down || !IsLocalChunk(chunk_id))
        GenerateGridEdges(annulus_id, chunk_id, cell_id, q);

    if (std::get<1>(annulus) >= clique_thres_ && std::max(TotalGridSizeForAnnulus(annulus_id), config_.k) > 1) {
        // Continue right
        if (current_min_phi_ < min_cell_phi
            || (OutOfBounds(current_min_phi_) && !(std::abs(min_cell_phi - 0.0) < cell_eps_))) {
            SInt next_chunk_id = chunk_id;
            if (cell_id == 0)
                next_chunk_id = (chunk_id + config_.k - 1) % config_.k;
            SInt next_cell_id = (cell_id + GridSizeForAnnulus(annulus_id) - 1) % GridSizeForAnnulus(annulus_id);
            // if (rank_ == 2)
            //   std::cout << "go right " << next_chunk_id << " " << annulus_id << " " << next_cell_id << std::endl;
            GenerateVertices(annulus_id, next_chunk_id, next_cell_id);
            QueryRightNeighbor(
                annulus_id, next_chunk_id, next_cell_id, q, std::abs(min_cell_phi - 0.0) < cell_eps_, search_down);
        }

        // Continue left
        if (current_max_phi_ > max_cell_phi
            || (OutOfBounds(current_max_phi_) && !(std::abs(max_cell_phi - 2 * M_PI) < cell_eps_))) {
            SInt next_chunk_id = chunk_id;
            if (cell_id == GridSizeForAnnulus(annulus_id) - 1)
                next_chunk_id = (chunk_id + config_.k + 1) % config_.k;
            SInt next_cell_id = (cell_id + GridSizeForAnnulus(annulus_id) + 1) % GridSizeForAnnulus(annulus_id);
            // if (rank_ == 2)
            //   std::cout << "go left " << next_chunk_id << " " << annulus_id << " " << next_cell_id << std::endl;
            GenerateVertices(annulus_id, next_chunk_id, next_cell_id);
            QueryLeftNeighbor(
                annulus_id, next_chunk_id, next_cell_id, q, std::abs(max_cell_phi - 2 * M_PI) < cell_eps_, search_down);
        }
    }

    // Continue with next annuli_
    SInt next_annulus;
    if (search_down)
        next_annulus = annulus_id + 1;
    else
        next_annulus = annulus_id - 1;

    if (next_annulus >= total_annuli_ || (LONG)next_annulus < 0)
        return;

    // Find next cell
    auto&   chunk         = chunks_[chunk_id];
    LPFloat min_chunk_phi = std::get<1>(chunk);
    LPFloat max_chunk_phi = std::get<2>(chunk);
    LPFloat grid_phi      = (max_chunk_phi - min_chunk_phi) / GridSizeForAnnulus(next_annulus);
    SInt    next_cell_id  = floor((std::get<0>(q) - min_chunk_phi) / grid_phi);

    Query(next_annulus, chunk_id, next_cell_id, q, search_down);
}

void Hyperbolic::QueryRightNeighbor(
    const SInt annulus_id, SInt chunk_id, SInt cell_id, const Vertex& q, bool phase, bool search_down) {
    while (true) {
        // Boundaries
        if (phase && current_min_phi_ < 0.0)
            current_min_phi_ += 2 * M_PI;
        if (phase && OutOfBounds(current_min_phi_))
            return;
        // || std::get<1>(annuli_[annulus_id]) < clique_thres_))

        auto&   cell           = cells_[ComputeGlobalCellId(annulus_id, chunk_id, cell_id)];
        LPFloat min_cell_phi   = std::get<1>(cell);
        right_processed_chunk_ = chunk_id;
        right_processed_cell_  = cell_id;

        // Iterate over cell
        if ((search_down && IsLocalChunk(chunk_id) && min_cell_phi > std::get<0>(q)) || !IsLocalChunk(chunk_id))
            GenerateGridEdges(annulus_id, chunk_id, cell_id, q);

        phase = phase || std::abs(min_cell_phi - 0.0) < cell_eps_;
        if (current_min_phi_ < min_cell_phi || OutOfBounds(current_min_phi_)) {
            SInt next_chunk_id = chunk_id;
            if (cell_id == 0)
                next_chunk_id = (chunk_id + config_.k - 1) % config_.k;
            SInt next_cell_id = (cell_id + GridSizeForAnnulus(annulus_id) - 1) % GridSizeForAnnulus(annulus_id);
            // if (rank_ == 2)
            //   std::cout << "go right " << next_chunk_id << " " << annulus_id << " " << next_cell_id << std::endl;
            GenerateVertices(annulus_id, next_chunk_id, next_cell_id);
            cell_id  = next_cell_id;
            chunk_id = next_chunk_id;
            continue;
        }
        return;
    }
}

void Hyperbolic::QueryLeftNeighbor(
    const SInt annulus_id, SInt chunk_id, SInt cell_id, const Vertex& q, bool phase, bool search_down) {
    while (true) {
        // Boundaries
        if (phase && current_max_phi_ >= 2 * M_PI)
            current_max_phi_ -= 2 * M_PI;
        if (phase && OutOfBounds(current_max_phi_))
            return;
        // || std::get<1>(annuli_[annulus_id]) < clique_thres_))
        if (chunk_id == right_processed_chunk_ && cell_id == right_processed_cell_)
            return;

        auto&   cell         = cells_[ComputeGlobalCellId(annulus_id, chunk_id, cell_id)];
        LPFloat max_cell_phi = std::get<2>(cell);

        // Iterate over cell
        if (search_down || !IsLocalChunk(chunk_id))
            GenerateGridEdges(annulus_id, chunk_id, cell_id, q);

        phase = phase || std::abs(max_cell_phi - 2 * M_PI) < cell_eps_;
        if (current_max_phi_ > max_cell_phi || OutOfBounds(current_max_phi_)) {
            SInt next_chunk_id = chunk_id;
            if (cell_id == GridSizeForAnnulus(annulus_id) - 1)
                next_chunk_id = (chunk_id + config_.k + 1) % config_.k;
            SInt next_cell_id = (cell_id + GridSizeForAnnulus(annulus_id) + 1) % GridSizeForAnnulus(annulus_id);
            // if (rank_ == 2)
            //   std::cout << "go left " << next_chunk_id << " " << annulus_id << " " << next_cell_id << std::endl;
            GenerateVertices(annulus_id, next_chunk_id, next_cell_id);
            cell_id  = next_cell_id;
            chunk_id = next_chunk_id;
            continue;
        }
        return;
    }
}

void Hyperbolic::GenerateGridEdges(const SInt annulus_id, const SInt chunk_id, const SInt cell_id, const Vertex& q) {
    // Check if vertices not generated
    SInt global_cell_id = ComputeGlobalCellId(annulus_id, chunk_id, cell_id);
    GenerateVertices(annulus_id, chunk_id, cell_id);

    // Gather vertices
    const std::vector<Vertex>& cell_vertices = vertices_[global_cell_id];
    // Same cell
    if (current_annulus_ == annulus_id && current_chunk_ == chunk_id && current_cell_ == cell_id) {
        for (SInt j = 0; j < cell_vertices.size(); ++j) {
            const Vertex& v = cell_vertices[j];
            // Skip if larger angle or same angle and larger radius
            if (std::get<0>(v) > std::get<0>(q)
                || (std::abs(std::get<0>(v) - std::get<0>(q)) < point_eps_ && std::get<1>(v) < std::get<1>(q)))
                continue;
            if (std::abs(std::get<1>(v) - std::get<1>(q)) < point_eps_
                && std::abs(std::get<0>(v) - std::get<0>(q)) < point_eps_)
                continue;
            // Generate edge
            if (PGGeometry::HyperbolicDistance(q, v) <= pdm_target_r_) {
                PushEdge(std::get<5>(q), std::get<5>(v));
                PushEdge(std::get<5>(v), std::get<5>(q));
            }
        }
    }
    // Different cells
    else {
        for (SInt j = 0; j < cell_vertices.size(); ++j) {
            const Vertex& v = cell_vertices[j];
            if (PGGeometry::HyperbolicDistance(q, v) <= pdm_target_r_) {
                PushEdge(std::get<5>(q), std::get<5>(v));
                if (IsLocalChunk(chunk_id)) {
                    PushEdge(std::get<5>(v), std::get<5>(q));
                }
            }
        }
    }
}

inline std::pair<LPFloat, LPFloat>
Hyperbolic::GetBoundaryPhis(const LPFloat boundary_phi, const LPFloat boundary_r, const SInt annulus_id) const {
    auto&   boundary    = boundaries_[annulus_id];
    LPFloat cosh_min_r  = std::get<0>(boundary);
    LPFloat sinh_min_r  = std::get<1>(boundary);
    LPFloat diff        = acos((cosh(boundary_r) * cosh_min_r - cosh_target_r_) / (sinh(boundary_r) * sinh_min_r));
    LPFloat lower_bound = boundary_phi - diff;
    LPFloat upper_bound = boundary_phi + diff;
    LPFloat min_phi     = std::min(lower_bound, upper_bound);
    LPFloat max_phi     = std::max(lower_bound, upper_bound);
    return std::make_pair(min_phi, max_phi);
}

inline bool Hyperbolic::OutOfBounds(const LPFloat num) const {
    return (std::isnan(num) || num < -2 * M_PI || num > 2 * M_PI);
}

inline SInt Hyperbolic::ComputeGlobalChunkId(const SInt annulus, const SInt chunk) const {
    // return annulus * config_.k + chunk;
    return chunk * total_annuli_ + annulus;
}

inline SInt Hyperbolic::ComputeGlobalCellId(const SInt annulus, const SInt chunk, const SInt cell) {
    return global_cell_ids_[annulus] + chunk * GridSizeForAnnulus(annulus) + cell;
}

SInt Hyperbolic::TotalGridSizeForAnnulus(const SInt annulus_id) {
    if (cells_per_annulus_[annulus_id] != std::numeric_limits<SInt>::max())
        return cells_per_annulus_[annulus_id];
    LPFloat min_r = annulus_id * target_r_ / total_annuli_;
    LPFloat max_r = (annulus_id + 1) * target_r_ / total_annuli_;
    LPFloat ring_area =
        PGGeometry::RadiusToHyperbolicArea(alpha_ * max_r) - PGGeometry::RadiusToHyperbolicArea(alpha_ * min_r);
    LPFloat total_area = PGGeometry::RadiusToHyperbolicArea(alpha_ * target_r_);

    SInt exp_points = config_.n * (LPFloat)ring_area / total_area;
    SInt cells      = exp_points / config_.hyp_base;

    SInt result                    = std::max<SInt>(1, cells);
    cells_per_annulus_[annulus_id] = result;
    return result;
}

inline SInt Hyperbolic::GridSizeForAnnulus(const SInt annulus_id) {
    return std::max<SInt>(1, TotalGridSizeForAnnulus(annulus_id) / size_);
}

inline bool Hyperbolic::IsLocalChunk(const SInt chunk_id) const {
    return (chunk_id >= local_chunk_start_ && chunk_id < local_chunk_end_);
}
} // namespace kagen
