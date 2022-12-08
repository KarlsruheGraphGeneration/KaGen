#include "kagen/generators/image/image_mesh.h"

#include "kagen/definitions.h"
#include "kagen/generators/generator.h"

namespace kagen {
std::unique_ptr<Generator> ImageMeshFactory::Create(const PGeneratorConfig& config, PEID rank, PEID size) const {
    return std::make_unique<ImageMesh>(config, rank, size);
}

ImageMesh::ImageMesh(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size) {}

void ImageMesh::GenerateImpl() {}
} // namespace kagen

