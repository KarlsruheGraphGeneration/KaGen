#pragma once

#include "kagen/context.h"
#include "kagen/generators/generator.h"

namespace kagen {
class RMATFactory : public GeneratorFactory {
public:
    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const override;
};

class RMAT : public Generator {
public:
    RMAT(const PGeneratorConfig& config, PEID rank, PEID size);

protected:
    void GenerateImpl() override;

private:
    const PGeneratorConfig& config_;
    PEID                    rank_, size_;
};
} // namespace kagen
