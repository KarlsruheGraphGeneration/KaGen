#include "kagen/generators/brain/brain.h"

#include "kagen/context.h"
#include "kagen/kagen.h"

#include <algorithm/AlgorithmEnum.h>
#include <mpi.h>

#include "algorithm/BarnesHutInternal/BarnesHutCell.h"
#include "algorithm/BarnesHutInternal/BarnesHutInvertedCell.h"
#include "algorithm/Kernel/Kernel.h"
#include "algorithm/NaiveInternal/NaiveCell.h"
#include "neurons/CalciumCalculator.h"
#include "neurons/helper/SynapseDeletionFinder.h"
#include "neurons/input/BackgroundActivityCalculators.h"
#include "neurons/input/FiredStatusCommunicationMap.h"
#include "neurons/input/SynapticInputCalculator.h"
#include "neurons/input/SynapticInputCalculators.h"
#include "neurons/models/NeuronModels.h"
#include "neurons/models/SynapticElements.h"
#include "sim/Essentials.h"
#include "sim/Simulation.h"
#include "sim/random/SubdomainFromNeuronPerRank.h"
#include "structure/Partition.h"
#include <memory>
namespace kagen {
PGeneratorConfig
BrainFactory::NormalizeParameters(PGeneratorConfig config, PEID, const PEID size, const bool /*output*/) const {
    config.brain.num_neurons_per_rank = (config.n + size - 1) / size;
    return config;
}
std::unique_ptr<Generator>
BrainFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<BrainGenerator>(config, rank, size);
}

BrainGenerator::BrainGenerator(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : _sim(std::make_unique<Essentials>(), std::shared_ptr<Partition>(nullptr)),
      _config(config.brain),
      _rank(rank) {
    MPIWrapper::init(1, nullptr);
    AlgorithmEnum algorithm;
    switch (_config.algorithm) {
        case BrainSynapseCreationModel::NAIVE:
            algorithm = AlgorithmEnum::Naive;
            break;
        case BrainSynapseCreationModel::BARNES_HUT:
            algorithm = AlgorithmEnum::BarnesHut;
            break;
        case BrainSynapseCreationModel::BARNES_HUT_INVERTED:
            algorithm = AlgorithmEnum::BarnesHutInverted;
            break;
    }
    SetVertexRange(rank * _config.num_neurons_per_rank, (rank + 1) * _config.num_neurons_per_rank);
    // AlgorithmEnum chosen_algorithm = AlgorithmEnum::BarnesHut;

    // RelearnTypes::step_type simulation_steps{1001};

    // unsigned int random_seed{0};
    // app.add_option("-r,--random-seed", random_seed, "Random seed. Default: 0.");

    // RelearnTypes::number_neurons_type number_neurons_per_rank{1000};

    double fraction_excitatory_neurons{1.0};

    double um_per_neuron{1.0};

    // double gaussian_sigma{ GaussianDistributionKernel::default_sigma };
    // double gaussian_sigma{2};

    // double gaussian_mu{ GaussianDistributionKernel::default_mu };
    // double gaussian_mu{ 0 };

    // double synapses = 5.0;
    double synaptic_elements_init_lb = _config.synapses_lb;

    double synaptic_elements_init_ub = _config.synapses_ub;

    RelearnException::check(
        synaptic_elements_init_lb >= 0.0, "The minimum number of vacant synaptic elements must not be negative");
    RelearnException::check(
        synaptic_elements_init_ub >= synaptic_elements_init_lb,
        "The minimum number of vacant synaptic elements must not be larger than the maximum number");

    omp_set_num_threads(1);

    std::size_t current_seed = 0;
    boost::hash_combine(current_seed, rank);
    boost::hash_combine(current_seed, config.seed);

    RandomHolder::seed_all(current_seed);
    LogFiles::disable = true;

    auto essentials = std::make_unique<Essentials>();

    MPIWrapper::barrier();

    // Timers::start(TimerRegion::INITIALIZATION);

    auto prepare_algorithm = [&]() -> void {
        // Set the correct kernel and initialize the MPIWrapper to return the correct type
        if (algorithm == AlgorithmEnum::BarnesHut) {
            MPIWrapper::init_buffer_octree<BarnesHutCell>();
        } else if (algorithm == AlgorithmEnum::BarnesHutInverted) {
            MPIWrapper::init_buffer_octree<BarnesHutInvertedCell>();
        } else {
            RelearnException::check(
                algorithm == AlgorithmEnum::Naive, "An algorithm was chosen that is not supported");
            MPIWrapper::init_buffer_octree<NaiveCell>();
        }

        GaussianDistributionKernel::set_sigma(_config.gaussian_sigma);
        GaussianDistributionKernel::set_mu(_config.gaussian_mu);
    };
    prepare_algorithm();

    auto partition = std::make_shared<Partition>(size, MPIRank{rank});

    auto subdomain = std::make_unique<SubdomainFromNeuronPerRank>(
        _config.num_neurons_per_rank, fraction_excitatory_neurons, um_per_neuron, partition);

    auto background_activity_calculator = std::make_unique<NullBackgroundActivityCalculator>();

    auto stimulus_calculator = std::make_unique<Stimulus>();

    auto construct_input = [&]() -> std::unique_ptr<SynapticInputCalculator> {
        auto fired_status_communicator = std::make_unique<FiredStatusCommunicationMap>(MPIWrapper::get_num_ranks());
        return std::make_unique<LinearSynapticInputCalculator>(
            SynapticInputCalculator::default_conductance, std::move(fired_status_communicator));
    };
    auto input_calculator = construct_input();

    auto neuron_model = std::make_unique<models::PoissonModel>(
        NeuronModel::default_h, std::move(input_calculator), std::move(background_activity_calculator),
        std::move(stimulus_calculator), models::PoissonModel::default_x_0, models::PoissonModel::default_tau_x,
        models::PoissonModel::default_refractory_period);

    auto construct_calcium_calculator = [&]() -> std::unique_ptr<CalciumCalculator> {
        auto calcium_calculator = std::make_unique<CalciumCalculator>();

        auto initial_calcium_calculator = [](MPIRank /*mpi_rank*/, NeuronID::value_type /*neuron_id*/) {
            return 0.0;
        };
        calcium_calculator->set_initial_calcium_calculator(std::move(initial_calcium_calculator));

        auto target_calcium_calculator = [](MPIRank /*mpi_rank*/, NeuronID::value_type /*neuron_id*/) {
            return CalciumCalculator::default_C_target;
        };
        calcium_calculator->set_target_calcium_calculator(std::move(target_calcium_calculator));

        return calcium_calculator;
    };
    auto calcium_calculator = construct_calcium_calculator();

    auto axons_model = std::make_shared<SynapticElements>(
        ElementType::Axon, SynapticElements::default_eta_Axons, SynapticElements::default_nu,
        SynapticElements::default_vacant_retract_ratio, synaptic_elements_init_lb, synaptic_elements_init_ub);

    auto excitatory_dendrites_model = std::make_shared<SynapticElements>(
        ElementType::Dendrite, SynapticElements::default_eta_Dendrites_exc, SynapticElements::default_nu,
        SynapticElements::default_vacant_retract_ratio, synaptic_elements_init_lb, synaptic_elements_init_ub);

    auto inhibitory_dendrites_model = std::make_shared<SynapticElements>(
        ElementType::Dendrite, SynapticElements::default_eta_Dendrites_inh, SynapticElements::default_nu,
        SynapticElements::default_vacant_retract_ratio, synaptic_elements_init_lb, synaptic_elements_init_ub);

    auto synapse_deletion_finder = std::make_unique<RandomSynapseDeletionFinder>();

    _sim = Simulation(std::move(essentials), partition);
    _sim.set_neuron_model(std::move(neuron_model));
    _sim.set_calcium_calculator(std::move(calcium_calculator));
    _sim.set_synapse_deletion_finder(std::move(synapse_deletion_finder));
    _sim.set_axons(std::move(axons_model));
    _sim.set_dendrites_ex(std::move(excitatory_dendrites_model));
    _sim.set_dendrites_in(std::move(inhibitory_dendrites_model));

    if (is_barnes_hut(algorithm)) {
        _sim.set_acceptance_criterion_for_barnes_hut(Constants::bh_default_theta);
    }

    _sim.set_algorithm(algorithm);
    _sim.set_subdomain_assignment(std::move(subdomain));

    /**********************************************************************************/

    // The barrier ensures that every rank finished its local stores.
    // Otherwise, a "fast" rank might try to read from the RMA window of another
    // rank which has not finished (or even begun) its local stores
    MPIWrapper::barrier(); // TODO(future) Really needed?

    // Lock local RMA memory for local stores
    MPIWrapper::lock_window(MPIWindow::Window::Octree, MPIRank{rank}, MPI_Locktype::Exclusive);

    _sim.initialize();

    // Unlock local RMA memory and make local stores visible in public window copy
    MPIWrapper::unlock_window(MPIWindow::Window::Octree, MPIRank{rank});
}

void BrainGenerator::GenerateEdgeList() {
    _sim.simulate(_config.simulation_steps);

    MPIWrapper::barrier();

    _sim.finalize();
    auto graph        = _sim.get_network_graph();
    auto to_global_id = [&](auto neuron_id) {
        PEID                 rank;
        NeuronID::value_type local_id;
        if constexpr (std::is_same_v<std::remove_cv_t<decltype(neuron_id)>, NeuronID>) {
            rank     = this->_rank;
            local_id = neuron_id.get_neuron_id();
        } else if constexpr (std::is_same_v<std::remove_cv_t<decltype(neuron_id)>, RankNeuronId>) {
            rank     = neuron_id.get_rank().get_rank();
            local_id = neuron_id.get_neuron_id().get_neuron_id();
        } else {
            static_assert(std::is_same_v<std::remove_cv_t<decltype(neuron_id)>, RankNeuronId>, "Invalid Neuron type");
        }

        auto local_offset = rank * _config.num_neurons_per_rank;
        return local_offset + local_id;
    };

    for (auto const& local_neuron_id: NeuronID::range(graph->get_number_neurons())) {
        std::unordered_set<SInt> neighbors;
        auto                     push_edges = [&](auto const& edge_list) {
            auto const [plastic_edges, static_edges] = edge_list;
            for (auto const& [head, _weight]: plastic_edges) {
                auto head_global = to_global_id(head);
                if (neighbors.insert(head_global).second) {
                    PushEdge(to_global_id(local_neuron_id), head_global);
                }
            }
            for (auto const& [head, _weight]: static_edges) {
                auto head_global = to_global_id(head);
                if (neighbors.insert(head_global).second) {
                    PushEdge(to_global_id(local_neuron_id), head_global);
                }
            }
        };
        push_edges(graph->get_local_out_edges(local_neuron_id));
        push_edges(graph->get_distant_out_edges(local_neuron_id));
        push_edges(graph->get_local_in_edges(local_neuron_id));
        push_edges(graph->get_distant_in_edges(local_neuron_id));
    }
    MPIWrapper::finalize(false);
}

} // namespace kagen
