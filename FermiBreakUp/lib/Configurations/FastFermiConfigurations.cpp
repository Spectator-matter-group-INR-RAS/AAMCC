//
// Created by Artem Novikov on 18.02.2024.
//

#include "Utilities/ConfigurationProperties.h"
#include "Utilities/Randomizer.h"

#include "FastFermiConfigurations.h"

FastFermiConfigurations::FastFermiConfigurations(NucleiData nuclei_data, FermiFloat total_energy) {
  FastFermiConfigurations::GenerateSplits(nuclei_data, total_energy);
}

VFermiConfigurations& FastFermiConfigurations::GenerateSplits(NucleiData nuclei_data, FermiFloat total_energy) {
  last_nuclei_= nuclei_data;
  auto max_fragments_count = FermiUInt(nuclei_data.mass_number);
  auto& cache = cached_configurations_[nuclei_data];
  if (cache.empty()) {
    cache.reserve(100);

    for (uint32_t particle_count = 2; particle_count <= max_fragments_count; ++particle_count) {
      for (auto& split : FermiSplit(nuclei_data, particle_count)) {
        cache.emplace_back(std::move(split));  /// split is moved!
      }
    }
  }

  configurations_.clear();
  weights_.clear();

  FermiFloat total_weight = 0;
  for (size_t i = 0; i < cache.size(); ++i) {
    auto split_weight =
        ConfigurationProperties::DecayProbability(cache[i], nuclei_data.mass_number, total_energy);

    if (split_weight != 0) {
      total_weight += split_weight;

      weights_.push_back(split_weight);
      configurations_.emplace_back(i);
    }
  }

  std::transform(weights_.begin(), weights_.end(),
                 weights_.begin(), std::bind(std::divides<>(), std::placeholders::_1, total_weight));

  return *this;
}

std::optional<FragmentVector> FastFermiConfigurations::ChooseSplit() {
  if (configurations_.empty()) {
    return {};
  }

  FermiFloat wheel_result = Randomizer::UniformRealDistribution();
  FermiFloat accumulated_weight = 0;

  for (size_t i = 0; i < weights_.size(); ++i) {
    accumulated_weight += weights_[i];

    if (accumulated_weight >= wheel_result) {
      return cached_configurations_[last_nuclei_][configurations_[i]];
    }
  }

  throw std::runtime_error("No split chosen, something went wrong!");
}

std::unique_ptr<VFermiConfigurations> FastFermiConfigurations::Clone() const {
  return std::make_unique<FastFermiConfigurations>(*this);
}
