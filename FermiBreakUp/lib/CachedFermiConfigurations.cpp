//
// Created by Artem Novikov on 24.05.2023.
//

#include "CachedFermiConfigurations.h"
#include "Utilities/ConfigurationProperties.h"
#include "Utilities/Randomizer.h"

CachedFermiConfigurations::CachedFermiConfigurations(NucleiData nuclei_data, FermiFloat total_energy) {
  CachedFermiConfigurations::GenerateSplits(nuclei_data, total_energy);
}

VFermiConfigurations& CachedFermiConfigurations::GenerateSplits(NucleiData nuclei_data, FermiFloat total_energy) {
  auto max_fragments_count = FermiUInt(nuclei_data.mass_number);
  if (nuclei_data != last_nuclei_) {
    last_nuclei_ = nuclei_data;
    cached_configurations_.clear();

    for (uint32_t particle_count = 2; particle_count <= max_fragments_count; particle_count++) {
      for (auto& split : FermiSplit(nuclei_data, particle_count)) {
        cached_configurations_.emplace_back(std::move(split));  /// split is moved!
      }
    }
  }

  configurations_.clear();
  weights_.clear();

  FermiFloat total_weight = 0;
  for(size_t i = 0; i < cached_configurations_.size(); ++i) {
    auto split_weight = ConfigurationProperties::DecayProbability(cached_configurations_[i], nuclei_data.mass_number, total_energy);

    if (split_weight != 0) {
      total_weight += split_weight;

      weights_.push_back(split_weight);
      configurations_.emplace_back(i);
    }
  }

  std::transform(weights_.begin(), weights_.end(),
                 weights_.begin(), std::bind(std::divides<FermiFloat>(), std::placeholders::_1, total_weight));

  return *this;
}

std::optional<FragmentVector> CachedFermiConfigurations::ChooseSplit() {
  if (configurations_.empty()) {
    return {};
  }

  FermiFloat wheel_result = Randomizer::UniformRealDistribution();
  FermiFloat accumulated_weight = 0;

  for (size_t i = 0; i < weights_.size(); ++i) {
    accumulated_weight += weights_[i];

    if (accumulated_weight >= wheel_result) {
      return cached_configurations_[configurations_[i]];
    }
  }

  throw std::runtime_error("No split chosen, something went wrong!");
}
