//
// Created by Artem Novikov on 19.02.2023.
//

#include <functional>

#include "FermiConfigurations.h"
#include "Utilities/ConfigurationProperties.h"
#include "Utilities/Randomizer.h"

FermiConfigurations::FermiConfigurations(NucleiData nuclei_data, FermiFloat total_energy) {
  FermiConfigurations::GenerateSplits(nuclei_data, total_energy);
}

VFermiConfigurations& FermiConfigurations::GenerateSplits(NucleiData nuclei_data, FermiFloat total_energy) {
  configurations_.clear();
  weights_.clear();

  /// let's split nucleus into 2,...,A fragments
  FermiFloat total_weight = 0;
  auto max_fragments_count = FermiUInt(nuclei_data.mass_number);

  for (uint32_t particle_count = 2; particle_count <= max_fragments_count; particle_count++) {
    /// initialize configuration for k fragments
    for (auto& split : FermiSplit(nuclei_data, particle_count)) {
      /// non-normalized statistical weight for given channel with k fragments
      auto split_weight = ConfigurationProperties::DecayProbability(split, nuclei_data.mass_number, total_energy);
      if (split_weight != 0) {
        total_weight += split_weight;

        weights_.push_back(split_weight);
        configurations_.emplace_back(std::move(split));  /// split is moved!
      }
    }
  }

  /// let's normalize statistical weights of channels
  std::transform(weights_.begin(), weights_.end(),
                 weights_.begin(), std::bind(std::divides(), std::placeholders::_1, total_weight));

  return *this;
}

std::optional<FragmentVector> FermiConfigurations::ChooseSplit() {
  if (configurations_.empty()) {
    return {};
  }

  FermiFloat wheel_result = Randomizer::UniformRealDistribution();
  FermiFloat accumulated_weight = 0;

  for (size_t i = 0; i < weights_.size(); ++i) {
    accumulated_weight += weights_[i];

    if (accumulated_weight >= wheel_result) {
      return configurations_[i];
    }
  }

  throw std::runtime_error("No split chosen, something went wrong!");
}
