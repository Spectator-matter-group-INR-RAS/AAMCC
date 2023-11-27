//
// Created by Artem Novikov on 21.02.2023.
//

#include <numeric>
#include <algorithm>
#include <functional>

#include "Decay.h"
#include "Utilities/Randomizer.h"

static const size_t MaxTryCount = 1e6;

std::vector<LorentzVector> Decay::CalculateDecay(const LorentzVector& momentum,
                                                 const std::vector<FermiFloat>& fragments_mass) const {
  std::vector<LorentzVector> result;
  result.reserve(fragments_mass.size());
  auto total_energy = momentum.m();

  /// 2 bodies case is faster
  if (fragments_mass.size() == 2) {
    auto momentum4 = TwoBodyDecay(total_energy, fragments_mass[0], fragments_mass[1]);
    result.push_back(momentum4.first);
    result.push_back(momentum4.second);
    return result;
  }

  /// N body case
  FermiFloat total_mass = std::accumulate(fragments_mass.begin(), fragments_mass.end(), FermiFloat(0));
  FermiFloat max_energy = total_energy - total_mass + fragments_mass[0];
  FermiFloat max_weight = CalculateMaxWeight(fragments_mass, max_energy);

  size_t try_count = 0;
  FermiFloat configuration_weight = 0;
  std::vector<FermiFloat> momentum_magnitudes(fragments_mass.size());
  while (configuration_weight < Randomizer::UniformRealDistribution() * max_weight) {
    auto virtual_masses = CalculateVirtualMasses(fragments_mass, total_energy - total_mass);

    configuration_weight = CalculateMomentumMagnitudes(momentum_magnitudes, fragments_mass, virtual_masses);

    if (try_count++ > MaxTryCount) {
      throw std::runtime_error("FermiPhaseSpaceDecay::Decay: Cannot determine decay kinematics");
    }
  }

  ParticleMomentum pivot_momentum = Randomizer::IsotropicVector(momentum_magnitudes[0]);
  result.emplace_back(pivot_momentum, std::sqrt(pivot_momentum.mag2() + std::pow(fragments_mass[0], 2)));  // 0
  result.emplace_back(-pivot_momentum, std::sqrt(pivot_momentum.mag2() + std::pow(fragments_mass[1], 2))); // 1

  for (size_t i = 2; i < fragments_mass.size(); ++i) {
    pivot_momentum = Randomizer::IsotropicVector(momentum_magnitudes[i - 1]);
    result.emplace_back(pivot_momentum, std::sqrt(pivot_momentum.mag2() + std::pow(fragments_mass[i], 2)));
    Vector3 boost_vector = (-1.0) * result.back().boostVector();
    /// boost already created particles
    for (auto& fragment_momentum: result) {
      fragment_momentum.boost(boost_vector);
    }
  }

  return result;
}

std::vector<FermiFloat> Decay::CalculateVirtualMasses(const std::vector<FermiFloat>& masses, FermiFloat energy) {
  auto probability_distribution = Randomizer::ProbabilityDistribution(masses.size());

  /// Calculate virtual masses
  std::vector<FermiFloat> virtual_masses(masses.size());
  virtual_masses[0] = 0;
  std::partial_sum(masses.begin(), masses.end(), virtual_masses.begin());

  std::transform(probability_distribution.begin(),
                 probability_distribution.end(),
                 probability_distribution.begin(),
                 std::bind(std::multiplies<FermiFloat>(), std::placeholders::_1, energy));

  std::transform(probability_distribution.begin(),
                 probability_distribution.end(),
                 virtual_masses.begin(),
                 virtual_masses.begin(),
                 std::plus<FermiFloat>());
  return virtual_masses;
}

FermiFloat Decay::CalculateMomentumMagnitudes(std::vector<FermiFloat>& daughter_momentum,
                                              const std::vector<FermiFloat>& masses,
                                              const std::vector<FermiFloat>& virtual_masses) {
  size_t sz = daughter_momentum.size();
  FermiFloat weight = 1;
  for (size_t i = 0; i < sz - 1; ++i) {
    daughter_momentum[i] = TwoBodyMomentum(virtual_masses[i + 1], virtual_masses[i], masses[i + 1]);
    weight *= daughter_momentum[i];

    if (daughter_momentum[i] < 0) {
      std::cerr << "FermiPhaseSpaceDecay::Decay: Daughter momentum less than zero\n";
      return 0;
    }
  }
  daughter_momentum[sz - 1] = TwoBodyMomentum(virtual_masses[sz - 2], masses[sz - 2], masses[sz - 1]);
  return weight;
}

FermiFloat Decay::CalculateMaxWeight(const std::vector<FermiFloat>& masses, FermiFloat max_energy) {
  FermiFloat min_energy = 0;
  FermiFloat max_weight = 1.0;
  for (size_t i = 0; i < masses.size(); ++i) {
    max_energy += masses[i];
    min_energy += masses[i - 1];
    max_weight *= TwoBodyMomentum(max_energy, min_energy, masses[i]);
  }
  return max_weight;
}



