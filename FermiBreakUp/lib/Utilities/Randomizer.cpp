//
// Created by Artem Novikov on 19.02.2023.
//

#include <random>
#include <CLHEP/Units/PhysicalConstants.h>

#include "Randomizer.h"

FermiFloat Randomizer::UniformRealDistribution() {
  static std::random_device device;
  static std::mt19937 generator(device());
  static std::uniform_real_distribution<FermiFloat> distribution(0, 1);
  return distribution(generator);
}

FermiFloat Randomizer::NormalDistribution(FermiFloat mean, FermiFloat deviation) {
  static std::random_device device;
  static std::mt19937 generator(device());
  static std::normal_distribution<FermiFloat> distribution(0, 1);
  return mean + distribution(generator) * deviation;
}

ParticleMomentum Randomizer::IsotropicVector(FermiFloat magnitude) {
  auto cos = 1.0 - 2.0 * UniformRealDistribution();
  auto sin = std::sqrt(1.0 - std::pow(cos, 2));
  auto phi = CLHEP::twopi * UniformRealDistribution();
  ParticleMomentum momentum(magnitude * std::cos(phi) * sin,
                            magnitude * std::sin(phi) * sin,
                            magnitude * cos);
//  ParticleMomentum momentum(NormalDistribution(), NormalDistribution(), NormalDistribution());
//  momentum = momentum * (magnitude) / momentum.mag();
  return momentum;
}

std::vector<FermiFloat> Randomizer::ProbabilityDistribution(size_t point_count) {
  /// Sample uniform random numbers in increasing order
  std::vector<FermiFloat> probability_distribution;
  probability_distribution.reserve(point_count);

  probability_distribution.push_back(0);
  std::generate_n(std::back_inserter(probability_distribution), point_count - 2, Randomizer::UniformRealDistribution);
  probability_distribution.push_back(1);

  std::sort(probability_distribution.begin(), probability_distribution.end());

  return probability_distribution;
}
