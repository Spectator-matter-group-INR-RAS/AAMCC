//
// Created by Artem Novikov on 21.02.2023.
//

#include <numeric>
#include <functional>
#include <CLHEP/Random/RandGamma.h>

#include "KopylovDecay.h"
#include "Utilities/Randomizer.h"

std::vector<LorentzVector> KopylovDecay::CalculateDecay(const LorentzVector& momentum,
                                                        const std::vector<FermiFloat>& fragments_mass) const {
  std::vector<LorentzVector> result(fragments_mass.size());

  FermiFloat parent_mass = momentum.m();
  FermiFloat total_fragments_mass = std::accumulate(fragments_mass.begin(), fragments_mass.end(), FermiFloat(0));
  FermiFloat mu = total_fragments_mass;
  FermiFloat mass = parent_mass;
  FermiFloat kinetic_energy = parent_mass - total_fragments_mass;

  auto momentum_rest_lab = LorentzVector(0, 0, 0, parent_mass);
  for (size_t i = fragments_mass.size() - 1; i > 0; --i) {
    mu -= fragments_mass[i];
    kinetic_energy *= i > 1 ? BetaKopylov(i) : 0;
    auto rest_mass = mu + kinetic_energy;

    auto momentum_magnitude_fragments_cm = TwoBodyMomentum(mass, fragments_mass[i], rest_mass);
    if (momentum_magnitude_fragments_cm < 0) {
      throw std::runtime_error("FermiPhaseSpaceDecay::KopylovNBodyDecay: Error sampling fragments momenta!!");
    }

    ParticleMomentum random_momentum = Randomizer::IsotropicVector(momentum_magnitude_fragments_cm);
    auto momentum_fragments_cm = LorentzVector(random_momentum,
                                               std::sqrt(random_momentum.mag2() + std::pow(fragments_mass[i], 2)));
    auto momentum_rest_cm = LorentzVector(-random_momentum, std::sqrt(random_momentum.mag2() + std::pow(rest_mass, 2)));

    /// change framework
    Vector3 boost_vector = momentum_rest_lab.boostVector();

    momentum_rest_lab = LorentzVector(momentum_rest_cm).boost(boost_vector);
    auto momentum_fragments_lab = LorentzVector(momentum_fragments_cm).boost(boost_vector);

    result[i] = momentum_fragments_lab;

    mass = rest_mass;
  }

  result[0] = momentum_rest_lab;

  return result;
}

FermiFloat KopylovDecay::BetaKopylov(size_t k) {
  /// Notice that alpha > beta always
  static const FermiFloat beta = 1.5;
  FermiFloat alpha = 1.5 * static_cast<FermiFloat>(k - 1);
  FermiFloat y1 = CLHEP::RandGamma::shoot(alpha, 1);
  FermiFloat y2 = CLHEP::RandGamma::shoot(beta, 1);

  return y1 / (y1 + y2);
}
