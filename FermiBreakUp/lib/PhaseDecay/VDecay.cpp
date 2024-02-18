//
// Created by Artem Novikov on 21.02.2023.
//

#include "VDecay.h"
#include "Utilities/Randomizer.h"

VDecay::~VDecay() = default;

FermiFloat VDecay::TwoBodyMomentum(FermiFloat total_energy, FermiFloat mass1, FermiFloat mass2) {
  FermiFloat sq_probability =
      (total_energy + mass1 + mass2) * (total_energy + mass1 - mass2) * (total_energy - mass1 + mass2)
          * (total_energy - mass1 - mass2) / (4.0 * std::pow(total_energy, 2));

  if (sq_probability <= 0) { return -1; }

  return std::sqrt(sq_probability);
}

std::pair<LorentzVector, LorentzVector> VDecay::TwoBodyDecay(FermiFloat total_energy,
                                                             FermiFloat mass1,
                                                             FermiFloat mass2) {
  auto psqr = TwoBodyMomentum(total_energy, mass1, mass2);
  ParticleMomentum momentum = Randomizer::IsotropicVector(std::sqrt(psqr));

  std::pair<LorentzVector, LorentzVector> momentum4 = {LorentzVector(momentum, std::sqrt(psqr + std::pow(mass1, 2))),
                                                       LorentzVector(-momentum, std::sqrt(psqr + std::pow(mass2, 2)))};

  return momentum4;
}