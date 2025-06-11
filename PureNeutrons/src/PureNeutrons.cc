//
// Created by Artem Novikov on 25.11.2023.
//

#include "G4FermiPhaseSpaceDecay.hh"

#include "../include/PureNeutrons.hh"

PureNeutrons::PureNeutrons(G4double neutron_mass) : neutron_mass_(neutron_mass) {}

G4FragmentVector PureNeutrons::BreakItUp(const G4Fragment& fragment) const {
  G4FermiPhaseSpaceDecay phase_decay;

  auto fragments = G4FragmentVector();
  fragments.reserve(fragment.GetA_asInt());

  auto boost_vector = fragment.GetMomentum().boostVector();
  auto fragments_momentum = std::move(*phase_decay.Decay(fragment.GetMomentum().m(),
                                                         std::vector<G4double>(fragment.GetA_asInt(), neutron_mass_)));
  for (auto momentum_ptr: fragments_momentum) {
    fragments.push_back(new G4Fragment(1, 0, momentum_ptr->boost(boost_vector)));
  }

  return fragments;
}
