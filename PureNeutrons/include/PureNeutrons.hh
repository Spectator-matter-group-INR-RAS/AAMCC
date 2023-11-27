//
// Created by Artem Novikov on 25.11.2023.
//

#ifndef GRATE_PURENEUTRONS_SRC_PURENEUTRONS_H_
#define GRATE_PURENEUTRONS_SRC_PURENEUTRONS_H_

#include "G4Fragment.hh"

class PureNeutrons {
 public:
  PureNeutrons(G4double neutron_mass = CLHEP::neutron_mass_c2);

  PureNeutrons(PureNeutrons&&) = default;

  PureNeutrons(const PureNeutrons&) = delete;

  PureNeutrons& operator=(const PureNeutrons&) = delete;

  PureNeutrons& operator=(PureNeutrons&&) = default;

  G4FragmentVector BreakItUp(const G4Fragment& fragment) const;

 private:
  G4double neutron_mass_;
};

#endif //GRATE_PURENEUTRONS_SRC_PURENEUTRONS_H_
