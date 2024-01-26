//
// Created by Artem Novikov on 17.02.2023.
//

#include "FermiFragment.h"
#include "Utilities/NucleiProperties/NucleiProperties.h"

FermiFragment::FermiFragment(MassNumber mass_number, ChargeNumber charge_number,
                             int polarization, FermiFloat excitation_energy)
    : mass_number_(mass_number), charge_number_(charge_number),
      polarization(polarization), excitation_energy_(excitation_energy) {}

MassNumber FermiFragment::GetA() const {
  return GetMassNumber();
}

MassNumber FermiFragment::GetMassNumber() const {
  return mass_number_;
}

ChargeNumber FermiFragment::GetZ() const {
  return GetChargeNumber();
}

ChargeNumber FermiFragment::GetChargeNumber() const {
  return charge_number_;
}

int32_t FermiFragment::GetPolarization() const {
  return polarization;
}

FermiFloat FermiFragment::GetExcitationEnergy() const {
  return excitation_energy_;
}

FermiFloat FermiFragment::GetFragmentMass() const {
  return properties::NucleiProperties().GetNuclearMass(mass_number_, charge_number_);
}

FermiFloat FermiFragment::GetTotalEnergy() const {
  return GetFragmentMass() + GetExcitationEnergy();
}

FermiFragment::~FermiFragment() = default;
