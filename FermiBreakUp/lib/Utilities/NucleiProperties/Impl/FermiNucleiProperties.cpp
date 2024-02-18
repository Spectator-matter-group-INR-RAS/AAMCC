//
// Created by Artem Novikov on 09.03.2023.
//

#include "FermiNucleiProperties.h"
#include "../DataStorage/DefaultNuclearMass.h"

namespace properties {

FermiNucleiProperties::FermiNucleiProperties() : FermiNucleiProperties(DefaultNuclearMass()) {}

FermiFloat FermiNucleiProperties::GetNuclearMass(MassNumber mass_number, ChargeNumber charge_number) const {
  auto it = nuclei_mass_.find(NucleiData{mass_number, charge_number});
  if (it != nuclei_mass_.end()) {
    return it->second;
  }

  #ifdef _DEBUG
  std::cerr << "Unknown particle A = " + std::to_string(mass_number) + ", Z = " + std::to_string(charge_number) << '\n';
  #endif
  return NuclearMass(mass_number, charge_number);
}

bool FermiNucleiProperties::IsStable(MassNumber mass_number, ChargeNumber charge_number) const {
  if (IsInvalidNuclei(mass_number, charge_number)) {
    PrintInvalidNuclei(mass_number, charge_number);
    return false;
  }

  return nuclei_mass_.contains(NucleiData{mass_number, charge_number});
}

void FermiNucleiProperties::AddMass(MassNumber mass_number, ChargeNumber charge_number, FermiFloat mass) {
  return AddMass(NucleiData{mass_number, charge_number}, mass);
}

void FermiNucleiProperties::AddMass(NucleiData nuclei_data, FermiFloat mass) {
  nuclei_mass_[nuclei_data] = mass;
}

} // namespace properties
