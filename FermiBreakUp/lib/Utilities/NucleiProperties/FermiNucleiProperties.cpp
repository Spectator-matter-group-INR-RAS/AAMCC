//
// Created by Artem Novikov on 09.03.2023.
//

#include "FermiNucleiProperties.h"
#include "Builder/DefaultBuilder.h"

namespace properties {

std::unique_ptr<FermiNucleiProperties::MassMap> FermiNucleiProperties::nuclei_mass_ = nullptr;

FermiNucleiProperties::FermiNucleiProperties() {
  if (nuclei_mass_ == nullptr) {
    Build(DefaultBuilder());
  }
}

FermiNucleiProperties::FermiNucleiProperties(const VFermiPropertiesBuilder& builder) {
  Build(builder);
}

FermiFloat FermiNucleiProperties::GetNuclearMass(MassNumber mass_number, ChargeNumber charge_number) const {
  auto it = nuclei_mass_->find(NucleiData{mass_number, charge_number});
  if (it != nuclei_mass_->end()) {
    return it->second;
  }

//  throw std::runtime_error("Unknown particle");
  std::cerr << "Unknown particle A = " + std::to_string(mass_number) + ", Z = " + std::to_string(charge_number) << '\n';
  return NuclearMass(mass_number, charge_number);
}

bool FermiNucleiProperties::IsStable(MassNumber mass_number, ChargeNumber charge_number) const {
  if (IsInvalidNuclei(mass_number, charge_number)) {
    std::cerr << "Unsupported values for A = " << mass_number << " and Z = " << charge_number << std::endl;
    return false;
  }

  return nuclei_mass_->find(NucleiData{mass_number, charge_number}) != nuclei_mass_->end();
}

void FermiNucleiProperties::Build(const VFermiPropertiesBuilder& builder) {
  nuclei_mass_ = std::make_unique<MassMap>();
  builder.BuildTable(*nuclei_mass_);
}

} // namespace properties
