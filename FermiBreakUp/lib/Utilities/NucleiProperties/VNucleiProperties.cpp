//
// Created by Artem Novikov on 09.03.2023.
//

#include <CLHEP/Units/PhysicalConstants.h>

#include "VNucleiProperties.h"

namespace properties {

VNucleiProperties::~VNucleiProperties() = default;

FermiFloat VNucleiProperties::AtomicMass(MassNumber mass_number, ChargeNumber charge_number) {
  static const FermiFloat hydrogen_mass_excess = 7.28897;
  static const FermiFloat neutron_mass_excess = 8.07132;

  FermiFloat mass = FermiFloat(FermiUInt(mass_number) - FermiUInt(charge_number)) * neutron_mass_excess
      + FermiFloat(charge_number) * hydrogen_mass_excess - BindingEnergy(mass_number, charge_number)
      + FermiFloat(mass_number) * CLHEP::amu_c2;

  return mass;
}

FermiFloat VNucleiProperties::NuclearMass(MassNumber mass_number, ChargeNumber charge_number) {
  FermiFloat mass = AtomicMass(mass_number, charge_number);
  /// atomic mass is converted to nuclear mass according formula in AME03
  mass -= FermiFloat(charge_number) * CLHEP::electron_mass_c2;
  mass += (14.4381 * std::pow(charge_number, 2.39) + 1.55468 * 1e-6 * std::pow(charge_number, 5.35)) * CLHEP::eV;

  return mass;
}

FermiFloat VNucleiProperties::BindingEnergy(MassNumber mass_number, ChargeNumber charge_number) {
  /// Weitzsaecker's Mass formula
  FermiInt N_pairing = (FermiInt(mass_number) - FermiInt(charge_number)) % 2;                      /// pairing
  FermiInt Z_pairing = FermiInt(charge_number) % 2;

  FermiFloat binding =
      -15.67 * FermiFloat(mass_number)                                                           /// nuclear volume
          + 17.23 * std::pow(mass_number, 2. / 3.)                                                    /// surface energy
          + 93.15 * std::pow(FermiFloat(mass_number) / 2. - FermiFloat(charge_number), 2)
              / FermiFloat(mass_number)                                                                   /// asymmetry
          + 0.6984523 * std::pow(FermiUInt(charge_number), 2) * std::pow(mass_number, -1. / 3.);      /// coulomb

  if (N_pairing == Z_pairing) {
    binding += (N_pairing + Z_pairing - 1) * 12.0 / std::sqrt(FermiFloat(mass_number));               /// pairing
  }

  return -binding * CLHEP::MeV;
}

bool VNucleiProperties::IsInvalidNuclei(MassNumber mass_number, ChargeNumber charge_number) {
  return mass_number < 1_m || charge_number < 0_c || FermiUInt(charge_number) > FermiUInt(mass_number);
}

void VNucleiProperties::PrintInvalidNuclei(MassNumber mass_number, ChargeNumber charge_number) {
  std::cerr << "Unsupported values for A = " << mass_number << " and Z = " << charge_number << std::endl;
}

} // namespace properties
