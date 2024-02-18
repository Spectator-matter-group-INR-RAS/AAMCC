//
// Created by Artem Novikov on 09.03.2023.
//

#include <CLHEP/Units/PhysicalConstants.h>

#include "NativeNucleiProperties.h"
#include "TableValues/NucleiPropertiesTable.h"
#include "TableValues/NucleiPropertiesTableAME12.h"

namespace properties {

const FermiFloat NativeNucleiProperties::proton_mass = CLHEP::proton_mass_c2;

const FermiFloat NativeNucleiProperties::neutron_mass = CLHEP::neutron_mass_c2;

const FermiFloat NativeNucleiProperties::deuteron_mass = 1.875613 * CLHEP::GeV;

const FermiFloat NativeNucleiProperties::triton_mass = 2.808921 * CLHEP::GeV;

const FermiFloat NativeNucleiProperties::alpha_mass = 3.727379 * CLHEP::GeV;

const FermiFloat NativeNucleiProperties::He3_mass = 2.808391 * CLHEP::GeV;

FermiFloat NativeNucleiProperties::GetNuclearMass(MassNumber mass_number, ChargeNumber charge_number) const {
  if (IsInvalidNuclei(mass_number, charge_number)) {
    std::cerr << "NucleiProperties::GetNuclearMass: Wrong values for A = " << mass_number
              << " and Z = " << charge_number << std::endl;
    return 0;
  }

  if (charge_number <= 2_c) {
    if (charge_number == 1_c && mass_number == 1_m) {
      return proton_mass;
    }
    if (charge_number == 0_c && mass_number == 1_m) {
      return neutron_mass;
    }
    if (charge_number == 1_c && mass_number == 2_m) {
      return deuteron_mass;
    }
    if (charge_number == 1_c && mass_number == 3_m) {
      return triton_mass;
    }
    if (charge_number == 2_c && mass_number == 4_m) {
      return alpha_mass;
    }
    if (charge_number == 2_c && mass_number == 3_m) {
      return He3_mass;
    }
  }

  NucleiPropertiesTableAME12 table_AME12;
  if (table_AME12.ContainsParticle(mass_number, charge_number)) {
    /// AME 03 table
    return table_AME12.GetNuclearMass(mass_number, charge_number);
  }

  NucleiPropertiesTable table_theoretical;
  if (table_theoretical.ContainsParticle(mass_number, charge_number)) {
    /// Theoretical table
    return table_theoretical.GetNuclearMass(mass_number, charge_number);
  }

  /// Raw calculation
  return NuclearMass(mass_number, charge_number);
}

bool NativeNucleiProperties::IsStable(MassNumber mass_number, ChargeNumber charge_number) const {
  if (IsInvalidNuclei(mass_number, charge_number)) {
    std::cerr << "NucleiProperties::IsInStableTable: Wrong values for A = "
              << mass_number << " and Z = " << charge_number << std::endl;
    return false;
  }

  return NucleiPropertiesTable().ContainsParticle(mass_number, charge_number);
}

FermiFloat NativeNucleiProperties::AtomicMass(MassNumber mass_number, ChargeNumber charge_number) {
  NucleiPropertiesTableAME12 table_AME12;
  const FermiFloat hydrogen_mass_excess = table_AME12.GetMassExcess(1_m, 1_c);
  const FermiFloat neutron_mass_excess  = table_AME12.GetMassExcess(1_m, 0_c);

  FermiFloat mass = FermiFloat(FermiUInt(mass_number) - FermiUInt(charge_number)) * neutron_mass_excess
      + FermiFloat(charge_number) * hydrogen_mass_excess - BindingEnergy(mass_number, charge_number)
      + FermiFloat(mass_number) * CLHEP::amu_c2;

  return mass;
}

FermiFloat NativeNucleiProperties::NuclearMass(MassNumber mass_number, ChargeNumber charge_number) {
  if (IsInvalidNuclei(mass_number, charge_number)) {
    std::cerr << "NativeNucleiProperties::NuclearMass: Wrong values for A = "
              << mass_number << " and Z = " << charge_number << std::endl;
    return 0;
  }

  FermiFloat mass = AtomicMass(mass_number, charge_number);
  /// atomic mass is converted to nuclear mass according formula in AME03
  mass -= FermiFloat(charge_number) * CLHEP::electron_mass_c2;
  mass += (14.4381 * std::pow(charge_number, 2.39) + 1.55468 * 1e-6 * std::pow(charge_number, 5.35)) * CLHEP::eV;

  return mass;
}

FermiFloat NativeNucleiProperties::BindingEnergy(MassNumber mass_number, ChargeNumber charge_number) {
  /// Weitzsaecker's Mass formula
  FermiUInt N_pairing = (FermiUInt(mass_number) - FermiUInt(charge_number)) % 2;                      /// pairing
  FermiUInt Z_pairing = FermiUInt(charge_number) % 2;

  FermiFloat binding =
          - 15.67 * FermiFloat(mass_number)                                                           /// nuclear volume
          + 17.23 * std::pow(mass_number, 2. / 3.)                                                    /// surface energy
          + 93.15 * ((FermiFloat(mass_number) / 2. - FermiFloat(charge_number))
          * (FermiFloat(mass_number) / 2. - FermiFloat(charge_number))) / FermiFloat(mass_number)     /// asymmetry
          + 0.6984523 * std::pow(FermiUInt(charge_number), 2) * std::pow(mass_number, -1. / 3.);      /// coulomb

  if (N_pairing == Z_pairing) {
    binding += (N_pairing + Z_pairing - 1) * 12.0 / std::sqrt(FermiFloat(mass_number));               /// pairing
  }

  return -binding * CLHEP::MeV;
}

bool NativeNucleiProperties::IsInvalidNuclei(MassNumber mass_number, ChargeNumber charge_number) {
  return mass_number < 1_m || charge_number < 0_c || FermiUInt(charge_number) > FermiUInt(mass_number);
}

} // namespace properties
