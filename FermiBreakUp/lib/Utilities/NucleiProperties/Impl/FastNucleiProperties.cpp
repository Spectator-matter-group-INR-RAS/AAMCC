//
// Created by Artem Novikov on 18.02.2024.
//

#include "../DataStorage/DefaultNuclearMass.h"

#include "FastNucleiProperties.h"

namespace properties {

FastNucleiProperties::FastNucleiProperties() : FastNucleiProperties(DefaultNuclearMass()) {}

FermiFloat FastNucleiProperties::GetNuclearMass(MassNumber mass_number, ChargeNumber charge_number) const {
#if __cplusplus >= 202002L
  if (FermiInt(mass_number) < FermiInt(charge_number)) [[unlikely]] {
#else
  if (FermiInt(mass_number) < FermiInt(charge_number)) {
#endif
    throw std::runtime_error(
        "invalid nuclei A = " + std::to_string(mass_number) + ", Z = " + std::to_string(charge_number));
  }

  auto slot = GetSlot(mass_number, charge_number);
  if (slot < nuclei_masses_.size() && nuclei_masses_[slot].is_cached) {
    return nuclei_masses_[slot].mass;
  }

#ifdef _DEBUG
  std::cerr << "Unknown particle A = " + std::to_string(mass_number) + ", Z = " + std::to_string(charge_number) << '\n';
#endif
  return NuclearMass(mass_number, charge_number);
}

FermiFloat FastNucleiProperties::GetNuclearMass(MassNumber mass_number, ChargeNumber charge_number) {
#if __cplusplus >= 202002L
  if (FermiInt(mass_number) < FermiInt(charge_number)) [[unlikely]] {
#else
  if (FermiInt(mass_number) < FermiInt(charge_number)) {
#endif
    throw std::runtime_error(
        "invalid nuclei A = " + std::to_string(mass_number) + ", Z = " + std::to_string(charge_number));
  }

  auto slot = GetSlot(mass_number, charge_number);
  if (slot < nuclei_masses_.size() && nuclei_masses_[slot].is_cached) {
    return nuclei_masses_[slot].mass;
  }

#ifdef _DEBUG
  std::cerr << "Unknown particle A = " + std::to_string(mass_number) + ", Z = " + std::to_string(charge_number) << '\n';
#endif
  if (slot >= nuclei_masses_.size()) {
    nuclei_masses_.resize(slot + 1);
  }
  nuclei_masses_[slot] = MassData{.mass=NuclearMass(mass_number, charge_number), .is_valid=false, .is_cached=true};
  return nuclei_masses_[slot].mass;
}

bool FastNucleiProperties::IsStable(MassNumber mass_number, ChargeNumber charge_number) const {
#if __cplusplus >= 202002L
  if (IsInvalidNuclei(mass_number, charge_number)) [[unlikely]] {
#else
  if (IsInvalidNuclei(mass_number, charge_number)) {
#endif
    PrintInvalidNuclei(mass_number, charge_number);
    return false;
  }

  auto slot = GetSlot(mass_number, charge_number);
  if (slot >= nuclei_masses_.size() || !nuclei_masses_[slot].is_valid) {
    return false;
  }

  return true;
}

void FastNucleiProperties::AddMass(MassNumber mass_number, ChargeNumber charge_number, FermiFloat mass) {
  auto slot = GetSlot(mass_number, charge_number);
  if (slot >= nuclei_masses_.size()) {
    nuclei_masses_.resize(slot + 1);
  }
  nuclei_masses_[slot] = MassData{.mass=mass, .is_valid=true, .is_cached=true};
}

void FastNucleiProperties::AddMass(NucleiData nuclei_data, FermiFloat mass) {
  return AddMass(nuclei_data.mass_number, nuclei_data.charge_number, mass);
}

size_t FastNucleiProperties::GetSlot(MassNumber mass_number, ChargeNumber charge_number) {
  auto mass = FermiInt(mass_number);
  auto charge = FermiInt(charge_number);
  return (mass * (mass + 1)) / 2 + charge;
}

} // properties
