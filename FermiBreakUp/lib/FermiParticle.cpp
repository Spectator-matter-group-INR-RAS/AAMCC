//
// Created by Artem Novikov on 17.02.2023.
//

#include <iomanip>
#include <CLHEP/Units/PhysicalConstants.h>

#include "Utilities/NucleiProperties/NucleiProperties.h"
#include "FermiParticle.h"

FermiParticle::FermiParticle(MassNumber mass_number, ChargeNumber charge_number, const LorentzVector& momentum)
    : mass_number_(mass_number), charge_number_(charge_number), momentum_(momentum) {
  CalculateGroundStateMass();
  CalculateExcitationEnergy();
}

NucleiData FermiParticle::GetNucleiData() const {
  return {mass_number_, charge_number_};
}


MassNumber FermiParticle::GetMassNumber() const {
  return mass_number_;
}

ChargeNumber FermiParticle::GetChargeNumber() const {
  return charge_number_;
}

const LorentzVector& FermiParticle::GetMomentum() const {
  return momentum_;
}

FermiFloat FermiParticle::GetExcitationEnergy() const {
  return excitation_energy_;
}

FermiFloat FermiParticle::GetGroundStateMass() const {
  return ground_state_mass_;
}

FermiFloat FermiParticle::GetBindingEnergy() const {
  return (FermiUInt(mass_number_) - FermiUInt(charge_number_)) * CLHEP::neutron_mass_c2
      + FermiFloat(charge_number_) * CLHEP::proton_mass_c2 - ground_state_mass_;
}

bool FermiParticle::IsStable() const {
  return excitation_energy_ <= 0;
}

void FermiParticle::SetMassAndCharge(MassNumber mass_number, ChargeNumber charge_number) {
  mass_number_ = mass_number;
  charge_number_ = charge_number;
  CalculateGroundStateMass();
}

void FermiParticle::SetMomentum(const LorentzVector& momentum) {
  momentum_ = momentum;
  CalculateExcitationEnergy();
}

void FermiParticle::CalculateExcitationEnergy() {
  excitation_energy_ = momentum_.mag() - ground_state_mass_;
  if (excitation_energy_ < 0) {
//    if (excitation_energy_ < -10 * CLHEP::eV) {
//      throw std::runtime_error("Excitation Energy is negative");
//    } /// TODO sometimes is raised in collaboration with other models
    excitation_energy_ = 0;
  }
}

void FermiParticle::CalculateGroundStateMass() {
  ground_state_mass_ = properties::NucleiProperties().GetNuclearMass(mass_number_, charge_number_);
}

std::ostream& operator<<(std::ostream& out, const FermiParticle& particle) {
  auto old_float_field = out.flags();
  out.setf(std::ios::floatfield);

  out << "Fragment: A = " << std::setw(3) << particle.GetMassNumber()
      << ", Z = " << std::setw(3) << particle.GetChargeNumber();
  out.setf(std::ios::scientific, std::ios::floatfield);

  // Store user's precision setting and reset to (3) here: back-compatibility
  std::streamsize old_user_precision = out.precision();

  out << std::setprecision(3)
      << ", U = " << particle.GetExcitationEnergy() / CLHEP::MeV
      << " MeV  IsStable= " << particle.IsStable() << std::endl
      << "          P = ("
      << particle.GetMomentum().x() / CLHEP::MeV << ","
      << particle.GetMomentum().y() / CLHEP::MeV << ","
      << particle.GetMomentum().z() / CLHEP::MeV
      << ") MeV   E = "
      << particle.GetMomentum().t() / CLHEP::MeV << " MeV"
      << std::endl;

  /// What about Angular momentum???

  out.setf(old_float_field, std::ios::floatfield);
  out.precision(old_user_precision);

  return out;
}
