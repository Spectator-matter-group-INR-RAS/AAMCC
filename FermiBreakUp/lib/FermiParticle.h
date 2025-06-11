//
// Created by Artem Novikov on 17.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_PARTICLE_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_PARTICLE_H_

#include <vector>

#include "Utilities/DataTypes.h"

class FermiParticle;

using ParticleVector = std::vector<FermiParticle>;

class FermiParticle {
 public:
  FermiParticle() = delete;

  FermiParticle(const FermiParticle&) = default;

  FermiParticle& operator=(const FermiParticle&) = default;

  FermiParticle(MassNumber mass_number, ChargeNumber charge_number, const LorentzVector& momentum);

  NucleiData GetNucleiData() const;

  MassNumber GetMassNumber() const;

  ChargeNumber GetChargeNumber() const;

  const LorentzVector& GetMomentum() const;

  FermiFloat GetExcitationEnergy() const;

  FermiFloat GetGroundStateMass() const;

  FermiFloat GetBindingEnergy() const;

  bool IsStable() const;

  void SetMassAndCharge(MassNumber mass_number, ChargeNumber charge_number);

  void SetMomentum(const LorentzVector& momentum);

  ~FermiParticle() = default;

 private:
   void CalculateExcitationEnergy();

   void CalculateGroundStateMass();

  MassNumber mass_number_;
  ChargeNumber charge_number_;
  LorentzVector momentum_;

  FermiFloat excitation_energy_ = 0;
  FermiFloat ground_state_mass_ = 0;
};

std::ostream& operator<<(std::ostream&, const FermiParticle&);

#endif //FERMIBREAKUP_MYFERMIBREAKUP_PARTICLE_H_
