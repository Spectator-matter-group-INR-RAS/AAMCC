//
// Created by Artem Novikov on 17.02.2023.
//

#include "StableFermiFragment.h"

ParticleVector StableFermiFragment::GetFragment(const LorentzVector& momentum) const {
  return {FermiParticle(GetMassNumber(), GetChargeNumber(), momentum)};
}
