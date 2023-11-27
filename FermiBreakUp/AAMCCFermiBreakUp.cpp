//
// Created by Artem Novikov on 20.05.2023.
//

#include <vector>

#include "G4Fragment.hh"
#include "AAMCCFermiBreakUp.h"
#include "FermiBreakUp.h"

AAMCCFermiBreakUp::AAMCCFermiBreakUp() : fermi_model_(std::make_unique<FermiBreakUp>()) {}

AAMCCFermiBreakUp::AAMCCFermiBreakUp(std::unique_ptr<VFermiBreakUp>&& model) : fermi_model_(std::move(model)) {}

void AAMCCFermiBreakUp::BreakFragment(G4FragmentVector* fragments_ptr, G4Fragment* fragment) {
  auto results = fermi_model_->BreakItUp(FermiParticle(MassNumber(fragment->GetA_asInt()),
                                                      ChargeNumber(fragment->GetZ_asInt()),
                                                      fragment->GetMomentum()));

  for (auto& particle : results) {
    fragments_ptr->push_back(new G4Fragment(G4int(particle.GetMassNumber()),
                                            G4int(particle.GetChargeNumber()),
                                            G4LorentzVector(particle.GetMomentum())));
  }
}

void AAMCCFermiBreakUp::Initialise() {}

G4bool AAMCCFermiBreakUp::IsFermiPossible(G4int Z, G4int A, [[maybe_unused]] G4double excitation_energy) {
  return Z < 9 && A < 19;
}

G4bool AAMCCFermiBreakUp::IsApplicable(G4int Z, G4int A, G4double excitation_energy) const {
  return IsFermiPossible(Z, A, excitation_energy);
}
