//
// Created by Artem Novikov on 20.05.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_AAMCCFERMIBREAKUP_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_AAMCCFERMIBREAKUP_H_

#include <memory>

#include "G4VFermiBreakUp.hh"
#include "VFermiBreakUp.h"

class AAMCCFermiBreakUp : public G4VFermiBreakUp {
 public:
  AAMCCFermiBreakUp();

  AAMCCFermiBreakUp(std::unique_ptr<VFermiBreakUp>&& model);

  void Initialise() override;

  void BreakFragment(G4FragmentVector* fragments_ptr, G4Fragment* fragment) override;

  static G4bool IsFermiPossible(G4int Z, G4int A, G4double excitation_energy);

  G4bool IsApplicable(G4int Z, G4int A, G4double excitation_energy) const override;

 private:
  std::unique_ptr<VFermiBreakUp> fermi_model_;
};

#endif //FERMIBREAKUP_MYFERMIBREAKUP_AAMCCFERMIBREAKUP_H_
