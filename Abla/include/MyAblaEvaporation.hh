//
// Created by Artem Novikov on 25.11.2023.
//

#ifndef GRATE_ABLA_INCLUDE_MYABLAEVAPORATION_H_
#define GRATE_ABLA_INCLUDE_MYABLAEVAPORATION_H_

#include <memory>

#include "globals.hh"

#include "G4VPreCompoundModel.hh"
#include "G4ReactionProduct.hh"
#include "G4Fragment.hh"
#include "G4HadFinalState.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4Abla.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"

class MyAblaEvaporation {
 public:
  MyAblaEvaporation();

  MyAblaEvaporation(const MyAblaEvaporation&) = delete;

  MyAblaEvaporation(MyAblaEvaporation&&) = default;

  MyAblaEvaporation& operator=(const MyAblaEvaporation&) = delete;

  MyAblaEvaporation& operator=(MyAblaEvaporation&&) = default;

  ~MyAblaEvaporation() = default;

  G4FragmentVector BreakItUp(const G4Fragment& fragment);

  void SetFreezeOutT(G4double T) {
    T_freeze_out = T;
    abla_model_->SetFreezeOutT(T);
  };

  G4double GetFreezeOutT() const { return T_freeze_out; }

 private:
  G4VarNtp abla_result_;
  G4Volant volant_;
  std::unique_ptr<G4Abla> abla_model_;

  G4int event_counter_ = 0;

  G4double T_freeze_out = 1e100;
};

#endif //GRATE_ABLA_INCLUDE_MYABLAEVAPORATION_H_
