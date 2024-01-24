//
// Created by Artem Novikov on 25.11.2023.
//

#include "../include/AblaEvaporation.hh"

AblaEvaporation::AblaEvaporation() : abla_model_(std::make_unique<G4Abla>(&volant_, &abla_result_)) {
  abla_model_->initEvapora();
  abla_model_->SetParameters();
  abla_model_->SetFreezeOutT(T_freeze_out);
}

G4FragmentVector AblaEvaporation::BreakItUp(const G4Fragment& fragment) {
  volant_.clear();
  abla_result_.clear();

  ++event_counter_;

  abla_model_->DeexcitationAblaxx(fragment.GetA_asInt(),
                                  fragment.GetZ_asInt(),
                                  fragment.GetExcitationEnergy(),
                                  fragment.GetAngularMomentum().mag() / CLHEP::hbar_Planck,
                                  fragment.GetMomentum().x(),
                                  fragment.GetMomentum().y(),
                                  fragment.GetMomentum().z(),
                                  event_counter_);

  G4FragmentVector results;
  results.reserve(abla_result_.ntrack);

  for (auto i = 0; i < abla_result_.ntrack; ++i) {
    auto A = abla_result_.avv[i];
    auto Z = abla_result_.zvv[i];
    auto additional_energy = (A == 0 && Z == 0) ? 0 : G4NucleiProperties::GetNuclearMass(A, Z);
    results.push_back(new G4Fragment(A, Z, G4LorentzVector(abla_result_.pxlab[i],
                                                           abla_result_.pylab[i],
                                                           abla_result_.pzlab[i],
                                                           abla_result_.enerj[i] + additional_energy)));
  }

  return results;
}
