//
// Created by Artem Novikov on 25.11.2023.
//

#ifndef GRATE_INCLUDE_MYDEEXCITATIONHANDLER_H_
#define GRATE_INCLUDE_MYDEEXCITATIONHANDLER_H_

#include <memory>
#include "G4ParticleTypes.hh"

#include "PureNeutrons/include/PureNeutrons.hh"

#include "ExcitationHandler.hh"

class DeexcitationHandler : public ExcitationHandler {
 public:
  DeexcitationHandler();

  DeexcitationHandler(const DeexcitationHandler&) = delete;

  DeexcitationHandler(DeexcitationHandler&&) = default;

  DeexcitationHandler& operator=(const DeexcitationHandler&) = delete;

  DeexcitationHandler& operator=(DeexcitationHandler&&) = default;

  std::vector<G4ReactionProduct> G4BreakItUp(const G4Fragment &fragment);

  std::vector<G4ReactionProduct> BreakUpPureNeutrons(const G4Fragment &fragment);

  std::vector<G4ReactionProduct> AAMCCBreakItUp(const G4Fragment &fragment);

  std::vector<G4ReactionProduct> BreakUp(const G4Fragment &fragment, const G4String& modelName);

  /// parameters setters

  ExcitationHandler& SetPureNeutrons(std::unique_ptr<PureNeutrons>&& model = DefaultPureNeutrons()) {
    pure_neutrons_ = std::move(model);
    return *this;
  }

  template <class F>
  ExcitationHandler& SetPureNeutronsCondition(F&& f) {
    pure_neutrons_ = std::forward<F>(f);
    return *this;
  }

  /// parameters getters
  std::unique_ptr<PureNeutrons>& GetPureNeutrons() { return pure_neutrons_; }

  const std::unique_ptr<PureNeutrons>& GetPureNeutrons() const { return pure_neutrons_; }

  Condition& GetPureNeutronsCondition() { return pure_neutrons_condition_; }

  const Condition& GetPureNeutronsCondition() const { return pure_neutrons_condition_; }

 private:
  static std::unique_ptr<PureNeutrons> DefaultPureNeutrons();

  static Condition DefaultPureNeutronsCondition();

  std::unique_ptr<PureNeutrons> pure_neutrons_;

  Condition pure_neutrons_condition_;
};

#endif //GRATE_INCLUDE_DEEXCITATIONHANDLER_H_
