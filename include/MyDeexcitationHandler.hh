//
// Created by Artem Novikov on 25.11.2023.
//

#ifndef GRATE_INCLUDE_MYMyDEEXCITATIONHANDLER_H_
#define GRATE_INCLUDE_MYMyDEEXCITATIONHANDLER_H_

#include <memory>
#include "G4ParticleTypes.hh"

#include "PureNeutrons/include/PureNeutrons.hh"

#include "ExcitationHandler.hh"

class MyDeexcitationHandler : public ExcitationHandler {
 public:
  MyDeexcitationHandler();

  MyDeexcitationHandler(const MyDeexcitationHandler&) = delete;

  MyDeexcitationHandler(MyDeexcitationHandler&&) = default;

  MyDeexcitationHandler& operator=(const MyDeexcitationHandler&) = delete;

  MyDeexcitationHandler& operator=(MyDeexcitationHandler&&) = default;

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

#endif //GRATE_INCLUDE_MYDEEXCITATIONHANDLER_H_
