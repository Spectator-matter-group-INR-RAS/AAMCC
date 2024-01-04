//
// Created by Artem Novikov on 25.11.2023.
//

#ifndef GRATE_INCLUDE_MYDEEXCITATIONHANDLER_H_
#define GRATE_INCLUDE_MYDEEXCITATIONHANDLER_H_

#include <memory>
#include "G4ParticleTypes.hh"

#include "Abla/include/MyAblaEvaporation.hh"
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

  std::vector<G4ReactionProduct> AblaBreakItUp(const G4Fragment &fragment);

  std::vector<G4ReactionProduct> BreakUp(const G4Fragment &fragment, const G4String& modelName);

  /// parameters setters
  ExcitationHandler& SetAblaEvaporation(std::unique_ptr<MyAblaEvaporation>&& model = DefaultAblaEvaporation()) {
    abla_evaporation_ = std::move(model);
    return *this;
  }

  ExcitationHandler& SetPureNeutrons(std::unique_ptr<PureNeutrons>&& model = DefaultPureNeutrons()) {
    pure_neutrons_ = std::move(model);
    return *this;
  }

  template <class F>
  ExcitationHandler& SetAblaEvaporationCondition(F&& f) {
    abla_evaporation_ = std::forward<F>(f);
    return *this;
  }

  template <class F>
  ExcitationHandler& SetPureNeutronsCondition(F&& f) {
    pure_neutrons_ = std::forward<F>(f);
    return *this;
  }

  /// parameters getters
  std::unique_ptr<MyAblaEvaporation>& GetAblaEvaporation() { return abla_evaporation_; }

  const std::unique_ptr<MyAblaEvaporation>& GetAblaEvaporation() const { return abla_evaporation_; }

  std::unique_ptr<PureNeutrons>& GetPureNeutrons() { return pure_neutrons_; }

  const std::unique_ptr<PureNeutrons>& GetPureNeutrons() const { return pure_neutrons_; }

  Condition& GetAblaEvaporationCondition() { return abla_condition_; }

  const Condition& GetAblaEvaporationCondition() const { return abla_condition_; }

  Condition& GetPureNeutronsCondition() { return pure_neutrons_condition_; }

  const Condition& GetPureNeutronsCondition() const { return pure_neutrons_condition_; }

 private:
  static std::unique_ptr<MyAblaEvaporation> DefaultAblaEvaporation();

  static std::unique_ptr<PureNeutrons> DefaultPureNeutrons();

  static Condition DefaultAblaEvaporationCondition();

  static Condition DefaultPureNeutronsCondition();

  std::unique_ptr<PureNeutrons> pure_neutrons_;
  std::unique_ptr<MyAblaEvaporation> abla_evaporation_;

  Condition pure_neutrons_condition_;
  Condition abla_condition_;
};

#endif //GRATE_INCLUDE_MYDEEXCITATIONHANDLER_H_
