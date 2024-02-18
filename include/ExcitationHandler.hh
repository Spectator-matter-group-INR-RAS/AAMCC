//
// Created by Artem Novikov on 17.05.2023.
//

#ifndef FERMIBREAKUP_HANDLER_EXCITATIONHANDLER_H_
#define FERMIBREAKUP_HANDLER_EXCITATIONHANDLER_H_

#include <functional>
#include <memory>
#include <vector>
#include <queue>

#include "G4Fragment.hh"
#include "G4ReactionProductVector.hh"
#include "G4IonTable.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4Fragment.hh"
#include "G4NistManager.hh"

#include "G4VMultiFragmentation.hh"
#include "G4VEvaporation.hh"
#include "G4VFermiBreakUp.hh"

class ExcitationHandler {
 public:
  using Condition = std::function<bool(const G4Fragment&)>;

  ExcitationHandler();

  ExcitationHandler(const ExcitationHandler&) = delete;

  ExcitationHandler(ExcitationHandler&&) = default;

  ~ExcitationHandler();

  ExcitationHandler& operator=(const ExcitationHandler&) = delete;

  ExcitationHandler& operator=(ExcitationHandler&&) = default;

  std::vector<G4ReactionProduct> BreakItUp(const G4Fragment& fragment);

  /// parameters setters
  ExcitationHandler& SetMultiFragmentation(std::unique_ptr<G4VMultiFragmentation>&& model = DefaultMultiFragmentation()) {
    multi_fragmentation_model_ = std::move(model);
    return *this;
  }

  ExcitationHandler& SetFermiBreakUp(std::unique_ptr<G4VFermiBreakUp>&& model = DefaultFermiBreakUp(),
                                     bool set_models = true) {
    fermi_break_up_model_ = std::move(model);
    if (set_models) {
      evaporation_model_->SetFermiBreakUp(fermi_break_up_model_.get());
    }
    return *this;
  }

  ExcitationHandler& SetEvaporation(std::unique_ptr<G4VEvaporation>&& model = DefaultEvaporation(),
                                    bool set_models = true) {
    evaporation_model_ = std::move(model);
    if (set_models) {
      evaporation_model_->SetPhotonEvaporation(photon_evaporation_model_.get());
      evaporation_model_->SetFermiBreakUp(fermi_break_up_model_.get());
    }
    return *this;
  }

  ExcitationHandler& SetPhotonEvaporation(std::unique_ptr<G4VEvaporationChannel>&& model = DefaultPhotonEvaporation(),
                                          bool set_models = true) {
    photon_evaporation_model_ = std::move(model);
    if (set_models) {
      evaporation_model_->SetPhotonEvaporation(photon_evaporation_model_.get());
    }
    return *this;
  }

  template <class F>
  ExcitationHandler& SetMultiFragmentationCondition(F&& f) {
    multi_fragmentation_condition_ = std::forward<F>(f);
    return *this;
  }

  ExcitationHandler& SetMultiFragmentationCondition() {
    return SetMultiFragmentationCondition(DefaultMultiFragmentationCondition());
  }

  template <class F>
  ExcitationHandler& SetFermiBreakUpCondition(F&& f) {
    fermi_condition_ = std::forward<F>(f);
    return *this;
  }

  ExcitationHandler& SetFermiBreakUpCondition() {
    return SetFermiBreakUpCondition(DefaultFermiBreakUpCondition());
  }

  template <class F>
  ExcitationHandler& SetEvaporationCondition(F&& f) {
    evaporation_condition_ = std::forward<F>(f);
    return *this;
  }

  ExcitationHandler& SetEvaporationCondition() {
    return SetEvaporationCondition(DefaultEvaporationCondition());
  }

  template <class F>
  ExcitationHandler& SetPhotonEvaporationCondition(F&& f) {
    photon_evaporation_condition_ = std::forward<F>(f);
    return *this;
  }

  ExcitationHandler& SetPhotonEvaporationCondition() {
    return SetPhotonEvaporationCondition(DefaultPhotonEvaporationCondition());
  }

  ExcitationHandler& SetStableThreshold(Float threshold) {
    stable_threshold_ = threshold;
    return *this;
  }

  /// parameters getters
  std::unique_ptr<G4VMultiFragmentation>& GetMultiFragmentation() { return multi_fragmentation_model_; }

  const std::unique_ptr<G4VMultiFragmentation>& GetMultiFragmentation() const { return multi_fragmentation_model_; }

  std::unique_ptr<G4VFermiBreakUp>& GetFermiBreakUp() { return fermi_break_up_model_; }

  const std::unique_ptr<G4VFermiBreakUp>& GetFermiBreakUp() const { return fermi_break_up_model_; }

  std::unique_ptr<G4VEvaporation>& GetEvaporation() { return evaporation_model_; }

  const std::unique_ptr<G4VEvaporation>& GetEvaporation() const { return evaporation_model_; }

  Condition& GetMultiFragmentationCondition() { return multi_fragmentation_condition_; }

  const Condition& GetMultiFragmentationCondition() const { return multi_fragmentation_condition_; }

  Condition& GetFermiBreakUpCondition() { return fermi_condition_; }

  const Condition& GetFermiBreakUpCondition() const { return fermi_condition_; }

  Condition& GetEvaporationCondition() { return evaporation_condition_; }

  const Condition& GetEvaporationCondition() const { return evaporation_condition_; }

  Condition& GetPhotonEvaporationCondition() { return photon_evaporation_condition_; }

  const Condition& GetPhotonEvaporationCondition() const { return photon_evaporation_condition_; }

  Float GetStableThreshold() const { return stable_threshold_; }

 protected:
  static G4ParticleDefinition* SpecialParticleDefinition(const G4Fragment& fragment);

  static void EvaporationError(const G4Fragment& fragment, const G4Fragment& current_fragment, size_t iter);

  /// default models and conditions
  static std::unique_ptr<G4VMultiFragmentation> DefaultMultiFragmentation();

  static std::unique_ptr<G4VFermiBreakUp> DefaultFermiBreakUp();

  static std::unique_ptr<G4VEvaporation> DefaultEvaporation();

  static std::unique_ptr<G4VEvaporationChannel> DefaultPhotonEvaporation();

  static Condition DefaultMultiFragmentationCondition();

  static Condition DefaultFermiBreakUpCondition();

  static Condition DefaultEvaporationCondition();

  static Condition DefaultPhotonEvaporationCondition();

  static void CleanUp(G4FragmentVector& v, std::queue<G4Fragment*>& q1, std::queue<G4Fragment*>& q2);

  bool IsGroundState(const G4Fragment& fragment) const;

  bool IsStable(const G4Fragment& fragment, const G4NistManager* nist) const;

  void ApplyMultiFragmentation(std::unique_ptr<G4Fragment> fragment, G4FragmentVector& results,
                               std::queue<G4Fragment*>& next_stage);

  void ApplyFermiBreakUp(std::unique_ptr<G4Fragment> fragment, G4FragmentVector& results,
                         std::queue<G4Fragment*>& next_stage);

  void ApplyEvaporation(std::unique_ptr<G4Fragment> fragment, G4FragmentVector& results,
                        std::queue<G4Fragment*>& next_stage);

  void ApplyPhotonEvaporation(std::unique_ptr<G4Fragment> fragment, G4FragmentVector& results);

  void GroupFragments(const G4FragmentVector& fragments, G4FragmentVector& results,
                      std::queue<G4Fragment*>& next_stage);

  std::vector<G4ReactionProduct> ConvertResults(const G4FragmentVector& results);

  std::unique_ptr<G4VMultiFragmentation> multi_fragmentation_model_;
  std::unique_ptr<G4VFermiBreakUp> fermi_break_up_model_;
  std::unique_ptr<G4VEvaporationChannel> photon_evaporation_model_;
  std::unique_ptr<G4VEvaporation> evaporation_model_;

  Condition multi_fragmentation_condition_;
  Condition fermi_condition_;
  Condition photon_evaporation_condition_;
  Condition evaporation_condition_;

  Float stable_threshold_ = 0;

  static const size_t EvaporationIterationThreshold;

  static const char* ErrorNoModel;
};

#endif //FERMIBREAKUP_HANDLER_EXCITATIONHANDLER_H_
