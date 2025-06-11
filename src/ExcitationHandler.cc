//
// Created by Artem Novikov on 17.05.2023.
//

#include <cmath>
#include <CLHEP/Units/PhysicalConstants.h>
#include <G4BosonConstructor.hh>
#include <G4LeptonConstructor.hh>
#include <G4MesonConstructor.hh>
#include <G4BaryonConstructor.hh>
#include <G4IonConstructor.hh>
#include <G4ProcessManager.hh>
#include <G4StateManager.hh>
#include <G4RunManager.hh>

#include "G4LorentzVector.hh"
#include "G4NistManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4Ions.hh"
#include "G4Electron.hh"

#include "G4Evaporation.hh"
#include "G4PhotonEvaporation.hh"
#include "G4StatMF.hh"
#include "FermiBreakUp/AAMCCFermiBreakUp.hh"
#include "ExcitationHandler.hh"

const size_t ExcitationHandler::EvaporationIterationThreshold = 1e3;

const char* ExcitationHandler::ErrorNoModel = "no model was applied, check conditions";

ExcitationHandler::ExcitationHandler()
    : multi_fragmentation_model_(DefaultMultiFragmentation()),
      fermi_break_up_model_(DefaultFermiBreakUp()),
      evaporation_model_(DefaultEvaporation()),
      multi_fragmentation_condition_(DefaultMultiFragmentationCondition()),
      fermi_condition_(DefaultFermiBreakUpCondition()),
      evaporation_condition_(DefaultEvaporationCondition()),
      photon_evaporation_condition_(DefaultPhotonEvaporationCondition()) {
  G4BosonConstructor pCBos;
  pCBos.ConstructParticle();

  G4LeptonConstructor pCLept;
  pCLept.ConstructParticle();

  G4MesonConstructor pCMes;
  pCMes.ConstructParticle();

  G4BaryonConstructor pCBar;
  pCBar.ConstructParticle();

  G4IonConstructor pCIon;
  pCIon.ConstructParticle();

  G4GenericIon* gion = G4GenericIon::GenericIon();
  gion->SetProcessManager(new G4ProcessManager(gion));

  G4StateManager::GetStateManager()->SetNewState(G4State_Init); // To let create ions
  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  G4IonTable* ion_table = partTable->GetIonTable();
  partTable->SetReadiness();
  ion_table->CreateAllIon();
  ion_table->CreateAllIsomer();
}

void ExcitationHandler::CleanUp(G4FragmentVector& v, std::queue<G4Fragment*>& q1, std::queue<G4Fragment*>& q2) {
  for (auto ptr : v) {
    delete ptr;
  }
  while (!q1.empty()) {
    delete q1.front();
    q1.pop();
  }
  while (!q2.empty()) {
    delete q2.front();
    q2.pop();
  }
}

std::vector<G4ReactionProduct> ExcitationHandler::BreakItUp(const G4Fragment& fragment) {
  auto nist = G4NistManager::Instance();
  G4FragmentVector results;
  std::queue<G4Fragment*> evaporation_queue;
  std::queue<G4Fragment*> photon_evaporation_queue;

  /// In case A <= 1 the fragment will not perform any nucleon emission
  auto initial_fragment_ptr = std::make_unique<G4Fragment>(fragment);
  if (IsStable(fragment, nist)) {
    results.push_back(initial_fragment_ptr.release());
  } else {
    if (multi_fragmentation_condition_(fragment)) {
      ApplyMultiFragmentation(std::move(initial_fragment_ptr), results, evaporation_queue);
    } else {
      evaporation_queue.push(initial_fragment_ptr.release());
    }

    for (size_t iteration_count = 0; !evaporation_queue.empty(); ++iteration_count) {
      auto fragment_ptr = std::unique_ptr<G4Fragment>(evaporation_queue.front());
      evaporation_queue.pop();

      /// infinite loop
      if (iteration_count == EvaporationIterationThreshold) {
        /// exception safety
        CleanUp(results, evaporation_queue, photon_evaporation_queue);

        EvaporationError(fragment, *fragment_ptr, iteration_count);
        /// process is dead
      }

      /// FermiBreakUp part
      if (fermi_condition_(*fragment_ptr)) {
        ApplyFermiBreakUp(std::move(fragment_ptr), results, photon_evaporation_queue);
        continue;
      }

      /// Evaporation part
      if (evaporation_condition_(*fragment_ptr)) {
        ApplyEvaporation(std::move(fragment_ptr), results, evaporation_queue);
        continue;
      }

      /// exception safety
      CleanUp(results, evaporation_queue, photon_evaporation_queue);
      throw std::runtime_error(ErrorNoModel);
    }

    /// Photon Evaporation part
    while (!photon_evaporation_queue.empty()) {
      auto fragment_ptr = std::unique_ptr<G4Fragment>(photon_evaporation_queue.front());
      photon_evaporation_queue.pop();

      if (photon_evaporation_condition_(*fragment_ptr)) {
        ApplyPhotonEvaporation(std::move(fragment_ptr), results);
        continue;
      }

      /// exception safety
      CleanUp(results, evaporation_queue, photon_evaporation_queue);
      throw std::runtime_error(ErrorNoModel);
    }
  }

  auto reaction_products = ConvertResults(results);

  for (auto ptr : results) {
    delete ptr;
  }

  return reaction_products;
}

std::unique_ptr<G4VMultiFragmentation> ExcitationHandler::DefaultMultiFragmentation() {
  return std::make_unique<G4StatMF>();
}

std::unique_ptr<G4VFermiBreakUp> ExcitationHandler::DefaultFermiBreakUp() {
  return std::make_unique<AAMCCFermiBreakUp>();
}

std::unique_ptr<G4VEvaporation> ExcitationHandler::DefaultEvaporation() {
  auto evaporation = std::make_unique<G4Evaporation>(DefaultPhotonEvaporation().release());
  evaporation->SetFermiBreakUp(DefaultFermiBreakUp().release());
  return evaporation;
}

std::unique_ptr<G4VEvaporationChannel> ExcitationHandler::DefaultPhotonEvaporation() {
  return std::make_unique<G4PhotonEvaporation>();
}

ExcitationHandler::Condition ExcitationHandler::DefaultMultiFragmentationCondition() {
  static const G4int max_A = 19;
  static const G4int max_Z = 9;

  return [](const G4Fragment& fragment) {
    auto A = fragment.GetA_asInt();
    auto Z = fragment.GetZ_asInt();
    auto Ex = fragment.GetExcitationEnergy();
    if (A < max_A && Z < max_Z) { return false; }
    G4double lower_bound_transition_MF = 3 * CLHEP::MeV;
    G4double upper_bound_transition_MF = 5 * CLHEP::MeV;
    G4double aE = 1 / (2. * (upper_bound_transition_MF - lower_bound_transition_MF));
    G4double E0 = (upper_bound_transition_MF + lower_bound_transition_MF) / 2.;
    G4double w = G4RandFlat::shoot();
    G4double transF = 0.5 * std::tanh((Ex / A - E0) / aE) + 0.5;

    if (Ex < lower_bound_transition_MF * A) { return false; }

    if (w < transF && Ex < upper_bound_transition_MF * A) { return true; }

    if (w > transF && Ex < upper_bound_transition_MF * A) { return false; }

    if (Ex > upper_bound_transition_MF * A) { return true; }

    return false;
  };
}

ExcitationHandler::Condition ExcitationHandler::DefaultFermiBreakUpCondition() {
  return [](const G4Fragment& fragment) {
    return AAMCCFermiBreakUp::IsFermiPossible(fragment.GetZ_asInt(),
                                              fragment.GetA_asInt(),
                                              fragment.GetExcitationEnergy());
  };
}

ExcitationHandler::Condition ExcitationHandler::DefaultEvaporationCondition() {
  return [](const G4Fragment&) { return true; };
}

ExcitationHandler::Condition ExcitationHandler::DefaultPhotonEvaporationCondition() {
  return [](const G4Fragment&) { return true; };
}

bool ExcitationHandler::IsGroundState(const G4Fragment& fragment) const {
  return fragment.GetExcitationEnergy() < stable_threshold_;
}

bool ExcitationHandler::IsStable(const G4Fragment& fragment, const G4NistManager* nist) const {
  return fragment.GetA_asInt() <= 1
      || (IsGroundState(fragment) && nist->GetIsotopeAbundance(fragment.GetZ_asInt(), fragment.GetA_asInt()) > 0);
}

void ExcitationHandler::ApplyMultiFragmentation(std::unique_ptr<G4Fragment> fragment,
                                                G4FragmentVector& results,
                                                std::queue<G4Fragment*>& next_stage) {
  auto fragments = std::unique_ptr<G4FragmentVector>(multi_fragmentation_model_->BreakItUp(*fragment));
  if (fragments == nullptr || fragments->size() <= 1) {
    next_stage.push(fragment.release());
    return;
  }

  GroupFragments(*fragments, results, next_stage);
}

void ExcitationHandler::ApplyFermiBreakUp(std::unique_ptr<G4Fragment> fragment,
                                          G4FragmentVector& results,
                                          std::queue<G4Fragment*>& next_stage) {
  G4FragmentVector fragments;
  fermi_break_up_model_->BreakFragment(&fragments, fragment.get());
  // auto fragments = std::unique_ptr<G4FragmentVector>(fermi_break_up_model_->BreakItUp(fragment.get()))

  if (fragments.size() == 1) {
    next_stage.emplace(fragment.release());
    return;
  }

  GroupFragments(fragments, results, next_stage);
}

void ExcitationHandler::ApplyEvaporation(std::unique_ptr<G4Fragment> fragment,
                                         G4FragmentVector& results,
                                         std::queue<G4Fragment*>& next_stage) {
  G4FragmentVector fragments;
  evaporation_model_->BreakFragment(&fragments, fragment.get());
  // auto fragments = std::unique_ptr<G4FragmentVector>(evaporation_model_->BreakItUp(fragment.get()))

  /// because evaporation adjusts it
  auto fragment_ptr = fragment.release();
  if (fragments.empty() || fragments.back() != fragment_ptr) {
    fragments.emplace_back(fragment_ptr);
  }

  if (fragments.size() == 1) {
    results.emplace_back(fragment_ptr);
    return;
  }

  GroupFragments(fragments, results, next_stage);
}

void ExcitationHandler::ApplyPhotonEvaporation(std::unique_ptr<G4Fragment> fragment, G4FragmentVector& results) {
  /// photon de-excitation only for hot fragments
  if (!IsGroundState(*fragment)) {
    G4FragmentVector fragments;

    evaporation_model_->GetPhotonEvaporation()->BreakUpChain(&fragments, fragment.get());

    for (auto fragment_ptr : fragments) {
      results.emplace_back(fragment_ptr);
    }
  }

  /// primary fragment is kept
  results.emplace_back(fragment.release());
}

void ExcitationHandler::GroupFragments(const G4FragmentVector& fragments,
                                       G4FragmentVector& results,
                                       std::queue<G4Fragment*>& next_stage) {
  auto nist = G4NistManager::Instance();

  for (auto fragment_ptr : fragments) { /// fragment pointers is moved to unique and will be deleted later
    /// gamma, p, n or stable nuclei
    if (IsStable(*fragment_ptr, nist)) {
      results.emplace_back(fragment_ptr);
    } else {
      next_stage.emplace(fragment_ptr);
    }
  }
}

constexpr G4int HashParticle(G4int A, G4int Z) { return A * 1000 + Z; }

G4ParticleDefinition* ExcitationHandler::SpecialParticleDefinition(const G4Fragment& fragment) {
  switch (HashParticle(fragment.GetA_asInt(), fragment.GetZ_asInt())) {
    case HashParticle(0, 0): {
      return G4Gamma::GammaDefinition();
    }

    case HashParticle(0, -1): {
      return G4Electron::ElectronDefinition();
    }

    case HashParticle(-1, 1): {
      return G4PionPlus::PionPlus();
    }

    case HashParticle(-1, -1): {
      return G4PionMinus::PionMinus();
    }

    case HashParticle(-1, 0): {
      return G4PionZero::PionZero();
    }

    case HashParticle(1, 0): {
      return G4Neutron::NeutronDefinition();
    }

    case HashParticle(1, 1): {
      return G4Proton::ProtonDefinition();
    }

    case HashParticle(2, 1): {
      return G4Deuteron::DeuteronDefinition();
    }

    case HashParticle(3, 1): {
      return G4Triton::TritonDefinition();
    }

    case HashParticle(3, 2): {
      return G4He3::He3Definition();
    }

    case HashParticle(4, 2): {
      return G4Alpha::AlphaDefinition();
    }
  }

  return nullptr;
}

std::vector<G4ReactionProduct> ExcitationHandler::ConvertResults(const G4FragmentVector& results) {
  std::vector<G4ReactionProduct> reaction_products;
  reaction_products.reserve(results.size());
  auto ion_table = G4ParticleTable::GetParticleTable()->GetIonTable();

  for (auto& fragment_ptr : results) {
    auto fragment_definition = SpecialParticleDefinition(*fragment_ptr);
    if (fragment_definition == nullptr) {
      auto excitation_energy = fragment_ptr->GetExcitationEnergy();
      auto level = fragment_ptr->GetFloatingLevelNumber();
      if (IsGroundState(*fragment_ptr)) {
        excitation_energy = 0;
        level = 0;
      }
      fragment_definition = ion_table->GetIon(fragment_ptr->GetZ_asInt(), fragment_ptr->GetA_asInt(),
                                              excitation_energy, G4Ions::FloatLevelBase(level));
    }
    /// fragment wasn't found, ground state is created
    if (fragment_definition == nullptr) {
      fragment_definition = ion_table->GetIon(fragment_ptr->GetZ_asInt(), fragment_ptr->GetA_asInt(), 0, noFloat, 0);
      if (fragment_definition == nullptr) {
        throw std::runtime_error("ion table isn't created");
      }
      G4double ion_mass = fragment_definition->GetPDGMass();
      if (fragment_ptr->GetMomentum().e() <= ion_mass) {
        fragment_ptr->SetMomentum(G4LorentzVector(ion_mass));
      } else {
        auto momentum = fragment_ptr->GetMomentum();
        G4double momentum_modulus = std::sqrt((momentum.e() - ion_mass) * (momentum.e() + ion_mass));
        momentum.setVect(momentum.vect().unit() * momentum_modulus);
        fragment_ptr->SetMomentum(momentum);
      }
    }

    reaction_products.emplace_back(fragment_definition);
    reaction_products.back().SetMomentum(fragment_ptr->GetMomentum().vect());
    reaction_products.back().SetTotalEnergy((fragment_ptr->GetMomentum()).e());
    reaction_products.back().SetFormationTime(fragment_ptr->GetCreationTime());
  }

  return reaction_products;
}

void ExcitationHandler::EvaporationError(const G4Fragment& fragment, const G4Fragment& current_fragment, size_t iter) {
  G4ExceptionDescription ed;
  ed << "Infinite loop in the de-excitation module: " << iter
     << " iterations \n"
     << "      Initial fragment: \n" << fragment
     << "\n      Current fragment: \n" << current_fragment;
  G4Exception("ExcitationHandler::BreakItUp", "", FatalException,
              ed, "Stop execution");
}
