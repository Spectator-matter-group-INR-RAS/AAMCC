//
// Created by Artem Novikov on 25.11.2023.
//

#include <G4NistManager.hh>

#include "MyDeexcitationHandler.hh"

MyDeexcitationHandler::MyDeexcitationHandler() :
    pure_neutrons_(DefaultPureNeutrons()),
    abla_evaporation_(DefaultAblaEvaporation()),
    pure_neutrons_condition_(DefaultPureNeutronsCondition()),
    abla_condition_(DefaultAblaEvaporationCondition()) {}

std::vector<G4ReactionProduct> MyDeexcitationHandler::G4BreakItUp(const G4Fragment& fragment) {
  if (pure_neutrons_condition_(fragment)) {
    return BreakUpPureNeutrons(fragment);
  }
  return ExcitationHandler::BreakItUp(fragment);
}

std::vector<G4ReactionProduct> MyDeexcitationHandler::BreakUpPureNeutrons(const G4Fragment& fragment) {
  auto fragments = pure_neutrons_->BreakItUp(fragment);
  auto reaction_products = ConvertResults(fragments);

  for (auto fragment_ptr : fragments) {
    delete fragment_ptr;
  }
  return reaction_products;
}

std::vector<G4ReactionProduct> MyDeexcitationHandler::AAMCCBreakItUp(const G4Fragment& fragment) {
  if (pure_neutrons_condition_(fragment)) {
    return BreakUpPureNeutrons(fragment);
  }

  if (abla_evaporation_->GetFreezeOutT() < 0) {
    abla_evaporation_->SetFreezeOutT(1e100);
  }

  auto nist = G4NistManager::Instance();
  G4FragmentVector results;
  std::queue<G4Fragment*> evaporation_queue;
  std::queue<G4Fragment*> secondary_evaporation_queue;

  auto initial_fragment_ptr = std::make_unique<G4Fragment>(fragment);
  if (IsStable(fragment, nist)) {
    results.push_back(initial_fragment_ptr.release());
  } else {
    auto is_multi_fragmentation = multi_fragmentation_condition_(fragment);
    if (is_multi_fragmentation) {
      ApplyMultiFragmentation(std::move(initial_fragment_ptr), results, evaporation_queue);
    } else if (fermi_condition_(fragment)) {
      ApplyFermiBreakUp(std::move(initial_fragment_ptr), results, evaporation_queue);
    } else if (abla_condition_(fragment)) {
      auto fragments = abla_evaporation_->BreakItUp(fragment);
      results.insert(results.end(), fragments.begin(), fragments.end());
    } else {
      throw std::runtime_error(ErrorNoModel);
    }

    for (size_t iteration_count = 0; !evaporation_queue.empty(); ++iteration_count) {
      auto fragment_ptr = std::unique_ptr<G4Fragment>(evaporation_queue.front());
      evaporation_queue.pop();

      /// infinite loop
      if (iteration_count == EvaporationIterationThreshold) {
        /// exception safety
        CleanUp(results, evaporation_queue, secondary_evaporation_queue);

        EvaporationError(fragment, *fragment_ptr, iteration_count);
        /// process is dead
      }

      if (fermi_condition_(*fragment_ptr) && is_multi_fragmentation) {
        ApplyFermiBreakUp(std::move(fragment_ptr), results, secondary_evaporation_queue);
        continue;
      }

      if (abla_condition_(*fragment_ptr)) {
        auto fragments = abla_evaporation_->BreakItUp(*fragment_ptr);
        results.insert(results.end(), fragments.begin(), fragments.end());
        continue;
      }

      /// exception safety
      CleanUp(results, evaporation_queue, secondary_evaporation_queue);
      throw std::runtime_error(ErrorNoModel);
    }
  }

  while (!secondary_evaporation_queue.empty()) {
    auto fragment_ptr = std::unique_ptr<G4Fragment>(secondary_evaporation_queue.front());
    secondary_evaporation_queue.pop();
    auto fragments = abla_evaporation_->BreakItUp(*fragment_ptr);
    results.insert(results.end(), fragments.begin(), fragments.end());
  }

  auto reaction_products = ConvertResults(results);

  CleanUp(results, evaporation_queue, secondary_evaporation_queue);

  return reaction_products;
}

std::vector<G4ReactionProduct> MyDeexcitationHandler::AblaBreakItUp(const G4Fragment& fragment) {
  if (abla_evaporation_->GetFreezeOutT() > 0) {
    abla_evaporation_->SetFreezeOutT(-6.5); // Let ABLAXX decide freeze-out T
  }
  auto fragments = abla_evaporation_->BreakItUp(fragment);
  auto reaction_products = ConvertResults(fragments);

  for (auto fragment_ptr : fragments) {
    delete fragment_ptr;
  }

  return reaction_products;
}

enum class Models {
  G4 = 0,
  ABLAXX = 1,
  AAMCC = 2,
  MIX = 3,
  NONE = 4
};

Models HashModelName(const G4String& name) {
  if (name == "G4") {
    return Models::G4;
  }

  if (name == "ABLAXX") {
    return Models::ABLAXX;
  }

  if (name == "AAMCC") {
    return Models::AAMCC;
  }

  if (name == "MIX") {
    return Models::MIX;
  }

  return Models::NONE;
}

std::vector<G4ReactionProduct> MyDeexcitationHandler::BreakUp(const G4Fragment& fragment,
                                                              const G4String& modelName) {
  switch (HashModelName(modelName)) {
    case Models::G4:return G4BreakItUp(fragment);

    case Models::ABLAXX:return AblaBreakItUp(fragment);

    case Models::AAMCC:return AAMCCBreakItUp(fragment);

    case Models::MIX: {
      const int OptionsCount = static_cast<int>(Models::MIX);
      auto model = Models(rand() % OptionsCount);
      if (model == Models::G4) return G4BreakItUp(fragment);
      if (model == Models::ABLAXX) return AblaBreakItUp(fragment);
      return AAMCCBreakItUp(fragment);
    }
    default: {
      std::cout << "Wrong model name " << modelName << ". G4, ABLAXX, AAMCC or MIX are available \n";
      return AAMCCBreakItUp(fragment);
    }
  }
}

std::unique_ptr<MyAblaEvaporation> MyDeexcitationHandler::DefaultAblaEvaporation() {
  return std::make_unique<MyAblaEvaporation>();
}

std::unique_ptr<PureNeutrons> MyDeexcitationHandler::DefaultPureNeutrons() {
  return std::make_unique<PureNeutrons>();
}

ExcitationHandler::Condition MyDeexcitationHandler::DefaultAblaEvaporationCondition() {
  return [](const G4Fragment&) -> bool {
    return true;
  };
}

ExcitationHandler::Condition MyDeexcitationHandler::DefaultPureNeutronsCondition() {
  return [](const G4Fragment& fragment) -> bool {
    return fragment.GetZ_asInt() == 0 && fragment.GetA_asInt() > 0;
  };
}
