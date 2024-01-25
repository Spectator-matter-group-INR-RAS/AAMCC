//
// Created by Artem Novikov on 25.11.2023.
//

#include <G4NistManager.hh>

#include "MyDeexcitationHandler.hh"

MyDeexcitationHandler::MyDeexcitationHandler() :
    pure_neutrons_(DefaultPureNeutrons()),
    pure_neutrons_condition_(DefaultPureNeutronsCondition()) {}

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

  auto nist = G4NistManager::Instance();
  G4FragmentVector results;
  std::queue<G4Fragment*> secondary_decay;

  auto initial_fragment_ptr = std::make_unique<G4Fragment>(fragment);
  if (IsStable(fragment, nist)) {
    results.push_back(initial_fragment_ptr.get());
    return ConvertResults(results);
  }

  auto is_multi_fragmentation = multi_fragmentation_condition_(fragment);
  if (is_multi_fragmentation) {
    ApplyMultiFragmentation(std::move(initial_fragment_ptr), results, secondary_decay);
  } else if (fermi_condition_(fragment)) {
    ApplyFermiBreakUp(std::move(initial_fragment_ptr), results, secondary_decay);
  } else if (evaporation_condition_(fragment)) {
    G4FragmentVector fragments;
    evaporation_model_->BreakFragment(&fragments, initial_fragment_ptr.get());
    results.insert(results.end(), fragments.begin(), fragments.end());
  } else {
    throw std::runtime_error(ErrorNoModel);
  }

  auto reaction_products = ConvertResults(results);

  /// prevent secondary MF
  Condition mf_condition;
  if (is_multi_fragmentation) {
    mf_condition = std::move(multi_fragmentation_condition_);
    multi_fragmentation_condition_ = [](const G4Fragment& frag) { return false; };
  }

  try {
    while (!secondary_decay.empty()) {
      auto fragment_ptr = std::unique_ptr<G4Fragment>(secondary_decay.back());
      secondary_decay.pop();

      auto fragment_reaction_products = BreakItUp(*fragment_ptr);
      for (auto& product : fragment_reaction_products) {
        reaction_products.emplace_back(std::move(product));
      }
    }
  } catch(...) {
    /// exception safety
    if (is_multi_fragmentation) {
      multi_fragmentation_condition_ = std::move(mf_condition);
    }
    for (auto& fragment_ptr : results) {
      delete fragment_ptr;
    }
    throw;
  }

  /// change MF back
  if (is_multi_fragmentation) {
    multi_fragmentation_condition_ = std::move(mf_condition);
  }

  for (auto& fragment_ptr : results) {
    delete fragment_ptr;
  }

  return reaction_products;
}

enum class Models {
  G4 = 0,
  AAMCC = 1,
  NONE = 2
};

Models HashModelName(const G4String& name) {
  if (name == "G4") {
    return Models::G4;
  }

  if (name == "AAMCC") {
    return Models::AAMCC;
  }

  return Models::NONE;
}

std::vector<G4ReactionProduct> MyDeexcitationHandler::BreakUp(const G4Fragment& fragment,
                                                              const G4String& modelName) {
  switch (HashModelName(modelName)) {
    case Models::G4:return G4BreakItUp(fragment);

    case Models::AAMCC:return AAMCCBreakItUp(fragment);

    default: {
      std::cout << "Wrong model name " << modelName << ". G4, ABLAXX, AAMCC or MIX are available \n";
      return AAMCCBreakItUp(fragment);
    }
  }
}

std::unique_ptr<PureNeutrons> MyDeexcitationHandler::DefaultPureNeutrons() {
  return std::make_unique<PureNeutrons>();
}

ExcitationHandler::Condition MyDeexcitationHandler::DefaultPureNeutronsCondition() {
  return [](const G4Fragment& fragment) -> bool {
    return fragment.GetZ_asInt() == 0 && fragment.GetA_asInt() > 0;
  };
}
