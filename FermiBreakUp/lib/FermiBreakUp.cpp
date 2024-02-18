//
// Created by Artem Novikov on 21.02.2023.
//

#include "FermiBreakUp.h"
#include "Utilities/DataTypes.h"

FermiBreakUp::FermiBreakUp() : fermi_configurations_(DefaultConfigurations()) {}

FermiBreakUp::FermiBreakUp(std::unique_ptr<VFermiConfigurations>&& configurations)
    : fermi_configurations_(std::move(configurations)) {}

ParticleVector FermiBreakUp::BreakItUp(const FermiParticle& nucleus) {
  /// CHECK that Excitation Energy > 0
  if (nucleus.GetExcitationEnergy() < 0) {
    return {nucleus};
  }

  /// Total energy of nucleus in nucleus rest frame
  FermiFloat total_energy = nucleus.GetMomentum().m();

  /// Split the nucleus
  auto fragment_split = fermi_configurations_->GenerateSplits(nucleus.GetNucleiData(), total_energy).ChooseSplit();
  if (!fragment_split.has_value()) {
    return {nucleus};
  }

  return ConvertToParticles(nucleus, fragment_split.value());
}

std::unique_ptr<VFermiConfigurations> FermiBreakUp::DefaultConfigurations() {
  return std::make_unique<FermiConfigurations>();
}
