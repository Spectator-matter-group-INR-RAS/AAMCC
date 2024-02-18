//
// Created by Artem Novikov on 21.05.2023.
//

#include "VFermiBreakUp.h"
#include "PhaseDecay/FermiPhaseSpaceDecay.h"

ParticleVector VFermiBreakUp::ConvertToParticles(const FermiParticle& source_nucleus, const FragmentVector& split) {
  ParticleVector particle_split;
  particle_split.reserve(2 * split.size());

  std::vector<FermiFloat> split_masses;
  split_masses.reserve(split.size());
  for (auto fragment_ptr : split) {
    split_masses.push_back(fragment_ptr->GetTotalEnergy());
  }

  FermiPhaseSpaceDecay phase_sampler;
  auto particles_momentum = phase_sampler.CalculateDecay(source_nucleus.GetMomentum(), split_masses);

  auto boost_vector = source_nucleus.GetMomentum().boostVector();

  /// Go back to the Lab Frame
  for (size_t fragment_idx = 0; fragment_idx < split.size(); ++fragment_idx) {
    ParticleVector fragment_particles = split[fragment_idx]->GetFragments(
        particles_momentum[fragment_idx].boost(boost_vector));

    particle_split.insert(particle_split.end(), std::make_move_iterator(fragment_particles.begin()),
                          std::make_move_iterator(fragment_particles.end()));
  }

  return particle_split;
}