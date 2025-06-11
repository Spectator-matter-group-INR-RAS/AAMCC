//
// Created by Artem Novikov on 19.03.2023.
//

#include <gtest/gtest.h>
#include <CLHEP/Units/PhysicalConstants.h>

#include "Utilities/NucleiProperties/NucleiProperties.h"
#include "Utilities/Randomizer.h"

#include "FermiBreakUp.h"
#include "Configurations/CachedFermiConfigurations.h"
#include "Configurations/FastFermiConfigurations.h"
#include "Configurations/FermiConfigurations.h"

using namespace properties;

float CalculateFragmentCount(MassNumber mass, ChargeNumber charge, const Vector3& vec,
                             FermiFloat energy_per_nucleon, size_t tests,
                             std::unique_ptr<VFermiConfigurations>&& configurations = FermiBreakUp::DefaultConfigurations()) {
  auto model = FermiBreakUp(std::move(configurations));
  size_t parts_counter = 0;
  auto energy = energy_per_nucleon * FermiFloat(mass);
  auto total_energy = std::sqrt(std::pow(NucleiProperties()->GetNuclearMass(mass, charge) + energy, 2) + vec.mag2());
  auto mom = LorentzVector(vec, total_energy);
  auto particle = FermiParticle(mass, charge, mom);
  for (size_t i = 0; i < tests; ++i) {
    auto particles = model.BreakItUp(particle);
    parts_counter += particles.size();
  }
  return float(parts_counter) / float(tests);
}

class ConfigurationsFixture : public ::testing::TestWithParam<VFermiConfigurations*> {
 protected:
  VFermiConfigurations* configurations;
};

INSTANTIATE_TEST_SUITE_P(
    FermiBreakUpTests,
    ConfigurationsFixture,
    ::testing::Values(
        new FermiConfigurations(),
        new CachedFermiConfigurations(),
        new FastFermiConfigurations()
    ));

TEST_P(ConfigurationsFixture, CarbonTest) {
  size_t runs = 1e3;
  auto config = ConfigurationsFixture::GetParam()->Clone();

  /// Carbons shouldn't break up [0, 0.7]
  ASSERT_NEAR(1, CalculateFragmentCount(12_m, 6_c, {0, 0, 0}, 0 * CLHEP::MeV, runs, ConfigurationsFixture::GetParam()->Clone()), 0.01);
  ASSERT_NEAR(1, CalculateFragmentCount(12_m, 6_c, {0, 0, 0}, 0.25 * CLHEP::MeV, runs, ConfigurationsFixture::GetParam()->Clone()), 0.01);
  ASSERT_NEAR(1, CalculateFragmentCount(12_m, 6_c, {0, 0, 0}, 0.5 * CLHEP::MeV, runs, ConfigurationsFixture::GetParam()->Clone()), 0.01);
  ASSERT_NEAR(1, CalculateFragmentCount(12_m, 6_c, {0, 0, 0}, 0.7 * CLHEP::MeV, runs, ConfigurationsFixture::GetParam()->Clone()), 0.01);

  /// Carbons should break into 3 parts [1, 1.4]
  ASSERT_NEAR(3, CalculateFragmentCount(12_m, 6_c, {0, 0, 0}, 1 * CLHEP::MeV, runs, ConfigurationsFixture::GetParam()->Clone()), 0.01);
  ASSERT_NEAR(3, CalculateFragmentCount(12_m, 6_c, {0, 0, 0}, 1.25 * CLHEP::MeV, runs, ConfigurationsFixture::GetParam()->Clone()), 0.01);
  ASSERT_NEAR(3, CalculateFragmentCount(12_m, 6_c, {0, 0, 0}, 1.4 * CLHEP::MeV, runs, ConfigurationsFixture::GetParam()->Clone()), 0.01);


  /// Carbons should break into less than 3 parts [1.5, 3.5]
  ASSERT_LE(CalculateFragmentCount(12_m, 6_c, {0, 0, 0}, 1.5 * CLHEP::MeV, runs, ConfigurationsFixture::GetParam()->Clone()), 3);
  ASSERT_LE(CalculateFragmentCount(12_m, 6_c, {0, 0, 0}, 2 * CLHEP::MeV, runs, ConfigurationsFixture::GetParam()->Clone()), 3);
  ASSERT_LE(CalculateFragmentCount(12_m, 6_c, {0, 0, 0}, 3 * CLHEP::MeV, runs, ConfigurationsFixture::GetParam()->Clone()), 3);
  ASSERT_LE(CalculateFragmentCount(12_m, 6_c, {0, 0, 0}, 3.5 * CLHEP::MeV, runs, ConfigurationsFixture::GetParam()->Clone()), 3);

  /// Carbons should break into more than 3 parts [5, ...]
  ASSERT_GE(CalculateFragmentCount(12_m, 6_c, {0, 0, 0}, 5 * CLHEP::MeV, runs, ConfigurationsFixture::GetParam()->Clone()), 3);
  ASSERT_GE(CalculateFragmentCount(12_m, 6_c, {0, 0, 0}, 7 * CLHEP::MeV, runs, ConfigurationsFixture::GetParam()->Clone()), 3);
  ASSERT_GE(CalculateFragmentCount(12_m, 6_c, {0, 0, 0}, 9 * CLHEP::MeV, runs, ConfigurationsFixture::GetParam()->Clone()), 3);
  ASSERT_GE(CalculateFragmentCount(12_m, 6_c, {0, 0, 0}, 20 * CLHEP::MeV, runs, ConfigurationsFixture::GetParam()->Clone()), 3);
}

TEST_P(ConfigurationsFixture, UnstableNucleiTest) {
  size_t runs = 1e3;

  /// protons should break up
  for (int i = 2; i <= 16; ++i) {
    ASSERT_NEAR(i, CalculateFragmentCount(MassNumber(i), ChargeNumber(i), {0, 0, 0}, 0, runs, ConfigurationsFixture::GetParam()->Clone()), 0.01);
  }

  /// neutron should break up
  for (int i = 2; i <= 16; ++i) {
    ASSERT_NEAR(i, CalculateFragmentCount(MassNumber(i), 0_c, {0, 0, 0}, 0, runs, ConfigurationsFixture::GetParam()->Clone()), 0.01);
  }
}

TEST_P(ConfigurationsFixture, MomentumConservation) {
  auto model = FermiBreakUp(ConfigurationsFixture::GetParam()->Clone());
  int seed = 5;
  srand(seed);
  int tries = 500;
  size_t runs = 1e3;
  for (int t = 0; t < tries; ++t) {
    MassNumber mass(rand() % 16 + 1);
    ChargeNumber charge(rand() % (int(mass) + 1));
    FermiFloat energy = (rand() % 1000) * CLHEP::MeV * FermiFloat(mass);
    auto vec = Randomizer::IsotropicVector() * (rand() % 1000) * CLHEP::MeV;
    auto
        total_energy = std::sqrt(std::pow(NucleiProperties()->GetNuclearMass(mass, charge) + energy, 2) + vec.mag2());
    auto mom = LorentzVector(vec, total_energy);
    auto particle = FermiParticle(mass, charge, mom);
    for (size_t i = 0; i < runs; ++i) {
      LorentzVector sum(0, 0, 0, 0);
      auto particles = model.BreakItUp(particle);
      for (auto& fragment : particles) {
        sum += fragment.GetMomentum();
      }
      ASSERT_NEAR(sum.x(), mom.x(), 1e-5);
      ASSERT_NEAR(sum.y(), mom.y(), 1e-5);
      ASSERT_NEAR(sum.z(), mom.z(), 1e-5);
      ASSERT_NEAR(sum.e(), mom.e(), 1e-5);
    }
  }
}

TEST_P(ConfigurationsFixture, BaryonAndChargeConservation) {
  auto model = FermiBreakUp(ConfigurationsFixture::GetParam()->Clone());
  int seed = 5;
  srand(seed);
  int tries = 500;
  size_t runs = 1e3;
  for (int t = 0; t < tries; ++t) {
    MassNumber mass(rand() % 16 + 1);
    ChargeNumber charge(rand() % (int(mass) + 1));
    FermiFloat energy = (rand() % 1000) * CLHEP::MeV * FermiFloat(mass);
    auto vec = Randomizer::IsotropicVector() * (rand() % 1000) * CLHEP::MeV;
    auto
        total_energy = std::sqrt(std::pow(NucleiProperties()->GetNuclearMass(mass, charge) + energy, 2) + vec.mag2());
    auto mom = LorentzVector(vec, total_energy);
    auto particle = FermiParticle(mass, charge, mom);
    for (size_t i = 0; i < runs; ++i) {
      MassNumber fragments_mass_sum(0);
      ChargeNumber fragments_charge_sum(0);
      auto particles = model.BreakItUp(particle);
      for (auto& fragment : particles) {
        fragments_mass_sum = MassNumber(FermiUInt(fragments_mass_sum) + FermiUInt(fragment.GetMassNumber()));
        fragments_charge_sum = ChargeNumber(FermiUInt(fragments_charge_sum) + FermiUInt(fragment.GetChargeNumber()));
      }
      ASSERT_EQ(fragments_mass_sum, mass);
      ASSERT_EQ(fragments_charge_sum, charge);
    }
  }
}
