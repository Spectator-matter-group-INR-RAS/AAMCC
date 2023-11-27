//
// Created by Artem Novikov on 26.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_TABLEVALUES_NUCLEIPROPERTIESTABLEAME12_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_TABLEVALUES_NUCLEIPROPERTIESTABLEAME12_H_

#include <vector>

#include "Utilities/DataTypes.h"

class NucleiPropertiesTableAME12 {
 public:
  NucleiPropertiesTableAME12();

  FermiFloat GetNuclearMass(MassNumber mass_number, ChargeNumber charge_number) const;

  FermiFloat GetMassExcess(MassNumber mass_number, ChargeNumber charge_number) const;

  FermiFloat GetBindingEnergy(MassNumber mass_number, ChargeNumber charge_number) const;

  FermiFloat GetBetaDecayEnergy(MassNumber mass_number, ChargeNumber charge_number) const;

  bool ContainsParticle(MassNumber mass_number, ChargeNumber charge_number) const;

  FermiUInt GetMaxMass() const;

  FermiUInt GetMaxCharge() const;

  ~NucleiPropertiesTableAME12() = default;

 private:
  size_t GetIndex(MassNumber mass_number, ChargeNumber charge_number) const;

  FermiFloat GetAtomicMass(MassNumber mass_number, ChargeNumber charge_number) const;

  bool VerifyNuclei(MassNumber mass_number, ChargeNumber charge_number) const;

  static const MassNumber MaxMassNumber;

  static const ChargeNumber MaxChargeNumber;

  static const size_t ParticleCount;

  static const std::vector<FermiFloat> MassExcess;

  static const std::vector<FermiFloat> BetaEnergy;

  static const std::vector<size_t> IndexCharge;

  static const std::vector<size_t> IndexMass;

  static const size_t ShortTableCount;

  static const std::vector<size_t> ShortTable;

  static std::vector<FermiFloat>* electron_mass_;
};

std::ostream& operator<<(std::ostream& out, const NucleiPropertiesTableAME12& table);

#endif //FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_TABLEVALUES_NUCLEIPROPERTIESTABLEAME12_H_
