//
// Created by Artem Novikov on 18.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_TABLEVALUES_NUCLEIPROPERTIESTABLE_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_TABLEVALUES_NUCLEIPROPERTIESTABLE_H_

#include <vector>

#include "Utilities/DataTypes.h"

class NucleiPropertiesTable {
 public:
  NucleiPropertiesTable();

  FermiFloat GetNuclearMass(MassNumber mass_number, ChargeNumber charge_number) const;

  bool ContainsParticle(MassNumber mass_number, ChargeNumber charge_number) const;

  FermiUInt GetMaxMass() const;

  FermiUInt GetMaxCharge() const;

  ~NucleiPropertiesTable() = default;

 private:
  static FermiFloat ElectronicBindingEnergy(ChargeNumber charge_number);

  size_t GetIndex(MassNumber mass_number, ChargeNumber charge_number) const;

  FermiFloat GetAtomicMass(MassNumber mass_number, ChargeNumber charge_number) const;

  bool VerifyNuclei(MassNumber mass_number, ChargeNumber charge_number) const;

  static const MassNumber MaxMassNumber;

  static const ChargeNumber MaxChargeNumber;

  static const size_t ParticleCount;

  static const std::vector<FermiFloat> AtomicMassExcess;

  static const std::vector<size_t> IndexCharge;

  static const std::vector<size_t> IndexMass;

  static const size_t ShortTableCount;

  static const std::vector<size_t> ShortTable;
};

std::ostream& operator<<(std::ostream& out, const NucleiPropertiesTable& table);

#endif //FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_TABLEVALUES_NUCLEIPROPERTIESTABLE_H_
