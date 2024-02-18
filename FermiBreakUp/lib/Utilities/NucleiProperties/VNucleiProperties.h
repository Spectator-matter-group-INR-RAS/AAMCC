//
// Created by Artem Novikov on 09.03.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_VNUCLEIPROPERTIES_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_VNUCLEIPROPERTIES_H_

#include "Utilities/DataTypes.h"

namespace properties {

class VNucleiProperties {
 public:
  virtual FermiFloat GetNuclearMass(MassNumber mass_number, ChargeNumber charge_number) const = 0;

  virtual bool IsStable(MassNumber mass_number, ChargeNumber charge_number) const = 0;

  virtual ~VNucleiProperties() = 0;

 protected:
  static FermiFloat AtomicMass(MassNumber mass_number, ChargeNumber charge_number);

  static FermiFloat NuclearMass(MassNumber mass_number, ChargeNumber charge_number);

  static FermiFloat BindingEnergy(MassNumber mass_number, ChargeNumber charge_number);

  static void PrintInvalidNuclei(MassNumber mass_number, ChargeNumber charge_number);

  static bool IsInvalidNuclei(MassNumber mass_number, ChargeNumber charge_number);
};

} // namespace properties

#endif //FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_VNUCLEIPROPERTIES_H_
