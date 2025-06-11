//
// Created by Artem Novikov on 09.03.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_NATIVENUCLEIPROPERTIES_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_NATIVENUCLEIPROPERTIES_H_

#include "../VNucleiProperties.h"

namespace properties {

class NativeNucleiProperties : public VNucleiProperties {
 public:
  [[nodiscard]] FermiFloat GetNuclearMass(MassNumber mass_number, ChargeNumber charge_number) const override;

  [[nodiscard]] bool IsStable(MassNumber mass_number, ChargeNumber charge_number) const override;

 private:
  static FermiFloat AtomicMass(MassNumber mass_number, ChargeNumber charge_number);

  static FermiFloat NuclearMass(MassNumber mass_number, ChargeNumber charge_number);

  static FermiFloat BindingEnergy(MassNumber mass_number, ChargeNumber charge_number);

  static bool IsInvalidNuclei(MassNumber mass_number, ChargeNumber charge_number);

  static const FermiFloat proton_mass;
  static const FermiFloat neutron_mass;
  static const FermiFloat deuteron_mass;
  static const FermiFloat triton_mass;
  static const FermiFloat alpha_mass;
  static const FermiFloat He3_mass;
};

} // namespace properties

#endif //FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_NATIVENUCLEIPROPERTIES_H_
