//
// Created by Artem Novikov on 09.03.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_FERMINUCLEIPROPERTIES_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_FERMINUCLEIPROPERTIES_H_

#include <map>
#include <memory>

#include "VNucleiProperties.h"

namespace properties {

class VFermiPropertiesBuilder;

class FermiNucleiProperties : public VNucleiProperties {
 public:
  using MassMap = std::map<NucleiData, FermiFloat>;

  FermiNucleiProperties();

  FermiNucleiProperties(const VFermiPropertiesBuilder& builder);

  [[nodiscard]] FermiFloat GetNuclearMass(MassNumber mass_number, ChargeNumber charge_number) const override;

  [[nodiscard]] bool IsStable(MassNumber mass_number, ChargeNumber charge_number) const override;

 private:
  static void Build(const VFermiPropertiesBuilder& builder);

  static std::unique_ptr<MassMap> nuclei_mass_;
};

} // namespace properties

#endif //FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_FERMINUCLEIPROPERTIES_H_
