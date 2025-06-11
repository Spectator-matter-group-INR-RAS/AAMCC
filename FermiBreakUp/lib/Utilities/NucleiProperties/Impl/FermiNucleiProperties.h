//
// Created by Artem Novikov on 09.03.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_FERMINUCLEIPROPERTIES_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_FERMINUCLEIPROPERTIES_H_

#include <map>
#include <memory>

#include "../VNucleiProperties.h"

namespace properties {

class FermiNucleiProperties : public VNucleiProperties {
 public:
  FermiNucleiProperties();

  template <typename DataSource>
  FermiNucleiProperties(const DataSource& data_source);

  template <typename Iter>
  FermiNucleiProperties(Iter begin, Iter end);

  [[nodiscard]] FermiFloat GetNuclearMass(MassNumber mass_number, ChargeNumber charge_number) const override;

  [[nodiscard]] bool IsStable(MassNumber mass_number, ChargeNumber charge_number) const override;

  void AddMass(MassNumber mass_number, ChargeNumber charge_number, FermiFloat mass);

  void AddMass(NucleiData nuclei_data, FermiFloat mass);

 private:
  std::map<NucleiData, FermiFloat> nuclei_mass_;
};

template <typename DataSource>
FermiNucleiProperties::FermiNucleiProperties(const DataSource& data_source)
    : FermiNucleiProperties(data_source.begin(), data_source.end()) {}

template <typename Iter>
FermiNucleiProperties::FermiNucleiProperties(Iter begin, Iter end) {
  static_assert(std::is_same_v<typename Iter::value_type, std::pair<const NucleiData, FermiFloat>>, "invalid iterator");
  for (auto it = begin; it != end; ++it) {
    AddMass(it->first, it->second);
  }
}

} // namespace properties

#endif //FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_FERMINUCLEIPROPERTIES_H_
