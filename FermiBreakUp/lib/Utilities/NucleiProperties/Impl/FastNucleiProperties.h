//
// Created by Artem Novikov on 18.02.2024.
//

#ifndef HANDLER_MYFERMIBREAKUP_UTILITIES_NUCLEIPROPERTIES_IMPL_FERMIFASTNUCLEIPROPERTIES_H_
#define HANDLER_MYFERMIBREAKUP_UTILITIES_NUCLEIPROPERTIES_IMPL_FERMIFASTNUCLEIPROPERTIES_H_

#include <vector>

#include "Utilities/DataTypes.h"
#include "../VNucleiProperties.h"

namespace properties {

class FastNucleiProperties : public VNucleiProperties {
 public:
  FastNucleiProperties();

  template <typename DataSource>
  FastNucleiProperties(const DataSource& data_source);

  template <typename Iter>
  FastNucleiProperties(Iter begin, Iter end);

  [[nodiscard]] FermiFloat GetNuclearMass(MassNumber mass_number, ChargeNumber charge_number) const override;

  [[nodiscard]] FermiFloat GetNuclearMass(MassNumber mass_number, ChargeNumber charge_number);

  [[nodiscard]] bool IsStable(MassNumber mass_number, ChargeNumber charge_number) const override;

  void AddMass(MassNumber mass_number, ChargeNumber charge_number, FermiFloat mass);

  void AddMass(NucleiData nuclei_data, FermiFloat mass);

  ~FastNucleiProperties() = default;

 private:
  struct MassData {
    FermiFloat mass;
    bool is_valid = false;
    bool is_cached = false;
  };

  [[nodiscard]] static size_t GetSlot(MassNumber mass_number, ChargeNumber charge_number);

  std::vector<MassData> nuclei_masses_;
};

template <typename DataSource>
FastNucleiProperties::FastNucleiProperties(const DataSource& data_source)
    : FastNucleiProperties(data_source.begin(), data_source.end()) {}

template <typename Iter>
FastNucleiProperties::FastNucleiProperties(Iter begin, Iter end) {
  static_assert(std::is_same_v<typename Iter::value_type, std::pair<const NucleiData, FermiFloat>>, "invalid iterator");
  for (auto it = begin; it != end; ++it) {
    AddMass(it->first, it->second);
  }
}

} // namespace properties

#endif //HANDLER_MYFERMIBREAKUP_UTILITIES_NUCLEIPROPERTIES_IMPL_FERMIFASTNUCLEIPROPERTIES_H_
