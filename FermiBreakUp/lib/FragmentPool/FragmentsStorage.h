//
// Created by Artem Novikov on 30.01.2024.
//

#ifndef HANDLER_MYFERMIBREAKUP_FRAGMENTPOOL_STORAGE_FRAGMENTSSTORAGE_H_
#define HANDLER_MYFERMIBREAKUP_FRAGMENTPOOL_STORAGE_FRAGMENTSSTORAGE_H_

#include <vector>

#include "Utilities/DataTypes.h"
#include "Fragments/FermiFragment.h"

namespace pool {

class FragmentsStorage {
 private:
  using Container = std::vector<const FermiFragment*>;

 public:
  using RangeIterators = std::pair<Container::const_iterator, Container::const_iterator>;

  FragmentsStorage();

  template <typename DataSource>
  FragmentsStorage(const DataSource& data_source);

  template <typename Iter>
  FragmentsStorage(Iter begin, Iter end);

  [[nodiscard]] size_t Count(MassNumber mass_number, ChargeNumber charge_number) const;

  [[nodiscard]] size_t Count(NucleiData nuclei) const;

  [[nodiscard]] RangeIterators GetFragments(MassNumber mass_number, ChargeNumber charge_number) const;

  [[nodiscard]] RangeIterators GetFragments(NucleiData nuclei) const;

  void AddFragment(const FermiFragment& fragment);

 private:
  [[nodiscard]] static size_t GetSlot(MassNumber mass_number, ChargeNumber charge_number);

  std::vector<Container> fragments_;
};

} // namespace pool

#endif //HANDLER_MYFERMIBREAKUP_FRAGMENTPOOL_STORAGE_FRAGMENTSSTORAGE_H_
