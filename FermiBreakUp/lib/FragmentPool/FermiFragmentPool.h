//
// Created by Artem Novikov on 18.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_FERMIFRAGMENTPOOL_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_FERMIFRAGMENTPOOL_H_

#include <memory>
#include <map>

#include "Utilities/DataTypes.h"
#include "Fragments/FermiFragment.h"

namespace pool {

class VPoolBuilder;

class FermiFragmentPool {
 private:
  using RangeIterators = std::pair<std::multimap<NucleiData, const FermiFragment*>::const_iterator,
                                   std::multimap<NucleiData, const FermiFragment*>::const_iterator>;
 public:
  using Map = std::multimap<NucleiData, const FermiFragment*>;

  FermiFragmentPool();

  FermiFragmentPool(const VPoolBuilder& builder);

  size_t Count(MassNumber mass_number, ChargeNumber charge_number) const;

  size_t Count(NucleiData nuclei) const;

  RangeIterators GetFragments(MassNumber mass_number, ChargeNumber charge_number) const;

  RangeIterators GetFragments(NucleiData nuclei) const;

 private:
  void Build(const VPoolBuilder& builder);

  static std::unique_ptr<Map> fragments_pool_;
};

} // namespace pool

#endif //FERMIBREAKUP_MYFERMIBREAKUP_FERMIFRAGMENTPOOL_H_
