//
// Created by Artem Novikov on 18.02.2023.
//

#include "FermiFragmentPool.h"

#include "Fragments/StableFermiFragment.h"
#include "Fragments/B9FermiFragment.h"
#include "Fragments/Be8FermiFragment.h"
#include "Fragments/He5FermiFragment.h"
#include "Fragments/Li5FermiFragment.h"

#include "Builder/DefaultBuilder.h"

using namespace CLHEP;

namespace pool {

std::unique_ptr<FermiFragmentPool::Map> FermiFragmentPool::fragments_pool_ = nullptr;

FermiFragmentPool::FermiFragmentPool() {
  if (fragments_pool_ == nullptr) {
    Build(DefaultBuilder());
  }
}

FermiFragmentPool::FermiFragmentPool(const VPoolBuilder& builder) {
  Build(builder);
}

size_t FermiFragmentPool::Count(MassNumber mass_number, ChargeNumber charge_number) const {
  return Count({mass_number, charge_number});
}

size_t FermiFragmentPool::Count(NucleiData nuclei) const {
  return fragments_pool_->count(nuclei);
}

FermiFragmentPool::RangeIterators FermiFragmentPool::GetFragments(MassNumber mass_number,
                                                                  ChargeNumber charge_number) const {
  return GetFragments({mass_number, charge_number});
}

FermiFragmentPool::RangeIterators FermiFragmentPool::GetFragments(NucleiData nuclei) const {
  return fragments_pool_->equal_range(nuclei);
}

void FermiFragmentPool::Build(const VPoolBuilder& builder) {
  fragments_pool_ = std::make_unique<Map>();
  builder.Build(*fragments_pool_);
}

} // namespace pool
