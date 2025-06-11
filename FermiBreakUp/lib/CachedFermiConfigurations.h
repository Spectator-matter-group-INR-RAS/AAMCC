//
// Created by Artem Novikov on 24.05.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_CACHEDFERMICONFIGURATIONS_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_CACHEDFERMICONFIGURATIONS_H_

#include "FermiSplit.h"
#include "VFermiConfigurations.h"

class CachedFermiConfigurations : public VFermiConfigurations {
 public:
  CachedFermiConfigurations() = default;

  CachedFermiConfigurations(NucleiData nuclei_data, FermiFloat total_energy);

  VFermiConfigurations& GenerateSplits(NucleiData nuclei_data, FermiFloat total_energy) override;

  std::optional<FragmentVector> ChooseSplit() override;

 private:
  std::vector<size_t> configurations_;
  std::vector<FermiFloat> weights_;

  std::vector<FragmentVector> cached_configurations_;
  NucleiData last_nuclei_{};
};

#endif //FERMIBREAKUP_MYFERMIBREAKUP_CACHEDFERMICONFIGURATIONS_H_
