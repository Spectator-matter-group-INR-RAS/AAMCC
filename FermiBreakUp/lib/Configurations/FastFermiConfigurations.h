//
// Created by Artem Novikov on 18.02.2024.
//

#ifndef HANDLER_MYFERMIBREAKUP_FASTFERMICONFIGURATION_H_
#define HANDLER_MYFERMIBREAKUP_FASTFERMICONFIGURATION_H_

#include <unordered_map>

#include "FermiSplit.h"
#include "VFermiConfigurations.h"

class FastFermiConfigurations : public VFermiConfigurations {
 public:
  FastFermiConfigurations() = default;

  FastFermiConfigurations(NucleiData nuclei_data, FermiFloat total_energy);

  VFermiConfigurations& GenerateSplits(NucleiData nuclei_data, FermiFloat total_energy) override;

  std::optional<FragmentVector> ChooseSplit() override;

  [[nodiscard]] std::unique_ptr<VFermiConfigurations> Clone() const override;

 private:
  std::vector<size_t> configurations_;
  std::vector<FermiFloat> weights_;

  std::unordered_map<NucleiData, std::vector<FragmentVector>> cached_configurations_;
  NucleiData last_nuclei_{};
};

#endif //HANDLER_MYFERMIBREAKUP_FASTFERMICONFIGURATION_H_
