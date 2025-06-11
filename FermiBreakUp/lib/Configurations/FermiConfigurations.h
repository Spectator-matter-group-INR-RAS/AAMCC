//
// Created by Artem Novikov on 19.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_FERMICONFIGURATIONS_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_FERMICONFIGURATIONS_H_

#include <vector>

#include "FermiSplit.h"
#include "VFermiConfigurations.h"

class FermiConfigurations : public VFermiConfigurations {
 public:
  FermiConfigurations() = default;

  FermiConfigurations(NucleiData nuclei_data, FermiFloat total_energy);

  VFermiConfigurations& GenerateSplits(NucleiData nuclei_data, FermiFloat total_energy) override;

  std::optional<FragmentVector> ChooseSplit() override;

  [[nodiscard]] std::unique_ptr<VFermiConfigurations> Clone() const override;

 private:
  std::vector<FragmentVector> configurations_;
  std::vector<FermiFloat> weights_;
};

#endif //FERMIBREAKUP_MYFERMIBREAKUP_FERMICONFIGURATIONS_H_
