//
// Created by Artem Novikov on 21.05.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_VFERMICONFIGURATIONS_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_VFERMICONFIGURATIONS_H_

#include <optional>
#include <memory>

#include "Utilities/DataTypes.h"
#include "FermiSplit.h"

class VFermiConfigurations {
 public:
  virtual VFermiConfigurations& GenerateSplits(NucleiData nuclei_data, FermiFloat total_energy) = 0;

  virtual std::optional<FragmentVector> ChooseSplit() = 0;

  virtual std::unique_ptr<VFermiConfigurations> Clone() const = 0;

  virtual ~VFermiConfigurations() = default;
};

#endif //FERMIBREAKUP_MYFERMIBREAKUP_VFERMICONFIGURATIONS_H_
