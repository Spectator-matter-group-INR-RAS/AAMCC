//
// Created by Artem Novikov on 21.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_VDECAY_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_VDECAY_H_

#include <vector>
#include "Utilities/DataTypes.h"

class VDecay {
 public:
  VDecay() = default;

  virtual std::vector<LorentzVector> CalculateDecay(const LorentzVector& momentum,
                                            const std::vector<FermiFloat>& fragments_mass) const = 0;

  virtual ~VDecay() = 0;

 protected:
  static FermiFloat TwoBodyMomentum(FermiFloat total_energy, FermiFloat mass1, FermiFloat mass2);

  static std::pair<LorentzVector, LorentzVector> TwoBodyDecay(
      FermiFloat total_energy, FermiFloat mass1, FermiFloat mass2);
};

#endif //FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_VDECAY_H_
