//
// Created by Artem Novikov on 21.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_DECAY_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_DECAY_H_

#include "VDecay.h"

class Decay : public VDecay {
 public:
  std::vector<LorentzVector> CalculateDecay(const LorentzVector& momentum,
                                            const std::vector<FermiFloat>& fragments_mass) const override;

 private:
  static std::vector<FermiFloat> CalculateVirtualMasses(const std::vector<FermiFloat>& masses, FermiFloat energy);

  /// modifies vector for a optimization purposes
  static FermiFloat CalculateMomentumMagnitudes(std::vector<FermiFloat>& daughter_momentum,
                                        const std::vector<FermiFloat>& masses,
                                        const std::vector<FermiFloat>& virtual_masses);

  static FermiFloat CalculateMaxWeight(const std::vector<FermiFloat>& masses, FermiFloat max_energy);
};

#endif //FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_DECAY_H_
