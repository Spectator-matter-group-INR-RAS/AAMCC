//
// Created by Artem Novikov on 21.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_KOPYLOVDECAY_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_KOPYLOVDECAY_H_

#include "VDecay.h"

class KopylovDecay : public VDecay{
 public:
  std::vector<LorentzVector> CalculateDecay(const LorentzVector& momentum,
                                            const std::vector<FermiFloat>& fragments_mass) const override;

 private:
  static FermiFloat BetaKopylov(size_t k);
};

#endif //FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_KOPYLOVDECAY_H_
