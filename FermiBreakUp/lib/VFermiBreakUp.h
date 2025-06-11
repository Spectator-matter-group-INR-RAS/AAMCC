//
// Created by Artem Novikov on 21.05.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_VFERMIBREAKUP_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_VFERMIBREAKUP_H_

#include "FermiParticle.h"
#include "FermiSplit.h"

class VFermiBreakUp {
 public:
  virtual ParticleVector BreakItUp(const FermiParticle& nucleus) = 0;

  virtual ~VFermiBreakUp() = default;

 protected:
  static ParticleVector ConvertToParticles(const FermiParticle& source_nucleus, const FragmentVector& split);
};

#endif //FERMIBREAKUP_MYFERMIBREAKUP_VFERMIBREAKUP_H_
