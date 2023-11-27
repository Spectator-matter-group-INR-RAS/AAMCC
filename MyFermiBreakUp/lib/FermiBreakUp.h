//
// Created by Artem Novikov on 21.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_FERMIBREAKUP_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_FERMIBREAKUP_H_

#include <memory>

#include "FermiParticle.h"
#include "FermiConfigurations.h"
#include "VFermiBreakUp.h"
#include "VFermiConfigurations.h"

class FermiBreakUp : public VFermiBreakUp {
 public:
  FermiBreakUp();

  FermiBreakUp(std::unique_ptr<VFermiConfigurations>&& configurations);

  ParticleVector BreakItUp(const FermiParticle& nucleus) override;

  ~FermiBreakUp() = default;

 private:
  static std::unique_ptr<VFermiConfigurations> DefaultConfigurations();

  std::unique_ptr<VFermiConfigurations> fermi_configurations_;
};

#endif //FERMIBREAKUP_MYFERMIBREAKUP_FERMIBREAKUP_H_
