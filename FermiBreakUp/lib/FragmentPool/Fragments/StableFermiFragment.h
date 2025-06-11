//
// Created by Artem Novikov on 17.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTS_STABLEFERMIFRAGMENT_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTS_STABLEFERMIFRAGMENT_H_

#include "FermiFragment.h"

class StableFermiFragment : public FermiFragment {
 public:
  using FermiFragment::FermiFragment;

  StableFermiFragment() = delete;

  StableFermiFragment(const StableFermiFragment&) = delete;

  StableFermiFragment& operator=(const StableFermiFragment&) = delete;

  ParticleVector GetFragment(const LorentzVector& momentum) const override;
};

#endif //FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTS_STABLEFERMIFRAGMENT_H_
