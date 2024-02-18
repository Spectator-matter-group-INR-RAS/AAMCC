//
// Created by Artem Novikov on 17.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTS_UNSTABLEFERMIFRAGMENT_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTS_UNSTABLEFERMIFRAGMENT_H_

#include <vector>

#include "FermiFragment.h"

class UnstableFermiFragment : public FermiFragment {
 public:
  using FermiFragment::FermiFragment;

  UnstableFermiFragment() = delete;

  UnstableFermiFragment(const UnstableFermiFragment&) = delete;

  const UnstableFermiFragment& operator=(const UnstableFermiFragment&) = delete;

  ParticleVector GetFragments(const LorentzVector& momentum) const override;

 protected:
  void Build(const std::vector<NucleiData>& decay_data);

  void Build(std::vector<NucleiData>&& decay_data);

  void FillMasses();

 private:
  std::vector<NucleiData> decay_data_;
  std::vector<FermiFloat> masses_;
};

#endif //FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTS_UNSTABLEFERMIFRAGMENT_H_
