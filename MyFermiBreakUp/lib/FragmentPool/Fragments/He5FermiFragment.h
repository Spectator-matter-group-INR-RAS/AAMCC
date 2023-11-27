//
// Created by Artem Novikov on 18.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTS_HE5FERMIFRAGMENT_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTS_HE5FERMIFRAGMENT_H_

#include "UnstableFermiFragment.h"

class He5FermiFragment : public UnstableFermiFragment {
 public:
  He5FermiFragment(MassNumber mass_number, ChargeNumber charge_number, int polarization, FermiFloat excitation_energy);

  He5FermiFragment() = delete;

  He5FermiFragment(const He5FermiFragment&) = delete;

  He5FermiFragment& operator=(const He5FermiFragment&) = delete;
};

#endif //FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTS_HE5FERMIFRAGMENT_H_
