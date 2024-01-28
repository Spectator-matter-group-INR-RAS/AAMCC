//
// Created by Artem Novikov on 18.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTS_BE8FERMIFRAGMENT_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTS_BE8FERMIFRAGMENT_H_

#include "UnstableFermiFragment.h"

class Be8FermiFragment : public UnstableFermiFragment {
 public:
  Be8FermiFragment(MassNumber mass_number, ChargeNumber charge_number, int polarization, FermiFloat excitation_energy);

  Be8FermiFragment() = delete;

  Be8FermiFragment(const Be8FermiFragment&) = delete;

  Be8FermiFragment& operator=(const Be8FermiFragment&) = delete;
};

#endif //FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTS_BE8FERMIFRAGMENT_H_
