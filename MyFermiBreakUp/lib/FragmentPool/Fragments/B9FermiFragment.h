//
// Created by Artem Novikov on 18.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTS_B9FERMIFRAGMENT_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTS_B9FERMIFRAGMENT_H_

#include "UnstableFermiFragment.h"

class B9FermiFragment : public UnstableFermiFragment {
 public:
  B9FermiFragment(MassNumber mass_number, ChargeNumber charge_number, int polarization, FermiFloat excitation_energy);

  B9FermiFragment() = delete;

  B9FermiFragment(const B9FermiFragment&) = delete;

  B9FermiFragment& operator=(const B9FermiFragment&) = delete;
};

#endif //FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTS_B9FERMIFRAGMENT_H_
