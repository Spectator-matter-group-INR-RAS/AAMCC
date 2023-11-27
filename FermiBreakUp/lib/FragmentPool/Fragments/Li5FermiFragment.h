//
// Created by Artem Novikov on 18.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTS_LI5FERMIFRAGMENT_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTS_LI5FERMIFRAGMENT_H_

#include "UnstableFermiFragment.h"

class Li5FermiFragment : public UnstableFermiFragment {
 public:
  Li5FermiFragment(MassNumber mass_number, ChargeNumber charge_number, int polarization, FermiFloat excitation_energy);

  Li5FermiFragment() = delete;

  Li5FermiFragment(const Li5FermiFragment&) = delete;

  Li5FermiFragment& operator=(const Li5FermiFragment&) = delete;
};

#endif //FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTS_LI5FERMIFRAGMENT_H_
