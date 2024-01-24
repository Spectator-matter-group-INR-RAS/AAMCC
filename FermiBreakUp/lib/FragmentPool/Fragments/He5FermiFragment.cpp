//
// Created by Artem Novikov on 18.02.2023.
//

#include "He5FermiFragment.h"

He5FermiFragment::He5FermiFragment(MassNumber mass_number,
                                   ChargeNumber charge_number,
                                   int polarization,
                                   FermiFloat excitation_energy)
    : UnstableFermiFragment(mass_number, charge_number, polarization, excitation_energy) {
  // He5 ----> alpha + neutron

  Build({NucleiData{4_m, 2_c},
         NucleiData{1_m, 0_c}});
}
