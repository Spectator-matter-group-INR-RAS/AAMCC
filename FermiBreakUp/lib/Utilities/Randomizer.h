//
// Created by Artem Novikov on 19.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_RANDOMIZER_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_RANDOMIZER_H_

#include <vector>
#include <algorithm>

#include "DataTypes.h"

class Randomizer {
 public:
  static FermiFloat UniformRealDistribution();

  static FermiFloat NormalDistribution(FermiFloat mean = 0, FermiFloat deviation = 1);

  static ParticleMomentum IsotropicVector(FermiFloat magnitude = 1);

  static std::vector<FermiFloat> ProbabilityDistribution(size_t point_count);
};

#endif //FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_RANDOMIZER_H_
