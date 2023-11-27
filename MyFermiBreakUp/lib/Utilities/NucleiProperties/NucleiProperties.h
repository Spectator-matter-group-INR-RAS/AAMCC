//
// Created by Artem Novikov on 18.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_NUCLEIPROPERTIES_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_NUCLEIPROPERTIES_H_

#include "Utilities/DataTypes.h"
#include "VNucleiProperties.h"
#include "FermiNucleiProperties.h"
//#include "NativeNucleiProperties.h"

namespace properties {

using NucleiProperties = FermiNucleiProperties;

static_assert(std::is_base_of<VNucleiProperties, NucleiProperties>::value, "Incorrect Nuclei properties class");

} // namespace properties

#endif //FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_NUCLEIPROPERTIES_H_
