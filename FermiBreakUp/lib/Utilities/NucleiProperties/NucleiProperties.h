//
// Created by Artem Novikov on 18.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_NUCLEIPROPERTIES_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_NUCLEIPROPERTIES_H_

#include "Utilities/DataTypes.h"
#include "Utilities/Singleton.h"

#include "VNucleiProperties.h"
#include "Impl/FermiNucleiProperties.h"
#include "Impl/NativeNucleiProperties.h"
#include "Impl/FastNucleiProperties.h"

namespace properties {

/// it is possible to use polymorphism here, but this way no inline table access is possible and it is a bottleneck
using NucleiProperties = Singleton<FastNucleiProperties>;

static_assert(std::is_base_of_v<VNucleiProperties, std::remove_reference_t<decltype(NucleiProperties::Instance())>>,
              "Incorrect Nuclei properties class");

} // namespace properties

#endif //FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_NUCLEIPROPERTIES_H_
