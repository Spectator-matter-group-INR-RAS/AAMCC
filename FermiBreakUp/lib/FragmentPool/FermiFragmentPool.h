//
// Created by Artem Novikov on 18.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_FERMIFRAGMENTPOOL_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_FERMIFRAGMENTPOOL_H_

#include <memory>
#include <map>

#include "Utilities/DataTypes.h"
#include "Utilities/Singleton.h"
#include "Fragments/FermiFragment.h"
#include "FragmentsStorage.h"

namespace pool {

using FermiFragmentPool = Singleton<FragmentsStorage>;

} // namespace pool

#endif //FERMIBREAKUP_MYFERMIBREAKUP_FERMIFRAGMENTPOOL_H_
