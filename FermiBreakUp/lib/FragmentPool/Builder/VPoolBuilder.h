//
// Created by Artem Novikov on 25.05.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTPOOL_VPOOLBUILDER_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTPOOL_VPOOLBUILDER_H_

#include "../FermiFragmentPool.h"

namespace pool {

class VPoolBuilder {
 public:
  virtual void Build(FermiFragmentPool::Map& map) const = 0;

  virtual ~VPoolBuilder() = default;
};

} // namespace pool

#endif //FERMIBREAKUP_MYFERMIBREAKUP_FRAGMENTPOOL_VPOOLBUILDER_H_
