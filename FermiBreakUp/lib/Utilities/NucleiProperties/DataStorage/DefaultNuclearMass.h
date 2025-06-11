//
// Created by Artem Novikov on 18.02.2024.
//

#ifndef HANDLER_MYFERMIBREAKUP_UTILITIES_NUCLEIPROPERTIES_DATASTORAGE_DEFAULTNUCLEARMASS_H_
#define HANDLER_MYFERMIBREAKUP_UTILITIES_NUCLEIPROPERTIES_DATASTORAGE_DEFAULTNUCLEARMASS_H_

#include <vector>

#include "Utilities/DataTypes.h"

namespace properties {

class DefaultNuclearMass {
 private:
  using Container = std::vector<std::pair<const NucleiData, FermiFloat>>;

 public:
  DefaultNuclearMass();

  ~DefaultNuclearMass() = default;

  Container::iterator begin() { return masses_.begin(); }

  [[nodiscard]] Container::const_iterator begin() const { return masses_.begin(); }

  [[nodiscard]] Container::const_iterator cbegin() const { return masses_.cbegin(); }

  Container::iterator end() { return masses_.end(); }

  [[nodiscard]] Container::const_iterator end() const { return masses_.end(); }

  Container::const_iterator cend() const { return masses_.cend(); }

 private:
  Container masses_;
};

} // properties

#endif //HANDLER_MYFERMIBREAKUP_UTILITIES_NUCLEIPROPERTIES_DATASTORAGE_DEFAULTNUCLEARMASS_H_
