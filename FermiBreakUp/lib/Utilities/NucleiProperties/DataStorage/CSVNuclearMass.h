//
// Created by Artem Novikov on 18.02.2024.
//

#ifndef HANDLER_MYFERMIBREAKUP_UTILITIES_NUCLEIPROPERTIES_DATASTORAGE_CSVNUCLEARMASS_H_
#define HANDLER_MYFERMIBREAKUP_UTILITIES_NUCLEIPROPERTIES_DATASTORAGE_CSVNUCLEARMASS_H_

#include <map>

#include "Utilities/DataTypes.h"

namespace properties {

class CSVNuclearMass {
 private:
  using Container = std::map<NucleiData, FermiFloat>;

 public:
  CSVNuclearMass(const std::string& csv_filename,
                 const std::string& mass_number_name = "A",
                 const std::string& charge_number_name = "Z",
                 const std::string& mass_name = "mass");

  ~CSVNuclearMass() = default;

  Container::iterator begin() { return masses_.begin(); }

  [[nodiscard]] Container::const_iterator begin() const { return masses_.begin(); }

  [[nodiscard]] Container::const_iterator cbegin() const { return masses_.cbegin(); }

  Container::iterator end() { return masses_.end(); }

  [[nodiscard]] Container::const_iterator end() const { return masses_.end(); }

  [[nodiscard]] Container::const_iterator cend() const { return masses_.cend(); }

 private:
  Container masses_;
};

} // namespace properties

#endif //HANDLER_MYFERMIBREAKUP_UTILITIES_NUCLEIPROPERTIES_DATASTORAGE_CSVNUCLEARMASS_H_
