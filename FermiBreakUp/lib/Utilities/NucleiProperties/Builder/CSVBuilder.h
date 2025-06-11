//
// Created by Artem Novikov on 21.05.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_NUCLEIPROPERTIES_CSVBUILDER_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_NUCLEIPROPERTIES_CSVBUILDER_H_

#include "VFermiPropertiesBuilder.h"

namespace properties {

class CSVBuilder : public VFermiPropertiesBuilder {
 public:
  CSVBuilder(const std::string& csv_filename, const std::string& mass_number_name = "A",
             const std::string& charge_number_name = "Z", const std::string& mass_name = "mass");

  void BuildTable(FermiNucleiProperties::MassMap& data) const override;

 private:
  std::string csv_filename_;
  std::string mass_number_name_;
  std::string charge_number_name_;
  std::string mass_name_;
};

} // namespace properties

#endif //FERMIBREAKUP_MYFERMIBREAKUP_UTILITIES_NUCLEIPROPERTIES_CSVBUILDER_H_
