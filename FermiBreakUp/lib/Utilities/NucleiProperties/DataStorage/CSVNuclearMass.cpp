//
// Created by Artem Novikov on 18.02.2024.
//

#include <fstream>

#include "CSVNuclearMass.h"

namespace properties {

CSVNuclearMass::CSVNuclearMass(const std::string& csv_filename,
                               const std::string& mass_number_name,
                               const std::string& charge_number_name,
                               const std::string& mass_name) {
  std::ifstream in(csv_filename);

  if (!in.is_open()) {
    throw std::runtime_error("invalid CSV file path");
  }

  int mass_number_idx = -1, charge_number_idx = -1, mass_idx = -1;

  std::string line;
  in >> line;
  line += ',';
  int column_idx = 0;
  uint start = 0;
  auto comma = line.find(',', start);
  while (comma != std::string::npos) {
    if (line.substr(start, comma - start) == mass_number_name) {
      mass_number_idx = column_idx;
    }

    if (line.substr(start, comma - start) == charge_number_name) {
      charge_number_idx = column_idx;
    }

    if (line.substr(start, comma - start) == mass_name) {
      mass_idx = column_idx;
    }

    ++column_idx;
    start = comma + 1;
    comma = line.find(',', start);
  }
  int columns_count = column_idx;

  if (mass_number_idx == -1) {
    throw std::runtime_error("no mass number found");
  }

  if (charge_number_idx == -1) {
    throw std::runtime_error("no charge number found");
  }

  if (mass_idx == -1) {
    throw std::runtime_error("no nuclei mass found");
  }

  MassNumber m = 0_m;
  ChargeNumber c = 0_c;
  FermiFloat mass;
  while (in >> line) {
    line += ',';
    column_idx = 0;
    start = 0;
    comma = line.find(',', start);
    while (comma != std::string::npos) {
      if (column_idx == mass_number_idx) {
        m = MassNumber(std::stoi(line.substr(start, comma - start)));
      }

      if (column_idx == charge_number_idx) {
        c = ChargeNumber(std::stoi(line.substr(start, comma - start)));
      }

      if (column_idx == mass_idx) {
        mass = FermiFloat(std::stod(line.substr(start, comma - start)));
      }

      ++column_idx;
      start = comma + 1;
      comma = line.find(',', start);
    }

    if (columns_count != column_idx) {
      throw std::runtime_error("invalid row format: " + line);
    }

    masses_.emplace(NucleiData{m, c}, mass);
  }
}

} // namespace properties
