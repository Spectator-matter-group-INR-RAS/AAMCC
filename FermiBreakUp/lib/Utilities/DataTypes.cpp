//
// Created by Artem Novikov on 07.03.2023.
//

#include "DataTypes.h"

ChargeNumber operator ""_c(unsigned long long charge) {
  return ChargeNumber(charge);
}

namespace std {
std::string to_string(ChargeNumber charge) {
  return std::to_string(FermiUInt(charge));
}
}

std::ostream& operator<<(std::ostream& out, ChargeNumber charge) {
  out << FermiUInt(charge);
  return out;
}

std::istream& operator>>(std::istream& in, ChargeNumber& charge) {
  FermiUInt val;
  in >> val;
  charge = ChargeNumber(val);
  return in;
}

MassNumber operator ""_m(unsigned long long mass) {
  return MassNumber(mass);
}

namespace std {
std::string to_string(MassNumber mass) {
  return std::to_string(FermiUInt(mass));
}
}

std::ostream& operator<<(std::ostream& out, const MassNumber& mass) {
  out << FermiUInt(mass);
  return out;
}

std::istream& operator>>(std::istream& in, MassNumber& mass) {
  FermiUInt val;
  in >> val;
  mass = MassNumber(val);
  return in;
}