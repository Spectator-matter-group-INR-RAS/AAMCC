//
// Created by Artem Novikov on 07.02.2023.
//

#ifndef FERMIBREAKUP_MYFERMIBREAKUP_DATATYPES_H_
#define FERMIBREAKUP_MYFERMIBREAKUP_DATATYPES_H_

#include <CLHEP/Vector/LorentzVector.h>

using FermiInt = int;

using FermiUInt = uint32_t;

using FermiFloat = double;

using LorentzVector = CLHEP::HepLorentzVector;

using Vector3 = CLHEP::Hep3Vector;

using ParticleMomentum = Vector3;

class MassNumber {
 public:
  using ValueType = FermiUInt;

  MassNumber() = default;

  explicit MassNumber(ValueType mass) : mass_(mass) {}

  MassNumber(const MassNumber& other) = default;

  MassNumber(MassNumber&& other) = default;

  MassNumber& operator=(const MassNumber& other) = default;

  MassNumber& operator=(MassNumber&& other) = default;

  operator FermiUInt() const { return mass_; }

  operator FermiInt() const { return mass_; }

  operator FermiFloat() const { return mass_; }

  bool operator<(const MassNumber& other) const { return mass_ < other.mass_; }

  bool operator>(const MassNumber& other) const { return mass_ > other.mass_; }

  bool operator<=(const MassNumber& other) const { return mass_ <= other.mass_; }

  bool operator>=(const MassNumber& other) const { return mass_ >= other.mass_; }

  bool operator==(const MassNumber& other) const { return mass_ == other.mass_; }

  bool operator!=(const MassNumber& other) const { return mass_ != other.mass_; }

 private:
  uint32_t mass_;
};

MassNumber operator ""_m(unsigned long long mass);

namespace std {
std::string to_string(MassNumber mass);
}

std::ostream& operator<<(std::ostream& out, const MassNumber& mass);

std::istream& operator>>(std::istream& in, MassNumber& mass);

class ChargeNumber {
 public:
  using ValueType = FermiUInt;

  ChargeNumber() = default;

  explicit ChargeNumber(ValueType charge) : charge_(charge) {}

  ChargeNumber(const ChargeNumber& other) = default;

  ChargeNumber(ChargeNumber&& other) = default;

  ChargeNumber& operator=(const ChargeNumber& other) = default;

  ChargeNumber& operator=(ChargeNumber&& other) = default;

  operator FermiUInt() const { return charge_; }

  operator FermiInt() const { return charge_; }

  operator FermiFloat() const { return charge_; }

  bool operator<(const ChargeNumber& other) const { return charge_ < other.charge_; }

  bool operator>(const ChargeNumber& other) const { return charge_ > other.charge_; }

  bool operator<=(const ChargeNumber& other) const { return charge_ <= other.charge_; }

  bool operator>=(const ChargeNumber& other) const { return charge_ >= other.charge_; }

  bool operator==(const ChargeNumber& other) const { return charge_ == other.charge_; }

  bool operator!=(const ChargeNumber& other) const { return charge_ != other.charge_; }

 private:
  uint32_t charge_;
};

ChargeNumber operator ""_c(unsigned long long charge);

namespace std {
std::string to_string(ChargeNumber charge);
}

std::ostream& operator<<(std::ostream& out, ChargeNumber charge);

std::istream& operator>>(std::istream& in, ChargeNumber& charge);

struct NucleiData {
  MassNumber mass_number;
  ChargeNumber charge_number;

  bool operator<(const NucleiData& other) const {
    return mass_number < other.mass_number
        || (mass_number == other.mass_number && charge_number < other.charge_number);
  }

  bool operator==(const NucleiData& other) const {
    return mass_number == other.mass_number && charge_number == other.charge_number;
  }

  bool operator!=(const NucleiData& other) const {
    return mass_number != other.mass_number || charge_number != other.charge_number;
  }
};

#endif //FERMIBREAKUP_MYFERMIBREAKUP_DATATYPES_H_
