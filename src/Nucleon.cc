#include "Nucleon.hh"


namespace aamcc {
    void Nucleon::Clean() {
        x = 0;
        y = 0;
        z = 0;
        isospin = false;
        isParticipant = false;
        Nucl = "";
    }

    Nucleon::Nucleon() {
    }

    int NucleonVector::GetZ(std::string Nucl) {
        int tZ = 0;
        for (long unsigned int k = 0; k < (*this).size(); ++k) {
            if (!at(k).isParticipant && at(k).Nucl == Nucl && at(k).isospin) { tZ++; }
        }
        return tZ;
    }

    int NucleonVector::GetA(std::string Nucl) {
        int tA = 0;
        for (long unsigned int k = 0; k < (*this).size(); ++k) {
            if (!at(k).isParticipant && at(k).Nucl == Nucl) { tA++; }
        }
        return tA;
    }

    int NucleonVector::GetTotA(std::string Nucl) {
        int tA = 0;
        for (long unsigned int k = 0; k < (*this).size(); ++k) { if (at(k).Nucl == Nucl) { tA++; }}
        return tA;
    }

    int NucleonVector::GetTotZ(std::string Nucl) {
        int tZ = 0;
        for (long unsigned int k = 0; k < (*this).size(); ++k) { if (at(k).Nucl == Nucl && at(k).isospin) { tZ++; }}
        return tZ;
    }

    NucleonVector *NucleonVector::GetNucleons(std::string side) {
        NucleonVector *out = new NucleonVector();
        if (side == "A" || side == "B") {
            for (int k = 0; k < this->size(); ++k) { if (this->at(k).Nucl == side) { out->push_back(this->at(k)); }}
        } else { out = this; }
        return out;
    }
}   //aamcc