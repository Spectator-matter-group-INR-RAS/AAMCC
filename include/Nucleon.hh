#include <iostream>
#include <vector>
#if !defined(NUCLEON_HH)
#define NUCLEON_HH

namespace aamcc {

    class Nucleon {
    public:
        Nucleon();

        ~Nucleon() = default;

        double x;
        double y;
        double z;
        bool isospin; // 1 - proton
        bool isParticipant; // 1 - participant
        std::string Nucl;

        void Clean();

        double GetX() const { return x; }

        double GetY() const { return y; }

        double GetZ() const { return z; }
    };

    class NucleonVector : public std::vector<Nucleon> {
    public:
        int GetA(std::string Nucl);

        int GetZ(std::string Nucl);

        int GetTotA(std::string Nucl);

        int GetTotZ(std::string Nucl);

        NucleonVector *GetNucleons(std::string side);
    };
}   //aamcc
#endif