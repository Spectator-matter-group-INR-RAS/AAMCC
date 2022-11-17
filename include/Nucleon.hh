#include <iostream>
#include <vector>
#include <math.h>

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

   inline int IsotopeToPDG(int Z, int A){
        if(Z == 0 && A == 1){return 2112;}
        if(Z == 1 && A == 1){return 2212;}

        int zHund = Z/100;
        int zDoz = Z/10;
        int aHund = A/100;
        int aDoz = A/10;

        int pdgpid = std::pow(10,9)+std::pow(10,6)*zHund+std::pow(10,5)*(zDoz-10*zHund) + std::pow(10,4)*(Z - zDoz*10) + std::pow(10,3)*aHund + 100*(aDoz - aHund*10) + 10*(A - aDoz*10);

        return pdgpid;
    }
}   //aamcc
#endif