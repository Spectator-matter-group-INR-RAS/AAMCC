#include "VCollisionReader.hh"
#include "../TGlauber/TGlauNucleon.hh"
#include "TObjArray.h"
class GlauberCollisionReader : public VCollisionReader{

public:
    GlauberCollisionReader();
    ~GlauberCollisionReader() = default;
    void Read(TObjArray* nucleons_in);
    NucleonVector GetNucleons();
    inline NucleonVector GetNucleons(TObjArray* nucleons_in) {this->Read(nucleons_in); return GetNucleons();};

private:
    NucleonVector nucleonVector;
    TObjArray* nucleons;
};