#ifndef VCollisionReader_h
#define VCollisionReader_h 1
#include "VCollisionReader.hh"
#endif


#include "../TGlauber/TGlauNucleon.hh"
#include "TObjArray.h"
class GlauberCollisionReader : public VCollisionReader{

public:
    GlauberCollisionReader();
    ~GlauberCollisionReader() = default;
    void Read(TObjArray* nucleons_in);
    NucleonVector GetNucleons();

private:
    NucleonVector nucleonVector;
    TObjArray* nucleons;
};