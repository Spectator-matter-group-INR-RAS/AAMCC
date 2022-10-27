
#include "VCollisionReader.hh"
#include "AAMCConstants.hh"

#include "../TGlauber/TGlauNucleon.hh"
#include "TObjArray.h"
class GlauberCollisionReader : public VCollisionReader{

public:
    GlauberCollisionReader() = default;
    ~GlauberCollisionReader() = default;
    void Read(TObjArray* nucleons_in);
    AAMCCinput operator()() final;
    inline AAMCCinput GetAfInput(TObjArray* nucleons_in){this->Read(nucleons_in); return (*this)();};

private:
    AAMCCinput data;
    TObjArray* nucleons;
};
