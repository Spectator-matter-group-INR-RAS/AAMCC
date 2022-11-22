#include "VCollisionReader.hh"

#include "../TGlauber/TGlauNucleon.hh"
#include "TObjArray.h"
class GlauberCollisionReader final : public VCollisionReader{

public:
    GlauberCollisionReader() = default;
    ~GlauberCollisionReader() final  = default;
    // No proper memory management due to legacy code. ROOT library has its own garbage collector.
    // Using delete explicitly or smart pointers results in segmentation violation.
    void Read(TObjArray*);
    inline AAMCCinput operator()() final {return data;}
    inline AAMCCinput GetNucleons(TObjArray* nucleons_in){this->Read(nucleons_in); return (*this)();}
private:
    AAMCCinput data;
};
