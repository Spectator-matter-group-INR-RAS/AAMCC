#ifndef AAMCC_h
#define AAMCC_h 1


#include "G4SystemOfUnits.hh"

#include "G4ReactionProductVector.hh"
#include "InitialConditions.hh"
#include "ExcitationEnergy.hh"
#include "DeexcitationHandler.hh"
#include "Nucleon.hh"
#include "AAMCConstants.hh"
#include "GRATEPhysicsList.hh"
#include "../TGlauber/TGlauberMC.hh"
#include "../TGlauber/TGlauNucleon.hh"
#include "../TGlauber/TGlauNucleus.hh"
#include "TVector3.h"
#include "TObjArray.h"
#include "TObject.h"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4Threading.hh"

#include "G4IonTable.hh"
#include "G4GenericIon.hh"
#include "G4Ions.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclearLevelData.hh"

#include "../include/GMSTClustering.hh"

#include <fstream>
#include <cassert>
#include <memory>
#include <functional>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

//#ifndef VCollisionReader_h
//#define VCollisionReader_h 1
#include "VCollisionReader.hh"
//#endif
#include "URun.h"

class AAMCC {
public:
    explicit AAMCC();
    ~AAMCC();

    template <class CollisionGeneratorOut>
           G4ReactionProductVector* Do(CollisionGeneratorOut * CollisionOut);

    inline G4double GetZpfA(AAMCCEvent ev){G4double Zpf = 0; for(auto iZ = ev.ChargeOnSideA.begin(); iZ != ev.ChargeOnSideA.end(); ++iZ){Zpf += *iZ;}; return Zpf;};
    inline G4double GetZpfB(AAMCCEvent ev){G4double Zpf = 0; for(auto iZ = ev.ChargeOnSideB.begin(); iZ != ev.ChargeOnSideB.end(); ++iZ){Zpf += *iZ;}; return Zpf;};
    inline G4double GetApfA(AAMCCEvent ev){G4double Apf = 0; for(auto iA = ev.MassOnSideA.begin(); iA != ev.MassOnSideA.end(); ++iA){Apf += *iA;}; return Apf;};
    inline G4double GetApfB(AAMCCEvent ev){G4double Apf = 0; for(auto iA = ev.MassOnSideB.begin(); iA != ev.MassOnSideB.end(); ++iA){Apf += *iA;}; return Apf;};
    inline G4double GetMultA(AAMCCEvent ev){return ev.ChargeOnSideA.size();};
    inline G4double GetMultB(AAMCCEvent ev){return ev.ChargeOnSideB.size();};


private:
    ExcitationEnergy* ExEnA;
    ExcitationEnergy* ExEnB;

};

#endif