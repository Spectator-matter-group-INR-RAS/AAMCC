#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4SystemOfUnits.hh"

#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4ExcitationHandler.hh"
#include "G4NucleiProperties.hh"
#include "G4Evaporation.hh"
#include "GRATEmanager.hh"
#include "InitialConditions.hh"
#include "ExcitationEnergy.hh"
#include "DeexcitationHandler.hh"
//#include "Nucleon.hh"
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

#include "G4UImanager.hh"
#include "G4IonTable.hh"
#include "G4GenericIon.hh"
#include "G4Ions.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclearLevelData.hh"

#include "../include/GMSTClustering.hh"

#include <fstream>
#include <cassert>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "VCollisionReader.hh"

class AAMCC {
public:
    explicit AAMCC(GRATEmanager* histoManager);
    ~AAMCC();

    template <class CollisionGeneratorOut>
           G4ReactionProductVector* Do(CollisionGeneratorOut * CollisionOut);

private:
    ExcitationEnergy* ExEnA;
    ExcitationEnergy* ExEnB;

};
