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

#ifndef VCollisionReader_h
#define VCollisionReader_h 1
#include "VCollisionReader.hh"
#endif

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

struct AAMCCEvent{

    std::vector<G4float> MassOnSideA;
    std::vector<G4float> MassOnSideB;
    std::vector<G4float> ChargeOnSideA;
    std::vector<G4float> ChargeOnSideB;
    std::vector<G4double> pXonSideA;
    std::vector<G4double> pYonSideA;
    std::vector<G4double> pZonSideA;
    std::vector<G4double> pXonSideB;
    std::vector<G4double> pYonSideB;
    std::vector<G4double> pZonSideB;
    std::vector<G4double> pseudorapidity_A;
    std::vector<G4double> pseudorapidity_B;

    G4float b;
    G4float ExEn;
    G4int id;
    G4int Nhard;
    G4int Ncoll;
    G4int Ncollpp;
    G4int Ncollpn;
    G4int Ncollnn;
    G4int Npart;
    G4int NpartA;
    G4int NpartB;

    G4double FermiMomA_x;
    G4double FermiMomA_y;
    G4double FermiMomA_z;
    G4double FermiMomB_x;
    G4double FermiMomB_y;
    G4double FermiMomB_z;

    G4float PhiRotA;
    G4float ThetaRotA;
    G4float PhiRotB;
    G4float ThetaRotB;
    G4float Ecc[10];

    G4int ClustNumA;
    G4int ClustNumB;
    G4double d_MstA;
    G4double d_MstB;
    std::vector<G4int> A_cl;
    std::vector<G4int> Z_cl;
    std::vector<G4int> Ab_cl;
    std::vector<G4int> Zb_cl;

};
