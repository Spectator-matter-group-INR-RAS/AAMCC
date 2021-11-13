#include "AAMCC.hh"

AAMCC::AAMCC(GRATEmanager *histoManager) {
   ExEnA = new ExcitationEnergy(histoManager->GetStatType(), histoManager->GetSourceA());
   ExEnB = new ExcitationEnergy(histoManager->GetStatType(), histoManager->GetSourceAb());

    if(histoManager->GetStatType()> 2){
        //Parameters for ALADIN parametrization
        G4double e_0=11.5*MeV;//was 8.13 MeV
        G4double sigma0 = 0.005; //was 0.01
        G4double c0 = 2; // From Bondorf 1995
        G4double sigmaE0 = 1*MeV;
        G4double b0 = 0.1;
        G4double Pe = 24*MeV;
        G4double Pm = 0.2;
        ExEnA->SetParametersALADIN(e_0,sigma0,c0);
        ExEnB->SetParametersALADIN(e_0,sigma0,c0);
        ExEnA->SetParametersParabolicApproximation(Pe, Pm, sigma0, b0, 0.01);
        ExEnB->SetParametersParabolicApproximation(Pe, Pm, sigma0, b0, 0.01);
        //ExEnA->SetParametersCorrectedALADIN(0.01,1000,sigma0,b0,0);
        //ExEnB->SetParametersCorrectedALADIN(0.01,1000,sigma0,b0,0);
    }
}

AAMCC::~AAMCC() = default;

template<class CollisionGeneratorOut>
G4ReactionProductVector* AAMCC::Do(CollisionGeneratorOut *CollisionOut) {
    __assert(std::is_base_of<VCollisionReader, CollisionGeneratorOut>::value, "type parameter of this method must be derived form VCollisionReader");
    NucleonVector nucleons = CollisionOut->GetNucleons();
    return nullptr;
}


