#ifndef ExcitationEnergy_h
#define ExcitationEnergy_h 1

#include "globals.hh"
#include "Randomize.hh"
#include "G4ExceptionHandler.hh"
#include "cmath"
#include <fstream>
#include <libgen.h>

class ExcitationEnergy {
public:
    ExcitationEnergy(G4int ExEnLabel_in, G4int initA_in);

    ~ExcitationEnergy();

public:
    G4double GetEnergy(G4int A);

    G4double GetEnergyEricson(G4int A);

    G4double GetEnergyGaimardSchmidt(G4int A);

    G4double GetEnergyALADIN(G4int A);

    G4double GetEnergyCorrectedALADIN(G4int A);

    G4double GetEnergyParabolicApproximation(G4int A);
   
    G4double GetEnergyDampEricson(G4int A);

    G4double GetEnergyHybridFit(G4int A);

    void SetInitNuclMass(G4int initA_in);

    void SetParametersEricson(G4double g0_in);

    void SetParametersGaimardSchmidt(G4double g0_in, G4double g1_in);

    void SetParametersALADIN(G4double e0_in, G4double sigma0_in, G4double b0_in);

    void SetParametersCorrectedALADIN(G4double e0_in,G4double c0_in, G4double sigma0_in, G4double b0_in, G4double b1_in);

    void SetParametersCorrectedALADINFromFile();

    void SetParametersParabolicApproximation(G4double Pe_in, G4double Pm_in, G4double sigmaP_in, G4double bP0_in, G4double bP1_in);

    void SetParametersHybridFit(G4double a0_in, G4double a1_in, G4double a2_in, G4double a3_in, G4double a4_in, G4double a5_in, G4double a6_in, G4double sigma1_in, G4double sigma2_in, G4double sigma3_in);

private:
     G4double g0;
     G4double g1;
     G4double e0;
     G4double sigma0;
     G4double b0;
     G4double c0;
     G4double b1;
     G4double Pe;
     G4double Pm;
     G4double sigmaP;
     G4double bP0;
     G4double bP1;
     G4double alphaSwitch;

    G4int initA;
    G4int ExEnLabel;

     G4double LowExEn;
     G4double UpExEn;
     G4double Ebound;

     G4double a0;
     G4double a1;
     G4double a2;
     G4double a3;
     G4double a4;
     G4double a5;
     G4double a6;
     G4double sigma1;
     G4double sigma2;
     G4double sigma3;

     std::ifstream ParamFile;
};
#endif
