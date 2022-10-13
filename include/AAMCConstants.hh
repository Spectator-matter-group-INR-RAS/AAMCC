//#ifndef AAMCConstants_h
//#define AAMCConstants_h 1

#pragma once
#include "G4SystemOfUnits.hh"
#include "Nucleon.hh"
static constexpr double nucleonAverMass = 0.93891875434*CLHEP::GeV;

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
    G4float ExEnA;
    G4float ExEnB;
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


struct AAMCCinput{
    NucleonVector nucleons;
    G4int ZinitA;
    G4int AinitA;
    G4int ZinitB;
    G4int AinitB;

    G4String SysA;
    G4String SysB;
    G4String fileName;

    G4double KinEnPerNucl;
    G4bool isCollider;
};

struct AAMCCrun{
    G4int ZinitA;
    G4int AinitA;
    G4int ZinitB;
    G4int AinitB;

    G4String SysA;
    G4String SysB;
    G4String fileName;

    G4double KinEnPerNucl;
    G4double SqrtSnn;
    G4double pzA;
    G4double pzB;
    G4bool isCollider;

    G4int iterations;

    G4double XsectNN;
    G4double XsectTot;
};


//#endif
