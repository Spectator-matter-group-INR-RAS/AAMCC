//#ifndef AAMCConstants_h
//#define AAMCConstants_h 1

#pragma once
#include "G4SystemOfUnits.hh"
#include "Nucleon.hh"
#include "G4Types.hh"
#include "G4String.hh"

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

    G4float b = 0;
    G4float ExEnA = 0;
    G4float ExEnB = 0;
    G4int id = 0;
    G4int Nhard = 0;
    G4int Ncoll = 0;
    G4int Ncollpp = 0;
    G4int Ncollpn = 0;
    G4int Ncollnn = 0;
    G4int Npart = 0;
    G4int NpartA = 0;
    G4int NpartB = 0;

    G4double FermiMomA_x = 0;
    G4double FermiMomA_y = 0;
    G4double FermiMomA_z = 0;
    G4double FermiMomB_x = 0;
    G4double FermiMomB_y = 0;
    G4double FermiMomB_z = 0;

    G4float PhiRotA = 0;
    G4float ThetaRotA = 0;
    G4float PhiRotB = 0;
    G4float ThetaRotB = 0;
    G4float Ecc[10]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    G4int ClustNumA = 0;
    G4int ClustNumB = 0;
    G4double d_MstA = 0;
    G4double d_MstB = 0;
    std::vector<G4int> A_cl{0};
    std::vector<G4int> Z_cl{0};
    std::vector<G4int> Ab_cl{0};
    std::vector<G4int> Zb_cl{0};

};


struct AAMCCinput{
    aamcc::NucleonVector nucleons;
    G4double FermiMomA_x = 0;
    G4double FermiMomA_y = 0;
    G4double FermiMomA_z = 0;
    G4double FermiMomB_x = 0;
    G4double FermiMomB_y = 0;
    G4double FermiMomB_z = 0;
};

struct AAMCCrun{
    G4int ZinitA = -1;
    G4int AinitA = -1;
    G4int ZinitB = -1;
    G4int AinitB = -1;

    G4String SysA = "";
    G4String SysB = "";
    G4String fileName = "";
    G4String fileRName = "";

    G4double KinEnPerNucl = -1.0;
    G4double SqrtSnn = -1.0;
    G4double pzA = -1.0;
    G4double pzB = -1.0;
    G4bool isCollider;

    G4int iterations = -1;

    G4double XsectNN = -1.0;
    G4double XsectTot = -1.0;

    G4String DeExModel = "";
    G4int ExExStatLabel = -1;

    G4double lowLimitB = -1.0; // MB if negative
    G4double upperLimitB = -2.0; // MB if upperLimitB < lowLimitB

    G4double CritDist = 2.7;

    G4bool InFileOrNot = false;

    G4bool isQMD = false;
};

//#endif
