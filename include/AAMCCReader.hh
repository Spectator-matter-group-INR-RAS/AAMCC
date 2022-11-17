#ifndef GRATE_AAMCCREADER_HH
#define GRATE_AAMCCREADER_HH

#include "VReader.hh"
#include "G4SystemOfUnits.hh"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include <mutex>
#include <memory>

class AAMCCReader: public VReader {
public:
    AAMCCReader();
    ~AAMCCReader();
    void operator()(AAMCCrun* run, InitialConditions* InCond);
private:

    G4bool IsInitFile = FALSE;
    G4int  AbrasionModelInt = -1;
    TString inputFileName;
    G4bool NucleusInputLabel = FALSE;

    G4double KinEn = -1.0;
    G4double SqrtSnn = -1.0;

    void ReadSecondStage(AAMCCrun* runData, InitialConditions* InCond);
    void ReadFirstStage(AAMCCrun* runData, InitialConditions* InCond);
};


#endif //GRATE_AAMCCREADER_HH
