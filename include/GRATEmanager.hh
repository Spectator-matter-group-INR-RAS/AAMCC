#ifndef GRATEmanager_h
#define GRATEmanager_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "TTree.h"
#include "TParticle.h"
#include <fstream>
#include <cmath>
#include <cassert>

#include "../TGlauber/TGlauberMC.hh"
#include "../TGlauber/TGlauNucleon.hh"
#include "../TGlauber/TGlauNucleus.hh"

#include "InitialConditions.hh"
#include "GMSTClustering.hh"

#include "AAMCConstants.hh"
#include "AAMCC.hh"
#include "VWriter.hh"

class TFile;
class TH1D;
class TH2D;

class GRATEmanager
{
  public:

  explicit GRATEmanager();
  explicit GRATEmanager(TFile* file, AAMCCrun (*getTheRunData)(TFile* file));
  explicit GRATEmanager(AAMCCrun (*getTheRunData)());
   ~GRATEmanager();

  public:

  TH1D* GetHisto(G4int id) {return histo[id];}; //zombie function waining to the migration of this->CalcNucleonDensity() to the WriteToFile()
  void CalcNucleonDensity(TObjArray* nucleons_pre, G4double b); //calling it will cause the seg fault
  void ToFile(AAMCCEvent* ev, NucleonVector* nucleons, std::function<void(AAMCCEvent*, AAMCCrun*, NucleonVector*)> toFile) {
      runData.XsectNN = XsectNN; runData.XsectTot = XsectTot;
      toFile(ev, &runData, nucleons);
  }
  /*
  void ToFile(AAMCCEvent* ev, NucleonVector* nucleons, VWriter* writer){
      runData.XsectNN = XsectNN; runData.XsectTot = XsectTot;
      (*writer)(ev, &runData, nucleons);
  }*/
  void WriteNucleonsCoordinatesInFile(GMSTClusterVector clusters_to_excit_A, GMSTClusterVector clusters_to_excit_B, G4double);

  inline G4String GetSysA() {return SysA;}
  inline G4String GetSysB() {return SysB;}
  inline G4String GetDeexModel() {return DeExModel;};
  inline G4int GetSourceZ() {return this->GetInitialContidions().GetSourceZ();}
  inline G4int GetSourceA() {return this->GetInitialContidions().GetSourceA();}
  inline G4int GetSourceZb() {return this->GetInitialContidions().GetSourceZb();}
  inline G4int GetSourceAb() {return this->GetInitialContidions().GetSourceAb();}
  inline G4int GetStatType() {return ExEnStatLabel;}
  inline G4int GetIterations()  {return iterations;};
  inline G4double GetXsectNN() {return this->GetInitialContidions().GetXsectNN();}
  inline G4double GetKinEn() {return KinEn;};
  inline G4double GetSqrtSnn() {return SqrtSnn;};
  inline G4double GetLowEn() {return lowLimitExEn;};
  inline G4double GetUpEn() {return upperLimitExEn;};
  inline G4double GetLowEnB() {return lowLimitExEnB;};
  inline G4double GetUpEnB() {return upperLimitExEnB;};
  inline G4double GetLowB() {return lowLimitB;};
  inline G4double GetUpB() {return upperLimitB;};
  inline G4bool   WriteMomentum() {return wM;};
  inline G4bool   WritePseudorapidity() {return wP;};
  inline InitialConditions GetInitialContidions() {return *InCond;};
  inline G4double GetCriticalDistance() {return CritDist;}
  inline G4double GetAngle() {return CLHEP::pi*angle/180;} //Left for the future development of polarized beams
  inline G4bool ToFileOrNot() {return InFileOrNot;}

  inline void SetXsectTot(G4double Xsect){XsectTot = Xsect; runData.XsectTot = XsectTot;};


  
  private:

  TH1D*  histo[20];
  TH2D*  histo2[10];


  G4int sourceZ;
  G4int sourceA;
  G4int sourceZb;
  G4int sourceAb;
  G4int iterations; 
  G4int ExEnStatLabel;
  G4bool NucleusInputLabel;
  G4bool IsCollider;
  G4bool InFileOrNot;

  G4String fileName;
  G4String fileType;
  G4String fileOpenPath;
  G4String DeExModel;

  G4String     fileFullName;
  G4String     SysA;
  G4String     SysB;
  G4int        compressionFactor;

  G4double XsectNN = -1;
  G4double XsectTot = 0;

  G4double KinEn;
  G4double SqrtSnn;
  G4double lowLimitExEn;
  G4double upperLimitExEn;
  G4double lowLimitExEnB;
  G4double upperLimitExEnB;
  G4double lowLimitB;
  G4double upperLimitB;
  G4int    binsExEn;
  G4int    eventsPerBin;

  G4double CritDist;
  G4double angle;

  G4bool   wM;
  G4bool   wP = false;

  InitialConditions* InCond = new InitialConditions();
  AAMCCEvent event;
  AAMCCrun runData;

  std::ifstream XsectFile;

};

#endif
