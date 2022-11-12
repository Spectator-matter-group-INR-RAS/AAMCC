#ifndef GRATEmanager_h
#define GRATEmanager_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "TTree.h"
#include "TParticle.h"
#include <fstream>
#include <cmath>

#include "../TGlauber/TGlauberMC.hh"
#include "../TGlauber/TGlauNucleon.hh"
#include "../TGlauber/TGlauNucleus.hh"

#include "InitialConditions.hh"
#include "GMSTClustering.hh"

#include "AAMCConstants.hh"
#include "AAMCC.hh"

class TFile;
class TH1D;
class TH2D;

class GRATEmanager
{
  public:

  explicit GRATEmanager();
   ~GRATEmanager();

  public:

  TH1D* GetHisto(G4int id) {return histo[id];}; //zombie function waining to the migration of this->CalcNucleonDensity() to the WriteToFile()
  void CalcNucleonDensity(TObjArray* nucleons_pre, G4double b); //calling it will cause the seg fault
  void ToFile(AAMCCEvent* ev, NucleonVector* nucleons, std::function<void(AAMCCEvent*, AAMCCrun*, NucleonVector*)> toFile) {
      runData.XsectNN = XsectNN; runData.XsectTot = XsectTot;
      toFile(ev, &runData, nucleons);
  }
  void WriteNucleonsCoordinatesInFile(GMSTClusterVector clusters_to_excit_A, GMSTClusterVector clusters_to_excit_B, G4double);

  inline G4String GetSysA() {return runData.SysA;}
  inline G4String GetSysB() {return runData.SysB;}
  inline G4String GetDeexModel() {return runData.DeExModel;};
  inline G4int GetSourceZ() {return this->GetInitialContidions().GetSourceZ();} // not used
  inline G4int GetSourceA() {return this->GetInitialContidions().GetSourceA();}
  inline G4int GetSourceZb() {return this->GetInitialContidions().GetSourceZb();} // not used
  inline G4int GetSourceAb() {return this->GetInitialContidions().GetSourceAb();}
  inline G4int GetStatType() {return runData.ExExStatLabel;}
  inline G4int GetIterations()  {return runData.iterations;};
  inline G4double GetXsectNN() {return this->GetInitialContidions().GetXsectNN();}
  inline G4double GetLowB() {return lowLimitB;};
  inline G4double GetUpB() {return upperLimitB;};
  inline InitialConditions GetInitialContidions() {return *InCond;};
  inline G4double GetCriticalDistance() {return CritDist;}
  inline G4double GetAngle() {return CLHEP::pi*angle/180;} //Left for the future development of polarized beams
  inline G4bool ToFileOrNot() {return InFileOrNot;}

  inline void SetXsectTot(G4double Xsect){XsectTot = Xsect; runData.XsectTot = XsectTot;};


  
  private:

  TH1D*  histo[20];

  G4bool NucleusInputLabel = FALSE;

  G4bool InFileOrNot = FALSE;
  G4bool IsInitFile = FALSE;
  G4int  AbrasionModelInt = -1;

  G4double XsectNN = -1.0;
  G4double XsectTot = -1.0;

  G4double KinEn;
  G4double SqrtSnn;

  G4double lowLimitB = -1.0; // MB if negative
  G4double upperLimitB = -2.0; // MB if upperLimitB < lowLimitB

  G4double CritDist;
  G4double angle;

  InitialConditions* InCond = new InitialConditions();
  AAMCCrun runData;

  TString inputFileName;

};

#endif
