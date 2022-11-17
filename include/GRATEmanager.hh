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
#include "AAMCCReader.hh"

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
  
  void CalcNucleonDensity(TObjArray* nucleons_pre, G4double b); //calling it will cause the seg fault (uses GetHisto)
  
  void ToFile(AAMCCEvent* ev, NucleonVector* nucleons, std::function<void(AAMCCEvent*, AAMCCrun*, NucleonVector*)> toFile) {
    toFile(ev, &runData, nucleons);
  }

  void ReadInitCond(std::function<void(AAMCCrun*, InitialConditions* InCond)> readInitCond){
    readInitCond(&runData, InCond);
  }

  void WriteNucleonsCoordinatesInFile(GMSTClusterVector clusters_to_excit_A, GMSTClusterVector clusters_to_excit_B, G4double);

  inline G4String GetSysA() {return runData.SysA;}
  inline G4String GetSysB() {return runData.SysB;}
  inline G4String GetDeexModel() {return runData.DeExModel;}
  inline G4int GetSourceZ() {return this->GetInitialContidions().GetSourceZ();} // not used
  inline G4int GetSourceA() {return this->GetInitialContidions().GetSourceA();}
  inline G4int GetSourceZb() {return this->GetInitialContidions().GetSourceZb();} // not used
  inline G4int GetSourceAb() {return this->GetInitialContidions().GetSourceAb();}
  inline G4int GetStatType() {return runData.ExExStatLabel;}
  inline G4int GetIterations()  {return runData.iterations;}
  inline G4double GetXsectNN() {return this->GetInitialContidions().GetXsectNN();}
  inline G4double GetLowB() {return runData.lowLimitB;}
  inline G4double GetUpB() {return runData.upperLimitB;}
  inline InitialConditions GetInitialContidions() {return *InCond;}
  inline G4double GetCriticalDistance() {return runData.CritDist;}
  inline G4double GetAngle() {return CLHEP::pi*angle/180;} // left for the future development of polarized beams
  inline G4bool ToFileOrNot() {return runData.InFileOrNot;}

  inline void SetXsectTot(G4double Xsect){runData.XsectTot = Xsect;}
  
  private:

  TH1D*  histo[20];

  G4double angle = 0.0; // left for future development of polarized beams

  InitialConditions* InCond = new InitialConditions();
  AAMCCrun runData;
};

#endif
