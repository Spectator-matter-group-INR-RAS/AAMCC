#include "GRATEmanager.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h" 
#include "TMath.h"


GRATEmanager::GRATEmanager()
  : sourceZ(0), sourceA(0), KinEn(-1.), SqrtSnn(-1.), XsectNN(-1.), lowLimitExEn( 0.), upperLimitExEn( 100.), binsExEn(1), eventsPerBin(1), ExEnStatLabel(-1), iterations(1), wM(0), wP(0), NucleusInputLabel(0), IsCollider(0), upperLimitB(-1), CritDist(2.7), angle(0), DeExModel("")
{  
  std::cout << "######### Abrasion-Ablation model using Glauber Monte Carlo and Geant4" <<std::endl;
  while(!NucleusInputLabel){
    std::cout << "Please enter colliding nucleus name (side A). U, U2, Pb, Pbrw, Pbpn, Pbpnrw, Au, Aurw, Au2, Au2rw, Xe, Ag, Br, Cu, Ca2 (with SRC), Ar, Al, O, O2 (with SRC), Oho (for HO param.), C, He4, He3, H3, d, p is available : ";
    std::cin >> SysA;
    NucleusInputLabel = InCond->SetSysA(SysA);
    sourceA = InCond->GetSourceA();
    sourceZ = InCond->GetSourceZ();
  }

  NucleusInputLabel = 0;

  while(!NucleusInputLabel){
    std::cout << "Please enter colliding nucleus name (side B). U, U2, Pb, Pbrw, Pbpn, Pbpnrw, Au, Aurw, Au2, Au2rw, Xe, Ag, Br, Cu, Ca2 (with SRC), Ar, Al, O, O2 (with SRC), Oho (for HO param.), C, He4, He3, H3, d, p is available : ";
    std::cin >> SysB;
    NucleusInputLabel = InCond->SetSysB(SysB);
    sourceAb = InCond->GetSourceAb();
    sourceZb = InCond->GetSourceZb();
  }

  std::cout<<"Input lower limit for impact parameter in fm (MB if negative) : ";
  std::cin >> lowLimitB;

  if(lowLimitB > -0.000001){
    std::cout<<"Input upper limit for impact parameter in fm (MB if negative) : ";
    std::cin >> upperLimitB;
  }
  else if(0){
    lowLimitB = 0;
    upperLimitB = 20;
  }

  std::cout<<"Do you want to calculate collisions for collider or for fixed target geometry (1 for collider, 0 for fixed target) : ";
  std::cin >> IsCollider;

  while ( KinEn<280.0*MeV/GeV && SqrtSnn <0.) {
    if (!IsCollider) {
      std::cout << "Please enter kinetic enegy of projectile nucleus (per nucleon in GeV) : ";
      std::cin >> KinEn;
      if(KinEn<280.0*MeV/GeV){
        std::cout << "AAMCC works at kinetic energies above 280A MeV, please input higher energy!" << std::endl;
      }
    }
    else {
      std::cout << "Please enter s^1/2 of colliding nuclei (per nucleon in GeV) : ";
      std::cin >> SqrtSnn;
      G4double NuclMass = nucleonAverMass/CLHEP::GeV;
      G4double CollKinCheck = (SqrtSnn/2.0 - NuclMass);
      G4double KinEnAtFixTargetCheck = (2.0*(CollKinCheck + NuclMass)*(CollKinCheck + NuclMass)/(NuclMass*NuclMass) - 1.0)*NuclMass;
      if(KinEnAtFixTargetCheck < 280.0*MeV/GeV){
        std::cout << "AAMCC works at kinetic energies above 280A MeV, please input higher energy!" << std::endl;
        SqrtSnn = -1.0;
      }
    }
  }

  InCond->SetCollider(IsCollider);
  if(IsCollider){
    InCond->SetKinematics(SqrtSnn);
    SqrtSnn = InCond->GetSqrtSnn();
    KinEn = InCond->GetKinEnergy();
  }
  else{
    InCond->SetKinematics(KinEn);
    SqrtSnn = InCond->GetSqrtSnn();
    KinEn = InCond->GetKinEnergy();
  }

  sourceA=InCond->GetSourceA();
  sourceAb=InCond->GetSourceAb();
  lowLimitExEnB = lowLimitExEn;
  upperLimitExEnB = upperLimitExEn;
  lowLimitExEn *= sourceA;
  lowLimitExEn *= MeV;
  upperLimitExEn *=sourceA;
  upperLimitExEn *=MeV;
  lowLimitExEnB *= sourceAb;
  lowLimitExEnB *= MeV;
  upperLimitExEnB *=sourceAb;
  upperLimitExEnB *=MeV;

  while ((ExEnStatLabel < 0) || (ExEnStatLabel > 7) || (upperLimitExEn < lowLimitExEn) ) {
    std::cout << "Please choose the level density function to be used: 1 - Ericson, 2 - Gaimard-Schmidt, 3 - ALADIN parametrization, 4 - Hybrid of 1 and 3, 7 - Fit of Hybrid : ";
    std::cin >> ExEnStatLabel;
  }

  while(CritDist < 0) {
    std::cout<<"Please enter critical distance (in fm) : ";
    std::cin >> CritDist;
  }

  while(DeExModel.empty()){
      std::cout<<"Choose a model for fragment deexcitation. G4, ABLAXX, AAMCC or MIX (random mix of G4, ABLAXX and AAMCC) options are available: ";
      std::cin>>DeExModel;
  }

  std::cout<<"Write coordinates of nucleons in the text file or not (one event)? (1 - yes, 0 - no): ";
  std::cin >> InFileOrNot;

  /*
  if(!InFileOrNot){
    std::cout<<"Write momentum of each fragment?  (1 - yes, 0 - no) ";
    std::cin>>wM;

    std::cout<<"Write pseudorapidity of each fragments?  (1 - yes, 0 - no) ";
    std::cin>>wP;
  }
*/

  while ( ((iterations<2) || (iterations>10000000)) && !InFileOrNot ) {
    std::cout<<"Please enter number of events to be generated: ";
    std::cin >> iterations;
  }

  std::cout << "Please enter the file name to write histograms (.root will be supplied): ";
  std::cin >> fileName; 

  XsectNN = InCond->GetXsectNN();

  runData.ZinitA = sourceZ; runData.ZinitB = sourceZb; runData.AinitA = sourceA; runData.AinitB = sourceAb; runData.SysA = SysA; runData.SysB = SysB; runData.fileName = fileName;
  runData.KinEnPerNucl = InCond->GetKinEnergyPerNucl(); runData.pzA = InCond->GetPzA(); runData.pzB = InCond->GetPzB(); runData.SqrtSnn = InCond->GetSqrtSnn(); runData.XsectNN = XsectNN;
  runData.isCollider = IsCollider; runData.iterations = iterations; runData.DeExModel = DeExModel; runData.ExExStatLabel = ExEnStatLabel;

  //for (G4int j=0; j<20; j++) histo[j] = 0;
  //for (G4int l=0; l<10; l++) histo2[l] = 0;
}

GRATEmanager::GRATEmanager(TFile* file, AAMCCrun (*getTheRunData)(TFile* file))
        : sourceZ(0), sourceA(0), KinEn(-1.), SqrtSnn(-1.), XsectNN(-1.), lowLimitExEn( 0.), upperLimitExEn( 100.), binsExEn(1), eventsPerBin(1), ExEnStatLabel(-1), iterations(1), wM(0), wP(0), NucleusInputLabel(0), IsCollider(0), upperLimitB(-1), CritDist(2.7), angle(0), DeExModel("")
{
    runData = getTheRunData(file);
    InCond->SetConditions(runData);
}

GRATEmanager::GRATEmanager(AAMCCrun (*getTheRunData)())
        : sourceZ(0), sourceA(0), KinEn(-1.), SqrtSnn(-1.), XsectNN(-1.), lowLimitExEn( 0.), upperLimitExEn( 100.), binsExEn(1), eventsPerBin(1), ExEnStatLabel(-1), iterations(1), wM(0), wP(0), NucleusInputLabel(0), IsCollider(0), upperLimitB(-1), CritDist(2.7), angle(0), DeExModel("")
{
    runData = getTheRunData();
    InCond->SetConditions(runData);
}


GRATEmanager::~GRATEmanager()
{
    delete InCond;
}


void GRATEmanager::CalcNucleonDensity(TObjArray* nucleons_pre, G4double b)
{
  G4float X_pre, Y_pre, Z_pre, R_pre;
  for (G4int iter = 0; iter < nucleons_pre->GetEntries(); iter++) {
    TGlauNucleon* nucleon_pre = (TGlauNucleon*)(nucleons_pre->At(iter));
    if (nucleon_pre->IsInNucleusA()){X_pre = nucleon_pre->GetX() + b/2;}
    if (nucleon_pre->IsInNucleusB()){X_pre = nucleon_pre->GetX() - b/2;}
    Y_pre = nucleon_pre->GetY();
    Z_pre = nucleon_pre->GetZ();
    R_pre = std::sqrt(pow(X_pre,2) + pow(Y_pre,2) + pow(Z_pre,2));
    if (nucleon_pre->IsProton() && nucleon_pre->IsInNucleusA()) { (*this).GetHisto(9)->Fill(R_pre); }
    if (nucleon_pre->IsNeutron() && nucleon_pre->IsInNucleusA()) { (*this).GetHisto(8)->Fill(R_pre); }
    if (nucleon_pre->IsProton() && nucleon_pre->IsInNucleusB()) { (*this).GetHisto(11)->Fill(R_pre); }
    if (nucleon_pre->IsNeutron() && nucleon_pre->IsInNucleusB()) { (*this).GetHisto(10)->Fill(R_pre); }
  }
}

void GRATEmanager::WriteNucleonsCoordinatesInFile(GMSTClusterVector clusters_to_excit_A, GMSTClusterVector clusters_to_excit_B, G4double b)
{
  ofstream Event("Event.txt");

  Event<<"Clusters on side A: \n";

  for(G4int i = 0; i < clusters_to_excit_A.size(); ++i) {
    Event<<"Nucleons from "<<i+1<<"-th cluster\n";
    for(G4int iCoord = 0; iCoord<(clusters_to_excit_A[i]).GetCoordinates().size(); iCoord++){
      Event<<"X = "<<(clusters_to_excit_A[i]).GetCoordinates().at(iCoord).X()<<" Y = "<<(clusters_to_excit_A[i]).GetCoordinates().at(iCoord).Y()<<" Z = "<<(clusters_to_excit_A[i]).GetCoordinates().at(iCoord).Z()<<"\n";
    }
  }

  Event<<"\n";  Event<<"Clusters on side B: \n";

  for(G4int i = 0; i < clusters_to_excit_B.size(); ++i) {
    Event<<"Nucleons from "<<i+1<<"-th cluster\n";
    for(G4int iCoord = 0; iCoord<(clusters_to_excit_B[i]).GetCoordinates().size(); iCoord++){
      Event<<"X = "<<(clusters_to_excit_B[i]).GetCoordinates().at(iCoord).X()<<" Y = "<<(clusters_to_excit_B[i]).GetCoordinates().at(iCoord).Y()<<" Z = "<<(clusters_to_excit_B[i]).GetCoordinates().at(iCoord).Z()<<"\n";
    }
  }
  Event<<"\nb = "<<b<<" fm \n";
}