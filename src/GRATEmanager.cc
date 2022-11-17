#include "GRATEmanager.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h" 
#include "TMath.h"


GRATEmanager::GRATEmanager(){}

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
