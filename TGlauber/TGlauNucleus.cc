/*
 $Id: runglauber.C 173 2018-06-07 19:47:47Z loizides $
 -------------------------------------------------------------------------------------
 Latest documentation: https://arxiv.org/abs/1710.07098
 -------------------------------------------------------------------------------------
 To run the code, you need to have the ROOT (http://root.cern.ch/drupal/)
 environment. On the root prompt, then enter
 root [0] gSystem->Load("libMathMore")
 root [1] .L runglauber_X.Y.C+
 (where X.Y denotes the version number).
 If you do not have libMathMore comment out "#define HAVE_MATHMORE" below.
 See the documentation for more information.
 -------------------------------------------------------------------------------------
 v3.1:
  Fixes related to spherical nuclei, as well as consistent set of reweighted profiles 
  for Cu, Au and Xe, see https://arxiv.org/abs/1710.07098v2
 -------------------------------------------------------------------------------------
 v3.0:
  Major update to include separate profile for protons and neutrons, placement of nucleon 
  dof on lattice, as well as reweighted profiles for recentering, 
  see https://arxiv.org/abs/1710.07098v1
 -------------------------------------------------------------------------------------
 v2.6:
  Includes runAndCalcDens macro, as well as definition for Al, and fixes beta4 for Si2,
  see https://arxiv.org/abs/1408.2549v8
 -------------------------------------------------------------------------------------
 v2.5:
  Include core/corona determination in Npart, and if requested for area from mc and eccentricity,
  as well as various Xe parameterizations including deformation,
  see https://arxiv.org/abs/1408.2549v7
 -------------------------------------------------------------------------------------
 v2.4: 
  Minor update to include Xenon and fix of the TGlauberMC::Draw function, 
  see https://arxiv.org/abs/1408.2549v4
 -------------------------------------------------------------------------------------
 v2.3: 
  Small bugfixes, see https://arxiv.org/abs/1408.2549v3
 -------------------------------------------------------------------------------------
 v2.2:
  Minor update to provide higher harmonic eccentricities up to n=5, and the average
  nucleon--nucleon impact parameter (bNN) in tree output. 
 -------------------------------------------------------------------------------------
 v2.1: 
  Minor update to include more proton pdfs, see https://arxiv.org/abs/1408.2549v2
 -------------------------------------------------------------------------------------
 v2.0: 
  First major update with inclusion of Tritium, Helium-3, and Uranium, as well as the 
  treatment of deformed nuclei and Glauber-Gribov fluctuations of the proton in p+A 
  collisions, see https://arxiv.org/abs/1408.2549v1
 -------------------------------------------------------------------------------------
 v1.1: 
  First public release of the PHOBOS MC Glauber, see https://arxiv.org/abs/0805.4411
 -------------------------------------------------------------------------------------

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>
*/
#include "TGlauNucleon.hh"
#include "TGlauNucleus.hh"


// //---------------------------------------------------------------------------------
// ClassImp(TGlauNucleus)
// //---------------------------------------------------------------------------------
TGlauNucleus::TGlauNucleus(const char* iname, Int_t iN, Double_t iR, Double_t ia, Double_t iw, TF1* ifunc) : 
  TNamed(iname,""),
  fN(iN),fR(iR),fA(ia),fW(iw),fR2(0),fA2(0),fW2(0),fBeta2(0),fBeta4(0),
  fMinDist(0.4),fNodeDist(0.0),fSmearing(0.0),fRecenter(1),fLattice(0),fSmax(99),
  fF(0),fTrials(0),fNonSmeared(0),fFunc1(ifunc),fFunc2(0),fFunc3(0),fNucleons(0),
  fPhiRot(0),fThetaRot(0),fNucCounter(-1),fIsUsed(0),fMaxR(14)
{
  if (fN==0) {
    cout << "Setting up nucleus " << iname << endl;
    Lookup(iname);
  }
}

TGlauNucleus::~TGlauNucleus()
{
  if (fIsUsed)
    delete fIsUsed;
  if (fNucleons)
    delete fNucleons;
  delete fFunc1;
  delete fFunc2;
  delete fFunc3;
}

void TGlauNucleus::Draw(Double_t xs, Int_t colp, Int_t cols)
{
  Double_t r = 0.5*TMath::Sqrt(xs/TMath::Pi()/10.);
  TEllipse en;
  en.SetLineStyle(1);
  en.SetLineWidth(1);
  en.SetFillStyle(1001);
  for (Int_t i = 0; i<fNucleons->GetEntries(); ++i) {
    TGlauNucleon* gn = (TGlauNucleon*) fNucleons->At(i);
    if (!gn->IsSpectator()) {
      en.SetFillColor(colp);
      en.DrawEllipse(gn->GetX(),gn->GetY(),r,r,0,360,0,"");
    } else {
      en.SetFillColor(cols);
      en.SetFillStyle(1001);
      en.DrawEllipse(gn->GetX(),gn->GetY(),r,r,0,360,0,"");
    }
  }
}

void TGlauNucleus::Lookup(const char* name)
{
  SetName(name);
  TString tmp(name);
  Double_t r0=0, r1=0, r2=0;

  if      (TString(name) == "p")       {fN = 1;   fR = 0.234;      fA = 0;      fW =  0;       fF = 0;  fZ=1;}
  else if (TString(name) == "pg")      {fN = 1;   fR = 0.514;      fA = 0;      fW =  0;       fF = 9;  fZ=1;} 
  else if (TString(name) == "pdg")     {fN = 1;   fR = 1;          fA = 0;      fW =  0;       fF = 10; fZ=1;} // from arXiv:1101.5953
  else if (TString(name) == "dpf")     {fN = 2;   fR = 0.01;       fA = 0.5882; fW =  0;       fF = 1;  fZ=1;} // deuteron 2pf (tuned to Hulthen)
  else if (TString(name) == "dh")      {fN = 2;   fR = 0.2283;     fA = 1.1765; fW =  0;       fF = 3;  fZ=1;} // deuteron Hulthen free
  else if (TString(name) == "d")       {fN = 2;   fR = 0.2283;     fA = 1.1765; fW =  0;       fF = 4;  fZ=1;} // deuteron Hulthen constrained
  else if (TString(name) == "He3")     {fN = 3;   fR = 0.00;       fA = 0.0000; fW =  0;       fF = 6;  fZ=1;} // read configurations from file
  else if (TString(name) == "H3")      {fN = 3;   fR = 0.00;       fA = 0.0000; fW =  0;       fF = 6;  fZ=2;} // read configurations from file
  else if (TString(name) == "He4")     {fN = 4;   fR = 0.00;       fA = 0.0000; fW =  0;       fF = 6;  fZ=2;} // read configurations from file
  else if (TString(name) == "C")       {fN = 12;  fR = 2.608;      fA = 0.513;  fW = -0.051;   fF = 6;  fZ=6;} // read configurations from file  
  else if (TString(name) == "O")       {fN = 16;  fR = 2.608;      fA = 0.513;  fW = -0.051;   fF = 6;  fZ=8;} // read configurations from file
  else if (TString(name) == "O2")    {fN = 16;  fR = 2.608;      fA = 0.513;  fW = -0.051;   fF = 6;  fZ=8;} // read configurations from file
  else if (TString(name) == "Opar")    {fN = 16;  fR = 2.608;      fA = 0.513;  fW = -0.051;   fF = 1;  fZ=8;} // WS parameterization
  else if (TString(name) == "Oho")     {fN = 16;  fR = 1.833;      fA = 1.544;  fW =  0;       fF = 15; fZ=8;} // Harmonic oscillator parameterization
  else if (TString(name) == "Al")      {fN = 27;  fR = 3.34;       fA = 0.580;  fW = 0.0;      fF = 8;  fZ=13; fBeta2=-0.448; fBeta4=0.239;}
  else if (TString(name) == "Si")      {fN = 28;  fR = 3.34;       fA = 0.580;  fW = -0.233;   fF = 1;  fZ=14;}
  else if (TString(name) == "Si2")     {fN = 28;  fR = 3.34;       fA = 0.580;  fW =  0;       fF = 8;  fZ=14; fBeta2=-0.478; fBeta4=0.250;}
  else if (TString(name) == "S")       {fN = 32;  fR = 2.54;       fA = 2.191;  fW =  0.16;    fF = 2;  fZ=16;}
  else if (TString(name) == "Ar")      {fN = 40;  fR = 3.53;       fA = 0.542;  fW =  0;       fF = 1;  fZ=18;}
  else if (TString(name) == "Ca")      {fN = 40;  fR = 3.766;      fA = 0.586;  fW = -0.161;   fF = 1;  fZ=20;}
  else if (TString(name) == "Ca2")     {fN = 40;  fR = 3.766;      fA = 0.586;  fW = -0.161;   fF = 6;  fZ=20;} // read configuration from file
  else if (TString(name) == "Ni")      {fN = 58;  fR = 4.309;      fA = 0.517;  fW = -0.1308;  fF = 1;  fZ=28;}
  else if (TString(name) == "Cu")      {fN = 63;  fR = 4.20;       fA = 0.596;  fW =  0;       fF = 1;  fZ=29;}
  else if (TString(name) == "Curw ")   {fN = 63;  fR = 4.20;       fA = 0.596;  fW =  0;       fF = 12; fZ=29; r0=1.00898; r1=-0.000790403; r2=-0.000389897;} 
  else if (TString(name) == "Cu2")     {fN = 63;  fR = 4.20;       fA = 0.596;  fW =  0;       fF = 8;  fZ=29; fBeta2=0.162; fBeta4=-0.006;}  
  else if (TString(name) == "Cu2rw")   {fN = 63;  fR = 4.20;       fA = 0.596;  fW =  0;       fF = 14; fZ=29; fBeta2=0.162; fBeta4=-0.006; r0=1.01269; r1=-0.00298083; r2=-9.97222e-05;}  
  else if (TString(name) == "CuHN")    {fN = 63;  fR = 4.28;       fA = 0.5;    fW =  0;       fF = 1;  fZ=29;} // from arXiv:0904.4080v1
  else if (TString(name) == "Br")      {fN = 79;  fR = 4.1629;     fA = 0.56;   fW =  0;       fF = 1;  fZ=35;} // from the ceiling
  else if (TString(name) == "Ag")      {fN = 109; fR = 4.5638;     fA = 0.53;   fW =  0;       fF = 1;  fZ=47;} // from the ceiling
  else if (TString(name) == "Xe")      {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 1;  fZ=54;} // adapted from arXiv:1703.04278
  else if (TString(name) == "Xes")     {fN = 129; fR = 5.42;       fA = 0.57;   fW =  0;       fF = 1;  fZ=54;} // scale from Sb (Antimony, A=122, r=5.32) by 1.019 = (129/122)**0.333
  else if (TString(name) == "Xe2")     {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 8;  fZ=54; fBeta2=0.161; fBeta4=-0.003;} // adapted from arXiv:1703.04278 and Z. Physik (1974) 270: 113
  else if (TString(name) == "Xe2a")    {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 8;  fZ=54; fBeta2=0.18; fBeta4=0;} // ALICE parameters (see public note from 2018 at https://cds.cern.ch/collection/ALICE%20Public%20Notes?ln=en)
  else if (TString(name) == "Xerw")    {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 12; fZ=54; r0=1.00911; r1=-0.000722999; r2=-0.0002663;}
  else if (TString(name) == "Xesrw")   {fN = 129; fR = 5.42;       fA = 0.57;   fW =  0;       fF = 12; fZ=54; r0=1.0096; r1=-0.000874123; r2=-0.000256708;}
  else if (TString(name) == "Xe2arw")  {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 14; fZ=54; fBeta2=0.18; fBeta4=0; r0=1.01246; r1=-0.0024851; r2=-5.72464e-05;} 
  else if (TString(name) == "W")       {fN = 186; fR = 6.58;       fA = 0.480;  fW =  0;       fF = 1;  fZ=74;}
  else if (TString(name) == "Au")      {fN = 197; fR = 6.38;       fA = 0.535;  fW =  0;       fF = 1;  fZ=79;}
  else if (TString(name) == "Aurw")    {fN = 197; fR = 6.38;       fA = 0.535;  fW =  0;       fF = 12; fZ=79; r0=1.00899; r1=-0.000590908; r2=-0.000210598;}
  else if (TString(name) == "Au2")     {fN = 197; fR = 6.38;       fA = 0.535;  fW =  0;       fF = 8;  fZ=79; fBeta2=-0.131; fBeta4=-0.031; }
  else if (TString(name) == "Au2rw")   {fN = 197; fR = 6.38;       fA = 0.535;  fW =  0;       fF = 14; fZ=79; fBeta2=-0.131; fBeta4=-0.031; r0=1.01261; r1=-0.00225517; r2=-3.71513e-05;}
  else if (TString(name) == "AuHN")    {fN = 197; fR = 6.42;       fA = 0.44;   fW =  0;       fF = 1;  fZ=79;} // from arXiv:0904.4080v1
  else if (TString(name) == "Pb")      {fN = 208; fR = 6.62;       fA = 0.546;  fW =  0;       fF = 1;  fZ=82;}
  else if (TString(name) == "Pbrw")    {fN = 208; fR = 6.62;       fA = 0.546;  fW =  0;       fF = 12; fZ=82; r0=1.00863; r1=-0.00044808; r2=-0.000205872;} //only Pb 207 was tested but should be the same for 208
  else if (TString(name) == "Pb*")     {fN = 208; fR = 6.624;      fA = 0.549;  fW =  0;       fF = 1;  fZ=82;}
  else if (TString(name) == "PbHN")    {fN = 208; fR = 6.65;       fA = 0.460;  fW =  0;       fF = 1;  fZ=82;}
  else if (TString(name) == "Pbpn")    {fN = 208; fR = 6.68;       fA = 0.447;  fW =  0;       fF = 11; fZ=82; fR2=6.69; fA2=0.56; fW2=0;}
  else if (TString(name) == "Pbpnrw")  {fN = 208; fR = 6.68;       fA = 0.447;  fW =  0;       fF = 13; fZ=82; fR2=6.69; fA2=0.56; fW2=0;}
  // Uranium description taken from Heinz & Kuhlman, nucl-th/0411054.  In this code, fR is defined as 6.8*0.91, fW=6.8*0.26
  else if (TString(name) == "U")       {fN = 238; fR = 6.188;      fA = 0.54;   fW =  1.77;    fF = 5;  fZ=92;}  
  else if (TString(name) == "U2")      {fN = 238; fR = 6.67;       fA = 0.44;   fW =  0;       fF = 8;  fZ=92; fBeta2=0.280; fBeta4=0.093;}
  else {
    cout << "Could not find nucleus " << name << endl;
    return;
  }

  switch (fF) {
    case 0: // Proton exp
      fFunc1 = new TF1(name,"x*x*exp(-x/[0])",0,5);
      fFunc1->SetParameter(0,fR);
      break;
    case 1: // 3pF
      fFunc1 = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW);
      break;
    case 2: // 3pG
      fFunc1 = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x**2-[0]**2)/[1]**2))",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW);
      break;
    case 3: // Hulthen (see nucl-ex/0603010)
    case 4: // same but constrain the neutron opposite to the proton event-by-event
      fFunc1 = new TF1(name,"x*x*([0]*[1]*([0]+[1]))/(2*pi*(pow([0]-[1],2)))*pow((exp(-[0]*x)-exp(-[1]*x))/x,2)",0,fMaxR);
      fFunc1->SetParameters(fR,fA);
      break;
    case 5: // Ellipsoid (Uranium)
      fFunc1 = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,fMaxR);
      fFunc1->SetParameters(fR,fA,0); // same as 3pF but setting W to zero
      break;
    case 6: // He3/H3
      fFunc1 = 0; // read in file instead
      break;
    case 7: // Deformed nuclei, box method
#ifndef HAVE_MATHMORE
      cerr << "Need libMathMore.so for deformed nuclei" << endl;
      gSystem->Exit(123);
#endif
     fFunc1 = 0; // no func: only need beta parameters and use uniform box distribution
      break;
    case 8: // Deformed nuclei, TF2 method
      fFunc3 = new TF2(name,"x*x*TMath::Sin(y)/(1+exp((x-[0]*(1+[2]*0.315*(3*pow(cos(y),2)-1.0)+[3]*0.105*(35*pow(cos(y),4)-30*pow(cos(y),2)+3)))/[1]))",0,fMaxR,0.0,TMath::Pi());
      fFunc3->SetNpx(120);
      fFunc3->SetNpy(120);
      fFunc3->SetParameters(fR,fA,fBeta2,fBeta4);
      break;
    case 9: // Proton gaus
      fFunc1 = new TF1(name,"x*x*exp(-x*x/[0]/[0]/2)",0,5);
      fFunc1->SetParameter(0,fR);
      break;
    case 10: // Proton dgaus
      fFunc1 = new TF1(name,"x*x*((1-[0])/[1]^3*exp(-x*x/[1]/[1])+[0]/(0.4*[1])^3*exp(-x*x/(0.4*[1])^2))",0,5);
      fFunc1->SetParameter(0,0.5);
      fFunc1->SetParameter(1,fR);
      break;
    case 11: // 3pF for proton and neutrons
      fFunc1 = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW);
      fFunc2 = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,fMaxR);
      fFunc2->SetParameters(fR2,fA2,fW2);
      break;
    case 12: // reweighted
      fFunc1 = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))/([3]+[4]*x+[5]*x^2)",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW,r0,r1,r2); 
      fRecenter=1;
      fSmax=0.1;
      break;
    case 13: // Pb for proton and neutrons reweighted
      fFunc1 = new TF1(Form("%s_prot",name),"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))/([3]+[4]*x+[5]*x^2)",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW,1.00866,-0.000461484,-0.000203571);
      fFunc2 = new TF1(Form("%s_neut",name),"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))/([3]+[4]*x+[5]*x^2)",0,fMaxR);
      fFunc2->SetParameters(fR2,fA2,fW2,1.00866,-0.000461484,-0.000203571);
      fRecenter=1;
      fSmax=0.1;
      break;
    case 14: // Deformed nuclei, TF2 method, reweighted
      fFunc3 = new TF2(name,"x*x*TMath::Sin(y)/(1+exp((x-[0]*(1+[2]*0.315*(3*pow(cos(y),2)-1.0)+[3]*0.105*(35*pow(cos(y),4)-30*pow(cos(y),2)+3)))/[1]))/([4]+[5]*x+[6]*x^2)",0,fMaxR,0.0,TMath::Pi());
      fFunc3->SetNpx(120);
      fFunc3->SetNpy(120);
      fFunc3->SetParameters(fR,fA,fBeta2,fBeta4,r0,r1,r2);
      fRecenter=1;
      fSmax=0.1;
      break;
    case 15: // harmonic oscillator model 
      fFunc1 = new TF1(name,"x^2*(1+[0]*(x^2/[1]^2))*exp(-x^2/[1]^2)",0,fMaxR);
      fFunc1->SetParameters(fR,fA);
      break;
    default:
      cerr << "Could not find function type " << fF << endl;
  }
  return;
}

void TGlauNucleus::SetA(Double_t ia, Double_t ia2)
{
  fA  = ia;
  fA2 = ia2;
  switch (fF) {
    case 1:  // 3pF
    case 12: // 3pF with pol2 normalization
    case 2:  // 3pG
    case 5:  // Ellipsoid (Uranium)
      fFunc1->SetParameter(1,fA);
      break;
    case 8:
      fFunc3->SetParameter(1,fA);
      break;
    case 11: //p&n
      fFunc1->SetParameter(1,fA);//proton
      fFunc2->SetParameter(1,fA2);//neutron
      break;
    default:
      cout << "Error: fA not needed for function " << fF <<endl;
  }
}

void TGlauNucleus::SetBeta(Double_t b2, Double_t b4) 
{
  fBeta2=b2; 
  fBeta4=b4;      
  if (fFunc3) {
    fFunc3->SetParameter(2,fBeta2);
    fFunc3->SetParameter(3,fBeta4);
  }
}

void TGlauNucleus::SetR(Double_t ir, Double_t ir2)
{
  fR  = ir;
  fR2 = ir2;
  switch (fF) {
    case 0:  // Proton exp
    case 9:  // Proton gaus
    case 1:  // 3pF
    case 12: // 3pF with pol2 normalization
    case 2:  // 3pG
    case 5:  // Ellipsoid (Uranium)
      fFunc1->SetParameter(0,fR);
      break;
    case 8:
      fFunc3->SetParameter(0,fR);
      break;
    case 10: // Proton
      fFunc1->SetParameter(1,fR);
      break;
    case 11: // p&n
      fFunc1->SetParameter(0,fR);//proton
      fFunc2->SetParameter(0,fR2);//neutron
      break;
    default:
      cout << "Error: fR not needed for function " << fF <<endl;
  }
}

void TGlauNucleus::SetW(Double_t iw)
{
  fW = iw;
  switch (fF) {
    case 1: // 3pF
    case 2: // 3pG
      fFunc1->SetParameter(2,fW);
      break;
    default:
      cout << "Error: fW not needed for function " << fF <<endl;
  }
}

Bool_t TGlauNucleus::TestMinDist(Int_t n, Double_t x, Double_t y, Double_t z) const
{
  if (fMinDist<=0)
    return kTRUE;
  const Double_t md2 = fMinDist*fMinDist; 
  for (Int_t j = 0; j<n; ++j) {
    TGlauNucleon *other=(TGlauNucleon*)fNucleons->At(j);
    Double_t xo=other->GetX();
    Double_t yo=other->GetY();
    Double_t zo=other->GetZ();
    Double_t dist2 = (x-xo)*(x-xo)+
                     (y-yo)*(y-yo)+
                     (z-zo)*(z-zo);
    if (dist2<md2) {
      return kFALSE;
    }
  }
  return kTRUE;
}

TVector3 &TGlauNucleus::ThrowNucleons(Double_t xshift)
{

  if (fNucleons==0) {
    fNucleons=new TObjArray(fN);
    fNucleons->SetOwner();
    for (Int_t i=0; i<fN; ++i) {
      TGlauNucleon *nucleon=new TGlauNucleon(); 
      nucleon->SetType(0);
      if (i<fZ) 
        nucleon->SetType(1);
      fNucleons->Add(nucleon); 
    }
  } 
  if (1) { //randomize p and n in nucleus
    for (Int_t i=0,iz=0; i<fN; ++i) {
      TGlauNucleon *nucleon=(TGlauNucleon*)fNucleons->At(i);
      Double_t frac=double(fZ-iz)/double(fN-i);
      Double_t rn=gRandom->Uniform(0,1);
      if (rn<frac) {
        nucleon->SetType(1);
        ++iz;
      } else {
        nucleon->SetType(0);
      }
    }
  }

 cmscheck: /* start over here in case shift was too large */

  fTrials = 0;
  fNonSmeared = 0;
  fPhiRot = gRandom->Rndm()*2*TMath::Pi();
  const Double_t cosThetaRot = 2*gRandom->Rndm()-1;
  fThetaRot = TMath::ACos(cosThetaRot);
  fXRot = gRandom->Rndm()*2*TMath::Pi();
  fYRot = gRandom->Rndm()*2*TMath::Pi();
  fZRot = gRandom->Rndm()*2*TMath::Pi();

  const Bool_t hulthen = (fF==3||fF==4);
  TString tmpname(GetName());
  Bool_t nucleonsfromfile = false;
  if ((tmpname=="He3") || (tmpname=="H3") ||
      (tmpname=="He4") || (tmpname=="C")   || 
      (tmpname=="O") || (tmpname=="O2") || (tmpname=="Ca2")){nucleonsfromfile = true;}
  
  if (fN==1) { //special treatment for proton
    Double_t r = fFunc1->GetRandom();
    Double_t phi = gRandom->Rndm() * 2 * TMath::Pi();
    Double_t ctheta = 2*gRandom->Rndm() - 1;
    Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);
    TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(0));
    nucleon->Reset();
    nucleon->SetXYZ(r * stheta * TMath::Cos(phi),
        r * stheta * TMath::Sin(phi),
        r * ctheta);
    fTrials = 1;

  } else if (fN==2 && hulthen) { //special treatment for Hulten

    Double_t r = fFunc1->GetRandom()/2;
    Double_t phi = gRandom->Rndm() * 2 * TMath::Pi();
    Double_t ctheta = 2*gRandom->Rndm() - 1;
    Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);

    TGlauNucleon *nucleon1=(TGlauNucleon*)(fNucleons->At(0));
    TGlauNucleon *nucleon2=(TGlauNucleon*)(fNucleons->At(1));
    nucleon1->Reset();
    nucleon1->SetXYZ(r * stheta * TMath::Cos(phi),
                     r * stheta * TMath::Sin(phi),
                     r * ctheta);
    nucleon2->Reset();
    if (fF==4) { // place opposite of 1
      nucleon2->SetXYZ(-nucleon1->GetX(),
           -nucleon1->GetY(),
           -nucleon1->GetZ());
    } else {
      r = fFunc1->GetRandom()/2;
      phi = gRandom->Rndm() * 2 * TMath::Pi();
      ctheta = 2*gRandom->Rndm() - 1;
      stheta = TMath::Sqrt(1-ctheta*ctheta);
      nucleon2->SetXYZ(r * stheta * TMath::Cos(phi),
           r * stheta * TMath::Sin(phi),
           r * ctheta);
    }
    fTrials = 1;

  } else if (fN > 2 && fN < 41 && nucleonsfromfile) {

    // if the first call, then read in the file configurations
    if (fNucCounter == -1) {
        std::string filepath(__FILE__);
        std::string locfilename(basename(__FILE__));
        filepath.erase(filepath.length() - locfilename.length(), locfilename.length());
        // read in the ascii file into the array and step through the counter
      char filename[100] = "foo.dat";
      if (tmpname=="He3") {
        filepath += "he3_plaintext.dat";
        //sprintf(filename,filepath.c_str());
      } else if (tmpname=="H3") {
        filepath += "h3_plaintext.dat";
        //sprintf(filename,filepath.c_str());
      } else if (tmpname=="He4") {
        filepath += "he4_plaintext.dat";
        //sprintf(filename,"../TGlauber/he4_plaintext.dat");
      } else if (tmpname=="C") {
        filepath += "carbon_plaintext.dat";
        //sprintf(filename,"../TGlauber/carbon_plaintext.dat");
      } else if (tmpname=="O") {
          filepath += "oxygen_plaintext.dat";
        //sprintf(filename,"../TGlauber/oxygen_plaintext.dat");
      } else if (tmpname=="O2") {
        filepath += "o16_alv.dat";
        //sprintf(filename,"../TGlauber/o16_alv.dat");
      } else if (tmpname=="Ca2") {
        filepath += "ca40_alv.dat";
        //sprintf(filename,"../TGlauber/ca40_alv.dat");
      }
      sprintf(filename,filepath.c_str());
      cout << "Reading in " << filename << " for nucleon configurations with fN = " << fN << endl;
      ifstream myfile;
      myfile.open(filename);
      if (!myfile) {
        cout << "ERROR:  no file for nucleon configurations found with name = " << filename << endl;
        gSystem->Exit(123);
      }

      Int_t inputcounter = 0;
      while (myfile) {
        //if (inputcounter > 5999) break;
          if (fNucCounter > 9999 && (tmpname=="O2" || tmpname=="Ca2") ) {break;}
          else if(fNucCounter > 5999) {break;}
        Double_t foo;
  if (fN == 3) {
    myfile >> fNucArr[inputcounter][0][0] >> fNucArr[inputcounter][0][1] >> fNucArr[inputcounter][0][2]
     >> fNucArr[inputcounter][1][0] >> fNucArr[inputcounter][1][1] >> fNucArr[inputcounter][1][2]
     >> fNucArr[inputcounter][2][0] >> fNucArr[inputcounter][2][1] >> fNucArr[inputcounter][2][2]
     >> foo >> foo >> foo >> foo;
  } else if (fN == 4) {
    // no extra data with isospin information at the end of the nucleon configurations
    myfile >> fNucArr[inputcounter][0][0] >> fNucArr[inputcounter][0][1] >> fNucArr[inputcounter][0][2]
     >> fNucArr[inputcounter][1][0] >> fNucArr[inputcounter][1][1] >> fNucArr[inputcounter][1][2]
     >> fNucArr[inputcounter][2][0] >> fNucArr[inputcounter][2][1] >> fNucArr[inputcounter][2][2]
     >> fNucArr[inputcounter][3][0] >> fNucArr[inputcounter][3][1] >> fNucArr[inputcounter][3][2];
  } else if (fN == 12) {
    // no extra data with isospin information at the end of the nucleon configurations
    // two extra words at the beginning --> foo foo
    myfile >> foo >> foo 
           >> fNucArr[inputcounter][0][0] >> fNucArr[inputcounter][0][1] >> fNucArr[inputcounter][0][2]
     >> fNucArr[inputcounter][1][0] >> fNucArr[inputcounter][1][1] >> fNucArr[inputcounter][1][2]
     >> fNucArr[inputcounter][2][0] >> fNucArr[inputcounter][2][1] >> fNucArr[inputcounter][2][2]
     >> fNucArr[inputcounter][3][0] >> fNucArr[inputcounter][3][1] >> fNucArr[inputcounter][3][2]
     >> fNucArr[inputcounter][4][0] >> fNucArr[inputcounter][4][1] >> fNucArr[inputcounter][4][2]
     >> fNucArr[inputcounter][5][0] >> fNucArr[inputcounter][5][1] >> fNucArr[inputcounter][5][2]
     >> fNucArr[inputcounter][6][0] >> fNucArr[inputcounter][6][1] >> fNucArr[inputcounter][6][2]
     >> fNucArr[inputcounter][7][0] >> fNucArr[inputcounter][7][1] >> fNucArr[inputcounter][7][2]
     >> fNucArr[inputcounter][8][0] >> fNucArr[inputcounter][8][1] >> fNucArr[inputcounter][8][2]
     >> fNucArr[inputcounter][9][0] >> fNucArr[inputcounter][9][1] >> fNucArr[inputcounter][9][2]
     >> fNucArr[inputcounter][10][0] >> fNucArr[inputcounter][10][1] >> fNucArr[inputcounter][10][2]
     >> fNucArr[inputcounter][11][0] >> fNucArr[inputcounter][11][1] >> fNucArr[inputcounter][11][2];
  } else if (fN == 16 && tmpname=="O") {
    // no extra data with isospin information at the end of the nucleon configurations
    myfile >> fNucArr[inputcounter][0][0] >> fNucArr[inputcounter][0][1] >> fNucArr[inputcounter][0][2]
     >> fNucArr[inputcounter][1][0] >> fNucArr[inputcounter][1][1] >> fNucArr[inputcounter][1][2]
     >> fNucArr[inputcounter][2][0] >> fNucArr[inputcounter][2][1] >> fNucArr[inputcounter][2][2]
     >> fNucArr[inputcounter][3][0] >> fNucArr[inputcounter][3][1] >> fNucArr[inputcounter][3][2]
     >> fNucArr[inputcounter][4][0] >> fNucArr[inputcounter][4][1] >> fNucArr[inputcounter][4][2]
     >> fNucArr[inputcounter][5][0] >> fNucArr[inputcounter][5][1] >> fNucArr[inputcounter][5][2]
     >> fNucArr[inputcounter][6][0] >> fNucArr[inputcounter][6][1] >> fNucArr[inputcounter][6][2]
     >> fNucArr[inputcounter][7][0] >> fNucArr[inputcounter][7][1] >> fNucArr[inputcounter][7][2]
     >> fNucArr[inputcounter][8][0] >> fNucArr[inputcounter][8][1] >> fNucArr[inputcounter][8][2]
     >> fNucArr[inputcounter][9][0] >> fNucArr[inputcounter][9][1] >> fNucArr[inputcounter][9][2]
     >> fNucArr[inputcounter][10][0] >> fNucArr[inputcounter][10][1] >> fNucArr[inputcounter][10][2]
     >> fNucArr[inputcounter][11][0] >> fNucArr[inputcounter][11][1] >> fNucArr[inputcounter][11][2]
     >> fNucArr[inputcounter][12][0] >> fNucArr[inputcounter][12][1] >> fNucArr[inputcounter][12][2]
     >> fNucArr[inputcounter][13][0] >> fNucArr[inputcounter][13][1] >> fNucArr[inputcounter][13][2]
     >> fNucArr[inputcounter][14][0] >> fNucArr[inputcounter][14][1] >> fNucArr[inputcounter][14][2]
     >> fNucArr[inputcounter][15][0] >> fNucArr[inputcounter][15][1] >> fNucArr[inputcounter][15][2];
  }
  else if (fN == 16 && tmpname=="O2") {
      // read from file by M. Alvioli et al. Phys. Lett. B680 (2009) 225
      myfile >> fNucArrAlv[inputcounter][0][0] >> fNucArrAlv[inputcounter][0][1] >> fNucArrAlv[inputcounter][0][2] >> fNucArrAlv[inputcounter][0][3];
      myfile >> fNucArrAlv[inputcounter][1][0] >> fNucArrAlv[inputcounter][1][1] >> fNucArrAlv[inputcounter][1][2] >> fNucArrAlv[inputcounter][1][3];
      myfile >> fNucArrAlv[inputcounter][2][0] >> fNucArrAlv[inputcounter][2][1] >> fNucArrAlv[inputcounter][2][2] >> fNucArrAlv[inputcounter][2][3];
      myfile >> fNucArrAlv[inputcounter][3][0] >> fNucArrAlv[inputcounter][3][1] >> fNucArrAlv[inputcounter][3][2] >> fNucArrAlv[inputcounter][3][3];
      myfile >> fNucArrAlv[inputcounter][4][0] >> fNucArrAlv[inputcounter][4][1] >> fNucArrAlv[inputcounter][4][2] >> fNucArrAlv[inputcounter][4][3];
      myfile >> fNucArrAlv[inputcounter][5][0] >> fNucArrAlv[inputcounter][5][1] >> fNucArrAlv[inputcounter][5][2] >> fNucArrAlv[inputcounter][5][3];
      myfile >> fNucArrAlv[inputcounter][6][0] >> fNucArrAlv[inputcounter][6][1] >> fNucArrAlv[inputcounter][6][2] >> fNucArrAlv[inputcounter][6][3];
      myfile >> fNucArrAlv[inputcounter][7][0] >> fNucArrAlv[inputcounter][7][1] >> fNucArrAlv[inputcounter][7][2] >> fNucArrAlv[inputcounter][7][3];
      myfile >> fNucArrAlv[inputcounter][8][0] >> fNucArrAlv[inputcounter][8][1] >> fNucArrAlv[inputcounter][8][2] >> fNucArrAlv[inputcounter][8][3];
      myfile >> fNucArrAlv[inputcounter][9][0] >> fNucArrAlv[inputcounter][9][1] >> fNucArrAlv[inputcounter][9][2] >> fNucArrAlv[inputcounter][9][3];
      myfile >> fNucArrAlv[inputcounter][10][0] >> fNucArrAlv[inputcounter][10][1] >> fNucArrAlv[inputcounter][10][2] >> fNucArrAlv[inputcounter][10][3];
      myfile >> fNucArrAlv[inputcounter][11][0] >> fNucArrAlv[inputcounter][11][1] >> fNucArrAlv[inputcounter][11][2] >> fNucArrAlv[inputcounter][11][3];
      myfile >> fNucArrAlv[inputcounter][12][0] >> fNucArrAlv[inputcounter][12][1] >> fNucArrAlv[inputcounter][12][2] >> fNucArrAlv[inputcounter][12][3];
      myfile >> fNucArrAlv[inputcounter][13][0] >> fNucArrAlv[inputcounter][13][1] >> fNucArrAlv[inputcounter][13][2] >> fNucArrAlv[inputcounter][13][3];
      myfile >> fNucArrAlv[inputcounter][14][0] >> fNucArrAlv[inputcounter][14][1] >> fNucArrAlv[inputcounter][14][2] >> fNucArrAlv[inputcounter][14][3];
      myfile >> fNucArrAlv[inputcounter][15][0] >> fNucArrAlv[inputcounter][15][1] >> fNucArrAlv[inputcounter][15][2] >> fNucArrAlv[inputcounter][15][3];
  }
  else if (fN == 40 && tmpname=="Ca2") {
      // read from file by M. Alvioli et al. Phys. Lett. B680 (2009) 225
    for(int nn = 0; nn < 40; nn++){myfile >> fNucArrAlv[inputcounter][nn][0] >> fNucArrAlv[inputcounter][nn][1] >> fNucArrAlv[inputcounter][nn][2] >> fNucArrAlv[inputcounter][nn][3];}
  }

        ++inputcounter;
      }
      myfile.close();
      fNucCounter=0;
    } // done reading in the file the first time

    if (fNucCounter > 9999 && (tmpname=="O2" || tmpname=="Ca2")) {fNucCounter = 0;}
    else if(fNucCounter > 5999) {fNucCounter = 0;}

    // change to loop over fN nucleons!
    for (Int_t i = 0; i<fN; ++i) {
      TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
      nucleon->Reset();
      if ((tmpname=="O2" )||(tmpname=="Ca2" )){
      nucleon->SetXYZ(fNucArrAlv[fNucCounter][i][0],
          fNucArrAlv[fNucCounter][i][1],
          fNucArrAlv[fNucCounter][i][2]);
      nucleon->SetType(fNucArrAlv[fNucCounter][i][3]);
      nucleon->RotateXYZ(fPhiRot,fThetaRot);}
      else{
          nucleon->SetXYZ(fNucArr[fNucCounter][i][0],
                          fNucArr[fNucCounter][i][1],
                          fNucArr[fNucCounter][i][2]);
          nucleon->RotateXYZ(fPhiRot,fThetaRot);
      }
    }

    ++fNucCounter;
    fTrials = 1;
  } else { // all other nuclei 

    const Double_t startingEdge  = 20; // throw nucleons within a cube of this size (fm)
    const Double_t startingEdgeX = startingEdge + fNodeDist*gRandom->Rndm() - 0.5*fNodeDist;
    const Double_t startingEdgeY = startingEdge + fNodeDist*gRandom->Rndm() - 0.5*fNodeDist;
    const Double_t startingEdgeZ = startingEdge + fNodeDist*gRandom->Rndm() - 0.5*fNodeDist;
    const Int_t nslots = 2*startingEdge/fNodeDist+1;
    if (fNodeDist>0) {
      if (fMinDist>fNodeDist) {
        cout << "Minimum distance (nucleon hard core diameter) [" 
          << fMinDist << "] cannot be larger than the nodal spacing of the grid [" 
          << fNodeDist << "]." << endl;
        cout << "Quitting...." << endl;
        gSystem->Exit(123);
      }
      if (!fIsUsed)
        fIsUsed = new TBits(nslots*nslots*nslots);
      else
        fIsUsed->ResetAllBits();
    }
    for (Int_t i = 0; i<fN; ++i) {
      TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
      nucleon->Reset();
      while (1) {
        ++fTrials;
        Bool_t nucleon_inside = 0;
        Double_t x=999, xsmeared=999;
        Double_t y=999, ysmeared=999;
        Double_t z=999, zsmeared=999;
        if (fF==5||fF==7) { // the extended way, throw in a box and test the weight
          while (!nucleon_inside) {
            x = (fR*2)*(gRandom->Rndm() * 2 - 1);
            y = (fR*2)*(gRandom->Rndm() * 2 - 1);
            z = (fR*2)*(gRandom->Rndm() * 2 - 1);
            Double_t r = TMath::Sqrt(x*x+y*y);
            Double_t theta = TMath::ATan2(r,z);
            Double_t R = TMath::Sqrt(x*x+y*y+z*z);
            Double_t Rtheta = fR;
            if (fF==5)
              Rtheta= fR + fW*TMath::Cos(theta)*TMath::Cos(theta);
            if (fF==7)
#ifdef HAVE_MATHMORE
              Rtheta = fR*(1+fBeta2*ROOT::Math::sph_legendre(2,0,theta)+fBeta4*ROOT::Math::sph_legendre(4,0,theta));
#else
            cerr << "Should not end here because you do not have libMathMore" << endl;
#endif
            Double_t prob = 1/(1+TMath::Exp((R-Rtheta)/fA));
            if (gRandom->Rndm()<prob) 
              nucleon_inside=1;
          }
        } else if ((fF==8) || (fF==14)) { // use TF2
          Double_t r;
          Double_t theta;
          fFunc3->GetRandom2(r,theta);
          Double_t phi = 2*TMath::Pi()*gRandom->Rndm();
          x = r * TMath::Sin(phi) * TMath::Sin(theta);
          y = r * TMath::Cos(phi) * TMath::Sin(theta);
          z = r *                   TMath::Cos(theta);
        } else { // all other types
          TF1 *ff = fFunc1;
          if ((fFunc2) && (nucleon->GetType()==0))
            ff = fFunc2;
          if (fNodeDist<=0) { // "continuous" mode
            Double_t r = ff->GetRandom();
            Double_t phi = 2*TMath::Pi()*gRandom->Rndm();
            Double_t ctheta = 2*gRandom->Rndm() - 1 ;
            Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);
            x = r * stheta * TMath::Cos(phi);
            y = r * stheta * TMath::Sin(phi);
            z = r * ctheta;
          } else { // "grid/lattice" mode
            Int_t iNode = Int_t((2*startingEdge/fNodeDist)*gRandom->Rndm());
            Int_t jNode = Int_t((2*startingEdge/fNodeDist)*gRandom->Rndm());
            Int_t kNode = Int_t((2*startingEdge/fNodeDist)*gRandom->Rndm());
            Int_t index=iNode*nslots*nslots+jNode*nslots+kNode;
            if (fIsUsed->TestBitNumber(index))
              continue;
            if (fLattice==1) {       // Primitive cubic system (PCS) -> https://en.wikipedia.org/wiki/Cubic_crystal_system
              x = fNodeDist*(iNode) - startingEdgeX;
              y = fNodeDist*(jNode) - startingEdgeY;
              z = fNodeDist*(kNode) - startingEdgeZ;
            } else if (fLattice==2) { //Body centered cubic (BCC) -> http://mathworld.wolfram.com/CubicClosePacking.html
              x = 0.5*fNodeDist*(-iNode+jNode+kNode) - 0.5*startingEdgeX;
              y = 0.5*fNodeDist*(+iNode-jNode+kNode) - 0.5*startingEdgeY;
              z = 0.5*fNodeDist*(+iNode+jNode-kNode) - 0.5*startingEdgeZ;
            } else if (fLattice==3) { //Face Centered Cubic (FCC) -> http://mathworld.wolfram.com/CubicClosePacking.html
              x = 0.5*fNodeDist*(jNode+kNode) - startingEdgeX;
              y = 0.5*fNodeDist*(iNode+kNode) - startingEdgeY;
              z = 0.5*fNodeDist*(iNode+jNode) - startingEdgeZ;
            } else {                  //Hexagonal close packing (HCP) -> https://en.wikipedia.org/wiki/Close-packing_of_equal_spheres
              x = 0.5*fNodeDist*(2*iNode+((jNode+kNode)%2))          - startingEdgeX;
              y = 0.5*fNodeDist*(TMath::Sqrt(3)*(jNode+(kNode%2)/3)) - startingEdgeY;
              z = 0.5*fNodeDist*(kNode*2*TMath::Sqrt(6)/3)           - startingEdgeZ;
            }
            const Double_t r2 = x*x + y*y + z*z;
            const Double_t r  = TMath::Sqrt(r2);
      if ((r>fMaxR)||(r2*gRandom->Rndm()>ff->Eval(r)))
        continue;
            if (fSmearing>0.0) {
              Int_t nAttemptsToSmear = 0;
              while (1) {
                xsmeared = x*gRandom->Gaus(1.0,fSmearing);
                ysmeared = y*gRandom->Gaus(1.0,fSmearing);
                zsmeared = z*gRandom->Gaus(1.0,fSmearing);
                nAttemptsToSmear++;
                if (TestMinDist(i,xsmeared,ysmeared,zsmeared)) {
                  x = xsmeared;
                  y = ysmeared;
                  z = zsmeared;
                  break;
                }
                if (nAttemptsToSmear>=99) {
                  cerr << "Could not place on this node :: [" << x <<","<< y <<","<< z <<"] r = " << TMath::Sqrt(x*x+y*y+z*z) << " fm; "
                    << "Node (" << iNode << "," << jNode << "," << kNode << ") not smeared !!!" << endl;
                  ++fNonSmeared;
                  break;
                }
              }
            }
            fIsUsed->SetBitNumber(index);
          } /* end "grid/lattice mode" */
  }
  nucleon->SetXYZ(x,y,z);
  if (fF==5||fF==7||fF==8||fF==14) 
    nucleon->RotateXYZ(fPhiRot,fThetaRot); // Uranium etc.
  if (fNodeDist>0) {
    nucleon->RotateXYZ_3D(fXRot,fYRot,fZRot);
    break;
  }
  if (TestMinDist(i,x,y,z))
    break;
      }
    }
  }    

  // calculate center of mass
  Double_t sumx=0;       
  Double_t sumy=0;       
  Double_t sumz=0;       
  for (Int_t i = 0; i<fN; ++i) {
    TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
    sumx += nucleon->GetX();
    sumy += nucleon->GetY();
    sumz += nucleon->GetZ();
  }
  sumx = sumx/fN;
  sumy = sumy/fN;
  sumz = sumz/fN;

  static TVector3 finalShift;
  finalShift.SetXYZ(sumx,sumy,sumz);
  if (finalShift.Mag()>fSmax)
    goto cmscheck;
  Double_t fsumx = 0;
  Double_t fsumy = 0;
  Double_t fsumz = 0;
  if (fRecenter==1) {
    fsumx = sumx;
    fsumy = sumy;
    fsumz = sumz;
  } else if (fRecenter==2) {
    TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(fN-1));
    Double_t x = nucleon->GetX() - fN*sumx;
    Double_t y = nucleon->GetY() - fN*sumy;
    Double_t z = nucleon->GetZ() - fN*sumz;
    nucleon->SetXYZ(x,y,z);
  } else if ((fRecenter==3)||(fRecenter==4)) {
    TVector3 zVec;
    zVec.SetXYZ(0,0,1);
    TVector3 shiftVec;
    shiftVec.SetXYZ(sumx,sumy,sumz);
    TVector3 orthVec;
    orthVec = shiftVec.Cross(zVec);
    TRotation myRot;
    myRot.Rotate(shiftVec.Angle(zVec),orthVec);
    TVector3 myNuc;
    for (Int_t i = 0; i<fN; ++i) {
      TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
      myNuc.SetXYZ(nucleon->GetX(),nucleon->GetY(),nucleon->GetZ());
      myNuc.Transform(myRot);
      nucleon->SetXYZ(myNuc.X(), myNuc.Y(), myNuc.Z());
    }
    if (fRecenter==3)
      fsumz = shiftVec.Mag();
  }

  // recenter and shift
  for (Int_t i = 0; i<fN; ++i) {
    TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
    nucleon->SetXYZ(nucleon->GetX()-fsumx + xshift,
        nucleon->GetY()-fsumy,
        nucleon->GetZ()-fsumz);
  }

  return finalShift;

}
