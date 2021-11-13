#ifndef TGlauNucleus_h
#define TGlauNucleus_h 1

#define HAVE_MATHMORE

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TBits.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TH2.h>
#include <TLine.h>
#include <TMath.h>
#include <TNamed.h>
#include <TNtuple.h>
#include <TObjArray.h>
#include <TRandom.h>
#include <TRotation.h>
#include <TString.h>
#include <TSystem.h>
#include <TVector3.h>
#include <Math/SpecFuncMathMore.h>
using namespace std;
#endif

#if !defined(__CINT__) || defined(__MAKECINT__)
#define _runglauber_ 3
#endif


//---------------------------------------------------------------------------------
class TGlauNucleus : public TNamed
{
  private:
    Int_t      fN;                   //Number of nucleons
    Int_t      fZ;                   //Number of protons
    Double_t   fR;                   //Parameters of function
    Double_t   fA;                   //Parameters of function (fA+fZ=fN)
    Double_t   fW;                   //Parameters of function
    Double_t   fR2;                  //Parameters of function (for p and n separately)
    Double_t   fA2;                  //Parameters of function (for p and n separately)
    Double_t   fW2;                  //Parameters of function (for p and n separately)
    Double_t   fBeta2;               //Beta2 (deformed nuclei) 
    Double_t   fBeta4;               //Beta4 (deformed nuclei) 
    Double_t   fMinDist;             //Minimum separation distance
    Double_t   fNodeDist;            //Average node distance (set to <=0 if you do not want the "crystal lattice")
    Double_t   fSmearing;            //Node smearing (relevant if fNodeDist>0)
    Int_t      fRecenter;            //=1 by default (0=no recentering, 1=recenter all, 2=recenter displacing only one nucleon, 3=recenter by rotation)
    Int_t      fLattice;             //=0 use HCP by default (1=PCS, 2=BCC, 3=FCC)
    Double_t   fSmax;                //Maximum magnitude of cms shift tolerated (99, ie all by default) 
    Int_t      fF;                   //Type of radial distribution
    Int_t      fTrials;              //Store trials needed to complete nucleus
    Int_t      fNonSmeared;          //Store number of non-smeared-node nucleons
    TF1*       fFunc1;               //!Probability density function rho(r)
    TF1*       fFunc2;               //!Probability density function rho(r) -> if set 1 is for p, 2 is for n
    TF2*       fFunc3;               //!Probability density function rho(r,theta) for deformed nuclei
    TObjArray* fNucleons;            //!Array of nucleons
    Double_t   fPhiRot;              //!Angle phi for nucleus
    Double_t   fThetaRot;            //!Angle theta for nucleus
    Double_t   fXRot;                //!Angle around X axis for nucleus
    Double_t   fYRot;                //!Angle around Y axis for nucleus
    Double_t   fZRot;                //!Angle around Z axis for nucleus
    Double_t   fNucArr[6000][20][3]; //!Array of events (max 6000), up to 20 nucleons (only for small nuclei), 3 coordinates
    Double_t   fNucArrAlv[16000][40][4]; //!Array of events (max 16000), up to 40 nucleons (only for small nuclei), 3 coordinates + isospin
    Int_t      fNucCounter;          //!Event counter
    TBits     *fIsUsed;              //!Bits for lattice use  
    Double_t   fMaxR;                //!maximum radius (15fm)
    void       Lookup(const char* name);
    Bool_t     TestMinDist(Int_t n, Double_t x, Double_t y, Double_t z) const;

  public:
    TGlauNucleus(const char* iname="Pb", Int_t iN=0, Double_t iR=0, Double_t ia=0, Double_t iw=0, TF1* ifunc=0);
    virtual ~TGlauNucleus();
    using      TObject::Draw;
    void       Draw(Double_t xs, Int_t colp, Int_t cols);
    Double_t   GetA()             const {return fA;}
    TF1*       GetFunc1()         const {return GetFuncP();}
    TF1*       GetFunc2()         const {return GetFuncN();}
    TF2*       GetFunc3()         const {return GetFuncDef();}
    TF1*       GetFuncP()         const {return fFunc1;}
    TF1*       GetFuncN()         const {return fFunc2;}
    TF2*       GetFuncDef()       const {return fFunc3;}
    Double_t   GetMinDist()       const {return fMinDist;}
    Int_t      GetN()             const {return fN;}
    Double_t   GetNodeDist()      const {return fNodeDist;}
    TObjArray *GetNucleons()      const {return fNucleons;}
    Int_t      GetRecenter()      const {return fRecenter;}
    Double_t   GetR()             const {return fR;}
    Double_t   GetPhiRot()        const {return fPhiRot;}
    Double_t   GetThetaRot()      const {return fThetaRot;}
    Int_t      GetTrials()        const {return fTrials;}
    Int_t      GetNonSmeared()    const {return fNonSmeared;}
    Double_t   GetShiftMax()      const {return fSmax;}
    Double_t   GetW()             const {return fW;}
    Double_t   GetXRot()          const {return fXRot;}
    Double_t   GetYRot()          const {return fYRot;}
    Double_t   GetZRot()          const {return fZRot;}
    void       SetA(Double_t ia, Double_t ia2=-1);
    void       SetBeta(Double_t b2, Double_t b4); 
    void       SetLattice(Int_t i)               {fLattice=i;}
    void       SetMinDist(Double_t min)          {fMinDist=min;}
    void       SetN(Int_t in)                    {fN=in;}
    void       SetNodeDist(Double_t nd)          {fNodeDist=nd;}
    void       SetR(Double_t ir, Double_t ir2=-1);
    void       SetRecenter(Int_t b)              {fRecenter=b;}
    void       SetShiftMax(Double_t s)           {fSmax=s;}
    void       SetSmearing(Double_t s)           {fSmearing=s;}
    void       SetW(Double_t iw);
    TVector3  &ThrowNucleons(Double_t xshift=0.);
#if !defined (__CINT__) || defined (__MAKECINT__)
      //ClassDef(TGlauNucleus,6)
#endif   
};

#endif
