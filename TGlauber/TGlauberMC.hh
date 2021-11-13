#ifndef TGlauberMC_h
#define TGlauberMC_h 1

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
#include <TRandom1.h>
#include <TRotation.h>
#include <TString.h>
#include <TSystem.h>
#include <TVector3.h>
#include "TGlauNucleon.hh"
#include "TGlauNucleus.hh"
#ifdef HAVE_MATHMORE
 #include <Math/SpecFuncMathMore.h>
#endif
using namespace std;
#endif

#if !defined(__CINT__) || defined(__MAKECINT__)
#define _runglauber_ 3
#endif

//---------------------------------------------------------------------------------
TF1 *getNNProf(Double_t snn=67.6, Double_t omega=0.4, Double_t G=1);

//---------------------------------------------------------------------------------
class TGlauberMC : public TNamed
{
  public:
    class Event {
      public:

        Float_t Npart;       //Number of wounded (participating) nucleons in current event
        Float_t Ncoll;       //Number of binary collisions in current event
        Float_t Nhard;       //Number of hard collisions in current event (based on fHardFrac)
        Float_t Nvoid;       //Number of collisions without nucleon-nucleon interaction
        Float_t B;           //[0,0,16] Impact parameter (b)
        Float_t BNN;         //[0,0,16] Average NN impact parameter
        Float_t Ncollpp;     //Ncoll pp
        Float_t Ncollpn;     //Ncoll pn
        Float_t Ncollnn;     //Ncoll nn
        Float_t VarX;        //[0,0,16] Variance of x of wounded nucleons
        Float_t VarY;        //[0,0,16] Variance of y of wounded nucleons
        Float_t VarXY;       //[0,0,16] Covariance of x and y of wounded nucleons
        Float_t NpartA;      //Number of wounded (participating) nucleons in Nucleus A
        Float_t NpartB;      //Number of wounded (participating) nucleons in Nucleus B
        Float_t Npart0;      //Number of singly-wounded (participating) nucleons
        Float_t AreaW;       //[0,0,16] area defined by width of participants
        Float_t Psi1;        //[0,0,16] psi1
        Float_t Ecc1;        //[0,0,16] eps1
        Float_t Psi2;        //[0,0,16] psi2
        Float_t Ecc2;        //[0,0,16] eps2
        Float_t Psi3;        //[0,0,16] psi3
        Float_t Ecc3;        //[0,0,16] eps3
        Float_t Psi4;        //[0,0,16] psi4
        Float_t Ecc4;        //[0,0,16] eps4
        Float_t Psi5;        //[0,0,16] psi5
        Float_t Ecc5;        //[0,0,16] eps5
        Float_t AreaA;       //[0,0,16] area defined by "and" of participants
        Float_t AreaO;       //[0,0,16] area defined by "or" of participants
        Float_t X0;          //[0,0,16] production point in x
        Float_t Y0;          //[0,0,16] production point in y
        Float_t Phi0;        //[0,0,16] direction in phi
        Float_t Length;      //[0,0,16] length in phi0
        Float_t MeanX;       //[0,0,16] <x> of wounded nucleons
        Float_t MeanY;       //[0,0,16] <y> of wounded nucleons
        Float_t MeanX2;      //[0,0,16] <x^2> of wounded nucleons
        Float_t MeanY2;      //[0,0,16] <y^2> of wounded nucleons
        Float_t MeanXY;      //[0,0,16] <xy> of wounded nucleons
        Float_t MeanXSystem; //[0,0,16] <x> of all nucleons
        Float_t MeanYSystem; //[0,0,16] <y> of all nucleons  
        Float_t MeanXA;      //[0,0,16] <x> of nucleons in nucleus A
        Float_t MeanYA;      //[0,0,16] <y> of nucleons in nucleus A
        Float_t MeanXB;      //[0,0,16] <x> of nucleons in nucleus B
        Float_t MeanYB;      //[0,0,16] <y> of nucleons in nucleus B
        Float_t PhiA;        //[0,0,16] phi angle nucleus A
        Float_t ThetaA;      //[0,0,16] theta angle nucleus B
        Float_t PhiB;        //[0,0,16] phi angle nucleus B
        Float_t ThetaB;      //[0,0,16] theta angle nucleus B
        void    Reset()      {Npart=0;Ncoll=0;Nhard=0;B=0;BNN=0;Ncollpp=0;Ncollpn=0;Ncollnn=0;VarX=0;VarY=0;VarXY=0;NpartA=0;NpartB=0;Npart0=0;AreaW=0;
                              Psi1=0;Ecc1=0;Psi2=0;Ecc2=0;Psi3=0;Ecc3=0;Psi4=0;Ecc4=0;Psi5=0;Ecc5=0;
                              AreaA=0;AreaO=0;X0=0;Y0=0;Phi0=0;Length=0;
                              MeanX=0;MeanY=0;MeanX2=0;MeanY2=0;MeanXY=0;MeanXSystem=0;MeanYSystem=0;MeanXA=0;MeanYA=0;MeanXB=0;MeanYB=0;
                              PhiA=0;ThetaA=0;PhiB=0;ThetaB=0;} // order must match that given in vars below
       // ClassDef(TGlauberMC::Event, 1)
    };

  protected:
    TGlauNucleus  fANucleus;       //Nucleus A
    TGlauNucleus  fBNucleus;       //Nucleus B
    Double_t      fXSect;          //Nucleon-nucleon cross section
    Double_t      fXSectOmega;     //StdDev of Nucleon-nucleon cross section
    Double_t      fXSectLambda;    //Jacobian from tot to inelastic (Strikman)
    Double_t      fXSectEvent;     //Event value of Nucleon-nucleon cross section
    TObjArray*    fNucleonsA;      //Array of nucleons in nucleus A
    TObjArray*    fNucleonsB;      //Array of nucleons in nucleus B
    TObjArray*    fNucleons;       //Array which joins Nucleus A & B
    Int_t         fAN;             //Number of nucleons in nucleus A
    Int_t         fBN;             //Number of nucleons in nucleus B
    TNtuple*      fNt;             //Ntuple for results (created, but not deleted)
    Int_t         fEvents;         //Number of events with at least one collision
    Int_t         fTotalEvents;    //All events within selected impact parameter range
    Double_t      fBmin;           //Minimum impact parameter to be generated
    Double_t      fBmax;           //Maximum impact parameter to be generated
    Double_t      fHardFrac;       //Fraction of cross section used for Nhard (def=0.65)
    Int_t         fDetail;         //Detail to store (99=all by default)
    Bool_t        fCalcArea;       //If true calculate overlap area via grid (slow, off by default)
    Bool_t        fCalcLength;     //If true calculate path length (slow, off by default)
    Bool_t        fDoCore;         //If true calculate area and eccentricy only for core participants (off by default)
    Int_t         fMaxNpartFound;  //Largest value of Npart obtained
    Double_t      fPsiN[10];       //Psi N
    Double_t      fEccN[10];       //Ecc N
    Double_t      f2Cx;            //Two-component x
    TF1          *fPTot;           //Cross section distribution
    TF1          *fNNProf;         //NN profile (hard-sphere == 0 by default)
    Event         fEv;             //Glauber event (results of calculation stored in tree)
    Bool_t        fBC[999][999];   //Array to record binary collision
    Bool_t        CalcResults(Double_t bgen);
    Bool_t        CalcEvent(Double_t bgen);

  public:
    TGlauberMC(const char* NA = "Pb", const char* NB = "Pb", Double_t xsect = 42, Double_t xsectsigma=0, ULong_t seed = 65539);
  virtual            ~TGlauberMC() {delete fNt; fNt=0;}

    Double_t            CalcDens(TF1 &prof, Double_t xval, Double_t yval) const;
    void                Draw(Option_t* option="");
    Double_t            GetB()                 const {return fEv.B;}
    Double_t            GetBNN()               const {return fEv.BNN;}
    Double_t            GetBmax()              const {return fBmax;}
    Double_t            GetBmin()              const {return fBmin;}
    Double_t            GetEcc(Int_t i=2)      const {return fEccN[i];}
    Double_t            GetHardFrac()          const {return fHardFrac;}
    Double_t            GetMeanX()             const {return fEv.MeanX;}
    Double_t            GetMeanXParts()        const {return fEv.MeanX;}
    Double_t            GetMeanXSystem()       const {return fEv.MeanXSystem;}
    Double_t            GetMeanY()             const {return fEv.MeanY;}
    Double_t            GetMeanYParts()        const {return fEv.MeanY;}
    Double_t            GetMeanYSystem()       const {return fEv.MeanYSystem;}
    Double_t            GetPsi(Int_t i=2)      const {return fPsiN[i];}
    Double_t            GetSx2()               const {return fEv.VarX;}    
    Double_t            GetSxy()               const {return fEv.VarXY;}    
    Double_t            GetSy2()               const {return fEv.VarY;}    
    Double_t            GetTotXSect()          const;
    Double_t            GetTotXSectErr()       const;
    Double_t            GetXSectEvent()        const {return fXSectEvent;}
    Int_t               GetNcoll()             const {return fEv.Ncoll;}
    Int_t               GetNcollnn()           const {return fEv.Ncollnn;}
    Int_t               GetNcollpn()           const {return fEv.Ncollpn;}
    Int_t               GetNcollpp()           const {return fEv.Ncollpp;}
    Int_t               GetNpart()             const {return fEv.Npart;}
    Int_t               GetNpart0()            const {return fEv.Npart0;}
    Int_t               GetNpartA()            const {return fEv.NpartA;}
    Int_t               GetNpartB()            const {return fEv.NpartB;}
    Int_t               GetNpartFound()        const {return fMaxNpartFound;}
    Int_t               GetNvoid()             const {return fEv.Nvoid;}
    Int_t               GetNhard()             const {return fEv.Nhard;}
    TF1*                GetXSectDist()         const {return fPTot;}
    TGlauNucleus*       GetNucleusA()                {return &fANucleus;}
    TGlauNucleus*       GetNucleusB()                {return &fBNucleus;}
    TNtuple*            GetNtuple()            const {return fNt;}
    TObjArray          *GetNucleons();
    const Event        &GetEvent()             const {return fEv;}
    const Event        *GetEvent()                   {return &fEv;}
    const TGlauNucleus* GetNucleusA()          const {return &fANucleus;}
    const TGlauNucleus* GetNucleusB()          const {return &fBNucleus;}
    Bool_t              IsBC(Int_t i, Int_t j) const {return fBC[i][j];}
    Bool_t              NextEvent(Double_t bgen=-1);
    void                Reset()                      {delete fNt; fNt=0; }
    Bool_t              ReadNextEvent(Bool_t calc=1, const char *fname=0);       
    void                Run(Int_t nevents,Double_t b=-1);
    void                Set2Cx(Double_t x)           {f2Cx = x;}
    void                SetBmax(Double_t bmax)       {fBmax = bmax;}
    void                SetBmin(Double_t bmin)       {fBmin = bmin;}
    void                SetCalcArea(Bool_t b)        {fCalcArea = b;}
    void                SetCalcCore(Bool_t b)        {fDoCore = b;}
    void                SetCalcLength(Bool_t b)      {fCalcLength = b;}
    void                SetDetail(Int_t d)           {fDetail = d;}
    void                SetHardFrac(Double_t f)      {fHardFrac=f;}
    void                SetLattice(Int_t i)          {fANucleus.SetLattice(i); fBNucleus.SetLattice(i);}
    void                SetMinDistance(Double_t d)   {fANucleus.SetMinDist(d); fBNucleus.SetMinDist(d);}
    void                SetNNProf(TF1 *f1)           {fNNProf = f1;}
    void                SetNodeDistance(Double_t d)  {fANucleus.SetNodeDist(d); fBNucleus.SetNodeDist(d);}
    void                SetRecenter(Int_t b)         {fANucleus.SetRecenter(b); fBNucleus.SetRecenter(b);}
    void                SetShiftMax(Double_t s)      {fANucleus.SetShiftMax(s); fBNucleus.SetShiftMax(s);}
    void                SetSmearing(Double_t s)      {fANucleus.SetSmearing(s); fBNucleus.SetSmearing(s);}
    const char         *Str()                  const {return Form("gmc-%s%s-snn%.1f-md%.1f-nd%.1f-rc%d-smax%.1f",fANucleus.GetName(),fBNucleus.GetName(),fXSect,fBNucleus.GetMinDist(),fBNucleus.GetNodeDist(),fBNucleus.GetRecenter(),fBNucleus.GetShiftMax());}
    static void         PrintVersion()               {cout << "TGlauberMC " << Version() << endl;}
    static const char  *Version()                    {return "v3.1a";}
#if !defined (__CINT__) || defined (__MAKECINT__)
      //ClassDef(TGlauberMC,6)
#endif
};

#endif
