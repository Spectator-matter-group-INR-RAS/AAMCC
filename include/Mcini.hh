/* This is a collection of mcini data classes.
 * Taken from https://github.com/OlegGolosov/mcini to add mcini filetype support
 */

#ifndef GRATE_MCINI_HH
#define GRATE_MCINI_HH

#include "TNamed.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TObject.h"
#include "TParticle.h"

namespace mcini {

    typedef unsigned short IdType;

    enum eNucleonCollisionTypes {
        kNoCollision = 0,
        kElasticWithInitialNucleon,
        kElasticWithProducedParticle,
        kInelasticWithInitialNucleon,
        kInelasticWithProducedParticle,
        kNnucleonCollisionTypes
    };

    class Nucleon : public TObject {

    public:
        Nucleon() : id(0), pdgId(-1), momentum(0, 0, 0, 999), position(0, 0, 0, 999) {}
        Nucleon(int pdgId, const TLorentzVector &momentum, const TLorentzVector &position, unsigned short collisionType)
                : id(0), pdgId(pdgId), momentum(momentum), position(position), collisionType(collisionType) {}
        Nucleon(IdType id,
                int pdgId,
                const TLorentzVector &momentum,
                const TLorentzVector &position,
                unsigned short collisionType,
                const std::vector<IdType> &collidedNucleonIndices)
                : id(id),
                  pdgId(pdgId),
                  momentum(momentum),
                  position(position),
                  collisionType(collisionType),
                  collidedNucleonIndices(collidedNucleonIndices) {}
        IdType getId() const {
            return id;
        }
        void setId(IdType id) {
            Nucleon::id = id;
        }
        int getPdgId() const {
            return pdgId;
        }
        void setPdgId(int pdgId) {
            Nucleon::pdgId = pdgId;
        }
        const TLorentzVector &getMomentum() const {
            return momentum;
        }
        void setMomentum(const TLorentzVector &momentum) {
            Nucleon::momentum = momentum;
        }
        const TLorentzVector &getPosition() const {
            return position;
        }
        void setPosition(const TLorentzVector &position) {
            Nucleon::position = position;
        }
        unsigned short getCollisionType() const {
            return collisionType;
        }
        void setCollisionType(unsigned short collisionType) {
            Nucleon::collisionType = collisionType;
        }
        const std::vector<IdType> &getCollidedNucleonIndices() const {
            return collidedNucleonIndices;
        }
        void setCollidedNucleonIndices(const std::vector<IdType> &collidedNucleonIndices) {
            Nucleon::collidedNucleonIndices = collidedNucleonIndices;
        }
        void addCollidedNucleonIndex(int index) {
            collidedNucleonIndices.push_back(index);
        }

        void Clear(Option_t * = "") {
            id = 0;
            pdgId = -1;
            momentum.SetXYZT(0, 0, 0, 999);
            position.SetXYZT(0, 0, 0, 999);
            collisionType = kNoCollision;
            collidedNucleonIndices.clear();
        }

    private:
        friend class BaseConverter;
        IdType id;

        int pdgId;

        TLorentzVector momentum;
        TLorentzVector position;

        unsigned short collisionType = kNoCollision;
        std::vector<IdType> collidedNucleonIndices;
    };



    class UParticle : public TObject {

    private:
        Int_t      fIndex;        // index of this particle
        Int_t      fPdg;          // PDG code
        Int_t      fStatus;       // Status
        Int_t      fParent;       // Index of parent
        Int_t      fParentDecay;  // Parent decay index
        Int_t      fMate;         // index of last collision partner
        Int_t      fDecay;        // decay index (-1 if not decayed)
        Int_t      fChild[2];     // index of first and last child
        Double32_t fPx;           // px (GeV)
        Double32_t fPy;           // py (GeV)
        Double32_t fPz;           // pz (GeV)
        Double32_t fE;            // Energy (GeV)
        Double32_t fX;            // x (fm)
        Double32_t fY;            // y (fm)
        Double32_t fZ;            // z (fm)
        Double32_t fT;            // t (fm)
        Double32_t fWeight;       // weight

    public:
        UParticle();
        UParticle(Int_t index, Int_t pdg, Int_t status,
                  Int_t parent, Int_t parentDecay,
                  Int_t mate, Int_t decay, Int_t child[2],
                  Double_t px, Double_t py, Double_t pz, Double_t e,
                  Double_t x, Double_t y, Double_t z, Double_t t,
                  Double_t weight);
        UParticle(Int_t index, Int_t pdg, Int_t status,
                  Int_t parent, Int_t parentDecay,
                  Int_t mate, Int_t decay, Int_t child[2],
                  TLorentzVector mom, TLorentzVector pos,
                  Double_t weight);
        UParticle(const UParticle& right);
        UParticle(const TParticle& right);
        virtual ~UParticle();
        const UParticle& operator =  (const UParticle& right);
        const UParticle& operator =  (const TParticle& right);
        Bool_t     operator == (const UParticle& right) const;
        void Print(Option_t* = "") const;
        inline Int_t    GetIndex()       const {return fIndex;}
        inline Int_t    GetPdg()         const {return fPdg;}
        inline Int_t    GetStatus()      const {return fStatus;}
        inline Int_t    GetParent()      const {return fParent;}
        inline Int_t    GetParentDecay() const {return fParentDecay;}
        inline Int_t    GetMate()        const {return fMate;}
        inline Int_t    GetDecay()       const {return fDecay;}
        inline Int_t    GetFirstChild()  const {return fChild[0];}
        inline Int_t    GetLastChild()   const {return fChild[1];}
        inline Double_t Px()             const {return fPx;}
        inline Double_t Py()             const {return fPy;}
        inline Double_t Pz()             const {return fPz;}
        inline Double_t E()              const {return fE;}
        inline TLorentzVector GetMomentum() const {return TLorentzVector(fPx,fPy,fPz,fE);}
        inline void Momentum(TLorentzVector& mom) const {mom.SetPxPyPzE(fPx,fPy,fPz,fE);}
        inline Double_t X()              const {return fX;}
        inline Double_t Y()              const {return fY;}
        inline Double_t Z()              const {return fZ;}
        inline Double_t T()              const {return fT;}
        inline TLorentzVector GetPosition() const {return TLorentzVector(fX,fY,fZ,fT);}
        inline void Position(TLorentzVector& pos) const {pos.SetXYZT(fX,fY,fZ,fT);}
        inline Double_t GetWeight()      const {return fWeight;}
        inline void SetIndex      (Int_t index)       {fIndex = index;}
        inline void SetPdg        (Int_t pdg)         {fPdg = pdg;}
        inline void SetStatus     (Int_t status)      {fStatus = status;}
        inline void SetParent     (Int_t parent)      {fParent = parent;}
        inline void SetParentDecay(Int_t parentDecay) {fParentDecay = parentDecay;}
        inline void SetMate       (Int_t mate)        {fMate = mate;}
        inline void SetDecay      (Int_t decay)       {fDecay = decay;}
        inline void SetChild      (Int_t child[2])    {fChild[0]=child[0]; fChild[1]=child[1];}
        inline void SetFirstChild (Int_t child)       {fChild[0] = child;}
        inline void SetLastChild  (Int_t child)       {fChild[1] = child;}
        inline void SetPx         (Double_t px)       {fPx = px;}
        inline void SetPy         (Double_t py)       {fPy = py;}
        inline void SetPz         (Double_t pz)       {fPz = pz;}
        inline void SetE          (Double_t e)        {fE = e;}
        inline void SetMomentum(Double_t px, Double_t py, Double_t pz, Double_t e)
        {fPx = px; fPy = py; fPz = pz; fE = e;}
        inline void SetMomentum(TLorentzVector mom) {fPx=mom.Px(); fPy=mom.Py(); fPz=mom.Pz(); fE=mom.E();}
        inline void SetX          (Double_t x)        {fX = x;}
        inline void SetY          (Double_t y)        {fY = y;}
        inline void SetZ          (Double_t z)        {fZ = z;}
        inline void SetT          (Double_t t)        {fT = t;}
        inline void SetPosition(Double_t x, Double_t y, Double_t z, Double_t t)
        {fX = x; fY = y; fZ = z; fT = t;}
        inline void SetPosition(TLorentzVector pos) {fX=pos.X(); fY=pos.Y(); fZ=pos.Z(); fT=pos.T();}
        inline void SetWeight     (Double_t weight)   {fWeight = weight;}

    };


#include "TNamed.h"
#include "TString.h"


    class URun : public TNamed {

    private:
        TString    fGenerator;     // Generator description
        TString    fComment;       // Run comment
        TString    fDecayer;       // Decayer description
        Int_t      fAProj;         // Projectile mass number
        Int_t      fZProj;         // Projectile charge
        Double32_t fPProj;         // Projectile momentum per nucleon (GeV)
        Int_t      fATarg;         // Target mass number
        Int_t      fZTarg;         // Target charge
        Double32_t fPTarg;         // Target momentum per nucleon (GeV)
        Double32_t fBMin;          // Minimum impact parameter
        Double32_t fBMax;          // Maximum impact parameter
        Int_t      fBWeight;       // Impact parameter weighting
        // 0 for geometrical weights (bdb)
        // 1 for flat distribution
        Double32_t fPhiMin;        // Event plane minimum angle (rad)
        Double32_t fPhiMax;        // Event plane maximum angle (rad)
        Double32_t fSigma;         // Cross-section (mb)
        Int_t      fNEvents;       // Requested number of events

    public:
        URun();
        URun(const char* generator, const char* comment, Int_t aProj,
             Int_t zProj, Double_t pProj, Int_t aTarg, Int_t zTarg,
             Double_t pTarg, Double_t bMin, Double_t bMax, Int_t bWeight,
             Double_t phiMin, Double_t phiMax, Double_t sigma, Int_t nEvents);
        virtual ~URun();
        void Print(Option_t* = "") const;
        void GetGenerator(TString& generator) {generator = fGenerator;}
        void GetComment(TString& comment)     {comment = fComment;}
        void GetDecayer(TString& decayer)     {decayer = fDecayer;}
        inline Int_t       GetAProj()   const {return fAProj;}
        inline Int_t       GetZProj()   const {return fZProj;}
        inline Double_t    GetPProj()   const {return fPProj;}
        inline Int_t       GetATarg()   const {return fATarg;}
        inline Int_t       GetZTarg()   const {return fZTarg;}
        inline Double_t    GetPTarg()   const {return fPTarg;}
        inline Double_t    GetBMin()    const {return fBMin;}
        inline Double_t    GetBMax()    const {return fBMax;}
        inline Int_t       GetBWeight() const {return fBWeight;}
        inline Double_t    GetPhiMax()  const {return fPhiMax;}
        inline Double_t    GetPhiMin()  const {return fPhiMin;}
        inline Double_t    GetSigma()   const {return fSigma;}
        inline Int_t       GetNEvents() const {return fNEvents;}
        Double_t    GetSqrtS();
        Double_t    GetNNSqrtS();
        Double_t    GetProjectileEnergy();
        Double_t    GetTargetEnergy();
        Double_t    GetBetaCM();
        Double_t    GetGammaCM();
        inline void SetNEvents(Int_t nEvents)   {fNEvents=nEvents;}
        inline void SetPProj  (Double_t pProj)  {fPProj=pProj;}
        inline void SetPTarg  (Double_t pTarg)  {fPTarg=pTarg;}
        inline void SetDecayer(TString decayer) {fDecayer=decayer;}
    };


//--------------------------------------------------------------------
    inline UParticle::UParticle()
            : TObject(),
              fIndex(0),
              fPdg(0),
              fStatus(0),
              fParent(0),
              fParentDecay(0),
              fMate(0),
              fDecay(0),
              fChild(),
              fPx(0.),
              fPy(0.),
              fPz(0.),
              fE(0.),
              fX(0.),
              fY(0.),
              fZ(0.),
              fT(0.),
              fWeight(0.)
    {
        fChild[0] = 0;
        fChild[1] = 0;
    }
//--------------------------------------------------------------------


//--------------------------------------------------------------------
    inline UParticle::UParticle(Int_t index, Int_t pdg, Int_t status,
                         Int_t parent, Int_t parentDecay,
                         Int_t mate, Int_t decay, Int_t child[2],
                         Double_t px, Double_t py, Double_t pz, Double_t e,
                         Double_t x, Double_t y, Double_t z, Double_t t,
                         Double_t weight)
            : TObject(),
              fIndex(index),
              fPdg(pdg),
              fStatus(status),
              fParent(parent),
              fParentDecay(parentDecay),
              fMate(mate),
              fDecay(decay),
              fPx(px),
              fPy(py),
              fPz(pz),
              fE(e),
              fX(x),
              fY(y),
              fZ(z),
              fT(t),
              fWeight(weight)
    {
        fChild[0]    = child[0];
        fChild[1]    = child[1];
    }
//--------------------------------------------------------------------


//--------------------------------------------------------------------
    inline UParticle::UParticle(Int_t index, Int_t pdg, Int_t status,
                         Int_t parent, Int_t parentDecay,
                         Int_t mate, Int_t decay, Int_t child[2],
                         TLorentzVector mom, TLorentzVector pos,
                         Double_t weight)
            : TObject(),
              fIndex(index),
              fPdg(pdg),
              fStatus(status),
              fParent(parent),
              fParentDecay(parentDecay),
              fMate(mate),
              fDecay(decay),
              fPx(mom.Px()),
              fPy(mom.Py()),
              fPz(mom.Pz()),
              fE(mom.E()),
              fX(pos.X()),
              fY(pos.Y()),
              fZ(pos.Z()),
              fT(pos.T()),
              fWeight(weight)
    {
        fChild[0]    = child[0];
        fChild[1]    = child[1];
    }
//--------------------------------------------------------------------


//--------------------------------------------------------------------
    inline UParticle::UParticle(const UParticle& right)
            : TObject(right),
              fIndex(right.fIndex),
              fPdg(right.fPdg),
              fStatus(right.fStatus),
              fParent(right.fParent),
              fParentDecay(right.fParentDecay),
              fMate(right.fMate),
              fDecay(right.fDecay),
              fPx(right.fPx),
              fPy(right.fPy),
              fPz(right.fPz),
              fE(right.fE),
              fX(right.fX),
              fY(right.fY),
              fZ(right.fZ),
              fT(right.fT),
              fWeight(right.fWeight)
    {
    }
//--------------------------------------------------------------------


//--------------------------------------------------------------------
    inline UParticle::UParticle(const TParticle &right)
            : TObject(),
              fIndex(0),
              fPdg(right.GetPdgCode()),
              fStatus(right.GetStatusCode()),
              fParent(right.GetFirstMother()),
              fParentDecay(0),
              fMate(0),
              fDecay(0),
              fPx(right.Px()),
              fPy(right.Py()),
              fPz(right.Pz()),
              fE(right.Energy()),
              fX(right.Vx()),
              fY(right.Vy()),
              fZ(right.Vz()),
              fT(right.T()),
              fWeight(right.GetWeight())
    {
        fChild[0] = right.GetFirstDaughter();
        fChild[1] = right.GetLastDaughter();
    }
//--------------------------------------------------------------------


//--------------------------------------------------------------------
    inline UParticle::~UParticle()
    {
        // Destructor
    }
//--------------------------------------------------------------------


//--------------------------------------------------------------------
    inline const UParticle& UParticle::operator = (const UParticle& right)
    {
        // Assignment operator
        TObject::operator=(right);
        fIndex       = right.fIndex;
        fPdg         = right.fPdg;
        fStatus      = right.fStatus;
        fParent      = right.fParent;
        fParentDecay = right.fParentDecay;
        fMate        = right.fMate;
        fDecay       = right.fDecay;
        fChild[0]    = right.fChild[0];
        fChild[1]    = right.fChild[1];
        fPx          = right.fPx;
        fPy          = right.fPy;
        fPz          = right.fPz;
        fE           = right.fE;
        fX           = right.fX;
        fY           = right.fY;
        fZ           = right.fZ;
        fT           = right.fT;
        fWeight      = right.fWeight;
        return (*this);
    }
//--------------------------------------------------------------------


//--------------------------------------------------------------------
    inline const UParticle& UParticle::operator = (const TParticle &right)
    {
        // Assignment operator from the TParticle
        fIndex = 0;
        fPdg = right.GetPdgCode();
        fStatus = right.GetStatusCode();
        fParent = right.GetFirstMother();
        fParentDecay = 0;
        fMate = 0;
        fDecay = 0;
        fChild[0] = right.GetFirstDaughter();
        fChild[1] = right.GetLastDaughter();
        fPx = right.Px();
        fPy = right.Py();
        fPz = right.Pz();
        fE = right.Energy();
        fX = right.Vx();
        fY = right.Vy();
        fZ = right.Vz();
        fT = right.T();
        fWeight = right.GetWeight();
        return (*this);
    }
//--------------------------------------------------------------------


//--------------------------------------------------------------------
    inline Bool_t UParticle::operator == (const UParticle& right) const
    {
        // If equal operator
        return (
                fIndex       == right.fIndex &&
                fPdg         == right.fPdg &&
                fStatus      == right.fStatus &&
                fParent      == right.fParent &&
                fParentDecay == right.fParentDecay &&
                fMate        == right.fMate &&
                fDecay       == right.fDecay &&
                fChild[0]    == right.fChild[0] &&
                fChild[1]    == right.fChild[1] &&
                ((TMath::Abs((fPx-right.fPx)/fPx)<0.0001) ||
                 (TMath::Abs(fPx)<1e-16&&TMath::Abs(right.fPx)<1e-16)) &&
                ((TMath::Abs((fPy-right.fPy)/fPy)<0.0001) ||
                 (TMath::Abs(fPy)<1e-16&&TMath::Abs(right.fPy)<1e-16)) &&
                ((TMath::Abs((fPz-right.fPz)/fPz)<0.0001) ||
                 (TMath::Abs(fPz)<1e-16&&TMath::Abs(right.fPz)<1e-16)) &&
                ((TMath::Abs((fE-right.fE)/fE)<0.0001) ||
                 (TMath::Abs(fE)<1e-16&&TMath::Abs(right.fE)<1e-16)) &&
                ((TMath::Abs((fX-right.fX)/fX)<0.0001) ||
                 (TMath::Abs(fX)<1e-16&&TMath::Abs(right.fX)<1e-16)) &&
                ((TMath::Abs((fY-right.fY)/fY)<0.0001) ||
                 (TMath::Abs(fY)<1e-16&&TMath::Abs(right.fY)<1e-16)) &&
                ((TMath::Abs((fZ-right.fZ)/fZ)<0.0001) ||
                 (TMath::Abs(fZ)<1e-16&&TMath::Abs(right.fZ)<1e-16)) &&
                ((TMath::Abs((fT-right.fT)/fT)<0.0001) ||
                 (TMath::Abs(fT)<1e-16&&TMath::Abs(right.fT)<1e-16)) &&
                ((TMath::Abs((fWeight-right.fWeight)/fWeight)<0.0001) ||
                 (TMath::Abs(fWeight)<1e-16&&TMath::Abs(right.fWeight)<1e-16))
        );
    }
//--------------------------------------------------------------------


//--------------------------------------------------------------------
    inline void UParticle::Print(Option_t* /*option*/) const
    {
        // Print the data members to the standard output
        cout << "------------------------------------------------" << endl
             << "-I-                 Particle                 -I-" << endl
             << "Index                       : " << fIndex << endl
             << "PDG code                    : " << fPdg << endl
             << "Status code                 : " << fStatus << endl
             << "Parent index                : " << fParent << endl
             << "Parent decay index          : " << fParentDecay << endl
             << "Last collision partner      : " << fMate << endl
             << "Decay index                 : " << fDecay << endl
             << "First child index           : " << fChild[0] << endl
             << "Last child index            : " << fChild[1] << endl
             << "Momentum (px, py, pz) (GeV) : (" << fPx << ", " << fPy << ", " << fPz << ")" << endl
             << "Energy (GeV)                : " << fE << endl
             << "Position (x, y, z) (fm)     : (" << fX << ", " << fY << ", " << fZ << ")" << endl
             << "Creation time (fm)          : " << fT << endl
             << "Weight                      : " << fWeight << endl
             << "------------------------------------------------" << endl;
    }
//--------------------------------------------------------------------


    inline URun::URun()
            : TNamed("run","Run Header"),
              fGenerator(""),
              fComment(""),
              fDecayer(""),
              fAProj(0),
              fZProj(0),
              fPProj(0.),
              fATarg(0),
              fZTarg(0),
              fPTarg(0.),
              fBMin(0.),
              fBMax(0.),
              fBWeight(0),
              fPhiMin(0.),
              fPhiMax(0.),
              fSigma(0.),
              fNEvents(0)
    {
    }
//--------------------------------------------------------------------



//--------------------------------------------------------------------
    inline URun::URun(const char* generator, const char* comment, Int_t aProj,
               Int_t zProj, Double_t pProj, Int_t aTarg, Int_t zTarg,
               Double_t pTarg, Double_t bMin, Double_t bMax, Int_t bWeight,
               Double_t phiMin, Double_t phiMax, Double_t sigma,
               Int_t nEvents)
            : TNamed("run", "Run Header"),
              fGenerator(generator),
              fComment(comment),
              fDecayer(""),
              fAProj(aProj),
              fZProj(zProj),
              fPProj(pProj),
              fATarg(aTarg),
              fZTarg(zTarg),
              fPTarg(pTarg),
              fBMin(bMin),
              fBMax(bMax),
              fBWeight(bWeight),
              fPhiMin(phiMin),
              fPhiMax(phiMax),
              fSigma(sigma),
              fNEvents(nEvents)
    {
    }
//--------------------------------------------------------------------



//--------------------------------------------------------------------
    inline URun::~URun()
    {
        // Destructor
    }
//--------------------------------------------------------------------



//--------------------------------------------------------------------
    inline void URun::Print(Option_t* /*option*/) const
    {
        // Print all data members to the standard output
        cout << "--------------------------------------------------" << endl
             << "-I-                 Run Header                 -I-" << endl
             << "Generator                     : " << fGenerator << endl
             << "Comment                       : " << fComment << endl
             << "Decayer                       : " << fDecayer << endl
             << "Projectile mass               : " << fAProj << endl
             << "Projectile charge             : " << fZProj << endl
             << "Projectile momentum (AGeV/c)  : " << fPProj << endl
             << "Target mass                   : " << fATarg << endl
             << "Target charge                 : " << fZTarg << endl
             << "Target momentum (AGeV/c)      : " << fPTarg << endl
             << "Minimal impact parameter (fm) : " << fBMin << endl
             << "Maximal impact parameter (fm) : " << fBMax << endl
             << "Impact parameter weightning   : " << fBWeight << endl
             << "Minimal azimuthal angle (rad) : " << fPhiMin << endl
             << "Maximal azimuthal angle (rad) : " << fPhiMax << endl
             << "Cross-section (mb)            : " << fSigma << endl
             << "Requested number of events    : " << fNEvents << endl
             << "--------------------------------------------------" << endl;
    }
//--------------------------------------------------------------------



//--------------------------------------------------------------------
    inline Double_t URun::GetProjectileEnergy()
    {
        // Get the projectile energy
        Double_t mProt = 0.938272029;
        Double_t mNeut = 0.939565360;
        Double_t mPion = 0.13957018;
        Double_t eProj = 0.;
        if ( fAProj > 0 )          // nucleus
            eProj = fZProj  * TMath::Sqrt( fPProj*fPProj + mProt*mProt )
                    + (fAProj - fZProj) * TMath::Sqrt( fPProj*fPProj + mNeut*mNeut );
        else if ( fAProj ==  0 )   // photon
            eProj = fPProj;
        else if ( fAProj == -1 )   // pion
            eProj = TMath::Sqrt( fPProj*fPProj + mPion*mPion );
        else cout << "Warning:: URun: Projectile mass " << fAProj
                  << " not valid! " << endl;
        return eProj;
    }
//--------------------------------------------------------------------



//--------------------------------------------------------------------
    inline Double_t URun::GetTargetEnergy()
    {
        // Get the target energy
        Double_t mProt = 0.938272029;
        Double_t mNeut = 0.939565360;
        Double_t mPion = 0.13957018;
        Double_t eTarg = 0.;
        if ( fATarg > 0 )            // nucleus
            eTarg = fZTarg  * TMath::Sqrt( fPTarg*fPTarg + mProt*mProt )
                    + (fATarg - fZTarg) * TMath::Sqrt( fPTarg*fPTarg + mNeut*mNeut );
        else if ( fAProj ==  0 )     // photon
            eTarg = fPTarg;
        else if ( fAProj == -1 )     // pion
            eTarg = TMath::Sqrt( fPTarg*fPTarg + mPion*mPion );
        else cout << "Warning:: URun: Target mass " << fATarg
                  << " not valid! " << endl;
        return eTarg;
    }
//--------------------------------------------------------------------



//--------------------------------------------------------------------
    inline Double_t URun::GetNNSqrtS()
    {
        // Get the cm energy
        Double_t eSum = TMath::Sqrt( fPTarg*fPTarg + 0.938272029*0.938272029 ) + TMath::Sqrt( fPProj*fPProj + 0.938272029*0.938272029 );
        Double_t pSum = Double_t(fPProj + fPTarg);
        Double_t ecm = TMath::Sqrt( eSum*eSum - pSum*pSum );
        return ecm;
    }
//--------------------------------------------------------------------

    inline Double_t URun::GetSqrtS()
    {
        // Get the cm energy
        Double_t eSum = GetProjectileEnergy() + GetTargetEnergy();
        Double_t pSum = Double_t(fAProj) * fPProj + Double_t(fATarg) * fPTarg;
        Double_t ecm = TMath::Sqrt( eSum*eSum - pSum*pSum );
        return ecm;
    }
//--------------------------------------------------------------------



//--------------------------------------------------------------------
    inline Double_t URun::GetBetaCM()
    {
        // Get cm velocity
        Double_t eSum = GetProjectileEnergy() + GetTargetEnergy();
        Double_t pSum = Double_t(fAProj) * fPProj + Double_t(fATarg) * fPTarg;
        return pSum / eSum;
    }
//--------------------------------------------------------------------



//--------------------------------------------------------------------
    inline Double_t URun::GetGammaCM()
    {
        // Get cm gamma factor
        Double_t betaCM = GetBetaCM();
        return 1. / TMath::Sqrt( 1. - betaCM*betaCM );
    }
}

#endif //GRATE_MCINI_HH
