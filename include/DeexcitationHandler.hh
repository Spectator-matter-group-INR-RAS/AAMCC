
#include <list>
#include <math.h>

 #include "G4ExcitationHandler.hh"
 #include "G4SystemOfUnits.hh"
 #include "G4LorentzVector.hh"
 #include "G4NistManager.hh"
 #include "G4ParticleTable.hh"
 #include "G4ParticleTypes.hh"
 #include "G4Ions.hh"
 #include "G4Evaporation.hh"
 #include "G4StatMF.hh"
 #include "G4PhotonEvaporation.hh"
 #include "G4Pow.hh"
 #include "G4FermiPhaseSpaceDecay.hh"
 #include "../FermiBreakUp/G4FermiBreakUp.hh"
 #include "../Abla/include/AblaEvaporation.hh"
 #include "../Multifragmentation/include/G4StatMF.hh"

 class DeexcitationHandler: public G4ExcitationHandler {
 public:
     DeexcitationHandler();
    ~DeexcitationHandler();

    G4ReactionProductVector* BreakUp(const G4Fragment &theInitialFragment);
    G4ReactionProductVector* BreakUpPureNeutrons(const G4Fragment &theInitialFragment);
    G4ReactionProductVector* futureBreakItUp(const G4Fragment &theInitialFragment);
    inline void SetMaxAforPureNeutronFragments(G4int in_A) {MaxAforFermiBreakUpForPureNeutronFragments = in_A;};
    inline void SetMaxAforFermiBreakUp(G4int in_A) {MaxAforFermiBreakUp = in_A;};
    inline void SetMaxZforFermiBreakUp(G4int in_Z) {MaxZforFermiBreakUp = in_Z;};
    inline void SetMinExForFermiBreakUp(G4double in_Ex) {minExForFBU = in_Ex;};
    inline void SetExForMF(G4double in_lowEx, G4double in_upEx) {
        lowBoundTransitionForMF = in_lowEx; upBoundTransitionForMF = in_upEx;
        aE = 1/(2.*(upBoundTransitionForMF - lowBoundTransitionForMF));
        E0 = (upBoundTransitionForMF + lowBoundTransitionForMF)/2.;
    };
    inline void SetMinEx(G4double in_minEx) {minEx = in_minEx;};
 private:
    G4int MaxAforFermiBreakUpForPureNeutronFragments = 200;
    G4int MaxAforFermiBreakUp = 19;
    G4int MaxZforFermiBreakUp =  9;
    G4double minExForFBU = 0.1*MeV;
    G4double mn = 939.5731*MeV;//TODO switch nucelon mass to Geant4 nucleon mass. ALSO FOR FERMI MOMENTUM.
    G4double lowBoundTransitionForMF = 3*MeV;
    G4double upBoundTransitionForMF = 5*MeV;
    G4double minEx = 0.001*MeV;

    G4double aE = 1;
    G4double E0 = 0;

    G4Fragment* toFragment(G4ReactionProduct* product);
    G4ReactionProduct* toReactionProduct(G4Fragment* fragment);
    G4ParticleDefinition* toParticleDefinition(G4int A, G4int Z) const;
    bool isMultifragmentation(G4double A, G4double Z, G4double Ex);
    bool isFermiBreakUp(G4double A, G4double Z, G4double Ex);
    bool isDecay(G4double A, G4double Z, G4double Ex);
    bool isDecayOfPureNeutrons(G4double A, G4double Z);

    G4FermiPhaseSpaceDecay PhaseSpaceDecay;
    G4FermiBreakUp FermiBreakUp;
    AblaEvaporation ablaEvaporation;
    G4StatMF theMultifragmentation;

    G4NistManager* nist;
};
