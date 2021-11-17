
#include <list>

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
 private:
    G4int MaxAforFermiBreakUpForPureNeutronFragments = 0;
    G4int MaxAforFermiBreakUp = 19;
    G4int MaxZforFermiBreakUp =  9;
    G4double mn = 939.5731*MeV;
    G4Fragment* toFragment(G4ReactionProduct* product);
    G4ReactionProduct* toReactionProduct(G4Fragment* fragment);
    G4ParticleDefinition* toParticleDefinition(G4int A, G4int Z) const;

    G4FermiPhaseSpaceDecay PhaseSpaceDecay;
    G4FermiBreakUp* FermiBreakUp;
};
