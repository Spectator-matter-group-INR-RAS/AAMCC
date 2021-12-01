
#include "globals.hh"

#include "G4VPreCompoundModel.hh"
#include "G4ReactionProduct.hh"
#include "G4Fragment.hh"
#include "G4HadFinalState.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4Abla.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"

class AblaEvaporation {
    public:
     AblaEvaporation();
     ~AblaEvaporation();

     G4ReactionProductVector *DeExcite(const G4Fragment &aFragment);

     inline void SetFreezeOutT(G4double T){T_freeze_out = T; if(theABLAModel) theABLAModel->SetFreezeOutT(T);};
     inline G4double GetFreezeOutT() {return T_freeze_out;}

    private:
        G4VarNtp *ablaResult;
        G4Volant *volant;
        G4Abla *theABLAModel;
        G4long eventNumber;

        G4double T_freeze_out = 1e100;

        /// \brief Convert an Abla particle to a G4ReactionProduct
        G4ReactionProduct *toG4Particle(G4int A, G4int Z , G4double kinE, G4double px, G4double py, G4double pz) const;

        /// \brief Convert A and Z to a G4ParticleDefinition
        G4ParticleDefinition *toG4ParticleDefinition (G4int A, G4int Z) const;

    };
