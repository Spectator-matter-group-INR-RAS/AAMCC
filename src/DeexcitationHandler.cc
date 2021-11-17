#include "DeexcitationHandler.hh"


DeexcitationHandler::DeexcitationHandler(){
}

DeexcitationHandler::~DeexcitationHandler() = default;

G4ReactionProductVector *DeexcitationHandler::BreakUp(const G4Fragment &theInitialFragment) {
    if(theInitialFragment.GetA()>MaxAforFermiBreakUpForPureNeutronFragments && theInitialFragment.GetZ()==0){return BreakUpPureNeutrons(theInitialFragment);}
   else { return BreakItUp(theInitialFragment);} //parent class method
} 

G4ReactionProductVector *DeexcitationHandler::BreakUpPureNeutrons(const G4Fragment &theInitialFragment) {
    G4ParticleDefinition *neutron = G4Neutron::NeutronDefinition();
    G4ReactionProductVector *outVec = new G4ReactionProductVector();
    std::vector<G4double> mr;
    mr.reserve(theInitialFragment.GetA_asInt());
    G4double M = theInitialFragment.GetMomentum().m();
    Hep3Vector momInitPerNucleon = theInitialFragment.GetMomentum().vect() / theInitialFragment.GetA();

        for (G4int k = 0; k < theInitialFragment.GetA_asInt(); k++) {
            mr.push_back(mn);
        }
        std::vector<G4LorentzVector *> *mom = PhaseSpaceDecay.Decay(M, mr);
        for (G4int k = 0; k < theInitialFragment.GetA_asInt(); k++) {
            G4ReactionProduct *New = new G4ReactionProduct(neutron);
            New->SetMomentum(mom->at(k)->vect() + momInitPerNucleon);
            //Needs an additional boost since now;
            double etot = mn + New->GetMomentum().mag() * c_light;
            New->SetTotalEnergy(etot);
            outVec->push_back(New);
        }
        return outVec;
    }

G4ReactionProductVector *DeexcitationHandler::futureBreakItUp(const G4Fragment &theInitialFragment) {
    return nullptr;
}

G4Fragment *DeexcitationHandler::toFragment(G4ReactionProduct *product) {
    G4LorentzVector fragLorentzVector(product->GetTotalEnergy(), product->GetMomentum());
    G4Fragment* newFrag = new G4Fragment(fragLorentzVector, product->GetDefinition());
    return newFrag;
}

G4ReactionProduct *DeexcitationHandler::toReactionProduct(G4Fragment* fragment) {
    const G4ParticleDefinition *def = fragment->GetParticleDefinition();
    if(def == 0) {return 0;}
    G4ReactionProduct* newProduct = new G4ReactionProduct(def);
    const G4LorentzVector LV = fragment->GetMomentum();
    newProduct->SetMomentum(LV.px(), LV.py(), LV.pz());
    newProduct->SetTotalEnergy(LV.e());
    newProduct->SetFormationTime(fragment->GetCreationTime());
    return newProduct;
}

G4ParticleDefinition *DeexcitationHandler::toParticleDefinition(G4int A, G4int Z) const{ if     (A == 1 && Z == 1)  return G4Proton::Proton();
    else if(A == 1 && Z == 0)  return G4Neutron::Neutron();
    else if(A == -1 && Z == 1)  return G4PionPlus::PionPlus();
    else if(A == -1 && Z == -1) return G4PionMinus::PionMinus();
    else if(A == -1 && Z == 0)  return G4PionZero::PionZero();
    else if(A == 0 && Z == 0)  return G4Gamma::Gamma();
    else if(A == 2 && Z == 1)  return G4Deuteron::Deuteron();
    else if(A == 3 && Z == 1)  return G4Triton::Triton();
    else if(A == 3 && Z == 2)  return G4He3::He3();
    else if(A == 4 && Z == 2)  return G4Alpha::Alpha();
    else if(A > 0 && Z > 0 && A >= Z) { // Returns ground state ion definition
        return G4IonTable::GetIonTable()->GetIon(Z, A, 0);
    } else { // Error, unrecognized particle
        G4cout << "Can't convert particle with A=" << A << ", Z=" << Z << " to G4ParticleDefinition, trouble ahead" << G4endl;
        return 0;
    }
};

