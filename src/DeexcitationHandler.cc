#include "DeexcitationHandler.hh"


DeexcitationHandler::DeexcitationHandler(){
    aE = 1/(2.*(upBoundTransitionForMF - lowBoundTransitionForMF));
    E0 = (upBoundTransitionForMF + lowBoundTransitionForMF)/2.;
    nist = G4NistManager::Instance();
}

DeexcitationHandler::~DeexcitationHandler() = default;

G4ReactionProductVector *DeexcitationHandler::BreakUp(const G4Fragment &theInitialFragment) {
    //std::cout<<"isMF = "<<isMultifragmentation(theInitialFragment.GetA(), theInitialFragment.GetZ(), theInitialFragment.GetExcitationEnergy())<<" A = "<<theInitialFragment.GetA()<<"\n";
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
    G4FragmentVector* tempResult = new G4FragmentVector();
    G4ReactionProductVector* theResult = new G4ReactionProductVector();
    G4ReactionProductVector* ablaTempResult = new G4ReactionProductVector();
    G4FragmentVector* toDecayVector = new G4FragmentVector();

    G4double Ain = theInitialFragment.GetA();
    G4double Zin = theInitialFragment.GetZ();
    G4double exEn = theInitialFragment.GetExcitationEnergy();

    if(isDecayOfPureNeutrons(Ain, Zin)){return BreakUpPureNeutrons(theInitialFragment);}

    if(isDecay(Ain, Zin, exEn)) {
        if (isMultifragmentation(Ain,Zin,exEn)) {
            tempResult = theMultifragmentation.BreakItUp(theInitialFragment);
        } else if (isFermiBreakUp(Ain,Zin,exEn)) {
            tempResult = FermiBreakUp.BreakItUp(theInitialFragment);
        } else {
            theResult = ablaEvaporation.DeExcite(theInitialFragment);
            return theResult;
        }
    } else{theResult->push_back(toReactionProduct(const_cast<G4Fragment *>(&theInitialFragment))); return theResult;}

    // TODO find a method to reduce code duplication
    // - probably some other conditions for cycle finish will help

    for (G4FragmentVector::iterator j = tempResult->begin(); j != tempResult->end(); ++j) {
        if (!isDecay((*j)->GetA(), (*j)->GetZ(), (*j)->GetExcitationEnergy())) {
             theResult->push_back(toReactionProduct((*j)));
        }
        else { toDecayVector->push_back((*j)); }
    }
    tempResult->clear();

    bool flagFBU = true;

    while (toDecayVector->size() != 0) {


        G4FragmentVector::iterator j = toDecayVector->begin();
        while (j != toDecayVector->end()) {
            if (isFermiBreakUp((*j)->GetA(), (*j)->GetZ(), (*j)->GetExcitationEnergy()) && flagFBU) {
                tempResult = FermiBreakUp.BreakItUp(*(*j));
                if(tempResult->size() == 1) flagFBU = false;
            } else {
                ablaTempResult = ablaEvaporation.DeExcite(*(*j));
                theResult->insert(theResult->end(), ablaTempResult->begin(), ablaTempResult->end());
            }
            j = toDecayVector->erase(j);
        }
        //std::cout<<"1) toDecayVector->size() "<<toDecayVector->size()<<std::endl;
        for (G4FragmentVector::iterator j = tempResult->begin(); j != tempResult->end(); ++j) {
            if (!isDecay((*j)->GetA(), (*j)->GetZ(), (*j)->GetExcitationEnergy())) {
                theResult->push_back(toReactionProduct((*j)));
            }
            else { toDecayVector->push_back((*j)); }
        }

        //if(toDecayVector->size() != 0) std::cout<<"2) toDecayVector->at(0) "<<toDecayVector->at(0)->GetA()<<std::endl;
        tempResult->clear();

    }

   //std::cout<<"theResult->size() "<<theResult->size()<<std::endl;
    return theResult;
}

G4Fragment *DeexcitationHandler::toFragment(G4ReactionProduct *product) {
    G4LorentzVector fragLorentzVector(product->GetTotalEnergy(), product->GetMomentum());
    G4Fragment* newFrag = new G4Fragment(fragLorentzVector, product->GetDefinition());
    return newFrag;
}

G4ReactionProduct *DeexcitationHandler::toReactionProduct(G4Fragment* fragment) {
    const G4ParticleDefinition *def = toParticleDefinition(fragment->GetA(), fragment->GetZ());
    if(def == 0) {return 0;}
    G4ReactionProduct* newProduct = new G4ReactionProduct(def);
    const G4LorentzVector LV = fragment->GetMomentum();
    newProduct->SetMomentum(LV.px(), LV.py(), LV.pz());
    newProduct->SetTotalEnergy(LV.e());
    newProduct->SetFormationTime(fragment->GetCreationTime());
    return newProduct;
}

G4ParticleDefinition *DeexcitationHandler::toParticleDefinition(G4int A, G4int Z) const{
    if     (A == 1 && Z == 1)  return G4Proton::Proton();
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
}

bool DeexcitationHandler::isMultifragmentation(G4double A, G4double Z, G4double Ex) {
    if(isFermiBreakUp(A,Z,Ex)){return false;}
    G4double w = G4RandFlat::shoot();
    G4double transF = 0.5*tanh((Ex/A - E0)/aE) + 0.5;
    if (Ex < lowBoundTransitionForMF*A) {return false;}
    else if (w < transF && Ex < upBoundTransitionForMF*A){return true;}
    else if (w > transF && Ex < upBoundTransitionForMF*A){return false;}
    else if (Ex > upBoundTransitionForMF*A) {return true;}
    else {return false;}
}

bool DeexcitationHandler::isFermiBreakUp(G4double A, G4double Z, G4double Ex) {
    if(A < MaxAforFermiBreakUp && Z < MaxZforFermiBreakUp && Ex > minExForFBU*A){return true;}
    else{return false;}
}

bool DeexcitationHandler::isDecay(G4double A, G4double Z, G4double Ex) {
    if(A > 1 && Ex > minEx){return true;}
    else if (nist->GetIsotopeAbundance(A,Z) > 0 ){return true;}
    else{return false;}
}

bool DeexcitationHandler::isDecayOfPureNeutrons(G4double A, G4double Z) {
    if(A<MaxAforFermiBreakUpForPureNeutronFragments && Z==0) return true;
    else return false;
}

