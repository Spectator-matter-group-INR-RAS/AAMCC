#include "../include/AblaEvaporation.hh"

AblaEvaporation::AblaEvaporation() :
        //G4VPreCompoundModel(NULL, "ABLA"),
        ablaResult(new G4VarNtp),
        volant(new G4Volant),
        theABLAModel(new G4Abla(volant, ablaResult)),
        eventNumber(0)
{
    theABLAModel->initEvapora();
    theABLAModel->SetParameters();
    theABLAModel->SetFreezeOutT(T_freeze_out);
}

AblaEvaporation::~AblaEvaporation() {
    delete volant;
    delete ablaResult;
    delete theABLAModel;
}

G4ReactionProductVector *AblaEvaporation::DeExcite(const G4Fragment &aFragment) {
    volant->clear();
    ablaResult->clear();

    const G4int ARem = aFragment.GetA_asInt();
    const G4int ZRem = aFragment.GetZ_asInt();
    const G4double eStarRem = aFragment.GetExcitationEnergy() / MeV;
    const G4double jRem = aFragment.GetAngularMomentum().mag() / CLHEP::hbar_Planck;
    const G4LorentzVector &pRem = aFragment.GetMomentum();
    const G4double pxRem = pRem.x() / MeV;
    const G4double pyRem = pRem.y() / MeV;
    const G4double pzRem = pRem.z() / MeV;

    eventNumber++;

    theABLAModel->DeexcitationAblaxx(ARem, ZRem,
                                     eStarRem,
                                     jRem,
                                     pxRem,
                                     pyRem,
                                     pzRem,
                                     eventNumber);

    G4ReactionProductVector *result = new G4ReactionProductVector;

    for(int j = 0; j < ablaResult->ntrack; ++j) {
        //std::cout<<" volant->pcv["<<j<<"] = "<<volant->pcv[j]<<"\n";
        G4ReactionProduct *product = toG4Particle(ablaResult->avv[j],
                                                  ablaResult->zvv[j],
                                                  ablaResult->enerj[j],
                                                  ablaResult->pxlab[j],
                                                  ablaResult->pylab[j],
                                                  ablaResult->pzlab[j]);
                                                  /*ablaResult->plab[j]*std::sin(ablaResult->tetlab[j]*CLHEP::pi/180.0)*std::cos(ablaResult->philab[j]*CLHEP::pi/180.0),
                                                  ablaResult->plab[j]*std::sin(ablaResult->tetlab[j]*CLHEP::pi/180.0)*std::sin(ablaResult->philab[j]*CLHEP::pi/180.0),
                                                  ablaResult->plab[j]*std::cos(ablaResult->tetlab[j]*CLHEP::pi/180.0));*/
        if(product)
            result->push_back(product);
    }
    return result;
}

G4ParticleDefinition *AblaEvaporation::toG4ParticleDefinition(G4int A, G4int Z) const {
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

G4ReactionProduct *AblaEvaporation::toG4Particle(G4int A, G4int Z,
                                                 G4double kinE,
                                                 G4double px,
                                                 G4double py, G4double pz) const {
    G4ParticleDefinition *def = toG4ParticleDefinition(A, Z);
    if(def == 0) { // Check if we have a valid particle definition
        return 0;
    }
    G4ReactionProduct* newProduct = new G4ReactionProduct(def);
    newProduct->SetMomentum(px,py,pz);
    if(def->GetParticleName() != "gamma") newProduct->SetTotalEnergy(G4NucleiProperties::GetNuclearMass(A,Z) + kinE);
    else{newProduct->SetTotalEnergy(kinE);}
    newProduct->SetFormationTime(0);
    return newProduct;
}

