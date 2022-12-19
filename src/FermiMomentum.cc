#include <utility>

#include "FermiMomentum.hh"


FermiMomentum::FermiMomentum(std::shared_ptr<AAMCCinput> input_in, const std::string& model_in) : engine(nullptr),
    randGauss(nullptr), randFlat(nullptr), randFlatphi(nullptr), input(std::move(input_in)) {
    if(model_in == "M") modelInt = 0;
    if(model_in == "G") modelInt = 1;
    if(model_in == "V") modelInt = 2;
    if(model_in == "VAL") modelInt = -1;

    if (modelInt >= 0) {
        engine = std::unique_ptr<CLHEP::RanecuEngine>(new CLHEP::RanecuEngine());
        randGauss = std::unique_ptr<CLHEP::RandGauss>(new CLHEP::RandGauss(*engine, 0, 1));
        randFlatphi = std::unique_ptr<CLHEP::RandFlat>(new CLHEP::RandFlat(*engine, 6.28318530718));
        randFlat = std::unique_ptr<CLHEP::RandFlat>(new CLHEP::RandFlat(*engine, 1.));
    }
}

vect3 FermiMomentum::GetMomentum(std::string side) {
    //std::cout<<"Called for side "<<side<<std::endl;
    vect3 out = {0, 0, 0};
        A = input->nucleons.GetTotA(side);
        Spec = input->nucleons.GetA(side);
        Npart = A - Spec;

        if(side == "A"){
            if ((pF.px == 0 && pF.py == 0 && pF.pz == 0) || modelInt == -1){
                switch (modelInt) {
                    case 0: out = GetMorrisey(); break;
                    case 1: out = GetGoldhaber(); break;
                    case 2: out = GetVanBiber(); break;
                    case -1: out = {input->FermiMomA_x, input->FermiMomA_y, input->FermiMomA_z}; break;
                    default: out = GetMorrisey(); break;
                }
                pF = out;
                if(modelInt > 0) out.pz = gammafacA * (out.pz + betafacA * pow(out.mag2() + pow (nucleonAverMass * (input->nucleons.GetA(side)), 2), 0.5));
                return out;
            }
            else {out = pF; pF = {0, 0, 0};
                if(modelInt > 0) out.pz = gammafacA * (out.pz + betafacA * pow(out.mag2() + pow (nucleonAverMass * (input->nucleons.GetA(side)), 2), 0.5));
                return out;}
        }
        else if(side == "B") {
            if ((pF.px == 0 && pF.py == 0 && pF.pz == 0) || modelInt == -1){
                switch (modelInt) {
                    case 0: out = GetMorrisey(); break;
                    case 1: out = GetGoldhaber(); break;
                    case 2: out = GetVanBiber(); break;
                    case -1: out = {input->FermiMomB_x, input->FermiMomB_y, input->FermiMomB_z}; break;
                    default: out = GetMorrisey(); break;
                }
                pF = {-out.px, -out.py, -out.pz};
                if(modelInt > 0)  out.pz = gammafacB * (out.pz - betafacB * pow(out.mag2() + pow (nucleonAverMass * (input->nucleons.GetA(side)), 2), 0.5));
                return out;
            }
            else {out = {-pF.px, -pF.py, -pF.pz}; pF = {0, 0, 0};
                if(modelInt > 0)   out.pz = gammafacB * (out.pz - betafacB * pow(out.mag2() + pow (nucleonAverMass * (input->nucleons.GetA(side)), 2), 0.5));
                return out;}


        }
        else {std::cout<< "Error, side argument is invalid. Please, use A or B as an argument"; return out;}
}

vect3 FermiMomentum::GetGoldhaber() {
    double GoldhaberDev = GoldhaberDev0 * pow(double(Spec * Npart) / double(A - 1), 0.5);
    double p = GoldhaberDev*randGauss->fire();
    double phi = randFlatphi->fire();
	double ctheta = 2 * randFlat->fire() - 1;
	double stheta = std::sqrt(1 - ctheta * ctheta);
	return { p * std::cos(phi) * stheta, p * std::sin(phi) * stheta, p * ctheta };
}

vect3 FermiMomentum::GetMorrisey() {
    double MorisseyDev = MorisseyDev0 * pow(Npart,0.5);
    double p = randGauss->fire()*MorisseyDev; //*Spec
    double phi = randFlatphi->fire();
	double ctheta = 2 * randFlat->fire() - 1;
	double stheta = std::sqrt(1 - ctheta * ctheta);
	return { p * std::cos(phi) * stheta, p * std::sin(phi) * stheta, p * ctheta };
}

vect3 FermiMomentum::GetVanBiber() {
    double MorisseyDev = MorisseyDev0 * pow(Npart,0.5);
    double VBDev = VBDev1*pow(Spec*(Spec-1)/(A*(A-1)),0.5);
    double p = randGauss->fire()*std::sqrt(MorisseyDev*MorisseyDev + VBDev*VBDev);
    double phi = randFlatphi->fire();
	double ctheta = 2 * randFlat->fire() - 1;
	double stheta = std::sqrt(1 - ctheta * ctheta);
	return { p * std::cos(phi) * stheta, p * std::sin(phi) * stheta, p * ctheta };
}

void FermiMomentum::SetPzPerNucleon(double pZA, double pZB) {
    pZperNucleonA = pZA; pZperNucleonB = pZB;
    gammafacA = pow(1 + pow((pZperNucleonA) / (nucleonAverMass), 2), 0.5);
    gammafacB = pow(1 + pow((pZperNucleonB) / (nucleonAverMass), 2), 0.5);
    betafacA = pow(1 - 1 / (gammafacA * gammafacA), 0.5); betafacB = pow(1 - 1 / (gammafacB * gammafacB), 0.5);
}

CLHEP::Hep3Vector FermiMomentum::GetBoost(std::string side) {
    if(input->nucleons.GetA(side) != 0){
        vect3 p = GetMomentum(side);
    double  E =std::sqrt( p.mag2()/(input->nucleons.GetA(side)*input->nucleons.GetA(side)) + nucleonAverMass*nucleonAverMass);
    CLHEP::HepLorentzVector  futureBoost(p.px/input->nucleons.GetA(side), p.py/input->nucleons.GetA(side), p.pz/input->nucleons.GetA(side), E);
    if(!futureBoost.isTimelike()) std::cout<<side<<" "<<E<<" nucl = "<<input->nucleons.GetA(side)<<std::endl;
    return futureBoost.boostVector();}
    else { return CLHEP::Hep3Vector(0,0,0);}
}

CLHEP::Hep3Vector FermiMomentum::toHep3Vector(vect3 vectorIn) {
    return {vectorIn.px, vectorIn.py, vectorIn.pz};
}

CLHEP::HepLorentzVector FermiMomentum::toHepLorentzVector(vect3 vectorIn , std::string side) {
    if(input->nucleons.GetA(side) != 0){
        vect3 p = vectorIn;
        double  E =std::sqrt( p.mag2()/(input->nucleons.GetA(side)*input->nucleons.GetA(side)) + nucleonAverMass*nucleonAverMass);
        return CLHEP::HepLorentzVector(p.px, p.py, p.pz, E*input->nucleons.GetA(side));
    } else{return {0,0,0, 0};}
}

CLHEP::Hep3Vector FermiMomentum::GetMomentumHep3(std::string side) {
    return toHep3Vector(GetMomentum(side));
}

CLHEP::HepLorentzVector FermiMomentum::GetLorentzVector(std::string side) {
    return toHepLorentzVector(GetMomentum(side),side);
}

