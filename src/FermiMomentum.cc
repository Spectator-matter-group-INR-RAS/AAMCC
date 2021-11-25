#include "FermiMomentum.hh"

FermiMomentum::FermiMomentum( NucleonVector* nucleons_in,std::string model_in) {
    if(model_in == "M") modelInt = 0;
    if(model_in == "G") modelInt = 1;
    if(model_in == "V") modelInt = 2;

    engine = new CLHEP::RanecuEngine();
    randGauss = new CLHEP::RandGauss(engine,0,1);
    randFlatphi = new CLHEP::RandFlat(engine,6.28318530718);
	randFlat = new CLHEP::RandFlat(engine, 1.);
    nucleons = nucleons_in;
}

vect3 FermiMomentum::GetMomentum(std::string side) {
    //std::cout<<"Called for side "<<side<<std::endl;
    vect3 out = {0, 0, 0};
        A = nucleons->GetTotA(side);
        Spec = nucleons->GetA(side);
        Npart = A - Spec;

        if(side == "A"){
            if (pF.px == 0 && pF.py == 0 && pF.pz == 0){
                switch (modelInt) {
                    case 0: out = GetMorrisey(); break;
                    case 1: out = GetGoldhaber(); break;
                    case 2: out = GetVanBiber(); break;
                    default: out = GetMorrisey(); break;
                }
                pF = out;
                out.pz = gammafacA * (out.pz + betafacA * pow(out.mag2() + pow (nucleonAverMass * (nucleons->GetA(side)), 2), 0.5));
                return out;
            }
            else {out = pF; pF = {0, 0, 0};
                out.pz = gammafacA * (out.pz + betafacA * pow(out.mag2() + pow (nucleonAverMass * (nucleons->GetA(side)), 2), 0.5));
                return out;}
        }
        else if(side == "B") {
                if (pF.px == 0 && pF.py == 0 && pF.pz == 0){
                    switch (modelInt) {
                        case 0: out = GetMorrisey(); break;
                        case 1: out = GetGoldhaber(); break;
                        case 2: out = GetVanBiber(); break;
                        default: out = GetMorrisey(); break;
                    }
                    pF = {-out.px, -out.py, -out.pz};
                    out.pz = gammafacB * (out.pz - betafacB * pow(out.mag2() + pow (nucleonAverMass * (nucleons->GetA(side)), 2), 0.5));
                    return out;
                }
                else {out = {-pF.px, -pF.py, -pF.pz}; pF = {0, 0, 0};
                    out.pz = gammafacB * (out.pz - betafacB * pow(out.mag2() + pow (nucleonAverMass * (nucleons->GetA(side)), 2), 0.5));
                    return out;}
        }
        else {std::cout<< "Error, side argument is invalid. Please, use A or B as an argument"; return out;}

    return out;
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
    if(nucleons->GetA(side) != 0){
        vect3 p = GetMomentum(side);
    double  E =std::sqrt( p.mag2()/(nucleons->GetA(side)*nucleons->GetA(side)) + nucleonAverMass*nucleonAverMass);
    CLHEP::HepLorentzVector  futureBoost(p.px/nucleons->GetA(side), p.py/nucleons->GetA(side), p.pz/nucleons->GetA(side), E);
    if(!futureBoost.isTimelike()) std::cout<<side<<" "<<E<<" nucl = "<<nucleons->GetA(side)<<std::endl;
    return futureBoost.boostVector();}
    else { return CLHEP::Hep3Vector(0,0,0);}
}
