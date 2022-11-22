#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "Nucleon.hh"
#include "G4LorentzVector.hh"
#include "AAMCConstants.hh"

struct vect3{
    double px;
    double py;
	double pz;
    double  mag(){
        return std::sqrt(px*px + py*py + pz*pz);
    };
    double mag2(){
        return px*px + py*py + pz*pz;
    }
};

class FermiMomentum{
public:
     FermiMomentum(std::shared_ptr<AAMCCinput> input, const std::string& model_in = "M");
    ~FermiMomentum() = default;
    vect3 GetMomentum(std::string side);
    CLHEP::Hep3Vector GetBoost(std::string side);
    CLHEP::Hep3Vector GetMomentumHep3(std::string side);
    CLHEP::HepLorentzVector GetLorentzVector(std::string side);

    void SetPzPerNucleon(double  pZA, double pZB);

private:
    //Goldhaber model parameter
    double GoldhaberDev0 = 193*MeV;
    //Morissey model parameter
    double MorisseyDev0 = 150*MeV;
    //VanBiber model parameter
    double VBDev1 = 146*MeV;

    double Npart;
    double Spec;
    double A;

    double pZperNucleonA;
    double pZperNucleonB;
    double gammafacA;
    double gammafacB;
    double betafacA;
    double betafacB;

    std::unique_ptr<CLHEP::RanecuEngine> engine;
    std::unique_ptr<CLHEP::RandGauss> randGauss;
    std::unique_ptr<CLHEP::RandFlat> randFlat;
    std::unique_ptr<CLHEP::RandFlat> randFlatphi;

    int modelInt = 0;
    std::shared_ptr<AAMCCinput> input;

    vect3 pF = {0, 0, 0};

    vect3 GetGoldhaber();
    vect3 GetMorrisey();
    vect3 GetVanBiber();

    CLHEP::Hep3Vector toHep3Vector(vect3 vectorIn);
    CLHEP::HepLorentzVector toHepLorentzVector(vect3 vectorIn, std::string side);

};
