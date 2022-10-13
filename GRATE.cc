// Igor Pshenichnov 06.08.2008 + Alexandr Svetlichnyi at present
// ROOT 6 or higher and Geant4installation is reqired
// Geant4 v10.4 -- v9.6 is supported.

#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4SystemOfUnits.hh"
  
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4ExcitationHandler.hh"
#include "G4NucleiProperties.hh"
#include "G4Evaporation.hh"

#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4Threading.hh"

#include "G4IonTable.hh"
#include "G4GenericIon.hh"
#include "G4Ions.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclearLevelData.hh"

#include "TVector3.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

#include "GRATEmanager.hh"
/*#include "InitialConditions.hh"
#include "ExcitationEnergy.hh"
//#include "DeexcitationHandler.hh"
#include "GRATEPhysicsList.hh"
#include "TGlauber/TGlauberMC.hh"
#include "TGlauber/TGlauNucleon.hh"
#include "TGlauber/TGlauNucleus.hh"

#include "include/GMSTClustering.hh"
*/
#include "GlauberCollisionReader.hh"
#include "FermiMomentum.hh"
#include "AAMCC.hh"
#include "WriteToFile.hh"

#include <fstream>


int main()
{
    //Seting parameters for Deexcitation
    G4NuclearLevelData* fLevelData = G4NuclearLevelData::GetInstance(); 
    G4DeexPrecoParameters* fParam = fLevelData->GetParameters();
    //fParam->SetDeexModelType(0);
    fParam->SetMinExPerNucleounForMF(4.*MeV); 
    fParam->SetDeexChannelsType(fEvaporation); 

    G4StateManager* fStateManager = G4StateManager::GetStateManager();

    // G4Random::setTheEngine(new CLHEP::RanluxEngine());

    //auto seed1 = 0x7ffce81312f0;
    
    //CLHEP::HepRandom::setTheSeed(seed1);

    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    CLHEP::HepRandom::setTheSeed((unsigned)clock());
    // CLHEP::HepRandom::setTheSeeds(CLHEP::HepRandom::getTheSeeds());

    CLHEP::RandFlat randFlat(new CLHEP::RanecuEngine);

    G4RunManager * runManager = new G4RunManager;
    runManager->SetUserInitialization(new GRATEPhysicsList);
    G4BosonConstructor pCBos;
    pCBos.ConstructParticle();

    G4LeptonConstructor pCLept;
    pCLept.ConstructParticle();

    G4MesonConstructor pCMes;
    pCMes.ConstructParticle();

    G4BaryonConstructor pCBar;
    pCBar.ConstructParticle();

    G4IonConstructor pCIon;
    pCIon.ConstructParticle();

    G4GenericIon* gion = G4GenericIon::GenericIon();
    gion->SetProcessManager(new G4ProcessManager(gion));

    G4StateManager::GetStateManager()->SetNewState(G4State_Init); // To let create ions
    G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
    G4IonTable* ions = partTable->GetIonTable();
    partTable->SetReadiness();
    ions->CreateAllIon();
    ions->CreateAllIsomer();

    // The user will be asked for the nuclear name to simulate it's collisions.
    GRATEmanager histoManager;
    AAMCCEvent event;
    void WriteToFile(AAMCCEvent*, AAMCCrun*, NucleonVector*); //Function that writes the data to the file

    //Get Z and A of nuclei
    G4int sourceA = histoManager.GetInitialContidions().GetSourceA();
    G4int sourceAb = histoManager.GetInitialContidions().GetSourceAb();
    G4double rA_plus_rB = 1.3*(pow(G4double(sourceA),1./3.)+pow(G4double(sourceAb),1./3.));

    //START OF THE MODELLING
    //Setting up GMST
    GMSTClustering* clusters = new GMSTClustering(histoManager.GetCriticalDistance(),sourceA,sourceAb);
    clusters->SetCD(histoManager.GetCriticalDistance());

    //Setting up ExcitationHandler
    DeexcitationHandler* handlerNew = new DeexcitationHandler();
    //handlerNew->SetMinEForMultiFrag(3*MeV);
    handlerNew->SetMaxAandZForFermiBreakUp(19, 9);
    handlerNew->SetMinExcitation(1e-4*MeV);
    handlerNew->SetMaxAforFermiBreakUp(19);
    handlerNew->SetMaxZforFermiBreakUp(9);
    handlerNew->SetMaxAforPureNeutronFragments(200);
    handlerNew->SetMinExForFermiBreakUp(0.01*MeV);
    handlerNew->SetExForMF(3*MeV, 5*MeV);
    //Setting up Glauber code
    //histoManager.CalcXsectNN();
    G4float omega = -1;
    G4float signn = histoManager.GetXsectNN();
    auto seed = static_cast<unsigned long int>(*CLHEP::HepRandom::getTheSeeds()); //setting the same seed to TGlauber

    TGlauberMC *mcg=new TGlauberMC(histoManager.GetInitialContidions().GetSysA(),histoManager.GetInitialContidions().GetSysB(),signn,omega,seed);
    mcg->SetMinDistance(0.4);
    //mcg->SetNodeDistance(0);
    mcg->SetCalcLength(0);
    mcg->SetCalcArea(0);
    mcg->SetCalcCore(0);
    mcg->SetDetail(99);

    if((histoManager.GetUpB() > 0 && histoManager.GetLowB() >= 0) && (histoManager.GetUpB() > histoManager.GetLowB()) && histoManager.GetUpB() < rA_plus_rB+7.5){
        mcg->SetBmin(histoManager.GetLowB());
        mcg->SetBmax(histoManager.GetUpB());
    } 
    else {std::cout << " \n------>Calculation will be MB \n"<<std::endl;}
 
    //Setting Excitation Energy
    ExcitationEnergy* ExEnA = new ExcitationEnergy(histoManager.GetStatType(), sourceA);
    ExcitationEnergy* ExEnB = new ExcitationEnergy(histoManager.GetStatType(), sourceAb);

    if(histoManager.GetStatType()> 2){
        //Parameters for ALADIN parametrization
        G4double e_0=11.5*MeV;//was 8.13 MeV
        G4double sigma0 = 0.005; //was 0.01
        G4double b0 = 2; // From Bondorf 1995
        G4double sigmaE0 = 1*MeV;
        G4double c0 = 0.1;
        G4double Pe = 24*MeV;
        G4double Pm = 0.2;
        ExEnA->SetParametersALADIN(e_0, sigma0, b0);
        ExEnB->SetParametersALADIN(e_0, sigma0, b0);
        ExEnA->SetParametersParabolicApproximation(Pe, Pm, sigma0, c0, 0.01);
        ExEnB->SetParametersParabolicApproximation(Pe, Pm, sigma0, c0, 0.01);
        //ExEnA->SetParametersCorrectedALADIN(0.01,1000,sigma0,c0,0);
        //ExEnB->SetParametersCorrectedALADIN(0.01,1000,sigma0,c0,0);
        ExEnA->SetParametersHybridFit(11.46648905*MeV, -1.84830078*MeV,  -58.53674677*MeV,  284.66431513*MeV, -637.51406293*MeV,  652.80324427*MeV, -251.28205381*MeV, 0.4*MeV, 0.5, 0.2);
        ExEnB->SetParametersHybridFit(11.46648905*MeV, -1.84830078*MeV,  -58.53674677*MeV,  284.66431513*MeV, -637.51406293*MeV,  652.80324427*MeV, -251.28205381*MeV, 0.4*MeV, 0.5, 0.2);
    }

    NucleonVector nV;
    GlauberCollisionReader reader;
    FermiMomentum FermiMom(&nV, "G");
    FermiMom.SetPzPerNucleon(histoManager.GetInitialContidions().GetPzA()/ sourceA, histoManager.GetInitialContidions().GetPzB() / sourceAb);

    for(G4int count=0;count<histoManager.GetIterations() ;count++){
        event.id = count;
        //An event generated by GlauberMC is here.
        mcg->Run(1);

        //histoManager.CalcNucleonDensity(mcg->GetNucleons(), mcg->GetB()); //newData
        //Side A $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        TGlauNucleus *nucA   = mcg->GetNucleusA(); 
        TGlauNucleus *nucB = mcg->GetNucleusB();
        event.NpartA = mcg->GetNpartA();
        event.NpartB = mcg->GetNpartB();
        event.Npart = mcg->GetNpart();
        event.b = mcg->GetB();
        event.PhiRotA = nucA->GetPhiRot();
        event.ThetaRotA = nucA->GetThetaRot();
        event.PhiRotB = nucB->GetPhiRot();
        event.ThetaRotB = nucB->GetThetaRot();
        for(int k = 0; k<10; k++){event.Ecc[k] = mcg->GetEcc(k);}
        event.Nhard = mcg->GetNhard();
        event.Ncoll = mcg->GetNcoll();
        event.Ncollpp = mcg->GetNcollpp();
        event.Ncollpn = mcg->GetNcollpn();
        event.Ncollnn = mcg->GetNcollnn();
        TObjArray* nucleons=mcg->GetNucleons();

        G4int A = 0;
        G4int Z = 0;
        G4int Ab = 0;
        G4int Zb = 0;

        nV = reader.GetNucleons(nucleons);
        Z = nV.GetZ("A"); A = nV.GetA("A"); Zb = nV.GetZ("B"); Ab = nV.GetA("B");

        if(!(A == 0 && Ab ==0)){
            G4int thisEventNumFragments = 0;
            std::cout.setf(std::ios::scientific, std::ios::floatfield);
            CLHEP::RandGauss randGauss(0, 1);

            // Excitation energy is calculated only for unstable clusters
            G4double energy_A = ExEnA->GetEnergy(A);
            G4double energy_B = ExEnB->GetEnergy(Ab);
            event.ExEnA = energy_A / G4double(A);
            event.ExEnB = energy_B/G4double(Ab);
            //histoManager.GetHisto2(1)->Fill(event.ExEnA, G4double(A) / sourceA); //newData


            CLHEP::HepLorentzVector Fermi4MomA = FermiMom.GetLorentzVector("A");
            CLHEP::Hep3Vector boostA = Fermi4MomA.boostVector();
            event.FermiMomA_x = Fermi4MomA.px();
            event.FermiMomA_y = Fermi4MomA.py();
            event.FermiMomA_z = Fermi4MomA.pz();
            CLHEP::HepLorentzVector Fermi4MomB = FermiMom.GetLorentzVector("B");;
            CLHEP::Hep3Vector boostB = Fermi4MomB.boostVector();
            event.FermiMomB_x = Fermi4MomB.px();
            event.FermiMomB_y = Fermi4MomB.py();
            event.FermiMomB_z = Fermi4MomB.pz();

            //if(sourceA - A == 1 ) histoManager.GetHisto2(8)->Fill(event.FermiMomA_x, event.FermiMomA_y);//newData


            std::vector<G4FragmentVector> MstClustersVector = clusters->GetClusters(&nV, energy_A, energy_B, boostA, boostB); //d = const if energy is negative

            event.d_MstA = clusters->GetCD("A");
            event.d_MstB = clusters->GetCD("B");

            //histoManager.GetHisto2(7)->Fill(event.ExEnA , event.d_MstA);//newData

            for(G4int i = 0; i < MstClustersVector.at(0).size(); ++i) {

                G4Fragment aFragment = (*MstClustersVector.at(0).at(i));
                G4LorentzVector p4 = aFragment.GetMomentum();
                //if((aFragment.GetMomentum().m() - G4NucleiProperties::GetNuclearMass(aFragment.GetA(), aFragment.GetZ()) - ExEnA*aFragment.GetA() !=0) && aFragment.GetA() != 1){std::cout<<"dE_x = "<<(aFragment.GetMomentum().m() - G4NucleiProperties::GetNuclearMass(aFragment.GetA(), aFragment.GetZ()))/aFragment.GetA() - ExEnA<<"\n";}
                event.A_cl.push_back(aFragment.GetA());
                event.Z_cl.push_back(aFragment.GetZ());
                G4double eta_A = 0;
                //if(abs(p4.px()) < 1) std::cout<<G4double(clfrag_A)<<" "<<G4double(sourceA)<<"\n";


                // HANDLER
                G4ReactionProductVector * theProduct = handlerNew->BreakUp(aFragment, histoManager.GetDeexModel());

                thisEventNumFragments = theProduct->size();

                //histoManager.GetHisto(1)->Fill(thisEventNumFragments);//newData

                for (G4ReactionProductVector::iterator iVector = theProduct->begin(); iVector != theProduct->end(); ++iVector) {
                    G4double thisFragmentZ = 0;
                    G4double thisFragmentA = 0;

                    const G4ParticleDefinition *pd = (*iVector)->GetDefinition();

                    G4String particleEmitted = pd->GetParticleName();

                    if (particleEmitted != "gamma" && particleEmitted != "e-" && particleEmitted != "e+") {
                        thisFragmentZ = pd->GetAtomicNumber();
                        thisFragmentA = pd->GetAtomicMass();
                        if (pd->GetAtomicMass() == 0) { G4cout << "ERROR, pn = " << pd->GetParticleName() << G4endl;}
                        event.MassOnSideA.push_back(thisFragmentA);
                        event.ChargeOnSideA.push_back(thisFragmentZ);

                        G4double eeA = (*iVector)->GetTotalEnergy();
                        G4LorentzVector product_p4((*iVector)->GetMomentum().x(), (*iVector)->GetMomentum().y(),
                        (*iVector)->GetMomentum().z(), eeA);
                        G4double pXonA = product_p4.x() / MeV;
                        G4double pYonA = product_p4.y() / MeV;
                        G4double pZonA = product_p4.z() / MeV;
                        p4 = p4 - product_p4;

                        eta_A = 0.5 * log((std::sqrt(pXonA * pXonA + pYonA * pYonA + pZonA * pZonA) + pZonA) /
                        (std::sqrt(pXonA * pXonA + pYonA * pYonA + pZonA * pZonA) - pZonA));

                        event.pXonSideA.push_back(pXonA); //if(abs(pXonA) < 1) std::cout<<"thisFragmentZ = "<<thisFragmentZ<<"\n";
                        event.pYonSideA.push_back(pYonA);
                        event.pZonSideA.push_back(pZonA);
                        event.pseudorapidity_A.push_back(eta_A);

                        //newData
                        /*
                        if (thisFragmentZ == 0) {
                            histoManager.GetHisto2(3)->Fill(pXonA, pYonA);
                            histoManager.GetHisto(2)->Fill(pZonA);
                        } else if (thisFragmentZ == 1 && thisFragmentA == 1) {
                            histoManager.GetHisto2(4)->Fill(pXonA, pYonA);
                            histoManager.GetHisto(3)->Fill(pZonA);
                        } else if (thisFragmentZ < 20 && thisFragmentZ > 2) {
                            histoManager.GetHisto2(5)->Fill(pXonA, pYonA);
                            histoManager.GetHisto(4)->Fill(pZonA);
                        } else {
                            histoManager.GetHisto2(6)->Fill(pXonA, pYonA);
                            histoManager.GetHisto(5)->Fill(pZonA);
                        }*/
                    }


                    //newData  histoManager.GetHisto(6)->Fill(thisFragmentZ);
                    //newData histoManager.GetHisto(7)->Fill(thisFragmentA);
                    //newData  histoManager.GetHisto2(2)->Fill(thisFragmentZ, thisFragmentA);
                    delete (*iVector);
                }
                //if(p4.mag() > 0.01) std::cout<<"p4.mag() = "<<p4.mag()<<" b = "<<b<<"\n";
                delete theProduct;
            }

            event.ClustNumA = event.A_cl.size();

//Side B $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

            for(G4int i = 0; i < MstClustersVector.at(1).size(); ++i) {
                G4Fragment aFragmentB = (*MstClustersVector.at(1).at(i));
                G4LorentzVector p4b = aFragmentB.GetMomentum();
                event.Ab_cl.push_back(aFragmentB.GetA());
                event. Zb_cl.push_back(aFragmentB.GetZ());
                G4double eta_B = 0;

                // HANDLER
                G4ReactionProductVector *theProductB = handlerNew->BreakUp(aFragmentB, histoManager.GetDeexModel());

                for (G4ReactionProductVector::iterator kVector = theProductB->begin(); kVector != theProductB->end(); ++kVector) {
                    G4int thisFragmentZb = 0;
                    G4int thisFragmentAb = 0;

                    const G4ParticleDefinition *pdB = (*kVector)->GetDefinition();

                    G4String particleEmittedB = pdB->GetParticleName();

                    if (particleEmittedB != "gamma" && particleEmittedB != "e-" && particleEmittedB != "e+") {
                        thisFragmentZb = pdB->GetAtomicNumber();
                        thisFragmentAb = pdB->GetAtomicMass();
                        event.MassOnSideB.push_back(thisFragmentAb);
                        event.ChargeOnSideB.push_back(thisFragmentZb);

                        G4double eeB = (*kVector)->GetTotalEnergy();
                        G4LorentzVector product_p4b((*kVector)->GetMomentum().x(), (*kVector)->GetMomentum().y(), (*kVector)->GetMomentum().z(), eeB);
                        G4double pXonB = product_p4b.x() / MeV;
                        G4double pYonB = product_p4b.y() / MeV;
                        G4double pZonB = product_p4b.z() / MeV;
                        p4b = p4b - product_p4b;

                        eta_B = 0.5 * log((std::sqrt(pXonB * pXonB + pYonB * pYonB + pZonB * pZonB) + pZonB) /
                        (std::sqrt(pXonB * pXonB + pYonB * pYonB + pZonB * pZonB) - pZonB));
                        event.pXonSideB.push_back(pXonB);
                        event.pYonSideB.push_back(pYonB);
                        event.pZonSideB.push_back(pZonB);
                        event.pseudorapidity_B.push_back(eta_B);
                    }
                    //newData  histoManager.GetHisto(0)->Fill(thisFragmentZb);
                    delete (*kVector);
                }
                delete theProductB;
            }

            event.ClustNumB = event.Ab_cl.size();


           histoManager.SetXsectTot(mcg->GetTotXSect());
           histoManager.ToFile(&event, &nV, &WriteToFile);

            event.A_cl.clear();
            event.Z_cl.clear();
            event.Ab_cl.clear();
            event.Zb_cl.clear();

            event.MassOnSideA.clear();
            event.MassOnSideB.clear();
            event.ChargeOnSideB.clear();
            event.ChargeOnSideA.clear();
            event.pXonSideA.clear();
            event.pXonSideB.clear();
            event.pYonSideA.clear();
            event.pYonSideB.clear();
            event.pZonSideA.clear();
            event.pZonSideB.clear();
            event.pseudorapidity_A.clear();
            event.pseudorapidity_B.clear();

            // Events calc info update
            if (!G4bool(count % 100)) { G4cout << "Program is working," << count << " events calculated    \r" << std::flush; }

        }
    }

    G4cout<<"----> collided "<<histoManager.GetIterations()<<" nuclei "<<histoManager.GetSysA()<< " with " << histoManager.GetSysB() <<" at N-N x-section "<<signn<<" mb"<<G4endl;

    if(!histoManager.ToFileOrNot()){
        G4cout<<"----> total x-sect = "<<mcg->GetTotXSect()<< " +- " << mcg->GetTotXSectErr() <<" b";
    }
    else
    {
        G4cout<<"----> Only 1 event, no tot x-sect";
    }

    delete runManager;
    delete clusters;
    delete handlerNew;
    delete mcg;
    delete ExEnA;
    delete ExEnB;
    return 0;
}
