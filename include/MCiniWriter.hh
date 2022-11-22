#ifndef GRATE_MCINIWRITER_HH
#define GRATE_MCINIWRITER_HH

#include "VWriter.hh"
#include "G4SystemOfUnits.hh"
#include "TTree.h"
#include "TFile.h"
#include "URun.h"
#include "UEvent.h"
#include "UParticle.h"
#include "EventInitialState.h"
#include <mutex>
#include <memory>

class MCINIWriter : public VWriter{
public:
    MCINIWriter();
    ~MCINIWriter();
    void operator()(AAMCCEvent* ev, AAMCCrun* run, aamcc::NucleonVector* nucleons);
    void SetBMin(double bMin_in){bMin = bMin_in;}
    void SetBMax(double bMax_in){bMax = bMax_in;}
    void SetPhiMin(double phi){phiMin = phi;}
    void SetPhiMax(double phi){phiMax = phi;}
private:
    std::shared_ptr<TTree> tMCini;
    std::shared_ptr<URun> runMCini;

    double_t bMin = 0;
    double_t bMax = 0;
    Int_t distParam = 0;
    double_t phiMin = 0;
    double_t phiMax = 2*3.141592;

    AAMCCEvent event;
    AAMCCrun runData;

    std::unique_ptr<UEvent>uevent  = std::unique_ptr<UEvent>(new UEvent());
    std::unique_ptr<EventInitialState> iniState = std::unique_ptr<EventInitialState>(new EventInitialState());

    void InitTree(){
        std::shared_ptr<TTree> tmc(new TTree("events", "AAMCC"));
        tMCini = tmc;
        tMCini->SetDirectory(0);
        tMCini->Branch("event","UEvent",uevent.get());
        tMCini->Branch("iniState", "EventInitialState", iniState.get());
    }

    std::unique_ptr<URun> GetURun(TTree* tTree, AAMCCrun run){
        std::unique_ptr<URun> uRun( new URun("AAMCC", "", run.AinitA, run.ZinitA, run.pzA, run.AinitB, run.ZinitB, run.pzB, bMin,bMax, distParam, phiMin,phiMax, run.XsectTot, tTree->GetEntries()));
        return uRun;
    }

    std::vector<Nucleon> NucleonsToNucleons(aamcc::NucleonVector nucleonVector){
        std::vector <Nucleon> nucleons;

        return nucleons;
    }

    void FillTree(AAMCCEvent* ev){
        uevent->Clear();
        iniState->Clear();
        uevent->SetParameters(ev->id, ev->b, ev->PhiRotA, 1, 0,0);

        int partid = 0;
        Int_t child[2] = { 0,0 };

        for(int k = 0; k<int(ev->MassOnSideA.size()); ++k){
            ++partid;
            double energy = std::sqrt(std::pow(std::pow(10,-3)*ev->pXonSideA.at(k),2) + std::pow(std::pow(10,-3)*ev->pYonSideA.at(k),2) + std::pow(std::pow(10,-3)*ev->pZonSideA.at(k),2) + std::pow(ev->MassOnSideA.at(k)*nucleonAverMass,2));
            int pdg = aamcc::IsotopeToPDG(ev->ChargeOnSideA.at(k), ev->MassOnSideA.at(k));
            uevent->AddParticle(partid, pdg, 0,0,0,0,0,child,std::pow(10,-3)*ev->pXonSideA.at(k),std::pow(10,-3)*ev->pYonSideA.at(k),std::pow(10,-3)*ev->pZonSideA.at(k), energy, 0,0,0,1,1.);
        }
        for(int k = 0; k<int(ev->MassOnSideB.size()); ++k){
            ++partid;
            double energy = std::sqrt(std::pow(std::pow(10,-3)*ev->pXonSideB.at(k),2) + std::pow(std::pow(10,-3)*ev->pYonSideB.at(k),2) + std::pow(std::pow(10,-3)*ev->pZonSideB.at(k),2) + std::pow(ev->MassOnSideB.at(k)*nucleonAverMass,2));
            int pdg = aamcc::IsotopeToPDG(ev->ChargeOnSideB.at(k), ev->MassOnSideB.at(k));
            uevent->AddParticle(partid, pdg, 0,0,0,0,0,child,std::pow(10,-3)*ev->pXonSideB.at(k),std::pow(10,-3)*ev->pYonSideB.at(k),std::pow(10,-3)*ev->pZonSideB.at(k), energy, 0,0,0,1,1.);
        }
        iniState->setNColl(ev->Ncoll);
        iniState->setNPart(ev->Npart);
        tMCini->Fill();
    }

};


#endif //GRATE_MCINIWRITER_HH
