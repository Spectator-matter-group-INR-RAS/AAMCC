//
// Created by apbus_amp_k on 12.11.22.
//

#include "G4LorentzVector.hh"

#include "MCIniReader.hh"

#include "TClonesArray.h"

MCIniReader::MCIniReader(MCIniReader && rha) noexcept  :
    ftree(rha.ftree), curr(nullptr), curr_st(nullptr), asum(rha.asum), iter(rha.iter), treen(rha.treen) {
    rha.ftree = nullptr;
    ftree->SetBranchAddress("event", &curr);
    //ftree->SetBranchAddress("iniState", &curr_st); //TODO:make it back when iniState is back
}

MCIniReader& MCIniReader::operator=(MCIniReader && rha) noexcept {
    if (this == &rha)
        return *this;
    ftree = rha.ftree;
    rha.ftree = nullptr;
    curr = nullptr;
    rha.curr = nullptr;
    curr_st = nullptr;
    rha.curr_st = nullptr;
    asum = rha.asum;
    iter = rha.iter;
    treen = rha.treen;
    ftree->SetBranchAddress("event", &curr);
    //ftree->SetBranchAddress("iniState", &curr_st); //TODO:make it back when iniState is back
}

MCIniReader::MCIniReader(const std::unique_ptr<TFile>& tfile) : curr(nullptr), curr_st(nullptr), iter(0) {
    tfile->GetObject("events", ftree);          // Legacy weirdness. ptr has to be lvalue thus preventing smart pointer usage
    treen = ftree->GetEntries();
    ftree->SetBranchAddress("event", &curr);
    //ftree->SetBranchAddress("iniState", &curr_st); //TODO:make it back when iniState is back
    URun* run_data;
    tfile->GetObject("run", run_data);// legacy
    aa = run_data->GetAProj();
    ab = run_data->GetATarg();
    asum = run_data->GetATarg() + run_data->GetAProj();
}

AAMCCinput MCIniReader::operator()() {
//    if (iter >= treen)
//        throw std::out_of_range("No more events in file!");
    ftree->GetEntry(iter);
    iter++;
    AAMCCinput cache;
    auto arr = curr->GetParticleList();
    int size = curr->GetNpa();
    for (int i = 0; i < size; i++) {
        auto particle = (UParticle*) arr->At(i);
        //new nucleon to write into cache->nucleons
        //std::cout<<particle->GetIndex()<< "\r" <<std::flush; //errors with indexing leads to out of range
        //if(particle->GetIndex() > asum || (particle->GetIndex() < asum && particle->GetIndex() != 0 ? curr_st->getNucleon(particle->GetIndex()).getCollisionType() : 1000)  > 0) //TODO: Add collision type check nucleon->collisionType == 0 (no collision), 1 (el. collision with initial nucleus), 2 (el. collision with produced particle), 3 (nonel. with init nucl), 4 (nonel. with produced)
        if((particle->GetStatus() != 0 || particle->GetParent() != 0))
            continue;
        aamcc::Nucleon nucl;
        switch (particle->GetPdg()) {
            case 2212:                  //proton pdg code
                nucl.isospin = true;
                break;
            case 2112:                  //neutron pdg code
                nucl.isospin = false;
                break;
            default:
                continue;               //break is unnecessary
        }
        nucl.isParticipant = false;
        nucl.x = particle->X();
        nucl.y = particle->Y();         //coords
        nucl.z = particle->Z();
       if(particle->Pz() > 0) {        //proj
        //if(particle->GetIndex() < aa) {        //proj
            nucl.Nucl = "A";
            cache.FermiMomA_x += particle->Px() * GeV;  //MCIni stores energy/impulses in GeV
            cache.FermiMomA_y += particle->Py() * GeV;
            cache.FermiMomA_z += particle->Pz() * GeV;
        }
        else {                          //targ
            nucl.Nucl = "B";
            cache.FermiMomB_x += particle->Px() * GeV;
            cache.FermiMomB_y += particle->Py() * GeV;
            cache.FermiMomB_z += particle->Pz() * GeV;
        }
        cache.nucleons.push_back(std::move(nucl));     //write the nucleon data
    }
    return cache;
}
