#include "GlauberCollisionReader.hh"

void GlauberCollisionReader::Read(TObjArray* nucleons_in) {
    data = make_unique<AAMCCinput>();
    aamcc::Nucleon nucl;
    for(int iArray = 0; iArray < nucleons_in->GetEntries(); iArray++){
        auto nucleon = (TGlauNucleon*)(nucleons_in->At(iArray));        //legacy
        nucl.x = nucleon->GetX();
        nucl.y = nucleon->GetY();
        nucl.z = nucleon->GetZ();
        nucl.isospin = nucleon->IsProton();
        nucl.isParticipant = nucleon->IsWounded();
        if(nucleon->IsInNucleusA()){nucl.Nucl = "A";} else {nucl.Nucl = "B";};
        data->nucleons.push_back(nucl);
        nucl.Clean();
    }
}

