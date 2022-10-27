#include "GlauberCollisionReader.hh"

AAMCCinput GlauberCollisionReader::operator()() {
    return data;
}

void GlauberCollisionReader::Read(TObjArray *nucleons_in) {
    data.nucleons.clear();
    Nucleon nucl;
    nucleons = nucleons_in;
    for(int iArray = 0; iArray < nucleons->GetEntries(); iArray++){
        auto *nucleon=(TGlauNucleon*)(nucleons->At(iArray));
        nucl.x = nucleon->GetX();
        nucl.y = nucleon->GetY();
        nucl.z = nucleon->GetZ();
        nucl.isospin = nucleon->IsProton();
        nucl.isParticipant = nucleon->IsWounded();
        if(nucleon->IsInNucleusA()){nucl.Nucl = "A";} else {nucl.Nucl = "B";};
        data.nucleons.push_back(nucl);
        nucl.Clean();
    }
}

