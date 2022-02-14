#include "InitialConditions.hh"

InitialConditions::InitialConditions() {

}

InitialConditions::~InitialConditions() {

}

G4bool InitialConditions::SetSysA(G4String SysA_in) {
    if(SysA_in == "Pb"){ sourceA = 208; sourceZ = 82; SysA = SysA_in; SysA+="*";}
    else if(SysA_in == "Pbpn"){sourceA = 208; sourceZ = 82; SysA = SysA_in;}
    else if(SysA_in == "Pbrw"){sourceA = 208; sourceZ = 82; SysA = SysA_in;}
    else if(SysA_in == "Pbpnrw"){sourceA = 208; sourceZ = 82; SysA = SysA_in;}
    else if(SysA_in == "Cu"){sourceA = 64; sourceZ = 29; SysA = SysA_in; SysA+="2";}
    else if(SysA_in == "O"){sourceA = 16; sourceZ = 8; SysA = SysA_in;}
    else if(SysA_in == "O2"){sourceA = 16; sourceZ = 8; SysA = SysA_in;}
    else if(SysA_in == "Oho"){sourceA = 16; sourceZ = 8; SysA = SysA_in;}
    else if(SysA_in == "Au"){sourceA = 197; sourceZ = 79; SysA = SysA_in;}
    else if(SysA_in == "Aurw"){sourceA = 197; sourceZ = 79; SysA = SysA_in;}
    else if(SysA_in == "Au2"){sourceA = 197; sourceZ = 79; SysA = SysA_in;}
    else if(SysA_in == "Au2rw"){sourceA = 197; sourceZ = 79; SysA = SysA_in;}
    else if(SysA_in == "Ag"){sourceA = 109; sourceZ = 47; SysA = SysA_in;}
    else if(SysA_in == "Br"){sourceA = 79; sourceZ = 35; SysA = SysA_in;}
    else if(SysA_in == "Xe"){sourceA = 129; sourceZ = 54; SysA = SysA_in;}
    else if(SysA_in == "Ar"){sourceA = 40; sourceZ = 18; SysA = SysA_in;}
    else if(SysA_in == "Ca2"){sourceA = 40; sourceZ = 20; SysA = SysA_in;}
    //else if(SysA_in == "Ni"){sourceA = 58; sourceZ = 28; SysA = SysA_in;} //issue on GlauberMC side
    else if(SysA_in == "C"){sourceA = 12; sourceZ = 6; SysA = SysA_in;}
    else if(SysA_in == "Al"){sourceA = 27; sourceZ = 13; SysA = SysA_in;}
    else if(SysA_in == "U"){sourceA = 238; sourceZ = 92; SysA = SysA_in;}
    else if(SysA_in == "U2"){sourceA = 238; sourceZ = 92; SysA = SysA_in;}
    else if(SysA_in == "He4"){sourceA = 4; sourceZ = 2; SysA = SysA_in;}
    else if(SysA_in == "He3"){sourceA = 3; sourceZ = 2; SysA = SysA_in;}
    else if(SysA_in == "H3"){sourceA = 3; sourceZ = 1; SysA = SysA_in;}
    else if(SysA_in == "d"){sourceA = 2; sourceZ = 1; SysA = SysA_in;}
    else if(SysA_in == "p"){sourceA = 1; sourceZ = 1; SysA = SysA_in;}
    else{ G4Exception("Nucleus input in GRATE", "GRATE-0", JustWarning, "There is no matched nucleus in GRATE");
        return 0;
    }

    return 1;
}

G4bool InitialConditions::SetSysB(G4String SysB_in) {
    if(SysB_in == "Pb"){sourceAb = 208; sourceZb = 82; SysB = SysB_in; SysB+="*";}
    else if(SysB_in == "Pbpn"){sourceAb = 208; sourceZb = 82; SysB = SysB_in;}
    else if(SysB_in == "Pbrw"){sourceAb = 208; sourceZb = 82; SysB = SysB_in;}
    else if(SysB_in == "Pbpnrw"){sourceAb = 208; sourceZb = 82; SysB = SysB_in;}
    else if(SysB_in == "Cu"){sourceAb = 64; sourceZb = 29; SysB = SysB_in; SysB+="2";}
    else if(SysB_in == "O") {sourceAb = 16; sourceZb = 8; SysB = SysB_in;}
    else if(SysB_in == "O2") {sourceAb = 16; sourceZb = 8; SysB = SysB_in;}
    else if(SysB_in == "Oho") {sourceAb = 16; sourceZb = 8; SysB = SysB_in;}
    else if(SysB_in == "Au"){sourceAb = 197; sourceZb = 79; SysB = SysB_in;}
    else if(SysB_in == "Aurw"){sourceAb = 197; sourceZb = 79; SysB = SysB_in;}
    else if(SysB_in == "Au2"){sourceAb = 197; sourceZb = 79; SysB = SysB_in;}
    else if(SysB_in == "Au2rw"){sourceAb = 197; sourceZb = 79; SysB = SysB_in;}
    else if(SysB_in == "Ag"){sourceAb = 109; sourceZb = 47; SysB = SysB_in;}
    else if(SysB_in == "Br"){sourceAb = 79; sourceZb = 35; SysB = SysB_in;}
    else if(SysB_in == "Xe"){sourceAb = 129; sourceZb = 54; SysB = SysB_in;}
    else if(SysB_in == "Ar"){sourceAb = 40; sourceZb = 18; SysB = SysB_in;}
    else if(SysB_in == "Ca2") {sourceAb = 40; sourceZb = 20; SysB = SysB_in;}
    //else if(SysB_in == "Ni") {sourceAb = 58; sourceZb = 28; SysB = SysB_in;} //issue on GlauberMC side
    else if(SysB_in == "C") {sourceAb = 12; sourceZb = 6; SysB = SysB_in;}
    else if(SysB_in == "Al"){sourceAb = 27; sourceZb = 13; SysB = SysB_in;}
    else if(SysB_in == "U") {sourceAb = 238; sourceZb = 92; SysB = SysB_in;}
    else if(SysB_in == "U2"){sourceAb = 238; sourceZb = 92; SysB = SysB_in;}
    else if(SysB_in == "He4"){sourceAb = 4; sourceZb = 2; SysB = SysB_in;}
    else if(SysB_in == "He3"){sourceAb = 3; sourceZb = 2; SysB = SysB_in;}
    else if(SysB_in == "H3"){sourceAb = 3; sourceZb = 1; SysB = SysB_in;}
    else if(SysB_in == "d"){sourceAb = 2; sourceZb = 1; SysB = SysB_in;}
    else if(SysB_in == "p"){sourceAb = 1; sourceZb = 1; SysB = SysB_in;}
    else{ G4Exception("NuclPzAeus input in GRATE", "GRATE-0", JustWarning, "There is no matched nucleus in GRATE");
        return 0;
    }

    return 1;
}

void InitialConditions::SetKinematics(G4double Energy_in) {

    if(IsCollider){
        PzA =    pow(Energy_in*Energy_in*0.25 - nucleonAverMass*nucleonAverMass,0.5);
        PzB = -1*pow(Energy_in*Energy_in*0.25 - nucleonAverMass*nucleonAverMass,0.5);
        KinEn = (Energy_in/2.0 - nucleonAverMass);
        SqrtSnn = Energy_in;
    } else{
        PzA = pow(Energy_in*(Energy_in+2*nucleonAverMass),0.5);
        PzB = 0;
        KinEn = Energy_in;
        SqrtSnn = pow(2*nucleonAverMass*nucleonAverMass+2*Energy_in*nucleonAverMass, 0.5);
    }
    SqrtSnn*= GeV;
    KinEn*= sourceA*GeV;
    PzA*= sourceA*GeV;
    PzB*= sourceAb*GeV;
}


