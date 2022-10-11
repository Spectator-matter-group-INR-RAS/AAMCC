
#ifndef InitialConditions_h
#define InitialConditions_h 1

#include "globals.hh"
#include "cmath"
#include "G4SystemOfUnits.hh"
#include <fstream>
#include "AAMCConstants.hh"

class InitialConditions {
public:
	 InitialConditions();
    ~InitialConditions();
     InitialConditions(G4double KinEn_in, G4String SysA_in, G4String SysB_in, G4bool IsCollider_in);
public:

    G4bool SetSysA(G4String SysA_in);
    G4bool SetSysB(G4String SysB_in);
    inline  void   SetCollider(G4bool IsCollider_in) {IsCollider = IsCollider_in;};
    void   SetKinEn(G4double KinEn_in);

    inline G4double GetKinEnergy(){return KinEn;};
    inline G4double GetPzA()      {return PzA;};
    inline G4double GetPzB()      {return PzB;};
    inline G4bool   GetCollider() {return  IsCollider;};

    inline G4int    GetSourceA()  {return sourceA;};
    inline G4int    GetSourceAb() {return sourceAb;};
    inline G4int    GetSourceZ()  {return sourceZ;};
    inline G4int    GetSourceZb() {return sourceZb;};

    inline G4String GetSysA()     {return  SysA;};
    inline G4String GetSysB()     {return  SysB;};

private:
    G4double KinEn;
    G4double PzA;
    G4double PzB;
    G4bool   IsCollider = 0;

    G4int sourceA;
    G4int sourceAb;
    G4int sourceZ;
    G4int sourceZb;

    G4String SysA;
    G4String SysB;

};

#endif