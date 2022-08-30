#ifndef InitialConditions_h
#define InitialConditions_h 1


#include "AAMCC.hh"
#include "globals.hh"
#include "cmath"
#include "G4SystemOfUnits.hh"
#include <fstream>
#include "AAMCConstants.hh"
#include "Nucleon.hh"

class InitialConditions {
public:
	 InitialConditions();
    ~InitialConditions();
public:

    G4bool SetSysA(G4String SysA_in);
    G4bool SetSysB(G4String SysB_in);

    inline  void   SetCollider(G4bool IsCollider_in) {IsCollider = IsCollider_in;};
    void   SetKinematics(G4double Energy_in);

    void SetConditions(AAMCCinput cond_in);

    inline G4double GetKinEnergy(){return KinEn;};
    inline G4double GetSqrtSnn(){return SqrtSnn;};
    inline G4double GetPzA()      {return PzA;};
    inline G4double GetPzB()      {return PzB;};
    inline G4bool   GetCollider() {return  IsCollider;};

    inline G4int    GetSourceA()  {return sourceA;};
    inline G4int    GetSourceAb() {return sourceAb;};
    inline G4int    GetSourceZ()  {return sourceZ;};
    inline G4int    GetSourceZb() {return sourceZb;};

    inline G4String GetSysA()     {return  SysA;};
    inline G4String GetSysB()     {return  SysB;};

public:

    G4double GetXsectNN(); //in barn

private:
    inline void    SetSourceA(G4int A_in)  {sourceA = A_in;};
    inline void    SetSourceAb(G4int A_in) {sourceAb = A_in;};
    inline void    SetSourceZ(G4int Z_in)  {sourceZ = Z_in;};
    inline void    SetSourceZb(G4int Z_in) {sourceZb = Z_in;};

    G4double KinEn;
    G4double SqrtSnn;
    G4double PzA;
    G4double PzB;
    G4bool   IsCollider{false};

    G4int sourceA;
    G4int sourceAb;
    G4int sourceZ;
    G4int sourceZb;

    G4String SysA;
    G4String SysB;

    G4double XsectNN = -1;
};



#endif