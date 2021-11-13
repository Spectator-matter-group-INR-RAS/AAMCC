#include "GaimardSchmidt.hh"

G4double GaimardSchmidt(G4double E,G4double EvaporationEnergy, G4int a,G4int A){ 
  G4double g0=16.;
  G4double g1=0.7;
  G4int RemovedNucleons = A-a;
  G4double res = 0.;
  G4double testRes = 0.;

  if(E>0.5*EvaporationEnergy*RemovedNucleons){ //Energy restriction for easier caclulations
    return res; 
  }
  else{

   //First term calculation at H-S formula
   G4double term=g0;

   for(G4int k = 2; k < RemovedNucleons-1; k++){

     term*=(g0*E)/G4double(k*(k-1));
  
    }

    res=term;

    testRes = term;
    //Rest sum calculating with recursive form of H-S formula
    for(G4int m = 1; m<RemovedNucleons; m++){
     term *=(g1*E/g0)*((-1)*G4double(RemovedNucleons-m+1)/(G4double(m)*G4double(RemovedNucleons+m-1)));
    //G4cout<<testRes<<G4endl; //NaN is result of inf-inf becouse term becomes inf at some step.
    res+=term;
    } 
    
   if(res!=res){res=0;} //- костыль!!!
   if(res<0){res=0;}

   return res; 

 }
}
 	
