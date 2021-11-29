#include "Ericson.hh"

G4double Ericson(G4double E,G4double EvaporationEnergy, G4int a,G4int A){ 
  G4double g0=16; //16 was in Shidenberger - not influence calculations at all cause its freeze out in distr normalization. 
  G4int RemovedNucleons = A-a;

if(E>EvaporationEnergy*(A-a)){
   G4double s=0;  
  return s; 
}
else{

G4double s=g0;

for(G4int k = 2; k < RemovedNucleons; k++){

s*=(g0*E)/G4double(k*(k-1));

}


 return s; 
}
 }

