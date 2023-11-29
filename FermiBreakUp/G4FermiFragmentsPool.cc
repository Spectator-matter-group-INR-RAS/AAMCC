//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//


#include "G4FermiFragmentsPool.hh"

G4bool G4FermiFragmentsPool::MapIsEmpty(true);


std::multimap<const std::pair<G4int,G4int>, const G4VFermiFragment* , std::less<const std::pair<G4int,G4int> > >  &
G4FermiFragmentsPool::GetMap()
{
  static std::vector<const G4VFermiFragment * > fragment_pool;
  //                                                             A  Z  Pol  ExcitE
  static const G4StableFermiFragment Fragment00(  1, 0,  2,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment00);
  static const G4StableFermiFragment Fragment01(  1, 1,  2,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment01);
  static const G4StableFermiFragment Fragment02(  2, 1,  3,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment02);
  static const G4StableFermiFragment Fragment03(  3, 1,  2,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment03);
  static const G4StableFermiFragment Fragment04(  3, 2,  2,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment04);
  static const G4StableFermiFragment Fragment05(  4, 2,  1,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment05);
  static const G4He5FermiFragment    Fragment06(  5, 2,  4, 16.76*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment06);// He5
  static const G4Li5FermiFragment    Fragment07(  5, 3,  4, 16.66*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment07);// Li5
  static const G4StableFermiFragment Fragment08(  6, 2,  1,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment08);
  static const G4StableFermiFragment Fragment09(  6, 3,  3,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment09);
  
  static const G4StableFermiFragment Fragment10(  6, 3,  1,  3.56*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment10);
  static const G4StableFermiFragment Fragment11(  7, 3,  4,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment11);
  static const G4StableFermiFragment Fragment12(  7, 3,  2,  0.48*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment12);
  static const G4StableFermiFragment Fragment13(  7, 4,  4,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment13);
  static const G4StableFermiFragment Fragment14(  7, 4,  2,  0.43*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment14);
  static const G4StableFermiFragment Fragment15(  8, 3,  5,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment15);
  static const G4StableFermiFragment Fragment16(  8, 3,  3,  0.98*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment16);
  static const G4Be8FermiFragment    Fragment17(  8, 4,  1,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment17); // Be8
  static const G4StableFermiFragment Fragment18(  9, 4,  4,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment18);
  static const G4B9FermiFragment     Fragment19(  9, 5,  4,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment19); // B9
  
  static const G4StableFermiFragment Fragment20( 10, 4,  1,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment20);
  static const G4StableFermiFragment Fragment21( 10, 4,  5,  3.37*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment21);
  static const G4StableFermiFragment Fragment22( 10, 4,  8,  5.96*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment22);
  static const G4StableFermiFragment Fragment23( 10, 4,  1,  6.18*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment23);
  static const G4StableFermiFragment Fragment24( 10, 4,  5,  6.26*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment24);
  static const G4StableFermiFragment Fragment25( 10, 5,  7,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment25);
  static const G4StableFermiFragment Fragment26( 10, 5,  3,  0.72*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment26);
  static const G4StableFermiFragment Fragment27( 10, 5,  1,  1.74*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment27);
  static const G4StableFermiFragment Fragment28( 10, 5,  3,  2.15*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment28);
  static const G4StableFermiFragment Fragment29( 10, 5,  5,  3.59*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment29);
  
  static const G4StableFermiFragment Fragment30( 10, 6,  3,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment30);
  static const G4StableFermiFragment Fragment31( 10, 6,  5,  3.35*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment31);
  static const G4StableFermiFragment Fragment32( 11, 5,  4,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment32);
  static const G4StableFermiFragment Fragment33( 11, 5,  2,  2.13*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment33);
  static const G4StableFermiFragment Fragment34( 11, 5,  6,  4.44*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment34);
  static const G4StableFermiFragment Fragment35( 11, 5,  4,  5.02*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment35);
  static const G4StableFermiFragment Fragment36( 11, 5, 10,  6.76*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment36);
  static const G4StableFermiFragment Fragment37( 11, 5,  6,  7.29*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment37);
  static const G4StableFermiFragment Fragment38( 11, 5,  4,  7.98*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment38);
  static const G4StableFermiFragment Fragment39( 11, 5,  6,  8.56*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment39);
  
  static const G4StableFermiFragment Fragment40( 11, 6,  4,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment40);
  static const G4StableFermiFragment Fragment41( 11, 6,  2,  2.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment41);
  static const G4StableFermiFragment Fragment42( 11, 6,  6,  4.32*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment42);
  static const G4StableFermiFragment Fragment43( 11, 6,  4,  4.80*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment43);
  static const G4StableFermiFragment Fragment44( 11, 6,  2,  6.34*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment44);
  static const G4StableFermiFragment Fragment45( 11, 6,  8,  6.48*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment45);
  static const G4StableFermiFragment Fragment46( 11, 6,  6,  6.90*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment46);
  static const G4StableFermiFragment Fragment47( 11, 6,  4,  7.50*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment47);
  static const G4StableFermiFragment Fragment48( 11, 6,  4,  8.10*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment48);
  static const G4StableFermiFragment Fragment49( 11, 6,  6,  8.42*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment49);
  
  static const G4StableFermiFragment Fragment50( 11, 6,  8,  8.66*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment50);
  static const G4StableFermiFragment Fragment51( 12, 5,  3,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment51);
  static const G4StableFermiFragment Fragment52( 12, 5,  5,  0.95*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment52);
  static const G4StableFermiFragment Fragment53( 12, 5,  5,  1.67*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment53);
  static const G4StableFermiFragment Fragment54( 12, 5,  4,  2.65*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment54);
  static const G4StableFermiFragment Fragment55( 12, 6,  1,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment55);
  static const G4StableFermiFragment Fragment56( 12, 6,  5,  4.44*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment56);
  static const G4StableFermiFragment Fragment57( 13, 6,  2,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment57);
  static const G4StableFermiFragment Fragment58( 13, 6,  2,  3.09*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment58);
  static const G4StableFermiFragment Fragment59( 13, 6,  4,  3.68*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment59);
  
  static const G4StableFermiFragment Fragment60( 13, 6,  6,  3.85*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment60);
  static const G4StableFermiFragment Fragment61( 13, 7,  2,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment61);
  static const G4StableFermiFragment Fragment62( 14, 6,  1,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment62);
  static const G4StableFermiFragment Fragment63( 14, 6,  3,  6.09*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment63);
  static const G4StableFermiFragment Fragment64( 14, 6,  8,  6.69*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment64);
  static const G4StableFermiFragment Fragment65( 14, 6,  6,  6.96*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment65);
  static const G4StableFermiFragment Fragment66( 14, 6,  5,  7.34*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment66);
  static const G4StableFermiFragment Fragment67( 14, 7,  3,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment67);
  static const G4StableFermiFragment Fragment68( 14, 7,  1,  2.31*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment68);
  static const G4StableFermiFragment Fragment69( 14, 7,  3,  3.95*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment69);
  
  static const G4StableFermiFragment Fragment70( 14, 7,  1,  4.92*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment70);
  static const G4StableFermiFragment Fragment71( 14, 7,  5,  5.11*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment71);
  static const G4StableFermiFragment Fragment72( 14, 7,  3,  5.69*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment72);
  static const G4StableFermiFragment Fragment73( 14, 7,  7,  5.83*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment73);
  static const G4StableFermiFragment Fragment74( 14, 7,  3,  6.20*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment74);
  static const G4StableFermiFragment Fragment75( 14, 7,  7,  6.44*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment75);
  static const G4StableFermiFragment Fragment76( 14, 7,  5,  7.03*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment76);
  static const G4StableFermiFragment Fragment77( 15, 7,  2,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment77);
  static const G4StableFermiFragment Fragment78( 15, 7,  8,  5.28*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment78);
  static const G4StableFermiFragment Fragment79( 15, 7,  4,  6.32*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment79);
  
  static const G4StableFermiFragment Fragment80( 15, 7, 10,  7.22*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment80);
  static const G4StableFermiFragment Fragment81( 15, 7,  8,  7.57*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment81);
  static const G4StableFermiFragment Fragment82( 15, 7,  2,  8.31*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment82);
  static const G4StableFermiFragment Fragment83( 15, 7,  4,  8.57*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment83);
  static const G4StableFermiFragment Fragment84( 15, 7, 14,  9.15*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment84);
  static const G4StableFermiFragment Fragment85( 15, 7, 14,  9.79*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment85);
  static const G4StableFermiFragment Fragment86( 15, 7,  8, 10.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment86);
  static const G4StableFermiFragment Fragment87( 15, 8,  2,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment87);
  static const G4StableFermiFragment Fragment88( 15, 8,  8,  5.22*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment88);
  static const G4StableFermiFragment Fragment89( 15, 8,  4,  6.18*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment89);
  
  static const G4StableFermiFragment Fragment90( 15, 8, 10,  6.83*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment90);
  static const G4StableFermiFragment Fragment91( 15, 8,  8,  7.28*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment91);
  static const G4StableFermiFragment Fragment92( 16, 7,  5,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment92);
  static const G4StableFermiFragment Fragment93( 16, 7,  1,  0.12*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment93);
  static const G4StableFermiFragment Fragment94( 16, 7,  7,  0.30*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment94);
  static const G4StableFermiFragment Fragment95( 16, 7,  3,  0.40*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment95);
  static const G4StableFermiFragment Fragment96( 16, 8,  1,  0.00*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment96);
  static const G4StableFermiFragment Fragment97( 16, 8,  8,  6.10*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment97);
  static const G4StableFermiFragment Fragment98( 16, 8,  5,  6.92*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment98);
  static const G4StableFermiFragment Fragment99( 16, 8,  3,  7.12*MeV ); if(MapIsEmpty) fragment_pool.push_back(&Fragment99);

  static std::multimap<const std::pair<G4int,G4int>, const G4VFermiFragment* , std::less<const std::pair<G4int,G4int> > >  
    theMapOfFragments;

  if (MapIsEmpty) 
    {
      for(size_t i=0; i<fragment_pool.size(); i++)
      {
        theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(fragment_pool[i]->GetA(),fragment_pool[i]->GetZ()),fragment_pool[i]));
      }
      MapIsEmpty = false;
    }

  return theMapOfFragments;
}


G4FermiFragmentsPool::G4FermiFragmentsPool()
{
}



G4FermiFragmentsPool::~G4FermiFragmentsPool()
{
}





G4FermiFragmentsPool::G4FermiFragmentsPool(const G4FermiFragmentsPool&)
{
  // It is meant to not be accesable
}

const G4FermiFragmentsPool & G4FermiFragmentsPool::operator=(const G4FermiFragmentsPool& )
{
  // It is meant to not be accesable
  return *this;
}

G4bool G4FermiFragmentsPool::operator==(const G4FermiFragmentsPool&) const
{
  // It is meant to not be accesable
  return false;
}

G4bool G4FermiFragmentsPool::operator!=(const G4FermiFragmentsPool&) const
{
  // It is meant to not be accesable
  return true;
}
