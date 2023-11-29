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
// $Id: G4Be8FermiFragment.cc,v 1.7 2006/06/29 20:12:46 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4Be8FermiFragment.hh"
#include "G4IonTable.hh"
#include "G4HadronicException.hh"
#include "G4PhysicalConstants.hh"

G4Be8FermiFragment::G4Be8FermiFragment()
{
}

G4Be8FermiFragment::G4Be8FermiFragment(const G4Be8FermiFragment &) : G4UnstableFermiFragment()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4Be8FermiFragment::copy_constructor meant to not be accessable");
}


G4Be8FermiFragment::~G4Be8FermiFragment()
{
}


const G4Be8FermiFragment & G4Be8FermiFragment::operator=(const G4Be8FermiFragment &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4Be8FermiFragment::operator= meant to not be accessable");
    return *this;
}


G4bool G4Be8FermiFragment::operator==(const G4Be8FermiFragment &) const
{
    return false;
}

G4bool G4Be8FermiFragment::operator!=(const G4Be8FermiFragment &) const
{
    return true;
}



G4Be8FermiFragment::G4Be8FermiFragment(const G4int anA, const G4int aZ, const G4int Pol, const G4double ExE)
  : G4UnstableFermiFragment(anA,aZ,Pol,ExE)
{
  // Be8 ----> alpha + alpha 
  G4double alpha_mass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(2,4);
  
  G4double be8_mass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(4,8);
  //std::cout << "Hello from G4Be8FermiFragment, Q = " << (be8_mass - 2.*alpha_mass)/keV <<" keV"<< '\n';
  

  Masses.push_back(alpha_mass);
  Masses.push_back(alpha_mass);
  
  AtomNum.push_back(4);
  AtomNum.push_back(4);
  
  Charges.push_back(2);
  Charges.push_back(2);
 
}
