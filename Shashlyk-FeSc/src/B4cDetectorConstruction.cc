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
// $Id$
// 
/// \file B4cDetectorConstruction.cc
/// \brief Implementation of the B4cDetectorConstruction class
/*
***************************************MAP FOR TILE NUMBERS************************************

	[00][10][20][30][40][50]
	[01][11][21][31][41][51]
	[02][12][22][32][42][52]
	[03][13][23][33][43][53]
	[04][14][24][34][44][54]
	[05][15][25][35][45][55]

*/
#include "B4cDetectorConstruction.hh"
#include "B4cCalorimeterSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4GenericMessenger.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <stdio.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cDetectorConstruction::B4cDetectorConstruction()
 : G4VUserDetectorConstruction(),
   fMessenger(0),
   fMagField(0),
   fCheckOverlaps(true)
{
  // Define /B4/det commands using generic messenger class
  fMessenger 
    = new G4GenericMessenger(this, "/B4/det/", "Detector construction control");

  // Define /B4/det/setMagField command
  G4GenericMessenger::Command& setMagFieldCmd
    = fMessenger->DeclareMethod("setMagField", 
                                &B4cDetectorConstruction::SetMagField, 
                                "Define magnetic field value (in X direction");
  setMagFieldCmd.SetUnitCategory("Magnetic flux density");                                
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cDetectorConstruction::~B4cDetectorConstruction()
{ 
  delete fMagField;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4cDetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cDetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;
  nistManager->FindOrBuildMaterial("G4_Pb", fromIsotopes);
  nistManager->FindOrBuildMaterial("G4_Fe",fromIsotopes);
  nistManager->FindOrBuildMaterial("G4_POLYSTYRENE",fromIsotopes);
  
  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4cDetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  G4int nofLayers = 66;//66
  G4int nofLayers2 = 36;
  G4double absoThickness = 1.5*mm;
  G4double gapThickness =  4.*mm;
  G4double calorSizeXY  = 600.*mm;
G4double hadcell=10.*cm;
G4double emcell=5.0*cm;
  G4double absoThickness2 = 20. *mm;
  G4double gapThickness2 = 3. *mm;//used to be 6.0 mm
G4double num_of_towers=12;
  G4double layerThickness = absoThickness + gapThickness;
  G4double layerThickness2 = absoThickness2+gapThickness2;
  G4double calorThickness = nofLayers * layerThickness;
  G4double calorThickness2 =nofLayers2*layerThickness2;
  G4double worldSizeXY = 10.0*(5.0* calorSizeXY);
  G4double worldSizeZ  = 5.0 *( 4*calorThickness2+calorThickness); 
  
  // Get materials
  G4Material* defaultMaterial = G4Material::GetMaterial("Galactic");
  G4Material* absorberMaterial = G4Material::GetMaterial("G4_Pb");
  G4Material* gapMaterial = G4Material::GetMaterial("liquidArgon");
  G4Material* gapMaterial2 =G4Material::GetMaterial("G4_POLYSTYRENE");
  G4Material* absorberMaterial2= G4Material::GetMaterial("G4_Fe");
  
  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
    G4cerr << "Cannot retrieve materials already defined. " << G4endl;
    G4cerr << "Exiting application " << G4endl;
    exit(1);
  }  
G4MaterialPropertiesTable* mtpt = new G4MaterialPropertiesTable();
gapMaterial2->SetMaterialPropertiesTable(mtpt);
gapMaterial2->GetIonisation()->SetBirksConstant(.2*mm/MeV);  
  //     
  // World
  //
  G4VSolid* worldS 
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size
                         
  G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                               
  // Calorimeter
  //
G4LogicalVolume* calorLV[6][6];

G4LogicalVolume* calor2LV[6][6]; 
G4LogicalVolume* calor3LV[6][6]; 
char bu[200];
char bu2[200];
char bu3[200];

for(int i=0;i<6;i++){

	for(int j=0;j<6;j++){
sprintf(bu,"Calorimeter%d%d",i,j);
sprintf(bu2,"calor%d%d",i,j);
sprintf(bu3,"Calorimeter2%d%d",i,j);
  G4VSolid* calorimeterS
    = new G4Box("Calorimeter",     // its name
                 hadcell/2, hadcell/2, calorThickness2/2); // its size
           
  calorLV[i][j]
    = new G4LogicalVolume(
                 calorimeterS,     // its solid
                 defaultMaterial,  // its material
                 bu);   // its name

calor2LV[i][j] = new G4LogicalVolume(calorimeterS,defaultMaterial,bu3);
//calor3LV[i][j] = new G4LogicalVolume(calorimeterS,defaultMaterial,bu2);
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(-25.*cm+i*10.*cm,25.*cm-j*10.*cm,calorThickness/2+calorThickness2/2),  // at (0,0,0)
                 calorLV[i][j],          // its logical volume                         
                 bu,    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
	new G4PVPlacement(0,G4ThreeVector(-25.*cm+i*10.*cm,25.*cm-j*10.*cm,3*calorThickness2/2+calorThickness/2),calor2LV[i][j],bu3,worldLV,false,i+j,fCheckOverlaps);
//	new G4PVPlacement(0,G4ThreeVector(-25.*cm+i*10.*cm,25.*cm-j*10.*cm,5*calorThickness2/2+calorThickness/2),calor3LV[i][j],bu2,worldLV,false,i+j,fCheckOverlaps);
	}

}

G4LogicalVolume* calorEM[13][13];
G4int copynono=0;
for(int i=0;i<num_of_towers;i++){
	for(int j=0;j<num_of_towers;j++){
		sprintf(bu,"CalorimeterEM%d%d",i,j);
	//	sprintf(bu2,"calorem%d%d",i,j);

		  G4VSolid* calorimeterEM
    = new G4Box("CalorimeterEM_",     // its name
                 emcell/2, emcell/2, calorThickness/2);
calorEM[i][j] = new G4LogicalVolume(calorimeterEM,defaultMaterial,bu);
//	copynono++;
//G4cout<< -27.5+5*i<<"meow "<<27.5-5*j<<G4endl;
        new G4PVPlacement(0,G4ThreeVector(-27.5*cm+i*5.*cm,27.5*cm-j*5.*cm,0.),calorEM[i][j],bu,worldLV,false,copynono,fCheckOverlaps);
	copynono++;
	}
}

  
  //                                 
  // Layer
  //

G4LogicalVolume* layerLV[6][6];
G4LogicalVolume* layer2LV[6][6];

G4LogicalVolume* layer3LV[6][6];

for(int i=0;i<6;i++){
	
	for(int j=0;j<6;j++){
sprintf(bu,"Layer%d%d",i,j);
sprintf(bu2,"Layer2%d%d",i,j);
sprintf(bu3,"Layer3%d%d",i,j);
  G4VSolid* layerS 
    = new G4Box("Layer",           // its name
                 hadcell/2, hadcell/2, layerThickness2/2); //its size
                         
  layerLV[i][j]
    = new G4LogicalVolume(
                 layerS,           // its solid
                 defaultMaterial,  // its material
                 bu);         // its name
layer2LV[i][j]=new G4LogicalVolume(layerS,defaultMaterial,bu2);
  new G4PVReplica(
                 "Layer",          // its name
                 layerLV[i][j],          // its logical volume
                 calorLV[i][j],          // its mother
                 kZAxis,           // axis of replication
                 nofLayers2,        // number of replica
                 layerThickness2);  // witdth of replica
//layer3LV[i][j]=new G4LogicalVolume(layerS,defaultMaterial,bu3);
//new G4PVReplica("Layer3",layer3LV[i][j],calor3LV[i][j],kZAxis,nofLayers2,layerThickness2);
new G4PVReplica("Layer2",layer2LV[i][j],calor2LV[i][j],kZAxis,nofLayers2,layerThickness2);
	}
}

G4LogicalVolume* layerEMLV[13][13];
for(int i=0;i<num_of_towers;i++){
	for(int j=0;j<num_of_towers;j++){
		sprintf(bu,"LayerEM%d%d",i,j);

		  G4VSolid* layerEM 
    = new G4Box("LayerEM__",           // its name
                 emcell/2, emcell/2, layerThickness/2); //its size
layerEMLV[i][j]=new G4LogicalVolume(layerEM,defaultMaterial,bu);
new G4PVReplica("LayerEM_",layerEMLV[i][j],calorEM[i][j],kZAxis,nofLayers,layerThickness);

	}

}


  //                               
  // Absorber
  //
G4LogicalVolume* absorberLV[6][6];
G4LogicalVolume* absorberLV2[6][6];

G4LogicalVolume* absorberLV3[6][6];

for(int i=0;i<6;i++){
	for(int j=0;j<6;j++){
sprintf(bu,"Abso%d%d",i,j);
sprintf(bu2,"Abso2%d%d",i,j);

sprintf(bu3,"Abso2%d%d",i,j);

  G4VSolid* absorberS 
    = new G4Box("Abso",            // its name
                 hadcell/2, hadcell/2, absoThickness2/2); // its size
                         
  absorberLV[i][j]
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial2, // its material
                 bu);          // its name
                            
absorberLV2[i][j] = new G4LogicalVolume(absorberS,absorberMaterial2,bu2);
//absorberLV3[i][j] = new G4LogicalVolume(absorberS,absorberMaterial2,bu3);

       
   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -gapThickness2/2), // its position
                 absorberLV[i][j],       // its logical volume                         
                 "Abso",           // its name
                 layerLV[i][j],          // its mother  volume
                 false,            // no boolean operation
                 i+j,                // copy number
                 fCheckOverlaps);  // checking overlaps 
new G4PVPlacement(0,G4ThreeVector(0.,0.,-gapThickness2/2),absorberLV2[i][j],"Abso2",layer2LV[i][j],false,i+j,fCheckOverlaps);
//new G4PVPlacement(0,G4ThreeVector(0.,0.,-gapThickness2/2),absorberLV3[i][j],"Abs32",layer3LV[i][j],false,i+j,fCheckOverlaps);

	}
}

G4LogicalVolume* absorberEMLV[13][13];
copynono=0;
for(int i=0;i<num_of_towers;i++){
	for(int j=0;j<num_of_towers;j++){
		sprintf(bu,"AbsoEM%d%d",i,j);
  G4VSolid* absorberEM
    = new G4Box("AbsoEM",            // its name
                 emcell/2, emcell/2, absoThickness/2); // its size
absorberEMLV[i][j] = new G4LogicalVolume(absorberEM,absorberMaterial,bu);
new G4PVPlacement(0,G4ThreeVector(0.,0.,-gapThickness/2),absorberEMLV[i][j],"ABsoem",layerEMLV[i][j],false,copynono,fCheckOverlaps);
copynono++;
	}

}

  //                               
  // Gap
  //
char buffer [200];
G4LogicalVolume* gapLV[6][6];
G4LogicalVolume* gapLV2[6][6];

G4LogicalVolume* gapLV3[6][6];

for(int i=0;i<6;i++){


	for(int j=0;j<6;j++){
	sprintf(buffer,"gap%d%d",i,j);
	sprintf(bu2,"gap2%d%d",i,j);
        sprintf(bu3,"gap3%d%d",i,j);

  G4VSolid* gapS 
    = new G4Box("Gap",             // its name
                 hadcell/2, hadcell/2, gapThickness2/2); // its size

                         
  gapLV[i][j]
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial2,      // its material
                 buffer);           // its name
gapLV2[i][j] = new G4LogicalVolume(gapS,gapMaterial2,bu2);        
//gapLV3[i][j] = new G4LogicalVolume(gapS,gapMaterial2,bu3);
                           
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., absoThickness2/2), // its position
                 gapLV[i][j],            // its logical volume                         
                 buffer,            // its name
                 layerLV[i][j],          // its mother  volume
                 false,            // no boolean operation
                 i+j,                // copy number
                 fCheckOverlaps);  // checking overlaps 
	new G4PVPlacement(0,G4ThreeVector(0.,0., absoThickness2/2),gapLV2[i][j],bu2,layer2LV[i][j],false,i+j,fCheckOverlaps);
// new G4PVPlacement(0,G4ThreeVector(0.,0., absoThickness2/2),gapLV3[i][j],bu3,layer3LV[i][j],false,i+j,fCheckOverlaps);

	}
}

G4LogicalVolume* gapEMLV[13][13];
copynono=0;
for(int i=0;i<num_of_towers;i++){
	for(int j=0;j<num_of_towers;j++){
        sprintf(bu3,"gapEM%d%d",i,j);

  G4VSolid* gapEM
    = new G4Box("GapEM_",             // its name
                 emcell/2, emcell/2, gapThickness/2); // its size
gapEMLV[i][j] = new G4LogicalVolume(gapEM,gapMaterial2,bu3);
        new G4PVPlacement(0,G4ThreeVector(0.,0., absoThickness/2),gapEMLV[i][j],bu3,layerEMLV[i][j],false,copynono,fCheckOverlaps);
copynono++;

	}
}

  //
  // print parameters
  //
 /* G4cout << "\n------------------------------------------------------------"
         << "\n---> The calorimeter is " << nofLayers << " layers of: [ "
         << absoThickness/mm << "mm of " << absorberMaterial->GetName() 
         << " + "
         << gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " 
         << "\n------------------------------------------------------------\n";
  */
  
  // 
  // Sensitive detectors
  //

B4cCalorimeterSD* absoSD[6][6];
B4cCalorimeterSD* gapSD[6][6];
B4cCalorimeterSD* abso2SD[6][6];
B4cCalorimeterSD* gap2SD[6][6];
B4cCalorimeterSD* abso3SD[6][6];
B4cCalorimeterSD* gap3SD[6][6];

char buffer2 [200];
char buffer3 [200];
char buffer4 [200];
char buffer5 [200];
char bu4[200];
char bu5[200];




for(int i=0;i<6;i++){
	for(int j=0;j<6;j++){
	sprintf(buffer2,"AbsorberHitsCollection%d%d",i,j);
	sprintf(buffer3,"AbsorberSD%d%d",i,j);
	sprintf(bu,"Absorber2HitsCollection%d%d",i,j);
	sprintf(bu2,"Absorber2SD%d%d",i,j);
        sprintf(bu4,"Absorber3HitsCollection%d%d",i,j);
        sprintf(bu5,"Absorber3SD%d%d",i,j);
  absoSD[i][j] 
    = new B4cCalorimeterSD(buffer3, buffer2, nofLayers2);
  G4SDManager::GetSDMpointer()->AddNewDetector(absoSD[i][j] );
  absorberLV[i][j]->SetSensitiveDetector(absoSD[i][j]);
abso2SD[i][j] = new B4cCalorimeterSD(bu2,bu, nofLayers2);
G4SDManager::GetSDMpointer()->AddNewDetector(abso2SD[i][j]);
absorberLV2[i][j]->SetSensitiveDetector(abso2SD[i][j]);

//abso3SD[i][j] = new B4cCalorimeterSD(bu5,bu4, nofLayers2);
//G4SDManager::GetSDMpointer()->AddNewDetector(abso3SD[i][j]);
//absorberLV3[i][j]->SetSensitiveDetector(abso3SD[i][j]);


sprintf(buffer4,"GapHitsCollection%d%d",i,j);
sprintf(buffer5,"GapSD%d%d",i,j);
sprintf(bu4,"Gap2HitsCollection%d%d",i,j);
sprintf(bu5,"Gap2SD%d%d",i,j);
sprintf(bu2,"Gap3HitsCollection%d%d",i,j);
sprintf(bu3,"Gap3SD%d%d",i,j);

  gapSD[i][j] 
    = new B4cCalorimeterSD(buffer5, buffer4, nofLayers2);
  G4SDManager::GetSDMpointer()->AddNewDetector(gapSD[i][j] );
  gapLV[i][j]->SetSensitiveDetector(gapSD[i][j]);
gap2SD[i][j] = new B4cCalorimeterSD(bu5,bu4,nofLayers2);
G4SDManager::GetSDMpointer()->AddNewDetector(gap2SD[i][j]);
gapLV2[i][j]->SetSensitiveDetector(gap2SD[i][j]);

//gap3SD[i][j] = new B4cCalorimeterSD(bu3,bu2,nofLayers2);
//G4SDManager::GetSDMpointer()->AddNewDetector(gap3SD[i][j]);
//gapLV3[i][j]->SetSensitiveDetector(gap3SD[i][j]);
	}
}

B4cCalorimeterSD* absoEMSD[13][13];
B4cCalorimeterSD* gapEMSD[13][13];
copynono=0;
for(int i=0;i<num_of_towers;i++){
	for(int j=0;j<num_of_towers;j++){
		        sprintf(buffer2,"AbsorberEMHitsCollection%02d%02d",i,j);
        sprintf(buffer3,"AbsorberEMSD%02d%02d",i,j);
  absoEMSD[i][j]
    = new B4cCalorimeterSD(buffer3, buffer2, nofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(absoEMSD[i][j]);
  absorberEMLV[i][j]->SetSensitiveDetector(absoEMSD[i][j]);
sprintf(buffer4,"GapEMHitsCollection%02d%02d",i,j);
sprintf(buffer5,"GapEMSD%02d%02d",i,j);
copynono++;
  gapEMSD[i][j]
    = new B4cCalorimeterSD(buffer5, buffer4, nofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gapEMSD[i][j]);
  gapEMLV[i][j]->SetSensitiveDetector(gapEMSD[i][j]);

	}

}

  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());//G4VisAttributes::Invisible

  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1,0,1));
  simpleBoxVisAtt->SetVisibility(true);
for(int i=0;i<6;i++){
	for(int j=0;j<6;j++){
	//G4VisAttributes* simpleBoxVisAtt = new G4VisAttributes(G4Colour(i,j,1));
	//simpleBoxVisAtt->SetVisibility(true);
 // calorLV[i][j]->SetVisAttributes(simpleBoxVisAtt);
}
}
  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cDetectorConstruction::SetMagField(G4double fieldValue)
{
  // Apply a global uniform magnetic field along X axis
  G4FieldManager* fieldManager
    = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  // Delete the existing magnetic field
  if ( fMagField )  delete fMagField; 

  if ( fieldValue != 0. ) {
    // create a new one if not null
    fMagField 
      = new G4UniformMagField(G4ThreeVector(fieldValue, 0., 0.));
      
    fieldManager->SetDetectorField(fMagField);
    fieldManager->CreateChordFinder(fMagField);
  } 
  else {
    fMagField = 0;
    fieldManager->SetDetectorField(fMagField);
  }
}
