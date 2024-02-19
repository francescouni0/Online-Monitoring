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
/// \file electromagnetic/TestEm11/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4RunManager.hh"
#include <iomanip>
//Fadd
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSDoseDeposit.hh"
#include "G4Tubs.hh" // Include G4Tubs for a tube
#include "G4Sphere.hh" // Include G4Sphere for a sphere
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "Detector.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4VisAttributes.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),fDefaultMaterial(0),fPhysiWorld(0),
 fDetectorMessenger(0)
{
  // default parameter values of the absorbers
  fNbOfAbsor = 1;
  fAbsorThickness[0] = 0*mm;        //dummy, for initialization   
  fAbsorThickness[1] = 1*mm;  
  fAbsorSizeYZ       = 1.*mm;
  for (G4int iAbs=0; iAbs<kMaxAbsor; iAbs++) {
    fNbOfDivisions[iAbs]  = 1;
  }  
  ComputeParameters();

  // materials
  DefineMaterials();
  SetAbsorMaterial(1,"G4_PLEXIGLASS");

  // create commands for interactive definition of the calorimeter
  fDetectorMessenger = new DetectorMessenger(this);






}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{ 
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  G4NistManager* man = G4NistManager::Instance();
  
  man->FindOrBuildMaterial("G4_Al");
  man->FindOrBuildMaterial("G4_Si");
  man->FindOrBuildMaterial("G4_Fe");
  man->FindOrBuildMaterial("G4_Cu");  
  man->FindOrBuildMaterial("G4_Ge");
  man->FindOrBuildMaterial("G4_Mo");
  man->FindOrBuildMaterial("G4_Ta");
  man->FindOrBuildMaterial("G4_W");
  man->FindOrBuildMaterial("G4_Au");
  man->FindOrBuildMaterial("G4_Pb");  
  man->FindOrBuildMaterial("G4_PbWO4");
  man->FindOrBuildMaterial("G4_SODIUM_IODIDE");
  man->FindOrBuildMaterial("G4_PLEXIGLASS"); 
  man->FindOrBuildMaterial("G4_AIR");
  man->FindOrBuildMaterial("G4_WATER");
  
  G4Element* H = man->FindOrBuildElement("H"); 
  G4Element* O = man->FindOrBuildElement("O");
  G4Element* C = man->FindOrBuildElement("C"); //Fadd
  G4Element* Na = man->FindOrBuildElement("Na");
  G4Element* I = man->FindOrBuildElement("I");

  
  G4Material* H2O = 
  new G4Material("Water", 1.000*g/cm3, 2);
  H2O->AddElement(H, 2);
  H2O->AddElement(O, 1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);

  G4Material* plexiglass = //Fadd
  new G4Material("Plexiglass", 1.19 * g/cm3, 3);
  plexiglass->AddElement(C, 5);
  plexiglass->AddElement(H, 8);
  plexiglass->AddElement(O, 2);
  plexiglass->GetIonisation()->SetMeanExcitationEnergy(74.0*eV);

  G4MaterialPropertiesTable* plexiglasMPT = new G4MaterialPropertiesTable();

  // Set refractive index for Plexiglas (example values)
  G4double plexiglasRIndex[] = { 1.49, 1.49, 1.49 };
  G4double photonEnergy[] = { 2.0 * eV, 3.0 * eV, 4.0 * eV };
  plexiglasMPT->AddProperty("RINDEX", photonEnergy, plexiglasRIndex, 3);
  
  // Assign the material properties table to Plexiglas
  plexiglass->SetMaterialPropertiesTable(plexiglasMPT);

////

  
  // Assign the material properties table to Plexiglas

  G4Material* NaI = //Fadd
  new G4Material("NaI", 3.67 * g/cm3, 2);
  NaI->AddElement(Na, 1);
  NaI->AddElement(I, 1);  
  G4MaterialPropertiesTable* NaI_MPT = new G4MaterialPropertiesTable();
  G4double energy[2]={1.239841939*eV/0.9,1.239841939*eV/0.2};

  G4double rindexNaI[2]={1.78,1.78};
  G4double fraction[2]={1.0,1.0};

  //NaI_MPT->AddProperty("RINDEX", energy, rindexNaI,2);
  //NaI_MPT->AddProperty("SCINTILLATIONCOMPONENT1", energy,fraction,2);
  //NaI_MPT->AddConstProperty("SCINTILLATIONYIELD", 38./keV);
  //NaI_MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 250*ns);
  //NaI_MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  //G4double NaI_absLength[] = {380.*cm, 380.*cm}; //sbagliata
  //NaI_MPT->AddProperty("ABSLENGTH", energy, NaI_absLength, 2);
  //NaI->SetMaterialPropertiesTable(NaI_MPT);

  G4double density     = universe_mean_density;    //from PhysicalConstants.h
  G4double pressure    = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  G4Material* Galactic =   
  new G4Material("Galactic", 1., 1.008*g/mole, density,
                             kStateGas,temperature,pressure);
                             
  fDefaultMaterial = Galactic;
  
  //  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ComputeParameters()
{
  // Compute total thickness of absorbers
  fAbsorSizeX = 0.;
  for (G4int iAbs=1; iAbs<=fNbOfAbsor; iAbs++) {
    fAbsorSizeX += fAbsorThickness[iAbs];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // complete the Calor parameters definition
  ComputeParameters();

  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //
  // World
  //
G4double worldSizeX = fAbsorSizeX*50;  // Make the world volume twice as big as the absorber
G4double worldSizeYZ = fAbsorSizeYZ*20; 

G4Box* solidWorld =
  new G4Box("World",
             worldSizeX/2, worldSizeYZ/2, worldSizeYZ/2);     //size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,              //solid
                        fDefaultMaterial,        //material
                        "World");                //name

  fPhysiWorld = 
    new G4PVPlacement(0,                        //no rotation
                        G4ThreeVector(),        //at (0,0,0)
                      logicWorld,               //logical volume
                      "World",                  //name
                       0,                       //mother volume
                       false,                   //no boolean operation
                       0);                      //copy number

  //FADD Compton
//  G4NistManager* nist = G4NistManager::Instance();
//  G4Material* plastic=nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"); 
//
//  
//  ///////////////////7
//  G4Rotate3D rotZscatt(180*deg, G4ThreeVector(0,1,0));
//  
//  G4Translate3D transXscatt(G4ThreeVector(0*cm ,0*cm,-20*cm));  //prende la metà della dimensione -4
//  G4Transform3D transforscatt=(rotZscatt)*(transXscatt);
//
//  // Create the tube (hollow cylindrical shape) with inner and outer radii and height
//  G4Box* solidScatterer= new G4Box("Scatterer", 15*cm, 10*cm, 1*cm);
//  G4LogicalVolume* logicScatt =
//      new G4LogicalVolume(solidScatterer,
//                          plastic, // Air material
//                          "Scatterer");
//  new G4PVPlacement(transforscatt,
//                    logicScatt,
//                    "Scatterer",
//                    logicWorld,
//                    false,
//                    0);
////
//piombino
//  G4NistManager* nist = G4NistManager::Instance();
//  G4Material* lead = nist->FindOrBuildMaterial("G4_Pb");
//  G4double radius = 0.13* cm;
//  G4Sphere* solidSphere = new G4Sphere("Sphere", 0, radius, 0, 2 * M_PI, 0, M_PI);
//  G4LogicalVolume* logicSphere =
//      new G4LogicalVolume(solidSphere,
//                          lead, // Air material
//                          "Sphere");
//   // Set visual attributes for visualization
//  G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(0, 0, 1)); // blue color
//  VisAtt->SetVisibility(true);
//  logicSphere->SetVisAttributes(VisAtt);
//  
//  G4ThreeVector spherePosition = G4ThreeVector(-13.3*mm, 0, 0); // Position of the sphere inside the world
//  new G4PVPlacement(0,
//                    spherePosition,
//                    logicSphere,
//                    "Sphere",
//                    logicWorld,
//                    false,
//                    0);
////piombino2
//
//  G4Sphere* solidSphere1 = new G4Sphere("Sphere1", 0, radius, 0, 2 * M_PI, 0, M_PI);
//  G4LogicalVolume* logicSphere1 =
//      new G4LogicalVolume(solidSphere1,
//                          lead, // Air material
//                          "Sphere1");
//   // Set visual attributes for visualization
//  VisAtt->SetVisibility(true);
//  logicSphere->SetVisAttributes(VisAtt);
//  
//  G4ThreeVector spherePosition1 = G4ThreeVector(-15.3*mm, 0, 0); // Position of the sphere inside the world
//  new G4PVPlacement(0,
//                    spherePosition1,
//                    logicSphere1,
//                    "Sphere1",
//                    logicWorld,
//                    false,
//                    0);
//  G4Box* solidBox = new G4Box("Box",1*cm, 1*cm, 0.5*cm);
//  G4NistManager* nist = G4NistManager::Instance();
//  G4Material* lead = nist->FindOrBuildMaterial("G4_Pb");
//  G4LogicalVolume* logicBox =
//      new G4LogicalVolume(solidBox,
//                          lead, // Air material
//                          "Box");
//  G4ThreeVector boxposition= G4ThreeVector(-20*mm, 0, 0);
//  new G4PVPlacement(rotation,
//                    boxposition,
//                    logicBox,
//                    "Box",
//                    logicWorld,
//                    false,
//                    0);
   // Set visual attributes for visualization

//Cono
//  // Set cone dimensions
//  G4double innerRadius1 = 0.0 * mm; // Inner radius of the cone at the small end
//  G4double outerRadius1 = 2.5 * mm; // Outer radius of the cone at the small end
//  G4double innerRadius2 = 0.0 * mm; // Inner radius of the cone at the large end
//  G4double outerRadius2 = 5.0 * mm; // Outer radius of the cone at the large end
//  G4double height = 10 * mm; // Height of the cone
//  // Create solid cone
//  G4Cons* solidCone = new G4Cons("solidCone", innerRadius1, outerRadius1, innerRadius2, outerRadius2, height/2.0 , 0.0, 360.0 * deg);
//  // Create logical volume for the cone
//  G4LogicalVolume* logicalCone = new G4LogicalVolume(solidCone, tubemat, "logicalCone");
//  // Set visual attributes for visualization
//  G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(0, 0, 1)); // blue color
//  VisAtt->SetVisibility(true);
//  logicalCone->SetVisAttributes(VisAtt);
//
//  G4ThreeVector conePosition = G4ThreeVector(-30.6*cm, 0, 0); // Position of the sphere inside the world
//  G4RotationMatrix* rotation = new G4RotationMatrix();
//  rotation->rotateY(-90.0 * deg);
//  new G4PVPlacement(rotation,
//                    conePosition,
//                    logicalCone,
//                    "Cone",
//                    logicWorld,
//                    false,
//                    0);

  //CONE 2
//    // Set cone dimensions
//  G4double innerRadius3 = 0.0 * mm; // Inner radius of the cone at the small end
//  G4double outerRadius3 = 5.0 * mm; // Outer radius of the cone at the small end
//  G4double innerRadius4 = 0.0 * mm; // Inner radius of the cone at the large end
//  G4double outerRadius4 = 5.0 * mm; // Outer radius of the cone at the large end
//  G4double height2 = 0.5 * mm; // Height of the cone
//
//  // Create solid cone
//  G4Cons* solidCone2 = new G4Cons("solidCone2", innerRadius3, outerRadius3, innerRadius4, outerRadius4, height2/2.0 , 0.0, 360.0 * deg);
//  // Create logical volume for the cone
//  G4LogicalVolume* logicalCone2 = new G4LogicalVolume(solidCone2, tubemat, "logicalCone2");
//  // Set visual attributes for visualization
//  G4VisAttributes* VisAtt2 = new G4VisAttributes(G4Colour(0.2, 1, 1)); // blue color
//  VisAtt2->SetVisibility(true);
//  logicalCone2->SetVisAttributes(VisAtt2);
//
//  G4ThreeVector conePosition2 = G4ThreeVector(-30*cm, 0, 0); // Position of the sphere inside the world
//  new G4PVPlacement(rotation,
//                    conePosition2,
//                    logicalCone2,
//                    "Cone2",
//                    logicWorld,
//                    false,
//                    0);



  //COLLIMATORE buchi
  //  // Define collimator dimensions
  //  G4double collimatorX = 10.0 * cm;
  //  G4double collimatorY = 12.0 * cm;
  //  G4double collimatorZ = 1.0 * cm;
  //  G4double separationBetweenCollimators = 40.0 * cm; // Adjust separation between collimators as needed
  //  G4double holeRadius = 6 * mm; // Adjust hole radius as needed
  //  G4double holeSeparation = 1.5*cm; // Adjust hole separation as needed
//
  //  // Get collimator material (e.g., lead)
  //  G4NistManager* nist = G4NistManager::Instance();
  //  G4Material* lead = nist->FindOrBuildMaterial("G4_Pb");
//
  //  // Create collimator solid shape (box in this case)
  //  G4VSolid* collimatorSolid = new G4Box("Collimator_Solid", collimatorX / 2, collimatorY / 2, collimatorZ / 2);
//
  //  // Create a single hole solid (cylindrical hole)
  //  G4Tubs* holeSolid = new G4Tubs("Hole_Solid", 0.0, holeRadius, collimatorZ / 2 + 1.0 * cm, 0.0, 360.0 * degree);
//
  //  // Calculate the number of holes that can fit in the collimator
  //  int numHolesX = static_cast<int>((collimatorX - 2 * holeSeparation) / holeSeparation);
  //  int numHolesY = static_cast<int>((collimatorY - 2 * holeSeparation) / holeSeparation);
//
  //  // Create the holes pattern by subtracting hole solids from the collimator solid
  //  G4VSolid* collimatorWithHolesSolid = collimatorSolid;
  //  for (int i = 0; i < numHolesX; ++i) {
  //      for (int j = 0; j < numHolesY; ++j) {
  //          G4ThreeVector holePosition((-collimatorX / 2) + holeSeparation + i * holeSeparation, (-collimatorY / 2) + holeSeparation + j * holeSeparation, 0.0);
  //          G4SubtractionSolid* collimatorMinusHole = new G4SubtractionSolid("Collimator_Minus_Hole", collimatorWithHolesSolid, holeSolid, 0, holePosition);
  //          collimatorWithHolesSolid = collimatorMinusHole;
  //      }
  //  }
//
  //  // Create logical volume for the first collimator
  //  G4LogicalVolume* firstCollimatorLogical = new G4LogicalVolume(collimatorWithHolesSolid, lead, "FirstCollimator_Logical");
//
  //  // Position the first collimator
  //  G4ThreeVector firstCollimatorPosition = G4ThreeVector(0.0, 0.0, -20.0 * cm); // Example position 20 cm in front of the detector
  //  G4RotationMatrix* noRotation = new G4RotationMatrix();
  //  G4Transform3D firstCollimatorTransform = G4Transform3D(*noRotation, firstCollimatorPosition);
//
  //  // Place the first collimator in the detector setup
  //  new G4PVPlacement(firstCollimatorTransform, firstCollimatorLogical, "FirstCollimator_Physical", logicWorld, false, 0);
//
  //  // Create a second collimator with holes
  //  G4VSolid* secondCollimatorWithHolesSolid = collimatorWithHolesSolid;
//
  //  // Create logical volume for the second collimator
  //  G4LogicalVolume* secondCollimatorLogical = new G4LogicalVolume(secondCollimatorWithHolesSolid, lead, "SecondCollimator_Logical");
//
  //  // Position the second collimator
  //  G4ThreeVector secondCollimatorPosition = G4ThreeVector(0.0, 0.0, -20*cm+separationBetweenCollimators); // Adjust position for separation
  //  G4Transform3D secondCollimatorTransform = G4Transform3D(*noRotation, secondCollimatorPosition);
//
  //  // Place the second collimator in the detector setup
  //  new G4PVPlacement(secondCollimatorTransform, secondCollimatorLogical, "SecondCollimator_Physical", logicWorld, false, 1);
//
  //  // Create red visualization attribute
  //  G4VisAttributes* redVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // Red color
//
  //  // Set visualization attributes for the first collimator
  //  firstCollimatorLogical->SetVisAttributes(redVisAtt);
//
  //  // Set visualization attributes for the second collimator
  //  secondCollimatorLogical->SetVisAttributes(redVisAtt);


//Septa1 & 2

  // Define septa dimensions
//  G4double septaX = 3 * cm;
//  G4double septaY = 67.5 * cm;
//  G4double septaZ = 20 * cm;
//
//  G4Material* lead = nist->FindOrBuildMaterial("G4_Pb");
//  G4VSolid* collimatorSolid = new G4Box("Collimator_Solid", septaX / 2, septaY / 2, septaZ / 2);
//
//  //  // Create logical volume for the first collimator
//
//  G4LogicalVolume* firstCollimatorLogical = new G4LogicalVolume(collimatorSolid, lead, "FirstCollimator_Logical");
//  // Position the first collimator
//  G4ThreeVector firstCollimatorPosition = G4ThreeVector(10 * cm, 0.0, 44*cm); // Example position 20 cm in front of the detector
//  G4RotationMatrix* Rotation = new G4RotationMatrix();
//  Rotation->rotateX(90.0 * deg);
//  G4Transform3D firstCollimatorTransform = G4Transform3D(*Rotation, firstCollimatorPosition);
//  // Place the first collimator in the detector setup
//  new G4PVPlacement(firstCollimatorTransform, firstCollimatorLogical, "FirstCollimator_Physical", logicWorld, false, 0);
//  // Create a second collimator with holes
//  G4VSolid* secondCollimatorSolid = collimatorSolid;
//  // Create logical volume for the second collimator
//  G4LogicalVolume* secondCollimatorLogical = new G4LogicalVolume(secondCollimatorSolid, lead, "SecondCollimator_Logical");
//  // Position the second collimator
//  G4ThreeVector secondCollimatorPosition = G4ThreeVector(10*cm, 0.0, -44*cm); // Adjust position for separation
//  G4Transform3D secondCollimatorTransform = G4Transform3D(*Rotation, secondCollimatorPosition);
//  // Place the second collimator in the detector setup
//  new G4PVPlacement(secondCollimatorTransform, secondCollimatorLogical, "SecondCollimator_Physical", logicWorld, false, 1);
//  // Create red visualization attribute
//  G4VisAttributes* redVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // Red color
//  // Set visualization attributes for the first collimator
//  firstCollimatorLogical->SetVisAttributes(redVisAtt);
//  // Set visualization attributes for the second collimator
//  secondCollimatorLogical->SetVisAttributes(redVisAtt);
//
//
//  //Septa 3&4
//  // Define septa dimensions
//
//
//  G4VSolid* thirdcollimatorSolid = new G4Box("ThirdCollimator_Solid", septaX / 2, septaY / 2, septaZ / 2);
//
//  //  // Create logical volume for the first collimator
//
//  G4LogicalVolume* thirdCollimatorLogical = new G4LogicalVolume(thirdcollimatorSolid, lead, "ThirdCollimator_Logical");
//  // Position the first collimator
//  G4ThreeVector thirdCollimatorPosition = G4ThreeVector(16 * cm, 0.0, 44*cm); // Example position 20 cm in front of the detector
//  G4Transform3D thirdCollimatorTransform = G4Transform3D(*Rotation, thirdCollimatorPosition);
//  // Place the first collimator in the detector setup
//  new G4PVPlacement(thirdCollimatorTransform, thirdCollimatorLogical, "ThirdCollimator_Physical", logicWorld, false, 0);
//  // Create a second collimator with holes
//  G4VSolid* fourthCollimatorSolid = collimatorSolid;
//  // Create logical volume for the fourth collimator
//  G4LogicalVolume* fourthCollimatorLogical = new G4LogicalVolume(fourthCollimatorSolid, lead, "fourthCollimator_Logical");
//  // Position the fourth collimator
//  G4ThreeVector fourthCollimatorPosition = G4ThreeVector(16*cm, 0.0, -44*cm); // Adjust position for separation
//  G4Transform3D fourthCollimatorTransform = G4Transform3D(*Rotation, fourthCollimatorPosition);
//  // Place the fourth collimator in the detector setup
//  new G4PVPlacement(fourthCollimatorTransform, fourthCollimatorLogical, "fourthCollimator_Physical", logicWorld, false, 1);
//  // Create red visualization attribute
//  G4VisAttributes* cVisAtt = new G4VisAttributes(G4Colour(0.5, 0.0, 0.25)); // Red color
//  // Set visualization attributes for the first collimator
//  thirdCollimatorLogical->SetVisAttributes(cVisAtt);
//  // Set visualization attributes for the fourth collimator
//  fourthCollimatorLogical->SetVisAttributes(cVisAtt);

  



  //
  // Absorbers
  //
  fXfront[0] = -0.5*fAbsorSizeX;
  //
  for (G4int k=1; k<=fNbOfAbsor; k++) {
    G4Material* material = fAbsorMaterial[k];
    G4String matname = material->GetName();
      
    G4Box* solidAbsor =
      new G4Box(matname,fAbsorThickness[k]/2,fAbsorSizeYZ/2,fAbsorSizeYZ/2);

    G4LogicalVolume* logicAbsor =
      new G4LogicalVolume(solidAbsor,           // solid
                          material,             // material
                          matname);             // name
                                     
    fXfront[k] = fXfront[k-1] + fAbsorThickness[k-1];    
    G4double xcenter = fXfront[k]+0.5*fAbsorThickness[k];
    G4ThreeVector position = G4ThreeVector(0.,0.,0.);
  
      new G4PVPlacement(0,                     //no rotation
                              G4ThreeVector(0.,0.,0.),        //position
                        logicAbsor,            //logical volume        
                        matname,               //name
                        logicWorld,            //mother
                        false,                 //no boulean operat
                        k);                    //copy number
    
    // divisions, if any
    //
    G4double LayerThickness = fAbsorThickness[k]/fNbOfDivisions[k];
    G4Box* solidLayer =   
      new G4Box(matname,LayerThickness/2,fAbsorSizeYZ/2,fAbsorSizeYZ/2);
                       
    G4LogicalVolume* logicLayer =
      new G4LogicalVolume(solidLayer,      //solid
                          material,        //material
                          matname);        //name
                                       
      new G4PVReplica(matname,             //name
                            logicLayer,    //logical volume
                            logicAbsor,    //mother
                      kXAxis,              //axis of replication
                      fNbOfDivisions[k],   //number of replica
                      LayerThickness);     //witdth of replica    
  }
  //FADD

  G4NistManager* man = G4NistManager::Instance();

  


  //solidScintillator= new G4Box("Scintillator", 5*mm, 5*mm, 15*cm);

  logicScintillator= new G4LogicalVolume(solidScintillator, man->FindOrBuildMaterial("G4_Pb"), "logicalScintillator");

  fscoringVolume = logicScintillator;

  solidDetector= new G4Box("solidDetector", 15*cm, 10*cm, 1*mm); //prende il doppio della dimensione

  logicDetector= new G4LogicalVolume(solidDetector, man->FindOrBuildMaterial("G4_SODIUM_IODIDE"),"logicDetector");

  for(G4int i=0;i<1;i++)
  {
    for(G4int j=0;j<2;j++)
    {
      G4Rotate3D rotZ(j*180*deg, G4ThreeVector(0,1,0));
      //fondo
      //G4Rotate3D rotZ(90*deg, G4ThreeVector(0,1,0));

      //j*180*deg
      //j*22.5*deg rotaaizone
      //G4Translate3D transXscint(G4ThreeVector(-40*cm + i*20*mm,0.*mm,5./tan(1.8/2*deg)*mm + 5.*mm ));

      G4Translate3D transXdet(G4ThreeVector(0*cm ,0*cm,30*cm));  //prende la metà della dimensione -45

      //-40*cm + i*10*mm Direzione lungo asse x
      //5./tan(22.5/2*deg)*cm + 5.*cm  lungo z

      //G4Transform3D transformScint=(rotZ)*(transXscint);
      G4Transform3D transformdet=(rotZ)*(transXdet);


      //physiScintillator = new G4PVPlacement(transformScint, logicScintillator, "physScintillator", logicWorld, false, 0);

      physiDetector = new G4PVPlacement(transformdet, logicDetector, "physDetector", logicWorld, false, j+i*200);



    }
  }

  PrintParameters();

  //always return the physical World
  //
  return fPhysiWorld;
}

//FADD PET

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//void DetectorConstruction::ConstructScintillator(){
//  G4NistManager* man = G4NistManager::Instance();
//  man->FindOrBuildMaterial("G4_PbWO4");
//  
//
//  solidScintillator= new G4Box("Scintillator", 5*mm, 5*mm, 6*mm);
//
//  logicScintillator= new G4LogicalVolume(solidScintillator, man->FindOrBuildMaterial("G4_SODIUM_IODIDE"), "logicalScintillator");
//
//  fscoringVolume = logicScintillator;
//
//  for(G4int i=0;i<6;i++)
//  {
//    for(G4int j=0;j<16;j++)
//    {
//      G4Rotate3D rotZ(j*22.5*deg, G4ThreeVector(1,0,0));
//      G4Translate3D transXscint(G4ThreeVector(-60*cm + i*15*mm,0.*mm,5./tan(22.5/2*deg)*mm + 5.*mm +50*cm));
//
//      G4Transform3D transformScint=(rotZ)*(transXscint);
//
//      physiScintillator = new G4PVPlacement(transformScint, logicScintillator, "physScintillator", logicWorld, false, 0,true);
//
//
//    }
//  }
//
//
//}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void DetectorConstruction::PrintParameters()
{
  G4cout << "\n-------------------------------------------------------------"
         << "\n ---> The Absorber is " << fNbOfAbsor << " layers of:";
  for (G4int i=1; i<=fNbOfAbsor; i++)
     {
      G4cout << "\n \t" << std::setw(12) << fAbsorMaterial[i]->GetName() <<": "
              << std::setw(6) << G4BestUnit(fAbsorThickness[i],"Length")
              << "  divided in " << fNbOfDivisions[i] << " slices";
     }
  G4cout << "\n-------------------------------------------------------------\n"
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfAbsor(G4int ival)
{
  // set the number of Absorbers
  //
  if (ival < 1 || ival > (kMaxAbsor-1))
    { G4cout << "\n ---> warning from SetfNbOfAbsor: "
             << ival << " must be at least 1 and and most " << kMaxAbsor-1
             << ". Command refused" << G4endl;
      return;
    }
  fNbOfAbsor = ival;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorMaterial(G4int iabs,const G4String& material)
{
  // search the material by its name
  //
  if (iabs > fNbOfAbsor || iabs <= 0)
    { G4cout << "\n --->warning from SetfAbsorMaterial: absor number "
             << iabs << " out of range. Command refused" << G4endl;
      return;
    }

  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(material);
  if (pttoMaterial) {
      fAbsorMaterial[iabs] = pttoMaterial;
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
      G4cout << "\n " << pttoMaterial << G4endl; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorThickness(G4int iabs,G4double val)
{
  // change Absorber thickness
  //
  if (iabs > fNbOfAbsor || iabs <= 0)
    { G4cout << "\n --->warning from SetfAbsorThickness: absor number "
             << iabs << " out of range. Command refused" << G4endl;
      return;
    }
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetfAbsorThickness: thickness "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  fAbsorThickness[iabs] = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorSizeYZ(G4double val)
{
  // change the transverse size
  //
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetfAbsorSizeYZ: thickness "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  fAbsorSizeYZ = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfDivisions(G4int iabs, G4int ival)
{
  // set the number of divisions
  //
  if (iabs > fNbOfAbsor || iabs < 1)
    { G4cout << "\n --->warning from SetNbOfDivisions: absor number "
             << iabs << " out of range. Command refused" << G4endl;
      return;
    }
      
  if (ival < 1)
    { G4cout << "\n --->warning from SetNbOfDivisions: "
             << ival << " must be at least 1. Command refused" << G4endl;
      return;
    }
  fNbOfDivisions[iabs] = ival;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
    if ( fFieldMessenger.Get() == 0 ) {
        // Create global magnetic field messenger.
        // Uniform magnetic field is then created automatically if
        // the field value is not zero.
        G4ThreeVector fieldValue = G4ThreeVector();
        G4GlobalMagFieldMessenger* msg =
        new G4GlobalMagFieldMessenger(fieldValue);
        //msg->SetVerboseLevel(1);
        G4AutoDelete::Register(msg);
        fFieldMessenger.Put( msg );
    }

    SensitiveDetector* det =  new SensitiveDetector("SensDet");

    logicDetector->SetSensitiveDetector(det);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

