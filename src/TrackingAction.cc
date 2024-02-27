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
/// \file electromagnetic/TestEm11/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
#include "Run.hh"
#include "HistoManager.hh"

#include "G4Track.hh"
#include "G4RunManager.hh"

//Fadd
#include "G4SystemOfUnits.hh"
#include "TrackingAction.hh"
#include "DetectorConstruction.hh"
#include "Run.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include <fstream>




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* det)
:G4UserTrackingAction(),fDetector(det)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{//Fadd per positroni
// Check if the track corresponds to a positron that has annihilated

if (track->GetDefinition()->GetParticleName() == "e+") {
    if (track->GetTrackStatus() == fStopAndKill) {
        G4VPhysicalVolume* volume = track->GetVolume();
        //G4cout << volume->GetName() << G4endl;

        if (volume && volume->GetName() == "Plexiglass") {
            // Retrieve the position of annihilation
            G4ThreeVector annihilationPos = track->GetPosition();
            annihilationPos = annihilationPos ;
            // Save the position in a CSV file
            std::ofstream outputFile("annihilation_positions.csv", std::ios_base::app);
            if (outputFile.is_open()) {
                // Write the position to the file
                outputFile << annihilationPos.x()/mm << "," << annihilationPos.y()/mm << ","
                           << annihilationPos.z()/mm << std::endl;
                outputFile.close();
            } else {
                G4cout << "Unable to open the output file." << G4endl;
            }
        }
    }
}








 //G4int trackID = track->GetTrackID();
 //
 //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 //Run* run 
 //  = static_cast<Run*>(
 //      G4RunManager::GetRunManager()->GetNonConstCurrentRun());
 //       
 ////track length of primary particle or charged secondaries
 ////
 //G4double tracklen = track->GetTrackLength();
 //if (trackID == 1) {
 //   run->AddTrackLength(tracklen);
 //   analysisManager->FillH1(3, tracklen);
 //} else if (track->GetDefinition()->GetPDGCharge() != 0.)
 //   analysisManager->FillH1(6, tracklen);
 //          
 ////extract projected range of primary particle
 ////
 //if (trackID == 1) {
 //  G4double x = track->GetPosition().x() + 0.5*fDetector->GetAbsorSizeX();
 //  run->AddProjRange(x);
 //  analysisManager->FillH1(5, x);
 //}
 //           
 ////mean step size of primary particle
 ////
 //if (trackID == 1) {
 //  G4int nbOfSteps = track->GetCurrentStepNumber();
 //  G4double stepSize = tracklen/nbOfSteps;
 //  run->AddStepSize(nbOfSteps,stepSize);
 //}
 //           
 ////status of primary particle : absorbed, transmited, reflected ?
 ////
 //if (trackID == 1) {
 // G4int status = 0;
 // if (!track->GetNextVolume()) {
 //   if (track->GetMomentumDirection().x() > 0.) status = 1;
 //   else                                        status = 2;
 // }
 // run->AddTrackStatus(status);
 //}   





 
//FADD
//   // Declare variables outside of the if-block
//    G4String particleName = track->GetDefinition()->GetParticleName();
//    G4double energy = track->GetKineticEnergy();
//
//    // Check if the particle is leaving the volume
//    //NOTA BENE: A QUANTO PARE TIENE CONTO DEL WORLD VOLUME!!!!
//    if (!track->GetNextVolume()) {
//        if (particleName == "gamma") {
//            // This is a gamma particle
//            // Retrieve the position of exit
//            G4ThreeVector exitPos = track->GetPosition();
//
//            // Retrieve the momentum direction of the particle
//            G4ThreeVector momentumDir = track->GetMomentumDirection();
//
//            // Retrieve the volume at the point of exit
//            G4LogicalVolume* exitVolume = track->GetVolume()->GetLogicalVolume();
//
//            // Retrieve the normal vector of the surface at the point of exit
//            G4ThreeVector normalVector = exitVolume->GetSolid()->SurfaceNormal(exitPos);
//
//
//            // Define the reference axis (e.g., the x-axis)
//            G4ThreeVector xAxis(1.0, 0.0, 0.0);
//            
//            // Calculate the angle between the momentum direction and the reference axis
//            G4double angle = momentumDir.angle(normalVector) * 180.0 / CLHEP::pi; // Convert to degrees
//            
//            // Determine if the gamma is moving forward or backward
//            G4bool isForward = (momentumDir.dot(xAxis) > 0.0);
//            
//            // Assign a minus sign to the angle if isForward is false
//            if (!isForward) {
//                angle = -angle;
//            }
//            // Save the energy, exit position, and angle in the same CSV file for gamma particles
//            std::ofstream gammaDataFile("gamma_data.csv", std::ios_base::app);
//            if (gammaDataFile.is_open()) {
//                gammaDataFile << energy / MeV << "," << exitPos.x() / mm << ","
//                              << exitPos.y() / mm << "," << exitPos.z() / mm << ","
//                              << angle << std::endl;
//                gammaDataFile.close();
//            } else {
//                G4cout << "Unable to open the gamma data output file." << G4endl;
//            }
//        } else if (particleName == "e+") {
//            // This is a positron
//            // Save positron energy to a CSV file
//            std::ofstream positronFile("positrons.csv", std::ios_base::app);
//            if (positronFile.is_open()) {
//                positronFile << energy / MeV << std::endl;
//                positronFile.close();
//            }
//        }else if (particleName == "e+") {
//            // This is a positron
//            // Save positron energy to a CSV file
//            std::ofstream positronFile("positrons.csv", std::ios_base::app);
//            if (positronFile.is_open()) {
//                positronFile << energy / MeV << std::endl;
//                positronFile.close();
//            }
//        } else if (particleName == "e-") {
//            // This is an electron
//            // Save electron energy to a CSV file
//            std::ofstream electronFile("electrons.csv", std::ios_base::app);
//            if (electronFile.is_open()) {
//                electronFile << energy / MeV << std::endl;
//                electronFile.close();
//            }
//        } else if (particleName == "neutron") {
//            // This is a neutron
//            // Save neutron energy to a CSV file
//            std::ofstream neutronFile("neutrons.csv", std::ios_base::app);
//            if (neutronFile.is_open()) {
//                neutronFile << energy / MeV << std::endl;
//                neutronFile.close();
//            }
//        }  else if (particleName == "GenericIon") {
//            // This is a neutron
//            // Save neutron energy to a CSV file
//            std::ofstream GenericIonFile("GenericIon.csv", std::ios_base::app);
//            if (GenericIonFile.is_open()) {
//                GenericIonFile << energy / MeV << std::endl;
//                GenericIonFile.close();
//            }
//        }
//    }


 G4String particleType = track->GetDefinition()->GetParticleName();

   // // Check if the track is entering the "Tube" volume
   // if (track->GetVolume()->GetName() == "Tube") {
   //     G4double energy = track->GetKineticEnergy();
   //     G4ThreeVector position = track->GetPosition();
//
   //     // Generate a separate CSV file for each particle type
   //     G4String fileName = particleType + "_data_tube.csv";
//
   //     // Save the particle type, energy, and position to the corresponding CSV file
   //     std::ofstream outputFile(fileName, std::ios_base::app);
   //     if (outputFile.is_open()) {
   //         outputFile << energy / MeV << ","
   //                    << position.x() / mm << "," << position.y() / mm << ","
   //                    << position.z() / mm << std::endl;
   //         outputFile.close();
   //     } else {
   //         G4cout << "Unable to open the output file." << G4endl;
   //     }
   // }
//
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

