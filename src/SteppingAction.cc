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
/// \file electromagnetic/TestEm11/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "Run.hh"
#include "HistoManager.hh"

#include "G4SteppingManager.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"

//FADD
#include "G4Neutron.hh"
#include "G4Electron.hh"
#include "G4SystemOfUnits.hh"
#include "DetectorConstruction.hh"
#include "G4Step.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Isotope.hh"

     size_t c11Counter = 0;  


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* event)
:G4UserSteppingAction(),fDetector(det), fEventAction(event)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
//FADD
//    // Get the particle's name
//    G4String particleName = step->GetTrack()->GetDefinition()->GetParticleName();
//
//    // Check if the particle is entering the "Tube" volume
//    G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetPhysicalVolume();
//    if (volume && volume->GetName() == "Tube") {
//        // Retrieve the kinetic energy and position
//        G4double kineticEnergy = step->GetPreStepPoint()->GetKineticEnergy();
//        G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
//
//        // Save the data to a file named after the particle name
//        std::string fileName = particleName + "_entering_particles.csv";
//        std::ofstream outputFile(fileName, std::ios_base::app);
//        if (outputFile.is_open()) {
//            outputFile << kineticEnergy / MeV << "," << position.x() / mm << ","
//                       << position.y() / mm << "," << position.z() / mm << std::endl;
//            outputFile.close();
//        } else {
//            G4cout << "Unable to open the output file." << G4endl;
//        }
//    }


 //G4double edep = step->GetTotalEnergyDeposit();
 //if (edep <= 0.) return;
 //
 //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();     
 //
 ////longitudinal profile of deposited energy
 ////randomize point of energy deposotion
 ////
 //G4StepPoint* prePoint  = step->GetPreStepPoint();
 //G4StepPoint* postPoint = step->GetPostStepPoint(); 
 //G4ThreeVector P1 = prePoint ->GetPosition();
 //G4ThreeVector P2 = postPoint->GetPosition();
 //G4ThreeVector point = P1 + G4UniformRand()*(P2 - P1);
 //if (step->GetTrack()->GetDefinition()->GetPDGCharge() == 0.) point = P2;
 //G4double x = point.x();
 //G4double xshifted = x + 0.5*fDetector->GetAbsorSizeX();
 //analysisManager->FillH1(1, xshifted, edep);
//
 ////"normalized" histogram
 //// 
 //Run* run
 //  = static_cast<Run*>(
 //      G4RunManager::GetRunManager()->GetNonConstCurrentRun());
 //G4int iabs = prePoint->GetTouchableHandle()->GetCopyNumber(1);
 //G4double csdaRange  = run->GetCsdaRange(iabs);
 //if (csdaRange > 0.) { 
 //  G4double density = fDetector->GetAbsorMaterial(iabs)->GetDensity();
 //  G4double xfront  = fDetector->GetXfront(iabs);
 //  G4double xfrontNorm = run->GetXfrontNorm(iabs);
 //  G4double xnorm = xfrontNorm + (x - xfront)/csdaRange;
 //  analysisManager->FillH1(8, xnorm, edep/(csdaRange*density));
 //}
 ////
 ////total energy deposit in absorber
 ////
 //fEventAction->AddEdep(iabs, edep);
 
 //step size of primary particle or charged secondaries
 //
 //FADD
//   G4Track* track = step->GetTrack();
//    G4ParticleDefinition* particleDef = track->GetDefinition();
//    
//    // Check if the particle is leaving the volume
//    if (track->GetNextVolume() == nullptr) {
//        // Extract information about the particle
//        G4String particleName = particleDef->GetParticleName();
//        
//        // Check if the particle is a photon
//        if (particleName == "gamma" || particleName == "opticalphoton") {
//            // This is a photon
//            G4double energy = track->GetKineticEnergy();
//            // Now you can process or store this information as needed
//            G4cout << "Photon Exiting: Energy: " << energy << " MeV" << G4endl;
//        } else {
//            // This is another particle (e.g., electron)
//            G4double energy = track->GetKineticEnergy();
//            // Now you can process or store this information as needed
//            G4cout << "Particle Exiting: " << particleName << ", Energy: " << energy << " MeV" << G4endl;
//        }
//    }
//FADD track neutroni
    if (step->GetTrack()->GetDefinition() == G4Electron::ElectronDefinition()) {
        // Check if new neutrons are created in this step
        const G4TrackVector* secondary = step->GetSecondary();
        for (size_t i = 0; i < secondary->size(); ++i) {
            if ((*secondary)[i]->GetDefinition() == G4Neutron::NeutronDefinition()) {
                G4double neutronEnergy = (*secondary)[i]->GetKineticEnergy() / CLHEP::MeV;
                G4cout << "New Neutron Produced: Energy = " << neutronEnergy << " MeV" << G4endl;
            }
        }
    }

G4Track* track = step->GetTrack();
G4String particleName = track->GetDynamicParticle()->GetDefinition()->GetParticleName();
// Check if the particle is a gamma photon
if (particleName == "gamma")
{
    const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();

    // Check if the process is the photonuclear process
    if (process->GetProcessName() == "photonNuclear")
    {
        // Access the particles produced in the step
        const std::vector<const G4Track*>* secondaryTracks = step->GetSecondaryInCurrentStep();
        if (secondaryTracks)
        {
            for (size_t i = 0; i < secondaryTracks->size(); ++i)
            {
                G4String producedParticleName = (*secondaryTracks)[i]->GetDynamicParticle()->GetDefinition()->GetParticleName();

                // Check if the produced particle is C11
                if (producedParticleName == "C13")
                {
                    // Increment the counter for C11
                    c11Counter++;

                    // Store the position of the C11 particle
                    G4ThreeVector c11Position = (*secondaryTracks)[i]->GetPosition();
                    c11Positions.push_back(c11Position);

                    // You can now process the produced particle information as needed
                    // Print, store, or analyze the information as required
                    G4cout << "Photonuclear Process, Produced Particle: " << producedParticleName << G4endl;
                    G4cout << "Count: " << c11Counter << G4endl;
                }
            }
        }
    }
}

//std::ofstream outputFile("C11_positions.csv"); // Open the file for writing
//
//if (outputFile.is_open()) {
//    // Write the header if the file is just created
//    outputFile << "X,Y,Z\n";
//
//    // Loop through the stored C11 positions and write them to the file
//    for (size_t i = 0; i < c11Positions.size(); ++i) {
//        G4ThreeVector c11Position = c11Positions[i];
//        outputFile << c11Position.x() << "," << c11Position.y() << "," << c11Position.z() << "\n";
//    }
//
//    // Close the file after writing
//    outputFile.close();
//} else {
//    // Handle the case where the file could not be opened
//    G4cout << "Error: Unable to open the file for writing." << G4endl;
//}


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

