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
/// \file electromagnetic/TestEm11/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "Run.hh"
//#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"

#include "G4TrackingManager.hh"
#include "G4VProcess.hh"
#include <fstream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<G4ThreeVector>& EventAction::GetC11Positions()
{
    return c11Positions;
}

EventAction::EventAction(DetectorConstruction* det)
:G4UserEventAction(), fDetector(det)
{
    

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  //energy deposited per event
  for (G4int k=0; k<kMaxAbsor; k++) { fEdepAbsor[k] = 0.0; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event )
{
  //get Run
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  
  //plot energy deposited per event
  //
  G4double TotalEdep(0.);
  for (G4int k=1; k<=fDetector->GetNbOfAbsor(); k++) {
    if (fEdepAbsor[k] > 0.) {
      run->AddEdep(k,fEdepAbsor[k]);
      //G4AnalysisManager::Instance()->FillH1(10 + k, fEdepAbsor[k]);
      TotalEdep += fEdepAbsor[k];
    }
  }
  
  if (TotalEdep > 0.) {
    run->AddTotEdep(TotalEdep);  
    //G4AnalysisManager::Instance()->FillH1(2,TotalEdep);
  }
  
  //FADD

      //FADD: Save C11 positions
    std::ofstream outputFile("C11_positions_eventwise.csv", std::ios_base::app); // Open the file for appending

    if (outputFile.is_open())
    {
        // Loop through the stored C11 positions and write them to the file
        for (size_t i = 0; i < c11Positions.size(); ++i)
        {
            G4ThreeVector c11Position = c11Positions[i];
            outputFile << c11Position.x() << "," << c11Position.y() << "," << c11Position.z() << "\n";
        }

        // Close the file after writing
        outputFile.close();

        // Clear the vector for the next event
        c11Positions.clear();
    }
    else
    {
        // Handle the case where the file could not be opened
        G4cout << "Error: Unable to open the file for writing C11 positions." << G4endl;
    }
  
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

