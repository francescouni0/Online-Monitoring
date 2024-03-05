#include "Detector.hh"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4TrackingManager.hh"

SensitiveDetector::SensitiveDetector(G4String name)
: G4VSensitiveDetector(name)
{
    G4String HCname;
    collectionName.insert(HCname="detectorCollection");
}

SensitiveDetector::~SensitiveDetector() {}

G4bool SensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory* /*ROhist*/)
{
    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();


    G4Track* track = aStep->GetTrack();
    if (IsGeneratedInScoringVolume(aStep)) {
       G4String particleName = track->GetDefinition()->GetParticleName();
       G4String processName;
       
       if (track->GetCreatorProcess() != nullptr) {
           processName = track->GetCreatorProcess()->GetProcessName();
       } else {
           processName = "Unknown Process";
       }    
       G4ThreeVector pos = preStepPoint->GetPosition();
       G4double initialEnergy = track->GetVertexKineticEnergy() / MeV;
   
       // Calculate cylindrical coordinates
   
           // Calculate TOF
       G4double TOF = preStepPoint->GetGlobalTime(); // This is the time of flight
   
       const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
   
       G4int copyNo = touchable->GetVolume()->GetCopyNo();
   
       G4VPhysicalVolume* physicalVolume = touchable->GetVolume();
   
       G4ThreeVector localPos = physicalVolume->GetTranslation();
   
       G4double x = pos.x();
       G4double y = pos.y();
       G4double z = pos.z();
       G4int evt= G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

       G4int parentID = track->GetParentID();




   
       G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
       analysisManager->FillNtupleDColumn(2, x);
       analysisManager->FillNtupleDColumn(3, y);
       analysisManager->FillNtupleDColumn(4, z);
       analysisManager->FillNtupleDColumn(5, initialEnergy);
       analysisManager->FillNtupleDColumn(6, TOF); // TOF information
       analysisManager->FillNtupleSColumn(7, particleName);
       analysisManager->FillNtupleSColumn(8, processName);
       analysisManager->FillNtupleIColumn(9, evt);
       analysisManager->FillNtupleIColumn(10, parentID);
       analysisManager->AddNtupleRow(0);
   
       // Set the track status to stop and kill
       track->SetTrackStatus(fStopAndKill);

    }

    return true;
}

G4bool SensitiveDetector::IsGeneratedInScoringVolume(G4Step* aStep) {
    G4Track* track = aStep->GetTrack();
    
    // Get the creation point of the particle
    G4ThreeVector creationPoint = track->GetVertexPosition();
    
    // Check if the creation point is within the detector volume
    bool isInDetectorVolume =
        (creationPoint.x() >= -300 * cm) && (creationPoint.x() <=300 * cm) &&  // Half-length in x direction
        (creationPoint.y() >= -30 * cm) && (creationPoint.y() <= 30 * cm) &&  // Half-length in y direction
        (creationPoint.z() >= -30 * cm) && (creationPoint.z() <= 30 * cm);    // Half-length in z direction    
    return isInDetectorVolume;
}