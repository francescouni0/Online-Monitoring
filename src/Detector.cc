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
   
       G4double x = localPos.x();
       G4double y = localPos.y();
       G4double z = localPos.z();
       G4int evt= G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

       G4int parentID = track->GetParentID();
   
       G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
       analysisManager->FillNtupleDColumn(0, x);
       analysisManager->FillNtupleDColumn(1, y);
       analysisManager->FillNtupleDColumn(2, z);
       analysisManager->FillNtupleDColumn(3, initialEnergy);
       analysisManager->FillNtupleDColumn(4, TOF); // TOF information
       analysisManager->FillNtupleSColumn(5, particleName);
       analysisManager->FillNtupleSColumn(6, processName);
       analysisManager->FillNtupleIColumn(7, evt);
       analysisManager->FillNtupleIColumn(8, parentID);
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
        (creationPoint.y() >= -20 * cm) && (creationPoint.y() <= 20 * cm) &&  // Half-length in y direction
        (creationPoint.z() >= -20 * cm) && (creationPoint.z() <= 20 * cm);    // Half-length in z direction    
    return isInDetectorVolume;
}