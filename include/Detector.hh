#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"

class SensitiveDetector : public G4VSensitiveDetector {
public:
    SensitiveDetector(G4String);
    virtual ~SensitiveDetector();

    //virtual void Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent);
    //virtual void EndOfEvent(G4HCofThisEvent* hitsCollectionOfThisEvent);

private:
    G4int fHitsCollectionID;
    virtual G4bool ProcessHits(G4Step* , G4TouchableHistory* );
    virtual G4bool IsGeneratedInScoringVolume(G4Step* aStep);

};

#endif // DETECTOR_HH
