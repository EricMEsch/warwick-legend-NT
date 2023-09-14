#ifndef WLGDCrystalSD_h
#define WLGDCrystalSD_h 1

#include "G4VSensitiveDetector.hh"

#include "WLGDCrystalHit.hh"

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// WLGDCrystal sensitive detector class
///
/// The hits are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step. A hit is created with each step with non zero
/// energy deposit.

class WLGDCrystalSD : public G4VSensitiveDetector
{
public:
  WLGDCrystalSD(const G4String& name, const G4String& hitsCollectionName,
                G4String setupName);
  virtual ~WLGDCrystalSD();

  // methods from base class
  virtual void   Initialize(G4HCofThisEvent* hitCollection);
  virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
  virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

private:
  WLGDCrystalHitsCollection* fHitsCollection;
  G4String                   fGeometryName;
  std::vector<G4int>         fListOfGe77IDs;
};

#endif
