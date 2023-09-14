#include "WLGDCrystalSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLGDCrystalSD::WLGDCrystalSD(const G4String& name, const G4String& hitsCollectionName,
                             G4String setupName)
: G4VSensitiveDetector(name)
, fHitsCollection(NULL)
, fGeometryName("baseline")
{
  collectionName.insert(hitsCollectionName);
  fGeometryName = setupName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLGDCrystalSD::~WLGDCrystalSD() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLGDCrystalSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection =
    new WLGDCrystalHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce

  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool WLGDCrystalSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();

  if(edep == 0.)
    return false;

  // particle filter on Ge-77
  auto iZ = aStep->GetTrack()->GetDefinition()->GetAtomicNumber();
  auto iA = aStep->GetTrack()->GetDefinition()->GetAtomicMass();
  if(!(iZ == 32 && iA == 77))
    return false;

  if(std::find(fListOfGe77IDs.begin(), fListOfGe77IDs.end(), aStep->GetTrack()->GetParentID()) != fListOfGe77IDs.end())
    return false;
  
  fListOfGe77IDs.push_back(aStep->GetTrack()->GetTrackID());


  WLGDCrystalHit* newHit = new WLGDCrystalHit();

  newHit->SetTID(aStep->GetTrack()->GetTrackID());
  newHit->SetTime(aStep->GetTrack()->GetGlobalTime());
  newHit->SetWeight(aStep->GetTrack()->GetWeight());
  newHit->SetEdep(edep);
  newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
  newHit->SetExcitedState(aStep->GetTrack()->GetDefinition()->GetPDGEncoding()%10);

  // -- added the output for the re-entrance tube and the detector number 
  G4double tmp_x = aStep->GetTrack()->GetPosition().getX();
  G4double tmp_y = aStep->GetTrack()->GetPosition().getY();
  G4int    whichReentranceTube;
  if(abs(tmp_x) > abs(tmp_y) && tmp_x > 0)
    whichReentranceTube = 0;
  if(abs(tmp_x) > abs(tmp_y) && tmp_x < 0)
    whichReentranceTube = 2;
  if(abs(tmp_x) < abs(tmp_y) && tmp_y > 0)
    whichReentranceTube = 1;
  if(abs(tmp_x) < abs(tmp_y) && tmp_y < 0)
    whichReentranceTube = 3;
  if(fGeometryName == "hallA")
    whichReentranceTube = 0;
  newHit->SetWhichReentranceTube(whichReentranceTube);
  if(fGeometryName == "baseline_large_reentrance_tube"){
    newHit->SetWhichDetector(
      aStep->GetPostStepPoint()->GetTouchable()->GetVolume(1)->GetCopyNo());
      whichReentranceTube = 0;
  }
  else
    newHit->SetWhichDetector(
      aStep->GetPostStepPoint()->GetTouchable()->GetVolume(1)->GetCopyNo() +
      whichReentranceTube * 96);
  fHitsCollection->insert(newHit);

  G4Track* nonConstTrack = const_cast<G4Track*>(aStep->GetTrack());
  //nonConstTrack->SetKineticEnergy(0);
  nonConstTrack->SetTrackStatus(fStopAndKill);

  const auto secondaries = aStep->GetSecondaryInCurrentStep();
  for (const auto& secondary : *secondaries) {
      G4Track* nonConstTrack = const_cast<G4Track*>(secondary);
      //nonConstTrack->SetKineticEnergy(0);
      nonConstTrack->SetTrackStatus(fStopAndKill);
  }

  return true;
}

void WLGDCrystalSD::EndOfEvent(G4HCofThisEvent*)
{
  fListOfGe77IDs.clear();
  if(verboseLevel > 1)
  {
    G4int nofHits = fHitsCollection->entries();
    G4cout << G4endl << "-------->Hits Collection: in this event they are " << nofHits
           << " hits in the crystal: " << G4endl;
  }
}
