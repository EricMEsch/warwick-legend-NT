#ifndef WLGDCrystalHit_h
#define WLGDCrystalHit_h 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"

/// Crystal hit class
///
/// It defines data members to store the energy deposit,
/// and position in a selected volume:

class WLGDCrystalHit : public G4VHit
{
public:
  WLGDCrystalHit();
  WLGDCrystalHit(const WLGDCrystalHit&);
  virtual ~WLGDCrystalHit();

  // operators
  const WLGDCrystalHit& operator=(const WLGDCrystalHit&);
  G4bool                operator==(const WLGDCrystalHit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  // methods from base class
  virtual void Draw();
  virtual void Print();

  // Set methods
  void SetTID(G4int tid) { fTid = tid; };
  void SetTime(G4double ti) { fTime = ti; };
  void SetWeight(G4double we) { fWeight = we; };
  void SetEdep(G4double de) { fEdep = de; };
  void SetPos(G4ThreeVector xyz) { fPos = xyz; };
  void SetWhichReentranceTube(G4int n) { fwhichReentranceTuber = n; };
  void SetWhichDetector(G4int n) { fwhichDetector = n; };
  void SetExcitedState(G4int n) { fExcitedState = n; };

  // Get methods
  G4int         GetTID() const { return fTid; };
  G4double      GetTime() const { return fTime; };
  G4double      GetWeight() const { return fWeight; };
  G4double      GetEdep() const { return fEdep; };
  G4ThreeVector GetPos() const { return fPos; };
  G4int         GetWhichReentranceTube() const { return fwhichReentranceTuber; };
  G4int         GetWhichDetector() const { return fwhichDetector; };
  G4int         GetExcitedState() const {return fExcitedState;};

private:
  G4int         fTid;
  // -- the global time of a hit
  G4double      fTime;
  // -- the current weight of the incoming particle
  G4double      fWeight;
  // -- the energy deposition of the hit
  G4double      fEdep;
  // -- the position of the hit
  G4ThreeVector fPos;
  // -- the number of the re-entrance tube in which the hit occured in
  G4int         fwhichReentranceTuber;
  // -- the detector number in which the event occured in
  G4int         fwhichDetector;
  G4int         fExcitedState;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<WLGDCrystalHit> WLGDCrystalHitsCollection;

extern G4ThreadLocal G4Allocator<WLGDCrystalHit>* WLGDCrystalHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* WLGDCrystalHit::operator new(size_t)
{
  if(!WLGDCrystalHitAllocator)
    WLGDCrystalHitAllocator = new G4Allocator<WLGDCrystalHit>;
  return (void*) WLGDCrystalHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void WLGDCrystalHit::operator delete(void* hit)
{
  WLGDCrystalHitAllocator->FreeSingle((WLGDCrystalHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
