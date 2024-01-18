#ifndef WLGDTrackingAction_h
#define WLGDTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "WLGDEventAction.hh"
#include "WLGDRunAction.hh"
#include <random>
#include <map>

class WLGDTrackingAction : public G4UserTrackingAction
{
public:
  WLGDTrackingAction();
  WLGDTrackingAction(WLGDEventAction* event, WLGDRunAction* run)
  {
    fEventAction = event;
    fRunAction   = run;
    WLGDTrackingAction();
  }
  virtual ~WLGDTrackingAction(){};

  virtual void PreUserTrackingAction(const G4Track*);
  virtual void PostUserTrackingAction(const  G4Track*);

  void KillSecundaries( const G4Track* aTrack);
  G4ThreeVector sample_direction();
  void ReplaceNeutronCaptureGammas( const G4Track* aTrack, int AtomicMass);

private:
  double           tmp_neutronXpos, tmp_neutronYpos, tmp_neutronZpos;
  double           tmp_neutronXmom, tmp_neutronYmom, tmp_neutronZmom;
  double           tmp_neutronTime;
  double           tmp_MuonXpos, tmp_MuonYpos, tmp_MuonZpos;
  double           tmp_MuonXmom, tmp_MuonYmom, tmp_MuonZmom;
  double           tmp_MuonEnergy;
  WLGDEventAction* fEventAction;
  WLGDRunAction*   fRunAction;
  std::random_device rd;
  std::ranlux24      generator;
};

#endif
