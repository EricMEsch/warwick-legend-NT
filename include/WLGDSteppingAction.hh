//
// Created by moritz on 3/4/21.
//

#ifndef WARWICK_LEGEND_WLGDSTEPPINGACTION_HH
#define WARWICK_LEGEND_WLGDSTEPPINGACTION_HH

#include "G4GenericMessenger.hh"
#include "G4UserSteppingAction.hh"
#include "WLGDDetectorConstruction.hh"
#include "WLGDRunAction.hh"
#include "G4HadronicProcess.hh"
#include "G4BiasingProcessInterface.hh"
#include "globals.hh"

class WLGDTrackingAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class WLGDSteppingAction : public G4UserSteppingAction
{
public:
  WLGDSteppingAction() = default;
  WLGDSteppingAction(WLGDEventAction* event, WLGDRunAction* run,
                     WLGDDetectorConstruction* det);
  virtual ~WLGDSteppingAction() = default;

  virtual void UserSteppingAction(const G4Step*);
  void         GetDepositionInfo(G4int answer);
  void         GetIndividualDepositionInfo(G4int answer);
  void         AllowForLongTimeEmissionReadout(G4int answer);
  void         DefineCommands();
  void         OutputCaptureIntegerFromName(G4String targetName);

private:
  WLGDRunAction*            fRunAction;
  WLGDEventAction*          fEventAction;
  WLGDDetectorConstruction* fDetectorConstruction;
  G4GenericMessenger*       fStepMessenger;
  G4int                     fDepositionInfo                  = 0;
  G4int                     fIndividualDepositionInfo        = 0;
  G4int                     fAllowForLongTimeEmissionReadout = 0;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
