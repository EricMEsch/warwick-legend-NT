#include "WLGDActionInitialization.hh"
#include "WLGDEventAction.hh"
#include "WLGDPrimaryGeneratorAction.hh"
#include "WLGDRunAction.hh"
#include "WLGDStackingAction.hh"
#include "WLGDSteppingAction.hh"
#include "WLGDTrackingAction.hh"
#include "WLGDPetersGammaCascadeReader.hh"

WLGDActionInitialization::WLGDActionInitialization(WLGDDetectorConstruction* det,
                                                   G4String                  name)
: G4VUserActionInitialization()
, fDet(det)
, foutname(std::move(name))
{}

WLGDActionInitialization::~WLGDActionInitialization() = default;

void WLGDActionInitialization::BuildForMaster() const
{
  auto event = new WLGDEventAction;
  SetUserAction(new WLGDRunAction(event, foutname));
}

void WLGDActionInitialization::Build() const
{
  // forward detector
  // Edit: Moritz Neuberger
  // added connection between the event action and run/tracking/stepping action classes
  // added connection between run action and tracking/stepping action class
  // Yeah, I know its not pretty but for my first G4 project, I think it's ok. Sorry for the mess. If you have any problems, just write me on Slack
  SetUserAction(new WLGDPrimaryGeneratorAction(fDet));
  auto event = new WLGDEventAction();
  SetUserAction(event);
  auto run = new WLGDRunAction(event, foutname);
  SetUserAction(run);
  SetUserAction(new WLGDTrackingAction(event, run));
  SetUserAction(new WLGDSteppingAction(event, run, fDet));
  SetUserAction(new WLGDStackingAction);
  WLGDPetersGammaCascadeReader::GetInstance();
}
