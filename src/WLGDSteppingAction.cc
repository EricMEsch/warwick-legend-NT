//
// Created by moritz on 3/4/21.
//

#include <iostream>

using namespace std;
#include "WLGDRunAction.hh"
#include "WLGDSteppingAction.hh"
#include "WLGDTrackingAction.hh"

#include "G4SystemOfUnits.hh"

#include "G4RunManager.hh"

WLGDSteppingAction::WLGDSteppingAction(WLGDEventAction* event, WLGDRunAction* run,
                                       WLGDDetectorConstruction* det)
{
  fEventAction          = event;
  fRunAction            = run;
  fDetectorConstruction = det;
  DefineCommands();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLGDSteppingAction::UserSteppingAction(const G4Step* aStep)
{
#define MostOuterRadiusTracking 0
  // Edit: 2021/03/05 by Moritz Neuberger
  // Adding tracking of amount of neutrons crossing the detectors
#if MostOuterRadiusTracking == 1
  if(aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "neutron")
  {
    G4double tmp_x = aStep->GetTrack()->GetPosition().getX() / m;
    G4double tmp_y = aStep->GetTrack()->GetPosition().getY() / m;
    fEventAction->UpdateMostOuterRadius(sqrt(tmp_x * tmp_x + tmp_y * tmp_y));
    if(aStep->GetTrack()->GetNextVolume())
    {
      auto physVol1 = aStep->GetTrack()->GetVolume();
      auto physVol2 = aStep->GetTrack()->GetNextVolume();
      if(physVol1->GetName() != "Ge_phys" && physVol2->GetName() == "Ge_phys")
      {
        if(fRunAction->getWriteOutGeneralNeutronInfo() == 1)
          fRunAction->increaseNumberOfCrossingNeutrons();
      }
    }
  }
#endif

  if(fRunAction->getReadMuonCrossingWLSR())
  {
    if((aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "mu-" || aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "mu+"))
    {
      if(aStep->GetTrack()->GetNextVolume())
      {
        auto physVol1 = aStep->GetTrack()->GetVolume();
        auto physVol2 = aStep->GetTrack()->GetNextVolume();
        if(physVol2->GetName() == "WLSR_LAr_physical")
          fEventAction->Add_Muon_WLSR_Edep(aStep->GetTotalEnergyDeposit() / eV);
        if((physVol1->GetName() != "WLSR_LAr_physical" && physVol2->GetName() == "WLSR_LAr_physical") || (physVol1->GetName() == "WLSR_LAr_physical" && physVol2->GetName() != "WLSR_LAr_physical"))
          fEventAction->Add_Muon_WLSR_intersect(aStep->GetPostStepPoint()->GetPosition().getX() / m,aStep->GetPostStepPoint()->GetPosition().getY() / m,aStep->GetPostStepPoint()->GetPosition().getZ() / m);
      }
    }
  }

  // Edit: 2021/04/07 by Moritz Neuberger
  // Adding total energy deposition inside LAr

  if(fDepositionInfo == 1)
  {
    if(aStep->GetTrack()->GetVolume()->GetName() == "World_phys")
      return;


    if((aStep->GetPostStepPoint()
          ->GetTouchable()
          ->GetVolume(0)
          ->GetLogicalVolume()
          ->GetName() == "Lar_log") ||
       aStep->GetPostStepPoint()
           ->GetTouchable()
           ->GetVolume(0)
           ->GetLogicalVolume()
           ->GetName() == "ULar_log" ||
       aStep->GetPostStepPoint()
           ->GetTouchable()
           ->GetVolume(0)
           ->GetLogicalVolume()
           ->GetName() == "Ge_log" ||
       aStep->GetPostStepPoint()
           ->GetTouchable()
           ->GetVolume(0)
           ->GetLogicalVolume()
           ->GetName() == "Water_log")
    {
      if(aStep->GetTotalEnergyDeposit() > 0)
      {

        // get position of deposition
        G4double tmp_x = aStep->GetTrack()->GetPosition().getX();
        G4double tmp_y = aStep->GetTrack()->GetPosition().getY();
        G4double tmp_z = aStep->GetTrack()->GetPosition().getZ();


        // get id of which reentrance tube
        /*
        G4int whichReentranceTube;
        if(abs(tmp_x) > abs(tmp_y) && tmp_x > 0)
          whichReentranceTube = 0;
        if(abs(tmp_x) > abs(tmp_y) && tmp_x < 0)
          whichReentranceTube = 2;
        if(abs(tmp_x) < abs(tmp_y) && tmp_y > 0)
          whichReentranceTube = 1;
        if(abs(tmp_x) < abs(tmp_y) && tmp_y < 0)
          whichReentranceTube = 3;

        if(fDetectorConstruction->GetGeometryName() == "hallA" ||
           aStep->GetPostStepPoint()
               ->GetTouchable()
               ->GetVolume(0)
               ->GetLogicalVolume()
               ->GetName() == "Lar_log" ||
           aStep->GetPostStepPoint()
               ->GetTouchable()
               ->GetVolume(0)
               ->GetLogicalVolume()
               ->GetName() == "Water_log")
          whichReentranceTube = 0;
        */

        // calculate total energy deposition in water tank for muon veto
        if(aStep->GetPostStepPoint()
             ->GetTouchable()
             ->GetVolume(0)
             ->GetLogicalVolume()
             ->GetName() == "Water_log")
        {
          // if(aStep->GetPostStepPoint()->GetGlobalTime() / us < 10.)
          fEventAction->IncreaseEdepWater_prompt(aStep->GetTotalEnergyDeposit() / eV);
          // else if(aStep->GetPostStepPoint()->GetGlobalTime() / ms < 1.)
          // fEventAction->IncreaseEdepWater_delayed(aStep->GetTotalEnergyDeposit() / eV);
          return;
        }
  
        
        G4int detector_number;
        if(fDetectorConstruction->GetGeometryName() == "baseline_large_reentrance_tube"){
          detector_number = aStep->GetPostStepPoint()->GetTouchable()->GetVolume(1)->GetCopyNo();
          whichReentranceTube = 0;
        }
        else{
          detector_number = aStep->GetPostStepPoint()->GetTouchable()->GetVolume(1)->GetCopyNo() + whichReentranceTube * 96;
        }

        // LAr veto

        G4int whichVolume = -1;
        if(aStep->GetPostStepPoint()
               ->GetTouchable()
               ->GetVolume(0)
               ->GetLogicalVolume()
               ->GetName() == "ULar_log" ||
           (fDetectorConstruction->GetGeometryName() == "hallA" &&
            aStep->GetTrack()->GetLogicalVolumeAtVertex()->GetName() == "Lar_log"))
        {
          whichVolume = 0;
          if(aStep->GetPostStepPoint()->GetGlobalTime() / us < 10.)
          {
            fEventAction->IncreaseLArEnergyDeposition(aStep->GetTotalEnergyDeposit() / eV,
                                                      whichReentranceTube);
          }  // prompt
          else
          {
            if(aStep->GetPostStepPoint()->GetGlobalTime() / ms < 1.)
            {
              fEventAction->IncreaseLArEnergyDeposition_delayed(
                aStep->GetTotalEnergyDeposit() / eV, whichReentranceTube);
            }  // delayed
            if(aStep->GetPostStepPoint()->GetGlobalTime() / s < 1. &&
               fRunAction->getWriteOutAdvancedMultiplicity())
            {
              fEventAction->IncreaseLArEnergyDeposition_delayed_long(
                aStep->GetTotalEnergyDeposit() / eV, whichReentranceTube);
            }  // long delayed
          }
          if(aStep->GetPostStepPoint()->GetGlobalTime() / s > 1. &&
             fRunAction->getWriteOutAdvancedMultiplicity())
          {
            fEventAction->IncreaseLArEnergyDeposition_after_delayed(
              aStep->GetTotalEnergyDeposit() / eV, whichReentranceTube);
          }  // after delayed
        }

        // Ge energy



        if(aStep->GetPostStepPoint()
             ->GetTouchable()
             ->GetVolume(0)
             ->GetLogicalVolume()
             ->GetName() == "Ge_log")
        {
          whichVolume = 1;
          if(aStep->GetPostStepPoint()->GetGlobalTime() / us < 10.)
          {
            if(aStep->GetPostStepPoint()
                 ->GetTouchable()
                 ->GetVolume(1)
                 ->GetLogicalVolume()
                 ->GetName() == "Layer_log")
            {

              fEventAction->IncreaseEdepPerDetector(
                detector_number,
                aStep->GetTotalEnergyDeposit() / eV);
              if(fRunAction->getWriteOutAdvancedMultiplicity())
              {
                if(fEventAction->GetIDListOfGdSiblingParticles().count(
                     aStep->GetTrack()->GetParentID()))
                  fEventAction->IncreaseEdepPerDetector_prompt_onlyGd(
                    detector_number,
                    aStep->GetTotalEnergyDeposit() / eV);
                else
                  fEventAction->IncreaseEdepPerDetector_prompt_woGd(
                    detector_number,
                    aStep->GetTotalEnergyDeposit() / eV);
              }  // w/ and w/o Gd info (redundant)
            }
            else
              G4cout
                << "Trying to access Layer_log for the prompt multiplicity but it is "
                << aStep->GetPostStepPoint()
                     ->GetTouchable()
                     ->GetVolume(1)
                     ->GetLogicalVolume()
                     ->GetName()
                << G4endl;
          }  // prompt
          else
          {
            if(aStep->GetPostStepPoint()->GetGlobalTime() / ms < 1.)
            {
              if(aStep->GetPostStepPoint()
                   ->GetTouchable()
                   ->GetVolume(1)
                   ->GetLogicalVolume()
                   ->GetName() == "Layer_log")
              {
                fEventAction->IncreaseEdepPerDetector_delayed(
                  detector_number,
                  aStep->GetTotalEnergyDeposit() / eV);

                if(fRunAction->getWriteOutAdvancedMultiplicity())
                {
                  if(fEventAction->GetIDListOfGdSiblingParticles().count(
                       aStep->GetTrack()->GetParentID()))
                    fEventAction->IncreaseEdepPerDetector_delayed_onlyGd(
                      detector_number,
                      aStep->GetTotalEnergyDeposit() / eV);
                  else
                    fEventAction->IncreaseEdepPerDetector_delayed_woGd(
                      detector_number,
                      aStep->GetTotalEnergyDeposit() / eV);
                }
              }
              else
                G4cout
                  << "Trying to access Layer_log for the delayed multiplicity but it is "
                  << aStep->GetPostStepPoint()
                       ->GetTouchable()
                       ->GetVolume(1)
                       ->GetLogicalVolume()
                       ->GetName()
                  << G4endl;
            }  // delayed
            if(aStep->GetPostStepPoint()->GetGlobalTime() / s < 1. &&
               fRunAction->getWriteOutAdvancedMultiplicity())
            {
              fEventAction->IncreaseEdepPerDetector_delayed_long(
                detector_number,
                aStep->GetTotalEnergyDeposit() / eV);
            }  // long delayed
          }
          if(aStep->GetPostStepPoint()->GetGlobalTime() / s > 1. &&
             fRunAction->getWriteOutAdvancedMultiplicity())
          {
            fEventAction->IncreaseGeEnergyDeposition_after_delayed(
              aStep->GetTotalEnergyDeposit() / eV, whichReentranceTube);
          }  // after delayed
        }




        if(fRunAction->getIndividualGeDepositionInfo())
        {
          if(fEventAction->GetIDListOfGe77().count(aStep->GetTrack()->GetParentID()))
          {
            

            fEventAction->AddGe77mGammaEmission_timing(
              aStep->GetPostStepPoint()->GetGlobalTime() / s);
            fEventAction->AddGe77mGammaEmission_x(
              aStep->GetPostStepPoint()->GetGlobalTime() / m);
            fEventAction->AddGe77mGammaEmission_y(
              aStep->GetPostStepPoint()->GetGlobalTime() / m);
            fEventAction->AddGe77mGammaEmission_z(
              aStep->GetPostStepPoint()->GetGlobalTime() / m);
            fEventAction->AddGe77mGammaEmission_edep(aStep->GetTotalEnergyDeposit() / eV);
            fEventAction->AddGe77mGammaEmission_id(aStep->GetTrack()->GetTrackID());
            fEventAction->AddGe77mGammaEmission_type(
              aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding());
            fEventAction->AddGe77mGammaEmission_whichGe77(
              aStep->GetTrack()->GetParentID());
            int whichVolume = -3;
            if(aStep->GetPostStepPoint()
                 ->GetTouchable()
                 ->GetVolume()
                 ->GetLogicalVolume()
                 ->GetName() == "Lar_log")
              whichVolume = -2;
            if(aStep->GetPostStepPoint()
                 ->GetTouchable()
                 ->GetVolume(0)
                 ->GetLogicalVolume()
                 ->GetName() == "ULar_log")
              whichVolume = -1;
            if(aStep->GetPostStepPoint()
                 ->GetTouchable()
                 ->GetVolume(0)
                 ->GetLogicalVolume()
                 ->GetName() == "Ge_log")
              whichVolume =
                detector_number;
            fEventAction->AddGe77mGammaEmission_whichVolume(whichVolume);
          }

          if(fEventAction->GetIDListOfGe77SiblingParticles().count(
               aStep->GetTrack()->GetParentID()))
          {
            fEventAction->AddGe77Siblings_timing(
              aStep->GetPostStepPoint()->GetGlobalTime() / s);
            fEventAction->AddGe77Siblings_x(
              aStep->GetPostStepPoint()->GetPosition().getX() / m);
            fEventAction->AddGe77Siblings_y(
              aStep->GetPostStepPoint()->GetPosition().getY() / m);
            fEventAction->AddGe77Siblings_z(
              aStep->GetPostStepPoint()->GetPosition().getZ() / m);
            fEventAction->AddGe77Siblings_edep(aStep->GetTotalEnergyDeposit() / eV);
            fEventAction->AddGe77Siblings_id(aStep->GetTrack()->GetTrackID());
            fEventAction->AddGe77Siblings_type(
              aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding());
            int whichVolume = -3;
            if(aStep->GetPostStepPoint()
                 ->GetTouchable()
                 ->GetVolume()
                 ->GetLogicalVolume()
                 ->GetName() == "Lar_log")
              whichVolume = -2;
            if(aStep->GetPostStepPoint()
                 ->GetTouchable()
                 ->GetVolume(0)
                 ->GetLogicalVolume()
                 ->GetName() == "ULar_log")
              whichVolume = -1;
            if(aStep->GetPostStepPoint()
                 ->GetTouchable()
                 ->GetVolume(0)
                 ->GetLogicalVolume()
                 ->GetName() == "Ge_log")
              whichVolume =
                detector_number;
            fEventAction->AddGe77Siblings_whichVolume(whichVolume);
          }
        }  // individual deposition of sibling or secundary Ge77 interactions

        if(fRunAction->getIndividualGdDepositionInfo())
        {
          if(fEventAction->GetIDListOfGdSiblingParticles().count(
               aStep->GetTrack()->GetParentID()))
          {
            // G4cout << " ______________________________________________ " << G4endl <<
            // aStep->GetTrack()->GetTrackID() << " " <<
            // aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding() << " "  <<
            // aStep->GetTrack()->GetParentID() << G4endl << "
            // ______________________________________________ " << G4endl;
            fEventAction->AddGdSiblings_timing(
              aStep->GetPostStepPoint()->GetGlobalTime() / s);
            fEventAction->AddGdSiblings_x(
              aStep->GetPostStepPoint()->GetPosition().getX() / m);
            fEventAction->AddGdSiblings_y(
              aStep->GetPostStepPoint()->GetPosition().getY() / m);
            fEventAction->AddGdSiblings_z(
              aStep->GetPostStepPoint()->GetPosition().getZ() / m);
            fEventAction->AddGdSiblings_edep(aStep->GetTotalEnergyDeposit() / eV);
            fEventAction->AddGdSiblings_id(aStep->GetTrack()->GetTrackID());
            fEventAction->AddGdSiblings_type(
              aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding());
            int whichVolume = -4;
            if(aStep->GetPostStepPoint()
                 ->GetTouchable()
                 ->GetVolume()
                 ->GetLogicalVolume()
                 ->GetName() == "Water_log")
              whichVolume = -3;
            if(aStep->GetPostStepPoint()
                 ->GetTouchable()
                 ->GetVolume()
                 ->GetLogicalVolume()
                 ->GetName() == "Lar_log")
              whichVolume = -2;
            if(aStep->GetPostStepPoint()
                 ->GetTouchable()
                 ->GetVolume(0)
                 ->GetLogicalVolume()
                 ->GetName() == "ULar_log")
              whichVolume = -1;
            if(aStep->GetPostStepPoint()
                 ->GetTouchable()
                 ->GetVolume(0)
                 ->GetLogicalVolume()
                 ->GetName() == "Ge_log")
              whichVolume =
                detector_number;
            fEventAction->AddGdSiblings_whichVolume(whichVolume);
          }
        }  // individual Gd interactions


        if(aStep->GetPostStepPoint()->GetGlobalTime() / s > 1 &&
           fAllowForLongTimeEmissionReadout == 0)
          return;  // skip all interactiosn >1s

        if(fIndividualDepositionInfo == 0 &&
           fRunAction->getIndividualGeDepositionInfo() == 0)
          return;

        if(fRunAction->getIndividualGeDepositionInfo() == 1 &&
           aStep->GetPostStepPoint()
               ->GetTouchable()
               ->GetVolume(1)
               ->GetLogicalVolume()
               ->GetName() != "Layer_log")
          return;
        // all individual interactions
        {
          fEventAction->AddIndividualEnergyDeposition_Timing(
            aStep->GetPostStepPoint()->GetGlobalTime() / (1000 * ns));
          fEventAction->AddIndividualEnergyDeposition_Energy(
            aStep->GetTotalEnergyDeposit() / eV);
          fEventAction->AddIndividualEnergyDeposition_ReentranceTube(whichReentranceTube);
          fEventAction->AddIndividualEnergyDeposition_Position_x(tmp_x / m);
          fEventAction->AddIndividualEnergyDeposition_Position_y(tmp_y / m);
          fEventAction->AddIndividualEnergyDeposition_Position_z(tmp_z / m);
          fEventAction->AddIndividualEnergyDeposition_LArOrGe(whichVolume);
          fEventAction->AddIndividualEnergyDeposition_ID(aStep->GetTrack()->GetTrackID());
          fEventAction->AddIndividualEnergyDeposition_ParentID(aStep->GetTrack()->GetParentID());
          fEventAction->AddIndividualEnergyDeposition_Type(
            aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding());
          int tmp = -1;
          if(aStep->GetPostStepPoint()
               ->GetTouchable()
               ->GetVolume(1)
               ->GetLogicalVolume()
               ->GetName() == "Layer_log")
            tmp = detector_number;
          fEventAction->AddIndividualEnergyDeposition_DetectorNumber(tmp);
          
          int VolumeNumber = -1;
          if(aStep->GetPostStepPoint()
                 ->GetTouchable()
                 ->GetVolume(0)
                 ->GetLogicalVolume()
                 ->GetName() == "Ge_log")
            VolumeNumber = 0;
          
          if(aStep->GetPostStepPoint()
                 ->GetTouchable()
                 ->GetVolume(0)
                 ->GetLogicalVolume()
                 ->GetName() == "ULar_log")
            VolumeNumber = 1;

          if(aStep->GetPostStepPoint()
                 ->GetTouchable()
                 ->GetVolume(0)
                 ->GetLogicalVolume()
                 ->GetName() == "Lar_log")
            VolumeNumber = 2;

          fEventAction->AddIndividualEnergyDeposition_VolumeNumber(VolumeNumber);
        }

      }
    }
  }
  
  // EDIT: by Eric Esch
  // Neutron Capture for Neutrontagger
  G4VProcess* process = const_cast<G4VProcess*>(aStep->GetPostStepPoint()->GetProcessDefinedStep());
  if(aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "neutron")
  {
    if(aStep->GetTrack()->GetTrackStatus() == fStopAndKill)
    { 
      if(process->GetProcessName() == "biasWrapper(nCapture)" || process->GetProcessName() == "nCapture")
      {
        // Check if the process is a biasing wrapper
        G4BiasingProcessInterface* biasingInterface = dynamic_cast<G4BiasingProcessInterface*>(process);
        G4HadronicProcess* hproc = nullptr;
        if(biasingInterface)
        {
          // Get the underlying process from the biasing wrapper
          G4VProcess* underlyingProcess = biasingInterface->GetWrappedProcess();
          // Check if the underlying process is a hadronic process
          hproc = dynamic_cast<G4HadronicProcess*>(underlyingProcess);
        }
        else
        {
          // If it's not a biasing wrapper, check if it's a direct hadronic process
          hproc = dynamic_cast<G4HadronicProcess*>(process);
        }
        //Get Target Isotope name and information
        const G4Isotope* target = NULL;
        if (hproc) target = hproc->GetTargetIsotope();
        G4String targetName = "XXXX";
        if (target) 
        {
          fEventAction->Add_NCaptureIsotopes(target->GetName());
          G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();
          fEventAction->Add_NCaptureX(position.getX());
          fEventAction->Add_NCaptureY(position.getY());
          fEventAction->Add_NCaptureZ(position.getZ());
          fEventAction->Add_NCaptureEnergy(aStep->GetTotalEnergyDeposit() / eV);
          fEventAction->Add_NCaptureTime(aStep->GetPostStepPoint()->GetGlobalTime() / ns);
        }
      }
    }
  }
  

}

void WLGDSteppingAction::GetDepositionInfo(G4int answer) { fDepositionInfo = answer; }
void WLGDSteppingAction::GetIndividualDepositionInfo(G4int answer)
{
  fIndividualDepositionInfo = answer;
}
void WLGDSteppingAction::AllowForLongTimeEmissionReadout(G4int answer)
{
  fAllowForLongTimeEmissionReadout = answer;
}

void WLGDSteppingAction::DefineCommands()
{
  // Define geometry command directory using generic messenger class
  fStepMessenger = new G4GenericMessenger(this, "/WLGD/step/",
                                          "Commands for controlling stepping action");

  // switch command
  fStepMessenger
    ->DeclareMethod("getDepositionInfo", &WLGDSteppingAction::GetDepositionInfo)
    .SetGuidance(
      "Set whether to obtain energy deposition information inside reentrance tubes")
    .SetGuidance("0 = don't")
    .SetGuidance("1 = do")
    .SetCandidates("0 1")
    .SetDefaultValue("0");

  fStepMessenger
    ->DeclareMethod("getIndividualDepositionInfo",
                    &WLGDSteppingAction::GetIndividualDepositionInfo)
    .SetGuidance("Set whether to obtain individual energy deposition information inside "
                 "reentrance tubes")
    .SetGuidance("0 = don't")
    .SetGuidance("1 = do")
    .SetCandidates("0 1")
    .SetDefaultValue("0");

  fStepMessenger
    ->DeclareMethod("AllowForLongTimeEmissionReadout",
                    &WLGDSteppingAction::AllowForLongTimeEmissionReadout)
    .SetGuidance("Set whether to individual deposition information can be gathered from "
                 "longer than 1sec after start of event.")
    .SetGuidance("0 = don't")
    .SetGuidance("1 = do")
    .SetCandidates("0 1")
    .SetDefaultValue("0");
}
