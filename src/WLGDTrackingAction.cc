#include "WLGDTrackingAction.hh"
// #include "WLGDTrackInformation.hh"
#include "WLGDTrajectory.hh"
#include "WLGDPetersGammaCascadeReader.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4TrackingManager.hh"
#include "G4UnitsTable.hh"
#include "G4IonTable.hh"
#include "G4Gamma.hh"
#include "G4StackManager.hh"
#include <map>

WLGDTrackingAction::WLGDTrackingAction(){
  generator.seed(rd());  // set a random seed
}

void WLGDTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // Create trajectory for track if requested
  if(fpTrackingManager->GetStoreTrajectory() > 0)
  {
    fpTrackingManager->SetTrajectory(new WLGDTrajectory(aTrack));
  }

  // add Ge77 events to ListOfGe77
  if(aTrack->GetParticleDefinition()->GetAtomicMass() == 77 &&
     aTrack->GetParticleDefinition()->GetPDGCharge() == 32)
  {
    fEventAction->AddIDListOfGe77(aTrack->GetTrackID());
  }

  // Adding tracking of initial muons
  if(aTrack->GetParticleDefinition()->GetParticleName() == "mu-" ||
     aTrack->GetParticleDefinition()->GetParticleName() == "mu+")
  {
    auto tmp_vector = aTrack->GetVertexPosition();
    tmp_MuonXpos    = tmp_vector.getX() / m;
    tmp_MuonYpos    = tmp_vector.getY() / m;
    tmp_MuonZpos    = tmp_vector.getZ() / m;
    tmp_vector      = aTrack->GetMomentumDirection();
    tmp_MuonXmom    = tmp_vector.getX();
    tmp_MuonYmom    = tmp_vector.getY();
    tmp_MuonZmom    = tmp_vector.getZ();
    tmp_MuonZmom    = tmp_vector.getZ();
    tmp_MuonEnergy  = aTrack->GetKineticEnergy();
  }

  // Edit: 2021/03/30 by Moritz Neuberger
  // Adding tracking of neutrons being later captured by Ge-76 as well as general produced
  // in LAr
  if(aTrack->GetParticleDefinition()->GetParticleName() == "neutron")
  {
    {
      auto tmp_vector = aTrack->GetVertexPosition();
      tmp_neutronXpos = tmp_vector.getX() / m;
      tmp_neutronYpos = tmp_vector.getY() / m;
      tmp_neutronZpos = tmp_vector.getZ() / m;
      tmp_vector      = aTrack->GetMomentumDirection();
      tmp_neutronXmom = tmp_vector.getX();
      tmp_neutronYmom = tmp_vector.getY();
      tmp_neutronZmom = tmp_vector.getZ();
    }  // readout of neutron production info

    // if all neutron info is readout, sent to output
    if(fRunAction->getWriteOutAllNeutronInfoRoot() == 1)
    {
      fEventAction->AddNeutronxLoc(tmp_neutronXpos);
      fEventAction->AddNeutronyLoc(tmp_neutronYpos);
      fEventAction->AddNeutronzLoc(tmp_neutronZpos);
      fEventAction->AddNeutronxMom(tmp_neutronXmom);
      fEventAction->AddNeutronyMom(tmp_neutronYmom);
      fEventAction->AddNeutronzMom(tmp_neutronZmom);
      fEventAction->AddNeutronTime(tmp_neutronTime);
      // fEventAction->AddNeutronEkin(aTrack->GetKineticEnergy() / eV);
    }

    // initial value for furthest position of neutron away from center (for testing, can
    // be removed)
    fEventAction->SetMostOuterRadius(
      sqrt(tmp_neutronXpos * tmp_neutronXpos + tmp_neutronYpos * tmp_neutronYpos));

    // increase number of neutrons produced
    if(fRunAction->getWriteOutGeneralNeutronInfo() == 1)
      fEventAction->IncreaseByOne_NeutronInEvent();

    // if production of neutrons in Ge76 nC is recorded, channel to output
    if(fRunAction->getWriteOutNeutronProductionInfo() == 1)
    {
      fRunAction->increaseTotalNumberOfNeutronsInLAr();
      fRunAction->addEventNumber(
        G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID());
      fRunAction->addCoordinatsToFile(tmp_neutronXpos, tmp_neutronYpos, tmp_neutronZpos);
      fRunAction->addMomentumToFile(tmp_neutronXmom, tmp_neutronYmom, tmp_neutronZmom);
      fRunAction->addEnergyToFile(aTrack->GetKineticEnergy() / eV);
      fRunAction->addParentParticleType(
        fEventAction->neutronProducerMap.find(aTrack->GetParentID())->second);
    }
  }
}

void WLGDTrackingAction::KillSecundaries(const G4Track* aTrack){
  const auto secondaries = aTrack->GetStep()->GetSecondaryInCurrentStep();
  for (const auto& secondary : *secondaries) {
      G4Track* nonConstTrack = const_cast<G4Track*>(secondary);
      //nonConstTrack->SetKineticEnergy(0);
      nonConstTrack->SetTrackStatus(fStopAndKill);
  }

  //G4Step* nonConstStep = const_cast<G4Step*>(aTrack->GetStep());
  //nonConstStep->DeleteSecondaryVector();
  //G4TrackVector* new_secondary = nonConstStep->NewSecondaryVector();
  //nonConstStep->SetSecondary(new_secondary);
}

G4ThreeVector WLGDTrackingAction::sample_direction(){
  std::uniform_real_distribution<> rndm(0.0, 1.0);
  double phi = CLHEP::twopi * rndm(generator);
  double theta = CLHEP::pi * rndm(generator);
  G4double costheta = std::cos(theta);  // get a random number
  G4double sintheta = std::sin(theta);
  G4double sinphi = std::sin(phi);
  G4double cosphi = std::cos(phi);

  G4double      px = -sintheta * cosphi;
  G4double      py = -sintheta * sinphi;
  G4double      pz = -costheta;  // default downwards: pz = -1.0

  G4ThreeVector mom_direction(px,py,pz);
  return mom_direction;
}

void WLGDTrackingAction::ReplaceNeutronCaptureGammas(const G4Track* aTrack){


  G4StackManager* stackManager = G4EventManager::GetEventManager()->GetStackManager();

  G4double ekin = aTrack->GetKineticEnergy()/keV;
  WLGDPetersGammaCascadeReader* pgcr = WLGDPetersGammaCascadeReader::GetInstance();
  GammaCascadeLine input = pgcr->GetNextEntry(ekin);

  G4ThreeVector location = aTrack->GetStep()->GetPostStepPoint()->GetPosition();

  G4IonTable* theTable = G4IonTable::GetIonTable();
  G4ParticleDefinition * particleDef_nuc = theTable->GetIon(static_cast<G4int>(32), static_cast<G4int>(77), (double)(input.em * keV));
  G4DynamicParticle* particle_nuc = new G4DynamicParticle(particleDef_nuc, sample_direction());
  G4Track* secondary_nuc = new G4Track(particle_nuc, aTrack->GetGlobalTime(), location);
  secondary_nuc->SetParentID(aTrack->GetTrackID());
  secondary_nuc->SetTrackID(aTrack->GetTrackID()+1);
  secondary_nuc->SetTouchableHandle(aTrack->GetTouchableHandle());
  secondary_nuc->SetWeight(1.0);
  //G4cout << input.em << " - " << secondary_nuc->GetTrackID() << G4endl;
  stackManager->PushOneTrack(secondary_nuc);
  //G4Track* nonConstTrack = const_cast<G4Track*>(aTrack);
  //nonConstTrack->AddSecondary(secondary_nuc);
  //G4Step* nonConstStep = const_cast<G4Step*>(aTrack->GetStep());
  //G4ParticleChange* nonConstParticleChange = const_cast<G4ParticleChange*>(aTrack->GetParticleChange());
  //nonConstParticleChange->AddSecondary(secondary_nuc);
  //nonConstStep->GetfSecondary()->push_back(secondary_nuc);
  //aTrack->AddSecondary(secondary_nuc);    

  for(int i = 0; i < input.m; i++){

    G4ParticleDefinition* particleDef_gamma = G4Gamma::Gamma();
    G4DynamicParticle* particle_gamma = new G4DynamicParticle(particleDef_gamma, input.eg[i]*keV,sample_direction() );
    G4Track* secondary_gamma = new G4Track(particle_gamma, aTrack->GetGlobalTime(), location);
    secondary_gamma->SetParentID(aTrack->GetTrackID());
    secondary_gamma->SetTouchableHandle(aTrack->GetTouchableHandle());
    secondary_gamma->SetWeight(1.0);
    secondary_gamma->SetKineticEnergy(input.eg[i]*keV);
    secondary_gamma->SetTrackID(aTrack->GetTrackID()+2+i);

    //G4cout << input.eg[i] << " - " << secondary_gamma->GetTrackID() << " - " << secondary_gamma->GetTotalEnergy ()  << G4endl;
    stackManager->PushOneTrack(secondary_gamma);
    
    //nonConstStep->GetfSecondary()->push_back(secondary_gamma);
    //nonConstParticleChange->AddSecondary(secondary_gamma);
    //nonConstTrack->AddSecondary(secondary_gamma);
    //particleChange->AddSecondary(secondary_gamma);
    //aTrack->AddSecondary(secondary_gamma);
  }
}

void WLGDTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{

  // for tracking of particles creatd in Gd interactions
  if(fRunAction->getIndividualGdDepositionInfo())
  {
    if(fEventAction->GetIDListOfGdSiblingParticles().count(aTrack->GetParentID()))
    {
      fEventAction->AddIDListOfGdSiblingParticles(aTrack->GetTrackID());
    }
  }

  // for tracking of sibling and secundary particles of Ge76 nC
  if(fEventAction->GetIDListOfGe77SiblingParticles().count(aTrack->GetParentID()) &&
     fRunAction->getIndividualGeDepositionInfo())
  {
    fEventAction->AddIDListOfGe77SiblingParticles(aTrack->GetTrackID());
  }
  if(fEventAction->GetIDListOfGe77().count(aTrack->GetParentID()) &&
     fRunAction->getIndividualGeDepositionInfo())
  {
    fEventAction->AddIDListOfGe77(aTrack->GetTrackID());
  }

  // write to output muon information
  if(aTrack->GetTrackID() == 1)
  {
    fEventAction->AddMuonxLoc(tmp_MuonXpos);
    fEventAction->AddMuonyLoc(tmp_MuonYpos);
    fEventAction->AddMuonzLoc(tmp_MuonZpos);
    fEventAction->AddMuonxMom(tmp_MuonXmom);
    fEventAction->AddMuonyMom(tmp_MuonYmom);
    fEventAction->AddMuonzMom(tmp_MuonZmom);
    fEventAction->AddMuonEnergy(tmp_MuonEnergy);
  }

  if(fRunAction->getWriteOutNeutronProductionInfo() == 1)
  {
    G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
    if(secondaries != nullptr)
    {
      size_t nSeco = secondaries->size();
      if(nSeco > 0)
      {
        for(size_t i = 0; i < nSeco; i++)
        {
          // Edit: 2021/03/30 by Moritz Neuberger
          // Adding map of parent particles that create neutrons used above
          if((*secondaries)[i]->GetParticleDefinition()->GetParticleName() == "neutron")
          {
            fEventAction->neutronProducerMap.insert(
              std::make_pair((int) aTrack->GetTrackID(),
                             (int) aTrack->GetParticleDefinition()->GetPDGEncoding()));
          }
        }
      }
    }
  }

  // For Ge77m IC readout
  if(aTrack->GetParticleDefinition()->GetPDGEncoding() == 1000320771)
  {
    fEventAction->SetisMetastable(1);
    int NumberOfSecundaries = aTrack->GetStep()->GetSecondaryInCurrentStep()->size();
    for(int i = 0; i < NumberOfSecundaries; i++)
    {
      if(aTrack->GetStep()
             ->GetSecondaryInCurrentStep()
             ->at(i)
             ->GetParticleDefinition()
             ->GetParticleName() == "gamma" &&
         abs(aTrack->GetStep()->GetSecondaryInCurrentStep()->at(i)->GetTotalEnergy() /
               eV -
             160e3) < 1e3)
        fEventAction->SetisIC(1);
    }
  }

  // Edit: 2021/03/30 by Moritz Neuberger

  // Adding tracking of nC on different nuclei
  if(aTrack->GetParticleDefinition()->GetParticleName() == "neutron")
  {
    if(aTrack->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() ==
       "biasWrapper(nCapture)")  // altered name necessary due to biasing
    {
      int NumberOfSecundaries = aTrack->GetStep()->GetSecondaryInCurrentStep()->size();
      bool is_neutron_capture_on_Ge76 = false;
      for(int i = 0; i < NumberOfSecundaries; i++){
        if(aTrack->GetStep()
               ->GetSecondaryInCurrentStep()
               ->at(i)
               ->GetParticleDefinition()
               ->GetAtomicMass() == 77 &&
           aTrack->GetStep()
               ->GetSecondaryInCurrentStep()
               ->at(i)
               ->GetParticleDefinition()
               ->GetPDGCharge() == 32)
        {
          is_neutron_capture_on_Ge76 = true;
        }  // sending info of neutron producing Ge77 to output
      }
      if(is_neutron_capture_on_Ge76){
        double tmp_x, tmp_y, tmp_z;
        tmp_x = aTrack->GetStep()->GetPostStepPoint()->GetPosition().getX() / m;
        tmp_y = aTrack->GetStep()->GetPostStepPoint()->GetPosition().getY() / m;
        tmp_z = aTrack->GetStep()->GetPostStepPoint()->GetPosition().getZ() / m;
        fEventAction->AddEkin(aTrack->GetStep()->GetPreStepPoint()->GetKineticEnergy() /
                              eV);
        fEventAction->AddNeutronxLoc(tmp_neutronXpos);
        fEventAction->AddNeutronyLoc(tmp_neutronYpos);
        fEventAction->AddNeutronzLoc(tmp_neutronZpos);
        fEventAction->AddNeutronxMom(tmp_neutronXmom);
        fEventAction->AddNeutronyMom(tmp_neutronYmom);
        fEventAction->AddNeutronzMom(tmp_neutronZmom);
        fEventAction->AddNeutronTime(tmp_neutronTime);
        fEventAction->WriteMostOuterRadius();
        fEventAction->AddIDListOfGe77SiblingParticles(aTrack->GetTrackID());

        if(fRunAction->getPetersGammaCascadeModel()){
          KillSecundaries(aTrack);
          ReplaceNeutronCaptureGammas(aTrack);
        }

        for(int j = 0; j < aTrack->GetStep()->GetSecondaryInCurrentStep()->size(); j++){            
          fEventAction->AddNeutronCaptureSiblings_edep(aTrack->GetStep()->GetSecondaryInCurrentStep()->at(j)->GetKineticEnergy() /eV);
          fEventAction->AddNeutronCaptureSiblings_type(aTrack->GetStep()->GetSecondaryInCurrentStep()->at(j)->GetParticleDefinition()->GetPDGEncoding());
        }
      }
      
      for(int i = 0; i < NumberOfSecundaries; i++)
      {
        if(aTrack->GetStep()
            ->GetSecondaryInCurrentStep()
            ->at(i)
            ->GetParticleDefinition()
            ->GetPDGCharge() == 18)
        {
          fEventAction->AddnCAr_timing(
            aTrack->GetStep()->GetPostStepPoint()->GetGlobalTime() / s);
          fEventAction->AddnCAr_x(
            aTrack->GetStep()->GetPostStepPoint()->GetPosition().getX() / m);
          fEventAction->AddnCAr_y(
            aTrack->GetStep()->GetPostStepPoint()->GetPosition().getY() / m);
          fEventAction->AddnCAr_z(
            aTrack->GetStep()->GetPostStepPoint()->GetPosition().getZ() / m);
          fEventAction->AddnCAr_A(aTrack->GetStep()
                                    ->GetSecondaryInCurrentStep()
                                    ->at(i)
                                    ->GetParticleDefinition()
                                    ->GetAtomicMass());
        }  // send info of nC on Ar
        else
        {
          if(aTrack->GetStep()
              ->GetSecondaryInCurrentStep()
              ->at(i)
              ->GetParticleDefinition()
              ->GetPDGCharge() == 64)
          {
            fEventAction->AddnCGd_timing(
              aTrack->GetStep()->GetPostStepPoint()->GetGlobalTime() / s);
            fEventAction->AddnCGd_x(
              aTrack->GetStep()->GetPostStepPoint()->GetPosition().getX() / m);
            fEventAction->AddnCGd_y(
              aTrack->GetStep()->GetPostStepPoint()->GetPosition().getY() / m);
            fEventAction->AddnCGd_z(
              aTrack->GetStep()->GetPostStepPoint()->GetPosition().getZ() / m);
            fEventAction->AddnCGd_A(aTrack->GetStep()
                                      ->GetSecondaryInCurrentStep()
                                      ->at(i)
                                      ->GetParticleDefinition()
                                      ->GetAtomicMass());
            fEventAction->AddIDListOfGdSiblingParticles(aTrack->GetTrackID());
          }  // send info of nC on Gd
          else
          {
            fEventAction->AddnCOther_timing(
              aTrack->GetStep()->GetPostStepPoint()->GetGlobalTime() / s);
            fEventAction->AddnCOther_x(
              aTrack->GetStep()->GetPostStepPoint()->GetPosition().getX() / m);
            fEventAction->AddnCOther_y(
              aTrack->GetStep()->GetPostStepPoint()->GetPosition().getY() / m);
            fEventAction->AddnCOther_z(
              aTrack->GetStep()->GetPostStepPoint()->GetPosition().getZ() / m);
            fEventAction->AddnCOther_A(aTrack->GetStep()
                                        ->GetSecondaryInCurrentStep()
                                        ->at(i)
                                        ->GetParticleDefinition()
                                        ->GetAtomicMass());
            fEventAction->AddnCOther_ZC(aTrack->GetStep()
                                          ->GetSecondaryInCurrentStep()
                                          ->at(i)
                                          ->GetParticleDefinition()
                                          ->GetPDGCharge());
          }  // send info of all other nC
        }
      }
    }
  }
  if(fEventAction->isAllProductions()){
    int NumberOfSecundaries = aTrack->GetStep()->GetSecondaryInCurrentStep()->size();
    for(int i = 0; i < NumberOfSecundaries; i++)
    {
      if(aTrack->GetStep()
                ->GetSecondaryInCurrentStep()
                ->at(i)
                ->GetParticleDefinition()
                ->GetAtomicMass() > 0){
                    fEventAction->Addprod_timing(
                      aTrack->GetStep()->GetPostStepPoint()->GetGlobalTime() / s);
                    fEventAction->Addprod_x(
                      aTrack->GetStep()->GetPostStepPoint()->GetPosition().getX() / m);
                    fEventAction->Addprod_y(
                      aTrack->GetStep()->GetPostStepPoint()->GetPosition().getY() / m);
                    fEventAction->Addprod_z(
                      aTrack->GetStep()->GetPostStepPoint()->GetPosition().getZ() / m);
                    fEventAction->Addprod_A(aTrack->GetStep()
                                ->GetSecondaryInCurrentStep()
                                ->at(i)
                                ->GetParticleDefinition()
                                ->GetAtomicMass());
                    fEventAction->Addprod_ZC(aTrack->GetStep()
                                ->GetSecondaryInCurrentStep()
                                ->at(i)
                                ->GetParticleDefinition()
                                ->GetPDGCharge());
                    fEventAction->Addprod_parentType(aTrack->GetParticleDefinition()
                                ->GetPDGEncoding());
                }
    }
  }
}
