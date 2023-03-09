#include "WLGDEventAction.hh"
#include "WLGDTrajectory.hh"
#include "G4AnalysisManager.hh"

#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include "WLGDCrystalSD.hh"

#include "Randomize.hh"
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <vector>

WLGDEventAction::WLGDEventAction() { DefineCommands(); }

WLGDCrystalHitsCollection* WLGDEventAction::GetHitsCollection(G4int          hcID,
                                                              const G4Event* event) const
{
  auto hitsCollection =
    static_cast<WLGDCrystalHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));

  if(!hitsCollection)
  {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID;
    G4Exception("WLGDEventAction::GetHitsCollection()", "MyCode0001", FatalException,
                msg);
  }

  return hitsCollection;
}

void WLGDEventAction::makeMap()
{
  // make a lookup table to translate Logical Volume name
  // to a unique int for storage in the ntuple.
  lookup["Cavern_log"]   = 0;
  lookup["Hall_log"]     = 1;
  lookup["Tank_log"]     = 2;
  lookup["Water_log"]    = 3;
  lookup["Cout_log"]     = 4;
  lookup["Cvac_log"]     = 5;
  lookup["Cinn_log"]     = 6;
  lookup["Lar_log"]      = 7;
  lookup["Lid_log"]      = 8;
  lookup["Bot_log"]      = 9;
  lookup["Copper_log"]   = 10;
  lookup["ULar_log"]     = 11;
  lookup["Ge_log"]       = 12;
  lookup["Pu_log"]       = 13;
  lookup["Membrane_log"] = 14;
}

G4int WLGDEventAction::GeomID(G4String name)
{
  auto it = lookup.find(name);

  if(it == lookup.end())
  {
    G4ExceptionDescription msg;
    msg << "Name  " << name << " Inot in look up table";
    G4Exception("WLGDEventAction::GeomID()", "MyCode0005", FatalException, msg);
  }
  return it->second;
}

void WLGDEventAction::BeginOfEventAction(const G4Event*
                                         /*event*/)
{

  // -- clearing all the vectors so that every new event saves also the info of the previous ones and our disk do not explode

  nGe77.clear();
  edep.clear();
  ekin.clear();
  htrid.clear();
  thit.clear();
  whit.clear();
  xloc.clear();
  yloc.clear();
  zloc.clear();

  Muonxloc.clear();
  Muonyloc.clear();
  Muonzloc.clear();
  Muonxmom.clear();
  Muonymom.clear();
  Muonzmom.clear();
  Muonenergy.clear();

  neutronxloc.clear();
  neutronyloc.clear();
  neutronzloc.clear();
  neutronxmom.clear();
  neutronymom.clear();
  neutronzmom.clear();
  neutrontime.clear();

  tmp_mostOuterRadius = 0;
  out_mostOuterRadius.clear();

  // clear trajectory data
  trjpdg.clear();
  trjnpts.clear();
  nameid.clear();
  trjxvtx.clear();
  trjyvtx.clear();
  trjzvtx.clear();
  trjxpos.clear();
  trjypos.clear();
  trjzpos.clear();

  ReentranceTube.clear();
  DetectorNumber.clear();

  Multiplicity_prompt.clear();
  Multiplicity_prompt.push_back(0);
  Multiplicity_prompt.push_back(0);
  Multiplicity_prompt.push_back(0);
  Multiplicity_prompt.push_back(0);

  Multiplicity_delayed.clear();
  Multiplicity_delayed.push_back(0);
  Multiplicity_delayed.push_back(0);
  Multiplicity_delayed.push_back(0);
  Multiplicity_delayed.push_back(0);

  Multiplicity_delayed_long.clear();
  Multiplicity_delayed_long.push_back(0);
  Multiplicity_delayed_long.push_back(0);
  Multiplicity_delayed_long.push_back(0);
  Multiplicity_delayed_long.push_back(0);

  v_EdepPerDetector_prompt.clear();

  v_EdepPerDetector_delayed.clear();

  v_EdepPerDetector_delayed_long.clear();

  v_NDetector_prompt.clear();

  v_NDetector_delayed.clear();

  v_NDetector_delayed_long.clear();

  Multiplicity_prompt_woGd.clear();
  Multiplicity_prompt_woGd.push_back(0);
  Multiplicity_prompt_woGd.push_back(0);
  Multiplicity_prompt_woGd.push_back(0);
  Multiplicity_prompt_woGd.push_back(0);

  Multiplicity_delayed_woGd.clear();
  Multiplicity_delayed_woGd.push_back(0);
  Multiplicity_delayed_woGd.push_back(0);
  Multiplicity_delayed_woGd.push_back(0);
  Multiplicity_delayed_woGd.push_back(0);

  v_EdepPerDetector_prompt_woGd.clear();

  v_EdepPerDetector_delayed_woGd.clear();

  v_NDetector_prompt_woGd.clear();

  v_NDetector_delayed_woGd.clear();

  Multiplicity_prompt_onlyGd.clear();
  Multiplicity_prompt_onlyGd.push_back(0);
  Multiplicity_prompt_onlyGd.push_back(0);
  Multiplicity_prompt_onlyGd.push_back(0);
  Multiplicity_prompt_onlyGd.push_back(0);

  Multiplicity_delayed_onlyGd.clear();
  Multiplicity_delayed_onlyGd.push_back(0);
  Multiplicity_delayed_onlyGd.push_back(0);
  Multiplicity_delayed_onlyGd.push_back(0);
  Multiplicity_delayed_onlyGd.push_back(0);

  v_EdepPerDetector_prompt_onlyGd.clear();

  v_EdepPerDetector_delayed_onlyGd.clear();

  v_NDetector_prompt_onlyGd.clear();

  v_NDetector_delayed_onlyGd.clear();

  NumberOfNeutronsProducedInEvent.clear();
  NumberOfNeutronsProducedInEvent.push_back(0);

  v_EdepWater_prompt.clear();
  v_EdepWater_prompt.push_back(0);

  v_EdepWater_delayed.clear();
  v_EdepWater_delayed.push_back(0);

  v_MuonVeto_flag.clear();
  v_MuonVeto_flag.push_back(0);

  v_isIC.clear();
  v_isIC.push_back(0);

  v_isMetastable.clear();
  v_isMetastable.push_back(0);

  v_nCAr_timing.clear();
  v_nCAr_x.clear();
  v_nCAr_y.clear();
  v_nCAr_z.clear();
  v_nCAr_A.clear();

  v_nCGd_timing.clear();
  v_nCGd_x.clear();
  v_nCGd_y.clear();
  v_nCGd_z.clear();
  v_nCGd_A.clear();

  v_nCOther_timing.clear();
  v_nCOther_x.clear();
  v_nCOther_y.clear();
  v_nCOther_z.clear();
  v_nCOther_A.clear();
  v_nCOther_ZC.clear();

  v_prod_timing.clear();
  v_prod_x.clear();
  v_prod_y.clear();
  v_prod_z.clear();
  v_prod_A.clear();
  v_prod_ZC.clear();
  v_prod_parentType.clear();

  TotalEnergyDepositionInLAr_prompt.clear();
  TotalEnergyDepositionInLAr_prompt.push_back(0);
  TotalEnergyDepositionInLAr_prompt.push_back(0);
  TotalEnergyDepositionInLAr_prompt.push_back(0);
  TotalEnergyDepositionInLAr_prompt.push_back(0);

  TotalEnergyDepositionInGe_prompt.clear();
  TotalEnergyDepositionInGe_prompt.push_back(0);
  TotalEnergyDepositionInGe_prompt.push_back(0);
  TotalEnergyDepositionInGe_prompt.push_back(0);
  TotalEnergyDepositionInGe_prompt.push_back(0);

  TotalEnergyDepositionInLAr_delayed.clear();
  TotalEnergyDepositionInLAr_delayed.push_back(0);
  TotalEnergyDepositionInLAr_delayed.push_back(0);
  TotalEnergyDepositionInLAr_delayed.push_back(0);
  TotalEnergyDepositionInLAr_delayed.push_back(0);

  TotalEnergyDepositionInGe_delayed.clear();
  TotalEnergyDepositionInGe_delayed.push_back(0);
  TotalEnergyDepositionInGe_delayed.push_back(0);
  TotalEnergyDepositionInGe_delayed.push_back(0);
  TotalEnergyDepositionInGe_delayed.push_back(0);

  TotalEnergyDepositionInLAr_delayed_long.clear();
  TotalEnergyDepositionInLAr_delayed_long.push_back(0);
  TotalEnergyDepositionInLAr_delayed_long.push_back(0);
  TotalEnergyDepositionInLAr_delayed_long.push_back(0);
  TotalEnergyDepositionInLAr_delayed_long.push_back(0);

  TotalEnergyDepositionInGe_delayed_long.clear();
  TotalEnergyDepositionInGe_delayed_long.push_back(0);
  TotalEnergyDepositionInGe_delayed_long.push_back(0);
  TotalEnergyDepositionInGe_delayed_long.push_back(0);
  TotalEnergyDepositionInGe_delayed_long.push_back(0);

  TotalEnergyDepositionInLAr_after_delayed.clear();
  TotalEnergyDepositionInLAr_after_delayed.push_back(0);
  TotalEnergyDepositionInLAr_after_delayed.push_back(0);
  TotalEnergyDepositionInLAr_after_delayed.push_back(0);
  TotalEnergyDepositionInLAr_after_delayed.push_back(0);

  TotalEnergyDepositionInGe_after_delayed.clear();
  TotalEnergyDepositionInGe_after_delayed.push_back(0);
  TotalEnergyDepositionInGe_after_delayed.push_back(0);
  TotalEnergyDepositionInGe_after_delayed.push_back(0);
  TotalEnergyDepositionInGe_after_delayed.push_back(0);

  IndividualEnergyDeposition_Timing.clear();
  IndividualEnergyDeposition_LArOrGe.clear();
  IndividualEnergyDeposition_Energy.clear();
  IndividualEnergyDeposition_Position_x.clear();
  IndividualEnergyDeposition_Position_y.clear();
  IndividualEnergyDeposition_Position_z.clear();
  IndividualEnergyDeposition_ReentranceTube.clear();
  IndividualEnergyDeposition_ID.clear();
  IndividualEnergyDeposition_Type.clear();
  IndividualEnergyDeposition_DetectorNumber.clear();

  EdepPerDetector.clear();
  EdepPerDetector_delayed.clear();
  EdepPerDetector_delayed_long.clear();

  EdepPerDetector_prompt_woGd.clear();
  EdepPerDetector_delayed_woGd.clear();

  EdepPerDetector_prompt_onlyGd.clear();
  EdepPerDetector_delayed_onlyGd.clear();

  IDListOfGe77.clear();
  IDListOfGe77SiblingParticles.clear();
  IDListOfGdSiblingParticles.clear();

  v_Ge77Siblings_timing.clear();
  v_Ge77Siblings_x.clear();
  v_Ge77Siblings_y.clear();
  v_Ge77Siblings_z.clear();
  v_Ge77Siblings_edep.clear();
  v_Ge77Siblings_id.clear();
  v_Ge77Siblings_type.clear();
  v_Ge77Siblings_whichVolume.clear();

  v_NeutronCaptureSiblings_edep.clear();
  v_NeutronCaptureSiblings_type.clear();

  v_GdSiblings_timing.clear();
  v_GdSiblings_x.clear();
  v_GdSiblings_y.clear();
  v_GdSiblings_z.clear();
  v_GdSiblings_edep.clear();
  v_GdSiblings_id.clear();
  v_GdSiblings_type.clear();
  v_GdSiblings_whichVolume.clear();

  v_Ge77mGammaEmission_timing.clear();
  v_Ge77mGammaEmission_x.clear();
  v_Ge77mGammaEmission_y.clear();
  v_Ge77mGammaEmission_z.clear();
  v_Ge77mGammaEmission_edep.clear();
  v_Ge77mGammaEmission_id.clear();
  v_Ge77mGammaEmission_type.clear();
  v_Ge77mGammaEmission_whichVolume.clear();
  v_Ge77mGammaEmission_whichGe77.clear();

  v_Muon_WLSR_intersect_x.clear();
  v_Muon_WLSR_intersect_y.clear();
  v_Muon_WLSR_intersect_z.clear();
  v_Muon_WLSR_Edep.clear();
  v_Muon_WLSR_Edep.push_back(0);

  makeMap();
}

void WLGDEventAction::EndOfEventAction(const G4Event* event)
{
  // -- writing all the information out

  G4int eventID = event->GetEventID();
  if(((int) eventID) % 100 == 0)
    G4cout << ">>> Event: " << eventID << G4endl;

  // Get crystal hits collections IDs
  if(fHID < 0)
    fHID = G4SDManager::GetSDMpointer()->GetCollectionID("CrystalHitsCollection");

  // Get entries from hits collections
  //
  auto CrysHC = GetHitsCollection(fHID, event);

  nGe77.push_back(ekin.size());

  if(CrysHC->entries() <= 0 && fAllEvents == 0)
  {
    return;  // no action on no hit
  }

  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // fill Hits output from SD
  G4int nofHits = CrysHC->entries();

  for(G4int i = 0; i < nofHits; i++)
  {
    auto hh = (*CrysHC)[i];

    htrid.push_back(hh->GetTID());
    whit.push_back(hh->GetWeight());
    thit.push_back(hh->GetTime() / G4Analysis::GetUnitValue("ns"));
    edep.push_back(hh->GetEdep() / G4Analysis::GetUnitValue("MeV"));
    xloc.push_back((hh->GetPos()).x() / G4Analysis::GetUnitValue("m"));
    yloc.push_back((hh->GetPos()).y() / G4Analysis::GetUnitValue("m"));
    zloc.push_back((hh->GetPos()).z() / G4Analysis::GetUnitValue("m"));
    ReentranceTube.push_back(hh->GetWhichReentranceTube());
    DetectorNumber.push_back(hh->GetWhichDetector());
  }

  for(auto const& x : EdepPerDetector)
  {
    if(x.second < 1e4)
      continue;
    int tmp_i = (int) (x.first / 96);
    Multiplicity_prompt[tmp_i] += 1;
    v_NDetector_prompt.push_back(x.first);
    v_EdepPerDetector_prompt.push_back(x.second);
    IncreaseGeEnergyDeposition(x.second, tmp_i);
  }

  for(auto const& x : EdepPerDetector_delayed)
  {
    if(x.second < 1e4)
      continue;
    int tmp_i = (int) (x.first / 96);
    Multiplicity_delayed[tmp_i] += 1;
    v_NDetector_delayed.push_back(x.first);
    v_EdepPerDetector_delayed.push_back(x.second);
    IncreaseGeEnergyDeposition_delayed(x.second, tmp_i);
  }

  for(auto const& x : EdepPerDetector_delayed_long)
  {
    if(x.second < 1e4)
      continue;
    int tmp_i = (int) (x.first / 96);
    Multiplicity_delayed_long[tmp_i] += 1;
    v_NDetector_delayed_long.push_back(x.first);
    v_EdepPerDetector_delayed_long.push_back(x.second);
    IncreaseGeEnergyDeposition_delayed_long(x.second, tmp_i);
  }

  for(auto const& x : EdepPerDetector_prompt_woGd)
  {
    v_NDetector_prompt_woGd.push_back(x.first);
    v_EdepPerDetector_prompt_woGd.push_back(x.second);
  }
  for(auto const& x : EdepPerDetector_delayed_woGd)
  {
    v_NDetector_delayed_woGd.push_back(x.first);
    v_EdepPerDetector_delayed_woGd.push_back(x.second);
  }
  for(auto const& x : EdepPerDetector_prompt_onlyGd)
  {
    v_NDetector_prompt_onlyGd.push_back(x.first);
    v_EdepPerDetector_prompt_onlyGd.push_back(x.second);
  }
  for(auto const& x : EdepPerDetector_delayed_onlyGd)
  {
    v_NDetector_delayed_onlyGd.push_back(x.first);
    v_EdepPerDetector_delayed_onlyGd.push_back(x.second);
  }

  if(v_EdepWater_prompt[0] > 120e6)
    v_MuonVeto_flag[0] = 1;

  // fill trajectory data if available
  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  G4int                  n_trajectories =
    (trajectoryContainer == nullptr) ? 0 : trajectoryContainer->entries();

  if(n_trajectories > 0)
  {
    // temporary full storage
    std::vector<G4int>    temptid, temppid, temppdg, tempnpts;
    std::vector<G4String> tempname;
    std::vector<G4double> tempxvtx, tempyvtx, tempzvtx;
    std::vector<G4double> tempxpos, tempypos, tempzpos;

    for(G4int i = 0; i < n_trajectories; i++)
    {
      WLGDTrajectory* trj = (WLGDTrajectory*) ((*(event->GetTrajectoryContainer()))[i]);
      temptid.push_back(trj->GetTrackID());
      temppid.push_back(trj->GetParentID());
      temppdg.push_back(trj->GetPDGEncoding());
      tempname.push_back(trj->GetVertexName());
      tempxvtx.push_back((trj->GetVertex()).x());
      tempyvtx.push_back((trj->GetVertex()).y());
      tempzvtx.push_back((trj->GetVertex()).z());
      tempnpts.push_back(trj->GetPointEntries());
      for(int nn = 0; nn < trj->GetPointEntries(); ++nn)
      {
        tempxpos.push_back((trj->GetPoint(nn)->GetPosition()).x());
        tempypos.push_back((trj->GetPoint(nn)->GetPosition()).y());
        tempzpos.push_back((trj->GetPoint(nn)->GetPosition()).z());
      }
    }

    // store filtered trajectories only
    for(const int& item : htrid)
    {
      std::vector<int> res = FilterTrajectories(item, temptid, temppid);
      for(int& idx : res)
      {
        trjpdg.push_back(temppdg.at(idx));
        nameid.push_back(GeomID(tempname.at(idx)));
        trjxvtx.push_back(tempxvtx.at(idx));
        trjyvtx.push_back(tempyvtx.at(idx));
        trjzvtx.push_back(tempzvtx.at(idx));
        trjnpts.push_back(tempnpts.at(idx));
        int start = std::accumulate(tempnpts.begin(), tempnpts.begin() + idx, 0);
        for(int i = start; i < (start + tempnpts.at(idx)); ++i)
        {
          trjxpos.push_back(tempxpos.at(i));
          trjypos.push_back(tempypos.at(i));
          trjzpos.push_back(tempzpos.at(i));
        }
      }
    }
    temptid.clear();
    temppid.clear();
    temppdg.clear();
    tempnpts.clear();
    tempname.clear();
    tempxvtx.clear();
    tempyvtx.clear();
    tempzvtx.clear();
    tempxpos.clear();
    tempypos.clear();
    tempzpos.clear();
  }
  // fill the ntuple
  analysisManager->AddNtupleRow();

  // printing
  //  G4int eventID = event->GetEventID();
  /*  G4cout << ">>> Event: " << eventID << G4endl;
    G4cout << "    " << edep.size() << " hits stored in this event." << G4endl;
    G4cout << "    " << trjpdg.size() << " trajectories stored in this event." <<
    G4endl;*/
}

void WLGDEventAction::SaveAllEvents(G4int answer) { fAllEvents = answer; }
void WLGDEventAction::SaveAllProductions(G4int answer) { fAllProductions = answer; }

void WLGDEventAction::DefineCommands()
{
  // Define geometry command directory using generic messenger class
  fEventMessenger =
    new G4GenericMessenger(this, "/WLGD/event/", "Commands for controlling event action");

  // switch command
  fEventMessenger->DeclareMethod("saveAllEvents", &WLGDEventAction::SaveAllEvents)
    .SetGuidance("Set whether to save not only Ge77 but all events")
    .SetGuidance("0 = only Ge77 events")
    .SetGuidance("1 = all events")
    .SetCandidates("0 1")
    .SetDefaultValue("0");

  fEventMessenger->DeclareMethod("saveAllProductions", &WLGDEventAction::SaveAllProductions)
    .SetGuidance("Set whether all productions of nucleis are saved or not")
    .SetGuidance("0 = no additional output")
    .SetGuidance("1 = productions are saved")
    .SetCandidates("0 1")
    .SetDefaultValue("0");
}
