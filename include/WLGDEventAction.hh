#ifndef WLGDEventAction_h
#define WLGDEventAction_h 1

#include <algorithm>
#include <numeric>
#include <vector>

#include "WLGDCrystalHit.hh"

#include "G4GenericMessenger.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"
#include <map>
#include <set>

/// Event action class
///

class WLGDEventAction : public G4UserEventAction
{
public:
  WLGDEventAction();
  virtual ~WLGDEventAction() = default;

  virtual void BeginOfEventAction(const G4Event* event);
  virtual void EndOfEventAction(const G4Event* event);

  // ---------------------------------
  // -- lots of getters and setters of variabls explaind below
  // ---------------------------------

  // to create columns for Ntuple
  std::vector<G4int>&    GetNGe77() { return nGe77; }
  std::vector<G4int>&    GetEvntID() { return evntID; }
  std::vector<G4int>&    GetHitTID() { return htrid; }
  std::vector<G4double>& GetHitEdep() { return edep; }
  std::vector<G4double>& GetHitEkin() { return ekin; }
  std::vector<G4double>& GetHitTime() { return thit; }
  std::vector<G4double>& GetHitWeight() { return whit; }
  std::vector<G4double>& GetHitxLoc() { return xloc; }
  std::vector<G4double>& GetHityLoc() { return yloc; }
  std::vector<G4double>& GetHitzLoc() { return zloc; }

  std::vector<G4double>& GetMuonxLoc() { return Muonxloc; }
  std::vector<G4double>& GetMuonyLoc() { return Muonyloc; }
  std::vector<G4double>& GetMuonzLoc() { return Muonzloc; }
  std::vector<G4double>& GetMuonxMom() { return Muonxmom; }
  std::vector<G4double>& GetMuonyMom() { return Muonymom; }
  std::vector<G4double>& GetMuonzMom() { return Muonzmom; }
  std::vector<G4double>& GetMuonEnergy() { return Muonenergy; }

  std::vector<G4double>& GetNeutronxLoc() { return neutronxloc; }
  std::vector<G4double>& GetNeutronyLoc() { return neutronyloc; }
  std::vector<G4double>& GetNeutronzLoc() { return neutronzloc; }
  std::vector<G4double>& GetNeutronxMom() { return neutronxmom; }
  std::vector<G4double>& GetNeutronyMom() { return neutronymom; }
  std::vector<G4double>& GetNeutronzMom() { return neutronzmom; }
  std::vector<G4double>& GetNeutronTime() { return neutrontime; }
  std::vector<G4int>&    GetNumberOfNeutronsInEvent()
  {
    return NumberOfNeutronsProducedInEvent;
  }
  std::vector<G4double>& GetNeutronsMostOuterRadius() { return out_mostOuterRadius; }

  std::vector<G4double>& GetNeutronxTrack() { return out_neutron_track_x; }
  std::vector<G4double>& GetNeutronyTrack() { return out_neutron_track_y; }
  std::vector<G4double>& GetNeutronzTrack() { return out_neutron_track_z; }

  std::vector<G4double>& GetLArEnergyDeposition()
  {
    return TotalEnergyDepositionInLAr_prompt;
  }
  std::vector<G4double>& GetGeEnergyDeposition()
  {
    return TotalEnergyDepositionInGe_prompt;
  }
  std::vector<G4double>& GetLArEnergyDeposition_delayed()
  {
    return TotalEnergyDepositionInLAr_delayed;
  }
  std::vector<G4double>& GetGeEnergyDeposition_delayed()
  {
    return TotalEnergyDepositionInGe_delayed;
  }
  std::vector<G4double>& GetLArEnergyDeposition_delayed_long()
  {
    return TotalEnergyDepositionInLAr_delayed_long;
  }
  std::vector<G4double>& GetGeEnergyDeposition_delayed_long()
  {
    return TotalEnergyDepositionInGe_delayed_long;
  }
  std::vector<G4double>& GetLArEnergyDeposition_after_delayed()
  {
    return TotalEnergyDepositionInLAr_after_delayed;
  }
  std::vector<G4double>& GetGeEnergyDeposition_after_delayed()
  {
    return TotalEnergyDepositionInGe_after_delayed;
  }
  std::vector<G4int>& GetReentranceTube() { return ReentranceTube; }
  std::vector<G4int>& GetDetectorNumber() { return DetectorNumber; }

  std::vector<G4double>& GetIndividualEnergyDeposition_Timing()
  {
    return IndividualEnergyDeposition_Timing;
  }
  std::vector<G4double>& GetIndividualEnergyDeposition_Energy()
  {
    return IndividualEnergyDeposition_Energy;
  }
  std::vector<G4double>& GetIndividualEnergyDeposition_Position_x()
  {
    return IndividualEnergyDeposition_Position_x;
  }
  std::vector<G4double>& GetIndividualEnergyDeposition_Position_y()
  {
    return IndividualEnergyDeposition_Position_y;
  }
  std::vector<G4double>& GetIndividualEnergyDeposition_Position_z()
  {
    return IndividualEnergyDeposition_Position_z;
  }
  std::vector<G4int>& GetIndividualEnergyDeposition_ReentranceTube()
  {
    return IndividualEnergyDeposition_ReentranceTube;
  }
  std::vector<G4int>& GetIndividualEnergyDeposition_LArOrGe()
  {
    return IndividualEnergyDeposition_LArOrGe;
  }
  std::vector<G4int>& GetIndividualEnergyDeposition_ID()
  {
    return IndividualEnergyDeposition_ID;
  }
  std::vector<G4int>& GetIndividualEnergyDeposition_ParentID()
  {
    return IndividualEnergyDeposition_ParentID;
  }
  std::vector<G4int>& GetIndividualEnergyDeposition_VolumeNumber()
  {
    return IndividualEnergyDeposition_VolumeNumber;
  }
  std::vector<G4int>& GetIndividualEnergyDeposition_Type()
  {
    return IndividualEnergyDeposition_Type;
  }
  std::vector<G4int>& GetIndividualEnergyDeposition_DetectorNumber()
  {
    return IndividualEnergyDeposition_DetectorNumber;
  }

  std::vector<G4int>& GetMultiplicity_prompt() { return Multiplicity_prompt; }
  std::vector<G4int>& GetMultiplicity_delayed() { return Multiplicity_delayed; }
  std::vector<G4int>& GetMultiplicity_delayed_long() { return Multiplicity_delayed_long; }
  std::vector<G4double>& GetEdepPerDetector_prompt() { return v_EdepPerDetector_prompt; }
  std::vector<G4double>& GetEdepPerDetector_delayed()
  {
    return v_EdepPerDetector_delayed;
  }
  std::vector<G4double>& GetEdepPerDetector_delayed_long()
  {
    return v_EdepPerDetector_delayed_long;
  }
  std::vector<G4double>& GetNDetector_prompt() { return v_NDetector_prompt; }
  std::vector<G4double>& GetNDetector_delayed() { return v_NDetector_delayed; }
  std::vector<G4double>& GetNDetector_delayed_long() { return v_NDetector_delayed_long; }

  std::vector<G4int>& GetMultiplicity_prompt_woGd() { return Multiplicity_prompt_woGd; }
  std::vector<G4int>& GetMultiplicity_delayed_woGd() { return Multiplicity_delayed_woGd; }
  std::vector<G4double>& GetEdepPerDetector_prompt_woGd()
  {
    return v_EdepPerDetector_prompt_woGd;
  }
  std::vector<G4double>& GetEdepPerDetector_delayed_woGd()
  {
    return v_EdepPerDetector_delayed_woGd;
  }
  std::vector<G4double>& GetNDetector_prompt_woGd() { return v_NDetector_prompt_woGd; }
  std::vector<G4double>& GetNDetector_delayed_woGd() { return v_NDetector_delayed_woGd; }

  std::vector<G4int>& GetMultiplicity_prompt_onlyGd()
  {
    return Multiplicity_prompt_onlyGd;
  }
  std::vector<G4int>& GetMultiplicity_delayed_onlyGd()
  {
    return Multiplicity_delayed_onlyGd;
  }
  std::vector<G4double>& GetEdepPerDetector_prompt_onlyGd()
  {
    return v_EdepPerDetector_prompt_onlyGd;
  }
  std::vector<G4double>& GetEdepPerDetector_delayed_onlyGd()
  {
    return v_EdepPerDetector_delayed_onlyGd;
  }
  std::vector<G4double>& GetNDetector_prompt_onlyGd()
  {
    return v_NDetector_prompt_onlyGd;
  }
  std::vector<G4double>& GetNDetector_delayed_onlyGd()
  {
    return v_NDetector_delayed_onlyGd;
  }

  std::vector<G4double>& GetEdepWater_prompt() { return v_EdepWater_prompt; }
  std::vector<G4double>& GetEdepWater_delayed() { return v_EdepWater_delayed; }
  std::vector<G4int>&    GetMuonVeto_flag() { return v_MuonVeto_flag; }
  std::vector<G4int>&    GetisIC() { return v_isIC; }
  std::vector<G4int>&    GetisMetastable() { return v_isMetastable; }

  // tajectory methods
  std::vector<G4int>&    GetTrjPDG() { return trjpdg; }
  std::vector<G4int>&    GetTrjEntries() { return trjnpts; }
  std::vector<G4int>&    GetNameID() { return nameid; }
  std::vector<G4double>& GetTrjXVtx() { return trjxvtx; }
  std::vector<G4double>& GetTrjYVtx() { return trjyvtx; }
  std::vector<G4double>& GetTrjZVtx() { return trjzvtx; }
  std::vector<G4double>& GetTrjXPos() { return trjxpos; }
  std::vector<G4double>& GetTrjYPos() { return trjypos; }
  std::vector<G4double>& GetTrjZPos() { return trjzpos; }

  std::vector<G4double>& GetnCAr_timing() { return v_nCAr_timing; }
  std::vector<G4double>& GetnCAr_x() { return v_nCAr_x; }
  std::vector<G4double>& GetnCAr_y() { return v_nCAr_y; }
  std::vector<G4double>& GetnCAr_z() { return v_nCAr_z; }
  std::vector<G4int>&    GetnCAr_A() { return v_nCAr_A; }

  std::vector<G4double>& GetnCGd_timing() { return v_nCGd_timing; }
  std::vector<G4double>& GetnCGd_x() { return v_nCGd_x; }
  std::vector<G4double>& GetnCGd_y() { return v_nCGd_y; }
  std::vector<G4double>& GetnCGd_z() { return v_nCGd_z; }
  std::vector<G4int>&    GetnCGd_A() { return v_nCGd_A; }

  std::vector<G4double>& GetnCOther_timing() { return v_nCOther_timing; }
  std::vector<G4double>& GetnCOther_x() { return v_nCOther_x; }
  std::vector<G4double>& GetnCOther_y() { return v_nCOther_y; }
  std::vector<G4double>& GetnCOther_z() { return v_nCOther_z; }
  std::vector<G4int>&    GetnCOther_A() { return v_nCOther_A; }
  std::vector<G4int>&    GetnCOther_ZC() { return v_nCOther_ZC; }

  std::vector<G4double>& Getprod_timing() { return v_prod_timing; }
  std::vector<G4double>& Getprod_x() { return v_prod_x; }
  std::vector<G4double>& Getprod_y() { return v_prod_y; }
  std::vector<G4double>& Getprod_z() { return v_prod_z; }
  std::vector<G4int>&    Getprod_A() { return v_prod_A; }
  std::vector<G4int>&    Getprod_ZC() { return v_prod_ZC; }
  std::vector<G4int>&    Getprod_parentType() { return v_prod_parentType;}

  std::vector<G4double>& GetGe77mGammaEmission_timing()
  {
    return v_Ge77mGammaEmission_timing;
  }
  std::vector<G4double>& GetGe77mGammaEmission_x() { return v_Ge77mGammaEmission_x; }
  std::vector<G4double>& GetGe77mGammaEmission_y() { return v_Ge77mGammaEmission_y; }
  std::vector<G4double>& GetGe77mGammaEmission_z() { return v_Ge77mGammaEmission_z; }
  std::vector<G4double>& GetGe77mGammaEmission_edep()
  {
    return v_Ge77mGammaEmission_edep;
  }
  std::vector<G4int>& GetGe77mGammaEmission_id() { return v_Ge77mGammaEmission_id; }
  std::vector<G4int>& GetGe77mGammaEmission_type() { return v_Ge77mGammaEmission_type; }
  std::vector<G4int>& GetGe77mGammaEmission_whichVolume()
  {
    return v_Ge77mGammaEmission_whichVolume;
  }
  std::vector<G4int>& GetGe77mGammaEmission_whichGe77()
  {
    return v_Ge77mGammaEmission_whichGe77;
  }

  std::vector<G4double>& GetGe77Siblings_timing() { return v_Ge77Siblings_timing; }
  std::vector<G4double>& GetGe77Siblings_x() { return v_Ge77Siblings_x; }
  std::vector<G4double>& GetGe77Siblings_y() { return v_Ge77Siblings_y; }
  std::vector<G4double>& GetGe77Siblings_z() { return v_Ge77Siblings_z; }
  std::vector<G4double>& GetGe77Siblings_edep() { return v_Ge77Siblings_edep; }
  std::vector<G4int>&    GetGe77Siblings_id() { return v_Ge77Siblings_id; }
  std::vector<G4int>&    GetGe77Siblings_type() { return v_Ge77Siblings_type; }
  std::vector<G4int>& GetGe77Siblings_whichVolume() { return v_Ge77Siblings_whichVolume; }

  std::vector<G4double>& GetGdSiblings_timing() { return v_GdSiblings_timing; }
  std::vector<G4double>& GetGdSiblings_x() { return v_GdSiblings_x; }
  std::vector<G4double>& GetGdSiblings_y() { return v_GdSiblings_y; }
  std::vector<G4double>& GetGdSiblings_z() { return v_GdSiblings_z; }
  std::vector<G4double>& GetGdSiblings_edep() { return v_GdSiblings_edep; }
  std::vector<G4int>&    GetGdSiblings_id() { return v_GdSiblings_id; }
  std::vector<G4int>&    GetGdSiblings_type() { return v_GdSiblings_type; }
  std::vector<G4int>&    GetGdSiblings_whichVolume() { return v_GdSiblings_whichVolume; }


  std::vector<G4double>& GetNeutronCaptureSiblings_edep() { return v_NeutronCaptureSiblings_edep; }
  std::vector<G4int>&    GetNeutronCaptureSiblings_type() { return v_NeutronCaptureSiblings_type; }

  std::vector<G4double>& Get_Muon_WLSR_intersect_x(){ return v_Muon_WLSR_intersect_x; }
  std::vector<G4double>& Get_Muon_WLSR_intersect_y(){ return v_Muon_WLSR_intersect_y; }
  std::vector<G4double>& Get_Muon_WLSR_intersect_z(){ return v_Muon_WLSR_intersect_z; }
  std::vector<G4double>& Get_Muon_WLSR_Edep(){ return v_Muon_WLSR_Edep; }

  void Add_Muon_WLSR_intersect(G4double x, G4double y, G4double z){
    v_Muon_WLSR_intersect_x.push_back(x);
    v_Muon_WLSR_intersect_y.push_back(y);
    v_Muon_WLSR_intersect_z.push_back(z);
  }
  void Add_Muon_WLSR_Edep(double Edep){
    v_Muon_WLSR_Edep[0] += Edep;
  }

  void SetisIC(int ans) { v_isIC[0] = ans; }
  void SetisMetastable(int ans) { v_isMetastable[0] = ans; }

  void AddGdSiblings_timing(double timing) { v_GdSiblings_timing.push_back(timing); }
  void AddGdSiblings_x(double x) { v_GdSiblings_x.push_back(x); }
  void AddGdSiblings_y(double y) { v_GdSiblings_y.push_back(y); }
  void AddGdSiblings_z(double z) { v_GdSiblings_z.push_back(z); }
  void AddGdSiblings_edep(double edep) { v_GdSiblings_edep.push_back(edep); }
  void AddGdSiblings_id(int id) { v_GdSiblings_id.push_back(id); }
  void AddGdSiblings_type(int type) { v_GdSiblings_type.push_back(type); }
  void AddGdSiblings_whichVolume(int volume)
  {
    v_GdSiblings_whichVolume.push_back(volume);
  }

  void AddGe77Siblings_timing(double timing) { v_Ge77Siblings_timing.push_back(timing); }
  void AddGe77Siblings_x(double x) { v_Ge77Siblings_x.push_back(x); }
  void AddGe77Siblings_y(double y) { v_Ge77Siblings_y.push_back(y); }
  void AddGe77Siblings_z(double z) { v_Ge77Siblings_z.push_back(z); }
  void AddGe77Siblings_edep(double edep) { v_Ge77Siblings_edep.push_back(edep); }
  void AddGe77Siblings_id(int id) { v_Ge77Siblings_id.push_back(id); }
  void AddGe77Siblings_type(int type) { v_Ge77Siblings_type.push_back(type); }
  void AddGe77Siblings_whichVolume(int volume)
  {
    v_Ge77Siblings_whichVolume.push_back(volume);
  }

  void AddNeutronCaptureSiblings_edep(double edep) { v_NeutronCaptureSiblings_edep.push_back(edep); }
  void AddNeutronCaptureSiblings_type(int type) { v_NeutronCaptureSiblings_type.push_back(type); }

  void AddnCAr_timing(G4double time) { v_nCAr_timing.push_back(time); }
  void AddnCAr_x(G4double x) { v_nCAr_x.push_back(x); }
  void AddnCAr_y(G4double y) { v_nCAr_y.push_back(y); }
  void AddnCAr_z(G4double z) { v_nCAr_z.push_back(z); }
  void AddnCAr_A(G4int A) { v_nCAr_A.push_back(A); }

  void AddnCGd_timing(G4double time) { v_nCGd_timing.push_back(time); }
  void AddnCGd_x(G4double x) { v_nCGd_x.push_back(x); }
  void AddnCGd_y(G4double y) { v_nCGd_y.push_back(y); }
  void AddnCGd_z(G4double z) { v_nCGd_z.push_back(z); }
  void AddnCGd_A(G4int A) { v_nCGd_A.push_back(A); }

  void AddnCOther_timing(G4double time) { v_nCOther_timing.push_back(time); }
  void AddnCOther_x(G4double x) { v_nCOther_x.push_back(x); }
  void AddnCOther_y(G4double y) { v_nCOther_y.push_back(y); }
  void AddnCOther_z(G4double z) { v_nCOther_z.push_back(z); }
  void AddnCOther_A(G4int A) { v_nCOther_A.push_back(A); }
  void AddnCOther_ZC(G4int A) { v_nCOther_ZC.push_back(A); }

  void Addprod_timing(G4double time) { v_prod_timing.push_back(time); }
  void Addprod_x(G4double x) { v_prod_x.push_back(x); }
  void Addprod_y(G4double y) { v_prod_y.push_back(y); }
  void Addprod_z(G4double z) { v_prod_z.push_back(z); }
  void Addprod_A(G4int A) { v_prod_A.push_back(A); }
  void Addprod_ZC(G4int A) { v_prod_ZC.push_back(A); }
  void Addprod_parentType(G4int A){ v_prod_parentType.push_back(A); }


  void AddGe77mGammaEmission_timing(G4double time)
  {
    v_Ge77mGammaEmission_timing.push_back(time);
  }
  void AddGe77mGammaEmission_x(G4double x) { v_Ge77mGammaEmission_x.push_back(x); }
  void AddGe77mGammaEmission_y(G4double y) { v_Ge77mGammaEmission_y.push_back(y); }
  void AddGe77mGammaEmission_z(G4double z) { v_Ge77mGammaEmission_z.push_back(z); }
  void AddGe77mGammaEmission_edep(double edep)
  {
    v_Ge77mGammaEmission_edep.push_back(edep);
  }
  void AddGe77mGammaEmission_id(int id) { v_Ge77mGammaEmission_id.push_back(id); }
  void AddGe77mGammaEmission_type(int type) { v_Ge77mGammaEmission_type.push_back(type); }
  void AddGe77mGammaEmission_whichVolume(int volume)
  {
    v_Ge77mGammaEmission_whichVolume.push_back(volume);
  }
  void AddGe77mGammaEmission_whichGe77(int volume)
  {
    v_Ge77mGammaEmission_whichGe77.push_back(volume);
  }

  void AddMuonxLoc(G4double posx) { Muonxloc.push_back(posx); }
  void AddMuonyLoc(G4double posy) { Muonyloc.push_back(posy); }
  void AddMuonzLoc(G4double posz) { Muonzloc.push_back(posz); }
  void AddMuonxMom(G4double posx) { Muonxmom.push_back(posx); }
  void AddMuonyMom(G4double posy) { Muonymom.push_back(posy); }
  void AddMuonzMom(G4double posz) { Muonzmom.push_back(posz); }
  void AddMuonEnergy(G4double energy) { Muonenergy.push_back(energy); }

  void AddEkin(G4double kin) { ekin.push_back(kin); }
  void AddNeutronxLoc(G4double posx) { neutronxloc.push_back(posx); }
  void AddNeutronyLoc(G4double posy) { neutronyloc.push_back(posy); }
  void AddNeutronzLoc(G4double posz) { neutronzloc.push_back(posz); }
  void AddNeutronxMom(G4double posx) { neutronxmom.push_back(posx); }
  void AddNeutronyMom(G4double posy) { neutronymom.push_back(posy); }
  void AddNeutronzMom(G4double posz) { neutronzmom.push_back(posz); }
  void AddNeutronTime(G4double time) { neutrontime.push_back(time); }

  void AddIndividualEnergyDeposition_Timing(G4double time)
  {
    IndividualEnergyDeposition_Timing.push_back(time);
  }
  void AddIndividualEnergyDeposition_Energy(G4double energy)
  {
    IndividualEnergyDeposition_Energy.push_back(energy);
  }
  void AddIndividualEnergyDeposition_Position_x(G4double posx)
  {
    IndividualEnergyDeposition_Position_x.push_back(posx);
  }
  void AddIndividualEnergyDeposition_Position_y(G4double posy)
  {
    IndividualEnergyDeposition_Position_y.push_back(posy);
  }
  void AddIndividualEnergyDeposition_Position_z(G4double posz)
  {
    IndividualEnergyDeposition_Position_z.push_back(posz);
  }
  void AddIndividualEnergyDeposition_ReentranceTube(G4int n)
  {
    IndividualEnergyDeposition_ReentranceTube.push_back(n);
  }
  void AddIndividualEnergyDeposition_LArOrGe(G4int n)
  {
    IndividualEnergyDeposition_LArOrGe.push_back(n);
  }
  void AddIndividualEnergyDeposition_ID(G4int n)
  {
    IndividualEnergyDeposition_ID.push_back(n);
  }
  void AddIndividualEnergyDeposition_ParentID(G4int n)
  {
    IndividualEnergyDeposition_ParentID.push_back(n);
  }
  void AddIndividualEnergyDeposition_VolumeNumber(G4int n)
  {
    IndividualEnergyDeposition_VolumeNumber.push_back(n);
  }
  void AddIndividualEnergyDeposition_Type(G4int n)
  {
    IndividualEnergyDeposition_Type.push_back(n);
  }
  void AddIndividualEnergyDeposition_DetectorNumber(G4int n)
  {
    IndividualEnergyDeposition_DetectorNumber.push_back(n);
  }

  void            AddIDListOfGe77(G4int ID) { IDListOfGe77.insert(ID); }
  std::set<G4int> GetIDListOfGe77() { return IDListOfGe77; }

  void AddIDListOfGe77SiblingParticles(G4int ID)
  {
    IDListOfGe77SiblingParticles.insert(ID);
  }
  void RemoveIDListOfGe77SiblingParticles(G4int ID)
  {
    IDListOfGe77SiblingParticles.erase(ID);
  }
  std::set<G4int> GetIDListOfGe77SiblingParticles()
  {
    return IDListOfGe77SiblingParticles;
  }

  void AddIDListOfGdSiblingParticles(G4int ID) { IDListOfGdSiblingParticles.insert(ID); }
  void RemoveIDListOfGdSiblingParticles(G4int ID)
  {
    IDListOfGdSiblingParticles.erase(ID);
  }
  std::set<G4int> GetIDListOfGdSiblingParticles() { return IDListOfGdSiblingParticles; }

  std::map<int, int> neutronProducerMap;

  void IncreaseByOne_NeutronInEvent() { NumberOfNeutronsProducedInEvent[0] += 1; }
  void IncreaseLArEnergyDeposition(G4double Edep, G4int whichReEntranceTube)
  {
    TotalEnergyDepositionInLAr_prompt[whichReEntranceTube] += Edep;
  }  
  void IncreaseGeEnergyDeposition(G4double Edep, G4int whichReEntranceTube)
  {
    TotalEnergyDepositionInGe_prompt[whichReEntranceTube] += Edep;
  }  
  void IncreaseLArEnergyDeposition_delayed(G4double Edep, G4int whichReEntranceTube)
  {
    TotalEnergyDepositionInLAr_delayed[whichReEntranceTube] += Edep;
  }  
  void IncreaseGeEnergyDeposition_delayed(G4double Edep, G4int whichReEntranceTube)
  {
    TotalEnergyDepositionInGe_delayed[whichReEntranceTube] += Edep;
  }  
  void IncreaseLArEnergyDeposition_delayed_long(G4double Edep, G4int whichReEntranceTube)
  {
    TotalEnergyDepositionInLAr_delayed_long[whichReEntranceTube] += Edep;
  }  
  void IncreaseGeEnergyDeposition_delayed_long(G4double Edep, G4int whichReEntranceTube)
  {
    TotalEnergyDepositionInGe_delayed_long[whichReEntranceTube] += Edep;
  }  
  void IncreaseLArEnergyDeposition_after_delayed(G4double Edep, G4int whichReEntranceTube)
  {
    TotalEnergyDepositionInLAr_after_delayed[whichReEntranceTube] += Edep;
  }  
  void IncreaseGeEnergyDeposition_after_delayed(G4double Edep, G4int whichReEntranceTube)
  {
    TotalEnergyDepositionInGe_after_delayed[whichReEntranceTube] += Edep;
  }  
  void IncreaseEdepPerDetector(G4int copyNumber, G4double Edep)
  {
    EdepPerDetector[copyNumber] = EdepPerDetector[copyNumber] + Edep;
  }
  void IncreaseEdepPerDetector_delayed(G4int copyNumber, G4double Edep)
  {
    EdepPerDetector_delayed[copyNumber] = EdepPerDetector_delayed[copyNumber] + Edep;
  }
  void IncreaseEdepPerDetector_delayed_long(G4int copyNumber, G4double Edep)
  {
    EdepPerDetector_delayed_long[copyNumber] =
      EdepPerDetector_delayed_long[copyNumber] + Edep;
  }
  void IncreaseEdepWater_prompt(G4double Edep) { v_EdepWater_prompt[0] += Edep; }
  void IncreaseEdepWater_delayed(G4double Edep) { v_EdepWater_delayed[0] += Edep; }

  void IncreaseEdepPerDetector_prompt_woGd(G4int copyNumber, G4double Edep)
  {
    EdepPerDetector_prompt_woGd[copyNumber] =
      EdepPerDetector_prompt_woGd[copyNumber] + Edep;
  }
  void IncreaseEdepPerDetector_delayed_woGd(G4int copyNumber, G4double Edep)
  {
    EdepPerDetector_delayed_woGd[copyNumber] =
      EdepPerDetector_delayed_woGd[copyNumber] + Edep;
  }

  void IncreaseEdepPerDetector_prompt_onlyGd(G4int copyNumber, G4double Edep)
  {
    EdepPerDetector_prompt_onlyGd[copyNumber] =
      EdepPerDetector_prompt_onlyGd[copyNumber] + Edep;
  }
  void IncreaseEdepPerDetector_delayed_onlyGd(G4int copyNumber, G4double Edep)
  {
    EdepPerDetector_delayed_onlyGd[copyNumber] =
      EdepPerDetector_delayed_onlyGd[copyNumber] + Edep;
  }

  // Getters and setter for Neutron Capture
  // For Neutron tagger
  void Add_NCaptureIsotopes(G4int name)
  {
    NCaptureIsotopes.push_back(name);
  }

  void Add_NCaptureX(G4double pos)
  {
    NCaptureX.push_back(pos);
  }

  void Add_NCaptureY(G4double pos)
  {
    NCaptureY.push_back(pos);
  }

  void Add_NCaptureZ(G4double pos)
  {
    NCaptureZ.push_back(pos);
  }

  void Add_NCaptureEnergy(G4double energy)
  {
    NCaptureEnergy.push_back(energy);
  }

  void Add_NCaptureTime(G4double time)
  {
    NCaptureTime.push_back(time);
  }

  std::vector<G4int>&       GetNCaptureAmount()   
  {
    NCaptureAmount.push_back(NCaptureIsotopes.size());
    return NCaptureAmount; 
  }
  std::vector<G4int>&    GetNCaptureIsotopes() { return NCaptureIsotopes; }
  std::vector<G4double>&    GetNCaptureX() { return NCaptureX; }
  std::vector<G4double>&    GetNCaptureY() { return NCaptureY; }
  std::vector<G4double>&    GetNCaptureZ() { return NCaptureZ; }
  std::vector<G4double>&    GetNCaptureEnergy() { return NCaptureEnergy; }
  std::vector<G4double>&    GetNCaptureTime() { return NCaptureTime; }

  // The default is that only events with Ge77 production are saved. This function relaxes the condition so that all events are saved.
  void SaveAllEvents(G4int answer);
  void SaveAllProductions(G4int answer);

  void DefineCommands();

  void SaveNeutronTrack()
  {
    out_neutron_track_x = tmp_neutron_track_x;
    out_neutron_track_y = tmp_neutron_track_y;
    out_neutron_track_z = tmp_neutron_track_z;
  }

  void ClearNeutronTrack()
  {
    tmp_neutron_track_x.clear();
    tmp_neutron_track_y.clear();
    tmp_neutron_track_z.clear();
  }

  void PushNeutronTrack_x(G4double x) { tmp_neutron_track_x.push_back(x); }
  void PushNeutronTrack_y(G4double y) { tmp_neutron_track_y.push_back(y); }
  void PushNeutronTrack_z(G4double z) { tmp_neutron_track_z.push_back(z); }

  void SetMostOuterRadius(double r) { tmp_mostOuterRadius = r; }
  void UpdateMostOuterRadius(double r)
  {
    if(tmp_mostOuterRadius < r)
      tmp_mostOuterRadius = r;
  }
  void WriteMostOuterRadius() { out_mostOuterRadius.push_back(tmp_mostOuterRadius); }
  G4int isAllProductions(){return fAllProductions;};

private:
  // methods
  WLGDCrystalHitsCollection* GetHitsCollection(G4int hcID, const G4Event* event) const;
  G4int                      GeomID(G4String name);
  void                       makeMap();
  G4GenericMessenger*        fEventMessenger;
  G4int                      fAllEvents = 0;
  G4int                      fAllProductions = 0;

  //! Brief description
  /*!
   * This filter method retrieves the track history of a Germanium hit.
   * Starting with the hit id, the identical track id is found in
   * the trajectory container with vectors extracted separately,
   * here the track id and the corresponding track parent id. The parent
   * id then serves as the new id to look for, all the way back to the
   * primary particle. That way, only the relevant track history from
   * start to finish is stored.
   *
   * \param[in] item hit id to look for as starting point of history
   * \param[in] tid vector of all track id's stored in event
   * \param[in] pid corresponding vector of all parent id's for all tracks in event.
   */
  std::vector<int> FilterTrajectories(int item, const std::vector<G4int>& tid,
                                      const std::vector<G4int>& pid)
  {
    int              idx  = 0;
    int              pidx = 0;
    std::vector<int> result;
    auto             it = std::find(tid.begin(), tid.end(), item);

    while(it != tid.end())  // find all links in the chain
    {
      idx = (it - tid.begin());  // location of id
      result.push_back(idx);

      pidx = pid.at(idx);  // next to look for
      it   = std::find(tid.begin(), tid.end(), pidx);
    }

    return result;
  }

  // data members
  // hit data

  // -- a lot of output variables

  // - variables concerning the Ge77 production
  std::vector<G4int>    nGe77;
  std::vector<G4int>    evntID;
  G4int                 fHID = -1;
  std::vector<G4int>    htrid;
  std::vector<G4int>    ReentranceTube;
  std::vector<G4int>    DetectorNumber;
  std::vector<G4double> edep;
  std::vector<G4double> ekin;
  std::vector<G4double> thit;
  std::vector<G4double> whit;
  std::vector<G4double> xloc;
  std::vector<G4double> yloc;
  std::vector<G4double> zloc;

  // - production location and momentum at production of the captured neutron
  std::vector<G4double> neutronxloc;
  std::vector<G4double> neutronyloc;
  std::vector<G4double> neutronzloc;
  std::vector<G4double> neutronxmom;
  std::vector<G4double> neutronymom;
  std::vector<G4double> neutronzmom;
  std::vector<G4double> neutrontime;

  // - redundant code for previously reading out the neutron track
  std::vector<G4double> out_neutron_track_x;
  std::vector<G4double> out_neutron_track_y;
  std::vector<G4double> out_neutron_track_z;
  std::vector<G4double> tmp_neutron_track_x;
  std::vector<G4double> tmp_neutron_track_y;
  std::vector<G4double> tmp_neutron_track_z;
  std::vector<G4int>    NumberOfNeutronsProducedInEvent;

  // - information about the initial muon
  std::vector<G4double> Muonxloc;
  std::vector<G4double> Muonyloc;
  std::vector<G4double> Muonzloc;
  std::vector<G4double> Muonxmom;
  std::vector<G4double> Muonymom;
  std::vector<G4double> Muonzmom;
  std::vector<G4double> Muonenergy;

  // -- additional data for other produced particles
  // - production location, timing and mass of nuclei produced in neutron capture in Ar
  std::vector<G4double> v_nCAr_timing;
  std::vector<G4double> v_nCAr_x;
  std::vector<G4double> v_nCAr_y;
  std::vector<G4double> v_nCAr_z;
  std::vector<G4int>    v_nCAr_A;

  // - production location, timing and mass of nuclei produced in neutron capture in Gd
  std::vector<G4double> v_nCGd_timing;
  std::vector<G4double> v_nCGd_x;
  std::vector<G4double> v_nCGd_y;
  std::vector<G4double> v_nCGd_z;
  std::vector<G4int>    v_nCGd_A;

  // - production location, timing and mass of nuclei produced in neutron capture in other elements
  std::vector<G4double> v_nCOther_timing;
  std::vector<G4double> v_nCOther_x;
  std::vector<G4double> v_nCOther_y;
  std::vector<G4double> v_nCOther_z;
  std::vector<G4int>    v_nCOther_A;
  std::vector<G4int>    v_nCOther_ZC;

  // - production location, timing and mass of nuclei produced in a general matter
  std::vector<G4double> v_prod_timing;
  std::vector<G4double> v_prod_x;
  std::vector<G4double> v_prod_y;
  std::vector<G4double> v_prod_z;
  std::vector<G4int>    v_prod_A;
  std::vector<G4int>    v_prod_ZC;
  std::vector<G4int>    v_prod_parentType;

  // - output of the subsidiary particles after an IC of a Ge77m
  std::vector<G4double> v_Ge77mGammaEmission_timing;
  std::vector<G4double> v_Ge77mGammaEmission_x;
  std::vector<G4double> v_Ge77mGammaEmission_y;
  std::vector<G4double> v_Ge77mGammaEmission_z;
  std::vector<G4double> v_Ge77mGammaEmission_edep;
  std::vector<G4int>    v_Ge77mGammaEmission_id;
  std::vector<G4int>    v_Ge77mGammaEmission_type;
  std::vector<G4int>    v_Ge77mGammaEmission_whichVolume;
  std::vector<G4int>    v_Ge77mGammaEmission_whichGe77;

  // - output of the total energy deposited in the different time frames (prompt < 10µs, delayed > 10µs and <1mus, after delayed > 1ms)
  std::vector<G4double> TotalEnergyDepositionInLAr_prompt;
  std::vector<G4double> TotalEnergyDepositionInLAr_delayed;

  std::vector<G4double> TotalEnergyDepositionInLAr_delayed_long;
  std::vector<G4double> TotalEnergyDepositionInLAr_after_delayed;

  std::vector<G4double> TotalEnergyDepositionInGe_prompt;
  std::vector<G4double> TotalEnergyDepositionInGe_delayed;
  std::vector<G4double> TotalEnergyDepositionInGe_delayed_long;
  std::vector<G4double> TotalEnergyDepositionInGe_after_delayed;

  // - all energy depositions in the cryostat (should only be put out for small number of events (<10000))
  std::vector<G4double> IndividualEnergyDeposition_Timing;
  std::vector<G4double> IndividualEnergyDeposition_Energy;
  std::vector<G4double> IndividualEnergyDeposition_Position_x;
  std::vector<G4double> IndividualEnergyDeposition_Position_y;
  std::vector<G4double> IndividualEnergyDeposition_Position_z;
  std::vector<G4int>    IndividualEnergyDeposition_ReentranceTube;
  std::vector<G4int>    IndividualEnergyDeposition_LArOrGe;
  std::vector<G4int>    IndividualEnergyDeposition_ID;
  std::vector<G4int>    IndividualEnergyDeposition_ParentID;
  std::vector<G4int> IndividualEnergyDeposition_VolumeNumber;
  std::vector<G4int>    IndividualEnergyDeposition_Type;
  std::vector<G4int>    IndividualEnergyDeposition_DetectorNumber;

  // -- output for the Ge77m vetoing in the different time frames (prompt < 10µs, delayed > 10µs and <1mus, after delayed > 1ms)
  // - number of detectors triggered per re-entrance tube  
  std::vector<G4int>        Multiplicity_prompt;
  std::vector<G4int>        Multiplicity_delayed;
  std::vector<G4int>        Multiplicity_delayed_long;
  // - energy deposited per detector
  std::vector<G4double>     v_EdepPerDetector_prompt;
  std::vector<G4double>     v_EdepPerDetector_delayed;
  std::vector<G4double>     v_EdepPerDetector_delayed_long;
  // - the id numbers of the triggered detectors
  std::vector<G4double>     v_NDetector_prompt;
  std::vector<G4double>     v_NDetector_delayed;
  std::vector<G4double>     v_NDetector_delayed_long;
  // - now at consideration a bit redudant, but it seems like it is a map with the id of the detector pointing at the energy deposited in it
  std::map<G4int, G4double> EdepPerDetector;
  std::map<G4int, G4double> EdepPerDetector_delayed;
  std::map<G4int, G4double> EdepPerDetector_delayed_long;

  // -- some additional variables, mainly for small investigations, but generally redundant and should not be used
  std::vector<G4int>        Multiplicity_prompt_woGd;
  std::vector<G4int>        Multiplicity_delayed_woGd;
  std::vector<G4double>     v_EdepPerDetector_prompt_woGd;
  std::vector<G4double>     v_EdepPerDetector_delayed_woGd;
  std::vector<G4double>     v_NDetector_prompt_woGd;
  std::vector<G4double>     v_NDetector_delayed_woGd;
  std::map<G4int, G4double> EdepPerDetector_prompt_woGd;
  std::map<G4int, G4double> EdepPerDetector_delayed_woGd;
  std::vector<G4int>        Multiplicity_prompt_onlyGd;
  std::vector<G4int>        Multiplicity_delayed_onlyGd;
  std::vector<G4double>     v_EdepPerDetector_prompt_onlyGd;
  std::vector<G4double>     v_EdepPerDetector_delayed_onlyGd;
  std::vector<G4double>     v_NDetector_prompt_onlyGd;
  std::vector<G4double>     v_NDetector_delayed_onlyGd;
  std::map<G4int, G4double> EdepPerDetector_prompt_onlyGd;
  std::map<G4int, G4double> EdepPerDetector_delayed_onlyGd;

  // -- some general information
  // - energy deposited in the water tank (prompt, delayed)
  std::vector<G4double> v_EdepWater_prompt;
  std::vector<G4double> v_EdepWater_delayed;
  // - if prompt more than 120GeV deposited in the water, raise the muon veto flag
  std::vector<G4int>    v_MuonVeto_flag;
  // - if the Ge77 performs an internal conversion, raise flag
  std::vector<G4int>    v_isIC;
  // - if the Ge77 is produced in the meta-stable state, raise flag (one should consider that with muon primaries this will never be the case)
  std::vector<G4int>    v_isMetastable;



  // -- the particle id list of the particles of interest 
  std::set<G4int> IDListOfGe77;
  std::set<G4int> IDListOfGdSiblingParticles;
  std::set<G4int> IDListOfGe77SiblingParticles;

  // -- information about the sibling particles of the neutron capture on Ge76
  std::vector<G4double> v_Ge77Siblings_timing;
  std::vector<G4double> v_Ge77Siblings_x;
  std::vector<G4double> v_Ge77Siblings_y;
  std::vector<G4double> v_Ge77Siblings_z;
  std::vector<G4double> v_Ge77Siblings_edep;
  std::vector<G4int>    v_Ge77Siblings_id;
  std::vector<G4int>    v_Ge77Siblings_type;
  std::vector<G4int>    v_Ge77Siblings_whichVolume;

  // -- information about the sibling particles of the neutron capture on Gd
  std::vector<G4double> v_GdSiblings_timing;
  std::vector<G4double> v_GdSiblings_x;
  std::vector<G4double> v_GdSiblings_y;
  std::vector<G4double> v_GdSiblings_z;
  std::vector<G4double> v_GdSiblings_edep;
  std::vector<G4int>    v_GdSiblings_id;
  std::vector<G4int>    v_GdSiblings_type;
  std::vector<G4int>    v_GdSiblings_whichVolume;
  
  std::vector<G4double> v_NeutronCaptureSiblings_edep;
  std::vector<G4int>    v_NeutronCaptureSiblings_type;

  // trajectory data
  std::vector<G4int>        trjpdg;
  std::vector<G4int>        trjnpts;
  std::vector<G4int>        nameid;
  std::vector<G4double>     trjxvtx;
  std::vector<G4double>     trjyvtx;
  std::vector<G4double>     trjzvtx;
  std::vector<G4double>     trjxpos;
  std::vector<G4double>     trjypos;
  std::vector<G4double>     trjzpos;
  std::map<G4String, G4int> lookup;

  G4double              tmp_mostOuterRadius;
  std::vector<G4double> out_mostOuterRadius;


  std::vector<G4double> v_Muon_WLSR_intersect_x;
  std::vector<G4double> v_Muon_WLSR_intersect_y;
  std::vector<G4double> v_Muon_WLSR_intersect_z;
  std::vector<G4double> v_Muon_WLSR_Edep;

  // EDIT: by Eric Esch
  // Neutron Capture for Neutrontagger
  std::vector<G4int>    NCaptureAmount;
  std::vector<G4int> NCaptureIsotopes;
  std::vector<G4double> NCaptureX;
  std::vector<G4double> NCaptureY;
  std::vector<G4double> NCaptureZ;
  std::vector<G4double> NCaptureEnergy;
  std::vector<G4double> NCaptureTime;
};

#endif
