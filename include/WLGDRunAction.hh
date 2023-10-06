#ifndef WLGDRunAction_h
#define WLGDRunAction_h 1

#include "G4GenericMessenger.hh"
#include "G4UserRunAction.hh"
#include "globals.hh"

class WLGDEventAction;
class G4Run;

#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
/// Run action class
///

class WLGDRunAction : public G4UserRunAction
{
public:
  WLGDRunAction(WLGDEventAction* eventAction, G4String name);
  virtual ~WLGDRunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);
  // -- some high level output about the number of neutrons produced (can be used to identify whether you use the correct primary algorithm for example)
  void         increaseNumberOfCrossingNeutrons() { fNumberOfCrossingNeutrons++; }
  void         increaseTotalNumberOfNeutronsInLAr() { fTotalNumberOfNeutronsInLAr++; }
  G4int        getNumberOfNeutronsInLAr() { return fTotalNumberOfNeutronsInLAr; }

  // -- some redudant output for the production of neutrons inside the cryostat (not only Ge77 producing) should be removed in the future
  void         addEventNumber(G4int i)
  {
    vector_eventNumber.push_back(i);
  }  
  void addCoordinatsToFile(G4double x, G4double y, G4double z)
  {
    vector_x_dir.push_back(x);
    vector_y_dir.push_back(y);
    vector_z_dir.push_back(z);
  } 
  void addMomentumToFile(G4double x, G4double y, G4double z)
  {
    vector_x_mom.push_back(x);
    vector_y_mom.push_back(y);
    vector_z_mom.push_back(z);
  }  
  void addEnergyToFile(G4double ekin)
  {
    vector_energy.push_back(ekin);
  } 
  void addParentParticleType(G4int PDG_code)
  {
    vector_parentParticleType.push_back(PDG_code);
  }
  
  G4int getWriteOutNeutronProductionInfo() { return fWriteOutNeutronProductionInfo; }
  G4int getWriteOutGeneralNeutronInfo() { return fWriteOutGeneralNeutronInfo; }
  G4int getWriteOutAdvancedMultiplicity() { return fWriteOutAdvancedMultiplicity; }
  G4int getWriteOutAllNeutronInfoRoot() { return fWriteOutAllNeutronInfoRoot; }
  G4int getIndividualGeDepositionInfo() { return fIndividualGeDepositionInfo; }
  G4int getIndividualGdDepositionInfo() { return fIndividualGdDepositionInfo; }
  G4int getReadMuonCrossingWLSR()       { return fReadMuonCrossingWLSR; }
  G4int getNeutronCaptureSiblings()     { return fNeutronCaptureSiblings; }
  G4int getPetersGammaCascadeModel()    { return fPetersGammaCascadeModel; }

private:
  void DefineCommands();

private:
  G4GenericMessenger*   fMessenger;
  WLGDEventAction*      fEventAction;  // have event information for run
  G4String              fout;          // output file name
  G4int                 fNumberOfCrossingNeutrons;
  G4int                 fTotalNumberOfNeutronsInLAr;
  G4String              fout2, fout3;
  ofstream              outputStream_2;
  G4int                 fWriteOutNeutronProductionInfo = 0; 
  G4int                 fWriteOutNeutronProductionInfo_only_on_nc_on_Ge76 = 0; 
  G4int                 fWriteOutGeneralNeutronInfo    = 0;
  G4int                 fWriteOutAdvancedMultiplicity  = 0;
  G4int                 fWriteOutAllNeutronInfoRoot    = 0;
  G4int                 fIndividualGeDepositionInfo    = 0;
  G4int                 fIndividualGdDepositionInfo    = 0;
  G4int                 fReadMuonCrossingWLSR          = 0;
  G4int                 fNeutronCaptureSiblings        = 0;
  G4int                 fPetersGammaCascadeModel       = 0;
  std::vector<G4int>    vector_eventNumber;
  std::vector<G4double> vector_x_dir;
  std::vector<G4double> vector_y_dir;
  std::vector<G4double> vector_z_dir;
  std::vector<G4double> vector_x_mom;
  std::vector<G4double> vector_y_mom;
  std::vector<G4double> vector_z_mom;
  std::vector<G4double> vector_energy;
  std::vector<G4int>    vector_parentParticleType;
  void                  SetWriteOutNeutronProductionInfo(G4int answer);
  void                  SetWriteOutNeutronProductionInfo_only_on_nc_on_Ge76(G4int answer);
  void                  SetWriteOutGeneralNeutronInfo(G4int answer);
  void                  SetWriteOutAdvancedMultiplicity(G4int answer);
  void                  SetWriteOutAllNeutronInfoRoot(G4int answer);
  void                  SetIndividualGeDepositionInfo(G4int answer);
  void                  SetIndividualGdDepositionInfo(G4int answer);
  void                  SetReadMuonCrossingWLSR(G4int answer);
  void                  SetNeutronCaptureSiblings(G4int answer);
  void                  SetPetersGammaCascadeModel(G4int answer);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
