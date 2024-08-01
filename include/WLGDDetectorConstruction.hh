#ifndef WLGDDetectorConstruction_h
#define WLGDDetectorConstruction_h 1

#include "G4Cache.hh"
#include "G4GenericMessenger.hh"
#include "G4Material.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class WLGDCrystalSD;

class WLGDDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  WLGDDetectorConstruction();
  ~WLGDDetectorConstruction();

public:
  virtual G4VPhysicalVolume* Construct();
  virtual void               ConstructSDandField();
  G4String GetGeometryName() { return fGeometryName; }

  // -- general size of the experiment 
  G4double GetWorldSizeZ() { return fvertexZ; }  
  G4double GetWorldExtent() { return fmaxrad; }  

  // -- get geometry for turbine like structure (currently necessary for the particle generation inside the turbine structure)
  G4int    GetBoratedType() { return fWithBoratedPET; }
  G4double GetBoratedTurbineRadius() { return fBoratedTurbineRadius; }
  G4double GetBoratedTurbineLength() { return fBoratedTurbineLength; }
  G4double GetBoratedTurbineAngle() { return fBoratedTurbineAngle; }
  G4double GetBoratedTurbineWidth() { return fBoratedTurbineWidth; }
  G4double GetBoratedTurbineHeight() { return fBoratedTurbineHeight; }
  G4double GetBoratedTurbinezPosition() { return fBoratedTurbinezPosition; }
  G4int    GetBoratedTurbinezNPanels() { return fNPanels; }
  // -- a getter to know whether Gd water is used
  G4int isSetWithGdWater() { return fWithGdWater; }
  G4int isSetWithWoWater() { return fWithWoWater; }

  
  void  SetGeometry(const G4String& name);
  void  ExportGeometry(const G4String& file);
  
  // -- setter to adjust the cross-section biasing factors
  void  SetNeutronBiasFactor(G4double nf);
  void  SetMuonBiasFactor(G4double mf);
  void  SetNeutronYieldBias(G4double mf);

  // -- setters to adjust the geometry via macros

  // - w/wo cupper tubes (options: 0:[no], 1:yes)
  void  SetWithoutCupperTubes(G4int answer);

  // - select geometry of Neutron Moderators (options: 0: [no moderators], 1: around re-entrance tubes, 2: in turbine mode, 3: in large hollow tube mode)
  void  SetNeutronModerator(G4int answer);

  // - what material should the moderators be made of (options: [BoratedPE], PolyEthylene, PMMA)
  void  SetMaterial(G4String answer);

  // - w/wo Gd in the water (options: 0:[no], 1:yes)
  void  SetGdWater(G4int answer);
  void  SetWoWater(G4int answer);

  // - setters to adjust the size and radius of the turbines and tubes (I know confusingly named)
  void  SetTurbineAndTubeRadius(G4double radius);
  void  SetTurbineAndTubeLength(G4double length);
  void  SetTurbineAndTubeAngle(G4double deg);
  void  SetTurbineAndTubeWidth(G4double width);
  void  SetTurbineAndTubeHeight(G4double height);
  void  SetTurbineAndTubezPosition(G4double zPosition);
  void  SetTurbineAndTubeNPanels(G4double nPanels);

  // - set the concentration of Xe and He3 in the LAr
  void  SetXeConc(G4double nf);
  void  SetHe3Conc(G4double nf);
  void  SetMaGeMaterial(G4int answer);

  // - set the size of the cryostat
  void  SetOuterCryostatRadius(G4double rad);
  void  SetCryostatHeight(G4double height);
  
  // - adjust between original simulation position of detectors or updates ones
  void  SetPositionOfDetectors(G4String name);

private:
  void DefineCommands();
  void DefineMaterials();

  G4VPhysicalVolume* SetupBaseline();
  G4VPhysicalVolume* SetupAlternative();
  G4VPhysicalVolume* SetupHallA();

  G4GenericMessenger*     fDetectorMessenger       = nullptr;
  G4GenericMessenger*     fBiasMessenger           = nullptr;
  G4GenericMessenger*     fMaterialMessenger       = nullptr;
  G4double                fvertexZ                 = -1.0;
  G4double                fmaxrad                  = -1.0;
  G4String                fGeometryName            = "baseline";
  G4String                fDetectorPosition        = "baseline";
  G4double                fNeutronBias             = 1.0;
  G4double                fMuonBias                = 1.0;
  G4double                fNeutronYieldBias        = 1.0;
  G4Cache<WLGDCrystalSD*> fSD                      = nullptr;
  G4double                fXeConc                  = 0.0;
  G4double                fHe3Conc                 = 0.0;
  G4double                fCryostatOuterRadius     = 350.0;
  G4double                fCryostatHeight          = 350.0;
  G4double                fBoratedTurbineRadius    = 200.0;
  G4double                fBoratedTurbineLength    = 50.0;
  G4double                fBoratedTurbineAngle     = 45.0;
  G4double                fBoratedTurbineWidth     = 5.0;
  G4double                fBoratedTurbineHeight    = 600.;
  G4double                fBoratedTurbinezPosition = 0.;
  G4int                   fNPanels;
  G4int                   fBoratedTurbineNPanels = 0;
  G4int                   fWithOutCupperTubes    = 0;
  G4int                   fWithBoratedPET        = 0;
  G4String                fSetMaterial           = "";
  G4int                   fWithGdWater           = 0;
  G4int                   fWithWoWater           = 0;
  G4int                   fMaGeMaterial          = 0;
  G4Material*             CombinedArXeHe3;
  G4Material*             water;  //NOTE THAT THIS IS GADOLINIUM LOADED!!
  G4Material*             larMat;

  // For the Optical Parameters
  std::vector<G4double> fenergySmall;
  // To reduce computing time:
  // isolated volumes that do not have detectors inside should not be optical
  // Therefore create extra optical active water material to discern
  G4Material*             fOpticalWater;
};

#endif
