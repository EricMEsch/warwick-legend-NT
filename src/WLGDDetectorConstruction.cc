#include "WLGDDetectorConstruction.hh"

#include <cmath>
#include <set>

#include "G4RunManager.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4GDMLParser.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4Tubs.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4SDManager.hh"
#include "WLGDCrystalSD.hh"

#include "WLGDBiasMultiParticleChangeCrossSection.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"


WLGDDetectorConstruction::WLGDDetectorConstruction()
{
  DefineCommands();
  DefineMaterials();
}

WLGDDetectorConstruction::~WLGDDetectorConstruction()
{
  delete fDetectorMessenger;
  delete fBiasMessenger;
}

auto WLGDDetectorConstruction::Construct() -> G4VPhysicalVolume*
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  if(fGeometryName == "baseline" || fGeometryName == "baseline_smaller" || fGeometryName == "baseline_large_reentrance_tube" || fGeometryName == "baseline_large_reentrance_tube_4m_cryo")
  {
    return SetupBaseline();
  }
  else if(fGeometryName == "hallA" || fGeometryName == "hallA_wo_ge" || fGeometryName == "hallA_only_WLSR")
  {
    return SetupHallA();
  }

  return SetupAlternative();
}

void WLGDDetectorConstruction::DefineMaterials()
{
  G4NistManager* nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Galactic");
  nistManager->FindOrBuildMaterial("G4_lAr");
  nistManager->FindOrBuildMaterial("G4_AIR");
  nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  nistManager->FindOrBuildMaterial("G4_Cu");
  nistManager->FindOrBuildMaterial("G4_WATER");

  auto* C  = new G4Element("Carbon", "C", 6., 12.011 * g / mole);
  auto* O  = new G4Element("Oxygen", "O", 8., 16.00 * g / mole);
  auto* Ca = new G4Element("Calcium", "Ca", 20., 40.08 * g / mole);
  auto* Mg = new G4Element("Magnesium", "Mg", 12., 24.31 * g / mole);

  // Standard Rock definition, similar to Gran Sasso rock
  // with density from PDG report
  auto* stdRock = new G4Material("StdRock", 2.65 * g / cm3, 4);
  stdRock->AddElement(O, 52.0 * perCent);
  stdRock->AddElement(Ca, 27.0 * perCent);
  stdRock->AddElement(C, 12.0 * perCent);
  stdRock->AddElement(Mg, 9.0 * perCent);

  auto* H     = new G4Element("Hydrogen", "H", 1., 1.00794 * g / mole);
  auto* N     = new G4Element("Nitrogen", "N", 7., 14.00 * g / mole);
  auto* puMat = new G4Material("polyurethane", 0.3 * g / cm3, 4);  // high density foam
  puMat->AddElement(H, 16);
  puMat->AddElement(O, 2);
  puMat->AddElement(C, 8);
  puMat->AddElement(N, 2);

  auto*    B10           = new G4Isotope("B10", 5, 10, 10 * g / mole);
  auto*    B11           = new G4Isotope("B11", 5, 11, 11 * g / mole);
  auto*    SpecialB      = new G4Element("SpecialB", "SpecialB", 2);
  G4double B_MassRatio   = 0.2;
  G4double B_NumberRatio = 1 / (1 + (1 - B_MassRatio) / B_MassRatio * 10. / 11.);
  SpecialB->AddIsotope(B10, B_NumberRatio);
  SpecialB->AddIsotope(B11, 1 - B_NumberRatio);

  // Borated Polyethylene according to GEM TN-92-172 TART Calculations of Neutron
  // Attenuation and Neutron-induced Photons on 5 % and 20 % Borated Polyethylene Slabs
  auto* BoratedPET =
    new G4Material("BoratedPET", 0.95 * g / cm3, 4);  // high density foam
  BoratedPET->AddElement(H, 0.116);
  BoratedPET->AddElement(C, 0.612);
  BoratedPET->AddElement(SpecialB, 0.05);
  BoratedPET->AddElement(O, 0.222);

  // Estimated using the number of elements per chain elements  (C_5 O_2 H_8)
  auto* PMMA = new G4Material("PMMA", 1.18 * g / cm3, 3);
  PMMA->AddElement(H, 0.08);
  PMMA->AddElement(C, 0.60);
  PMMA->AddElement(O, 0.32);

  // Estimated using the number of elements per molecule (C_2 H_4)
  auto* PolyEthylene = new G4Material("PolyEthylene", 0.95 * g / cm3, 2);
  PolyEthylene->AddElement(H, 0.142);
  PolyEthylene->AddElement(C, 0.857);

  G4Element* elGd = new G4Element("Gadolinium", "Gd", 64, 157.25 * g / mole);
  G4Element* elS  = new G4Element("Sulfur", "S", 16., 32.066 * g / mole);
  G4cout << elGd << G4endl;

  G4double density = 3.01 * g / cm3;  // https://www.sigmaaldrich.com/catalog/product/aldrich/203300?lang=de&region=DE
                                      // @room temp
  G4Material* gadoliniumSulfate =
    new G4Material("GadoliniumSulfate", density, 3);  // Gd2(SO4)3
  gadoliniumSulfate->AddElement(elGd, 2);
  gadoliniumSulfate->AddElement(elS, 3);
  gadoliniumSulfate->AddElement(O, 12);

  G4Material* purewater = G4Material::GetMaterial(
    "G4_WATER");  // EDIT: changed water to purewater & use it to create "special" water
  water = new G4Material("GdLoadedWater", 1.000000 * g / cm3, 2);
  water->AddMaterial(purewater, 1. - 0.002);
  water->AddMaterial(gadoliniumSulfate, 0.002);

  G4Material* gadoliniumOxide =
  new G4Material("GadoliniumOxide", 7.407 * g / cm3, 2);  // Gd2O3
  gadoliniumOxide->AddElement(elGd, 2);
  gadoliniumOxide->AddElement(O, 3);

  // For Water Neutron tagger
  G4Material* GdPMMA = new G4Material("GdLoadedPMMA", 1.18 * g / cm3, 2);
  GdPMMA->AddMaterial(PMMA, 1. - 0.002);
  GdPMMA->AddMaterial(gadoliniumOxide, 0.002);

  // enriched Germanium from isotopes
  auto* Ge_74 = new G4Isotope("Ge74", 32, 74, 74.0 * g / mole);
  auto* Ge_76 = new G4Isotope("Ge76", 32, 76, 76.0 * g / mole);

  G4Isotope* Ge_70 = new G4Isotope("Ge70",  32, 70, 69.92*g/mole);
  G4Isotope* Ge_72 = new G4Isotope("Ge72",  32, 72, 71.92*g/mole);
  G4Isotope* Ge_73 = new G4Isotope("Ge73",  32, 73, 73.0*g/mole);

  G4Element* eGe = new G4Element("enriched Germanium", "enrGe", 2);

  eGe->AddIsotope(Ge_76, 87. * perCent);
  eGe->AddIsotope(Ge_74, 13. * perCent);
  density      = 5.323 * g / cm3;

  if(fMaGeMaterial){
    density = 5.56 * g / cm3;
    eGe = new G4Element("enriched Germanium", "enrGe", 5);  
    eGe->AddIsotope(Ge_70,0.0*perCent);
    eGe->AddIsotope(Ge_72,0.1*perCent);
    eGe->AddIsotope(Ge_73,0.2*perCent);
    eGe->AddIsotope(Ge_74,13.1*perCent);
    eGe->AddIsotope(Ge_76,86.6*perCent);
  }

  auto* roiMat = new G4Material("enrGe", density, 1);
  roiMat->AddElement(eGe, 1);

  // Edit: 2020/02/17 by Moritz Neuberger
  // Added new def. of LAr in order to add doping with Xe and He-3

  G4double dLAr  = 1.393 * g / cm3;
  G4double dLXe  = 3.02 * g / cm3;
  G4double dLHe3 = 0.059 * g / cm3;

  G4double fArConc = 1 - (fXeConc + fHe3Conc);

  G4double dComb = 1 / ((fArConc / dLAr) + (fXeConc / dLXe) + (fHe3Conc / dLHe3));

  G4cout << "___________________________________________" << G4endl;
  G4cout << "Mass ratios of cryostat:" << G4endl;
  G4cout << "LAr:   " << fArConc << G4endl;
  G4cout << "LXe:   " << fXeConc << G4endl;
  G4cout << "LHe3:   " << fHe3Conc << G4endl;
  G4cout << "dComb:   " << dComb << G4endl;
  G4cout << "___________________________________________" << G4endl;

  
  //auto* eLAr = new G4Element("LAr", "Ar", 18., 39.95 * g / mole);
  larMat = G4Material::GetMaterial("G4_lAr");
  auto*      eLXe   = new G4Element("LXe", "Xe", 54., 131.29 * g / mole);
  auto*      eHe3   = new G4Element("He3", "He3", 1);
  G4Isotope* iHe3   = new G4Isotope("He3", 2, 3);
  eHe3->AddIsotope(iHe3, 1);

  if(fMaGeMaterial){
    double densityLAr = 1.396  * g / cm3; 
    G4Element* larEl = new G4Element("LAr", "Ar", 18., 39.95 * g / mole);
    larMat = new G4Material("LAr", densityLAr, 1);
    larMat->AddElement(larEl, 1);
  }

  CombinedArXeHe3 =
    new G4Material("CombinedArXeHe3", dComb, 3, kStateLiquid, 87. * kelvin);
  CombinedArXeHe3->AddMaterial(larMat, fArConc);
  CombinedArXeHe3->AddElement(eHe3, fHe3Conc);
  CombinedArXeHe3->AddElement(eLXe, fXeConc);

  G4cout << CombinedArXeHe3 << G4endl;

  // Optical Parameters:

  fenergySmall = { 1.239841939*eV/0.6,1.239841939*eV/0.1 };

  std::vector<G4double> energyAbs = {
  1.239841939*eV/0.6, 1.239841939*eV/0.55, 1.239841939*eV/0.50, 1.239841939*eV/0.45,  1.239841939*eV/0.40,  
  1.239841939*eV/0.35, 1.239841939*eV/0.30, 1.239841939*eV/0.25, 1.239841939*eV/0.20, 1.239841939*eV/0.19, 1.239841939*eV/0.1
  }; 

  std::vector<G4double> rindexWater = {1.33,1.33};

  std::vector<G4double> absH2O =  {
  10*m, 20*m, 50*m, 100*m,  100*m,  
  100*m, 90*m, 20*m, 1*m, 0.001*mm, 0.0001*mm // See https://www.researchgate.net/publication/307856024_Ultraviolet_250-550_nm_absorption_spectrum_of_pure_water
  };

  fOpticalWater = new G4Material("WaterWithOptical", 1.0*g/cm3, 2);
  fOpticalWater->AddElement(nistManager->FindOrBuildElement("H"), 2);
  fOpticalWater->AddElement(nistManager->FindOrBuildElement("O"), 1);

  G4MaterialPropertiesTable *mptH2O = new G4MaterialPropertiesTable();
  mptH2O->AddProperty("RINDEX", fenergySmall, rindexWater);
  mptH2O->AddProperty("ABSLENGTH", energyAbs, absH2O);
  fOpticalWater->SetMaterialPropertiesTable(mptH2O);

}

void WLGDDetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // Only need to construct the (per-thread) SD once
  if(!fSD.Get())
  {
    G4String       crystalSDname = "CrystalSD";
    WLGDCrystalSD* aCrystalSD =
      new WLGDCrystalSD(crystalSDname, "CrystalHitsCollection", fGeometryName);
    fSD.Put(aCrystalSD);

    // Also only add it once to the SD manager!
    G4SDManager::GetSDMpointer()->AddNewDetector(fSD.Get());

    SetSensitiveDetector("Ge_log", fSD.Get());

    // ----------------------------------------------
    // -- operator creation and attachment to volume:
    // ----------------------------------------------
    G4LogicalVolumeStore* volumeStore = G4LogicalVolumeStore::GetInstance();
    /*
        // -- Attach neutron XS biasing to Germanium -> enhance nCapture
        auto *biasnXS = new WLGDBiasMultiParticleChangeCrossSection();
        biasnXS->SetNeutronFactor(fNeutronBias);
        biasnXS->SetMuonFactor(fMuonBias);
        biasnXS->SetNeutronYieldFactor(fNeutronYieldBias);
        G4cout << " >>> Detector: set neutron bias to " << fNeutronBias << G4endl;
        biasnXS->AddParticle("neutron");
        G4LogicalVolume *logicGe = volumeStore->GetVolume("Ge_log");
        biasnXS->AttachTo(logicGe);
        G4LogicalVolume *logicLar = volumeStore->GetVolume("Lar_log");
        biasnXS->AttachTo(logicLar);
        G4LogicalVolume *logicTank = volumeStore->GetVolume("Tank_log");
        biasnXS->AttachTo(logicTank);
        G4LogicalVolume *logicCavern = volumeStore->GetVolume("Cavern_log");
        biasnXS->AttachTo(logicCavern);
        G4LogicalVolume *logicHall = volumeStore->GetVolume("Hall_log");
        biasnXS->AttachTo(logicHall);
    */
    // -- Attach neutron XS biasing to Germanium -> enhance nCapture
    auto* biasnXS = new WLGDBiasMultiParticleChangeCrossSection();
    biasnXS->SetNeutronFactor(fNeutronBias);
    biasnXS->SetMuonFactor(fMuonBias);
    G4cout << " >>> Detector: set neutron bias to " << fNeutronBias << G4endl;
    biasnXS->AddParticle("neutron");
    G4LogicalVolume* logicGe = volumeStore->GetVolume("Ge_log");
    biasnXS->AttachTo(logicGe);

    // -- Attach muon XS biasing to all required volumes consistently
    auto* biasmuXS = new WLGDBiasMultiParticleChangeCrossSection();
    biasmuXS->SetNeutronFactor(fNeutronBias);
    biasmuXS->SetMuonFactor(fMuonBias);
    biasmuXS->SetNeutronYieldFactor(fNeutronYieldBias);
    G4cout << " >>> Detector: set muon bias to " << fMuonBias << G4endl;
    biasmuXS->AddParticle("mu-");

    if(fNeutronYieldBias != 1)
    {
      biasmuXS->AddParticle("neutron");
      biasmuXS->AddParticle("pi+");
      biasmuXS->AddParticle("pi-");
      biasmuXS->AddParticle("gamma");
      biasmuXS->AddParticle("kaon-");
      biasmuXS->AddParticle("proton");
    }

    // G4LogicalVolume* logicGe = volumeStore->GetVolume("Ge_log");
    // biasmuXS->AttachTo(logicGe);
    G4LogicalVolume* logicCavern = volumeStore->GetVolume("Cavern_log");
    biasmuXS->AttachTo(logicCavern);
    G4LogicalVolume* logicHall = volumeStore->GetVolume("Hall_log");
    biasmuXS->AttachTo(logicHall);
    G4LogicalVolume* logicTank = volumeStore->GetVolume("Tank_log");
    biasmuXS->AttachTo(logicTank);
    G4LogicalVolume* logicLar = volumeStore->GetVolume("Lar_log");
    biasmuXS->AttachTo(logicLar);

    // non hallA have these volumes
    if(fGeometryName != "hallA")
    {
      G4LogicalVolume* logicCu = volumeStore->GetVolume("Copper_log");
      biasmuXS->AttachTo(logicCu);
      G4LogicalVolume* logicULar = volumeStore->GetVolume("ULar_log");
      biasmuXS->AttachTo(logicULar);
    }

    // Baseline also has a water volume and cryostat
    if(fGeometryName == "baseline" || fGeometryName == "hallA" || fGeometryName == "hallA_wo_ge" || fGeometryName == "hallA_only_WLSR" ||
       fGeometryName == "baseline_smaller" || fGeometryName == "baseline_large_reentrance_tube" || fGeometryName == "baseline_large_reentrance_tube_4m_cryo")
    {
      G4LogicalVolume* logicWater = volumeStore->GetVolume("Water_log");
      biasmuXS->AttachTo(logicWater);
      G4LogicalVolume* logicCout = volumeStore->GetVolume("Cout_log");
      biasmuXS->AttachTo(logicCout);
      G4LogicalVolume* logicCinn = volumeStore->GetVolume("Cinn_log");
      biasmuXS->AttachTo(logicCinn);
      G4LogicalVolume* logicCLid = volumeStore->GetVolume("Lid_log");
      biasmuXS->AttachTo(logicCLid);
      G4LogicalVolume* logicCBot = volumeStore->GetVolume("Bot_log");
      biasmuXS->AttachTo(logicCBot);
    }
    // Alternative has the membrane and insulator
    else if(fGeometryName == "alternative")
    {
      G4LogicalVolume* logicPu = volumeStore->GetVolume("Pu_log");
      biasmuXS->AttachTo(logicPu);
      G4LogicalVolume* logicMembrane = volumeStore->GetVolume("Membrane_log");
      biasmuXS->AttachTo(logicMembrane);
    }
  }
  else
  {
    G4cout << " >>> fSD has entry. Repeated call." << G4endl;
  }
}

auto WLGDDetectorConstruction::SetupAlternative() -> G4VPhysicalVolume*
{
  // Get materials
  auto* worldMaterial = G4Material::GetMaterial("G4_Galactic");
  //auto* larMat        = G4Material::GetMaterial("G4_lAr");
  auto* airMat        = G4Material::GetMaterial("G4_AIR");
  auto* steelMat      = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  auto* copperMat     = G4Material::GetMaterial("G4_Cu");
  auto* stdRock       = G4Material::GetMaterial("StdRock");
  auto* puMat         = G4Material::GetMaterial("polyurethane");
  auto* roiMat        = G4Material::GetMaterial("enrGe");
  auto* larMat_alt    = G4Material::GetMaterial("CombinedArXeHe3");

  if(fXeConc != 0 || fHe3Conc != 0)
    larMat = larMat_alt;
  // size parameter, unit [cm]
  // cavern
  G4double stone     = 100.0;  // Hall wall thickness 1 m
  G4double hallhside = 850.0;  // Hall cube side 17 m
  // cryostat
  G4double tankhside  = 650;   // cryostat cube side 13 m
  G4double outerwall  = 1.2;   // outer SS wall thickness
  G4double insulation = 80.0;  // polyurethane foam
  G4double innerwall  = 0.12;  // inner SS membrane
  // copper tubes with Germanium ROI
  G4double copper    = 0.35;    // tube thickness 3.5 mm
  G4double curad     = 40.0;    // copper tube diam 80 cm
  G4double cuhheight = 334.34;  // copper tube height 7 m inside cryostat
  G4double cushift   = 234.34;  // shift cu tube inside cryostat to top
  G4double ringrad   = 100.0;   // cu tube placement ring radius
  // Ge cylinder for 250 kg at 5.32 g/cm3
  G4double roiradius = 30.0;  // string radius curad - Ge radius - gap

  G4double gerad          = 4.0;                    // Ge radius
  G4double gehheight      = 5.0;                    // full height 10 cm
  G4double gegap          = 3.0;                    // gap between Ge 3cm
  G4double layerthickness = gegap + 2 * gehheight;  // 13 cm total
  G4int    nofLayers      = 8;   // 8 Ge + 7 gaps = 1010 mm string height
  G4int    nofStrings     = 12;  // 12 strings  of 8 Ge each

  // total
  G4double offset =
    hallhside - tankhside;  // shift cavern floor to keep detector centre at origin
  G4double worldside = hallhside + stone + offset + 0.1;  // larger than rest
  G4double larside =
    tankhside - outerwall - insulation - innerwall;  // cube side of LAr volume

  fvertexZ = (worldside - stone - 0.1) * cm;  // max vertex height
  fmaxrad  = hallhside * cm;                  // max vertex circle radius

  // Volumes for this geometry

  //
  // World
  //
  auto* worldSolid = new G4Box("World", worldside * cm, worldside * cm, worldside * cm);
  auto* fWorldLogical  = new G4LogicalVolume(worldSolid, worldMaterial, "World_log");
  auto* fWorldPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fWorldLogical,
                                           "World_phys", nullptr, false, 0);

  //
  // Cavern
  //
  auto* cavernSolid    = new G4Box("Cavern", (hallhside + stone) * cm,
                                   (hallhside + stone) * cm, (hallhside + stone) * cm);
  auto* fCavernLogical = new G4LogicalVolume(cavernSolid, stdRock, "Cavern_log");
  auto* fCavernPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., offset * cm), fCavernLogical,
                      "Cavern_phys", fWorldLogical, false, 0);

  //
  // Hall
  //
  auto* hallSolid = new G4Box("Cavern", hallhside * cm, hallhside * cm, hallhside * cm);
  auto* fHallLogical = new G4LogicalVolume(hallSolid, airMat, "Hall_log");
  auto* fHallPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -stone * cm), fHallLogical,
                      "Hall_phys", fCavernLogical, false, 0, true);

  //
  // Tank
  //
  auto* tankSolid    = new G4Box("Tank", tankhside * cm, tankhside * cm, tankhside * cm);
  auto* fTankLogical = new G4LogicalVolume(tankSolid, steelMat, "Tank_log");
  auto* fTankPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -stone * cm), fTankLogical,
                      "Tank_phys", fHallLogical, false, 0, true);

  //
  // Insulator
  //
  auto* puSolid     = new G4Box("Insulator", (tankhside - outerwall) * cm,
                                (tankhside - outerwall) * cm, (tankhside - outerwall) * cm);
  auto* fPuLogical  = new G4LogicalVolume(puSolid, puMat, "Pu_log");
  auto* fPuPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fPuLogical, "Pu_phys",
                                        fTankLogical, false, 0, true);

  //
  // Membrane
  //
  auto* membraneSolid = new G4Box("Membrane", (tankhside - outerwall - insulation) * cm,
                                  (tankhside - outerwall - insulation) * cm,
                                  (tankhside - outerwall - insulation) * cm);
  auto* fMembraneLogical = new G4LogicalVolume(membraneSolid, steelMat, "Membrane_log");
  auto* fMembranePhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(), fMembraneLogical, "Membrane_phys",
                      fPuLogical, false, 0, true);

  //
  // LAr
  //
  auto* larSolid     = new G4Box("LAr", larside * cm, larside * cm, larside * cm);
  auto* fLarLogical  = new G4LogicalVolume(larSolid, larMat, "Lar_log");
  auto* fLarPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fLarLogical,
                                         "Lar_phys", fMembraneLogical, false, 0, true);

  //
  // copper tubes, hollow cylinder shell
  //
  auto* copperSolid = new G4Tubs("Copper", (curad - copper) * cm, curad * cm,
                                 cuhheight * cm, 0.0, CLHEP::twopi);

  //
  // ULAr bath, solid cylinder
  //
  auto* ularSolid = new G4Tubs("ULar", 0.0 * cm, (curad - copper) * cm, cuhheight * cm,
                               0.0, CLHEP::twopi);

  // tower; logical volumes
  auto* fCopperLogical = new G4LogicalVolume(copperSolid, copperMat, "Copper_log");
  auto* fUlarLogical   = new G4LogicalVolume(ularSolid, larMat, "ULar_log");

  //
  // Germanium, solid cylinder
  //
  // layers in tower
  auto* layerSolid = new G4Tubs("LayerSolid", 0.0 * cm, gerad * cm,
                                (gehheight + gegap / 2.0) * cm, 0.0, CLHEP::twopi);

  auto* fLayerLogical = new G4LogicalVolume(layerSolid, larMat, "Layer_log");

  // fill one layer
  auto* geSolid =
    new G4Tubs("ROI", 0.0 * cm, gerad * cm, gehheight * cm, 0.0, CLHEP::twopi);

  auto* fGeLogical = new G4LogicalVolume(geSolid, roiMat, "Ge_log");
  new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, -gegap / 2.0 * cm), fGeLogical,
                    "Ge_phys", fLayerLogical, false, 0, true);

  auto* gapSolid =
    new G4Tubs("Gap", 0.0 * cm, gerad * cm, gegap / 2.0 * cm, 0.0, CLHEP::twopi);

  auto* fGapLogical = new G4LogicalVolume(gapSolid, larMat, "Gap_log");
  new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, gehheight * cm), fGapLogical,
                    "Gap_phys", fLayerLogical, false, 0, true);

  // place layers as mother volume with unique copy number
  G4double step = (gehheight + gegap / 2) * cm;
  G4double xpos;
  G4double ypos;
  G4double angle = CLHEP::twopi / nofStrings;

  // layer logical into ULarlogical
  for(G4int j = 0; j < nofStrings; j++)
  {
    xpos = roiradius * cm * std::cos(j * angle);
    ypos = roiradius * cm * std::sin(j * angle);
    for(G4int i = 0; i < nofLayers; i++)
    {
      new G4PVPlacement(
        nullptr,
        G4ThreeVector(xpos, ypos,
                      -step + (nofLayers / 2 * layerthickness - i * layerthickness) * cm),
        fLayerLogical, "Layer_phys", fUlarLogical, false, i + j * nofLayers, true);
    }
  }

  // placements
  new G4PVPlacement(nullptr, G4ThreeVector(ringrad * cm, 0., cushift * cm),
                    fCopperLogical, "Copper_phys", fLarLogical, false, 0, true);

  new G4PVPlacement(nullptr, G4ThreeVector(ringrad * cm, 0., cushift * cm), fUlarLogical,
                    "ULar_phys", fLarLogical, false, 0, true);

  // tower 2
  new G4PVPlacement(nullptr, G4ThreeVector(0., ringrad * cm, cushift * cm),
                    fCopperLogical, "Copper_phys2", fLarLogical, false, 1, true);

  new G4PVPlacement(nullptr, G4ThreeVector(0., ringrad * cm, cushift * cm), fUlarLogical,
                    "ULar_phys2", fLarLogical, false, 1, true);

  // tower 3
  new G4PVPlacement(nullptr, G4ThreeVector(-ringrad * cm, 0., cushift * cm),
                    fCopperLogical, "Copper_phys3", fLarLogical, false, 2, true);

  new G4PVPlacement(nullptr, G4ThreeVector(-ringrad * cm, 0., cushift * cm), fUlarLogical,
                    "ULar_phys3", fLarLogical, false, 2, true);

  // tower 4
  new G4PVPlacement(nullptr, G4ThreeVector(0., -ringrad * cm, cushift * cm),
                    fCopperLogical, "Copper_phys4", fLarLogical, false, 3, true);

  new G4PVPlacement(nullptr, G4ThreeVector(0., -ringrad * cm, cushift * cm), fUlarLogical,
                    "ULar_phys4", fLarLogical, false, 3, true);

  //
  // Visualization attributes
  //
 // fWorldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());

  auto* redVisAtt = new G4VisAttributes(G4Colour::Red());
  redVisAtt->SetVisibility(true);
  auto* greyVisAtt = new G4VisAttributes(G4Colour::Grey());
  greyVisAtt->SetVisibility(true);
  auto* greenVisAtt = new G4VisAttributes(G4Colour::Green());
  greenVisAtt->SetVisibility(true);
  auto* blueVisAtt = new G4VisAttributes(G4Colour::Blue());
  blueVisAtt->SetVisibility(true);

  fCavernLogical->SetVisAttributes(redVisAtt);
  fHallLogical->SetVisAttributes(greyVisAtt);
  fTankLogical->SetVisAttributes(blueVisAtt);
  fPuLogical->SetVisAttributes(greyVisAtt);
  fMembraneLogical->SetVisAttributes(blueVisAtt);
  fLarLogical->SetVisAttributes(greyVisAtt);
  fCopperLogical->SetVisAttributes(greenVisAtt);
  fUlarLogical->SetVisAttributes(greyVisAtt);
  fGeLogical->SetVisAttributes(redVisAtt);

  return fWorldPhysical;
}

auto WLGDDetectorConstruction::SetupBaseline() -> G4VPhysicalVolume*
{
  // Get materials
  auto* worldMaterial = G4Material::GetMaterial("G4_Galactic");
  // auto* larMat        = G4Material::GetMaterial("G4_lAr");
  auto* airMat   = G4Material::GetMaterial("G4_AIR");
  auto* waterMat = G4Material::GetMaterial("G4_WATER");
  if(fWithGdWater == 1)
    waterMat = water;
  if(fWithWoWater == 1)
    waterMat = airMat;
  auto* steelMat      = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  auto* copperMat     = G4Material::GetMaterial("G4_Cu");
  auto* stdRock       = G4Material::GetMaterial("StdRock");
  auto* roiMat        = G4Material::GetMaterial("enrGe");
  auto* BoratedPETMat = G4Material::GetMaterial("BoratedPET");
  auto* GdLoadedPMMA = G4Material::GetMaterial("GdLoadedPMMA");
  if(fSetMaterial == "PolyEthylene")
    BoratedPETMat = G4Material::GetMaterial("PolyEthylene");
  if(fSetMaterial == "PMMA")
    BoratedPETMat = G4Material::GetMaterial("PMMA");
  //auto* larMat /*_alt*/ = G4Material::GetMaterial("CombinedArXeHe3");
  // if(fXeConc != 0 || fHe3Conc != 0) BoratedPET
  
  if(!((fXeConc + fHe3Conc) == 0))
    larMat = CombinedArXeHe3;
  //larMat = G4Material::GetMaterial("G4_lAr");

  // G4cout << larMat << G4endl;

  // Edit: 2021/03/30 by Moritz Neuberger
  // Adjusted size of stone (1m->5m) s.t. MUSUN cuboid lies inside.
  // Also adjusted relative geometry relations s.t. they are independent of the stone
  // size.

  // constants
  // size parameter, unit [cm]
  G4double offset = 200.0;  // shift cavern floor to keep detector centre at origin
  G4double offset_2 =
    100.0;  // shift s.t. cavern, hall and tank are in line for different stone sizes
  G4double offset_3 = 100;  // 69.5; // shift to get to the baseline lowered position
  if(fDetectorPosition == "original")
    offset_3 = 0;
  if(fGeometryName == "baseline_smaller"  || fGeometryName == "baseline_large_reentrance_tube_4m_cryo")
    offset_3 = 37.5;

  // cavern
  G4double stone       = 500.0;  // Hall wall thickness 5 m
  G4double hallrad     = 600.0;  // Hall cylinder diam 12 m
  G4double hallhheight = 850.0;  // Hall cylinder height 17 m
  // water tank
  G4double tankwalltop = 0.6;  // water tank thickness at top 6 mm
  G4double tankwallbot = 0.8;  // water tank thickness at bottom 8 mm
  G4double tankrad     = 550;  // water tank diam 11 m
  G4double tankhheight = 650;  // water tank height 13 m
  // cryostat
  G4double cryowall   = 3.0;                   // cryostat wall thickness from GERDA
  G4double vacgap     = 1.0;                   // vacuum gap between walls
  G4double cryrad     = fCryostatOuterRadius;  // 350.0;  // cryostat diam 7 m
  G4double cryhheight = fCryostatHeight;       // 350.0;  // cryostat height 7 m
  // For Optical foil on scaffolding
  G4double ScaffoldRadius = cryrad + 80;

  if(fGeometryName == "baseline_large_reentrance_tube")
  {
    vacgap     = 50.0;  
  }

  if(fGeometryName == "baseline_smaller" || fGeometryName == "baseline_large_reentrance_tube_4m_cryo")
  {
    cryrad     = 200.0;  // cryostat diam 4 m
    cryhheight = 225.0;  // cryostat height 4.5 m
  }
  // Borated PET tubes around copper tubes
  G4double BoratedPETouterrad = 5.0;  // tube thickness 5 cm
  // copper tubes with Germanium ROI
  G4double copper    = 0.35;  // tube thickness 3.5 mm
  G4double curad     = 40.0;  // copper tube diam 80 cm
  G4double cuhheight = (400 - (350 - fCryostatHeight)) / 2.;  // 200.0;  // copper tube height 4 m inside cryostat
  G4double cushift = fCryostatHeight - cuhheight;  // 150.0;  // shift cu tube inside cryostat to top
  if(fGeometryName == "baseline_smaller" || fGeometryName == "baseline_large_reentrance_tube_4m_cryo")
  {
    cuhheight = 137.5;  // cupper height 2.25 m
    cushift   = 87.5;   // shift
  }
  if(fGeometryName == "baseline_large_reentrance_tube" || fGeometryName == "baseline_large_reentrance_tube_4m_cryo")
  {
    curad     = 95.0;  
  }
  G4double ringrad = 100.0;  // cu tube placement ring radius
  // Ge cylinder 2.67 kg at 5.32 g/cm3
  G4double roiradius = 30.0;  // string radius curad - Ge radius - gap
  // total mass 1026.86 kg in 4 towers, each with 8 Ge stacked in 12 strings
  G4double gerad          = 4.0;                    // Ge radius
  G4double gehheight      = 5.0;                    // full height 10 cm
  G4double gegap          = 3.0;                    // gap between Ge 3cm
  G4double layerthickness = gegap + 2 * gehheight;  // 13 cm total
  //    G4int nofLayers = 8;   // 8 Ge + 7 gaps = 1010 mm string height
  //    G4int nofStrings = 12;  // 12 strings  of 8 Ge each
  G4int nofLayers  = 7;   // 8 Ge + 7 gaps = 1010 mm string height
  G4int nofStrings = 14;  // 12 strings  of 8 Ge each

  fvertexZ = (hallhheight + offset) * cm;
  fmaxrad  = hallrad * cm;

  // Volumes for this geometry

  //
  // World
  //
  auto* worldSolid =
    new G4Tubs("World", 0.0 * cm, (hallrad + stone + 0.1) * cm,
               (hallhheight + stone + offset + 0.1) * cm, 0.0, CLHEP::twopi);
  auto* fWorldLogical  = new G4LogicalVolume(worldSolid, worldMaterial, "World_log");
  auto* fWorldPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fWorldLogical,
                                           "World_phys", nullptr, false, 0);

  //
  // Cavern
  //
  auto* cavernSolid    = new G4Tubs("Cavern", 0.0 * cm, (hallrad + stone) * cm,
                                    (hallhheight + stone) * cm, 0.0, CLHEP::twopi);
  auto* fCavernLogical = new G4LogicalVolume(cavernSolid, stdRock, "Cavern_log");
  auto* fCavernPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., offset * cm), fCavernLogical,
                      "Cavern_phys", fWorldLogical, false, 0);

  //
  // Hall
  //
  auto* hallSolid =
    new G4Tubs("Hall", 0.0 * cm, hallrad * cm, hallhheight * cm, 0.0, CLHEP::twopi);
  auto* fHallLogical = new G4LogicalVolume(hallSolid, airMat, "Hall_log");
  auto* fHallPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., /*-stone*/ -offset_2 * cm),
                      fHallLogical, "Hall_phys", fCavernLogical, false, 0, true);

  //
  // Tank
  //
  auto* tankSolid =
    new G4Cons("Tank", 0.0 * cm, (tankrad + tankwallbot) * cm, 0.0 * cm,
               (tankrad + tankwalltop) * cm, tankhheight * cm, 0.0, CLHEP::twopi);
  auto* fTankLogical  = new G4LogicalVolume(tankSolid, steelMat, "Tank_log");
  auto* fTankPhysical = new G4PVPlacement(
    nullptr,
    G4ThreeVector(0., 0., /*-stone -offset_2*/ -(hallhheight - tankhheight) * cm),
    fTankLogical, "Tank_phys", fHallLogical, false, 0, true);

  //
  // Water
  //
  auto* OuterWaterSolid     = new G4Tubs("OuterWater_sol", 0.0 * cm, tankrad * cm,
                                    (tankhheight - tankwallbot) * cm, 0.0, CLHEP::twopi);
  auto* OuterWaterLogical  = new G4LogicalVolume(OuterWaterSolid, waterMat, "OuterWater_log");
  auto* OuterWaterPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), OuterWaterLogical,
                                           "Water_phys", fTankLogical, false, 0, true);

  //
  // Seperate outer water volume with reflective foil on Scaffolding from inner optical Water volume
  //                                         
  auto* ScaffoldSolid = new G4Tubs("Scaffold_sol", 0.0 * cm, (ScaffoldRadius + 0.1) * cm,
                                    ScaffoldRadius * cm, 0.0, CLHEP::twopi);
  auto* ScaffoldLogical  = new G4LogicalVolume(ScaffoldSolid, steelMat, "Scaffold_log");
  auto* ScaffoldPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), ScaffoldLogical,
                                         "Scaffold_phys", OuterWaterLogical, false, 0, true);
  
  //
  // Now add inner optical water volume
  //
  auto* InnerWaterSolid = new G4Tubs("InnerWater_sol", 0.0 * cm, ScaffoldRadius * cm,
                                      ScaffoldRadius * cm, 0.0, CLHEP::twopi);
  auto* InnerWaterLogical  = new G4LogicalVolume(InnerWaterSolid, fOpticalWater, "InnerWater_log");
  auto* InnerWaterPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), InnerWaterLogical,
                                           "InnerWater_phys", ScaffoldLogical, false, 0, true);

  //
  //Place PMTs. 500 PMTs to be exact
  //
  G4double PMTradius = 12.7;
  auto* PMTsolid = new G4Tubs("PMT_sol", 0.0 * cm, PMTradius * cm,
                                   0.5 * cm, 0.0, CLHEP::twopi); // 10 inch PMTs, just make them as 1cm surface for now
  auto* PMTLogical = new G4LogicalVolume(PMTsolid, steelMat, "PMT_log");

  G4int PMTperRow = 45;      // For now
  G4int Rows = 11;           // Current setup. Will probably change
  G4double VerticalPMTDistance = 50.;

  // Now calculate the transformation and place all the PMTs. This is where the fun begins
  // First spin the PMT to face x-axis
  G4Rotate3D rotateY(90*deg, G4ThreeVector(0,1,0));
  for(G4int j = 0; j < Rows; j++)
  {
      for(G4int i = 0; i < PMTperRow; i++)
      {
        // Move PMT to the correct radius and z-position
        G4Translate3D transX(G4ThreeVector((ScaffoldRadius - 1.)*cm,0., 
                            (-5*(VerticalPMTDistance + PMTradius) + j*(VerticalPMTDistance + PMTradius))*cm));
        // Rotate PMT position (phi in cylindrical coordinates). Will always face center like this cause whole object is rotated.
        G4Rotate3D rotateZ(i*(360./(PMTperRow))*deg,G4ThreeVector(0,0,1));
        G4Transform3D transformPMT = rotateZ * transX * rotateY; // Matrix multiplication from right to left
        G4int DetectorNumber = (j+1) * 1000 + (i+1); // Detectornumber in the format of row column e.g. (1003) row 1 column 3
        G4String PMTname =  "PMT" + std::to_string(DetectorNumber);

        G4cout << "i: " << i << "j: " << j << "Name: " << PMTname << G4endl;
        auto* PMTPhysical = new G4PVPlacement(transformPMT, PMTLogical, PMTname, InnerWaterLogical, false, 0, true);
      }
  }

  //
  // PMMA around cryostat for Neutrontagger
  //
  auto* PMMAsolid =
    new G4Tubs("PMMA_sol", 0.0 * cm, (cryrad + 5) * cm, (cryhheight + 5) * cm, 0.0, CLHEP::twopi);
  auto* fPMMALogical  = new G4LogicalVolume(PMMAsolid, GdLoadedPMMA, "PMMA_log");
  auto* fPMMAPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fPMMALogical,
                                          "PMMA_phys", InnerWaterLogical, false, 0, true);
  //
  // outer cryostat
  //
  auto* coutSolid =
    new G4Tubs("Cout", 0.0 * cm, cryrad * cm, cryhheight * cm, 0.0, CLHEP::twopi);
  auto* fCoutLogical  = new G4LogicalVolume(coutSolid, steelMat, "Cout_log");
  auto* fCoutPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fCoutLogical,
                                          "Cout_phys", fPMMALogical, false, 0, true);

  //
  // vacuum gap
  //
  auto* cvacSolid     = new G4Tubs("Cvac", 0.0 * cm, (cryrad - cryowall) * cm,
                                   (cryhheight - cryowall) * cm, 0.0, CLHEP::twopi);
  auto* fCvacLogical  = new G4LogicalVolume(cvacSolid, worldMaterial, "Cvac_log");
  auto* fCvacPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fCvacLogical,
                                          "Cvac_phys", fCoutLogical, false, 0, true);

  //
  // inner cryostat
  //
  auto* cinnSolid     = new G4Tubs("Cinn", 0.0 * cm, (cryrad - cryowall - vacgap) * cm,
                                   (cryhheight - cryowall - vacgap) * cm, 0.0, CLHEP::twopi);
  auto* fCinnLogical  = new G4LogicalVolume(cinnSolid, steelMat, "Cinn_log");
  auto* fCinnPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fCinnLogical,
                                          "Cinn_phys", fCvacLogical, false, 0, true);

  //
  // LAr bath
  //
  auto* larSolid     = new G4Tubs("LAr", 0.0 * cm, (cryrad - 2 * cryowall - vacgap) * cm,
                                  (cryhheight - 2 * cryowall - vacgap) * cm, 0.0, CLHEP::twopi);
  auto* fLarLogical  = new G4LogicalVolume(larSolid, larMat, "Lar_log");
  auto* fLarPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fLarLogical,
                                         "Lar_phys", fCinnLogical, false, 0, true);

  //
  // cryostat Lid
  //
  auto* lidSolid =
    new G4Tubs("Lid", 0.0 * cm, cryrad * cm, cryowall / 2.0 * cm, 0.0, CLHEP::twopi);
  auto* fLidLogical = new G4LogicalVolume(lidSolid, steelMat, "Lid_log");
  auto* fLidPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., (cryhheight + cryowall / 2.0) * cm),
                      fLidLogical, "Lid_phys", fPMMALogical, false, 0, true);
  auto* fBotLogical = new G4LogicalVolume(lidSolid, steelMat, "Bot_log");
  auto* fBotPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -(cryhheight + cryowall / 2.0) * cm),
                      fBotLogical, "Bot_phys", fPMMALogical, false, 0, true);

  //
  // copper tubes, hollow cylinder shell
  //


  // Box variables
  G4double b_width  = fBoratedTurbineWidth / 2. * cm;   // 2.5 * cm;
  G4double b_length = fBoratedTurbineLength / 2. * cm;  // 0.25 * m;
  G4double b_height = fBoratedTurbineHeight / 2. * cm;  // fCryostatHeight * cm - 0.5 * m;

 /* auto* boratedPETSolid_Tube =
    new G4Tubs("BoratedPET", curad * cm, (BoratedPETouterrad + curad) * cm,
               cuhheight * cm, 0.0, CLHEP::twopi);
 */

  auto* boratedPETSolid_Tube =
    new G4Tubs("BoratedPET", curad * cm, (2*b_width + curad) * cm,
               cuhheight * cm, 0.0, CLHEP::twopi);

  auto* boratedPETSolid_Box = new G4Box("BoratedPET", b_length, b_width, b_height);

  //
  // copper tubes, hollow cylinder shell
  auto* copperSolid = new G4Tubs("Copper", (curad - copper) * cm, curad * cm,
                                 cuhheight * cm, 0.0, CLHEP::twopi);
  //
  // ULAr bath, solid cylinder
  auto* ularSolid = new G4Tubs("ULar", 0.0 * cm, (curad - copper) * cm,
                               (cuhheight - copper) * cm, 0.0, CLHEP::twopi);

  // tower; logical volumes
  auto* fBoratedPETLogical_Tube =
    new G4LogicalVolume(boratedPETSolid_Tube, BoratedPETMat, "BoratedPET_Logical");
  auto* fBoratedPETLogical_Box =
    new G4LogicalVolume(boratedPETSolid_Box, BoratedPETMat, "BoratedPET_Logical");
  auto* fCopperLogical = new G4LogicalVolume(copperSolid, copperMat, "Copper_log");
  auto* fUlarLogical   = new G4LogicalVolume(ularSolid, larMat, "ULar_log");

  //
  // Germanium, solid cylinder
  //
  // layers in tower
  auto* layerSolid = new G4Tubs("LayerSolid", 0.0 * cm, gerad * cm,
                                (gehheight + gegap / 2.0) * cm, 0.0, CLHEP::twopi);

  auto* fLayerLogical = new G4LogicalVolume(layerSolid, larMat, "Layer_log");

  // fill one layer
  auto* geSolid =
    new G4Tubs("ROI", 0.0 * cm, gerad * cm, gehheight * cm, 0.0, CLHEP::twopi);

  auto* fGeLogical = new G4LogicalVolume(geSolid, roiMat, "Ge_log");
  new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, -gegap / 2.0 * cm), fGeLogical,
                    "Ge_phys", fLayerLogical, false, 0, true);

  auto* gapSolid =
    new G4Tubs("Gap", 0.0 * cm, gerad * cm, gegap / 2.0 * cm, 0.0, CLHEP::twopi);

  auto* fGapLogical = new G4LogicalVolume(gapSolid, larMat, "Gap_log");
  new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, gehheight * cm), fGapLogical,
                    "Gap_phys", fLayerLogical, false, 0, true);

  // place layers as mother volume with unique copy number
  G4double step = (gehheight + gegap / 2) * cm;
  G4double xpos;
  G4double ypos;
  G4double angle = CLHEP::twopi / nofStrings;

  // layer logical into ULarlogical
  if(fGeometryName == "baseline_large_reentrance_tube" || fGeometryName == "baseline_large_reentrance_tube_4m_cryo"){
    G4double length = 16.75 * cm;

    G4double vec_main_x = 1/2.;
    G4double vec_main_y = sqrt(1 - vec_main_x*vec_main_x);
    G4double vec_left_x = -1/2.;
    G4double vec_left_y = sqrt(1 - vec_left_x*vec_left_x);
    G4double vec_right_x = 1;
    G4double vec_right_y = 0;

    int N = 0;

    for(G4int i = -4; i <= 4; i++){

      for(G4int j = 0; j <= 4 && j <= 4 - i; j++){

        // leave out calibration ports
        if( ( i == -3 && j == 0 ) || ( i == 0 && j == 0 ) || ( i == 3 && j == 0) || ( i == -3 && j == 3) || ( i == 0 && j == 3) ) continue;

        xpos = length * ( vec_main_x * i + vec_left_x * j);
        ypos = length * ( vec_main_y * i + vec_left_y * j) ;

        G4cout << "i: " << i << " j: " << j << " | x: " << xpos << " - y: " << ypos << G4endl;

        for(G4int k = 0; k < nofLayers; k++)
        {
          // Cube coordinates
          int x_coor = i;
          int y_coor = j;
          int z_coor = -(i + j);
          int coordinate = N;//(x_coor<0)*1e6 + abs(x_coor)*1e5 + (y_coor<0)*1e4 + abs(y_coor)*1e3 + (z_coor<0)*1e2 + abs(z_coor)*1e1 + k;
          G4cout << "coordinate: " << coordinate << G4endl;
          new G4PVPlacement(
            nullptr,
            G4ThreeVector(xpos, ypos,
                          -step + (nofLayers / 2 * layerthickness - k * layerthickness) * cm -
                            offset_3 * cm),
            fLayerLogical, "Layer_phys", fUlarLogical, false, coordinate, true);
            N++;
        }
      }

      for(G4int j = 1; j <= 4 && j <= 4 - i; j++){

        // leave out calibration ports
        if( ( i == -3 && j == 3 ) || ( i == 0 && j == 3 ) ) continue;

        xpos = length * ( vec_main_x * i + vec_right_x * j);
        ypos = length * ( vec_main_y * i + vec_right_y * j) ;

        for(G4int k = 0; k < nofLayers; k++)
        {
          // Cube coordinates
          int x_coor = i + j;
          int y_coor = -j;
          int z_coor = -i;
          int coordinate = N;//(x_coor<0)*1e6 + abs(x_coor)*1e5 + (y_coor<0)*1e4 + abs(y_coor)*1e3 + (z_coor<0)*1e2 + abs(z_coor)*1e1 + k;
          G4cout << "coordinate: " << coordinate << G4endl;
          new G4PVPlacement(
            nullptr,
            G4ThreeVector(xpos, ypos,
                          -step + (nofLayers / 2 * layerthickness - k * layerthickness) * cm -
                            offset_3 * cm),
            fLayerLogical, "Layer_phys", fUlarLogical, false, coordinate, true);
            N++;
        }
      }
    }
  }
  else{
    for(G4int j = 0; j < nofStrings; j++)
    {
      xpos = roiradius * cm * std::cos(j * angle);
      ypos = roiradius * cm * std::sin(j * angle);

      for(G4int i = 0; i < nofLayers; i++)
      {
        new G4PVPlacement(
          nullptr,
          G4ThreeVector(xpos, ypos,
                        -step + (nofLayers / 2 * layerthickness - i * layerthickness) * cm -
                          offset_3 * cm),
          fLayerLogical, "Layer_phys", fUlarLogical, false, i + j * nofLayers, true);
      }
    }
  }

  if(fGeometryName == "baseline")
  {
    // placements
    if(fWithBoratedPET == 1)
      new G4PVPlacement(nullptr, G4ThreeVector(ringrad * cm, 0., cushift * cm),
                        fBoratedPETLogical_Tube, "BoratedPET_phys", fLarLogical, false, 0,
                        true);

    if(fWithOutCupperTubes == 0)
      new G4PVPlacement(nullptr, G4ThreeVector(ringrad * cm, 0., cushift * cm),
                        fCopperLogical, "Copper_phys", fLarLogical, false, 0, true);

    new G4PVPlacement(nullptr, G4ThreeVector(ringrad * cm, 0., cushift * cm),
                      fUlarLogical, "ULar_phys", fLarLogical, false, 0, true);

    // tower 2
    if(fWithBoratedPET == 1)
      new G4PVPlacement(nullptr, G4ThreeVector(0., ringrad * cm, cushift * cm),
                        fBoratedPETLogical_Tube, "BoratedPET_phys2", fLarLogical, false,
                        1, true);

    if(fWithOutCupperTubes == 0)
      new G4PVPlacement(nullptr, G4ThreeVector(0., ringrad * cm, cushift * cm),
                        fCopperLogical, "Copper_phys2", fLarLogical, false, 1, true);

    new G4PVPlacement(nullptr, G4ThreeVector(0., ringrad * cm, cushift * cm),
                      fUlarLogical, "ULar_phys2", fLarLogical, false, 1, true);

    // tower 3
    if(fWithBoratedPET == 1)
      new G4PVPlacement(nullptr, G4ThreeVector(-ringrad * cm, 0., cushift * cm),
                        fBoratedPETLogical_Tube, "BoratedPET_phys3", fLarLogical, false,
                        2, true);

    if(fWithOutCupperTubes == 0)
      new G4PVPlacement(nullptr, G4ThreeVector(-ringrad * cm, 0., cushift * cm),
                        fCopperLogical, "Copper_phys3", fLarLogical, false, 2, true);

    new G4PVPlacement(nullptr, G4ThreeVector(-ringrad * cm, 0., cushift * cm),
                      fUlarLogical, "ULar_phys3", fLarLogical, false, 2, true);

    // tower 4
    if(fWithBoratedPET == 1)
      new G4PVPlacement(nullptr, G4ThreeVector(0., -ringrad * cm, cushift * cm),
                        fBoratedPETLogical_Tube, "BoratedPET_phys4", fLarLogical, false,
                        3, true);

    if(fWithOutCupperTubes == 0)
      new G4PVPlacement(nullptr, G4ThreeVector(0., -ringrad * cm, cushift * cm),
                        fCopperLogical, "Copper_phys4", fLarLogical, false, 3, true);

    new G4PVPlacement(nullptr, G4ThreeVector(0., -ringrad * cm, cushift * cm),
                      fUlarLogical, "ULar_phys4", fLarLogical, false, 3, true);
  }

  if(fGeometryName == "baseline_smaller" || fGeometryName == "baseline_large_reentrance_tube" || fGeometryName == "baseline_large_reentrance_tube_4m_cryo")
  {
    // placements
    if(fWithBoratedPET == 1)
      new G4PVPlacement(nullptr, G4ThreeVector(0., 0., cushift * cm),
                        fBoratedPETLogical_Tube, "BoratedPET_phys", fLarLogical, false, 0,
                        true);

    if(fWithOutCupperTubes == 0)
      new G4PVPlacement(nullptr, G4ThreeVector(0., 0., cushift * cm), fCopperLogical,
                        "Copper_phys", fLarLogical, false, 0, true);

    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., cushift * cm), fUlarLogical,
                      "ULar_phys", fLarLogical, false, 0, true);
  }

#define whichGeometry 1

#if whichGeometry == 0
  if(fWithBoratedPET == 2)
  {
    int               NPanels        = 28;
    double            radiusOfPanels = 2 * m;
    double            anglePanel     = 360 / 28. * deg;
    G4double          zpos           = 0 * cm;
    G4RotationMatrix* rotMat;
    for(G4int j = 0; j < NPanels; j++)
    {
      xpos   = radiusOfPanels * std::cos(j * anglePanel);
      ypos   = radiusOfPanels * std::sin(j * anglePanel);
      rotMat = new G4RotationMatrix;
      rotMat->rotateZ(-(j + 1) * anglePanel + 90 * deg);
      new G4PVPlacement(rotMat, G4ThreeVector(xpos, ypos, zpos), fBoratedPETLogical_Box,
                        "BoratedPET_phys", fLarLogical, false, j, true);
    }
  }
#endif
#if whichGeometry == 1

  if(fWithBoratedPET == 2 || fWithBoratedPET == 4)
  {
    G4double densityOfBPE   = 0.95;
    double   radiusOfPanels = fBoratedTurbineRadius * cm;
    double   constantAngle  = fBoratedTurbineAngle * deg;  // 45 * deg;
    int      NPanels;

    if(fBoratedTurbineNPanels == 0)
      NPanels = ceil(2 * 3.14159265 * radiusOfPanels / cm /
                     (0.95 * b_length / cm * cos(constantAngle)));
    else
      NPanels = fBoratedTurbineNPanels;

    fNPanels          = NPanels;
    double anglePanel = 360. / NPanels * deg;

    G4double totalVolume =
      NPanels * 2 * b_length / cm * 2 * b_width / cm * 2 * b_height / cm;

    G4cout << "Total Mass of B-PE: " << totalVolume * densityOfBPE << G4endl;

    G4double          zpos = 0 * cm;
    G4RotationMatrix* rotMat;

    for(G4int j = 0; j < NPanels; j++)
    {
      xpos   = radiusOfPanels * std::cos(j * anglePanel);
      ypos   = radiusOfPanels * std::sin(j * anglePanel);
      rotMat = new G4RotationMatrix;
      rotMat->rotateZ(-j * anglePanel + 90 * deg + constantAngle);
      new G4PVPlacement(rotMat, G4ThreeVector(xpos, ypos, zpos), fBoratedPETLogical_Box,
                        "BoratedPET_phys", fLarLogical, false, j, true);
    }

    if(fWithBoratedPET == 4)
    {
      boratedPETSolid_Tube =
        new G4Tubs("BoratedPET", 0, (fBoratedTurbineRadius * cm + b_width * 2), b_width,
                   0.0, CLHEP::twopi);
      fBoratedPETLogical_Tube =
        new G4LogicalVolume(boratedPETSolid_Tube, BoratedPETMat, "BoratedPET_Logical");

      new G4PVPlacement(
        nullptr, G4ThreeVector(0, 0, fBoratedTurbinezPosition * cm - b_height - b_width),
        fBoratedPETLogical_Tube, "BoratedPET_phys", fLarLogical, false, 0, true);
      new G4PVPlacement(
        nullptr, G4ThreeVector(0, 0, fBoratedTurbinezPosition * cm + b_height + b_width),
        fBoratedPETLogical_Tube, "BoratedPET_phys", fLarLogical, false, 0, true);
    }
  }
#endif

  if(fWithBoratedPET == 3)
  {
    G4double densityOfBPE   = 0.95;
    double   radiusOfPanels = fBoratedTurbineRadius * cm;

    boratedPETSolid_Tube =
      new G4Tubs("BoratedPET", fBoratedTurbineRadius * cm,
                 (fBoratedTurbineRadius * cm + b_width * 2), b_height, 0.0, CLHEP::twopi);
    fBoratedPETLogical_Tube =
      new G4LogicalVolume(boratedPETSolid_Tube, BoratedPETMat, "BoratedPET_Logical");

    G4cout << "Total Mass of B-PE: "
           << 3.141592653589 * b_height / cm *
                (pow(fBoratedTurbineRadius + b_width / cm * 2, 2) -
                 pow(fBoratedTurbineRadius, 2)) *
                densityOfBPE
           << G4endl;

    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, fBoratedTurbinezPosition * cm),
                      fBoratedPETLogical_Tube, "BoratedPET_phys", fLarLogical, false, 0,
                      true);

    boratedPETSolid_Tube =
      new G4Tubs("BoratedPET", 0, (fBoratedTurbineRadius * cm + b_width * 2), b_width,
                 0.0, CLHEP::twopi);
    fBoratedPETLogical_Tube =
      new G4LogicalVolume(boratedPETSolid_Tube, BoratedPETMat, "BoratedPET_Logical");

    new G4PVPlacement(
      nullptr, G4ThreeVector(0, 0, fBoratedTurbinezPosition * cm - b_height - b_width),
      fBoratedPETLogical_Tube, "BoratedPET_phys", fLarLogical, false, 0, true);
    new G4PVPlacement(
      nullptr, G4ThreeVector(0, 0, fBoratedTurbinezPosition * cm + b_height + b_width),
      fBoratedPETLogical_Tube, "BoratedPET_phys", fLarLogical, false, 0, true);
  }

  //
  // Visualization attributes
  //
  //fWorldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4Color testColor(0., 109 / 225., 119 / 225.);
  auto*   testVisAtt = new G4VisAttributes(testColor);
  testVisAtt->SetVisibility(true);

  G4Color testColor2(131 / 255., 197 / 225., 190 / 225.);
  auto*   testVisAtt2 = new G4VisAttributes(testColor2);
  testVisAtt2->SetVisibility(true);

  G4Color testColor3(226 / 255., 149 / 225., 120 / 225.);
  auto*   testVisAtt3 = new G4VisAttributes(testColor3);
  testVisAtt3->SetVisibility(true);

  G4Color testColor4(255 / 255., 221 / 225., 210 / 225.);
  auto*   testVisAtt4 = new G4VisAttributes(testColor4);
  testVisAtt4->SetVisibility(true);


  G4Color testColor_copper(1., 190 / 225., 102 / 225.); //255, 190, 182
  auto*   testVisAtt_copper = new G4VisAttributes(testColor_copper);
  testVisAtt_copper->SetVisibility(true);

  G4Color testColor_LAr(153 / 255., 193 / 225., 151 / 225.); //153, 193, 151
  auto*   testVisAtt_LAr = new G4VisAttributes(testColor_LAr);
  testVisAtt_LAr->SetVisibility(true);

  G4Color testColor_Ge(193 / 255., 193 / 225., 193 / 225.); //193, 193, 193
  auto*   testVisAtt_Ge = new G4VisAttributes(testColor_Ge);
  testVisAtt_Ge->SetVisibility(true);

  G4Color testColor_water(170 / 255., 191 / 225., 219 / 225.); //170, 191, 219
  auto*   testVisAtt_water = new G4VisAttributes(testColor_water);
  testVisAtt_water->SetVisibility(true);


  auto* redVisAtt = new G4VisAttributes(G4Colour::Red());
  redVisAtt->SetVisibility(true);
  auto* whiteVisAtt = new G4VisAttributes(G4Colour::White());
  whiteVisAtt->SetVisibility(true);
  auto* orangeVisAtt = new G4VisAttributes(G4Colour::Brown());
  orangeVisAtt->SetVisibility(true);

  auto* greyVisAtt = new G4VisAttributes(G4Colour::Grey());
  greyVisAtt->SetVisibility(true);
  auto* greenVisAtt = new G4VisAttributes(G4Colour::Green());
  greenVisAtt->SetVisibility(true);
  auto* blueVisAtt = new G4VisAttributes(G4Colour::Blue());
  blueVisAtt->SetVisibility(true);

  fCavernLogical->SetVisAttributes(greyVisAtt);
  fHallLogical->SetVisAttributes(whiteVisAtt);
  fTankLogical->SetVisAttributes(greyVisAtt);
  //fWaterLogical->SetVisAttributes(testVisAtt);
  OuterWaterLogical->SetVisAttributes(testVisAtt_water);
  fLarLogical->SetVisAttributes(testVisAtt2);
  fCoutLogical->SetVisAttributes(greyVisAtt);
  fCvacLogical->SetVisAttributes(greyVisAtt);
  fCinnLogical->SetVisAttributes(greyVisAtt);
  fLidLogical->SetVisAttributes(greyVisAtt);
  fBotLogical->SetVisAttributes(greyVisAtt);
  //fCopperLogical->SetVisAttributes(testVisAtt4);
  fCopperLogical->SetVisAttributes(testVisAtt_copper);
  //fUlarLogical->SetVisAttributes(testVisAtt2);
  fUlarLogical->SetVisAttributes(testVisAtt_LAr);
  //fGapLogical->SetVisAttributes(testVisAtt2);
  fGapLogical->SetVisAttributes(testVisAtt_LAr);
  //fGeLogical->SetVisAttributes(testVisAtt3);
  fGeLogical->SetVisAttributes(testVisAtt_Ge);
  fBoratedPETLogical_Tube->SetVisAttributes(testVisAtt4);
  fBoratedPETLogical_Box->SetVisAttributes(testVisAtt4);

  // Set Optical Boundaries
  G4OpticalSurface *OpSurface = new G4OpticalSurface("Teflon");
  G4LogicalBorderSurface *FoilCloseToCryostat = new G4LogicalBorderSurface("FoilCloseToCryostat", InnerWaterPhysical, fPMMAPhysical , OpSurface);
  G4LogicalBorderSurface *FoilFarFromCryostat = new G4LogicalBorderSurface("FoilFarFromCryostat", InnerWaterPhysical, ScaffoldPhysical , OpSurface);
  OpSurface->SetType(dielectric_dielectric);
  OpSurface->SetModel(unified);
  OpSurface->SetFinish(groundfrontpainted);

  G4MaterialPropertiesTable *SMPT = new G4MaterialPropertiesTable();
  std::vector<G4double> reflectivity = {0.95,0.95};
  std::vector<G4double> transmission = {0.,0.};
  SMPT->AddProperty("REFLECTIVITY", fenergySmall, reflectivity);
  SMPT->AddProperty("TRANSMITTANCE", fenergySmall, transmission);
  OpSurface->SetMaterialPropertiesTable(SMPT);

  return fWorldPhysical;
}

auto WLGDDetectorConstruction::SetupHallA() -> G4VPhysicalVolume*
{
  // Full copy of baseline set up but smaller as a Gerda mock-up.
  // Get materials
  auto* worldMaterial = G4Material::GetMaterial("G4_Galactic");
  //auto* larMat        = G4Material::GetMaterial("G4_lAr");
  auto* airMat        = G4Material::GetMaterial("G4_AIR");
  auto* waterMat      = G4Material::GetMaterial("G4_WATER");
  if(fWithGdWater == 1)
    waterMat = water;
  if(fWithWoWater == 1)
    waterMat = airMat;
  auto* steelMat   = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  auto* copperMat  = G4Material::GetMaterial("G4_Cu");
  auto* stdRock    = G4Material::GetMaterial("StdRock");
  auto* roiMat     = G4Material::GetMaterial("enrGe");
  auto* larMat_alt = G4Material::GetMaterial("CombinedArXeHe3");

  if(fXeConc != 0 || fHe3Conc != 0)
    larMat = larMat_alt;
  // constants
  // size parameter, unit [cm]
  G4double offset = 250.0;  // shift cavern floor to keep detector centre at origin
  // cavern
  G4double stone       = 400.0;  // Hall wall thickness 1 m
  G4double hallrad     = 700.0;  // Hall cylinder diam 16 m
  G4double hallhheight = 650.0;  // Hall cylinder height 13 m
  // water tank
  G4double tankwalltop = 0.6;  // water tank thickness at top 6 mm
  G4double tankwallbot = 0.8;  // water tank thickness at bottom 8 mm
  G4double tankrad     = 500;  // water tank diam 10 m
  G4double tankhheight = 400;  // water tank height 8 m
  // cryostat
  G4double cryowall   = 1.5;    // cryostat wall thickness from GERDA
  G4double vacgap     = 1.0;    // vacuum gap between walls
  G4double cryrad     = 200.0;  // cryostat diam 4 m
  G4double cryhheight = 225.0;  // cryostat height 4.5 m
  // Ge cylinder for 35.6 kg at 5.32 g/cm3
  G4double roiradius      = 15.0;                     // detector region diam 30 cm
  G4double roihalfheight  = 20.0;                     // detector region height 40 cm
  G4double gerad          = 3.7;                      // BEGe radius
  G4double gehheight      = 1.5;                      // full height 3 cm
  G4double begegap        = 2.0;                      // gap between BEGe 2cm
  G4double layerthickness = begegap + 2 * gehheight;  // 5 cm total
  G4int    nofLayers      = 8;

  fvertexZ = (hallhheight + offset) * cm;
  fmaxrad  = hallrad * cm;  // 8 m radius

  // Volumes for this geometry

  //
  // World
  //
  auto* worldSolid =
    new G4Tubs("World", 0.0 * cm, (hallrad + stone + 0.1) * cm,
               (hallhheight + stone + offset + 0.1) * cm, 0.0, CLHEP::twopi);
  auto* fWorldLogical  = new G4LogicalVolume(worldSolid, worldMaterial, "World_log");
  auto* fWorldPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fWorldLogical,
                                           "World_phys", nullptr, false, 0);

  //
  // Cavern
  //
  auto* cavernSolid    = new G4Tubs("Cavern", 0.0 * cm, (hallrad + stone) * cm,
                                    (hallhheight + stone) * cm, 0.0, CLHEP::twopi);
  auto* fCavernLogical = new G4LogicalVolume(cavernSolid, stdRock, "Cavern_log");
  auto* fCavernPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., offset * cm), fCavernLogical,
                      "Cavern_phys", fWorldLogical, false, 0);

  //
  // Hall
  //
  auto* hallSolid =
    new G4Tubs("Hall", 0.0 * cm, hallrad * cm, hallhheight * cm, 0.0, CLHEP::twopi);
  auto* fHallLogical  = new G4LogicalVolume(hallSolid, airMat, "Hall_log");
  auto* fHallPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fHallLogical,
                                          "Hall_phys", fCavernLogical, false, 0, true);

  //
  // Tank
  //
  auto* tankSolid =
    new G4Cons("Tank", 0.0 * cm, (tankrad + tankwallbot) * cm, 0.0 * cm,
               (tankrad + tankwalltop) * cm, tankhheight * cm, 0.0, CLHEP::twopi);
  auto* fTankLogical = new G4LogicalVolume(tankSolid, steelMat, "Tank_log");
  auto* fTankPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -offset * cm), fTankLogical,
                      "Tank_phys", fHallLogical, false, 0, true);

  //
  // Water
  //
  auto* waterSolid     = new G4Tubs("Water", 0.0 * cm, tankrad * cm,
                                    (tankhheight - tankwallbot) * cm, 0.0, CLHEP::twopi);
  auto* fWaterLogical  = new G4LogicalVolume(waterSolid, waterMat, "Water_log");
  auto* fWaterPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fWaterLogical,
                                           "Water_phys", fTankLogical, false, 0, true);

  //
  // outer cryostat
  //
  auto* coutSolid =
    new G4Tubs("Cout", 0.0 * cm, cryrad * cm, cryhheight * cm, 0.0, CLHEP::twopi);
  auto* fCoutLogical  = new G4LogicalVolume(coutSolid, steelMat, "Cout_log");
  auto* fCoutPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fCoutLogical,
                                          "Cout_phys", fWaterLogical, false, 0, true);

  //
  // vacuum gap
  //
  auto* cvacSolid     = new G4Tubs("Cvac", 0.0 * cm, (cryrad - cryowall) * cm,
                                   (cryhheight - cryowall) * cm, 0.0, CLHEP::twopi);
  auto* fCvacLogical  = new G4LogicalVolume(cvacSolid, worldMaterial, "Cvac_log");
  auto* fCvacPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fCvacLogical,
                                          "Cvac_phys", fCoutLogical, false, 0, true);

  //
  // inner cryostat
  //
  auto* cinnSolid     = new G4Tubs("Cinn", 0.0 * cm, (cryrad - cryowall - vacgap) * cm,
                                   (cryhheight - cryowall - vacgap) * cm, 0.0, CLHEP::twopi);
  auto* fCinnLogical  = new G4LogicalVolume(cinnSolid, steelMat, "Cinn_log");
  auto* fCinnPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fCinnLogical,
                                          "Cinn_phys", fCvacLogical, false, 0, true);

  //
  // LAr bath
  //
  auto* larSolid     = new G4Tubs("LAr", 0.0 * cm, (cryrad - 2 * cryowall - vacgap) * cm,
                                  (cryhheight - 2 * cryowall - vacgap) * cm, 0.0, CLHEP::twopi);
  auto* fLarLogical  = new G4LogicalVolume(larSolid, larMat, "Lar_log");
  auto* fLarPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fLarLogical,
                                         "Lar_phys", fCinnLogical, false, 0, true);

/*
  //
  // cryostat Lid
  //
  auto* lidSolid =
    new G4Tubs("Lid", 0.0 * cm, cryrad * cm, cryowall / 2.0 * cm, 0.0, CLHEP::twopi);
  auto* fLidLogical = new G4LogicalVolume(lidSolid, steelMat, "Lid_log");
  auto* fLidPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., (cryhheight + cryowall / 2.0) * cm),
                      fLidLogical, "Lid_phys", fWaterLogical, false, 0, true);
  auto* fBotLogical = new G4LogicalVolume(lidSolid, steelMat, "Bot_log");
  auto* fBotPhysical =
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -(cryhheight + cryowall / 2.0) * cm),
                      fBotLogical, "Bot_phys", fWaterLogical, false, 0, true);
*/

  //
  // String tower
  //
  // auto* towerSolid =
  //   new G4Tubs("String", 0.0 * cm, gerad * cm, roihalfheight * cm, 0.0, CLHEP::twopi);

  // auto* fTowerLogical  = new G4LogicalVolume(towerSolid, larMat, "Tower_log");

  // layers in tower
  auto* layerSolid = new G4Tubs("LayerSolid", 0.0 * cm, gerad * cm,
                                (gehheight + begegap / 2.0) * cm, 0.0, CLHEP::twopi);

  auto* fLayerLogical = new G4LogicalVolume(layerSolid, larMat, "Layer_log");

  // fill one layer
  auto* geSolid =
    new G4Tubs("ROI", 0.0 * cm, gerad * cm, gehheight * cm, 0.0, CLHEP::twopi);

  auto* fGeLogical = new G4LogicalVolume(geSolid, roiMat, "Ge_log");
  new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, -begegap / 2.0 * cm), fGeLogical,
                    "Ge_phys", fLayerLogical, false, 0, true);

  auto* gapSolid =
    new G4Tubs("Gap", 0.0 * cm, gerad * cm, begegap / 2.0 * cm, 0.0, CLHEP::twopi);

  auto* fGapLogical = new G4LogicalVolume(gapSolid, larMat, "Gap_log");
  new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, gehheight * cm), fGapLogical,
                    "Gap_phys", fLayerLogical, false, 0, true);

  // place layers as mother volume with unique copy number
  G4double step = (gehheight + begegap / 2) * cm;
  G4double xpos;
  G4double ypos;
  G4double angle = CLHEP::twopi / 6.0;

  if(fGeometryName != "hallA_wo_ge" && fGeometryName != "hallA_only_WLSR"){
    for(G4int j = 0; j < 6; j++)
    {
      xpos = roiradius * cm * std::cos(j * angle);
      ypos = roiradius * cm * std::sin(j * angle);
      for(G4int i = 0; i < nofLayers; i++)
      {
        new G4PVPlacement(
          nullptr,
          G4ThreeVector(xpos, ypos,
                        -step + (nofLayers / 2 * layerthickness - i * layerthickness) * cm),
          fLayerLogical, "Layer_phys", fLarLogical, false, i + j * nofLayers, true);
      }
    }
  }

  // WLSR volume - its a workaround to access the volume of the WLSR w/o implementing it properly
  if(fGeometryName == "hallA_only_WLSR"){
    // cavern
    G4double WLSR_r       = 70.0;  
    G4double WLSR_h       = 300.0/2.; 

    auto WLSR_LAr_solid = new G4Tubs("WLSR_LAr_solid", 0.0 * cm, WLSR_r * cm, WLSR_h * cm, 0.0, CLHEP::twopi);  
    auto* fWLSR_LAr_logical = new G4LogicalVolume(WLSR_LAr_solid, larMat, "WLSR_LAr_logical");
    new G4PVPlacement(nullptr, G4ThreeVector(), fWLSR_LAr_logical,
                    "WLSR_LAr_physical", fLarLogical, false, 0, true);
  }

  //
  // Visualization attributes
  //
  //fWorldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());

  auto* redVisAtt = new G4VisAttributes(G4Colour::Red());
  redVisAtt->SetVisibility(true);
  auto* greyVisAtt = new G4VisAttributes(G4Colour::Grey());
  greyVisAtt->SetVisibility(true);
  auto* greenVisAtt = new G4VisAttributes(G4Colour::Green());
  greenVisAtt->SetVisibility(true);
  auto* blueVisAtt = new G4VisAttributes(G4Colour::Blue());
  blueVisAtt->SetVisibility(true);

  fCavernLogical->SetVisAttributes(redVisAtt);
  fHallLogical->SetVisAttributes(greyVisAtt);
  fTankLogical->SetVisAttributes(greenVisAtt);
  fWaterLogical->SetVisAttributes(greyVisAtt);
  fLarLogical->SetVisAttributes(greyVisAtt);
  fCoutLogical->SetVisAttributes(blueVisAtt);
  fCvacLogical->SetVisAttributes(greyVisAtt);
  fCinnLogical->SetVisAttributes(blueVisAtt);
  //fLidLogical->SetVisAttributes(blueVisAtt);
  //fBotLogical->SetVisAttributes(blueVisAtt);
  //fLayerLogical->SetVisAttributes(G4VisAttributes::Invisible);
  fGapLogical->SetVisAttributes(greyVisAtt);
  fGeLogical->SetVisAttributes(redVisAtt);
  return fWorldPhysical;
}
void WLGDDetectorConstruction::SetPositionOfDetectors(G4String name)
{
  std::set<G4String> knownGeometries = { "baseline", "original" };
  if(knownGeometries.count(name) == 0)
  {
    G4Exception("WLGDDetectorConstruction::SetGeometry", "WLGD0001", JustWarning,
                ("Invalid geometry setup name '" + name + "'").c_str());
    return;
  }

  fDetectorPosition = name;
  // Reinit wiping out stores
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void WLGDDetectorConstruction::SetGeometry(const G4String& name)
{
  std::set<G4String> knownGeometries = { "baseline", "baseline_smaller", "baseline_large_reentrance_tube", "alternative",
                                         "hallA", "hallA_wo_ge", "hallA_only_WLSR", "baseline_large_reentrance_tube_4m_cryo" };
  if(knownGeometries.count(name) == 0)
  {
    G4Exception("WLGDDetectorConstruction::SetGeometry", "WLGD0001", JustWarning,
                ("Invalid geometry setup name '" + name + "'").c_str());
    return;
  }

  fGeometryName = name;
  // Reinit wiping out stores
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void WLGDDetectorConstruction::ExportGeometry(const G4String& file)
{
  G4GDMLParser parser;
  parser.Write(file);
}

void WLGDDetectorConstruction::SetNeutronBiasFactor(G4double nf) { fNeutronBias = nf; }

void WLGDDetectorConstruction::SetMuonBiasFactor(G4double mf) { fMuonBias = mf; }

void WLGDDetectorConstruction::SetNeutronYieldBias(G4double nf)
{
  fNeutronYieldBias = nf;
}

// Additional settings for adjusting the detector geometry

// changing the concentration of Xe in the LAr
void WLGDDetectorConstruction::SetXeConc(G4double nf)
{
  fXeConc = nf * 1e-3;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

// changing the concentration of He-3 in the LAr
void WLGDDetectorConstruction::SetHe3Conc(G4double nf)
{
  fHe3Conc = nf * 1e-3;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

// changing the radius if the cryostat
void WLGDDetectorConstruction::SetOuterCryostatRadius(G4double rad)
{
  fCryostatOuterRadius = rad;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

// changing the height of the cryostat
void WLGDDetectorConstruction::SetCryostatHeight(G4double height)
{
  fCryostatHeight = height;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

// option to remove the cupper tubes to see what impact it has on the production rate
void WLGDDetectorConstruction::SetWithoutCupperTubes(G4int answer)
{
  fWithOutCupperTubes = answer;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

// option to include borated PE in the setup (1: tubes around the re-entrance tubes, 2:
// trubine structure)
void WLGDDetectorConstruction::SetNeutronModerator(G4int answer)
{
  fWithBoratedPET = answer;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

// option to include borated PE in the setup (1: tubes around the re-entrance tubes, 2:
// trubine structure)
void WLGDDetectorConstruction::SetMaterial(G4String answer)
{
  fSetMaterial = answer;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

// option to set the radius of the turbine structure
void WLGDDetectorConstruction::SetTurbineAndTubeRadius(G4double radius)
{
  fBoratedTurbineRadius = radius;
}

// option to set the length of the turbine structure
void WLGDDetectorConstruction::SetTurbineAndTubeLength(G4double length)
{
  fBoratedTurbineLength = length;
}

// option to set the width of the turbine structure
void WLGDDetectorConstruction::SetTurbineAndTubeWidth(G4double width)
{
  fBoratedTurbineWidth = width;
}

// option to set the angle of the turbine structure
void WLGDDetectorConstruction::SetTurbineAndTubeAngle(G4double deg)
{
  fBoratedTurbineAngle = deg;
}

// option to set the height of the turbine structure
void WLGDDetectorConstruction::SetTurbineAndTubeHeight(G4double height)
{
  fBoratedTurbineHeight = height;
}

// option to set the zPosition of the turbine structure
void WLGDDetectorConstruction::SetTurbineAndTubezPosition(G4double zPosition)
{
  fBoratedTurbinezPosition = zPosition;
}

// option to set the zPosition of the turbine structure
void WLGDDetectorConstruction::SetTurbineAndTubeNPanels(G4double nPanels)
{
  fBoratedTurbineNPanels = nPanels;
}

// option to change the water to gadolinium weighted water in the water tank
void WLGDDetectorConstruction::SetGdWater(G4int answer)
{
  fWithGdWater = answer;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void WLGDDetectorConstruction::SetWoWater(G4int answer)
{
  fWithWoWater = answer;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void WLGDDetectorConstruction::SetMaGeMaterial(G4int answer)
{
  fMaGeMaterial = answer;
  WLGDDetectorConstruction::DefineMaterials();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void WLGDDetectorConstruction::DefineCommands()
{
  // Define geometry command directory using generic messenger class
  fDetectorMessenger = new G4GenericMessenger(this, "/WLGD/detector/",
                                              "Commands for controlling detector setup");

  // switch command
  fDetectorMessenger
    ->DeclareMethod("setPositionOfDetectors",
                    &WLGDDetectorConstruction::SetPositionOfDetectors)
    .SetGuidance("Set geometry model of cavern and detector")
    .SetGuidance("original = original Warwick-Legend Position")
    .SetGuidance("baseline = lowered position")
    .SetCandidates("original baseline")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);

  // switch command
  fDetectorMessenger->DeclareMethod("setGeometry", &WLGDDetectorConstruction::SetGeometry)
    .SetGuidance("Set geometry model of cavern and detector")
    .SetGuidance("baseline = NEEDS DESCRIPTION")
    .SetGuidance("baseline_smaller = Gerda cryostat with only one module")
    .SetGuidance("baseline_large_reentrance_tube = large single re-entrance tube")
    .SetGuidance("baseline_large_reentrance_tube_4m_cryo")
    .SetGuidance("alternative = NEEDS DESCRIPTION")
    .SetGuidance("hallA = Gerda mock-up for validation.")
    .SetGuidance("hallA_wo_ge = Gerda mock-up w/o HPGe")
    .SetGuidance("hallA_only_WLSR = Gerda w/o detectors but with a workaround WLSR volume")
    .SetCandidates("baseline baseline_smaller baseline_large_reentrance_tube alternative hallA hallA_wo_ge hallA_only_WLSR baseline_large_reentrance_tube_4m_cryo")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);

  // GDML Export
  fDetectorMessenger
    ->DeclareMethod("exportGeometry", &WLGDDetectorConstruction::ExportGeometry)
    .SetGuidance("Export current geometry to a GDML file")
    .SetParameterName("filename", false)
    .SetDefaultValue("wlgd.gdml")
    .SetStates(G4State_Idle)
    .SetToBeBroadcasted(false);

  // Edit: 2021/02/17 by Moritz Neuberger
  // Adding options to adjust the concentration of Xe and He-3 in the LAr to test change
  // in Ge-77 production

  // changing the concentration of Xe in the LAr
  fDetectorMessenger->DeclareMethod("XeConc", &WLGDDetectorConstruction::SetXeConc)
    .SetGuidance("Set concentration of Xe in the LAr [mg/g]")
    .SetDefaultValue("0.0")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);

  // changing the concentration of He-3 in the LAr
  fDetectorMessenger->DeclareMethod("He3Conc", &WLGDDetectorConstruction::SetHe3Conc)
    .SetGuidance("Set concentration of He3 in the LAr [mg/g]")
    .SetDefaultValue("0.0")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);

  // Edit: 2021/03/30 by Moritz Neuberger
  // Mainly for the baseline geometry

  // changing the radius if the cryostat
  fDetectorMessenger
    ->DeclareMethod("Cryostat_Radius_Outer",
                    &WLGDDetectorConstruction::SetOuterCryostatRadius)
    .SetGuidance("Set the outer radius of the cryostat [cm]")
    .SetDefaultValue("350.0")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);

  // changing the height of the cryostat
  fDetectorMessenger
    ->DeclareMethod("Cryostat_Height", &WLGDDetectorConstruction::SetCryostatHeight)
    .SetGuidance("Set the height of the cryostat [cm]")
    .SetDefaultValue("350.0")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);

  // option to remove the cupper tubes to see what impact it has on the production rate
  fDetectorMessenger
    ->DeclareMethod("Without_Cupper_Tubes",
                    &WLGDDetectorConstruction::SetWithoutCupperTubes)
    .SetGuidance("Set whether to include cupper tubes or not")
    .SetGuidance("0 = with cupper tubes")
    .SetGuidance("1 = without cupper tubes")
    .SetCandidates("0 1")
    .SetDefaultValue("0");

  // option to include borated PE in the setup (1: tubes around the re-entrance tubes, 2:
  // trubine structure)
  fDetectorMessenger
    ->DeclareMethod("With_NeutronModerators",
                    &WLGDDetectorConstruction::SetNeutronModerator)
    .SetGuidance("Set whether to include Neutron Moderators or not")
    .SetGuidance("0 = without Neutron Moderators")
    .SetGuidance("1 = with Neutron Moderators around Re-Entrance tubes")
    .SetGuidance("2 = with Neutron Moderators in turbine mode")
    .SetGuidance("3 = with Neutron Moderators in large tub")
    .SetGuidance("4 = with Neutron Moderators in turbine mode with lids")
    .SetCandidates("0 1 2 3 4")
    .SetDefaultValue("0");

  // option to include borated PE in the setup (1: tubes around the re-entrance tubes, 2:
  // trubine structure)
  fDetectorMessenger
    ->DeclareMethod("Which_Material", &WLGDDetectorConstruction::SetMaterial)
    .SetGuidance("Set which material should be used instead of Boarated PE")
    .SetGuidance("BoratedPE = normal case")
    .SetGuidance("PolyEthylene = without Boron")
    .SetGuidance("PMMA = instead using PMMA")
    .SetCandidates("BoratedPE PolyEthylene PMMA")
    .SetDefaultValue("BoratedPE");

  // option to set the radius of the turbine structure
  fDetectorMessenger
    ->DeclareMethod("TurbineAndTube_Radius",
                    &WLGDDetectorConstruction::SetTurbineAndTubeRadius)
    .SetGuidance("Set the radius on which the borated PE pannels are aligned on [cm]")
    .SetDefaultValue("200.0")
    .SetToBeBroadcasted(false);

  // option to set the radius of the turbine structure
  fDetectorMessenger
    ->DeclareMethod("TurbineAndTube_Length",
                    &WLGDDetectorConstruction::SetTurbineAndTubeLength)
    .SetGuidance("Set the length of a panel of borated PE [cm]")
    .SetDefaultValue("50.0")
    .SetToBeBroadcasted(false);

  // option to set the radius of the turbine structure
  fDetectorMessenger
    ->DeclareMethod("TurbineAndTube_Angle",
                    &WLGDDetectorConstruction::SetTurbineAndTubeAngle)
    .SetGuidance("Set the angle on which the borated PE pannels are aligned on [deg]")
    .SetDefaultValue("45.0")
    .SetToBeBroadcasted(false);

  // option to set the radius of the turbine structure
  fDetectorMessenger
    ->DeclareMethod("TurbineAndTube_Width",
                    &WLGDDetectorConstruction::SetTurbineAndTubeWidth)
    .SetGuidance("Set the width of the borated PE pannels [cm]")
    .SetDefaultValue("5.0")
    .SetToBeBroadcasted(false);

  // option to set the radius of the turbine structure
  fDetectorMessenger
    ->DeclareMethod("TurbineAndTube_Height",
                    &WLGDDetectorConstruction::SetTurbineAndTubeHeight)
    .SetGuidance("Set the height of the borated PE pannels [cm] (this is total height)")
    .SetDefaultValue("400")
    .SetToBeBroadcasted(false);

  // option to set the radius of the turbine structure
  fDetectorMessenger
    ->DeclareMethod("TurbineAndTube_zPosition",
                    &WLGDDetectorConstruction::SetTurbineAndTubezPosition)
    .SetGuidance("Set the zPosition of the borated PE pannels [cm]")
    .SetDefaultValue("0")
    .SetToBeBroadcasted(false);

  // option to set the number of panels of the turbine structure
  fDetectorMessenger
    ->DeclareMethod("TurbineAndTube_NPanels",
                    &WLGDDetectorConstruction::SetTurbineAndTubeNPanels)
    .SetGuidance("Set the number of panels of the borated PE pannels [cm]")
    .SetDefaultValue("0")
    .SetToBeBroadcasted(false);

  // option to change the water to gadolinium weighted water in the water tank
  fDetectorMessenger
    ->DeclareMethod("With_Gd_Water", &WLGDDetectorConstruction::SetGdWater)
    .SetGuidance("Set whether to include Gd water or not")
    .SetGuidance("0 = without Gd water")
    .SetGuidance("1 = with Gd water")
    .SetCandidates("0 1")
    .SetDefaultValue("0");
  // option to change the water to gadolinium weighted water in the water tank
  fDetectorMessenger
    ->DeclareMethod("Without_Water", &WLGDDetectorConstruction::SetWoWater)
    .SetGuidance("Set whether to include water or not")
    .SetGuidance("0 = with water")
    .SetGuidance("1 = without water")
    .SetCandidates("0 1")
    .SetDefaultValue("0");


  // Define bias operator command directory using generic messenger class
  fBiasMessenger =
    new G4GenericMessenger(this, "/WLGD/bias/", "Commands for controlling bias factors");

  // switch commands
  fBiasMessenger
    ->DeclareMethod("setNeutronBias", &WLGDDetectorConstruction::SetNeutronBiasFactor)
    .SetGuidance("Set Bias factor for neutron capture process.")
    .SetDefaultValue("1.0")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);
  fBiasMessenger
    ->DeclareMethod("setMuonBias", &WLGDDetectorConstruction::SetMuonBiasFactor)
    .SetGuidance("Set Bias factor for muon nuclear process.")
    .SetDefaultValue("1.0")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);
  fBiasMessenger
    ->DeclareMethod("setNeutronYieldBias", &WLGDDetectorConstruction::SetNeutronYieldBias)
    .SetGuidance("Set Bias factor for neutron yield process.")
    .SetDefaultValue("1.0")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);

  // Define bias operator command directory using generic messenger class
  fMaterialMessenger =
    new G4GenericMessenger(this, "/WLGD/material/", "Commands for controlling the material choice");  
    
  fMaterialMessenger
    ->DeclareMethod("setMaGeMaterial", &WLGDDetectorConstruction::SetMaGeMaterial)
    .SetGuidance("Set LAr and Ge material to MaGe settings")
    .SetDefaultValue("0")
    .SetStates(G4State_PreInit)
    .SetToBeBroadcasted(false);
}
