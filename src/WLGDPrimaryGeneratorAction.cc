// us
#include "WLGDPrimaryGeneratorAction.hh"
#include "WLGDPetersGammaCascadeReader.hh"
#include "WLGDDetectorConstruction.hh"

// geant
#include "G4Event.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// std
#include <fstream>
#include <random>
#include <set>
/*#include "TH1F.h"
#include "TFile.h"*/
//#include "TH1.h"

// G4String WLGDPrimaryGeneratorAction::fFileName;
// std::ifstream* WLGDPrimaryGeneratorAction::fInputFile;

vector<int> WLGDPrimaryGeneratorAction::WLGD_get_coord_from_index(int index) {
    int n = 0;
    for (int i = -4; i <= 4; i++) {
        for (int j = 0; j <= std::min(4, 4 - i); j++) {
            if ((abs(i) == 3 || i == 0) && (j == 0 || j == 3)) {
                continue;
            }
            for (int k = 0; k < 7; k++) {
                if (n == index) {
                  std::vector<int> output{ i, j, -(i + j), k }; 
                  return output;
                } 
                n++;
            }
        }
        for (int j = 1; j <= std::min(4, 4 - i); j++) {
            if ((abs(i) == 3 || i == 0) && (j == 0 || j == 3)) {
                continue;
            }
            for (int k = 0; k < 7; k++) {
                if (n == index) {
                  std::vector<int> output{ i, j, -(i + j), k }; 
                  return output;
                } 
                n++;
            }
        }
    }
    std::vector<int> output{ 0,0,0 }; 
    return output;
}

vector<double> WLGDPrimaryGeneratorAction::WLGD_get_carthesian_coordinates(int index){
  vector<int> p = WLGD_get_coord_from_index(index);

  double offset = -20 * cm;
  double length = 16.75 * cm;
  double vec_main_x = 1./2.;
  double vec_main_y = std::sqrt(1. - vec_main_x*vec_main_x);
  double vec_left_x = -1./2.;
  double vec_left_y = std::sqrt(1. - vec_left_x*vec_left_x);
  double vec_right_x = 1;
  double vec_right_y = 0;
  int nofLayers = 7;
  double gehheight = 5.0 * cm;
  double gegap = 3.0 * cm;
  double layerthickness = gegap + 2 * gehheight;
  double step = (gehheight + gegap / 2);

  std::vector<double> output(3, 0.);
  if (p[1] >= 0) {
      output[0] = length * ( vec_main_x * p[0] + vec_left_x * p[1]);
      output[1] = length * ( vec_main_y * p[0] + vec_left_y * p[1]);
      output[2] = offset - p[3] * layerthickness;
  }
  if (p[1] < 0) {
      output[0] = length * ( vec_main_x * p[0] - vec_right_x * p[1]);
      output[1] = length * ( vec_main_y * p[0] - vec_right_y * p[1]);
      output[2] = offset - p[3] * layerthickness;
  }
  return output;
}

vector<double> WLGDPrimaryGeneratorAction::WLGD_generate_random_position_in_detectors(){
  double gehheight = 5.0 * cm;
  double gerad     = 4.0 * cm; 
  std::uniform_real_distribution<>   rndm(0.0, 1.0);  
  double z = gehheight * (1. - 2. * rndm(generator));
  double r = gerad * rndm(generator);
  double phi = CLHEP::twopi * rndm(generator);

  double x = r * std::cos(phi);
  double y = r * std::sin(phi);

  vector<double> output{x,y,z};
  return output;
}



WLGDPrimaryGeneratorAction::WLGDPrimaryGeneratorAction(WLGDDetectorConstruction* det)
: G4VUserPrimaryGeneratorAction()
, fDetector(det)
, fParticleGun(nullptr)
, fMessenger(nullptr)
, fDepth(0.0)
, fGenerator("Musun")
, fZShift(200.0 * cm)
{
  generator.seed(rd());  // set a random seed

  G4int nofParticles = 1;
  fParticleGun       = new G4ParticleGun(nofParticles);

  auto particleTable = G4ParticleTable::GetParticleTable();

  // default particle kinematics
  fParticleGun->SetParticleDefinition(particleTable->FindParticle("mu-"));

  // define commands for this class
  DefineCommands();
  fFileName = "";
}

WLGDPrimaryGeneratorAction::~WLGDPrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fMessenger;
  if(fInputFile.is_open())
    fInputFile.close();
}

// -- for the Musun method, input files have to be provided

void WLGDPrimaryGeneratorAction::OpenFile()
{
  fInputFile.open(fFileName, std::ifstream::in);
  if(!(fInputFile.is_open()))
  {  
     G4cerr << "Musung file not valid! Name: " << fFileName << G4endl;
  }
}


void WLGDPrimaryGeneratorAction::ChangeFileName(G4String newFile)
{
  if(fFileName != newFile)  // check if the new file is equal to the other
  {
    if(fInputFile.is_open())
      fInputFile.close();  // close the old file
    G4cout << "opening file: " << newFile << G4endl;
    fFileName = newFile;
    OpenFile();  // open the new one
  }
}

// -- depending on the name of the generator given, a different method is used to generate the primaries

void WLGDPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  if(fGenerator == "SimpleNeutronGun"){
    G4ThreeVector momentumDir(1, 0, 0);
    fParticleGun->SetParticleMomentumDirection(momentumDir);
    fParticleGun->SetParticleEnergy(neutron_ekin * eV);
    fParticleGun->SetParticlePosition(G4ThreeVector(coord_x*cm, coord_y*cm, coord_z*cm));
    auto particleTable = G4ParticleTable::GetParticleTable();
    fParticleGun->SetParticleDefinition(particleTable->FindParticle("neutron"));
    fParticleGun->GeneratePrimaryVertex(event);
  }
  if(fGenerator == "MeiAndHume")
  {
    using pld_type = std::piecewise_linear_distribution<double>;

    int    nw             = 100;     // number of bins
    double lower_bound    = 1.0;     // energy interval lower bound [GeV]
    double upper_bound    = 3000.0;  // upper bound [GeV]
    double nearhorizontal = 1.0e-5;
    double fullcosangle   = 1.0;

    // custom probability distributions
    pld_type ed(nw, lower_bound, upper_bound, MuEnergy(fDepth));
    pld_type cosd(nw, nearhorizontal, fullcosangle, MuAngle(fDepth));

    // momentum vector
    G4double costheta = cosd(generator);  // get a random number
    G4double sintheta = std::sqrt(1. - costheta * costheta);

    std::uniform_real_distribution<> rndm(0.0, 1.0);   // azimuth angle
    G4double phi    = CLHEP::twopi * rndm(generator);  // random uniform number
    G4double sinphi = std::sin(phi);
    G4double cosphi = std::cos(phi);

    G4double      px = -sintheta * cosphi;
    G4double      py = -sintheta * sinphi;
    G4double      pz = -costheta;  // default downwards: pz = -1.0
    G4ThreeVector momentumDir(px, py, pz);
    fParticleGun->SetParticleMomentumDirection(momentumDir);
    // G4cout << "Momentum direction Primary: " << momentumDir << G4endl;

    G4double ekin = ed(generator);  // get random number
    ekin *= GeV;
    fParticleGun->SetParticleEnergy(ekin);

    // position, top of world, sample circle uniformly
    G4double zvertex = fDetector->GetWorldSizeZ();  // inline on WLGDDetectorConstruction
    G4double radius  = fDetector->GetWorldExtent() * rndm(generator);  // fraction of max
    phi              = CLHEP::twopi * rndm(generator);  // another random angle
    G4double vx      = radius * std::cos(phi);
    G4double vy      = radius * std::sin(phi);

    fParticleGun->SetParticlePosition(G4ThreeVector(vx, vy, zvertex - 1.0 * cm));

    fParticleGun->GeneratePrimaryVertex(event);
  }
  if(fGenerator == "Musun")
  {
    G4int    nEvent = 0;
    G4double time   = 0.0;
    G4double energy = 0.0 * MeV;
    G4double px, py, pz;
    G4double theta, phi;
    G4double x = 0, y = 0, z = 0;
    G4int    particleID = 0;

    fInputFile >> nEvent >> particleID >> energy >> x >> y >> z >> theta >> phi;

     //G4cout  << nEvent << " " << x << " " << y << " " << z << G4endl;
    if(fInputFile.eof())
    {
      fInputFile.close();
      G4cerr << "File over: not enough events! Debugoutput" << G4endl;
      G4Exception("WLGDPrimaryGeneratorAction::GeneratePrimaryVertex()", "err001",
                  FatalException, "Exit Warwick");
      return;
    }

    G4double particle_time = time * s;
    energy                 = energy * GeV;
    theta                  = theta * rad;
    phi                    = phi * rad;
    x                      = x * cm;
    y                      = y * cm;
    z                      = fZShift + (z * cm);

    //   G4cout << "Primary coordinates: " << position/m << " m" << G4endl;
    //   G4cout << "Primary coordinates: " << x/cm << " " <<  y/cm << " " << z/cm << " "
    //   << G4endl; G4cout << "Primary energy: " << energy/GeV << " GeV" << G4endl; G4cout
    //   << "Theta: " << theta/deg << " deg; Phi: " << phi/deg << " deg" << G4endl;

    G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();

    G4String particleName = " ";

    if(particleID == 10)
      particleName = "mu+";
    else
      particleName = "mu-";

    G4double theMass     = theParticleTable->FindParticle(particleName)->GetPDGMass();
    G4double totMomentum = std::sqrt(energy * energy + 2 * theMass * energy);
    pz                   = -1 * std::cos(theta);
    px                   = std::sin(theta) * cos(phi);
    py                   = std::sin(theta) * sin(phi);
    G4ThreeVector momentumDir(px, py, pz);

    fParticleGun->SetParticleMomentumDirection(momentumDir);

    fParticleGun->SetParticleEnergy(energy);

    fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));

    fParticleGun->GeneratePrimaryVertex(event);
  }

  if(fGenerator == "Musun_alternative")
  {
    G4int    nEvent = 0;
    G4double time   = 0.0;
    G4double energy = 0.0 * MeV;
    G4double px, py, pz;
    G4double theta, phi;
    G4double x = 0, y = 0, z = 0;
    G4int    particleID = 0;

    fInputFile >> nEvent >> particleID >> energy >> x >> y >> z >> px >> py >> pz; // >> theta >> phi;

     //G4cout  << nEvent << " " << x << " " << y << " " << z << G4endl;
    if(fInputFile.eof())
    {
      fInputFile.close();
      G4cerr << "File over: not enough events! Debugoutput" << G4endl;
      G4Exception("WLGDPrimaryGeneratorAction::GeneratePrimaryVertex()", "err001",
                  FatalException, "Exit Warwick");
      return;
    }

    G4double particle_time = time * s;
    energy                 = energy * GeV;
    theta                  = theta * rad;
    phi                    = phi * rad;
    x                      = x * cm;
    y                      = y * cm;
    z                      = fZShift + (z * cm);

    //   G4cout << "Primary coordinates: " << position/m << " m" << G4endl;
    //   G4cout << "Primary coordinates: " << x/cm << " " <<  y/cm << " " << z/cm << " "
    //   << G4endl; G4cout << "Primary energy: " << energy/GeV << " GeV" << G4endl; G4cout
    //   << "Theta: " << theta/deg << " deg; Phi: " << phi/deg << " deg" << G4endl;

    G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();

    G4String particleName = " ";

    if(particleID == 10)
      particleName = "mu+";
    else
      particleName = "mu-";

    G4double theMass     = theParticleTable->FindParticle(particleName)->GetPDGMass();
    G4double totMomentum = std::sqrt(energy * energy + 2 * theMass * energy);
    //pz                   = -1 * std::cos(theta);
    //px                   = std::sin(theta) * cos(phi);
    //py                   = std::sin(theta) * sin(phi);
    G4ThreeVector momentumDir(px, py, pz);

    fParticleGun->SetParticleMomentumDirection(momentumDir);

    fParticleGun->SetParticleEnergy(energy);

    fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));

    fParticleGun->GeneratePrimaryVertex(event);
  }
  if(fGenerator == "Ge77m" || fGenerator == "Ge77andGe77m" || fGenerator == "GammaCascade")
  {
    G4double cushift        = 150.;
    G4double roiradius      = 30.0;  // string radius curad - Ge radius - gap
    G4double gerad          = 4.0;   // Ge radius
    G4int    nofLayers      = 7;     // 8;   // 8 Ge + 7 gaps = 1010 mm string height
    G4int    nofStrings     = 14;    // 12 strings  of 8 Ge each
    G4double gehheight      = 5.0;   // full height 10 cm
    G4double gegap          = 3.0;   // gap between Ge 3cm
    G4double step           = (gehheight + gegap / 2) * cm;
    G4double layerthickness = gegap + 2 * gehheight;  // 13 cm total
    G4double angle          = CLHEP::twopi / nofStrings;

    std::uniform_int_distribution<int> distribution(0,
                                                    4 * (int) (nofLayers * nofStrings));
    std::uniform_real_distribution<>   rndm(0.0, 1.0);

    G4int    detectorNumber      = distribution(generator);
    G4int    whichReentranceTube = detectorNumber / ((int) (nofLayers * nofStrings));
    G4int    z_number            = detectorNumber % nofLayers;
    G4int    string_number = (detectorNumber % (nofLayers * nofStrings)) / nofLayers;
    G4double radius        = gerad * cm * rndm(generator);
    //G4double phi           = CLHEP::twopi * rndm(generator);

    double offset_x, offset_y;
    if(whichReentranceTube == 0)
    {
      offset_x = 1 * m;
      offset_y = 0 * m;
    }
    if(whichReentranceTube == 1)
    {
      offset_x = 0 * m;
      offset_y = 1 * m;
    }
    if(whichReentranceTube == 2)
    {
      offset_x = -1 * m;
      offset_y = 0 * m;
    }
    if(whichReentranceTube == 3)
    {
      offset_x = 0 * m;
      offset_y = -1 * m;
    }

    G4double particle_time = 0 * s;
    G4double energy        = 0 * GeV;
    G4double theta         = 0 * rad;
    G4double phi           = 0 * rad;

/*
    G4double x =
      radius * cos(phi) + roiradius * cm * std::cos(string_number * angle) + offset_x;
    G4double y =
      radius * sin(phi) + roiradius * cm * std::sin(string_number * angle) + offset_y;
    */
    double vec_main_x = 1/2.;
    double vec_main_y = sqrt(1 - vec_main_x*vec_main_x);
    double vec_left_x = -1/2.;
    double vec_left_y = sqrt(1 - vec_left_x*vec_left_x);
    double vec_right_x = 1;
    double vec_right_y = 0;

    std::uniform_int_distribution<>   rndm_detector_index(0,378);

    vector<double> coor = WLGD_get_carthesian_coordinates(rndm_detector_index(generator));
    vector<double> pos_in_det = WLGD_generate_random_position_in_detectors();

    G4double x = coor[0] + pos_in_det[0];
    G4double y = coor[1] + pos_in_det[1];
    G4double z = coor[2] + pos_in_det[2];
    
 /*   cushift * cm - step +
                 (nofLayers / 2 * layerthickness - z_number * layerthickness) * cm +
                 gehheight * cm * (1 - 2 * rndm(generator));*/

    //   G4cout << "Primary coordinates: " << position/m << " m" << G4endl;
    //   G4cout << "Primary coordinates: " << x/cm << " " <<  y/cm << " " << z/cm << " "
    //   << G4endl; G4cout << "Primary energy: " << energy/GeV << " GeV" << G4endl; G4cout
    //   << "Theta: " << theta/deg << " deg; Phi: " << phi/deg << " deg" << G4endl;

    G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
    G4double         theMass = theParticleTable->GetIonTable()->GetIonMass(32, 77, 0, 1);

    if(fGenerator != "GammaCascade")
    {
      if(fGenerator == "Ge77andGe77m")
      {
        std::uniform_int_distribution<int> distribution_2(0, 1);
        G4int                              whichGe77State = distribution_2(generator);
        if(whichGe77State == 0)
        {
          fParticleGun->SetParticleDefinition(
            theParticleTable->GetIonTable()->GetIon(32, 77, 0 * keV));
          theMass = theParticleTable->GetIonTable()->GetIonMass(32, 77, 0, 0);
        }
        else
        {
          fParticleGun->SetParticleDefinition(
            theParticleTable->GetIonTable()->GetIon(32, 77, 159.71 * keV));
          theMass = theParticleTable->GetIonTable()->GetIonMass(32, 77, 0, 1);
        }
      }

      if(fGenerator == "Ge77m")
      {
        fParticleGun->SetParticleDefinition(
          theParticleTable->GetIonTable()->GetIon(32, 77, 159.71 * keV));
        theMass = theParticleTable->GetIonTable()->GetIonMass(32, 77, 0, 1);
      }
      G4double      totMomentum = std::sqrt(energy * energy + 2 * theMass * energy);
      G4double      pz          = -1 * std::cos(theta);
      G4double      px          = std::sin(theta) * cos(phi);
      G4double      py          = std::sin(theta) * sin(phi);
      G4ThreeVector momentumDir(px, py, pz);

      fParticleGun->SetParticleMomentumDirection(momentumDir);

      fParticleGun->SetParticleEnergy(energy);

      fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));

      fParticleGun->GeneratePrimaryVertex(event);
    }
    else{
      fParticleGun->SetParticleDefinition(theParticleTable->FindParticle("gamma"));
      WLGDPetersGammaCascadeReader* pgcr = WLGDPetersGammaCascadeReader::GetInstance();
      GammaCascadeLine input = pgcr->GetNextEntry(pgc_neutron_ekin);
      for(int i = 0; i < input.m; i++){
        G4double theta = rndm(generator) * 180. * deg;
        G4double phi   = rndm(generator) * 360. * deg;

        G4double      pz          = -1 * std::cos(theta);
        G4double      px          = std::sin(theta) * cos(phi);
        G4double      py          = std::sin(theta) * sin(phi);
        G4ThreeVector momentumDir(px, py, pz);
        fParticleGun->SetParticleMomentumDirection(momentumDir);
        fParticleGun->SetParticleEnergy(input.eg[i]*keV);
        fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
        fParticleGun->GeneratePrimaryVertex(event);
      }
    }
  }
  if(fGenerator == "ModeratorNeutrons")
  {
    if(neutronEnergySpectrumInBPE == 0)
    {
      string   filelistName = "../data/resultingSpectrum.txt";
      ifstream filestream;
      filestream.open(filelistName.c_str());
      vector<double> x_val;
      vector<double> y_val;
      double         tmp_x, tmp_y;
      while(filestream >> tmp_x >> tmp_y)
      {
        x_val.push_back(tmp_x);
        y_val.push_back(tmp_y);
      }
      neutronEnergySpectrumInBPE = new std::piecewise_linear_distribution<double>(
        x_val.begin(), x_val.end(), y_val.begin());
    }

    G4int type = fDetector->GetBoratedType();
    if(type == 0)
      throw std::runtime_error(
        std::string("Do not use BoratedPENeutrons generator without using Neutron Moderators! ):"));

    std::uniform_int_distribution<int> distribution(0, 4);
    std::uniform_real_distribution<>   rndm(0.0, 1.0);

    G4double ran_x, ran_y, ran_z;

    // - depending on the different types of moderator design
    if(type == 1)
    {
      G4double curad              = 40.0;
      G4double BoratedPETouterrad = 5.0;
      G4double cuhheight          = 400.0 / 2.;

      G4int    whichReentranceTube = distribution(generator);
      G4double offset_x, offset_y;
      if(whichReentranceTube == 0)
      {
        offset_x = 1 * m;
        offset_y = 0 * m;
      }
      if(whichReentranceTube == 1)
      {
        offset_x = 0 * m;
        offset_y = 1 * m;
      }
      if(whichReentranceTube == 2)
      {
        offset_x = -1 * m;
        offset_y = 0 * m;
      }
      if(whichReentranceTube == 3)
      {
        offset_x = 0 * m;
        offset_y = -1 * m;
      }

      G4double ran_rad = curad * cm + BoratedPETouterrad * cm * rndm(generator);
      G4double ran_phi = 360 * deg * rndm(generator);

      ran_x = ran_rad * sin(ran_phi) + offset_x;
      ran_y = ran_rad * cos(ran_phi) + offset_y;
      ran_z = cuhheight * cm * (1 - 2 * rndm(generator));
    }

    if(type == 2)
    {
      G4int                              BPE_N = fDetector->GetBoratedTurbinezNPanels();
      double                             anglePanel = 360. / BPE_N * deg;
      std::uniform_int_distribution<int> distribution_2(0, BPE_N);
      G4int                              whichPanel = distribution_2(generator);

      G4double BPE_rad  = fDetector->GetBoratedTurbineRadius();
      G4double offset_x = BPE_rad * cm * std::cos(whichPanel * anglePanel);
      G4double offset_y = BPE_rad * cm * std::sin(whichPanel * anglePanel);

      G4double BPE_wid  = fDetector->GetBoratedTurbineWidth();
      G4double BPE_len  = fDetector->GetBoratedTurbineLength();
      G4double BPE_hei  = fDetector->GetBoratedTurbineHeight();
      G4double BPE_ang  = fDetector->GetBoratedTurbineAngle();
      G4double BPE_zPos = fDetector->GetBoratedTurbinezPosition() * cm - 100 * cm;

      G4double tmp_x = BPE_wid / 2. * cm * (1 - 2 * rndm(generator));
      G4double tmp_y = BPE_len / 2. * cm * (1 - 2 * rndm(generator));

      G4double tmp_ang = whichPanel * anglePanel + BPE_ang * deg;

      ran_x = (tmp_x * cos(tmp_ang) + tmp_y * sin(tmp_ang)) + offset_x;
      ran_y = (tmp_y * cos(tmp_ang) - tmp_x * sin(tmp_ang)) + offset_y;
      ran_z = BPE_hei / 2. * cm * (1 - 2 * rndm(generator)) + BPE_zPos;
    }

    if(type == 3)
    {
      G4double BPE_rad  = fDetector->GetBoratedTurbineRadius();
      G4double BPE_wid  = fDetector->GetBoratedTurbineWidth();
      G4double BPE_hei  = fDetector->GetBoratedTurbineHeight() / 2.;
      G4double BPE_zPos = fDetector->GetBoratedTurbinezPosition() * cm - 100 * cm;

      G4double volume_cyl =
        3.1415926535 * BPE_hei * 2 * (pow(BPE_rad + BPE_wid, 2) - pow(BPE_rad, 2));
      G4double volume_top = 3.1415926535 * BPE_wid * pow(BPE_rad + BPE_wid, 2);

      G4double prob_cyl = volume_cyl / (volume_cyl + 2 * volume_top);
      G4double prob_top = (1 - prob_cyl) / 2.;

      std::discrete_distribution<> distribution_2({ prob_cyl, prob_top, prob_top });

      G4int where = distribution_2(generator);

      if(where == 0)
      {
        G4double ran_rad = BPE_rad * cm + BPE_wid * cm * rndm(generator);
        G4double ran_phi = 360 * deg * rndm(generator);
        ran_x            = ran_rad * sin(ran_phi);
        ran_y            = ran_rad * cos(ran_phi);
        ran_z            = BPE_hei * (1 - 2 * rndm(generator));
      }
      if(where > 0)
      {
        G4double ran_rad = BPE_rad * cm * rndm(generator);
        G4double ran_phi = 360 * deg * rndm(generator);
        ran_x            = ran_rad * sin(ran_phi);
        ran_y            = ran_rad * cos(ran_phi);
        ran_z            = BPE_wid * (1 - 2 * rndm(generator));
        if(where == 1)
          ran_z += BPE_hei;
        if(where == 2)
          ran_z -= BPE_hei;
      }
    }

    G4double particle_time = 0 * s;
    G4double energy        = (*neutronEnergySpectrumInBPE)(generator) *keV;
    //G4cout << energy << G4endl;
    G4double theta = rndm(generator) * 180. * deg;
    G4double phi   = rndm(generator) * 360. * deg;
    G4double x     = ran_x;
    G4double y     = ran_y;
    G4double z     = ran_z;

    G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
    fParticleGun->SetParticleDefinition(theParticleTable->FindParticle("neutron"));
    G4double theMass = theParticleTable->FindParticle("neutron")->GetPDGMass();

    G4double      totMomentum = std::sqrt(energy * energy + 2 * theMass * energy);
    G4double      pz          = -1 * std::cos(theta);
    G4double      px          = std::sin(theta) * cos(phi);
    G4double      py          = std::sin(theta) * sin(phi);
    G4ThreeVector momentumDir(px, py, pz);

    fParticleGun->SetParticleMomentumDirection(momentumDir);

    fParticleGun->SetParticleEnergy(energy);

    fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));

    fParticleGun->GeneratePrimaryVertex(event);
  }
  if(fGenerator == "ExternalNeutrons")
  {
    if(neutronEnergySpectrumFromOutside == 0)
    {
      string   filelistName = "../data/FluxOverEnergy.txt";
      ifstream filestream;
      filestream.open(filelistName.c_str());
      vector<double> x_val;
      vector<double> y_val;
      double         tmp_x, tmp_y;
      while(filestream >> tmp_x >> tmp_y)
      {
        x_val.push_back(tmp_x);
        y_val.push_back(tmp_y);
        G4cout << tmp_x << " " << tmp_y << G4endl;
      }
      neutronEnergySpectrumFromOutside = new std::piecewise_linear_distribution<double>(
        x_val.begin(), x_val.end(), y_val.begin());
    }

    std::uniform_real_distribution<> rndm(0.0, 1.0);  // azimuth angle

    G4double WaterTankHeight = (650 + 0.8) * cm;
    G4double WaterTankRadius = (550 + 0.6) * cm; 
    // area: 3201372.3 cm^2
    // n per sec: 6.04 /s ec
    G4double Offset          = (200 - 100 - (850 - 650)) * cm;

    G4double area_cyl =
      2 * CLHEP::twopi * WaterTankRadius * 2 * WaterTankHeight;  // 3.1415926535*BPE_hei*2*(pow(BPE_rad+BPE_wid,2)
                                                                 // - pow(BPE_rad,2));
    G4double area_top = CLHEP::twopi / 2. * WaterTankRadius *
                        WaterTankRadius;  // 3.1415926535*BPE_wid*pow(BPE_rad+BPE_wid,2);

    G4double prob_cyl = area_cyl / (area_cyl + 2 * area_top);
    G4double prob_top = (1 - prob_cyl) / 2.;

    std::discrete_distribution<> distribution_2({ prob_cyl, prob_top, prob_top });

    G4int    where = distribution_2(generator);
    G4double px, py, pz, pos_x, pos_y, pos_z;
    if(where == 0)
    {
      G4double pos_phi    = CLHEP::twopi * rndm(generator);
      G4double pos_height = WaterTankHeight * (1 - 2 * rndm(generator));

      pos_x = WaterTankRadius * cos(pos_phi);
      pos_y = WaterTankRadius * sin(pos_phi);
      pos_z = pos_height;

      G4double mom_phi   = CLHEP::twopi / 4. * (3 - 2 * rndm(generator)) + pos_phi;
      G4double mom_theta = CLHEP::twopi / 2. * rndm(generator);

      px = sin(mom_theta) * cos(mom_phi);
      py = sin(mom_theta) * sin(mom_phi);
      pz = cos(mom_theta);

      //            G4cout << pos_phi << " " << pos_height << " | " << mom_phi << " " <<
      //            mom_theta << " | " << sin(mom_theta) << " " << cos(mom_phi) << G4endl;
    }

    if(where > 0)
    {
      G4double pos_phi = CLHEP::twopi * rndm(generator);
      G4double pos_height;
      if(where == 1)
        pos_height = WaterTankHeight;
      if(where == 2)
        pos_height = -WaterTankHeight;
      G4double pos_rad = WaterTankRadius * rndm(generator);

      pos_x = pos_rad * cos(pos_phi);
      pos_y = pos_rad * sin(pos_phi);
      pos_z = pos_height;

      G4double mom_phi   = CLHEP::twopi * rndm(generator);
      G4double mom_theta = -CLHEP::twopi / 4. * rndm(generator);
      if(where == 2)
        mom_theta += CLHEP::twopi / 4.;

      px = sin(mom_theta) * cos(mom_phi);
      py = sin(mom_theta) * sin(mom_phi);
      pz = cos(mom_theta);

      //          G4cout << pos_phi << " " << pos_height << " " << pos_rad << " | " <<
      //          mom_phi << " " << mom_theta << " | " << sin(mom_theta) << " " <<
      //          cos(mom_phi) << G4endl;
    }

    //    G4cout << " ------------ " << pos_x << " " << pos_y << " " << pos_z << " | " <<
    //    px << " " << py << " " << pz << G4endl;

    G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
    fParticleGun->SetParticleDefinition(theParticleTable->FindParticle("neutron"));
    G4ThreeVector momentumDir(px, py, pz);
    fParticleGun->SetParticleMomentumDirection(momentumDir);
    G4double ekin = (*neutronEnergySpectrumFromOutside)(
      generator);  // get random number (*neutronEnergySpectrumInBPE)(generator) * keV
    ekin *= MeV;
    fParticleGun->SetParticleEnergy(ekin);
    fParticleGun->SetParticlePosition(G4ThreeVector(pos_x, pos_y, pos_z + Offset));
    fParticleGun->GeneratePrimaryVertex(event);
  }
}

void WLGDPrimaryGeneratorAction::SetGenerator(const G4String& name)
{
  std::set<G4String> knownGenerators = {
    "MeiAndHume",        "Musun",           "Ge77m", "Ge77andGe77m",
    "ModeratorNeutrons", "ExternalNeutrons", "Musun_alternative", "SimpleNeutronGun", "GammaCascade"
  };
  if(knownGenerators.count(name) == 0)
  {
    G4Exception("WLGDPrimaryGeneratorAction::SetGenerator", "WLGD0101", JustWarning,
                ("Invalid generator name '" + name + "'").c_str());
    return;
  }
  fGenerator = name;
}

void WLGDPrimaryGeneratorAction::SetSimpleNeutronGun_coord_x(const G4double& x)
{
  coord_x = x;
}
void WLGDPrimaryGeneratorAction::SetSimpleNeutronGun_coord_y(const G4double& y)
{
  coord_y = y;
}
void WLGDPrimaryGeneratorAction::SetSimpleNeutronGun_coord_z(const G4double& z)
{
  coord_z = z;
}
void WLGDPrimaryGeneratorAction::SetSimpleNeutronGun_ekin(const G4double& ekin)
{
  neutron_ekin = ekin;
}
void WLGDPrimaryGeneratorAction::PetersGammaCascade_ekin(const G4double& ekin)
{
  pgc_neutron_ekin = ekin;
}

void WLGDPrimaryGeneratorAction::shortcutToChangeFileName(const G4String& newFile)
{
  // G4cout <<
  // "___________________________________________________________________________________-"
  // << G4endl; G4cout << "MUSUN FileName:    " << newFile << G4endl;
  ChangeFileName(newFile);
}

void WLGDPrimaryGeneratorAction::DefineCommands()
{
  // Define /WLGD/generator command directory using generic messenger class
  fMessenger =
    new G4GenericMessenger(this, "/WLGD/generator/", "Primary generator control");

  // depth command
  auto& depthCmd = fMessenger->DeclareProperty("depth", fDepth,
                                               "Underground laboratory depth [km.w.e.].");
  depthCmd.SetParameterName("d", true);
  depthCmd.SetRange("d>=0.");
  depthCmd.SetDefaultValue("0.");

  // musun file command
  auto& musunfileCmd =
    fMessenger
      ->DeclareMethod("setMUSUNFile",
                      &WLGDPrimaryGeneratorAction::shortcutToChangeFileName)
      .SetGuidance("Set MUSUN file name")
      .SetParameterName("filename", false)
      .SetDefaultValue("./musun_gs_100M.dat");
  // generator command
  // switch command
  fMessenger->DeclareMethod("setGenerator", &WLGDPrimaryGeneratorAction::SetGenerator)
    .SetGuidance("Set generator model of primary muons")
    .SetGuidance("SimpleNeutronGun = generate neutrons with zero energy at a certain location")
    .SetGuidance("MeiAndHume = WW standard case")
    .SetGuidance("Musun = Used in previous MaGe simulation")
    .SetGuidance("Musun_alternative = Alternative Musun input")
    .SetGuidance("Ge77m = generate Ge77m inside the HPGe detectors")
    .SetGuidance("Ge77andGe77m = generate 50% Ge77, 50% Ge77m inside the HPGe detectors")
    .SetGuidance("GammaCascade = generate Peters gamma cascade inside the HPGe detectors")
    .SetGuidance("ModeratorNeutrons = generate neutrons inside the neutron moderators")
    .SetGuidance("ExternalNeutrons = generate neutrons from outside the water tank")
    .SetCandidates(
      "MeiAndHume Musun Musun_alternative Ge77m Ge77andGe77m ModeratorNeutrons ExternalNeutrons SimpleNeutronGun GammaCascade");

  fMessenger->DeclareMethod("SimpleNeutronGun_coord_x", &WLGDPrimaryGeneratorAction::SetSimpleNeutronGun_coord_x)    
    .SetGuidance("Set the x coordinate for the neutron gun")
    .SetDefaultValue("0");

  fMessenger->DeclareMethod("SimpleNeutronGun_coord_y", &WLGDPrimaryGeneratorAction::SetSimpleNeutronGun_coord_y)    
    .SetGuidance("Set the y coordinate for the neutron gun")
    .SetDefaultValue("0");

  fMessenger->DeclareMethod("SimpleNeutronGun_coord_z", &WLGDPrimaryGeneratorAction::SetSimpleNeutronGun_coord_z)    
    .SetGuidance("Set the z coordinate for the neutron gun")
    .SetDefaultValue("0");
    
  fMessenger->DeclareMethod("SimpleNeutronGun_ekin", &WLGDPrimaryGeneratorAction::SetSimpleNeutronGun_ekin)    
    .SetGuidance("Set the ekin of the neutron")
    .SetDefaultValue("0");

  fMessenger->DeclareMethod("PetersGammaCascade_ekin", &WLGDPrimaryGeneratorAction::PetersGammaCascade_ekin)    
    .SetGuidance("Set the ekin of the neutron")
    .SetDefaultValue("1e-3");


}
