// us
#include "WLGDPrimaryGeneratorAction.hh"
#include "WLGDDetectorConstruction.hh"

// geant
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// std
#include <set>
#include <fstream>
#include <random>


//G4String WLGDPrimaryGeneratorAction::fFileName;
//std::ifstream* WLGDPrimaryGeneratorAction::fInputFile;

WLGDPrimaryGeneratorAction::WLGDPrimaryGeneratorAction(WLGDDetectorConstruction* det)
        : G4VUserPrimaryGeneratorAction()
        , fDetector(det)
        , fParticleGun(nullptr)
        , fMessenger(nullptr)
        , fDepth(0.0)
        , fGenerator("Musun")
        , fZShift(200.0*cm)
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
    if (fInputFile.is_open()) fInputFile.close();
}

void WLGDPrimaryGeneratorAction::OpenFile()
{

    fInputFile.open(fFileName,std::ifstream::in);

    if (!(fInputFile.is_open())) {//not open correctly

        //  G4cerr << "File not valid!" << G4endl;
    }
    //else
    //  G4cout << "Was opened!" << G4endl;
}

void WLGDPrimaryGeneratorAction::ChangeFileName(G4String newFile)
{
    if (fFileName != newFile) //check if the new file is equal to the other
    {
        if (fInputFile.is_open()) fInputFile.close(); //close the old file
        G4cout << "opening file: " << newFile << G4endl;
        fFileName = newFile;
        OpenFile(); //open the new one

    }
}


void WLGDPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
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

        G4int nEvent=0;
        G4double time=0.0;
        G4double energy = 0.0*MeV;
        G4double px,py,pz;
        G4double theta,phi;
        G4double x = 0, y = 0, z = 0;
        G4int particleID = 0;

        fInputFile >> nEvent >> particleID >> energy >> x >> y >> z >> theta >> phi;

        // G4cout  << nEvent << " " << x << " " << y << " " << z << G4endl;
        if (fInputFile.eof())
        {
            fInputFile.close();
            G4cerr << "File over: not enough events!" << G4endl;
            G4Exception("WLGDPrimaryGeneratorAction::GeneratePrimaryVertex()", "err001", FatalException, "Exit Warwick");
            return;
        }

        G4double particle_time = time*s;
        energy = energy*GeV;
        theta = theta*rad;
        phi = phi*rad;
        x = x * cm;
        y = y * cm;
        z = fZShift + (z * cm);

        //   G4cout << "Primary coordinates: " << position/m << " m" << G4endl;
        //   G4cout << "Primary coordinates: " << x/cm << " " <<  y/cm << " " << z/cm << " " << G4endl;
        //   G4cout << "Primary energy: " << energy/GeV << " GeV" << G4endl;
        //   G4cout << "Theta: " << theta/deg << " deg; Phi: " << phi/deg << " deg" << G4endl;

        G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();

        G4String particleName = " ";

        if (particleID == 10) particleName = "mu+";
        else particleName = "mu-";

        G4double theMass = theParticleTable->FindParticle(particleName)->GetPDGMass();
        G4double totMomentum = std::sqrt(energy*energy+2*theMass*energy);
        pz = -1*std::cos(theta);
        px = std::sin(theta)*cos(phi);
        py = std::sin(theta)*sin(phi);
        G4ThreeVector momentumDir(px, py, pz);

        fParticleGun->SetParticleMomentumDirection(momentumDir);

        fParticleGun->SetParticleEnergy(energy);

        fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));

        fParticleGun->GeneratePrimaryVertex(event);

    }
    if(fGenerator == "Ge77m")
    {




        G4double cushift	     = 150.;
        G4double roiradius     = 30.0;   // string radius curad - Ge radius - gap
        G4double gerad          = 4.0;                      // Ge radius
        G4int    nofLayers      = 8;   // 8 Ge + 7 gaps = 1010 mm string height
        G4int    nofStrings     = 12;  // 12 strings  of 8 Ge each
        G4double gehheight      = 5.0;                      // full height 10 cm
        G4double gegap          = 3.0;                      // gap between Ge 3cm
        G4double step = (gehheight + gegap / 2) * cm;
        G4double layerthickness = gegap + 2 * gehheight;    // 13 cm total
        G4double angle = CLHEP::twopi / nofStrings;

        std::uniform_int_distribution<int> distribution(0,4*(int)(nofLayers*nofStrings));
        std::uniform_real_distribution<> rndm(0.0, 1.0);

        G4int detectorNumber = distribution(generator);
        G4int whichReentranceTube = detectorNumber/((int)(nofLayers*nofStrings));
        G4int z_number = detectorNumber%nofLayers;
        G4int string_number = (detectorNumber%(nofLayers*nofStrings))/nofLayers;
        G4double radius = gerad * cm * rndm(generator);
        G4double ṕhi = CLHEP::twopi * rndm(generator);
        G4double height = gehheight * cm * (1 - 2*rndm(generator));

        double offset_x, offset_y;
        if(whichReentranceTube == 0) {offset_x = 1*m; offset_y = 0*m;}
        if(whichReentranceTube == 1) {offset_x = 0*m; offset_y = 1*m;}
        if(whichReentranceTube == 2) {offset_x = -1*m; offset_y = 0*m;}
        if(whichReentranceTube == 3) {offset_x = 0*m; offset_y = -1*m;}

        G4double particle_time = 0 * s;
        G4double energy = 0 * GeV;
        G4double theta = 0 * rad;
        G4double phi = 0 * rad;
        G4double x = radius*cos(ṕhi) + roiradius * cm * std::cos(string_number * angle) + offset_x;
        G4double y = radius*sin(ṕhi) + roiradius * cm * std::sin(string_number * angle) + offset_y;
        G4double z = cushift*cm - step + (nofLayers / 2 * layerthickness - z_number * layerthickness)*cm; ;

        //   G4cout << "Primary coordinates: " << position/m << " m" << G4endl;
        //   G4cout << "Primary coordinates: " << x/cm << " " <<  y/cm << " " << z/cm << " " << G4endl;
        //   G4cout << "Primary energy: " << energy/GeV << " GeV" << G4endl;
        //   G4cout << "Theta: " << theta/deg << " deg; Phi: " << phi/deg << " deg" << G4endl;

        G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();

        G4int particleID = 1000320771;


        fParticleGun->SetParticleDefinition(theParticleTable->GetIonTable()->FindIon(32,77,1));


        G4double theMass = theParticleTable->GetIonTable()->FindIon(32,77,1)->GetPDGMass();
        G4double totMomentum = std::sqrt(energy*energy+2*theMass*energy);
        G4double pz = -1*std::cos(theta);
        G4double px = std::sin(theta)*cos(phi);
        G4double py = std::sin(theta)*sin(phi);
        G4ThreeVector momentumDir(px, py, pz);

        fParticleGun->SetParticleMomentumDirection(momentumDir);

        fParticleGun->SetParticleEnergy(energy);

        fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));

        fParticleGun->GeneratePrimaryVertex(event);

    }
}

void WLGDPrimaryGeneratorAction::SetGenerator(const G4String& name)
{
//  G4cout << "__________________________________________xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx__________________________________________" << G4endl << name << G4endl << "__________________________________________xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx__________________________________________" << G4endl;


    std::set<G4String> knownGenerators = { "MeiAndHume", "Musun" , "Ge77m"};
    if(knownGenerators.count(name) == 0)
    {
        G4Exception("WLGDPrimaryGeneratorAction::SetGenerator", "WLGD0101", JustWarning,
                    ("Invalid generator name '" + name + "'").c_str());
        return;
    }
    fGenerator = name;
}

void WLGDPrimaryGeneratorAction::shortcutToChangeFileName(const G4String& newFile){
    // G4cout << "___________________________________________________________________________________-" << G4endl;
    // G4cout << "MUSUN FileName:    " << newFile << G4endl;
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
    auto& musunfileCmd = fMessenger->DeclareMethod("setMUSUNFile",&WLGDPrimaryGeneratorAction::shortcutToChangeFileName)
            .SetGuidance("Set MUSUN file name")
            .SetParameterName("filename", false)
            .SetDefaultValue("./musun_gs_100M.dat");
    // generator command
    // switch command
    fMessenger->DeclareMethod("setGenerator", &WLGDPrimaryGeneratorAction::SetGenerator)
            .SetGuidance("Set generator model of primary muons")
            .SetGuidance("MeiAndHume = WW standard case")
            .SetGuidance("Musun = Used in previous MaGe simulation")
            .SetGuidance("Ge77m = generate Ge77m inside the HPGe detectors")
            .SetCandidates("MeiAndHume Musun Ge77m" );
}

