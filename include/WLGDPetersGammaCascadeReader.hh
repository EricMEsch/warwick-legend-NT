#ifndef WLGDPetersGammaCascadeReader_h
#define WLGDPetersGammaCascadeReader_h 1

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "globals.hh"
#include "G4GenericMessenger.hh"
#include <random>


using namespace std;
struct GammaCascadeLine {
    int en; // neutron energy [keV]
    int ex; // excitation energy [keV]
    int m; // multiplicity of gamma cascade
    int em; // missing energy [keV]
    std::vector<int> eg; // energies of photons [keV]
};


class WLGDPetersGammaCascadeReader {
public:
    static WLGDPetersGammaCascadeReader* GetInstance();
    ~WLGDPetersGammaCascadeReader();
    //WLGDPetersGammaCascadeReader& operator=(const WLGDPetersGammaCascadeReader&) = delete;

    void ParseFileList(const std::string& file_list);
    void CloseFiles();
    void OpenFiles();
    int GetIndexFromEnergy(double Ekin);
    void DefineCommands();

    GammaCascadeLine GetNextEntry(double Ekin);

private:

    static WLGDPetersGammaCascadeReader* instance;
    WLGDPetersGammaCascadeReader();
    std::vector<std::string> file_names;
    //std::vector<std::unique_ptr<std::ifstream>> files;
    std::ifstream* files[100];
    std::vector<double> lower_edge;
    std::vector<double> upper_edge;
    G4GenericMessenger* fMessenger;
    G4int fGammaCascadeRandomStartLocation = 0;
  std::random_device rd;
  std::ranlux24      generator;

    void SetGammaCascadeFileList(const G4String& file_name);
    void SetGammaCascadeRandomStartLocation(const int answer);
};



#endif
