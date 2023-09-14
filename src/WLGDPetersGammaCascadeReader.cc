#include "WLGDPetersGammaCascadeReader.hh"

WLGDPetersGammaCascadeReader* WLGDPetersGammaCascadeReader::instance = nullptr;

WLGDPetersGammaCascadeReader* WLGDPetersGammaCascadeReader::GetInstance() {
    if (instance == nullptr) {
        instance = new WLGDPetersGammaCascadeReader();
    }
    return instance;
}

void WLGDPetersGammaCascadeReader::ParseFileList(const std::string& file_list){
    CloseFiles();
    std::ifstream file_list_file = std::ifstream(file_list);
    file_names.clear();
    lower_edge.clear();
    upper_edge.clear();
    if(file_list_file.is_open()){
        std::string line;
        std::string tmp_file_name;
        double tmp_lower_edge, tmp_upper_edge;
        G4cout << ">>> Reading Peters gamma cascade file list:" << G4endl;
        while(std::getline(file_list_file, line)){
            istringstream ss(line);
            ss >> tmp_file_name >> tmp_lower_edge >> tmp_upper_edge;
            file_names.push_back(tmp_file_name);
            lower_edge.push_back(tmp_lower_edge);
            upper_edge.push_back(tmp_upper_edge);
            G4cout << tmp_file_name << ": " << tmp_lower_edge << " - " << tmp_upper_edge << G4endl;
        }
        G4cout << ">>> Finish" << G4endl;
    }
    else{
            throw std::runtime_error("Failed to open file list " + file_list);
    }
}

void WLGDPetersGammaCascadeReader::CloseFiles(){
    for (int i = 0; i < file_names.size(); i++) {
        if (files[i]->is_open()) {
            G4cout << "Trying to close: " << files[i] << G4endl;
            files[i]->close();
        }
    }
}

void WLGDPetersGammaCascadeReader::OpenFiles() {
    G4cout << ">>> fGammaCascadeRandomStartLocation: " << fGammaCascadeRandomStartLocation << G4endl;
    G4cout << ">>> Open Peters gamma cascade file list:" << G4endl;

  std::uniform_real_distribution<> rndm(0.0, 1.0);
    for (int i = 0; i < file_names.size(); i++) {
        G4cout << "File " << i << ": " << file_names[i] << G4endl;
        files[i] = new ifstream(file_names[i]);
    
        std::string line;
        int header_length = 0;
        do{
            std::getline(*(files[i]), line);
            header_length++;
        } while(line[0] == '%' || (line.find("version") != std::string::npos));

        if(fGammaCascadeRandomStartLocation){

            int n_entries_in_file = 0;
            while(std::getline(*(files[i]), line))
                n_entries_in_file++;

            files[i]->clear(); // clear EOF flag
            files[i]->seekg(0, std::ios::beg); // move to beginning of file

            int start_location = (int) (n_entries_in_file * rndm(generator) + header_length);

            G4cout << ">>> Random start location: " << start_location << G4endl;
            for(int j = 0; j < start_location; j++)
                std::getline(*(files[i]), line);

        }
    }

    G4cout << ">>> Finish" << G4endl;
}

int WLGDPetersGammaCascadeReader::GetIndexFromEnergy(double Ekin){
    for(int i = 0; i < lower_edge.size(); i++)
        if(Ekin > lower_edge[i] && Ekin < upper_edge[i])
            return i;
    return 0; // If no falling in a defined range, take first file.
}

GammaCascadeLine WLGDPetersGammaCascadeReader::GetNextEntry(double Ekin) {
    // determine the file to read from based on energy
    int file_index = GetIndexFromEnergy(Ekin);
    // read next line from file
    std::string line;
    do{
    if (!std::getline(*(files[file_index]), line)) {
        // if end-of-file is reached, reset file and read first line
        G4cout << "Reopening file " << file_names[file_index] << G4endl;
        files[file_index]->clear(); // clear EOF flag
        files[file_index]->seekg(0, std::ios::beg); // move to beginning of file
        if (!std::getline(*(files[file_index]), line)) {
            throw std::runtime_error("Failed to read next line from file");
        }
    }
    } while(line[0] == '%' || (line.find("version") != std::string::npos));

    // parse line and return as struct
    GammaCascadeLine gamma_cascade;
    std::istringstream iss(line);
    iss >> gamma_cascade.en >> gamma_cascade.ex >> gamma_cascade.m >> gamma_cascade.em;
    gamma_cascade.eg.reserve(gamma_cascade.m);
    int eg_value;
    for (int i = 0; i < gamma_cascade.m; i++) {
        if (!(iss >> eg_value)) {
            throw std::runtime_error("Failed to read photon energy from line");
        }
        gamma_cascade.eg.push_back(eg_value);
    }
    return gamma_cascade;
}

WLGDPetersGammaCascadeReader::WLGDPetersGammaCascadeReader() {
    generator.seed(rd());  // set a random seed
    DefineCommands();
}

WLGDPetersGammaCascadeReader::~WLGDPetersGammaCascadeReader() {
    CloseFiles();
}

void WLGDPetersGammaCascadeReader::SetGammaCascadeFileList(const G4String& file_name)
{
  ParseFileList(file_name);
  OpenFiles();
}
void WLGDPetersGammaCascadeReader::SetGammaCascadeRandomStartLocation(const int answer)
{
    fGammaCascadeRandomStartLocation = answer;
    G4cout << ">>> setting fGammaCascadeRandomStartLocation to: " << fGammaCascadeRandomStartLocation << G4endl;
    CloseFiles();
    OpenFiles();
}

void WLGDPetersGammaCascadeReader::DefineCommands()
{
  fMessenger =
    new G4GenericMessenger(this, "/WLGD/PetersGammaCascade/", "Controll Peters gamma cascade model");

  auto& GammaCascadeFileListCmd =
    fMessenger
        ->DeclareMethod("SetGammaCascadeFileList",
                    &WLGDPetersGammaCascadeReader::SetGammaCascadeFileList)
        .SetGuidance("Set the path to the file list of all gamma cascade input files")
        .SetParameterName("filename", false)
        .SetDefaultValue("/home/ga63zay/atlas_neuberger/L1K_simulations/peters_gamma_cascade/file_list.txt");


  auto& GammaCascadeRandomStartLocationCmd =
    fMessenger
        ->DeclareMethod("SetGammaCascadeRandomStartLocation",
                    &WLGDPetersGammaCascadeReader::SetGammaCascadeRandomStartLocation)
        .SetGuidance("Set the whether the start location in the gamma cascade file is random or not")
        .SetGuidance("0 = don't")
        .SetGuidance("1 = do")
        .SetCandidates("0 1")
        .SetDefaultValue("0");

}
