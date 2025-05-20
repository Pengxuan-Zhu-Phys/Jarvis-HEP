#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"
#include <iostream>
#include <string>
#include <memory>  // Include for std::shared_ptr and std::make_shared

using namespace Pythia8;

class NoPtCutOnSUSY : public UserHooks {
public:
    NoPtCutOnSUSY(const std::vector<int>& susyIds) : susyIds_(susyIds) {}

    virtual bool doVetoPT(int iPos, const Event& event) {
        for (int id : susyIds_) {
            if (event[iPos].id() == id) return true;
        }
        return false;
    }

private:
    std::vector<int> susyIds_;
};

int main(int argc, char* argv[]) {
    // Variables to store command line arguments
    std::string slhaFile, outputFile;
    int numberOfEvents = 1000;  // Default number of events
    int randomSeed = 0;  // Default seed, 0 means random

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-i" && i + 1 < argc) {
            slhaFile = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            outputFile = argv[++i];
        } else if (arg == "-n" && i + 1 < argc) {
            numberOfEvents = std::stoi(argv[++i]);
        } else if (arg == "-s" && i + 1 < argc) {
            randomSeed = std::stoi(argv[++i]);
        } else {
            std::cerr << "Usage: " << argv[0] << " -i <slha_file> -o <output_file> -n <number_of_events> -s <random_seed>" << std::endl;
            return 1;
        }
    }

    // Ensure both input and output files are provided
    if (slhaFile.empty() || outputFile.empty()) {
        std::cerr << "Error: Both SLHA input file and HepMC output file must be specified." << std::endl;
        return 1;
    }

    // Initialize Pythia
    Pythia pythia;

    // Set the SLHA file from the command line
    pythia.readString("SLHA:file = " + slhaFile);

    // Enable SUSY process: chargino-neutralino production
    pythia.readString("SUSY:all = off");  // First turn off all SUSY processes
    pythia.readString("SUSY:qqbar2chi+-chi0 = on");  // Enable qqbar to chargino-neutralino
    pythia.readString("1000022:mayDecay = false"); // Neutralino chi_10
    pythia.readString("1000023:mayDecay = false"); // Neutralino chi_20
    pythia.readString("1000025:mayDecay = false"); // Neutralino chi_30
    pythia.readString("1000035:mayDecay = false"); // Neutralino chi_40
    pythia.readString("1000024:mayDecay = false"); // Chargino chi_1+
    pythia.readString("1000037:mayDecay = false"); // Chargino chi_2+
    // Set center of mass energy (e.g., 13 TeV for LHC)

    pythia.readString("Beams:idA = 2212;");  // PDG code for proton
    pythia.readString("Beams:idB = 2212;");  // PDG code for proton
    pythia.readString("Beams:eCM = 13000.");
    // pythia.readString("HardQCD:all = on");   // 举例：启动所有硬QCD过程
    // pythia.readString("PartonLevel:All = off");
    pythia.readString("HadronLevel:All = off");

    pythia.readString("Init:showChangedSettings = on");
    pythia.readString("Next:numberShowInfo = 1"); // 显示每个事件的信息
    pythia.readString("Next:numberShowProcess = 1"); // 显示过程信息
    pythia.readString("Next:numberShowEvent = 1"); // 显示事件详情

    // pythia.readString("PhaseSpace:pTHatMinDiverge = 0.5");

    pythia.settings.forceParm("PhaseSpace:pTHatMinDiverge", 1e-9);
    // pythia.readString("PhaseSpace:pTHatMin = 0.0");
    // pythia.readString("PhaseSpace:pTHatMax = 50.0");
    std::vector<int> susyParticleIds = {1000022, 1000023, 1000025, 1000035, 1000024, 1000037};
    auto noPtCutOnSUSY = std::make_shared<NoPtCutOnSUSY>(susyParticleIds);
    pythia.setUserHooksPtr(noPtCutOnSUSY);

    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(randomSeed));  // 0 表示每次都自动选择一个种子

    // Initialize the Pythia simulation
    pythia.init();

    // Create an instance to convert Pythia events to HepMC3 format
    HepMC3::Pythia8ToHepMC3 toHepMC3;

    // Create a HepMC3 event output file with the specified output file name
    HepMC3::WriterAscii hepmcOutput(outputFile);

    // Generate events
    for (int iEvent = 0; iEvent < numberOfEvents; ++iEvent) {
        if (!pythia.next()) continue;

        // Create a HepMC event
        HepMC3::GenEvent hepmcEvent(HepMC3::Units::GEV, HepMC3::Units::MM);

        // Fill the HepMC event with data from Pythia
        toHepMC3.fill_next_event(pythia, &hepmcEvent);

        // Write the HepMC event to file
        hepmcOutput.write_event(hepmcEvent);
    }

    // Close the output file
    hepmcOutput.close();

    // Print statistics
    pythia.stat();
    return 0;
}
