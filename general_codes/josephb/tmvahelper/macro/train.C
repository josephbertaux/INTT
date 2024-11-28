#ifndef TRAIN_C
#define TRAIN_C

#include "config.C"

#include <tmvahelper/TMVAHelper.h>
R__LOAD_LIBRARY(libtmvahelper.so)

#include <filesystem>

void
train (
) {
	// Helper
	TMVAHelper tmva_helper;
	tmva_helper.read_branches(config::branches);
	tmva_helper.read_training(config::training);

	// Initialize factory and dataloader
	TFile* factory_file = TFile::Open(config::factory_file_name.c_str(), "RECREATE");
	if (!factory_file) {
		std::cerr << "file: " << config::factory_file_name << std::endl;
		return EXIT_FAILURE;
	}

	TMVA::Factory* factory = new TMVA::Factory (
		"factory", factory_file,
		"!V:!Silent:AnalysisType=Classification"
	);
	TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataloader");

	// Add variables and cuts
	tmva_helper.branch(dataloader);
	// This method also adds the cuts that have been added to the tmva_helper instance
	// And also adds cuts that protect against NaN values in the input TTrees

	// Sideband cut for training
	dataloader->AddCut(config::get_sideband_cut(), "Background");

	// Add input files
	for (auto const& entry : std::filesystem::directory_iterator{config::data_dir}) {
		if (!entry.is_regular_file()) continue;

		std::string filename = entry.path().filename();
		if (filename.find(config::channel) == std::string::npos) continue;
		if (filename.find("sig_KFP") == std::string::npos) continue;

		TTree* tree = tmva_helper.get_tree(entry.path().string(), "DecayTree");
		if (!tree) {
			std::cerr << "file: " << entry.path() << std::endl;
			continue;
		}
		dataloader->AddSignalTree(tree);
	}

	for (auto const& entry : std::filesystem::directory_iterator{config::data_dir}) {
		if (!entry.is_regular_file()) continue;

		std::string filename = entry.path().filename();
		if (filename.find(config::channel) == std::string::npos) continue;
		if (filename.find("bak_KFP") == std::string::npos) continue;

		TTree* tree = tmva_helper.get_tree(entry.path().string(), "DecayTree");
		if (!tree) {
			std::cerr << "file: " << entry.path() << std::endl;
			continue;
		}
		dataloader->AddBackgroundTree(tree);
	}

	// Train
	factory->BookMethod(dataloader, config::method_type, config::method_name.c_str(), config::method_options.c_str());
	factory->TrainAllMethods();
	factory->TestAllMethods();
	factory->EvaluateAllMethods();

	// Get Cut value from the training
	Double_t sig_cut, max_sig;
	auto method = dynamic_cast<TMVA::MethodBase*>(factory->GetMethod("dataloader", config::method_name.c_str()));
	sig_cut = method->GetMaximumSignificance(config::ratio, 1.0, max_sig);

	// Cleanup
	factory_file->Close();
	delete dataloader;
	delete factory;

	config::cut_val = sig_cut;
	std::cout << "cut: " << sig_cut << " significance: " << max_sig << std::endl;
	// std::cout << config::get_sideband_cut().GetTitle() << std::endl;
}

#endif//TRAIN_C

