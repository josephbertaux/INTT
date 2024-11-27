#ifndef TRAIN_C
#define TRAIN_C

#include "config.C"

#include <tmvahelper/TMVAHelper.h>
R__LOAD_LIBRARY(libtmvahelper.so)

void
train (
	std::vector<std::string> const& signal_files = {
		"outputKFParticle_D0_Kpi_0.root",
		"outputKFParticle_D0_Kpi_1.root",
		"outputKFParticle_D0_Kpi_2.root",
		"outputKFParticle_D0_Kpi_3.root",
		"outputKFParticle_D0_Kpi_4.root",
		"outputKFParticle_D0_Kpi_5.root",
		"outputKFParticle_D0_Kpi_6.root",
		"outputKFParticle_D0_Kpi_7.root",
		"outputKFParticle_D0_Kpi_8.root",
		"outputKFParticle_D0_Kpi_9.root",
	},
	std::vector<std::string> const& background_files = {
		"outputMinBiasKFParticle_D2Kpi_0.root",
		"outputMinBiasKFParticle_D2Kpi_1.root",
		"outputMinBiasKFParticle_D2Kpi_2.root",
		"outputMinBiasKFParticle_D2Kpi_3.root",
		"outputMinBiasKFParticle_D2Kpi_4.root",
		"outputMinBiasKFParticle_D2Kpi_5.root",
		"outputMinBiasKFParticle_D2Kpi_6.root",
		"outputMinBiasKFParticle_D2Kpi_7.root",
		"outputMinBiasKFParticle_D2Kpi_8.root",
		"outputMinBiasKFParticle_D2Kpi_9.root",
	}
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
	for (auto const& signal_file : signal_files) {
		TTree* tree = tmva_helper.get_tree(config::data_dir + "/" + signal_file, "DecayTree");
		if (!tree) {
			std::cerr << "file: " << config::data_dir + "/" + signal_file << std::endl;
			continue;
		}
		dataloader->AddSignalTree(tree);
	}

	for (auto const& background_file : background_files) {
		TTree* tree = tmva_helper.get_tree(config::data_dir + "/" + background_file, "DecayTree");
		if (!tree) {
			std::cerr << "file: " << config::data_dir + "/" + background_file << std::endl;
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

