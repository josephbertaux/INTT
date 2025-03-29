#ifndef FILTER_BACKGROUND_C
#define FILTER_BACKGROUND_C

#include "config.C"
#include "filter_signal.C"
#include "fit.C"

#include <tmvahelper/TMVAHelper.h>
R__LOAD_LIBRARY(libtmvahelper.so)

#include <filesystem>
#include <boost/format.hpp>

void
filter_background (
	std::string const& data_dir
) {
	filter_signal(data_dir);
	fit();

	std::vector<std::string> background_cuts = config::cuts;
	background_cuts.push_back(config::get_sideband_cut().GetTitle());

	// Helper
	TMVAHelper tmva_helper;
	tmva_helper.read_branches(config::branches);
	tmva_helper.read_training(config::training);
	tmva_helper.read_cuts(background_cuts);

	TFile* background_file = TFile::Open("background.root", "RECREATE");
	TTree* background_tree = new TTree("DecayTree", "DecayTree");
	background_tree->SetDirectory(background_file);
	tmva_helper.make_branches(background_tree);

	Long64_t counts{0};
	Long64_t files{0};
	for (auto const& entry : std::filesystem::directory_iterator{data_dir}) {
		if (!entry.is_regular_file()) continue;

		std::string filename = entry.path().filename();
		if (filename.find(config::channel) == std::string::npos) continue;
		if (filename.find("bak_KFP") == std::string::npos) continue;

		TTree* tree = tmva_helper.get_tree(entry.path().string(), "DecayTree");
		if (!tree || tmva_helper.branch(tree)) {
			std::cerr << entry.path().c_str() << std::endl;
			continue;
		}
		++files;

		for (Long64_t n = 0, N = tree->GetEntriesFast(); n < N; ++n) {
			tree->GetEntry(n);

			if (tmva_helper.eval()) continue;

			background_tree->Fill();
			++counts;
		}
	}

	background_file->cd();
	background_tree->Write();
	background_file->Write();
	background_file->Close();

	std::cout << "finished" << std::endl;
	std::cout << "counts: " << counts << std::endl;
	std::cout << "files:  " << files << std::endl;
}

#endif//FILTER_BACKGROUND_C
