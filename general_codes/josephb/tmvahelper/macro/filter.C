#ifndef FILTER_C
#define FILTER_C

#include "config.C"
#include <sPhenixStyle.C>

#include <tmvahelper/TMVAHelper.h>
R__LOAD_LIBRARY(libtmvahelper.so)

#include <filesystem>
#include <boost/format.hpp>

void
filter (
	std::string const& data_dir
) {
	// Helper
	TMVAHelper tmva_helper;
	tmva_helper.read_branches(config::branches);
	tmva_helper.read_training(config::training);
	tmva_helper.read_cuts(config::signal_cuts);

	int num_real{0};
	int total{0};
	for (auto const& entry : std::filesystem::directory_iterator{data_dir}) {
		if (!entry.is_regular_file()) continue;

		std::string filename = entry.path().filename();
		if (filename.find(config::channel) == std::string::npos) continue;
		if (filename.find("sig_KFP") == std::string::npos) continue;

		TTree* tree = tmva_helper.get_tree(entry.path().string(), "DecayTree");
		if (!tree || tmva_helper.branch(tree)) {
			std::cerr << entry.path().c_str() << std::endl;
			continue;
		}

		int parent_index[3];
		std::vector<int>** true_track_history_PDG_ID = new std::vector<int>*[3];
		std::vector<int>** true_track_history_px = new std::vector<int>*[3];

		for (int i = 0; i < 3; ++i) {
			std::string name;

			true_track_history_PDG_ID[i] = new std::vector<int>;
			true_track_history_px[i] = new std::vector<int>;

			name = (boost::format("track_%d_true_track_history_PDG_ID") % (i + 1)).str();
			tree->SetBranchAddress(name.c_str(), &(true_track_history_PDG_ID[i]));

			name = (boost::format("track_%d_true_track_history_px") % (i + 1)).str();
			tree->SetBranchAddress(name.c_str(), &(true_track_history_px[i]));
		}

		for (Long64_t n = 0, N = tree->GetEntriesFast(); n < N; ++n) {
			tree->GetEntry(n);
			++total;

			bool should_continue = false;
			for (int i = 0; i < 3; ++i ) {
				if (!true_track_history_PDG_ID[i]) should_continue = true;
				if (!true_track_history_px[i]) should_continue = true;
			}
			if (should_continue) continue;

			for (int i = 0; i < 3; ++i) {
				for(parent_index[i] = 0; parent_index[i] < true_track_history_PDG_ID[i]->size(); ++parent_index[i]) {
					if (abs(true_track_history_PDG_ID[i]->at(parent_index[i])) == 4122) break;
				}
				if (parent_index[i] == true_track_history_PDG_ID[i]->size()) should_continue = true;
			}
			if (should_continue) continue;

			for (int i = 0; i < 2; ++i) {
				for (int j = i; j < 3; ++j) {
					if (abs(true_track_history_px[i]->at(parent_index[i]) - true_track_history_px[j]->at(parent_index[j])) < 1E-4) continue; 
					should_continue = true;
					break;
				}
				if (should_continue) break;
			}
			if (should_continue) continue;

			++num_real;
		}

		for (int i = 0; i < 3; ++i) {
			delete true_track_history_PDG_ID[i];
			delete true_track_history_px[i];
		}
		delete[] true_track_history_PDG_ID;
		delete[] true_track_history_px;

	}

	std::cout
		<< "num_real: " << num_real << "\n"
		<< "total:    " << total << "\n"
		<< std::flush;

	// Draw
	// SetsPhenixStyle();
	// TCanvas cnvs(
	// 	(config::channel + "_sig_fit_cnvs").c_str(),
	// 	(config::channel + "_sig_fit_cnvs").c_str(),
	// 	600, 800
	// );
	// cnvs.cd();

}

#endif//FILTER_C
