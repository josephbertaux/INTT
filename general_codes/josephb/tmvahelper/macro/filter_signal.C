#ifndef FILTER_SIGNAL_C
#define FILTER_SIGNAL_C

#include "config.C"

#include <tmvahelper/TMVAHelper.h>
R__LOAD_LIBRARY(libtmvahelper.so)

#include <filesystem>
#include <boost/format.hpp>

void
filter_signal (
	std::string const& data_dir
) {
	// Helper
	TMVAHelper tmva_helper;
	tmva_helper.read_branches(config::branches);
	tmva_helper.read_training(config::training);
	tmva_helper.read_cuts(config::cuts);

	TFile* signal_file = TFile::Open("signal.root", "RECREATE");
	TTree* signal_tree = new TTree("DecayTree", "DecayTree");
	signal_tree->SetDirectory(signal_file);
	tmva_helper.make_branches(signal_tree);

	Long64_t counts{0};
	Long64_t files{0};
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
		++files;

		Float_t* mass = static_cast<Float_t*>(tmva_helper.get_branch(config::mass_branch));

		int parent_index[3];
		int true_track_ID[3];
		int track_PDG_ID[3];

		std::vector<int>** true_track_history_PDG_ID = new std::vector<int>*[3];

		std::vector<int>** true_track_history_px = new std::vector<int>*[3];
		std::vector<int>** true_track_history_py = new std::vector<int>*[3];
		std::vector<int>** true_track_history_pz = new std::vector<int>*[3];


		for (int i = 0; i < 3; ++i) {
			std::string name;

			true_track_history_PDG_ID[i] = new std::vector<int>;

			true_track_history_px[i] = new std::vector<int>;
			true_track_history_py[i] = new std::vector<int>;
			true_track_history_pz[i] = new std::vector<int>;

			name = (boost::format("track_%d_true_track_history_PDG_ID") % (i + 1)).str();
			tree->SetBranchAddress(name.c_str(), &(true_track_history_PDG_ID[i]));

			name = (boost::format("track_%d_true_track_history_px") % (i + 1)).str();
			tree->SetBranchAddress(name.c_str(), &(true_track_history_px[i]));
			name = (boost::format("track_%d_true_track_history_py") % (i + 1)).str();
			tree->SetBranchAddress(name.c_str(), &(true_track_history_py[i]));
			name = (boost::format("track_%d_true_track_history_pz") % (i + 1)).str();
			tree->SetBranchAddress(name.c_str(), &(true_track_history_pz[i]));

			name = (boost::format("track_%d_true_ID") % (i + 1)).str();
			tree->SetBranchAddress(name.c_str(), &(true_track_ID[i]));

			name = (boost::format("track_%d_PDG_ID") % (i + 1)).str();
			tree->SetBranchAddress(name.c_str(), &(track_PDG_ID[i]));
		}

		for (Long64_t n = 0, N = tree->GetEntriesFast(); n < N; ++n) {
			tree->GetEntry(n);

			if (tmva_helper.eval()) continue;

			bool should_continue = false;
			for (int i = 0; i < 3; ++i ) {
				if (!true_track_history_PDG_ID[i]) should_continue = true;
				if (!true_track_history_px[i]) should_continue = true;
				if (!true_track_history_py[i]) should_continue = true;
				if (!true_track_history_pz[i]) should_continue = true;

				if (track_PDG_ID[i] != true_track_ID[i]) should_continue = true;
			}
			// if (should_continue) {std::cout << __LINE__ << std::endl; continue;}
			if (should_continue) continue;

			for (int i = 0; i < 3; ++i) {
				for(parent_index[i] = 0; parent_index[i] < true_track_history_PDG_ID[i]->size(); ++parent_index[i]) {
					if (abs(true_track_history_PDG_ID[i]->at(parent_index[i])) == 4122) break;
				}
				if (parent_index[i] == true_track_history_PDG_ID[i]->size()) should_continue = true;
			}
			// if (should_continue) {std::cout << __LINE__ << std::endl; continue;}
			if (should_continue) continue;

			for (int i = 0; i < 2; ++i) {
				for (int j = i; j < 3; ++j) {
					if (abs(true_track_history_px[i]->at(parent_index[i]) - true_track_history_px[j]->at(parent_index[j])) > 1E-4) should_continue = true; 
					if (abs(true_track_history_py[i]->at(parent_index[i]) - true_track_history_py[j]->at(parent_index[j])) > 1E-4) should_continue = true; 
					if (abs(true_track_history_pz[i]->at(parent_index[i]) - true_track_history_pz[j]->at(parent_index[j])) > 1E-4) should_continue = true; 
					if (should_continue) break;
				}
				if (should_continue) break;
			}
			// if (should_continue) {std::cout << __LINE__ << std::endl; continue;}
			if (should_continue) continue;

			signal_tree->Fill();
			++counts;
		}

		for (int i = 0; i < 3; ++i) {
			delete true_track_history_PDG_ID[i];
			delete true_track_history_px[i];
			delete true_track_history_py[i];
			delete true_track_history_pz[i];
		}
		delete[] true_track_history_PDG_ID;
		delete[] true_track_history_px;
		delete[] true_track_history_py;
		delete[] true_track_history_pz;
	}

	signal_file->cd();
	signal_tree->Write();
	signal_file->Write();
	signal_file->Close();

	std::cout << "finished" << std::endl;
	std::cout << "counts: " << counts << std::endl;
	std::cout << "files:  " << files << std::endl;
}

#endif//FILTER_SIGNAL_C
