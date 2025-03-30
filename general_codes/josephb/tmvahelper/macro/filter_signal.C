#ifndef FILTER_SIGNAL_C
#define FILTER_SIGNAL_C

#include "config.C"

#include <tmvahelper/TmvaHelper.h>
#include <tmvahelper/PidFilter.h>
R__LOAD_LIBRARY(libtmvahelper.so)

#include <filesystem>
#include <boost/format.hpp>

void
filter_signal (
	std::string const& data_dir
) {
	// Helper
	TmvaHelper tmva_helper;
	tmva_helper.read_branches(config::branches);
	tmva_helper.read_training(config::training);
	tmva_helper.read_cuts(config::cuts);

	PidFilter pid_filter;
	pid_filter.set_mother_pdg_id(config::particle_trigger);
	pid_filter.set_num_daughters(3);

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
		if (!tree || tmva_helper.branch(tree) || pid_filter.branch(tree)) {
			std::cerr << entry.path().c_str() << std::endl;
			break;
		}
		++files;

		for (Long64_t n = 0, N = tree->GetEntriesFast(); n < N; ++n) {
			tree->GetEntry(n);

			if (tmva_helper.eval()) continue;
			if (pid_filter.eval()) continue;

			signal_tree->Fill();
			++counts;
		}

		tree->GetDirectory()->GetFile()->Close();
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
