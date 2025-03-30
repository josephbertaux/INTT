#ifndef FILTER_BACKGROUND_C
#define FILTER_BACKGROUND_C

#include "config.C"

#include <tmvahelper/TmvaHelper.h>
R__LOAD_LIBRARY(libtmvahelper.so)

#include <filesystem>
#include <boost/format.hpp>

void
filter_background (
	std::string const& data_dir,
	std::string const& fit_path =
		"/sphenix/tg/tg01/hf/jbertaux/dEdx_fits/dedx_fitparam.root"
) {
	// Get sideband cut from signal roofit
	TFile* file = TFile::Open("signal_fit.root", "READ");
	if (!file) {std::cerr << __LINE__ << std::endl; return;}

	RooWorkspace* w = dynamic_cast<RooWorkspace*>(file->Get("w"));
	if (!w) {std::cerr << __LINE__ << std::endl; return;}
	w->Print();
	RooArgSet args = w->allVars();

	Float_t sigma_eff{0};
	Float_t total{0};
	Float_t mean = dynamic_cast<RooRealVar&>(args["mean"]).getValV();
	for (int i = 0; i < 2; ++i) {
		std::string name;

		name = (boost::format("sigma_%d") % i).str();
		Float_t sigma = dynamic_cast<RooRealVar&>(args[name.c_str()]).getValV();
		std::cout << name << " " << sigma << std::endl;

		name = (boost::format("coeff_%d") % i).str();
		Float_t coeff = dynamic_cast<RooRealVar&>(args[name.c_str()]).getValV();
		std::cout << name << " " << coeff << std::endl;

		sigma_eff += sigma * sigma * coeff;
		total += coeff;
	}
	sigma_eff /= total;
	sigma_eff = sqrt(sigma_eff);
	std::cout << sigma_eff << std::endl;

	std::string sideband_cut = (
		boost::format("(%f < %s && %s < %f) || (%f < %s && %s < %f)")
		% (mean - 6.0 * sigma_eff) % config::mass_branch % config::mass_branch % (mean - 3.0 * sigma_eff)
		% (mean + 3.0 * sigma_eff) % config::mass_branch % config::mass_branch % (mean + 6.0 * sigma_eff)
	).str();

	std::vector<std::string> background_cuts = config::cuts;
	background_cuts.push_back(sideband_cut);
	std::cout << sideband_cut << std::endl;

	// Helper
	TmvaHelper tmva_helper;
	tmva_helper.read_branches(config::branches);
	tmva_helper.read_training(config::training);
	tmva_helper.read_cuts(background_cuts);

	TFile* background_file = TFile::Open("background.root", "RECREATE");
	TTree* background_tree = new TTree("DecayTree", "DecayTree");
	background_tree->SetDirectory(background_file);
	tmva_helper.make_branches(background_tree);

	// dEdx
	Float_t track_dEdx_fit[3];
	for (int i = 0; i < 3; ++i) {
		std::string name;
		name = (boost::format("track_%d_dEdx_fit") % (i + 1)).str();
		background_tree->Branch(name.c_str(), &track_dEdx_fit[i]);
	}

	TFile* fit_file = TFile::Open(fit_path.c_str(), "READ");
	if (!fit_file || !fit_file->IsOpen()) {
		std::cerr
			<< "\t" << fit_path << "\n"
			<< std::flush;
		return;
	}

	std::map<int, TF1*> fit_map;
	fit_map.insert({  211, dynamic_cast<TF1*>(fit_file->Get("f_piband"))});
	fit_map.insert({ -211, dynamic_cast<TF1*>(fit_file->Get("f_piminus_band"))});
	fit_map.insert({  321, dynamic_cast<TF1*>(fit_file->Get("f_Kband"))});
	fit_map.insert({ -321, dynamic_cast<TF1*>(fit_file->Get("f_Kminus_band"))});
	fit_map.insert({ 2212, dynamic_cast<TF1*>(fit_file->Get("f_pband"))});
	fit_map.insert({-2212, dynamic_cast<TF1*>(fit_file->Get("f_pbar_band"))});
	for (auto const& [pid, func_ptr] : fit_map) {
		if (func_ptr) continue;
		std::cerr
			<< "\t" << pid << "\n"
			<< std::flush;
		return;
	}

	Int_t* track_id[3];
	Float_t* track_p[3];

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
			return;
		}
		++files;

		for (int i = 0; i < 3; ++i) {
			std::string name;

			name = (boost::format("track_%d_PDG_ID") % (i + 1)).str();
			track_id[i] = static_cast<Int_t*>(tmva_helper.get_branch(name));
			if (!track_id[i]) {
				std::cerr
					<< "\t" << name << "\n"
					<< std::flush;
				return;
			}

			name = (boost::format("track_%d_p") % (i + 1)).str();
			track_p[i] = static_cast<Float_t*>(tmva_helper.get_branch(name));
			if (!track_p[i]) {
				std::cerr
					<< "\t" << name << "\n"
					<< std::flush;
				return;
			}
		}

		for (Long64_t n = 0, N = tree->GetEntriesFast(); n < N; ++n) {
			tree->GetEntry(n);

			if (tmva_helper.eval()) continue;

			bool should_continue = false;
			for (int i = 0; i < 3; ++i) {
				auto itr = fit_map.find(abs(*track_id[i]));
				if (itr == fit_map.end()) {
					should_continue = true;
					break;
				}
	
				track_dEdx_fit[i] = itr->second->Eval(abs(*track_p[i]));
			}
			if (should_continue) continue;

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
