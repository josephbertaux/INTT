#ifndef FIT_C
#define FIT_C

#include "config.C"

#include <tmvahelper/TMVAHelper.h>
R__LOAD_LIBRARY(libtmvahelper.so)

#include <filesystem>

void
fit (
) {
	// Helper
	TMVAHelper tmva_helper;
	tmva_helper.read_branches(config::branches);
	tmva_helper.read_training(config::training);
	// tmva_helper.read_cuts(config::pT_cuts[0]); // CHANGE ME

	// Welford online algorithm
	Double_t num = 0, avg = 0, err = 0;
	std::map<Double_t, Int_t> pdf;
	for (auto const& entry : std::filesystem::directory_iterator{config::data_dir}) {
		if (!entry.is_regular_file()) continue;

		std::string filename = entry.path().filename();
		if (filename.find(config::channel) == std::string::npos) continue;
		if (filename.find("sig_KFP") == std::string::npos) continue;

		TTree* tree = tmva_helper.get_tree(entry.path().string(), "DecayTree");
		if (!tree || tmva_helper.branch(tree)) {
			std::cerr << entry.path().c_str() << std::endl;
			continue;
		}

		tree->SetBranchStatus("*", 0);
		tree->SetBranchStatus(config::mass_branch.c_str(), 1);
		Float_t* mass = static_cast<Float_t*>(tmva_helper.get_branch(config::mass_branch));
		for (Int_t n = 0, N = tree->GetEntriesFast(); n < N; ++n) {
			tree->GetEntry(n);

			// Welford's online algorithm
			++num;
			double del_1 = *mass - avg;
			avg += del_1 / num;
			double del_2 = *mass - avg;
			err += del_2 * del_1;

			++pdf[*mass];
		}
	}
	err = sqrt(err / num);

	config::mean =  avg;
	config::sigma = err;

	// Freedman-Diaconis rule
	Double_t bin_width = 3.49 * err / pow(num, 0.3333);
	config::num_bins = (config::max_mass - config::min_mass) / bin_width;

	std::cout << "count: " << num << std::endl;
	std::cout << "mean:  " << config::mean  << std::endl;
	std::cout << "sigma: " << config::sigma << std::endl;
	std::cout << "nbins: " << config::num_bins << std::endl;

	// fill hist
	TH1D* fit_hist = new TH1D (
		"mass_fit_hist", "mass_fit_hist",
		config::num_bins, config::min_mass, config::max_mass
	);
	fit_hist->SetLineColor(kBlue);
	for (auto const& [mass_val, count] : pdf) {
		Int_t bin = fit_hist->FindBin(mass_val);
		fit_hist->SetBinContent(bin, fit_hist->GetBinContent(bin) + 1);
	}

	// fit hist
	TF1* fit_func = new TF1 (
		(config::channel + "_mass_fit").c_str(),
		"gausn(0) + gausn(3)",
		config::min_mass, config::max_mass
	);
	fit_func->SetLineColor(kRed);

	fit_func->SetParameter(0, 0.5 * num * bin_width);
	fit_func->SetParameter(1, config::mean);
	fit_func->SetParameter(2, 0.5 * config::sigma);

	fit_func->SetParameter(3, 0.5 * num * bin_width);
	fit_func->SetParameter(4, config::mean);
	fit_func->SetParameter(5, 2.0 * config::sigma);

	// Draw
	fit_hist->Fit(fit_func, "L");
	fit_hist->Draw();
	fit_func->Draw("same");

	// fit_hist->Fit(fit_func, "L");

	// fit_hist->Fit("gausn", "+L");
	// TF1* fit_func = dynamic_cast<TF1*>(fit_hist->GetListOfFunctions()->FindObject("gausn"));
	// if (!fit_func) {
	// 	std::cerr << "func" << std::endl;
	// 	return;
	// }

}

#endif//FIT_C
