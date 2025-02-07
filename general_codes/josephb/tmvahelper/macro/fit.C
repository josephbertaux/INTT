#ifndef FIT_C
#define FIT_C

#include "config.C"
#include <sPhenixStyle.C>

#include <tmvahelper/TMVAHelper.h>
R__LOAD_LIBRARY(libtmvahelper.so)

#include <filesystem>

void
fit (
	std::string const& data_dir
) {
	// Helper
	TMVAHelper tmva_helper;
	tmva_helper.read_branches(config::branches);
	tmva_helper.read_training(config::training);
	tmva_helper.read_cuts(config::signal_cuts);

	// Pass 1 for stats
	Double_t num = 0, avg = 0, err = 0;
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

		Float_t* mass = static_cast<Float_t*>(tmva_helper.get_branch(config::mass_branch));
		for (Int_t n = 0, N = tree->GetEntriesFast(); n < N; ++n) {
			tree->GetEntry(n);
			if (tmva_helper.eval()) continue;

			// Welford's online algorithm
			++num;
			double del_1 = *mass - avg;
			avg += del_1 / num;
			double del_2 = *mass - avg;
			err += del_2 * del_1;
		}
	}
	err = sqrt(err / num);

	// Freedman-Diaconis rule
	Double_t bin_width = 3.49 * err / pow(num, 0.3333);
	Double_t num_bins = (config::max_mass - config::min_mass) / bin_width;

	// Pass 2 for histogram
	TH1D fit_hist (
		"mass_fit_hist", (config::channel + ";GeV;Counts").c_str(),
		num_bins, config::min_mass, config::max_mass
	);
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

		tree->SetBranchStatus("*", 0);
		tree->SetBranchStatus(config::mass_branch.c_str(), 1);
		Float_t* mass = static_cast<Float_t*>(tmva_helper.get_branch(config::mass_branch));
		for (Int_t n = 0, N = tree->GetEntriesFast(); n < N; ++n) {
			tree->GetEntry(n);
			if (tmva_helper.eval()) continue;

			// Fill
			fit_hist.Fill(*mass);
		}
	}

	// Fit
	TF1 fit_func(
		(config::channel + "_mass_fit").c_str(),
		"gausn(0) + gausn(3)",
		config::min_mass, config::max_mass
	);

	fit_func.SetParameter(0, 0.5 * num * bin_width);
	fit_func.SetParameter(1, avg);
	fit_func.SetParameter(2, 0.5 * err);
	fit_func.SetParameter(3, 0.5 * num * bin_width);
	fit_func.SetParameter(4, avg);
	fit_func.SetParameter(5, 2.0 * err);

	fit_hist.Fit(&fit_func, "L0");
	Double_t norm = fit_func.GetParameter(0) + fit_func.GetParameter(3);

	config::num0 = fit_func.GetParameter(0) / norm;
	config::sig0 = fit_func.GetParameter(2);
	config::num1 = fit_func.GetParameter(3) / norm;
	config::sig1 = fit_func.GetParameter(5);
	config::mean = config::num0 * fit_func.GetParameter(1) + config::num1 * fit_func.GetParameter(4);

	err = sqrt(config::num0 * config::sig0 * config::sig0 + config::num1 * config::sig1 * config::sig1);

	std::cout
		<< "mean: " << config::mean << "\n"
		<< "\tnum0: " << config::num0 << "\tsig0: " << config::sig0 << "\n"
		<< "\tnum1: " << config::num1 << "\tsig1: " << config::sig1 << "\n"
		<< "min: " << config::mean - 6.0 * err << "\n"
		<< "max: " << config::mean + 6.0 * err << "\n"
		<< std::flush;

	// Draw
	SetsPhenixStyle();
	TCanvas cnvs(
		(config::channel + "_sig_fit_cnvs").c_str(),
		(config::channel + "_sig_fit_cnvs").c_str(),
		600, 800
	);
	cnvs.cd();

	fit_hist.SetLineColor(kBlue);
	fit_func.SetLineColor(kRed);

	fit_hist.Draw();
	fit_func.Draw("same");

	TLine line;
	line.DrawLine(config::mean - 6.0 * err, 0.0, config::mean - 6.0 * err, 1.1 * fit_hist.GetBinContent(fit_hist.GetMaximumBin()));
	line.DrawLine(config::mean - 3.0 * err, 0.0, config::mean - 3.0 * err, 1.1 * fit_hist.GetBinContent(fit_hist.GetMaximumBin()));
	line.DrawLine(config::mean + 3.0 * err, 0.0, config::mean + 3.0 * err, 1.1 * fit_hist.GetBinContent(fit_hist.GetMaximumBin()));
	line.DrawLine(config::mean + 6.0 * err, 0.0, config::mean + 6.0 * err, 1.1 * fit_hist.GetBinContent(fit_hist.GetMaximumBin()));

	cnvs.Update();
	cnvs.Show();
	cnvs.SaveAs((config::channel + "_sig_fit.png").c_str());
}

#endif//FIT_C
