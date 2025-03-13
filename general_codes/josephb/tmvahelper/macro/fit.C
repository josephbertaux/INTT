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
	tmva_helper.read_cuts(config::cuts);

	TTree* signal_tree = TMVAHelper::get_tree("signal.root", "DecayTree");
	if (!signal_tree) {
		std::cerr << "expected file 'signal.root' not present" << std::endl;
		return;
	}
	tmva_helper.branch(signal_tree);

	// Pass 1 for stats
	Double_t num = 0, avg = 0, err = 0, min = std::numeric_limits<Float_t>::max(), max = -std::numeric_limits<Float_t>::max();
	Float_t* mass = static_cast<Float_t*>(tmva_helper.get_branch(config::mass_branch));
	for (Int_t n = 0, N = signal_tree->GetEntriesFast(); n < N; ++n) {
		signal_tree->GetEntry(n);
		if (tmva_helper.eval()) continue;

		if (*mass < min) min = *mass;
		if (max < *mass) max = *mass;

		// Welford's online algorithm
		++num;
		double del_1 = *mass - avg;
		avg += del_1 / num;
		double del_2 = *mass - avg;
		err += del_2 * del_1;
	}
	err = sqrt(err / num);

	// Freedman-Diaconis rule
	Double_t bin_width = 3.49 * err / pow(num, 0.3333); // Freedman-Diaconis rule
	Double_t num_bins = (max - min) / bin_width;

	// Pass 2 for histogram
	TH1D fit_hist (
		"mass_fit_hist", (config::channel + ";GeV;Counts").c_str(),
		num_bins, min, max
	);
	for (Int_t n = 0, N = signal_tree->GetEntriesFast(); n < N; ++n) {
		signal_tree->GetEntry(n);
		if (tmva_helper.eval()) continue;

		// Fill
		fit_hist.Fill(*mass);
	}

	// Fit
	TF1 fit_func(
		(config::channel + "_mass_fit").c_str(), "gausn(0)",
		min, max
	);

	fit_func.SetParameter(0, num * bin_width);
	fit_func.SetParameter(1, avg);
	fit_func.SetParameter(2, err);

	fit_hist.Fit(&fit_func, "LR0");
	config::mean = fit_func.GetParameter(1);
	config::sig  = fit_func.GetParameter(2);

	std::cout
		<< "mean: " << config::mean << "\n"
		<< "sig:  " << config::sig  << "\n"
		<< std::flush;

	std::cout << config::get_sideband_cut().GetTitle() << std::endl;

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
	line.DrawLine(config::mean - 6.0 * config::sig, 0.0, config::mean - 6.0 * config::sig, 1.1 * fit_hist.GetBinContent(fit_hist.GetMaximumBin()));
	line.DrawLine(config::mean - 3.0 * config::sig, 0.0, config::mean - 3.0 * config::sig, 1.1 * fit_hist.GetBinContent(fit_hist.GetMaximumBin()));
	line.DrawLine(config::mean + 3.0 * config::sig, 0.0, config::mean + 3.0 * config::sig, 1.1 * fit_hist.GetBinContent(fit_hist.GetMaximumBin()));
	line.DrawLine(config::mean + 6.0 * config::sig, 0.0, config::mean + 6.0 * config::sig, 1.1 * fit_hist.GetBinContent(fit_hist.GetMaximumBin()));

	cnvs.Update();
	cnvs.Show();
	cnvs.SaveAs((config::channel + "_sig_fit.png").c_str());
}

#endif//FIT_C
