#ifndef FIT_C
#define FIT_C

#include "config.C"

#include <tmvahelper/TMVAHelper.h>
R__LOAD_LIBRARY(libtmvahelper.so)

void
fit (
	std::vector<std::string> const& signal_files = {
		"outputKFP_D0_Kpi_0.root",
		"outputKFP_D0_Kpi_1.root",
		"outputKFP_D0_Kpi_2.root",
		"outputKFP_D0_Kpi_3.root",
		"outputKFP_D0_Kpi_4.root",
		"outputKFP_D0_Kpi_5.root",
		"outputKFP_D0_Kpi_6.root",
		"outputKFP_D0_Kpi_7.root",
		"outputKFP_D0_Kpi_8.root",
		"outputKFP_D0_Kpi_9.root",
	}
) {
	// Helper
	TMVAHelper tmva_helper;
	tmva_helper.read_branches(config::branches);
	tmva_helper.read_training(config::training);
	// tmva_helper.read_cuts(config::pT_cuts[0]); // CHANGE ME

	Long64_t pdf_size = 0;
	std::map<Float_t, Long64_t> pdf;
	for (auto const& signal_file : signal_files) {
		TTree* tree = tmva_helper.get_tree(config::data_dir + "/" + signal_file, "DecayTree");
		if (!tree || tmva_helper.branch(tree)) {
			std::cerr << (config::data_dir + "/" + signal_file) << std::endl;
			continue;
		}

		Float_t* mass = tmva_helper.get_branch(config::mass_branch);
		for (Int_t n = 0, N = tree->GetEntriesFast(); n < N; ++n) {
			tree->GetEntry(n);
			++pdf[*mass];
			++pdf_size;
		}
	}

	// pdf...
	Long64_t counts = 0;
	Float_t quartiles[5] = {};

	for (auto const& [mass_val, count] : pdf) {
		counts += count;
		for (int i = 0; i < 5; ++i)
			if (counts < 0.25 * i * pdf_size) quartiles[i] = mass_val;
	}

	// Freedman-Diaconis rule
	Float_t bin_width = 2.59 * (quartiles[3] - quartiles[1]) / pow(pdf_size, 0.3333);
	Float_t lower = quartiles[2] - 2.5 * (quartiles[3] - quartiles[1]);
	Float_t upper = quartiles[2] + 2.5 * (quartiles[3] - quartiles[1]);
	Int_t num_bins = (upper - lower) / bin_width;
	for (auto quartile : quartiles) {
		std::cout << quartile << std::endl;
	}

	// fill hist
	TH1D* fit_hist = new TH1D (
		"mass_fit_hist", "mass_fit_hist",
		num_bins, lower, upper
	);
	for (auto const& [mass_val, count] : pdf) {
		fit_hist->Fill(mass_val);
	}

	// fit hist
	fit_hist->Fit("gausn", "+");
	TF1* fit_func = dynamic_cast<TF1*>(fit_hist->GetListOfFunctions()->FindObject("gausn"));
	if (!fit_func) {
		std::cerr << "func" << std::endl;
		return;
	}

	config::mean =  fit_func->GetParameter(1);
	config::sigma = fit_func->GetParameter(2);
}

#endif//FIT_C
