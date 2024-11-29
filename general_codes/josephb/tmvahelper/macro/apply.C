#ifndef APPLY_C
#define APPLY_C

#include "config.C"

#include <tmvahelper/TMVAHelper.h>
R__LOAD_LIBRARY(libtmvahelper.so)

void
apply (
	std::vector<std::string> const& input_files = {
		// "outputMinBiasKFParticle_Lc_pKpi_0.root",
		"outputMinBiasKFParticle_D2Kpi_0.root",
		"outputMinBiasKFParticle_D2Kpi_1.root",
		"outputMinBiasKFParticle_D2Kpi_2.root",
		"outputMinBiasKFParticle_D2Kpi_3.root",
		"outputMinBiasKFParticle_D2Kpi_4.root",
		"outputMinBiasKFParticle_D2Kpi_5.root",
		"outputMinBiasKFParticle_D2Kpi_6.root",
	}
) {
	// Helper
	TMVAHelper tmva_helper;
	tmva_helper.read_branches(config::branches);
	tmva_helper.read_training(config::training);

	// Initialize reader
	TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");
	tmva_helper.branch(reader);
	reader->BookMVA (
		config::method_name,
		(boost::format("dataloader/weights/factory_%s.weights.xml") % config::method_name.c_str()).str().c_str()
	);

	Long64_t ref_pdf_size = 0, pdf_size = 0;
	std::map<Float_t, Long64_t> ref_pdf, pdf;
	for (auto const& input_file : input_files) {
		TTree* tree = tmva_helper.get_tree(config::data_dir + "/" + input_file, "DecayTree");
		if (!tree || tmva_helper.branch(tree)) {
			std::cerr << (config::data_dir + "/" + input_file) << std::endl;
			continue;
		}

		Float_t* mass = tmva_helper.get_branch(config::mass_branch);
		for (Int_t n = 0, N = tree->GetEntriesFast(); n < N; ++n) {
			tree->GetEntry(n);

			if (tmva_helper.eval()) continue;

			++ref_pdf[*mass];
			++ref_pdf_size;
			if (reader->EvaluateMVA(config::method_name.c_str()) < config::cut_val) continue;

			++pdf[*mass];
			++pdf_size;
		}
	}

	std::cout << "pdf_size: " << pdf_size << std::endl;
	std::cout << "ref_pdf_size: " << ref_pdf_size << std::endl;

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
		"mass_hist", "mass_hist",
		num_bins, lower, upper
	);
	fit_hist->SetLineColor(kRed);
	for (auto const& [mass_val, count] : pdf) {
		int bin = fit_hist->FindBin(mass_val);
		fit_hist->AddBinContent(bin, count);
	}

	// reference
	TH1D* ref_fit_hist = new TH1D (
		"ref_mass_hist", "ref_mass_hist",
		num_bins, lower, upper
	);
	ref_fit_hist->SetLineColor(kBlue);
	for (auto const& [mass_val, count] : ref_pdf) {
		int bin = ref_fit_hist->FindBin(mass_val);
		ref_fit_hist->AddBinContent(bin, count);
	}

	TCanvas* cnvs = new TCanvas (
		"cnvs", "cnvs", 800, 600
	);
	cnvs->cd();
	// cnvs->SetLogy();

	ref_fit_hist->Draw();
	fit_hist->Draw("same");

	cnvs->Update();
	cnvs->SaveAs("png/cnvs.png");
	// delete cnvs;
}

#endif//APPLY_C

