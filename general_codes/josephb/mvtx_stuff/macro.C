#ifndef MACRO_C
#define MACRO_C

#include <sPhenixStyle.C>

Double_t phi_min = -1.0 * TMath::Pi() / 2.0;
int multiplicity_cutoff = 10000;
int chip_hits_cutoff = 500;
int clus_cutoff = 20;
Double_t width_cutoff = 0.2;

template <class T> 
void draw_canvas (int, int, std::vector<T*>);

std::vector<std::set<Int_t>>* g_clusters{nullptr};
std::set<Int_t>* g_clus{nullptr};
TH1* g_hist{nullptr};

Double_t
prob_func (
	Double_t* vars,
	Double_t* pars
) {
	Int_t bin = g_hist->FindBin(vars[0], vars[1]);
	if (g_clus->find(bin) == g_clus->end()) return 0;

	// bool is_in_big_enough_clus = false;
	// for (auto const& clus : *g_clusters) {
	// 	if (clus.find(bin) == clus.end()) continue;
	// 	is_in_big_enough_clus = true;
	// 	break;
	// }
	// if (!is_in_big_enough_clus) return 0;

	Double_t dist = pars[0] + pars[1] * vars[0] - vars[1];
	dist /= sqrt(1.0 + pars[1] * pars[1]);
	return pars[3] * TMath::Gaus(dist, 0.0, pars[2]);
}

template <class T> 
void
make_cluster (
	int start_bin,
	std::map<Int_t, bool>& added,
	std::set<Int_t>& cluster,
	T* hist
) {
	if (added[start_bin]) return;

	if (!hist->GetBinContent(start_bin)) return;

	cluster.insert(start_bin);
	added[start_bin] = true;

	Int_t bin_x, bin_y, bin_z;
	hist->GetBinXYZ(start_bin, bin_x, bin_y, bin_z);

	// Z coords are unnused

	if (bin_x == 0) return; // underflow
	if (bin_y == 0) return; // underflow

	if (bin_x == hist->GetNbinsX() + 1) return; // overflow
	if (bin_y == hist->GetNbinsY() + 1) return; // overflow

	for (int dx = -1; dx < 2; dx += 1) {
		for (int dy = -1; dy < 2; dy += 1) {
			start_bin = hist->GetBin(bin_x + dx, bin_y + dy, bin_z);
			make_cluster (start_bin, added, cluster, hist);
		}
	}
}

void
macro (
	std::string const& file_name =
		"/gpfs/mnt/gpfs02/sphenix/user/cdean/public/mvtx_standalone_cluster/macros/outputHits_00054271.root",
	int runnumber = 54271
) {
	SetsPhenixStyle();
	gROOT->SetBatch();

	TFile* file = TFile::Open(file_name.c_str(), "READ");
	if (!file) {
		std::cerr << "file" << std::endl;
		return;
	}

	TTree* tree = dynamic_cast<TTree*>(file->Get("Hits"));
	if (!tree) {
		std::cerr << "tree" << std::endl;
		return;
	}

	Int_t event =                     0;                        tree->SetBranchAddress("event",          &event);
	Int_t layer =                     0;                        tree->SetBranchAddress("layer",          &layer);
	Int_t stave =                     0;                        tree->SetBranchAddress("stave",          &stave);
	Int_t chip =                      0;                        tree->SetBranchAddress("chip",           &chip);
	Int_t chip_hits =                 0;                        tree->SetBranchAddress("chip_hits",      &chip_hits);

	Float_t chip_occupancy =          0;                        tree->SetBranchAddress("chip_occupancy", &chip_occupancy);

	std::vector<Float_t>* localX =    new std::vector<Float_t>; tree->SetBranchAddress("localX",         &localX);
	std::vector<Float_t>* localY =    new std::vector<Float_t>; tree->SetBranchAddress("localY",         &localY);
	std::vector<Float_t>* globalX =   new std::vector<Float_t>; tree->SetBranchAddress("globalX",        &globalX);
	std::vector<Float_t>* globalY =   new std::vector<Float_t>; tree->SetBranchAddress("globalY",        &globalY);
	std::vector<Float_t>* globalZ =   new std::vector<Float_t>; tree->SetBranchAddress("globalZ",        &globalZ);
	std::vector<Float_t>* globalPhi = new std::vector<Float_t>; tree->SetBranchAddress("globalPhi",      &globalPhi);

	int n_bins_phi = 100, n_bins_z = 100;
	std::map<int, std::vector<TH1*>> prof_map;

	int prev_event = 0;
	int multiplicity = 0;

	tree->SetBranchStatus("*", 0);
	tree->SetBranchStatus("event", 1);
	tree->SetBranchStatus("chip_hits", 1);
	for (int n = 0, N = tree->GetEntriesFast(); n < N; ++n) {
		tree->GetEntry(n);

		if (prev_event == event) {
			// multiplicity += chip_hits;
			if (multiplicity < chip_hits) multiplicity = chip_hits;
			continue;
		}

		if (multiplicity < multiplicity_cutoff) {
			prev_event = event;
			multiplicity = chip_hits;
			continue;
		}

		std::cout << "\t" << prev_event << std::endl;

		prof_map[prev_event] = {

			new TH2D (Form("event_%d_layer_%d_hist", event, 0), "MVTX Layer 0; Layer 0    Z (cm);Layer 0    Phi (Radians)",
					n_bins_z, -13.5, 13.5, n_bins_phi, phi_min, phi_min + 2.0 * TMath::Pi()),
			new TH2D (Form("event_%d_layer_%d_hist", event, 1), "MVTX Layer 1; Layer 1    Z (cm);Layer 1    Phi (Radians)",
					n_bins_z, -13.5, 13.5, n_bins_phi, phi_min, phi_min + 2.0 * TMath::Pi()),
			new TH2D (Form("event_%d_layer_%d_hist", event, 2), "MVTX Layer 2; Layer 2    Z (cm);Layer 2    Phi (Radians)",
					n_bins_z, -13.5, 13.5, n_bins_phi, phi_min, phi_min + 2.0 * TMath::Pi()),

			new TProfile (Form("event_%d_layer_%d_prof", event, 0), "MVTX Layer 0; Layer 0    Z (cm);Layer 0    Phi (Radians)",
					n_bins_z, -13.5, 13.5, phi_min, phi_min + 2.0 * TMath::Pi()),
			new TProfile (Form("event_%d_layer_%d_prof", event, 1), "MVTX Layer 1; Layer 1    Z (cm);Layer 1    Phi (Radians)",
					n_bins_z, -13.5, 13.5, phi_min, phi_min + 2.0 * TMath::Pi()),
			new TProfile (Form("event_%d_layer_%d_prof", event, 2), "MVTX Layer 2; Layer 2    Z (cm);Layer 2    Phi (Radians)",
					n_bins_z, -13.5, 13.5, phi_min, phi_min + 2.0 * TMath::Pi()),

			new TH2D (Form("event_%d_layer_%d_hist_entire", event, 0), "MVTX Layer 0; Layer 0    Z (cm);Layer 0    Phi (Radians)",
					n_bins_z, -13.5, 13.5, n_bins_phi, phi_min, phi_min + 2.0 * TMath::Pi()),
			new TH2D (Form("event_%d_layer_%d_hist_entire", event, 1), "MVTX Layer 1; Layer 1    Z (cm);Layer 1    Phi (Radians)",
					n_bins_z, -13.5, 13.5, n_bins_phi, phi_min, phi_min + 2.0 * TMath::Pi()),
			new TH2D (Form("event_%d_layer_%d_hist_entire", event, 2), "MVTX Layer 2; Layer 2    Z (cm);Layer 2    Phi (Radians)",
					n_bins_z, -13.5, 13.5, n_bins_phi, phi_min, phi_min + 2.0 * TMath::Pi()),
		};

		for (int i = 3; i < 6; ++i) {
			prof_map[prev_event][i]->Sumw2();
		}

		prev_event = event;
		multiplicity = chip_hits;
	}

	std::cout << prof_map.size() << std::endl;

	tree->SetBranchStatus("layer", 1);
	tree->SetBranchStatus("globalZ", 1);
	tree->SetBranchStatus("globalPhi", 1);
	for (int n = 0, N = tree->GetEntriesFast(); n < N; ++n) {
		tree->GetEntry(n);

		auto itr = prof_map.find(event);
		if (itr == prof_map.end()) continue;

		for (int h = 0; h < chip_hits; ++h) {

			double phi = globalPhi->at(h);
			while (phi < phi_min) {
				phi += 2.0 * TMath::Pi();
			}

			itr->second[layer + 6]->Fill(globalZ->at(h), phi);

			if (chip_hits < chip_hits_cutoff) continue;
			itr->second[layer]->Fill(globalZ->at(h), phi);
			itr->second[layer + 3]->Fill(globalZ->at(h), phi);
		}
	}

	// int i = 0;
	for (auto& [evt, vec] : prof_map) {
		// if (++i < 2) continue;
		draw_canvas(runnumber, evt, vec);
		// break;
	}

	gSystem->Exit(0);
}

template <class T>
void
draw_canvas (
	int runnumber,
	int event,
	std::vector<T*> hist
) {
	gStyle->SetOptStat(0);

	TCanvas* cnvs = new TCanvas(Form("cnvs%d", event), Form("cnvs%d", event), 1600, 900);
	cnvs->Range(0.0, 0.0, 1.0, 1.0);
	cnvs->cd();
	cnvs->Draw();

	Double_t max = 0.0;
	for (int i = 0; i < 3; ++i) {
		Double_t d = hist[i]->GetBinContent(hist[i]->GetMaximumBin());
		if (d < max) continue;
		max = d;
	}

	for (auto& hist_ptr : hist) {
		hist_ptr->SetTitleSize(0.06);
		hist_ptr->SetTitleOffset(0.2);

		hist_ptr->GetXaxis()->SetTitleSize(0.06);
		hist_ptr->GetXaxis()->SetTitleOffset(0.8);
		hist_ptr->GetXaxis()->CenterTitle(kTRUE);

		hist_ptr->GetYaxis()->SetTitleSize(0.06);
		hist_ptr->GetYaxis()->SetTitleOffset(0.4);
		hist_ptr->GetYaxis()->CenterTitle(kTRUE);

		hist_ptr->GetZaxis()->SetRangeUser(1, max);
	}

	for (int i = 0; i < 3; ++i) {
		cnvs->cd();
		TPad* pad = new TPad (
			Form("left_pad%d", i), Form("left_pad%d", i),
			0.0, 0.3 * i + 0.0,
			0.5, 0.3 * i + 0.3
		);
		pad->SetFillStyle(4000);
		pad->Range(0.0, 0.0, 1.0, 1.0);
		pad->SetBottomMargin(0.15);
		pad->SetTopMargin(0.0);
		pad->SetLeftMargin(0.05);
		pad->SetRightMargin(0.1);
		pad->SetLogz();
		pad->Draw();

		pad->cd();

		hist[i + 6]->Draw("COL");

	}

	for (int i = 0; i < 3; ++i) {
		cnvs->cd();
		TPad* pad = new TPad (
			Form("right_pad%d", i), Form("right_pad%d", i),
			0.5, 0.3 * i + 0.0,
			1.0, 0.3 * i + 0.3
		);
		pad->SetFillStyle(4000);
		pad->Range(0.0, 0.0, 1.0, 1.0);
		pad->SetBottomMargin(0.15);
		pad->SetTopMargin(0.0);
		pad->SetLeftMargin(0.05);
		pad->SetRightMargin(0.1);
		pad->SetLogz();
		pad->Draw();

		pad->cd();

		hist[i]->Draw("COLZ");
		// hist[i + 3]->Draw("same");

		if (!hist[i]->GetEntries() || !hist[i+3]->GetEntries()) {
			std::cout << "TH2D:  " << event << " " << i << " " << hist[i + 0]->GetEntries() << std::endl;
			std::cout << "TProf: " << event << " " << i << " " << hist[i + 3]->GetEntries() << std::endl;
			continue;
		}

		if (!hist[i]->GetEntries() || !hist[i+3]->GetEntries()) {
			std::cout << "TH2D:  " << event << " " << i << " " << hist[i + 0]->GetEntries() << std::endl;
			std::cout << "TProf: " << event << " " << i << " " << hist[i + 3]->GetEntries() << std::endl;
			continue;
		}

		// clustering

		std::map<Int_t, bool> added;
		std::vector<std::set<Int_t>> clusters;
		for (int n = 0; n < hist[i]->GetNcells(); ++n) {
			if (added[n]) continue;

			clusters.push_back({});
			make_cluster (n, added, clusters.back(), hist[i]);
		}

		std::set<Int_t>& max_clus = clusters.back();
		for (auto const& clus : clusters) {
			if(max_clus.size() < clus.size()) max_clus = clus;
		}

		// clusters.erase(std::remove_if (
		// 	clusters.begin(), clusters.end(),
		// 	[](std::set<Int_t> const& clus) {return clus.size() < clus_cutoff;}
		// ), clusters.end());

		// max_clus = clusters.back();
		// for (auto const& clus : clusters) {
		// 	if(max_clus.size() < clus.size()) max_clus = clus;
		// }

		for (int n = 0; n < hist[i]->GetNcells(); ++n) {

			// Set bins not in these clusters to 0
			// bool is_in_big_enough_cluster = false;
			// for (auto const& clus : clusters) {
			// 	if (clus.find(n) != clus.end()) continue;
			// 	is_in_big_enough_cluster = true;
			// 	break;
			// }

			// if (is_in_big_enough_cluster) continue;

			// Set bins not in max_clus to 0
			if (max_clus.find(n) != max_clus.end()) continue;
			hist[i]->SetBinContent(n, 0);
		}

		TF1* line = new TF1("line", "[0] + [1] * x", -13.5, 13.5);
		line->SetParameter(0, 0.0);
		line->SetParameter(1, 0.0);
		hist[i + 3]->Fit(line, "WQN0");

		g_clusters = &clusters;
		g_clus = &max_clus;
		g_hist = hist[i];

		TF2* fit_func = new TF2("fit_func", prob_func, -13.5, 13.5, phi_min, phi_min + 2.0 * TMath::Pi(), 4);
		fit_func->SetParameter(0, line->GetParameter(0));
		fit_func->SetParameter(1, line->GetParameter(1));
		fit_func->SetParameter(2, 0.1);
		fit_func->SetParameter(3, chip_hits_cutoff);

		fit_func->SetParLimits(2, 0.0, width_cutoff);
		hist[i]->Fit(fit_func, "LQN0");
		// fit_func->Draw("cont1 same");

		line->SetParameter(0, fit_func->GetParameter(0));
		line->SetParameter(1, fit_func->GetParameter(1));
		line->Draw("same");

		// Double_t y_min = phi_min;
		// Double_t y_max = phi_min + 2.0 * TMath::Pi();
		Double_t y_min = fit_func->GetParameter(0) - 1.0; // 3.0 * fit_func->GetParameter(2) - 13.5 * abs(fit_func->GetParameter(1));
		Double_t y_max = fit_func->GetParameter(0) + 1.0; // 3.0 * fit_func->GetParameter(2) + 13.5 * abs(fit_func->GetParameter(1));

		hist[i + 0]->GetYaxis()->SetRangeUser(y_min, y_max);
		hist[i + 3]->GetYaxis()->SetRangeUser(y_min, y_max);

		Double_t pos = (y_max + y_min) * 0.5 + (y_max - y_min) * 3.0 / 8.0;
		TPaveText* text = new TPaveText (
			- 2, pos - (y_max - y_min) * 1.0 / 8.0,
			+12, pos + (y_max - y_min) * 1.0 / 8.0
		);
		text->SetTextSize(0.06);
		text->SetTextAlign(12);

		text->AddText (
			Form("m = %.2E #pm %.2E rad/cm", fit_func->GetParameter(1), fit_func->GetParError(1))
		);
		text->AddText (
			Form("#sigma = %.2E #pm %.2E rad (limit %.2E)", fit_func->GetParameter(2), fit_func->GetParError(2), width_cutoff)
		);

		Double_t error   = 3.0 * fit_func->GetParameter(2) / 13.5;
		Double_t err_err = 3.0 * fit_func->GetParError(2) / 13.5;
		text->AddText (
			Form("3#sigma / (13.5 cm) = %.2E #pm %.2E rad / cm", error, err_err)
		);

		text->Draw();
	}

	bool wtf = true;
	for (int i = 0; i < 3; ++i) {
		if (!hist[i]->GetEntries()) continue;
		wtf = false;
	}

	if (wtf) {
		std::cout << event << ": skipping" << std::endl;
		delete cnvs;
		return;
	}

	cnvs->cd();
	TPad* title_pad = new TPad (
		"title_pad", "title_pad",
		0.0, 0.9, 1.0, 1.0
	);
	title_pad->SetFillStyle(4000);
	title_pad->Range(0.0, 0.0, 1.0, 1.0);
	title_pad->Draw();

	title_pad->cd();
	TText title_text;
	title_text.SetTextAlign(22);
	title_text.SetTextSize(0.25);
	title_text.DrawText(0.5, 0.80, Form("Run %08d MVTX Z-Phi Occupancy For Event %08d", runnumber, event));
	title_text.SetTextSize(0.15);
	title_text.DrawText(0.5, 0.60, "Fit PDF(z, phi) as being Gaussian in distance from line phi = m * z + b");
	title_text.DrawText(0.5, 0.40, "Likelihood fit using only non-empty bins");

	title_text.DrawText(0.25, 0.25, "Entire event");
	title_text.DrawText(0.75, 0.25, "Only fitted bins; zoomed phi range");

	cnvs->Update();
	cnvs->SaveAs(Form("mvtx_occupancy_run%08d_event%08d.png", runnumber, event));
	cnvs->Draw();
	cnvs->Close();
	delete cnvs;
}

#endif//MACRO_C

