#ifndef MACRO_C
#define MACRO_C

#include <sPhenixStyle.C>

template <class T> 
void draw_canvas (int, int, std::vector<T*>);

void
macro (
	std::string const& file_name =
		"/gpfs/mnt/gpfs02/sphenix/user/cdean/public/mvtx_standalone_cluster/macros/outputHits_00054271.root"
	int runnumber = 54271,
	int multiplicity_cutoff = 10000
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

	int n_bins_phi = 60, n_bins_z = 60;
	std::map<int, std::vector<TProfile*>> prof_map;
	std::vector<TH2D*> hist = {
		new TH2D ("hist0", "MVTX Layer 0; Layer 0    Z (cm);Layer 0    Phi (Radians)", n_bins_z, -13.5, 13.5, n_bins_phi, -3.1416, 3.1416),
		new TH2D ("hist1", "MVTX Layer 1; Layer 1    Z (cm);Layer 1    Phi (Radians)", n_bins_z, -13.5, 13.5, n_bins_phi, -3.1416, 3.1416),
		new TH2D ("hist2", "MVTX Layer 2; Layer 2    Z (cm);Layer 2    Phi (Radians)", n_bins_z, -13.5, 13.5, n_bins_phi, -3.1416, 3.1416),
	};

	tree->SetBranchStatus("*", 0);
	tree->SetBranchStatus("event", 1);
	tree->SetBranchStatus("layer", 1);
	tree->SetBranchStatus("chip_hits", 1);
	for (int n = 0, N = tree->GetEntriesFast(); n < N; ++n) {
		tree->GetEntry(n);
		if (chip_hits < multiplicity_cutoff) continue;
		if (prof_map.find(event) != prof_map.end()) continue;
		prof_map[event] = {
			new TProfile (Form("event_%d_layer_%d", event, 0), "MVTX Layer 0; Layer 0    Z (cm);Layer 0    Phi (Radians)", n_bins_z, -13.5, 13.5, -3.1416, 3.1416),
			new TProfile (Form("event_%d_layer_%d", event, 1), "MVTX Layer 1; Layer 1    Z (cm);Layer 1    Phi (Radians)", n_bins_z, -13.5, 13.5, -3.1416, 3.1416),
			new TProfile (Form("event_%d_layer_%d", event, 2), "MVTX Layer 2; Layer 2    Z (cm);Layer 2    Phi (Radians)", n_bins_z, -13.5, 13.5, -3.1416, 3.1416),
		};
	}

	tree->SetBranchStatus("globalZ", 1);
	tree->SetBranchStatus("globalPhi", 1);
	for (int n = 0, N = tree->GetEntriesFast(); n < N; ++n) {
		tree->GetEntry(n);
		if (prof_map.find(event) == prof_map.end()) continue;
		for (int h = 0; h < chip_hits; ++h) {
			hist[layer]->Fill(globalZ->at(h), globalPhi->at(h));
			prof_map[event][layer]->Fill(globalZ->at(h), globalPhi->at(h));
		}
	}

	draw_canvas(runnumber, -1, hist);
	for (auto const& [evt, prof] : prof_map) {
		draw_canvas(runnumber, evt, prof);
	}
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

	for (int i = 0; i < 3; ++i) {
		cnvs->cd();
		TPad* pad = new TPad (
			Form("pad%d", i), Form("pad%d", i),
			0.0, 0.3 * i + 0.0,
			1.0, 0.3 * i + 0.3
		);
		pad->SetFillStyle(4000);
		pad->Range(0.0, 0.0, 1.0, 1.0);
		pad->SetBottomMargin(0.15);
		pad->SetTopMargin(0.0);
		pad->SetLeftMargin(0.03);
		pad->SetRightMargin(0.1);
		pad->SetLogz();
		pad->Draw();

		pad->cd();
		hist[i]->SetTitleSize(0.06);
		hist[i]->SetTitleOffset(0.2);

		hist[i]->GetXaxis()->SetTitleSize(0.06);
		hist[i]->GetXaxis()->SetTitleOffset(0.8);
		hist[i]->GetXaxis()->CenterTitle(kTRUE);

		hist[i]->GetYaxis()->SetTitleSize(0.06);
		hist[i]->GetYaxis()->SetTitleOffset(0.2);
		hist[i]->GetYaxis()->CenterTitle(kTRUE);
		hist[i]->Draw("COLZ");
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
	title_text.SetTextSize(0.3);
	if (event) {
		title_text.DrawText(0.5, 0.60, Form("Run %08d MVTX Z-Phi Occupancy For Event %d", runnumber, event));
		title_text.SetTextSize(0.2);
		title_text.DrawText(0.5, 0.40, Form("(a chip's multiplicity exceeded %d in this event)", multiplicity_cutoff));
	} else {
		title_text.DrawText(0.5, 0.60, Form("Run %08d MVTX Z-Phi High-Multiplicity Event Occupancy", runnumber));
		title_text.SetTextSize(0.2);
		title_text.DrawText(0.5, 0.40, Form("(only and all hits from events where a chip's multiplicity exceeded %d)", multiplicity_cutoff));
	}

	cnvs->Update();
	cnvs->SaveAs(Form("plot_event%d.png", event));
	cnvs->Close();
	delete cnvs;
}

#endif//MACRO_C

