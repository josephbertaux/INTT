#ifndef DRAW_SOMETHING_C
#define DRAW_SOMETHING_C

#include <filesystem>

void
DrawSomething (
	std::string const& file_name = "dat/clusters_seeds_196clusters_tracks_41992-23.root_resid.root"
) {
	if (file_name.empty()) {
		std::cerr << __LINE__ << std::endl;
		return;
	}

	if (!std::filesystem::exists(file_name)) {
		std::cerr << __LINE__ << std::endl;
		return;
	}

	TFile* file = TFile::Open(file_name.c_str(), "READ");
	if (!file) {
		std::cerr << __LINE__ << std::endl;
		return;
	}

	TTree* residual_tree = dynamic_cast<TTree*>(file->Get("residualtree"));
	TTree* cluster_tree  = dynamic_cast<TTree*>(file->Get("clustertree"));
	TTree* vertex_tree   = dynamic_cast<TTree*>(file->Get("vertextree"));
	TTree* event_tree    = dynamic_cast<TTree*>(file->Get("eventtree"));

	if (!residual_tree || !cluster_tree || !vertex_tree || !event_tree) {
		std::cerr << __LINE__ << std::endl;
		return;
	}

	Int_t ntracks = 0;
	Float_t vx = 0.0;
	Float_t vy = 0.0;
	Float_t vz = 0.0;
	std::vector<Float_t>* gx = new std::vector<Float_t>;
	std::vector<Float_t>* gy = new std::vector<Float_t>;
	std::vector<Float_t>* gz = new std::vector<Float_t>;

	vertex_tree->SetBranchAddress("ntracks", &ntracks);
	vertex_tree->SetBranchAddress("vx", &vx);
	vertex_tree->SetBranchAddress("vy", &vy);
	vertex_tree->SetBranchAddress("vz", &vz);
	vertex_tree->SetBranchAddress("gx", &gx);
	vertex_tree->SetBranchAddress("gy", &gy);
	vertex_tree->SetBranchAddress("gz", &gz);

	for (Int_t n = 0, N = vertex_tree->GetEntriesFast(); n < N; ++n) {
		vertex_tree->GetEntry(n);

		// if (n != 314) continue;

		// std::cout << std::endl;
		// std::cout << "entry:    " << n << std::endl;
		// std::cout << "ntracks:  " << ntracks << std::endl;
		// std::cout << "gx->size: " << gx->size() << std::endl;
		// std::cout << "gy->size: " << gy->size() << std::endl;
		// std::cout << "gz->size: " << gz->size() << std::endl;
		// std::cout << "vx:       " << vx << std::endl;
		// std::cout << "vy:       " << vy << std::endl;
		// std::cout << std::endl;

		if (gx->size() != gy->size()) continue;
		int len = gx->size();

		Double_t* vertex_x = new Double_t[1];
		Double_t* vertex_y = new Double_t[1];
		vertex_x[0] = vx;
		vertex_y[0] = vy;

		Double_t* x = new Double_t[len];
		Double_t* y = new Double_t[len];
		for (int i = 0; i < len; ++i) {
			x[i] = gx->at(i);
			y[i] = gy->at(i);
		}

		TCanvas* cnvs = new TCanvas ("cnvs", "cnvs", 1600, 900);
		cnvs->cd();

		TPad* pad = new TPad ("pad", "pad", 0.0, 0.0, 1.0, 1.0);
		pad->Range(0.0, 0.0, 1.0, 1.0);
		pad->SetFillStyle(4000);
		pad->Draw();
		pad->cd();

		TMultiGraph* mg = new TMultiGraph();
		mg->SetTitle(Form("event_%05d", int{n}));

		TGraph* graph = new TGraph (len, x, y);
		graph->SetMarkerColor(kBlack);
		graph->SetMarkerStyle(20);
		graph->SetMarkerSize(0.5);
		mg->Add(graph, "p");

		TGraph* vertex = new TGraph (1, vertex_x, vertex_y);
		vertex->SetMarkerColor(kRed);
		vertex->SetMarkerStyle(20);
		vertex->SetMarkerSize(0.5);
		mg->Add(vertex, "p");

		mg->GetXaxis()->SetLimits(-12.0, 12.0);
		mg->GetXaxis()->SetRangeUser(-12.0, 12.0);
		mg->GetYaxis()->SetRangeUser(-12.0, 12.0);
		mg->Draw("a");

		cnvs->Update();
		cnvs->Show();

		cnvs->SaveAs(Form("png/plot_%05d.png", int{n}));

		delete mg;
		delete cnvs;
		delete[] vertex_x;
		delete[] vertex_y;
		delete[] x;
		delete[] y;
	}

	// ...
}

#endif//DRAW_SOMETHING_C
