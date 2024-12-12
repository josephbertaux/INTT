#ifndef CONFIG_C
#define CONFIG_C

#include <boost/format.hpp>

namespace config {

	// gen
	std::string const channel            = "D0_Kpi"; // "Lc_pKpi";
	std::string const pythia_config_file = "steering_cards/pythia8_D0_Kpi.cfg"; // "steering_cards/pythia8_Lc_pKpi.cfg";
	std::string const evtgen_config_file = "dec_files/D0_Kpi.DEC"; // "dec_files/Lc_pKpi.DEC";
	std::string const decay_descriptor   = "[D0 -> K^- pi^+]cc"; // "[Lambda_c+ -> proton^+ K^- pi^+]cc";
	int         const particle_trigger   = 421; // 4122;

	// fit
	Double_t const min_mass = 1.5; // 1.1; // 1.8;
	Double_t const max_mass = 2.2; // 2.6; // 2.7;

	Double_t mean = 1.85834e+00;
	Double_t num0 = 0.401421;
	Double_t sig0 = 1.43871e-02;
	Double_t num1 = 0.598579;
	Double_t sig1 = 1.48976e-01;
	Int_t num_bins = 100;

	// training/application
	Double_t const ratio = 0.1;
	std::string const factory_file_name = "factory/factory.root";

	TMVA::Types::EMVA const method_type = TMVA::Types::kBDT;
	std::string       const method_name = "BDT";
	std::string       const method_options =
		"!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20";
	Double_t cut_val = 0.0520139;

	// TMVA::Types::EMVA const method_type = TMVA::Types::kMLP;
	// std::string       const method_name = "MLP";
	// std::string       const method_options =
	// 	"!H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator";
	// Double_t cut_val = 0.0;

	std::string const particle_name = "D0"; // "Lambda_cplus";
	std::string const mass_branch = particle_name + "_mass";
	std::vector<std::string> const branches = {
		mass_branch, particle_name + "_pT",

		particle_name + "_decayLength", particle_name + "_decayLengthErr",

		particle_name +  "_x", particle_name +  "_y",
		particle_name + "_px", particle_name + "_py",

		"track_1_pT", "track_1_pTErr",
		"track_2_pT", "track_2_pTErr",
		// "track_3_pT", "track_3_pTErr",

		"track_1_PDG_ID/I",
		"track_2_PDG_ID/I",
		// "track_3_PDG_ID/I",
	};
	
	std::vector<std::string> const training = {
		// decay_length_significance
		std::string{"decay_len_sig:=log(abs("}
			+ particle_name + "_decayLength / " + particle_name + "_decayLengthErr"
			+ "))",

		// pointing angle
		std::string{"alpha:=acos("}
			+        "(" + particle_name + "_px * " + particle_name +  "_x + " + particle_name + "_py * " + particle_name +  "_y)"
			+ " / sqrt(" + particle_name +  "_x * " + particle_name +  "_x + " + particle_name +  "_y * " + particle_name +  "_y)"
			+ " / sqrt(" + particle_name + "_px * " + particle_name + "_px + " + particle_name + "_py * " + particle_name + "_py)"
			+ ")",

		// kaon pT significance
		std::string{"kaon_pT_sig:="}
			+ "(abs(track_1_PDG_ID) == 321) * (track_1_pT / track_1_pTErr) + "
			+ "(abs(track_2_PDG_ID) == 321) * (track_2_pT / track_2_pTErr)",
			// ...

		// pion pT significance
		std::string{"pion_pT_sig:="}
			+ "(abs(track_1_PDG_ID) == 211) * (track_1_pT / track_1_pTErr) + "
			+ "(abs(track_2_PDG_ID) == 211) * (track_2_pT / track_2_pTErr)",
			// ...

		// proton pT significance
		// ...
	};
	
	TCut get_sideband_cut () {
		Double_t err = sqrt(num0 * sig0 * sig0 + num1 * sig1 * sig1);
		std::string cut = std::string{"("}
			+ std::to_string(mean - 6.0 * err) + " < " + mass_branch + " && " + mass_branch + " < " + std::to_string(mean - 3.0 * err)
			+ ") || ("
			+ std::to_string(mean + 3.0 * err) + " < " + mass_branch + " && " + mass_branch + " < " + std::to_string(mean + 6.0 * err)
			+ ")";
		return TCut(cut.c_str());
	}

	Double_t get_bin_width (
		Double_t N
	) {
		Double_t err = sqrt(num0 * sig0 * sig0 + num1 * sig1 * sig1);
		Double_t bin_width = 3.49 * err / pow(N, 0.3333); // Freedman-Diaconis rule
		return (max_mass - min_mass) / bin_width;
	}
};

#endif//CONFIG_C

