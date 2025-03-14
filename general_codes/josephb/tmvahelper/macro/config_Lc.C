#ifndef CONFIG_C
#define CONFIG_C

#include <sPhenixStyle.C>
#include <boost/format.hpp>

namespace config {

	// gen
	std::string const channel            = "Lc_pKpi";
	std::string const pythia_config_file = "steering_cards/pythia8_Lc_pKpi.cfg";
	std::string const evtgen_config_file = "dec_files/Lc_pKpi.DEC";
	std::string const decay_descriptor   = "[Lambda_c+ -> proton^+ K^- pi^+]cc";
	int         const particle_trigger   = 4122;

	// fit
	Double_t const min_mass = 1.1; // 1.8; // 
	Double_t const max_mass = 3.7; // 2.6; // 

	Double_t mean = 2.27653e+00;
	Double_t sig  = 2.00000e-01;
	Int_t num_bins = 20;

	// training/application
	// Double_t const ratio = 0.015;
	Double_t const ratio = 2.0E-3;
	std::string const factory_file_name = "factory/factory.root";

	TMVA::Types::EMVA const method_type = TMVA::Types::kBDT;
	std::string       const method_name = "BDT";
	std::string       const method_options =
		"!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20";
	// cut: 0.312493 significance: 0.00769632
	Double_t cut_val =  0.312493;

	// TMVA::Types::EMVA const method_type = TMVA::Types::kMLP;
	// std::string       const method_name = "MLP";
	// std::string       const method_options =
	// 	"!H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator";
	// Double_t cut_val = 0.0;

	std::string const particle_name = "Lambda_cplus";
	std::string const mass_branch = particle_name + "_mass";

	std::vector<std::string> const branches = {
		mass_branch, particle_name + "_pT",

		particle_name + "_vertex_volume",
		particle_name + "_DIRA",

		particle_name + "_decayLength", particle_name + "_decayLengthErr",

		// particle_name +  "_x", particle_name +  "_y",
		// particle_name + "_px", particle_name + "_py",

		"track_1_pT", "track_1_pTErr",
		"track_2_pT", "track_2_pTErr",
		"track_3_pT", "track_3_pTErr",

		// "track_1_PDG_ID/I",
		// "track_2_PDG_ID/I",
		// "track_3_PDG_ID/I",
	};
	
	std::vector<std::string> const training = {
		// decay_length_significance
		std::string{"decay_len_sig:=log(abs("} + particle_name + "_decayLength / " + particle_name + "_decayLengthErr" + "))",

		// pointing angle
		std::string{"alpha:="} + "acos(" + particle_name + "_DIRA" + ")",

		// vertex volume
		std::string{"vertex_volume:="} + particle_name + "_vertex_volume",

		// std::string{"alpha:=acos("}
		// 	+        "(" + particle_name + "_px * " + particle_name +  "_x + " + particle_name + "_py * " + particle_name +  "_y)"
		// 	+ " / sqrt(" + particle_name +  "_x * " + particle_name +  "_x + " + particle_name +  "_y * " + particle_name +  "_y)"
		// 	+ " / sqrt(" + particle_name + "_px * " + particle_name + "_px + " + particle_name + "_py * " + particle_name + "_py)"
		// 	+ ")",

		std::string{"track_1_pT_sig:=(track_1_pT / track_1_pTErr)"},
		std::string{"track_2_pT_sig:=(track_2_pT / track_2_pTErr)"},
		std::string{"track_3_pT_sig:=(track_3_pT / track_3_pTErr)"},

		// // kaon pT significance
		// std::string{"kaon_pT_sig:="}
		// 	+ "(abs(track_1_PDG_ID) == 321) * (track_1_pT / track_1_pTErr) + "
		// 	+ "(abs(track_2_PDG_ID) == 321) * (track_2_pT / track_2_pTErr) + "
		// 	+ "(abs(track_3_PDG_ID) == 321) * (track_3_pT / track_3_pTErr)",
		// 	// ...

		// // pion pT significance
		// std::string{"pion_pT_sig:="}
		// 	+ "(abs(track_1_PDG_ID) == 211) * (track_1_pT / track_1_pTErr) + "
		// 	+ "(abs(track_2_PDG_ID) == 211) * (track_2_pT / track_2_pTErr) + "
		// 	+ "(abs(track_3_PDG_ID) == 211) * (track_3_pT / track_3_pTErr)",
		// 	// ...

		// // proton pT significance
		// std::string{"proton_pT_sig:="}
		// 	+ "(abs(track_1_PDG_ID) == 2212) * (track_1_pT / track_1_pTErr) + "
		// 	+ "(abs(track_2_PDG_ID) == 2212) * (track_2_pT / track_2_pTErr) + "
		// 	+ "(abs(track_3_PDG_ID) == 2212) * (track_3_pT / track_3_pTErr)",
		// 	// ...
		// // ...
	};
	
	std::vector<std::string> const cuts = {
		// "2.5 < " + particle_name + "_pT && " + particle_name + "_pT < 3.5",
	};

	TCut get_sideband_cut () {
		std::string cut = std::string{"("}
			+ std::to_string(mean - 6.0 * sig) + " < " + mass_branch + " && " + mass_branch + " < " + std::to_string(mean - 3.0 * sig)
			+ ") || ("
			+ std::to_string(mean + 3.0 * sig) + " < " + mass_branch + " && " + mass_branch + " < " + std::to_string(mean + 6.0 * sig)
			+ ")";
		return TCut(cut.c_str());
	}

	Double_t get_bin_width (
		Double_t num
	) {
		Double_t bin_width = 3.49 * sig / pow(num, 0.3333); // Freedman-Diaconis rule
		return (max_mass - min_mass) / bin_width;
	}

};

#endif//CONFIG_C

