#ifndef DEFS_C
#define DEFS_C

#include <boost/format.hpp>

namespace defs {
	std::string const data_dir =          "/sphenix/user/jbertaux/Repositories/analysis/HF-Particle/KFParticle_sPHENIX/hf_generator/dat";
	std::string const factory_file_name = "factory/factory.root";

	Double_t ratio = 0.1;
	Double_t mean =  1.864e+00;
	Double_t sigma = 1.7e-01;

	TMVA::Types::EMVA const method_type = TMVA::Types::kBDT;
	std::string       const method_name = "BDT";
	std::string       const method_options =
		"!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20";

	// TMVA::Types::EMVA const method_type = TMVA::Types::kMLP;
	// std::string       const method_name = "MLP";
	// std::string       const method_options =
	// 	"!H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator";

	std::string mass_branch = "D0_mass";
	std::vector<std::string> const branches = {
		mass_branch,
		"D0_pT",
		"D0_decayLength",
		"D0_decayLengthErr",
		"D0_DIRA",
		"D0_FDchi2",
		"D0_IP",
		"D0_IPchi2",
		"D0_IPErr",
		"D0_IP_xy",
	};
	
	std::vector<std::string> const training = {
		"LogDLS := log(abs(D0_decayLength / D0_decayLengthErr))",
		"D0_DIRA",
		"D0_FDchi2",
		"D0_IP",
		"D0_IPchi2",
		"D0_IPErr",
		"D0_IP_xy",
	};
	
	std::vector<std::vector<std::string>> const pT_cuts = {
		{"2.0 < D0_pT", "D0_pT < 5.0"},
	};

	// std::vector<std::vector<std::string>> const cent_cuts = {
	// 	{"2.0 < D0_pT", "D0_pT < 5.0"},
	// };

	TCut get_sideband_cut () {
		std::string mass = "D0_mass";
		boost::format cut("(%f < %s && %s < %f) || (%f < %s && %s < %f)");
		return TCut((cut
				% (mean - 6.0 * sigma) % mass.c_str() % mass.c_str() % (mean - 3.0 * sigma)
				% (mean + 3.0 * sigma) % mass.c_str() % mass.c_str() % (mean + 6.0 * sigma)).str().c_str());
	}
	//...
};

#endif//DEFS_C

