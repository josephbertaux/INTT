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

	std::string const particle_name = "D0"; // "Lambda_cplus";
	std::string const mass_branch = particle_name + "_mass";
	std::vector<std::string> const branches = {
		mass_branch,

		particle_name + "_decayLength", particle_name + "_decayLengthErr",

		particle_name + "_px", particle_name + "_py",
		"track_1_px", "track_1_py",
		"track_2_px", "track_2_py",
		// "track_3_px", "track_3_py",

		"track_1_pT", "track_1_pTErr",
		"track_2_pT", "track_2_pTErr",
		// "track_3_pT", "track_3_pTErr",
	};
	
	std::vector<std::string> const training = {
		std::string{"LogDLS := log(abs("} + particle_name + "_decayLength / " + particle_name + "_decayLengthErr))",

		std::string{"acos("} +
			  "((track_1_px + track_2_px) * " + particle_name + "_px" +
			" + (track_1_py + track_2_py) * " + particle_name + "_py)" +
			"/ sqrt((track_1_px + track_2_px) * (track_1_px + track_2_px)" +
			    " + (track_1_py + track_2_py) * (track_1_py + track_2_py))" +
			"/ sqrt(" + particle_name + "_px * " + particle_name + "_px + " + particle_name + "_py * " + particle_name + "_py))",

		// std::string{"acos("} +
		// 	  "((track_1_px + track_2_px + track_3_px) * " + particle_name + "_px" +
		// 	" + (track_1_py + track_2_py + track_3_py) * " + particle_name + "_py)" +
		// 	"/ sqrt((track_1_px + track_2_px + track_3_px) * (track_1_px + track_2_px + track_3_px)" +
		// 	    " + (track_1_py + track_2_py + track_3_py) * (track_1_py + track_2_py + track_3_py))" +
		// 	"/ sqrt(" + particle_name + "_px * " + particle_name + "_px + " + particle_name + "_py * " + particle_name + "_py))",

		"track_1_pT / track_1_pTErr",
		"track_2_pT / track_2_pTErr",
		// "track_3_pT / track_3_pTErr",
	};
	
	TCut get_sideband_cut () {
		boost::format cut("(%f < %s && %s < %f) || (%f < %s && %s < %f)");
		return TCut((cut
				% (mean - 6.0 * sigma) % mass_branch.c_str() % mass_branch.c_str() % (mean - 3.0 * sigma)
				% (mean + 3.0 * sigma) % mass_branch.c_str() % mass_branch.c_str() % (mean + 6.0 * sigma)).str().c_str());
	}
};

#endif//DEFS_C

