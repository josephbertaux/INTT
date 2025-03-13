#ifndef CHECK_DST_C
#define CHECK_DST_C

#include "config.C"

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllServer.h>
R__LOAD_LIBRARY(libfun4all.so)

#include <tmvahelper/KFAnalyzer.h>
R__LOAD_LIBRARY(libtmvahelper.so)

#include <phool/recoConsts.h>
R__LOAD_LIBRARY(libphool.so)

#include <filesystem>

void
check_dst (
	std::string const& data_dir
) {

	recoConsts *rc = recoConsts::instance();
	// rc->set_StringFlag("CDB_GLOBALTAG",CDB::global_tag);
	// rc->set_uint64Flag("TIMESTAMP",CDB::timestamp);
	rc->set_IntFlag("RUNNUMBER",1);

	Fun4AllServer *se = Fun4AllServer::instance();
	se->Verbosity(1);

	for (auto const& entry : std::filesystem::directory_iterator{data_dir}) {
		if (!entry.is_regular_file()) continue;

		std::string filename = entry.path().filename();
		if (filename.find(config::channel) == std::string::npos) continue;
		if (filename.find("sig_DST") == std::string::npos) continue;

		std::cout << filename << std::endl;

		// Input
		Fun4AllDstInputManager* in = new Fun4AllDstInputManager("DSTIN", filename);
		se->registerInputManager(in);
		break;
	}

	// KFAnalyzer* kf_ana = new KFAnalyzer;
	// se->registerSubsystem(kf_ana);

	se->run(1);
	se->Print("NODETREE");

	se->End();
	gSystem->Exit(0);
}

#endif//CHECK_DST_C
