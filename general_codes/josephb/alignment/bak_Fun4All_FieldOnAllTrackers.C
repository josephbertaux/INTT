#include <fun4all/Fun4AllUtils.h>

#include <G4_ActsGeom.C>
#include <G4_Global.C>
#include <G4_Magnet.C>
#include <G4_Mbd.C>
#include <GlobalVariables.C>
#include <QA.C>
#include <Trkr_Clustering.C>
#include <Trkr_LaserClustering.C>
#include <Trkr_Reco.C>
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>

#include <ffamodules/CDBInterface.h>
R__LOAD_LIBRARY(libffamodules.so)

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
R__LOAD_LIBRARY(libfun4all.so)

#include <phool/recoConsts.h>
R__LOAD_LIBRARY(libphool.so)

#include <cdbobjects/CDBTTree.h>
R__LOAD_LIBRARY(libcdbobjects.so)

#include <tpccalib/PHTpcResiduals.h>

#include <trackingqa/InttClusterQA.h>
#include <trackingqa/MicromegasClusterQA.h>
#include <trackingqa/MvtxClusterQA.h>
R__LOAD_LIBRARY(libtrackingqa.so)

#include <tpcqa/TpcRawHitQA.h>
R__LOAD_LIBRARY(libtpcqa.so)

#include <trackingdiagnostics/SiliconOnlyTrackResiduals.h>
#include <trackingdiagnostics/TrkrNtuplizer.h>
#include <trackingdiagnostics/KshortReconstruction.h>
R__LOAD_LIBRARY(libTrackingDiagnostics.so)

#include <trackermillepedealignment/AlignmentDefs.h>
#include <trackermillepedealignment/HelicalFitter.h>
R__LOAD_LIBRARY(libtrackeralign.so)

#include <stdio.h>

#include <iostream>
#include <fstream>
#include <filesystem>

R__LOAD_LIBRARY(libmvtx.so)
R__LOAD_LIBRARY(libintt.so)

void Fun4All_FieldOnAllTrackers(
	int line_number = 1,
	std::string const& list_file = "segments.list",
	std::string outfilename = "dat/clusters_seeds_",
	int const nEvents = 100000,
	// int const nEvents = 10,
	bool convertSeeds = true
) {
	outfilename += std::to_string(line_number);

	std::string inputtpcRawHitFile;
	std::ifstream segments_file(list_file, std::ios_base::in);
	if (!segments_file.good()) {
		std::cerr << list_file << std::endl;
		return;
	}

	for (; std::getline(segments_file, inputtpcRawHitFile); --line_number) {
		if (line_number < 0) break;
	}

	if (inputtpcRawHitFile.empty()) return;
	std::cout << "using inputtpcRawHitFile: " << inputtpcRawHitFile << std::endl;

	TRACKING::pp_mode = true;

	G4TRACKING::convert_seeds_to_svtxtracks = convertSeeds;
	std::cout << "Converting to seeds =	" << G4TRACKING::convert_seeds_to_svtxtracks << std::endl;
	std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(std::filesystem::path{inputtpcRawHitFile}.filename());
	int runnumber = runseg.first;
	int segment = runseg.second;

	TpcReadoutInit( runnumber );
	std::cout << " run: " << runnumber
	          << " samples: " << TRACKING::reco_tpc_maxtime_sample
	          << " pre: " << TRACKING::reco_tpc_time_presample
	          << " vdrift: " << G4TPC::tpc_drift_velocity_reco
	          << std::endl;
	// distortion calibration mode
	G4TRACKING::SC_CALIBMODE = false;

	ACTSGEOM::mvtxMisalignment = 100.;
	ACTSGEOM::inttMisalignment = 100.;
	ACTSGEOM::tpotMisalignment = 100.;
	TString outfile = outfilename + "clusters_tracks_" + runnumber + "-" + segment + ".root";
	std::string theOutfile = outfile.Data();
	auto se = Fun4AllServer::instance();
	se->Verbosity(1);
	auto rc = recoConsts::instance();
	rc->set_IntFlag("RUNNUMBER", runnumber);
	rc->set_IntFlag("RUNSEGMENT", segment);

	Enable::CDB = true;
	rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
	rc->set_uint64Flag("TIMESTAMP", runnumber);
	std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");

	Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
	ingeo->AddFile(geofile);
	se->registerInputManager(ingeo);

	//Flag for running the tpc hit unpacker with zero suppression on
	G4MAGNET::magfield_rescale = 1;
	TrackingInit();

	auto hitsin = new Fun4AllDstInputManager("InputManager");
	hitsin->fileopen(inputtpcRawHitFile);
	// hitsin->AddFile(inputMbd);
	se->registerInputManager(hitsin);

	Mvtx_HitUnpacking();
	Intt_HitUnpacking();

	Mvtx_Clustering();
	Intt_Clustering();

	// Silicon Seeding
	auto silicon_Seeding = new PHActsSiliconSeeding;
	silicon_Seeding->Verbosity(999);
	silicon_Seeding->setinttRPhiSearchWindow(1.0);
	silicon_Seeding->setinttZSearchWindow(7.0); 
	silicon_Seeding->seedAnalysis(false);
	se->registerSubsystem(silicon_Seeding);

	auto merger = new PHSiliconSeedMerger;
	merger->Verbosity(0);
	se->registerSubsystem(merger);


	// Either converts seeds to tracks with a straight line/helix fit
	// or run the full Acts track kalman filter fit
	if (G4TRACKING::convert_seeds_to_svtxtracks) {
		auto converter = new TrackSeedTrackMapConverter;
		// Default set to full SvtxTrackSeeds. Can be set to
		//SiliconTrackSeedContainer or TpcTrackSeedContainer
		//converter->setTrackSeedName("SvtxTrackSeedContainer");
		converter->setTrackSeedName("SiliconTrackSeedContainer");
		converter->setFieldMap(G4MAGNET::magfield_tracking);
		converter->Verbosity(999);
		se->registerSubsystem(converter);
	} else {
		auto deltazcorr = new PHTpcDeltaZCorrection;
		deltazcorr->Verbosity(999);
		se->registerSubsystem(deltazcorr);

		// perform final track fit with ACTS
		auto actsFit = new PHActsTrkFitter;
		actsFit->Verbosity(0);
		actsFit->commissioning(G4TRACKING::use_alignment);
		// in calibration mode, fit only Silicons and Micromegas hits
		actsFit->fitSiliconMMs(G4TRACKING::SC_CALIBMODE);
		actsFit->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);
		actsFit->set_pp_mode(TRACKING::pp_mode);
		actsFit->set_use_clustermover(true);	// default is true for now
		actsFit->useActsEvaluator(false);
		actsFit->useOutlierFinder(false);
		actsFit->setFieldMap(G4MAGNET::magfield_tracking);
		se->registerSubsystem(actsFit);

		auto cleaner = new PHTrackCleaner();
		cleaner->Verbosity(1);
		se->registerSubsystem(cleaner);
	}

	auto finder = new PHSimpleVertexFinder;
	finder->Verbosity(0);
	finder->setDcaCut(0.5);
	finder->setTrackPtCut(-99999.);
	finder->setBeamLineCut(1);
	//finder->setTrackQualityCut(150);
	finder->setTrackQualityCut(1000000000);
	finder->setNmvtxRequired(3);
	//	finder->setOutlierPairCut(1);
	finder->setOutlierPairCut(0.1);
	se->registerSubsystem(finder);

	TString residoutfile = theOutfile + "_resid.root";
	std::string residstring(residoutfile.Data());

	auto resid = new SiliconOnlyTrackResiduals();
	resid->outfileName(residstring);
	resid->Verbosity(0);
	resid->alignment(true);

	// adjust track map name
	if(G4TRACKING::SC_CALIBMODE && !G4TRACKING::convert_seeds_to_svtxtracks)
	{
		resid->trackmapName("SvtxSiliconMMTrackMap");
	}

	resid->clusterTree();
	resid->hitTree();
	resid->convertSeeds(G4TRACKING::convert_seeds_to_svtxtracks);
	//resid->set_rejectLaserEvent(true);
	// resid->linefitAll();	// default isTPC only if not set
	se->registerSubsystem(resid);

	std::string hfbinstring = outfilename+"helical_out_"+std::to_string(runnumber)+"-"+std::to_string(segment)+".bin";
	std::string hfsteerstring = outfilename+"helical_steer_"+std::to_string(runnumber)+"-"+std::to_string(segment)+".txt";
	std::string hfntpstring = outfilename+"helical_ntuple_"+std::to_string(runnumber)+"-"+std::to_string(segment)+".root";
	std::cout << hfbinstring << "	" << hfsteerstring << "	" << hfntpstring << std::endl;

	auto hf = new HelicalFitter();
	hf->Verbosity(999);
	hf->set_silicon_track_map_name("SiliconTrackSeedContainer");
	hf->set_datafile_name(hfbinstring);
	hf->set_steeringfile_name(hfsteerstring);
	hf->set_mvtx_grouping(AlignmentDefs::mvtxGrp::snsr);
	hf->set_intt_grouping(AlignmentDefs::inttGrp::lad);
	hf->set_tpc_grouping(AlignmentDefs::tpcGrp::tp);
	hf->set_layer_param_fixed(0, 0);
	hf->set_layer_param_fixed(0, 1);
	hf->set_layer_param_fixed(0, 2);
	hf->set_layer_param_fixed(1, 0);
	hf->set_layer_param_fixed(1, 1);
	hf->set_layer_param_fixed(1, 2);
	hf->set_layer_param_fixed(2, 0);
	hf->set_layer_param_fixed(2, 1);
	hf->set_layer_param_fixed(2, 2);	
	hf->set_layer_param_fixed(3, 0);
	hf->set_layer_param_fixed(3, 1);
	hf->set_layer_param_fixed(3, 2);	
	hf->set_layer_param_fixed(4, 0);
	hf->set_layer_param_fixed(4, 1);
	hf->set_layer_param_fixed(4, 2);
	hf->set_layer_param_fixed(5, 0);
	hf->set_layer_param_fixed(5, 1);
	hf->set_layer_param_fixed(5, 2);
	hf->set_layer_param_fixed(6, 0);
	hf->set_layer_param_fixed(6, 1);
	hf->set_layer_param_fixed(6, 2);
	hf->set_intt_layer_fixed(3);
	hf->set_intt_layer_fixed(4);
	hf->set_use_event_vertex(true);
	hf->set_vertex_param_fixed(0);
	hf->set_vertex_param_fixed(1);
	hf->set_ntuplefile_name(hfntpstring);
	hf->set_fitted_subsystems(true, false, false);	// silicon, tpc, all
	se->registerSubsystem(hf);

	se->run(nEvents);
	se->End();
	se->PrintTimer();

	delete se;
	std::cout << "Finished" << std::endl;
	gSystem->Exit(0);
}
