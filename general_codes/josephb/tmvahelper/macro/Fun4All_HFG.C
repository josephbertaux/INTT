#ifndef MACRO_FUN4ALLG4SPHENIX_C
#define MACRO_FUN4ALLG4SPHENIX_C

#include <G4_Input.C>
#include <G4_Global.C>
#include <G4Setup_sPHENIX.C>

#include <Trkr_RecoInit.C>
#include <Trkr_Clustering.C>
#include <Trkr_TruthTables.C>
#include <Trkr_Reco.C>

#include <phpythia8/PHPy8ParticleTrigger.h>

#include <decayfinder/DecayFinder.h>
R__LOAD_LIBRARY(libdecayfinder.so)

#include <hftrackefficiency/HFTrackEfficiency.h>
R__LOAD_LIBRARY(libhftrackefficiency.so)

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wundefined-internal"
#include <kfparticle_sphenix/KFParticle_sPHENIX.h>
#pragma GCC diagnostic pop
R__LOAD_LIBRARY(libkfparticle_sphenix.so)

#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>
R__LOAD_LIBRARY(libffamodules.so)
#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
R__LOAD_LIBRARY(libfun4all.so)

struct decay_config_s {
	std::string channel;
	std::string pythia_config_file;
	std::string evtgen_config_file;
	std::string decay_descriptor;
	int particle_trigger;
};

decay_config_s bs2jpsiks0 {
	.channel            = "bs2jpsiks0",
	.pythia_config_file = "steeringCards/pythia8_bs2jpsiks0.cfg",
	.evtgen_config_file = "decFiles/Bs2JpsiKS0.DEC",
	.decay_descriptor   = "[B_s0 -> {J/psi -> e^+ e^-} {K_S0 -> pi^+ pi^-}]cc",
	.particle_trigger   = 531
};

decay_config_s D0_Kpi {
	.channel            = "D0_Kpi",
	.pythia_config_file = "steeringCards/pythia8_D2Kpi.cfg",
	.evtgen_config_file = "decFiles/D2Kpi.DEC",
	.decay_descriptor   = "[D0 -> K^- pi^+]cc",
	.particle_trigger   = 421
};

decay_config_s Lc_pKpi {
	.channel            = "Lc_pKpi",
	.pythia_config_file = "steeringCards/pythia8_Lc_pKpi.cfg",
	.evtgen_config_file = "decFiles/Lc_pKpi.DEC",
	.decay_descriptor   = "[Lambda_c+ -> proton^+ K^- pi^+]cc",
	.particle_trigger   = 4122
};

int Fun4All_HFG (
	std::string processID = "0",
	int nEvents = 2e3
) {
	decay_config_s decay_config = D0_Kpi;
	// decay_config_s decay_config = Lc_pKpi;

	//F4A setup
	
	Fun4AllServer *se = Fun4AllServer::instance();
	se->Verbosity(1);

	PHRandomSeed::Verbosity(1);
	recoConsts *rc = recoConsts::instance();

	//Generator setup

	Input::PYTHIA8 = true;
	PYTHIA8::config_file     = decay_config.pythia_config_file;
	EVTGENDECAYER::DecayFile = decay_config.evtgen_config_file;

	Input::BEAM_CONFIGURATION = Input::pp_COLLISION;

	InputInit();

	PHPy8ParticleTrigger * p8_hf_signal_trigger = new PHPy8ParticleTrigger();
	p8_hf_signal_trigger->AddParticles( decay_config.particle_trigger);
	p8_hf_signal_trigger->AddParticles(-decay_config.particle_trigger);

	p8_hf_signal_trigger->SetPtLow(1.);
	p8_hf_signal_trigger->SetEtaHighLow(1.3, -1.3); // sample a rapidity range higher than the sPHENIX tracking pseudorapidity
	p8_hf_signal_trigger->SetStableParticleOnly(false); // process unstable particles that include quarks
	p8_hf_signal_trigger->PrintConfig();
	INPUTGENERATOR::Pythia8->register_trigger(p8_hf_signal_trigger);
	INPUTGENERATOR::Pythia8->set_trigger_OR();

	Input::ApplysPHENIXBeamParameter(INPUTGENERATOR::Pythia8);

	InputRegister();

	//CDB flags and such

	rc->set_IntFlag("RUNNUMBER",1);

	SyncReco *sync = new SyncReco();
	se->registerSubsystem(sync);

	HeadReco *head = new HeadReco();
	se->registerSubsystem(head);

	FlagHandler *flag = new FlagHandler();
	se->registerSubsystem(flag);


	DecayFinder *myFinder = new DecayFinder("myFinder");
	myFinder->Verbosity(INT_MAX);
	myFinder->setDecayDescriptor(decay_config.decay_descriptor);
	myFinder->saveDST(1);
	myFinder->allowPi0(1);
	myFinder->allowPhotons(1);
	myFinder->triggerOnDecay(1);
	myFinder->setPTmin(0.16); //Note: sPHENIX min pT is 0.2 GeV for tracking
	myFinder->setEtaRange(-1.2, 1.2); //Note: sPHENIX acceptance is |eta| <= 1.1
	myFinder->useDecaySpecificEtaRange(false);
	se->registerSubsystem(myFinder);	

	//Simulation setup

	Enable::MBDFAKE = true;
	Enable::PIPE = true;
	Enable::PIPE_ABSORBER = true;
	Enable::MVTX = true;
	Enable::INTT = true;
	Enable::TPC = true;
	Enable::MICROMEGAS = true;

	//Tracking setup

	Enable::CDB = true;
	rc->set_StringFlag("CDB_GLOBALTAG",CDB::global_tag);
	rc->set_uint64Flag("TIMESTAMP",CDB::timestamp);

	MagnetInit();
	MagnetFieldInit();

	G4Init();
	G4Setup();

	Mbd_Reco();
	Mvtx_Cells();
	Intt_Cells();
	TPC_Cells();
	Micromegas_Cells();

	TrackingInit();

	Mvtx_Clustering();
	Intt_Clustering();
	TPC_Clustering();
	Micromegas_Clustering();

	Tracking_Reco();
	Global_Reco();

	build_truthreco_tables();

	//HF stuff
	HFTrackEfficiency *myTrackEff = new HFTrackEfficiency("myTrackEff");
	myTrackEff->Verbosity(INT_MAX);
	myTrackEff->setDFNodeName("myFinder");
	myTrackEff->triggerOnDecay(1);
	myTrackEff->writeSelectedTrackMap(true);
	myTrackEff->writeOutputFile(true);
	std::string outputHFEffFile = "./dat/outputHFTrackEff_" + decay_config.channel + "_" + processID + ".root";
	myTrackEff->setOutputFileName(outputHFEffFile);
	se->registerSubsystem(myTrackEff);

	//KFParticle stuff
	KFParticle_sPHENIX* myKFParticle = new KFParticle_sPHENIX("myKFParticle");
	myKFParticle->setDecayDescriptor(decay_config.decay_descriptor);
	myKFParticle->setTrackMapNodeName("HFSelected_SvtxTrackMap");

	myKFParticle->constrainToPrimaryVertex(true);
	myKFParticle->setMotherIPchi2(999.0);
	myKFParticle->setFlightDistancechi2(-1.);
	myKFParticle->setMinDIRA(-1.1);
	myKFParticle->setDecayLengthRange(-1*FLT_MAX, FLT_MAX);

	//Track parameters
	myKFParticle->setMinimumTrackPT(0.0);
	myKFParticle->setMinimumTrackIPchi2(-1.0);
	myKFParticle->setMinimumTrackIP(-1.0);
	myKFParticle->setMaximumTrackchi2nDOF(999.0);

	//Vertex parameters
	myKFParticle->setMaximumVertexchi2nDOF(999.0);
	myKFParticle->setMaximumDaughterDCA(999.0);

	//Parent parameters
	myKFParticle->setMotherPT(0);
	myKFParticle->setMinimumMass(1.50);
	myKFParticle->setMaximumMass(3.50);
	myKFParticle->setMaximumMotherVertexVolume(999.0);
	myKFParticle->saveDST();
	std::string outputKFParticleFile = "./dat/outputKFParticle_" + decay_config.channel + "_" + processID + ".root";
	myKFParticle->setOutputName(outputKFParticleFile);
	se->registerSubsystem(myKFParticle);

	//Output file handling

	string FullOutFile = "./dat/" + decay_config.channel + "_DST_" + processID + ".root";
	Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", FullOutFile);
	out->StripNode("G4HIT_PIPE");
	out->StripNode("G4HIT_SVTXSUPPORT");
	out->StripNode("PHG4INEVENT");
	out->StripNode("Sync");
	out->StripNode("myFinder_DecayMap");
	out->StripNode("G4HIT_PIPE");
	out->StripNode("G4HIT_MVTX");
	out->StripNode("G4HIT_INTT");
	out->StripNode("G4HIT_TPC");
	out->StripNode("G4HIT_MICROMEGAS");
	out->StripNode("TRKR_HITSET");
	out->StripNode("TRKR_HITTRUTHASSOC");
	out->StripNode("TRKR_CLUSTER");
	out->StripNode("TRKR_CLUSTERHITASSOC");
	out->StripNode("TRKR_CLUSTERCROSSINGASSOC");
	out->StripNode("TRAINING_HITSET");
	out->StripNode("TRKR_TRUTHTRACKCONTAINER");
	out->StripNode("TRKR_TRUTHCLUSTERCONTAINER");
	out->StripNode("alignmentTransformationContainer");
	out->StripNode("alignmentTransformationContainerTransient");
	out->StripNode("SiliconTrackSeedContainer");
	out->StripNode("TpcTrackSeedContainer");
	out->StripNode("SvtxTrackSeedContainer");
	out->StripNode("ActsTrajectories");
	// out->StripNode("SvtxTrackMap");
	out->StripNode("SvtxAlignmentStateMap");
	out->SaveRunNode(0);
	se->registerOutputManager(out);

	se->run(nEvents);
	se->End();
	gSystem->Exit(0);

	return 0;
}

#endif
