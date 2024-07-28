#include "intt_decoding_check.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainer.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>

#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <array>
#include <iostream>
#include <limits>

intt_decoding_check::intt_decoding_check (
	const std::string& name
) : SubsysReco(name) {
	// do nothing
}

intt_decoding_check::~intt_decoding_check (
) {
}

int
intt_decoding_check::Init (
	PHCompositeNode*
) {
	return Fun4AllReturnCodes::EVENT_OK;
}

int
intt_decoding_check::InitRun (
	PHCompositeNode* top_node
) {
	m_evts = 0;
	for(int i = 0; i < 8; ++i) {
	   	m_begin_bco[i] = std::numeric_limits<bco_t>::max();
	   	m_end_bco[i] = std::numeric_limits<bco_t>::max();
	   	m_bco_evt[i] = 0;
	}

	for(auto& bco_map : m_bco_map) {
		bco_map.clear();
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int
intt_decoding_check::process_event (
	PHCompositeNode* top_node
) {
	/*
	Event* evt = findNode::getClass<Event>(top_node, m_prdf_node_name);
	if(!evt) {
		__CERR_HERE__
		fprintf(stderr, "\tCould not get node:\n");
		fprintf(stderr, "\t\t%s\n", m_prdf_node_name.c_str());
		fprintf(stderr, "\tExiting\n");

		gSystem->Exit(1);
		exit(1);
	}

	for(int pid = 3001; pid < 3009; ++pid) {
		Packet* pkt = evt->getPacket(pid);
		if(!pkt) continue;

		std::cout << "has packet: " << pid << std::endl;
	}
	*/

	InttRawHitContainer* raw_hit_container = findNode::getClass<InttRawHitContainer>(top_node, m_intt_node_name);
	if(!raw_hit_container) {
		__CERR_HERE__
		fprintf(stderr, "\tCould not get node:\n");
		fprintf(stderr, "\t\t%s\n", m_intt_node_name.c_str());
		fprintf(stderr, "\tExiting\n");

		gSystem->Exit(1);
		exit(1);
	}

	for(int n = 0, N = raw_hit_container->get_nhits(); n < N; ++n) {
		InttRawHit* hit{nullptr};
		if(!(hit = raw_hit_container->get_hit(n))) continue;

		int which_intt = hit->get_packetid() - 3001;
		bco_t bco = hit->get_bco();

		bco_t bco_projection = bco & m_BCO_PROJECTION;
		++m_bco_map[which_intt][bco_projection];

		if(m_bco_less(bco, m_begin_bco[which_intt]) || m_begin_bco[which_intt] == std::numeric_limits<bco_t>::max()) {
			m_begin_bco[which_intt] = bco;
		}

		if(m_bco_less(m_end_bco[which_intt], bco) || m_begin_bco[which_intt] == std::numeric_limits<bco_t>::max()) {
			m_end_bco[which_intt] = bco;
			m_bco_evt[which_intt] = m_evts;
		}
	}

	if(m_verbosity && !(m_evts % m_evt_per_cout)) {
		std::cout << m_evts << std::endl;
	}
	++m_evts;

	return Fun4AllReturnCodes::EVENT_OK;
}

int
intt_decoding_check::ResetEvent (
	PHCompositeNode*
) {
	return Fun4AllReturnCodes::EVENT_OK;
}

int
intt_decoding_check::EndRun (
	const int runnumber
) {
	for(int i = 0; i < 8; ++i) {
		std::cout << "intt" << i << " recieved hits through RCDAQ event " << m_bco_evt[i] << std::endl;
	}

	std::cout << "bco list:" << std::endl;
	for(int i = 0; i < 8; ++i) {
		std::cout << "intt" << i << std::endl;
		for(auto const& [bco, count] : m_bco_map[i]) {
			std::cout << std::dec << bco << " " << count << std::endl;
		}
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int
intt_decoding_check::End (
	PHCompositeNode*
) {
	return Fun4AllReturnCodes::EVENT_OK;
}

int
intt_decoding_check::Reset (
	PHCompositeNode*
) {
	return Fun4AllReturnCodes::EVENT_OK;
}

void
intt_decoding_check::Print (
	std::string const&
) const {
	// do nothing
}

bool
intt_decoding_check::bco_comparator::operator() (
	intt_decoding_check::bco_t const& lhs,
	intt_decoding_check::bco_t const& rhs
) const {
	return ((rhs % m_BCO_MAX) - (lhs % m_BCO_MAX) + 2 * m_BCO_MAX) % m_BCO_MAX < ((lhs % m_BCO_MAX) - (rhs % m_BCO_MAX) + 2 * m_BCO_MAX) % m_BCO_MAX;
}
