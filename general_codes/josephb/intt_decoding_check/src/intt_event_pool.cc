#include "intt_event_pool.h"

#include <Event/Event.h>
#include <phool/phool.h>

#include <stdlib.h>

#include <iostream>
#include <filesystem>
#include <fstream>
#include <set>

// bool
// bco_comparator::operator() (
// 	bco_t const& lhs,
// 	bco_t const& rhs
// ) {
// 	return ((rhs % BCO_MAX) - (lhs % BCO_MAX) + 2 * BCO_MAX) % BCO_MAX < ((lhs % BCO_MAX) - (rhs % BCO_MAX) + 2 * BCO_MAX) % BCO_MAX;
// }

intt_event_pool::intt_event_pool (
) {
}

intt_event_pool::~intt_event_pool (
) {
	delete m_file_evt_itr;
}

int
intt_event_pool::open_list (
	std::string const& list_file
) {
	if(list_file.empty()) {
		std::cerr << PHWHERE << "\n"
		          << "\tArgument 'list_file' is empty string" << std::endl;
		return EXIT_FAILURE;
	}

	if(!std::filesystem::exists(list_file)) {
		std::cerr << PHWHERE << "\n"
		          << "\tFile '" << list_file << "' not found" << std::endl;
		return EXIT_FAILURE;
	}

	std::ifstream list(list_file);
	for(std::string line; std::getline(list, line);) {
		open(line);
	}

	return EXIT_SUCCESS;
}

int
intt_event_pool::open (
	std::string const& file_name
) {
	if(file_name.empty()) {
		std::cerr << PHWHERE << "\n"
		          << "\tArgument 'file_name' is empty string" << std::endl;
		return EXIT_FAILURE;
	}

	if(!std::filesystem::exists(file_name)) {
		std::cerr << PHWHERE << "\n"
		          << "\tFile '" << file_name << "' not found" << std::endl;
		return EXIT_FAILURE;
	}

	m_files_to_read.push_back(file_name.c_str());
	if(m_verbosity) {
		std::cout << PHWHERE << "\n"
		          << "\tAppended file for reading: " << file_name << std::endl;
	}

	if(!m_file_evt_itr) {
		open_next();
	}

	return EXIT_SUCCESS;
}

int
intt_event_pool::open_next (
) {
	if(!m_files_to_read.size()) {
		if(m_verbosity) {
			std::cout << PHWHERE << "\n"
			          << "\tFinished reading" << std::endl;
		}
		return EXIT_FAILURE;
	}

	if(m_verbosity) {
		std::cout << PHWHERE << "\n"
		          << "\tOpening next file: " << *m_files_to_read.begin() << std::endl;
	}

	m_current_file = *m_files_to_read.begin();
	m_files_to_read.erase(m_files_to_read.begin(), std::next(m_files_to_read.begin()));

	delete m_file_evt_itr;
	m_file_evt_itr = new fileEventiterator(m_current_file.c_str());

	return EXIT_SUCCESS;
}

int
intt_event_pool::set_output_file (
	std::string const& file_name
) {
	delete m_file;
	m_file = TFile::Open(file_name.c_str(), "RECREATE");

	if(!m_file) {
		std::cerr << PHWHERE << "\n"
		          << "\tCould not (re)create file: " << file_name << std::endl;
		return EXIT_FAILURE;
	}

	m_file->cd();
	delete m_tree;
	m_tree = new TTree("tree", "tree");
	m_tree->SetDirectory(m_file);

	m_file->cd();
	delete m_hist;
	m_hist = new TH1I("hist", "hist", 511, -255.5, 255.5);
	m_hist->SetDirectory(m_file);

	return EXIT_SUCCESS;
}

int
intt_event_pool::write_output_file (
) {
	if(!m_file) {
		std::cerr << PHWHERE << "\n"
		          << "\tMember 'm_file' is null at call\n"
		          << "\tWas this preceded by a successful call to 'set_output_file'?" << std::endl;
		return EXIT_FAILURE;
	}

	if(!m_tree) {
		std::cerr << PHWHERE << "\n"
		          << "\tMember 'm_tree' is null at call\n"
		          << "\t(This state should be unreachable)" << std::endl;
		return EXIT_FAILURE;
	} else {
		Long64_t head; m_tree->Branch("head", &head);
		Long64_t min;  m_tree->Branch("min",  &min);
		Long64_t max;  m_tree->Branch("max",  &max);
		Long64_t num;  m_tree->Branch("num",  &num);

		for(auto const& [bco_head, bco_bin_struct] : m_bco_map) {
			head = (Long64_t)bco_head;
			min =  (Long64_t)bco_bin_struct.min;
			max =  (Long64_t)bco_bin_struct.max;
			num =  (Long64_t)bco_bin_struct.count;

			m_tree->Fill();
		}

		m_tree->Write();
	}

	if(!m_hist) {
		std::cerr << PHWHERE << "\n"
		          << "\tMember 'm_hist' is null at call\n"
		          << "\t(This state should be unreachable)" << std::endl;
		return EXIT_FAILURE;
	} else {
		m_hist->Write();
	}

	m_file->Write();
	m_file->Close();

	return EXIT_SUCCESS;
}

int
intt_event_pool::next (
) {
	if(!m_file_evt_itr) {
		std::cerr << PHWHERE << "\n"
		          << "\tMember 'm_file_evt_itr' is null at call\n"
		          << "\tWas this preceded by a successful call to 'open' or 'open_next'?" << std::endl;
		return EXIT_FAILURE;
	}

	Event* evt = m_file_evt_itr->getNextEvent();
	while(!evt) {
		if(open_next()) {
			return EXIT_FAILURE;
		}
		evt = m_file_evt_itr->getNextEvent();
	}

	if(1 < m_verbosity && !(m_evts % m_evt_per_cout)) {
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "file: " << m_current_file << std::endl;
		std::cout << "event: " << std::dec << m_evts << std::endl;
	}

	std::set<bco_t> bco_set;
	for(int pid = 3001; pid < 3009; ++pid) {
		Packet* pkt = evt->getPacket(pid);
		if(!pkt)continue;

		if(1 < m_verbosity && !(m_evts % m_evt_per_cout)) {
			pkt->identify();
			std::cout << "num bcos: " << std::dec << pkt->iValue(0, "NR_BCOS") << std::endl;
			std::cout << "num hits: " << std::dec << pkt->iValue(0, "NR_HITS") << std::endl;
		}

		// for(int n = 0, N = pkt->iValue(0, "NR_HITS"); n < N; ++n) {
		// 	int   fee      =  pkt->iValue(n, "FEE");
		// 	int   chip     = (pkt->iValue(n, "CHIP_ID") + 25) % 26;
		// 	int   channel  =  pkt->iValue(n, "CHANNEL_ID");
		// 	int   fphx_bco =  pkt->iValue(n, "FPHX_BCO");
		// 	int   adc      =  pkt->iValue(n, "ADC");
		// 	bco_t bco      =  pkt->lValue(n, "BCO");
		// 	// ...
		// }

		for(int n = 0, N = pkt->iValue(0, "NR_BCOS"); n < N; ++n) {
		 	bco_t bco = pkt->lValue(n, "BCOLIST");
			bco_set.insert(bco);

			bco_t bco_key = bco & m_BCO_PROJ;

			if(bco < m_bco_map[bco_key].min) {
				m_bco_map[bco_key].min = bco;
			}

			if(m_bco_map[bco_key].max < bco) {
				m_bco_map[bco_key].max = bco;
			}

			++m_bco_map[bco_key].count;
		}

		delete pkt;
	}

	for(auto const& bco : bco_set) {
		if(m_prev_bco == std::numeric_limits<bco_t>::max()) {
			m_prev_bco = bco;
			continue;
		}

		if(m_prev_bco == bco) {
			continue;
		}

		signed long long diff = (signed long long)bco - (signed long long)m_prev_bco;
		if(diff < -255) {
			diff = -255;
		} else if (255 < diff) {
			diff = 255;
		}
		m_hist->Fill(diff);

		m_prev_bco = bco;
	}

	if(1 < m_verbosity && !(m_evts % m_evt_per_cout)) {
		std::cout << "prev bco: 0x" << std::hex << m_prev_bco << std::endl;
	}

	delete evt;
	++m_evts;

	return EXIT_SUCCESS;
}


