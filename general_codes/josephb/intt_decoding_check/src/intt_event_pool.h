#ifndef INTT_EVENT_POOL_H
#define INTT_EVENT_POOL_H

#include "intt_pool.h"

#include <Event/fileEventiterator.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>

#include <limits>
#include <map>
#include <string>
#include <vector>

typedef unsigned long long bco_t;
// struct bco_comparator {
// 	static const bco_t BCO_MAX{bco_t{1}<<40};
// 	bool operator()(bco_t const&, bco_t const&);
// };

class intt_event_pool {
public:
	intt_event_pool();
	~intt_event_pool();

	int open(std::string const&);
	int open_list(std::string const&);

	int verbosity() {return m_verbosity;}
	int verbosity(int const& v) {return m_verbosity = v;}

	int events_per_cout() {return m_evt_per_cout;}
	int events_per_cout(int const& v) {return m_evt_per_cout = 0 < v ? v : 1;}

	int set_output_file(std::string const&);
	int write_output_file();

	int next();

private:
	int open_next();

	static bco_t constexpr m_BCO_PROJ = (bco_t{1}<<40) - (bco_t{1}<<30);

	fileEventiterator* m_file_evt_itr{nullptr};
	intt_pool m_pool;

	std::string m_current_file;
	std::vector<std::string> m_files_to_read;

	int m_verbosity{0};
	unsigned long long m_evts{0};
	unsigned long long m_evt_per_cout{1};

	struct bco_bin_s {
		bco_t min{std::numeric_limits<bco_t>::max()};
		bco_t max{0};
		std::size_t count{0};
	};
	std::map<bco_t, bco_bin_s> m_bco_map;
	bco_t m_prev_bco{std::numeric_limits<bco_t>::max()};

	// bco_comparator m_bco_less;

	TFile* m_file{nullptr};
	TTree* m_tree{nullptr};
	TH1I* m_hist{nullptr};
};

#endif//INTT_EVENT_POOL_H
