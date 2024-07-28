#ifndef INTT_DECODING_CHECK_H
#define INTT_DECODING_CHECK_H

#include <cstdio>
#include <cstdint>
#define __COUT_HERE__ fprintf(stdout, "%s @ %s:%d\n", __PRETTY_FUNCTION__, __FILE__, __LINE__);
#define __CERR_HERE__ fprintf(stderr, "%s @ %s:%d\n", __PRETTY_FUNCTION__, __FILE__, __LINE__);

#include <map>
#include <string>
#include <vector>

#include <fun4all/SubsysReco.h>

#include <RtypesCore.h>

class PHCompositeNode;


class intt_decoding_check : public SubsysReco {
public:
	intt_decoding_check(std::string const& = "intt_decoding_check");
	~intt_decoding_check() override;

	typedef unsigned long long bco_t;
	struct bco_comparator {
		static const bco_t m_BCO_MAX{bco_t{1}<<40};
		bool operator()(bco_t const&, bco_t const&) const;
	};

	int Init(PHCompositeNode*) override;
	int InitRun(PHCompositeNode*) override;

	int process_event(PHCompositeNode*) override;
	int ResetEvent(PHCompositeNode*) override;

	int EndRun(const int runnumber) override;
	int End(PHCompositeNode*) override;

	int Verbosity(int const& verbosity) {return m_verbosity = verbosity;}
	int Verbosity() {return m_verbosity;}

	int EvtPerCout(int const& evt_per_cout) {return m_evt_per_cout = (0 < evt_per_cout) ? evt_per_cout : 1;}
	int EvtPerCout() {return m_evt_per_cout;}

	int Reset(PHCompositeNode*) override;

	void Print(std::string const& = "ALL") const override;

private:
	// inline std::string static const m_prdf_node_name = "PRDF";
	inline std::string static const m_intt_node_name = "INTTRAWHIT";

	int m_verbosity{0};
	int m_evt_per_cout = 100000;
	unsigned long long m_evts{0}; // RC DAQ event counter

	static constexpr bco_t m_BCO_PROJECTION = (bco_t{1}<<40) - (bco_t{1}<<24);
	bco_t m_begin_bco[8] = {};
	bco_t m_end_bco[8] = {};
	unsigned long long m_bco_evt[8] = {0};

	std::map<bco_t, std::size_t, bco_comparator> m_bco_map[8];

	bco_comparator m_bco_less;
};

#endif//INTT_DECODING_CHECK_H
