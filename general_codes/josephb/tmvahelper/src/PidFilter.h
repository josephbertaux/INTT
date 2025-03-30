#ifndef PID_FILTER_H
#define PID_FILTER_H

#include <map>
#include <string>
#include <vector>

class TTree;

// For filtering combinatorial background from KFP signal data

class PidFilter {
public:
	PidFilter() = default;
	~PidFilter() {clean();}

	int branch(TTree*);
	int eval();

	void set_mother_pdg_id(int);
	void set_num_daughters(int);

private:
	void clean();

	TTree* m_tree{};

	int m_mother_pdg_id{0};
	int m_num_daughters{0};

	int* m_true_id{};
	int* m_pdg_id{};

	std::vector<int>** m_true_track_history_pdg_id{};

	std::vector<float>** m_true_track_history_px{};
	std::vector<float>** m_true_track_history_py{};
	std::vector<float>** m_true_track_history_pz{};
};

#endif//PID_FILTER_H
