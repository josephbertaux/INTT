#include "PidFilter.h"

#include <TTree.h>

#include <cstdlib>
#include <iostream>

#include <boost/format.hpp>

int
PidFilter::branch (
	TTree* tree
) {
	if (!tree) {
		std::cerr
			<< __PRETTY_FUNCTION__ << " @ " << __FILE__ << ":" << __LINE__ << "\n"
			<< "\tArgument 'tree' is null\n"
			<< std::flush;
		return EXIT_FAILURE;
	}

	std::vector<std::string> missing_branches;
	for (int n = 0; n < m_num_daughters; ++n) {
		std::string name;

		name = (boost::format("track_%d_true_ID") % (n + 1)).str();
		if (tree->SetBranchAddress(name.c_str(), &(m_true_id[n])) != TTree::kMatch) missing_branches.push_back(name);

		name = (boost::format("track_%d_PDG_ID") % (n + 1)).str();
		if (tree->SetBranchAddress(name.c_str(), &(m_pdg_id[n])) != TTree::kMatch) missing_branches.push_back(name);

		name = (boost::format("track_%d_true_track_history_PDG_ID") % (n + 1)).str();
		if (tree->SetBranchAddress(name.c_str(), &(m_true_track_history_pdg_id[n])) != TTree::kMatch) missing_branches.push_back(name);
		name = (boost::format("track_%d_true_track_history_px") % (n + 1)).str();
		if (tree->SetBranchAddress(name.c_str(), &(m_true_track_history_px[n])) != TTree::kMatch) missing_branches.push_back(name);
		name = (boost::format("track_%d_true_track_history_py") % (n + 1)).str();
		if (tree->SetBranchAddress(name.c_str(), &(m_true_track_history_py[n])) != TTree::kMatch) missing_branches.push_back(name);
		name = (boost::format("track_%d_true_track_history_pz") % (n + 1)).str();
		if (tree->SetBranchAddress(name.c_str(), &(m_true_track_history_pz[n])) != TTree::kMatch) missing_branches.push_back(name);
	}

	if (missing_branches.size()) {
		std::cerr
			<< __PRETTY_FUNCTION__ << " @ " << __FILE__ << ":" << __LINE__ << "\n"
			<< "\tFailed to set branch addresses\n";
		for (auto const& name : missing_branches) std::cerr
			<< "\t" << name << std::endl;
		std::cerr
			<< std::flush;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

int
PidFilter::eval (
) {
	// Require every branch is as advertized
	for (int n = 0; n < m_num_daughters; ++n) {
		if (m_true_id[n] != m_pdg_id[n]) return EXIT_FAILURE;
	}

	// Require that all daughters came from the same mother
	float px{0}, py{0}, pz{0};
	for (int n = 0; n < m_num_daughters; ++n) {
		// find the mother in the track history of this daughter
		std::size_t mother_index = 0;
		for (auto const& pdg_id : *m_true_track_history_pdg_id[n]) {
			if (abs(pdg_id) == m_mother_pdg_id) break;
			++mother_index;
		}

		// No mother found
		if (mother_index == m_true_track_history_pdg_id[n]->size()) return EXIT_FAILURE;

		// mother momentum
		float px_n = m_true_track_history_px[n]->at(mother_index);
		float py_n = m_true_track_history_py[n]->at(mother_index);
		float pz_n = m_true_track_history_pz[n]->at(mother_index);

		// Successive passes compare
		if (n) {
			if(px != px_n) return EXIT_FAILURE;
			if(py != py_n) return EXIT_FAILURE;
			if(pz != pz_n) return EXIT_FAILURE;
		}

		px = px_n;
		py = py_n;
		pz = pz_n;
	}

	return EXIT_SUCCESS;
}

void
PidFilter::set_mother_pdg_id (
	int mother_pdg_id
) {
	m_mother_pdg_id = abs(mother_pdg_id);
}

void
PidFilter::set_num_daughters (
	int num_daughters
) {
	clean();
	m_num_daughters = num_daughters;

	m_true_id = new int[num_daughters];
	m_pdg_id = new int[num_daughters];

	m_true_track_history_pdg_id = new std::vector<int>*[num_daughters];
	m_true_track_history_px = new std::vector<float>*[num_daughters];
	m_true_track_history_py = new std::vector<float>*[num_daughters];
	m_true_track_history_pz = new std::vector<float>*[num_daughters];

	for (int n = 0; n < m_num_daughters; ++n) {
		m_true_track_history_pdg_id[n] = new std::vector<int>;
		m_true_track_history_px[n] = new std::vector<float>;
		m_true_track_history_py[n] = new std::vector<float>;
		m_true_track_history_pz[n] = new std::vector<float>;
	}
}

void
PidFilter::clean (
) {
	for (int n = 0; n < m_num_daughters; ++n) {
		delete m_true_track_history_pdg_id[n];
		delete m_true_track_history_px[n];
		delete m_true_track_history_py[n];
		delete m_true_track_history_pz[n];
	}

	delete[] m_true_track_history_pdg_id;
	delete[] m_true_track_history_px;
	delete[] m_true_track_history_py;
	delete[] m_true_track_history_pz;

	delete[] m_true_id;
	delete[] m_pdg_id;

	m_num_daughters = 0;
}
