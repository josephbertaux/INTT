#include "TMVAHelper.h"

#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <TFile.h>

#include <iostream>
#include <fstream>

#include <boost/format.hpp>

TTree*
TMVAHelper::get_tree (
	std::string const& file_name,
	std::string const& tree_name
) {
	if (file_name.empty()) {
		std::cerr
			<< __PRETTY_FUNCTION__ << " @ " << __FILE__ << ":" << __LINE__ << "\n"
			<< std::flush;
		return nullptr;
	}
	TFile* file = TFile::Open(file_name.c_str(), "READ");
	if (!file) {
		std::cerr
			<< __PRETTY_FUNCTION__ << " @ " << __FILE__ << ":" << __LINE__ << "\n"
			<< "\tfile: " << file_name << "\n"
			<< std::flush;
		return nullptr;
	}

	if (tree_name.empty()) {
		std::cerr
			<< __PRETTY_FUNCTION__ << " @ " << __FILE__ << ":" << __LINE__ << "\n"
			<< std::flush;
		return nullptr;
	}
	TTree* tree = dynamic_cast<TTree*>(file->Get(tree_name.c_str()));
	if (!tree) {
		std::cerr
			<< __PRETTY_FUNCTION__ << " @ " << __FILE__ << ":" << __LINE__ << "\n"
			<< "\tfile: " << file_name << "\n"
			<< "\ttree: " << tree_name << "\n"
			<< std::flush;
		return nullptr;
	}

	return tree;
}

int
TMVAHelper::read_file (
	std::string const& file_name,
	std::vector<std::string>& names
) {
	std::ifstream file(file_name, std::ios_base::in);
	if (!file.good()) {
		std::cerr
			<< __PRETTY_FUNCTION__ << " @ " << __FILE__ << ":" << __LINE__ << "\n"
			<< "file: " << file_name << "\n"
			<< std::flush;
		return EXIT_FAILURE;
	}

	for (std::string line; std::getline(file, line);) {
		for (std::size_t pos; (pos = line.find("#")) != std::string::npos;) {
			line = line.substr(0, pos);
		}

		if (line.empty()) continue;

		names.push_back(line);
	}

	return EXIT_SUCCESS;
}

void
TMVAHelper::read_branches (
	std::vector<std::string> const& branches_names
) {
	for (auto const& name : branches_names) {
		m_branches_names.push_back(name);
	}

	init_branches();
}

void
TMVAHelper::read_training (
	std::vector<std::string> const& training_names
) {
	for (auto const& name : training_names) {
		m_training_names.push_back(name);
	}

	init_training();
}

void
TMVAHelper::read_cuts (
	std::vector<std::string> const& cut_names
) {
	for (auto const& name : cut_names) {
		m_cuts_names.push_back(name);
	}

	init_cuts();
}

int
TMVAHelper::read_branches (
	std::string const& branches_file_name
) {
	if (read_file(branches_file_name, m_branches_names)) return EXIT_FAILURE;

	init_branches();
	return EXIT_SUCCESS;
}

int
TMVAHelper::read_training (
	std::string const& training_file_name
) {
	if (read_file(training_file_name, m_training_names)) return EXIT_FAILURE;

	init_training();
	return EXIT_SUCCESS;
}

int
TMVAHelper::read_cuts (
	std::string const& cuts_file_name
) {
	if (read_file(cuts_file_name, m_cuts_names)) return EXIT_FAILURE;

	init_cuts();
	return EXIT_SUCCESS;
}

void
TMVAHelper::init_branches (
) {
	m_branches_map_i.clear();
	m_branches_map_f.clear();

	m_branches_args.Clear();
	for (auto const& name : m_branches_names) {

		std::size_t pos = name.find("/");
		std::string n = pos == std::string::npos ? name : name.substr(0,  pos);
		std::string t = pos == std::string::npos ? name : name.substr(pos + 1);

		if (pos == std::string::npos || t == "F") { // Float_t by default
			m_branches_map_f[n] = 0.0;
		} else if (t == "I") { // Int
			m_branches_map_i[n] = 0;
		} else {
			std::cerr
				<< __PRETTY_FUNCTION__ << " @ " << __FILE__ << ":" << __LINE__ << "\n"
				<< "\tBranch " << n << " specifies untreated type '" << t << "'\n"
				<< "\t\t" << name << "\n"
				<< std::flush;
			continue;
		}

		m_branches_args.addOwned ( *new RooRealVar (
			n.c_str(), n.c_str(), 0.0,
			-std::numeric_limits<Float_t>::max(), std::numeric_limits<Float_t>::max()
		) );
	}
}

void
TMVAHelper::init_training (
) {
	m_training_map.clear();
	m_training_args.Clear();
	for (auto const& name : m_training_names) {

		std::size_t pos = name.find(":=");
		std::string f = pos == std::string::npos ? name : name.substr(pos + 2);

		m_training_map[name] = 0.0;
		m_training_args.addOwned ( *new RooFormulaVar (
			name.c_str(), f.c_str(),
			m_branches_args, kFALSE
		) );
	}
}

void
TMVAHelper::init_cuts (
) {
	m_cuts_map.clear();
	m_cuts_args.Clear();
	for (auto const& name : m_cuts_names) {
		m_cuts_map[name] = 0.0;
		m_cuts_args.addOwned ( *new RooFormulaVar (
			name.c_str(), name.c_str(),
			m_branches_args, kFALSE
		) );
	}
}

int
TMVAHelper::branch (
	TTree* tree
) {
	int rv = EXIT_SUCCESS;
	for (auto& [name, val] : m_branches_map_i) {
		if (!tree->GetBranch(name.c_str())) {
			std::cerr
				<< __FILE__ << ":" << __LINE__ << "\n"
				<< "\tbranch name: " << name << "\n"
				<< std::flush;
			rv = EXIT_FAILURE;
			continue;
		}
		tree->SetBranchAddress(name.c_str(), &val);
	}

	for (auto& [name, val] : m_branches_map_f) {
		if (!tree->GetBranch(name.c_str())) {
			std::cerr
				<< __FILE__ << ":" << __LINE__ << "\n"
				<< "\tbranch name: " << name << "\n"
				<< std::flush;
			rv = EXIT_FAILURE;
			continue;
		}
		tree->SetBranchAddress(name.c_str(), &val);
	}

	return rv;
}

void
TMVAHelper::branch (
	TMVA::DataLoader* dataloader
) const {
	for (auto const& name : m_training_names) {
		dataloader->AddVariable(name.c_str());
	}

	boost::format no_nan("%s == %s");
	for (auto const& name : m_branches_names) {
		std::size_t pos = name.find("/");
		std::string n = pos == std::string::npos ? name : name.substr(0,  pos);

		TCut cut = (no_nan % n % n).str().c_str();
		dataloader->AddCut(cut, "Signal");
		dataloader->AddCut(cut, "Background");
	}

	for (auto const& name : m_cuts_names) {
		TCut cut = name.c_str();
		dataloader->AddCut(cut, "Signal");
		dataloader->AddCut(cut, "Background");
	}
}

void
TMVAHelper::branch (
	TMVA::Reader* reader
) {
	for (auto& name : m_training_names) {
		if (m_training_map.find(name) == m_training_map.end()) continue;
		reader->AddVariable(name.c_str(), &(m_training_map[name]));
	}
}

void*
TMVAHelper::get_branch (
	std::string const& name
) {
	if (m_branches_map_i.find(name) != m_branches_map_i.end()) return static_cast<void*>(&m_branches_map_i[name]);
	if (m_branches_map_f.find(name) != m_branches_map_f.end()) return static_cast<void*>(&m_branches_map_f[name]);
	return nullptr;
}

int
TMVAHelper::eval (
) {
	for (auto const& [name, val] : m_branches_map_i) {
		// if (!(val == val)) return EXIT_FAILURE; // int never represents NaN
		dynamic_cast<RooRealVar&>(m_branches_args[name]).setVal(val);
	}
	for (auto const& [name, val] : m_branches_map_f) {
		if (!(val == val)) return EXIT_FAILURE; // IEEE NaN filtering
		dynamic_cast<RooRealVar&>(m_branches_args[name]).setVal(val);
	}

	for (auto& [name, val] : m_training_map) {
		val = dynamic_cast<RooFormulaVar&>(m_training_args[name]).getValV();
		if (!(val == val)) return EXIT_FAILURE; // IEEE NaN filtering
	}

	for (auto& [name, val] : m_cuts_map) {
		val = dynamic_cast<RooFormulaVar&>(m_cuts_args[name]).getValV();
		if (!(val == val)) return EXIT_FAILURE; // IEEE NaN filtering
		if (val == 0) return EXIT_FAILURE; // Doesn't pass cut criteria
	}

	return EXIT_SUCCESS;
}

void
TMVAHelper::show (
) const {
	std::cout << __PRETTY_FUNCTION__ << " @ " << __FILE__ << ":" << __LINE__ << std::endl;

	for (auto const& [name, val] : m_branches_map_f) {
		std::cout << "\t" << name << ": " << val << std::endl;
	}
	for (auto const& [name, val] : m_branches_map_i) {
		std::cout << "\t" << name << ": " << val << std::endl;
	}

	for (auto const& [name, val] : m_training_map) {
		std::cout << "\t" << name << ": " << val << std::endl;
	}
	std::cout << std::endl;

	for (auto const& [name, val] : m_cuts_map) {
		std::cout << "\t" << name << ": " << val << std::endl;
	}
}

