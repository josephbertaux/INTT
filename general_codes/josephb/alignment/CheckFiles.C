#ifndef CHECK_FILES
#define CHECK_FILES

#include <filesystem>

void
CheckFiles (
) {
	for (auto const& dir_entry : std::filesystem::directory_iterator{"dat"}) {
		if (!dir_entry.is_regular_file()) continue;
		std::string file_name = dir_entry.path().filename();

		// Only ROOT files
		if (file_name.find(".root") == std::string::npos) continue;

		// Get only helical fitter related output for now
		if (file_name.find("helical") == std::string::npos) continue;

		TFile* file = TFile::Open(std::string{dir_entry.path()}.c_str(), "READ");
		if (!file) continue;

		TNtuple* ntp       = dynamic_cast<TNtuple*>(file->Get("ntp"      ));
		TNtuple* track_ntp = dynamic_cast<TNtuple*>(file->Get("track_ntp"));

		if (!ntp || !track_ntp) {
			std::cout << file_name << " problem" << std::endl;
			continue;
		}

		if (ntp->GetEntriesFast() == 0 && track_ntp->GetEntriesFast() == 0) continue;
		std::cout << file_name << " non-empty" << std::endl;

		file->Close();
	}
}

#endif//CHECK_FILES
