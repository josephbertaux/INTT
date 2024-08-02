#ifndef CHECK_C
#define CHECK_C

#include <intt_decoding_check/bco_ana.h>
#include <intt_decoding_check/bco_ana_v2.h>
#include <intt_decoding_check/intt_event_pool.h>
R__LOAD_LIBRARY(libintt_decoding_check.so)

#include <phool/recoConsts.h>
R__LOAD_LIBRARY(libphool.so)

#include <filesystem>
#include <iostream>
#include <string>

void
macro (
	int which_intt,
	int run_num,
	int num_evt,
	std::string const& intt_format,
	std::string const& out_format
) {
	int verbosity = 2;
	int evt_per_cout = 100000;

	char buff[256];
	bco_ana_v2 ana;

	snprintf(buff, sizeof(buff), intt_format.c_str(), run_num, which_intt);
	if(!std::filesystem::exists(buff)) {
		std::cerr << "List file: '" << buff << "' does not exist" << std::endl;

		gSystem->Exit(1);
		exit(1);
	}

	ana.verbosity(verbosity);
	ana.events_per_cout(evt_per_cout);
	ana.open_list(buff);

	snprintf(buff, sizeof(buff), out_format.c_str(), run_num, which_intt);
	ana.set_output_file(buff);

	if(num_evt) {
		for(int n = 0; n < num_evt; ++n) {
			if(ana.next()) break;
		}
	} else {
		for(;;) {
			if(ana.next()) break;
		}
	}

	ana.write_output_file();

	std::cout << "\n\n" << std::endl;
	std::cout << "Macro done" << std::endl;

	gSystem->Exit(0);
	exit(0);
}

#endif//MACRO_C
