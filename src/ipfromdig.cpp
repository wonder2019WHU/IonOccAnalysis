#include<getopt.h>
#include<unistd.h>
#include<iostream>
#include "IonPrf.h"
using namespace std;

int main(int argc, char* argv[]) {
	//About the option
	int opt;
	string argument;
	struct option longopts[] = {
		{"dirent", 0, NULL, 'd'},
		{"regex", 1, NULL, 'r'},
		{"qcscore", 1, NULL, 'q'},
		{"nmf2", 0, NULL, 'n'},
		{"hmf2", 0, NULL, 'h'},
		{"parameter", 0, NULL, 'p'},
		{0, 0, 0, 0}
	};
	bool is_dir = false, is_parN=false, is_parH=false;
	string reg_str(".*dig*");
	int qc_core = 90;

	while ((opt = getopt_long(argc, argv, ":dnhpr::q::", longopts, NULL)) != -1) {
		switch (opt) {
		case 'd':
			is_dir = true;
			break;
		case 'r':
			reg_str.assign(optarg);
			break;
		case 'q':
			qc_core=File::str2num<int>(string(optarg));
			break;
		case 'n':
			is_parN = true;
			break;
		case 'h':
			is_parH = true;
			break;
		case 'p':
			is_parN = true; is_parH = true;
			break;
			//******
		case ':':
			cout << "need a value:\t" << char(optopt) << endl;
			return 0;
		case '?':
			cout << "unknown option:\t" << char(optopt) << endl;
			return 0;
		}
	}

	//About the other command arguments
	bool is_autoname;
	vector<string> files_in, files_out;
	vector<string>::iterator it;
	switch (argc - optind) {
	case 0://no other argument
		cout << "Not enough argument\t\t" << "ipfromdig:\t[options] srcfile [dstfile]" << endl;
		return 0;
	case 1://only one parameter:srcfile,so we output to the std
		is_autoname = false;
		break;
	case 2://have two parameter
		is_autoname = true;

		break;
	default:
		cout << "Too much argument\t\t" << endl;
		return 0;

	}
	//Judge if the argument is a directory
	if (is_dir) {
		File::browseDir(string(argv[optind]), files_in);
	}
	else {
		files_in.push_back(string(argv[optind]));
	}
	//Judge if the outname need be automatically generated
	if (is_autoname) {
		argument.assign(argv[optind + 1]);
		File::split(argument, string(";"), files_out);
	}
	else {
		File::addPrefix(files_in, files_out, std::string("report.dig"));
	}
	//Convert the digisonde file(s)
	if (files_in.empty()) {
		cout << "Not find file" << endl;
	}
	else {
		for (it = files_in.begin(); it != files_in.end(); it++) {
			//write the data
			IonPoint::ConvertByDig(*it, files_out[it - files_in.begin()], qc_core, is_parN, is_parH);
		}
	}
	
	cout << "Success converting digisonde" << endl;
	return 0;
}