#include<getopt.h>
#include<iostream>
#include<unistd.h>
#include "IonPrf.h"
using namespace std;

int main(int argc, char* argv[]) {
	int opt;
	string argument;
	bool is_dir = false;
	string reg_str(".*[nN][cC]$");
	AziIndex msg;
	struct option longopts[] = {
		{ "dirent", 0, NULL, 'd' },
		{"regex", 1, NULL, 'r'},
		{"amf2", 2, NULL, 'a'},
		//*****
		{ 0, 0, 0, 0 }
	};
	while ((opt = getopt_long(argc, argv, ":dr::", longopts, NULL)) != -1) {
		switch (opt) {
		case 'd':
			is_dir = true;
			break;
		case 'r':
			reg_str.assign(optarg);
			break;
		case 'a':
			if(optarg) {
				argument.assign(optarg);
				msg.set(argument);
			}
			break;
		//*********
		case ':':
			cout << "need a value:\t" << char(optopt) << endl;
			return 0;
		case '?':
			cout << "unknown option:\t" << char(optopt) << endl;
			return 0;
		}
	
	}

	//other parameter
	bool is_outdir;
	string file_in;
	string file_out;
	vector<string> myfiles;
	vector<string>::iterator it;
	IonPrf myfile;

	switch (argc - optind) {
	case 0://no other argument
		cout << "Not enough argument\t\t" << "ip2rep:\t[options] srcfile [dstfile]" << endl;
		return 0;
	case 1://only one parameter:we set the default outfilename
		is_outdir = false;
		break;
	case 2://have two parameter
		is_outdir = true;
		file_out.assign(argv[optind + 1]);
		if (!File::is_existDir(file_out)) {
			cout << "Not Exist Dir:\t\t" << file_out << endl;
			return 0;
		}
		break;
	default:
		cout << "Too much argument\t\t" << endl;
		return 0;
	}
	file_in.assign(argv[optind]);

	if (is_dir) {
		File::browseDir(file_in, myfiles);
	}
	else {
		myfiles.push_back(file_in);
	}

	//read and write
	int count_all = myfiles.size();
	int count_rc = 0;
	if (myfiles.empty()) {
		cout << "Not find such file\t\t" << endl;
		return 0;
	}
	else {
		for (it = myfiles.begin(); it != myfiles.end(); it++) {
			//write the data
			if(! myfile.ReadIonPrf(*it, msg)) {
				count_rc++;
				continue;
			}
			
			if (is_outdir) {
				myfile.WriteInfo(myfile.generate_filename_txt(file_out));
			}
			else {
				myfile.WriteInfo();
			}
			//cout << "Finish:\t" << *it << endl;
		}
	}
	cout << "***Finish writing information:\t" 
	     <<count_all<<"(all), "
		 <<count_rc<<"(rc), "
		 <<"***"<<endl;
	return 0;
}
