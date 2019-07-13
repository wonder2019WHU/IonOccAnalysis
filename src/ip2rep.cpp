#include<getopt.h>
#include<iostream>
#include<unistd.h>
#include "IonPrf.h"
using namespace std;

int main(int argc, char* argv[]) {
	int opt;
	string argument;
	IonIndex qc;
	AziIndex msg;
	qc.set_switch_all(false); qc.set_switch_invalid(true); qc.set_switch_qcfail(true);
	bool is_dir = false, is_files = false;
	string reg_str(".*[nN][cC]$");
	string rom_str("COSMIC");
	string copath("#");
	struct option longopts[] = {
		{ "md", 2, NULL, 'm' },
		{ "sigma", 2, NULL, 's' },
		{ "gk", 0, NULL, 'g' },
		{ "lk", 2, NULL, 'l' },
		{ "whole", 2, NULL, 'w' },
		{"dirent", 0, NULL, 'd'},
		{"ignore", 0, NULL, 'i'},
		{"quality", 2, NULL, 'q'},
		{"regex", 1, NULL, 'r'},
		{"files", 0, NULL, 'f'},
		{"amf2", 2, NULL, 'a'},
		//*****
		{ 0, 0, 0, 0 }
	};
	while ((opt = getopt_long(argc, argv, ":m::s::gl::w::diq::r::fa::", longopts, NULL)) != -1) {
		switch (opt) {
		case 'm':
			qc.set_switch_md(true);//open the md-switch
			if (optarg) { //have parameter
				argument.assign(optarg);//get the argument
				qc.set_par_win(File::str2num<int>(argument));//write the parameter
			}
			// not have parameter, we set the defaul(winsize = 9)
			break;
		case 's':
			qc.set_switch_sigma(true);//open the sigma-switch
			if (optarg) { //have parameter
				argument.assign(optarg);//get the argument
				qc.set_par_win(File::str2num<int>(argument));
			}
			// not have parameter, we set the defaul(winsize = 9)
			break;
		case 'g':
			qc.set_switch_gk(true);//open the gk-switch
			break;
		case 'l':
			qc.set_switch_lk(true);//open the lk-switch
			if (optarg) { //have parameter
				argument.assign(optarg);//get the argument
				qc.set_par_lri(argument);
			}
			// not have parameter, we set the defaul(localregress = [420, 490])
			break;
		case 'w':
			qc.set_switch_all(true);//open the all switch
			if (optarg) { //have parameter
				argument.assign(optarg);
				qc.set_par_all(argument);
			}
			// not have parameter, we set the defaul(winsize = 9, localregress = [420, 490])
			break;
		case 'd':
			is_dir = true;
			break;
		case 'i':
			qc.set_switch_invalid(false);//ignore invalid
			break;
		case 'q':
			qc.set_switch_invalid(false);//ignore invalid
			qc.set_switch_qcfail(false);//ignore qcfail
			if (optarg) { //have parameter
				argument.assign(optarg);
				qc.set_par_qc(argument);
			}
			break;
			// not have parameter, we set the defaul(md = [0,0.1], sigma = [0, 0.05], gk = (<0), lk = (<0))
		case 'r':
			reg_str.assign(optarg);
			break;
		case 'f':
			is_files = true;
			break;
		case 'a':
			qc.set_switch_parA(true);//open the parA-switch
			if (optarg) {//have parameter
				argument.assign(optarg);
				msg.set(argument);
			}
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
	//get the other argument
	IonPrf myfile;
	string rep;
	vector<string> myfiles;
	vector<string>::iterator it;
	ofstream out;
	bool is_outfile;

	switch (argc - optind) {
	case 0://no other argument
		cout << "Not enough argument\t\t" << "ip2rep:\t[options] srcfile [dstfile]" << endl;
		return 0;
	case 1://only one parameter:srcfile,so we output to the std
		is_outfile = false;
		break;
	case 2://have two parameter
		is_outfile = true;
		out.open(string(argv[optind + 1]));
		break;
	default:
		cout << "Too much argument\t\t" << endl;
		return 0;
	}

	//Judge if the argument is a directory
	if (is_dir) {
		File::browseDir(string(argv[optind]), myfiles);
	}
	else {
		if(is_files)
			File::browseFiles(string(argv[optind]), myfiles);
		else		
			myfiles.push_back(string(argv[optind]));
	}
	cout<<"***Finish Files Search***"<<endl;
	//Read the ionprf file(s)
	int count_all = myfiles.size();
	int count_rc = 0, count_qc = 0;
	if (myfiles.empty()) {
		cout << "Not find .nc file" << endl;
	}
	else {
		for (it = myfiles.begin(); it != myfiles.end(); it++) {			
			//write the data
			//read control
			if(! myfile.ReadIonPrf(*it, msg)) {
				count_rc ++;
				continue;
			}
			//quality control
			if(is_outfile) {
				if(!myfile.WriteReport(qc, out))
					count_qc++;
			}
			else {
				if(!myfile.WriteReport(qc, cout))
					count_qc++;
			}
			cout<<"Finish: "<<*it<<endl;
		}
	
	}
	cout << "***Finish writing report:\t" 
	     <<count_all<<"(all), "
	     <<count_rc<<"(rc), "
		 <<count_qc<<"(qc), "
		 <<"***"<<endl;

	return 0;
}
