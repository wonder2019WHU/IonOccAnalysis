#include<getopt.h>
#include<iostream>
#include<unistd.h>
#include "IonPrf.h"
using namespace std;

int main(int argc, char* argv[]) {
	int opt;
	string argument;
	IonGrid alws;
	IonGridSample alw;
	bool is_ig = false, is_igsp = true, is_corr = false, is_icp = false;
	bool is_file_corr = false, is_file_icp = false;
	ofstream out_corr, out_icp;
	string str_stra("AVERAGE");
	MatchIndex info;
	struct option longopts[] = {
		{ "grid", 2, NULL, 'g' },
		//grid="t_start t_end [ept_s][ept_e][rvs] t_step; 
		//lat_start lat_end [ept_s][ept_e][rvs] lat_step; 
		//lon_start lon_end [ept_s][ept_e][rvs] lon_step"
		{ "sample", 2, NULL, 's'},//sample="alw_t  alw_lat  alw_lon"
		{ "corrcoef", 2, NULL, 'c'},//corrcoef="file"
		{ "icp", 2, NULL, 'p'},//icp="file"
		{"strategy", 1, NULL, 't'},//strategy="type"
		{"nmf2maxdiff", 1, NULL, 'n'},//nmf2maxdiff="maxdiff"
		{"hmf2maxdiff", 1, NULL, 'h'},//hmf2maxdiff="maxdiff"
		{"amf2maxdiff", 1, NULL, 'a'},//amf2maxdiff="maxdiff"
		{"amf2", 1, NULL, 'A'},//amf2 in "start~end"
		{"Rnmf2maxdiff", 1, NULL, 'N'},//Relative nmf2maxdiff="maxdiff",%
		{"Rhmf2maxdiff", 1, NULL, 'H'},//Relative hmf2maxdiff="maxdiff",%
		{"abs", 0, NULL, 'b'},//amf2 is set abs
		//******
		{ 0, 0, 0, 0 }
	};
	while ((opt = getopt_long(argc, argv, ":g::s::c::p::t::n::h::p::a::A::N::H::b", longopts, NULL)) != -1) {
		switch (opt) {
		case 'g':
			is_ig = true; is_igsp = false; //The [s] and [g] can only use once, the latter will overwrite the former
			if (optarg) {
				argument.assign(optarg);
				alws.parse(argument);
			}
			break;
		case 's':
			is_igsp = true; is_ig = false; //The [s] and [g] can only use once, the latter will overwrite the former
			if (optarg) {
				argument.assign(optarg);
				alw.parse(argument);
			}
			break;
		case 'c':
			is_corr = true;
			if (optarg) {
				argument.assign(optarg);
				out_corr.open(argument);
				is_file_corr = true;
			}	
			break;
		case 'p':
			is_icp = true;
			if (optarg) {
				argument.assign(optarg);
				out_icp.open(argument);
				is_file_icp = true;
			}
			break;
		case 't':
			str_stra.assign(optarg);
			break;
		case 'n':
			argument.assign(optarg);
			info.set_qcdN(argument);
			break;
		case 'h':
			argument.assign(optarg);
			info.set_qcdH(argument);
			break;
		case 'a':
			argument.assign(optarg);
			info.set_qcdA(argument);
			break;
		case 'A':
			argument.assign(optarg);
			info.set_qcA(argument);
			break;
		case 'N':
			argument.assign(optarg);
			info.set_qcRdN(argument);
			break;
		case 'H':
			argument.assign(optarg);
			info.set_qcRdH(argument);
			break;
		case 'b':
			info.set_is_abs(true);
			break;
			//*********************
		case ':':
			cout << "need a value:\t" << char(optopt) << endl;
			return 0;
		case '?':
			cout << "unknown option:\t" << char(optopt) << endl;
			return 0;
		}	
	}

	//other parameter
	string file_svy, file_ref;
	switch (argc - optind) {
	case 0://no other argument
	case 1://only one parameter
		cout << "Not enough argument\t\t" << "ipmatch:\t[options] svyfile reffile" << endl;
		return 0;
	case 2://have two parameter
		file_svy.assign(argv[optind]); file_ref.assign(argv[optind + 1]);
		break;
	default:
		cout << "Too much argument\t\t" << endl;
		return 0;
	}
	//calculation
	vector<IonGroup> dst;
	IonMatch im;
	vector<IonMatch> ims;
	vector<ICP> dst_svy, dst_ref;
	STRATEGY _stra = File::str2STRATEGY(str_stra);
	if (is_igsp) {  //we set a IonGridSample
		if (!is_icp && !is_corr) { //none parameter
			cout << "You must set at least one match option\t\t" << "ipmatch:\t[options] svyfile reffile" << endl;
			return 0;
		}
		IonPoint::Match(file_svy, file_ref, alw, dst, info);
		if (is_icp) { 
		    IonGroup::toICP(dst, dst_svy, dst_ref, _stra);//cout<<dst_svy.size()<<" , "<<dst_ref.size()<<endl;
			is_file_icp ? (IonPoint::WriteMatchICP(out_icp, dst_svy, dst_ref)) 
						: (IonPoint::WriteMatchICP(cout, dst_svy, dst_ref));
		}
		if (is_corr) { 
			IonPoint::CorrCoefficient(dst, im, _stra);
			is_file_corr ? (IonPoint::WriteMatchCorrCoef(out_corr, alw, im))
						: (IonPoint::WriteMatchCorrCoef(cout, alw, im));
		}

	}

	if (is_ig) { //we set a IonGrid
		if (!is_icp && !is_corr) { //none parameter
			cout << "You must set at least one match option\t\t" << "ipmatch:\t[options] svyfile reffile" << endl;
			return 0;
		}
		if (is_icp) {
			cout << "When [g] and [p] at the same time, you must set parameter for [c]" << endl;
			return 0;
		}
		if (is_corr) {
			IonPoint::Match(file_svy, file_ref, alws, ims, info, _stra);
			is_file_corr ? (IonPoint::WriteMatchCorrCoef(out_corr, alws, ims))
						: (IonPoint::WriteMatchCorrCoef(cout, alws, ims));
		
		}
	
	}
	
	int count_group = dst.size();
	int count_couple = dst_svy.size();
	int count_window = (is_igsp)?(1):(ims.size());
	cout << "***finish matching:\t"
	     <<count_group<<"(group), "
		 <<count_couple<<"(couple), "
		 <<count_window<<"(window), "
	     <<"***" << endl;
	return 0;
}
