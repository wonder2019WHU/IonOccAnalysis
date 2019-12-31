//IonPrf.cpp
//******************************************************************
#include "IonPrf.h"
using namespace std;

//*******************************************************************
//Global function
const string ICP::epfFlag = string("EPF");
//*******************************************************************
//Member function
//****DigStationList class****
const vector<DigStation> DigStationList::stas = {
	DigStation("WU430", 30.6, 114.4),
	DigStation("ZH466", 66.8, 123.4),
	DigStation("PSJ5J", -51.6, -57.9),
	DigStation("LV12P", -28.5, 21.2),
	DigStation("BVJ03", 2.8, -60.7),
	DigStation("AT138", 38.0, 23.5),
	DigStation("AL945", 45.1, -83.6),
	DigStation("AU930", 30.4, -97.7),
	DigStation("BP440", 40.3, 116.2),
	DigStation("SA418", 18.3, 109.6),
	DigStation("LM42B", -21.8, 114.0),
	DigStation("FZA0M", -3.8, -38.0),
	DigStation("YA462", 62.0, 129.6),
	DigStation("MH453", 53.5, 112.3),
	DigStation("RL052", 51.5, -0.6),
	DigStation("AS00Q", -8.0, -14.0),
	DigStation("CAJ2M", -22.7, -45.0),
	DigStation("TR169", 69.9, 19.2),
	DigStation("PQ052", 50.0, 14.6),
	DigStation("PA836", 34.7, -120.6),
	DigStation("DW41K", -12.5, 130.9),
	DigStation("PRJ18", 18.5, -67.1),
	DigStation("BC840", 40.0, -105.3),
	DigStation("EB040", 40.8, 0.5)
};
//****Time class****
const int Time::DOY_FDM[13] = {//nonleap year
	1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366
};
//****File class*****
void File::browseDir(string src_dir, vector<string>&dst_files) {
	//get the dir by char*
	const char *dir = src_dir.c_str();
	//open the dir
	DIR *dp;
	if ((dp = opendir(dir)) == NULL) {
		cout << "cannot open dir: " << string(dir) << endl;
		return;
	}

	//read the fir
	struct dirent *entry;
	struct stat statbuf;
	string filename, prefix_path;	
	prefix_path = src_dir.at(src_dir.size() - 1) == '/' ? (src_dir) : (src_dir + "/");

	chdir(dir);
	while ((entry = readdir(dp)) != NULL) {
		lstat(entry->d_name, &statbuf);
		filename = string(entry->d_name);
		if (S_ISDIR(statbuf.st_mode)) { // if diretory, browse it
			if (filename == "." || filename == "..") { //skip "." and ".."
				continue;
			}
			browseDir(prefix_path + filename, dst_files);
		}
		if (S_ISREG(statbuf.st_mode)) { //if regular file, judge it
				dst_files.push_back(prefix_path + filename);
		}
	}
	chdir("..");
	//close the dir
	closedir(dp);
}
void File::browseFiles(string src_file, vector<string>&dst_files) {
	//open the file lsit
	ifstream read(src_file);
	//read the filename by line
	string line;
	while(getline(read, line))
		dst_files.push_back(line);
	//close the file list
	read.close();
};
//****IonPoint class****
void IonPoint::ConvertByDig(string file_in, string file_out, int qc_core, bool isParN, bool isParH)
{
	File dig(file_in);
	string dig_path = dig.get_filepath(), dig_file = dig.get_filename();
	string dig_station, dig_str;
	//IonPoint::split_filename(file_in, dig_path, dig_file);//get the path and file of the digisonde
	ifstream in_dig(file_in);
	ofstream out_dig(file_out);
	IonPoint dig_data;
	//get the info from the title
	DigStation sta;
	dig_station = dig_file.substr(0, 5);
	if (!DigStationList::search_sta(dig_station, sta)) { //can not find
		cout << "Station is not in the known lists" << endl;
		return;
	}
	dig_data.loc.Lat = sta.Lat; dig_data.loc.Lon = sta.Lon; dig_data.par.profile = file_in;
	//get the info from the context
	stringstream ss; char delimiter;
	string str_ymd, str_dd, str_hms;
	int core; Index parN(0,20), parH(200,500);
	while (getline(in_dig, dig_str)) {
		if (IonPoint::is_valid_dig(dig_str)) {
			//split the all part
			ss << dig_str;
			ss >> str_ymd >> str_dd >> str_hms >> core >> dig_data.par.NmF2 >> dig_data.par.hmF2;
			dig_data.par.f0F2_to_NmF2();//f0F2 to NmF2
			ss.str(""); ss.clear(); //clear the stringstream
			if(core < qc_core)
				continue;
			if(isParN && !parN.is_include(dig_data.par.NmF2))
				continue;
			if(isParH && !parH.is_include(dig_data.par.hmF2))
				continue;
			//split the time part
			//**ymd
			ss << str_ymd;
			ss >> dig_data.t.Year;
			ss.str(""); ss.clear();//clear the stringstream
			//**dd
			ss << str_dd;
			ss >> delimiter >> dig_data.t.Doy;
			ss.str(""); ss.clear();//clear the stringstream
			//**hms
			ss << str_hms;
			ss >> dig_data.t.Hour >> delimiter >> dig_data.t.Min;
			ss.str(""); ss.clear();//clear the stringstream
			//write to file
			out_dig << dig_data << endl;
		}
	}
	//Free
	in_dig.close(); out_dig.close();
}
void IonPoint::ConvertByBool(vector<double>&src, vector<double>&dst, vector<bool>&_bool) {
	//calculation
	vector<double> tem;
	transform(src.begin(), src.end(), _bool.begin(), back_inserter(tem), std::multiplies<double>());
	//value
	dst.empty() ? (void(0 == 0)) : (dst.clear());
	remove_copy_if(tem.begin(), tem.end(), back_inserter(dst), bind2nd(std::equal_to<double>(), double(0)));
}
int IonPoint::BiSearchInsideLoc(vector<IonPoint>&data, IonPoint&p, IonGridSample&alw) {
	//Binary search the inside location
	vector<IonPoint>::iterator it;
	int start = 0, length = data.size();
	double interval = 0, alw_t = alw.t;

	while (length >= 1) {
		if (p.t.is_appless(p.t, (data.begin() + start)->t, alw_t) ||
			p.t.is_appgreater(p.t, (data.begin() + start + length - 1)->t, alw_t)) {
			return -1; //means can't find
		}
		it = data.begin() + start + length / 2;
		interval = p.t.cal_interval(p.t, it->t);
		if (abs(interval) <= alw_t) {
			return it - data.begin();
		}
		else {
			if (interval > alw_t) {
				start += (length + 1) / 2; length /= 2;
			}
			else {
				start = start; length /= 2;
			}
		}
	}
	return -1;

}

int IonPoint::BiSearch(vector<IonPoint>&data, IonPoint&p, int startloc, IonGridSample&alw, vector<IonGroup>&dst,
	MatchIndex&info) {
	double alw_t = alw.t, alw_lat = alw.lat, alw_lon = alw.lon;
	int mid;
	if (startloc == -1)
		mid = BiSearchInsideLoc(data, p, alw);
	else
		mid = startloc;
	//**find the mid(not precisely,just inside) location of the possible area
	if (mid == -1)
		return -1;
	else {
		int after = mid, before = mid - 1;
		vector<IonPoint>::iterator it = data.begin();
		IonGroup _ig;
		//after the mid
		while ((it + after) != data.end() && !Time::is_appless(p.t, (it + after)->t, alw_t)) {
			if (Coord::is_appequal(p.loc, (it + after)->loc, alw_lat, alw_lon) 
					&& info.is_include((it + after)->par.NmF2-p.par.NmF2, 
										(it + after)->par.hmF2-p.par.hmF2,
										(it + after)->par.AmF2, p.par.AmF2, p.par.NmF2, p.par.hmF2)) {
				_ig.set_cmp(*(it + after));
			}
			after++;
		}
		//before the mid
		while (before >= 0 && !Time::is_appgreater(p.t, (it + before)->t, alw_t)) {
			if (Coord::is_appequal(p.loc, (it + before)->loc, alw_lat, alw_lon)
					&& info.is_include((it + before)->par.NmF2-p.par.NmF2, 
										(it + before)->par.hmF2-p.par.hmF2,
										(it + before)->par.AmF2, p.par.AmF2, p.par.NmF2, p.par.hmF2)) {
				_ig.set_cmp(*(it + before));
			}
			before--;
		}
		int N = _ig.get_cmpSize();
		if (N != 0) {
			_ig.set_ref(p);
			dst.push_back(_ig);
			//cout << "Mathced:\t" << p.toString() << endl;
			return after;
		}
		else {
			//cout << "Not Matched:\t" << p.toString() << endl;
			//return mid;
			return -1;
		}
	}
}
void IonPoint::Match(vector<IonPoint>&svy, vector<IonPoint>&ref, IonGridSample&alw, vector<IonGroup>&dst,
	MatchIndex&info) {
	//Initial the vital parameter
	vector<IonPoint>::iterator it_ref(ref.begin());
	int startloc = BiSearchInsideLoc(svy, ref[0], alw);
	//For each Ionpoint in the ref
	do{
		if (it_ref == ref.end())
			break;
		startloc = BiSearch(svy, ref[it_ref - ref.begin()], startloc, alw, dst, info);
		//cout << "Finish Match\t" << it_ref->toString() << endl;
		it_ref++;
	} while (true);
}

void IonPoint::Match(string file_svy, string file_ref, IonGridSample&alw, vector<IonGroup>&dst,
	MatchIndex&info) {
	//Preprocess
	vector<IonPoint> svys, refs;
	IonPoint::MatchPreprocess(file_svy, file_ref, svys, refs);
	//Match
	cout << "****Match Start*****" << endl;
	IonPoint::Match(svys, refs, alw, dst, info);
	cout << "****Match Finish****" << endl;
}

//**********************************************************************************************************
void IonPoint::Match(vector<IonPoint>&svy, vector<IonPoint>&ref, IonGridSample&alw, IonMatch&dst,
	MatchIndex&info, STRATEGY _stra) {
	vector<IonGroup> _ig;
	IonPoint::Match(svy, ref, alw, _ig, info);
	IonPoint::CorrCoefficient(_ig, dst, _stra);
}
void IonPoint::Match(vector<IonPoint>&svy, vector<IonPoint>&ref, IonGrid&alws, vector<IonMatch>&dst,
	MatchIndex&info, STRATEGY _stra) {
	vector<IonGridSample> igs;
	vector<IonGridSample>::iterator it_igs;
	IonMatch im;
	//Get the iongridsamples
	alws.SampleGrid(igs);
	//calculation
	for (it_igs = igs.begin(); it_igs != igs.end(); it_igs++) {
		IonPoint::Match(svy, ref, *it_igs, im, info, _stra);
		dst.push_back(im);
		cout << "Finish\t" << it_igs->toString() <<" " <<im<<endl;
	}
}
void IonPoint::Match(string file_svy, string file_ref, IonGridSample&alw, IonMatch&dst,
	MatchIndex&info, STRATEGY _stra) {
	//Preprocess
	vector<IonPoint> svys, refs;
	IonPoint::MatchPreprocess(file_svy, file_ref, svys, refs);
	//Match
	cout << "****Match Start*****" << endl;
	IonPoint::Match(svys, refs, alw, dst, info, _stra);
	cout << "****Match Finish****" << endl;
}
void IonPoint::Match(string file_svy, string file_ref, IonGrid&alws, vector<IonMatch>&dst,
	MatchIndex&info, STRATEGY _stra) {
	//Preprocess
	vector<IonPoint> svys, refs;
	IonPoint::MatchPreprocess(file_svy, file_ref, svys, refs);
	//Match
	cout << "****Match Start*****" << endl;
	IonPoint::Match(svys, refs, alws, dst, info, _stra);
	cout << "****Match Finish****" << endl;
}
//**********************************************************************************************************

void IonPoint::MatchPreprocess(string file_svy, string file_ref, vector<IonPoint>&svys, vector<IonPoint>&refs) {
	//Construct svy
	cout << "****Svy data is loading*****" << endl;
	IonPoint::ReadIonPointFile(file_svy, svys);
	cout << "****Svy data is ok*****" << endl;
	//Construct the file_ref
	cout << "****Ref data is loading*****" << endl;
	IonPoint::ReadIonPointFile(file_ref, refs);
	cout << "****Ref data is ok*****" << endl;
	//Sort the svy and ref by up-order
	cout << "***All the data is sorting***" << endl;
	std::sort(svys.begin(), svys.end(), IonPoint::less_t);
	std::sort(refs.begin(), refs.end(), IonPoint::less_t);
	cout << "***Sorting is ok***" << endl;
}
void IonPoint::MoveAverage(vector<double>&src, vector<double>&dst, int _win) {
	assert(int(src.size()) >= _win);
	vector<double> _dst;
	vector<double>::iterator it_src;
	for (it_src = src.begin(); it_src != src.end() - _win + 1; it_src++)
		_dst.push_back((accumulate(it_src, it_src + _win, double(0))) / _win);
	// it_src+win(becasuse there is end,end is excess the filed by one elem)
	dst = _dst;
}
void IonPoint::LineRegression(vector<double>&src_x, vector<double>&src_y, double&dst_k, double&dst_b) {
	//The regression line is "y=kx+b"
	int N = src_x.size();
	if (N < 2 || src_x.size() != src_y.size()) {
		dst_k = std::numeric_limits<double>::quiet_NaN();
		dst_b = std::numeric_limits<double>::quiet_NaN();
		return;
	}
	//Calculate the mean
	double mean_x = accumulate(src_x.begin(), src_x.end(), double(0)) / N;
	double mean_y = accumulate(src_y.begin(), src_y.end(), double(0)) / N;
	//Calculate the sum_mul and sum_squ
	double sxy = inner_product(src_x.begin(), src_x.end(), src_y.begin(), double(0));
	double sxx = inner_product(src_x.begin(), src_x.end(), src_x.begin(), double(0));
	//Calculate the dst_a and dst_b
	dst_k = (sxy - N*mean_x*mean_y) / (sxx - N*mean_x*mean_x);
	dst_b = mean_y - dst_k*mean_x;

}
void IonPoint::CorrCoefficient(vector<double>&src_x, vector<double>&src_y, double&dst_r) {
	//inspect the size
	int N = src_x.size();
	if (N < 2 || src_x.size() != src_y.size()) {
		dst_r = std::numeric_limits<double>::quiet_NaN();
		return;
	}
	//calculate the sx,sy,sxy,sxx,syy
	double sx, sy, sxy, sxx, syy = 0;
	sx = accumulate(src_x.begin(), src_x.end(), double(0));
	sy = accumulate(src_y.begin(), src_y.end(), double(0));
	sxy = inner_product(src_x.begin(), src_x.end(), src_y.begin(), double(0));
	sxx = inner_product(src_x.begin(), src_x.end(), src_x.begin(), double(0));
	syy = inner_product(src_y.begin(), src_y.end(), src_y.begin(), double(0));
	//calculate the correlation coefficient
	dst_r = (sxy - sx*sy / N) / sqrt((sxx - sx*sx / N)*(syy - sy*sy / N));
}

void IonPoint::Cal_Bias(vector<double>&src_x, vector<double>&src_y, double&dst_am, double&dst_av, double&dst_rm, double&dst_rv) {
	//src_x is svy, src_y is ref
	int N = src_x.size();
	if (N < 1 || src_x.size() != src_y.size()) {
		dst_am = std::numeric_limits<double>::quiet_NaN();dst_av = std::numeric_limits<double>::quiet_NaN();
		dst_rm = std::numeric_limits<double>::quiet_NaN(); dst_rv = std::numeric_limits<double>::quiet_NaN();
	}
	//calculate
	vector<double> diff_ab, diff_re;
	diff_ab.resize(N); diff_re.resize(N);
	//**absolutely diff
	std::transform(src_x.begin(), src_x.end(), src_y.begin(), diff_ab.begin(), std::minus<double>());
	//**relatively diff
	std::transform(diff_ab.begin(), diff_ab.end(), src_y.begin(), diff_re.begin(), std::divides<double>());
	//**mean
	dst_am = std::accumulate(diff_ab.begin(), diff_ab.end(), double(0)) / N;
	dst_rm = std::accumulate(diff_re.begin(), diff_re.end(), double(0)) / N;
	//**variance
	std::transform(diff_ab.begin(), diff_ab.end(), diff_ab.begin(), std::bind2nd(std::minus<double>(), dst_am));
	std::transform(diff_re.begin(), diff_re.end(), diff_re.begin(), std::bind2nd(std::minus<double>(), dst_rm));
	dst_av = std::inner_product(diff_ab.begin(), diff_ab.end(), diff_ab.begin(), double(0)) / N;
	dst_rv = std::inner_product(diff_re.begin(), diff_re.end(), diff_re.begin(), double(0)) / N;
}
void IonPoint::Cal_Bias(vector<double>&src_x, vector<double>&src_y, IonBias&ib) {
	double _am, _av, _rm, _rv;
	IonPoint::Cal_Bias(src_x, src_y, _am, _av, _rm, _rv);
	ib.set(_am, _av, _rm, _rv);
}
void IonPoint::CorrCoefficient(vector<IonGroup>&src, IonMatch&dst_r, STRATEGY _stra) {
	//transform
	vector<ICP> src_x, src_y;
	IonGroup::toICP(src, src_x, src_y, _stra);
	//original code
	vector<double> x, y;
	vector<ICP>::iterator it_x, it_y;

	double rn, rh;
	IonBias bn, bh;
	int n = src_x.size();
	x.resize(n); y.resize(n);
	//corrcoef and bias for nmf2
	for (it_x = src_x.begin(), it_y = src_y.begin(); it_x != src_x.end(), it_y != src_y.end(); it_x++, it_y++) {
		x[it_x - src_x.begin()] = it_x->NmF2; y[it_y - src_y.begin()] = it_y->NmF2;
	}
	IonPoint::CorrCoefficient(x, y, rn);
	IonPoint::Cal_Bias(x, y, bn);
	//corrcoef and bias for hmf2
	for (it_x = src_x.begin(), it_y = src_y.begin(); it_x != src_x.end(), it_y != src_y.end(); it_x++, it_y++) {
		x[it_x - src_x.begin()] = it_x->hmF2; y[it_y - src_y.begin()] = it_y->hmF2;
	}
	IonPoint::CorrCoefficient(x, y, rh);
	IonPoint::Cal_Bias(x, y, bh);
	//Construct IonMatch
	dst_r.set(rn, rh, n, bn, bh);
}

//****IonPrf class****
bool IonPrf::ReadIonPrf(string _ionprf, AziIndex&msg) {
	//Add**********************************************************
	ROM _rom = File::file2ROM(_ionprf);
	switch(_rom) {
		case IONPRF: return ReadIonPrfByIONPRF(_ionprf, msg);
		case FY3C: return ReadIonPrfByFY3C(_ionprf, msg);
		default: cout<<setiosflags(ios::left)<<setw(35)<<"Error .nc format"<<endl;return false;
	}
	
}
bool  IonPrf::ReadIonPrfByIONPRF(string _ionprf, AziIndex&msg)
{
	//File address parse
	this->datafile.addressParse(_ionprf);
	//Open  the ionprf file
	netCDF::NcFile _file(_ionprf, netCDF::NcFile::read);
	if(_file.isNull()) {
		cout<<setiosflags(ios::left)<<setw(35)<<"can not open the .nc file:  "<<_ionprf<<endl;
		return false;
	}
	//Get the dimemsion
	if(_file.getDimCount() < 1) {
		cout<<setiosflags(ios::left)<<setw(35)<<"Not Enough dimension:\t"<<_ionprf<<endl;
		return false;
	}
	int Dim = _file.getDim("MSL_alt").getSize();
	//Get the data
	double*_Dens = new double[Dim], *_Alt = new double[Dim], *_Azi = new double[Dim];
	double*_Lat = new double[Dim], *_Lon = new double[Dim];
	_file.getVar("MSL_alt").getVar(_Alt); _file.getVar("ELEC_dens").getVar(_Dens); _file.getVar("OCC_azi").getVar(_Azi);
	_file.getVar("GEO_lat").getVar(_Lat); _file.getVar("GEO_lon").getVar(_Lon);
	//Construt the vetcor 
	//Array2Vector(_Alt, Dim, IonPrf::alt); Array2Vector(_Dens, Dim, IonPrf::dens);//Construct the vector 
	//**1 if the vector not empty, clear it
	this->alt.resize(Dim);this->dens.resize(Dim);
	//**2 value the dens and alt
	copy(_Alt, _Alt + Dim, this->alt.begin()); 
	copy(_Dens, _Dens + Dim, this->dens.begin());
	//**3 change units to x10^5el/cm^3
	transform(dens.begin(), dens.end(), dens.begin(), bind2nd(divides<double>(), double(1e5)));
	//Mining the time 
	string _filename = this->datafile.get_filenames();
	int pos_s = _filename.find_first_of("."), pos_e = _filename.find_last_of(".");
	_filename = _filename.substr(pos_s + 1, pos_e - pos_s + 1 - 2);
	IonPrf::ip.t.ReadTime(_filename, string("."));
	//Calculate the ICP and the location corresponding
	int max_loc; double max;
	vector<double>::iterator it = max_element(dens.begin(), dens.end());
	max = *it; max_loc = it - dens.begin();
	IonPrf::ip.par.set(max, IonPrf::alt[max_loc], OccAzi::checkAmF2(_Azi[max_loc], msg.is_abs), _ionprf);
	IonPrf::ip.loc.set(_Lat[max_loc], _Lon[max_loc]);
	//Free the memory
	delete[] _Dens; delete[] _Alt; delete[] _Azi;delete[] _Lat; delete[] _Lon;
	
	return true;
}
bool IonPrf::ReadIonPrfByFY3C(string _ionprf, AziIndex&msg) {
		//File address parse
	this->datafile.addressParse(_ionprf);
	//Open  the ionprf file
	netCDF::NcFile _file(_ionprf.c_str(), netCDF::NcFile::read);
	if(_file.isNull()) {
		cout<<setiosflags(ios::left)<<setw(35)<<"can not open .nc file:  "<<_ionprf<<endl;
		return false;
	}
	//Get the dimemsion
	//int num= _file.num_dims();
	if(_file.getDimCount() < 2) {
		cout<<setiosflags(ios::left)<<setw(35)<<"Not Enough dimension:  "<<_ionprf<<endl;
		return false;
	}
	int Dim = _file.getDim("dim_lev2a").getSize();
	if(Dim < 10){
		cout<<setiosflags(ios::left)<<setw(35)<<"Not Enough dimension:  "<<_ionprf<<endl;
		return false;
	}
	//Get the data
	double* _Dens = new double[Dim], *_Alt = new double[Dim];
	_file.getVar("MSL_alt").getVar(_Alt); 
	_file.getVar("elec_Dens").getVar(_Dens);
	//Construt the vetcor 
	//Array2Vector(_Alt, Dim, IonPrf::alt); Array2Vector(_Dens, Dim, IonPrf::dens);//Construct the vector 
	//**1 if the vector not empty, clear it
	this->alt.resize(Dim);this->dens.resize(Dim);
	//**2 value the dens and alt
	copy(_Alt, _Alt + Dim, this->alt.begin()); 
	copy(_Dens, _Dens + Dim, this->dens.begin());
	//**3 change units to x10^5el/cm^3
	transform(dens.begin(), dens.end(), dens.begin(), bind2nd(divides<double>(), double(1e5)));
	//Mining the time 
	string _time = "", _timestr = "";
	_file.getAtt("Year").getValues(_timestr);_time += ( _timestr + string("  "));
	_file.getAtt("Doy").getValues(_timestr);_time += (_timestr + string("  "));
	_file.getAtt("Hour").getValues(_timestr);_time += (_timestr + string("  "));
	_file.getAtt("Minute").getValues(_timestr);_time += (_timestr + string("  "));
	IonPrf::ip.t.ReadTime(_time);
	//Mining the coordinate
	double _lat, _lon;
	string _locstr = "";
	_file.getAtt("Lat").getValues(_locstr);_lat = File::str2num<double>(_locstr);
	_file.getAtt("Lon").getValues(_locstr);_lon = File::str2num<double>(_locstr);
	if(abs(_lat) > 90 || abs(_lon) > 180) {
		cout<<setiosflags(ios::left)<<setw(35)<<"Location over the boundary:  "<<_ionprf<<endl;
		return false;
	}
	IonPrf::ip.loc.set(_lat, _lon);
	//Calculate the ICP and the location corresponding
	int max_loc; double max;
	vector<double>::iterator it = max_element(dens.begin(), dens.end());
	max = *it; max_loc = it - dens.begin();
	OccAzi oa(_ionprf, IonPrf::alt[max_loc], msg);
	IonPrf::ip.par.set(max, IonPrf::alt[max_loc], oa.getAmF2(), _ionprf);

	//Free the memory
	delete[] _Dens; delete[] _Alt; 
	return true;
	
}
string IonPrf::WriteReport(IonIndex&qc) {
	//Basic string
	string basic_str = this->ip.toString();
	//Index string
	string index_str;
	double num_md = 0, num_sigma = 0, num_gk = 0, num_lk = 0;
	ios_base::fmtflags flags = ios::fixed | ios::internal;
	vector<double> index_num;
	//md and sigma
	if (qc.md.get_enable() || qc.sigma.get_enable()) {
		this->cal_fluctuation(num_md, num_sigma,qc);
		if (qc.md.get_enable()) {
			index_str += File::num2str<double>(num_md, flags);
			index_num.push_back(num_md);
		}
		if (qc.sigma.get_enable()) {
			index_str += File::num2str<double>(num_sigma, flags);
			index_num.push_back(num_sigma);
		}
	}
	//gk and lk
	flags = ios::fixed | ios::internal | ios::showpos;
	if (qc.gk.get_enable() || qc.lk.get_enable()) {
		this->cal_regression(num_gk, num_lk,qc);
		if (qc.gk.get_enable()) {
			index_str += File::num2str<double>(num_gk, flags);
			index_num.push_back(num_gk);

		}
		if (qc.lk.get_enable()) {
			index_str += File::num2str<double>(num_lk, flags);
			index_num.push_back(num_lk);
		}
	}
	//parN, parh, parA--becauase we have calculate before,so we just put it into the index_num
	if (qc.parN.get_enable()) {
		index_num.push_back(this->ip.par.get_NmF2());
	}
	if (qc.parh.get_enable()) {
		index_num.push_back(this->ip.par.get_hmF2());
	}
	if(qc.parA.get_enable()) {
		index_num.push_back(this->ip.par.get_AmF2());
	}
	
	//qc-score (1-true, 0-false)
	string qc_str = qc.is_include(index_num)?("  1"):("  0");
	
	//judge whether ignore invalid value
	string result = basic_str + index_str + qc_str;
	if (!qc.is_accept_invalid && File::is_include_invalid(result)) { //include the invalid and not accept the invalid
		cout <<setiosflags(ios::left)<<setw(35)<< "Find NaN or Inf in this file:  " << result<<endl;//display the error
		result.clear();
		return result;
	}
	if (!qc.is_accept_qcfail && !qc.is_include(index_num)) { //find the qcfail and not accept the qcfail
		cout <<setiosflags(ios::left)<<setw(35)<< "Find Fail Index in this file:  "<< result << endl;//display the error
		result.clear();
		return result;
	}
	return result;
}
bool IonPrf::WriteReport(IonIndex&qc, ostream&out) {
	string rep_str = this->WriteReport(qc);
	if(!rep_str.empty()) {//not empty can be write
		out << rep_str << endl;
		return true;
	} 
	return false;
		
}
void IonPrf::WriteInfo(string file_out) {
	ofstream out(file_out);
	//Head
	//**Time
	out << setiosflags(ios::left)
		<< setw(60) << this->ip.t.toString()
		<< setw(20) << string("Time (Y-D-H-M)")
		<< resetiosflags(ios::left) << endl;
	//**Coord
	out << setiosflags(ios::left)
		<< setw(60) << this->ip.loc.toString()
		<< setw(20) << string("Loc (Lat-Lon)")
		<< resetiosflags(ios::left) << endl;
	//**ICP
	out << setiosflags(ios::left)
		<< setw(60) << this->ip.par.toSimpleString()
		<< setw(20) << string("ICP (NmF2-hmF2-AmF2)")
		<< resetiosflags(ios::left) << endl;
	//**Directory
	out << setiosflags(ios::left)
		<< setw(60) << this->datafile.get_filepath()
		<< setw(20) << string("Directory")
		<< resetiosflags(ios::left) << endl;
	//**FileName
	out << setiosflags(ios::left)
		<< setw(60) << this->datafile.get_filenames()
		<< setw(20) << string("Filenames")
		<< resetiosflags(ios::left) << endl;
	//**DataInfo
	out << setiosflags(ios::left)
		<< setw(60) << string("Density (x10^5el/cm^-3)    Altitude (Km)")
		<< setw(20) << string("DataFormat")
		<< resetiosflags(ios::left) << endl;
	//**End of Head
	out << setiosflags(ios::left)
		<< setw(60) << string("")
		<< setw(20) << string("End of Head")
		<< resetiosflags(ios::left) << endl;;
	//Data
	vector<double>::iterator it_dens(dens.begin()), it_alt(alt.begin());
	for (; it_dens != dens.end(), it_alt != alt.end(); it_dens++, it_alt++) {
		out << setiosflags(ios::fixed | ios::internal | ios::showpos)
			<< setprecision(6) << setfill('0')
			<< setw(10) << *it_dens << "  "
			<< setw(11) << *it_alt << " "
			<< resetiosflags(ios::fixed | ios::internal | ios::showpos) << endl;
	}

	out.close();

}

void IonPrf::cal_fluctuation(double&dst_MD, double&dst_sigma, IonIndex&qc)
{
	assert(!IonPrf::is_empty());//IonPrf must not be null
	vector<double> mov_ave, single_ssd, single_srd;
	//calculate the ssd
	double ssd = 0, srd = 0; int _win = qc.winsize;
	IonPoint::MoveAverage(IonPrf::dens, mov_ave, _win);
	transform(mov_ave.begin(), mov_ave.end(), dens.begin() + _win / 2,
		back_inserter(single_ssd), IonPoint::op_square_diff);//get the single ssd
	transform(mov_ave.begin(), mov_ave.end(), dens.begin() + _win / 2,
		back_inserter(single_srd), IonPoint::op_relat_diff);//get the single srd
	ssd = accumulate(single_ssd.begin(), single_ssd.end(), double(0));//get the ssd
	srd = accumulate(single_srd.begin(), single_srd.end(), double(0));//get the srd

	//calculate the MD and sigma
	dst_sigma = sqrt(ssd / mov_ave.size()) / IonPrf::ip.par.get_NmF2();
	dst_MD = srd / mov_ave.size();
	
	//limit the extremely value
	if(dst_sigma > 999)dst_sigma = 999;
	if(dst_MD > 999)dst_MD = 999;

}
void IonPrf::cal_regression(double&dst_gk, double&dst_lk, IonIndex&qc)
{
	double tem_b = 0, ls = qc.localregress.get_min(), le =qc.localregress.get_max();
	vector<double> _x, _y; 
	vector<bool> _bool;
	//Calculate the dst_gk
	transform(IonPrf::alt.begin(), IonPrf::alt.end(),back_inserter(_bool), 
			bind2nd(std::greater_equal<double>(), IonPrf::ip.par.get_hmF2()));
	IonPoint::ConvertByBool(dens, _y, _bool);
	IonPoint::ConvertByBool(alt, _x, _bool);
	IonPoint::LineRegression(_x, _y, dst_gk, tem_b);
	//Calculate the dst_lk
	vector<bool> _bool1, _bool2;
	transform(IonPrf::alt.begin(), IonPrf::alt.end(), back_inserter(_bool1),
			bind2nd(std::greater_equal<double>(), ls));
	transform(IonPrf::alt.begin(), IonPrf::alt.end(), back_inserter(_bool2),
			bind2nd(std::less_equal<double>(), le));
	transform(_bool1.begin(), _bool1.end(), _bool2.begin(), _bool.begin(), std::logical_and<bool>());
	IonPoint::ConvertByBool(dens, _y, _bool);
	IonPoint::ConvertByBool(alt, _x, _bool);
	IonPoint::LineRegression(_x, _y, dst_lk, tem_b);

}

//class OccAzi
string OccAzi::creIonPhs(std::string _copath, std::string _subdir) {
	ROM _rom = File::file2ROM(this->datafile.get_filename());
	switch(_rom) {
		case IONPRF: return creIonPhsByIONPRF(_copath, _subdir);
		case FY3C: return creIonPhsByFY3C(_copath, _subdir);
		default:return string("");
	}
};

string OccAzi::creIonPhsByIONPRF(std::string _copath, std::string _subdir) {
	//split the info
	int pos = this->datafile.get_filename().find_first_of("_");
	string base = this->datafile.get_filenames().substr(pos+1);
	
	//check the path
	if(_copath[_copath.size()-1] != '/')
		_copath += string("/");
	string prefix = (_copath==string("#"))?(this->datafile.get_filepath()+string("/")):(_copath);
	
	//check the subpdir
	//PreProcess	
        vector<string> info;
        File::split(this->datafile.get_filename(), string("."), info);
        string year_str = info[1];
        string doy_str = info[2];
        //SubDir
        string dir;
	switch(_subdir[_subdir.size()-1]) {	
	  case '*': //Standard SubDir
                dir = prefix + year_str + string("/") + year_str + string(".") + doy_str + string("/");
                break;
          case '|': //Half Standard SubDir
                dir = prefix + year_str + string(".") + doy_str + string("/");
                break; 
          case '#': //None SubDir
                dir = prefix;
                break;
          default: //Default to None SubDir
                dir = prefix;
                break;
	}
	
	//epoch ionphs filename
	return dir + string("ionPhs_") + base;
};
string OccAzi::creIonPhsByFY3C(std::string _copath, std::string _subdir) {
	//split the info
	vector<string> info;
	File::split(this->datafile.get_filename(), string("_"), info);
	string time_str = (info[7] + string("_") + info[8] + string("_"));
	string sat_str = (string("IE") + info[9].substr(0,3) + string("_"));
	
	//check the path
	if(_copath[_copath.size()-1] != '/')
		_copath += string("/");
	string prefix = (_copath==string("#"))?(this->datafile.get_filepath()+string("/")):(_copath);
	
	//check the subpdir
	//PreProcess
	string year_str = info[7].substr(0,4);
        int year = File::str2num<int>(year_str);
        int mon = File::str2num<int>(info[7].substr(4,2));
        int day = File::str2num<int>(info[7].substr(6,2));		
        int doy = Time::cal_doy(year, mon, day);
        std::ios_base::fmtflags flags(ios::right);
        string doy_str = (File::num2str<int>(doy, flags, 3, 0, '0')).substr(0, 3);
        //SubDir
        string dir;
	switch(_subdir[_subdir.size()-1]) {
          case '*': //Standard SubDir
		dir = prefix + year_str + string("/") + year_str + string(".") + doy_str + string("/");
                break;
          case '|': //Half Sttandard SubDir
                dir = prefix + year_str + string(".") + doy_str + string("/");
                break;
          case '#': //None SubDir
                dir = prefix;
                break;
          default: //Default to None SubDir
                dir = prefix;
                break;
	}
	//epoch ionphs filename
	return dir + string("FY3C_GNOSX_GBAL_L1_") + time_str + sat_str + string("MS.NC");
};


void OccAzi::set(string _ionprf, double _hmF2, AziIndex&msg) {
	//get ionprf
	this->datafile.addressParse(_ionprf);
	//get ionphs
	string _ionphs = this->creIonPhs(msg.copath, msg.subdir);
	if(!File::is_existFile(_ionphs)) {
		this->AmF2 = 999; return;
	}
	//read ionphs
	double EARTH_R = 6378.137;
	netCDF::NcFile _file(_ionphs, netCDF::NcFile::read);
	if(_file.isNull() || _file.getDimCount() < 1) {
		this->AmF2 = 999; return;
	}
		
	int _dim = _file.getDim("time").getSize();
	double *xLeo = new double[_dim], *yLeo = new double[_dim], *zLeo = new double[_dim];
	double *xGps = new double[_dim], *yGps = new double[_dim], *zGps = new double[_dim];
		
	_file.getVar("xLeo").getVar(xLeo); _file.getVar("yLeo").getVar(yLeo); _file.getVar("zLeo").getVar(zLeo);
	_file.getVar("xGps").getVar(xGps); _file.getVar("yGps").getVar(yGps); _file.getVar("zGps").getVar(zGps);
	
	Cartesian rL, rG, rH;
	Polar pL, pG, pH;
	double dx = 0, dy = 0;
	vector<double> occ_height;
	vector<double> occ_azi;
	
	for(int i = 0; i < _dim; i++) {
		rL = Cartesian(xLeo[i], yLeo[i], zLeo[i]);
		rG = Cartesian(xGps[i], yGps[i], zGps[i]);
		rH = rL.triHeight(rG);
		if(rH.cross(rL).calDot(rH.cross(rG)) > 0)
			continue;
		occ_height.push_back(rH.getR()-EARTH_R);
		pL = rL.toPolar(); pG = rG.toPolar(); pH = rH.toPolar();
		pL.setLambda(Coordinate::azimuthZS2ZM(pL.getLambda()));
		pG.setLambda(Coordinate::azimuthZS2ZM(pG.getLambda()));
		pH.setLambda(Coordinate::azimuthZS2ZM(pH.getLambda()));
		//method3-sphere azimuth-rH~rG
		double p = SIND(pG.getLambda()-pH.getLambda())*COSD(pG.getPhi());
		double q = COSD(pH.getPhi())*SIND(pG.getPhi())-SIND(pH.getPhi())*COSD(pG.getPhi())*COSD(pG.getLambda()-pH.getLambda());
		occ_azi.push_back(Coordinate::calAzimuth(p, q, false, false, false));
	}	
	delete[] xLeo, yLeo, zLeo;
	delete[] xGps, yGps, zGps;
	
	//get AmF2
	double min_height = *min_element(occ_height.begin(), occ_height.end());
	double max_height = *max_element(occ_height.begin(), occ_height.end());
	if(_hmF2 < min_height || _hmF2 > max_height) {
		this->AmF2 = 999; return;
	}
	transform(occ_height.begin(), occ_height.end(), occ_height.begin(),
		bind(minus<double>(), placeholders::_1, _hmF2));
	transform(occ_height.begin(), occ_height.end(), occ_height.begin(),
		occ_height.begin(), multiplies<double>());
	int loc = min_element(occ_height.begin(), occ_height.end())-occ_height.begin();
	this->AmF2 = OccAzi::checkAmF2(occ_azi[loc], msg.is_abs);
}





