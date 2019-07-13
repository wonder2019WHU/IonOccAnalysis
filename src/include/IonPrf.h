//IonPrf.h
//******************************************************************
//Author
//Time
//Function
//Instruction
//******************************************************************
#ifndef IONPRF_H_INCLUDE
#define IONPRF_H_INCLUDE
#include<cstdlib>
#include<cmath>
#include<iostream>
#include<fstream>
#include<sstream>
#include<assert.h>
#include<string>
#include<iomanip>
#include<vector>
#include<iterator>
#include<algorithm>
#include<functional>
#include<numeric>
#include<unistd.h>
#include<sys/stat.h>
#include<dirent.h>
#include<netcdf>
#include"Coordinate.h"
//#include"netcdfcpp.h"
//******************************************************************
enum ROM {IONPRF, FY3C, ERRORTYPE};
enum STRATEGY {AVERAGE, EQUALITY, OPTIMAL, ERRORSTRA};
class IonPoint;
class IonPrf;
class IonGroup;
//****First class
class File {
protected:
	std::string filepath;
	std::string filename;
	std::string fileformat;
public:
	File() {};
	File(std::string address) { addressParse(address); }
	~File() {};
	void addressParse(std::string address) {
		//Initial format
		int pos = address.find_last_of("."); //The position of "."
		this->fileformat = address.substr(pos + 1, address.size() - pos);
		address = address.substr(0, pos);
		//Initial filename
		pos = address.find_last_of("/"); //The position of "/"(dirent seperator)
		this->filename = address.substr(pos + 1, address.size() - pos);
		address = address.substr(0, pos);
		//Initial filepath
		this->filepath = address;
	}
	const std::string get_filepath()const { return filepath; }
	const std::string get_filename()const { return filename; }
	const std::string get_fileformat()const { return fileformat; }
	const std::string get_filenames()const { return filename + "." + fileformat;}
	const std::string get_address()const {
		return filepath + std::string("/") + filename + std::string(".") + fileformat;
	}
	std::string toString() { std::ostringstream ss; ss << *this; return ss.str(); }
	std::string replace_filepath(std::string src_dir) {
		using namespace std;
		string prefix_path = src_dir.at(src_dir.size() - 1) == '/' ? (src_dir) : (src_dir + "/");
		return prefix_path + this->filename + string(".") + fileformat;
	}
	//state judgement
	bool is_empty() { return filepath.empty() || filename.empty() || fileformat.empty(); }
	//***using "||" is for only when all the string is not empty, we can use the class formally
	//Declare frined
	friend std::ostream& operator << (std::ostream&os, const File&p) {  //not member funtion
		os << p.filepath << " , " << p.filename << " , " << p.fileformat << " ";
		return os;
	}
	//static function
	static bool is_existFile(std::string src_file) {
		std::ifstream read(src_file);
		bool dst = (read)?(true):(false);
		read.close();
		return dst;
	}
	static bool is_existDir(std::string src_dir) {
		//get the dir by char*
		const char *dir = src_dir.c_str();
		DIR*dp = opendir(dir);
		bool is_exist = (dp == NULL) ? (false) : (true);
		closedir(dp);
		return is_exist;
	}
	static void browseDir(std::string src_dir, std::vector<std::string>&dst_files);
	static void browseFiles(std::string src_file, std::vector<std::string>&dst_files);
	static void split(std::string src, std::string delimiters, std::vector<std::string>& dst) {
		dst.empty() ? (void(0 == 0)) : (dst.clear());
		//split the src by delimiters
		if ("" == src || "" == delimiters)
			dst.push_back(src);//get the whole src when src or delimites is null
		else {
			std::string _src = src + delimiters;
			int pos = _src.find(delimiters);
			while (_src.npos != pos) { //When find the pos
				dst.push_back(_src.substr(0, pos));
				_src = _src.substr(pos + 1, _src.size() - pos);
				pos = _src.find(delimiters);
			}
		}
	}
	static void addPrefix(std::vector<std::string>&file_src, std::vector<std::string>&file_dst, std::string prefix) {
		using namespace std;
		vector<string>::iterator it(file_src.begin());
		File p;
		for (; it != file_src.end(); it++) {
			p.addressParse(*it);
			p.filename = prefix + p.filename;
			file_dst.push_back(p.get_address());
		}
	}
	template<class T>
	static void str2num(std::vector<std::string>&strs, std::vector<T>&dst) {
		using namespace std;
		dst.empty() ? (void(0 == 0)) : (dst.clear());
		vector<string>::iterator it;
		stringstream ss; T num = 0;
		for (it = strs.begin(); it != strs.end(); it++) {
			ss << *it; ss >> num; //value
			dst.push_back(num);//push number
			ss.str(""); ss.clear();//initial the stream
		}
	}
	template<class T>
	static T str2num(std::string str) { std::stringstream ss; T num; ss << str; ss >> num; return num; }
	template<class T>
	static std::string num2str(T num) { std::stringstream ss; ss << num; return ss.str(); }
	template<class T>
	static std::string num2str(T num, std::ios_base::fmtflags flags,
		int w = 7, int prec = 3, char fill = '0') {
		//transform to str by self-define way
		using namespace std;
		stringstream ss;
		ss << setiosflags(flags)
			<< setw(w) << setprecision(prec) << setfill(fill)
			<< num << " "
			<< resetiosflags(flags);
		return ss.str();
	}
	static ROM file2ROM(const std::string&str) {
		using namespace std;
		File p(str);
		int pos = p.get_filename().find_first_of("_");
		string _rom = p.get_filename().substr(0, pos);
		if(_rom == string("ionPrf"))
			return IONPRF;
		if(_rom == std::string("FY3C"))
			return FY3C;
		
		return ERRORTYPE;
	}
	static STRATEGY str2STRATEGY(const std::string&str) {
		using namespace std;
		if(str == string("AVERAGE"))
			return AVERAGE;
		if(str == string("EQUALITY"))
			return EQUALITY;
		if(str == string("OPTIMAL"))
			return OPTIMAL;
		return ERRORSTRA;
	}
	static bool is_num(std::string str) {
		std::stringstream os(str); double num;
		return os >> num ? true : false;
	};
	static bool is_include_nan(std::string str) {
		return (str.find("1.#QO") != std::string::npos ||
			str.find("nan") != std::string::npos) ?  (true) : (false);
	}
	static bool is_include_inf(std::string str) {
		return (str.find("1.#IO") != std::string::npos ||
			str.find("inf") != std::string::npos) ? (true) : (false);
	}
	static bool is_include_invalid(std::string str) {
		return is_include_nan(str) || is_include_inf(str);
	}
};
class IntervalPoint {
protected:
	double value;
	bool except;
public:
	IntervalPoint() :value(0), except(false){};
	IntervalPoint(double _value, bool _ept = false) :value(_value), except(_ept){};
	~IntervalPoint(){};
	double get_value()const{ return value; }
	bool get_state()const{ return except; }
	void set(double _value, bool _ept = false) {
		this->value = _value;
		this->except = _ept;
	}
	bool operator <= (double num) { //if num == value ,but the "except" is true, we still return false
		return (value < num) || (value == num && (!except));
	}
	bool operator >= (double num) { //if num == value ,but the "except" is true, we still return false
		return (value > num) || (value == num && (!except));
	}
	bool operator == (double num) { //if num == value ,but the "except" is true, we still return false
		return (value == num) && (!except);
	}
	//static function
	static bool myxor(bool p, bool q) {
		return (p&&q) || (!(p || q));
	}
};
class Interval {
protected:
	IntervalPoint start;
	IntervalPoint end;
	bool reverse;
public:
	Interval() :reverse(false){};
	Interval(double _min, double _max, bool ept_min = false, bool ept_max = false, bool _rvs = false) :
		start(_min, ept_min), end(_max, ept_max), reverse(_rvs){};
	~Interval(){};
	virtual void set(double _min, double _max, bool ept_min = false, bool ept_max = false, bool _rvs = false) {
		this->start.set(_min, ept_min);
		this->end.set(_max, ept_max);
		this->reverse = _rvs;
	}
	virtual void parse(std::string str) {
		//str format！！"[start][end][ept_s][ept_e][rvs]"
		//1--true, 0--false
		//seprator！！" "(space)
		std::stringstream ss(str);
		Interval::parse(ss);

	}
	virtual void parse(std::stringstream&ss) {
		std::string par_str;
		double _min = 0, _max = 0;
		bool ept_min = false, ept_max = false, _rvs = false;
		//par1:start
		if (ss >> par_str)
			_min = File::str2num<double>(par_str);
		//par2:end
		if (ss >> par_str)
			_max = File::str2num<double>(par_str);
		//par3:ept_s
		if (ss >> par_str)
			ept_min = File::str2num<bool>(par_str);
		//par4:ept_e
		if (ss >> par_str)
			ept_max = File::str2num<bool>(par_str);
		//par5:rvs
		if (ss >> par_str)
			_rvs = File::str2num<bool>(par_str);
		//set
		this->set(_min, _max, ept_min, ept_max, _rvs);
	}
	double get_min()const{ return start.get_value(); }
	double get_max()const{ return end.get_value(); }
	virtual bool is_include(double num) {
		return IntervalPoint::myxor((start <= num && end >= num), (!reverse));
	}
};
class Index :public Interval{
protected:
	bool enable;
public:
	Index() :enable(false){};
	Index(double _min, double _max, bool ept_min = false, bool ept_max = false, bool rvs = false,
		bool eb = true) :enable(eb), Interval(_min, _max, ept_min, ept_max, rvs){};
	void parse(std::string str) {
		//str format！！"[start][end][ept_s][ept_e][rvs][enable]"
		//1--true, 0--false
		//seprator！！" "(space)
		std::stringstream ss(str);
		Index::parse(ss);
	}
	void parse(std::stringstream&ss) {
		std::string enable_str;
		bool _eb = true;
		Interval::parse(ss);
		if (ss >> enable_str)
			_eb = File::str2num<bool>(enable_str);
		enable = _eb;
	}
	void set(double _min, double _max, bool ept_min = false, bool ept_max = false, bool rvs = false,
		bool eb = true) {
		this->enable = eb;
		Interval::set(_min, _max, ept_min, ept_max, rvs);
	}
	void set_enable(bool eb){ this->enable = eb; }
	void open_enable() { this->enable = true; }
	void close_enable() { this->enable = false; }
	bool get_enable()const { return this->enable; }
	bool is_include(double num) {
		return (!enable) || Interval::is_include(num);
	}
	~Index(){};
};
class Grid :public Interval {
protected:
	double step;
public:
	Grid() :step(0){};
	Grid(double _min, double _step, double _max, bool _epts = false, bool _epte = false, bool _rvs = false)
		:step(_step), Interval(_min, _max, _epts, _epte, _rvs){};
	~Grid(){};
	void parse(std::stringstream&ss) {
		std::string step_str;
		double _step = 1.0;
		Interval::parse(ss);
		if (ss >> step_str)
			_step = File::str2num<double>(step_str);
		this->step = _step;

	}
	void parse(std::string str) {
		//str format！！"[start][end][ept_s][ept_e][rvs][step]"
		//1--true, 0--false
		//seprator！！" "(space)
		std::stringstream ss(str);
		Grid::parse(ss);
	}
	void set(double _min, double _max, double _step, bool ept_min = false, bool ept_max = false, bool _rvs = false) {
		this->step = _step;
		Interval::set(_min, _max, ept_min, ept_max, _rvs);
	}
	double get_begin() const{
		int size = get_size();
		if (size <= 0)
			return std::numeric_limits<double>::quiet_NaN();
		return  start.get_state() ? (start.get_value() + step) : (start.get_value());
	}
	int get_size()const {
		if (this->reverse)
			return -1;
		bool ept_min = start.get_state(), ept_max = end.get_state();
		int size = int((end.get_value() - start.get_value()) / step + 1);
		int count_min = ept_min ? (-1) : (0);
		bool is_eq = abs(start.get_value() + step*(size - 1) - end.get_value()) <= std::numeric_limits<double>::epsilon();
		int count_max = ept_max && is_eq ? (-1) : (0);
		return size + count_min + count_max;
	}
	double operator[] (int i) const{ //from zero to start
		double begin = Grid::get_begin();
		if (i < 0 || i >= Grid::get_size()) {
			return std::numeric_limits<double>::quiet_NaN();
		}
		return begin + step*i;
	}
};
class IonGridSample {
protected:
	double t;
	double lat;
	double lon;
public:
	IonGridSample() :t(0), lat(0), lon(0) {};
	IonGridSample(double _t, double _lat, double _lon) :t(_t), lat(_lat), lon(_lon){};
	~IonGridSample(){};
	void set(double _t, double _lat, double _lon) {
		this->t = _t;
		this->lat = _lat;
		this->lon = _lon;
	}
	void parse(std::string str) {
		using namespace std;
		stringstream ss(str);
		double _t = 0, _lat = 0, _lon = 0;
		ss >> _t >> _lat >> _lon;
		this->set(_t, _lat, _lon);
	
	}
	std::string toString() { std::ostringstream os; os << *this; return os.str(); }
	//Declare a friend
	friend IonPoint;
	friend IonPrf;
	friend std::ostream& operator<< (std::ostream&os, const IonGridSample&igs) {
		using namespace std;
		os << setiosflags(ios::fixed)
			<< setprecision(1) << setfill('0')
			<< setw(4) << igs.t << " "
			<< setw(4) << igs.lat << " "
			<< setw(4) << igs.lon << " "
			<< resetiosflags(ios::fixed);
		return os;
	}
};
class IonGrid {
protected:
	Grid win_t;
	Grid win_lat;
	Grid win_lon;
	bool is_bind_bl;
public:
	IonGrid() { Initial(); };
	~IonGrid(){};
	void Initial() {
		win_t.set(0, 15, 1);
		win_lat.set(0, 10, 1);
		win_lon.set(0, 10, 1);
		is_bind_bl = false;
	}
	void parse(std::string str) {
		using namespace std;
		vector<string> strs;
		File::split(str, string(";"), strs);
		vector<string>::iterator it(strs.begin());
		if (it != strs.end() && !it->empty()) {
			win_t.parse(*it++);
		}
		if (it != strs.end() && !it->empty()) {
			win_lat.parse(*it++); is_bind_bl = true;
		}
		if (it != strs.end() && !it->empty()) {
			win_lon.parse(*it++); is_bind_bl = false;
		}
	}
	void SampleGrid(std::vector<IonGridSample>&dst) {
		int t_i = 0, lat_i = 0, lon_i = 0;
		int t_size = win_t.get_size(), lat_size = win_lat.get_size(), lon_size = win_lon.get_size();
		double t_v = 0, lat_v = 0, lon_v = 0;
		if (is_bind_bl) {
			for (t_i = 0; t_i < t_size; t_i++) {
				for (lat_i = 0; lat_i < lat_size; lat_i++) {
					t_v = win_t[t_i]; lat_v = win_lat[lat_i]; lon_v = lat_v;
					dst.push_back(IonGridSample(t_v, lat_v, lon_v));
				}
			}
		}
		else {
			for (t_i = 0; t_i < t_size; t_i++) {
				for (lat_i = 0; lat_i < lat_size; lat_i++) {
					for (lon_i = 0; lon_i < lon_size; lon_i++) {
						t_v = win_t[t_i]; lat_v = win_lat[lat_i]; lon_v = win_lon[lon_i];
						dst.push_back(IonGridSample(t_v, lat_v, lon_v));
					}
				}
			}
		}

	}
	//************************************
	//Declare a friend
	friend IonPoint;
	friend IonPrf;
};
class IonMatch;
class IonBias{
protected:
	double AbMean;
	double AbVar;
	double ReMean;
	double ReVar;
public:
	IonBias() :AbMean(0), AbVar(0), ReMean(0), ReVar(0){};
	~IonBias(){};
	IonBias(double _am, double _av, double _rm, double _rv) :AbMean(_am), AbVar(_av), ReMean(_rm), ReVar(_rv){};
	IonBias(const IonBias&ib) :AbMean(ib.AbMean), AbVar(ib.AbVar), ReMean(ib.ReMean), ReVar(ib.ReVar){};
	void set(double _am, double _av, double _rm, double _rv) {
		this->AbMean = _am; this->AbVar = _av;
		this->ReMean = _rm; this->ReVar = _rv;
	}
	void set(const IonBias&ib) {
		this->AbMean = ib.AbMean; this->AbVar = ib.AbVar;
		this->ReMean = ib.ReMean; this->ReVar = ib.ReVar;
	}
	std::string toString(){ std::ostringstream os; os << *this; return os.str(); }
	//************************************
	//Declare a friend
	friend IonMatch;
	friend std::ostream& operator<< (std::ostream&os, const IonBias&ib) {
		using namespace std;
		os << setiosflags(ios::fixed | ios::internal | ios::showpos)
			<< setprecision(3) << setfill('0')
			<< setw(6) << ib.AbMean << "  "
			<< setw(6) << ib.AbVar << "  "
			<< setw(6) << ib.ReMean << "  "
			<< setw(6) << ib.ReVar << "  "
			<< resetiosflags(ios::fixed | ios::internal | ios::showpos);
		return os;
	}
};
class IonMatch {
protected:
	double R_nmF2;
	double R_hmF2;
	IonBias B_nmF2;
	IonBias B_hmF2;
	int N;
public:
	IonMatch() :R_nmF2(0), R_hmF2(0), N(0){};
	~IonMatch(){};
	IonMatch(double _rn, double _rh, int _n,
		     double _amn, double _avn, double _rmn, double _rvn,
			 double _amh, double _avh, double _rmh, double _rvh) 
			 :R_nmF2(_rn), R_hmF2(_rh), N(_n), B_nmF2(_amn, _avn, _rmn, _rvn), B_hmF2(_amh, _avh, _rmh, _rvh){};
	IonMatch(double _rn, double _rh, int _n, const IonBias&ibn, const IonBias&ibh)
			 :R_nmF2(_rn), R_hmF2(_rh), N(_n), B_nmF2(ibn), B_hmF2(ibh){};
	void set(double _rn, double _rh, int _n,
		double _amn, double _avn, double _rmn, double _rvn,
		double _amh, double _avh, double _rmh, double _rvh) {
		this->R_nmF2 = _rn;
		this->R_hmF2 = _rh;
		this->N = _n;
		this->B_nmF2.set(_amn, _avn, _rmn, _rvn);
		this->B_hmF2.set(_amh, _avh, _rmh, _rvh);
	}
	void set(double _rn, double _rh, int _n, const IonBias&ibn, const IonBias&ibh) {
		this->R_nmF2 = _rn;
		this->R_hmF2 = _rh;
		this->N = _n;
		this->B_nmF2.set(ibn);
		this->B_hmF2.set(ibh);
	}
	std::string toString(){ std::ostringstream os; os << *this; return os.str(); }
	//Declare a friend
	friend IonPoint;
	friend IonPrf;
	friend std::ostream& operator<< (std::ostream&os, const IonMatch&im) {
		using namespace std;
		os << setiosflags(ios::fixed | ios::internal | ios::showpos)
			<< setprecision(0) << setfill('0')
			<< setw(5) << im.N << "  "
			<< setprecision(3)
			<< setw(6) << im.R_nmF2 << "  "
			<< setw(6) << im.R_hmF2 << "  "
			<< im.B_nmF2 << " "
			<< im.B_hmF2 << " "
			<< resetiosflags(ios::fixed | ios::internal | ios::showpos);
		return os;
	}
};
class MatchIndex {
protected:
	Index qcdN;
	Index qcdH;
	Index qcdA;
	Index qcA;
	Index qcRdN;
	Index qcRdH;
	bool is_abs;
	bool is_ept_svy;//a1
	bool is_ept_ref;//a2
public:
	MatchIndex():is_abs(false), is_ept_svy(false), is_ept_ref(false){};
	void set_qcdN(Index p){this->qcdN = p;}
	void set_qcdH(Index p){this->qcdH = p;}
	void set_qcdA(Index p){this->qcdA = p;}
	void set_qcA(Index p, bool _svy, bool _ref){
		this->qcA = p; 
		this->is_ept_svy = _svy; this->is_ept_ref = _ref;
	}
	void set_qcdN(std::string&str) {
		double maxdiffN = File::str2num<double>(str);
		this->qcdN.set(-maxdiffN, maxdiffN);
	}
	void set_qcdH(std::string&str) {
		double maxdiffH = File::str2num<double>(str);
		this->qcdH.set(-maxdiffH, maxdiffH);
	}
	void set_qcdA(std::string&str) {
		double maxdiffA = File::str2num<double>(str);
		this->qcdA.set(-maxdiffA, maxdiffA);
	}
	void set_qcA(std::string&str) {
		using namespace std;
		vector<string> pars;
		File::split(str, string(";"), pars);
		int num = pars.size();
		if(num > 0) {
			this->qcA.parse(pars[0]);
		}
		if(num > 1) {
			this->is_ept_svy=(File::str2num<int>(pars[1])!=0);
		}
		if(num > 2) {
			this->is_ept_ref=(File::str2num<int>(pars[2])!=0);
		}
	}
	void set_qcRdN(std::string&str) {
		double maxdiffN = File::str2num<double>(str);
		this->qcRdN.set(-maxdiffN, maxdiffN);
	}
	void set_qcRdH(std::string&str) {
		double maxdiffH = File::str2num<double>(str);
		this->qcRdH.set(-maxdiffH, maxdiffH);
	}
	void set_is_abs(bool _eb){this->is_abs = _eb;}
	bool is_include(double _dn, double _dh, double _a1, double _a2, double _n2, double _h2) {
		bool flag1 = qcdN.is_include(_dn);
		bool flag2 = qcdH.is_include(_dh);
		double _A1 = (is_abs)?(_a1<0?(_a1+180):(_a1)):(_a1);
		double _A2 = (is_abs)?(_a2<0?(_a2+180):(_a2)):(_a2);
		double _dA = _A1 - _A2;
		bool flag3 = qcdA.is_include(_dA)
					&& (!qcdA.get_enable() || _A1 < 360)
					&& (!qcdA.get_enable() || _A2 < 360);
		bool flag4 = is_ept_svy || qcA.is_include(_A1);
		bool flag5 = is_ept_ref || qcA.is_include(_A2);
		bool flag6 = qcRdN.is_include(_dn/_n2*100);//%
		bool flag7 = qcRdH.is_include(_dh/_h2*100);//%
		return flag1 && flag2 && flag3 && flag4 && flag5 && flag6 && flag7;
	}
	friend IonPoint;
};
//****Futher class
class Coord {
protected:
	double Lat;//latitude--degree
	double Lon;//lontitude--degree
public:
	Coord() :Lat(0), Lon(0){};
	Coord(double _lat, double _lon) :Lat(_lat), Lon(_lon){};
	~Coord(){};
	void set(double _lat, double _lon) { Lat = _lat; Lon = _lon; }
	std::string toString(){ std::ostringstream os; os << *this; return os.str(); };
	bool is_appequal(const Coord&p, double alw_lat, double alw_lon) {
		return Coord::is_appequal(*this, p, alw_lat, alw_lon);
	}
	bool is_equal(const Coord&p) {
		double eps = std::numeric_limits<double>::epsilon();
		return Coord::is_appequal(*this, p, eps, eps);
	}
	//Declare frined
	friend IonPoint;
	friend std::ostream& operator<<(std::ostream&os, const Coord&cor)
	{
		using namespace std;//Only in this function
		os << setiosflags(ios::fixed | ios::internal | ios::showpos)
			<< setprecision(3) << setfill('0')
			<< setw(7) << cor.Lat << " "
			<< setw(8) << cor.Lon << " "
			<< resetiosflags(ios::fixed | ios::internal | ios::showpos);

		return os;
	}
	//static function
	static bool is_appequal(const Coord&p1, const Coord&p2, double alw_lat, double alw_lon) {
		return (abs(p1.Lat - p2.Lat) <= alw_lat) &&
			(abs(p1.Lon - p2.Lon) <= alw_lon);
	}
	static bool is_equal(const Coord&p1, const Coord&p2) {
		double eps = std::numeric_limits<double>::epsilon();
		return is_appequal(p1, p2, eps, eps);
	}

};
class Time {
protected:
	int Year, Doy, Hour, Min;
public:
	static const int DOY_FDM[13];//DOY of First Day in a Month (nonleap year)
	Time() :Year(0), Doy(0), Hour(0), Min(0){};
	Time(int _y, int _d, int _h, int _m) :Year(_y), Doy(_d), Hour(_h), Min(_m){};
	~Time(){};
	void set(int _y, int _d, int _h, int _m) { Year = _y; Doy = _d; Hour = _h; Min = _m; };
	void ReadTime(std::string&dataline) { //the delimiter is [space]
		std::stringstream ss(dataline);
		ss >> Year >> Doy >> Hour >> Min;
		ss.str(""); ss.clear();
	}
	void ReadTime(std::string&dataline, std::string delimiter) {
		using namespace std;
		vector<string> strs;
		vector<int> nums;
		stringstream ss;
		File::split(dataline, delimiter, strs);
		File::str2num<int>(strs, nums);
		this->Year = nums[0]; this->Doy = nums[1];
		this->Hour = nums[2]; this->Min = nums[3];
	}
	std::string toString(){ std::ostringstream os; os << *this; return os.str(); };
	double get_JD()const {
		int y = this->Year, m = 1, d = 1;
		double h = this->Hour + this->Min / 60.0;
		y = y - 1; m = m + 12;//when m<=2,we must deal some details
		double JD = int(365.25*y) + int(30.6001*(m + 1)) + d + h / 24 + 1720981.5;
		JD += this->Doy - 1;

		return JD;
	};
	double get_MJD()const{ return this->get_JD() - 2400000.5; };
	//Declare frined
	friend IonPoint;
	friend std::ostream &operator<<(std::ostream&os, const Time&t)
	{
		using namespace std;//Only in this function
		os << setfill('0')
			<< setw(4) << t.Year << " "
			<< setw(3) << t.Doy << " "
			<< setw(2) << t.Hour << " "
			<< setw(2) << t.Min << " ";
		return os;
	}
	//static function
	static bool is_leap(int year) { return (year%100 == 0)?(year%400 == 0):(year%4 == 0);}
	static bool is_appequal(const Time&t1, const Time&t2, double alw) {
		return abs(cal_interval(t1, t2)) <= alw;

	}
	static bool is_equal(const Time&t1, const Time&t2) {
		return t1.Year == t2.Year&&t1.Doy == t2.Doy&&t1.Hour == t2.Hour&&t1.Min == t2.Min;
	}
	static bool is_appless(const Time&t1, const Time&t2, double alw) {
		return cal_interval(t1, t2) < -alw;
	};
	static bool is_appgreater(const Time&t1, const Time&t2, double alw) {
		return cal_interval(t1, t2) > alw;
	};
	static double cal_interval(const Time&t1, const Time&t2) {
		return (t1.get_MJD() - t2.get_MJD()) * 1440;
	};//interval units: min;
	static int cal_doy(int year, int mon, int day) {
		int flag = (Time::is_leap(year) && mon >= 3)?(1):(0);
		int doy = flag + Time::DOY_FDM[mon-1] + (day-1);
		return doy;
	};

};
class ICP {
protected:
	double NmF2;//max electron density in F2--x10^5el/cm^-3
	double hmF2;//corresponding height of NmF2--Km
	double AmF2;//corresponding occultation plane azimuth of NmF2--degree, [-80,180)
	std::string profile;
public:
	static const std::string epfFlag;
	ICP() :NmF2(0), hmF2(0), AmF2(999),profile(epfFlag){};//999 represents unknown value
	ICP(double _nmf2, double _hmF2, double _amF2)
		:NmF2(_nmf2), hmF2(_hmF2), AmF2(_amF2), profile(epfFlag){};
	ICP(double _nmf2, double _hmF2, double _amF2, std::string _profile)
		:NmF2(_nmf2), hmF2(_hmF2), AmF2(_amF2), profile(_profile){};
	////if there is [f0f2], the best way is use <foF2_to_NmF2> at first in <ICP>construct function
	~ICP(){};
	const double get_NmF2()const { return NmF2; };
	const double get_hmF2()const { return hmF2; };
	const double get_AmF2()const { return AmF2; };
	const std::string get_profile()const { return profile;}
	void set(double _nmf2, double _hmf2, double _amF2) { 
		NmF2 = _nmf2; hmF2 = _hmf2; AmF2 = _amF2;
	};
	void set(double _nmf2, double _hmf2, double _amF2, std::string _profile) { 
		NmF2 = _nmf2; hmF2 = _hmf2; AmF2 = _amF2; profile = _profile;
	}
	std::string toString(){ std::ostringstream os; os << *this; return os.str(); };
	std::string toSimpleString() {
		using namespace std;
		ostringstream os;
		os << setiosflags(ios::fixed | ios::internal | ios::showpos) 
			<< setprecision(3) << setfill('0')
			<< setw(8) << this->NmF2 << " "
			<< setw(8) << this->hmF2 << " "
			<< setw(8) << this->AmF2 << " "
			<< resetiosflags(ios::fixed | ios::internal | ios::showpos);
		return os.str();
	}
	void f0F2_to_NmF2(){ this->NmF2 = ICP::f0F2_to_NmF2(this->NmF2); };
	bool hasProfile() { return this->profile != epfFlag;}
	//Declare frined
	friend IonPoint;
	friend std::ostream& operator<<(std::ostream&os, const ICP&icp)
	{
		using namespace std;
		os << setiosflags(ios::fixed | ios::internal | ios::showpos) 
			<< setprecision(3) << setfill('0')
			<< setw(8) << icp.NmF2 << " "
			<< setw(8) << icp.hmF2 << " "
			<< setw(8) << icp.AmF2 << " "
			<< resetiosflags(ios::fixed | ios::internal | ios::showpos)
			<< setiosflags(ios::left) << setfill(' ')
			<< setw(35)<< icp.profile<<" "
			<< resetiosflags(ios::left);
		return os;
	}
	//static function
	static double f0F2_to_NmF2(double _f0F2){ return 1.24*1e-1*_f0F2*_f0F2; }
};

//****Final class
//DigStation
class DigStation :public Coord {
protected:
	std::string Name;
public:
	DigStation() :Name(""){};
	DigStation(std::string _name, double _lat, double _lon) :Name(_name), Coord(_lat, _lon){};
	~DigStation(){};
	const std::string&get_name()const{ return Name; }
};
class DigStationList {
protected:
	DigStationList(){};
	~DigStationList(){};//Stop instantiation
public:
	static const std::vector<DigStation> stas;
	static bool search_sta(std::string&name,DigStation&dst) {
		std::vector<DigStation>::const_iterator  it_stas = stas.begin();
		////we can also use "auto" when c++11 or newer
		while (it_stas != stas.end()) {
			if (it_stas->get_name() == name) {
				dst = *it_stas;
				return true;
			}
			it_stas++;
		}
		return false;
	}
};
//IonPoint
class IonPoint {
protected:
	Coord loc;
	Time t;
	ICP par;
	//Inner function
	static double op_square_diff(double _src1, double _src2){ return (_src1 - _src2)*(_src1 - _src2); };
	static double op_relat_diff(double _base, double _src){ return abs((_src - _base) / _base); };
public:
	IonPoint(){};
	IonPoint(std::string&dataline){ ReadIonPoint(dataline); };
	IonPoint(const IonPoint&ip) {this->set(ip);}
	~IonPoint(){};
	void set(const IonPoint&ip) {
		this->loc.set(ip.loc.Lat, ip.loc.Lon);
		this->t.set(ip.t.Year, ip.t.Doy, ip.t.Hour, ip.t.Min);
		this->par.set(ip.par.NmF2, ip.par.hmF2, ip.par.AmF2, ip.par.profile);
	}
	void ReadIonPoint(std::string&dataline) {
		std::stringstream isos(dataline);
		isos >> t.Year >> t.Doy >> t.Hour >> t.Min
			 >> loc.Lat >> loc.Lon
			 >> par.NmF2 >> par.hmF2 >> par.AmF2;
		isos >> par.profile;
		par.profile = par.profile.empty()?(ICP::epfFlag):(par.profile);
	};  //get information from the string
	std::string toString(){ std::stringstream os; os << *this; return os.str(); };
	//Declare the friend
	friend IonPrf;
	friend IonGroup;
	friend std::ostream&operator<<(std::ostream&os, const IonPoint&ip)
	{
		os << ip.t << " "
			<< ip.loc << " "
			<< ip.par << " ";
		return os;
	}
	//static function
	static void ConvertByDig(std::string file_in, std::string file_out, int qc_core, bool isParN, bool isParH);
	static bool is_valid_dig(std::string str) {
		return File::is_num(str) && str.find("---") == std::string::npos;
	}
	static void ConvertByBool(std::vector<double>&src, std::vector<double>&dst, std::vector<bool>&_bool);
	static bool less_t(const IonPoint&ip1, const IonPoint&ip2) {
		return ip1.t.get_JD() < ip2.t.get_JD();
	}
	static void ReadIonPointFile(std::string file, std::vector<IonPoint>&ips) {
		//Open 
		std::ifstream in(file);
		IonPoint ip;
		std::string str;
		//Construct 
		while (std::getline(in, str)) {
			ip.ReadIonPoint(str);//get the IonPoint
			ips.push_back(ip);
		}
		//Close 
		in.close();
		//std::copy(ips.begin(), ips.end(), std::ostream_iterator<IonPoint>(std::cout,"\n"));
		//std::cout<<"****************************"<<std::endl;
	}
	static void WriteMatchICP(std::ostream&out, std::vector<ICP>&svy, std::vector<ICP>&ref) {
		using namespace std;
		//The head str
		//string head_str = "svy_nmF2\t svy_hmF2\t svy_profile\t ref_nmF2\t ref_hmF2\t ref_profile\t";
		//out << head_str << endl;
		//The data
		vector<ICP>::iterator it_svy(svy.begin()), it_ref(ref.begin());
		for (; it_svy != svy.end(), it_ref != ref.end(); it_svy++, it_ref++) {
			out << it_svy->toString() << " " << it_ref->toString() << endl;
		}
	}
	static void WriteMatchCorrCoef(std::ostream&out, IonGrid&alws, std::vector<IonMatch>&ims) {
		using namespace std;
		//The head str
		//string head_str = "Win_t\t Win_lat\t Win_lon\t Number\t R_nmF2\t R_hmF2\t [am av rm rv][nmF2 hmF2]\t";
		//out << head_str << endl;
		//The data
		vector<IonGridSample> igs;
		alws.SampleGrid(igs);

		vector<IonGridSample>::iterator it_igs(igs.begin());
		vector<IonMatch>::iterator it_ims(ims.begin());
		for (; it_igs != igs.end(), it_ims != ims.end(); it_igs++, it_ims++) {
			out << it_igs->toString() << " " << it_ims->toString() << endl;
		}
	
	}
	static void WriteMatchCorrCoef(std::ostream&out, IonGridSample&alw, IonMatch&im) {
		using namespace std;
		//The head str
		//string head_str = "Win_t\t Win_lat\t Win_lon\t Number\t R_nmF2\t R_hmF2\t [am av rm rv][nmF2 hmF2]\t";
		//out << head_str << endl;
		//The data
		out << alw << " " << im << endl;

	}
	static int BiSearchInsideLoc(std::vector<IonPoint>&data, IonPoint&p, IonGridSample&alw);
	static int BiSearch(std::vector<IonPoint>&data, IonPoint&p, int startloc, IonGridSample&alw, std::vector<IonGroup>&dst,
		MatchIndex&info);
	static void Match(std::vector<IonPoint>&svy, std::vector<IonPoint>&ref, IonGridSample&alw, std::vector<IonGroup>&dst, 
		MatchIndex&info);
	static void Match(std::string file_svy, std::string file_ref, IonGridSample&alw, std::vector<IonGroup>&dst, 
		MatchIndex&info);
	
	static void Match(std::vector<IonPoint>&svy, std::vector<IonPoint>&ref, IonGridSample&alw, IonMatch&dst,
		MatchIndex&info, STRATEGY _stra = AVERAGE);
	static void Match(std::vector<IonPoint>&svy, std::vector<IonPoint>&ref, IonGrid&alws, std::vector<IonMatch>&dst,
		MatchIndex&info, STRATEGY _stra = AVERAGE);
	static void Match(std::string file_svy, std::string file_ref, IonGridSample&alw, IonMatch&dst,
		MatchIndex&info, STRATEGY _stra = AVERAGE);
	static void Match(std::string file_svy, std::string file_ref, IonGrid&alws, std::vector<IonMatch>&dst,
		MatchIndex&info, STRATEGY _stra = AVERAGE);
	
	static void MatchPreprocess(std::string file_svy, std::string file_ref, std::vector<IonPoint>&svys, std::vector<IonPoint>&refs);
	static void MoveAverage(std::vector<double>&src, std::vector<double>&dst, int _win = 9);
	//**y=kx+b
	static void LineRegression(std::vector<double>&src_x, 
				std::vector<double>&src_y, double&dst_k, double&dst_b);
	//**coreltion coeficient
	static void CorrCoefficient(std::vector<double>&src_x, std::vector<double>&src_y, double&dst_r);
	static void Cal_Bias(std::vector<double>&src_x, std::vector<double>&src_y, double&dst_am, double&dst_av, double&dst_rm, double&dst_rv);//x-svy, y-ref
	static void Cal_Bias(std::vector<double>&src_x, std::vector<double>&src_y, IonBias&ib);
	static void CorrCoefficient(std::vector<IonGroup>&src, IonMatch&dst_r, STRATEGY _stra = AVERAGE);
    

};
class IonGroup {
protected:
    IonPoint ref;
	std::vector<IonPoint> cmp;
public:
    IonGroup(){};
	~IonGroup(){};
	void set_ref(IonPoint&ip) {ref.set(ip);}
	void set_cmp(IonPoint&ip) {cmp.push_back(IonPoint(ip));}
	int get_cmpSize() {return this->cmp.size();}
	ICP get_cmp_ave() {
		using namespace std;
		double NmF2 = 0, hmF2 = 0, AmF2 = 0;
		ICP dst;
		vector<IonPoint>::iterator it = cmp.begin();
		while(it != cmp.end()) {
			NmF2 += it->par.get_NmF2();
			hmF2 += it->par.get_hmF2();
			AmF2 += it->par.get_AmF2();
			it++;
		}
		NmF2 /= cmp.size(); hmF2 /= cmp.size();  AmF2 /= cmp.size();
		dst.set(NmF2, hmF2, AmF2);
		return dst;
	}
	void get_cmp_eqa(std::vector<ICP>&_cmp, std::vector<ICP>&_ref) {
		using namespace std;
		vector<IonPoint>::iterator it = cmp.begin();
		while(it != cmp.end()) {
			_ref.push_back(this->ref.par);
			_cmp.push_back(it->par);
			it++;
		}
	}

	static void toICP(std::vector<IonGroup>&_src, std::vector<ICP>&_cmp, std::vector<ICP>&_ref, STRATEGY _stra = AVERAGE) {
		using namespace std;
		switch(_stra) {
			case AVERAGE: toICP_Ave(_src, _cmp, _ref);break;
			case EQUALITY: toICP_Eqa(_src, _cmp, _ref);break;
			default: cout<<"Error Strategy Type"<<endl;return;
		}
	}
	static void toICP_Ave(std::vector<IonGroup>&_src, std::vector<ICP>&_cmp, std::vector<ICP>&_ref) {
		using namespace std;
		vector<IonGroup>::iterator it = _src.begin();
		while(it != _src.end()) {
			_ref.push_back(it->ref.par);
			_cmp.push_back(it->get_cmp_ave());
			it++;
		}
	}
	static void toICP_Eqa(std::vector<IonGroup>&_src, std::vector<ICP>&_cmp, std::vector<ICP>&_ref) {
		using namespace std;
		vector<IonGroup>::iterator it = _src.begin();
		while(it != _src.end()) {
			it->get_cmp_eqa(_cmp , _ref);
			it++;
		}
	}
	friend class IonPoint;
};
class IonIndex {
protected:
	Index md;
	Index sigma;
	Index gk;
	Index lk;
	Index parN;
	Index parh;
	Index parA;
	int winsize;
	Interval localregress;
	bool is_accept_invalid;
	bool is_accept_qcfail;
public:
	IonIndex(){ Initial(); };
	~IonIndex(){};
	void Initial() {
		md.set(0, 0.1);//[0,0.1]
		sigma.set(0, 0.05);//[0, 0.05]
		gk.set(-1*std::numeric_limits<double>::infinity(), 0, false, true);//<0
		lk.set(-1*std::numeric_limits<double>::infinity(), 0, false, true);//<0
		parN.set(0, 20);//[0,20]x10^5el/cm^3
		parh.set(200, 500);//[200,500]km
		parA.set(-180,1000);//[-180,180]degree, and also accept 999
		winsize = 9;
		localregress.set(420, 490);//[420, 490]km
		is_accept_invalid = true;
		is_accept_qcfail = true;
	}
	std::string get_switchstate() {
		std::string switch_str;
		if (md.get_enable())
			switch_str += "MD\t";
		if (sigma.get_enable())
			switch_str += "SIGMA\t";
		if (gk.get_enable())
			switch_str += "GK\t";
		if (lk.get_enable())
			switch_str += "LK\t";
		if (parN.get_enable())
			switch_str += "parN\t";
		if (parh.get_enable())
			switch_str += "parh\t";
		if (parA.get_enable())
			switch_str += "parA\t";
		if (switch_str.empty())
			switch_str = "NULL\t";
		return switch_str;
	}
	void set_switch(bool _md, bool _sigma, 
					bool _gk, bool _lk,
					bool _parN = true, bool _parh = true, bool _parA = true,
					bool _iv = true, bool _qf = true) {
		this->md.set_enable(_md); this->sigma.set_enable(_sigma);
		this->gk.set_enable(_gk); this->lk.set_enable(_lk);
		this->parN.set_enable(_parN); this->parh.set_enable(_parh); this->parA.set_enable(_parA);
		this->is_accept_invalid = _iv; this->is_accept_qcfail = _qf;
	}
	void set_switch_all(bool _eb) {
		set_switch(_eb, _eb, _eb, _eb, _eb, _eb, _eb, _eb, _eb);
	}
	void set_switch_md(bool _md) { this->md.set_enable(_md); }
	void set_switch_sigma(bool _sigma) { this->sigma.set_enable(_sigma); }
	void set_switch_gk(bool _gk) { this->gk.set_enable(_gk); }
	void set_switch_lk(bool _lk) { this->lk.set_enable(_lk); }
	void set_switch_parN(bool _parN) { this->parN.set_enable(_parN); }
	void set_switch_parh(bool _parh) { this->parh.set_enable(_parh); }
	void set_switch_parA(bool _parA) { this->parA.set_enable(_parA); }
	void set_switch_invalid(bool _iv) { this->is_accept_invalid = _iv; }
	void set_switch_qcfail(bool _qf) { this->is_accept_qcfail = _qf; }
	void set_par_win(int _win) { this->winsize = _win; };
	void set_par_lri(double _s, double _e, bool ept_s = false, bool ept_e = false, bool _rvs = false) { 
		this->localregress.set(_s, _e, ept_s, ept_e, _rvs); 
	}
	void set_par_lri(std::string&str) {
		this->localregress.parse(str);
	}
	void set_par_all(std::string&str) {
		using namespace std;
		stringstream ss(str);
		string md_str, lri_str;
		//get the md_str
		ss >> md_str;
		this->winsize = File::str2num<int>(md_str);
		//get the lri_str
		lri_str = str.substr(md_str.size(), str.size());
		this->localregress.parse(lri_str);
	}
	void set_par_qc(std::string&str) {
		using namespace std;
		vector<string> par_qc;
		File::split(str, string(";"), par_qc);
		vector<string>::iterator it = par_qc.begin();
		if (it != par_qc.end() && !it->empty() && md.get_enable()) {
			md.parse(*it);
		}
		it++;

		if (it != par_qc.end() && !it->empty() && sigma.get_enable()) {
			sigma.parse(*it);
		}
		it++;

		if (it != par_qc.end() && !it->empty() && gk.get_enable()) {
			gk.parse(*it);
		}
		it++;

		if (it != par_qc.end() && !it->empty() && lk.get_enable()) {
			lk.parse(*it);
		}
		it++;

		if (it != par_qc.end() && !it->empty() && parN.get_enable()) {
			parN.parse(*it);
		}
		it++;

		if (it != par_qc.end() && !it->empty() && parh.get_enable()) {
			parh.parse(*it);
		}
		it++;
		
		if (it != par_qc.end() && !it->empty() && parA.get_enable()) {
			parA.parse(*it);
		}
		it++;
	}
	bool is_include(std::vector<double>&src) {
		using namespace std;
		vector<double>::iterator it = src.begin();
		if (md.get_enable() && it != src.end() && !md.is_include(*it++)) { 
			//"++" is because only if we callback md.is_include, this index number is used, can not use again
			return false;//just find one unquality index, we can return false
		}
		if (sigma.get_enable() && it != src.end() && !sigma.is_include(*it++)) {
			return false;
		}
		if (gk.get_enable() && it != src.end() && !gk.is_include(*it++)) {
			return false;
		}
		if (lk.get_enable() && it != src.end() && !lk.is_include(*it++)) {
			return false;
		}
		if (parN.get_enable() && it != src.end() && !parN.is_include(*it++)) {
			return false;
		}
		if (parh.get_enable() && it != src.end() && !parh.is_include(*it++)) {
			return false;
		}
		if (parA.get_enable() && it != src.end() && !parA.is_include(*it++)) {
			return false;
		}
		return true;
	}
	//Declare a friend
	friend IonPrf;
};

struct AziIndex {
	std::string copath;
	std::string subdir;
	bool is_abs;
	AziIndex():copath("#"),subdir("*"),is_abs(false){};
	AziIndex(std::string&str) {this->set(str);}
	void set(std::string&str) {
		using namespace std;
		this->copath = "#"; this->subdir = "*"; this->is_abs = false;
		vector<string> pars;
		File::split(str, string(";"), pars);
		int num = pars.size();
		if(num > 0)
			this->copath = pars[0];
		if(num > 1)
			this->subdir = pars[1];
		if(num > 2)
			this->is_abs = (File::str2num<int>(pars[2]) != 0);
	};
};

//IonPrf
class IonPrf {
protected:
	std::vector<double> dens;//Density data
	std::vector<double> alt;//Altitude data
	IonPoint ip;//ionpoint data
	File datafile;//Ionprf file information
public:
	//Constuctor and destructor
	IonPrf(){};
	IonPrf(std::string _ionprf, AziIndex&msg){ 
		ReadIonPrf(_ionprf, msg); 
	};
	~IonPrf(){};
	bool ReadIonPrf(std::string _ionprf, AziIndex&msg);
	bool ReadIonPrfByIONPRF(std::string _ionprf, AziIndex&msg);
	bool ReadIonPrfByFY3C(std::string _ionprf, AziIndex&msg);
	//Write Report
	std::string WriteReport(IonIndex&qc);
	bool WriteReport(IonIndex&qc, std::ostream&out);
	//Write Info
	std::string generate_filename_txt(){ return datafile.get_address() + std::string(".txt"); }
	std::string generate_filename_txt(std::string src_dir) { return datafile.replace_filepath(src_dir) + std::string(".txt"); }
	void WriteInfo(std::string file_out);
	void WriteInfo() {
		std::string file_out = this->generate_filename_txt();
		IonPrf::WriteInfo(file_out);
	}
	//Status judgement
	bool is_empty(){ return datafile.is_empty(); };
	//Information mining
	void cal_fluctuation(double&dst_MD, double&dst_sigma, IonIndex&qc);
	void cal_regression(double&dst_gk, double&dst_lk, IonIndex&qc);
	//static function


};

class OccAzi {
protected:
	double AmF2;
	File datafile;//Ionprf file information
	//innner function
	std::string creIonPhs(std::string _copath, std::string _subdir);
	std::string creIonPhsByIONPRF(std::string _copath, std::string _subdir);
	std::string creIonPhsByFY3C(std::string _copath, std::string _subdir);

public:	
	//con-get-set
	//**con
	OccAzi():AmF2(999){};
	OccAzi(std::string _ionprf, double _hmF2, AziIndex&msg){ 
		this->set(_ionprf, _hmF2, msg);
	}
	//**set
	void set(std::string _ionprf, double _hmF2, AziIndex&msg);
	//**get
	const double getAmF2()const {return AmF2;}
	
	static double checkAmF2(double _amF2, bool is_abs){
		if(is_abs && _amF2 < 0)_amF2+=180;
		return _amF2;
	};
};




#endif//IONPRF_H_INCLUDE
///~