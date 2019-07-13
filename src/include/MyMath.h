//Matrix.h
//20190321  V0
#ifndef MyMath_H_INCLUDE
#define MyMath_H_INCLUDE

#include<iostream>
#include<numeric>
#include<functional>
#include<vector>
#include<iterator>
#include<algorithm>
#include<utility>
#include<cmath>
#include<sstream>
#include<string>
//****************************************************************
//**************const parameter****************
const double PI = 3.1415926;
const double EPS = 1e-10;
enum FRAMEAXIS {AXISX, AXISY, AXISZ};
//*************static function*****************
//str2num
template<class T>
static T str2num(std::string str) { std::stringstream ss; T num; ss << str; ss >> num; return num; }
template<class T>
static std::string num2str(T num) { std::stringstream ss; ss << num; return ss.str(); }

//double2int
static int INT(double value) { return (value>0)?(int(value + EPS)):(int(abs(value)+EPS)*-1);};
static double FRAC(double value) { return value - INT(value);}
static void INT_FRAC(double value, int&dst_int, double&dst_frac){ 
	dst_int = INT(value); dst_frac = value - dst_int;};
static int FLOOR(double value) { return floor(value+EPS);}

//deg-arc transform
static double DEG2ARC(double _deg) { return _deg/180.0*PI;};
static double ARC2DEG(double _arc) { return _arc/PI*180.0;};

//Trigonometric function
static double SIND(double _deg) { return sin(DEG2ARC(_deg));};
static double COSD(double _deg) { return cos(DEG2ARC(_deg));};
static double TAND(double _deg) { return tan(DEG2ARC(_deg));};
static double DASIN(double _arc) { return ARC2DEG(asin(_arc));};
static double DACOS(double _arc) { return ARC2DEG(acos(_arc));};
static double DATAN(double _arc) { return ARC2DEG(atan(_arc));};

//value compare
template<typename T> static T MAX(T _p1, T _p2){ return _p1<_p2?(_p2):(_p1);}
template<typename T> static T MIN(T _p1, T _p2){ return _p1<_p2?(_p1):(_p2);}
template<typename T> static T ABS(T _p){ return _p>0?(_p):(-_p);}
template<typename T> static bool APPR(T _a, T _b, T _ale){ return ABS<T>(_a-_b) <= _ale;}
template<typename T> static bool EQUAL(T _a, T _b){ return ABS<T>(_a-_b) <= EPS;}
template<typename T> static bool NEQUAL(T _a, T _b){ return ABS<T>(_a-_b) > EPS;}
template<typename T> static bool LEQUAL(T _a, T _b) { return EQUAL<T>(_a, _b) || _a < _b;}
template<typename T> static bool GEQUAL(T _a, T _b) { return EQUAL<T>(_a, _b) || _a > _b;}
template<typename T> static bool LESS(T _a, T _b) { return NEQUAL<T>(_a, _b) && _a < _b;}
template<typename T> static bool GREAT(T _a, T _b) { return NEQUAL<T>(_a, _b) && _a > _b;}
template<typename T> static int SIGN(T _value){ return (EQUAL<T>(_value, 0))?(0):(_value>0)?(1):(-1);}
template<typename T> static int SIGNPOS(T _value){ return (GEQUAL<T>(_value, 0))?(1):(-1);}

//interval with k-based	
static int KBasedInterval(int value, int k, int&dst_newvalue, bool check_pos = false) {
	//dst_newvalue: remain number
	//return: carry number
	//tips: (1)calculate complement fro value in k-based interval when check_pos is true
	//      (2)value, dst_newvalue can be the same variable
	int carry = (value<0 && value%k!=0)?(value/k-1):(value/k);
	//use "carry" is for the condition where dst_newvalue is value
	dst_newvalue = (check_pos || value>=0)?(value-k*carry):(value-k*(carry+1));
	return carry;
};
static int KBasedInterval(double value, int k, double&dst_newvalue, bool check_pos = false) {
	//dst_newvalue: remain number
	//return: carry number
	//tips: (1)calculate complement fro value in k-based interval when check_pos is true
	//      (2)value, dst_newvalue can be the same variable
	int carry = FLOOR(value/k);
	//use "carry" is for the condition where dst_newvalue is value
	dst_newvalue = (check_pos || !(value<0))?(value-k*carry):(value-k*(carry+1));
	return carry;
};
//large number calculation

//equation
template<typename FunObject>
static double BinaryEquation(FunObject object, double xs, double xe,
					double y = 0, double eps = 1e-3) {
	//FunObject: y=f(x), continous
	//xs: start of the interval
	//xe: end of the interval
	//y: value of fx
	//eps: allow error for fx and interval of (xe-xs)
	//return: x which meet f(x)=y
	//***********************************************
	if(xe<xs)return std::numeric_limits<double>::quiet_NaN();
	double fxs = object(xs), fxe = object(xe);
 	bool flag = (fxs-y)*(fxe-y)<0;
	if(abs(fxs-y)<eps)return xs;
	if(abs(fxe-y)<eps)return xe;
	if(!flag)return std::numeric_limits<double>::quiet_NaN();
	//**x in the (xs, xe)
	double x = 0, fx = 0;
	do {
		x = 0.5*(xs + xe); fx = object(x);//cout<<x<<"  "<<fx<<endl;
		flag = (fx-y)*(object(xe)-y)<0;
		xs = flag?(x):(xs); xe = flag?(xe):(x);
	}while(abs(fx-y)>eps && abs(xe-xs)>eps);
	return x;
}

template<typename FunObject>
static double IterEquation(FunObject object, 
				double x0 = 0, double eps = 1e-3) {
	//FunObject: y=f(x), continous
	//x0: initial value of x
	//eps: allow error for x
	//return: x which meet f(x)=x
	//***********************************************
	double x1 = x0;
	do {
		x0 = x1;
		x1 = object(x0);
	}while(abs(x0-x1)>eps);
	return x1;
}

//interged
static double AbelTransformObj(double xmin, double xmax, double ymin, double ymax, double x0) {
	//xmin, xmax: minimum and maximum of integral variable
	//ymin, ymax: minimum and maximum of f(xmin), f(xmax)
	//x0: const number in sqrt(x*x-a0*a0) partial formula
	//           xmax        f(x)
	// return:  S     ------------------- dx
	//           xmin   sqrt(x*x-a0*a0)
	//tips: (1)xmin and max must be very close
	//***************************************************
	double k, b, ds, ds1, ds2, dln = 0;
	//coefficent
	k = (ymax-ymin)/(xmax-xmin);//slope
	b = (xmax*ymin - xmin*ymax)/(xmax-xmin);//intercept
	//integration
	ds1 = sqrt((xmax+x0)*(xmax-x0));
	ds2 = sqrt((xmin+x0)*(xmin-x0));
	ds = ds1 - ds2;
	dln = log(abs((xmax + ds1)/(xmin + ds2)));
	//accumulate = sum(coefficent*integration)
	return k*ds + b*dln;
}

//triangle
static double TriHeight(double a, double b, double theta) {
	//a: distance of one line
	//b: distance of another line
	//theta: angle of line-a and line-b, [0, 180), degree
	//return the height of the opposite line
	//***********************************************
	double TanAlpha = a*SIND(theta)/(b-a*COSD(theta));
	return b*TanAlpha/sqrt(1+TanAlpha*TanAlpha);
}

static double TriBase(double a, double b, double theta) {
	//a: distance of one line
	//b: distance of another line
	//theta: angle of line-a and line-b, degree, (0, 180)
	//return the base of the opposite line
	//***********************************************
	return sqrt(a*a+b*b-2*a*b*COSD(theta));
}


//peak search
template<typename RandomAccessIter, typename T>
static RandomAccessIter SearchPeakLowerBound(RandomAccessIter begin, RandomAccessIter end, const T& threshold) {
	//begin: begin of the sequence, [begin, end)
	//end: end of the sequence, [begin, end)
	//threshold: threshold for peak search (suppress noise)
	//return: the lower boundary for peak range
	//tips: (1)return and return+1 can calculate the precise intersection point
	//***********************************************
	while(begin != end-1) {
		if(*begin<threshold && *(begin+1)>=threshold)
			return begin;
		begin++;
	}
	return end-1;
}
template<typename RandomAccessIter, typename T, typename Compare>
static RandomAccessIter SearchPeakLowerBound(RandomAccessIter begin, RandomAccessIter end, const T& threshold, Compare comp) {
	//begin: begin of the sequence, [begin, end)
	//end: end of the sequence, [begin, end)
	//threshold: threshold for peak search (suppress noise)
	//comp: define <
	//return: the lower boundary for peak range
	//tips: (1)return and return+1 can calculate the precise intersection point
	//***********************************************
	while(begin != end-1) {
		if(comp(*begin, threshold) && !comp(*(begin+1), threshold))
			return begin;
		begin++;
	}
	return end-1;
}

template<typename RandomAccessIter, typename T>
static RandomAccessIter SearchPeakUpperBound(RandomAccessIter begin, RandomAccessIter end, const T& threshold) {
	//begin: begin of the sequence, [begin, end)
	//end: end of the sequence, [begin, end)
	//threshold: threshold for peak search (suppress noise)
	//return: the upper boundary for peak range
	//tips: (1)return and return-1 can calculate the precise intersection point
	//***********************************************
	while(begin != end-1) {
		if(*begin>=threshold && *(begin+1)<threshold)
			return begin+1;
		begin++;
	}
	return end;
}
template<typename RandomAccessIter, typename T, typename Compare>
static RandomAccessIter SearchPeakUpperBound(RandomAccessIter begin, RandomAccessIter end, const T& threshold, Compare comp) {
	//begin: begin of the sequence, [begin, end)
	//end: end of the sequence, [begin, end)
	//threshold: threshold for peak search (suppress noise)
	//comp: define <
	//return: the upper boundary for peak range
	//tips: (1)return and return-1 can calculate the precise intersection point
	//***********************************************
	while(begin != end-1) {
		if(!comp(*begin, threshold) && comp(*(begin+1), threshold))
			return begin+1;
		begin++;
	}
	return end;
}

template<typename RandomAccessIter, typename T>
static std::pair<RandomAccessIter, RandomAccessIter> SearchPeakRange(RandomAccessIter begin, RandomAccessIter end,  const T& threshold) {
	//begin: begin of the sequence, [begin, end)
	//end: end of the sequence, [begin, end)
	//threshold: threshold for peak search (suppress noise)
	//return: the lower and upper boundary for peak range
	//tips: (1)return.first and return.first+1 can calculate the precise intersection point
	//		(2)return.second and return.second-1 can calculate the precise intersection point
	//      (3)max of [return.first, return.second) is the peak value
	//***********************************************
	RandomAccessIter it = SearchPeakLowerBound(begin, end, threshold);
	return make_pair(it, SearchPeakUpperBound(it, end, threshold));
}

template<typename RandomAccessIter, typename T, typename Compare>
static std::pair<RandomAccessIter, RandomAccessIter> SearchPeakRange(RandomAccessIter begin, RandomAccessIter end,  const T& threshold, Compare comp) {
	//begin: begin of the sequence, [begin, end)
	//end: end of the sequence, [begin, end)
	//threshold: threshold for peak search (suppress noise)
	//comp: define <
	//return: the lower and upper boundary for peak range
	//tips: (1)return.first and return.first+1 can calculate the precise intersection point
	//		(2)return.second and return.second-1 can calculate the precise intersection point
	//      (3)max of [return.first, return.second) is the peak value
	//***********************************************
	RandomAccessIter it = SearchPeakLowerBound(begin, end, threshold, comp);
	return make_pair(it, SearchPeakUpperBound(it, end, threshold, comp));
}

//statistic
template<typename ForwardIter>
static double Average(ForwardIter begin, ForwardIter end) {
	//begin: sequence begin
	//end: sequence end
	//return:   1
	//         ---(X1+X2+...+Xn)
	//          n
	//***********************************************
	int length = end - begin;
	return std::accumulate(begin, end, double(0)) / length;
}

template<typename ForwardIter>
static double StandardDeviation(ForwardIter begin, ForwardIter end) {
	//begin: sequence begin
	//end: sequence end
	//return:   1  
	//   sqrt( ----((X1-u)^2+(X2-u)^2+...+(Xn-u)^2) )
	//          n-1
	//tips: (1)u is the average
	//***********************************************
	int length = end - begin;
	double average = Average(begin, end);
	double sum = 0;
	while(begin != end) {
		sum += ((*begin-average)*(*begin-average));
		begin++;
	}
	return sqrt(sum/(length-1));
}

//move calculation
template<typename RandomAccessIter, typename OutputIter, typename BinaryFunction>
static OutputIter MoveCalculation(RandomAccessIter begin, RandomAccessIter end, OutputIter dst, int _win, BinaryFunction f) {
	//begin: input sequence begin
	//end: input sequence end
	//dst: output sequence end
	//_win: moving window
	//f: binary function(iterator, _win)
	//return: output sequence end
	//tips: (1)_win must >=1
	//      (2)length=end-begin must >= _win
	//***********************************************
	if(end-begin < _win || _win < 1)return dst;
	while(begin != end-_win+1) {
		*dst++ = f(begin, _win);
		begin++;//must reach the end, so end is begin+_win, not begin+_win-1
	}
	return dst;
}


template<typename RandomAccessIter, typename OutputIter>
static OutputIter MoveAverage(RandomAccessIter begin, RandomAccessIter end, OutputIter dst, int _win) {
	//begin: input sequence begin
	//end: input sequence end
	//dst: output sequence end
	//_win: moving-average of _win points
	//return: output sequence end
	//tips: (1)_win must >=1
	//      (2)length=end-begin must >= _win
	//      (3)suggest that _win is odd number 
	//***********************************************
	if(end-begin < _win || _win < 1)return dst;
	while(begin != end-_win+1) {
		*dst++ = Average(begin, begin+_win);
		begin++;//must reach the end, so end is begin+_win, not begin+_win-1
	}
	return dst;
}
template<typename RandomAccessIter1, typename RandomAccessIter2, typename OutputIter1, typename OutputIter2>
static std::pair<OutputIter1, OutputIter2> MoveAverage(RandomAccessIter1 begin1, RandomAccessIter1 end1, RandomAccessIter2 begin2, 
	OutputIter1 dst1, OutputIter2 dst2, int _win) {
	//begin1: input sequence begin, to operation
	//end1: input sequence end, to operation
	//begin2: input sequence begin, to keep
	//dst1: output sequence end, to operation
	//dst2: output sequence end, to keep
	//_win: moving-average of _win points
	//return: output sequence end
	//tips: (1)_win must >=1
	//      (2)length=end1-begin1 must >= _win
	//      (3)suggest that _win is odd number 
	//***********************************************
	if(end1-begin1 < _win || _win < 1)return make_pair(dst1, dst2);
	while(begin1 != end1-_win+1) {
		*dst1++ = Average(begin1, begin1+_win);
		*dst2++ = *(begin2 + _win/2);
		begin1++;begin2++;//must reach the end, so end is begin+_win, not begin+_win-1
	}
	return make_pair(dst1, dst2);
}


template<typename RandomAccessIter, typename OutputIter>
static OutputIter MoveNormalize(RandomAccessIter begin, RandomAccessIter end, OutputIter dst, int _win) {
	//begin: input sequence begin
	//end: input sequence end
	//dst: output sequence end
	//_win: moving-average of _win points
	//return: output sequence end
	//tips: (1)_win must >=1
	//      (2)length=end-begin must >= _win
	//      (3)suggest that _win is odd number 
	//***********************************************
	if(end-begin < _win || _win < 1)return dst;
	while(begin != end-_win+1) {
		*dst++ = *(begin + _win/2) / Average(begin, begin+_win);
		begin++;//must reach the end, so end is begin+_win, not begin+_win-1
	}
	return dst;
}
template<typename RandomAccessIter1, typename RandomAccessIter2, typename OutputIter1, typename OutputIter2>
static std::pair<OutputIter1, OutputIter2> MoveNormalize(RandomAccessIter1 begin1, RandomAccessIter1 end1, RandomAccessIter2 begin2, 
	OutputIter1 dst1, OutputIter2 dst2, int _win) {
	//begin1: input sequence begin, to operation
	//end1: input sequence end, to operation
	//begin2: input sequence begin, to keep
	//dst1: output sequence end, to operation
	//dst2: output sequence end, to keep
	//_win: moving-average of _win points
	//return: output sequence end
	//tips: (1)_win must >=1
	//      (2)length=end1-begin1 must >= _win
	//      (3)suggest that _win is odd number 
	//***********************************************
	if(end1-begin1 < _win || _win < 1)return make_pair(dst1, dst2);
	while(begin1 != end1-_win+1) {
		*dst1++ = *(begin1 + _win/2) / Average(begin1, begin1+_win);
		*dst2++ = *(begin2 + _win/2);
		begin1++; begin2++;//must reach the end, so end is begin+_win, not begin+_win-1
	}
	return make_pair(dst1, dst2);
}


template<typename RandomAccessIter, typename OutputIter>
static OutputIter MoveStandardDeviation(RandomAccessIter begin, RandomAccessIter end, OutputIter dst, int _win) {
	//begin: input sequence begin
	//end: input sequence end
	//dst: output sequence end
	//_win: moving-average of _win points
	//return: output sequence end
	//tips: (1)_win must >=1
	//      (2)length=end-begin must >= _win
	//      (3)suggest that _win is odd number 
	//***********************************************
	if(end-begin < _win || _win < 1)return dst;
	while(begin != end-_win+1) {
		*dst++ = StandardDeviation(begin, begin+_win);
		begin++;//must reach the end, so end is begin+_win, not begin+_win-1
	}
	return dst;
}

template<typename RandomAccessIter1, typename RandomAccessIter2, typename OutputIter1, typename OutputIter2>
static std::pair<OutputIter1, OutputIter2> MoveStandardDeviation(RandomAccessIter1 begin1, RandomAccessIter1 end1, RandomAccessIter2 begin2, 
	OutputIter1 dst1, OutputIter2 dst2, int _win) {
	//begin1: input sequence begin, to operation
	//end1: input sequence end, to operation
	//begin2: input sequence begin, to keep
	//dst1: output sequence end, to operation
	//dst2: output sequence end, to keep
	//_win: moving-average of _win points
	//return: output sequence end
	//tips: (1)_win must >=1
	//      (2)length=end1-begin1 must >= _win
	//      (3)suggest that _win is odd number 
	//***********************************************
	if(end1-begin1 < _win || _win < 1)return make_pair(dst1, dst2);
	while(begin1 != end1-_win+1) {
		*dst1++ = StandardDeviation(begin1, begin1+_win);
		*dst2++ = *(begin2 + _win/2);
		begin1++; begin2++;//must reach the end, so end is begin+_win, not begin+_win-1
	}
	return make_pair(dst1, dst2);
}


#endif //MyMath_H_INCLUDE
//~
