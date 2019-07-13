//Coordinate.cpp
//~~~~~~~~~~~~~~~~
#include"Coordinate.h"
using namespace std;

//class Coordinate
double Coordinate::calAzimuth(double _x, double _y, bool _ax, bool _rvs, bool _zs) {
	//_x: x-horizontal, _y: y-vertical
	//_ax: reference axis is x-horizontal when true
	//_rvs: positive oriention is reverse when true
	//_zs: zero degree is the start when true
	//pb: reference base axis, nb: not reference base axis
	//return: degree
	double pb = (_ax)?(_x):(_y), nb = (_ax)?(_y):(_x);
	if(abs(nb)<EPS)return pb<0?(180):(0);
	double base = DATAN(abs(nb/pb)), dst = 0;//atan will return base=90 when pb=0
	switch(Coordinate::calQuadrant(pb, nb)) {
		case 1:dst = base;break;
		case 2:dst = 180 - base;break;
		case 3:dst = (_zs)?(180 + base):(base - 180); break;
		case 4:dst = (_zs)?(360 - base):(-base);break;
	}
	_rvs = (_ax)?(_rvs):(!_rvs);
	dst = (_rvs)?(dst):(Coordinate::azimuthReverse(dst, _zs));
	return dst;
}
double Coordinate::calElevation(double _x, double _y, bool _zs) {
	//_x: x-horizontal, _y: y-vertical
	//_zs: zero degree is the start when true
	//return: degree
	double _r = sqrt(_x*_x + _y*_y);
	return (_zs)?(DACOS(_y/_r)):(DASIN(_y/_r));
}

//class Cartesian
Polar Cartesian::toPolar()const {
	double _r = this->getR();
	double _lambda = Coordinate::calAzimuthCartesian(this->x, this->y);//degree
	double _phi = DASIN(this->z/_r);//degree
	return Polar(_r, _lambda, _phi);
};
Cylinder Cartesian::toCylinder()const {
	double _rho = this->getRho();
	double _lambda = Coordinate::calAzimuthCartesian(this->x, this->y);//degree
	return Cylinder(_rho, _lambda, this->z);
};
ElpCentric Cartesian::toElpCentric(const Ellipse&p)const {
	double _ca = Coordinate::calAzimuthCartesian(this->x, this->y);
	return ElpCentric(_ca);
};
ElpFocal Cartesian::toElpFocal(const Ellipse&p)const {
	double c = p.getC();
	double _fa = Coordinate::calAzimuthCartesian(this->x-c, this->y);
	return ElpFocal(_fa);
};
ElpECentric Cartesian::toElpECentric(const Ellipse&p)const {
	double a = p.getA(), b = p.getB();
	double _ea = Coordinate::calAzimuthCartesian(this->x/a, this->y/b);
	return ElpECentric(_ea);
};
ElpNormal Cartesian::toElpNormal(const Ellipse&p)const {
	double a = p.getA(), b = p.getB(), e2 = p.getE2();
	double _na = Coordinate::calAzimuthCartesian(this->x/a*sqrt(1-e2), this->y/b);
	return ElpNormal(_na);
};
Geodetic Cartesian::toGeodetic(const Ellipsoid&p)const {
	//L
	double _l = Coordinate::calAzimuth(this->x, this->y, true, true, false);//degree, [-180,180)
	//B, H
	double r0 = this->getRho(), N = 0, _b, _h;
	double cs = p.getCs(), e2 = p.getE2(), es2 = p.getEs2();
	if(abs(r0)<EPS) {//polar point
		_b = (this->z>0)?(90):(-90); _h = this->z - p.getB();
		return Geodetic(_b, _l, _h);
	}
	double ti = this->z/r0, pp = cs*e2/r0, k = 1 + es2, t0;
	do {
		t0 = ti;
		ti = this->z/r0 + pp*t0/sqrt(k+t0*t0);
	} while(abs(ti-t0)>EPS);
    _b = DATAN(ti);//degree
	N = Geodetic::calN(_b, p); _h = this->z*SIND(_b) + r0*COSD(_b) - N + N*e2*SIND(_b)*SIND(_b);
	return Geodetic(_b, _l, _h);
}

//class Polar
Cartesian Polar::toCartesian()const {
	double _x = this->r*COSD(this->phi)*COSD(this->lambda);
	double _y = this->r*COSD(this->phi)*SIND(this->lambda);
	double _z = this->r*SIND(this->phi);
	return Cartesian(_x, _y, _z);
};

//class cylinder
Cartesian Cylinder::toCartesian()const {
	double _x = this->rho*COSD(this->lambda);
	double _y = this->rho*SIND(this->lambda);
	return Cartesian(_x, _y, this->z);
};

//***************************Point in Ellipse*******************************//
//class Ellipse
void Ellipse::set(double _p1, double _p2, Format _format) {
		switch(_format) {
		case A_B:
			//input parameters
			this->a = _p1; this->b = _p2;
			//other parameters
			this->initial_alpha(); this->initial_e(); this->initial_es();
			break;
		case A_ALPHA:
			//input parameters
			this->a = _p1; this->alpha = _p2;
			//other parameters
			this->b = a*(1-alpha);
			this->initial_e(); this->initial_es();
			break;
		case A_E:
			//input parameters
			this->a = _p1; this->e = _p2;
			//other parameters
			this->b = a*sqrt(1-e*e);
			this->initial_alpha(); this->initial_es();
			break;
		case A_ES:
			//input parameters
			this->a = _p1; this->es = _p2;
			//other parameters
			this->b = a/sqrt(1+es*es);
			this->initial_alpha(); this->initial_e();
			break;
		case A_E2:
			//input parameters
			this->a = _p1; this->e = sqrt(_p2);
			//other parameters
			this->b = a*sqrt(1-_p2);
			this->initial_alpha(); this->initial_es();
			break;
		case A_ES2:
			//input parameters
			this->a = _p1; this->es = sqrt(_p2);
			//other parameters
			this->b = a/sqrt(1+_p2);
			this->initial_alpha(); this->initial_e();
			break;
		case A_C:
			//input parameters
			this->a = _p1; this->b = sqrt(a*a-_p2*_p2);
			//other parameters
			this->initial_alpha(); this->initial_e(); this->initial_es();
			break;
	}
}

//class ElpNormal
Cartesian ElpNormal::toCartesian(const Ellipse&p)const {
	double a = p.getA(), b = p.getB();
	double W = this->calW(p), V = this->calV(p);
	double _x = a*COSD(this->na)/W;
	double _y = b*SIND(this->na)/V;
	return Cartesian(_x, _y);
};
double ElpNormal::calW(double _b, const Ellipse&p) {
	return sqrt(1-p.getE2()*SIND(_b)*SIND(_b));
};
double ElpNormal::calV(double _b, const Ellipse&p) {
	return sqrt(1+p.getEs2()*COSD(_b)*COSD(_b));
};
double ElpNormal::calEta2(double _b, const Ellipse&p) {
	return p.getEs2()*COSD(_b)*COSD(_b);
};
double ElpNormal::calU(double _b, const Ellipse&p) {
	return DATAN(TAND(_b)*p.getB()/p.getA());
};
double ElpNormal::calPhi(double _b, const Ellipse&p) {
	return DATAN(TAND(_b)*(1-p.getE2()));
};

//class ElpCentric
Cartesian ElpCentric::toCartesian(const Ellipse&p)const {
	double a = p.getA(), b = p.getB(), tc = TAND(this->ca);
	double a2 = a*a, b2 = b*b, tc2 = tc*tc;
	double _x = a*b/sqrt(b2+a2*tc2)*SIGNPOS<double>(COSD(this->ca));
	double _y = _x*tc;
	return Cartesian(_x, _y);
};
double ElpCentric::calB(double _phi, const Ellipse&p) {
	return DATAN(TAND(_phi)/(1-p.getE2()));
};

//class ElpECentric
Cartesian ElpECentric::toCartesian(const Ellipse&p)const {
	double a = p.getA(), b = p.getB();
	double _x = a*COSD(this->ea);
	double _y = b*SIND(this->ea);
	return Cartesian(_x, _y);
};
double ElpECentric::calB(double _u, const Ellipse&p) {
	return DATAN(TAND(_u)*p.getA()/p.getB());
};
double ElpECentric::calFocalR(double _ea, const Ellipse&p) {
	return p.getA()*(1-p.getE()*COSD(_ea));
};
//class ElpFocal
Cartesian ElpFocal::toCartesian(const Ellipse&p)const {
	double r = this->calFocalR(p), c = p.getC();
	double _x = r*COSD(this->fa) + c;
	double _y = r*SIND(this->fa);
	return Cartesian(_x, _y);
};
double ElpFocal::calFocalR(const Ellipse&p)const {
	double a = p.getA(), e = p.getE(), e2 = e*e;
	return a*(1-e2)/(1+e*COSD(this->fa));
};

//class Ellipsoid
const Ellipsoid Ellipsoid::WGS84(6378137.0, 0.00669437999013, Ellipsoid::A_E2);
const Ellipsoid Ellipsoid::CGCS2000(6378137.0, 0.00669438002290, Ellipsoid::A_E2);
const Ellipsoid Ellipsoid::KRASSOVSKY(6378245.0, 0.006693421622966, Ellipsoid::A_E2);
const Ellipsoid Ellipsoid::IE75(6378140.0, 0.006694384999588, Ellipsoid::A_E2);
const Ellipsoid Ellipsoid::EARTH_SPHERIOD(6378137.0, 0, Ellipsoid::A_E2);
//class Geodetic
Cartesian Geodetic::toCartesian(const Ellipsoid&p)const {
	double e2 = p.getE2(), N = Geodetic::calN(this->b, p);
	double _x = (N+this->h)*COSD(this->b)*COSD(this->l);
	double _y = (N+this->h)*COSD(this->b)*SIND(this->l);
	double _z = (N*(1-e2)+this->h)*SIND(this->b);
	return Cartesian(_x, _y, _z);
};
double Geodetic::calM(double _b, const Ellipsoid&p) {
	double W = Geodetic::calW(_b, p);
	return  p.getA()*(1-p.getE2())/(pow(W,3));
};
double Geodetic::calN(double _b, const Ellipsoid&p) {
	double W = Geodetic::calW(_b, p);
	return p.getA()/W;
};