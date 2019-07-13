//Coordinate.h
//~~~~~~~~~~~~~~~~
#ifndef COORDINATE_H_INCLUDE
#define COORDINATE_H_INCLUDE

#include<iostream>
#include<sstream>
#include<string>
#include<iomanip>
#include<vector>
#include<iterator>
#include<algorithm>
#include<functional>
#include<numeric>
#include<regex>
#include<cmath>
#include"Matrix.h"

//Point in Free Space
class Cartesian;//(x, y, z)
class Polar;//(r, lamda, phi)
class Cylinder;//(tho, lamda, z)
//Point in Ellipse
class Ellipse;
class ElpCentric;//(ca, <u>) 
class ElpFocal;//(fa)
class ElpECentric;//(ea, <phi>)
class ElpNormal;//(na, <b>)
//point in Ellipsoid
class Ellipsoid;
class Geodetic;//(b, l, h, <na>)  derived from ElpNormal

class Coordinate {
protected:
	static double azimuthReverse(double _alpha, bool _zs) { return (_zs)?(360 - _alpha):(-_alpha);}
	static int calQuadrant(double _pb, double _nb) { return (_pb>0)?((_nb>0)?(1):(4)):((_nb>0)?(2):(3));}
public:	
	//const parameters
	
	//static functions
	//**azimuth
	static double calAzimuth(double _x, double _y, bool _ax, bool _rvs, bool _zs);//x-horizontal, y-vertical
	static double calAzimuthCartesian(double _x, double _y) { return calAzimuth(_x, _y, true, true, true);}
	static double calAzimuthGauss(double _x, double _y) { return calAzimuth(_y, _x, false, false, false);}
	static double azimuthZS2ZM(double _alpha) { return (GEQUAL<double>(_alpha, 180))?(_alpha-360):(_alpha); };//[0,360) to [-180,180)
	static double azimuthZM2ZS(double _alpha) { return (GEQUAL<double>(_alpha, 0))?(_alpha):(360+_alpha);}//[-180,180) to [0,360)}
	static double calElevation(double _x, double _y, bool _zs);//x-horizontal, y-vertical
	static double elevationZS2ZM(double _alpha) { return 90 - _alpha;}//[0,180] to [-90,90]
	static double elevationZM2ZS(double _alpha) { return 90 - _alpha;}//[-90,90] to [0,180]
	static double azimuthZS2elevationZM(double _alpha) { return (_alpha>270)?(_alpha-360):((_alpha>90)?(180-_alpha):(_alpha));}//[0, 360) to [-90, 90]
	static double elevationZM2azimuthZS(double _alpha) { return azimuthZM2ZS(_alpha);}//[-90, 90] to [0, 360)
};

//*********************************Free Space****************************************//

class Cartesian {
protected:
	double x, y, z;
public:
	//con-set-get
	//**construct
	Cartesian():x(0),y(0),z(0){};
	Cartesian(double _x, double _y = 0, double _z = 0):x(_x), y(_y), z(_z){};
	Cartesian(const Cartesian&p):x(p.x), y(p.y), z(p.z){};
	Cartesian(const Matrix<double>&p){this->set(p);};
	//**set
	void set(double _x, double _y = 0, double _z = 0) {
		this->x = _x; this->y = _y; this->z = _z; 
	};
	void set(const Cartesian&p) {
		this->x = p.x; this->y = p.y; this->z = p.z; 
	}
	void set(const Matrix<double>&p) {
		this->x = p.getElem(1,1); 
		this->y = p.getElem(2,1); 
		this->z = p.getElem(3,1); 
	};
	void setX(double _x = 0) {this->x = _x;};
	void setY(double _y = 0) {this->y = _y;};
	void setZ(double _z = 0) {this->z = _z;};
	//**get
	const double getX()const {return this->x;};
	const double getY()const {return this->y;};
	const double getZ()const {return this->z;};
	const double getR()const {return sqrt(x*x+y*y+z*z);};
	const double getR2()const {return x*x+y*y+z*z;};
	const double getRho()const {return sqrt(x*x+y*y);};
	const double getRho2()const {return x*x+y*y;};
	const Matrix<double> getMat()const {//Matrix: 3*1
		double data[3] = {x, y, z};
		Matrix<double> dst(3, 1, data);
		return dst;
	}
	friend std::ostream& operator<< (std::ostream&os, const Cartesian&p) {
		using namespace std;
		char fillc = os.fill();
		os<<setiosflags(ios::fixed | ios::internal | ios::showpos)<<setfill('0')
		  <<setw(12)<<setprecision(3)<<p.x<<" "
		  <<setw(12)<<setprecision(3)<<p.y<<" "
		  <<setw(12)<<setprecision(3)<<p.z<<" "
		  <<resetiosflags(ios::fixed | ios::internal | ios::showpos)<<setfill(fillc);
		return os;
	}
	//**transform
	//****point in free sapce
	Polar toPolar()const;
	Cylinder toCylinder()const;
	//****Point in Ellipse
	ElpCentric toElpCentric(const Ellipse&p)const;
	ElpFocal toElpFocal(const Ellipse&p)const;
	ElpECentric toElpECentric(const Ellipse&p)const;
	ElpNormal toElpNormal(const Ellipse&p)const;
	//****point in Ellipsoid
	Geodetic toGeodetic(const Ellipsoid&p)const;
	//**calculation
	Cartesian operator+(const Cartesian&p)const {
		return Cartesian(this->x+p.x, this->y+p.y, this->z+p.z);
	}
	Cartesian operator-(const Cartesian&p)const {
		return Cartesian(this->x-p.x, this->y-p.y, this->z-p.z);
	}
	Cartesian operator*(const Cartesian&p)const {//cross product
		return Cartesian(y*p.z-z*p.y, z*p.x-x*p.z, x*p.y-y*p.x);
	}
	friend Cartesian operator+(const Cartesian&p, double _elem) {
		return Cartesian(p.getX()+_elem, p.getY()+_elem, p.getZ()+_elem);
	}
	friend Cartesian operator+(double _elem, const Cartesian&p) {
		return Cartesian(_elem+p.getX(), _elem+p.getY(), _elem+p.getZ());
	}
	friend Cartesian operator-(const Cartesian&p, double _elem) {
		return Cartesian(p.getX()-_elem, p.getY()-_elem, p.getZ()-_elem);
	}
	friend Cartesian operator-(double _elem, const Cartesian&p) {
		return Cartesian(_elem-p.getX(), _elem-p.getY(), _elem-p.getZ());
	}
	friend Cartesian operator*(const Cartesian&p, double _elem) {
		return Cartesian(p.getX()*_elem, p.getY()*_elem, p.getZ()*_elem);
	}
	friend Cartesian operator*(double _elem, const Cartesian&p) {
		return Cartesian(_elem*p.getX(), _elem*p.getY(), _elem*p.getZ());
	}
	Cartesian axisReverse(FRAMEAXIS _ai)const { 
		//axial reverse
		switch(_ai) {
			case AXISX:return Cartesian(-x, y, z);
			case AXISY:return Cartesian(x, -y, z);
			case AXISZ:return Cartesian(x, y, -z);
		}
	}
	Cartesian normlize()const { 
		double r = this->getR(); 
		return Cartesian(x/r, y/r, z/r);
	};
	Cartesian cross(const Cartesian&p)const {
		//cross product
		return Cartesian(y*p.z-z*p.y, z*p.x-x*p.z, x*p.y-y*p.x);
	}
	Cartesian project(const Cartesian&p)const {
		//project the [this] to the vector of [p]
		return this->getR()*this->calAngleCos(p)*p.normlize();
	}
	Cartesian project(const Cartesian&p1, const Cartesian&p2)const {
		//project the [this] to the plane of p1xxp2
		return *this - this->project(p1*p2);
	}
	Cartesian rotationX(double _alpha)const {
		//rotation, ref-axis equal X-axis
		return Cartesian(DOUBLEMAT::creRotX(-_alpha)*this->getMat());
		//[-] is because we move vector, not axis
	}
	Cartesian rotationY(double _alpha)const {
		//rotation, ref-axis equal Y-axis
		return Cartesian(DOUBLEMAT::creRotY(-_alpha)*this->getMat());
		//[-] is because we move vector, not axis
	}
	Cartesian rotationZ(double _alpha)const {
		//rotation, ref-axis equal Z-axis
		return Cartesian(DOUBLEMAT::creRotZ(-_alpha)*this->getMat());
		//[-] is because we move vector, not axis
	}	
	Cartesian rotation(double _alpha, FRAMEAXIS _axis)const {
		//rotation, ref-axis specified
		return Cartesian(DOUBLEMAT::creRot(-_alpha, _axis)*this->getMat());
		//[-] is because we move vector, not axis
	}
	Cartesian rotation(double _alpha, Cartesian&p)const {
		//rotation, arbitrary-axis, unit normal vector-p
		return Cartesian(DOUBLEMAT::creRot(-_alpha, p.x, p.y, p.z)*this->getMat());
		//[-] is because we move vector, not axis
	}
	Cartesian triHeight(const Cartesian&p)const {
		//height vector of triangle
		double a = this->getR(), b = p.getR(), c = this->calDot(p);
		double a2 = a*a, b2 = b*b;
		double coeff0 = a2 + b2 - 2*c, coeff1 = b2 -c, coeff2 = a2 -c;
		return coeff1/coeff0*(*this) + coeff2/coeff0*p;
	}
	double calDot(const Cartesian&p)const {
		//dot product
		return x*p.x+y*p.y+z*p.z;
	}
	double calDistance(const Cartesian&p)const {
		//distance
		return Cartesian(x-p.x, y-p.y, z-p.z).getR();
	}
	double calAngle(const Cartesian&p)const {
		//vector angle, [0, 180]degree
		return DACOS(this->calAngleCos(p));
	}
	double calAngleCos(const Cartesian&p)const {
		//cos of vector angle
		return this->calDot(p)/(this->getR()*p.getR());
	}
};

class Polar {
protected:
	double r;//distance
	double lambda, phi;//angle, degree, [0, 360)  [-90, 90)
public:
	//con-set-get
	//**construct
	Polar():r(0), lambda(0), phi(0){};
	Polar(double _r, double _lambda, double _phi = 0):r(_r), lambda(_lambda), phi(_phi){};
	Polar(const Polar&p):r(p.r), lambda(p.lambda), phi(p.phi){};
	//**set
	void set(double _r, double _lambda, double _phi = 0) {
		this->r = _r; this->lambda = _lambda; this->phi = _phi; 
	};
	void set(const Polar&p) {
		this->r = p.r; this->lambda = p.lambda; this->phi = p.phi; 
	};
	void setR(double _r = 0) {this->r = _r;};
	void setLambda(double _lambda = 0) {this->lambda = _lambda;};
	void setPhi(double _phi = 0) {this->phi = _phi;};
	//**get
	const double getR()const {return this->r;};
	const double getLambda()const {return this->lambda;};
	const double getPhi()const {return this->phi;};
	friend std::ostream& operator<< (std::ostream&os, const Polar&p) {
		using namespace std;
		char fillc = os.fill();
		os<<setiosflags(ios::fixed | ios::internal | ios::showpos)<<setfill('0')
		  <<setw(12)<<setprecision(3)<<p.r<<" "
		  <<setw(12)<<setprecision(7)<<p.lambda<<" "
		  <<setw(12)<<setprecision(7)<<p.phi<<" "
		  <<resetiosflags(ios::fixed | ios::internal | ios::showpos)<<setfill(fillc);
		return os;
	}
	//**transform
	Cartesian toCartesian()const;
};

class Cylinder {
protected:
	double rho, z;//distance
	double lambda;//angle, degree, [0, 360)
public:	
	//con-set-get
	//**construct
	Cylinder():rho(0), lambda(0), z(0){};
	Cylinder(double _rho, double _lambda, double _z = 0):rho(_rho), lambda(_lambda), z(_z){};
	Cylinder(const Cylinder&p):rho(p.rho), lambda(p.lambda), z(p.z){};
	//**set
	void set(double _rho, double _lambda, double _z = 0) {
		this->rho = _rho; this->lambda = _lambda; this->z = _z; 
	};
	void set(const Cylinder&p) {
		this->rho = p.rho; this->lambda = p.lambda; this->z = p.z; 
	};
	void setRho(double _rho = 0) {this->rho = _rho;};
	void setLambda(double _lambda = 0) {this->lambda = _lambda;};
	void setZ(double _z = 0) {this->z = _z;};
	//**get
	const double getRho()const {return this->rho;};
	const double getLambda()const {return this->lambda;};
	const double getZ()const {return this->z;};
	friend std::ostream& operator<< (std::ostream&os, const Cylinder&p) {
		using namespace std;
		char fillc = os.fill();
		os<<setiosflags(ios::fixed | ios::internal | ios::showpos)<<setfill('0')
		  <<setw(12)<<setprecision(3)<<p.rho<<" "
		  <<setw(12)<<setprecision(7)<<p.lambda<<" "
		  <<setw(12)<<setprecision(3)<<p.z<<" "
		  <<resetiosflags(ios::fixed | ios::internal | ios::showpos)<<setfill(fillc);
		return os;
	}
	//**transform
	Cartesian toCartesian()const;
};

//**************************************Point in Ellipse*************************//
class Ellipse {
protected:	
	//neccessary parameters
	double a;//semi-major axis
	double b;//semi-minor axis
	//additional parameters
	double alpha;//flattening
	double e;//first eccentricity
	double es;//second eccentricity
	//innner initial
	void initial_alpha(){ this->alpha= (a-b)/a; };
	void initial_e(){ this->e= sqrt(a*a-b*b)/a; };
	void initial_es(){ this->es= sqrt(a*a-b*b)/b; };
public:	
	//con-set-get
	enum Format {A_B, A_ALPHA, A_E, A_ES, A_E2, A_ES2, A_C};
	//**construct
	Ellipse():a(0), b(0), alpha(0), e(0), es(0){};
	Ellipse(double _p1, double _p2, Format _format){ this->set(_p1, _p2, _format); };
	Ellipse(const Ellipse&p){ this->set(p); };
	//**set
	void set(double _p1, double _p2, Format _format);
	void set(const Ellipse&p) {
		this->a = p.a; this->b = p.b; 
		this->alpha = p.alpha; this->e = p.e; this->es = p.es;
	}
	//**get
	const double getA()const {return this->a;};
	const double getB()const {return this->b;};
	const double getAlpha()const {return this->alpha;};
	const double getE()const {return this->e;};
	const double getEs()const {return this->es;};
	const double getC()const {return this->a*this->e;};
	const double getCs()const {return this->a*this->a/this->b;};//latus rectum
	const double getP()const {return this->getB2()/this->getC();};//distance between directrix and focal
	const double getA2()const {return this->a*this->a;};
	const double getB2()const {return this->b*this->b;};
	const double getE2()const {return this->e*this->e;};
	const double getEs2()const {return this->es*this->es;};
	const double getC2()const {return this->getC()*this->getC();};
};

class ElpNormal {
protected:
	//neccessary parameters
	double na;//degree, [0, 360)
	//additional parameters
	double b;//degree, [-90, 90], derived from  na
	//innner initial
	void initial_b() { this->b = Coordinate::azimuthZS2elevationZM(na);}
	void initial_na() { this->na = Coordinate::elevationZM2azimuthZS(b);}
public:
	//con-set-get
	//**construct
	ElpNormal(double _na = 0) {this->set(_na);};
	ElpNormal(bool, double _b = 0) {this->set(true, _b);}
	ElpNormal(const ElpNormal&p) {this->set(p);};
	//**set
	void set(double _na = 0) {
		this->na = _na; this->initial_b();
	};
	void set(bool, double _b = 0) {
		this->b = _b; this->initial_na();
	}
	void set(const ElpNormal&p) {
		this->na = p.na; this->b = p.b;
	};
	void setNA(double _na = 0) { this->set(_na);}
	void setB(double _b = 0) { this->set(true, _b);}
	//**get
	const double getNA()const {return this->na;};
	const double getB()const {return this->b;}
	friend std::ostream& operator<< (std::ostream&os, const ElpNormal&p) {
		using namespace std;
		char fillc = os.fill();
		os<<setiosflags(ios::fixed | ios::internal | ios::showpos)<<setfill('0')
		  <<setw(12)<<setprecision(7)<<p.na<<" "
		  <<resetiosflags(ios::fixed | ios::internal | ios::showpos)<<setfill(fillc);
		return os;
	}
	//transform
	Cartesian toCartesian(const Ellipse&p)const;
	//calculation
	static double calW(double _b, const Ellipse&p);//f(b)=f(na)
	static double calV(double _b, const Ellipse&p);//f(b)=f(na)
	static double calEta2(double _b, const Ellipse&p);//f(b)=f(na)
	static double calU(double _b, const Ellipse&p);//only f(b)
	static double calPhi(double _b, const Ellipse&p);//only f(b)
	double calW(const Ellipse&p)const { return calW(this->b, p);};//f(b)=f(na)
	double calV(const Ellipse&p)const { return calV(this->b, p);};//f(b)=f(na)
	double calEta2(const Ellipse&p)const { return calEta2(this->b, p);};//f(b)=f(na)
	double calU(const Ellipse&p)const { return calU(this->b, p);};//only f(b)
	double calPhi(const Ellipse&p)const { return calPhi(this->b, p);};//only f(b)
};

class ElpCentric {
protected:
	//neccessary parameters
	double ca;//degree, [0, 360)
	//additional parameters
	double phi;//degree, [-90, 90], derived from  ca
	//innner initial
	void initial_phi() { this->phi = Coordinate::azimuthZS2elevationZM(ca);}
	void initial_ca() { this->ca = Coordinate::elevationZM2azimuthZS(phi);}
public:
	//con-set-get
	//**construct
	ElpCentric(double _ca = 0) {this->set(_ca);};
	ElpCentric(bool, double _phi = 0) {this->set(true, _phi);}
	ElpCentric(const ElpCentric&p) {this->set(p);};
	//**set
	void set(double _ca = 0) {
		this->ca = _ca; this->initial_phi();
	};
	void set(bool, double _phi = 0) {
		this->phi = _phi; this->initial_ca();
	};
	void set(const ElpCentric&p) {
		this->ca = p.ca; this->phi = p.phi;
	};
	void setCA(double _ca = 0) { this->set(_ca);}
	void setPhi(double _phi = 0) { this->set(true, _phi);}
	//**get
	const double getCA()const {return this->ca;};
	const double getPhi()const {return this->phi;};
	friend std::ostream& operator<< (std::ostream&os, const ElpCentric&p) {
		using namespace std;
		char fillc = os.fill();
		os<<setiosflags(ios::fixed | ios::internal | ios::showpos)<<setfill('0')
		  <<setw(12)<<setprecision(7)<<p.ca<<" "
		  <<resetiosflags(ios::fixed | ios::internal | ios::showpos)<<setfill(fillc);
		return os;
	}
	//transform
	Cartesian toCartesian(const Ellipse&p)const;
	//calculation
	static double calB(double _phi, const Ellipse&p);//only f(phi)
	double calB(const Ellipse&p)const { return calB(this->phi, p);};//only f(phi)
};

class ElpECentric {
protected:
	//neccessary parameters
	double ea;//degree, [0, 360)
	//additional parameters
	double u;//degree, [-90, 90], derived from  ea
	//innner initial
	void initial_u() { this->u = Coordinate::azimuthZS2elevationZM(ea);}
	void initial_ea() { this->ea = Coordinate::elevationZM2azimuthZS(u);}
public:
	//con-set-get
	//**construct
	ElpECentric(double _ea = 0) {this->set(_ea);};
	ElpECentric(bool, double _u = 0) {this->set(true, _u);}
	ElpECentric(const ElpECentric&p) {this->set(p);};
	//**set
	void set(double _ea = 0) {
		this->ea = _ea; this->initial_u();
	};
	void set(bool, double _u = 0) {
		this->u = _u; this->initial_ea();
	}
	void set(const ElpECentric&p) {
		this->ea = p.ea; this->u = p.u;
	};
	void setEA(double _ea = 0) { this->set(_ea);}
	void setU(double _u = 0) { this->set(true, _u);}
	//**get
	const double getEA()const {return this->ea;};
	const double getU()const {return this->u;}
	friend std::ostream& operator<< (std::ostream&os, const ElpECentric&p) {
		using namespace std;
		char fillc = os.fill();
		os<<setiosflags(ios::fixed | ios::internal | ios::showpos)<<setfill('0')
		  <<setw(12)<<setprecision(7)<<p.ea<<" "
		  <<resetiosflags(ios::fixed | ios::internal | ios::showpos)<<setfill(fillc);
		return os;
	}
	//transform
	Cartesian toCartesian(const Ellipse&p)const;
	//calculation
	static double calFocalR(double _ea, const Ellipse&p);
	static double calB(double _u, const Ellipse&p);//only f(phi)
	double calB(const Ellipse&p)const { return calB(this->u, p);};//only f(phi)
	double calFocalR(const Ellipse&p)const { return calFocalR(this->ea, p);};//only f(ea)
};

class ElpFocal {
protected:
	double fa;//degree, [0, 360)
public:
	//con-set-get
	//**construct
	ElpFocal(double _fa = 0) {this->set(_fa);};
	ElpFocal(const ElpFocal&p) {this->set(p);};
	//**set
	void set(double _fa = 0) {
		this->fa = _fa;
	};
	void set(const ElpFocal&p) {
		this->fa = p.fa;
	};
	//**get
	const double getFA()const {return this->fa;};
	friend std::ostream& operator<< (std::ostream&os, const ElpFocal&p) {
		using namespace std;
		char fillc = os.fill();
		os<<setiosflags(ios::fixed | ios::internal | ios::showpos)<<setfill('0')
		  <<setw(12)<<setprecision(7)<<p.fa<<" "
		  <<resetiosflags(ios::fixed | ios::internal | ios::showpos)<<setfill(fillc);
		return os;
	}
	//transform
	Cartesian toCartesian(const Ellipse&p)const;
	//calculation
	double calFocalR(const Ellipse&p)const;
};

//**************************************Point in Ellipsoid*************************//
class Ellipsoid : public Ellipse {
protected:
public:
	static const Ellipsoid WGS84;
	static const Ellipsoid CGCS2000;
	static const Ellipsoid KRASSOVSKY;
	static const Ellipsoid IE75;
	static const Ellipsoid EARTH_SPHERIOD;
	//con-set-get
	//**construct
	Ellipsoid():Ellipse(){};
	Ellipsoid(double _p1, double _p2, Format _format):Ellipse(_p1, _p2, _format){};
	Ellipsoid(const Ellipsoid&p):Ellipse(p){};

};
class Geodetic : public ElpNormal {
protected:
	double l;//degree, [-180, 180]
	double h;//m
public:	
	//con-set-get
	//**construct
	Geodetic(double _b = 0, double _l = 0, double _h = 0) {this->set(_b, _l, _h);};
	Geodetic(const Geodetic&p){this->set(p);};
	//**set
	void set(double _b = 0, double _l = 0, double _h = 0) {
		this->ElpNormal::set(true, _b); 
		this->l = _l; this->h = _h;
	}
	void set(const Geodetic&p) {
		this->ElpNormal::set(p); 
		this->l = p.l; this->h = p.h;
	}
	void setL(double _l = 0) { this->l = _l;}
	void setH(double _h = 0) { this->h = _h;}
	//**get
	const double getL()const{ return this->l; }
	const double getH()const{ return this->h; }
	friend std::ostream& operator<< (std::ostream&os, const Geodetic&p) {
		using namespace std;
		char fillc = os.fill();
		os<<setiosflags(ios::fixed | ios::internal | ios::showpos)<<setfill('0')
		  <<setw(12)<<setprecision(7)<<p.b<<" "
		  <<setw(12)<<setprecision(7)<<p.l<<" "
		  <<setw(12)<<setprecision(3)<<p.h<<" "
		  <<resetiosflags(ios::fixed | ios::internal | ios::showpos)<<setfill(fillc);
		return os;
	}
	//transform
	Cartesian toCartesian(const Ellipsoid&p)const;
	//calculation
	static double calM(double _b, const Ellipsoid&p);//f(b)=f(na)
	static double calN(double _b, const Ellipsoid&p);//f(b)=f(na)
	double calM(const Ellipsoid&p)const { return calM(this->b, p);};//f(b)=f(na)
	double calN(const Ellipsoid&p)const { return calN(this->b, p);};//f(b)=f(na)

};
#endif //COORDINATE_H_INCLUDE
//~