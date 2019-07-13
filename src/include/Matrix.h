//Matrix.h
//20190331  V3.1
#ifndef Matrix_H_INCLUDE
#define Matrix_H_INCLUDE

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
#include"MyMath.h"

//****************************************************************
//Instruction
//1,All the  data_array ,vector start with 0 
//2,ALl the Matrix start with 1
//3,Store by row
//4,rows:row matrix; cols:col matrix; blocks:block matrix;
//  diags:diag matrix; ones: one-value matrix; zeros: zero-matrix; eyes: eye matrix
//5.diag:diag vector of a matrix; row: row value of a matrix; col: col value of a matrix
//  elem:someone value of a matrix
enum MATRIXERROR {MISMATCH, NOPOSITIVE, BEYOND, EMPTY, WRONGCAL, UNDEFINE, NOEXIST};

template<typename T>
class Matrix {
protected:
	int row;//number of rows
	int col;//number of columns
	std::vector<T> mat;//data of the Matrix
public:
	static void msgError(MATRIXERROR _ms);
	//con-set-get
	//**construct
	Matrix():row(0), col(0){};
	Matrix(int _row, int _col, T _elem){this->set(_row, _col, _elem);};
	Matrix(int _row, int _col, const std::vector<T>& _mat){this->set(_row, _col, _mat);};
	Matrix(int _row, int _col, T* _mat){this->set(_row, _col, _mat);};
	Matrix(int _row, int _col, T** _mat){this->set(_row, _col, _mat);};
	Matrix(int _row, int _col, std::istream&is){this->set(_row, _col, is);};
	Matrix(const Matrix<T>&p){this->set(p);}
	//Matrix(int _row, int _col, )
	//**set(setted empty when unmatched)
	void set(int _row, int _col, T _elem);
	void set(int _row, int _col, const std::vector<T>& _mat);
	void set(int _row, int _col, T* _mat);
	void set(int _row, int _col, T** _mat);
	void set(int _row, int _col, std::istream&is);
	void set(const Matrix<T>&p);
	void set(int _row, int _col);//reshape the row and col but keep mat
	T& operator()(int _i, int _j);//element can be modified
	Matrix<T>& operator=(const Matrix<T>&p){this->set(p.row, p.col, p.mat); return *this;};
	void setRows(const Matrix<T>&_rows, int _i);//replace rows-i
	void setCols(const Matrix<T>&_cols, int _j);//replace cols-j
	void setRows(const Matrix<T>&_rows, int _rs, int _re);//replace rows-rs and rows-re
	void setCols(const Matrix<T>&_cols, int _cs, int _ce);//replace cols-cs and cols-ce
	void setBlocks(const Matrix<T>&_blocks);//replace block from(1,1)
	void setBlocks(const Matrix<T>&_blocks, int _rs, int _cs);//repalce block from (_rs, _cs)
	void setDiags(const Matrix<T>&_diags);//replace diags
	//**get
	const int getRow()const { return this->row;}
	const int getCol()const { return this->col;}
	const int getIndex(int _i, int _j)const { return (_i - 1)*col + _j - 1;};
	const int getIndexInvRow(int _index)const{ return (col == 0)?(0):(_index/col+1);}
	const int getIndexInvCol(int _index)const{ return (col == 0)?(0):(_index%col+1);}
	const void getIndexInv(int _index, int&dsti, int&dstj)const { 
		dsti = (col == 0)?(0):(_index/col + 1); 
		dstj = (col == 0)?(0):(_index%col + 1);
	};
	const int getElem()const { return row*col;}//the number of elements
	const T getElem(int _i, int _j)const { return mat[getIndex(_i, _j)];}
	const Matrix<T> getRows(int _ri)const;
	const Matrix<T> getCols(int _cj)const;
	const Matrix<T> getRows(int _rs, int _re)const;
	const Matrix<T> getCols(int _cs, int _ce)const;
	const Matrix<T> getBlocks(int _rs, int _cs, int _re, int _ce)const;
	const Matrix<T> getBlocks(int _re, int _ce)const;
	const Matrix<T> getDiags()const;
	const void getMat(std::vector<T>&dst)const;
	const T getMax()const{ return *std::max_element(mat.begin(), mat.end());};
	const T getMin()const{ return *std::min_element(mat.begin(), mat.end());};
	const int getMaxRow()const{ return getIndexInvRow(std::max_element(mat.begin(), mat.end()) - mat.begin());};
	const int getMinRow()const{ return getIndexInvRow(std::min_element(mat.begin(), mat.end()) - mat.begin());};
	const int getMaxCol()const{ return getIndexInvCol(std::max_element(mat.begin(), mat.end()) - mat.begin());};
	const int getMinCol()const{ return getIndexInvCol(std::min_element(mat.begin(), mat.end()) - mat.begin());};
	const T getSum()const{ return std::accumulate(mat.begin(), mat.end(), T(0));}
	const T getProduct()const{ return std::accumulate(mat.begin(), mat.end(), T(1),std::multiplies<T>());}
	const T getSSQ()const{ return std::inner_product(mat.begin(), mat.end(), mat.begin(), T(0));}//Sum of SQuare
	const T getRSSQ()const{ return sqrt(std::inner_product(mat.begin(), mat.end(), mat.begin(), T(0)));}//Root Sum of SQuare
	const T getTr()const{ return this->getDiags().getSum();}
	const T getMax(int&dstr, int&dstc)const;
	const T getMin(int&dstr, int&dstc)const;
	const Matrix<T> getSimplestRows()const;
	const Matrix<T> getEchelonRows()const;
	const Matrix<T> getEchelonRows(int&flag)const;
	const int getNonZeroRow()const;
	const int getRank()const;
	const T getDet()const;
	const Matrix<T> getLPSubMatrix(int _p)const;//leading principle submatrix;p*p
	const T getLPMinor(int _p)const;//leading principle minor;det(p*p)
	const Matrix<T> getSumRows()const;
	const Matrix<T> getSumCols()const;
	const T getNormOne()const;
	const T getNormTwo()const;
	const T getNormInf()const;
	const T getCond()const{ return getNormInf()*invGJC().getNormInf();};
	const Matrix<T> getVTV()const { return tp()*(*this);};
	const Matrix<T> getVVT()const { return (*this)*tp();};
	friend std::ostream& operator<< (std::ostream&os, const Matrix<T>&p) {
		using namespace std;//Number(3),Decimal(3),Symbol(2)
		char fillc = os.fill();
		os<<setiosflags(ios::internal | ios::fixed | ios::showpos)<<setfill('0');
	  	for(int i = 1; i <= p.row; i++) {
			for(int j = 1; j <= p.col; j++)
				os<<setw(8)<<setprecision(3)<<p.getElem(i,j)<<" ";
			os<<endl;
		}
		os<<resetiosflags(ios::internal | ios::fixed | ios::showpos)<<setfill(fillc);
	  return os;
	}

	//calculation
	//**special matrix
	static Matrix<T> creOnes(int _index){ Matrix<T> dst(_index, _index, T(1)); return dst;};
	static Matrix<T> creOnes(int _row, int _col){ Matrix<T> dst(_row, _col, T(1)); return dst;};
	static Matrix<T> creZeros(int _index){ Matrix<T> dst(_index, _index, T(0)); return dst;};
	static Matrix<T> creZeros(int _row, int _col){ Matrix<T> dst(_row, _col, T(0)); return dst;};
	static Matrix<T> creDiags(const std::vector<T>&diags);
	static Matrix<T> creDiags(const Matrix<T>&diags);
	static Matrix<T> creEyes(int _index);
	static Matrix<T> creRotX(T _rx);//+ for inverse rotation, degree
	static Matrix<T> creRotY(T _ry);//+ for inverse rotation, degree
	static Matrix<T> creRotZ(T _rz);//+ for inverse rotation, degree
	static Matrix<T> creRot(T _r);//two-dimension
	static Matrix<T> creRot(T _r, FRAMEAXIS _a);//three-dimension, specified axis
	static Matrix<T> creRot(T _rx, T _ry, T _rz);//three-dimension, Rx*Ry*Rz*coordinate
	static Matrix<T> creRot(T _r, T _nx, T _ny, T _nz);//arbitrary axis of rotation, (nx, ny, nz)
	static Matrix<T> creVander(const std::vector<T>&vander);//Vandermonde matrix
	static Matrix<T> creStagger(const std::vector<T>&stagger, int _col);//stagger vector in matrix by cols-directrion
	static Matrix<T> creStaggers(const std::vector<T>&stagger, int _col);//S'*S in fast way
	//**plus, minus, multiply
	Matrix<T> operator+(const Matrix<T>&p)const;
	Matrix<T> operator-(const Matrix<T>&p)const;
	Matrix<T> operator*(const Matrix<T>&p)const;
	Matrix<T>& operator+=(const Matrix<T>&p);// cannot be const
	Matrix<T>& operator-=(const Matrix<T>&p);//cannot be const
	friend Matrix<T> operator+(const Matrix<T>&p, T _elem) {
		using namespace std;
		Matrix<T> dst(p.row, p.col, T(0));
		transform(p.mat.begin(), p.mat.end(), dst.mat.begin(), bind(plus<T>(),placeholders::_1, _elem));
		return dst;
	};
	friend Matrix<T> operator+(T _elem, const Matrix<T>&p) {
		using namespace std;
		Matrix<T> dst(p.row, p.col, T(0));
		transform(p.mat.begin(), p.mat.end(), dst.mat.begin(), bind(plus<T>(),placeholders::_1, _elem));
		return dst;
	};
	friend Matrix<T> operator-(const Matrix<T>&p, T _elem) {
		using namespace std;
		Matrix<T> dst(p.row, p.col, T(0));
		transform(p.mat.begin(), p.mat.end(), dst.mat.begin(), bind(minus<T>(),placeholders::_1, _elem));
		return dst;
	};
	friend Matrix<T> operator-(T _elem, const Matrix<T>&p) {
		using namespace std;
		Matrix<T> dst(p.row, p.col, T(0));
		transform(p.mat.begin(), p.mat.end(), dst.mat.begin(), bind(minus<T>(), _elem, placeholders::_1));
		return dst;
	};
	friend Matrix<T> operator*(const Matrix<T>&p, T _elem) {
		using namespace std;
		Matrix<T> dst(p.row, p.col, T(0));
		transform(p.mat.begin(), p.mat.end(), dst.mat.begin(), bind(multiplies<T>(), placeholders::_1, _elem));
		return dst;
	};
	friend Matrix<T> operator*(T _elem, const Matrix<T>&p) {
		using namespace std;
		Matrix<T> dst(p.row, p.col, T(0));
		transform(p.mat.begin(), p.mat.end(), dst.mat.begin(), bind(multiplies<T>(), _elem, placeholders::_1));
		return dst;
	};
	//**primary transform
	Matrix<T> exRows(int _ri, int _rj)const;//exchange row-i and row-j
	Matrix<T> exCols(int _ci, int _cj)const;//exchange col-i and col-j
	Matrix<T> mpRows(int _ri, T _elem)const;//multiply _elem for row-i
	Matrix<T> mpCols(int _ci, T _elem)const;//multiply _elem for col-i
	Matrix<T> plRows(int _ri, int _rj, T _elem)const;//plus _elem*row-j for row-i
	Matrix<T> plCols(int _ci, int _cj, T _elem)const;//plus _elem*row-j for col-i
	//**operation rows and cols
	Matrix<T> rpRows(const Matrix<T>&_rows, int _i)const;//replace rows-i
	Matrix<T> rpCols(const Matrix<T>&_cols, int _j)const;//replace cols-j
	Matrix<T> rpBlocks(const Matrix<T>&_blocks)const;//replace block
	Matrix<T> rpBlocks(const Matrix<T>&_blocks, int _rs, int _cs)const;//replace block
	Matrix<T> rpDiags(const Matrix<T>&_diags);//replace diags
	Matrix<T> addRows(const Matrix<T>&_rows)const;//add rows-last
	Matrix<T> addCols(const Matrix<T>&_cols)const;//add cols-last
	Matrix<T> substractRows(int num)const;//substract rows-last-num
	Matrix<T> substractCols(int num)const;//substract cols-last-num
	Matrix<T> insertRows(const Matrix<T>&_rows, int _i)const;//insert rows-i
	Matrix<T> insertCols(const Matrix<T>&_cols, int _j)const;//insert cols-j
	Matrix<T> deleteRows(int _i)const;//delete rows-i
	Matrix<T> deleteCols(int _j)const;//delete cols-j
	Matrix<T> deleteRows(int _rs, int _re)const;//delete rows-rs to rows-re
	Matrix<T> deleteCols(int _cs, int _ce)const;//delete cols-cs to cols-ce
	//**merge
	static Matrix<T> mergeRows(const Matrix<T>&p1, const Matrix<T>&p2);//merge by row-oriention
	static Matrix<T> mergeCols(const Matrix<T>&p1, const Matrix<T>&p2);//merge by col-oriention
	static Matrix<T> mergeBlocks(const Matrix<T>&p1, const Matrix<T>&p2);//merge by diag-oriention
	//**flip
	Matrix<T> flipRows()const;//flip by row-oriention
	Matrix<T> flipCols()const;//flip by row-oriention
	Matrix<T> flipDiags()const;//flip by row-oriention
	//**sort
	Matrix<T> downZero()const;//down the rows-zero to bottum
	//*equation solution
	Matrix<T> solveHomo()const;//solve Homogenous equation
	//**factorization
	void facLU(Matrix<T>&L, Matrix<T>&U)const;
	void facQR(Matrix<T>&Q, Matrix<T>&R)const;
	//**orthgonization
	Matrix<T> orthSMT()const;
	//**normlization
	Matrix<T> norm()const;
	Matrix<T> normRows()const;
	Matrix<T> normCols()const;
	//**eigenvalue and eigenvector
	Matrix<T> eigQR()const;//eigenvalue -- QR algorithm
	//**transpose
	Matrix<T> tp()const;
	//**inverse
	Matrix<T> invGJC()const;
	//**abs
	Matrix<T> absElems()const;
	//**inner_product
	T innerProduct(const Matrix<T>&p)const;
	static T innerProduct(const Matrix<T>&p1, const Matrix<T>&p2);
	//status judgement
	void clear()const{ this->row = 0; this->col = 0; this->mat.clear();}
	bool isInternal(int _i, int _j)const{ return _i > 0 && _i <= row && _j > 0 && _j <= col;}
	bool isInternalRow(int _i)const{ return _i > 0 && _i <= row;}
	bool isInternalCol(int _j)const{ return _j > 0 && _j <= col;}
	bool isInternalBlock(int _p)const{ return _p > 0 && _p <= MIN<int>(row, col);}
	bool isSquare()const{ return row == col;}
	bool isEmpty()const{ return row*col==0;}
	bool isSingle()const{ return row == 1 && col == 1;}
	bool isLessSize(const Matrix<T>&p)const{ return row < p.row && col < p.col;}
	bool isGreaterSize(const Matrix<T>&p)const{ return row > p.row && col > p.col;}
	bool isEqualSize(const Matrix<T>&p)const{ return row == p.row && col == p.col;}
	bool isEqual(const Matrix<T>&p)const{ return isEqualSize(p) && std::equal(mat.begin(), mat.end(), p.mat.begin());}
	bool isAppr(const Matrix<T>&p, T _ale)const{ return isEqualSize(p) && std::equal(mat.begin(), mat.end(), p.mat.begin(), 
		std::bind(APPR<T>, std::placeholders::_1, std::placeholders::_2, _ale));}
	bool isZero()const{ return isAppr(Matrix<T>(row, col, T(0)), 1e-6);}
	bool isMatchMultiply(const Matrix<T>&p)const{ return col == p.row;}//cur * p
	bool isMatchRows(const Matrix<T>&_rows)const{ return _rows.row == 1 && _rows.col == col;}//rows vector match
	bool isMatchCols(const Matrix<T>&_cols)const{ return _cols.col == 1 && _cols.row == row;}//cols vector match
	bool isInv()const{ return this->isSquare() && this->getRank() == row;}
	bool isVector()const{ return row == 1 || col == 1 && row*col != 0;}
	bool isUnitVector()const{ return isVector() && ABS<T>(getRSSQ()-T(1)) < EPS;}
	bool isUnitOneVector(int _i)const{ return isUnitVector() && ABS<T>(mat[_i-1]-T(1)) < EPS;}

};
//**********************class Matrix**************************
//static
template<typename T>
void Matrix<T>::msgError(MATRIXERROR _ms) {
	using namespace std;
	cout<<setiosflags(ios::left);
	switch(_ms) {
		case MISMATCH:
			cout<<setw(30)<<"dimension or size mismatch"<<endl;break;
		case NOPOSITIVE:
			cout<<setw(30)<<"dimension or size no positive"<<endl;break;
		case BEYOND:
			cout<<setw(30)<<"index visit beyond boundary"<<endl;break;
		case WRONGCAL:
			cout<<setw(30)<<"wrong caculation condition"<<endl;break;
		case UNDEFINE:
			cout<<setw(30)<<"undefined parameter type"<<endl;break;
		case NOEXIST:
			cout<<setw(30)<<"not exits this result"<<endl;break;
		default:
			cout<<setw(30)<<"error message"<<endl;break;
	}
	cout<<resetiosflags(ios::left);
}
//con-set-get
//**construct
//**set
template<typename T>
void Matrix<T>::set(int _row, int _col, T _elem) {
	this->row = _row; this->col = _col; this->mat.resize(_row*_col);
	std::fill(mat.begin(), mat.end(), _elem);
};

template<typename T>
void Matrix<T>::set(int _row, int _col, const std::vector<T>& _mat) {
	if(_row*_col <= 0) {
		msgError(NOPOSITIVE);exit(0);//underthinking
	}
	if(_row*_col != _mat.size()) {
		msgError(MISMATCH);exit(0);//underthinking
	}
	this->row = _row; this->col = _col;this->mat.resize(_row*_col);
	std::copy(_mat.begin(), _mat.end(), mat.begin());
};

template<typename T>
void Matrix<T>::set(int _row, int _col, T* _mat) {
	if(_row*_col <= 0) {
		msgError(NOPOSITIVE);exit(0);//underthinking
	}
	this->row = _row; this->col = _col;this->mat.resize(_row*_col);
	std::copy(_mat, _mat + (_row*_col), mat.begin());
};

template<typename T>
void Matrix<T>::set(int _row, int _col, T** _mat) {
	if(_row*_col <= 0) {
		msgError(NOPOSITIVE);exit(0);//underthinking
	}
	this->row = _row; this->col = _col;this->mat.resize(_row*_col);
	for(int i = 1; i <= _row; i++)
		for(int j = 1; j <= _col; j++)
			mat[getIndex(i,j)] = _mat[i-1][j-1]; 
};

template<typename T>
void Matrix<T>::set(int _row, int _col, std::istream&is) {
	if(_row*_col <= 0) {
		msgError(NOPOSITIVE);exit(0);//underthinking
	}
	this->row = _row; this->col = _col;this->mat.resize(_row*_col);
	for(int i = 1; i <= row*col; i++)
		is>>mat[i];
};

template<typename T>
void Matrix<T>::set(const Matrix<T>&p) {
	if(this == &p)return;
	this->row = p.row; this->col = p.col;this->mat.resize(p.row*p.col);
	std::copy(p.mat.begin(), p.mat.end(), mat.begin());
};

template<typename T>
void Matrix<T>::set(int _row, int _col) { //reshape the row and col but keep mat
	if(_row*_col != mat.size()) {
		msgError(MISMATCH);exit(0);//underthinking
	}
	this->row = _row; this->col = _col;
};

template<typename T>
T& Matrix<T>::operator()(int _i, int _j) { //element can be modified
	if(!isInternal(_i, _j)) {
		msgError(BEYOND);exit(0);//underthinking
	}
	return mat[getIndex(_i, _j)];
};

template<typename T>
void Matrix<T>::setRows(const Matrix<T>&_rows, int _i) {
	if(!isInternalRow(_i)) {
		msgError(BEYOND);exit(0);
	}
	if(_rows.row != 1 || _rows.col != col) {
		msgError(MISMATCH);exit(0);
	}
	for(int j = 1; j <= col; j++)
		mat[getIndex(_i, j)] = _rows.getElem(1, j);
};//replace rows-i

template<typename T>
void Matrix<T>::setCols(const Matrix<T>&_cols, int _j) {
	if(!isInternalCol(_j)) {
		msgError(BEYOND);exit(0);
	}
	if(_cols.row != row || _cols.col != 1) {
		msgError(MISMATCH);exit(0);
	}
	for(int i = 1; i <= row; i++)
		mat[getIndex(i, _j)] = _cols.getElem(i, 1);
};//replace cols-j

template<typename T>
void Matrix<T>::setRows(const Matrix<T>&_rows, int _rs, int _re) {
	if(!(isInternalRow(_rs) && isInternalRow(_re) && _re>=_rs)) {
		msgError(BEYOND);exit(0);
	}
	int nr = _re - _rs + 1;
	if(_rows.row != nr || _rows.col != col) {
		msgError(MISMATCH);exit(0);
	}
	for(int i = 1; i <= nr; i++)
		for(int j = 1; j <= col; j++)
			mat[getIndex(i+_rs-1, j)] = _rows.getElem(i, j);
};//replace rows-rs and rows-re

template<typename T>
void Matrix<T>::setCols(const Matrix<T>&_cols, int _cs, int _ce) {
	if(!(isInternalCol(_cs) && isInternalCol(_ce) && _ce>=_cs)) {
		msgError(BEYOND);exit(0);
	}
	int nc = _ce - _cs + 1;
	if(_cols.row != row || _cols.col != nc) {
		msgError(MISMATCH);exit(0);
	}
	for(int i = 1; i <= row; i++)
		for(int j = 1; j <= nc; j++)
			mat[getIndex(i, j+_cs-1)] = _cols.getElem(i, j);
};//replace cols-cs and cols-ce

template<typename T>
void Matrix<T>::setBlocks(const Matrix<T>&_blocks) {
	if(isLessSize(_blocks)) {
		msgError(BEYOND);exit(0);
	}
	if(isEqualSize(_blocks))return;
	for(int i = 1; i <= _blocks.row; i++)
		for(int j = 1; j <= _blocks.col; j++)
			mat[getIndex(i, j)] = _blocks.getElem(i, j);
};//replace block from(1,1)

template<typename T>
void Matrix<T>::setBlocks(const Matrix<T>&_blocks, int _rs, int _cs) {
	if(row > _blocks.row+_rs-1 || col > _blocks.col+_cs-1) {
		msgError(BEYOND);exit(0);
	}
	for(int i = 1; i <= _blocks.row; i++)
		for(int j = 1; j <= _blocks.col; j++)
			mat[getIndex(i+_rs-1, j+_cs-1)] = _blocks.getElem(i, j);
};//repalce block from (_rs, _cs)

template<typename T>
void Matrix<T>::setDiags(const Matrix<T>&_diags) {
	int n = MIN<int>(row, col);
	if(_diags.isVector() && n == _diags.getElem()) {
		for(int i = 1; i <= n; i++)
			mat[getIndex(i, i)] = _diags.mat[i-1];
		return;
	}
	if(isEqualSize(_diags)) {
		for(int i = 1; i <= n; i++)
			mat[getIndex(i, i)] = _diags.getElem(i, i);
		return;
	}
	msgError(MISMATCH);exit(0);
};//replace diags

//**get
template<typename T>
const Matrix<T> Matrix<T>::getRows(int _ri)const {
	Matrix<T> dst(1, col, T(0));;
	for(int j = 1; j <= col; j++)
		dst(1, j) = getElem(_ri, j);
	return dst;
};
template<typename T>
const Matrix<T> Matrix<T>::getCols(int _cj)const {
	Matrix<T> dst(row, 1, T(0));
	for(int i = 1; i <= row; i++)
		dst(i, 1) = getElem(i, _cj);
	return dst;
};
template<typename T>
const Matrix<T> Matrix<T>::getRows(int _rs, int _re)const {
	int nr = _re - _rs + 1;
	Matrix<T> dst(nr, col, T(0));
	for(int i = 1; i <= nr; i++)
		for(int j = 1; j <= col; j++)
			dst(i, j) = getElem(_rs+i-1, j);
	return dst;
};
template<typename T>
const Matrix<T> Matrix<T>::getCols(int _cs, int _ce)const {
	int nc = _ce - _cs + 1;
	Matrix<T> dst(row, nc, T(0));
	for(int i = 1; i <= row; i++)
		for(int j = 1; j <= nc; j++)
			dst(i, j) = getElem(i, _cs+j-1);
	return dst;
};
template<typename T>
const Matrix<T> Matrix<T>::getBlocks(int _rs, int _cs, int _re, int _ce)const {
	int rn = _re - _rs + 1, cn = _ce - _cs + 1;
	Matrix<T> dst(rn, cn, T(0));
	for(int i = 1; i <= rn; i++)
		for(int j =1; j <= cn; j++)
			dst(i, j) = getElem(_rs+i-1, _cs+j-1);
	return dst;
};
template<typename T>
const Matrix<T> Matrix<T>::getBlocks(int _re, int _ce)const {
	int rn = _re, cn = _ce;
	Matrix<T> dst(rn, cn, T(0));
	for(int i = 1; i <= rn; i++)
		for(int j =1; j <= cn; j++)
			dst(i, j) = getElem(i, j);
	return dst;
};
template<typename T>
const Matrix<T> Matrix<T>::getDiags()const {
	int nr = MIN<int>(row, col), nc = 1;
	Matrix<T> dst(nr, nc, T(0));
	for(int i = 1; i <= nr; i++)
		dst(i, 1) = getElem(i, i);
	return dst;
};
template<typename T>
const void Matrix<T>::getMat(std::vector<T>&dst)const {
	dst.resize(mat.size());
	copy(mat.begin(), mat.end(), dst.begin());
};
template<typename T>
const T Matrix<T>::getMax(int&dstr, int&dstc)const {
	int index = std::max_element(mat.begin(), mat.end()) - mat.begin();
	this->getIndexInv(index, dstr, dstc);
	return mat[index];
};
template<typename T>
const T Matrix<T>::getMin(int&dstr, int&dstc)const {
	int index = std::min_element(mat.begin(), mat.end()) - mat.begin();
	this->getIndexInv(index, dstr, dstc);
	return mat[index];
};
template<typename T>
const Matrix<T> Matrix<T>::getSimplestRows()const {
	Matrix<T> aug(*this);
	int dstr = 0; T elem;
	for(int i = 1; i <= col; i++) {
		//select main-element
		dstr = aug.getBlocks(i,i,row,i).absElems().getMaxRow() + i - 1;
		elem = aug.getElem(dstr, i);
		if(ABS<T>(elem) < EPS)continue;
		//exchange main-elem to (i,i)
		if(dstr != i)aug = aug.exRows(dstr, i);
		//divide main-elem
		for(int j = 1; j <= row; j++) {
			if(j == i)continue;	
			aug = aug.plRows(j, i, -aug.getElem(j, i)/elem);
		}
		aug = aug.mpRows(i, T(1)/elem);
	}
	return aug.downZero();
};
template<typename T>
const Matrix<T> Matrix<T>::getEchelonRows(int&flag)const {
	Matrix<T> aug(*this);flag = 1;
	int dstr = 0; T elem;
	for(int i = 1; i <= col; i++) {
		//select main-element
		dstr = aug.getBlocks(i,i,row,i).absElems().getMaxRow() + i - 1;
		elem = aug.getElem(dstr, i);
		if(ABS<T>(elem) < EPS)continue;
		//exchange main-elem to (i,i)
		if(dstr != i) {
			aug = aug.exRows(dstr, i); flag = -flag;
		}
		//divide main-elem
		for(int j = i+1; j <= row; j++)
			aug = aug.plRows(j, i, -aug.getElem(j, i)/elem);
	}
	return aug;
};
template<typename T>
const Matrix<T> Matrix<T>::getEchelonRows()const {
	Matrix<T> aug(*this);
	int dstr = 0; T elem;
	for(int i = 1; i <= col; i++) {
		//select main-element
		dstr = aug.getBlocks(i,i,row,i).absElems().getMaxRow() + i - 1;
		elem = aug.getElem(dstr, i);
		if(ABS<T>(elem) < EPS)continue;
		//exchange main-elem to (i,i)
		if(dstr != i)aug = aug.exRows(dstr, i);
		//divide main-elem
		for(int j = i+1; j <= row; j++)
			aug = aug.plRows(j, i, -aug.getElem(j, i)/elem);
	}
	return aug;
};
template<typename T>
const int Matrix<T>::getNonZeroRow()const {
	int count = 0;
	for(int i = 1; i <= row; i++)
		if(!this->getRows(i).isZero())count++;
	return count;
};
template<typename T>
const int Matrix<T>::getRank()const {
	return this->getSimplestRows().getNonZeroRow();
};
template<typename T>
const T Matrix<T>::getDet()const {
	if(!this->isSquare()) {
		msgError(WRONGCAL);exit(0);//underthinking
	}
	return this->getEchelonRows().getDiags().getProduct();
};
template<typename T>
const Matrix<T> Matrix<T>::getLPSubMatrix(int _p)const {
	if(!isInternalBlock(_p)) {
		msgError(BEYOND);exit(0);
	}
	return this->getBlocks(1, 1, _p, _p);
};//leading principle submatrix;p*p
template<typename T>
const T Matrix<T>::getLPMinor(int _p)const {
	return this->getLPSubMatrix(_p).getDet();
};//leading principle minor;det(p*p)
template<typename T>
const Matrix<T> Matrix<T>::getSumRows()const {
	Matrix<T> dst(row, 1, T(0));
	for(int i = 1; i <= row; i++)
		dst(i, 1) = this->getRows(i).getSum();
	return dst;
};//sum of each rows
template<typename T>
const Matrix<T> Matrix<T>::getSumCols()const {
	Matrix<T> dst(1, col, T(0));
	for(int j = 1; j <= col; j++)
		dst(1, j) = this->getCols(j).getSum();
	return dst;
};//sum of each cols
template<typename T>
const T Matrix<T>::getNormOne()const {
	return isVector()?(this->absElems().getSum())
		:(this->absElems().getSumCols().getMax());
};//norm-1 of vector or matrix 
template<typename T>
const T Matrix<T>::getNormTwo()const {
	return isVector()?(this->getRSSQ())
		:(sqrt(this->getVTV().eigQR().getMax()));
};//norm-2 of vector or matrix 
template<typename T>
const T Matrix<T>::getNormInf()const {
	return isVector()?(this->absElems().getMax())
		:(this->absElems().getSumRows().getMax());
};//norm-Inf of vector or matrix 

//calculation
//**special matrix
template<typename T>
Matrix<T> Matrix<T>::creDiags(const std::vector<T>&diags) {
	int _index = diags.size();
	Matrix<T> dst(_index, _index, T(0));
	for (int i = 1; i <= _index; i++)
		dst(i, i) = diags[i-1];
	return dst;
};
template<typename T>
Matrix<T> Matrix<T>::creDiags(const Matrix<T>&diags) {
	if(!diags.isVector()) {
		msgError(MISMATCH);exit(0);
	}
	int rc = diags.row*diags.col;
	Matrix<T> dst(rc, rc, T(0));
	for(int i = 1; i <= rc; i++)
		dst(i, i) = diags.mat[i-1];
	return dst;
};

template<typename T>
Matrix<T> Matrix<T>::creEyes(int _index) {
	Matrix<T> dst(_index, _index, T(0));
	for (int i = 1; i <= _index; i++)
		dst(i, i) = 1;
	return dst;
};

template<typename T>
Matrix<T> Matrix<T>::creRotX(T _rx) {//+ for inverse rotation, degree
	//axis moves
	//if you move vector not axis, transpose this mat
	T rx = _rx/180.0*PI;
	T data[9]={1, 0, 0,
			   0, cos(rx), sin(rx),
			   0, -sin(rx), cos(rx)};
	Matrix<T> dst(3,3,data);
	return dst;
};

template<typename T>
Matrix<T> Matrix<T>::creRotY(T _ry) {//+ for inverse rotation, degree
	//axis moves
	//if you move vector not axis, transpose this mat
	T ry = _ry/180.0*PI;
	T data[9]={cos(ry), 0, -sin(ry),
			   0, 1, 0,
			   sin(ry), 0, cos(ry)};
	Matrix<T> dst(3,3,data);
	return dst;
};

template<typename T>
Matrix<T> Matrix<T>::creRotZ(T _rz) {//+ for inverse rotation, degree
	//axis moves
	//if you move vector not axis, transpose this mat
	T rz = _rz/180.0*PI;
	T data[9]={cos(rz), sin(rz), 0,
			   -sin(rz), cos(rz), 0,
			   0, 0, 1};
	Matrix<T> dst(3,3,data);
	return dst;
};

template<typename T>
Matrix<T> Matrix<T>::creRot(T _r) {//+ for inverse rotation, degree, two-dimension
	//axis moves
	//if you move vector not axis, transpose this mat
	T r = _r/180.0*PI;
	T data[4]={cos(r), sin(r),
			   -sin(r), cos(r)};
	Matrix<T> dst(2, 2, data);
	return dst;
};

template<typename T>
Matrix<T> Matrix<T>::creRot(T _rx, T _ry, T _rz) {//three-dimension, Rx*Ry*Rz*coordinate
	//axis moves
	//if you move vector not axis, transpose this mat
	return creRotX(_rx)*creRotY(_ry)*creRotZ(_rz);
};

template<typename T>
Matrix<T> Matrix<T>::creRot(T _r, T _nx, T _ny, T _nz) {//arbitrary axis of rotation, unit normal vector(nx, ny, nz)
	//axis moves
	//if you move vector not axis, transpose this mat
	T r = _r/180.0*PI;
	T data[9]={_nx*_nx*(1-cos(r))+cos(r), _ny*_nx*(1-cos(r))+_nz*sin(r), _nz*_nx*(1-cos(r))-_ny*sin(r), 
			   _nx*_ny*(1-cos(r))-_nz*sin(r), _ny*_ny*(1-cos(r))+cos(r), _nz*_ny*(1-cos(r))+_nx*sin(r), 
			   _nx*_nz*(1-cos(r))+_ny*sin(r), _ny*_nz*(1-cos(r))-_nx*sin(r), _nz*_nz*(1-cos(r))+cos(r)};
	Matrix<T> dst(3,3,data);
	return dst;
};

template<typename T>
Matrix<T> Matrix<T>::creRot(T _r, FRAMEAXIS _a){//three-dimension, specified axis
	switch(_a) {
		case AXISX:return creRotX(_r);
		case AXISY:return creRotY(_r);
		case AXISZ:return creRotZ(_r);
		default:msgError(UNDEFINE);exit(0);
	}
};
template<typename T>
Matrix<T> Matrix<T>::creVander(const std::vector<T>&vander) {//Vandermonde matrix
	int rc = vander.size();
	if(rc == 0)return Matrix<T>();
	if(rc == 1)return Matrix<T>(1, 1, T(1));
	Matrix<T> dst(rc, rc, T(1));
	for(int i = 2; i <= rc; i++)
		for(int j =1; j <= rc; j++)
			dst(i, j) = pow(vander[j-1], i-1);
	return dst;
};
template<typename T>
Matrix<T> Matrix<T>::creStagger(const std::vector<T>&stagger, int _col) {//stagger vector in matrix by cols-directrion
	//(col-stagger.size()) x col
	int dc = stagger.size(), _row = _col - dc + 1;
	if(_row <= 0)return Matrix<T>();
	Matrix<T> dst(_row, _col, T(0));
	for(int i = 1; i <= _row; i++)
		for(int j = 1; j <= dc; j++)
			dst(i, i+j-1) = stagger[j-1];
	return dst;
};
template<typename T>
Matrix<T> Matrix<T>::creStaggers(const std::vector<T>&stagger, int _col) {//S'*S in fast way
	Matrix<T> S = Matrix<T>::creStagger(stagger, _col), ST = S.tp();
	Matrix<T> dst(ST.row, S.col, T(0));
	T sum = 0;
	for (int i = 1; i <= dst.row; i++) {
		for (int j = 1; j <= dst.col; j++) { 
			sum = 0;
			for (int k = 1; k <= stagger.size(); k++) {
				if(!ST.isInternalCol(i-k+1))continue;
				sum += ST.getElem(i, i-k+1)*S.getElem(i-k+1, j);
			}
			dst(i,j)=sum;//generate by the line
		}
	}
	return dst;
};

//**plus, minus, multiply
template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>&p)const {
	using namespace std;
	if(!this->isEqualSize(p)) {
		msgError(WRONGCAL);exit(0);//underthinking
	}
	Matrix<T> dst(row, col, T(0));
	transform(mat.begin(), mat.end(), p.mat.begin(), dst.mat.begin(), plus<T>());
	return dst;
};

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>&p)const {
	using namespace std;
	if(!this->isEqualSize(p)) {
		msgError(WRONGCAL);exit(0);//underthinking
	}
	Matrix<T> dst(row, col, T(0));
	transform(mat.begin(), mat.end(), p.mat.begin(), dst.mat.begin(), minus<T>());
	return dst;
};

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>&p)const {
	using namespace std;
	if(!this->isMatchMultiply(p)) {
		msgError(WRONGCAL);exit(0);//underthinking
	}
	Matrix<T> dst(row, p.col, T(0));
	T sum = 0;
	for (int i = 1; i <= row; i++) {
		for (int j = 1; j <= p.col; j++) { 
			sum = 0;
			for (int k = 1; k <= p.row; k++)
				sum += getElem(i, k)*p.getElem(k, j);
			dst(i,j)=sum;//generate by the line
		}
	}
	return dst;
};

template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>&p) {
	using namespace std;
	if(!this->isEqualSize(p)) {
		msgError(WRONGCAL);exit(0);//underthinking
	}
	transform(mat.begin(), mat.end(), p.mat.begin(), mat.begin(), plus<T>());
	return *this;
};

template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>&p) {
	using namespace std;
	if(!this->isEqualSize(p)) {
		msgError(WRONGCAL);exit(0);//underthinking
	}
	transform(mat.begin(), mat.end(), p.mat.begin(), mat.begin(), minus<T>());
	return *this;
};

//**normal transform
template<typename T>
Matrix<T> Matrix<T>::exRows(int _ri, int _rj)const {
	if(_ri == _rj)return *this;
	Matrix<T> dst(*this);
	for(int i = 1; i <= col; i++) {
		dst(_ri, i) = getElem(_rj, i);
		dst(_rj, i) = getElem(_ri, i);
	}
	return dst;
};//exchange row-i and row-j

template<typename T>
Matrix<T> Matrix<T>::exCols(int _ci, int _cj)const {
	if(_ci == _cj)return *this;
	Matrix<T> dst(*this);
	for(int i = 1; i <= row; i++) {
		dst(i, _ci) = getElem(i, _cj);
		dst(i, _cj) = getElem(i, _ci);
	}
	return dst;
};//exchange col-i and col-j

template<typename T>
Matrix<T> Matrix<T>::mpRows(int _ri, T _elem)const {
	Matrix<T> dst(*this);
	for(int i = 1; i <= col; i++)
		dst(_ri, i) = _elem*getElem(_ri, i);
	return dst;
};//multiply _elem for row-i

template<typename T>
Matrix<T> Matrix<T>::mpCols(int _ci, T _elem)const {
	Matrix<T> dst(*this);
	for(int i = 1; i <= row; i++)
		dst(i, _ci) = _elem*getElem(i, _ci);
	return dst;
};//multiply _elem for col-i

template<typename T>
Matrix<T> Matrix<T>::plRows(int _ri, int _rj, T _elem)const {
	Matrix<T> dst(*this);
	for(int i = 1; i <= col; i++)
		dst(_ri, i) = getElem(_ri, i) + _elem*getElem(_rj, i);
	return dst;
};//plus _elem*row-j for row-i

template<typename T>
Matrix<T> Matrix<T>::plCols(int _ci, int _cj, T _elem)const {
	Matrix<T> dst(*this);
	for(int i = 1; i <= col; i++)
		dst(i, _ci) = getElem(i, _ci) + _elem*getElem(i, _cj);
	return dst;
};//plus _elem*col-j for col-i

template<typename T>
Matrix<T> Matrix<T>::flipRows()const {
	Matrix<T> dst(row, col, T(0));
	for(int i = 1; i <= row; i++)
		for(int j = 1; j <= col; j++)
			dst(i, j) = getElem(row+1-i, j);
	return dst;
};//flip by row-oriention

template<typename T>
Matrix<T> Matrix<T>::flipCols()const {
	Matrix<T> dst(row, col, T(0));
	for(int i = 1; i <= row; i++)
		for(int j = 1; j <= col; j++)
			dst(i, j) = getElem(i, col+1-j);
	return dst;
};//flip by col-oriention

template<typename T>
Matrix<T> Matrix<T>::flipDiags()const {
	Matrix<T> dst(*this);
	int rc = MIN<T>(row, col);
	for(int i = 1; i <= rc; i++)
		for(int j = 1; j <= rc; j++)
			dst(i, j) = getElem(j, i);
	return dst;
};//flip by diag-oriention

//**sort
template<typename T>
Matrix<T> Matrix<T>::downZero()const {
	Matrix<T> dst(row, col, T(0));
	int count = 1;
	for(int i = 1; i <= row; i++) {
		if(this->getRows(i).isZero())continue;
		dst.setRows(this->getRows(i), count);
		count++;
	}
	return dst;
};//down the rows-zero to bottum

//*equation solution
template<typename T>
Matrix<T> Matrix<T>::solveHomo()const {
	int r = this->getRank(), count = 1;
	//zero-solution
	if(r == col)return Matrix<T>(col, 1, T(0));
	//not zero-solution
	Matrix<T> dst(col, col-r, T(0));
	Matrix<T> rsp = this->getSimplestRows();
	//**free x
	std::vector<int> fcs, nfcs, nfrs;
	for(int i = 1; i <= col; i++) {
		if(!rsp.getCols(i).isUnitOneVector(count)) {
			fcs.push_back(i);continue;
		}
		nfcs.push_back(i);nfrs.push_back(count);
		count++;
	}
	//**solution based on free x
	for(int i = 1; i <= col-r; i++) {
		dst(fcs[i-1] , i) = 1;
		for(int j = 1; j <= r; j++)
			dst(nfcs[j-1], i) = -rsp.getCols(fcs[i-1]).getElem(nfrs[j-1], 1);
	}
	return dst;
	
	
};//solve Homogenous equation

template<typename T>
Matrix<T> Matrix<T>::rpRows(const Matrix<T>&_rows, int _i)const {
	Matrix<T> dst(*this);
	dst.setRows(_rows, _i, _i+_rows.row-1);
	return dst;
};//row-i

template<typename T>
Matrix<T> Matrix<T>::rpCols(const Matrix<T>&_cols, int _j)const {
	Matrix<T> dst(*this);
	dst.setCols(_cols, _j, _j+_cols.col-1);
	return dst;
};//col-j

template<typename T>
Matrix<T> Matrix<T>::rpBlocks(const Matrix<T>&_blocks)const {
	if(isEqualSize(_blocks))return *this;
	//start (1, 1)
	Matrix<T> dst(*this);
	dst.setBlocks(_blocks);
	return dst;
};//block

template<typename T>
Matrix<T> Matrix<T>::rpBlocks(const Matrix<T>&_blocks, int _rs, int _cs)const {
	//start (_rs, _cs)
	Matrix<T> dst(*this);
	dst.setBlocks(_blocks, _rs, _cs);
	return dst;
};//block

template<typename T>
Matrix<T> Matrix<T>::rpDiags(const Matrix<T>&_diags) {
	Matrix<T> dst(*this);
	dst.setDiags(_diags);
	return dst;
};//replace diags

template<typename T>
Matrix<T> Matrix<T>::addRows(const Matrix<T>&_rows)const {
	int n = _rows.row;
	Matrix<T> dst(row+n, col, T(0));
	dst.setBlocks(*this);
	dst.setRows(_rows, row+1, row+n);
	return dst;
};//add rows-last

template<typename T>
Matrix<T> Matrix<T>::addCols(const Matrix<T>&_cols)const {
	int n = _cols.col;
	Matrix<T> dst(row, col+n, T(0));
	dst.setBlocks(*this);
	dst.setCols(_cols, col+1, col+n);
	return dst;
};//add cols-last

template<typename T>
Matrix<T> Matrix<T>::substractRows(int num)const {
	Matrix<T> dst(row-num, col, T(0));
	dst.setRows(this->getRows(1, row-num), 1, row-num);
	return dst;
};//substract rows-last-num

template<typename T>
Matrix<T> Matrix<T>::substractCols(int num)const {
	Matrix<T> dst(row, col-num, T(0));
	dst.setCols(this->getCols(1, col-num), 1, col-num);
	return dst;
};//substract cols-last-num

template<typename T>
Matrix<T> Matrix<T>::insertRows(const Matrix<T>&_rows, int _i)const {
	int n = _rows.row;
	Matrix<T> dst(row+n, col, T(0));
	dst.setBlocks(*this);
	dst.setRows(_rows, _i, _i+n-1);
	dst.setRows(this->getRows(_i, row), _i+n, row+n);
	return dst;
};//insert rows-i

template<typename T>
Matrix<T> Matrix<T>::insertCols(const Matrix<T>&_cols, int _j)const {
	int n = _cols.col;
	Matrix<T> dst(row, col+n, T(0));
	dst.setBlocks(*this);
	dst.setCols(_cols, _j, _j+n-1);
	dst.setCols(this->getCols(_j, col), _j+n, col+n);
	return dst;
};//insert cols-j

template<typename T>
Matrix<T> Matrix<T>::deleteRows(int _i)const {
	Matrix<T> dst(*this);
	dst.setRows(this->getRows(_i+1, row), _i, row-1);
	return dst.substractRows(1);
};//delete rows-i

template<typename T>
Matrix<T> Matrix<T>::deleteCols(int _j)const {
	Matrix<T> dst(*this);
	dst.setCols(this->getCols(_j+1, col), _j, col-1);
	return dst.substractCols(1);
};//delete cols-j

template<typename T>
Matrix<T> Matrix<T>::deleteRows(int _rs, int _re)const {
	Matrix<T> dst(*this);
	int n = _re - _rs + 1;
	dst.setRows(this->getRows(_re+1, row), _rs, row-n);
	return dst.substractRows(n);
};//delete rows-rs to rows-re

template<typename T>
Matrix<T> Matrix<T>::deleteCols(int _cs, int _ce)const {
	Matrix<T> dst(*this);
	int n = _ce - _cs + 1;
	dst.setCols(this->getCols(_ce+1, col), _cs, col-n);
	return dst.substractCols(n);
};//delete cols-cs to cols-ce

template<typename T>
Matrix<T> Matrix<T>::mergeRows(const Matrix<T>&p1, const Matrix<T>&p2) {
	int nr = p1.row + p2.row, nc = MAX<int>(p1.col, p2.col);
	Matrix<T> dst(nr, nc, T(0));
	//copy p1
	for(int i = 1; i <= p1.row; i++)
		for(int j = 1; j <= p1.col; j++)
			dst(i, j) = p1.getElem(i, j);
	//copy p2
	for(int i = 1; i <= p2.row; i++)
		for(int j = 1; j <= p2.col; j++)
			dst(i + p1.row, j) = p2.getElem(i, j);
	return dst;
};//merge by row-oriention

template<typename T>
Matrix<T> Matrix<T>::mergeCols(const Matrix<T>&p1, const Matrix<T>&p2) {
	int nr = MAX<int>(p1.row, p2.row), nc = p1.col + p2.col;
	Matrix<T> dst(nr, nc, T(0));
	//copy p1
	for(int i = 1; i <= p1.row; i++)
		for(int j = 1; j <= p1.col; j++)
			dst(i, j) = p1.getElem(i, j);
	//copy p2
	for(int i = 1; i <= p2.row; i++)
		for(int j = 1; j <= p2.col; j++)
			dst(i, j + p1.col) = p2.getElem(i, j);
	return dst;
};//merge by col-oriention

template<typename T>
Matrix<T> Matrix<T>::mergeBlocks(const Matrix<T>&p1, const Matrix<T>&p2) {
	int nr = p1.row + p2.row, nc = p1.col + p2.col;
	Matrix<T> dst(nr, nc, T(0));
	//copy p1
	for(int i = 1; i <= p1.row; i++)
		for(int j = 1; j <= p1.col; j++)
			dst(i, j) = p1.getElem(i, j);
	//copy p2
	for(int i = 1; i <= p2.row; i++)
		for(int j = 1; j <= p2.col; j++)
			dst(i + p1.row, j + p1.col) = p2.getElem(i, j);
	return dst;
};//merge by diag-oriention

//**matrix factorization
template<typename T>
void Matrix<T>::facLU(Matrix<T>&L, Matrix<T>&U)const {
	if(!this->isInv()) {
		msgError(WRONGCAL);exit(0);
	}
	//initial
	L.set(row, col, T(0)); U.set(row, col, T(0));
	//row-i-U and col-i-L
	T sumU, sumL;int n = row;
	for(int k = 1; k <= n; k++) {
		for(int j = k; j <= n; j++) {
			sumU = 0; sumL = 0;
			//Ukj
			for(int p = 1; p <= k-1; p++)
				sumU+=(L.getElem(k, p)*U.getElem(p, j));
			U(k, j) = getElem(k, j) - sumU;
			//Ljk
			if(j == k) { L(j, k) = T(1);continue; }
			for(int p = 1; p <= k-1; p++)
				sumL+=(L.getElem(j, p)*U.getElem(p, k));
			L(j, k) = (getElem(j, k) - sumL)/U.getElem(k,k);
 		}
	}
};
template<typename T>
void Matrix<T>::facQR(Matrix<T>&Q, Matrix<T>&R)const {
	if(!this->isInv()) {
		msgError(WRONGCAL);exit(0);
	}
	//Matrix-Q
	Q = this->orthSMT();
	//Matrix-R
	R = Q.tp()*(*this);
};
//**orthgonization
template<typename T>
Matrix<T> Matrix<T>::orthSMT()const {
	if(!this->isInv()) {
		msgError(WRONGCAL);exit(0);
	}
	Matrix<T> dst(*this), bi, ai, Bi;
	T coeff = 0;
	//Schimidt orthgonization
	for(int j = 2; j <= col; j++) {//b1 = a1
		ai = this->getCols(j); Bi.set(row, 1, T(0));
		for(int i = 1; i <= j-1; i++) {
			bi = dst.getCols(i);
			coeff = innerProduct(ai, bi)/bi.getSSQ();
			Bi += coeff*bi;
		}
		dst.setCols(ai-Bi, j);	
	}
	return dst.normCols();
};
//**normlization
template<typename T>
Matrix<T> Matrix<T>::norm()const {
	return (*this)*(1.0/this->getRSSQ());
};

template<typename T>
Matrix<T> Matrix<T>::normRows()const {
	Matrix<T> dst(row, col, T(0)), ri;
	for(int i = 1; i <= row; i++) {
		ri = this->getRows(i);
		dst = dst.rpRows(ri*(1.0/ri.getRSSQ()), i);
	}
	return dst;
};

template<typename T>
Matrix<T> Matrix<T>::normCols()const {
	Matrix<T> dst(row, col, T(0)), cj;
	for(int j = 1; j <= col; j++) {
		cj = this->getCols(j);
		dst = dst.rpCols(cj*(1.0/cj.getRSSQ()), j);
	}
	return dst;
};
//**eigenvalue and eigenvector
template<typename T>
Matrix<T> Matrix<T>::eigQR()const {
	if(!isSquare()) {
		msgError(WRONGCAL);exit(0);
	}
	if(isSingle())return Matrix<T>(1, 1, mat[0]);
	Matrix<T> Q, R, P1, P2(*this);
	do {
		P1 = P2;
		P1.facQR(Q, R);
		P2 = R*Q;
	}while(!P1.getDiags().isAppr(P2.getDiags(), 1e-6));
	return P2.getDiags();
};//eigenvalue -- QR algorithm

//**transpose
template<typename T>
Matrix<T> Matrix<T>::tp()const
{
	Matrix<T> dst(col, row, T(0));
	for (int i = 1; i <= col; i++)
		for (int j = 1; j <= row; j++)
			dst(i,j)=getElem(j, i);
	return dst;
};
//**inverse
template<typename T>
Matrix<T> Matrix<T>::invGJC()const {
	if(!this->isSquare()) {
		msgError(WRONGCAL);exit(0);//underthinking
	}
	//GJC: Gauss-Jordan by Column
	Matrix<T> aug(mergeCols(*this, creEyes(col)));
	int dstr = 0; T elem;
	for(int i = 1; i <= col; i++) {
		//select main-element
		dstr = aug.getBlocks(i,i,row,i).absElems().getMaxRow() + i - 1;
		elem = aug.getElem(dstr, i);
		if(ABS<T>(elem) < EPS) {
			msgError(NOEXIST);exit(0);
		}
		//exchange main-elem to (i,i)
		aug = aug.exRows(dstr, i);
		//divide main-elem
		for(int j = 1; j <= row; j++) {
			if(j == i)continue;	
			aug = aug.plRows(j, i, -aug.getElem(j, i)/elem);
		}
		aug = aug.mpRows(i, T(1)/elem);
	}
	//get the right-block which is inv-matrix
	return aug.getCols(col+1, aug.col);
};

//**abs
template<typename T>
Matrix<T> Matrix<T>::absElems()const {
	using namespace std;
	Matrix<T> dst(row, col, T(0));
	std::transform(mat.begin(), mat.end(), dst.mat.begin(), ptr_fun(ABS<T>));
	return dst;
};

//**inner_product
template<typename T>
T Matrix<T>::innerProduct(const Matrix<T>&p)const {
	if(!isEqualSize(p)) {
		msgError(WRONGCAL);exit(0);
	}
	T dst = std::inner_product(mat.begin(), mat.end(), p.mat.begin(), T(0));
	return dst;
};
template<typename T>
T Matrix<T>::innerProduct(const Matrix<T>&p1, const Matrix<T>&p2) {
	if(!p1.isEqualSize(p2)) {
		msgError(WRONGCAL);exit(0);
	}
	T dst = std::inner_product(p1.mat.begin(), p1.mat.end(), p2.mat.begin(), T(0));
	return dst;
};

//special class type
typedef Matrix<double> DOUBLEMAT;
typedef Matrix<int> INTMAT;


#endif //Matrix_H_INCLUDE
//~