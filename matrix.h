#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>

using namespace std;

template <class T>
class Matrix {
private:
	int nn;
	int mm;
	T **v;
public:
	Matrix(int n, int m);			// Zero-based array
	Matrix(int n, int m, const T &a);	//Initialize to constant
	Matrix(int n, int m, const T *a);	// Initialize to array
	Matrix & operator=(const Matrix &rhs);	//assignment
	inline T* operator[](const int i);	//subscripting: pointer to row i
	void print_matrix();
	inline int nrows() const;
	inline int ncols() const;
	~Matrix();
};

template <class T>
Matrix<T>::Matrix(int n, int m) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1;i<n;i++) v[i] = v[i-1] + m;
}

template <class T>
Matrix<T>::Matrix(int n, int m, const T &a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = a;
}

template <class T>
Matrix<T>::Matrix(int n, int m, const T *a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = *a++;
}

template <class T>
Matrix<T> & Matrix<T>::operator=(const Matrix<T> &rhs)
{
	if (this != &rhs) {
		int i,j,nel;
		if (nn != rhs.nn || mm != rhs.mm) {
			if (v != NULL) {
				delete[] (v[0]);
				delete[] (v);
			}
			nn=rhs.nn;
			mm=rhs.mm;
			v = nn>0 ? new T*[nn] : NULL;
			nel = mm*nn;
			if (v) v[0] = nel>0 ? new T[nel] : NULL;
			for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
		}
		for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
	}
	return *this;
}

template <class T>
inline T* Matrix<T>::operator[](const int i)	{return v[i];} //subscripting: pointer to row i

template <class T>
void Matrix<T>::print_matrix() {
	for (int i = 0; i < this->nrows();  ++i)  {      
		for (int j=0; j < this->ncols();  ++j) 		
			cout << (*this)[i][j] << "\t";
		cout << endl;
	}
}

template <class T>
inline int Matrix<T>::nrows() const
{
	return nn;
}

template <class T>
inline int Matrix<T>::ncols() const
{
	return mm;
}

template <class T>
Matrix<T>::~Matrix()
{
	if (v != NULL) {
		delete[] (v[0]);
		delete[] (v);
	}
}

#endif 