/* extra matrix iterface implementation
 * (C) 2011-2013 Felix Salfelder
 * License GPLv3 or later
 */
#ifndef MMATRIXEXTRA__
#define MMATRIXEXTRA__
#include "m_matrix.h"
#include <complex>

using std::complex;
/*-----------------------------------*/
template<>
inline BSMATRIX<double>::BSMATRIX(BSMATRIX<double>::REAL, const
		BSMATRIX<complex<double> >& m ) :
	_changed(new bool[m.size()+1]),
	_lownode(new unsigned[m.size()+1]),
	_space(new double[m._nzcount]),
	_rowptr(new double*[m._size+1]),
	_colptr(new double*[m._size+1]),
	_diaptr(new double*[m._size+1]),
	_nzcount(m._nzcount),
	_size(m._size),
	_zero(0.)
{
	memcpy(_changed,m._changed, (_size+1) * sizeof(bool) );
	memcpy(_lownode,m._lownode, (_size+1) * sizeof(unsigned) );

	_trash = abs (m._trash);
	_min_pivot = abs(m._min_pivot) ;
	_space = new double[_nzcount];
	for(unsigned i=0; i<_nzcount; i++ ){
		_space[i]=(m._space[i]).real();
	}

	for(unsigned i=0; i<_size+1; i++){
		_rowptr[i]= _space - intptr_t(m._space - m._rowptr[i]);
		_colptr[i]= _space - intptr_t(m._space - m._colptr[i]);
		_diaptr[i]= _space - intptr_t(m._space - m._diaptr[i]);
	}
}
/*-----------------------------------*/
template<>
inline BSMATRIX<double>::BSMATRIX(BSMATRIX<double>::IMAG, const
		BSMATRIX<complex<double> >& m ) :
	_changed(new bool[m.size()+1]),
	_lownode(new unsigned[m.size()+1]),
	_rowptr(new double*[m._size+1]),
	_colptr(new double*[m._size+1]),
	_diaptr(new double*[m._size+1]),
	_nzcount(m._nzcount),
	_size(m._size),
	_zero(0.)
{
	memcpy(_changed,m._changed, (_size+1) * sizeof(bool) );
	memcpy(_lownode,m._lownode, (_size+1) * sizeof(unsigned) );

	_trash = abs (m._trash);
	_min_pivot = abs(m._min_pivot) ;
	_space = new double[_nzcount];
	for(unsigned i=0; i<_nzcount; i++ ){
		_space[i]=(m._space[i]).imag();
	}

	for(unsigned i=0; i<_size+1; i++){
		_rowptr[i]= _space - intptr_t(m._space - m._rowptr[i]);
		_colptr[i]= _space - intptr_t(m._space - m._colptr[i]);
		_diaptr[i]= _space - intptr_t(m._space - m._diaptr[i]);
	}
}
/*-----------------------------------*/
template<>
inline BSMATRIX<double>::BSMATRIX(BSMATRIX<double>::SUM, const
		BSMATRIX<complex<double> >& m ) :
	_changed(new bool[m.size()+1]),
	_lownode(new unsigned[m.size()+1]),
	_space(new double[m._nzcount]),
	_rowptr(new double*[m._size+1]),
	_colptr(new double*[m._size+1]),
	_diaptr(new double*[m._size+1]),
	_nzcount(m._nzcount),
	_size(m._size),
	_zero(0.)
{
	memcpy(_changed,m._changed, (_size+1) * sizeof(bool) );
	memcpy(_lownode,m._lownode, (_size+1) * sizeof(unsigned) );

	_trash = abs (m._trash);
	_min_pivot = abs(m._min_pivot) ;
	_space = new double[_nzcount];
	for(unsigned i=0; i<_nzcount; i++ ){
		_space[i]= (m._space[i]).real() + (m._space[i]).imag() ;
	}

	for(unsigned i=0; i<_size+1; i++){
		_rowptr[i]= _space - intptr_t(m._space - m._rowptr[i]);
		_colptr[i]= _space - intptr_t(m._space - m._colptr[i]);
		_diaptr[i]= _space - intptr_t(m._space - m._diaptr[i]);
	}
}
/*-----------------------------------*/
template<class T>
inline T* BSMATRIX<T>::row(T* x, unsigned n)
{
	for(unsigned i=0; i<=_size; ++i)
	  x[i]=s(n,i);	
	return x;
}
/*-----------------------------------*/
template<class T>
inline T* BSMATRIX<T>::col(T* x, unsigned n)
{
	for(unsigned i=0; i<=_size; ++i)
	  x[i]=s(i,n);	
	return x;
}
/*--------------------------------------------------------------------------*/
template<class T>
inline T*  BSMATRIX<T>::rmul(T* b, const T* x)const
{ // probably possible to exploit sparse structure here
	assert(b!=x);
	for(unsigned i=1; i<=_size; ++i){
		b[i] = 0;
		for(unsigned j=1; j<=_size; ++j){
			b[i] += s(i,j) * x[j];
		}
	}
	return b;
}
/*--------------------------------------------------------------------------*/
// template<class T>
// inline BSMATRIX<T>& BSMATRIX<T>::operator-(const BSMATRIX<T>& s)
// {
// 	for(unsigned i=0; i<_nzcount; i++){
// 		_space[i]=-s._space[i];
// 	}
// 	return *this;
// }
/*--------------------------------------------------------------------------*/
// only works with similar layout.
template<class T>
inline BSMATRIX<T>& BSMATRIX<T>::operator=(const BSMATRIX<T>& s)
{
	for(unsigned i=0; i<_nzcount; i++){
		_space[i]=s._space[i];
	}
	return *this;
}
/*--------------------------------------------------------------------------*/
template<class T>
inline BSMATRIX<T>& BSMATRIX<T>::operator+=(const BSMATRIX<T>& s)
{
	for(unsigned i=0; i<_nzcount; i++){
		_space[i]+=s._space[i];
	}
	return *this;
}
/*--------------------------------------------------------------------------*/
template<class T>
inline BSMATRIX<T>& BSMATRIX<T>::operator*=(const T& s)
{
	for(unsigned i=0; i<_nzcount; i++){
		_space[i]*=s;
	}
	return *this;
}
/*--------------------------------------------------------------------------*/
template<class T>
inline BSMATRIX<T> operator*(const T s, const BSMATRIX<T>& m)
{
	BSMATRIX<T> M(m);
	return M *= s;
}
/*--------------------------------------------------------------------------*/
template<class T>
inline BSMATRIX<T> operator-(BSMATRIX<T> m)
{
	m *= -1;
	return m;
}
/*--------------------------------------------------------------------------*/
template<class T>
void BSMATRIX<T>::augment(T* r, T* /*c, notyet*/)
{
	_size++;
	T** colptr = new T*[_size+1];
	T** rowptr = new T*[_size+1];
	T** diaptr = new T*[_size+1];
	bool* changed = new bool[_size+1];
	unsigned* lownode = new unsigned[_size+1];

	memcpy(colptr,_colptr, (_size) * sizeof(T*) );
	memcpy(diaptr,_diaptr, (_size) * sizeof(T*) );
	memcpy(rowptr,_rowptr, (_size) * sizeof(T*) );
	memcpy(lownode,_lownode, (_size) * sizeof(unsigned) );
	memcpy(changed,_changed, (_size) * sizeof(bool) );

	delete[] _colptr; _colptr=colptr; colptr[_size] = r-1;
	delete[] _diaptr; _diaptr=diaptr; diaptr[_size] = r+_size-1;
	delete[] _rowptr; _rowptr=rowptr; rowptr[_size] = r+2*_size-1; // row is reversed
	delete[]_changed;_changed=changed;changed[_size]= 1;
	delete[]_lownode;_lownode=lownode;lownode[_size]= 1;
}
/*--------------------------------------------------------------------------*/
template<class T>
void BSMATRIX<T>::deaugment()
{
	_size--;
}
#endif
