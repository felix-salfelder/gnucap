/*                                -*- C++ -*-
 * Copyright (C) 2001 Albert Davis
 * Author: Albert Davis <aldavis@gnu.org>
 *
 * This file is part of "Gnucap", the Gnu Circuit Analysis Package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 *------------------------------------------------------------------
 * Sparse matrix package
 * Bump and spike - bordered block diagonal pattern
 * -------------------------------------------------
 * To use it (simple method) .....
 * 1. Declare it ...
 *	BSMATRIX<double> a;
 * 2. Then tell it what slots to allocate, in a loop ...
 *	for (all pairs i,j that might exist)
 *	   a.iwant(i,j);
 * 3. Then do the actual allocation ...
 *	a.allocate();
 * 4. If you want, you can add an offset to the diagonal ...
 *	a.dezero(gmin);
 * 5. You probably should set a minimum pivot value ...
 *	a.set_min_pivot(gmin);
 * 6. Fill in the matrix ...
 *	for (all pairs i,j you want to fill in)
 *	   a.m(i,j) += data;
 * 7. Then do the LU decomposition ...
 *	a.lu_decomp();
 * 8. Then get the solution by applying the right side vector (v)
 * and doing the fwd and back substitution.
 *	a.fbsub(v);
 * The solution is now in v.
 * -------------------------------------------------
 * To reuse it .....
 * Get rid of the old data ...
 *	a.zero().
 * Restore the offset, if you want ...
 *	a.dezero(gmin);
 * Then fill and solve as before. (steps 6-8)
 * -------------------------------------------------
 * In the above case, LU replaced the matrix and the solution replaced
 * the right side, losing the data.
 * To keep the matrix, and keep the right side ...
 * 1. Declare two matrices ...
 *	BSMATRIX<double> a;
 *	BSMATRIX<double> lu;
 * 2. as before
 * 2a. Say the lu has the same structure as a.
 *	lu.clone(a);
 * 3-5. Allocate both.
 * 6. Fill in "a" only.
 * 7. Do the LU decomposition, keeping a, with result in lu ...
 *	lu.lu_decomp(a, false);
 * 8. Do the f&B sub, keeping b (the right side), with result in x ...
 *	lu.fbsub(x, b);
 * -------------------------------------------------
 * To make a change to the matrix and re-solve ...
 * Apply the change to a ...
 *	for (all changes you want to make)
 *	   a.m(i,j) += delta;  a.set_changed(i).set_changed(j)
 * Then solve again .. step 7 above for a full solution,
 * or for an update ...
 *	lu.lu_decomp(a, true);
 * Then do the same step 8 as above.
 * -------------------------------------------------
 * some other functions that might be useful ....
 * reinit(newsize) -- change the size (and lose the contents and structure)
 * size()const -- return the size (# of rows & cols)
 * density()const -- return the matrix density, as a number between 0 and 1
 * sparsity()const -- 1-density
 * -------------------------------------------------
 * individual element access ...
 * 5 access functions are provided.
 * All return lvalues (references to the actual entry).
 * All take 2 args, row and column.  They differ in speed and generality.
 * Since they are used in equations, they all have 1 letter names.
 * s(r,c) -- "safe" -- most general case, with bounds checking
 *	if completely out of bounds, returns trash
 *	if in the matrix but out of band, returns a reference to zero
 *	Changing it will corrupt the matrix.
 *	All others have no bounds checking.
 * m(r,c) -- "matrix" -- known to be within the allocated band
 * r,c) -- "upper" -- known to be in the upper triangle. (c>=r)
 * l(r,c) -- "lower" -- known to be in the lower triangle. (r>=c)
 * d(r,c) -- "diagonal" -- known to be on the diagonal. (r==c)
 * Using s() will always work, but will be slower than the other methods.
 * You should use the most restricted function that you know is correct.
 * -------------------------------------------------
 * Other notes ...
 * The numbering starts at 1 (Fortran style).
 * When using "s" access, it is ok to load row and column 0, then ignore it.
 *   This may simplify load functions, at the expense of run speed.
 * "s" will let you change the value of zero,
 *   but you will find out about it later.
 */
//testing=script 2008.09.19
#ifndef M_MATRIX_H
#define M_MATRIX_H
/*--------------------------------------------------------------------------*/
#include "l_stlextra.h"
#include "io_.h"
#include <iostream>
#include <fstream>
#include <set>
/*--------------------------------------------------------------------------*/
enum needed_t  {nNO=0, nYES, nFILL};
/*--------------------------------------------------------------------------*/
class OMSTREAM;
template <class T>
class BSMATRIX {
  // hack.
  friend class BSMATRIX<double>; 
  friend class BSMATRIX<std::complex<double> >;
private:
  mutable bool* _changed;  // flag: this node changed value
  unsigned*	_lownode;  // lowest node connecting to this one
  T*	_space;		// ptr to actual memory space used
  T**	_rowptr;	// ptrs to col 0 of every row
  T**	_colptr;	// ptrs to row 0 of every col
  T**	_diaptr;	// ptrs to diagonal
  unsigned _nzcount;	// count of non-zero elements
  unsigned _size;	// # of rows and columns
  T	_zero;		// always 0 but not const
  T	_trash;		// depository for row and col 0, write only
  T	_min_pivot;	// minimum pivot value
  unsigned _bandwidth;  // (0=diagonal)
  std::vector<std::set<unsigned> > _adj;
  std::vector<std::set<unsigned> > _adj_ordered;
public:
  enum REAL {_REAL};
  enum IMAG {_IMAG};
  enum SUM {_SUM};
  BSMATRIX(const BSMATRIX<T>&);
private:
  explicit	BSMATRIX(REAL, const BSMATRIX<std::complex<double> >&);
  explicit	BSMATRIX(IMAG, const BSMATRIX<std::complex<double> >&);
  explicit	BSMATRIX(SUM, const BSMATRIX<std::complex<double> >&);
  void		uninit();
  void		init(unsigned s=0);
  T&		subtract_dot_product(uint_t r, uint_t c, uint_t d);
  T&		subtract_dot_product(uint_t r, uint_t c, uint_t d, const T& in);
  unsigned	lownode(unsigned i)const {return _lownode[i];}
  bool		is_changed(unsigned n)const {return _changed[n];}
  void		set_changed(uint_t n, bool x = true)const {_changed[n] = x;}
public:
  explicit	BSMATRIX(unsigned ss=0);
  		~BSMATRIX()		{uninit();}
  void		reinit(unsigned ss=0)	{uninit(); init(ss);}
  //void	clone(const BSMATRIX<T>&);
  BSMATRIX<T>*  copy()const // also copy contents
  { return new BSMATRIX(*this); }
  BSMATRIX<double> real()const
  { return( BSMATRIX<double>(BSMATRIX< double >::_REAL, *this) ) ; }
  BSMATRIX<double> imag()const
  { return( BSMATRIX<double>(BSMATRIX< double >::_IMAG, *this) ) ; }
  BSMATRIX<double> sum()const
  { return( BSMATRIX<double>(BSMATRIX< double >::_SUM, *this) ) ; }
  T* row(T*, unsigned);
  T* col(T*, unsigned);
  void		iwant(unsigned, unsigned);
  void		iwant(unsigned*); 
private:
  void          i_do_want(unsigned, unsigned);
  int* _rcm;
  void rcm(); // used by allocate.
public:
  void		unallocate();
  void		allocate(unsigned*nm=NULL);
  void		reallocate()		{unallocate(); allocate();}
  void		set_min_pivot(double x)	{_min_pivot = x;}
  void		zero();
  void		dezero(T& o);
  unsigned	size()const		{return _size;}
  double 	density();
  inline T*    rmul(T* b, const T* x)const; // returns A*x
  T 	d(unsigned r, unsigned  )const	{return *(_diaptr[r]);}
  T     s(unsigned r, unsigned c)const; 
  needed_t n(unsigned row, unsigned col)const;
private:
  T 	u(unsigned r, unsigned c)const	{ return _colptr[c][r];}
  T 	l(unsigned r, unsigned c)const	{ return *(_rowptr[r]-c);}
  T&	d(unsigned r, unsigned c);
  T&	u(unsigned r, unsigned c);
  T&	l(unsigned r, unsigned c);
  T&	m(unsigned r, unsigned c);
  T&	s(unsigned r, unsigned c);
public:
  template <class S,class X>
  friend S& operator<< ( S &o, const BSMATRIX<X>& m);

  BSMATRIX& operator*=( const T& s);
  BSMATRIX& operator+=( const BSMATRIX& s);
  BSMATRIX& operator=( const BSMATRIX& s);
//  BSMATRIX& operator-( BSMATRIX& s);
//  template <class X>
//  friend OMSTREAM& operator<< ( OMSTREAM &o, const BSMATRIX<X>& m);
  void		load_diagonal_point(uint_t i, T value);
  void		load_point(uint_t i, uint_t j, T value);
  void		load_couple(uint_t i, uint_t j, T value);
  void		load_symmetric(uint_t i, uint_t j, T value);
  void		load_asymmetric(uint_t r1, uint_t r2, uint_t c1, uint_t c2, T value);
  
  void		lu_decomp(const BSMATRIX<T>&, bool do_partial=false);
  void		lu_decomp();
  void		fbsub(T* v) const;
  void		fbsub(std::complex<T>* v) const;
  void		fbsub(T* x, const T* b, T* c ) const;
  T*		fbsub(T* x, const T* b ) const { fbsub(x,b,x); return x; }
  void		fbsubt(T* v) const;
public:
  void sink_forward(unsigned* nm);
  void sink_reverse(unsigned* nm);
public: //eigenvalue stuff
  void augment(T*, T*row=0);
  void deaugment();
};
/*--------------------------------------------------------------------------*/
// private implementations
/*--------------------------------------------------------------------------*/
template <class T>
void BSMATRIX<T>::uninit()
{
  unallocate();
  delete [] _lownode;
  _lownode = NULL;
  delete [] _changed;
  _changed = NULL;
}
/*--------------------------------------------------------------------------*/
template <class T>
void BSMATRIX<T>::init(unsigned ss)
{
  assert(!_lownode);
  assert(!_colptr);
  assert(!_rowptr);
  assert(!_diaptr);
  assert(!_space);

  assert(_zero == 0.);
  _min_pivot = _trash = 0.;
  _nzcount = 0;
  _size = ss;
  _lownode = new unsigned[size()+1];
  assert(_lownode);
  for (unsigned ii = 0;  ii <= size();  ++ii) {
    _lownode[ii] = ii;
  }
  _changed = new bool[size()+1];
  assert(_changed);
  for (unsigned ii = 0;  ii <= size();  ++ii) {
    set_changed(ii, false);
  }

  _adj.clear();
  _adj.resize(size()+1);
}
/*--------------------------------------------------------------------------*/
template <class T>
T& BSMATRIX<T>::subtract_dot_product(uint_t rr, uint_t cc, uint_t dd)
{
  assert(_lownode);
  uint_t kk = std::max(_lownode[rr], _lownode[cc]);
  assert(dd>=kk);
  unsigned len = unsigned(dd - kk);
  T& dot = m(rr, cc);
  if (len > 0) {
    T* row_ = &(l(rr,kk));
    T* col_ = &(u(kk,cc));
    /* for (ii = kk;   ii < dd;   ++ii) */
    for (unsigned ii = 0;   ii < len;   ++ii) {
      dot -= *(row_-ii) * col_[ii];
    }
  }else{
  }
  return dot;
}
/*--------------------------------------------------------------------------*/
template <class T>
T& BSMATRIX<T>::subtract_dot_product(uint_t rr, uint_t cc, uint_t dd, const T& in)
{
  assert(_lownode);
  unsigned kk = std::max(_lownode[rr], _lownode[cc]);
  int len = int(dd) - int(kk);
  T& dot = m(rr, cc);
  dot = in;
  if (len > 0) {
    T* row_ = &(l(rr,kk));
    T* col_ = &(u(kk,cc));
    /* for (ii = kk;   ii < dd;   ++ii) */
    for (int ii = 0;   ii < len;   ++ii) {
      dot -= *(row_-ii) * col_[ii];
    }
  }else{
  }
  return dot;
}
/*--------------------------------------------------------------------------*/
// public implementations
/*--------------------------------------------------------------------------*/
template <class T>
BSMATRIX<T>::BSMATRIX(unsigned ss)
  :_changed(NULL),
   _lownode(NULL),
   _space(NULL),
   _rowptr(NULL),
   _colptr(NULL),
   _diaptr(NULL),
   _nzcount(0),
   _size(ss),
   _zero(0.),
   _trash(0.),
   _min_pivot(0.),
   _bandwidth(0)
{
  init(ss);
}
/*--------------------------------------------------------------------------*/
#if 0
/* clone: copy to self the structure of another BSMATRIX
 * this does not copy the values stored in the matrix
 */
template <class T>
void BSMATRIX<T>::clone(const BSMATRIX<T> & aa)
{untested();
  reinit(aa.size());
  for (int ii = 0;  ii <= size();  ++ii) {untested();
    _lownode[ii] = aa.lownode(ii);
  }
}
#endif
/*--------------------------------------------------------------------------*/
/* iwant: indicate that "iwant" to allocate this spot in the matrix
 */
template <class T>
void BSMATRIX<T>::iwant(unsigned node1, unsigned node2)
{
  trace2("BSMATRIX::iwant", node1, node2);
  assert(_lownode);
  assert(node1 <= size());
  assert(node2 <= size());

//   _adj.resize( std::max(std::max(node1, node2), unsigned(_adj.size())));
  assert(node1 < _adj.size());
  assert(node2 < _adj.size());
  _adj[node1].insert(node2);
  _adj[node2].insert(node1);

}
/*--------------------------------------*/
template <class T>
void BSMATRIX<T>::iwant(unsigned* nm)
{
  _adj_ordered.clear();
  _adj_ordered.resize(_adj.size());
  for( unsigned j=0; j<size(); ++j ){
    for(std::set<unsigned>::const_iterator ii=_adj[j].begin();
        ii != _adj[j].end(); ++ii )
    {
      i_do_want(nm[j],nm[*ii]);
      _adj_ordered[nm[j]].insert(nm[*ii]);
      _adj_ordered[nm[*ii]].insert(nm[j]);
    }
  }
}
/*--------------------------------------*/
// i really want the nodes (former iwant)
template <class T>
void BSMATRIX<T>::i_do_want(unsigned node1, unsigned node2)
{
  if (node1 <= 0  ||  node2 <= 0) {
    // node 0 is ground, and doesn't count as a connection
    // negative is invalid, not used but still may be in a node list
  }else if (node1 < _lownode[node2]) {
    _lownode[node2]=node1;
  }else if (node2 < _lownode[node1]) {
    _lownode[node1]=node2;
  }else{
  }
}
/*--------------------------------------------------------------------------*/
template <class T>
void BSMATRIX<T>::unallocate()
{
  assert (_zero == 0.);
  delete [] _rowptr;
  delete [] _colptr;
  delete [] _diaptr;
  delete [] _space;

  _rowptr = _colptr = _diaptr = NULL;
  _space = NULL;
}
/*--------------------------------------------------------------------------*/
/* allocate: really get the space to work
 * hmm maybe pass nm here?
 */
template <class T>
void BSMATRIX<T>::allocate(unsigned* nm)
{
  if(nm)incomplete();
  assert(_lownode);
  assert(!_colptr);
  assert(!_rowptr);
  assert(!_diaptr);
  assert(!_space);

  _nzcount = 0;
  for (unsigned ii = 0;   ii <= size();   ++ii) {
    assert (ii >= _lownode[ii]);
    _nzcount += 2 * unsigned(ii - _lownode[ii]) + 1;
  }

  _colptr = new T*[size()+1];
  _rowptr = new T*[size()+1];
  _diaptr = new T*[size()+1];
  trace1("BSMATRIX::allocate", _nzcount);
  _space  = new T[_nzcount];

  assert(_colptr);
  assert(_rowptr);
  assert(_diaptr);

  zero();

  {
    T* point = _space;
    for (unsigned ii = 0; ii <= size(); ++ii) {
      _colptr[ii] = point - (int)_lownode[ii];
      _rowptr[ii] = _colptr[ii] + 2*ii;
      _diaptr[ii] = _colptr[ii] + ii;
      point += 2 * ((int)ii - (int)_lownode[ii]) + 1;
    }
  }
}
/*--------------------------------------------------------------------------*/
/* zero: wipe the whole array
 */
template <class T>
void BSMATRIX<T>::zero()
{
  assert(_space);
  assert(_zero == 0.);
  _trash = 0.;
  std::fill_n(_space, _nzcount, 0.);
}
/*--------------------------------------------------------------------------*/
/* dezero: make sure(?) the diagonal is non-zero
 */
template <class T>
void BSMATRIX<T>::dezero(T& offset)
{
  for (unsigned ii = 1;  ii <= size();  ++ii) {
    d(ii,ii) += offset;
  }
}
/*--------------------------------------------------------------------------*/
template <class T>
double BSMATRIX<T>::density()
{
  if (size() > 0) {
    assert(_lownode);
    _nzcount = 0;
    for (unsigned ii = 0;   ii <= size();   ++ii) {
      _nzcount += 2 * (ii - _lownode[ii]) + 1;
    }
    return static_cast<double>(_nzcount-1)/(static_cast<double>(size())*size());
  }else{
    return 0;
  }
}
/*--------------------------------------------------------------------------*/
/* d: fast matrix entry access, lvalue
 * It is known that the entry is valid and on the diagonal
 */
template <class T>
T& BSMATRIX<T>::d(unsigned r, unsigned c)
{
  assert(_diaptr);
  assert(r == c); USE(c);
  assert(r);
  assert(r <= _size);

  return *(_diaptr[r]);
}
/*--------------------------------------------------------------------------*/
/* u: fast matrix entry access
 * It is known that the entry is valid and in the upper triangle (or diagonal)
 */
template <class T>
T& BSMATRIX<T>::u(unsigned r, unsigned c)
{
  assert(_colptr);
  assert(_lownode);
  assert(r);
  assert(r <= c);
  //assert(c <= _size);
  assert(1 <= _lownode[c]);
  assert(_lownode[c] <= r);

  return _colptr[c][r];
}
/*--------------------------------------------------------------------------*/
/* l: fast matrix entry access
 * It is known that the entry is valid and in the lower triangle not diagonal
 */
template <class T>
T& BSMATRIX<T>::l(unsigned r, unsigned c)
{
  assert(_rowptr);
  assert(_lownode);
  assert(0 < c);
  assert(c < r);
  assert(r <= _size);
  assert(1 <= _lownode[r]);
  assert(_lownode[r] <= c);

  return *(_rowptr[r]-c);
}
/*--------------------------------------------------------------------------*/
/* m: semi-fast matrix entry access
 * It is known that the entry refers to a valid location
 * but it is not known whether lower, upper, or diagonal
 */
template <class T>
T& BSMATRIX<T>::m(unsigned r, unsigned c)
{
  if (c>=r) {
    return u(r,c);
  }else{
    return l(r,c);
  }
}
/*--------------------------------------------------------------------------*/
/* s: general matrix entry access.
 * It is known that the location is strictly in bounds,
 *   but it is not known whether the location actually exists.
 * If access is attempted to a non-allocated location, 
 *   it returns a reference to a shared zero variable.
 *   Writing to this zero is not prohibited,
 *   but will corrupt the matrix in a known and testable way.
 * If access is attempted to row 0 or column 0,
 *   it returns a reference to a shared trash variable.
 *   Writing to trash is allowed and encouraged,
 *   but reading it gives a number not useful for anything.
 */
template <class T>
T BSMATRIX<T>::s(unsigned row, unsigned col)const
{
  assert(_lownode);
  // assert(0 <= col);
  assert(col <= size());
  // assert(0 <= row);
  assert(row <= size());
  assert(_zero == 0.);

  if (col == row) {
    return d(row,col);
  }else if (col > row) {	/* above the diagonal */
    if (row == 0) {
      return _trash;
    }else if (row < _lownode[col]) {
      return _zero;
    }else{
      return (u(row,col));
    }
  }else{			/* below the diagonal */
    assert(col < row);
    if (col == 0) {
      return _trash;
    }else if (col < _lownode[row]) {
      return _zero;
    }else{
      return (l(row,col));
    }
  }
  unreachable();
}
/*--------------------------------------------------------------------------*/
// the entry is allocated (non-zero)
template <class T>
needed_t BSMATRIX<T>::n(unsigned row, unsigned col)const
{
  assert(_lownode);
  // assert(0 <= col);
  assert(col <= size());
  // assert(0 <= row);
  assert(row <= size());
  assert(_zero == 0.);

  if (_adj_ordered[row].find(col) != _adj_ordered[row].end()){
    return nYES;
  }

  if (col == row) {
    return nYES;
  }else if (col > row) {	/* above the diagonal */
    if (row == 0) {
      return nNO;
    }else if (row < _lownode[col]) {
      return nNO;
    }else{
      return nFILL;
    }
  }else{			/* below the diagonal */
    assert(col < row);
    if (col == 0) {
      return nNO;
    }else if (col < _lownode[row]) {
      return nNO;
    }else{
      return nFILL; 
    }
  }
  unreachable();
}
/*--------------------------------------------------------------------------*/
// rw access.
template <class T>
T& BSMATRIX<T>::s(unsigned row, unsigned col)
{
  assert(_lownode);
// assert(0 <= col);
  assert(col <= size());
//  assert(0 <= row);
  assert(row <= size());
  assert(_zero == 0.);

  if (col == row) {
    return d(row,col);
  }else if (col > row) {	/* above the diagonal */
    if (row == 0) {
      return _trash;
    }else if (row < _lownode[col]) {
      return _zero;
    }else{
      return (u(row,col));
    }
  }else{			/* below the diagonal */
    assert(col < row);
    if (col == 0) {
      return _trash;
    }else if (col < _lownode[row]) {
      return _zero;
    }else{
      return (l(row,col));
    }
  }
  unreachable();
}
/*--------------------------------------------------------------------------*/
template <class T>
void BSMATRIX<T>::load_point(uint_t i, uint_t j, T value)
{
  if (i > 0 && j > 0) {
    set_changed(j);
    set_changed(i);
    m(i,j) += value;
  }else{
  }
}
/*--------------------------------------------------------------------------*/
// load_point(i, i, value);
template <class T>
void BSMATRIX<T>::load_diagonal_point(uint_t i, T value)
{
  if (i > 0) {
    set_changed(i);
    d(i,i) += value;
  }else{untested();
  }
}
/*--------------------------------------------------------------------------*/
// load_point(i, j, -value);
// load_point(j, i, -value);
template <class T>
void BSMATRIX<T>::load_couple(uint_t i, uint_t j, T value)
{
  if (j > 0) {
    set_changed(j);
    if (i > 0) {
      set_changed(i);
      m(i,j) -= value;
      m(j,i) -= value;
    }else{
    }
  }else{untested();
  }
}
/*--------------------------------------------------------------------------*/
// load_point(i, i, value); or load_diagonal_point(i, value);
// load_point(j, j, value); or load_diagonal_point(j, value);
// load_point(i, j, -value);
// load_point(j, i, -value);
template <class T>
void BSMATRIX<T>::load_symmetric(uint_t i, uint_t j, T value)
{
  if (j > 0) {
    set_changed(j);
    d(j,j) += value;
    if (i > 0) {
      set_changed(i);
      d(i,i) += value;
      m(i,j) -= value;
      m(j,i) -= value;
    }else{
    }
  }else if (i > 0) {
    set_changed(i);
    d(i,i) += value;
  }else{
  }
}
/*--------------------------------------------------------------------------*/
// load_point(r1, c1, value);
// load_point(r2, c2, value);
// load_point(r1, c2, -value);
// load_point(r2, c1, -value);
template <class T>
void BSMATRIX<T>::load_asymmetric(uint_t r1,uint_t r2,uint_t c1,uint_t c2,T value)
{
  set_changed(c1);
  set_changed(c2);
  if (r1 > 0) {
    set_changed(r1);
    if (c1 > 0) {
      m(r1,c1) += value;
    }else{
    }
    if (c2 > 0) {
      m(r1,c2) -= value;
    }else{
    }
  }else{
  }
  if (r2 > 0) {
    set_changed(r2);
    if (c1 > 0) {
      m(r2,c1) -= value;
    }else{
    }
    if (c2 > 0) {
      m(r2,c2) += value;
    }else{
    }
  }else{
  }
}
/*--------------------------------------------------------------------------*/
template <class T>
void BSMATRIX<T>::lu_decomp(const BSMATRIX<T>& aa, bool do_partial)
{
  unsigned prop = 0;   /* change propagation indicator */
  assert(_lownode);
  assert(aa._lownode);
  assert(aa.size() == size());
  for (unsigned mm = 1;   mm <= size();   ++mm) {
    assert(aa.lownode(mm) == _lownode[mm]);
    unsigned bn = _lownode[mm];
    if (!do_partial  ||  aa.is_changed(mm)  ||  bn <= prop) {
      aa.set_changed(mm, false);
      prop = mm;
      if (bn < mm) {
	prop = mm;
	u(bn,mm) = aa.u(bn,mm) / d(bn,bn);
	for (unsigned ii = bn+1;  ii<mm;  ii++) {
	  /* u(ii,mm) = (aa.u(ii,mm) - dot(ii,mm,ii)) / d(ii,ii); */
	  subtract_dot_product(ii,mm,ii,aa.u(ii,mm)) /= d(ii,ii);
	}
	l(mm,bn) = aa.l(mm,bn);
	for (unsigned jj = bn+1;  jj<mm;  jj++) {
	  /* l(mm,jj) = aa.l(mm,jj) - dot(mm,jj,jj); */
	  subtract_dot_product(mm,jj,jj,aa.l(mm,jj));
	}
	{ /* jj == mm */
	  /* d(mm,mm) = aa.d(mm,mm) - dot(mm,mm,mm); then test */
	  if (subtract_dot_product(mm,mm,mm,aa.d(mm,mm)) == 0.) {
	    error(bWARNING, "open circuit: internal node %u\n", mm);
	    d(mm,mm) = _min_pivot;
	  }else{
	  }
	}
      }else{    /* bn == mm */
	d(mm,mm) = aa.d(mm,mm);
	if (d(mm,mm)==0.) {untested();
	  d(mm,mm) = _min_pivot;
	}else{
	}
      }
    }else{
    }
  }
}
/*--------------------------------------------------------------------------*/
template <class T>
void BSMATRIX<T>::lu_decomp()
{
  assert(_lownode);
  for (unsigned mm = 1;   mm <= size();   ++mm) {
    unsigned bn = _lownode[mm];
    if (bn < mm) {
      u(bn,mm) /= d(bn,bn);
      for (unsigned ii =bn+1;  ii<mm;  ii++) {
	/* (m(ii,mm) -= dot(ii,mm,ii)) /= d(ii,ii); */
	subtract_dot_product(ii,mm,ii) /= d(ii,ii);
      }
      for (unsigned jj = bn+1;  jj<mm;  jj++) {
	/* m(mm,jj) -= dot(mm,jj,jj); */
	subtract_dot_product(mm,jj,jj);
      }
      { /* jj == mm */
	/* m(mm,mm) -= dot(mm,mm,mm); then test */
	if (subtract_dot_product(mm,mm,mm) == 0.) {
	  error(bWARNING, "open circuit: internal node %u\n", mm);
	  d(mm,mm) = _min_pivot;
	}else{
	}
      }
    }else{    /* bn == mm */
      if (d(mm,mm)==0.) {
	d(mm,mm) = _min_pivot;
      }else{
      }
    }
  }
}
/*--------------------------------------------------------------------------*/
/* fbsub: forward and back sub, shared storage
 * v = right side vector, changed in place to solution vector
 */
template <class T>
void BSMATRIX<T>::fbsub(T* v) const
{
  assert(_lownode);
  assert(v);

  for (unsigned ii = 1; ii <= size(); ++ii) {	/* forward substitution */
    for (unsigned jj = _lownode[ii]; jj < ii; ++jj) {
      v[ii] -= l(ii,jj) * v[jj];
    }
    v[ii] /= d(ii,ii);
  }

  for (unsigned jj = size(); jj > 1; --jj) {		/* back substitution    */
    for (unsigned ii = _lownode[jj]; ii < jj; ++ii) {
      v[ii] -= u(ii,jj) * v[jj];
    }
  }
}
/*--------------------------------------------------------------------------*/
/* fbsub: forward and back sub, separate storage
 * x = solution vector
 * b = right side vector
 * c = intermediate vector after fwd sub
 */
template <class T>
void BSMATRIX<T>::fbsub(T* x, const T* b, T* c) const
{
  assert(_lownode);
  assert(x);
  assert(b);
  assert(c);

  {
    unsigned ii = 1;
    for (   ; ii <= size(); ++ii) {
      if (b[ii] != 0.) {
	break;
      }else{
      }
      c[ii] = 0.;
    }

    unsigned first_nz = ii;
    for (   ; ii <= size(); ++ii) {		/* forward substitution */
      unsigned low_node = std::max(_lownode[ii], first_nz);
      c[ii] = b[ii];
      for (unsigned jj = low_node; jj < ii; ++jj) {
	c[ii] -= l(ii,jj) * c[jj];
      }
      c[ii] /= d(ii,ii);
    }
  }

  notstd::copy_n(c, size()+1, x);

  for (unsigned jj = size(); jj > 1; --jj) {		/* back substitution    */
    for (unsigned ii = _lownode[jj]; ii < jj; ++ii) {
      x[ii] -= u(ii,jj) * x[jj];
    }
  }
  x[0] = 0.;
  // index starts at 1, but node 0 is ground
  // x[0]==0 eliminates a lot of "if" statements
}
/*--------------------------------------------------------------------------*/
template<class T>
BSMATRIX<T>::BSMATRIX(const BSMATRIX<T>& m) :
	_changed(new bool[m.size()+1]),
	_lownode(new unsigned[m.size()+1]),
	_space(new T[m._nzcount]),
	_rowptr(new T*[m._size+1]),
	_colptr(new T*[m._size+1]),
	_diaptr(new T*[m._size+1]),
	_nzcount(m._nzcount),
	_size(m._size),
	_zero(0),
	_trash(m._trash),
	_min_pivot(m._min_pivot)
{
	memcpy(_changed,m._changed, (_size+1) * sizeof(bool) );
	memcpy(_lownode,m._lownode, (_size+1) * sizeof(unsigned) );
	memcpy(_space,m._space, (_nzcount)  * sizeof(T));

	for(unsigned i=0; i<=_size; i++){
		_rowptr[i]= _space - intptr_t(m._space - m._rowptr[i]);
		_colptr[i]= _space - intptr_t(m._space - m._colptr[i]);
		_diaptr[i]= _space - intptr_t(m._space - m._diaptr[i]);
	}

}

/*--------------------------------------------------------------------------*/
// template <class X, class Y>
//Y& operator<< ( Y &o, const BSMATRIX<X>& m);

//template <class X>
//OMSTREAM& operator<< ( OMSTREAM &o, const BSMATRIX<X>& m);
template<>
void BSMATRIX<double>::sink_forward(unsigned* nm);
template<>
void BSMATRIX<double>::sink_reverse(unsigned* nm);
/*--------------------------------------------------------------------------*/
/* fbsubt: forward and back substitution with implicitly transposed matrix Ut Lt x = v
 * v = right side vector, changed in place to solution vector
 * GS:
 * this method s used to solve system A_t X = B (_t - transposed)
 * which corresponds to adjoint system
 * (LU)_t then transforms to U_t L_t
 *  added: Gennadiy Serdyuk <gserdyuk@gserdyuk.com>
 */
template <class T>
void BSMATRIX<T>::fbsubt(T* v) const
{
  assert(_lownode);
  assert(v);

  for (unsigned ii = 1; ii <= size(); ++ii) {	// forward substitution
    for (unsigned jj = _lownode[ii]; jj < ii; ++jj) {
      v[ii] -= u(jj,ii) * v [jj];
    }
  }

  for (unsigned jj = size(); jj > 1; --jj) {		// back substitution
	v[jj] /= d(jj,jj);
    for (unsigned ii = _lownode[jj]; ii < jj; ++ii) {
      v[ii] -= l(jj,ii) * v[jj];			
    }
  }
  v[1]/=d(1,1);
}
/*--------------------------------------------------------------------------*/
// apply real matrix to complex vector.
template<class T>
void BSMATRIX<T>::fbsub(std::complex<T>* v)const
{
  T part[size()+1];
  for(unsigned i=0; i<=size(); ++i){
    part[i] = v[i].real();
  }
  fbsub(part);
  for(unsigned i=0; i<=size(); ++i){
    v[i].real(part[i]);
    part[i] = v[i].imag();
  }
  fbsub(part);
  for(unsigned i=0; i<=size(); ++i){
    v[i].imag(part[i]);
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
#endif

// vim:ts=8:sw=2:noet:
