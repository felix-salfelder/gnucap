/*                               -*- C++ -*-
 * Copyright (C) 2014 Felix Salfelder
 * Author: Felix Salfelder
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
 * gsl supplementary header
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#if !defined(GNUCAP_GSL) && defined(HAVE_GSL_FIT_H)
#define GNUCAP_GSL

#ifdef NDEBUG
#define GSL_RANGE_CHECK_OFF
#endif

#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>

#include "m_matrix.h"

typedef int64_t index_t; // fixme: to_grid etc. to ssp header

inline double norm(const gsl_complex& x)
{
	return sqrt(GSL_REAL(x)*GSL_REAL(x)+GSL_IMAG(x)+GSL_IMAG(x));
}
/*--------------------------------------------------------------------------*/
inline double norm(const gsl_vector_complex& v)
{
	size_t n = v.size;
	double ret=0;
	for(unsigned i=0; i<n; ++i){
		gsl_complex x = gsl_vector_complex_get(&v,i);
		ret += GSL_REAL(x) * GSL_REAL(x);
		ret += GSL_IMAG(x) * GSL_IMAG(x);
	}
	return sqrt(ret);
}
/*--------------------------------------------------------------------------*/
inline gsl_vector& to_norm(gsl_vector& ret, const gsl_vector_complex& v)
{
	size_t n = v.size;
	assert(v.size==ret.size);
	for(unsigned i=0; i<n; ++i){
		const gsl_complex x = *gsl_vector_complex_const_ptr(&v,i);
		gsl_vector_set(&ret,i,norm(x));
	}
	return ret;
}
/*--------------------------------------------------------------------------*/
inline gsl_vector_complex& complex_array_to_gsl_vector_complex(gsl_vector_complex& v,COMPLEX* p, size_t n)
{
	for(unsigned i=0; i<n; ++i){
		gsl_complex x;
		GSL_REAL(x)=real(p[i]);
		GSL_IMAG(x)=imag(p[i]);
		gsl_vector_complex_set(&v,i,x);
	}
	return v;
}
/*--------------------------------------------------------------------------*/
//gsl_vector_complex& operator=(gsl_vector_complex& v, const gsl_vector& in)
inline gsl_vector_complex& set_vector (gsl_vector_complex& v, const gsl_vector& in)
{
	size_t n = in.size;
	assert(v.size==in.size);
	for(unsigned i=0; i<n; ++i){
		gsl_complex x;
		GSL_REAL(x) = gsl_vector_get(&in, i);
		GSL_IMAG(x) = 0;
		gsl_vector_complex_set(&v,i,x);
	}
	return v;
}
/*--------------------------------------------------------------------------*/
inline gsl_vector_complex& real_array_to_gsl_vector(gsl_vector_complex& v,double* p, size_t n)
{
	for(unsigned i=0; i<n; ++i){
		gsl_complex x;
		GSL_REAL(x) = p[i];
		GSL_IMAG(x) = 0;
		gsl_vector_complex_set(&v,i,x);
	}
	return v;
}
/*--------------------------------------------------------------------------*/
inline gsl_vector& array_to_gsl_vector(gsl_vector& v,double* p, size_t n)
{
	for(unsigned i=0; i<n; ++i){
		gsl_vector_set(&v,i,p[i]);
	}
	return v;
}
/*--------------------------------------------------------------------------*/
inline int sort_vector_index (gsl_permutation * , const gsl_vector * , size_t& )
{ incomplete();
	return 1;
}
/*--------------------------------------------------------------------------*/
inline int sort_vector_and_index (gsl_permutation * p, gsl_vector * v, size_t& n)
{
	int err;
	err = gsl_sort_vector_index (p, v);
	err+= gsl_permute_vector (p, v);

	for(unsigned i=0; i<n; ++i){
		if (!is_number(gsl_vector_get(v,i))){
			n = i;
			break;
		}
	}

	return err;
}
/*--------------------------------------------------------------------------*/
inline int mul(gsl_vector_complex& x, const gsl_vector& M)
{
	int err;
	gsl_vector_view xr = gsl_vector_complex_real (&x);
	gsl_vector_view xi = gsl_vector_complex_real (&x);
	err = gsl_vector_mul(&xr.vector,&M);
	err+= gsl_vector_mul(&xi.vector,&M);
	return err;
}
/*--------------------------------------------------------------------------*/
inline int mul(CBLAS_TRANSPOSE_t t, gsl_matrix M, const gsl_vector_complex& x, gsl_vector_complex& Mx)
{
	int err;
	gsl_vector_const_view xr = gsl_vector_complex_const_real (&x);
	gsl_vector_const_view xi = gsl_vector_complex_const_real (&x);
	gsl_vector_view Mxr = gsl_vector_complex_real (&Mx);
	gsl_vector_view Mxi = gsl_vector_complex_real (&Mx);
	err = gsl_blas_dgemv( t, 1., &M, &xr.vector, 0., &Mxr.vector );
	err+= gsl_blas_dgemv( t, 1., &M, &xi.vector, 0., &Mxi.vector );
	return err;
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
inline gsl_matrix& operator>>( const BSMATRIX<double>& in, gsl_matrix& out)
{
	size_t n = in.size();
	assert(out.size1 == n);
	assert(out.size2 == n);

	for(unsigned i=0; i<n; ++i){
		for(unsigned j=0; j<n; ++j){
			gsl_matrix_set(&out, i, j, in.s(i+1,j+1));
		}
	}
	return out;
}
/*--------------------------------------------------------------------------*/
inline double to_grid(double j, index_t& index, double mod)
{
	index = int64_t(floor((j+mod/2.)/mod));
	return double(index)*mod;
}
/*--------------------------------------------------------------------------*/
inline int64_t to_grid_index(double j, double mod)
{
	int index = int(floor((j+mod/2.)/mod));
	return index;
}
/*--------------------------------------------------------------------------*/
inline double to_grid(double j, double mod)
{
	int64_t index = to_grid_index(j,mod);
	return double(index)*mod;
}
/*--------------------------------------------------------------------------*/
inline gsl_vector_complex& to_grid( gsl_vector_complex& out, const gsl_vector& d, index_t* index)
{
	size_t n = d.size;
	assert(out.size==n);

	gsl_vector_view r = gsl_vector_complex_real(&out);
	gsl_vector_view i = gsl_vector_complex_imag(&out);

	for(unsigned ii=0; ii<n; ++ii) {
		double norm = gsl_vector_get(&d,ii);
		if(! std::isfinite(norm)){ itested();
			gsl_vector_set(&r.vector,ii,0.);
			gsl_vector_set(&i.vector,ii,0.);
		}else if(is_number(norm)){
			double r_ = *gsl_vector_const_ptr(&r.vector,ii);
			gsl_vector_set(&r.vector,ii,to_grid(r_,index[2*ii],norm));
			double i_ = *gsl_vector_const_ptr(&i.vector,ii);
			gsl_vector_set(&i.vector,ii,to_grid(i_,index[2*ii],norm));
		}else{ unreachable();
			//gsl_vector_set(&r.vector,ii,0.);
			//gsl_vector_set(&i.vector,ii,0.);
		}
	}
	return out;
}
/*--------------------------------------------------------------------------*/
inline gsl_vector_complex& to_grid( gsl_vector_complex_view out, const gsl_vector& d, index_t* index)
{
	return to_grid (out.vector,d, index);
}
/*--------------------------------------------------------------------------*/
inline gsl_vector_complex& to_grid( gsl_vector_complex& out, const gsl_vector& d, unsigned num=INT_MAX)
{
	size_t n = d.size;
	assert(out.size==n);

	gsl_vector_view r = gsl_vector_complex_real(&out);
	gsl_vector_view i = gsl_vector_complex_imag(&out);

	for(unsigned ii=0; ii<n && ii<num; ++ii) {
		double norm = gsl_vector_get(&d,ii);
		if(! std::isfinite(norm)){ itested();
			gsl_vector_set(&r.vector,ii,0.);
			gsl_vector_set(&i.vector,ii,0.);
		}else if(is_number(norm)){
			double r_ = *gsl_vector_const_ptr(&r.vector,ii);
			gsl_vector_set(&r.vector,ii,to_grid(r_,norm));
			double i_ = *gsl_vector_const_ptr(&i.vector,ii);
			gsl_vector_set(&i.vector,ii,to_grid(i_,norm));
		}else{ unreachable();
			//gsl_vector_set(&r.vector,ii,0.);
			//gsl_vector_set(&i.vector,ii,0.);
		}
	}
	return out;
}
/*--------------------------------------------------------------------------*/
inline gsl_vector_complex& to_grid( gsl_vector_complex& out, const gsl_vector_complex& d)
{ incomplete();
	size_t n = d.size;
	assert(out.size==n);

	gsl_vector_view r = gsl_vector_complex_real(&out);
	gsl_vector_view i = gsl_vector_complex_imag(&out);

	for(unsigned ii=0; ii<n; ++ii) {
		gsl_complex lambda = *gsl_vector_complex_const_ptr(&d,ii);
		double norm = sqrt(pow(GSL_IMAG(lambda),2.) + pow(GSL_REAL(lambda), 2.));

		if(! std::isfinite(norm)){ untested();
			gsl_vector_set(&r.vector,ii,0.);
			gsl_vector_set(&i.vector,ii,0.);
		} else if(is_number(norm)){
			double r_ = *gsl_vector_const_ptr(&r.vector,ii);
			gsl_vector_set(&r.vector,ii,to_grid(r_,norm));
			double i_ = *gsl_vector_const_ptr(&i.vector,ii);
			gsl_vector_set(&i.vector,ii,to_grid(i_,norm));
		}else{ unreachable();
			//gsl_vector_set(&r.vector,ii,0.);
			//gsl_vector_set(&i.vector,ii,0.);
		}
	}
	return out;
}
/*--------------------------------------------------------------------------*/
// hmmm
inline gsl_vector_complex& to_grid( gsl_vector_complex& out, const gsl_complex& d)
{ untested();
	size_t n = out.size;

	gsl_vector_view r = gsl_vector_complex_real(&out);
	gsl_vector_view i = gsl_vector_complex_imag(&out);

	for(unsigned ii=0; ii<n; ++ii) {
		if(! std::isfinite(GSL_REAL(d)) || !std::isfinite(GSL_IMAG(d))){ untested();
			gsl_vector_set(&r.vector,ii,0.);
			gsl_vector_set(&i.vector,ii,0.);
		} else if(is_number(GSL_REAL(d)) && is_number(GSL_IMAG(d)) ){
			double r_ = *gsl_vector_const_ptr(&r.vector,ii);
			gsl_vector_set(&r.vector,ii,to_grid(r_,GSL_REAL(d)));
			double i_ = *gsl_vector_const_ptr(&i.vector,ii);
			gsl_vector_set(&i.vector,ii,to_grid(i_,GSL_IMAG(d)));
		}else{ unreachable();
			//gsl_vector_set(&r.vector,ii,0.);
			//gsl_vector_set(&i.vector,ii,0.);
		}
	}
	return out;
}
/*--------------------------------------------------------------------------*/
inline gsl_vector& to_grid( gsl_vector& out, const double& d)
{
	size_t n = out.size;
	for(unsigned ii=0; ii<n; ++ii) {
		if(! std::isfinite(d) ){ untested();
			gsl_vector_set(&out,ii,0.);
		} else if(is_number(d) ) {
			double r_ = gsl_vector_get(&out,ii);
			gsl_vector_set(&out,ii,to_grid(r_,d));
		}else{ unreachable();
			//gsl_vector_set(&r.vector,ii,0.);
			//gsl_vector_set(&i.vector,ii,0.);
		}
	}
	return out;

}
/*--------------------------------------------------------------------------*/
inline gsl_vector_complex& to_grid( gsl_vector_complex& out, const double& d)
{
	return to_grid(out,gsl_complex_rect(d,d));
}
/*--------------------------------------------------------------------------*/
inline gsl_vector& operator/=( gsl_vector& out, const gsl_vector& d)
{
	size_t n = d.size;
	assert(out.size==n);

	for(unsigned ii=0; ii<n; ++ii) {
		if(0. != *gsl_vector_const_ptr(&d,ii)){
			gsl_vector_set(&out,ii,*gsl_vector_const_ptr(&out,ii)/ *gsl_vector_const_ptr(&d,ii));
		}else{
			gsl_vector_set(&out,ii,1./0.);
		}
	}
	return out;
}
/*--------------------------------------------------------------------------*/
inline gsl_vector_complex& operator/=( gsl_vector_complex& out, const gsl_vector& d)
{
	size_t n = d.size;
	assert(out.size==n);

	gsl_vector_view r = gsl_vector_complex_real(&out);
	gsl_vector_view i = gsl_vector_complex_imag(&out);

	for(unsigned ii=0; ii<n; ++ii) {
		if(0. != *gsl_vector_const_ptr(&d,ii)){
			gsl_vector_set(&r.vector,ii,*gsl_vector_const_ptr(&r.vector,ii)/ *gsl_vector_const_ptr(&d,ii));
			gsl_vector_set(&i.vector,ii,*gsl_vector_const_ptr(&i.vector,ii)/ *gsl_vector_const_ptr(&d,ii));
		}else{
			gsl_vector_set(&r.vector,ii,1./0.);
			gsl_vector_set(&i.vector,ii,0.);
		}
	}
	return out;
}
/*--------------------------------------------------------------------------*/
inline gsl_vector_complex& operator%=( gsl_vector_complex& out, const gsl_vector& d)
{
	size_t n = d.size;
	assert(out.size==n);

	gsl_vector_view r = gsl_vector_complex_real(&out);
	gsl_vector_view i = gsl_vector_complex_imag(&out);

	for(unsigned ii=0; ii<n; ++ii) {
		if(1 || (0. != *gsl_vector_const_ptr(&d,ii))){
			gsl_vector_set(&r.vector,ii,fmod(*gsl_vector_const_ptr(&r.vector,ii), *gsl_vector_const_ptr(&d,ii)));
			gsl_vector_set(&i.vector,ii,fmod(*gsl_vector_const_ptr(&i.vector,ii), *gsl_vector_const_ptr(&d,ii)));
		}else{
			gsl_vector_set(&r.vector,ii,1./0.);
			gsl_vector_set(&i.vector,ii,0.);
		}
	}
	return out;
}
/*--------------------------------------------------------------------------*/
inline int discretize(gsl_vector_complex &V, int64_t *index, double mod=1e-4)
{
	size_t n = V.size;
	for(unsigned j=0; j<n; ++j){
		gsl_complex* p = gsl_vector_complex_ptr (&V, j);

		GSL_SET_COMPLEX(
				p,
				to_grid( GSL_REAL(*p), index[2*j], mod),
				to_grid( GSL_IMAG(*p), index[2*j+1], mod)
				);
	}

	return 0;
}
/*--------------------------------------------------------------------------*/
inline int matrix_div_cols(gsl_matrix_complex &m, const gsl_vector& v)
{
	size_t n = v.size;
	int err = 0;
	gsl_vector_complex_view col;
	for(unsigned i=0; i<n; ++i){
		col = gsl_matrix_complex_row(&m, i);
		gsl_vector_view r_ = gsl_vector_complex_real(&col.vector);
		gsl_vector_view i_ = gsl_vector_complex_imag(&col.vector);

		err+= gsl_vector_div(&r_.vector,&v);
		err+= gsl_vector_div(&i_.vector,&v);
	}
	return err;
}
/*--------------------------------------------------------------------------*/
inline int matrix_scale_cols(gsl_matrix_complex &m, const gsl_vector& v)
{
	size_t n = v.size;
	int err = 0;
	gsl_vector_complex_view col;
	for(unsigned i=0; i<n; ++i){
		col = gsl_matrix_complex_row(&m, i);
		gsl_vector_view r_ = gsl_vector_complex_real(&col.vector);
		gsl_vector_view i_ = gsl_vector_complex_imag(&col.vector);

		err+= gsl_vector_mul(&r_.vector,&v);
		err+= gsl_vector_mul(&i_.vector,&v);
	}
	return err;
}
/*--------------------------------------------------------------------------*/
inline int discretize(gsl_matrix_complex &M, double mod=1e-4)
{
	size_t n = M.size1;
	for(unsigned i=0; i<n; ++i){
		for(unsigned j=0; j<n; ++j){
			gsl_complex* p = gsl_matrix_complex_ptr (&M, i, j);

			GSL_REAL(*p) = to_grid( GSL_REAL(*p), mod);
			GSL_IMAG(*p) = to_grid( GSL_IMAG(*p), mod);
		}
	}
	return 0;
}
/*--------------------------------------------------------------------------*/
inline gsl_vector_const_view rest(const gsl_vector& x)
{
	// x must have at least two cols.
	return gsl_vector_const_subvector(&x, 1, x.size-1);
}
/*--------------------------------------------------------------------------*/
inline gsl_vector_view rest(gsl_vector& x)
{
	// x must have at least two cols.
	return gsl_vector_subvector(&x, 1, x.size-1);
}
/*--------------------------------------------------------------------------*/
inline gsl_matrix_complex_view rest(gsl_matrix_complex& x)
{
	// x must have at least two cols.
	return gsl_matrix_complex_submatrix(&x, 0, 1, x.size1, x.size2-1);
}
/*--------------------------------------------------------------------------*/
inline gsl_vector_complex_view firstcol(gsl_matrix_complex& x)
{
	return gsl_matrix_complex_column(&x, 0);
}
/*--------------------------------------------------------------------------*/
inline gsl_vector_complex_view firstcol(gsl_matrix_complex_view x)
{
	return gsl_matrix_complex_column(&x.matrix, 0);
}
/*--------------------------------------------------------------------------*/
inline ostream& operator<<( ostream& s, const gsl_permutation& p)
{
	size_t n = p.size;

	for(unsigned i=0; i<n; ++i){
		if (i) s<< ", ";
		s << gsl_permutation_get(&p,i);
	}
	return s;
}
/*--------------------------------------------------------------------------*/
inline ostream& operator<<( ostream& s, const gsl_complex& z)
{
	s << "(";
	s << GSL_REAL(z);
	if(GSL_IMAG(z)!=0.){
		s << "+i*";
		s << GSL_IMAG(z);
	}
	s << ")";
	return s;
}
/*--------------------------------------------------------------------------*/
inline ostream& operator<<( ostream& s, const gsl_vector& v)
{
	size_t n = v.size;

	s << "(";
	for(unsigned i=0; i<n; ++i){
		if (i) s<< ", ";
		s << gsl_vector_get(&v,i);
	}
	s << ")";
	return s;
}
/*--------------------------------------------------------------------------*/
inline ostream& operator<<( ostream& s, gsl_vector_const_view v)
{
	return s << v.vector;
}
/*--------------------------------------------------------------------------*/
inline gsl_complex* gsl_vector_complex_ptr(gsl_vector_complex_view z0, size_t i)
{
	return gsl_vector_complex_ptr(&z0.vector,i);
}
/*--------------------------------------------------------------------------*/
inline gsl_complex gsl_vector_complex_get(gsl_vector_complex_const_view z0, size_t i)
{ itested();
	return gsl_vector_complex_get(&z0.vector,i);
}
/*--------------------------------------------------------------------------*/
inline gsl_complex gsl_vector_complex_get(gsl_vector_complex_view z0, size_t i)
{ itested();
	return gsl_vector_complex_get(&z0.vector,i);
}
/*--------------------------------------------------------------------------*/
inline double gsl_vector_get(gsl_vector_const_view z0, size_t i)
{ itested();
	return gsl_vector_get(&z0.vector,i);
}
/*--------------------------------------------------------------------------*/
inline double gsl_vector_get(gsl_vector_view z0, size_t i)
{ itested();
	return gsl_vector_get(&z0.vector,i);
}
/*--------------------------------------------------------------------------*/
inline int gsl_vector_complex_memcpy(gsl_vector_complex_view dest, gsl_vector_complex_const_view src)
{ itested();
	return gsl_vector_complex_memcpy(&dest.vector, &src.vector);
}
/*--------------------------------------------------------------------------*/
inline int gsl_vector_memcpy(gsl_vector_view dest, gsl_vector_const_view src)
{ untested();
	return gsl_vector_memcpy(&dest.vector, &src.vector);
}
/*--------------------------------------------------------------------------*/
inline int gsl_matrix_memcpy(gsl_matrix* dest, gsl_matrix_const_view src)
{ untested();
	return gsl_matrix_memcpy(dest, &src.matrix);
}
/*--------------------------------------------------------------------------*/
inline int gsl_matrix_memcpy(gsl_matrix_view dest, gsl_matrix_const_view src)
{ untested();
	return gsl_matrix_memcpy(&dest.matrix, &src.matrix);
}
/*--------------------------------------------------------------------------*/
inline int gsl_matrix_memcpy_real(gsl_matrix* out, const gsl_matrix_complex* in)
{
	assert(out->size1 == in->size1);
	assert(out->size2 == in->size2);
	int err=0;
	for(unsigned i=0; i<in->size1; ++i){ untested();
		gsl_vector_complex_const_view row = gsl_matrix_complex_const_row(in,i);
		gsl_vector_const_view r = gsl_vector_complex_const_real(&row.vector);
		gsl_vector_view o = gsl_matrix_row(out, i);
		err+= gsl_vector_memcpy(o,r);
	}
	return err;
}
/*--------------------------------------------------------------------------*/
inline int gsl_matrix_memcpy_real(gsl_matrix* dest, gsl_matrix_complex_const_view src)
{
	return gsl_matrix_memcpy_real(dest, &src.matrix);
}
/*--------------------------------------------------------------------------*/
inline ostream& operator<<( ostream& s, const gsl_vector_complex& v)
{
	size_t n = v.size;

	gsl_vector_const_view r = gsl_vector_complex_const_real(&v);
	gsl_vector_const_view i = gsl_vector_complex_const_imag(&v);

	s << "(";
	for(unsigned ii=0; ii<n; ++ii) {
		if (ii) s<< ", ";
		s << gsl_vector_get(&r.vector,ii);
		double tmp = gsl_vector_get(&i.vector,ii);
		if (0. != tmp) {
			s << "+i*" << tmp;
		}
	}
	s << ")";
	return s;
}
/*--------------------------------------------------------------------------*/
inline ostream& operator<<( ostream& s, const gsl_matrix& m)
{
	size_t I = m.size1;
	size_t J = m.size2;

	s << "(";
	for(unsigned i=0; i<I; ++i){
		if (i) s<< ",\n";
		s << "(";
		for(unsigned j=0; j<J; ++j){
			if (j) s<< ", ";
			s << gsl_matrix_get(&m,i,j);
		}
		s << ")";
	}
	s << ")";
	return s;
}
/*--------------------------------------------------------------------------*/
inline ostream& operator<<( ostream& s, const gsl_matrix_complex& m)
{
	size_t I = m.size1;
	size_t J = m.size2;

	s << "(";
	for(unsigned i=0; i<I; ++i){
		if (i) s<< ",\n";
		s << "(";
		for(unsigned j=0; j<J; ++j){
			if (j) s<< ", ";
			gsl_complex x = gsl_matrix_complex_get(&m,i,j);
			s << x.dat[0];
			if(0. != x.dat[1]) {
				s << "+i*" << x.dat[1];
			}
		}
		s << ")";
	}
	s << ")";
	return s;
}
/*--------------------------------------------------------------------------*/
inline ostream& operator<<( ostream& s, const gsl_vector_view& v)
{
	return s << v.vector;
}
/*--------------------------------------------------------------------------*/
inline ostream& operator<<( ostream& s, const gsl_vector_complex_view v){
	return s << v.vector;
}
/*--------------------------------------------------------------------------*/
inline ostream& operator<<( ostream& s, gsl_vector_complex_const_view v)
{
	return s << v.vector;
}
/*--------------------------------------------------------------------------*/
typedef int integer;
int gsl_blas_zgelss ( gsl_matrix_complex * A,
		gsl_matrix_complex * B,
		gsl_vector * S,
		double* RCOND, integer* RANK, integer* INFO );
/*--------------------------------------------------------------------------*/
inline int gsl_vector_sub(gsl_vector_complex_view a, const gsl_vector_complex* b)
{
	return gsl_vector_complex_sub(&a.vector, b);
}
/*--------------------------------------------------------------------------*/
inline int gsl_vector_sub(gsl_vector_complex_view a, gsl_vector_complex_view b)
{
	return gsl_vector_complex_sub(&a.vector, &b.vector);
}
/*--------------------------------------------------------------------------*/
inline int gsl_vector_sub(gsl_vector_complex_view a, gsl_vector_complex_const_view b)
{
	return gsl_vector_complex_sub(&a.vector, &b.vector);
}
/*--------------------------------------------------------------------------*/
inline int gsl_vector_mul(gsl_vector_complex_view a, gsl_vector_complex_view b)
{
	return gsl_vector_complex_mul(&a.vector, &b.vector);
}
/*--------------------------------------------------------------------------*/
inline int gsl_vector_scale(gsl_vector_view a, double s)
{
	return gsl_vector_scale(&a.vector, s);
}
/*--------------------------------------------------------------------------*/
inline int gsl_vector_mul(gsl_vector_complex_view a, gsl_vector_complex_const_view b)
{
	return gsl_vector_complex_mul(&a.vector, &b.vector);
}
/*--------------------------------------------------------------------------*/
inline int gsl_vector_sub(gsl_vector_view& a, gsl_vector_view& b)
{
	return gsl_vector_sub(&a.vector, &b.vector);
}
#endif
// vim works with any ts=sw
