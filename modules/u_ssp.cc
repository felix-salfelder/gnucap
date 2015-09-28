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
 * ssp implementation
 */
#include "u_ssp.h"
#include "u_opt.h"

/*--------------------------------------------------------------------------*/
SSP_TREE SSP_TREE::_root;
unsigned SSP_TREE::_num_charts;
unsigned SSP_TREE::_num_states;
unsigned SSP_TREE::_num_spls;
/*--------------------------------------------------------------------------*/
SSP_CHART& SSP_TREE::insert(gsl_matrix_complex& M, const gsl_vector& grid, unsigned n_inp)
{ itested();
	gsl_vector_complex_view f = firstcol(M);
	double g = gsl_vector_get(&grid,0);
	SSP_VECTOR v(f.vector, g);
	trace1("insert", f.vector);
	trace1("insert", v);
	SSP_TREE*& item = _list[v];

	if(M.size2 > 1){
		if(item == NULL){
			item = new SSP_TREE();
		}
		assert(!dynamic_cast<SSP_CHART*>(item));

		gsl_matrix_complex_view tail = rest(M);
		gsl_vector_const_view gridtail = rest(grid);
		return item->insert(tail.matrix, gridtail.vector, n_inp);
	}else{
		assert(M.size2==1);
		if(item == NULL){
			item = new SSP_CHART(n_inp);
		}else{
			assert(dynamic_cast<SSP_CHART*>(item));
			dynamic_cast<SSP_CHART*>(item)->_new = false;
		}
		assert(dynamic_cast<SSP_CHART*>(item)->n_inputs() == n_inp);
		return *(dynamic_cast<SSP_CHART*>(item));
	}
}
/*--------------------------------------------------------------------------*/
// FIXME: quantize here?
SSP_STATE& SSP_CHART::insert_state(const index_t* index, unsigned size)
{ itested();
	SSP_VECTOR v(index, size);
	SSP_TREE*& item = _list[v];
	if(item == NULL){ itested();
		trace0("new state (generic)");
		item = new SSP_STATE();
	}else{
		item->visit();
	}
	assert(dynamic_cast<SSP_STATE*>(item));
	return *(dynamic_cast<SSP_STATE*>(item));
}
/*--------------------------------------------------------------------------*/
// insert generic state.
SSP_STATE& SSP_CHART::insert()
{ itested();
	SSP_TREE*& item = _list[SSP_VECTOR()];
	if(item == NULL){ itested();
		trace0("new state (generic)");
		item = new SSP_STATE();
	}
	assert(dynamic_cast<SSP_STATE*>(item));
	return *(dynamic_cast<SSP_STATE*>(item));
}
/*--------------------------------------------------------------------------*/
SSP_STATE& SSP_CHART::insert(gsl_vector_complex& f, const double grid)
{ untested();
	SSP_VECTOR v(f, grid);
	SSP_TREE*& item = _list[v];
	if(item == NULL){ untested();
		trace1("new state at", v);
		item = new SSP_STATE();
	}
	assert(dynamic_cast<SSP_STATE*>(item));
	return *(dynamic_cast<SSP_STATE*>(item));
}
/*--------------------------------------------------------------------------*/
SSP_SPL& SSP_STATE::insert_spl(gsl_vector &op, const SSP_CHART* c, const double grid)
{ itested();
	SSP_VECTOR v(&op,grid);
	SSP_TREE*& item = _list[v];
	visit();
	if(item == NULL){ itested();
		trace1("new spl at", v);
		assert(c);
		item = new SSP_SPL(c);
	}else{ itested();
		SSP_SPL& s = *(prechecked_cast<SSP_SPL*>(item));
		s.visit();
	}
	assert(dynamic_cast<SSP_SPL*>(item));
	return *(dynamic_cast<SSP_SPL*>(item));
}
/*--------------------------------------------------------------------------*/
// hack.
SSP_TRANS SSP_SPL::insert_trans(const SSP_VECTOR& v, SSP_SPL* tgt)
{ untested();
	SSP_TREE*& item = _list[v];
	if(item == NULL){ itested();
		trace0("new trans");
		item = tgt;
	}
	iterator it = _list.find(v);
	return it;
}
/*--------------------------------------------------------------------------*/
SSP_TRANS SSP_SPL::insert_trans(const index_t* inp, SSP_SPL* tgt)
{ itested();
	SSP_VECTOR v(inp, _chart->n_inputs());
	SSP_TREE*& item = _list[v];
	if(item == NULL){ itested();
		trace0("new trans");
		item = tgt;
	}
	iterator it = _list.find(v);
	return it;
}
/*--------------------------------------------------------------------------*/
SSP_TRANS SSP_SPL::insert_trans(const gsl_vector* inp, SSP_SPL* tgt)
{ itested();
	SSP_VECTOR v(inp,1.);
	SSP_TREE*& item = _list[v];
	if(item == NULL){ itested();
		trace0("new trans");
		item = tgt;
	}
	iterator it = _list.find(v);
	return it;
}
/*--------------------------------------------------------------------------*/
SSP_SPL& SSP_TRANS::spl()const
{
	SSP_SPL* s = prechecked_cast<SSP_SPL*>(_t->second);
	assert(s);
	return *s;
}
/*--------------------------------------------------------------------------*/
//index_t* SSP_SPL::clone_input() const
//{
//	size_t size = _chart->n_inputs();
//	index_t* ret = new index_t[size];
//
//	memcpy(ret, , sizeof(index_t) * size);
//}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
SSP_VECTOR::SSP_VECTOR(const SSP_VECTOR& v) :
	_size(v._size),
	_scale(v._scale)
{
	_data = new int64_t[_size];
	memcpy(_data, v._data, sizeof(int64_t)*_size);
}
/*--------------------------------------------------------------------------*/
index_t* SSP_VECTOR::clone_data() const
{
	index_t* ret = new index_t[_size];
	memcpy(ret, _data, sizeof(index_t) * _size);
	return ret;
}
/*--------------------------------------------------------------------------*/
SSP_VECTOR::SSP_VECTOR(const index_t* data, const unsigned size, double scale) :
	_size(size),
	_scale(scale)
{
	_data = new int64_t[_size];
	for(unsigned i=0; i<_size; ++i){
		_data[i] = data[i];
	}
}
/*--------------------------------------------------------------------------*/
SSP_VECTOR::SSP_VECTOR(const double* data, const unsigned size) :
	_size(size),
	_scale(sqrt(OPT::abstol))
{
	_data = new int64_t[_size];
	for(unsigned i=0; i<_size; ++i){
		_data[i] = int64_t(floor(data[i]/sqrt(OPT::abstol) + sqrt(OPT::abstol)*.5));
	}
}
/*--------------------------------------------------------------------------*/
SSP_VECTOR::SSP_VECTOR(const gsl_vector* f, double grid) : // hmm, complex grid?
	_size((unsigned)f->size),
	_scale(grid)
{
	_data = new int64_t[_size];
	for(unsigned i=0; i<_size; ++i){
		double data = gsl_vector_get(f,i);
		_data[i] = int64_t(floor(data/sqrt(OPT::abstol) + sqrt(OPT::abstol)*.5));
	}
}
/*--------------------------------------------------------------------------*/
SSP_VECTOR::SSP_VECTOR(gsl_vector_complex& f, const double& grid) : // hmm, complex grid?
	_size((unsigned)f.size*2),
	_scale(grid)
{
	_data = new int64_t[_size];
	discretize(f, _data, grid);
}
/*--------------------------------------------------------------------------*/
bool SSP_VECTOR::operator==(const SSP_VECTOR&) const
{
	incomplete();
	return false;
}
/*--------------------------------------------------------------------------*/
ostream& SSP_TREE::printlist(ostream& s) const
{
	for(const_iterator i=_list.begin(); _list.end()!=i; ++i){ itested();
		s << (*i).first << "->" << hp(i->second) << " ";
	}
	return s;
}
/*--------------------------------------------------------------------------*/
bool SSP_VECTOR::operator<(const SSP_VECTOR&r) const
{
	if(_size != r._size){ itested();
		return _size < r._size;
	}

	for(unsigned i=0; i<_size; ++i){ itested();
		if(_data[i]<r._data[i]){ itested();
			return true;
		}else if(_data[i]>r._data[i]){
			return false;
		}
	}
	return _scale<r._scale;
}
/*--------------------------------------------------------------------------*/
// index_t& SSP_TRANS::index(const unsigned i)
// {
// 	return _t->first.index(i);
// }
/*--------------------------------------------------------------------------*/
// id of target sample.
unsigned SSP_TRANS::id() const
{
	const SSP_SPL* s = prechecked_cast<const SSP_SPL*>(_t->second);
	assert(s);
	return s->id();
}
/*--------------------------------------------------------------------------*/
double SSP_TRANS::operator[](const unsigned i)const
{
	double r = _t->first[i];
	return r;
}
/*--------------------------------------------------------------------------*/
