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
 * ssp
 */
#ifndef U_SSP__
#define U_SSP__

#include "u_gsl.h"
#include <map>

/*--------------------------------------------------------------------------*/
class SSP_VECTOR {
	public:
		explicit SSP_VECTOR() : _size(0), _data(NULL){itested();}
		SSP_VECTOR(const SSP_VECTOR&);
		SSP_VECTOR(gsl_vector_complex&, const double& grid);
		SSP_VECTOR(const gsl_vector*, double grid=0);
//		SSP_VECTOR(const int* index, const unsigned size, double scale=0.);
		SSP_VECTOR(const index_t* index, const unsigned size, double scale=1.);
		SSP_VECTOR(const double* data, const unsigned size);
		~SSP_VECTOR(){itested();delete[] _data;}
		bool operator<(const SSP_VECTOR&) const;
		bool operator==(const SSP_VECTOR&) const;

		ostream& print( ostream& s) const;

		double operator[](const size_t i)const{
			assert(i<_size); return double(_data[i])*_scale;
		}
		index_t& index(const size_t i){
			assert(i<_size); return _data[i];
		}
		index_t* clone_data() const;

	private:
		unsigned _size; // n, ht, #inputs
		index_t* _data; // int32_t is possibly sufficient.
		double _scale;
};
/*--------------------------------------------------------------------------*/
class SSP_CHART; // chart is a linear operation region
                 // contains coeffs == states
class SSP_STATE; // charge configuration
                 // holds quantized state
class SSP_SPL;   // sampling point.
                 // stores ddc values and probes.
class SSP_TRANS; // (tgt_input, dt) -> spl
/*--------------------------------------------------------------------------*/
class SSP_TREE {
	public:
		explicit SSP_TREE() : _new(true) {}
		virtual ~SSP_TREE() {}; // later?
		SSP_TREE(const SSP_TREE& t):_list(t._list){untested();}
		size_t size() const{return _list.size();}

	public:
		SSP_CHART& insert(gsl_matrix_complex& M, const gsl_vector& grid, unsigned n_inp); // complex grid?!

		// static size_t _size;
		static SSP_TREE _root;

		ostream& printlist(ostream&) const;
		bool is_new()const {return _new;}
		void visit() {_new=false;}
		static unsigned num_charts(){return _num_charts;}
		static unsigned num_states(){return _num_states;}
		static unsigned num_spls(){return _num_spls;}

	protected:
		typedef std::map<SSP_VECTOR, SSP_TREE*> list_t;
		typedef list_t::const_iterator const_iterator;
		typedef list_t::iterator iterator;
		list_t _list;

	protected: // status
		static unsigned _num_charts;
		static unsigned _num_states;
		static unsigned _num_spls;
	public:
	private:
		bool _new;

};
// typedef std::map<SSP_VECTOR, SSP_SPL>::const_reference SSP_TRANS;
/*--------------------------------------------------------------------------*/
class SSP_TRANS {
	public:
		typedef std::map<SSP_VECTOR, SSP_TREE*>::iterator trans_t;
		SSP_TRANS(const trans_t& t) : _t(t) {untested();}
		SSP_TRANS(trans_t& t) : _t(t) {untested();}

		double operator[](const unsigned) const;
		double dt() const;
		
		int index(const unsigned) const;
//		index_t& index(const unsigned);
		typedef std::pair<SSP_VECTOR, SSP_TREE*> p_;
		p_ p() const {return *_t;}
		operator p_() const {return *_t;}

		operator SSP_VECTOR(){return _t->first;}

		SSP_VECTOR input()const {return _t->first;}

		unsigned id() const;

		SSP_SPL& spl()const;
		operator SSP_SPL&(){return spl();}

	private:
		trans_t _t;
};
/*--------------------------------------------------------------------------*/
class SSP_CHART : public SSP_TREE {
		SSP_CHART() : _inputs(0) {unreachable();}
		SSP_CHART(const SSP_CHART&) : SSP_TREE(), _inputs(0) {unreachable();}
	public:
		SSP_CHART(unsigned inp) : _inputs(inp), _id(_num_charts) {
			_inputgrid = new double[inp];
			++_num_charts;
			trace0("new chart");
			_new=true;
			assert(_id < (1<<20));
			assert(_num_charts > 0);
		};
		~SSP_CHART() {assert(num_charts); --_num_charts;}
		SSP_STATE& insert_state();
		SSP_STATE& insert_state(gsl_vector_complex&, const double grid=1.);

		SSP_STATE& insert_state(const index_t* index, unsigned size);

		SSP_STATE& insert();
		SSP_STATE& insert(gsl_vector_complex&, const double grid=1.);

		// 		SSP_GRID* insert_grid();
		// 		SSP_GRID* insert_grid(gsl_vector& z);

		double grid(unsigned x)const{
			assert(x<_inputs);
			return _inputgrid[x];
		}
		unsigned n_inputs()const{return _inputs;}
		unsigned id()const {return _id;}

		double* _inputgrid; // HACK. private!
		bool _new; //hack. private
	private: // ad-hoc hack
		const unsigned _inputs;
	private:
		unsigned _id;
};
/*--------------------------------------------------------------------------*/
class SSP_STATE : public SSP_TREE {
	public:
		explicit SSP_STATE() {++_num_states;}
		~SSP_STATE() {--_num_states;}

	public:
		SSP_SPL& insert_spl(gsl_vector &op, const SSP_CHART*, const double grid=1.);
//		SSP_SPL& insert_spl(gsl_vector &op, const double grid=1.);
//		SSP_SPL& insert_spl(gsl_vector_view op, const double grid=1.){
//			return insert_spl(op.vector, grid);
//		}
		SSP_SPL& insert_spl(gsl_vector_view op, const SSP_CHART* c, const double grid=1.){
			return insert_spl(op.vector, c, grid);
		}
};
/*--------------------------------------------------------------------------*/
// class SSP_SPL : public SSP_BASE  ??
class SSP_SPL : public SSP_TREE {
	public:
		explicit SSP_SPL(const SSP_CHART* c) : SSP_TREE(),
		  	_id(_num_spls), _chart(c) {assert(c); ++_num_spls;}
		~SSP_SPL() {--_num_spls;}
		// SSP_SPL& insert_trans(gsl_vector *inputs, double time, SSP_SPL*);
		// std::map<SSP_VECTOR, SSP_SPL>::const_reference
		SSP_TRANS insert_trans(const index_t *inputs, SSP_SPL*);
		SSP_TRANS insert_trans(const gsl_vector *inputs, SSP_SPL*);
		SSP_TRANS insert_trans(const SSP_VECTOR& inputs, SSP_SPL*);
		unsigned id()const{return _id;}
		void set_old(){visit();}

		const SSP_CHART& chart()const{assert(_chart); return *_chart;}
		unsigned n_inputs()const{return _chart->n_inputs();}

//		index_t* clone_input() const;
	private:
		unsigned _id;
		const SSP_CHART* _chart;
};
/*--------------------------------------------------------------------------*/
inline ostream& operator<<( ostream& s, const SSP_TREE& v)
{
	return v.printlist(s);
}
/*--------------------------------------------------------------------------*/
inline ostream& operator<<( ostream& s, const SSP_VECTOR& v){
	return v.print(s);
}
/*--------------------------------------------------------------------------*/
inline ostream& SSP_VECTOR::print( ostream& s) const
{
	if (_scale!=1) {
		s << _scale << "*";
	}
		s << "(";
	for(unsigned i=0; i<_size; ++i){
		if (i) s<< ", ";
		s << _data[i];
		//     if(_data[2*i+1]){
		//       s << "+i*" << _data[2*i+1];
		//     }
	}
	return s << ")";
}
#endif
