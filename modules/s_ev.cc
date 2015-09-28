/* Copyright (C) 2014 Felix Salfelder
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
 * eigenvalue stuff
 */

#define ADD_VERSION

#include "config.h"

#ifdef HAVE_CLAPACK_H
// #include "u_atlas.h"
#endif


#ifndef HAVE_CLAPACK_H
#warning "untested"
#endif

#include "u_status.h"
#include "u_ssp.h"
#include "s_ev.h"
#include "u_gsl.h"
#include "u_prblst.h"
#include "u_cardst.h"
#include "e_elemnt.h"
#include "e_storag.h"
#include "s_ev.h"
#include "io_matrix.h"
#include "m_matrix_extra.h"
#include "e_aux.h"
using namespace std;
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
namespace { //
class EV : public EV_BASE { //
	public:
		explicit EV(): EV_BASE() {}
		~EV() {}
		void	do_it(CS&, CARD_LIST*);
	private:
		void	setup(CS&);
		explicit EV(const EV&): EV_BASE() {unreachable(); incomplete();}

	protected:
		void	sweep_recursive(int);

		void align_inf_eigenvectors(gsl_matrix_complex_view& EV_inf, size_t inputnumber=0);
		unsigned sens_inf_vectors(gsl_matrix_complex_view& EV);

	private:
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// EV_DISC?
void EV::do_it(CS& Cmd, CARD_LIST* Scope)
{
	trace0("EV::do_it(CS&, CARD_LIST*)");
	_scope = Scope;
	_sim->_time0 = 0.;
	_sim->set_command_ddc();
	_sim->_phase = p_INIT_DC;
	::status.ddc.reset().start();
	_sim->_temp_c = temp_c_in;
	_do_tran_step = 0;
	_dump_matrix = 0;
	command_base(Cmd);
	::status.ddc.stop();
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
void EV::setup(CS& Cmd)
{
	_cont = false;
	_trace = tNONE;
	_out = IO::mstdout;
	_out.reset(); //BUG// don't know why this is needed */
	bool ploton = IO::plotset  &&  plotlist().size() > 0;

	_n_inputs = 0;

	if (Cmd.more()) {
		for (_n_sweeps = 0; Cmd.more() && _n_sweeps < DCNEST; ++_n_sweeps) {
			if (Cmd.umatch("to")) {
				break;
			}
			CARD_LIST::fat_iterator ci = findbranch(Cmd, &CARD_LIST::card_list);
			if (!ci.is_end()) {
				if (ELEMENT* c = dynamic_cast<ELEMENT*>(*ci)) {
					_zap[_n_sweeps] = c;
					trace1("EV::setup", _zap[_n_sweeps]->long_label());
				} else { untested();
					throw Exception("dc/op: can't sweep " + (**ci).long_label() + '\n');
				}
				_inputno[_n_sweeps] = _n_inputs;
				++_n_inputs;
			}else if (Cmd.is_float()) {		// sweep the generator
				_zap[_n_sweeps] = NULL;
			}else if (Cmd.is_alpha()) {
				string pname;
				unsigned here = Cmd.cursor();
				Cmd >> pname;
				try {
					_param[_n_sweeps] = &(_scope->params()->find(pname));
				} catch(Exception) {
					Cmd.reset(here);
					// throw Exception("ddc: can't sweep " + pname + '\n');
				}
				_zap[_n_sweeps] = NULL;
			}else{ untested();
				// leave as it was .. repeat Cmd with no args
			}
			if (Cmd.match1("'\"({") || Cmd.is_float()) {
				_start[_n_sweeps] = "NA";
				_stop[_n_sweeps] = "NA";
				Cmd >> _start[_n_sweeps] >> _stop[_n_sweeps];
				_step[_n_sweeps] = 0.;
			}else{
				// leave it as it was .. repeat Cmd with no args
			}
			_sim->_genout = 0.;
			temp_c_in = OPT::temp_c;
			_sim->_temp_c = temp_c_in;
			options(Cmd,_n_sweeps);

		}
		unsigned here = Cmd.cursor();
		string output;
		bool newprobes = true;
		do{ itested();
			// FIXME: can we use the probe parser?!
			if( Cmd >> "v" ){
				output_t t;
				t.brh[1] = 0;
				trace1("SENS::setup, have v", Cmd.tail());

				int paren = Cmd.skip1b('(');
				Cmd >> output;
				t.label = "v("+output;

				CKT_NODE* node = dynamic_cast<CKT_NODE*>((*_scope).node(output));
				if(node)
					t.brh[0] = node;
				else{
					continue;
				}

				trace2("SENS::setup, have", output, Cmd.tail());
				if(!Cmd.match1(')')){
					output = Cmd.ctos(")");
					trace1("SENS:skip1b:setup, have2", output);
					t.label += "," + output;
					node = dynamic_cast<CKT_NODE*>((*_scope).node(output));
					if(node)
						t.brh[1] = node;
					else{
						Cmd.warn(bWARNING, "probelist: what's this?");
					}

				}else{
					trace0("SENS::setup no comma");
				}

				paren -= Cmd.skip1b(')');
				if (paren != 0) {untested();
					Cmd.warn(bWARNING, "need )");
				}else if (output.empty()) {untested();
					Cmd.warn(bWARNING, "probelist: what's this?");
				}else{
				}

				t.label+=")";
				if (newprobes){
					_output.clear();
					newprobes = false;
				}
				_output.push_back(t);
			}
			//try{
			//  _output.add_list(Cmd);
			//}
			//catch(Exception_Cant_Find)
			//{}
			ONE_OF
				|| Get(Cmd, "dm",	          &_dump_matrix)
				|| _out.outset(Cmd);
			;
		}while (Cmd.more() && !Cmd.stuck(&here));
	}else{ untested();
	}
	Cmd.check(bWARNING, "what's this?");

	// hat peter auskommentiert
	//_sim->_uic = _sim->_more_uic = true;

	IO::plotout = (ploton) ? IO::mstdout : OMSTREAM();
	initio(_out);

	assert(_n_sweeps > 0);
	for (int ii = 0;  ii < _n_sweeps;  ++ii) {
		trace1("DDC_setup", ii);
		_start[ii].e_val(0., _scope);
		fix_args(ii);

		if (_zap[ii]) {
			trace2("zap inc_probes" + _zap[ii]->long_label(), ii, _zap[ii]->probes() );
			_stash[ii] = _zap[ii];			// stash the std value
			_zap[ii]->inc_probes();			// we need to keep track of it

			STORAGE* s = dynamic_cast<STORAGE*>(_zap[ii]);
			if (s) {
				trace2("EV::setup ", _zap[ii]->long_label(), _zap[ii]->is_constant());
				_zap[ii]->set_constant(false);		   // so it will be updated
				_pushel[ii] = _zap[ii];	           // point to value to patch
				_uic_caplist.push_back(s);
				_sweepval[ii] = s->set__ic();
			}else{
				trace1("EV::setup, no STORAGE", ii);
				_zap[ii]->set_value(_zap[ii]->value(),0);  // zap out extensions
				_zap[ii]->set_constant(false);		   // so it will be updated
				_pushel[ii] = _zap[ii];	                   // element to patch
				_sweepval[ii] = _zap[ii]->set__value();
			}
		} else if (_param[ii]) {
			_sweepval[ii] = _param[ii]->pointer_hack();
		} else {
			trace1("EV::setup does not exist", ii);
			//_sweepval[ii] = 0;
			_sweepval[ii] = &_sim->_genout;
			_pushel[ii] = NULL; //&_sim->set_gen;			// point to value to patch
			// throw(Exception("something went wrong\n"));
		}
	}
	_sim->_freq = 0;

}
/*--------------------------------------------------------------------------*/
static EV p2;
static DISPATCHER<CMD>::INSTALL d2(&command_dispatcher, "evs", &p2);
/*--------------------------------------------------------------------------*/
}// namespace
static void register_status()
{
	static bool done;
	if(!done) { itested();
		new DISPATCHER<CKT_BASE>::INSTALL(&status_dispatcher, "ssp", &p2);
		done = true;
	}
}
/*--------------------------------------------------------------------------*/
static size_t n;
void EV_BASE::sweep()
{ itested();
	register_status();

	trace0("EV_BASE::sweep()");

	_sim->set_command_tran();
	head(_start[0], _stop[0], "chart  spl"); // here? hmmm
	//  _sim->_bypass_ok = false;
	_sim->set_inc_mode_bad();
	if (_cont) {untested();
		_sim->restore_voltages();
		CARD_LIST::card_list.tr_restore();
	}else{
		_sim->clear_limit();
		CARD_LIST::card_list.tr_begin();
	}

// 	allocate()?
	n = _sim->_total_nodes;
	// FIXME: allocate globally, where necessary.
	G_ = gsl_matrix_alloc(n,n);
	C_ = gsl_matrix_alloc(n,n);
	Q_ = gsl_matrix_alloc(n,n);
	Z_ = gsl_matrix_alloc(n,n);

	sweep_recursive(_n_sweeps);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
double	EV_BASE::temp_c_in = 0.;
/*--------------------------------------------------------------------------*/
	EV_BASE::EV_BASE()
:SIM(),
	_n_sweeps(1),
	_cont(false),
	_trace(tNONE)
{

	for (int ii = 0; ii < DCNEST; ++ii) {
		_loop[ii] = false;
		_reverse_in[ii] = false;
		_reverse[ii] = false;
		_step[ii]=0.;
		_linswp[ii]=true;
		_sweepval[ii]=&_sim->_genout;
		_zap[ii]=NULL;
		_stepmode[ii] = ONE_PT;
		_param[ii]=NULL;
	}

	temp_c_in=OPT::temp_c;
	_out=IO::mstdout;
}
/*--------------------------------------------------------------------------*/
void EV_BASE::finish(void)
{
	trace0("EV_BASE::finish");
	//  _sim->_phase = p_;
	for (int ii = 0;  ii < _n_sweeps;  ++ii) {
		if (_zap[ii]) { // component
			_stash[ii].restore();
			_zap[ii]->dec_probes();
			trace2("EV_BASE::finish dec_probes done", _zap[ii]->long_label(), _zap[ii]->probes());
			_zap[ii]->precalc_first();
			_zap[ii]->precalc_last();
		}else{
		}
	}
	_uic_caplist.clear();
}
/*--------------------------------------------------------------------------*/
void EV_BASE::fix_args(int Nest)
{
	trace1("EV_BASE::fix_args(int Nest)", Nest);

	_stop[Nest].e_val(_start[Nest], _scope);
	_step_in[Nest].e_val(0., _scope);
	_step[Nest] = _step_in[Nest];

	switch (_stepmode[Nest]) { untested();
		case ONE_PT:
		case LIN_STEP:
			_linswp[Nest] = true;
			break;
		case LIN_PTS:untested();
						 if (_step[Nest] <= 2.) {untested();
							 _step[Nest] = 2.;
						 }else{untested();
						 }
						 _linswp[Nest] = true;
						 break;
		case TIMES:untested();
					  if (_step[Nest] == 0.  &&  _start[Nest] != 0.) {untested();
						  _step[Nest] = _stop[Nest] / _start[Nest];
					  }else{untested();
					  }
					  _linswp[Nest] = false;
					  break;
		case OCTAVE:
					  if (_step[Nest] == 0.) {untested();
						  _step[Nest] = 1.;
					  }else{ untested();
					  }
					  _step[Nest] = pow(2.00000001, 1./_step[Nest]);
					  _linswp[Nest] = false;
					  break;
		case DECADE:
					  if (_step[Nest] == 0.) { untested();
						  _step[Nest] = 1.;
					  }else{ untested();
					  }
					  _step[Nest] = pow(10., 1./_step[Nest]);
					  _linswp[Nest] = false;
					  break;
	};

	if (_step[Nest] == 0.) {	// prohibit log sweep from 0
		_step[Nest] = _stop[Nest] - _start[Nest];
		_linswp[Nest] = true;
	}else{
	}
}
/*--------------------------------------------------------------------------*/
// process tail vector.
// find maximum coeff in v0prime_tail at piv, swap to 0, and
// exchange, rescale first basis vector to get z0_tail = (1,0,0 ...)
// maybe unnecessary?
int EV_BASE::process_tail(gsl_matrix_complex_view& relevant_EV)
{
	int err = 0;
	size_t n = _sim->_total_nodes;
	unsigned tailstart = unsigned(relevant_EV.matrix.size2);
	unsigned tailsize = unsigned(n-tailstart);
	trace2("preparing tail", tailstart, tailsize);
	if(tailstart < n) {
		gsl_matrix_complex_view EV_tail = gsl_matrix_complex_submatrix (EV_, 0, tailstart, n, tailsize);
		gsl_vector_complex_view v0prime_tail = gsl_vector_complex_subvector(z0_, tailstart, tailsize);
		size_t piv = gsl_blas_izamax(&v0prime_tail.vector);
		if(piv != 0){
			gsl_vector_complex_swap_elements(&v0prime_tail.vector,0,piv);
			gsl_matrix_complex_swap_columns(&EV_tail.matrix,0,piv);
		}
		gsl_complex pivot = gsl_vector_complex_get(&v0prime_tail.vector,0);
		if (0.!=GSL_REAL(pivot) || 0.!=GSL_IMAG(pivot)){
			gsl_vector_complex_view first = gsl_matrix_complex_column(EV_, tailstart);

			for(unsigned i=tailstart+1; i<n; ++i) {
				gsl_complex s = gsl_vector_complex_get(z0_,i);
				gsl_vector_complex_const_view b = gsl_matrix_complex_const_column(EV_, i);
				gsl_blas_zaxpy(gsl_complex_div(s,pivot), &b.vector, &first.vector);

				gsl_vector_complex_set(z0_,i,zero);
			}

			double norm = 1;
			gsl_complex* p = gsl_vector_complex_ptr(&v0prime_tail.vector,0);
			if(0){
				norm = 1/gsl_blas_dznrm2(&first.vector);
			}else if(1) {
				norm = gsl_complex_abs(*p);
			}
			trace1("normalize tail[0]", norm);
			err+= gsl_vector_complex_scale(&first.vector, gsl_complex_mul(one,gsl_complex_rect(norm,0.))); // inefficient.
			GSL_REAL(*p)/= norm;
			GSL_IMAG(*p)/= norm;
			trace1("tail", first);
		} else { itested();
			trace2("no pivot. possible zero op", tailstart, tailsize);
		}
	}else{
		trace2("not preparing tail. possible zero op", tailstart, tailsize);
	}
	return err;
}
/*--------------------------------------------------------------------------*/
// match _sweepval and dcsens coeffs.
inline void EV_BASE::compute_canonical_op(gsl_vector_complex_view EV_op, gsl_matrix_complex_const_view dcsens)
{
	trace1("compute_canonical_op", *z0_);
	for(unsigned i=0; i<_n_inputs; ++i) {
		gsl_complex* k_i = gsl_vector_complex_ptr(z0_, _n_states + i);
		gsl_complex delta = gsl_complex_sub(*k_i, gsl_complex_rect(*_sweepval[_inputno[i]],0.));
		trace3("inputs vs inputs", *k_i, *_sweepval[_inputno[i]], delta);

		*k_i = gsl_complex_rect(*_sweepval[_inputno[i]], 0.);

		// op -= k_i * delta;
		gsl_vector_complex_const_view s = gsl_matrix_complex_const_column(&dcsens.matrix,i);
		gsl_blas_zaxpy( delta, &s.vector, &EV_op.vector);
	}

//	for(unsigned i=0; i<_n_inputs; ++i) {
//		gsl_complex* op_k = gsl_vector_complex_ptr(z0_, _n_states + i);
//		trace2("inputs vs inputs", *op_k, *_sweepval[_inputno[i]]);
//	}
}
/*--------------------------------------------------------------------------*/
// foreach input
//   compute node voltage DC sensitivities, stash in columns of dc_sens
//   project to infspace
//   fill inf vectors, preserving linear independence.
//
//   append op to ensure v0 is in img. FIXME
//
inline unsigned EV_BASE::sens_inf_vectors(gsl_matrix_complex_view EV_inf, gsl_vector_view grid, gsl_matrix_complex* nodc_sens)
{
	size_t n = _sim->_total_nodes;
	assert(n==EV_inf.matrix.size1); // number of rows
	size_t m = EV_inf.matrix.size2;
	unsigned seek = 0;
	trace4("sens_inf_vectors", EV_inf.matrix.size1, _n_sweeps, n, EV_inf.matrix.size2);
	trace1("sens_inf_vectors", grid);

	for(unsigned i=0; i<unsigned(_n_sweeps); ++i){ itested();
		if(!_zap[i]){ untested();
			// sweeping over parameter?
			continue;
		}

		clear_sens();
		_zap[i]->sens_load("dummy"); // sens input
		complex_array_to_gsl_vector_complex(*sens_,_sim->_sens+1,n);
		trace3("sensstamp", _zap[i]->long_label(), *sens_, _step[i]);

		::status.back.start();
		_sim->_lu.fbsub(_sim->_sens);
		::status.back.stop();

		trace1("input sens\n", _sim->_aa);

		complex_array_to_gsl_vector_complex(*sens_,_sim->_sens+1,n);
		gsl_vector_complex_view nodcsens = gsl_matrix_complex_column(nodc_sens, seek);

		trace2("DC input sens", *sens_, *betabar_);

//		gsl_vector_complex_memcpy(&dcsens.vector, sens_);

		bool inf_sens = true;
		if(inf_sens){ incomplete();
			// hmm, no-DC sensitivity...
			gsl_linalg_complex_LU_svx(EVLU_, LU_perm, sens_); // too many, only need left ones...
			trace3("sensprime\n", *sens_, m,n);
			gsl_vector_complex_view sens_right = gsl_vector_complex_subvector(sens_, n-m, m);

			gsl_blas_zgemv( CblasNoTrans, one, &EV_inf.matrix, &sens_right.vector, zero, ztmp_ );

			if(1){ // sanitycheck
				mul(CblasTrans, *C_, *ztmp_, *sens_);
				trace1("zero", *sens_);
			}
			trace1("sens projection", *ztmp_);
			gsl_vector_complex_memcpy(&nodcsens.vector, ztmp_);
		}else{
		}

		complex_array_to_gsl_vector_complex(*sens_,_sim->_sens+1,n);

		if(1){ //emplace...
			gsl_vector_complex_view seek_ = gsl_matrix_complex_column(&EV_inf.matrix, seek);
			gsl_vector_complex_memcpy(&seek_.vector, sens_);
		}else{ // rewrite basis... inefficient
		  // BUG. refresh dc_sens!!
			gsl_vector_complex_const_view seek_ = gsl_matrix_complex_const_column(&EV_inf.matrix, seek);
			gsl_vector_complex_sub(sens_, &seek_.vector);
			if (0) for(size_t j=seek; j<m; ++j){
				trace1("replacing", j);
				gsl_vector_complex_view seek_ = gsl_matrix_complex_column(&EV_inf.matrix, j);
				gsl_vector_complex_add(&seek_.vector, sens_);

				// normalize
				gsl_complex len;
				GSL_IMAG(len) = 0.;
				GSL_REAL(len) = 1./norm(seek_.vector);
				gsl_vector_complex_scale(&seek_.vector, len);
			}
			// BUG. should use d z0dot/din as grid size!!
			// FIXME: use chart!!1
		}
		{
			double* g = gsl_vector_ptr(&grid.vector,seek);
			if (fabs(_step[i]) < fabs(*g)){ itested();
				*g = fabs(_step[i]);
		//		trace2("sens_inf_vectors", _step[i], c.grid(i)); // FIXME
			}
		}
		++seek;
	}
	trace1("done sens_inf_vectors\n", EV_inf.matrix);
	return seek;
}
/*--------------------------------------------------------------------------*/
gsl_complex operator-(const gsl_complex&x){
	return gsl_complex_negative(x);
}
/*--------------------------------------------------------------------------*/
SSP_STATE& EV_BASE::ev_solver(SSP_CHART*& chart)
{

	// FIXME: fetch more information from chart.
	//  - number of inputs
	//  - number of states.
	//  - grid stuff

	int err;
	size_t n = _sim->_total_nodes;

	trace2("========= ev_solver ========\n", _sim->_acx, SSP_TREE::num_charts());
	assert(SSP_TREE::num_charts() < (1<<20));
	//  gsl_eigen_gen_workspace* w = gsl_eigen_gen_alloc(n);
	gsl_eigen_genv_workspace* w = gsl_eigen_genv_alloc(n);
	gsl_vector *tmp_;
	gsl_vector_complex *v0_;
	gsl_vector *grid_;
	gsl_vector_view v0__;
	gsl_vector_complex *v0back_;

	z0_ = gsl_vector_complex_alloc(n);
	z0dot_ = gsl_vector_complex_alloc(n);
	v0back_ = gsl_vector_complex_alloc(n);

	BSMATRIX<double> G = _sim->_acx.real();
	BSMATRIX<double> C = _sim->_acx.imag();


	// P_ = gsl_matrix_alloc(n,n);
	EV_ = gsl_matrix_complex_alloc(n,n);
	EVLU_ = gsl_matrix_complex_alloc(n,n);
	LU_perm = gsl_permutation_alloc (n);

	alpha_ = gsl_vector_complex_alloc(n);
	beta_ = gsl_vector_alloc(n);
	betabar_ = gsl_vector_calloc(n);
	v0_ = gsl_vector_complex_alloc(n);
	sens_ = gsl_vector_complex_alloc(n);
	grid_ = gsl_vector_alloc(n);
	tmp_ = gsl_vector_alloc(n);
	ztmp_ = gsl_vector_complex_alloc(n);
	lambda_ = gsl_vector_complex_alloc(n); // fixme: only need maximum _n_states
	// _n_states may depend on operating point...
	G >> *G_;
	C >> *C_;
	v0__ = gsl_vector_view_array(_sim->_v0+1,n);
	trace1("v0", v0__.vector);

	trace1("A = G =\n", *G_);
	trace1("B = C =\n", *C_);

	// want \lambda C x = G x
	// or   1/\mu   C x = G x
	// solve beta G = alpha C
	// lambda_i = alpha/beta
	// mu_i     = beta/alpha

	err+= gsl_eigen_genv_QZ (G_, C_, alpha_, beta_, EV_, Q_, Z_, w);
	trace1("raw genv\n", *EV_);
	trace1("schur\n", *G_);
	trace1("schur\n", *C_);

	{ // check.
		err+= gsl_vector_complex_memcpy(lambda_,alpha_);
		*lambda_ /= *beta_;
	}
	trace2("solved alpha G = beta C\n", *alpha_, *beta_);
	trace1("lambda\n", *lambda_);

	C >> *C_;

	//gsl_eigen_gen_QZ (G_, C_, alpha_, beta_, Q_, Z_, w);

	if(0){ // hmm check
		gsl_matrix *test_;
		test_ = gsl_matrix_alloc(n,n);
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., Q_, Z_, 0., test_);
		trace1("QZ^T?", *test_);
	}
	//  { // hmm check
	//    gsl_matrix_complex *test_;
	//    test_ = gsl_matrix_complex_alloc(n,n);
	//    gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, one, EV_, EV_, zero, test_);
	//    trace1("identity?", *test_);
	//  }
	//
	to_norm(*grid_, *alpha_);
	gsl_vector_div(grid_, beta_);
	trace1("lambda tmp grid", *grid_);
	trace1("...", *lambda_);

	gsl_permutation *ev_perm = gsl_permutation_alloc (n);

	// order_by_lambda(ev_perm, grid_ find_dim)
	{ // order eigenvalues by size (starting with small)
		//
		// weed out huge lambdas.
		// replace 1st nan eigenvector by v0 ... hmm
		// apply order to EV_
		size_t fin_dim = (size_t)n;

		err+= sort_vector_and_index (ev_perm, grid_, fin_dim);
		trace2("evalue perm", *ev_perm, fin_dim);

		for(unsigned i=0; i<fin_dim; ++i){ itested();
			gsl_vector_set(betabar_,i,1);
		}

		trace1("beta", *betabar_);

		for(unsigned i=0; i<n; ++i){
			gsl_vector_complex_view v;
			gsl_vector_view vr;
			v = gsl_matrix_complex_row (EV_, i);
			err+= gsl_permute_vector_complex (ev_perm, &v.vector);
			vr = gsl_matrix_row (Q_, i);
			err+= gsl_permute_vector (ev_perm, &vr.vector);
			vr = gsl_matrix_row (Z_, i);
			err+= gsl_permute_vector (ev_perm, &vr.vector);
		}
		trace1("ordered\n", *EV_);
		if (0) {
			gsl_matrix_complex *test_;
			test_ = gsl_matrix_complex_alloc(n,n);
			gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, one, EV_, EV_, zero, test_);
			trace1("ordered, identity?", *test_);
		}
		// first n rows of EV_ (the non inf-EVs)


		real_array_to_gsl_vector(*z0_,_sim->_v0+1,n);
		trace1("v0", *z0_);
		//     {
		//       trace1("preimg1", *EV_);
		//       err+= gsl_blas_zgemv( CblasConjTrans, one, EV_, z0_, zero, ztmp_ );
		//       trace1("preimg1", *ztmp_);
		//     }

		gsl_matrix_complex_memcpy (EVLU_, EV_);

		int signum;
		gsl_linalg_complex_LU_decomp(EVLU_, LU_perm, &signum);

		// probably not necessary.
		// compute EV^-1 v0 below, with modified EV
		real_array_to_gsl_vector(*z0_,_sim->_v0+1,n);
		gsl_linalg_complex_LU_svx(EVLU_, LU_perm, z0_);
		trace2("preimg2  ", *z0_, signum);

		// reverse test
		err+= gsl_blas_zgemv( CblasNoTrans, one, EV_, z0_, zero, ztmp_ );
		trace1("v0 again?", *ztmp_);
		_n_states = (unsigned)fin_dim;
	}

	if(0){ //trying to figure out coeffs...
		err+= gsl_blas_dgemv( CblasTrans, 1, Z_, &v0__.vector, 0, tmp_ );
		err+= gsl_permute_vector (ev_perm, tmp_);
		trace1("preimg3 ", *tmp_);
		set_vector(*z0_,*tmp_);
		err+= gsl_blas_zgemv( CblasNoTrans, one, EV_, z0_, zero, ztmp_ );
		trace1("v0 again?", *ztmp_);
	}

	::status.lud.start();
	COMPLEX z = COMPLEX(OPT::gmin,0.);
	_sim->_acx.dezero(z);
	_sim->_acx.lu_decomp();
	::status.lud.stop();

	trace2("hmm", n, _n_inputs);
	real_array_to_gsl_vector(*v0_,_sim->_v0+1,n);

	// fin_dim:   dimension of noninf space.
	// input_dim: dimension of constant Z input space
	//            == number of inputs (?)

	unsigned input_dim;
	gsl_matrix_complex* nodc_sens_ = gsl_matrix_complex_alloc(n, _n_inputs);
	{ itested();
		trace1("pre input sens", *grid_);
		// compute input sensitivities.
		// add to basis ?
		gsl_matrix_complex_view EV_inf = gsl_matrix_complex_submatrix (EV_, 0, _n_states, n, n-_n_states);
		gsl_vector_view g = gsl_vector_subvector (grid_, _n_states, n-_n_states);
		trace3("infspace\n", EV_inf.matrix, EV_inf.matrix.size1, EV_inf.matrix.size2);
		input_dim = sens_inf_vectors(EV_inf, g, nodc_sens_);
		trace1("modified EV", input_dim);
		trace1("done sens_inf_vectors", *grid_);
		trace1("done sens_inf_vectors", *nodc_sens_);
	}
	assert(input_dim==unsigned(_n_inputs));

	trace1("tmp grid", *grid_);

	{  // recompute LU decomp. inefficient!
		// why?
		int signum;
		gsl_matrix_complex_memcpy (EVLU_, EV_);
		gsl_linalg_complex_LU_decomp(EVLU_, LU_perm, &signum);
	}

	{ // v0' = EV^-1 v0
		real_array_to_gsl_vector(*z0_,_sim->_v0+1,n);
		trace1("v0", *z0_);
		gsl_linalg_complex_LU_svx(EVLU_, LU_perm, z0_);
		trace1("", *z0_);
	}

	{
		gsl_matrix_complex_view relevant_EV;
		relevant_EV = gsl_matrix_complex_submatrix (EV_, 0, 0, n, _n_states + input_dim);
		trace1("unprocessed tail\n", *EV_);
		trace1("unprocessed tail", *z0_);
		process_tail(relevant_EV);
		trace1("processed tail\n", *EV_);
		trace1("processed tail", *z0_);
	}

//	gsl_vector_complex* z0limit_ = gsl_vector_complex_alloc(n);
#if 0
	if(_n_states) { // new grid?
		trace2("new grid", _n_states, input_dim);
		gsl_vector_complex_memcpy(z0dot_, z0_);
		gsl_vector_complex_view z0dot__ = gsl_vector_complex_subvector(z0dot_, 0, _n_states);
		do_ddc(z0dot__);
		gsl_vector_complex_memcpy(z0limit_, ztmp_);
	}else{ // no movement, already at limit
		gsl_vector_complex_memcpy(z0limit_, z0_);
	}
#endif

	trace1("precan op", *z0_);
	trace1("precan op\n" , *EV_);
	if(_n_states+_n_inputs < n) { itested();
		gsl_vector_complex_view op__ = gsl_matrix_complex_column(EV_, _n_states+_n_inputs);
		gsl_matrix_complex_const_view dcs = gsl_matrix_complex_const_submatrix (EV_, 0, _n_states, n, n -_n_states);
		trace2("can op EV", _n_states, _n_inputs);
		trace3("can op EV", _n_states, n, n-_n_states);
		compute_canonical_op(op__, dcs);
		trace1("can op", *z0_);
		trace1("can op", op__);
		trace1("can op\n", *EV_);
	}

	if(_n_states) { itested();
		gsl_vector_complex_view z0dot = gsl_vector_complex_subvector(z0dot_,0,_n_states);
		gsl_vector_complex_const_view z0 = gsl_vector_complex_const_subvector(z0_, 0, _n_states);
		compute_zdot(z0, z0dot);
		trace1("z0dot", z0dot);
	} else { untested();
	}

	if(0 && _n_states){
		gsl_vector_view g = gsl_vector_subvector (grid_, 0, _n_states);
		trace1("scale sens", *grid_);
		gsl_vector_scale(g, OPT::dtmin);
		trace1("scale sens", *grid_);
	}

	{
		gsl_matrix_complex_view relevant_EV;
		trace1("pre output sens", *grid_);
		trace2("pre output sens", _n_states, input_dim);
		relevant_EV = gsl_matrix_complex_submatrix (EV_, 0, 0, n, min(unsigned(n), _n_states + input_dim + 1)); // +1??
		output_sens_grid(relevant_EV, grid_, lambda_);
		trace1("done output sens", *grid_);
	}

	{  // recompute LU decomp for rhs. inefficient...
		int signum;
		gsl_matrix_complex_memcpy (EVLU_, EV_);
		gsl_linalg_complex_LU_decomp(EVLU_, LU_perm, &signum);
	}


	trace1("unscaled\n", *EV_);
	trace3("actual grid", *grid_, _n_states, input_dim);
	to_grid(*grid_, OPT::reltol * OPT::reltol);
	trace1("disc grid", *grid_);

	if(_n_states+_n_inputs < n) {
		// try to stabilize canonical op.
		// doesn't work
		gsl_vector_set(grid_,_n_states+_n_inputs,OPT::reltol);
	}

	chart = &SSP_TREE::_root.insert(*EV_, *grid_, _n_inputs);
	assert(chart);

	assert(chart->n_inputs() == _n_inputs);

	if(chart->is_new()){ itested();
		for(unsigned i=0; i<_n_inputs; ++i){
			chart->_inputgrid[i] = gsl_vector_get(grid_,i+_n_states);
		}
	}else{ itested();
		for(unsigned i=0; i<_n_inputs; ++i){
			assert(chart->grid(i) == gsl_vector_get(grid_,i+_n_states));
		}
	}

	assert(chart);
	trace1("discretized\n", *EV_);
	trace3("preminimize", _n_states, input_dim, *z0_);

	trace1("v0 raw  ", v0__.vector);
	{ // find z0, | v0 - EV_grid z0 | minimal.
		size_t numcols = _n_states+input_dim;
		bool have_op = numcols != n;
		double rcond = 0;
		integer rank = 0;
		gsl_complex opcoeff;
		integer info = 0;

		if(have_op){ untested();
			// hmm opcoeff could be one, always. but it's not (currently)
			opcoeff = gsl_vector_complex_get(z0_, _n_states+input_dim);
		}

		real_array_to_gsl_vector(*z0_,_sim->_v0+1,n);

		if(have_op){ untested();
			gsl_vector_complex_view op__ = gsl_matrix_complex_column(EV_, _n_states+_n_inputs);
			trace2("sub op", op__, opcoeff);
//			gsl_vector_complex_sub(z0_, &op__.vector);
			gsl_blas_zaxpy(-opcoeff, &op__.vector, z0_);
			trace1("subd op", *z0_);
		}

		gsl_matrix_complex_view v0pm = gsl_matrix_complex_view_vector (z0_, 1, n);
		trace1("zgelss", v0pm.matrix);

		// use EVLU as scratch space
		err+= gsl_matrix_complex_memcpy (EVLU_, EV_);
		err+= gsl_matrix_complex_transpose(EVLU_);
		gsl_matrix_complex_view zEVLU = gsl_matrix_complex_submatrix(EVLU_,0,0,numcols,n);
		trace1("zgelss", zEVLU.matrix);
		gsl_blas_zgelss (&zEVLU.matrix, &v0pm.matrix, tmp_, &rcond, &rank, &info );
		trace3("zgelss result", info, rank, *z0_);
		trace1("zgelss result", v0pm.matrix);
		if(have_op){ untested();
			gsl_vector_complex_set(z0_,_n_states+input_dim, opcoeff);
		}
		trace1("zgelss result+op", *z0_);
		err+= gsl_blas_zgemv( CblasNoTrans, one, EV_, z0_, zero, v0back_ );
		trace2("zgelss EV_*z0_", *v0back_, info);
	}

	trace1("hmmm", *z0_);
	trace1("hmmm", *grid_);

	for(unsigned i=0; i<_n_states; ++i){
		gsl_complex dotgrid = gsl_vector_complex_get(z0dot_,i);
		double on = 1./norm(dotgrid);
		double ln = gsl_vector_get(grid_, i);
		if(on<ln){ incomplete();
			// gsl_vector_set(grid_,i,on);
		}
	}

	{ // refine state grid. using nodc-sens
		for(unsigned i=0; i<_n_states; ++i){
			gsl_vector_complex_const_view state = gsl_matrix_complex_const_column(EV_, i);
			for(unsigned j=0; j<_n_inputs; ++j){
				gsl_vector_complex_const_view input = gsl_matrix_complex_const_column(EV_, _n_states + j);
				gsl_complex cpl;
				gsl_blas_zdotc (&state.vector, &input.vector, &cpl);

				trace2("nodc cpl", input, state);
				trace3("nodc cpl", i, j , cpl);
				trace1("nodc cpl", _step[j]);

				double cpln = norm(cpl);

				// set state grid to <= inputgrid/cpln
				if(cpln * gsl_vector_get(grid_,i) > gsl_vector_get(grid_,j+_n_inputs)){ itested();
					trace2("nodc reducing grid", i, gsl_vector_get(grid_,j+_n_inputs)/cpln);
					gsl_vector_set(grid_,i,gsl_vector_get(grid_,j+_n_inputs)/cpln);
				}else{
					trace0("nodc not reducing grid");
				}
			}
		}
	}


	// 	unsigned gridbase=2; not yet.
	index_t state_index[_n_states];
	if(_n_states){ itested();
		//		gsl_vector_complex_view z0dot = gsl_vector_complex_subvector(z0dot_,0,_n_states);
		gsl_vector_const_view tol = gsl_vector_const_subvector(grid_,0,_n_states);
		gsl_vector_complex_view z0 = gsl_vector_complex_subvector(z0_, 0, _n_states);
		gsl_vector_complex_const_view lambda = gsl_vector_complex_const_subvector(lambda_, 0, _n_states);
		trace2("quant tol", tol, lambda);
		trace1("quant input", z0);
		if (_quantize_states) { itested();
			quantize(z0, lambda, tol, state_index);
		}else{ untested();
			to_grid(z0, tol.vector, state_index);
		}
	}

	trace1("refined", *grid_);
	to_grid(*z0_, *grid_);
	trace1("discrete", *z0_);
	trace1("discretized\n", *EV_);
	err+= gsl_blas_zgemv( CblasNoTrans, one, EV_, z0_, zero, v0back_ );
	trace1("v0back", *v0back_);

	for(unsigned i=0; i<n; ++i){
		_sim->_v0[i+1] = GSL_REAL(*gsl_vector_complex_const_ptr(v0back_,i));
	}

	// FIXME: merge insert/discretize. use unit grid
	SSP_STATE* state;
	if(_n_states) { untested();
		gsl_vector_view r = gsl_vector_complex_real(z0_);
		gsl_vector_view i = gsl_vector_complex_imag(z0_);
		gsl_vector_div(&r.vector, grid_);
		gsl_vector_div(&i.vector, grid_);

		state = &chart->insert_state(state_index, _n_states);

	}else{ itested();
		state = &chart->insert();
	}
	USE(state);
	trace1("inserted", chart->size());

	trace1("========= ev_solver done ========\n", _sim->_acx);

	// BUG use global alloc/free the other stuff!
	gsl_eigen_genv_free(w);
	return *state;
}
/*--------------------------------------------------------------------------*/
void EV_BASE::options(CS& Cmd, int Nest)
{
	trace1("EV_BASE::options(CS&, int)", Nest);

	_quantize_states = true;
	_loop[Nest] = _reverse_in[Nest] = false;
	_sim->_uic = false;
	_tran_step = OPT::dtddc;
	unsigned here = Cmd.cursor();
	do{
		ONE_OF
			|| (Cmd.match1("'\"({")	&& ((Cmd >> _step_in[Nest]), (_stepmode[Nest] = LIN_STEP)))
			|| (Cmd.is_float()	&& ((Cmd >> _step_in[Nest]), (_stepmode[Nest] = LIN_STEP)))
			|| (Get(Cmd, "*",		  &_step_in[Nest]) && (_stepmode[Nest] = TIMES))
			|| (Get(Cmd, "+",		  &_step_in[Nest]) && (_stepmode[Nest] = LIN_STEP))
			|| (Get(Cmd, "by",	  &_step_in[Nest]) && (_stepmode[Nest] = LIN_STEP))
			|| (Get(Cmd, "trstep",	  &_tran_step))
			|| (Get(Cmd, "uic",	  &_uic))
			|| (Get(Cmd, "qu{antize}",	  &_quantize_states))
			|| (Get(Cmd, "step",	  &_step_in[Nest]) && (_stepmode[Nest] = LIN_STEP))
			|| (Get(Cmd, "d{ecade}",	  &_step_in[Nest]) && (_stepmode[Nest] = DECADE))
			|| (Get(Cmd, "ti{mes}",	  &_step_in[Nest]) && (_stepmode[Nest] = TIMES))
			|| (Get(Cmd, "lin",	  &_step_in[Nest]) && (_stepmode[Nest] = LIN_PTS))
			|| (Get(Cmd, "o{ctave}",	  &_step_in[Nest]) && (_stepmode[Nest] = OCTAVE))
			|| Get(Cmd, "c{ontinue}",   &_cont)
			|| Get(Cmd, "tr{s}",        &_do_tran_step)
			|| Get(Cmd, "sr",           &_slew_rate)
			|| Get(Cmd, "dm",           &_dump_matrix)
			|| Get(Cmd, "dt{emp}",	  &temp_c_in,   mOFFSET, OPT::temp_c)
			|| Get(Cmd, "lo{op}", 	  &_loop[Nest])
			|| Get(Cmd, "re{verse}",	  &_reverse_in[Nest])
			|| Get(Cmd, "te{mperature}",&temp_c_in)
			|| (Cmd.umatch("tr{ace} {=}") &&
					(ONE_OF
					 || Set(Cmd, "n{one}",      &_trace, tNONE)
					 || Set(Cmd, "o{ff}",       &_trace, tNONE)
					 || Set(Cmd, "w{arnings}",  &_trace, tUNDER)
					 || Set(Cmd, "i{terations}",&_trace, tITERATION)
					 || Set(Cmd, "v{erbose}",   &_trace, tVERBOSE)
					 || Cmd.warn(bWARNING, 
						 "need none, off, warnings, iterations, verbose")
					)
				)
			|| _out.outset(Cmd);
	}while (Cmd.more() && !Cmd.stuck(&here));

}
/*--------------------------------------------------------------------------*/
void EV_BASE::set_uic_caps_constant(bool x)
{
	for( vector<STORAGE*>::iterator i=_uic_caplist.begin(); i!=_uic_caplist.end(); ++i)
	{
		(*i)->set_constant(x);
	}
}
/*--------------------------------------------------------------------------*/
// void EV_BASE::compute_op(OPT::ITL itl)
void EV_BASE::get_op(OPT::ITL itl)
{ itested();
	_sim->_time0 = _sim->_dt0 = 0.0;
	_sim->tr_reset();
	_sim->_phase = p_INIT_DC; // BUG?

	_sim->_loadq.clear();
	_sim->_evalq->clear();
	_sim->_evalq_uc->clear();
	_sim->set_inc_mode_bad();

	TRACE oldtrace=_trace;
	_trace=tNONE;
	int converged = solve_with_homotopy(itl,_trace);
	_trace=oldtrace;

	size_t d = _sim->_total_nodes;


	if (!converged) {itested();
		error(bWARNING, "did not converge\n");
		throw(Exception("foobar"));
	}
	::status.accept.start();
	_sim->set_limit();

	if(_dump_matrix){
		_out << " ======================== \n";
		_out << "_i ( " << _sim->_i[1];
		for(unsigned a=2; a <= d; ++a){
			_out << " " <<  _sim->_i[a];
		}
		_out  << ") \n";
	}

	_sim->_loadq.clear();

	_sim->_uic = false;
	CARD_LIST::card_list.tr_accept(); // hmmm.
	CARD_LIST::card_list.do_tr();     // no-uic evaluates caps
	_sim->_phase = p_INIT_DC;
	_sim->_uic = false;
	_sim->_evalq->clear();
	_sim->_evalq_uc->clear();
	_sim->count_iterations(iTOTAL);
	::status.accept.stop();
	rhs_ = gsl_vector_alloc(d);
	array_to_gsl_vector(*rhs_,_sim->_i+1,d);

#if 0 // HACK, somehow use inc mode?!
	_sim->set_inc_mode_no();
	clear_arrays();
	CARD_LIST::card_list.tr_load();
#else
	while (!_sim->_loadq.empty()) {
		trace1("loading from q", _sim->_loadq.back()->long_label());
		_sim->_loadq.back()->tr_load();
		_sim->_loadq.pop_back();
	}
#endif
	rhs_ = gsl_vector_alloc(d);
	array_to_gsl_vector(*rhs_,_sim->_i+1,d);
	::status.lud.start();
	_sim->_lu.lu_decomp(_sim->_aa); // , bool(OPT::lubypass && _sim->is_inc_mode()));
	::status.lud.stop();


	trace1("fresh", *rhs_);
	trace1("fresh\n", _sim->_aa);
	_sim->_lu.fbsub(_sim->_i);
	array_to_gsl_vector(*rhs_,_sim->_i+1,d);
	trace1("G^-1 rhs", *rhs_);

	// ???
	finish_building_evalq();
	evaluate_models();
	//  CARD_LIST::card_list.tr_load();
	//   array_to_gsl_vector(*rhs_,_sim->_i+1,d);

	_sim->keep_voltages(); // vdc = v0

	_sim->_acx.reallocate(); // ?
	_sim->_jomega = COMPLEX(0., 1.);
	_sim->_phase = p_AC;
	CARD_LIST::card_list.ac_begin();

	if (_dump_matrix) {
		_out << "solved w/ht\n";
		_out << "i ( " << _sim->_i[1]; // K-put

		for(unsigned a=2; a <= _sim->_total_nodes; ++a){
			_out << " " <<  _sim->_i[a];
		}
		_out  << ") \n";
		_out << "v0 = ( " << _sim->_v0[1];
		for(unsigned a=2;a <= _sim->_total_nodes; ++a){
			_out << " " <<  _sim->_v0[a];
		}
		_out << ") \n";
	}

	_sim->_uic = false;
	_sim->set_inc_mode_yes();
	while (!_sim->_loadq.empty()) {
		trace1("loading from q", _sim->_loadq.back()->long_label());
		_sim->_loadq.back()->tr_load();
		_sim->_loadq.pop_back();
	}
}
/*--------------------------------------------------------------------------*/
void EV_BASE::get_op_()
{
	_sim->_time0 = 0.;
	_sim->set_inc_mode_no(); // uaargh
	_sim->_uic = false;
	_sim->_phase = p_RESTORE;
	size_t n = _sim->_total_nodes;

	gsl_vector_view rhs__ = gsl_vector_view_array(_sim->_i+1,n);
	trace1("get_op b4", rhs__);
	clear_arrays();

	// CARD_LIST::card_list.tr_accept(); // hmmm.
	CARD_LIST::card_list.tr_begin(); // hmmm.
	// finish_building_evalq();
	CARD_LIST::card_list.do_tr();

#if 1 // hack
	CARD_LIST::card_list.tr_load();
#else

	while (!_sim->_loadq.empty()) {
		trace1("get_op loading from q", _sim->_loadq.back()->long_label());
		_sim->_loadq.back()->tr_load();
		_sim->_loadq.pop_back();
	}
#endif

	assert(rhs_->size == n);
	array_to_gsl_vector(*rhs_,_sim->_i+1,n);
	::status.lud.start();
	_sim->_lu.lu_decomp(_sim->_aa); // , bool(OPT::lubypass && _sim->is_inc_mode()));
	::status.lud.stop();

	trace1("get_op fresh", *rhs_);
	trace1("get_op fresh\n", _sim->_aa);
	_sim->_lu.fbsub(_sim->_i);
	array_to_gsl_vector(*rhs_,_sim->_i+1,n);
	trace1("G^-1 rhs", *rhs_);

	{ // also fetch AC stuff
		_sim->_acx.zero();
		std::fill_n(_sim->_ac, n+1, 0.);
		::status.load.start();
		_sim->count_iterations(iTOTAL);
		CARD_LIST::card_list.do_ac();
		CARD_LIST::card_list.ac_load();
		::status.load.stop();
	}
}
/*--------------------------------------------------------------------------*/
void EV::sweep_recursive(int Nest)
{
	trace1("EV_BASE::sweep_recursive(int)", Nest);
	size_t d = _sim->_total_nodes;
	unsigned n = _sim->_total_nodes;
	--Nest;
	assert(Nest >= 0);
	assert(Nest < DCNEST);

	// double iddc[d];

	OPT::ITL itl = OPT::DCBIAS;

	first(Nest);
	do {
		trace1("EV_BASE::sweep_recursive loopstart", Nest);
		_sim->_temp_c = temp_c_in;
		if (Nest == 0) {
			_sim->_time0 = _sim->_dt0 = 0.0;
			_sim->tr_reset();
			// _sim->zero_currents();
			//
			// why not?
			_sim->_uic = true;
			_sim->_phase = p_INIT_DC;

			_sim->_loadq.clear();
			_sim->_evalq->clear();
			_sim->_evalq_uc->clear();
			trace0("hot");
			_sim->set_inc_mode_bad();

			_sim->_uic = true; // ?!
			get_op(itl);

			{ // fetch ACX
				trace0("AC::solve");
				_sim->_acx.zero();
				std::fill_n(_sim->_ac, _sim->_total_nodes+1, 0.);
				::status.load.start();
				_sim->count_iterations(iTOTAL);
				CARD_LIST::card_list.do_ac();
				CARD_LIST::card_list.ac_load();
				::status.load.stop();
			}

			// if verbose outdata(1./0.);
			SSP_CHART* chart = NULL;
			SSP_STATE& state = ev_solver(chart);
			gsl_vector_view v0__ = gsl_vector_view_array(_sim->_v0+1,n);
			SSP_SPL& spl = state.insert_spl(v0__, chart);
			USE(spl);
			USE(state);

			gsl_matrix_complex_const_view dcs = gsl_matrix_complex_const_submatrix (EV_, 0, _n_states, d, d -_n_states);
			USE(dcs);
			// state.insert_op(dcs);

			{ // some more AC stuff

				if(_dump_matrix){ untested();
					_out << "RS ( " << _Gu[0];
					for(unsigned a=1; a < d; ++a){
						_out << " " <<  _Gu[a];
					}
					_out  << ") \n";
				}

				// irgendwoher di/du holen...
				// C.fbsub( dv, Gu , dv );

			}

			if (spl.is_new()){
				EV_BASE::_out << chart->id();
				outdata(spl.id());
			}

			// here??
			//CARD_LIST::card_list.tr_regress(); incomplete(); // => do_tran_step ?
			itl = OPT::DCXFER;
		}else{
			sweep_recursive(Nest);
		}
	} while (next(Nest));
}
/*--------------------------------------------------------------------------*/
void EV_BASE::first(int Nest)
{
	trace2("EV_BASE::first", Nest, _start[Nest]);
	assert(Nest >= 0);
	assert(Nest < DCNEST);
	assert(_start);
	assert(_sweepval);
	// assert(_pushel[Nest]);

	(*_sweepval[Nest]) = _start[Nest];
	CARD_LIST::card_list.precalc_last();
	if (ELEMENT* c = dynamic_cast<ELEMENT*>(_zap[Nest])) {
		c->set_constant(false); // because of extra precalc_last
		// obsolete, once pointer hack is fixed
	}

	// here? (hack...)
	//if(_pushel[Nest])
	//   _pushel[Nest]->set_ic(_start[Nest]);
	_reverse[Nest] = false;
	if (_reverse_in[Nest]) {itested();
		while (next(Nest)) {itested();
			/* nothing */;
		}
		_reverse[Nest] = true;
		next(Nest);
	}else{
	}
	_sim->_phase = p_INIT_DC;
}
/*--------------------------------------------------------------------------*/
// EV_DISC
bool EV_BASE::next(int Nest)
{
	trace1("EV_BASE::next(int)", Nest);

	bool ok = false;
	if (_linswp[Nest]) {
		double fudge = _step[Nest] / 10.;
		if (_step[Nest] == 0.) {
			ok = false;
		}else{
			if (!_reverse[Nest]) {
				*(_sweepval[Nest]) += _step[Nest];
				fixzero(_sweepval[Nest], _step[Nest]);
				CARD_LIST::card_list.precalc_last();
				if (ELEMENT* c = dynamic_cast<ELEMENT*>(_zap[Nest])) {
					c->set_constant(false); // because of extra precalc_last
					// obsolete, once pointer hack is fixed
				}
				ok=in_order(_start[Nest]-fudge,(*_sweepval[Nest]),_stop[Nest]+fudge);
				//_pushel[Nest]->set_ic(_sweepval[Nest]);
				if (!ok  &&  _loop[Nest]) { untested();
					_reverse[Nest] = true;
				}else{
				}
			}else{ untested();
			}
			if (_reverse[Nest]) { untested();
				*(_sweepval[Nest]) -= _step[Nest];
				CARD_LIST::card_list.precalc_last();
				if (ELEMENT* c = dynamic_cast<ELEMENT*>(_zap[Nest])) {
					c->set_constant(false); // because of extra precalc_last
					// obsolete, once pointer hack is fixed
				}
				fixzero(_sweepval[Nest], _step[Nest]);
				ok=in_order(_start[Nest]-fudge,(*_sweepval[Nest]),_stop[Nest]+fudge);
				//_pushel[Nest]->set_ic(_sweepval[Nest]);
			}else{
			}
		}
	}else{ untested();
		double fudge = pow(_step[Nest], .1);
		if (_step[Nest] == 1.) {untested();
			ok = false;
		}else{ untested();
			if (!_reverse[Nest]) { untested();
				*(_sweepval[Nest]) *= _step[Nest];
				ok=in_order(_start[Nest]/fudge,*(_sweepval[Nest]),_stop[Nest]*fudge);
				//_pushel[Nest]->set_ic(_sweepval[Nest]);
				if (!ok  &&  _loop[Nest]) {untested();
					_reverse[Nest] = true;
				}else{ untested();
				}
			}else{ untested();
			}
			if (_reverse[Nest]) {untested();
				*(_sweepval[Nest]) /= _step[Nest];
				ok=in_order(_start[Nest]/fudge,*(_sweepval[Nest]),_stop[Nest]*fudge);
				//_pushel[Nest]->set_ic(_sweepval[Nest]);
			}else{ untested();
			}
		}
	}
	_sim->_phase = p_DC_SWEEP;
	return ok;
}
/*-----------------------------------------------------------*/
void EV_BASE::ac_snapshot()
{ itested(); // in sock
	trace0("EV_BASE::ac_snapshot");

	_sim->_acx.reallocate();
	_sim->_jomega = COMPLEX(0., 1.0);

	_sim->_phase = p_AC;

	// _sim->_acx.set_min_pivot(OPT::pivtol);
	//
	CARD_LIST::card_list.ac_begin();

	if(_dump_matrix){ untested();
		_out << "i ( " << _sim->_i[1];
		for(unsigned a=2; a <= _sim->_total_nodes; ++a){ untested();
			_out << " " <<  _sim->_i[a];
		}
		_out  << ") \n";
		_out << "v0 = ( " << _sim->_v0[1];
		for(unsigned a=2;a <= _sim->_total_nodes; ++a){ untested();
			_out << " " <<  _sim->_v0[a];
		}
		_out << ") \n";
	}

	// if verbose
	_sim->_uic = false;

	_sim->_acx.zero();
	std::fill_n(_sim->_ac, _sim->_total_nodes+1, 0.);

	::status.load.start();
	_sim->count_iterations(iTOTAL);
	CARD_LIST::card_list.do_ac();
	CARD_LIST::card_list.ac_load();
	::status.load.stop();
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// take some eigenvectors (alpha!=0, beta==0)
// rearrange considering input sensitivities.
// obsolete?
#if 0
void EV::align_inf_eigenvectors(gsl_matrix_complex_view& EV_inf, size_t i)
{
	if (!_zap[i]){ untested();
		return;
	}
	size_t n = _sim->_total_nodes;
	trace3("fixing\n", EV_inf.matrix, EV_inf.matrix.size1, EV_inf.matrix.size2);



	trace1("lu", _sim->_acx);

	assert(_zap[i]); // incomplete
	clear_sens();
	_zap[i]->sens_load("dummy");
	complex_array_to_gsl_vector_complex(*sens_,_sim->_sens+1,n);
	trace1("2 sensstamp", *sens_);

	::status.back.start();
	_sim->_acx.fbsub(_sim->_sens);
	::status.back.stop();

	complex_array_to_gsl_vector_complex(*sens_,_sim->_sens+1,n);
	trace2("sens", _zap[i]->long_label(), *sens_);

	gsl_vector_complex_view first_EV = firstcol(EV_inf);
	gsl_vector_complex_memcpy(&first_EV.vector, sens_);


	if(EV_inf.matrix.size2 > 1) { untested();
		EV_inf = gsl_matrix_complex_submatrix (&EV_inf.matrix, 0, 1, n, EV_inf.matrix.size2-1);
		align_inf_eigenvectors(EV_inf, i+1);
	}
}
#endif
/*--------------------------------------------------------------------------*/
// output sensitivities. compute grid.
// foreach output, use eigenvector -> output sensitivity
// FIXME: use actual probes
void EV_BASE::output_sens_grid(gsl_matrix_complex_view& relevant_EV, gsl_vector* grid_, const gsl_vector_complex* lambda)
{
	int err = 0;
	trace1("output sens\n", relevant_EV.matrix);
	trace1("output sens", *grid_);
	size_t n = _sim->_total_nodes;

	_output_iter = _output.begin();

	for(; _output_iter!=_output.end(); ++_output_iter){
		set_sens(_output_iter);
		gsl_vector_complex_view ztmp_left = gsl_vector_complex_subvector(ztmp_, 0, relevant_EV.matrix.size2 );

		complex_array_to_gsl_vector_complex(*sens_,_sim->_sens+1,n);
		trace2("sens stamp", _output_iter->label, *sens_);
		err+= gsl_blas_zgemv( CblasTrans, one, &relevant_EV.matrix, sens_, zero, &ztmp_left.vector );
		trace1("EV senstivities", *ztmp_);

		for(unsigned i=0; i<relevant_EV.matrix.size2; ++i){

			double ln = gsl_vector_get(grid_, i);
			gsl_complex lam = gsl_vector_complex_get(lambda, i);
			gsl_complex output_sens = gsl_vector_complex_get(ztmp_,i);
			trace3("grid", output_sens, lam, i);

			double n = 1;
			if(i<_n_states){
//				n = .1/sqrt(OPT::dtmin)/norm(lam);
//				trace3("grid", output_sens, lam, n);
				gsl_vector_set(grid_,i,1./0.); // reset. grid has been abused above
			}

			double on = norm(output_sens) * n;

			if (0.!=on){ itested();
				double outgrid = sqrt(OPT::reltol)/on;
//				outgrid = (OPT::reltol)/on;
				if(outgrid < fabs(ln)){ itested();
					trace3("reducing grid", i, outgrid, fabs(ln));
					gsl_vector_set(grid_,i,outgrid);
				}else{
					trace2("not reducing grid", outgrid, fabs(ln));
				}
			}
		}
	}
	trace1("output sens done", *grid_);
}
/*--------------------------------------------------------------------------*/
// quantize z0, stash index into index
// z0= \pm (2**(\pm index) - 1) * tol
// index = floor ( ln (z0 + reltol/ reltol) / ln(2) )
void EV_BASE::quantize(gsl_vector_complex_view z0, gsl_vector_complex_const_view /*lambda*/, gsl_vector_const_view tol, index_t* index)
{
	trace1("inp", z0);
	for(unsigned i=0; i<_n_states; ++i) {
		gsl_complex z = gsl_vector_complex_get(z0,i);
		double t_ = gsl_vector_get(tol,i);
		gsl_complex t = gsl_complex_rect(t_,0); // hmmm...
		int sign = 1-2* int(bool(signbit(GSL_REAL(z))));
		trace2("log discretized", z, sign);
		z = gsl_complex_mul_real(z, double(sign));

		z = gsl_complex_div(z, t);
		z = gsl_complex_add(z,gsl_complex_rect(1,0));
		z = gsl_complex_log(z);

		double ln2 = logf(2.);
		z = gsl_complex_add( gsl_complex_rect(ln2 * .5, 0. ),z) ;
		z = gsl_complex_div( z, gsl_complex_rect(ln2 ,0.));

		if (GSL_REAL(z) < 0) { untested();
			sign = 0;
		}

		index[i] = sign*int(floor(GSL_REAL(z)));
		GSL_SET_REAL(&z, sign * GSL_REAL(z));

		if (index[i]){
			z = gsl_complex_rect( sign*t_ * (pow(2,sign*index[i]-1)), 0 );
		} else {
			z = zero;
		}
		assert(fabs(norm(z))<10); // for now.
		gsl_vector_complex_set(&z0.vector,i,z);
	}
	trace2("quantized", z0, index[0]);
	// if err; throw...
}
/*--------------------------------------------------------------------------*/
string EV_BASE::status()const
{
	return string("state space: charts=")
		+ to_string(SSP_TREE::num_charts())
		+ ", states=" + to_string(SSP_TREE::num_states())
		+ ", samples=" + to_string(SSP_TREE::num_spls()) + "\n";
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
