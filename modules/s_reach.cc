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
 * reachability stuff
 */

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
#include "s_tr.h"
using namespace std;
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
namespace { //
/*--------------------------------------------------------------------------*/
class REACH : public EV_BASE, protected TRANSIENT
	{ //
	public:
		explicit REACH(): EV_BASE() {}
		~REACH() {}
		void	do_it(CS&, CARD_LIST*);
	private:
		void	setup(CS&);
		explicit REACH(const REACH&): EV_BASE(), TRANSIENT() {unreachable(); incomplete();}

		COMMON_COMPONENT* _zap_bm[DCNEST]; /* pointer to thing to sweep, dc command */
		void options_(CS&);

	protected:
		void sweep_recursive(int);

		void align_inf_eigenvectors(gsl_matrix_complex_view& EV_inf, size_t inputnumber=0);
		unsigned sens_inf_vectors(gsl_matrix_complex_view& REACH);


		void find_successors(index_t* swp, SSP_SPL& t, unsigned Nest, unsigned depth);
		void find_successors(const SSP_VECTOR& inp, SSP_SPL& t, unsigned=0, unsigned depth=0){
			index_t* swp = inp.clone_data();
			find_successors(swp, t, _n_inputs, depth);
			delete[] swp;
		}
		void find_successors(const gsl_vector* swp, SSP_SPL& t, unsigned Nest, unsigned depth){
			SSP_VECTOR v(swp);
			return find_successors(v, t, Nest, depth);
		}

		void tran_step(double dt);
	private:
		SSP_SPL* root_spl;
		unsigned _depth;
		unsigned _maxwidth;

	private:
		class SWEEP{
			public:
				SWEEP(REACH& p, SSP_TRANS& t, unsigned depth) :
				    _parent(p), _spl(t), _depth(depth) { itested();
					swp = t.input().clone_data();

				}
				~SWEEP(){}

				const SSP_CHART& chart()const{return _spl.chart();}
				unsigned n_inputs()const{return _spl.chart().n_inputs();}
				void find(unsigned Nest);
			private:
				double first(const SSP_SPL& c, index_t* swp, int Nest);
				bool next(const SSP_CHART& c, index_t* swp, int Nest) const;
				void do_it(int Nest, double dt, double bak);
				REACH &_parent;
				SSP_SPL& _spl;
				double* grid;
				index_t* swp;
				//index_t const* swp; ??
				unsigned _depth;
				index_t _last;
		};

};
/*--------------------------------------------------------------------------*/
void REACH::SWEEP::do_it(int Nest, double dt, double bak)
{
	size_t n = _sim->_total_nodes;
	gsl_vector_view v0__ = gsl_vector_view_array(_sim->_v0+1,n);
	_sim->tr_reset(); // required?
	assert(_sim->_time0 == 0.);

	_sim->restore_voltages();
	CARD_LIST::card_list.tr_restore(); // done by TRANSIENT::sweep?
	trace7("transition", _spl.id(), Nest, Nest, _parent._zap_bm[Nest]->param_name(10), _parent._zap_bm[Nest]->param_value(10), swp[Nest], dt);

	if(_parent.EV_BASE::_trace >= tVERBOSE){
		_parent.EV_BASE::outdata(-double(_spl.id()));
	}

	_parent.tran_step(dt);
	_parent.get_op_();
	double v0save[n];
	memcpy(v0save, _sim->_v0+1, n*sizeof(double));

	if(_parent.EV_BASE::_trace >= tITERATION){
		_parent.EV_BASE::outdata(1./0.);
	}

	trace1("presolver", v0__);
	SSP_CHART* chart = NULL;
	SSP_STATE& state = _parent.ev_solver(chart);
	trace1("solved", v0__);
	gsl_vector_view v0bak__ = gsl_vector_view_array(v0save,n);

	SSP_SPL& spl = state.insert_spl(v0__, chart);
	trace1("discretized", v0__);

	gsl_vector_sub(v0bak__, v0__);
	double error = gsl_blas_dnrm2(&v0bak__.vector);

	SSP_TRANS transition = _spl.insert_trans(swp, &spl);

	if (spl.is_new()) { itested();
		//				_parent.EV_BASE::_out << double(spl.id()) << double(depth) << double(chart->id());
		_parent.EV_BASE::_out << chart->id();
		_parent.EV_BASE::outdata(double(spl.id()));
	}

	_parent.EV_BASE::_out << _spl.chart().id() << "," << _spl.id() << " -> " << chart->id() << "," << spl.id();
	_parent.EV_BASE::_out << dt;
	for(unsigned i=0; i < n_inputs(); ++i){
		_parent.EV_BASE::_out << double(swp[i]) * chart->grid(i);
	}
	assert(fabs(error)<10);
	_parent.EV_BASE::_out << " err " << error << "\n";

	if (!spl.is_new()) { itested();
		trace2("notnew", spl.id(), _depth);
	}else if(!_depth){ untested();
		// border.
		// hmm push to stack?
	}else{ itested();
		trace2("isnew", spl.id(), swp);
		_sim->keep_voltages(true);
		SWEEP sweeper(_parent, transition, _depth-1);
		sweeper.find(spl.n_inputs());
		_sim->pop_voltages();
	}

	trace2("resetting iv", Nest, bak);
	_parent._zap_bm[Nest]->set_param_by_name("iv", to_string(bak));
}
/*--------------------------------------------------------------------------*/
// TODO: breadth-first search...?
void REACH::SWEEP::find(unsigned Nest)
{ itested();
	--Nest;

	trace3("find_successors", _spl.id(), Nest, swp[Nest]);
	double bak = first(_spl, swp, Nest);
	do { itested();
		if (Nest == 0) {

//			double dt = 1e-8;
//			do_it(Nest, dt, bak);
			do_it(Nest, 1e-2, bak);


		}else{ untested();
			find(Nest);
		}
	} while (next(_spl.chart(), swp, Nest));
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
void REACH::do_it(CS& Cmd, CARD_LIST* Scope)
{

	trace0("REACH::do_it(CS&, CARD_LIST*)");

	trace0("doing ddc");
	EV_BASE::_scope = Scope;
	_sim->_time0 = 0.;
	_sim->set_command_ddc();
	_sim->_phase = p_INIT_DC;
	::status.ddc.reset().start();
	_sim->_temp_c = temp_c_in;
	_do_tran_step = 0;
	EV_BASE::_dump_matrix = 0;
	EV_BASE::command_base(Cmd);
	::status.ddc.stop();
}
/*--------------------------------------------------------------------------*/
void REACH::setup(CS& Cmd)
{
	EV_BASE::_cont = false;
	EV_BASE::_trace = tNONE;
	EV_BASE::_out = IO::mstdout;
	EV_BASE::_out.reset(); //BUG// don't know why this is needed */
	bool ploton = IO::plotset  &&  EV_BASE::plotlist().size() > 0;
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
					c->tr_begin();
					trace2("REACH::setup", c->long_label(), c->value());
				} else { untested();
					throw Exception("dc/op: can't sweep " + (**ci).long_label() + '\n');
				}
				++_n_inputs;
			}else if (Cmd.is_float()) {		// sweep the generator
				_zap[_n_sweeps] = NULL;
			}else if (Cmd.is_alpha()) {
				string pname;
				unsigned here = Cmd.cursor();
				Cmd >> pname;
				try {
					_param[_n_sweeps] = &(EV_BASE::_scope->params()->find(pname));
				} catch(Exception) {
					Cmd.reset(here);
				}
				throw Exception("reach: can't sweep parameter " + pname + '\n');
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
			EV_BASE::options(Cmd,_n_sweeps);

		}
		unsigned here = Cmd.cursor();
		string output;
		bool newprobes = true;
		do{ itested();
			// FIXME: can we use the probe parser?!
			unsigned here2 = Cmd.cursor();
			if( Cmd >> "tran" ){ itested();
				Cmd.reset(here2);
				break;
			}else if( Cmd >> "v" ){
				output_t t;
				t.brh[1] = 0;
				trace1("SENS::setup, have v", Cmd.tail());

				int paren = Cmd.skip1b('(');
				Cmd >> output;
				t.label = "v("+output;

				CKT_NODE* node = dynamic_cast<CKT_NODE*>((*EV_BASE::_scope).node(output));
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
					node = dynamic_cast<CKT_NODE*>((*EV_BASE::_scope).node(output));
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
				|| Get(Cmd, "dm",	          &(EV_BASE::_dump_matrix))
				|| EV_BASE::_out.outset(Cmd);
			;
		}while (Cmd.more() && !Cmd.stuck(&here));
	}else{ untested();
	}

	TRANSIENT::_scope = EV_BASE::_scope;
	options_(Cmd);
	Cmd.check(bWARNING, "what's this?");

	IO::plotout = (ploton) ? IO::mstdout : OMSTREAM();
	initio(EV_BASE::_out);

	assert(_n_sweeps > 0);
	for (int ii = 0;  ii < _n_sweeps;  ++ii) { untested();
		trace1("REACH setup", ii);
		_start[ii].e_val(0., EV_BASE::_scope);
		fix_args(ii);

		if (_zap[ii]) {
			trace2("zap inc_probes" + _zap[ii]->long_label(), ii, _zap[ii]->value());
			_stash[ii] = _zap[ii];			// stash the std value
			_zap[ii]->inc_probes();			// we need to keep track of it

			// urghs. hack
			STORAGE* s = dynamic_cast<STORAGE*>(_zap[ii]);
			if (s) { unreachable(); incomplete();
				// trace2("REACH::setup ", _zap[ii]->long_label(), _zap[ii]->is_constant());
				// _zap[ii]->set_constant(false);		   // so it will be updated
				// _pushel[ii] = _zap[ii];	           // point to value to patch
				// _uic_caplist.push_back(s);
				// _sweepval[ii] = s->set__ic();
			}else{
				COMMON_COMPONENT* pulse = bm_dispatcher.clone("pulse");
				trace2("REACH::setup", ii, hp(pulse));
				pulse->set_param_by_name("iv", "0"); // BUG?.
				pulse->set_param_by_name("pv", "0"); // BUG?.
				_zap_bm[ii] = pulse;
				_zap[ii]->set_value(_zap[ii]->value(), pulse);
				assert(_zap[ii]->common() == pulse);
				_zap[ii]->set_constant(false);		       // so it will be updated
				_sweepval[ii] = _zap[ii]->set__value();
			}
		} else if (_param[ii]) {
			_sweepval[ii] = _param[ii]->pointer_hack();
		} else {
			trace1("REACH::setup does not exist", ii);
			//_sweepval[ii] = 0;
			_sweepval[ii] = &_sim->_genout;
			// throw(Exception("something went wrong\n"));
		}
	}
	_sim->_freq = 0;

	untested();
}
/*--------------------------------------------------------------------------*/
// unnested options
void REACH::options_(CS& Cmd)
{
	trace1("REACH::options(CS&, int)", Cmd.tail());
	_quantize_states = false;
	_depth = -1;
	_maxwidth = 10;

	_sim->_uic = false;
	unsigned here = Cmd.cursor();
	bool tr=false;
	do{ itested();
		ONE_OF
			||(Get(Cmd, "uic",            &_uic))
			|| Get(Cmd, "sr",             &_slew_rate)
			|| Get(Cmd, "dt{emp}",        &temp_c_in,   mOFFSET, OPT::temp_c)
			|| Get(Cmd, "te{mperature}",  &temp_c_in)
			|| Get(Cmd, "de{pth}",        &_depth)
			|| Get(Cmd, "max{width}",     &_maxwidth)
			|| (Get(Cmd, "qu{antize}",    &_quantize_states))
			|| (Cmd.umatch("tr{ace} {=}") &&
					(ONE_OF
					 || Set(Cmd, "n{one}",      &(EV_BASE::_trace), tNONE)
					 || Set(Cmd, "o{ff}",       &(EV_BASE::_trace), tNONE)
					 || Set(Cmd, "w{arnings}",  &(EV_BASE::_trace), tUNDER)
					 || Set(Cmd, "i{terations}",&(EV_BASE::_trace), tITERATION)
					 || Set(Cmd, "v{erbose}",   &(EV_BASE::_trace), tVERBOSE)
					 || Cmd.warn(bWARNING, 
						 "need none, off, warnings, iterations, verbose")
					)
				)
			// || ( Get(Cmd , "tran",      &tr) &&
			|| ( Cmd >> "tran" &&
					(trace1("TR", Cmd.tail()), true) &&
					(TRANSIENT::options(Cmd), tr=true)
				)

			|| (trace1("options_", Cmd.tail()), EV_BASE::_out.outset(Cmd));
	}while (Cmd.more() && !Cmd.stuck(&here) && !tr);

	if (!tr){untested();
		CS c(CS::_STRING,"");
		TRANSIENT::options(c);
	}
}
/*--------------------------------------------------------------------------*/
static REACH p2;
EV_BASE* p3 = dynamic_cast<EV_BASE*>(&p2);
static DISPATCHER<CMD>::INSTALL d2(&command_dispatcher, "reach", p3);
/*--------------------------------------------------------------------------*/
void REACH::sweep_recursive(int Nest)
{ untested();
	size_t n = _sim->_total_nodes;
	trace2("REACH::sweep_recursive(int)", Nest, _sim->_uic);
	//  unsigned n = _sim->_total_nodes;

	--Nest;
	assert(Nest >= 0);
	assert(Nest < DCNEST);

	// double iddc[d];

	OPT::ITL itl = OPT::DCBIAS;

	_sim->_uic = _uic;

	CARD_LIST::card_list.precalc_last();
	double inputvals[_n_inputs+1];

	for(unsigned i=0; i < _n_inputs; ++i) { untested();
		inputvals[i+1] = _zap[i]->value();
		trace2("setup pulse input", _zap[i]->long_label(), _zap[i]->value());
		_zap_bm[i]->set_param_by_name("iv", to_string(inputvals[i+1]));
		_zap_bm[i]->set_param_by_name("pv", to_string(inputvals[i+1]));
		_zap[i]->precalc_first();
		_zap[i]->precalc_last();
		_zap[i]->tr_begin();

		if (ELEMENT* c = dynamic_cast<ELEMENT*>(_zap[i])) { untested();
			c->set_constant(false);
		}
	}

	get_op(itl); // FIXME: only one root op.
	             // need to untangle input range from dc-sweep init points
	_sim->_uic = false;

	if(EV_BASE::_trace >= tVERBOSE){
		EV_BASE::outdata(1./0.);
	}

	SSP_CHART* chart = NULL;
	SSP_STATE& state = ev_solver(chart);
	assert(chart);
	assert(chart->id()<10);
	gsl_vector_view v0__ = gsl_vector_view_array(_sim->_v0+1,n);

	double probes[EV_BASE::printlist().size()];

	size_t k = 0;
	for (PROBELIST::const_iterator
			p=EV_BASE::printlist().begin();  p!=EV_BASE::printlist().end();  ++p) {
		double v= (*p)->value();
		probes[k++] = v;
	}
	USE(probes);

	trace1("to", v0__);
	untested();
	SSP_CHART* ch = chart; // new SSP_CHART(_n_inputs);
	SSP_SPL& spl = state.insert_spl(v0__, ch);
	trace1("ins", v0__);

	inputvals[_n_inputs] = 1./0.;
	for(unsigned i=0; i < _n_inputs; ++i){
		inputvals[i] = *(_sweepval[i]);
		ch->_inputgrid[i] = _step[i];
	}
	gsl_vector_view inp__ = gsl_vector_view_array(inputvals,_n_inputs+1);
	root_spl = &spl;
	SSP_TRANS transition = root_spl->insert_trans(&inp__.vector, &spl);

	SSP_VECTOR inp(inputvals,_n_inputs+1);

	if (spl.is_new()){
		EV_BASE::_out << chart->id();
		EV_BASE::outdata(spl.id());
	}

	if (!_depth){ untested();
	}else	if (spl.is_new()){ untested();
		trace1("isnew", v0__.vector);
		_sim->keep_voltages(true);
		SSP_VECTOR swp(transition.input());
		SWEEP sweeper(*this, transition, _depth-1);
		sweeper.find(spl.n_inputs());
		_sim->pop_voltages();
	}else{ unreachable(); // for now.
		trace1("notnew", v0__.vector);
	}
}
/*--------------------------------------------------------------------------*/
//double REACH::first(SSP_VECTOR& t, int Nest) const
double REACH::SWEEP::first(const SSP_SPL& s, index_t* swp, int Nest)
{ itested();
	assert(Nest >= 0);
	assert(Nest < DCNEST);
	const SSP_CHART& ch = s.chart();
	assert (unsigned(Nest) < ch.n_inputs());

	double start = -.2; // FIXME!
	double stop = -start;

	int width = _parent._maxwidth; // FIXME: must depend on SR

	double from = ch.grid(Nest) * double(swp[Nest]);

	index_t min_index = index_t(floor( (start+ch.grid(Nest)*.5) /(ch.grid(Nest))));
	index_t max_index = index_t(floor( (stop+ch.grid(Nest)*.5) /(ch.grid(Nest))));
	index_t here = swp[Nest];

	assert(min_index<=max_index);

	_last = min(max_index, here+width);
	swp[Nest] = max(min_index, here-width);

	double to = ch.grid(Nest) * double(swp[Nest]);
	trace5("REACH::first, pulse", ch.id(), Nest, ch.grid(Nest), from, to);
	trace3("REACH::first, pulse", ch.id(), swp[Nest], _last);
	trace2("REACH::first, pulse", min_index, max_index);

	_parent._zap_bm[Nest]->set_param_by_name("iv", to_string(from));
	_parent._zap_bm[Nest]->set_param_by_name("pv", to_string(to));

	*(_parent._sweepval[Nest]) = to; // needed by inputsens scaling

	return from;
}
/*--------------------------------------------------------------------------*/
// FIXME: individual precisions
// FIXME: individual inputno?
// FIXME: break upon leaving chart.
bool REACH::SWEEP::next(const SSP_CHART& c, index_t* swp, int Nest) const
{
	++swp[Nest];
	double p(c.grid(Nest));
	double q = double(swp[Nest]);
	assert(c.grid(Nest));
	double to = q*p;
	_parent._zap_bm[Nest]->set_param_by_name("pv", to_string(to));
	*(_parent._sweepval[Nest]) = to; // BUG? used by inputsens scaling

	trace5("EV_BASE::next", Nest, swp[Nest], p, to, _last);
	return swp[Nest] <= _last;
}
/*--------------------------------------------------------------------------*/
// take some eigenvectors (alpha!=0, beta==0)
// rearrange considering input sensitivities.
void REACH::align_inf_eigenvectors(gsl_matrix_complex_view& EV_inf, size_t i)
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

#if 0
	relevant_EV = gsl_matrix_complex_submatrix (EV_, 0, 0, n, fin_dim);
	gsl_vector_complex_view relevant_v0prime_ = gsl_vector_complex_subvector(v0prime_,0,fin_dim);
	trace1("preimg cut ", relevant_v0prime_.vector);
	gsl_matrix_complex_view relevant_EV = gsl_matrix_complex_submatrix (EV_, 0, 0, n, fin_dim);
	trace1("relevant EV\n", (relevant_EV.matrix));
	err+= gsl_blas_zgemv( CblasNoTrans, one, &relevant_EV.matrix, &relevant_v0prime_.vector, zero, ztmp_ );
	trace1("proj\n", *ztmp_);

	//extra vector: normalize(v0 - ztmp)
	err+= gsl_vector_complex_sub(v0_, ztmp_);
	double len = gsl_blas_dznrm2(v0_);
	trace1("res\n", *v0_);
	gsl_vector_complex_view v = gsl_matrix_complex_column (EV_, fin_dim);
	gsl_vector_complex_memcpy (&v.vector, v0_);

	//    gsl_vector_complex_set(relevant_v0prime.vector,0,len);
	gsl_complex c; GSL_REAL(c)=1./len; GSL_IMAG(c)=0;
	gsl_vector_complex_scale (&v.vector, c);
	fin_dim+= 1;
#endif

}
/*--------------------------------------------------------------------------*/
// FIXME: this is ad-hoc.
void REACH::tran_step(double dt)
{ itested();

	for(unsigned i=0; i < _n_inputs; ++i) {
		assert(_zap[i]);
		assert(_zap_bm[i]);
		_zap_bm[i]->set_param_by_name("rise", to_string(dt));
		_zap[i]->precalc_last();
		_zap[i]->tr_begin(); // hmmm?
	}

	size_t n = _sim->_total_nodes;
	_sim->set_command_tran();
	_tstart = 0.;
	TRANSIENT::_time1 = 0.;
	TRANSIENT::_cont_dc = true;
	TRANSIENT::_trace = tALLTIME;
	_sim->_time0 = 0.;
	_sim->_dtmin = dt/10.;
	_dtmax = dt/2.;
	_tstop = dt;
	_tstep = dt;
	_sim->_freezetime = true;
	TRANSIENT::sweep();
	_sim->_freezetime = false;
	gsl_vector_view v0__ = gsl_vector_view_array(_sim->_v0+1,n);
	trace1("swept", v0__);
}
/*--------------------------------------------------------------------------*/
}// namespace
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
