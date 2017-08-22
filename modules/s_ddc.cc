/** Copyright (C) 2010 Peter
 * Author: Peter, Felix
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
 * ddc analysis top
 */
#include "config.h"

#include "u_status.h"
#include "u_prblst.h"
#include "u_cardst.h"
#include "e_elemnt.h"
#include "e_storag.h"
#include "s__.h"
#include "s_ddc.h"
#include "io_matrix.h"
#include "m_matrix_extra.h"
#ifdef HAVE_CLAPACK_H
#include "u_atlas.h"
#endif
using namespace std;
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace { //
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
class DDC : public DDC_BASE { //
public:
  explicit DDC(): DDC_BASE() {}
  ~DDC() {}
  void	do_it(CS&, CARD_LIST*);
private:
  void	setup(CS&);
  explicit DDC(const DDC&): DDC_BASE() {unreachable(); incomplete();}
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
void DDC::do_it(CS& Cmd, CARD_LIST* Scope)
{

  trace0("DDC::do_it(CS&, CARD_LIST*)");

  trace0("doing ddc");
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
void DDC::setup(CS& Cmd)
{
  _cont = false;
  _trace = tNONE;
  _out = IO::mstdout;
  _out.reset(); //BUG// don't know why this is needed */
  bool ploton = IO::plotset  &&  plotlist().size() > 0;

  if (Cmd.more()) {
    for (_n_sweeps = 0; Cmd.more() && _n_sweeps < DCNEST; ++_n_sweeps) {
      CARD_LIST::fat_iterator ci = findbranch(Cmd, &CARD_LIST::card_list);
      if (!ci.is_end()) {
	if (ELEMENT* c = dynamic_cast<ELEMENT*>(*ci)) {
	  _zap[_n_sweeps] = c;
          trace1("DDC::setup", _zap[_n_sweeps]->long_label());
	} else { untested();
	  throw Exception("dc/op: can't sweep " + (**ci).long_label() + '\n');
	}
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

      // urghs. hack
      STORAGE* s = dynamic_cast<STORAGE*>(_zap[ii]);
      if (0) {
        trace2("DDC::setup ", _zap[ii]->long_label(), _zap[ii]->is_constant());
        _zap[ii]->set_constant(false);		   // so it will be updated
        _pushel[ii] = _zap[ii];	           // point to value to patch
        _uic_caplist.push_back(s);
        _sweepval[ii] = s->set__ic();
      }else{
        trace2("DDC::setup, no STORAGE", ii, _zap[ii]->long_label());
        _zap[ii]->set_value(_zap[ii]->value(),0);  // zap out extensions
        _zap[ii]->set_constant(false);		   // so it will be updated
        _pushel[ii] = _zap[ii];	                   // element to patch
        _sweepval[ii] = _zap[ii]->set__value();
      }
      //_sweepval[ii] = 0;	        
    } else if (_param[ii]) {
      _sweepval[ii] = _param[ii]->pointer_hack();
    } else {
      trace1("DDC::setup does not exist", ii);
      //_sweepval[ii] = 0;
      _sweepval[ii] = &_sim->_genout;
      _pushel[ii] = NULL; //&_sim->set_gen;			// point to value to patch
      // throw(Exception("something went wrong\n"));
    }
  }
  _sim->_freq = 0;

}
/*--------------------------------------------------------------------------*/
static DDC p2;
static DISPATCHER<CMD>::INSTALL d2(&command_dispatcher, "ddc", &p2);
}// namespace
/*--------------------------------------------------------------------------*/
double	DDC_BASE::temp_c_in = 0.;
/*--------------------------------------------------------------------------*/
DDC_BASE::DDC_BASE()
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
void DDC_BASE::finish(void)
{
  trace0("DDC_BASE::finish");

//  _sim->_phase = p_;

  for (int ii = 0;  ii < _n_sweeps;  ++ii) {
    if (_zap[ii]) { // component
      _stash[ii].restore();
      _zap[ii]->dec_probes();
      trace2("DDC_BASE::finish dec_probes done", _zap[ii]->long_label(), _zap[ii]->probes());
      _zap[ii]->precalc_first();
      _zap[ii]->precalc_last();
    }else{
    }
  }

  _uic_caplist.clear();
}
/*--------------------------------------------------------------------------*/
void DDC_BASE::do_tran_step()
{
  trace2("DDC_BASE::do_tran_step", _tran_step, _sim->iteration_tag());
  _sim->_phase = p_TRAN;
  SIM_MODE old_mode = _sim->_mode;
  //_sim->_bypass_ok = false;
  _sim->_mode = s_TRAN;
      _sim->_evalq->clear();
      _sim->_evalq_uc->clear();

#if 0
  _sim->restore_voltages();
//   assert(!_sim->_time0);
  // _sim->_phase = p_INIT_DC;
  _sim->_phase = p_RESTORE; // short cut differentiate
  CARD_LIST::card_list.tr_begin(); // breaks stuff. hmmm
    _sim->_cont = false;
#endif

  _sim->_phase = p_TRAN; //?

  _sim->_time0 = _sim->_dt0 = _tran_step;
  //_sim->_genout = gen();

  assert(!_sim->uic_now()); // bug?
  assert(_sim->analysis_is_tran()); // bug?
  _sim->count_iterations(iTOTAL);
  _sim->_loadq.clear();
  trace1("DDC_BASE::do_tran_step done_tr", _sim->iteration_tag());



  int tr_converged = solve(OPT::TRHIGH, _trace);
  _sim->count_iterations(iTOTAL);
  trace1("DDC_BASE::do_tran_step done_solve", _sim->iteration_tag());

  if (!tr_converged) { untested();
    error(bWARNING, "did not converge\n");
  }else{
  }
  ::status.accept.start();
  trace0("DDC_BASE::sweep_recursive solved a transient step");

  _sim->set_limit();
  CARD_LIST::card_list.tr_accept();
  trace0("DDC_BASE::sweep_recursive itr_accepted");

  ::status.accept.stop();

  _sim->_time0 = 0;
  _sim->_mode = old_mode;
  _sim->_phase = p_RESTORE;
  //_sim->restore_voltages(); ????
  _sim->keep_voltages(); //  vdc  = v0
  _sim->put_v1_to_v0(); // v0 = vt1
  trace2("DDC_BASE::do_tran_step regress", _tran_step, _sim->iteration_tag());
  CARD_LIST::card_list.tr_regress();
  trace2("DDC_BASE::do_tran_step done", _tran_step, _sim->iteration_tag());
}
/*=========================*/
void DDC_BASE::fix_args(int Nest)
{
  trace1("DDC_BASE::fix_args(int Nest)", Nest);

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
void DDC_BASE::options(CS& Cmd, int Nest)
{
  
  trace1("DDC_BASE::options(CS&, int)", Nest);

  _loop[Nest] = _reverse_in[Nest] = false;
  _sim->_uic = false;
  _old_solver = false;
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
      || (Get(Cmd, "step",	  &_step_in[Nest]) && (_stepmode[Nest] = LIN_STEP))
      || (Get(Cmd, "d{ecade}",	  &_step_in[Nest]) && (_stepmode[Nest] = DECADE))
      || (Get(Cmd, "ti{mes}",	  &_step_in[Nest]) && (_stepmode[Nest] = TIMES))
      || (Get(Cmd, "lin",	  &_step_in[Nest]) && (_stepmode[Nest] = LIN_PTS))
      || (Get(Cmd, "o{ctave}",	  &_step_in[Nest]) && (_stepmode[Nest] = OCTAVE))
      || Get(Cmd, "c{ontinue}",   &_cont)
      || Get(Cmd, "tr{s}",        &_do_tran_step)
      || Get(Cmd, "dm",           &_dump_matrix)
      || Get(Cmd, "old",          &_old_solver)
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
void DDC_BASE::sweep()
{

  trace0("DDC_BASE::sweep()");

  head(_start[0], _stop[0], " ");
//  _sim->_bypass_ok = false;
  _sim->set_inc_mode_bad();
  if (_cont) {untested();
    _sim->restore_voltages();
    CARD_LIST::card_list.tr_restore();
  }else{
    _sim->clear_limit();
    CARD_LIST::card_list.tr_begin();
  }

  unsigned d = _sim->_total_nodes; // 3
  U = new double[d*d];
  CU = new double[d*d];
  CUTCU = new double[d*d];

  unsigned d2 = 2 * d;
  trace1("", d2);
  // GLS: A * y = RS2
  A = new double[d2*d2];
  y = new double[d2];
  
  sweep_recursive(_n_sweeps);
}
/*--------------------------------------------------------------------------*/
void DDC_BASE::set_uic_caps_constant(bool x){
  for( vector<STORAGE*>::iterator i=_uic_caplist.begin(); i!=_uic_caplist.end(); ++i)
  {
    (*i)->set_constant(x);
  }
}
/*--------------------------------------------------------------------------*/
//here?? -> DDC
void DDC_BASE::sweep_recursive(int Nest)
{

  trace1("DDC_BASE::sweep_recursive(int)", Nest);

  unsigned d = _sim->_total_nodes; // 3
  //trace1("DDC_BASE::sweep_recursive",d);

  //trace1("DDC_BASE::sweep_recursive", Nest);
  --Nest;
  assert(Nest >= 0);
  assert(Nest < DCNEST);

  // double iddc[d];

  OPT::ITL itl = OPT::DCBIAS;
  
  first(Nest);
  do {
    trace1("DDC_BASE::sweep_recursive loopstart", Nest);
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
      int converged = solve_with_homotopy(itl,_trace);
      trace0("hot done");


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

      for(unsigned a=0; a <= _sim->_total_nodes; ++a){
        // iddc[a]=_sim->_i[a];
      }

//      _sim->_uic = false; // see DEV_CAP
      _sim->_uic = true; // see DEV_CAP
      CARD_LIST::card_list.tr_accept(); // make dv work (?)
      _sim->_phase = p_INIT_DC;
      _sim->_uic = false; // see DEV_CAP
      _sim->_evalq->clear();
      _sim->_evalq_uc->clear();
      _sim->_loadq.clear();
      _sim->count_iterations(iTOTAL);
      ::status.accept.stop();

      finish_building_evalq();
      evaluate_models();

      _sim->keep_voltages(); // vdc = v0

      _sim->_acx.reallocate();
      _sim->_jomega = COMPLEX(0., 1.0);
      // _sim->_mode=s_AC;
      // _sim->_acx.set_min_pivot(OPT::pivtol);
      //
//      set_uic_caps_constant(false);
      _sim->_phase = p_AC;

      {// AC::sweep
        CARD_LIST::card_list.ac_begin();
        //...
      }

      trace0("solved with homotopy");
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

      // if verbose
      _sim->_uic = false;
      _sim->set_inc_mode_yes();
      while (!_sim->_loadq.empty()) {
	trace1("loading from q", _sim->_loadq.back()->long_label());
	_sim->_loadq.back()->tr_load();
	_sim->_loadq.pop_back();
      }

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

      double Gul[_sim->_total_nodes+1];
      _Gu = Gul+1;
      double dv[_sim->_total_nodes+1];
      _dv = dv+1;

      if (_old_solver){
        old_solver();
      }else{
        block_solver();
      }

      if (_do_tran_step) {
	do_tran_step();
      } else {
      }

      { // some more AC stuff

	if(_dump_matrix){
	  _out << "RS ( " << _Gu[0];
	  for(unsigned a=1; a < d; ++a){
	    _out << " " <<  _Gu[a];
	  }
	  _out  << ") \n";
	}

	// irgendwoher di/du holen...
	// C.fbsub( dv, Gu , dv );

	for(unsigned a=0; a < d; ++a){
	  // stash hack, abuse _vt1, after transtep
	  _sim->_vt1[a+1] = _dv[a];
	}
      }

      //_sim->set_command(s_DC);
      outdata(*_sweepval[Nest], ofPRINT);
      //_sim->set_command(s_DDC);

      // here??
      //CARD_LIST::card_list.tr_regress(); incomplete(); // => do_tran_step ?
      itl = OPT::DCXFER;
    }else{
      sweep_recursive(Nest);
    }
  } while (next(Nest));
}
/*--------------------------------------------------------------------------*/
void DDC_BASE::first(int Nest)
{
  trace2("DDC_BASE::first", Nest, _start[Nest]);
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
bool DDC_BASE::next(int Nest)
{
  trace3("DDC_BASE::next(int)", Nest, _step[Nest], *_sweepval[Nest]);
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
	  c->set_constant(false);
	}
	trace1("DDC_BASE::next(int)", *_sweepval[Nest]);
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
void DDC_BASE::ac_snapshot()
{ itested(); // in sock
  trace0("DDC_BASE::ac_snapshot");

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
void DDC_BASE::old_solver()
{
#ifndef HAVE_CLAPACK_H
  error(bWARNING, "no atlas\n");
#else
  trace0("DDC_BASE::old_solver()");

  BSMATRIX<double> G = _sim->_acx.real();
  BSMATRIX<double> C = _sim->_acx.imag();

  unsigned d = _sim->_total_nodes;
  double RS[d];

  if(_dump_matrix){ untested();
    _out << "G\n" << G << "\n";
    _out << "C\n" << C << "\n";
  }

  double col[d+1];
  double CU[d*d];

  double* Gu = _Gu;
  // Gu = G * v0
  G.rmul(_Gu-1, _sim->_v0);
  G.lu_decomp();

  // U = G^{-1} C (column-major)
  for( unsigned i=0; i<d; ++i){
    C.col(col,1+i);
    double buf[d+1];

    // compute ith column of U
    G.fbsub(buf, col);
    for (unsigned k=0; k<d; ++k ){
      U[k+i*d] = buf[k+1];
    }
    if(_dump_matrix){ untested();
      _out << "U row " << i << ":  " << U[0 + i*d];
      for(unsigned a=1; a < d; ++a){ untested();
        _out << " " <<  U[a + i*d ];
      }
      _out  << ") \n";
    }
  }
  // CU=C*U (col maj)
  for(unsigned i=0; i < d; ++i){
    for(unsigned j=0; j < d; ++j){
      CU[i+j*d] = 0;
      for(unsigned k=0;k < d; ++k){
        CU[i+j*d] += ((const BSMATRIX<double>)C).s(i+1,k+1) * U[k+j*d];
      }
    }
  }

  /// old solver

  // U = G^{-1} C
  // want to solve C x = i - Gu, and G x = C y <=> x = U y
  // solving C U y = i - Gu

  double alpha=1;
  double wurk;
  int rank;
  int lwork=-1;
  int info=0;
  double rcond=0;
  int D = int(d);
  int one = 1;

  // does not work (not necessarily full rank)
  // clapack_dgesv(CblasColMajor, d, 1,
  //             CU, d, ipiv, X, d);

  // void dgels_(const char *trans, const int *M, const int *N, const int *nrhs,
  // double *A, const int *lda, double *b, const int *ldb, double *work, const
  // int * lwork, int *info);

  /*
     dgels_("No transpose", &D,&D,&one, CU, &D, RS, &D, &wurk, &lwork, &info );
     assert(!info);
     lwork = (int)wurk;
     double work[lwork];
     dgels_("No transpose", &D,&D,&one, CU, &D, RS, &D, work, &lwork, &info );
     */
  for(unsigned a = 0; a < d; ++a){
    RS[a] = - Gu[a] + _sim->_i[a+1] ;
    // X[a] = - Gu[a] + _sim->_i[a+1] ;
  }

  double S[d]; //singular values
  dgelss_(&D,&D,&one,CU, &D, RS, &D, S, &rcond, &rank, &wurk, &lwork, &info );
  assert(!info);
  lwork = (int)wurk;
  double work[lwork];
  dgelss_(&D,&D,&one,CU, &D, RS, &D, S, &rcond, &rank, work, &lwork, &info );
  assert(!info);

  assert(D == int(d));
  // result in RS

  if(_dump_matrix){ untested();
    _out << "after dgels " << lwork  << " RS( " << RS[0];
    for(unsigned a=1; a < d; ++a){ untested();
      _out << " " <<  RS[a];
    }
    _out  << ") \n";
  }

#ifndef HAVE_LIBGSLCBLAS
  incomplete();
  USE(alpha);
#else
  // dv = U * app
  cblas_dgemv(CblasColMajor,
      CblasNoTrans, D, D,
      1.0/alpha,
      U, D,
      RS, 1, 0,
      _dv,1);
#endif


  if(_dump_matrix){ untested();
    _out << "Gu ( " << Gu[0];
    for(unsigned a=1; a < d; ++a){ untested();
      _out << " " <<  Gu[a];
    }
    _out  << ") \n";
  }

  for(unsigned a=0; a < d; ++a){
    Gu[a] = - Gu[a] +  _sim->_i[a+1] ;
  }


  C.dezero( OPT::cmin ); 
  C.lu_decomp();

  if(_dump_matrix){ untested();
    _out << "G\n" << G << "\n";
    _out << "C\n" << C << "\n";
    //_out << "A\n" << A << "\n";
  }

//  _sim->_bypass_ok = false;

#endif
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
void DDC_BASE::block_solver()
{
#ifndef HAVE_CLAPACK_H
  error(bWARNING, "no atlas\n");
#else
  trace0("DDC_BASE::block_solver()");

  BSMATRIX<double> G = _sim->_acx.real();
  BSMATRIX<double> C = _sim->_acx.imag();

  if(_dump_matrix){
     _out << "==========G===========\n" << G << "\n";
     _out << "==========C===========\n" << C << "\n";
     _out << "=========RHS==========\n";
     string comma("");
     for (unsigned i=0; i<_sim->_total_nodes; ++i) {
       _out << comma << _sim->_i[i+1];
       comma = ", ";
     }
     _out << "\n======================\n";
  }

  unsigned d = _sim->_total_nodes;
  unsigned d2 = 2*d;
  double* Gu = _Gu;
  // Gu = G * v0
  G.rmul(_Gu-1, _sim->_v0);

  double RS2[d*2];
  double X[d];
  // fetch rhs
  for(unsigned a = 0; a < d; ++a){
    // RS[a] = - Gu[a] + _sim->_i[a+1] ;
    // X[a] = - Gu[a] + _sim->_i[a+1] ;
    RS2[a] = - Gu[a] + _sim->_i[a+1];
  }
  for(unsigned a = d; a < d*2; ++a){
    RS2[a] = 0;
  }

  if(_dump_matrix){
    _out << "RS2 ( " << RS2[0];
    for(unsigned a=1; a<d*2; a++){
      _out << " " << RS2[a];
    }
    _out << ") \n";
  }

  //_out << "Allokation von A\n";
  //_out << "d: " << d << "\n";

  // Das Array A allozieren
  for(unsigned i=0; i<d; i++){
    for(unsigned j=0; j<d; j++){
      A[i+j*d*2] = ((const BSMATRIX<double>)C).s(i+1,j+1);
    }
    for(unsigned j=d; j<d2; j++){
      A[i+j*d2] = 0;
    }
  }

  for(unsigned i=d; i<d2; i++){
    for(unsigned j=0; j<d; j++){
      A[i+j*d2] = ((const BSMATRIX<double>)G).s(i+1-d,j+1);
    }
    for(unsigned j=d; j<d2; j++){
      A[i+j*d2] = ((const BSMATRIX<double>)C).s(i+1-d,j+1-d);
    }
  }
  // Das Array A ausgeben
  for(unsigned i=0; i<d2; i++){
    if(_dump_matrix){
      _out << " A row "<< i << " ( " <<  A[i];
      for(unsigned a=1; a<d2; a++){
        _out << " " << A[i+a*d2];
      }
      _out << " ) \n";
    } 
  }


  // X = CU.trans * App
  //cblas_dgemv(CblasColMajor, 
  //    CblasTrans, d, d,
  //    1.0, 
  //    CU, d,
  //    RS, 1, 0,
  //    X,1);

  if(_dump_matrix){
    _out << "X ( " << X[0];
    for(unsigned a=1; a < d; ++a){
      _out << " " <<  X[a];
    }
    _out  << ") \n";
  }

  int D2 = int(d2);
  int ONE = 1;
  double RCOND = 0;
  int RANK;
  double WURK;
  int LWORK=-1;
  int INFO = 0;

  // der ansatz mit der 2x2 blockmatrix
  double S2[d2]; //singular values
  dgelss_( &D2, &D2, &ONE, A, &D2, RS2, &D2, S2, &RCOND, &RANK, &WURK, &LWORK, &INFO );
  assert(!INFO);
  LWORK = (int)WURK;
  double WORK[LWORK];
  dgelss_( &D2, &D2, &ONE, A, &D2, RS2, &D2, S2, &RCOND, &RANK, WORK, &LWORK, &INFO );
  //
  if(_dump_matrix){
    //        _out << "after dgels " << lwork  << " RS( " << RS[0];
    //        for(unsigned a=1; a < d; ++a){ untested();
    //          _out << " " <<  RS[a];
    //        }
    //        _out  << ") \n";

    _out << "after dgels "  << " RS2( " << RS2[0];
    for(unsigned a=1; a<d2; ++a){
      _out << " " << RS2[a];
    }
    _out << ") \n";
  }
  _sim->_bypass_ok = false;
  for(unsigned a=0; a < d; ++a){
    // stash hack
    _dv[a] = RS2[a];
  }

#endif
}
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet
