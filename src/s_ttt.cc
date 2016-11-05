/*                              -*- C++ -*-
 * Copyright (C) 2010-2012
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
 * Foundation, Inc., 53 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 *------------------------------------------------------------------
 * performs tt analysis
 *
 *	
 *	do this in a loop:
 *
 * - set T, dT
 *	- tt_begin/tt_continue
 *	-  tt_advance / regress
 *	-  do_tt
 *	-  TR::sweep (calls tr_accept)
 *	-  tr_stress_last
 *	-  tt_review
 *	-  tt_accept (if queued?)
 *
 */
#include "declare.h"	/* plclose, plclear, fft */
#include "u_prblst.h"
#include "s_tr.h"
#include "u_nodemap.h"
#include "e_node.h"
#include "u_sim_data.h"
#include "s__.h"
#include "u_status.h"
#include "u_function.h" // after?
#include "m_wave.h"
#include "u_prblst.h"
#include "ap.h"
#include "s_tt.h"


// #include "globals.h" ??
#include "e_adp.h"
#include "e_adplist.h"
using namespace std;
/*--------------------------------------------------------------------------*/
#define sanitycheck()  assert ( is_almost(  _sim->_Time0 , _Time1 + _sim->_dT0 ));\
                       assert ( _sim->_Time0 > _Time1 || _sim->tt_iteration_number()==0 );\
                       assert ( _sim->_dT0 == 0 || _sim->tt_iteration_number()!=0 );
/*--------------------------------------------------------------------------*/
namespace TT { //
/*--------------------------------------------------------------------------*/
OMSTREAM TTT::mstdout(stdout);
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
#define check_consistency_tt() {						\
	trace4("", __LINE__, newTime, almost_fixed_Time, fixed_Time);	\
	assert(almost_fixed_Time <= fixed_Time);				\
	assert(newTime <= fixed_Time);					\
	/*assert(newTime == fixed_Time || newTime <= fixed_Time -_dTmin);*/	\
	assert(newTime <= almost_fixed_Time);				\
	/*assert(newTime == almost_fixed_Time || newTime <= almost_fixed_Time - _dTmin);*/ \
	assert(newTime > _Time1);						\
	assert(newTime > refTime);						\
	assert(new_dT > 0.);						\
	assert(new_dT >= _dTmin);						\
	assert(newTime <= _Time_by_user_request);				\
	/*assert(newTime == _Time_by_user_request*/				\
	/*	   || newTime < _Time_by_user_request - _dTmin);	*/	\
}
/*--------------------------------------------------------------------------*/
int TTT::step_cause()const
{
	return ::status.control/100;
}
/*--------------------------------------------------------------------------*/
double behaviour_time()
{ untested();
	return 80000;
}
/*--------------------------------------------------------------------------*/
void TTT::set_step_cause(STEP_CAUSE C)
{
	switch (C) { //
		case scITER_A:untested();
		case scADT:untested();
		case scITER_R:
		case scINITIAL:
		case scSKIP:
		case scTE:
		case scAMBEVENT:
		case scEVENTQ:
		case scUSER:
		case scGROW:
		case scLAST:
			::status.control = 100*C;
			break;
		case scNO_ADVANCE: untested();
		case scZERO: untested();
		case scSMALL:
		case scREJECT:
			::status.control += 100*C;
			break;
	}
}
/*--------------------------------------------------------------------------*/
double TTT::get_new_dT()
{ untested();
	double factor=  (double) steps_total();
	double buf= 100 / CKT_BASE::tt_behaviour_rel * factor;
	return buf;
}
/*--------------------------------------------------------------------------*/
void TTT::rescale_behaviour()
{ untested();
	incomplete();
}
/*--------------------------------------------------------------------------*/
void TTT::first()
{
	trace2("TTT::first()", _Time1, _Tstart);
	_sim->force_tt_order(0);
//	_sim->zero_voltages();

	// ADP_NODE_LIST::adp_node_list.do_forall( &ADP_NODE::tt_clear );

	 // _Time1 = 0; BUG! a first after powerdown doesnt start at 0
	_sim->_dT0 = 0;
	_sim->_dT1 = 0;
	_sim->_dT2 = 0;
	_dT_by_adp = _tstop;

	if (_sim->_loadq.size()) {
		trace0("TTT::first loadq nonempty -- clearing");
		_sim->_loadq.clear() ; // why?
	} else {
	}

	// tell gnuplot that we start from the beginning.
//	TRANSIENT::_out << "* newline hack\n\n";

	assert(_sim->get_tt_order() == 0 );

 	_Time_by_user_request = _sim->_Time0 = _Tstart;

	if (_Tstop>0) assert(_Tstep>0);

	_sim->_tt_accepted = 0;
	_sim->_tt_rejects = 0;
	_sim->_tt_rejects_total = 0;
	// _stepno_tt = 0;

	assert(_sim->_Time0 >= 0 );

	advance_Time();

	_sim->set_command_tt();

	trace1("TTT::first", _sim->_Time0);

	if (TRANSIENT::_trace >= tDEBUG) { untested();
		TRANSIENT::outdata(_sim->_Time0, ofPRINT | ofSTORE);
	}else{
// 		outdata(_sim->_Time0);
	}
	::status.tran.reset().start();

	//CARD_LIST::card_list.tr_begin(); ///????
	trace0("TTT::first sweep");
	CARD_LIST::card_list.do_tt(); // before sweep_tr
	// BUG: outdata here maybe
	//
	print_head_tr(); // always (no review)
	TTT::sweep_tr();
	print_foot_tr();
	trace1("TTT::first sweep done", _sim->_Time0);
	// assert (_sim->_loadq.empty());

	::status.hidden_steps=1;
	++steps_total_tt;
	::status.review.stop();

	_sim->set_command_tt();
	set_step_cause(scUSER);

	assert(_sim->_Time0 >= 0 );
	assert( _sim->_mode  == s_TTT );
	_accepted_tt = true;

//	CARD_LIST::card_list.tt_review();
	TTT::review();
	assert(_time_by_ambiguous_event > _sim->_Time0);
	TTT::accept();
	outdata_tt(_sim->_Time0); // first.

	_Time_by_user_request = _sim->_last_Time = _sim->_Time0 + _sim->last_time();
	_sim->_dT0 = _sim->last_time();
	_Time1 = _Tstart;
} //first
/*--------------------------------------------------------------------------*/
void TTT::tt_begin()
{
	_sim->_dT0 = 0;
	_sim->_dT1 = 0;
	_sim->_dT2 = 0;
	_sim->tr_reset(); // assert  _sim->_last_time == 0;?
	trace0("TTT::tt_begin");
	_sim->_stepno = 0;
	_sim->_tt_uic = false;
	CARD_LIST::card_list.tt_begin();
	// advance_time() // must be done, as tt1 is used in do_tt...
	CARD_LIST::card_list.do_tt();
}
/*--------------------------------------------------------------*/
void TTT::tt_cont()
{
	// continue from externally set adp_node values...
	_sim->_dT0 = 0;
	_sim->_dT1 = 0;
	_sim->_dT2 = 0;
	_sim->_time0 = 0;
	trace0("TTT::cont");
	_sim->_stepno = 0;
	_sim->_tt_uic = true;
	CARD_LIST::card_list.tt_begin();
	CARD_LIST::card_list.do_tt();
}
/*--------------------------------------------------------------*/
void TTT::do_initial_dc()
{
	trace3("TTT::do_initial_dc", _sim->_Time0, _cont, _cont_dc);
	// set adp_nodes to initial values
	_sim->set_inc_mode_bad();
	CARD_LIST::card_list.tr_begin();
	_sim->_phase = p_INIT_DC;
	bool _converged = solve_with_homotopy(OPT::DCBIAS,TRANSIENT::_trace);
	assert(_converged); USE(_converged); // incomplete(); (why?)

//	 review();
	_accepted = true;
	TRANSIENT::accept();
	//
	trace0("TTT::do_initial_dc done");
	_sim->keep_voltages();
	assert(_sim->last_time() == 0.);
	_sim->set_inc_mode_yes();
	_sim->_mode = s_TTT;
}
/*--------------------------------------------------------------*/
void TTT::power_down(double until)
{
	_sim->_phase = p_PD;
	// do stress apply now, then tt advance and stress_last.

	_sim->force_tt_order(0);
	double time = until-_sim->_last_Time;
	trace2("TTT::power_down... ", _sim->_last_Time, until);

	/// only if _cont?
	//
	for(unsigned i=0; i<_sim->_adp_nodes;i++ ){
		trace2("TTT::power_down... ",i, _sim->_tt[i]);
	}
	for(unsigned i=0; i<_sim->_adp_nodes;i++ ) {
		trace2("TTT::power_down... ",i, _sim->_tr[i]);
	}

	_Time1 = _sim->_last_Time;
	_sim->_Time0 = _sim->_last_Time;
	assert(tt_iteration_number() == 0);
	assert(tt_iteration_number() == 0);

	_sim->_dT0 = 0;
	_sim->_dT1 = 0;
	_sim->_dT2 = 0;
	_sim->_time0 = 0;
	_sim->force_tt_order(0); assert(_sim->get_tt_order() == 0 );
	_time1 = 0.;

	advance_Time(); // fix last_iter time (dT0==0);
	assert(!_sim->_dT0);

	// bug. too early
	if(_trace>tDEBUG) { untested();
		outdata_b4(_sim->_Time0);
	}

	// start at age stored in tt.
	notstd::copy_n(_sim->_tt, _sim->_adp_nodes, _sim->_tt1);

	CARD_LIST::card_list.tt_advance();
//	CARD_LIST::card_list.do_tt(); // why?!

	for (uint_t ii = 1;  ii <= _sim->_total_nodes;  ++ii) {
		_sim->_nstat[_sim->_nm[ii]].set_last_change_time(0);
		_sim->_nstat[_sim->_nm[ii]].store_old_last_change_time();
		_sim->_nstat[_sim->_nm[ii]].set_final_time(0);
	}
	for (unsigned i=0; i<_sim->_adp_nodes; i++) {
		trace2("TTT::power_down... ",i, _sim->_tr[i]);
	}

	// NOW do the powerdown.

	_sim->zero_some_voltages();

	assert( _sim->_mode  == s_TTT );
	print_results_tt( _sim->_Time0 );

	trace0("TTT::power_down accept...");
	CARD_LIST::card_list.tr_accept();
	CARD_LIST::card_list.tr_stress_last(); // FIXME: merge into tt_review

	for (unsigned i=0; i<_sim->_adp_nodes; i++) {
		trace2("TTT::power_down... ",i, _sim->_tr[i]);
	}

	_sim->_dt0 = 0;
	_sim->_dT0 = time;
	_sim->_Time0 = until;

	advance_Time(); // check: move tr->tr1 correctly...
	_sim->force_tt_order(0);

	for(unsigned i=0; i<_sim->_adp_nodes;i++ ) {
		trace3("TTT::power_down advanced ",i, _sim->_tt[i], hp(_sim->_tt));
	}
	for(unsigned i=0; i<_sim->_adp_nodes;i++ ) {
		trace3("TTT::power_down advanced ",i, _sim->_tt1[i], hp(_sim->_tt1));
	}

	CARD_LIST::card_list.do_tt();

	trace0("TTT::power_down done");

	// go on with old voltages (good idea?)
	
	_sim->restore_voltages();
}
/*--------------------------------------------------------------------------*/
// FIXME:
// do a 2nd tt to check aging impact.
// ie. take some norm of node-voltage difference over time.
//
// printing (preliminary)
//   tALL,        /* all accepted time steps */
//   tGUESS    ,	/* show guess, (accepted only) */
//   tREJECTED  ,	/* show rejected time steps (+ guesses)
//   tITERATION ,	/* also show corrector iterations (not yet) */
//   tVERBOSE   ,	/* show more stuff? */
/*--------------------------------------------------------------------------*/
static void register_status();
/*--------------------------------------------------------------------------*/
void TTT::sweep()
{ untested();
	register_status();
	assert( _sim->_mode == s_TTT );
	assert( _sim->_Time0 >= 0 );
	_sim->_phase = p_NONE;
	_printed_steps = 0;
	::error(bLOG, "TTT sweep cont %d cont_tt %d new %d cont_dc %d\n", _cont, _cont_tt, _new, _cont_dc);
	::error(bLOG, "TTT sweep Tstart %f Tstop %f Tstep %f stepmode %d\n",
			(double)_Tstart, (double)_Tstop, (double)_Tstep, (int)_stepmode);
	assert(_sim->tt_iteration_number() == 0);

	//  sanitycheck();
	_sim->_tt_iter = 0;
	_sim->_tt_done = 0;

	head_tt(_tstart, _tstop, "TTime");

//	assert(!_new || !_cont_tt); ??

	if(_sim->_last_Time == 0. && !_new && !_cont_tt) {
		_new = 1;
	} else if (_sim->last_time() != 0. && !_cont_tt) {
//		_cont_dc = 1;
	} else {
	}

	if(_no_act) { untested();
		_out << "ttt " << string(_Tstart) << " " << string(_Tstop) << " " << string(_Tstep)
		     << " tran " << string(_tstop) << " " << string(_tstrobe) << "\n";
		_out << "stepmode " << _stepmode << "\n";
		return;
	}else{
	}

	_sim->force_tt_order(0); //hack

	if (_new) { untested();
		trace2("TTT::sweep tt_begin", _Tstart, _sim->_phase);
		_sim->_last_Time = _Tstart;
		tt_begin();
		// CARD_LIST::card_list.do_tt(); BUG: not idempotent (?!)
		set_step_cause(scUSER);
		assert(_sim->_mode == s_TTT);
	}else{ untested();
		tt_cont();
	}

	// bug: _new not correctly processed (?)

	if (_power_down) { untested();
		if (_new) {untested();
			print_results_tt(_sim->_Time0);
		}
		trace1("TTT::sweep pd until", _Tstop );
		power_down(_Tstop);
		print_results_tt(_Tstop);
		store_results(_Tstop);
		_sim->_last_Time = _Tstop;
		trace0("TTT::sweep pd done");
		return;
	}else if( _Tstop == _Tstart || double(_Tstop) == 0) {
		if (_new) {
			print_results_tt(_sim->_Time0);
		}
		trace0("TTT::sweep just apply");
		advance_Time();
		_sim->_tt_uic = 1;
		CARD_LIST::card_list.do_tt();
		_sim->_acceptq.clear();
		_sim->_tt_uic = 0;
		return;
	}else if (_Tstop == _Tstart) { untested();
		if (_new) {untested();
			print_results_tt(_sim->_Time0);
		}
		trace0("TTT::sweep just printing");
		outdata_b4(_sim->_Time0); // first output tt data
		return;
	}else if (!_cont_dc) {
		if (_new) {
			print_results_tt(_sim->_Time0);
		}
		// _cont_dc == skip initial dc
		trace0("TTT::sweep from 0");
		// print_head_tr(); // BUG
		do_initial_dc();

		// outdata_tt( _sim->_Time0); no.
		assert( _sim->_mode  == s_TTT );
		trace0("initial done");
		_cont = true;

		assert(_sim->_Time0 <= _Tstart);
		first();
		trace2("",_Time1,_Tstart);
		trace3("TTT::sweep first done", _sim->_Time0, _Time1, _Tstart);
		assert (_Time1 == _Tstart);
	}else if(_tstop==0. && _tstrobe==0. && (!_cont_tt)) { untested();
		if (_new) {untested();
			print_results_tt(_sim->_Time0);
		}
		print_head_tr();
		untested();
		do_initial_dc();
		return;
	}else if (_cont_tt) { untested();
		if (_new) {untested();
			print_results_tt(_sim->_Time0);
		}
		trace3("cont tt", _sim->last_Time(), _Tstart, _sim->_Time0);
		_sim->_Time0 = _sim->last_Time();
		_Tstart = _sim->last_Time();
		_sim->_dT0 = 0.;
		assert(_sim->_Time0 <= _Tstart);
		first();
	}else if (_cont_dc) {
		_sim->restore_voltages();
		if (_new) {
			print_results_tt(_sim->_Time0);
		}
		trace3("tt after tr?", _sim->last_Time(), _sim->last_time(), _sim->last_Time());
		_sim->_Time0 = _Tstart;
		_sim->_dT0 = _sim->last_time();
		//_sim->_dT1 = 0;
		_Time1 = _Tstart;
		assert(_sim->_Time0 <= _Tstart);
		first();
	} else { untested();
		error(bDANGER, "tt, unhandled. cont_tt %d cont_dc %d new %d %f %f\n",
				_cont_tt, _cont_dc, _new, _sim->last_time(), _sim->last_Time());
		unreachable();
	}

	// assert (_sim->_loadq.empty());
	assert(_sim->_Time0 == _Tstart);

	//_time_by_adp=NEVER;

	trace2("",_Time1,_Tstart);
	assert (_Time1 == _Tstart);
	assert( _sim->_mode  == s_TTT );

	while (next()) { untested();
		assert(step_cause());
		assert( _sim->_mode  == s_TTT );
		trace8( "TTT::sweep loop start ", _sim->_Time0, _Time1, _sim->_dT0,
				_accepted, _accepted_tt, tt_iteration_number(), _sim->_last_Time, _sim->last_time() );
		assert(_sim->_dT0 >= _sim->last_time());
		sanitycheck();

// 		_sim->_last_time = _tstop; // UARGH hack, should be implicit...
		trace2("TTT::sweep CARD::do_tt", _sim->_dT0, _sim->last_time());
		assert(  _sim->_dT0 - _sim->last_time() >= 0);

		_sim->_time0 = 0.;
		CARD_LIST::card_list.do_tt(); // before sweep_tr

		assert( _sim->_mode  == s_TTT );
		if(_trace >= tITERATION ) {
			outdata_b4(_sim->_Time0);
		}
		store_results_tt(_sim->_Time0); // first output tt data untested();
		assert( _sim->_mode  == s_TTT );
//		_sim->_time0 = _sim->last_time(); // time0 does not necessarily correspond to a
//		                                  // computed step.
//		assert( _sim->_time0 <= _sim->_dT0);

		assert(step_cause());
		TTT::sweep_tr();
		assert(step_cause());

		if (!_accepted) { incomplete();
			assert(_sim->_mode==s_TTT);
			_accepted_tt = false;
		} else {
			_accepted_tt = review();
		}

		trace4("", _sim->_Time0, _Time_by_user_request, _tstop, step_cause());
		// assert (_sim->_Time0 <= _Time_by_user_request + _tstop * 1.0001 );

		if (!_accepted_tt) {
			error(bDEBUG, "rejected step at %f, want %f, delta %E\n",
					_sim->_Time0, _Time_by_error_estimate,
					_sim->_Time0 - _Time_by_error_estimate);
			_sim->_tt_rejects++;
			_sim->_tt_rejects_total++;
			if (_trace >= tDEBUG) { untested();
				print_stored_results_tt(-_sim->_Time0);
			}

			if (_trace >= tREJECTED) {
				TRANSIENT::_out << "#r\n";
				assert(step_cause());
				outdata_tt(_sim->_Time0);
				assert( _sim->_mode  == s_TTT );
			}

		} else {
			// accepted_tt
			_sim->_time0 = _sim->last_time(); // time0 does not necessarily correspond to a
		                                     // computed step.
			_sim->keep_voltages(); // assume that v0 is still what we want to keep...
			assert(0.99999* _sim->last_time() <= _tstop);
			if (_trace == tGUESS) {
				// in verbose mode, stored data have been printed already...
				print_stored_results_tt(_sim->_Time0);
			}
			TTT::accept();
			assert( _sim->_mode  == s_TTT );
		}

		if (!_accepted_tt) {
//		} else if (_sim->_last_Time + 0.00001 > _Tstop) {
//			// good idea? better leave it to next() ...
//			_Time_by_user_request = _Tstop;
		} else if (_sim->_Time0 >= _Time_by_user_request // incomplete
				|| step_cause() == scUSER || step_cause() == scLAST) {
			trace3("TTT::sweep user step", _sim->_Time0, _sim->last_time(), _Time_by_user_request);
			assert(step_cause() == scUSER || step_cause() == scINITIAL || step_cause() == scLAST);
			outdata_tt(_sim->_Time0 );
			_printed_steps++;
			switch (_stepmode) { //
				case tts_LIN:
					_Time_by_user_request = _sim->_Time0 + _Tstep;
					// _out << "LIN" << _Time_by_user_request << " " << _Tstep << "\n";
					break;
				case tts_MUL:
					_Time_by_user_request = _sim->_Time0 + _sim->last_time() * pow((double)_Tstep, _printed_steps);
					break;
				default:
					unreachable();
			}
			assert(_Time_by_user_request>=_sim->_last_Time);
			// if(_Time_by_user_request > _Tstop - 2.*_tstop){ untested();
			// 	_Time_by_user_request = _Tstop - _tstop;
			// }else if(_Time_by_user_request > _Tstop - _tstop){ untested();
			// 	_Time_by_user_request = _Tstop - _tstop;
			// }else{
			// }
		} else if(_trace >= tALLTIME) { untested();
			assert(step_cause());
			outdata_tt(_sim->_Time0);
		} else { untested();
			trace1("TTT::sweep not printing", _sim->_Time0);
			TTT::store_results( _sim->_Time0 + _tstop );
		}

		if (!_accepted_tt) { untested();
		} else { untested();
			_Time_by_user_request = min(_Time_by_user_request, _Tstop-_sim->last_time());

			_sim->_last_Time = _sim->_Time0 + _sim->last_time();
			trace5("TTT::sweep bottom loop", _sim->_last_Time, _sim->_Time0, _Time_by_user_request, _stepmode, _Tstep );
			assert( _sim->_mode  == s_TTT );
		}
	}
	assert( _sim->_mode  == s_TTT );

//		_out << "* TTT::sweep =================== endof loop "<<_sim->_last_Time<<"\n";
	_sim->_Time0 = _sim->_last_Time;
	_sim->_dT0 = 0;

	advance_Time();

#ifndef NDEBUG
	/// invalidate history
	std::fill_n(_sim->_tr1, _sim->_adp_nodes, NAN);
	std::fill_n(_sim->_tt1, _sim->_adp_nodes, NAN);
#endif
	_sim->set_command_tt();
}
/*--------------------------------------------------------------------------*/
void TTT::sweep_tr() // tr sweep wrapper.
{
	_sim->_time0 = 0.;// for now?
	_sim->keep_voltages(true); // also sets _last_time=_time0
	unsigned hidden = ::status.hidden_steps;
	int ttcontrol = ::status.control;
	_sim->_mode = s_TRAN;
	//print_head_tr(); wrong. not b4 accept
	trace4("TTT::sweep_tr", storelist().size() , _sim->tt_iteration_number(), ttcontrol, _tstop);
	typedef std::map<std::string, WAVE > wavemap_t;
	for ( wavemap_t::iterator i = _sim->_waves["tran"].begin(); i != _sim->_waves["tran"].end(); ++i) {
		i->second.initialize();
	}
	_sim->_mode = s_TTT;

	trace1("TTT::sweep_tr() ", _inside_tt);
	_sim->_time0 = 0.0;

	try {
		_cont = true;
		trace4("TTT::sweep calling sweep", _sim->_time0, _sim->last_time(), _tstrobe, _sim->_phase);
		_inside_tt = true;
		assert(_sim->_mode == s_TTT);
		TRANSIENT::sweep();
		assert(_accepted);
		assert( _sim->last_time() <= _tstop + _sim->_dtmin );

	}catch (Exception& e) { untested();
		assert(_sim->_mode == s_TTT);
		error(bDANGER, "Sweep exception \"%s\" at %E, dT0 %E, step %i\n",
				e.message().c_str(), _sim->_Time0, _sim->_dT0, tt_iteration_number());
		assert(_sim->_dT0 > _sim->last_time()); // = means zero step
		                                       // if this fails, something is wrong.
		_accepted=_accepted_tt=false;
		::status.review.stop();
//		_sim->invalidate_tt();
		// throw(e); go on with smalller step?
		if (!tt_iteration_number()) { untested();
			_out << "* first one must not fail\n";
			throw(e);
		}
	}

	if(_sim->_tt_accepted && _agemode == amONCE){
		error(bNOERROR, "skipping tr_stress_last\n");
	}else{
		CARD_LIST::card_list.tr_stress_last(); // FIXME: merge into tt_review
	}

	_cont_dc = false;
	_sim->_mode = s_TTT;

	// restore tt status
	::status.control = ttcontrol;
	::status.hidden_steps = hidden;

	_sim->pop_voltages();
   _sim->_time0 = _sim->last_time();
	_sim->keep_voltages();
	trace4( "TTT::sweep_tr done", _sim->_Time0, _sim->_time0, _sim->last_time(), _tstop);
	assert(_sim->_last_time *0.9999 <= _tstop);
}
/*--------------------------------------------------------------------------*/
void TTT::accept()
{

	ADP_NODE_LIST::adp_node_list.do_forall( &ADP_NODE::tt_accept);

	//  FIXME:  _tt_acceptq -> corrector?
	if(_sim->_tt_accepted && _agemode == amONCE){
		error(bNOERROR, "skipping tt_accept\n");
	}else{
		while (!_sim->_tt_acceptq.empty()) {
			_sim->_tt_acceptq.back()->tt_accept();
			_sim->_tt_acceptq.pop_back();
		}
	}
	_sim->_tt_accepted++;
	//  _sim->_last_Time = _sim->_Time0+_tstop;
}
/*--------------------------------------------------------------------------*/
double TTT::time_by_voltages()
{ untested();
	// we started at vdc (unchanged), before keep, compare v0<->vdc

	double d=0;
	USE(d);
	for (unsigned i=1; i<=_sim->_total_nodes; ++i) { untested();
		double delta = fabs(_sim->_v0[i]-_sim->vdc()[i]);
		double sum = fabs(_sim->_v0[i])+fabs(_sim->vdc()[i]);
		USE(sum);

		double err = delta ;
		USE(err);
	}

	incomplete();
	return(0);

}
/*-----------------------------------------------------------*/
// need to check that boundaries are consistent.
bool TTT::conchk() const
{ untested();
	double a = OPT::abstol;
	double r = OPT::reltol;

	for (unsigned i=1; i<=_sim->_total_nodes; ++i) { untested();
		double n = fabs(_sim->_v0[i]);
		double o = fabs(_sim->vdc()[i]);

		if (std::abs(n-o) > (r * std::abs(n) + a))
			return false;
	}
	return true;
}
/*--------------------------------------------------------------------------*/
bool TTT::review()
{
	sanitycheck();

	assert(_sim->_Time0 >= _sim->_dT0);
	assert(_sim->_Time0 - _sim->_dT0 >=0 );

	TIME_PAIR Time_by = CARD_LIST::card_list.tt_review();
	_Time_by_error_estimate = Time_by._error_estimate;

	if(Time_by._event < _Time1 + _sim->last_time()){
		_time_by_ambiguous_event = _Time1 + _sim->last_time();
	}else{
		_time_by_ambiguous_event = Time_by._event;
	}

	trace6("TTT::review", _sim->_Time0, _sim->_dT0, _Time_by_error_estimate, _time_by_ambiguous_event, _sim->last_time(), _Time1);

	_dT_by_adp = (ADP_NODE_LIST::adp_node_list.tt_review())._event;

	//  double _dT_by_nodes = time_by_voltages();// 1.0 / (1e-20+voltage_distance());//

	// BUG. don't use _Time1 here.
	_time_by_adp = _Time1 + _dT_by_adp;

	if ((_dT_by_adp < _tstop) && (_tstop == _sim->_dT0 )) { untested();
		return true;
	}

	if( _time_by_adp < _sim->_Time0 ){ untested();
		//    _out << "* tt reject (adp) timestep " << _sim->_dT0 << " want dT " << _dT_by_adp << " \n";
		trace3( "TTT::review_tt: reject adptime", _sim->_Time0, _dT_by_adp , _sim->_dT0 );
	}

	bool acc = _Time_by_error_estimate + .01*_sim->_dtmin >= _sim->_Time0
           && _time_by_ambiguous_event >= _sim->_Time0;

	if(step_cause() == scLAST && !acc) {
		// endless loop?
		trace3( "TTT::review_tt: endless loop", _sim->_Time0, _sim->_dT0, _Tstop);
		acc = true;
	}else{
	}
	return acc;
}
/*--------------------------------------------------------------------------*/
void TTT::do_it(CS& Cmd, CARD_LIST* Scope)
{
	_scope = Scope;
	_sim->set_command_tt();
	_sim->set_label("tt");
//	_sim->_age = false; // a tr option
	::status.ttt.reset().start();
	::status.tt_tries = 0;
	command_base(Cmd); // s__init.cc
	assert(_sim->_mode==s_TTT);
	::status.ttt.stop();

	::status.four.stop(); // ?
	::status.total.stop(); // ?
}
/*--------------------------------------------------------------------------*/
void TTT::probeexpand()
{
	trace0("TTT::probeexpand");
	_sim->_mode = s_TRAN; // BUG: use phase.
	PROBELIST* transtore = &_probe_lists->store[s_TRAN];

	// append prior prints to store (aligned)
	for (PROBELIST::const_iterator p=printlist().begin();
			p!=printlist().end();  ++p) {
		PROBE* x=((*p)->clone());
		transtore->push_probe(x);
	}
	
	// additionally store during TR as user wishes. (if necessary)
	for (PROBELIST::const_iterator p=oldstore.begin();
			p!=oldstore.end();  ++p) {
		PROBE* x=((*p)->clone());
		transtore->merge_probe(x);
	}

	_sim->set_command_tt();

	PROBELIST measstore;
	for (PROBELIST::const_iterator
			p=printlist().begin();  p!=printlist().end();  ++p) {
		MEAS_PROBE* w = dynamic_cast<MEAS_PROBE*>(*p);
		if (w) {
			w->expand();
			string probe_name = w->probe_name;
			trace1("TTT::probeexpand", probe_name);
			CS* c = new CS(CS::_STRING,probe_name);
			measstore.add_list(*c);
			delete c;
		} else {
		}
	}

// merge some more measure nodes.
	for (PROBELIST::const_iterator p=measstore.begin();
			p!=measstore.end();  ++p) {
		PROBE* x=((*p)->clone());
		transtore->merge_probe(x);
	}
}
/*--------------------------------------------------------------------------*/
/* allocate:  allocate space for tt
*/
void TTT::allocate()
{
	int probes = printlist().size();
	assert(!_tt_store);
	trace1("TTT::allocate allocating probes ", probes);
	_tt_store = new double[probes]; // shadows the .print tt probes

	PROBELIST* transtore = &_probe_lists->store[s_TRAN];

	for (PROBELIST::const_iterator p=transtore->begin();
			p!=transtore->end(); ++p) {
		oldstore.push_probe((*p)->clone());
	}
	transtore->clear();

	probeexpand();

	assert(_sim->_mode==s_TTT);

	_sim->set_command_tran();

	_sim->_waves["tran"].clear();

	if(_wavep){
		delete[] _wavep;
	}else{
	}
	_wavep = new WAVE*[_probe_lists->store[s_TRAN].size()];

	unsigned ii = 0;
	for (PROBELIST::const_iterator
			p=storelist().begin();  p!=storelist().end();  ++p) {
		string l = (*p)->label();
		if (OPT::case_insensitive) {
			notstd::to_upper(&l);
		}
		_wavep[ii++] = &(_sim->_waves["tran"][l]);
	}
	_sim->set_command_tt();

	trace1("TTT::allocate allocated transtore waves", _probe_lists->store[s_TRAN].size());
}
/*--------------------------------------------------------------------------*/
/* unallocate:  unallocate space for tt
 * (merge further up?)
*/
void TTT::finish()
{
	// int probes = printlist().size();
	assert (_tt_store);
	trace2("TTT::unallocate freeing _tt_store", printlist().size(), storelist().size());
	// incomplete(); //unallocate at deconstruction...
	// make consecutive sims faster.
	delete[] _tt_store;
	_tt_store = NULL;

	_probe_lists->store[s_TRAN].clear();
	//FIXME: delete waves;

	//  if (_fdata_tt) { untested();
	//    for (int ii = 0;  ii < printlist().size();  ++ii) { untested();
	//      delete [] _fdata_tt[ii];
	//    }
	//    delete [] _fdata_tt;
	//    _fdata_tt = NULL;
	//  }else{itested();
	//  }
	//  TRANSIENT::finish?
	_sim->unalloc_vectors();
	_sim->_lu.unallocate();
	_sim->_aa.unallocate();

	PROBELIST* transtore = &_probe_lists->store[s_TRAN];

	transtore->clear();
	// insert previous store probes
	for (PROBELIST::const_iterator p=oldstore.begin();
			p!=oldstore.end(); ++p) {
		transtore->push_probe((*p)->clone());
	}

	// move waves back (hack, workaround)
//	trace1("TTT::unallocate deleting", _sim->_waves);
}
/*--------------------------------------------------------------------------*/
static TTT p10;
DISPATCHER<CMD>::INSTALL d10(&command_dispatcher, "twotimetran|ttr", &p10);
/*--------------------------------------------------------------------------*/
static void register_status()
{
	static bool done;

	if(!done) {
		new DISPATCHER<CKT_BASE>::INSTALL(&status_dispatcher, "tt", &p10);
		done = true;
	}
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
double behaviour_timestep()
{ untested();
	return 1000000;
}
/*--------------------------------------------------------------------------*/
bool TTT::next()
{
	if(_Time_by_user_request > _Tstop){ untested();
		return false;
	}else{ untested();
	}
	assert( _sim->_mode  == s_TTT );
	double new_dT;
	double new_Time0;
	double refTime = _Time1;
	_inside_tt = true;

	trace7( "TTT::next()", tt_iteration_number(), _sim->_Time0, _Time1,
			_sim->_last_Time, _accepted_tt, _Time_by_user_request, _Time_by_error_estimate);

	assert(_sim->_Time0 >=0 || tt_iteration_number() == 0);
	STEP_CAUSE new_control = scNO_ADVANCE;

	if (_sim->_dT0 == 0.) { unreachable();
		// first transient. is this reachable?
		new_dT = 0;
	}else if (_Time1 == _sim->_Time0) { untested();
		// second step...
		new_dT = _sim->last_time();
		new_control = scINITIAL;
	} else if (!_accepted_tt) {
		error(bTRACE, "rejected step %f... at %f\n", _sim->_dT0, _sim->_Time0);

		new_dT = 0.9 * ( _sim->_dT0 - _sim->last_time() ) + _sim->last_time(); // back at least this much.
		if ( _time_by_adp < _sim->_Time0 ){ untested();
			assert(_time_by_adp >= 0);
			new_dT = fmin(new_dT, _dT_by_adp);
		} else {
		}
		if ( _Time_by_error_estimate < _sim->_Time0 ){
			assert(_Time_by_error_estimate >= 0);
			new_dT = fmin(new_dT, _Time_by_error_estimate - _Time1);
			assert(new_dT>=0);
		} else {
		}
		if (!_accepted) {
			new_dT = fmin(new_dT,.5*_sim->_dT0);
		} else {
		}

		assert( is_number(new_dT));
		new_control = scREJECT;
	}else{
		// accepted step. calculating new_dT
		refTime = _sim->_Time0;

		new_dT = _Time_by_user_request - refTime;
		new_control = scUSER;

		if (new_dT > _dT_by_adp) { untested();
		}
		if (new_dT > _Time_by_error_estimate - refTime) {
			new_control = scTE;
			new_dT = _Time_by_error_estimate - refTime;
		}else{
		}
		if(OPT::ttstepgrow == 0) {
			// fixme: allow "inf" in u_opt
		}else if(new_dT > _sim->_dT0 * OPT::ttstepgrow) {
			new_control = scGROW;
			new_dT = _sim->_dT0 * OPT::ttstepgrow;
		} else {
		}

		// new_dT = max( new_dT, (double) _tstop );
		new_dT = max( new_dT, (double) _sim->last_time() );
		assert(new_dT >= 0.);

		if (new_dT < _dTmin) {
			error(bTRACE, "step too small %e %e adp %e\n", new_dT, _dTmin, _dT_by_adp );
		} else {
		}

		// last step handler, snap to edge

		assert( !OPT::ttstepgrow || new_dT < (OPT::ttstepgrow+1) * _sim->_dT0 || _sim->_dT1<=0 );

		++::status.hidden_steps_tt;
	} // accepted

	double newTime = refTime + new_dT;
	trace8("TTT::next", _time_by_ambiguous_event, _Time1, _sim->_last_Time, _tstop, _Tstart, newTime, _Time_by_user_request, _sim->last_time());
	assert(_time_by_ambiguous_event > _Time1);

	if (!_evt) { untested();
		new_dT = newTime - refTime; // bug?
		assert(new_dT >= 0.9999 * _sim->last_time());
	}else if (_time_by_ambiguous_event < newTime - _sim->_dTmin) {
		if (_time_by_ambiguous_event < _Time1 + 2*_sim->last_time()) {
			double minTime = _Time1 + 2*_sim->last_time();
			if (newTime - _sim->last_time() < minTime) {
				newTime = minTime;
				new_control = scAMBEVENT;
			}else{
			}
		}else{
			newTime = _time_by_ambiguous_event;
			new_control = scAMBEVENT;
		}
		trace5("ambevent", _time_by_ambiguous_event, _sim->_Time0, newTime, refTime, _accepted_tt);
		new_dT = newTime - refTime;
		assert(new_dT >= 0.);
	}else{
	}

	assert(step_cause());

	new_Time0 = _sim->_last_Time + new_dT - _sim->last_time();
	trace6("TTT::next ", new_dT, _tstop, new_dT, _Time_by_user_request, new_Time0, refTime);
	assert(new_Time0 > _Time1);
	_Time1 = refTime;
	assert(new_Time0 > _Time1);

	// hmm rethink. maybe don't shift forwards if not necessary....
	if ( (double)_Tstop - _sim->last_time() <= _sim->_last_Time) {
		//	new_control = scLAST;
	}else if ( ( (double)_Tstop - _sim->_last_Time < 2.00001 * _tstop
				|| new_dT + _sim->_last_Time > _Tstop ) ) {
		// snap forward...
		new_dT = _Tstop - _sim->_Time0 - _sim->last_time();
		assert(new_dT > 0.999*_sim->last_time()); //noise?
		new_Time0 = _sim->_last_Time + new_dT - _sim->last_time();
		new_control = scLAST;
	}else if(new_Time0 + _sim->last_time() >= _Time_by_user_request){
		new_control = scUSER;
	}else{
	}
	set_step_cause(new_control);

	if (_sim->_dT1 && new_dT < _sim->last_time() ) {
		new_dT = _sim->last_time();
		assert(new_dT > 0.999*_sim->last_time()); //noise?
		set_step_cause(scSMALL);
	}else{
	}

	bool another_step = true;

	if (! ( _Tstop - _sim->_last_Time >= _sim->last_time() )){
		another_step = false;
	}else if (! (  new_Time0 <= _Tstop - _sim->last_time()   )) {untested();
		another_step = false;
	}else if (! ( new_dT != 0 || !_sim->_dT0  )){ untested();
		another_step = false;
	}

	trace7("TTT::next", another_step, _sim->last_time(), new_Time0, new_dT, _sim->last_Time(), _Time1, step_cause());

	assert(new_Time0 >= new_dT);
	assert(new_Time0 > _Time1);

	if(_accepted_tt){ // fixme: obenhin, oder?
		_sim->_dT3 = _sim->_dT2;
		_sim->_dT2 = _sim->_dT1;
		_sim->_dT1 = _sim->_dT0;
	}

	if (!another_step) {
		trace6( "TTT::next no next @ Time0: " , _sim->_Time0,  _sim->_dT0, new_dT, _dTmin, _tstop, new_Time0 );
		trace5( "TTT::next no next @ Time0: " ,\
				_Time1, _Tstop, new_dT >= _dTmin ,  _Tstop - _Time1 >= _dTmin ,  new_Time0 + _tstop <= _Tstop );
		//    _sim->_last_Time = _Time1 + _tstop; // FIXME

		return (false); // bug
	} else {
		trace8("TTT::next another step ", _sim->_Time0, new_dT, _Time1, _Tstop,
				_sim->_dT0, new_Time0, _Time_by_user_request, _dTmin );
	}

	assert(new_dT > 0.999*_sim->last_time()); //noise?
	if (new_dT < _sim->last_time()) {
		// BUG
		new_dT = _sim->last_time();
	}else{
	}
	_sim->_dT0 = new_dT;
	_sim->_Time0 = new_Time0;

	assert( _sim->_Time0 >= _sim->_dT0*0.9999 );
	if (_sim->_Time0 < _sim->_dT0) {
		// BUG
		_sim->_Time0 = _sim->_dT0;
	}else{
	}

// 	assert(_sim->_dT0 >= _dTmin);
	assert(_sim->_dT0 + 0.001*_dTmin >= _sim->last_time());
	assert(_sim->_dT0 >= _sim->last_time());

	_time1 = 0.;

	assert( _sim->_Time0 < 1+ _Time1 + _sim->_dT0 );
	assert( _sim->_Time0 >= _sim->_dT0 );
	assert( _sim->_Time0 > 0  || tt_iteration_number() == 0);

	_sim->restore_voltages();

	trace0("praparing nodes nodes");

	for (uint_t ii = 1;  ii <= _sim->_total_nodes;  ++ii) {
		_sim->_nstat[_sim->_nm[ii]].set_last_change_time(0);
		_sim->_nstat[_sim->_nm[ii]].store_old_last_change_time();
		_sim->_nstat[_sim->_nm[ii]].set_final_time(0);
	}

	advance_Time();
	++steps_total_tt;
	++::status.hidden_steps;

	trace8( "TTT::next: exiting next: " , _sim->tt_iteration_number(),
			_sim->_Time0, _Time1, _sim->_dT0, _sim->_dT1 , _sim->_dT2, _sim->_last_Time, step_cause() );
	assert(_sim->_Time0 > _Time1);
	// assert(_sim->_dT0 >= _dTmin);
	assert(step_cause());
	return another_step;
} // next
/*--------------------------------------------------------------------------*/
/* SIM::head: print column headings and draw plot borders
*/

void TTT::head_tt(double start, double stop, const std::string& col1)
{
	trace2("TTT::head_tt", start, stop);
	assert(_sim->_mode==s_TTT);
	//PROBELIST* transtore = &_probe_lists->store[s_TRAN];

	print_tr_probe_number = printlist().size();
	{
		trace1("TTT::head tr ttt WAVE", storelist().size() );

		if (_wavep_tt) {
			delete[]_wavep_tt;
		} else {
		}
		_wavep_tt = new WAVE*[storelist().size()];
		unsigned ii = 0;
		for (PROBELIST::const_iterator
				p=storelist().begin();  p!=storelist().end();  ++p) {
			string l = (*p)->label();
			if (OPT::case_insensitive) {
				notstd::to_upper(&l);
			}
			_wavep_tt[ii] = &(_sim->_waves[_sim->_label][l]);
			_wavep_tt[ii++]->initialize();
		}
	}

	if (!plopen(start, stop, plotlist())) {
		// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		int width = std::min(OPT::numdgt+5, BIGBUFLEN-10);
		char format[20];
		//sprintf(format, "%%c%%-%u.%us", width, width);
		sprintf(format, "%%c%%-%us", width);
		_out.form(format, '#', col1.c_str());

		for (PROBELIST::const_iterator
				p=printlist().begin();  p!=printlist().end();  ++p) {
			_out.form(format, ' ', (*p)->label().c_str());
			//trace1("TTT::head_tt", (*p)->label() );
			if ( !(--print_tr_probe_number) ) break;
		}
		_out << '\n';
	}else{ untested();
	}
	_sim->_mode=s_TRAN;


	trace3("TTT::tt_head probe TRAN", printlist().size(), storelist().size(), oldstore.size());

//	_sim->_waves = new WAVE[storelist().size()];

	_sim->_mode=s_TTT;
	for (PROBELIST::const_iterator
			p=printlist().begin();  p!=printlist().end();  ++p) {
		_sim->_mode=s_TRAN;
		(*p)->precalc_last();
		_sim->_mode=s_TTT;
	}
}
/*--------------------------------------------------------------------------*/
/* SIM::head: initialize waves, override TRANSIENT::head()
*/
void TTT::head(double /* start */ , double /* stop */, const std::string& )
{
	trace0("TTT::head");
	assert(_sim->_mode==s_TTT);
	_sim->_mode=s_TRAN;
	if (_sim->tt_iteration_number() == 0) {
		TRANSIENT::_out.setfloatwidth(OPT::numdgt, OPT::numdgt+6);
		// print_head_tr();
	}
	trace3("TTT::head probe TTT", printlist().size(), storelist().size(), _sim->_phase);

	unsigned ii = 0;
	for (PROBELIST::const_iterator
			p=storelist().begin();  p!=storelist().end();  ++p) {
		string l = (*p)->label();
		if (OPT::case_insensitive) {
			notstd::to_upper(&l);
		}
		_wavep[ii] = &(_sim->_waves["tran"][l]);
		_wavep[ii++]->initialize();
	}

	_sim->_mode = s_TTT;
}
/*--------------------------------------------------------------------------*/
/* SIM::print_head: print column headers TR
*/
void TTT::print_head_tr()
{
	trace0("TTT::print_head_tr");

	//assert( _sim->_mode==s_TRAN );
	_sim->_mode = s_TRAN;
	SIM_MODE oldmode = _sim->_mode;

	if (!printlist().size()) {
		// incomplete(); // bug: what to do with empty printlists?
		                 // for now: suppress output...
		_sim->_mode = oldmode;
		return;
	}
	TRANSIENT::_out << "#TTime ";
	{
		// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		int width = std::min(OPT::numdgt+5, BIGBUFLEN-10);
		char format[20];
		sprintf(format, "%%c%%-%us", width);
		for (int i=0; i< OPT::numdgt; i++) {
			TRANSIENT::_out << " ";
		}
		TRANSIENT::_out << "time";
		for (int i=0; i< OPT::numdgt+1; i++) {
			TRANSIENT::_out << " ";
		}
		for (PROBELIST::const_iterator
				p=printlist().begin();  p!=printlist().end();  ++p) {
			TRANSIENT::_out.form(format, ' ', (*p)->label().c_str());
		}
		TRANSIENT::_out << '\n';
	}
	_sim->_mode = oldmode;
	TRANSIENT::_out.setfloatwidth(OPT::numdgt, OPT::numdgt+6);
}
/*--------------------------------------------------------------------------*/
void TTT::print_foot_tr()
{
	// make plotter happy. (only if _out!=TRANSIENT::_out?)
	_sim->set_command_tran();
	if (printlist().size()) {
		TRANSIENT::_out << '\n';
	}
}
/*--------------------------------------------------------------------------*/
// override TRANSIENT::outdata
// save things during TR.
void TTT::outdata(double time0, int x)
{ untested();
	assert(TRANSIENT::step_cause());
	assert(is_number(time0));
	_sim->_mode = s_TTT;
	::status.output.start();
	assert( _sim->_mode  == s_TTT );

	// SIM::alarm();
	_sim->_mode=s_TRAN;
	if ( _trace>=tDEBUG && ( x & ofPRINT ) ) { untested();
		TRANSIENT::print_results(time0);
	}else{ untested();
	}

	_sim->set_command_tran();
	if(printlist().size()==0){ untested();
	}else if (_sim->tt_iteration_number()==0) { untested();
		// will always accept 1st
		if( x & ofPRINT ){ untested();
			TRANSIENT::_out << (double)_sim->_Time0;
			TRANSIENT::print_results(time0);
		}else{untested();
		}
	} else { untested();
		// store_results(time0);
		if (TRANSIENT::_trace >= tDEBUG){ untested();
			TRANSIENT::_out << "*" << (double)_sim->_Time0;
			TRANSIENT::print_results(time0);
		}
	}
	_sim->set_command_tt();

	_sim->_mode=s_TRAN;
	// FIXME (only > 0)?
	if( x & ofPRINT ){ untested();
		TRANSIENT::store_results(time0);
	}else{ untested();
	}
	if( x & ofKEEP ){ untested();
		_sim->keep_voltages();
	}

	_sim->_mode = s_TTT;
	_sim->reset_iteration_counter(iPRINTSTEP);
	::status.hidden_steps = 0;
	::status.output.stop();
	assert( _sim->_mode  == s_TTT );
}
/*--------------------------------------------------------------------------*/
void TTT::outdata_b4(double time)
{ untested();
	assert( _sim->_mode  == s_TTT );
	::status.output.start();
	print_results_tt(time);
	::status.output.stop();
}
/*--------------------------------------------------------------------------*/
// print at end of timeframe. "now" is begin of timeframe...
void TTT::outdata_tt(double now)
{ untested();
	assert(step_cause());
	assert(_sim->_mode==s_TTT);
	if ( 0 && _accepted_tt && _sim->_dT0 && !is_almost (_sim->_dT0 + _sim->_last_Time, now + _sim->last_time() )) { untested();
		error(bWARNING, "EOF: %.9f, last_Time: %.9f, dT0: %f, now: %f\n", _sim->_Time0+_sim->last_time(), _sim->_last_Time, _sim->_dT0,now );
	}
	::status.output.start();
	assert( _sim->_mode == s_TTT );
	print_results_tr(now); //transient print?
	assert( _sim->_mode == s_TTT );
	print_results_tt(now + _sim->last_time());
	tt_alarm(_alarm, &_out);
	assert(_sim->_mode==s_TTT);
	if (!tt_iteration_number()) {
		TTT::store_results( now ); // store extra results.
	}
	TTT::store_results( now + _sim->last_time() ); // results are stored after TRAN
	assert(_sim->_mode==s_TTT);
	_sim->reset_iteration_counter(iPRINTSTEP);
	::status.hidden_steps = 0;
	::status.output.stop();
}
/*--------------------------------------------------------------------------*/
// fixme: merge back to s__out.cc
void TTT::tt_alarm(ALARM a, OMSTREAM* out)
{
	if (a==aNONE) {
		return;
	}
	if (!out) { untested();
		out = &_out;
	}
	bool abort=false;
	_out.setfloatwidth(OPT::numdgt, OPT::numdgt+6);
	for (PROBELIST::const_iterator
			p=alarmlist().begin();  p!=alarmlist().end();  ++p) {
		if (!(*p)->in_range()) {
			stringstream s;
			s << (*p)->label() << '=' << (*p)->value() << '\n';
			switch (a) {
				case aREDIR:
					cerr << s.str();
					break;
				case aABORT:
					abort = true;
				default:
					*out << string(s.str());
			}
		}
	}
	if (abort && _accepted_tt) {
		throw Exception_Alarm("values out of range\n");
	}
}
/*--------------------------------------------------------------------------*/
// print during tr. i.e. print from storelist
void TTT::print_results(double time)
{ untested();
	// deprecated call.
	assert(false);
	unreachable();
	print_results_tr(time);
}
/*--------------------------------------------------------------------------*/
void TTT::print_results_tr(double Time ) // Time is begin of frame
{
	SIM_MODE oldmode = _sim->_mode;
	trace4("TTT::print_results_tr()", tt_iteration_number(), Time, _sim->_mode, oldmode);
	USE(Time);
	_sim->set_command_tran();
	const WAVE* w = NULL;

	if ( printlist().size() == 0 ) {
		_sim->_mode = oldmode;
		return;
	}
	assert(_sim->_mode=s_TRAN);

	if (!IO::plotout.any() && _sim->tt_iteration_number() > 0 ) {
		// 1st step has already been printed (no need for review)
		TRANSIENT::_out.setfloatwidth(OPT::numdgt, OPT::numdgt+6);
		std::map<std::string, WAVE>::const_iterator wi = _sim->_waves["tran"].begin();
		assert(wi!=_sim->_waves["tran"].end());
		w = &wi->second;

		print_head_tr();

		int ii=0;
		// hmm correct?
		WAVE::const_iterator* myiterators = new WAVE::const_iterator[printlist().size()];

		for (PROBELIST::const_iterator p = printlist().begin();
				p!=printlist().end();  ++p) {
			myiterators[ii] = _wavep[ii]->begin();
			ii++;
		}
		for (WAVE::const_iterator i = w->begin(); i != w->end(); i++) {
			TRANSIENT::_out << _sim->_Time0;
			if( i->first > _sim->last_time()) { untested();
				break;
			}else{
			}
			TRANSIENT::_out << i->first;

			ii=0;
			for (PROBELIST::const_iterator
					p=printlist().begin();  p!=printlist().end();  ++p) {

				TRANSIENT::_out << myiterators[ii]->second;
				myiterators[ii]++;
				ii++;
			}

			TRANSIENT::_out << '\n';
		}
		delete[] myiterators;
		print_foot_tr();
	} else {
	}
	_sim->set_command_tt(); // FIXME
	_sim->_mode = oldmode;
	trace3("TTT::print_results_tr done", tt_iteration_number(), Time, _sim->_mode);
}
/*--------------------------------------------------------------------------*/
void TTT::print_stored_results_tt(double x)
{

	if ( printlist().size() ==0) { untested();
		trace0("no ttprint");
		return;
	} else {
	}
	assert( _sim->_mode  == s_TTT );
	if (!IO::plotout.any()) {
		assert(x != NOT_VALID);
		int i;
		_out << x;
		for (i=0; i <  printlist().size(); i++ ) {
			_out << _tt_store[i];
		}
		_out << '\n';
	} else { untested();
	}
}
/*--------------------------------------------------------------------------*/
// store things at begin of timeframe, so it may be printed
// after/if the step is accepted.
void TTT::store_results_tt(double x)
{ untested();
	trace0("TTT::store_results_tt()");
	if ( printlist().size() ==0) {
		return;
	}

	int i=0;
	assert( _sim->_mode  == s_TTT );
	if (!IO::plotout.any()) {
		assert(x != NOT_VALID); USE(x);
		for (PROBELIST::const_iterator
				p=printlist().begin();  p!=printlist().end();  ++p) {
			_tt_store[i++] =  (*p)->value();
		}
	} else { untested();
	}

}
/*--------------------------------------------------------------------------*/
void TTT::print_results_tt(double x)
{
	if( printlist().size() == 0 && _quiet ){ untested();
		return;
	}

	assert( _sim->_mode  == s_TTT );
	if (!IO::plotout.any()) {
		_out.setfloatwidth(OPT::numdgt, OPT::numdgt+6);
		assert(x != NOT_VALID);
		_out << x;
		for (PROBELIST::const_iterator
				p=printlist().begin();  p!=printlist().end();  ++p) {

			if (!(*p)) { unreachable(); // ??
				continue;
			}

			try {
				_out << (*p)->value();
			} catch( Exception_No_Match ) { untested();
				_out << "XXXX";
			}
		}
		_out << '\n';
	}else{ untested();
	}
}
/*--------------------------------------------------------------------------*/
string TTT::status()const
{
	return "twotime timesteps: accepted=" + ::to_string(steps_accepted())
		+ ", rejected=" + ::to_string(steps_rejected())
		+ ", total=" + ::to_string(steps_total()) + "\n";
}
/*--------------------------------------------------------------------------*/
/* SIM::store: store data in preparation for post processing
 * use storelist as index for waves
 */
void TTT::store_results(double time)
{ untested();
	trace3("TTT::store_results()", tt_iteration_number(), iteration_number(), time );
	assert(_sim->_mode==s_TTT);
	int ii = 0;
	for (PROBELIST::const_iterator
			p=storelist().begin();  p!=storelist().end();  ++p) {
		trace3("TTT::store_results", time, (*p)->label(), (*p)->value());
		_wavep_tt[ii++]->push(time, (*p)->value());
	}
}
/*--------------------------------------------------------------------------*/
void TTT::advance_Time(void)
{
	trace6("TTT::advance_Time", _sim->_tr[0], _sim->_tt[0], _sim->_Time0, _sim->_adp_nodes, _tstop, _sim->_time0);
	::status.tt_advance.start();

	static double last_iter_time;
	// _sim->_time0 = 0.;
	assert(_sim->last_time() == _sim->_time0);
	_sim->keep_voltages();
	assert(_sim->last_time()*0.9999 <= _tstop);
	if (_sim->_Time0 > 0) {
		if (_sim->_dT0 == 0) {
		} else if (_sim->_Time0 > last_iter_time && _accepted_tt ) {
			/* moving forward */
			_sim->_tt_iter++;
//			_sim->_tt_rejects = 0;
			_sim->update_tt_order();

//			trace2("TTT::advance_Time", _sim->_tr1[0], _sim->_tt1[0]);

			notstd::copy_n(_sim->_tr2, _sim->_adp_nodes, _sim->_tr3);
			notstd::copy_n(_sim->_tr1, _sim->_adp_nodes, _sim->_tr2);

			assert(is_number(_sim->_tt[0]) || !_sim->_adp_nodes);

			notstd::copy_n(_sim->_tr, _sim->_adp_nodes, _sim->_tr1);
			notstd::copy_n(_sim->_tt, _sim->_adp_nodes, _sim->_tt1);

			// invalidate ....
			// std::fill_n(_sim->_tr, _sim->_adp_nodes, NAN); // no. needed if ttsteporder==0
			// std::fill_n(_sim->_tt, _sim->_adp_nodes, NAN); // no. needed if ttsteporder==0

			trace3("TTT::advance_Time done", _sim->_tr[0], _sim->_tt[0], _sim->_Time0);
			trace2("TTT::advance_Time done", _sim->_tr1[0], _sim->_tt1[0]);

			// assert(is_number(_sim->_tr1[0]) || !_sim->_adp_nodes); // BUG?
			assert(is_number(_sim->_tt1[0]) || !_sim->_adp_nodes);

			ADP_NODE_LIST::adp_node_list.do_forall( &ADP_NODE::tt_advance );
			/// AFTER tr and tt have been shifted.
			CARD_LIST::card_list.tt_advance();

		} else {

			std::fill_n(_sim->_tr, _sim->_adp_nodes, NAN); // invalidate.
			std::fill_n(_sim->_tt, _sim->_adp_nodes, NAN); // invalidate.

			_sim->restore_voltages();
			_sim->_time0 = 0.;
			CARD_LIST::card_list.tt_regress();
			/* moving backward */
		}
	} else {
		trace0("TTT::advance_Time obsolete call?");
		ADP_NODE_LIST::adp_node_list.do_forall( &ADP_NODE::tt_advance );
		/// AFTER tr and tt have been shifted.
		CARD_LIST::card_list.tt_advance();
	}
	last_iter_time = _sim->_Time0;
	trace0("TTT::advance_Time() done");
	::status.tt_advance.stop();
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

} // namespace
// vim anything.
