/*$Id: s_tr_swp.cc 2016/09/22 al $ -*- C++ -*-
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
 * sweep time and simulate.  output results.
 * manage event queue
 */
#include "u_time_pair.h"
#include "u_sim_data.h"
#include "u_status.h"
#include "e_card.h" // debugging...
#include "declare.h"	/* gen */
#include "s_tr.h"
#include "e_adp.h" // hack. see below
using namespace std;
//#define ALT_CQ // alternatively clear queue (experimental)
/*--------------------------------------------------------------------------*/
//	void	TRANSIENT::sweep(void);
//	void	TRANSIENT::first(void);
//	bool	TRANSIENT::next(void);
//	void	TRANSIENT::accept(void);
//	void	TRANSIENT::reject(void);
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
namespace TR {
  static std::string step_cause[] = {
    "impossible",
    "user requested",
    "event queue",
    "command line \"skip\"",
    "convergence failure, reducing (itl4)",
    "slow convergence, holding (itl3)",
    "truncation error",
    "ambiguous event",
    "limit growth",
    "initial guess"
  };
}
/*--------------------------------------------------------------------------*/
void TRANSIENT::sweep()
{
  trace2("TRANSIENT::sweep", _cont, _cont_dc);
  _sim->_phase = p_INIT_DC;
  head(_tstart, _tstop, "Time");
  _sim->_bypass_ok = false;
  _sim->set_inc_mode_bad();
 
  if ( _print_only ) { untested();
    // better to be achieved with a singleton time interval...?
    _sim->_phase = p_RESTORE;
    _sim->restore_voltages();
    CARD_LIST::card_list.tr_restore();
      outdata(_sim->_time0, ofPRINT);
      return;

  } else if (_inside_tt && _cont_dc) {
    _sim->restore_voltages(); // required by some tr_begins.
                              // might be possible to move
                              // node dependency into tr_accept
    _sim->_cont = true; // keep coil from overwriting current
    _sim->_phase = p_RESTORE; // short cut differentiate
    CARD_LIST::card_list.tr_restore();
    CARD_LIST::card_list.do_tr();
    _sim->_cont = false;
    _sim->clear_limit();
    advance_time();
    _sim->reset_iteration_counter(iSTEP);
    CARD_LIST::card_list.do_tr();
  } else if (_cont_dc) {
    // continue from DC point.
    _sim->restore_voltages(); // required by some tr_begins.
                              // might be possible to move
                              // node dependency into tr_accept
    if(_inside_tt){ untested();
    }else if(_edge_detect & edYES){
      _sim->keep_voltages(true);
    }
    _sim->_cont = true; // keep coil from overwriting current
    _sim->_phase = p_RESTORE; // short cut differentiate
    CARD_LIST::card_list.tr_begin();
    _sim->_cont = false;
    _sim->clear_limit();
    advance_time();
    _sim->reset_iteration_counter(iSTEP);
    CARD_LIST::card_list.do_tr();

  }else if (_inside_tt) {
    trace0("TRANSIENT::sweep inside tt");
    assert(  _sim->_mode == s_TTT );

    _sim->restore_voltages();
    _sim->_phase = p_RESTORE;
//     CARD_LIST::card_list.tr_restore(); no. tt_advance takes care of it.
    CARD_LIST::card_list.do_tr();

  }else if (_cont) {
    // use the data from last time
    _sim->_phase = p_RESTORE;
    _sim->restore_voltages();
    CARD_LIST::card_list.tr_restore();
    if(_edge_detect & edYES){
      _sim->keep_voltages(true);
    }else{
    }
  }else{
    _sim->clear_limit();
    CARD_LIST::card_list.tr_begin();
  }
  
  first();
  _sim->_genout = gen();
  
  // assert (_sim->_loadq.empty());
  if ( /*_sim->more_uic_now() */ 0 ) {
    trace0("TRAN: more UIC now ");
    advance_time();
    _sim->zero_voltages();
    CARD_LIST::card_list.do_tr();
    while (!_sim->_late_evalq.empty()) {itested(); //BUG// encapsulation violation
      _sim->_late_evalq.front()->do_tr_last();
      _sim->_late_evalq.pop_front();
    }
//    _converged = true;
    _converged = solve_with_homotopy(OPT::DCBIAS,_trace);


  } else if (_sim->uic_now() || _inside_tt ) {
    trace3("TRANSIENT::sweep uic_now solve", _time1, _sim->_time0, _inside_tt);
    advance_time();
    trace0("TRANSIENT::sweep advanced");
    if (!_inside_tt) {
      _sim->zero_voltages(); // ?
    }
    CARD_LIST::card_list.do_tr();    //evaluate_models
    while (!_sim->_late_evalq.empty()) {untested(); //BUG// encapsulation violation
      _sim->_late_evalq.front()->do_tr_last();
      _sim->_late_evalq.pop_front();
    }
    _converged = true;
  } else if (_cont_dc) {
    // continue from DC point.
    // advance_time();
    assert(_sim->_late_evalq.empty());
    while (!_sim->_late_evalq.empty()) {itested(); //BUG// encapsulation violation
      _sim->_late_evalq.front()->do_tr_last();
      _sim->_late_evalq.pop_front();
    }
    _converged = true;
  } else if (_cont_dc && _sim->analysis_is_static()) { untested();
  } else {
    _converged = solve_with_homotopy(OPT::DCBIAS,_trace);
    if (!_converged) {
      error(bWARNING, "did not converge\n");
    }else{
    }
  }
  trace0("TRANSIENT::sweep review...");
  review(); 
  trace0("TRANSIENT::sweep reviewed");
  _accepted = true;
  accept();
  trace0("TRANSIENT::sweep accepted");

  {
    bool printnow = (_sim->_time0 == _tstart || _trace >= tALLTIME);
    int outflags = ofNONE;
    if (printnow) {
      outflags = ofPRINT | ofSTORE | ofKEEP;
    }else{
      outflags = ofSTORE;
    }
    outdata(_sim->_time0, outflags);
  }

  assert(_tstrobe >=OPT::dtmin ); // == wont work because of CAUSE
                                // we do only increase _time_by_user_request if
                                // CAUSE == user.
                                // the second step is always caused by initial guess...
                                // BUG?
  trace3("TRANSIENT::sweep entering loop", (STEP_CAUSE)step_cause(), _sim->_time0, _sim->last_time());
  _edge_break = false;
  _dt_by_edge1 = -NEVER;
  _dt_by_edge0 = -NEVER;
  while (next()) {
    trace3("TRANSIENT::sweep loop ... ", (STEP_CAUSE)step_cause(), _sim->_time0, _sim->last_time());
    _sim->_bypass_ok = false;
    _sim->_phase = p_TRAN;
    _sim->_genout = gen();
    _converged = solve(OPT::TRHIGH,_trace);

    _accepted = _converged && review();

    if (_accepted) {
      trace1("TRANSIENT::sweep accepted", (STEP_CAUSE)step_cause());
      assert(_converged);
      accept();
      if (step_cause() == scUSER) {
        trace3("up order", _sim->_time0-_sim->_dtmin, _time_by_user_request, _sim->_time0+_sim->_dtmin);
	assert( up_order(_sim->_time0-_sim->_dtmin, _time_by_user_request, _sim->_time0+_sim->_dtmin)
            || _sim->_time0 == _tstop );
	++(_sim->_stepno);
        trace1("TRANSIENT::sweep delivered req", _time_by_user_request);
        if (_time_by_user_request<_tstop){
          _time_by_user_request += _tstrobe;	/* advance user time */
          // _time_by_user_request = min(_time_by_user_request, (double)_tstop+_sim->_dtmin);
        } else {
          _time_by_user_request += _tstrobe;	/* advance user time */
        }
      }else{
      }
      assert(_sim->_time0 <= _time_by_user_request);
    }else{
      trace2("TRANSIENT::sweep NOT accepted", _sim->_time0, _trace);
      reject();
      assert(_time1 < _time_by_user_request);
    }
    {
      bool printnow =
	(_trace >= tREJECTED)
	|| (_accepted && (_trace >= tALLTIME
			  || step_cause() == scUSER
			  || (!_tstrobe.has_hard_value() && _sim->_time0+_sim->_dtmin > _tstart)));
      int outflags = ofNONE;
      if (printnow) {
	outflags = ofPRINT | ofSTORE | ofKEEP;
      }else if (_accepted) {untested();
	// ++::status.hidden_steps;
	outflags = ofSTORE;
      }else{untested();
      }
      outdata(_sim->_time0, outflags);
    }
    
    if (!_converged && OPT::quitconvfail) {untested();
      outdata(_sim->_time0, ofPRINT);
      throw Exception("convergence failure, giving up");
    }else{
    }
  }

//  _sim->_time0 = _sim->last_time(); // hack?
//                                    // this is required because keep_voltage
//                                    // calls set last_time to time0

  for (uint_t ii = _sim->_lu.size(); ii > 0; --ii) {
    _sim->_nstat[ii].set_discont(disNONE);
  }
  assert(!_sim->_nstat[0].discont());
}
/*--------------------------------------------------------------------------*/
void TRANSIENT::set_step_cause(STEP_CAUSE C)
{
  switch (C) {
  case scITER_A:untested();
  case scADT:untested();
  case scUSER:
  case scEVENTQ:
  case scSKIP:
  case scITER_R:
  case scTE:
  case scAMBEVENT:
  case scINITIAL:
    ::status.control = C;
    break;
  case scNO_ADVANCE:untested();
  case scZERO:untested();
  case scSMALL:
  case scREJECT:
    ::status.control += C;
    break;
  case scGROW:
  case scLAST:
        incomplete();
  }
}
/*--------------------------------------------------------------------------*/
int TRANSIENT::step_cause()const
{
  return ::status.control;
}
/*--------------------------------------------------------------------------*/
void TRANSIENT::first()
{
  /* usually, _sim->_time0, time1 == 0, from setup */
  assert(_sim->_time0 == _time1);
  // assert(_sim->_time0 <= _tstart); // oops?
  ::status.review.start();

  //_eq.Clear();					/* empty the queue */
  while (!_sim->_eq.empty()) {untested();
    _sim->_eq.pop();
  }
  _sim->_stepno = 0;

  //_time_by_user_request = _sim->_time0 + _tstrobe;	/* set next user step */
  //set_step_cause(scUSER);

  if (_sim->_time0 < _tstart) {			// skip until _tstart
    set_step_cause(scINITIAL);				// suppressed 
    _time_by_user_request = _tstart;			// set first strobe
  }else{					// no skip
    set_step_cause(scUSER);				// strobe here
    _time_by_user_request = _sim->_time0 + _tstrobe;	// set next strobe
  }

  ::status.hidden_steps = 0;
  ::status.review.stop();
}
/*--------------------------------------------------------------------------*/
#define check_consistency() {						\
    trace4("", __LINE__, newtime, almost_fixed_time, fixed_time);	\
    assert(almost_fixed_time <= fixed_time);				\
    assert(newtime <= fixed_time);					\
    /*assert(newtime == fixed_time || newtime <= fixed_time -_sim->_dtmin);*/	\
    assert(newtime <= almost_fixed_time);				\
    /*assert(newtime == almost_fixed_time || newtime <= almost_fixed_time - _sim->_dtmin);*/ \
    assert(newtime > _time1);						\
    assert(newtime > reftime);						\
    assert(new_dt > 0.);                                                \
    assert(new_dt >= _sim->_dtmin);                                     \
    /*assert(new_dt == 0 || new_dt / _sim->_dtmin > 0.999999);*/        \
   /* assert(newtime <= _time_by_user_request);		                */ \
    /*assert(newtime == _time_by_user_request*/				\
    /*    || newtime < _time_by_user_request - _sim->_dtmin);  */	\
  }
#define check_consistency2() {						\
    assert(newtime > _time1);						\
    assert(new_dt > 0.);						\
    assert(new_dt >= _sim->_dtmin);						\
    assert(newtime <= _time_by_user_request);				\
    /*assert(newtime == _time_by_user_request	*/			\
    /*	   || newtime < _time_by_user_request - _sim->_dtmin);*/		\
  }
/*--------------------------------------------------------------------------*/
/* next: go to next time step
 * Set _sim->_time0 to the next time step, store the old one in time1.
 * Try several methods.  Take the one that gives the shortest step.
 */
bool TRANSIENT::next()
{
  ::status.review.start();

  double old_dt = _sim->_time0 - _time1;
  assert(old_dt >= 0);

  
  double newtime = NEVER;
  double new_dt = NEVER;
  STEP_CAUSE new_control = scNO_ADVANCE;

  if (_sim->_time0 == _time1) {
    // initial step -- could be either t==0 or continue
    // for the first time, just guess
    // make it 100x smaller than expected
    new_dt = std::max(_dtmax/100., _sim->_dtmin);
    newtime = _sim->_time0 + new_dt;
    new_control = scINITIAL;

    if (_time_by_ambiguous_event < newtime) {
      newtime = _time_by_ambiguous_event;
      new_dt = _time_by_ambiguous_event - _time1;
      new_control = scAMBEVENT;
    }else{
    }

    // UGLY duplicate
    if (up_order(newtime-_sim->_dtmin, _time_by_user_request, newtime+_sim->_dtmin)) {
      new_control = scUSER;
    }
  }else if (!_converged) {
    new_dt = old_dt / OPT::trstepshrink;
    newtime = _time_by_iteration_count = _time1 + new_dt;
    new_control = scITER_R;
    //    fprintf(stderr,".");
    //    _trace=tDEBUG;
  }else{
  }
  {
    double reftime;
    if (_accepted) {
      reftime = _sim->_time0;
    }else{
      reftime = _time1;
      trace0("rejected");
    }
    trace2("TRANSIENT::next ", step_cause(), old_dt);
    trace3("TRANSIENT::next ", _time1, _sim->_time0, reftime);
    trace2("TRANSIENT::next ", _time_by_user_request, _sim->_dtmin);

    if (_time_by_user_request < newtime) {
      newtime = _time_by_user_request;
      new_dt = newtime - reftime;
      if (new_dt < _sim->_dtmin) { itested();
	// last step handler?
	new_dt = _sim->_dtmin;
	newtime = reftime + _sim->_dtmin;
      }else{itested();
      }
      new_control = scUSER;
    }else{
    }
    double fixed_time = newtime;
    double almost_fixed_time = newtime;
    trace2("TRANSIENT::next", _time_by_user_request, newtime);
    check_consistency();

    
    // event queue, events that absolutely will happen
    // exact time.  NOT ok to move or omit, even by _sim->_dtmin
    // some action is associated with it.
    if (!_sim->_eq.empty() && _sim->_eq.top() < newtime) {
      trace2("TRANSIENT eventq", newtime, _time1);
      newtime = _sim->_eq.top();
      assert( newtime >= _time1 );
      new_dt = newtime - reftime;
      if (new_dt < _sim->_dtmin) {untested();
	//new_dt = _sim->_dtmin;
	//newtime = reftime + new_dt;
      }else{
      }
      new_control = scEVENTQ;
      fixed_time = newtime;
      almost_fixed_time = newtime;
      trace2("checking", reftime, newtime);
      check_consistency();
    }else if ( !_sim->_eq.empty()  ) {
      trace2("TRANSIENT skipping non empty eq", _time1, _sim->_eq.top() );
    } else {
      trace0("TRANSIENT no events pending");
    }
    // device events that may not happen
    // not sure of exact time.  will be rescheduled if wrong.
    // ok to move by _sim->_dtmin.  time is not that accurate anyway.
    trace1("next ambevt", _time_by_ambiguous_event);
    if (_time_by_ambiguous_event < newtime - _sim->_dtmin) {  
      if (_time_by_ambiguous_event < _time1 + 2*_sim->_dtmin) {untested();
	double mintime = _time1 + 2*_sim->_dtmin;
	if (newtime - _sim->_dtmin < mintime) {untested();
	  newtime = mintime;
	  new_control = scAMBEVENT;
	}else{ untested();
	}
      }else{
	newtime = _time_by_ambiguous_event;
	new_control = scAMBEVENT;
      }
      new_dt = newtime - reftime;
      almost_fixed_time = newtime;
      check_consistency();
    }else{
    }
    
    // device error estimates
    if (_time_by_error_estimate < newtime - _sim->_dtmin) {
      newtime = _time_by_error_estimate;
      new_dt = newtime - reftime;
      new_control = scTE;
      check_consistency();
      trace3("TRANSIENT::next err", newtime, new_dt, _time_by_error_estimate);
    }else{
      trace3("TRANSIENT::next", newtime, new_dt, _time_by_error_estimate);
    }
    
    // skip parameter
    if (new_dt > _dtmax) {
      if (new_dt > _dtmax + _sim->_dtmin) {
	new_control = scSKIP;
      }else{
      }
      new_dt = _dtmax;
      newtime = reftime + new_dt;
      assert(newtime <= _time_by_user_request);
      check_consistency();
    }else{
    }

    // converged but with more iterations than we like
    if ((new_dt > (old_dt + _sim->_dtmin) * OPT::trstephold)
	&& _sim->exceeds_iteration_limit(OPT::TRLOW)) { untested();
      assert(_accepted);
      new_dt = old_dt * OPT::trstephold;
      newtime = reftime + new_dt;
      new_control = scITER_A;
      check_consistency();
    }else{
    }

    // limit growth
    if (_sim->analysis_is_tran_dynamic() && new_dt > old_dt * OPT::trstepgrow) { untested();
      new_dt = old_dt * OPT::trstepgrow;
      newtime = reftime + new_dt;
      new_control = scADT;
      check_consistency();
    }else{
    }

    // quantize
    if (newtime < almost_fixed_time) {
      assert(new_dt >= 0);
      if (newtime < _sim->_time0) {
	assert(reftime == _time1);
	assert(reftime < _sim->_time0); // not moving forward
	// try to pick a step that will end up repeating the rejected step
	// with an integer number of same size steps
	double target_dt = _sim->_time0 - reftime;
	assert(target_dt > new_dt);
	double steps = 1 + floor((target_dt - _sim->_dtmin) / new_dt);
	assert(steps > 0);
	new_dt = target_dt / steps;
	newtime = reftime + new_dt;
	check_consistency();
      }else if (newtime > reftime + old_dt*.8
	  && newtime < reftime + old_dt*1.5
	  && reftime + old_dt <= almost_fixed_time) {
	// new_dt is close enough to old_dt.
	// use old_dt, to avoid a step change.
        old_dt = max(_sim->_dtmin,old_dt);  // eliminate numerical noise to pass assertion
	assert(reftime == _sim->_time0); // moving forward
	assert(reftime > _time1);
	new_dt = old_dt;
	newtime = reftime + new_dt;
	if (newtime > almost_fixed_time) { untested();
	  new_control = scAMBEVENT;
	  newtime = almost_fixed_time;
	  new_dt = newtime - reftime;
	}else{
	}
        trace4("TRANSIENT::next quantized", new_dt, _sim->_dtmin, _sim->_dtmin-new_dt, new_control );
	check_consistency();
      }else{
	// There will be a step change.
	// Try to choose one that we will keep for a while.
	// Choose new_dt to be in integer fraction of target_dt.
	assert(reftime == _sim->_time0); // moving forward
	//assert(reftime > _time1); // _time1==_time0 on restart, ok
	double target_dt = fixed_time - reftime;
	assert(target_dt >= new_dt);
	double steps = 1 + floor((target_dt - _sim->_dtmin) / new_dt);
	assert(steps > 0);
	new_dt = target_dt / steps;
        trace3("TRANSIENT::next step change", reftime, new_dt, reftime+new_dt);
	newtime = reftime + new_dt;
	check_consistency();
      }
    }else{
      assert(newtime == almost_fixed_time);
    }

    // trap time step too small
    if (!_accepted && new_dt < _sim->_dtmin) { untested();
      new_dt = _sim->_dtmin;
      newtime = reftime + new_dt;
      new_control = scSMALL;
      check_consistency();
    }else{
    }

    // if all that makes it close to user_requested, make it official
    if (up_order(newtime-_sim->_dtmin, _time_by_user_request, newtime+_sim->_dtmin)) {
      //newtime = _time_by_user_request;
      //new_dt = newtime - reftime;
      new_control = scUSER;
    }
  } // end of else converged 

  set_step_cause(new_control);

  trace4("TRANSIENT::next got it i think", newtime, new_control, newtime-_sim->_time0, _time_by_user_request);
  
  /* check to be sure */
  if (newtime < _time1 + _sim->_dtmin) {
    /* It's really bad. */
    /* Reject the most recent step, back up as much as possible, */
    /* and creep along */
    trace3("TRANSIENT::next ", newtime, _time1, _sim->_dtmin );
    assert(!_accepted);
    assert(step_cause() < scREJECT);
    assert(step_cause() >= 0);
    error(bDANGER,"non-recoverable: " + TR::step_cause[step_cause()] + "\n");
    error(bDANGER, "newtime=%e  rejectedtime=%e  oldtime=%e  using=%e\n",
	  newtime, _sim->_time0, _time1, _time1 + _sim->_dtmin);
    newtime = _time1 + _sim->_dtmin;
    assert(newtime <= _time_by_user_request);
    set_step_cause(scSMALL);
    //check_consistency2();
    throw Exception("tried everything, still doesn't work, giving up");
  }else if (newtime < _sim->_time0) {
    /* Reject the most recent step. */
    /* We have faith that it will work with a smaller time step. */
    assert(!_accepted);
    assert(newtime >= _time1 + _sim->_dtmin);
    error(bLOG, "backwards time step\n");
    error(bLOG, "newtime=%e  rejectedtime=%e  oldtime=%e\n", newtime, _sim->_time0, _time1);
    set_step_cause(scREJECT);
    _sim->mark_inc_mode_bad();
    check_consistency2();
#if 0
  }else if (_sim->_time0 >= _time_by_user_request) {
    trace3("TRANSIENT::next: already there", _sim->_time0, _time1, newtime );
    _time1 = _sim->_time0;
    new_dt = 0;
    check_consistency2();
#endif
  }else if (newtime < _sim->_time0 + _sim->_dtmin) { untested();
    /* Another evaluation at the same time. */
    /* Keep the most recent step, but creep along. */
    assert(newtime > _sim->_time0 - _sim->_dtmin);
    error(bDANGER, "zero time step\n");
    error(bDANGER, "newtime=%e  rejectedtime=%e delta=%e time1=%e requested=%e dtmin=%e, control=%d\n",
        newtime, _sim->_time0, newtime-_sim->_time0, _time1, _time_by_user_request, _sim->_dtmin, step_cause());
    if (_accepted) {untested();
      _time1 = _sim->_time0;
    }else{untested();
      assert(_converged);
    }
    trace3( "TRANSIENT::next:", newtime, _time1, new_dt);
    check_consistency2();
    newtime = _sim->_time0 + _sim->_dtmin;
    if (newtime > _time_by_user_request) { untested();
      newtime = _time_by_user_request;
      set_step_cause(scUSER);
    }else{
    }
    assert (newtime<=_tstop);
    set_step_cause(scZERO);
    check_consistency2();
  }else{
    assert(_accepted);
    assert(newtime >= _sim->_time0 + _sim->_dtmin);
    /* All is OK.  Moving on. */
    /* Keep value of newtime */
    _time1 = _sim->_time0;
    check_consistency2();
  }
  _sim->_time0 = newtime;
  assert(_sim->_time0 <= _time_by_user_request);
  
  /* advance event queue (maybe) */
  /* We already looked at it.  Dump what's on top if we took it. */

#ifndef ALT_CQ

  while (!_sim->_eq.empty() && _sim->_eq.top() <= _sim->_time0) {
    trace1("eq", _sim->_eq.top());
    _sim->_eq.pop();
  }
  while (!_sim->_eq.empty() && _sim->_eq.top() < _sim->_time0 + _sim->_dtmin) {
    trace1("eq-extra", _sim->_eq.top());
    _sim->_eq.pop();
  }
#endif
  // event queue cleanup moved....
  //BUG// what if it is later rejected?  It's lost!
  // -> why not move to tr_advance? tr_accept? hmmm
  //

  if(_time1 < _tstop - _sim->_dtmin
      && _sim->_time0 > _tstop + _sim->_dtmin ) {
    _sim->_time0 = _tstop;
    check_consistency2();
  } else {
    check_consistency2();
  }
  assert(_sim->_time0 <= _time_by_user_request);

  check_consistency2();
  ++steps_total_;
  ::status.review.stop();
  bool ret= _sim->_time0 <= _tstop; // throw away last step if it helps.  + _sim->_dtmin;
  ret= _sim->_time0 < _tstop + _sim->_dtmin; // this once worked
  if(_accepted && _edge_break){
    trace3("accepted edge", _sim->_time0, _time1, _sim->_time0 - _time1);
    return false;
  }else if(_edge_break){
    trace5("unaccepted edge", _sim->_time0, _sim->v0dist_min(), _time1, _time_by_ambiguous_event, STEP_CAUSE(step_cause()));
  }

  return (ret);
}
/*--------------------------------------------------------------------------*/
bool TRANSIENT::review()
{
  ::status.review.start();
  _sim->count_iterations(iTOTAL);

  TIME_PAIR time_by = CARD_LIST::card_list.tr_review();
  _time_by_error_estimate = time_by._error_estimate;
  _time_by_ambiguous_event = time_by._event;
  _edge_break = false;

  double dt0 = _sim->_time0 - _time1;
  trace3("review times", _sim->_time0, dt0, _time1);

  // fixme: move to a component (?)... maybe later.
  if (_sim->_time0 == 0.){
  }else if (dt0 == 0. ){
  }else if (_sim->_time0 < _tstop * .2){
    // incomplete();
  }else if (_edge_detect & edEVT){
    double dist = 0;
    double dt_factor = _sim->v0dist_min(&dist);
    double dt_by_edge = dt_factor * dt0;
    // if(dt_by_edge < 0. && dtedge1 > 0){ untested();
    //   dt_by_edge = (dt_by_edge + dtedge1)*.5;
    // }

    double time_by_edge = dt_by_edge + _time1;
    trace4("edge", time_by_edge, _sim->_time0, dt_by_edge, dist);

    if(dist>10){
      time_by_edge = NEVER;
    }else if(fabs(dt_by_edge) > 10 * fabs(_dt_by_edge1)){
      time_by_edge = NEVER;
    }else if(_dt_by_edge1 < 0.){
      time_by_edge = NEVER;
    }else if(_dt_by_edge1 < dt0){
      time_by_edge = NEVER;
      // ?!
    }else if(dt_by_edge < 0. && _dt_by_edge1 > 0. ){
      trace3("negative", dt_by_edge, _sim->_time0, _dt_by_edge1);
      time_by_edge = _time1 + dt0 * .3;
      _edge_break = _edge_detect & edBREAK;
    }else if(dt_by_edge < 0. ){ untested();
      time_by_edge = _sim->_dt0 * .5 + _time1;
    }else if(dt_by_edge < dt0 && _dt_by_edge1 > dt0){
      trace4("got edge", dt_by_edge, _sim->_time0, _dt_by_edge1, dt_factor);
      _edge_break = _edge_detect & edBREAK; // edge between time0 and time1.
                                      // if it gets accepted, we are done
//      time_by_edge = dt_by_edge * .5 + _time1;
    }else if(dt_by_edge > dt0){
      // everything fine. edge is far ahead.
    }else if(dt_by_edge + _sim->_dtmin > dt0 && _dt_by_edge1 > dt0 ){ untested();
      trace4("close to edge", dt0, dt_by_edge, _sim->_time0, _dt_by_edge1);
      _edge_break = _edge_detect & edBREAK; // edge close enogh
    }else if(dt_by_edge < dt0){ untested();
      // past edge
      time_by_edge = NEVER;
    }

    if(!(_edge_detect & edEVT)){ untested();
    }else if(time_by_edge < _time_by_ambiguous_event){
      trace5("setting evt", time_by_edge, _sim->_time0, dt_by_edge, _time_by_ambiguous_event, dt_by_edge / _sim->_dtmin);
      _time_by_ambiguous_event = time_by_edge;
    }else{
    }

    _dt_by_edge0 = dt_by_edge;
  }

  // limit minimum time step
  // 2*_sim->_dtmin because _time[1] + _sim->_dtmin might be == _time[0].
  if (_time_by_ambiguous_event < _time1 + 2*_sim->_dtmin) {
    _time_by_ambiguous_event = _time1 + 2*_sim->_dtmin;
  }else{
  }
  // force advance when time too close to previous
  if (std::abs(_time_by_ambiguous_event - _sim->_time0) < 2*_sim->_dtmin) {
    _time_by_ambiguous_event = _sim->_time0 + 2*_sim->_dtmin;
  }else{
  }

  if (time_by._error_estimate < _time1 + 2*_sim->_dtmin) {itested();
    _time_by_error_estimate = _time1 + 2*_sim->_dtmin;
  }else{
    _time_by_error_estimate = time_by._error_estimate;
  }
  if (std::abs(_time_by_error_estimate - _sim->_time0) < 1.1*_sim->_dtmin) {
    _time_by_error_estimate = _sim->_time0 + 1.1*_sim->_dtmin;
  }else{
  }

  ::status.review.stop();

  return (_time_by_error_estimate > _sim->_time0  &&  _time_by_ambiguous_event > _sim->_time0);
}
/*--------------------------------------------------------------------------*/
void TRANSIENT::accept()
{
  ::status.accept.start();
  for (unsigned ii = _sim->_lu.size(); ii >= 1; --ii) {
      _sim->_nstat[ii].set_discont(disNONE);
  }

  _dt_by_edge1 = _dt_by_edge0;

#ifdef ALT_CQ_PRE
  while (!_sim->_eq.empty() && _sim->_eq.top() < _sim->_time0 + _sim->_dtmin) {itested();
    trace1("TRANSIENT::accept eq-pop-extra", _sim->_eq.top());
    _sim->_eq.pop();
  }
#endif

  trace2("TRANSIENT::accept dt0", _sim->_dt0, _sim->_time0 - _time1);
  _sim->_dt0 = _sim->_time0 - _time1;
  if(_sim->_dt0 <=0) assert (_sim->_stepno == 0);
  _sim->set_limit();
  if (OPT::traceload) { // traceload == "use queue"
    while (!_sim->_acceptq.empty()) {
      trace1("TRANSIENT::accept", _sim->_acceptq.back()->long_label());
      _sim->_acceptq.back()->tr_accept();
      _sim->_acceptq.pop_back();
    }
  }else{itested();
    _sim->_acceptq.clear();
    CARD_LIST::card_list.tr_accept();
  }
//  tmp hack don't know yet how to fix (always_q_for_accept?)
  ADP_LIST::adp_list.tr_accept();
  
  ++steps_accepted_;
  if( _sim->analysis_is_tt() || OPT::trage ) {
    trace0( "TRANSIENT::accept: done stressing cardlist");
    if ( OPT::trage ) {
      incomplete();
      CARD_LIST::card_list.do_forall( &CARD::do_tt );
    }
  }
# ifdef ALT_CQ
  while (!_sim->_eq.empty() && _sim->_eq.top() <= _sim->_time0) {
    trace1("eq", _sim->_eq.top());
    _sim->_eq.pop();
  }
# endif
  _sim->_nstat[0].set_discont(disNONE);
  ::status.accept.stop();
}
/*--------------------------------------------------------------------------*/
void TRANSIENT::reject()
{
  ::status.accept.start();
  _sim->_acceptq.clear();
  ++steps_rejected_;
  ::status.accept.stop();
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
