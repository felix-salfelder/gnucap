/*                              -*- C++ -*-
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
 * Transient analysis
 */
//testing=script,complete 2006.07.14
#ifndef S_TR_H
#define S_TR_H
#include "u_parameter.h"
#include "s__.h"
/*--------------------------------------------------------------------------*/
class TRANSIENT : public SIM {
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
		
protected:
  enum STEP_CAUSE {
    scUSER      =  1,	/* user requested				*/
    scEVENTQ    =  2,	/* an "event" from the queue			*/
    scSKIP      =  3,	/* effect of "skip" parameter			*/
    scITER_R    =  4,	/* iter count exceeds itl4 (reducing)		*/
    scITER_A    =  5,	/* iter count exceeds itl3 (holding)		*/
    scTE        =  6,	/* truncation error, or device stuff		*/
    scAMBEVENT	 =  7,	/* ambiguous event				*/
    scADT       =  8,	/* by iter count limited by max(rdt, 2*adt)	*/
    scINITIAL   =  9,	/* initial guess				*/
    scREJECT    = 10,	/* rejected previous time step			*/
    scZERO      = 11,	/* fixed zero time step				*/
    scSMALL     = 12,	/* time step too small				*/
    scNO_ADVANCE= 13,	/* after all that it still didn't advance	*/
    scGROW      = 14,	/* requested step too large, truncating */
    scLAST      = 15 	/* last step */
  };
  enum EDGE_DETECT {
    edNONE     =  0,
    edYES      =  1, // do edge computation
    edEVT      =  2, // set evt
    edBREAK    =  4	// stop at edge
  };
  template<typename T> friend T& operator<<( T& o, const STEP_CAUSE &c);
public:
  static const std::string STEP_CAUSE_label[];
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
public:
  explicit TRANSIENT():
    SIM(),
    _skip_in(1),
    _dtmax(0.),
    _cold(false),
    _cont(false),
    _cont_dc(false),
    _inside_tt(false),
    _trace(tNONE),
    _time_by_iteration_count(0.),
    _time_by_user_request(0.),
    _time_by_error_estimate(0.),
    _time_by_ambiguous_event(0.),
    _converged(false),
    _accepted(false),
	 _print_only(false)
  {
	  trace4( "TRANSIENT()",  steps_total_out_, _tstep , _tstop, _tstart );
  }
  ~TRANSIENT() {}
public:
  void	do_it(CS&, CARD_LIST* scope);
  std::string status()const;
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
private:		// s_tr_rev.cc
  bool	review();
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
protected:		// s_tr_set.cc
  void	setup(CS&);
protected:
  void	options(CS&);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
protected:		// s_tr_swp.cc
  void	sweep();
private:
  void	set_step_cause(STEP_CAUSE);
public:
  int stepno(){return _sim->_stepno;}		// count of visible (saved) steps
  int	step_cause()const;
  void	first();
  bool	next();
  void	accept();
  void	reject();
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
private:
  bool	is_step_rejected()const {return (step_cause() > scREJECT);}
  explicit TRANSIENT(const TRANSIENT&): SIM(),_skip_in(1) {unreachable(); incomplete();}
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
protected:
  PARAMETER<double> _tstart;	// sweep start time
  PARAMETER<double> _tstop;	// sweep stop time
  PARAMETER<double> _tstrobe;	// printed step size
  PARAMETER<double> _dtratio_in;// ratio of max/min dt
  PARAMETER<double> _dtmin_in;	// min internal step size
  PARAMETER<double> _dtmax_in;	// max internal step size (user)
  PARAMETER<int>    _skip_in;	// fixed step size: internal steps per external
  double _time1;		/* time at previous time step */
  double _dtmax;	// max internal step size (step / _skip)
  bool _cold;		// flag: start time=0, all voltages=0
  bool _cont;		// flag: continue from previous run
  bool _cont_dc;	// flag: skip initial dc.
//  bool _cont_tt;	// continue with adp_node values. (in TTT)
  bool _inside_tt;// hack: tell tr_swp to not reinit circuit
  //int _stepno;		// count of visible (saved) steps
  unsigned _edge_detect;
  bool _edge_break;
protected:
  TRACE _trace;		// enum: show extended diagnostics
private:
  double _time_by_iteration_count;
  double _time_by_user_request;
  double _time_by_error_estimate;
  double _dt_by_edge0; // hack.
  double _dt_by_edge1; // hack.
protected:
  double _time_by_ambiguous_event;
protected:
  bool _converged;
  bool _accepted;
  bool _print_only;
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
private:
  static int steps_accepted_;
  static int steps_rejected_;
  static int steps_total_;
protected:
  static uint_t steps_total_out_;
public:
  static int steps_accepted() {return steps_accepted_;}
  static int steps_rejected() {return steps_rejected_;}
  static int steps_total()    {return steps_total_;}
  static uint_t steps_total_out()    {return steps_total_out_;}
};

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 

/*--------------------------------------------------------------------------*/
template <class T>
inline T& operator<<( T& o, const TRANSIENT::STEP_CAUSE &c)
{
	if(c > TRANSIENT::scLAST){
		// FIXME: sometimes multiple causes sum up...
		return o << ::to_string(int(c));
	}
	return o << TRANSIENT::STEP_CAUSE_label[c];
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
#endif
// vim:ts=8:sw=2:noet:
