/*$Id: s__.h,v 1.9 2009-12-16 17:22:07 felix Exp $ -*- C++ -*-
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
 * base class for simulation methods
 */
//testing=script,complete 2006.07.14
#ifndef S___H
#define S___H
#include "u_opt.h"
#include "c_comand.h"
/*--------------------------------------------------------------------------*/
class CARD;
class CARD_LIST;
class CS;
class PROBELIST;
class COMPONENT;
class WAVE;
/*--------------------------------------------------------------------------*/
class SIM : public CMD {
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
protected:
  enum TRACE { // how much diagnostics to show
    tNONE      = 0,	/* no extended diagnostics			*/
    tUNDER     ,	/* show underlying analysis, important pts only	*/
    tALLTIME   ,	/* show every time step, including hidden 	*/
    tGUESS     ,	/* show guesses (from predictor) */
    tREJECTED  ,	/* show rejected time steps			*/
    tITERATION ,	/* show every iteration, including nonconverged	*/
    tVERBOSE   ,	/* show extended diagnostics			*/
    tMATRIX    = 16,	/* dump matrix.			*/
    tDEBUG     = 32
  };
  enum ALARM { // how much diagnostics to show
    aNONE      = 0,	/* ignore range violations */
    aREPORT    = 1,	/* report range violations */
    aREDIR     = 2,	/* redirect to stderr */
    aABORT     = 4	/* abort simulation */
  };
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  CARD_LIST* _scope;
  OMSTREAM   _out;		/* places to send the results		*/
  bool _dump_matrix;
public:
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
private:
  const std::string long_label()const {unreachable(); return "";}
private:
  virtual void	setup(CS&)	= 0;
  virtual void	sweep()		= 0;
  virtual void	finish()	{}
  virtual bool	is_step_rejected()const {return false;}

  explicit SIM(const SIM&):CMD(),_scope(NULL) {unreachable(); incomplete();}
protected:
  explicit SIM(): CMD(),_scope(NULL) {}
public:
  ~SIM();
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
protected:
  	 void	command_base(CS&);	/* s__init.cc */
	 bool common_options(CS& Cmd);
	 void	reset_timers();	
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
protected:
	 const PROBELIST& alarmlist()const;	/* s__out.cc */
	 const PROBELIST& plotlist()const;
	 const PROBELIST& printlist()const;
	 const PROBELIST& storelist()const;
  virtual void	outdata(double);
  virtual void	head(double,double,const std::string&);
  virtual void	print_results(double);
  virtual void	alarm();
  virtual void	store_results(double);
  virtual void	expect_results(double);
protected: // hack?
  WAVE** _wavep;
public:
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
protected:				/* s__solve.cc */
  bool	solve(OPT::ITL,TRACE);
  bool	solve_with_homotopy(OPT::ITL,TRACE);
  void	advance_time();
protected:
	void	finish_building_evalq();
	void	set_flags();
	void	clear_arrays();
	void	evaluate_models();
	void	set_damp();
	void	load_matrix();
	void	solve_equations(TRACE);
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
#endif
// vim:ts=8:sw=2:noet:
