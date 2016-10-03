/*                        -*- C++ -*-
 * Copyright (C) 2008 Albert Davis
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
 */
#include "u_parameter.h"
#include "m_wave.h"
#include "u_function.h"
/*--------------------------------------------------------------------------*/
namespace {
/*--------------------------------------------------------------------------*/
class MIN : public WAVE_FUNCTION {
  PARAMETER<double> before;
  PARAMETER<double> after;
  bool last;
  bool arg;
  bool derivative;

public:
  MIN(): 
    WAVE_FUNCTION(),
    before(BIGBIG),
    after(-BIGBIG),
    last(false),
    arg(false),
    derivative(false)
  {}
  virtual FUNCTION_BASE* clone()const { return new MIN(*this);}
  string label()const{return "min";}
  void expand(CS& Cmd, const CARD_LIST* Scope)
  {
    before = BIGBIG;
    after = -BIGBIG;
    last = false;
    arg = false;
    derivative = false;
    unsigned here = Cmd.cursor();
    Cmd >> probe_name;
    _w = find_wave(probe_name);

    if (!_w) {
      Cmd.reset(here);
    }else{
    }

    here = Cmd.cursor();
    do {
      ONE_OF
	|| Get(Cmd, "probe",  &probe_name)
	|| Get(Cmd, "before", &before)
	|| Get(Cmd, "after",  &after)
	|| Get(Cmd, "end",    &before)
	|| Get(Cmd, "begin",  &after)
	|| Set(Cmd, "deriv{ative}", &derivative, true)
	|| Set(Cmd, "arg",    &arg, true)
	|| Set(Cmd, "last",   &last, true)
	|| Set(Cmd, "first",  &last, false)
	;
    }while (Cmd.more() && !Cmd.stuck(&here));

    if (!_w) {
      _w = find_wave(probe_name);
    }else{
    }
    before.e_val(BIGBIG, Scope);
    after.e_val(-BIGBIG, Scope);
  }

  double wave_eval()const
  {
    trace1("MIN::wave_eval", probe_name);
    if (_w) {

      double time = (last) ? -BIGBIG : BIGBIG;
      double m = BIGBIG;
      WAVE::const_iterator begin = lower_bound(_w->begin(), _w->end(), DPAIR(after, -BIGBIG));
      WAVE::const_iterator end   = upper_bound(_w->begin(), _w->end(), DPAIR(before, BIGBIG));
      if (begin == end) return(NAN);

      double prev = 0;
      double t1 = begin->first;

      for (WAVE::const_iterator i = begin; i < end; ++i) {
        double val;
        if(derivative){ // very simple difference quotient.
          double dt = i->first - t1;
          t1 = i->first;
          if(!dt) {
            prev = i->second;
            continue;
          }
          val = (i->second - prev)/(dt);
          prev = i->second;
        }else{
          val = i->second;
        }

	double at = i->first;
        trace3("MIN::wave_eval", probe_name, at, val);
	if (val < m || (last && (val == m))) {
	  time = at;
	  m = val;
	}else{
	}
      }
      return ((arg) ? (time) : (m));
    }else{
      throw Exception_No_Match(probe_name);
    }
  }
} p2;
DISPATCHER<FUNCTION_BASE>::INSTALL d2(&measure_dispatcher, "min", &p2);
/*--------------------------------------------------------------------------*/
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
