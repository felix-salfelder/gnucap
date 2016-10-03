/*$Id: measure_at.cc,v 1.2 2009-12-13 17:55:01 felix Exp $ -*- C++ -*-
 * vim:ts=8:sw=2:et:
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
class MEASURE : public WAVE_FUNCTION {
  PARAMETER<double> x;
  bool derivative;
public:
  MEASURE() :
    WAVE_FUNCTION(),
    derivative(false)
  {}
  virtual FUNCTION_BASE* clone()const { return new MEASURE(*this);}
  string label() const{return "at";}
  void expand(CS& Cmd, const CARD_LIST* Scope){
    unsigned here = Cmd.cursor();
    Cmd >> probe_name;
    _w = find_wave(probe_name);

    if (!_w) {
      Cmd.reset(here);
    }else{
    }
    x = BIGBIG;
    derivative=false;

    here = Cmd.cursor();
    do {
      ONE_OF
        || Get(Cmd, "probe",	    &probe_name)
        || Get(Cmd, "x",	    &x)
        || Get(Cmd, "at",	    &x)
        || Get(Cmd, "deriv{ative}", &derivative)
        ;
    }while (Cmd.more() && !Cmd.stuck(&here));

    if (!_w) {
      _w = find_wave(probe_name);
    }else{
    }
    x.e_val(BIGBIG, Scope);
  }
/*--------------------------------------------------------------------------*/
  fun_t wave_eval()const
  {
    if (_w) {
      return to_fun_t((derivative) ? (_w->v_out(x).f1) : (_w->v_out(x).f0));
    }else{
      throw Exception_No_Match(probe_name);
    }
  }
} p1;
DISPATCHER<FUNCTION_BASE>::INSTALL d1(&measure_dispatcher, "at", &p1);
/*--------------------------------------------------------------------------*/
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
