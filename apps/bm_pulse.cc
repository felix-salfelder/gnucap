/*$Id: bm_pulse.cc,v 1.3 2009-12-13 17:55:01 felix Exp $ -*- C++ -*-
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
 * SPICE compatible PULSE
 */
//testing=script 2005.10.06
#include "e_elemnt.h"
#include "u_lang.h"
#include "bm.h"
#include "io_trace.h"
/*--------------------------------------------------------------------------*/
namespace {
/*--------------------------------------------------------------------------*/
using std::max;
using std::min;
/*--------------------------------------------------------------------------*/
const double _default_iv    (NOT_INPUT);
const double _default_pv    (NOT_INPUT);
const double _default_delay (0);
const double _default_rise  (0);
const double _default_fall  (0);
const double _default_width (BIGBIG);
const double _default_period(BIGBIG);
const double _default_t2    (NOT_INPUT);
const double _default_freq  (NOT_INPUT);
const double _default_duty  (NOT_INPUT); // duty cycle (whatever that means)
const double _default_area  (NOT_INPUT); // area
// const double _default_pwr  (0.5); // l2
const double _default_phase (0);
/*--------------------------------------------------------------------------*/
class EVAL_BM_PULSE : public EVAL_BM_ACTION_BASE {
private:
  PARAMETER<double> _iv_in;
  PARAMETER<double> _pv_in;
  PARAMETER<double> _delay;
  PARAMETER<double> _rise_in;
  PARAMETER<double> _fall_in;
  PARAMETER<double> _width;
  PARAMETER<double> _period_in;
  PARAMETER<double> _end;
  PARAMETER<double> _t2;
  PARAMETER<double> _freq;
  PARAMETER<double> _duty;
  PARAMETER<double> _area;
  PARAMETER<double> _phase;
  static std::map<string, PARA_BASE EVAL_BM_PULSE::*> _param_dict;

  explicit	EVAL_BM_PULSE(const EVAL_BM_PULSE& p);
public:
  explicit      EVAL_BM_PULSE(int c=0);
		~EVAL_BM_PULSE()	{}
  int param_count()const { untested();
    return 7 + EVAL_BM_ACTION_BASE::param_count();}
  string param_name(int i)const;
  string param_name(int i,int)const{return param_name(i);}
  string param_value(int)const;
  bool param_is_printable(int i)const;
private: // override vitrual
  bool		operator==(const COMMON_COMPONENT&)const;
  COMMON_COMPONENT* clone()const	{return new EVAL_BM_PULSE(*this);}
  void		print_common_obsolete_callback(OMSTREAM&, LANGUAGE*)const;

  void		precalc_last(const CARD_LIST*);
  void		tr_eval(ELEMENT*)const;
  TIME_PAIR	tr_review(COMPONENT*)const;
  std::string	name()const		{return "pulse";}
  bool		ac_too()const		{return false;}
  bool		parse_numlist(CS&);
  void      set_param_by_name(string Name, string Value);
  bool      parse_params_obsolete_callback(CS& cmd);
private:
  double _pv;
  double _iv;
  double _rise;
  double _fall;
  double _period;
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
EVAL_BM_PULSE::EVAL_BM_PULSE(int c)
  :EVAL_BM_ACTION_BASE(c),
   _iv_in(_default_iv),
   _pv_in(_default_pv),
   _delay(_default_delay),
   _rise_in(_default_rise),
   _fall_in(_default_fall),
   _width(_default_width),
   _period_in(_default_period),
   _end(NOT_VALID),
   _t2(_default_t2),
   _freq(_default_freq),
   _duty(_default_duty),
   _area(_default_area),
   _phase(_default_phase)
{
}
/*--------------------------------------------------------------------------*/
EVAL_BM_PULSE::EVAL_BM_PULSE(const EVAL_BM_PULSE& p)
  :EVAL_BM_ACTION_BASE(p),
   _iv_in(p._iv_in),
   _pv_in(p._pv_in),
   _delay(p._delay),
   _rise_in(p._rise_in),
   _fall_in(p._fall_in),
   _width(p._width),
   _period_in(p._period_in),
   _end(NOT_VALID),
   _t2(p._t2),
   _freq(p._freq),
   _duty(p._duty),
   _area(p._area),
   _phase(p._phase)
{
}
/*--------------------------------------------------------------------------*/
string EVAL_BM_PULSE::param_name(int i)const
{ untested();
  switch (param_count() - 1 - i) { untested();
  case 0:  return "iv";
  case 1:  return "pv";
  case 2:  return "delay";
  case 3:  return "rise";
  case 4:  return "fall";
  case 5:  return "width";
  case 6:  return "period";
  default: return EVAL_BM_ACTION_BASE::param_name(i);
  }
}
/*--------------------------------------------------------------------------*/
string EVAL_BM_PULSE::param_value(int i)const
{
  switch (param_count() - 1 - i) { untested();
  case 0:  return _iv_in.string();
  case 1:  return _pv_in.string();
  case 2:  return _delay.string();
  case 3:  return _rise_in.string();
  case 4:  return _fall_in.string();
  case 5:  return _width.string();
  case 6:  return _period_in.string();
  default: return EVAL_BM_ACTION_BASE::param_name(i);
  }
}
/*--------------------------------------------------------------------------*/
bool EVAL_BM_PULSE::param_is_printable(int i)const
{ untested();
  switch (param_count() - 1 - i) { untested();
  case 0:  return true;
  case 1:  return true;
  case 2:  return _delay.has_hard_value();
  case 3:  return true;
  case 4:  return true;
  case 5:  return true;
  case 6:  return _period_in.has_hard_value();
  default: return EVAL_BM_ACTION_BASE::param_is_printable(i);
  }
}
/*--------------------------------------------------------------------------*/
bool EVAL_BM_PULSE::operator==(const COMMON_COMPONENT& x)const
{
  const EVAL_BM_PULSE* p = dynamic_cast<const EVAL_BM_PULSE*>(&x);
  bool rv = p
    && _iv_in == p->_iv_in
    && _pv_in == p->_pv_in
    && _delay == p->_delay
    && _rise_in == p->_rise_in
    && _fall_in == p->_fall_in
    && _width == p->_width
    && _period_in == p->_period_in
    && _t2 == p->_t2
    && _freq == p->_freq
    && _duty == p->_duty
    && _area == p->_area
    && _phase == p->_phase
    && EVAL_BM_ACTION_BASE::operator==(x);
  return rv;
}
/*--------------------------------------------------------------------------*/
void EVAL_BM_PULSE::print_common_obsolete_callback(OMSTREAM& o, LANGUAGE* lang)const
{ incomplete();
  assert(lang);
  o << name();
  print_pair(o, lang, "iv", _iv_in);
  print_pair(o, lang, "pv", _pv_in);
  print_pair(o, lang, "delay", _delay);
  print_pair(o, lang, "rise", _rise_in);
  print_pair(o, lang, "fall", _fall_in);
  print_pair(o, lang, "width", _width);
  print_pair(o, lang, "period", _period_in);
  if (_t2.has_hard_value()) print_pair(o, lang, "t2", _t2);
  if (_freq.has_hard_value()) print_pair(o, lang, "freq", _freq);
  if (_duty.has_hard_value()) print_pair(o, lang, "duty", _duty);
  if (_area.has_hard_value()) print_pair(o, lang, "area", _area);
  if (_phase.has_hard_value()) print_pair(o, lang, "phase", _phase);
  EVAL_BM_ACTION_BASE::print_common_obsolete_callback(o, lang);
}
/*--------------------------------------------------------------------------*/
void EVAL_BM_PULSE::precalc_last(const CARD_LIST* Scope)
{
  assert(Scope);
  EVAL_BM_ACTION_BASE::precalc_last(Scope);
  _iv_in.e_val(_default_iv, Scope);
  _pv_in.e_val(_default_pv, Scope);
  _delay.e_val(_default_delay, Scope);
  _rise_in.e_val(_default_rise, Scope);
  _fall_in.e_val(_default_fall, Scope);
  _width.e_val(_default_width, Scope);
  _period_in.e_val(_default_period, Scope);
  _t2.e_val(_default_t2, Scope);
  _duty.e_val(_default_duty, Scope);
  _area.e_val(_default_area, Scope);
  _phase.e_val(_default_phase, Scope);
  _freq.e_val(_default_freq, Scope);

  _pv = _pv_in;
  _iv = _iv_in;
  _rise = _rise_in;
  _fall = _fall_in;
  _period = _period_in;

  trace6("pulse precalc_last", _duty, _width, _freq, CKT_BASE::_sim->_freq, _rise_in, _period_in);
  trace2("pulse precalc_last", _pv, _iv);
  if (_freq.has_hard_value()) {
    _width = NOT_INPUT; //hack
    _period = 1./_freq;
    if (_duty.has_good_value()) {

      _delay = - _phase * _period / M_TWO_PI; //  + _rise;
        //(_period*(1-_duty)-_rise)/2.0;

      if (_width.has_good_value()) {
        incomplete(); // probably worth an exception
      } else {
	_width = _duty*_period - (_rise+_fall)/2.0;
      }
    } else { untested();
    }
  }else if (_duty.has_hard_value() && _area.has_hard_value()) { untested();
    assert(0); // for now
    // throw Exception ... bla
  }else if (_duty.has_hard_value() || _area.has_hard_value()) {
    double eff = _duty;
    if (_area.has_hard_value()) {
      eff = _area;
    }

    assert(_rise_in); // for now
    assert(_fall_in); // for now
    _width = NOT_INPUT; //hack
    double freq = max(1.,CKT_BASE::_sim->_freq); // 1 in .dc
    _period = 1./freq;
    _delay = - _phase * _period / M_TWO_PI; //  + _rise;


    double factor = eff*_period/(_rise_in+_fall_in)*2.;
    double cofactor = (1-eff)*_period/(_rise_in+_fall_in)*2.;

    if (_area.has_hard_value()) {
      factor = sqrt(factor);
      cofactor = sqrt(cofactor);
      eff = _area;
    }

    _iv = _pv_in + (_iv_in - _pv_in) * min(1.,cofactor);
    _pv = _iv_in + (_pv_in - _iv_in) * min(1.,factor);

    _fall = _fall_in * min(1., min(cofactor, factor));
    _rise = _rise_in * min(1., min(cofactor, factor));

    if (eff*_period < (_rise+_fall)/2.) {
      _width = 0.;
    }else if ((1-eff)*_period < (_rise+_fall)/2.) {
      _width = _period - (_rise+_fall);
    }else{
      _width = eff*_period - (_rise+_fall)/2.;
    }

    trace6("pulse duty", _duty, _width, freq, CKT_BASE::_sim->_freq, _iv, _pv);
  }else{
  }

  if (_period_in == 0.) {
    trace2("pulse precalc_last zero", _period, _period_in);
    _period = _default_period;
  }else{
    trace2("pulse precalc_last nz", _period, _period_in);
  }

  if (_t2.has_good_value()){ incomplete();
    _width = _t2 - _delay;
  }
}
/*--------------------------------------------------------------------------*/
void EVAL_BM_PULSE::tr_eval(ELEMENT* d)const
{
  double eps = d->_sim->_dtmin * .01;
  double time = d->_sim->_time0;
  if (0 < _period && _period < BIGBIG) {
    //time = fmod(time,_period);
    if (time > _delay) {
      time = fmod(time - _delay, _period) + _delay;
    }else{
    }
  }else{
  }
  double ev = 0; // effective value
  if (time >= eps + _delay+_rise+_width+_fall) {
    /* past pulse	*/
    ev = _iv;
    if (_fall==0.) {
      d->_discont |= disFIRST;
    }
  }else if (time >= eps + _delay+_rise+_width) {
    /* falling 	*/
    double interp = (time - (_delay+_rise+_width)) / _fall;
	 assert(_pv_in != NOT_INPUT);
    ev = _pv + interp * (_iv - _pv);
  }else if (time >= eps + _delay + _rise) {
    /* pulse val 	*/
	 assert(_pv != NOT_INPUT);
    ev = _pv;
    if (_rise==0.) {
      d->_discont |= disFIRST;
    }
  }else if (time >= eps + _delay) {
    /* rising 	*/
    double interp = (time - _delay) / _rise;
    ev = _iv + interp * (_pv - _iv);
  }else{					/* init val	*/
    ev = _iv;
  }
  trace2("EVAL_BM_PULSE::tr_eval", ev, d->long_label());
  assert(is_number(ev));
  //d->q_accept();
  tr_finish_tdv(d, ev);
}
/*--------------------------------------------------------------------------*/
TIME_PAIR EVAL_BM_PULSE::tr_review(COMPONENT* d)const
{
  ELEMENT* e = prechecked_cast<ELEMENT*>(d);
  assert(e);

  double time = d->_sim->_time0;
  double eps = d->_sim->_dtmin * .01;
  time += eps; // hack to avoid duplicate events from numerical noise
  double raw_time = time;

  if (0 < _period && _period < BIGBIG) {
    if (time > _delay) {
      time = fmod(time - _delay, _period) + _delay;
    }else{
    }
  }else{
  }
  double time_offset = raw_time - time;

  if (time >= _delay+_rise+_width+_fall) {		/* past pulse	*/
    d->_time_by.min_event(_delay + _period + time_offset);
    if (d->_sim->_time0 < _delay + _rise + _width + _fall + eps + time_offset) {
      e->_discont |= disSECOND;
      d->q_accept();
    }
  }else if (time >= _delay+_rise+_width) {		/* falling 	*/
    if (d->_sim->_time0 < _delay + _rise + _width + eps + time_offset) {
      e->_discont |= disSECOND;
      d->q_accept();
    }
    d->_time_by.min_event(_delay + _rise + _width + _fall + time_offset);
  }else if (time >= _delay + _rise) {			/* pulse val 	*/
    if (d->_sim->_time0 < _delay + _rise + eps + time_offset) {
      e->_discont |= disSECOND;
      d->q_accept();
    }
    d->_time_by.min_event(_delay + _rise + _width + time_offset);
  }else if (time >= _delay) {				/* rising 	*/
    if (d->_sim->_time0 < _delay + eps + time_offset) {
      e->_discont |= disSECOND;
      d->q_accept();
    }
    d->_time_by.min_event(_delay + _rise + time_offset);
  }else{						/* init val	*/
    d->_time_by.min_event(_delay + time_offset);
  }

  return d->_time_by;
}
/*--------------------------------------------------------------------------*/
bool EVAL_BM_PULSE::parse_numlist(CS& cmd)
{
  unsigned start = cmd.cursor();
  unsigned here = cmd.cursor();
  for (PARAMETER<double>* i = &_iv_in;  i < &_end;  ++i) {
    PARAMETER<double> val(NOT_VALID);
    cmd >> val;
    if (cmd.stuck(&here)) {
      break;
    }else{
      *i = val;
    }
  }
  return cmd.gotit(start);
}
/*--------------------------------------------------------------------------*/
bool EVAL_BM_PULSE::parse_params_obsolete_callback(CS& cmd)
{ trace1("EVAL_BM_PULSE::parse_params_obsolete_callback", cmd.tail());
  return ONE_OF
    || Get(cmd, "iv", 	  &_iv_in)
    || Get(cmd, "pv", 	  &_pv_in)
    || Get(cmd, "delay",  &_delay)
    || Get(cmd, "rise",   &_rise_in)
    || Get(cmd, "fall",   &_fall_in)
    || Get(cmd, "width",  &_width)
    || Get(cmd, "period", &_period_in)
    || Get(cmd, "area",   &_area)
    || Get(cmd, "duty{cycle}", &_duty)
    || Get(cmd, "freq",   &_freq)
    || Get(cmd, "frequency",   &_freq)
    || Get(cmd, "phase",   &_phase)
    || EVAL_BM_ACTION_BASE::parse_params_obsolete_callback(cmd)
    ;
}
/*--------------------------------------------------------------------------*/
std::map<string, PARA_BASE EVAL_BM_PULSE::*> EVAL_BM_PULSE::_param_dict =
  boost::assign::map_list_of
    ("iv",    (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_iv_in)
    ("U1",    (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_iv_in)
    ("val0",  (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_iv_in)
    ("pv",    (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_pv_in)
    ("U1",    (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_pv_in)
    ("val1",  (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_pv_in)
    ("delay", (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_delay)
    ("T1",    (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_delay)
    ("rise",  (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_rise_in)
    ("Tr",    (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_rise_in)
    ("T2",    (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_t2)
    ("width", (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_width)
    ("fall",  (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_fall_in)
    ("Tf",    (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_fall_in)
    ("period",(PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_period_in)
    ("freq",  (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_freq)
    ("freq{ency}",(PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_freq)
    ("area",  (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_area)
    ("duty",  (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_duty)
    ("duty{cycle}",(PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_duty)
    ("ph",    (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_phase)
    ("phase", (PARA_BASE EVAL_BM_PULSE::*) &EVAL_BM_PULSE::_phase);
/*--------------------------------------------------------------------------*/
void EVAL_BM_PULSE::set_param_by_name(std::string Name, std::string Value)
{ untested();
  PARA_BASE EVAL_BM_PULSE::* x = (_param_dict[Name]);
  if(x) { itested();
    PARA_BASE* p = &(this->*x);
    *p = Value;
  }else{
    throw Exception_No_Match(Name);
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
EVAL_BM_PULSE p1(CC_STATIC);
DISPATCHER<COMMON_COMPONENT>::INSTALL d1(&bm_dispatcher, "pulse", &p1);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
