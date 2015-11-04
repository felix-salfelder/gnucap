/*$Id: bm_sin.cc,v 1.3 2009-12-13 17:55:01 felix Exp $ -*- C++ -*-
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
 * SPICE compatible SIN
 */
//testing=script,complete 2005.10.07
#include "e_elemnt.h"
#include "u_lang.h"
#include "l_denoise.h"
#include "bm.h"
/*--------------------------------------------------------------------------*/
namespace {
/*--------------------------------------------------------------------------*/
using std::map;
using std::string;
/*--------------------------------------------------------------------------*/
const double _default_offset	(0);
const double _default_amplitude	(1);
const double _default_frequency (NOT_INPUT);
const double _default_delay	(0);
const double _default_damping	(0);
const double _default_samples	(4);
const bool   _default_zero	(false);
const bool   _default_peak	(false);
/*--------------------------------------------------------------------------*/
class EVAL_BM_SIN : public EVAL_BM_ACTION_BASE {
private:
  PARAMETER<double> _offset;
  PARAMETER<double> _amplitude;
  PARAMETER<double> _frequency;
  PARAMETER<double> _delay;
  PARAMETER<double> _damping;
  PARAMETER<double> _end;
  PARAMETER<double> _samples;
  PARAMETER<bool>   _zero;
  PARAMETER<bool>   _peak;
  void set_param_by_name(string Name, string Value);

  mutable double _actual_frequency;
  static std::map<string, PARA_BASE EVAL_BM_SIN::*> _param_dict;
  explicit	EVAL_BM_SIN(const EVAL_BM_SIN& p);
public:
  explicit      EVAL_BM_SIN(int c=0);
		~EVAL_BM_SIN()		{}
  int param_count()const {return 8 + EVAL_BM_ACTION_BASE::param_count();}
  string param_name(int i)const;
  string param_name(int i,int)const{return param_name(i);}
  string param_value(int)const;
  bool param_is_printable(int i)const;
private: // override vitrual
  bool		operator==(const COMMON_COMPONENT&)const;
  COMMON_COMPONENT* clone()const	{return new EVAL_BM_SIN(*this);}
  void		print_common_obsolete_callback(OMSTREAM&, LANGUAGE*)const;

  void		precalc_first(const CARD_LIST*);
  void		tr_eval(ELEMENT*)const;
  TIME_PAIR	tr_review(COMPONENT*)const;
  std::string	name()const		{return "sin";}
  bool		ac_too()const		{return false;}
  bool		parse_numlist(CS&);
  bool		parse_params_obsolete_callback(CS&);
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
EVAL_BM_SIN::EVAL_BM_SIN(int c)
  :EVAL_BM_ACTION_BASE(c),
   _offset(_default_offset),
   _amplitude(_default_amplitude),
   _frequency(_default_frequency),
   _delay(_default_delay),
   _damping(_default_damping),
   _end(NOT_VALID),
   _samples(_default_samples),
   _zero(_default_zero),
   _peak(_default_peak),
   _actual_frequency(0)
{
}
/*--------------------------------------------------------------------------*/
EVAL_BM_SIN::EVAL_BM_SIN(const EVAL_BM_SIN& p)
  :EVAL_BM_ACTION_BASE(p),
   _offset(p._offset),
   _amplitude(p._amplitude),
   _frequency(p._frequency),
   _delay(p._delay),
   _damping(p._damping),
   _end(NOT_VALID),
   _samples(p._samples),
   _zero(p._zero),
   _peak(p._peak),
   _actual_frequency(p._actual_frequency)
{
}
/*--------------------------------------------------------------------------*/
bool EVAL_BM_SIN::operator==(const COMMON_COMPONENT& x)const
{
  const EVAL_BM_SIN* p = dynamic_cast<const EVAL_BM_SIN*>(&x);
  bool rv = p
    && _offset == p->_offset
    && _amplitude == p->_amplitude
    && _frequency == p->_frequency
    && _delay == p->_delay
    && _damping == p->_damping
    && _samples == p->_samples
    && _zero == p->_zero
    && _peak == p->_peak
    && EVAL_BM_ACTION_BASE::operator==(x);
  return rv;
}
/*--------------------------------------------------------------------------*/
void EVAL_BM_SIN::print_common_obsolete_callback(OMSTREAM& o, LANGUAGE* lang)const
{ untested();
  assert(lang);
  o << name();
  print_pair(o, lang, "offset",	   _offset);
  print_pair(o, lang, "amplitude", _amplitude);
  print_pair(o, lang, "frequency", _frequency);
  print_pair(o, lang, "delay",     _delay,   _delay.has_hard_value());
  print_pair(o, lang, "damping",   _damping, _damping.has_hard_value());
  print_pair(o, lang, "samples",   _samples, _samples.has_hard_value());
  print_pair(o, lang, "zero",      _zero,    _zero.has_hard_value());
  print_pair(o, lang, "peak",      _peak,    _peak.has_hard_value());
  EVAL_BM_ACTION_BASE::print_common_obsolete_callback(o, lang);
}
/*--------------------------------------------------------------------------*/
void EVAL_BM_SIN::precalc_first(const CARD_LIST* Scope)
{
  assert(Scope);
  EVAL_BM_ACTION_BASE::precalc_first(Scope);
  _offset.e_val(_default_offset, Scope);
  _amplitude.e_val(_default_amplitude, Scope);
  _frequency.e_val(_default_frequency, Scope);
  _delay.e_val(_default_delay, Scope);
  _damping.e_val(_default_damping, Scope);
  _samples.e_val(_default_samples, Scope);
  _zero.e_val(_default_zero, Scope);
  _peak.e_val(_default_peak, Scope);
}
/*--------------------------------------------------------------------------*/
void EVAL_BM_SIN::tr_eval(ELEMENT* d)const
{
  if (!_frequency.has_hard_value()) {
    _actual_frequency = d->_sim->_freq;
  }else{
    _actual_frequency = _frequency;
  }

  double reltime = ioffset(d->_sim->_time0);
  double ev = _offset;
  if (reltime > _delay) {
    double x=_amplitude*fixzero(sin(M_TWO_PI*_actual_frequency*(reltime-_delay)),1.);
    if (_damping != 0.) {
      x *= exp(-(reltime-_delay)*_damping);
    }else{
    }
    ev += x;
  }else{
  }
  tr_finish_tdv(d, ev);
}
/*--------------------------------------------------------------------------*/
TIME_PAIR EVAL_BM_SIN::tr_review(COMPONENT* d)const
{
  double eps = d->_sim->_dtmin * .01;
  double time = ioffset(d->_sim->_time0);
  double reltime = time + eps;
  ELEMENT* e = prechecked_cast<ELEMENT*>(d);
  assert(e);

  if (reltime > _delay) {
    if (_peak && _zero) {
      d->_time_by.min_event(floor(reltime * 4 * _actual_frequency + 1) / (4 * _actual_frequency));
    }else if (_peak) {
      d->_time_by.min_event((floor(reltime*2*_actual_frequency+.5) + .5) / (2*_actual_frequency));
    }else if (_zero) {
      d->_time_by.min_event(floor(reltime * 2 * _actual_frequency + 1) / (2 * _actual_frequency));
    }else{
    }
    d->_time_by.min_error_estimate(d->_sim->_time0 + 1. / (_samples * _actual_frequency));
    if (time < _delay+eps) {
      e->_discont |= disSECOND;
      d->q_accept();
    }
  }else{
    d->_time_by.min_event(_delay);
  }

  return d->_time_by;
}
/*--------------------------------------------------------------------------*/
bool EVAL_BM_SIN::parse_numlist(CS& cmd)
{
  unsigned start = cmd.cursor();
  unsigned here = cmd.cursor();
  for (PARAMETER<double>* i = &_offset;  i < &_end;  ++i) {
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
bool EVAL_BM_SIN::parse_params_obsolete_callback(CS& cmd)
{
  return ONE_OF
    || Get(cmd, "o{ffset}",	&_offset)
    || Get(cmd, "a{mplitude}",	&_amplitude)
    || Get(cmd, "f{requency}",	&_frequency)
    || Get(cmd, "de{lay}",	&_delay)
    || Get(cmd, "da{mping}",	&_damping)
    || Get(cmd, "sa{mples}",	&_samples)
    || Get(cmd, "ze{ro}", 	&_zero)
    || Get(cmd, "pe{ak}", 	&_peak)
    || EVAL_BM_ACTION_BASE::parse_params_obsolete_callback(cmd)
    ;
}
/*--------------------------------------------------------------------------*/
bool EVAL_BM_SIN::param_is_printable(int i)const
{
  switch (EVAL_BM_SIN::param_count() - 1 - i) { untested();
  case 0:  return (true);
  case 1:  return (true);
  case 2:  return (true);
  case 3:  return _delay.has_hard_value();
  case 4:  return _damping.has_hard_value();
  case 5:  return _samples.has_hard_value();
  case 6:  return _zero.has_hard_value();
  case 7:  return _peak.has_hard_value();
  default: return EVAL_BM_ACTION_BASE::param_is_printable(i);
  }
}
/*--------------------------------------------------------------------------*/
string EVAL_BM_SIN::param_name(int i)const
{
  switch (EVAL_BM_SIN::param_count() - 1 - i) { untested();
  case 0:  return "offset";
  case 1:  return "amplitude";
  case 2:  return "frequency";
  case 3:  return "delay";
  case 4:  return "damping";
  case 5:  return "samples";
  case 6:  return "zero";
  case 7:  return "peak";
  default: return EVAL_BM_ACTION_BASE::param_name(i);
  }
}
/*--------------------------------------------------------------------------*/
string EVAL_BM_SIN::param_value(int i)const
{
  switch (EVAL_BM_SIN::param_count() - 1 - i) { untested();
  case 0:  return _offset.string();
  case 1:  return _amplitude.string();
  case 2:  return _frequency.string();
  case 3:  return _delay.string();
  case 4:  return _damping.string();
  case 5:  return _samples.string();
  case 6:  return _zero.string();
  case 7:  return _peak.string();
  default: return EVAL_BM_ACTION_BASE::param_value(i);
  }
}
/*--------------------------------------------------------------------------*/
map<string, PARA_BASE EVAL_BM_SIN::*> EVAL_BM_SIN::_param_dict =
  boost::assign::map_list_of
    ("offset",   (PARA_BASE EVAL_BM_SIN::*) &EVAL_BM_SIN::_offset)
    ("o",        (PARA_BASE EVAL_BM_SIN::*) &EVAL_BM_SIN::_offset)
    ("amplitude",(PARA_BASE EVAL_BM_SIN::*) &EVAL_BM_SIN::_amplitude)
    ("a",        (PARA_BASE EVAL_BM_SIN::*) &EVAL_BM_SIN::_amplitude)
    ("frequency",(PARA_BASE EVAL_BM_SIN::*) &EVAL_BM_SIN::_frequency)
    ("f",        (PARA_BASE EVAL_BM_SIN::*) &EVAL_BM_SIN::_frequency)
    ("delay",    (PARA_BASE EVAL_BM_SIN::*) &EVAL_BM_SIN::_delay)
    ("de",       (PARA_BASE EVAL_BM_SIN::*) &EVAL_BM_SIN::_delay)
    ("damping",  (PARA_BASE EVAL_BM_SIN::*) &EVAL_BM_SIN::_damping)
    ("da",       (PARA_BASE EVAL_BM_SIN::*) &EVAL_BM_SIN::_damping)
    ("samples",  (PARA_BASE EVAL_BM_SIN::*) &EVAL_BM_SIN::_samples)
    ("sa",       (PARA_BASE EVAL_BM_SIN::*) &EVAL_BM_SIN::_samples)
    ("zero",     (PARA_BASE EVAL_BM_SIN::*) &EVAL_BM_SIN::_zero)
    ("ze",       (PARA_BASE EVAL_BM_SIN::*) &EVAL_BM_SIN::_zero)
    ("peak",     (PARA_BASE EVAL_BM_SIN::*) &EVAL_BM_SIN::_peak)
    ("pe",       (PARA_BASE EVAL_BM_SIN::*) &EVAL_BM_SIN::_peak);
/*--------------------------------------------------------------------------*/
void EVAL_BM_SIN::set_param_by_name(std::string Name, std::string Value)
{
  PARA_BASE EVAL_BM_SIN::* x = (_param_dict[Name]);
  if (x) {
    PARA_BASE* p = &(this->*x);
    *p = Value;
  }else{ untested();
    EVAL_BM_ACTION_BASE::set_param_by_name(Name, Value);
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
EVAL_BM_SIN p1(CC_STATIC);
DISPATCHER<COMMON_COMPONENT>::INSTALL d1(&bm_dispatcher, "sin|sine", &p1);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
