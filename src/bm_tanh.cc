/*                             -*- C++ -*-
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
 * behavioral modeling hyperbolic tangent
 */
//testing=script 2005.10.07
#include "u_lang.h"
#include "e_elemnt.h"
#include "bm.h"
/*--------------------------------------------------------------------------*/
namespace {
/*--------------------------------------------------------------------------*/
const double _default_gain (NOT_INPUT);
const double _default_limit (NOT_INPUT);
/*--------------------------------------------------------------------------*/
class EVAL_BM_TANH : public EVAL_BM_ACTION_BASE {
private:
  PARAMETER<double> _gain;
  PARAMETER<double> _limit;
  PARAMETER<double> _shear;
  explicit	EVAL_BM_TANH(const EVAL_BM_TANH& p);
public:
  explicit      EVAL_BM_TANH(int c=0);
		~EVAL_BM_TANH()		{}
private: // override virtual
  bool		operator==(const COMMON_COMPONENT&)const;
  COMMON_COMPONENT* clone()const	{return new EVAL_BM_TANH(*this);}
  void		print_common_obsolete_callback(OMSTREAM&, LANGUAGE*)const;

  void		precalc_last(const CARD_LIST*);
  void		tr_eval(ELEMENT*)const;
  std::string	name()const		{return "tanh";}
  bool		ac_too()const		{untested();return false;}
  bool		parse_numlist(CS&);
  bool		parse_params_obsolete_callback(CS&);
  bool		param_is_printable(int i)const;
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
EVAL_BM_TANH::EVAL_BM_TANH(int c)
  :EVAL_BM_ACTION_BASE(c),
   _gain(NOT_INPUT),
   _limit(NOT_INPUT),
   _shear(NOT_INPUT)
{
}
/*--------------------------------------------------------------------------*/
EVAL_BM_TANH::EVAL_BM_TANH(const EVAL_BM_TANH& p)
  :EVAL_BM_ACTION_BASE(p),
   _gain(p._gain),
   _limit(p._limit),
   _shear(p._shear)
{
}
/*--------------------------------------------------------------------------*/
bool EVAL_BM_TANH::operator==(const COMMON_COMPONENT& x)const
{
  const EVAL_BM_TANH* p = dynamic_cast<const EVAL_BM_TANH*>(&x);
  bool rv = p
    && _gain == p->_gain
    && _limit == p->_limit
    && _shear == p->_shear
    && EVAL_BM_ACTION_BASE::operator==(x);
  return rv;
}
/*--------------------------------------------------------------------------*/
void EVAL_BM_TANH::print_common_obsolete_callback(OMSTREAM& o, LANGUAGE* lang)const
{
  assert(lang);
  o << name();
  print_pair(o, lang, "gain", _gain);
  print_pair(o, lang, "limit", _limit);
  if(_shear.has_hard_value()) { untested();
	  print_pair(o, lang, "shear", _shear);
  }
  EVAL_BM_ACTION_BASE::print_common_obsolete_callback(o, lang);
}
/*--------------------------------------------------------------------------*/
void EVAL_BM_TANH::precalc_last(const CARD_LIST* Scope)
{
  assert(Scope);
  EVAL_BM_ACTION_BASE::precalc_last(Scope);
  _gain.e_val(_default_gain, Scope);
  _limit.e_val(_default_limit, Scope);
  _shear.e_val(0, Scope);
}
/*--------------------------------------------------------------------------*/
void EVAL_BM_TANH::tr_eval(ELEMENT* d)const
{
  double x = ioffset(d->_y[0].x);
  double aa = x * _gain/_limit;
  double f1, f0;
  if (aa > LOGBIGBIG) {
    f1 = _shear;
    f0 = _shear * x + _limit;
  }else if (aa < -LOGBIGBIG) {
    f1 = _shear;
    f0 = _shear * x -_limit;
  }else{
    double cosine = cosh(aa);
    f1 = _shear + _gain / (cosine*cosine);
    f0 = _shear * x + _limit * tanh(aa);
  }
  d->_y[0] = FPOLY1(x, f0, f1);
  tr_final_adjust(&(d->_y[0]), d->f_is_value());
}
/*--------------------------------------------------------------------------*/
bool EVAL_BM_TANH::parse_numlist(CS& cmd)
{
  unsigned here = cmd.cursor();
  PARAMETER<double> gain(NOT_VALID);
  PARAMETER<double> limit(NOT_VALID);
  cmd >> gain >> limit;
  if (cmd.gotit(here)) {
    _gain  = gain;
    _limit = limit;
    return true;
  }else{
    return false;
  }
}
/*--------------------------------------------------------------------------*/
bool EVAL_BM_TANH::parse_params_obsolete_callback(CS& cmd)
{
  return ONE_OF
    || Get(cmd, "gain",  &_gain)
    || Get(cmd, "limit", &_limit)
    || Get(cmd, "shear", &_shear)
    || EVAL_BM_ACTION_BASE::parse_params_obsolete_callback(cmd)
    ;
}
/*--------------------------------------------------------------------------*/
bool EVAL_BM_TANH::param_is_printable(int i)const
{
  switch (EVAL_BM_TANH::param_count() - 1 - i) { untested();
  case 0:  return (true);
  case 1:  return (true);
  case 2:  untested(); return _shear.has_hard_value();
  default: return EVAL_BM_ACTION_BASE::param_is_printable(i);
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
EVAL_BM_TANH p1(CC_STATIC);
DISPATCHER<COMMON_COMPONENT>::INSTALL d1(&bm_dispatcher, "tanh", &p1);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
