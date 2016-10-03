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
 * functions for fixed voltage sources
 * temporary kluge: it has resistance
 */
#include "e_elemnt.h"
#include "u_xprobe.h"
/*--------------------------------------------------------------------------*/
using namespace std;
/*--------------------------------------------------------------------------*/
namespace {
/*--------------------------------------------------------------------------*/
class DEV_VS : public ELEMENT {
private:
  explicit DEV_VS(const DEV_VS& p) :ELEMENT(p) {}
public:
  explicit DEV_VS()		:ELEMENT() {}
private: // override virtual
  char	   id_letter()const	{return 'V';}
  std::string value_name()const {return "dc";}
  // void set_param_by_name(string, string);
  string element_type()const {return "vsource";}
  uint_t	   max_nodes()const	{return 2;}
  uint_t	   min_nodes()const	{return 2;}
  uint_t	   matrix_nodes()const	{return 2;}
  uint_t	   net_nodes()const	{return 2;}
  bool	   is_source()const	{return true;}
  bool	   f_is_value()const	{return true;}
  bool	   has_iv_probe()const  {return true;}
  bool	   use_obsolete_callback_parse()const {return true;}
  CARD*	   clone()const		{return new DEV_VS(*this);}
  void     precalc_last();
  void	   tr_iwant_matrix()	{
	  tr_iwant_matrix_passive();}
  void	   tr_begin();
  bool	   do_tr();
  void	   tr_load()		{tr_load_shunt(); tr_load_source();}
  void	   tr_unload()		{untested();tr_unload_source();}
  hp_float_t   tr_involts()const	{return 0.;}
  hp_float_t   tr_involts_limited()const {unreachable(); return 0.;}
  void	   ac_iwant_matrix()	{ac_iwant_matrix_passive();}
  void	   ac_begin()	{_loss1 = _loss0 = 1./OPT::shortckt; _acg = _ev = 0.;}
  void	   do_ac();
  void	   ac_load()		{ac_load_shunt(); ac_load_source();}
  COMPLEX  ac_involts()const	{return 0.;}
  COMPLEX  ac_amps()const	{return (_acg + ac_outvolts()* (double)_loss0);}
  XPROBE sens_probe_ext(const std::string& x)const;
  void sens_load(const std::string&);

  std::string port_name(uint_t i)const {
    assert(i != INVALID_NODE);
    assert(i < 2);
    static std::string names[] = {"p", "n"};
    return names[i];
  }
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// no. ELEMENT must take care of this.
// sometimes value is stored in common...
// void DEV_VS::set_param_by_name(string Name, string Value)
// {
//   if (Umatch (Name,"dc|u")) { _value = Value; }
//   else{ ELEMENT::set_param_by_name(Name,Value); }
// }
/*--------------------------------------------------------------------------*/
void DEV_VS::precalc_last()
{
  trace2("DEV_VS::precalc_last()", value(), has_tr_eval());
  //ELEMENT::precalc_last();	//BUG// skip
  COMPONENT::precalc_last();
  set_constant(!has_tr_eval());
  set_converged(!has_tr_eval());
  set_constant(false);
}
/*--------------------------------------------------------------------------*/
void DEV_VS::tr_begin()
{
  trace0("DEV_VS::tr_begin()");
  ELEMENT::tr_begin();
  _y[0].x  = 0.;
  _y[0].f1 = value();
  _y1.f0 = _y[0].f0 = 0.;	//BUG// override
  _loss1 = _loss0 = 1./OPT::shortckt;
  _m0.x  = 0.;
  _m0.c0 = -_loss0 * _y[0].f1;
  _m0.c1 = 0.;
  _m1 = _m0;    
  if (!using_tr_eval()) {
    if (_n[OUT2].m_() == 0) {
      _sim->set_limit(value());
    }else if (_n[OUT1].m_() == 0) {
      _sim->set_limit(-value());
    }else{
      //BUG// don't set limit
    }
  }else{
  }
}
/*--------------------------------------------------------------------------*/
bool DEV_VS::do_tr()
{
  trace4("DEV_VS::do_tr", long_label(), value(), using_tr_eval(), tr_needs_eval());
  assert(_m0.x == 0.);
  if (using_tr_eval()) {
    _y[0].x = _sim->_time0;
    tr_eval();
    if (_n[OUT2].m_() == 0) {
      _sim->set_limit(_y[0].f1);
    }else if (_n[OUT1].m_() == 0) {
      _sim->set_limit(-_y[0].f1);
    }else{
      //BUG// don't set limit
    }
    store_values();
    q_load();
    _m0.c0 = -_loss0 * _y[0].f1;
    assert(_m0.c1 == 0.);
  }else{untested();
    assert(conchk(_loss0, 1./OPT::shortckt));
    assert(_y[0].x == 0.);
    assert(_y[0].f0 == 0.);
    assert(_y[0].f1 == value());
    assert(_m0.x == 0.);
    assert(conchk(_m0.c0, -_loss0 * _y[0].f1));
    assert(_m0.c1 == 0.);
    assert(_y1 == _y[0]);
    assert(converged());
  }

  trace3("", _y[0].x, _y[0].f0, _y[0].f1);

  return converged();
}
/*--------------------------------------------------------------------------*/
void DEV_VS::do_ac()
{
  if (using_ac_eval()) {
    ac_eval();
    _acg = -_loss0 * _ev;
  }else{itested();
    assert(_acg == 0.);
  }
}
/*--------------------------------------------------------------------------*/
XPROBE DEV_VS::sens_probe_ext(const std::string& x)const
{
  trace3("DEV_VS", x, value(), _loss0 );
  unsigned n1 = _n[OUT1].m_();
  unsigned n2 = _n[OUT2].m_();
  COMPLEX a = CKT_BASE::_sim->_sens[n1];
  COMPLEX b = CKT_BASE::_sim->_sens[n2];

  if (Umatch(x, "v ")) { untested();
    COMPLEX ddr = (a-b);
    return XPROBE( _loss0*ddr );
  }

  return ELEMENT::sens_probe_ext(x);
}
/*--------------------------------------------------------------------------*/
void DEV_VS::sens_load(const std::string&)
{ itested();
  unsigned n1 = _n[OUT1].m_();
  unsigned n2 = _n[OUT2].m_();
  CKT_BASE::_sim->_sens[n1] += _loss0;
  CKT_BASE::_sim->_sens[n2] -= _loss0;
}
/*--------------------------------------------------------------------------*/
DEV_VS p1;
DISPATCHER<CARD>::INSTALL d1(&device_dispatcher, "V|vsource", &p1);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
