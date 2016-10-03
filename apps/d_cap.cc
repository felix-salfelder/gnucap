/*$Id: d_cap.cc,v 1.10 2010-08-26 09:07:17 felix Exp $ -*- C++ -*-
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
 * capacitance devices:
 *	self-capacitance (C device)
 *	trans-capacitance (non-spice charge transfer device)
 *------------------------------------------------------------------
 * capacitor models
 * y.x = volts, y.f0 = coulombs, ev = y.f1 = farads
 * q = y history in time
 * i.x = volts, i.f0 = amps,	      i.f1 = mhos
 * m.x = volts, m.c0 = amps,    acg = m.c1 = mhos
 */
//testing=script 2006.07.17
#include "d_cap.h"
#include "e_storag.h"
namespace SOME_CAP_HACK {
/*--------------------------------------------------------------------------*/
class DEV_TRANSCAP : public DEV_CAPACITANCE {
private:
  explicit DEV_TRANSCAP(const DEV_TRANSCAP& p) :DEV_CAPACITANCE(p){}
public:
  explicit DEV_TRANSCAP()	:DEV_CAPACITANCE() {}
private: // override virtual
  char     id_letter()const	{untested();return '\0';}
  std::string value_name()const {untested(); return "c";}
  std::string dev_type()const	{return "tcap";}
  uint_t	   max_nodes()const	{return 4;}
  uint_t	   min_nodes()const	{return 4;}
  uint_t	   matrix_nodes()const	{return 4;}
  uint_t	   net_nodes()const	{return 4;}
  bool	   has_iv_probe()const  {untested(); return false;}
  bool	   f_is_value()const	{untested();return true;}
  CARD*	   clone()const		{return new DEV_TRANSCAP(*this);}
  void	   tr_iwant_matrix()	{tr_iwant_matrix_active();}
  void	   tr_load()		{tr_load_active();}
  hp_float_t   tr_involts()const	{return dn_diff(_n[IN1].v0(),_n[IN2].v0());}
  hp_float_t   tr_involts_limited()const {return volts_limited(_n[IN1],_n[IN2]);}
  void	    ac_iwant_matrix()	{ac_iwant_matrix_active();}
  void	    ac_load()		{untested(); ac_load_active();}
  std::string port_name(uint_t i)const {untested();
    assert(i != INVALID_NODE);
    assert(i < 4);
    static std::string names[] = {"p", "n", "ps", "ns"};
    return names[i];
  }
};
/*--------------------------------------------------------------------------*/
//BUG// doesn't model dynamic effects of control.
class DEV_VCCAP : public DEV_CAPACITANCE {
private:
  explicit DEV_VCCAP(const DEV_VCCAP& p) :DEV_CAPACITANCE(p) {}
public:
  explicit DEV_VCCAP()		:DEV_CAPACITANCE() {}
private: // override virtual
  char     id_letter()const	{untested();return '\0';}
  std::string value_name()const {untested(); return "c";}
  std::string dev_type()const	{return "vccap";}
  uint_t	   max_nodes()const	{return 4;}
  uint_t	   min_nodes()const	{return 4;}
  uint_t	   matrix_nodes()const	{return 4;}
  uint_t	   net_nodes()const	{return 4;}
  bool	   has_iv_probe()const  {untested(); return false;}
  bool	   f_is_value()const	{untested();return true;}
  CARD*	   clone()const		{return new DEV_VCCAP(*this);}
  void	   tr_iwant_matrix()	{tr_iwant_matrix_extended();}
  bool     do_tr();
  void uic_clanup();
  hp_float_t   tr_involts()const	{return dn_diff(_n[IN1].v0(),_n[IN2].v0());}
  hp_float_t   tr_involts_limited()const {return volts_limited(_n[IN1],_n[IN2]);}
  void	    ac_iwant_matrix()	{ac_iwant_matrix_extended();}

  std::string port_name(uint_t i)const {untested();
    assert(i != INVALID_NODE);
    assert(i < 4);
    static std::string names[] = {"p", "n", "ps", "ns"};
    return names[i];
  }

};
/*--------------------------------------------------------------------------*/
// obsolete
bool DEV_CAPACITANCE::has_ic() const
{
  if ( const EVAL_BM_ACTION_BASE* x = dynamic_cast<const EVAL_BM_ACTION_BASE*>(common()) ){
    if (x->_ic != NOT_INPUT) return 1;
  }
  return 0;
}
/*--------------------------------------------------------------------------*/
bool DEV_CAPACITANCE::do_tr()
{
  // FPOLY1* q=_y;
  //trace0(("DEV_CAPACITANCE::do_tr " + long_label()));
  //trace3(("DEV_CAPACITANCE::do_tr " + long_label()).c_str(), _y[0].f1, value(), tr_input() );
  if (using_tr_eval()) {
    trace1("DEV_CAPACITANCE::do_tr, tr_eval", has_tr_eval());
    _y[0].x = tr_input_limited();
    tr_eval();
  }else{
    trace1("DEV_CAPACITANCE::do_tr, no tr_eval", tr_input());
    _y[0].x = tr_input(); // tr_involts();
    assert(_y[0].f1 == value());
    _y[0].f0 = _y[0].x * _y[0].f1;
    assert(converged());
  }
  trace3("DEV_CAPACITANCE::do_tr ", _y[0].f1, _sim->uic_now(), is_constant() );
  store_values();
  q_load();

  if (_sim->uic_now() && has_ic()) {
    double G = 1./OPT::shortckt;
    assert(_time[0] == 0.);

    trace1("DEV_CAPACITANCE::do_tr: desired capacitance is ", _y[0].f1);
    trace1("DEV_CAPACITANCE::do_tr: desired charge is      ", _y[0].f0);
    trace1("DEV_CAPACITANCE::do_tr: desired voltage is     ", _y[0].x);

    // imitate voltage source... (d_vs.cc)
    _i[0] = FPOLY1( CPOLY1( 0., -_y[0].x * G,         G  ) ); 
//    trace2("2 quotienten", ( -_y[0].f0 / (OPT::shortckt)  ) /-_y[0].x  / (OPT::shortckt) ,
//			         (_y[0].f1 / (OPT::shortckt))/    1/(OPT::shortckt)  );

    // picked up by tr_load_*
    _loss0 = G;

  } else {
    trace4("q", _y[0].x, _y[0].f0, _y[0].f1, _sim->_phase);
    trace3("i", _i[0].x, _i[0].f0, _i[0].f1);
    _i[0] = differentiate(_y, _i, _time, _method_a);
    trace3("i", _i[0].x, _i[0].f0, _i[0].f1);
  }
  _m0 = CPOLY1(_i[0]);
  trace4("DEV_CAPACITANCE::do_tr done", long_label(), value(), _m0, _m1);
  return converged();
}
/*--------------------------------------------------------------------------*/
void DEV_CAPACITANCE::tr_accept()
{
  ELEMENT::tr_accept(); // commmon->accept..

  if (_loss0!=0) {
    // _loss0 should be irrelevant
    // sources might have been abused to enforce initial conditions
    // during dc analysis (otherwise this should be unreachable)

    if (0 && !_sim->uic_now()) { // s_ddc
      _m0.c0 = 0;
      trace2("DEV_CAPACITANCE::tr_accept ddc", _m0, _m1);
    } else {
      trace3("DEV_CAPACITANCE::tr_accept ddc", _m0, _m1, _sim->_uic);
      // vera hack. something about qdot
      // also used in dc uic + tr cont ...
      FPOLY1 m(_m0);

      m.x = tr_input();
      m.f0 = - _loss0 * tr_outvolts() - _m0.c0;
      m.f1 = 0;
      _m0 = CPOLY1(m);
      _i[0].f0=0;
      _i[0].f1=0;
    }

    _loss0 = 0.;
    set_converged(false);
  } else {
  }
}
/*--------------------------------------------------------------------------*/
void DEV_CAPACITANCE::do_ac()
{
  if (using_ac_eval()) {
    ac_eval();
  }else{
    assert(_ev == _y[0].f1);
    assert(has_tr_eval() || _ev == hp_float_t(value()));
  }
  trace1("DEV_CAPACITANCE::do_ac", _sim->_jomega);
  _acg =  (COMPLEX)_ev * _sim->_jomega;

}
/*--------------------------------------------------------------------------*/
double DEV_CAPACITANCE::tr_probe_num(const std::string& x)const
{
  if (Umatch(x, "q{cap} |ch{arge} ")) {
    return _y[0].f0;
  }else if (Umatch(x, "c{apacitance} ")) {
    return _y[0].f1;
  }else if (Umatch(x, "dcdt ")) {untested();
    return (_y[0].f1 - _y[1].f1) / _dt;
  }else if (Umatch(x, "dc ")) {untested();
    return (_y[0].f1 - _y[1].f1);
  }else if (Umatch(x, "dqdt ")) {
    return (_y[0].f0 - _y[1].f0) / _dt;
  }else if (Umatch(x, "dq ")) {
    return (_y[0].f0 - _y[1].f0);
  }else{
    return STORAGE::tr_probe_num(x);
  }
}
/*--------------------------------------------------------------------------*/
bool DEV_VCCAP::do_tr()
{
  _y[0].x = tr_input_limited();
  tr_eval();

  store_values();
  q_load();

  _y[0].x = tr_outvolts();
  _y[0].f1 = _y[0].f0;		 // self capacitance
  _y[0].f0 = _y[0].x * _y[0].f1; // charge

  _i[0] = differentiate(_y, _i, _time, _method_a);
  _m0.x  = _i[0].x;
  _m0.c1 = _i[0].f1;
  _m0.c0 = _i[0].f0 - _i[0].x * _i[0].f1;
  return converged();
}
/*--------------------------------------------------------------------------*/
void DEV_CAPACITANCE::tr_restore()
{

  if(_sim->_time0 < _time[0]){
    // recover from _freezetime ...

    _y[0].x = tr_input();
//    _i[0].x = tr_input();
//    _i[0].f0 = 0.; // amps
//    _i[0].f1 = 1./OPT::shortckt; // mhos
//    _i[1] = _i[0];
    if (using_tr_eval()) {
      assert(_y[0].f1 == value());
    }else{incomplete();
    }
    // BUG: compute _y[0].f1 from OP?
    _y[0].f0 = _y[0].x * _y[0].f1; // charge
    _y[1] = _y[0];
    _y1 = _y[0];
    _method_a = mINVALID;
    // _i[0] = differentiate(_y, _i, _time, _method_a);
    double G = 1./OPT::shortckt;
    _i[0] = FPOLY1( CPOLY1( 0., -_y[0].x * G,         G  ) ); 
    _m0 = CPOLY1(_i[0]);
    trace3("DEV_CAPACITANCE::tr_restore from freeze", _y[0], _i[0], _i[1]);
    STORAGE::tr_restore();
    _time[1] = 0.; /// hmm hack.
  }else{
    STORAGE::tr_restore();
  }

}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
DEV_CAPACITANCE p1;
DEV_TRANSCAP    p2;
DEV_VCCAP       p3;
DISPATCHER<CARD>::INSTALL
  d1(&device_dispatcher, "C|capacitor",	    &p1),
  d2(&device_dispatcher, "tcap|tcapacitor", &p2),
  d3(&device_dispatcher, "vccap",	    &p3);
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
