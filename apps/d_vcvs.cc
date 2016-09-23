/* Copyright (C) 2001 Albert Davis
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
 * functions for vcvs
 * temporary kludge: it has resistance
 */
#include "e_elemnt.h"
/*--------------------------------------------------------------------------*/
namespace {
/*--------------------------------------------------------------------------*/
class DEV_VCVS : public ELEMENT {
private:
  explicit DEV_VCVS(const DEV_VCVS& p) :ELEMENT(p) { }
public:
  explicit DEV_VCVS()		:ELEMENT() {} //{ _net_nodes=4; }
private: // override virtual
  char	   id_letter()const	{return 'E';}
  std::string value_name()const {return "gain";}
  std::string element_type()const {return "vcvs";}
  bool	   use_obsolete_callback_parse()const {return true;}
  CARD*	   clone()const		{return new DEV_VCVS(*this);}
  void     precalc_last();
  void     expand();
  void	   tr_iwant_matrix()	{tr_iwant_matrix_extended(); if(subdev())subdev()->tr_iwant_matrix();}
  void     tr_begin();
  bool	   do_tr();
  void	   tr_load();
  void	   tr_unload()		{untested();tr_unload_active();}
  hp_float_t tr_involts()const {return dn_diff(_n[IN1].v0(), _n[IN2].v0());}
  hp_float_t tr_involts_limited()const {return volts_limited(_n[IN1],_n[IN2]);}
  void	   ac_iwant_matrix()	{ac_iwant_matrix_extended(); if(subdev())subdev()->ac_iwant_matrix();}
  void	   ac_begin();
  void	   do_ac();
  void	   ac_load()		{ac_load_shunt(); ac_load_active();}
  COMPLEX  ac_involts()const	{return _n[IN1].vac() - _n[IN2].vac();}

  std::string port_name(uint_t i)const {
    assert(i !=INVALID_NODE);
    assert(i < 6);
    static std::string names[] = {"p", "n", "ps", "ns", "ps2", "ns2"};
    return names[i];
  }
  bool _master;
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
void DEV_VCVS::tr_load()
{
  if(subdev()){
    subdev()->tr_load();
  }
  tr_load_shunt();
  tr_load_active();
}
/*--------------------------------------------------------------------------*/
#if 0
unsigned DEV_VCVS::min_nodes() const
{
  if( const EVAL_BM_ACTION_BASE* e=dynamic_cast<const EVAL_BM_ACTION_BASE*>(common())){
    trace1("DEV_VCVS::min_nodes", e->input_order());
    return e->input_order() + 2;
  } else {
    trace0("DEV_VCVS::min_nodes 4");
  }
  return 4;
}
#endif
/*--------------------------------------------------------------------------*/
void DEV_VCVS::precalc_last()
{
  ELEMENT::precalc_last();
  set_constant(!has_tr_eval());
  set_converged(!has_tr_eval());
}
/*--------------------------------------------------------------------------*/
void DEV_VCVS::expand()
{
  trace3("DEV_VCVS::expand", long_label(), input_order(), hp(common()));
  assert(!subckt() || subckt()->size() < 2);
  ELEMENT::expand();
  _master = 1;
  if (has_common()) {
    if( const EVAL_BM_ACTION_BASE* e=dynamic_cast<const EVAL_BM_ACTION_BASE*>(common())){
      _master = input_order() == e->input_order();
    }
  }
}
/*--------------------------------------------------------------------------*/
void DEV_VCVS::tr_begin()
{
  trace1("DEV_VCVS::tr_begin", long_label());
  ELEMENT::tr_begin();
  double bm_ord = 1;
  if (has_common()) {
    if (const EVAL_BM_ACTION_BASE* e=dynamic_cast<const EVAL_BM_ACTION_BASE*>(common())) {
      bm_ord = e->input_order();
    }
  }
  _loss1 = _loss0 = bm_ord / OPT::shortckt;
  _m0.x  = 0.;
  _m0.c1 = -_loss0 * _y[0].f1 * bm_ord;
  _m0.c0 = 0.;
  _m1 = _m0;
  if(subdev()){
    subdev()->tr_begin();
  }
}
/*--------------------------------------------------------------------------*/
bool DEV_VCVS::do_tr()
{
  double bm_ord = 1;
  if (has_common()) {
    if (const EVAL_BM_ACTION_BASE* e=dynamic_cast<const EVAL_BM_ACTION_BASE*>(common())) {
      bm_ord = e->input_order();
    }
  }
  if (using_tr_eval()) {
    _y[0].x = _m0.x = tr_involts_limited();
    //_y[0].x = tr_input_limited();
    //assert(_y[0].x == _m0.x);
    tr_eval();
    assert(_y[0].f0 != LINEAR);
    store_values();
    trace5("DEV_VCVS::do_tr", long_label(), _loss0, _y[0].x, _y[0].f0, _y[0].f1);
    double oldf0=_y[0].f0;
    if (_master) {
      q_load();
    } else {
      _y[0].f0 = 0;
    }
    _m0 = CPOLY1(_y[0]);
    _y[0].f0=oldf0;
    _m0 *= -_loss0 * bm_ord;
  }else{
    assert(conchk(_loss0, 1./OPT::shortckt));
    assert(_y[0].f0 == LINEAR);
    assert(_y[0].f1 == value());
    assert(conchk(_m0.c1, -_loss0 * _y[0].f1));
    assert(_m0.c0 == 0.);
    assert(_y1 == _y[0]);
    assert(converged());
  }
  return converged();
}
/*--------------------------------------------------------------------------*/
void DEV_VCVS::ac_begin()
{
  _loss1 = _loss0 = 1./OPT::shortckt;
  _ev = _y[0].f1;
  _acg = -_loss0 * _ev;
  if(subdev()) subdev()->ac_begin();
}
/*--------------------------------------------------------------------------*/
void DEV_VCVS::do_ac()
{
  if (using_ac_eval()) {
    ac_eval();
    trace2("DEV_VCVS::do_ac", _loss0, _ev);
    _acg = -_loss0 * _ev;
  }else{
    assert(_ev == _y[0].f1);
    assert(has_tr_eval() || _ev == hp_float_t(value()));
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
DEV_VCVS p1;
DISPATCHER<CARD>::INSTALL d1(&device_dispatcher, "E|vcvs", &p1);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
