/*                               -*- C++ -*-
 * Copyright (C) 2001 Albert Davis
 *           (C) 2014 Felix Salfelder
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
 * inductor base class used by various inductor implementations
 */
#include "e_storag.h"

namespace INDUCTANCE {

class DEV_INDUCTANCE : public STORAGE {
protected:
  explicit DEV_INDUCTANCE(const DEV_INDUCTANCE& p) 
    :STORAGE(p), _c_model(p._c_model) {}
public:
  explicit DEV_INDUCTANCE()
    :STORAGE(), _c_model(false) {}
public: // override virtual
  char	   id_letter()const	{return 'L';}
  std::string value_name()const {return "l";}
  std::string dev_type()const	{return "inductor";}
  uint_t	   max_nodes()const	{return 2;}
  uint_t	   min_nodes()const	{return 2;}
  uint_t	   net_nodes()const	{return 2;}
  uint_t	   int_nodes()const     {return (!_c_model) ? 0 : 1;}
  uint_t	   matrix_nodes()const	{return net_nodes() + int_nodes();}

  bool	   has_inode()const	{return _c_model;}
  bool	   has_iv_probe()const  {return true;}
  bool	   use_obsolete_callback_parse()const {return true;}
  CARD*	   clone()const		{return new DEV_INDUCTANCE(*this);}
  void     expand();
  void	   tr_iwant_matrix();
  void     tr_begin();
  void     tr_restore();
  void tr_accept();
  bool	   do_tr();
  void	   tr_load();
  void	   tr_unload();
  double   tr_involts()const	{return tr_outvolts();}
  double   tr_input()const;
  double   tr_involts_limited()const {return tr_outvolts_limited();}
  double   tr_input_limited()const;
  double   tr_amps()const;
  double   tr_probe_num(const std::string&)const;
  void	   ac_iwant_matrix();
  void	   ac_begin()		{_loss1 = _loss0 = ((!_c_model) ? 0. : 1.); _ev = _y[0].f1;}
  void	   do_ac();
  void	   ac_load();
  COMPLEX  ac_involts()const	{return ac_outvolts();}
  COMPLEX  ac_amps()const;

  void set_param_by_name(string Name, string Value);

  std::string port_name(uint_t i)const {itested();
    assert(i != INVALID_NODE);
    assert(i < 2);
    static std::string names[] = {"p", "n"};
    return names[i];
  }
  bool _c_model;
private:
  bool has_ic() const;
};
/*--------------------------------------------------------------------------*/
inline void DEV_INDUCTANCE::expand()
{
  STORAGE::expand();
  if (_sim->is_first_expand()) {
    if (!_c_model) {
      _n[IN1].set_to_ground(this);
    }else{
      _n[IN1].new_model_node("." + long_label() + ".i", this);
    }
  }else{untested();
  }
}
/*--------------------------------------------------------------------------*/
inline double DEV_INDUCTANCE::tr_probe_num(const std::string& x)const
{
  if (Umatch(x, "flux ")) {
    return _y[0].f0;
  }else if (Umatch(x, "ind{uctance} |l ")) { itested();
    return _y[0].f1;
  }else if (Umatch(x, "dldt ")) { untested();
    return (_y[0].f1 - _y[1].f1) / _dt;
  }else if (Umatch(x, "dl ")) { untested();
    return (_y[0].f1 - _y[1].f1);
  }else if (Umatch(x, "dfdt ")) { untested();
    return (_y[0].f0 - _y[1].f0) / _dt;
  }else if (Umatch(x, "inode ")) {
    return has_inode();
  }else if (Umatch(x, "dflux ")) { untested();
    return (_y[0].f0 - _y[1].f0);
  }else{
    return STORAGE::tr_probe_num(x);
  }
}
/*--------------------------------------------------------------------------*/
void DEV_INDUCTANCE::tr_iwant_matrix()
{
  if (!_c_model) {
    tr_iwant_matrix_passive();
  }else{
    assert(matrix_nodes() == 3);
    
    assert(_n[OUT1].m_() != (uint_t) INVALID_NODE);
    assert(_n[OUT2].m_() != (uint_t) INVALID_NODE);
    assert(_n[IN1].m_() != (uint_t) INVALID_NODE);
    
    _sim->_aa.iwant(_n[OUT1].m_(),_n[IN1].m_());
    _sim->_aa.iwant(_n[OUT2].m_(),_n[IN1].m_());
    
    _sim->_lu.iwant(_n[OUT1].m_(),_n[IN1].m_());
    _sim->_lu.iwant(_n[OUT2].m_(),_n[IN1].m_());
  }
}
/*--------------------------------------------------------------------------*/
void DEV_INDUCTANCE::tr_begin()
{
  trace6("DEV_INDUCTANCE::tr_begin", _i[0].f0, long_label(), _sim->_cont, tr_involts(), _i[0].x, tr_outvolts());
  if (0&& _sim->_cont) {
    _i[0].f0 = tr_involts();
  }
  STORAGE::tr_begin();
  if (0 && _sim->_cont) {
 // i.x = amps,  i.f0 = volts,      i.f1 = ohms
 //  BUG: move to tr accept (if !uic?)
    _i[0].f0 = tr_involts();
    assert(_i[0].x != 0);
    _i[0].f1 = tr_involts() / _i[0].x;
    _m0.x = _i[0].f0;
    _m0.c0 = _i[0].x; // amps. iof?
    _m0.c1 = 0; // _i[0].x / tr_involts();
  }
  if (_c_model) {
    _loss1 = _loss0 = 1.;
  } else {
    _loss1 = _loss0 = 0.;
  }
}
/*--------------------------------------------------------------------------*/
void DEV_INDUCTANCE::tr_restore()
{

  if(_sim->_time0 >= _time[0]){
    // nothing.
    STORAGE::tr_restore();
  }else if(!_c_model){ incomplete();
    // impossible?
    STORAGE::tr_restore();
  }else{ // _sim->_time0 < _time[0]
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
  }
}
/*--------------------------------------------------------------------------*/
void DEV_INDUCTANCE::tr_accept()
{
  if (_sim->_uic) {
    _i[0].f0 = tr_involts();
  }
  STORAGE::tr_accept();
  if(_sim->_uic){
    _i[0].f0 = tr_involts();
    // assert(_i[0].x != 0);
    _i[0].f1 = tr_involts() / (_i[0].x+1e-20);
    _m0.x = _i[0].f0;
    _m0.c0 = _i[0].x; // amps. iof?
    _m0.c1 = 0; // _i[0].x / tr_involts();
    if (_c_model) {
      double idot = - tr_involts() / (_y[0].f1 + 1e-20);
      _m0.c0 = idot * _y[0].f1;
      trace2("coil",  _y[0].f1, _loss0);
    }
    set_converged(false);
  }
}
/*--------------------------------------------------------------------------*/
// quick hack. dont know how to do this in general.
inline void DEV_INDUCTANCE::set_param_by_name(std::string Name, std::string Value)
{
  trace2("DEV_INDUCTANCE::set_param_by_name", Name, Value);
  if (Umatch(Name, value_name()) && !has_common()) {
    set_value(Value);
  } else if (Umatch(Name, value_name())) { untested();
    ELEMENT::set_param_by_name(value_name(), value().string());
  } else if (has_common()) { untested();
    ELEMENT::set_param_by_name(Name, Value);
  } else {
    COMMON_COMPONENT* c = bm_dispatcher["eval_bm_value"]->clone();
    c->set_param_by_name("=", value().string());
    c->set_param_by_name(Name, Value);
    assert(c);
    attach_common(c);
    assert(has_common());
  }
  if (const EVAL_BM_ACTION_BASE* x = dynamic_cast<const EVAL_BM_ACTION_BASE*>(common())) {
    USE(x);
    trace2("DEV_INDUCTANCE::set_param_by_name", long_label(), x->_ic);
  }else{
  }
}
/*--------------------------------------------------------------------------*/
bool DEV_INDUCTANCE::has_ic() const
{
  if (const EVAL_BM_ACTION_BASE* x = dynamic_cast<const EVAL_BM_ACTION_BASE*>(common())) {
    if (x->_ic != NOT_INPUT) {
      return true;
    } else { untested();
      trace2("has_ic", long_label(), x->_ic);
      return false;
    }
  } else {
    return false;
  }
}
/*--------------------------------------------------------------------------*/
bool DEV_INDUCTANCE::do_tr()
{
  if (using_tr_eval()) {
    _y[0].x = tr_input_limited(); // _m0.c0 + _m0.c1 * x;
    tr_eval();
    if ((!_c_model) && (_y[0].f1 == 0.)) {untested();
      error(bDANGER, long_label() + ": short circuit,  L = 0\n");
      _y[0].f1 = OPT::shortckt;
      set_converged(conv_check());
    }else{
    }
  }else{
    _y[0].x = tr_input(); // _m0.c0 + _m0.c1 * x;
    assert(_y[0].f1 == value());
    _y[0].f0 = _y[0].x * _y[0].f1;
    assert(converged());
  }
  store_values();
  q_load();

  // i is really voltage ..
  if (_sim->uic_now() && has_ic()) {
    trace4("imitating cs", _y[0].x, value(), _y[0].f1, tr_involts());
    // imitate current src...
     // i.x = amps,  i.f0 = volts,      i.f1 = ohms
    _i[0].x = _y[0].x;
    _i[0].f1 = 0;
    if (!_c_model) {
      _m0.x = 0.;
      _m0.c0 = _y[0].x;
      _m0.c1 = 0.;
    } else {
      assert(_loss0==1);
    }
  } else {
    trace4("q", _y[1].x, _y[1].f0, _y[1].f1, _sim->_time0);
    trace4("q", _y[0].x, _y[0].f0, _y[0].f1, _sim->_time0);
    trace3("b4", _i[0].x, _i[0].f0, _i[0].f1);
    _i[0] = differentiate(_y, _i, _time, _method_a);
    trace3("d ", _i[0].x, _i[0].f0, _i[0].f1);
  }

  if (!_c_model) {
    _m0.x = NOT_VALID;
    _m0.c1 = 1 / ((_i[0].c1()==0) ? OPT::shortckt : _i[0].c1());
    if (_sim->uic_now() && has_ic()) {
      _m0.c0 = _y[0].x;
      _y[0].f0 = _y[0].x * _y[0].f1;
      _m0.c1 = 0;
    } else {
      _m0.c0 = -_i[0].c0() * _m0.c1;
    }
  } else {
    if (_sim->uic_now() && has_ic()) {
      // i.x = amps,  i.f0 = volts,      i.f1 = ohms
      _i[0].f0 = tr_involts();
      _i[0].f1 = 1. / OPT::shortckt;
      // m.x = volts, m.c0 = amps, acg = m.c1 = mhos
      _m0.x = 0; //_y[0].x;
      _m0.c1 = -_loss0 * _loss0 * _i[0].c1();
      _m0.c0 =  _loss0 * _loss0 * _i[0].c0(); //  + tr_involts(); // hmm
    } else {
      _m0.x = NOT_VALID;
      _m0.c1 = -_loss0 * _loss0 * _i[0].c1();
      _m0.c0 =  _loss0 * _loss0 * _i[0].c0();
    }
  }

  trace7("L::do_tr", _i[0].c0(), _i[0].c1(), _y[0].f1, tr_amps(), tr_involts(), _m0.c0, _m0.c1);
  return converged();
}
/*--------------------------------------------------------------------------*/
void DEV_INDUCTANCE::tr_load()
{
  if (!_c_model) {
    tr_load_passive();
  }else{
    trace6("DEV_INDUCTANCE::tr_load", _m0.c1, _m1.c1, _m0.c0, _m1.c0, _loss0, _loss1);
    tr_load_inode(); // load loss.
    tr_load_diagonal_point(_n[IN1], &_m0.c1, &_m1.c1);
    tr_load_source_point(_n[IN1], &_m0.c0, &_m1.c0); // load i
  }
}
/*--------------------------------------------------------------------------*/
void DEV_INDUCTANCE::tr_unload()
{untested();
  _loss0 = _m0.c0 = _m0.c1 = 0.;
  _sim->mark_inc_mode_bad();
  tr_load();
}
/*--------------------------------------------------------------------------*/
double DEV_INDUCTANCE::tr_input()const
{
  if (!_c_model) {
    return _m0.c0 + _m0.c1 * tr_involts();
  }else{
    return _n[IN1].v0();
  }
}
/*--------------------------------------------------------------------------*/
double DEV_INDUCTANCE::tr_input_limited()const
{
  if (!_c_model) {
    return _m0.c0 + _m0.c1 * tr_involts_limited();
  }else{
    return _n[IN1].v0();
  }
}
/*--------------------------------------------------------------------------*/
double DEV_INDUCTANCE::tr_amps()const
{
  if (!_c_model) {
    // m0c0 is "flux"
    return fixzero((_m0.c1 * tr_involts() + _m0.c0), _m0.c0);
  }else{
    return _loss0 * _n[IN1].v0();
  }
}
/*--------------------------------------------------------------------------*/
void DEV_INDUCTANCE::ac_iwant_matrix()
{
  if (!_c_model) {
    ac_iwant_matrix_passive();
  }else{
    assert(matrix_nodes() == 3);
    
    assert(_n[OUT1].m_() != INVALID_NODE);
    assert(_n[OUT2].m_() != INVALID_NODE);
    assert(_n[IN1].m_() != INVALID_NODE);
    
    _sim->_acx.iwant(_n[OUT1].m_(),_n[IN1].m_());
    _sim->_acx.iwant(_n[OUT2].m_(),_n[IN1].m_());
  }
}
/*--------------------------------------------------------------------------*/
void DEV_INDUCTANCE::do_ac()
{
  if (using_ac_eval()) {
    ac_eval();
  }else{
    assert(_ev == _y[0].f1);
    // uuh
//    assert(dynamic_cast<INDUCTANCE::DEV_MUTUAL_L*>(this) || has_tr_eval() || _ev == double(value()));
  }
  if (!_c_model) {
    if ((COMPLEX)_ev * _sim->_jomega == 0.) {untested();
      _acg = 1. / OPT::shortckt;
    }else{
      _acg = 1. / ((COMPLEX)_ev * _sim->_jomega);
    }
  }else{
    _acg = -(double)_loss0 *(double) _loss0 * (COMPLEX)_ev * _sim->_jomega;
  }
}
/*--------------------------------------------------------------------------*/
void DEV_INDUCTANCE::ac_load()
{
  if (!_c_model) {
    ac_load_passive();
  }else{
    ac_load_inode(); // 4x \pm loss.
    ac_load_diagonal_point(_n[IN1], _acg);
  }
}
/*--------------------------------------------------------------------------*/
COMPLEX DEV_INDUCTANCE::ac_amps()const
{
  if (!_c_model) {
    return (ac_involts() * _acg);
  }else{
    return (  (double)_loss0 * (COMPLEX)(_n[IN1].vac()) );
  }
}
/*--------------------------------------------------------------------------*/
} // INDUCTANCE
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
