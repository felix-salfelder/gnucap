/* * Copyright (C) 2011 Felix Salfelder
 * Author: Felix Salfelder
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
 * a vs using a branch node for its current.
 * still incomplete?. needed for sock
 */
//testing=script 2006.07.17
#include "e_elemnt.h"
#define BR 2
using namespace std;
/*--------------------------------------------------------------------------*/
namespace {
/*--------------------------------------------------------------------------*/
class DEV_BVS : public ELEMENT {
private:
  explicit DEV_BVS(const DEV_BVS& p) :ELEMENT(p) {}
public:
  explicit DEV_BVS()		:ELEMENT() {}
private: // state
  double _one0, _one1;
  double _source0, _source1;
private: // override virtual
  char	   id_letter()const	{return 'V';}
  std::string value_name()const {return "dc";}
  void set_param_by_name(string, string);
#if 0
  std::string dev_type()const	{return "vsb";}
#else
  std::string element_type()const       {return "vsource";}
#endif
  uint_t	   max_nodes()const	{return 2;}
  uint_t	   min_nodes()const	{return 2;}
  uint_t	   matrix_nodes()const	{return 3;}
  uint_t	   net_nodes()const	{return 2;}
  uint_t           int_nodes()const	{return 1;}
  uint_t           ext_nodes()const	{return 2;}
  bool	   is_source()const	{return true;}
  bool	   f_is_value()const	{return true;}
  bool	   has_iv_probe()const  {return true;}
  bool	   use_obsolete_callback_parse()const {return true;}
  CARD*	   clone()const		{return new DEV_BVS(*this);}
  void     precalc_last();
  void	   tr_iwant_matrix();//	{ tr_iwant_matrix_passive();}
  void	   tr_begin();
  bool	   do_tr();
  void	   tr_load();
  void	   tr_unload();
  hp_float_t   tr_involts()const	{return 0.;}
  hp_float_t   tr_involts_limited()const {unreachable(); return 0.;}
  void	   ac_iwant_matrix()	{ac_iwant_matrix_passive();}
  void	   ac_begin();//	{_loss1 = _loss0 = 1./OPT::shortckt; _acg = _ev = 0.;}
  void	   do_ac();
  void	   ac_load();//		{ac_load_shunt(); ac_load_source();}
  COMPLEX  ac_involts()const	{return 0.;}
  COMPLEX  ac_amps()const	{return (_acg + ac_outvolts()* (double)_loss0);}
  double tr_amps()const;

  void sens_load(const std::string&);

  std::string port_name(uint_t i)const { untested();
    assert(i != INVALID_NODE);
    assert(i < 3);
    static std::string names[] = {"p", "n", "b"};
    return names[i];
  }

  void expand();
};
/*--------------------------------------------------------------------------*/
void DEV_BVS::set_param_by_name(string Name, string Value)
{
  if (Umatch (Name,"dc|u")) { _value = Value; }
  else{ ELEMENT::set_param_by_name(Name,Value); }
}
/*--------------------------------------------------------------------------*/
inline void DEV_BVS::ac_begin()
{
  _acg = _ev = 0.;
}
/*--------------------------------------------------------------------------*/
inline void DEV_BVS::ac_load()
{
  trace0("DEV_BVS::ac_load");
  //_sim->_acx.load_symmetric(_n[OUT1].m_(), _n[OUT2].m_(), mfactor() * _loss0);
  double d=1;
  _sim->_acx.load_point(_n[OUT1].m_(), _n[BR].m_(), d);
  _sim->_acx.load_point(_n[OUT2].m_(), _n[BR].m_(), -d);
  _sim->_acx.load_point(_n[BR].m_(), _n[OUT1].m_(), d);
  _sim->_acx.load_point(_n[BR].m_(), _n[OUT2].m_(), -d);

  assert (_n[BR].m_() != 0);
  //  untested(); // does this do the right thing?
  /* load source */
  _n[BR].iac() += _acg; // mfactor() * _acg;

  trace3("DEV_BVS::ac_load source ", _acg, _n[BR].m_(),BR);

}
/*--------------------------------------------------------------------------*/
void DEV_BVS::sens_load(const std::string&)
{ untested(); incomplete();
  unsigned n1 = _n[BR].m_();
  CKT_BASE::_sim->_sens[n1] += 1.;
}
/*--------------------------------------------------------------------------*/
void DEV_BVS::expand()
{
  trace1("DEV_BVS::expand ", long_label() );
  if (!subckt()) {
    new_subckt(); // hmm probably not a good idea
  }else{ untested();
  }

  if (!(_n[BR].n_())) {
    _n[BR].new_model_node( "br", this);
  }
  trace2("DEV_BVS::expand ", long_label(),  _n[BR].m_() );
  COMPONENT::expand();

}
/*--------------------------------------------------------------------------*/
void DEV_BVS::tr_iwant_matrix()
{
/*
 *   +  -  b
 * +       1    u+
 * -      -1    u-   
 * b 1 -1  0    ib  val
 */


  assert(matrix_nodes() == 3);
  assert(is_device());
  //assert(!subckt()); ok for subckt to exist for logic

  if (_n[OUT1].m_() == INVALID_NODE) {
    std::cerr << "ELEMENT::tr_iwant_matrix_passive " << long_label() << "\n";
    exit(4);
  }
  assert(_n[OUT1].m_() != INVALID_NODE);
  assert(_n[OUT2].m_() != INVALID_NODE);
  assert(_n[BR].m_() != INVALID_NODE);

  trace3("DEV_BVS::tr_iwant_matrix ", _n[OUT1].m_(),_n[OUT2].m_(),_n[BR].m_());

  _sim->_aa.iwant(_n[OUT1].m_(),_n[BR].m_());
  _sim->_aa.iwant(_n[OUT2].m_(),_n[BR].m_());
  _sim->_lu.iwant(_n[OUT1].m_(),_n[BR].m_());
  _sim->_lu.iwant(_n[OUT2].m_(),_n[BR].m_());
  _sim->_acx.iwant(_n[OUT1].m_(),_n[BR].m_());
  _sim->_acx.iwant(_n[OUT2].m_(),_n[BR].m_());
//  ac_iwant_matrix_extended();
//  tr_iwant_matrix_extended();
}
/*--------------------------------------------------------------------------*/
inline void DEV_BVS::tr_load()
{
  trace7("DEV_BVS::tr_load", long_label(), _n[OUT1].m_(), _n[OUT2].m_(), _n[BR].m_(), _sim->_time0, _sim->iteration_number(), _sim->is_inc_mode());
  assert(_one0 == _one0);

/*----------shunt---------------------------------------------------------------*/
  double d = _one0 - _one1;
  if (!_sim->is_advance_or_first_iteration()) {
    _one0 = _one1 + d;
  }
  d = ((_sim->is_inc_mode()) ? d : _one0);

  if (d != 0.) {
    trace3("DEV_BVS::tr_load 4 times",d,_one0,_one1);
    _sim->_aa.load_point(_n[OUT1].m_(), _n[BR].m_(), d);
    _sim->_aa.load_point(_n[OUT2].m_(), _n[BR].m_(), -d);
    _sim->_aa.load_point(_n[BR].m_(), _n[OUT1].m_(), d);
    _sim->_aa.load_point(_n[BR].m_(), _n[OUT2].m_(), -d);
  }

  _one1 = _one0;

/*--load_source------------------------------------------------------------------------*/
  trace3("DEV_BVS::tr_load", long_label(), _m0.c0, _m1.c0);
  assert(_m0.c0 == _m0.c0);
  assert (_n[BR].m_() != 0);
  d = dampdiff(&_m0.c0, _m1.c0); // hmmm
  trace1("DEV_BVS::tr_load source",d);
  if (d != 0.) {
    _n[BR].i() += d;
  }else{
  }
  _m1 = _m0;
}
/*--------------------------------------------------------------------------*/
inline void DEV_BVS::tr_unload()
{
  cout << "bvs unload\n";
  _m0.c0 = 0; _source0 = 0;
  _one0 = 0;

  tr_load();
}
/*--------------------------------------------------------------------------*/
void DEV_BVS::precalc_last()
{
  trace0("DEV_BVS::precalc_last()");
  ELEMENT::precalc_last();
  set_constant(!has_tr_eval());
  set_converged(!has_tr_eval());
}
/*--------------------------------------------------------------------------*/
void DEV_BVS::tr_begin()
{
  ELEMENT::tr_begin();
  // _y1.f0 = _y[0].f0 = 0.; // override
  _m0.x  = 0.;
  _m0.c0 = value() ; //  -_loss0 * _y[0].f1;
  _m0.c1 = 0.;
  trace2("DEV_BVS::tr_begin", long_label(), value());
  _source0 = value();
  _source1 = 0;
  _one0 = 1;
  _one1 = 0;
  //  q_load(); auch in d_vs nich
  //  q_eval(); auch in d_vs nich
  // _m1 = _m0;    
  if (!using_tr_eval()) {
  }else{
  }
}
/*--------------------------------------------------------------------------*/
bool DEV_BVS::do_tr()
{

  assert(_m0.x == 0.);
  if (using_tr_eval()) {
    _y[0].x = _sim->_time0;
    tr_eval();
    trace2("DEV_BVS::do_tr, did tr_eval", _y[0], _y[0].f1);
    if (_n[OUT2].m_() == 0) {
      _sim->set_limit(_y[0].f1);
    }else if (_n[OUT1].m_() == 0) {
      _sim->set_limit(-_y[0].f1);
    }else{
      //BUG// don't set limit
    }
    store_values();
    q_load();
    _m0.c0 = _y[0].f1;
    assert(_m0.c1 == 0.);
  }else{itested();
    // assert(false); //??
    assert(_y[0].x == 0.);
    assert(_y[0].f1 == value());
    assert(_m0.x == 0.);
    assert(_m0.c1 == 0.);
    assert(_y1 == _y[0]);
    assert(converged());
  }
  trace1( "DEV_BVS::do_tr value", value() );
  q_load();
  return converged();
}
/*--------------------------------------------------------------------------*/
void DEV_BVS::do_ac()
{
  if (using_ac_eval()) { itested();
    ac_eval();
    _acg = _ev;
  }else{itested();
    assert(_acg == 0.);
    //  untested();
    _ev = value();
    _acg = _ev;
  }
}
/*--------------------------------------------------------------------------*/
double DEV_BVS::tr_amps()const
{ untested();
  trace2( "DEV_BVS::tr_amps", value(), _loss0 );
  unsigned n1 = _n[BR].m_();
  return (CKT_BASE::_sim->_v0[n1]);
}
/*--------------------------------------------------------------------------*/
DEV_BVS p1;
#ifdef BVS_TO_MODULE_TRANSITION
DISPATCHER<CARD>::INSTALL d1(&device_dispatcher, "vsb|branchvs", &p1);
#else
// override voltage source
DISPATCHER<CARD>::INSTALL d1(&device_dispatcher, "V|vsource", &p1);
#endif
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:et
