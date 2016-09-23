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
 * inductors
 * y.x = amps,  y.f0 = flux,  ev = y.f1 = henrys
 * q = y history in time
 * i.x = amps,  i.f0 = volts,      i.f1 = ohms
 * m.x = volts, m.c0 = amps, acg = m.c1 = mhos
 */
#include "e_subckt.h"
#include "e_ccsrc.h"
#include "e_storag.h"
#include "d_coil.h"
/*--------------------------------------------------------------------------*/
using std::string;
/*--------------------------------------------------------------------------*/
namespace INDUCTANCE {
/*--------------------------------------------------------------------------*/
class DEV_MUTUAL_L : public DEV_INDUCTANCE { //
private:
  std::string	  _output_label;
  DEV_INDUCTANCE* _output;
  std::string	  _input_label;
  DEV_INDUCTANCE* _input;
  double _lm;
  double _mf0_c0;	// matrix parameters, new
  double _mf1_c0;	// matrix parameters, 1 fill ago
  double _mr0_c0;	// matrix parameters, new
  double _mr1_c0;	// matrix parameters, 1 fill ago
  FPOLY1 _yf1;		// iteration parameters, 1 iter ago
  FPOLY1 _yf[OPT::_keep_time_steps];
  FPOLY1 _if[OPT::_keep_time_steps];
  FPOLY1 _yr1;		// iteration parameters, 1 iter ago
  FPOLY1 _yr[OPT::_keep_time_steps];
  FPOLY1 _ir[OPT::_keep_time_steps];
private:
  explicit	DEV_MUTUAL_L(const DEV_MUTUAL_L& p);
public:
  explicit	DEV_MUTUAL_L();
private: // override virtual
  char	   id_letter()const	{return 'K';}
  bool	   print_type_in_spice()const {return false;}
  std::string value_name()const {return "k";}
  std::string dev_type()const	{untested(); return "mutual_inductor";}
  uint_t	   max_nodes()const	{return 2;}
  uint_t	   min_nodes()const	{return 2;}
  uint_t	   matrix_nodes()const	{return 2;}
  uint_t	   net_nodes()const	{return 0;}
  uint_t	   num_current_ports()const {return 2;}
  bool	   has_iv_probe()const  {untested(); return false;}
  bool	   use_obsolete_callback_parse()const {return false;}
  CARD*	   clone()const		{return new DEV_MUTUAL_L(*this);}
  void     expand_first();
  void	   expand_last();
  void	   precalc_last();
  void     tr_iwant_matrix()	{tr_iwant_matrix_passive();}
  void     tr_begin();
  void     dc_advance();
  void     tr_advance();
  bool     do_tr()		{_sim->_late_evalq.push_back(this); return true;}
  bool     do_tr_last();
  void	   tr_load();
  TIME_PAIR tr_review()		{return TIME_PAIR(NEVER,NEVER);}
  void	   tr_unload();
  double   tr_input()const		{return tr_involts();}
  double   tr_input_limited()const	{untested(); return tr_involts_limited();}
  double   tr_amps()const		{untested(); return _loss0 * tr_outvolts();}
  double   tr_probe_num(const std::string&)const;

  void	   ac_iwant_matrix()	{ac_iwant_matrix_passive();}
  void	   ac_load();
  COMPLEX  ac_amps()const	{untested(); return (double)_loss0 * ac_outvolts();}

  void	   set_port_by_name(std::string& Name, std::string& Value)
		{untested(); COMPONENT::set_port_by_name(Name,Value);}
  void	   set_port_by_index(uint_t Index, std::string& Value)
		{set_current_port_by_index(Index, Value);}
  bool	   node_is_connected(uint_t i)const {
    switch (i) {
    case 0:  return _output_label != "";
    case 1:  return _input_label != "";
    default: unreachable(); return false;
    }
  }

  std::string port_name(uint_t)const {untested();
    return "";
  }
  std::string current_port_name(uint_t i)const {untested();
    assert(i != INVALID_NODE);
    assert(i < 2);
    static std::string names[] = {"l1", "l2"};
    return names[i];
  }
  const std::string current_port_value(uint_t i)const {
    switch (i) {
    case 0:  return _output_label;
    case 1:  return _input_label;
    default: unreachable(); return COMPONENT::current_port_value(i);
    }
  }
  void set_current_port_by_index(uint_t i, const std::string& s) {
    switch (i) {
    case 0:  _output_label = s;	break;
    case 1:  _input_label = s;	break;
    default: unreachable();	break;
    }
  }
};
/*--------------------------------------------------------------------------*/
class DEV_BRANCH_L : public DEV_INDUCTANCE { //
public:
  DEV_BRANCH_L() :DEV_INDUCTANCE() { _c_model = true; }
  DEV_BRANCH_L(const DEV_BRANCH_L&p) :DEV_INDUCTANCE(p) { _c_model = true; }
  CARD*	   clone()const		{return new DEV_BRANCH_L(*this);}
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
DEV_MUTUAL_L::DEV_MUTUAL_L()
  :DEV_INDUCTANCE(),
   _output_label(),
   _output(0),
   _input_label(),
   _input(0),
   _lm(NOT_INPUT),
   _mf0_c0(0.),
   _mf1_c0(0.),
   _mr0_c0(0.),
   _mr1_c0(0.)
{
  _c_model = true;
  assert(_yf[0].x == 0. && _yf[0].f0 == 0. && _yf[0].f1 == 0.);
  assert(_yf1 == _yf[0]);
  assert(_yr[0].x == 0. && _yr[0].f0 == 0. && _yr[0].f1 == 0.);
  assert(_yr1 == _yr[0]);
}
/*--------------------------------------------------------------------------*/
DEV_MUTUAL_L::DEV_MUTUAL_L(const DEV_MUTUAL_L& p)
  :DEV_INDUCTANCE(p),
   _output_label(p._output_label),
   _output(p._output),
   _input_label(p._input_label),
   _input(p._input),
   _lm(p._lm),
   _mf0_c0(0.),
   _mf1_c0(0.),
   _mr0_c0(0.),
   _mr1_c0(0.)
{
  _c_model = true;
  assert(_yf[0].x == 0. && _yf[0].f0 == 0. && _yf[0].f1 == 0.);
  assert(_yf1 == _yf[0]);
  assert(_yr[0].x == 0. && _yr[0].f0 == 0. && _yr[0].f1 == 0.);
  assert(_yr1 == _yr[0]);
}
/*--------------------------------------------------------------------------*/
void DEV_MUTUAL_L::expand_first()
{
  _output = dynamic_cast<DEV_INDUCTANCE*>(find_in_my_scope(_output_label));
  if (!_output) {
    throw Exception_Type_Mismatch(long_label(), _output_label, "inductor");
  }else{
    _output->_c_model = true;
  }

  _input = dynamic_cast<DEV_INDUCTANCE*>(find_in_my_scope(_input_label));
  if (!_input) {
    throw Exception_Type_Mismatch(long_label(), _input_label, "inductor");
  }else{
    _input->_c_model = true;
  }
}
/*--------------------------------------------------------------------------*/
void DEV_MUTUAL_L::expand_last()
{
  STORAGE::expand(); // skip DEV_INDUCTANCE
  if (_sim->is_first_expand()) {
    _n[OUT2] = _input->n_(IN1);
    _n[OUT1] = _output->n_(IN1);
  }else{untested();
  }
}
/*--------------------------------------------------------------------------*/
void DEV_MUTUAL_L::precalc_last()
{
  _output->precalc_last();
  _input->precalc_last();

  DEV_INDUCTANCE::precalc_last();

  double l1 = _output->value();
  double l2 = _input->value();
  _lm = value() * sqrt(l1 * l2);
  trace3(long_label().c_str(), l1, l2, _lm);

  if (_sim->has_op() == s_NONE) {
    assert(_y[0].x  == 0.);
    assert(_y[0].f0 == LINEAR);
    _y[0].f1 = -_lm; // override
    _yf[0] = _yr[0] = _y[0];
  }else{
  }
}
/*--------------------------------------------------------------------------*/
void DEV_MUTUAL_L::tr_begin()
{
  DEV_INDUCTANCE::tr_begin();
  assert(_y[0].x  == 0.);
  assert(_y[0].f0 == LINEAR);
  _y[0].f1 = -_lm; // override
  _y1 = _y[0];

  for (int i = 0;  i < OPT::_keep_time_steps;  ++i) {
    _if[i] = _ir[i] = FPOLY1(0., 0., 0.);
  }
  _mf1_c0 = _mf0_c0 = _mr1_c0 = _mr0_c0 = 0.;
}
/*--------------------------------------------------------------------------*/
void DEV_MUTUAL_L::dc_advance()
{
  STORAGE::dc_advance();

  for (int i = 1;  i < OPT::_keep_time_steps;  ++i) {
    _if[i] = _if[0];
    _ir[i] = _ir[0];
  }
}
/*--------------------------------------------------------------------------*/
void DEV_MUTUAL_L::tr_advance()
{
  STORAGE::tr_advance();
  
  for (int i=OPT::_keep_time_steps-1; i>0; --i) {
    _yf[i] = _yf[i-1];
    _yr[i] = _yr[i-1];
    _if[i] = _if[i-1];
    _ir[i] = _ir[i-1];
  }
}
/*--------------------------------------------------------------------------*/
bool DEV_MUTUAL_L::do_tr_last()
{
  double l1 = _output->_y[0].f1;
  double l2 = _input->_y[0].f1;
  _lm = value() * sqrt(l1 * l2);

  _y[0].x = _n[OUT1].v0() - _n[OUT2].v0(); // really current
  _y[0].f1 = -_lm;
  _y[0].f0 = _y[0].x * _y[0].f1;	   // flux = I * L
  trace3("", _y[0].x, _y[0].f0, _y[0].f1);
  store_values();
  _i[0] = differentiate(_y, _i, _time, _method_a);  // really voltage, v = df/dt
  trace3("", _i[0].x, _i[0].f0, _i[0].f1);
  _m0.x = NOT_VALID;
  _m0.c1 = -_loss0 * _loss0 * _i[0].c1();
  _m0.c0 = -_loss0 * _loss0 * _i[0].c0();
  trace3("", _m0.x, _m0.c0, _m0.c1);

  _yf[0].x = _n[OUT1].v0();
  _yf[0].f1 = -_lm;
  _yf[0].f0 = _yf[0].x * _yf[0].f1;
  trace3("", _yf[0].x, _yf[0].f0, _yf[0].f1);
  assert(_yf[0]==_yf[0]); // store_values();
  _yf1=_yf[0];		  // store_values();
  _if[0] = differentiate(_yf, _if, _time, _method_a);
  trace3("", _if[0].x, _if[0].f0, _if[0].f1);
  _mf0_c0 = -_loss0 * _loss0 * _if[0].c0();
  
  _yr[0].x = _n[OUT2].v0();
  _yr[0].f1 = -_lm;
  _yr[0].f0 = _yr[0].x * _yr[0].f1;
  trace3("", _yr[0].x, _yr[0].f0, _yr[0].f1);
  assert(_yr[0]==_yr[0]); // store_values();
  _yr1=_yr[0];		  // store_values();
  _ir[0] = differentiate(_yr, _ir, _time, _method_a);
  trace3("", _ir[0].x, _ir[0].f0, _ir[0].f1);
  _mr0_c0 = -_loss0 * _loss0 * _ir[0].c0();

  q_load();
  return true;
}
/*--------------------------------------------------------------------------*/
void DEV_MUTUAL_L::tr_load()
{
  tr_load_couple();
  tr_load_source();
  tr_load_source_point(_n[OUT2], &_mr0_c0, &_mr1_c0);
  tr_load_source_point(_n[OUT1], &_mf0_c0, &_mf1_c0);
}
/*--------------------------------------------------------------------------*/
void DEV_MUTUAL_L::tr_unload()
{untested();
  tr_unload_couple();
}
/*--------------------------------------------------------------------------*/
void DEV_MUTUAL_L::ac_load()
{
  ac_load_couple();
}
/*--------------------------------------------------------------------------*/
double DEV_MUTUAL_L::tr_probe_num(const std::string& x)const
{ untested();
  if (Umatch(x, "fflux ")) { untested();
    return _yf[0].f0;
  }else if (Umatch(x, "rflux ")) { untested();
    return _yr[0].f0;
  }else if (Umatch(x, "fiof{fset} ")) { untested();
    return _mf0_c0;
  }else if (Umatch(x, "riof{fset} ")) { untested();
    return _mr0_c0;
  }else{ untested();
    return DEV_INDUCTANCE::tr_probe_num(x);
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
DEV_MUTUAL_L   p1;
DEV_INDUCTANCE p2;
DEV_BRANCH_L p3;
DISPATCHER<CARD>::INSTALL
  d1(&device_dispatcher, "K|mutual_inductor", &p1),
  d2(&device_dispatcher, "L|inductor", &p2),
  d3(&device_dispatcher, "brl|branchl", &p3);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
