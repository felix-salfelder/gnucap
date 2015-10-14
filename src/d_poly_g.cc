/*                               -*- C++ -*-
 * Copyright (C) 2001 Albert Davis
 *               2015 Felix Salfelder
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
 * DEV_CPOLY_G
 * number of nodes = 2*n_ports
 * number of val, ov = n_ports+1
 * val[0] is the constant part, val[1] is self admittance,
 *   val[2+] are transadmittances, up to n_ports
 * node[0] and node[1] are the output.
 * node[2] up are inputs.
 * node[2*i] and node[2*i+1] correspond to val[i+1]
 */
#include <vector>
#include "u_parameter.h"
/*--------------------------------------------------------------------------*/
#include "e_elemnt.h"
#include "m_mvpoly.h"
/*--------------------------------------------------------------------------*/
namespace {
/*--------------------------------------------------------------------------*/
class COMMON_G_POLY_K : public EVAL_BM_ACTION_BASE {
public:
  explicit COMMON_G_POLY_K(int x) :
    EVAL_BM_ACTION_BASE(x),
    _poly(NULL),
    _n_ports(0)
  {}
  COMMON_G_POLY_K(const COMMON_G_POLY_K& p) :
    EVAL_BM_ACTION_BASE(p),
    _coeffs(p._coeffs),
    _poly(NULL),
    _n_ports(p._n_ports)
  {}
  COMMON_COMPONENT* clone()const{return new COMMON_G_POLY_K(*this);}
  bool operator==(const COMMON_COMPONENT&x)const { untested();
    const COMMON_G_POLY_K* p = dynamic_cast<const COMMON_G_POLY_K*>(&x);
    bool rv = p
      && _n_ports == p->_n_ports
      && _coeffs == p->_coeffs
      && EVAL_BM_ACTION_BASE::operator==(x);
    return rv;
  }

  bool use_obsolete_callback_parse()const {return false;}
  bool use_obsolete_callback_print()const {return false;}
  bool has_parse_params_obsolete_callback()const {return false;}
  std::string name()const{untested(); return "poly_k";}
  bool has_tr_eval()const{return true;}
  bool ac_too()const {untested();return false;}
  void set_param_by_name(std::string Name, std::string Value)
  { untested();
    if (Umatch(Name, "c{oeffs} ")) { untested();
      _coeffs = Value;
    }else{
      EVAL_BM_ACTION_BASE::set_param_by_name(Name, Value);
    }
  }
  void set_param_by_index(int i, std::string& s, int)
  { untested();
//    unsigned i = COMMON_G_POLY_K::param_count() - 1 - I;
    if(i < int(_coeff.size())) {
      //ok.
    }else{
      _coeff.resize(i+1);
    }
    _coeff[i] = s;
  }
  void set_nports(unsigned n_ports)
  {
    _n_ports = n_ports;
  }

  void precalc_first(const CARD_LIST* scope)
  { untested();
    trace2("COMMON_G_POLY_K::precfirst", _coeffs, _n_ports);
    _coeffs.e_val(vector<PARAMETER<double> >(), scope);
    if(_coeffs.has_hard_value()){ untested();
    }else if(_coeff.size()==0){ untested();
    }else if(_coeff[0].has_hard_value()){ untested();
      _coeffs = _coeff;
    }
    if(_poly){ untested();
    }else{untested();
      _poly = new MV_POLY<double>(vector<PARAMETER<double> >(_coeffs), _n_ports-1);
    }
  }
  void precalc_last(const CARD_LIST* scope)
  { untested();
    _coeffs.e_val(vector<PARAMETER<double> >(), scope);
    if(_poly){ untested();
    }else{untested();
      _poly = new MV_POLY<double>(vector<PARAMETER<double> >(_coeffs), _n_ports-1);
    }
  }
  void tr_eval(ELEMENT* e) const;
private:
  PARAMETER<vector<PARAMETER<double> > > _coeffs;
  vector<PARAMETER<double> > _coeff; // one by one...
  MV_POLY<double>* _poly;
  unsigned _n_ports;
};
/*--------------------------------------------------------------------------*/
class DEV_CPOLY_G : public ELEMENT {
  friend class COMMON_G_POLY_K;
protected:
  hp_float_t*  _values;
  hp_float_t*  _old_values;
  uint_t       _n_ports;
  double   _time;
  const double** _inputs;
protected:
  explicit DEV_CPOLY_G(const DEV_CPOLY_G& p);
public:
  explicit DEV_CPOLY_G(COMMON_COMPONENT* c);
  explicit DEV_CPOLY_G();
  ~DEV_CPOLY_G();
protected: // override virtual
  char	   id_letter()const	{unreachable(); return '\0';}
  std::string value_name()const	{return "p0";}
//  std::string dev_type()const	{untested(); return "cpoly_g";}
  uint_t	   max_nodes()const	{return 99; /*incomplete*/}
  uint_t	   min_nodes()const	{return net_nodes();}
  uint_t	   matrix_nodes()const	{return _n_ports*2;}
  uint_t	   net_nodes()const	{untested(); return _net_nodes;}
  uint_t	   ext_nodes()const	{return net_nodes();}
  CARD*	   clone()const		{return new DEV_CPOLY_G(*this);}
  void	   tr_iwant_matrix()	{tr_iwant_matrix_extended();}
  bool	   do_tr();
  void	   tr_load();
  void	   tr_unload();
  hp_float_t   tr_involts()const	{unreachable(); return NOT_VALID;}
  hp_float_t   tr_involts_limited()const {unreachable(); return NOT_VALID;}
  hp_float_t   tr_amps()const;
  void	   ac_iwant_matrix()	{ac_iwant_matrix_extended();}
  void	   ac_load();
  COMPLEX  ac_involts()const	{itested(); return NOT_VALID;}
  COMPLEX  ac_amps()const	{itested(); return NOT_VALID;}

  std::string port_name(uint_t n)const {itested();
    if(n==0){
      return "p";
    }else if(n==1){
      return "n";
    }else if(n%2){
      return "cn" + to_string(n/2-1);
    }else{
      return "cp" + to_string(n/2-1);
    }
  }
public:
  void set_parameters(const std::string& Label, CARD* Parent,
		      COMMON_COMPONENT* Common, double Value,
		      uint_t state_count, hp_float_t state[],
		      uint_t node_count, const node_t nodes[]);
  //		      const double* inputs[]=0);
  void set_param_by_index(int, std::string&, int);
  void set_param_by_name(const std::string, const std::string);
  void expand();
  void precalc_last();
protected:
  bool do_tr_con_chk_and_q();
};
/*--------------------------------------------------------------------------*/
void COMMON_G_POLY_K::tr_eval(ELEMENT* e) const
{
  DEV_CPOLY_G* d = prechecked_cast<DEV_CPOLY_G*>(e);
  for(unsigned i=2; i<=_n_ports; ++i) { itested();
    d->_values[i] = dn_diff(d->n_(2*i-2).v0(),d->n_(2*i-1).v0());
  }
  trace2("b4", d->_values[2], d->_values[3]);
  d->_values[1] = _poly->eval(d->_values+2);
  trace1("evald", d->_values[1]);
  trace2("deriv", d->_values[2], d->_values[3]);

  for(unsigned i=2; i<=_n_ports; ++i) { itested();
    d->_values[1] -= dn_diff(d->n_(2*i-2).v0(),d->n_(2*i-1).v0()) * d->_values[i];
  }

  trace3("after", d->_values[0], d->_values[1], d->_values[2]);
}
/*--------------------------------------------------------------------------*/
#if 0
class DEV_FPOLY_G : public DEV_CPOLY_G {
private:
  explicit DEV_FPOLY_G(const DEV_FPOLY_G& p)
    :DEV_CPOLY_G(p) {incomplete(); unreachable();}
public:
  explicit DEV_FPOLY_G() :DEV_CPOLY_G() {untested();}
private: // override virtual
  char	   id_letter()const	{unreachable(); return '\0';}
  std::string dev_type()const	{unreachable(); return "fpoly_g";}
  CARD*	   clone()const		{unreachable(); return new DEV_FPOLY_G(*this);}
  bool	   do_tr();
};
#endif
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
DEV_CPOLY_G::DEV_CPOLY_G(const DEV_CPOLY_G& p)
  :ELEMENT(p),
   _values(NULL),
   _old_values(NULL),
   _n_ports(p._n_ports),
   _time(NOT_VALID),
   _inputs(NULL)
{
  // not really a copy .. only valid to copy a default
  // too lazy to do it right, and that's all that is being used
  // to do it correctly requires a deep copy
  // just filling in defaults is better than a shallow copy, hence this:
  assert(!p._values);
  assert(!p._old_values);
  assert(p._n_ports == 0);
  assert(!p._inputs);
}
/*--------------------------------------------------------------------------*/
DEV_CPOLY_G::DEV_CPOLY_G()
  :ELEMENT(),
   _values(NULL),
   _old_values(NULL),
   _n_ports(0),
   _time(NOT_VALID),
   _inputs(NULL)
{
}
/*--------------------------------------------------------------------------*/
DEV_CPOLY_G::DEV_CPOLY_G(COMMON_COMPONENT* c)
  :ELEMENT(),
   _values(NULL),
   _old_values(NULL),
   _n_ports(0),
   _time(NOT_VALID),
   _inputs(NULL)
{ untested();
  attach_common(c);
}
/*--------------------------------------------------------------------------*/
DEV_CPOLY_G::~DEV_CPOLY_G()
{
  delete [] _old_values;
  if (net_nodes() > NODES_PER_BRANCH) {
    delete [] _n;
  }else{
    // it is part of a base class
  }
}
/*--------------------------------------------------------------------------*/
bool DEV_CPOLY_G::do_tr_con_chk_and_q()
{
  q_load();

  assert(_old_values);
  set_converged(conchk(_time, _sim->_time0));
  _time = _sim->_time0;
  for (uint_t i=0; converged() && i<=_n_ports; ++i) {
    set_converged(conchk(_old_values[i], _values[i]));
  }
  return converged();
}
/*--------------------------------------------------------------------------*/
bool DEV_CPOLY_G::do_tr()
{
  if (using_tr_eval()) { untested();
    assert(dynamic_cast<COMMON_G_POLY_K const*>(common())); // for now.
    tr_eval();
  }else{ untested();
  }
  assert(_values);
  _m0 = CPOLY1(0., _values[0], _values[1]);
  return do_tr_con_chk_and_q();
}
/*--------------------------------------------------------------------------*/
#if 0
bool DEV_FPOLY_G::do_tr()
{untested();
  assert(_values);
  double c0 = _values[0];
  if (_inputs) {untested();
    untested();
    for (int i=1; i<=_n_ports; ++i) {untested();
      c0 -= *(_inputs[i]) * _values[i];
      trace4("", i, *(_inputs[i]), _values[i], *(_inputs[i]) * _values[i]);
    }
  }else{untested();
    for (int i=1; i<=_n_ports; ++i) {untested();
      c0 -= volts_limited(_n[2*i-2],_n[2*i-1]) * _values[i];
      trace4("", i, volts_limited(_n[2*i-2],_n[2*i-1]), _values[i],
	     volts_limited(_n[2*i-2],_n[2*i-1]) * _values[i]);
    }
  }
  trace2("", _values[0], c0);
  _m0 = CPOLY1(0., c0, _values[1]);

  return do_tr_con_chk_and_q();
}
#endif
/*--------------------------------------------------------------------------*/
void DEV_CPOLY_G::tr_load()
{
  tr_load_passive();
  _old_values[0] = _values[0];
  _old_values[1] = _values[1];
  for (uint_t i=2; i<=_n_ports; ++i) {
    tr_load_extended(_n[OUT1], _n[OUT2], _n[2*i-2], _n[2*i-1], &(_values[i]), &(_old_values[i]));
  }
}
/*--------------------------------------------------------------------------*/
void DEV_CPOLY_G::tr_unload()
{
  std::fill_n(_values, _n_ports+1, 0.);
  _m0.c0 = _m0.c1 = 0.;
  _sim->mark_inc_mode_bad();
  tr_load();
}
/*--------------------------------------------------------------------------*/
hp_float_t DEV_CPOLY_G::tr_amps()const
{
  hp_float_t amps = _m0.c0;
  for (uint_t i=1; i<=_n_ports; ++i) {
    amps += dn_diff(_n[2*i-2].v0(),_n[2*i-1].v0()) * _values[i];
  }
  return amps;
}
/*--------------------------------------------------------------------------*/
void DEV_CPOLY_G::ac_load()
{
  _acg = _values[1];
  ac_load_passive();
  for (uint_t i=2; i<=_n_ports; ++i) {
    ac_load_extended(_n[OUT1], _n[OUT2], _n[2*i-2], _n[2*i-1], _values[i]);
  }
}
/*--------------------------------------------------------------------------*/
/* set: set parameters, used in model building
 */
void DEV_CPOLY_G::set_parameters(const std::string& Label, CARD *Owner,
				 COMMON_COMPONENT *Common, double Value,
				 uint_t n_states, hp_float_t states[],
				 uint_t n_nodes, const node_t nodes[])
  //				 const double* inputs[])
{
  bool first_time = (net_nodes() == 0);

  set_label(Label);
  set_owner(Owner);
  set_value(Value);
  attach_common(Common);

  if (first_time) {
    _n_ports = n_nodes/2; // sets num_nodes() = _n_ports*2
    assert(_n_ports+1 == n_states);

    assert(!_old_values);
    _old_values = new hp_float_t[n_states];

    if (matrix_nodes() > NODES_PER_BRANCH) {
      trace1("DEV_CPOLY_G::set_parameters, more nodes", matrix_nodes());
      // allocate a bigger node list
      _n = new node_t[matrix_nodes()];
    }else{
      // use the default node list, already set
    }      
  }else{
    assert(_n_ports == n_states-1);
    assert(_old_values);
    assert(net_nodes() == n_nodes);
    // assert could fail if changing the number of nodes after a run
  }

  //_inputs = inputs;
  _inputs = 0;
  _values = states;
  std::fill_n(_values, n_states, 0.);
  std::fill_n(_old_values, n_states, 0.);
  notstd::copy_n(nodes, net_nodes(), _n); // copy to _n
  assert(net_nodes() == _n_ports * 2);
}
/*--------------------------------------------------------------------------*/
void DEV_CPOLY_G::set_param_by_index(int pos, std::string& Value, int slot)
{
  trace2("DEV_CPOLY_G::set_param_by_index", pos, Value);
  if (common()){ untested();
    COMMON_COMPONENT* m = common()->clone();
    try{ untested();
      unsigned count=0; // incomplete;
      m->set_param_by_index(pos+count, Value, slot);
    }catch(...){ untested();
      delete m;
      throw;
    }
    attach_common(m);

  }else{ incomplete();
  }
}
/*--------------------------------------------------------------------------*/
void DEV_CPOLY_G::set_param_by_name(std::string Name, std::string Value)
{
  bool retry = false;
  if (!common()){ incomplete();
  }else if (Umatch(Name, "p0 ")) {
    incomplete();
  }else if (Umatch(Name, "nports ")) {
    if(_n_ports){ untested();
      throw Exception("only one nports allowed right now.");
    }
    _n_ports = atoi(Value.c_str());
    if (matrix_nodes() > NODES_PER_BRANCH) { untested();
      // allocate a bigger node list
      trace3("DEV_CPOLY_G more nodes", long_label(), matrix_nodes(), Value);
      _n = new node_t[matrix_nodes()];
    }else{ untested();
      trace3("DEV_CPOLY_G enough nodes", long_label(), matrix_nodes(), NODES_PER_BRANCH);
      // use the default node list, already set
    }

    // HACK
    assert(common());
    COMMON_COMPONENT* m = common()->clone();
    assert(m);
    COMMON_G_POLY_K* p = dynamic_cast<COMMON_G_POLY_K*>(m);
    assert(p);
    p->set_nports(_n_ports);
    attach_common(p);
    _values = new double[_n_ports*2];
    _values[0] = 0;
    _old_values = new double[_n_ports*2];
  }else if (Umatch(Name, "c{oeffs} ")) { untested();
    assert(common());
    COMMON_COMPONENT* m = common()->clone();
    try{ untested();
      m->set_param_by_name(Name,Value);
    }catch(Exception){untested();
      retry = true;
      delete m;
      m = NULL;
    }
    if(m){ untested();
      attach_common(m);
    }else{ untested();
    }
  }else{
   retry = true;
  }

  if(retry){ untested();
    ELEMENT::set_param_by_name(Name, Value);
  }else{ untested();
  }
}
/*--------------------------------------------------------------------------*/
void DEV_CPOLY_G::expand()
{
  if(_values){ untested();
    // hack: done by set_parameters
    // move to common somehow
  }else{ untested();
  }
}
/*--------------------------------------------------------------------------*/
void DEV_CPOLY_G::precalc_last()
{
  trace3("DEV_CPOLY_G::precalc_last", long_label(), net_nodes(), matrix_nodes());

  if(_values){ untested();
    // hack: done by set_parameters
    // move to common somehow
  }else{ untested();
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// raw device. used by modelgen
DEV_CPOLY_G p4;
DISPATCHER<CARD>::INSTALL d4(&device_dispatcher, "cpoly_g", &p4);

COMMON_G_POLY_K poly(CC_STATIC);
DISPATCHER<COMMON_COMPONENT>::INSTALL d1(&bm_dispatcher, "poly_k", &poly);

// not required, maybe. using bm_wrapper
// DEV_CPOLY_G p5(&poly);
// DISPATCHER<CARD>::INSTALL d5(&device_dispatcher, "g_poly_k", &p5);
}
/*--------------------------------------------------------------------------*/
// hack: put it here, for now.
vector<vector<unsigned> > MV_POLY_BASE::_spn_list=vector<vector<unsigned> >(1,vector<unsigned>(1,1));
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
