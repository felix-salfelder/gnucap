/*$Id: d_switch.cc,v 1.6 2010-07-09 12:14:22 felix Exp $ -*- C++ -*-
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
 * voltage (and current) controlled switch.
 * netlist syntax:
 * device:  Sxxxx n+ n- vc+ vc- mname <on>|<off> <model-card-args>
 * model:   .model mname SW <args>
 * current controlled switch
 * device:  Wxxxx n+ n- controlelement mname <on>|<off> <model-card-args>
 * model:   .model mname CSW <args>
 */
//testing=script 2006.06.14
#include "e_model.h"
#include "e_elemnt.h"
/*--------------------------------------------------------------------------*/
namespace {
/*--------------------------------------------------------------------------*/
using std::string;
/*--------------------------------------------------------------------------*/
  enum state_t {_UNKNOWN=0, _ON=1, _OFF=2};
/*--------------------------------------------------------------------------*/
class COMMON_SWITCH : public COMMON_COMPONENT {
private:
  explicit COMMON_SWITCH(const COMMON_SWITCH& p)
    :COMMON_COMPONENT(p), ic(p.ic), _ic(p._ic) {}
public:
  explicit COMMON_SWITCH(int c=0) 
    :COMMON_COMPONENT(c), ic(0), _ic(_UNKNOWN) {}
  bool operator==(const COMMON_COMPONENT&)const;
  COMMON_COMPONENT* clone()const {return new COMMON_SWITCH(*this);}
  std::string name()const	 {untested(); return "switch";}

  bool		param_is_printable(int)const;
  std::string	param_name(int)const;
  void		set_param_by_name(string, string);
  std::string	param_name(int,int)const;
  std::string	param_value(int)const;
  int param_count()const {return (1 + COMMON_COMPONENT::param_count());}
public:
  void precalc_last(const CARD_LIST* par_scope);
public:
  PARAMETER<int> ic;
public: // calculated params
  state_t	_ic;		/* initial state, belongs in common */
};
/*--------------------------------------------------------------------------*/
class SWITCH_BASE : public ELEMENT {
protected:
  explicit	SWITCH_BASE();
  explicit	SWITCH_BASE(const SWITCH_BASE& p);
protected: // override virtual
  std::string value_name()const	{return "";}
  std::string dev_type()const {assert(common()); return common()->modelname().c_str();}
  bool	   print_type_in_spice()const {return true;}
  uint_t	   tail_size()const	{return 1;}
  uint_t	   max_nodes()const	= 0;
  uint_t	   min_nodes()const	= 0;
  uint_t	   matrix_nodes()const	{return 2;}
  uint_t	   net_nodes()const	= 0;
  CARD*	   clone()const		= 0;
  void     expand();
  void     precalc_last();
  void	   tr_iwant_matrix()	{tr_iwant_matrix_passive();}
  void	   tr_begin();
  void     dc_advance();
  void     tr_advance();
  void     tr_regress();
  bool     tr_needs_eval()const {return _sim->analysis_is_static();} // also q by tr_advance
  bool	   do_tr();
  void	   tr_load()		{tr_load_passive();}
  TIME_PAIR tr_review();
  void	   tr_unload()		{untested(); tr_unload_passive();}
  hp_float_t   tr_involts()const {return tr_outvolts();}
  hp_float_t   tr_involts_limited()const {unreachable(); return tr_outvolts_limited();}
  void	   ac_iwant_matrix()	{ac_iwant_matrix_passive();}
  void	   ac_begin()		{_ev = _y[0].f1; _acg = _m0.c1;}
  void	   do_ac();
  void	   ac_load()		{ac_load_passive();}
  COMPLEX  ac_involts()const	{untested(); return ac_outvolts();}
  double tr_probe_num(const std::string& x)const;
  void     do_tt() { q_tt_accept(); }
  void     tt_advance();
  void     tt_regress();
  TIME_PAIR tt_review();
  void     tt_accept();
protected:
  const ELEMENT* _input;
private:
  double	_in[OPT::_keep_time_steps];
  double	_in_tt[2];
  state_t	_state[OPT::_keep_time_steps];
  state_t	_state_tt[2];
};
/*--------------------------------------------------------------------------*/
class DEV_VSWITCH : public SWITCH_BASE {
private:
  explicit  DEV_VSWITCH(const DEV_VSWITCH& p) :SWITCH_BASE(p) {}
public:
  explicit  DEV_VSWITCH()	:SWITCH_BASE() {}
private: // override virtual
  uint_t	    max_nodes()const	{return 4;}
  uint_t	    min_nodes()const	{return 4;}
  uint_t	    net_nodes()const	{return 4;}
  CARD*	    clone()const	{return new DEV_VSWITCH(*this);}
  char	    id_letter()const	{return 'S';}

  std::string port_name(uint_t i)const {itested();
    assert(i != INVALID_NODE);
    assert(i < 4);
    static std::string names[] = {"p", "n", "ps", "ns"};
    return names[i];
  }
};
/*--------------------------------------------------------------------------*/
class DEV_CSWITCH : public SWITCH_BASE {
private:
  explicit  DEV_CSWITCH(const DEV_CSWITCH& p) 
    :SWITCH_BASE(p), _input_label(p._input_label) {}
public:
  explicit  DEV_CSWITCH()	:SWITCH_BASE(), _input_label() {}
private: // override virtual
  uint_t	    max_nodes()const	{return 3;}
  uint_t	    ext_nodes()const	{return 2;}
  uint_t	    min_nodes()const	{return 3;}
  uint_t	    net_nodes()const	{return 2;}
  uint_t	    num_current_ports()const {return 1;}
  const std::string current_port_value(uint_t)const {return _input_label;};
  CARD*	    clone()const	{return new DEV_CSWITCH(*this);}
  void	    expand();
  char	    id_letter()const	{return 'W';}
  void	   set_port_by_name(std::string& Name, std::string& Value)
		{untested(); SWITCH_BASE::set_port_by_name(Name,Value);}
  void	   set_port_by_index(uint_t, std::string&);
  bool	   node_is_connected(uint_t)const;

  std::string port_name(uint_t i)const {itested();
    assert(i != INVALID_NODE);
    assert(i < 2);
    static std::string names[] = {"p", "n"};
    return names[i];
  }
  std::string current_port_name(uint_t i)const {
    assert(i != INVALID_NODE);
    assert(i < 1);
    static std::string names[] = {"in"};
    return names[i];
  }
private:
  std::string	 _input_label;
};
/*--------------------------------------------------------------------------*/
class MODEL_SWITCH : public MODEL_CARD {
private:
  explicit	MODEL_SWITCH(const MODEL_SWITCH& p);
public:
  explicit	MODEL_SWITCH(const SWITCH_BASE*);
private: // override virtual
  void		set_dev_type(const std::string& nt);
  std::string	dev_type()const;
  CARD*		clone()const	{return new MODEL_SWITCH(*this);}
  void		precalc_first();
//  void		set_param_by_name(string, string); // not yet
  void		set_param_by_index(int, std::string&, int);
  bool		param_is_printable(int)const;
  std::string	param_name(int)const;
  std::string	param_name(int,int)const;
  std::string	param_value(int)const;
  int param_count()const {return (6 + MODEL_CARD::param_count());}
public:
  PARAMETER<double> vt;		/* threshold voltage */
  PARAMETER<double> vh;		/* hysteresis voltage */
  PARAMETER<double> ron;	/* on resistance */
  PARAMETER<double> roff;	/* off resistance */
  PARAMETER<double> von;
  PARAMETER<double> voff;
  enum control_t {VOLTAGE, CURRENT};
  control_t type;	/* current or voltage controlled */
private:
  static double const _default_vt;
  static double const _default_vh;
  static double const _default_ron;
  static double const _default_roff;
};
/*--------------------------------------------------------------------------*/
static COMMON_SWITCH Default_SWITCH(CC_STATIC);
double const MODEL_SWITCH::_default_vt = 0.;
double const MODEL_SWITCH::_default_vh = 0.;
double const MODEL_SWITCH::_default_ron = 1.;
double const MODEL_SWITCH::_default_roff = 1e12;
/*--------------------------------------------------------------------------*/
bool COMMON_SWITCH::operator==(const COMMON_COMPONENT& x)const
{
  const COMMON_SWITCH* p = dynamic_cast<const COMMON_SWITCH*>(&x);
  return p
    && ic == p->ic
    && COMMON_COMPONENT::operator==(x);
}
/*--------------------------------------------------------------------------*/
bool COMMON_SWITCH::param_is_printable(int i)const
{
  switch (COMMON_SWITCH::param_count() - 1 - i) {
  case 0:  return (_ic == _ON || _ic == _OFF);
  default: return COMMON_COMPONENT::param_is_printable(i);
  }
}
/*--------------------------------------------------------------------------*/
void COMMON_SWITCH::set_param_by_name(string Name, string Value)
{
  if (Umatch (Name,"ic")) {
    ic = Value;
    trace4("COMMON_SWITCH::set_param_by_name", Name, Value, hp(this), ic);
  }else{ untested();
    return COMMON_COMPONENT::set_param_by_name(Name, Value);
  }
}
/*--------------------------------------------------------------------------*/
std::string COMMON_SWITCH::param_name(int i)const
{
  switch (COMMON_SWITCH::param_count() - 1 - i) {
  case 0:  return "ic";
  default: return COMMON_COMPONENT::param_name(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string COMMON_SWITCH::param_name(int i, int j)const
{itested();
  if (j == 0) {itested();
    return param_name(i);
  }else if (i >= COMMON_COMPONENT::param_count()) {itested();
    return "";
  }else{itested();
    return COMMON_COMPONENT::param_name(i, j);
  }
}
/*--------------------------------------------------------------------------*/
std::string COMMON_SWITCH::param_value(int i)const
{
  switch (COMMON_SWITCH::param_count() - 1 - i) {
  case 0:  return (_ic == _ON) ? "1" : "0";
  default: return COMMON_COMPONENT::param_value(i);
  }
}
/*--------------------------------------------------------------------------*/
void COMMON_SWITCH::precalc_last(const CARD_LIST* par_scope)
{
  assert(par_scope);
  COMMON_COMPONENT::precalc_last(par_scope);
  ic.e_val(0, par_scope);

  if (ic.has_hard_value()){
    _ic = (ic)?_ON:_OFF;
    trace2("COMMON_SWITCH::precalc_last", _ic, ic);
  }else{
    trace3("COMMON_SWITCH::precalc_last no hard value", _ic, ic, hp(this));
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
MODEL_SWITCH::MODEL_SWITCH(const SWITCH_BASE* p)
  :MODEL_CARD(p),
   vt(_default_vt),
   vh(_default_vh),
   ron(_default_ron),
   roff(_default_roff),
   von(_default_vt + _default_vh),
   voff(_default_vt - _default_vh),
   type(VOLTAGE)
{
}
/*--------------------------------------------------------------------------*/
MODEL_SWITCH::MODEL_SWITCH(const MODEL_SWITCH& p)
  :MODEL_CARD(p),
   vt(p.vt),
   vh(p.vh),
   ron(p.ron),
   roff(p.roff),
   von(p.von),
   voff(p.voff),
   type(p.type)
{
}
/*--------------------------------------------------------------------------*/
void MODEL_SWITCH::set_dev_type(const std::string& new_type)
{
  if (new_type == "sw") {
    type = VOLTAGE;
  }else if (new_type == "csw") {
    type = CURRENT;
  }else{
    MODEL_CARD::set_dev_type(new_type);
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_SWITCH::precalc_first()
{
  MODEL_CARD::precalc_first();

  const CARD_LIST* par_scope = scope();
  assert(par_scope);

  vt.e_val(_default_vt, par_scope);
  vh.e_val(_default_vh, par_scope);
  ron.e_val(_default_ron, par_scope);
  roff.e_val(_default_roff, par_scope);
  von.e_val(vt + vh, par_scope);
  voff.e_val(vt - vh, par_scope);
}
/*--------------------------------------------------------------------------*/
std::string MODEL_SWITCH::dev_type()const
{
  switch (type) {
  case VOLTAGE:
    return "sw";
    break;
  case CURRENT:
    return "csw";
    break;
  }
  unreachable();
  return "";
}
/*--------------------------------------------------------------------------*/
void MODEL_SWITCH::set_param_by_index(int i, std::string& value, int offset)
{
  switch (MODEL_SWITCH::param_count() - 1 - i) {
  case 0: vt = value; break;
  case 1: vh = value; break;
  case 2: von = value; break;
  case 3: voff = value; break;
  case 4: ron = value; break;
  case 5: roff = value; break;
  default: MODEL_CARD::set_param_by_index(i, value, offset); break;
  }
}
/*--------------------------------------------------------------------------*/
bool MODEL_SWITCH::param_is_printable(int i)const
{
  switch (MODEL_SWITCH::param_count() - 1 - i) {
  case 0: 
  case 1: 
  case 2: 
  case 3: 
  case 4: 
  case 5: return true;
  default: return MODEL_CARD::param_is_printable(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_SWITCH::param_name(int i)const
{
  string p = (type==VOLTAGE)?"v":"i";
  switch (MODEL_SWITCH::param_count() - 1 - i) {
    case 0: return p+"t";
    case 1: return p+"h";
    case 2: return p+"on";
    case 3: return p+"off";
    case 4: return "ron";
    case 5: return "roff";
    default: return MODEL_CARD::param_name(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_SWITCH::param_name(int i, int j)const
{
  if (j == 0) {
    return param_name(i);
  }else if (i >= MODEL_CARD::param_count()) {
    return "";
  }else{
    return MODEL_CARD::param_name(i, j);
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_SWITCH::param_value(int i)const
{
  switch (MODEL_SWITCH::param_count() - 1 - i) {
  case 0: return vt.string();
  case 1: return vh.string();
  case 2: return von.string();
  case 3: return voff.string();
  case 4: return ron.string();
  case 5: return roff.string();
  default: return MODEL_CARD::param_value(i);
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
SWITCH_BASE::SWITCH_BASE()
  :ELEMENT(),
   _input(NULL)
{
  attach_common(&Default_SWITCH);
  std::fill_n(_in, int(OPT::_keep_time_steps), 0.);
  std::fill_n(_state, int(OPT::_keep_time_steps), _UNKNOWN);
}
/*--------------------------------------------------------------------------*/
SWITCH_BASE::SWITCH_BASE(const SWITCH_BASE& p)
  :ELEMENT(p),
   _input(NULL)
{
  notstd::copy_n(p._in, int(OPT::_keep_time_steps), _in);
  notstd::copy_n(p._state, int(OPT::_keep_time_steps), _state);  
}
/*--------------------------------------------------------------------------*/
void SWITCH_BASE::expand()
{
  ELEMENT::expand();
  
  attach_model();
  const COMMON_SWITCH* c = dynamic_cast<const COMMON_SWITCH*>(common());
  assert(c);
  const MODEL_SWITCH* m = dynamic_cast<const MODEL_SWITCH*>(c->model());
  if (!m) {
    assert(has_common());
    throw Exception_Model_Type_Mismatch(long_label(), common()->modelname(),
					"switch (SW or CSW)");
  }else{
  }
}
/*--------------------------------------------------------------------------*/
void SWITCH_BASE::precalc_last()
{
  ELEMENT::precalc_last();
    
  if (_sim->is_first_expand()) {
    const COMMON_SWITCH* c = prechecked_cast<const COMMON_SWITCH*>(common());
    assert(c);
    const MODEL_SWITCH* m = prechecked_cast<const MODEL_SWITCH*>(c->model());
    assert(m);
    trace1("SWITCH_BASE::precalc_last", c->_ic);
    _y1.f1 = _y[0].f1 = (c->_ic == _ON) ? m->ron : m->roff;	// override, unknown is off
    
    assert(!is_constant()); // depends on input
    // converged?????
    
    _m0.c1 = 1./_y[0].f1;
    _m0.c0 = 0.;
    _m1 = _m0;
    _state[1] = _state[0] = c->_ic;
  }else{
  }

  assert(_loss0 == 0.);
  assert(_loss1 == 0.);
}
/*--------------------------------------------------------------------------*/
void SWITCH_BASE::tr_begin()
{
  ELEMENT::tr_begin();

  const COMMON_SWITCH* c = prechecked_cast<const COMMON_SWITCH*>(common());
  assert(c);
  const MODEL_SWITCH* m = prechecked_cast<const MODEL_SWITCH*>(c->model());
  assert(m);

  assert(_loss0 == 0.);
  assert(_loss1 == 0.);
  assert(_y[0].f0 == LINEAR);
  _y1.f1 = _y[0].f1 = ((c->_ic == _ON) ? m->ron : m->roff);  /* unknown is off */
  _m0.c1 = 1./_y[0].f1;
  assert(_m0.c0 == 0.);
  _m1 = _m0;
  if (c->_ic!=_UNKNOWN) {
//    _in[0] = m->vt + (3-2*c->_ic) * m->vh;
    q_load();
  }

  _state[1] = _state[0] = c->_ic;
  set_converged();
}
/*--------------------------------------------------------------------------*/
void SWITCH_BASE::dc_advance()
{
  ELEMENT::dc_advance();
  _state[1] = _state[0];
}
/*--------------------------------------------------------------------------*/
void SWITCH_BASE::tr_advance()
{
  ELEMENT::tr_advance();

  const COMMON_SWITCH* c = prechecked_cast<const COMMON_SWITCH*>(common());
  assert(c);
  const MODEL_SWITCH* m = prechecked_cast<const MODEL_SWITCH*>(c->model());
  assert(m);
  
  _state[1] = _state[0];
  _y[0].x = _in[1] = _in[0];
  
  if (_y[0].x >= m->von) {
    _state[0] = _ON;
  }else if (_y[0].x <= m->voff) {
    _state[0] = _OFF;
  }else{
    _state[0] = _state[1];
  }
  
  if (_state[1] != _state[0]) {
    _y[0].f1 = (_state[0] == _ON) ? m->ron : m->roff;	/* unknown is off */
    _m0.c1 = 1./_y[0].f1;
    q_eval();
  }else{
  }
  
  assert(_y[0].f1 == ((_state[0] == _ON) ? m->ron : m->roff));
  assert(_y[0].f0 == LINEAR);
  assert(_m0.c0 == 0.);
  set_converged();
}
/*--------------------------------------------------------------------------*/
void SWITCH_BASE::tr_regress()
{
  ELEMENT::tr_regress();

  const COMMON_SWITCH* c = prechecked_cast<const COMMON_SWITCH*>(common());
  assert(c);
  const MODEL_SWITCH* m = prechecked_cast<const MODEL_SWITCH*>(c->model());
  assert(m); USE(m);

  // heuristic hack
  if (OPT::disc>=dREJECT && _y[0].f1 <= 1){
    for (unsigned i=0;i < 2; ++i) {
      if (DISCONT d=_n[i]->discont()) {
        _n[1-i]->discont(d);
      } else {
      }
    }
  }

  assert(_y[0].f1 == ((_state[0] == _ON) ? m->ron : m->roff));
  assert(_y[0].f0 == LINEAR);
  // assert(_m0.c1 == 1./_y[0].f1); see do_tr (?)
  assert(_m0.c0 == 0.);
  set_converged();
}
/*--------------------------------------------------------------------------*/
TIME_PAIR SWITCH_BASE::tt_review()
{
  const COMMON_SWITCH* c = prechecked_cast<const COMMON_SWITCH*>(common());
  assert(c);
  const MODEL_SWITCH* m = prechecked_cast<const MODEL_SWITCH*>(c->model());
  assert(m);

  double ttin = _in_tt[0]; // pre tr
  _in_tt[0] = _in[0];
  state_t ttstate = _state_tt[0];
  _state_tt[0] = _state[0];

  _time_by.reset();

  double old_dt = _sim->_dT0;
  double old_dv = _in_tt[0] - _in_tt[1];
  double new_dt = NEVER;
  assert(_sim->last_time()>0);

  if (ttstate != _state[0]) {
     // back... catch exact switching time.
     _time_by.min_event(_sim->_Time0 - _sim->_dT0/2);
  }else if (_state[0] != _ON && ttstate != _ON && old_dv > 0) {
    double new_dv = m->von - ttin;
    new_dt = old_dt * new_dv / old_dv - _sim->_dT0;
    _time_by.min_event(_sim->_Time0 + new_dt + _sim->last_time() );
  }else if (_state[0] != _OFF && ttstate != _OFF && old_dv < 0) {
    double new_dv = m->voff - ttin;
    new_dt = old_dt * new_dv / old_dv;
    _time_by.min_event(_sim->_Time0 + new_dt + _sim->last_time() );
  }else{
    assert(_time_by._event == NEVER);
  }

  return _time_by;
}
/*--------------------------------------------------------------------------*/
void SWITCH_BASE::tt_accept()
{
  trace3("SWITCH_BASE::tt_acc", _sim->_Time0, _in[0], _in[1]);
  ELEMENT::tt_accept();
  return;
}
/*--------------------------------------------------------------------------*/
void SWITCH_BASE::tt_advance()
{
  ELEMENT::tt_advance();
  _state_tt[1] = _state_tt[0];
  _in_tt[1] = _in_tt[0];

  // _in_tt[0] == ??!
}
/*--------------------------------------------------------------------------*/
void SWITCH_BASE::tt_regress()
{
  ELEMENT::tt_regress();
  const COMMON_SWITCH* c = prechecked_cast<const COMMON_SWITCH*>(common());
  assert(c);
  const MODEL_SWITCH* m = prechecked_cast<const MODEL_SWITCH*>(c->model());
  assert(m);
  trace3("tt_regr", _sim->_Time0, _in_tt[0], _in_tt[1]);

  _in[0] = _in[1] = _in_tt[1] = _in_tt[1];
  _state[0] = _state[1] = _state_tt[0] = _state_tt[1];

  // necessary? tr_advance should take care(?);
  _y[0].f1 = (_state[0] == _ON) ? m->ron : m->roff;	/* unknown is off */
  if (_state[0] == _UNKNOWN) { untested();
  }

  _m0.c1 = 1./_y[0].f1;
  q_eval();
}
/*--------------------------------------------------------------------------*/
bool SWITCH_BASE::do_tr()
{
  const COMMON_SWITCH* c = prechecked_cast<const COMMON_SWITCH*>(common());
  assert(c);
  const MODEL_SWITCH* m = prechecked_cast<const MODEL_SWITCH*>(c->model());
  assert(m);

  if (_sim->analysis_is_static() ) {
    _y[0].x = (_input)			/* _y[0].x is controlling value */
      ? CARD::probe(_input,"I")		/* current controlled */
      : _n[IN1].v0() - _n[IN2].v0();	/* voltage controlled */
    
    state_t new_state;
    if (_sim->uic_now() && c->_ic!=_UNKNOWN) { // necessary?
      new_state = _state[1] = c->_ic;
      _state[0] = _state[1];
    } else if (_y[0].x > m->von) {
      new_state = _ON;
    }else if (_y[0].x < m->voff) {
      new_state = _OFF;
    }else{
      new_state = _state[1];
    }
    
    if (new_state != _state[0]) {
      _y[0].f1 = (new_state == _ON) ? m->ron : m->roff;	/* unknown is off */
      _state[0] = new_state;
      _m0.c1 = 1./_y[0].f1;
      //trace4("change", new_state, old_state, _y[0].f1, _m0.c1);
      q_load();
      store_values();
      set_not_converged();
    }else{
      trace3("no change", new_state, _y[0].f1, _m0.c1);
      set_converged();
    }
  }else{
    // Above isn't necessary because it was done in tr_advance,
    // and doesn't iterate.
    // I believe it is not necessary on restart, because it is stored.
    if (_state[0] != _state[1]) {
      q_load();
      store_values();
    }else{
      // gets here only on "nobypass"
    }
    assert(converged());
  }
  assert(_y[0].f1 == ((_state[0] == _ON) ? m->ron : m->roff));
  assert(_y[0].f0 == LINEAR);
  //assert(_m0.c1 == 1./_y[0].f1);
  assert(_m0.c0 == 0.);
  return converged();
}
/*--------------------------------------------------------------------------*/
TIME_PAIR SWITCH_BASE::tr_review()
{
  ELEMENT::tr_review();

  const COMMON_SWITCH* c = prechecked_cast<const COMMON_SWITCH*>(common());
  assert(c);
  const MODEL_SWITCH* m = prechecked_cast<const MODEL_SWITCH*>(c->model());
  assert(m);
  
  _in[0] = (_input)
    ? CARD::probe(_input,"I")
    : _n[IN1].v0() - _n[IN2].v0();

  if (_sim->uic_now() && c->_ic!=_UNKNOWN) {
    _in[0] = m->vt + (3-2*c->_ic) * m->vh;
  }

  double old_dt = _time[0] - _time[1];
  double old_dv = _in[0] - _in[1];
  double new_dt = NEVER;

  if (_state[0] != _ON  &&  old_dv > 0) {
    double new_dv = m->von - _in[1];
    new_dt = old_dt * new_dv / old_dv;
    _time_by.min_event(_time[1] + new_dt);
  }else if (_state[0] != _OFF  &&  old_dv < 0) {
    double new_dv = m->voff - _in[1];
    new_dt = old_dt * new_dv / old_dv;
    _time_by.min_event(_time[1] + new_dt);
  }else{
    assert(_time_by._event == NEVER);
  }
  // _time_by_event is the predicted switch time

  if (new_dt<_dt+ _sim->_dtmin) {
    _discont |= disSECOND;
    q_accept();
  }

  return _time_by;
}
/*--------------------------------------------------------------------------*/
void SWITCH_BASE::do_ac()
{
  assert(_ev  == _y[0].f1);
  assert(_acg == (double)_m0.c1);
}
/*--------------------------------------------------------------------------*/
double SWITCH_BASE::tr_probe_num(const std::string& x)const
{
  const COMMON_SWITCH* c = prechecked_cast<const COMMON_SWITCH*>(common());
  assert(c);
  if (Umatch(x, "state|s0 ")) {
    return _state[0];
  }else if (Umatch(x, "ic ")) {
    return c->_ic;
  }else if (Umatch(x, "s1 ")) {
    return _state[1];
  }
  return ELEMENT::tr_probe_num(x);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
void DEV_CSWITCH::expand()
{
  SWITCH_BASE::expand();
  _input = dynamic_cast<const ELEMENT*>(find_in_my_scope(_input_label));
  if (!_input) {
    throw Exception(long_label() + ": " + _input_label + " cannot be used as input");
  }else{
  }
}
/*--------------------------------------------------------------------------*/
void DEV_CSWITCH::set_port_by_index(uint_t Num, std::string& Value)
{
  if (Num == 2) {
    _input_label = Value;
  }else{
    SWITCH_BASE::set_port_by_index(Num, Value);
  }
}
/*--------------------------------------------------------------------------*/
bool DEV_CSWITCH::node_is_connected(uint_t i)const
{
  if (i == 2) {
    return _input_label != "";
  }else{
    return SWITCH_BASE::node_is_connected(i);
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
DEV_VSWITCH p2;
DEV_CSWITCH p3;
DISPATCHER<CARD>::INSTALL
  d2(&device_dispatcher, "S|vswitch", &p2),
  d3(&device_dispatcher, "W|cswitch|iswitch", &p3);

MODEL_SWITCH p1(&p2);
MODEL_SWITCH p4(&p3);
DISPATCHER<MODEL_CARD>::INSTALL
  d1(&model_dispatcher, "sw", &p1),
  d4(&model_dispatcher, "csw", &p4);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
