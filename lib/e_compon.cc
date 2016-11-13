/*$Id: e_compon.cc 2016/03/23 al $ -*- C++ -*-
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
 * Base class for elements of a circuit
 */
#include "u_lang.h"
#include "e_model.h"
#include "e_elemnt.h"
#include "io_trace.h"
#include "io_misc.h"

#include <typeinfo>
using namespace std;

/*--------------------------------------------------------------------------*/
using notstd::to_lower;
using std::string;
/*--------------------------------------------------------------------------*/
using notstd::to_lower;
using std::string;
/*--------------------------------------------------------------------------*/
COMMON_COMPONENT::COMMON_COMPONENT(const COMMON_COMPONENT& p)
  :CKT_BASE(p),
   _tnom_c(p._tnom_c),
   _dtemp(p._dtemp),
   _temp_c(p._temp_c),
   _mfactor(p._mfactor),
   _value(p._value),
   _modelname(p._modelname),
   _model(p._model),
   _attach_count(0)
{
  trace1(("COMMON_COMPONENT copy, modelname: "+p._modelname), hp(_model) );
}
/*--------------------------------------------------------------------------*/
COMMON_COMPONENT::COMMON_COMPONENT(int c)
  :CKT_BASE(),
   _tnom_c(NOT_INPUT),
   _dtemp(0),
   _temp_c(NOT_INPUT),
   _mfactor(1),
   _value(0),
   _modelname(""), // !
   _model(0),
   _attach_count(c)
{
  trace2("COMMON_COMPONENT::COMMON_COMPONENT"+ _modelname, c, hp(this));
}
/*--------------------------------------------------------------------------*/
COMMON_COMPONENT::~COMMON_COMPONENT()
{
  assert(_attach_count == 0 || _attach_count == CC_STATIC);
}
/*--------------------------------------------------------------------------*/
  // void	COMPONENT::attach_common(COMMON_COMPONENT*c) {
  //   COMMON_COMPONENT::attach_common(c,&_common);}
/*--------------------------------------------------------------------------*/
void COMMON_COMPONENT::attach_common(COMMON_COMPONENT*c, COMMON_COMPONENT**to)
{
  assert(to);
  if (c == *to) {
    // The new and old are the same object.  Do nothing.
  }else if (!c) { untested();
    // There is no new common.  probably a simple element
    detach_common(to);
  }else if (!*to) {
    // No old one, but have a new one.
    ++(c->_attach_count);
    *to = c;
    trace1("COMMON_COMPONENT::attach_common attached new ", hp(c) );
    // trace2( c->modelname(), c->attach_count(), hp(c));
  }else if (*c != **to) {
    // They are different, usually by edit.
    trace1("different... "+ c->name() + " " + c->modelname() , c->_attach_count );
    detach_common(to);
    ++(c->_attach_count);
    *to = c;
  }else if (c->_attach_count == 0) {
    // The new and old are identical.
    // Use the old one.
    // The new one is not used anywhere, so throw it away.
    trace1("COMMON_COMPONENT::attach_common, equal and c unused, deleting", hp(c));
    delete c;
  }else{ untested();
    trace0("COMMON_COMPONENT::attach_common same twice");
    untested();
    // The new and old are identical.
    // Use the old one.
    // The new one is also used somewhere else, so keep it.
  }
}
/*--------------------------------------------------------------------------*/
void COMMON_COMPONENT::detach_common(COMMON_COMPONENT** from)
{
  assert(from);
  if (*from) {
    assert((**from)._attach_count > 0);
    //assert((**from)._attach_count != CC_STATIC );
    --((**from)._attach_count);
    if ((**from)._attach_count == 0) {
      delete *from;
    }else{
    }
    *from = NULL;
  }else{
  }
}
/*--------------------------------------------------------------------------*/
void COMMON_COMPONENT::attach_model(const COMPONENT* d)const
{
  assert(d);
  trace1("COMMON_COMPONENT::attach_model to", d->short_label());

  _model = d->find_model(modelname());
  assert(_model);
  trace2("COMMON_COMPONENT::attach_model attached ", hp(_model), hp(this) );
}
/*--------------------------------------------------------------------------*/
void COMMON_COMPONENT::parse_modelname(CS& cmd)
{
  set_modelname(cmd.ctos(TOKENTERM));
  trace0(("COMMON_COMPONENT::parse_modelname " + _modelname).c_str() );
}
/*--------------------------------------------------------------------------*/
// called only by COMMON_COMPONENT::parse_obsolete
bool COMMON_COMPONENT::parse_param_list(CS& cmd)
{
  trace0(("COMMON_COMPONENT::parse_param_list " + cmd.tail()).c_str() );
  unsigned start = cmd.cursor();
  unsigned here = cmd.cursor();
  do{
    parse_params_obsolete_callback(cmd); //BUG//callback
  }while (cmd.more() && !cmd.stuck(&here));
  return cmd.gotit(start);
}
/*--------------------------------------------------------------------------*/
void COMMON_COMPONENT::parse_common_obsolete_callback(CS& cmd) //used
{
  trace1("COMMON_COMPONENT::parse_common_obsolete_callback ", cmd.tail() );
  if (cmd.skip1b('(')) {
    // start with a paren
    unsigned start = cmd.cursor();
    parse_param_list(cmd);
    if (cmd.gotit(start)) {		// ( params ( ....
      // named args before num list
      if (cmd.skip1b('(')) {		// ( params ( list ) params )
	parse_numlist(cmd);
	if (!cmd.skip1b(')')) {untested();
	  cmd.warn(bWARNING, "need )");
	}else{
	}
      }else{				// ( params list params )
	parse_numlist(cmd);		//BUG//
      }
      parse_param_list(cmd);
      if (!cmd.skip1b(')')) {untested();
	cmd.warn(bWARNING, "need )");
      }else{
      }
    }else{
      // no named args before num list
      // but there's a paren
      // not sure whether it belongs to all args or to num list
      if (cmd.skip1b('(')) {		// ( ( list ) params )
	// two parens
	parse_numlist(cmd);
	if (!cmd.skip1b(')')) {untested();
	  cmd.warn(bWARNING, "need )");
	}else{
	}
	parse_param_list(cmd);
	if (!cmd.skip1b(')')) {untested();
	  cmd.warn(bWARNING, "need )");
	}else{
	}
      }else{				// ( list ...
	// only one paren
	parse_numlist(cmd);
	if (cmd.skip1b(')')) {		// ( list ) params
	  // assume it belongs to num list
	  // and named params follow
	  parse_param_list(cmd);
	}else{				// ( list params )
	  // assume it belongs to all args
	  parse_param_list(cmd);
	  if (!cmd.skip1b(')')) {
	    cmd.warn(bWARNING, "need )");
	  }else{
	  }
	}
      }
    }
  }else{
    // does not start with a paren
    unsigned start = cmd.cursor();
    parse_param_list(cmd);
    if (cmd.gotit(start)) {
      if (cmd.skip1b('(')) {		// params ( list ) params
	parse_numlist(cmd);
	if (!cmd.skip1b(')')) {untested();
	  cmd.warn(bWARNING, "need )");
	}else{
	}
      }else if (!(cmd.is_alpha())) {	// params list params
	parse_numlist(cmd);
      }else{				// params   (only)
      }
    }else{				// list params
      assert(!(cmd.skip1b('(')));
      parse_numlist(cmd);
    }
    parse_param_list(cmd);
    if (cmd.skip1b(')')) {
      cmd.warn(bWARNING, start, "need (");
    }else{
    }
  }
}
/*--------------------------------------------------------------------------*/
void COMMON_COMPONENT::print_common_obsolete_callback(OMSTREAM& o, LANGUAGE* lang)const
{
  trace0(("COMMON_COMPONENT::print_common_obsolete_callback "+ _modelname).c_str());
  assert(lang);
  print_pair(o, lang, "tnom", _tnom_c,  _tnom_c.has_hard_value());
  print_pair(o, lang, "dtemp",_dtemp,   _dtemp.has_hard_value());
  print_pair(o, lang, "temp", _temp_c,  _temp_c.has_hard_value());
  print_pair(o, lang, "m",    _mfactor, _mfactor.has_hard_value());
}
/*--------------------------------------------------------------------------*/
void COMMON_COMPONENT::set_param_by_index(int i, std::string& value, int Offset)
{
  IString Value(value);
  switch (i) {
  case 0:untested();  _tnom_c = Value; break;
  case 1:untested();  _dtemp = Value; break;
  case 2:untested();  _temp_c = Value; break;
  case 3:  _mfactor = Value; break;
  default:untested();
			 throw Exception_Too_Many(unsigned(i), 3u, Offset); break;
  }
}
/*--------------------------------------------------------------------------*/
bool COMMON_COMPONENT::param_is_printable(int i)const
{
  switch (i) {
  case 0:  return _tnom_c.has_hard_value();
  case 1:  return _dtemp.has_hard_value();
  case 2:  return _temp_c.has_hard_value();
  case 3:  return _mfactor.has_hard_value();
  default:untested(); return false;
  }
}
/*--------------------------------------------------------------------------*/
std::string COMMON_COMPONENT::param_name(int i)const
{
  switch (i) {
  case 0:itested();  return "tnom";
  case 1:itested();  return "dtemp";
  case 2:itested();  return "temp";
  case 3:  return "m";
  default:untested(); return "";
  }
}
/*--------------------------------------------------------------------------*/
std::string COMMON_COMPONENT::param_name(int i, int j)const
{itested();
  return (j==0) ? param_name(i) : "";
}
/*--------------------------------------------------------------------------*/
std::string COMMON_COMPONENT::param_value(int i)const
{
  switch (i) {
  case 0:untested();  return _tnom_c.string();
  case 1:untested();  return _dtemp.string();
  case 2:untested();  return _temp_c.string();
  case 3:  return _mfactor.string();
  default:untested(); return "";
  }
}
/*--------------------------------------------------------------------------*/
void COMMON_COMPONENT::precalc_last(const CARD_LIST* Scope)
{
  assert(Scope);
  _tnom_c.e_val(OPT::tnom_c, Scope);
  _dtemp.e_val(0., Scope);
  _temp_c.e_val(_sim->_temp_c + _dtemp, Scope);
  _mfactor.e_val(1, Scope);
  _value.e_val(0, Scope);
}
/*--------------------------------------------------------------------------*/
void COMMON_COMPONENT::tt_commit(ELEMENT*x)const
{ untested();
  assert(_model);
  _model->do_tt_commit(x);
}
/*--------------------------------------------------------------------------*/
void COMMON_COMPONENT::tr_eval(ELEMENT*x)const
{ untested();
  assert(_model);

  //printf("typeid(_model): %s", typeid(_model).name());

  _model->tr_eval(x);
}
/*--------------------------------------------------------------------------*/
void COMMON_COMPONENT::ac_eval(ELEMENT*x)const
{untested(); // bug?
  assert(_model);
  _model->ac_eval(x);
}
/*--------------------------------------------------------------------------*/
bool COMMON_COMPONENT::operator==(const COMMON_COMPONENT& x)const
{
  return (_modelname == x._modelname
	  && _model == x._model
	  && _tnom_c == x._tnom_c
	  && _dtemp == x._dtemp
	  && _temp_c == x._temp_c
	  && _mfactor == x._mfactor
	  && _value == x._value);
}
/*--------------------------------------------------------------------------*/
map<string, PARA_BASE COMMON_COMPONENT::*> COMMON_COMPONENT::_param_dict
  = boost::assign::map_list_of
("tnom", (PARA_BASE COMMON_COMPONENT::*)  (&COMMON_COMPONENT::_tnom_c))
("dtemp", (PARA_BASE COMMON_COMPONENT::*)  (&COMMON_COMPONENT::_dtemp))
("temp", (PARA_BASE COMMON_COMPONENT::*)  (&COMMON_COMPONENT::_temp_c))
("m", (PARA_BASE COMMON_COMPONENT::*)  (&COMMON_COMPONENT::_mfactor));
/*--------------------------------------------------------------------------*/
void COMMON_COMPONENT::set_param_by_name(std::string Name, std::string Value)
{
  PARA_BASE COMMON_COMPONENT::* x = (OPT::case_insensitive)?
     (_param_dict[to_lower(Name)]) : (_param_dict[Name]);
  if(x) {
    PARA_BASE* p = &(this->*x);
    *p = Value;
    return;
  }

  if (has_parse_params_obsolete_callback()) {
    std::string args(Name + "=" + Value);
    CS cmd(CS::_STRING, args); //obsolete_callback
    trace3("COMMON_COMPONENT::set_param_by_name", Name, Value, cmd.fullstring());
    bool ok = parse_params_obsolete_callback(cmd); //BUG//callback
    if (!ok) {untested();
      throw Exception_No_Match(Name);
    }else{
    }
  }else if (Umatch(Name, "tnom")) { untested();
    _tnom_c = Value;
  }else if (Umatch(Name, "dtemp")) { untested();
    _dtemp = Value;
  }else if (Umatch(Name, "temp")) { untested();
    _temp_c = Value;
  }else if (Umatch(Name, "m")) {
    _mfactor = Value;
  }else{
    //BUG// ugly linear search
    for (int i = param_count() - 1;  i >= 0;  --i) {
      for (int j = 0;  param_name(i,j) != "";  ++j) {
	if (Umatch(Name, param_name(i,j) + ' ')) {
          cerr << typeid(this).name() << " linear search for " << Name << ": ";
          incomplete();
	  set_param_by_index(i, Value, 0/*offset*/);
	  return; //success
	}else{
	  //keep looking
	}
      }
    }
    throw Exception_No_Match(Name);
  }
}
/*--------------------------------------------------------------------------*/
//BUG// This is a kluge for the spice_wrapper, to disable virtual functions.
// It is called during expansion only.

void COMMON_COMPONENT::Set_param_by_name(std::string Name, std::string Value)
{untested();
  assert(!has_parse_params_obsolete_callback());

  //BUG// ugly linear search
  for (int i = COMMON_COMPONENT::param_count() - 1;  i >= 0;  --i) {untested();
    for (int j = 0;  COMMON_COMPONENT::param_name(i,j) != "";  ++j) {untested();
      if (Umatch(Name, COMMON_COMPONENT::param_name(i,j) + ' ')) {untested();
	COMMON_COMPONENT::set_param_by_index(i, Value, 0/*offset*/);
	return; //success
      }else{untested();
	//keep looking
      }
    }
  }
  throw Exception_No_Match(Name);
}
/*--------------------------------------------------------------------------*/
bool COMMON_COMPONENT::parse_numlist(CS&)
{
  return false;
}
/*--------------------------------------------------------------------------*/
bool COMMON_COMPONENT::parse_params_obsolete_callback(CS& cmd)
{
  return ONE_OF
    || Get(cmd, "tnom",   &_tnom_c)
    || Get(cmd, "dtemp",  &_dtemp)
    || Get(cmd, "temp",   &_temp_c)
    || Get(cmd, "m",	  &_mfactor)
    || Get(cmd, "mfactor",&_mfactor)
    ;
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
COMPONENT::COMPONENT()
  :CARD(),
   _common(0),
   _value(0),
   _mfactor(1),
   _mfactor_fixed(NOT_VALID),
   _converged(false),
   _q_for_eval(-1),
   _time_by(),
   _amps(0),
   _adp(0)
{
  trace1("COMPONENT::COMPONENT", long_label());
  if (_sim) {
    _sim->uninit();
  } else {
  }
}
/*--------------------------------------------------------------------------*/
COMPONENT::COMPONENT(const COMPONENT& p)
  :CARD(p),
   _common(0),
   _value(p._value),
   _mfactor(p._mfactor),
   _mfactor_fixed(p._mfactor_fixed),
   _converged(p._converged),
   _q_for_eval(-1),
   _time_by(p._time_by),
   _amps(0),
   _adp(0)
{
  if (_sim) {
    _sim->uninit();
  } else { untested();
  }
  trace0("attaching common");
  if(p._adp){
    attach_adp(p._adp->clone());
  }
  attach_common(p._common);
  assert(_common == p._common);
}
/*--------------------------------------------------------------------------*/
COMPONENT::~COMPONENT()
{
  if (_amps)     free (_amps);
  detach_common();
  if (_sim) {
    _sim->uninit();
  } else { itested();
  }
}
/*--------------------------------------------------------------------------*/
bool COMPONENT::node_is_grounded(uint_t i)const 
{
  assert(_n);
  assert(i != INVALID_NODE);
  assert(i < net_nodes());
  return _n[i].is_grounded();
}
/*--------------------------------------------------------------------------*/
bool COMPONENT::node_is_connected(uint_t i)const 
{
  assert(_n);
  assert(i != INVALID_NODE);
  if (i >= net_nodes()) {
    trace3( "COMPONENT::node_is_connected", i, net_nodes(), _net_nodes );
    assert(false);
  }
  return _n[i].is_connected();
}
/*--------------------------------------------------------------------------*/
void COMPONENT::set_port_by_name(std::string& int_name, std::string& ext_name)
{
  for (uint_t i=0; i<max_nodes(); ++i) {
    if (int_name == port_name(i)) {
      set_port_by_index(i, ext_name);
      return;
    }else{
    }
  }
  throw Exception_No_Match(int_name);
}
/*--------------------------------------------------------------------------*/
#ifdef DDO_TRACE
#include "e_cardlist.h"
#include "e_node.h"
#include "u_nodemap.h"
class NODE;
class NODE_MAP;
inline void trace_nodenames(const CARD_LIST* scope){
  trace0("CARD_LIST tracing nodenames");
  NODE_MAP* nm = scope->nodes();
  for (NODE_MAP::const_iterator ni = nm->begin(); ni != nm->end(); ++ni) {
    //NODE_BASE* n = (*ni).second;
    string label = (*ni).first;
    //trace2("CARD_LIST:... nodename ", label, n->user_number() );
  }
}
#else
inline void trace_nodenames(const CARD_LIST*){}
#endif
/*--------------------------------------------------------------------------*/
void COMPONENT::set_port_by_index(uint_t num, std::string& ext_name)
{
  if (num < max_nodes()) {
    _n[num].new_node(ext_name, this);
    trace1("COMPONENT::set_port_by_index", _n[num].t_());
    if (num+1 > _net_nodes) {
      // make the list bigger
      _net_nodes = num+1;
    }else{
      // it's already big enough, probably assigning out of order
    }
  }else{
    throw Exception_Too_Many(num+1, max_nodes(), 0/*offset*/);
  }
  trace_nodenames(scope());
}
/*--------------------------------------------------------------------------*/
void COMPONENT::set_port_to_ground(uint_t num)
{
  if (num < max_nodes()) {
    _n[num].set_to_ground(this);
    if (num+1 > _net_nodes) {
      _net_nodes = num+1;
    }else{untested();
    }
  }else{untested();
    throw Exception_Too_Many(num+1, max_nodes());
  }
}
/*--------------------------------------------------------------------------*/
void COMPONENT::set_dev_type(const std::string& new_type)
{
  trace1("COMPONENT::set_dev_type", new_type);
  if (common()) {
    if (new_type != dev_type()) {
      COMMON_COMPONENT* c = common()->clone();
      assert(c);
      c->set_modelname(new_type);
      attach_common(c);
    }else{
    }
  }else{
    CARD::set_dev_type(new_type);
  }
}
/*--------------------------------------------------------------------------*/
void COMPONENT::print_args_obsolete_callback(OMSTREAM& o, LANGUAGE* lang)const
{
  assert(lang);
  assert(has_common());
  trace0(("COMPONENT::print_args_obsolete_callback "+ short_label()).c_str());
  common()->print_common_obsolete_callback(o, lang);
  if(comment()!=""){ untested();
    o << " ; " << comment();
  }else{
  }
}
/*--------------------------------------------------------------------------*/
void COMPONENT::deflate_common()
{untested();
  unreachable();
  trace0("COMPONENT::deflate_common");
  if (has_common()) {untested();
    COMMON_COMPONENT* deflated_common = mutable_common()->deflate();
    if (deflated_common != common()) {untested();
      attach_common(deflated_common);
    }else{untested();
    }
  }else{untested();
    unreachable();
  }
}
/*--------------------------------------------------------------------------*/
void COMPONENT::expand()
{
  CARD::expand();
  if (has_common()) {
    COMMON_COMPONENT* new_common = common()->clone();
    new_common->expand(this);
    COMMON_COMPONENT* deflated_common = new_common->deflate();

    if ( deflated_common != common()) {
      attach_common(deflated_common);
    }else if ( deflated_common != new_common )  {
      delete new_common;
    }else{
      untested();
    }
  }else{
  }
  trace1("COMPONENT::expand done", long_label());
}
/*--------------------------------------------------------------------------*/
void COMPONENT::precalc_first()
{
  trace4("COMPONENT::precalc_first", long_label(), _value, *(scope()->params()), hp(mutable_common()));
  CARD::precalc_first();
  if (has_common()) {
    try {
      mutable_common()->precalc_first(scope());
    }catch (Exception_Precalc& e) {untested();
      trace0("COMPONENT::precalc_first exception...");
      error(bWARNING, long_label() + ": " + e.message());
    }
    _mfactor = common()->mfactor();
  }else{
  }

  //BUG//  _mfactor must be in precalc_first

  _mfactor.e_val(1, scope());
  trace2("COMPONENT::precalc_first", long_label(), double(_mfactor));
  if (const COMPONENT* o = prechecked_cast<const COMPONENT*>(owner())) {
    _mfactor_fixed = o->mfactor() * _mfactor;
  }else{
    _mfactor_fixed =  _mfactor;
  } 
  trace4("COMPONENT::precalc_first done", long_label(), _mfactor_fixed, _value, common()?common()->value():-1.);
}
/*--------------------------------------------------------------------------*/
void COMPONENT::precalc_last()
{
  CARD::precalc_last();
  if (has_common()) {
    trace1("COMPONENT::precalc_last", long_label());
    try {
      mutable_common()->precalc_last(scope());
    }catch (Exception_Precalc& e) {
      error(bWARNING, long_label() + ": " + e.message());
    }
  }else{
  }

  _value.e_val(0.,scope());
  trace4("COMPONENT::precalc_last done", long_label(), _mfactor_fixed, _value, common()?common()->value():-1.);
}
/*--------------------------------------------------------------------------*/
void COMPONENT::map_nodes()
{
  trace5("COMPONENT::map_nodes", long_label(), ext_nodes(), int_nodes(),
      max_nodes(), net_nodes());
  if(!is_device()){ unreachable();
    return;
  }
  //assert(min_nodes() <= net_nodes());
  assert(net_nodes() <= max_nodes());
  //assert(ext_nodes() + int_nodes() == matrix_nodes());

  for (uint_t ii = 0; ii < ext_nodes()+int_nodes(); ++ii) {
    unsigned oldm = _n[ii].m_(); USE(oldm);
    assert(_n[ii].t_() <= NODE::_sim->_total_nodes || _n[ii].t_()==INVALID_NODE || _n[ii].is_adp() );
    _n[ii].map();
    assert(_n[ii].m_() <= _sim->_total_nodes || _n[ii].m_() == INVALID_NODE || _n[ii].is_adp() );
    unsigned user_number = (_n[ii].n_())? _n[ii].n_()->user_number(): 0; USE(user_number);
    unsigned matrix_number = (_n[ii].n_())? _n[ii].n_()->matrix_number(): 0;
    IString node_label = (_n[ii].n_())? _n[ii].n_()->long_label(): "";
    // matrix_number not initialized yet.
    USE(user_number); USE(matrix_number);
    trace5("COMPONENT::map_nodes done", long_label(), ii, _n[ii].t_(), _n[ii].m_(), _n[ii].is_adp() );
    user_number = matrix_number = 0;
    node_label="";
  }

  if (subckt()) {
    subckt()->map_nodes();
  }else{
  }
}
/*--------------------------------------------------------------------------*/
void COMPONENT::tr_iwant_matrix()
{
  if (is_device()) {
    assert(matrix_nodes() == 0);
    if (subckt()) {
      subckt()->tr_iwant_matrix();
    }else{ // untested();
    }
  }else{
  }
}
/*--------------------------------------------------------------------------*/
void COMPONENT::ac_iwant_matrix()
{
  if (is_device()) {
    assert(matrix_nodes() == 0);
    if (subckt()) {
      subckt()->ac_iwant_matrix();
    }else{ // untested();
    }
  }else{
  }
}
/*--------------------------------------------------------------------------*/
/* set: set parameters, used in model building
 */
void COMPONENT::set_parameters(const std::string& Label, CARD *Owner,
			       COMMON_COMPONENT *Common, double Value,
			       uint_t , hp_float_t [],
			       uint_t node_count, const node_t Nodes[])
{
  trace6("COMPONENT::set_parameters", long_label(), hp(Common), hp(common()), max_nodes(), node_count, net_nodes());
  set_label(Label);
  set_owner(Owner);
  set_value(Value);
  attach_common(Common);
  _net_nodes = node_count;
  assert(node_count <= max_nodes());

  if (node_count > net_nodes()) {
    error(bDANGER, "net nodes problem in %s: only have %i, passed %i\n", long_label().c_str(), net_nodes(), node_count);
  }
  assert(node_count <= net_nodes());

  trace0("COMPONENT::set_parameters copy nodes...");
  notstd::copy_n(Nodes, node_count, _n);

  for(uint_t i=0; i<node_count;i++){
    trace6("COMPONENT::set_parameters", long_label(), i, Nodes[i].is_adp(), _n[i].is_adp(), _n[i].m_(), _n[i].t_());
  }
}
/*--------------------------------------------------------------------------*/
/* set_slave: force evaluation whenever the owner is evaluated.
 * and: deactivate q_eval...
 */
void COMPONENT::set_slave()
{
  mark_always_q_for_eval();
  if (subckt()) {
    subckt()->set_slave();
  }else{
  }
}
/*--------------------------------------------------------------------------*/
void COMPONENT::set_value(double v, COMMON_COMPONENT* c)
{
  if (c != _common) {
    detach_common();
    attach_common(c);
  }else{
  }
  set_value(v);
}
/*--------------------------------------------------------------------------*/
void COMPONENT::set_param_by_name(std::string Name, std::string Value)
{
  trace3("COMPONENT::set_param_by_name", Name, has_common(), long_label());
  if (has_common()) {
    COMMON_COMPONENT* c = common()->clone();
    assert(c);
    c->set_param_by_name(Name, Value);
    attach_common(c);
  }else{
    CARD::set_param_by_name(Name, Value);
  }
}
/*--------------------------------------------------------------------------*/
void COMPONENT::set_param_by_index(int i, std::string& value, int offset)
{
  IString Value(value);
  if (has_common()) {untested();
    COMMON_COMPONENT* c = common()->clone();
    assert(c);
    c->set_param_by_index(i, value, offset);
    attach_common(c);
  }else{
    switch (COMPONENT::param_count() - 1 - i) {
    case 0: _value = Value; break;
    case 1:untested(); _mfactor = Value; break;
    default:untested(); CARD::set_param_by_index(i, value, offset);
    }
  }
}
/*--------------------------------------------------------------------------*/
bool COMPONENT::param_is_printable(int i)const
{
  if (has_common()) {
    return common()->param_is_printable(i);
  }else{
    switch (COMPONENT::param_count() - 1 - i) {
    case 0:  return value().has_hard_value();
    case 1:  return _mfactor.has_hard_value();
    default:untested(); return CARD::param_is_printable(i);
    }
  }
}
/*--------------------------------------------------------------------------*/
std::string COMPONENT::param_name(int i)const
{
  if (has_common()) {
    return common()->param_name(i);
  }else{
    switch (COMPONENT::param_count() - 1 - i) {
    case 0:  return value_name();
    case 1:  return "m";
    default:untested(); return CARD::param_name(i);
    }
  }
}
/*--------------------------------------------------------------------------*/
std::string COMPONENT::param_name(int i, int j)const
{
  if (has_common()) {untested();
    return common()->param_name(i,j);
  }else{
    if (j == 0) {
      return param_name(i);
    }else if (i >= CARD::param_count()) {
      return "";
    }else{untested();
      return CARD::param_name(i,j);
    }
  }
}
/*--------------------------------------------------------------------------*/
std::string COMPONENT::param_value(int i)const
{
  if (has_common()) {
    return common()->param_value(i);
  }else{
    switch (COMPONENT::param_count() - 1 - i) {
    case 0:  return value().string();
    case 1:  return _mfactor.string();
    default:untested(); return CARD::param_value(i);
    }
  }
}
/*--------------------------------------------------------------------------*/
const std::string COMPONENT::port_value(uint_t i)const 
{
  assert(_n);
  assert(i !=INVALID_NODE);
  assert(i < net_nodes());
  return _n[i].short_label().to_string();
}
/*--------------------------------------------------------------------------*/
const std::string COMPONENT::current_port_value(uint_t)const 
{untested();
  unreachable();
  static std::string s;
  return s;
}
/*--------------------------------------------------------------------------*/
const MODEL_CARD* COMPONENT::find_model(const std::string& modelname)const
{
  trace2("COMPONENT::find_model", short_label(), modelname);
  if (modelname == "") {
    throw Exception(long_label() + ": missing args -- need model name");
    unreachable();
    return NULL;
  }else{
    const CARD* c = NULL;
    {
      int bin_count = 0;
      for (const CARD* Scope = this; Scope && !c; Scope = Scope->owner()) {
	// start here, looking out
	try {
	  c = Scope->find_in_my_scope(modelname);
          trace1("COMPONENT::find_model found model in scope", hp(c));
	}catch (Exception_Cant_Find& e1) {
          trace3("COMPONENT::find_model found no model in scope",
              modelname, Scope->long_label(), hp(Scope->scope()));
	  // didn't find plain model.  try binned models
	  bin_count = 0;
	  for (;;) {
	    // loop over binned models
	    std::string extended_name = modelname + '.' + ::to_string(++bin_count);
	    try {
	      c = Scope->find_in_my_scope(extended_name);
	    }catch (Exception_Cant_Find& e2) {
	      // that's all .. looked at all of them
	      c = NULL;
	      break;
	    }
	    const MODEL_CARD* m = dynamic_cast<const MODEL_CARD*>(c);
	    if (m && m->is_valid(this)) {
	      //matching name and correct bin
	      break;
	    }else{
              trace0("COMPONENT::find_model invalid");
	      // keep looking
	    }
	  }
	}
      }
      if (!c) {
        // this does not work!
        // const MODEL_CARD* m = dynamic_cast<const MODEL_CARD*>(c);
        // if (m && m->is_valid(this)) { untested();
        // }else
        if (bin_count <= 1) {
          if (model_dispatcher[modelname]) {
            trace1("COMPONENT::find_model there's a model... ", modelname);
          }
          trace1("COMPONENT::find_model Exception", modelname);
	  throw Exception_Cant_Find(long_label(), modelname);
	}else{
	  throw Exception(long_label() + ": no bins match: " + modelname);
	}
	unreachable();
      }else{
      }
    }
    // found something, what is it?
    assert(c);
    const MODEL_CARD* model = dynamic_cast<const MODEL_CARD*>(c);
    if (!model) {
      trace0("COMPONENT::find_model no model here");
      throw Exception_Type_Mismatch(long_label(), modelname, ".model");
    }else if (!model->is_valid(this)) {untested();
      error(bWARNING, long_label() + ", " + modelname
	   + "\nmodel and device parameters are incompatible, using anyway\n");
    }
    assert(model);
    return model;
  }
}
/*--------------------------------------------------------------------------*/
double COMPONENT::tr_probe_num(const std::string& x)const
{
  CS cmd(CS::_STRING, x);
  ADP_CARD* a=adp();
  if (cmd.umatch("v")) {
    int nn = cmd.ctoi();
    return (nn > 0 && nn <= int(net_nodes())) ? _n[nn-1].v0() : NOT_VALID;
  }else if (Umatch(x, "amps |amps0 ")) {
    return ( 17  );
  }else if (Umatch(x, "stress ") || Umatch(x, "hci ") ) {
    return(19);
    a->tr_probe_num(x);
  }else if (Umatch(x, "brel ")) {
    ///    std::cout << "rel\n";
    return ( tr_behaviour_rel  );
  }else if (Umatch(x, "bdel ")) {
    return ( tr_behaviour_del  );
  }else if (Umatch(x, "sv ")) { // some value (debugging)
         return (  tr_amps_diff_cur()  );
  }else if (Umatch(x, "next{time} ")) {
    return (_time_by._error_estimate < BIGBIG) ? _time_by._error_estimate : 0;
  }else if (Umatch(x, "error{time} ")) { itested();
    return (_time_by._error_estimate < BIGBIG) ? _time_by._error_estimate-_sim->_time0 : 0;
  }else if (Umatch(x, "timef{uture} ")) {
    return (_time_by._error_estimate < _time_by._event) 
      ? _time_by._error_estimate
      : _time_by._event;
  }else if (Umatch(x, "te{mperature} ")) {
    if(has_common())
      return common()->temp();
    return _sim->_temp_c;
  }else if (Umatch(x, "event{time} ")) {
    return (_time_by._event < BIGBIG) ? _time_by._event : 0;
  }else if (Umatch(x, "const " )) {
    return is_constant();
  }else if (Umatch(x, "conv{erged} ")) {
    return double( _converged );
#ifndef NDEBUG
  }else if (Umatch(x, "_m ")) {
    return double(hp( _common->model() ) );
  }else if (Umatch(x, "_c ")) {
    return double(  hp(_common) );
  }else if (Umatch(x, "common " )) {
    return hp(common());
#endif
  }

  return CARD::tr_probe_num(x);
}
/*--------------------------------------------------------------------------*/
void COMPONENT::tr_save_amps(int)
{
  return;
#if 0
  std::cerr << "COMPONENT::tr_save_amps " << short_label() << "\n";
  int j = net_nodes() - 1;
  assert(j==1);
  j=1;
  int k = j;
  hp_float_t tramps = 0;//tr_amps();
  hp_float_t* trampsp=&tramps;

  std::cerr << "saving _amps[ " << n << " ]. have " << net_nodes() << " nodes\n";


  if (_amps==0) _tr_amps_diff_cur = 0;
  if (_amps==0) _tr_amps_diff_max = 0; // maximum der delta i in diesem zeitschritt.

  while (j--> 0){
    if(_amps!=0) {
      _tr_amps_diff_cur = _amps[n*k + j] - trampsp[j];
      _tr_amps_diff_max = fmax( _tr_amps_diff_max, _tr_amps_diff_cur );
    }

    std::cerr << short_label() << ": saving _amps[ " << n*k+j << " ]" << _amps << " \n";

    hp_float_t tmp=trampsp[j];

    std::cerr << " have " << tmp << "\n";
    _amps_new[ n*k + j ]= tmp;
  }

  trace1("COMPONENT::tr_save_amps" ,  _tr_amps_diff_max);
  tr_behaviour_del = _tr_amps_diff_max;
#endif
}
/*--------------------------------------------------------------------------*/
void COMPONENT::tt_behaviour_update()
{
  trace1("COMPONENT::tt_behaviour_update" ,  tr_behaviour_del);

  tt_behaviour_del += tr_behaviour_del;
  tt_behaviour_rel += tr_behaviour_rel;

  //global bahaviour = maximum device bahaviour.
  if(tr_behaviour_del > CKT_BASE::tr_behaviour_del) 
  {
//    std::cerr << "COMPONENT::tr_behaviour " << _sim->_Time0  
//      << "dev:"  << short_label() << " " << tr_behaviour_del << "CKT_BASE: " << CKT_BASE::tr_behaviour_del << "\n";
  }

  CKT_BASE::tr_behaviour_del = fmax( CKT_BASE::tr_behaviour_del, tr_behaviour_del);
  CKT_BASE::tr_behaviour_rel = fmax( CKT_BASE::tr_behaviour_rel, tr_behaviour_rel);
  //std::cerr << "ELEMENT::tr_save_amps: " << short_label() << " " << "del: " << tr_behaviour_del << "rel: " <<  tr_behaviour_rel << "\n";

  assert (tr_behaviour_del >=0 );
}

/*--------------------------------------------------------------------------*/
/* q_eval: queue this device for evaluation on the next pass,
 * with a check against doing it twice.
 */
void COMPONENT::q_eval()
{
  trace1("COMPONENT::q_eval", long_label());
  if(!is_q_for_eval()) {
    mark_q_for_eval();
    _sim->_evalq_uc->push_back(this);
  }else{
  }
}
/*--------------------------------------------------------------------------*/
void COMPONENT::tr_queue_eval()
{
  if(tr_needs_eval()) {
    q_eval();
  }else{
  }
}
/*--------------------------------------------------------------------------*/
TIME_PAIR COMPONENT::tr_review()
{
  _time_by.reset();
  if(has_common()) {
    return _common->tr_review(this);
  }else{
    return _time_by;
  }
}
/*--------------------------------------------------------------------------*/
void COMPONENT::tr_accept()
{
  if(has_common()) {
    _common->tr_accept(this);
  }else{
  }
}
/*--------------------------------------------------------------------------*/
bool COMPONENT::use_obsolete_callback_parse()const
{
  if (has_common()) {
    return common()->use_obsolete_callback_parse();
  }else{untested();
    return false;
  }
}
/*--------------------------------------------------------------------------*/
bool COMPONENT::use_obsolete_callback_print()const
{
  if (has_common()) {
    return common()->use_obsolete_callback_print();
  }else{
    return false;
  }
}
/*--------------------------------------------------------------------------*/
void COMPONENT::obsolete_move_parameters_from_common(const COMMON_COMPONENT* dc)
{
  assert(dc);
  _value   = dc->value();
  _mfactor = dc->mfactor();
}
/*--------------------------------------------------------------------------*/
/* volts_limited: transient voltage, best approximation, with limiting
 */
double COMPONENT::volts_limited(const node_t & n1, const node_t & n2)
{
  bool limiting = false;

  double v1 = n1.v0();
  double v2 = n2.v0();

  trace3("COMPONENT::volts_limited", v1 ,v2, _sim->_time0);
  assert(v1 == v1);
  if (v1 < _sim->_vmin) {
    limiting = true;
    v1 = _sim->_vmin;
  }else if (v1 > _sim->_vmax) {
    limiting = true;
    v1 = _sim->_vmax;
  }

  assert(v2 == v2);
  if (v2 < _sim->_vmin) {
    limiting = true;
    v2 = _sim->_vmin;
  }else if (v2 > _sim->_vmax) {
    limiting = true;
    v2 = _sim->_vmax;
  }

  if (limiting) {
    _sim->_limiting = true;
    if (OPT::dampstrategy & dsRANGE) {
      _sim->_fulldamp = true;
      error(bTRACE, "range limit damp\n");
    }
    if (OPT::picky <= bTRACE) { itested();
      error(bNOERROR,"node limiting (n1,n2,dif) "
	    "was (%g %g %g) now (%g %g %g)\n",
	    n1.v0(), n2.v0(), n1.v0() - n2.v0(), v1, v2, v1-v2);
    }
  }

  return dn_diff(v1,v2);
}
/*--------------------------------------------------------------------------*/
double COMPONENT::tt_probe_num(const std::string& x)const
{
  if (Umatch(x, "bdel ")) {
	 return ( tt_behaviour_del  );
  }else if (Umatch(x, "brel ")) {
	 return ( tt_behaviour_rel  );
  }

  return  tr_probe_num(x);

}
/*--------------------------------------------------------------------------*/
void COMPONENT::tr_do_behaviour(){
	std::cerr << "COMPONENT::tr_do_behaviour() FIXME\n";
}


/*--------------------------------------------------------------------------*/
void COMPONENT::tt_accept()
{

  // not here...
  // double* tmp = _amps;
  //_amps = _amps_new;
  //_amps_new = tmp;

  tt_behaviour_del /= (_sim->_dT0);
  tt_behaviour_rel /= (_sim->_dT0);
}
/*--------------------------------------------------------------------------*/
void COMPONENT::attach_adp(ADP_CARD* a)
{
  if (!a){
    return;
  } else if(_adp){ untested();
    return;
  }
  _adp = a;
  ADP_LIST::adp_list.push_back( a );

}
/*--------------------------------------------------------------------------
 *
 * need to be careful!
 * attach_common will delete a new COMMON if its already there, so better not
 * enqueue it here (somethng like that)
COMMON_COMPONENT* COMMON_COMPONENT::deflate()
{
  for( list<const COMMON_COMPONENT*>::iterator i = _commons.begin();
      i != _commons.end(); ++i ){

    assert(*i);
    if (*this == **i){
      return const_cast<COMMON_COMPONENT*>( *i );
    }
  }
  _commons.push_back(this);
  return this;
}
------------------------*/
// ELEMENT? not yet.
double COMPONENT::tt_review_check_and_convert(double timestep)
{
  double dT = _sim->_dT0;
  double dTmin = _sim->last_time();
  double time_future;
  if (timestep == NEVER) {
    time_future = NEVER;
  }else{
    double atimestep = timestep;
    if (timestep >= dTmin) {
    }else{
      timestep = dTmin;
    }

    if (atimestep < (dT) * OPT::ttreject ) {
      error(bTRACE, "tt step rejected: %s\n", long_label().c_str());
      error(bTRACE, "new\n");
      error(bTRACE, "new=%f  required=%f\n", timestep, dT);
      time_future = _sim->_Time0 - dT + timestep;
      trace3("reject", timestep, dT, time_future);
    }else if (timestep - dTmin < dT - dTmin * .5) {
      time_future = _sim->_Time0 + timestep;
    }else{
      time_future = _sim->_Time0 + timestep;
      trace3("accept", timestep, dT, time_future);
    }
  }
  assert(time_future > 0.);
//  assert(time_future > _Time[1]);
  return time_future;
}
/*-------------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
