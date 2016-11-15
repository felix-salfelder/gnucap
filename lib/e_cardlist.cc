/*                                 -*- C++ -*-
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
 * Base class for "cards" in the circuit description file
 * This file contains functions that process a list of cards
 */
#include "u_time_pair.h"
#include "e_node.h"
#include "u_nodemap.h"
#include "e_model.h"
#include "io_misc.h"
/*--------------------------------------------------------------------------*/
#define trace_func_comp() trace2(__func__, (*ci)->long_label(), CKT_BASE::_sim->iteration_tag())
/*--------------------------------------------------------------------------*/
CARD_LIST::CARD_LIST(const CARD* owner, PARAM_LIST_MAP* p)
  :_parent(NULL),
   _nm(new NODE_MAP),
   _params(p),
   _language(NULL),
   _owner(owner),
   _origin(NULL)
{
  if (p) { // register owner, to disable param delete in ~
    assert(owner);
    _parent = owner->scope();
  }
}
/*--------------------------------------------------------------------------*/
CARD_LIST::CARD_LIST(const CARD* model, CARD* owner,
		     const CARD_LIST* scope, PARAM_LIST_BASE* p)
  :_parent(NULL),
   _nm(new NODE_MAP),
   _params(NULL),
   _language(NULL),
   _owner(owner),
   _origin(model->subckt())
{
  assert(model);
  assert(model->subckt());
  assert(owner);
  assert(!p || scope);

  attach_params(p, scope);
  shallow_copy(model->subckt());
  set_owner(owner);
  map_subckt_nodes(model, owner);
  // rewire_nodenames(model->subckt());
  // _origin = scope; // hack?

  trace2("CARD_LIST::CARD_LIST owner", owner->short_label(), _nm->how_many() );
  trace2("CARD_LIST::CARD_LIST owner", owner->short_label(), model->subckt()->nodes()->how_many() );
}
/*--------------------------------------------------------------------------*/
CARD_LIST::~CARD_LIST()
{
  erase_all();
  delete _nm;
  if (!_parent) {
    delete _params;
  }else{
  }
}
/*--------------------------------------------------------------------------*/
PARAM_LIST* CARD_LIST::params()
{
  if (_parent) {
    trace1("hmm", *(_parent->params()));
  }
  if (!_params) {
    assert(!_parent);
    _params = new PARAM_LIST;
  }else{
  }
  return _params;
}
/*--------------------------------------------------------------------------*/
PARAM_LIST* CARD_LIST::params()const
{
  if (_params) {
    return _params;
  }else if(_parent) { untested();
    return _parent->_params;
  }else{ //BUG//const
    static PARAM_LIST empty_params;
    return &empty_params;
  }
}
/*--------------------------------------------------------------------------*/
CARD_LIST::iterator CARD_LIST::find_again(const IString& short_name,
					  CARD_LIST::iterator Begin)
{
  trace0(("CARD_LIST::find_ name=" + short_name).c_str());
  return notstd::find_ptr(Begin, end(), short_name);
}
/*--------------------------------------------------------------------------*/
CARD_LIST::const_iterator CARD_LIST::find_again(const IString& short_name,
						CARD_LIST::const_iterator Begin)const
{
  return notstd::find_ptr(Begin, end(), short_name);
}
/*--------------------------------------------------------------------------*/
CARD_LIST& CARD_LIST::erase(iterator ci)
{
  assert(ci != end());
  delete *ci;
  _cl.erase(ci);
  return *this;
}
/*--------------------------------------------------------------------------*/
CARD_LIST& CARD_LIST::remove(CARD* c)
{ untested();
  _cl.remove(c);
  return *this;
}
/*--------------------------------------------------------------------------*/
CARD_LIST& CARD_LIST::erase(CARD* c)
{
  delete c;
  _cl.remove(c);
  return *this;
}
/*--------------------------------------------------------------------------*/
/* erase_all: empty the list, destroy contents
 * Beware: something else may be pointing to them, leaving dangling ptr.
 */
CARD_LIST& CARD_LIST::erase_all()
{
  while (!_cl.empty()) {
    delete _cl.back();
    _cl.pop_back();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
CARD_LIST& CARD_LIST::set_owner(CARD* owner)
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).set_owner(owner);
  }
  _owner = owner;
  return *this;
}
/*--------------------------------------------------------------------------*/
/* set_slave: set a whole circuit to "slave" mode.
 * Only useful for subckts.
 */
CARD_LIST& CARD_LIST::set_slave()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).set_slave();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/* expand: expand (flatten) a list of components (subckts)
 * Scan component list.  Expand each subckt: create actual elements
 * for flat representation to use for simulation.
 */
CARD_LIST& CARD_LIST::expand()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).precalc_first();
  }
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).expand_first();
  }
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).expand();
  }
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).expand_last();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
CARD_LIST& CARD_LIST::precalc_first()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).precalc_first();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
//CARD_LIST& CARD_LIST::tt_next()
//{
//  for (iterator ci=begin(); ci!=end(); ++ci) {
//    trace_func_comp();
//    (**ci).tt_next();
//  }
//  return *this;
//}
/*--------------------------------------------------------------------------*/
CARD_LIST& CARD_LIST::precalc_last()
{
  _eq.clear();
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).precalc_last();
    if(!OPT::prequeue){ itested();
    }else if((*ci)->is_constant()){ itested();
    }else{ itested();
      _eq.push_back(*ci);
    }
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/* map_nodes: create mapping between user node names and internal numbers
 */
CARD_LIST& CARD_LIST::map_nodes()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).map_nodes();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
unsigned CARD_LIST::adp_nodes()const{
  unsigned ret = nodes()->how_many_adp();
  for (const_iterator ci=begin(); ci!=end(); ++ci) {
    if ((*ci)->is_device()) {
      if ((*ci)->subckt()){
        ret+=(*ci)->subckt()->adp_nodes();
      }else{
        // untested();
      }
    }
  }
  return ret;
}
/*--------------------------------------------------------------------------*/
// FIXME: store nodecound in static CARDLIST::nodecount
// don't recalculate
unsigned CARD_LIST::total_nodes()const{
  unsigned ret = nodes()->how_many_ckt();
  for (const_iterator ci=begin(); ci!=end(); ++ci) {
    if ((*ci)->is_device()) {
      if ((*ci)->subckt()){
        ret+=(*ci)->subckt()->total_nodes();
      }else{
        // untested();
      }
    }
  }
  return ret;
}
/*--------------------------------------------------------------------------*/
//void CARD_LIST::init_node_count( unsigned* user_nodes, unsigned* subckt_nodes,
//    unsigned* _model_nodes, unsigned* _adp_nodes) const{
//
//  unsigned ret = nodes()->how_many_ckt();
//  for (const_iterator ci=begin(); ci!=end(); ++ci) {
//    if ((*ci)->is_device()) {
//      if ((*ci)->subckt()){
//        ret+=(*ci)->subckt()->total_nodes();
//      }else{
//        untested();
//      }
//    }
//  }
//
//}
/*--------------------------------------------------------------------------*/
NODE_BASE* CARD_LIST::node(IString s) const{
  trace1("CARD_LIST::node", s);

  NODE_MAP* NM = nodes();
  NODE_BASE* ret = (*NM)[s];
  trace1("CARD_LIST::node", hp(ret));
  return(ret);
}
/*--------------------------------------------------------------------------*/
/* tr_iwant_matrix: allocate solution matrix
 * also sets some flags for mixed-mode
 */
CARD_LIST& CARD_LIST::tr_iwant_matrix()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).tr_iwant_matrix();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/* tr_begin: first pass on a new transient simulation (initial DC)
 */
CARD_LIST& CARD_LIST::tr_begin()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).tr_begin();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/* tr_restore: first pass on restarting a transient simulation
 */
CARD_LIST& CARD_LIST::tr_restore()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).tr_restore();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/* keep_ic: latch node voltages into device state
 */
CARD_LIST& CARD_LIST::keep_ic()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).keep_ic();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/*
 */
CARD_LIST& CARD_LIST::tt_behaviour_commit()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
//    trace_func_comp();
    (**ci).tt_behaviour_commit();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/* dc_advance: first pass on a new step in a dc sweep
 */
CARD_LIST& CARD_LIST::dc_advance()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).dc_advance();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
CARD_LIST& CARD_LIST::do_tt() 
{
  for (const_iterator ci=begin(); ci!=end(); ++ci) {
    (*ci)->do_tt();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
const CARD_LIST& CARD_LIST::tr_stress_last(  ) const
{
  for (const_iterator ci=begin(); ci!=end(); ++ci) {
    (*ci)->tr_stress_last();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/* apply things to all cards
 */
const CARD_LIST& CARD_LIST::do_forall( void (CARD::*thing)( ) const  ) const
{
  for (const_iterator ci=begin(); ci!=end(); ++ci) {
    ((*ci)->*thing)( );
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/* apply things to all cards
 */
CARD_LIST& CARD_LIST::do_forall( void (CARD::*thing)( )  )
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    ((*ci)->*thing)( );
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
CARD_LIST& CARD_LIST::do_forall( void (CARD::*thing)( int ), int i  )
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    ((*ci)->*thing)( i );
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/* tr_advance: first pass on a new time step
 */
CARD_LIST& CARD_LIST::tr_advance()
{
  std::list<CARD*>* Q;
  if(!OPT::prequeue) { itested();
    Q = &_cl;
  }else{ itested();
    Q = &_eq;
  }
  for (iterator ci=Q->begin(); ci!=Q->end(); ++ci) {
    trace_func_comp();
    (**ci).tr_advance();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
CARD_LIST& CARD_LIST::tt_accept()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).tt_accept();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
CARD_LIST& CARD_LIST::tt_regress()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).tt_regress();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
CARD_LIST& CARD_LIST::tt_advance()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).tt_advance();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/* tr_regress: throw away the last result and try again, first pass on redo
 */
CARD_LIST& CARD_LIST::tr_regress()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).tr_regress();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/* tr_needs_eval: determine if anything needs to be evaluated
 */
bool CARD_LIST::tr_needs_eval()const
{
  for (const_iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    if ((**ci).tr_needs_eval()) {
      return true;
    }else{
    }
  }
  return false;
}
/*--------------------------------------------------------------------------*/
/* tr_queue_eval: build evaluator queue
 */
CARD_LIST& CARD_LIST::tr_queue_eval()
{
  std::list<CARD*>* Q;
  if(!OPT::prequeue) { itested();
    Q = &_cl;
  }else{ itested();
    Q = &_eq;
  }
  for (iterator ci=Q->begin(); ci!=Q->end(); ++ci) {
    trace_func_comp();
    assert(!OPT::prequeue || !(*ci)->is_constant());
    (**ci).tr_queue_eval();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/* tr_eval: evaluate a list of models
 * evaluates a list (or sublist), checks convergence, etc.
 * does not load the matrix
 * argument is the head of the netlist.
 * recursively called to evaluate subcircuits
 */
bool CARD_LIST::do_tr()
{
  bool isconverged = true;
  if (OPT::bypass) {
    for (iterator ci=begin(); ci!=end(); ++ci) {
      trace_func_comp();
      if ((**ci).tr_needs_eval()) {
	isconverged &= (**ci).do_tr();
      }else{
      }
    }
  }else{
    for (iterator ci=begin(); ci!=end(); ++ci) {
      trace_func_comp();
      isconverged &= (**ci).do_tr();
    }
  }
  return isconverged;
}
/*--------------------------------------------------------------------------*/
/* tr_load: load list of models to the matrix
 * recursively called to load subcircuits
 * Called only when either !OPT::traceload or !SIM::inc_mode
 */
CARD_LIST& CARD_LIST::tr_load()
{
  if (CKT_BASE::_sim->is_inc_mode()) {itested();
    assert(!OPT::traceload);
    for (iterator ci=begin(); ci!=end(); ++ci) {itested();
      trace_func_comp();
      CARD* brh = *ci;
      if (!brh->is_constant()) {itested();
	brh->tr_load();
      }else{itested();
      }
    }
  }else{
    for (iterator ci=begin(); ci!=end(); ++ci) {
      trace_func_comp();
      CARD* brh = *ci;
      brh->tr_load();
    }
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
CARD_LIST& CARD_LIST::tt_begin()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).tt_begin();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
TIME_PAIR CARD_LIST::tt_review()
{
  TIME_PAIR time_by;
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    time_by.min((**ci).tt_review());
  }
  return time_by;
}
/*--------------------------------------------------------------------------*/
TIME_PAIR CARD_LIST::tr_review()
{
  TIME_PAIR time_by(NEVER,NEVER);

  std::list<CARD*>* Q;
  if(!OPT::prequeue) { itested();
    Q = &_cl;
  }else{ itested();
    Q = &_eq;
  }
  for (iterator ci=Q->begin(); ci!=Q->end(); ++ci) {
    trace_func_comp();
    time_by.min((**ci).tr_review());
  }
  return time_by;
}
/*--------------------------------------------------------------------------*/
/* tr_accept: final acceptance of a time step, before moving on
 */
CARD_LIST& CARD_LIST::tr_accept()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).tr_accept();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/* tr_unload: remove a list of models from the matrix
 * recursively called to unload subcircuits
 */
CARD_LIST& CARD_LIST::tr_unload()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).tr_unload();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/* ac_iwant_matrix: allocate solution matrix
 */
CARD_LIST& CARD_LIST::ac_iwant_matrix()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).ac_iwant_matrix();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/* ac_begin: first pass on a new ac simulation
 */
CARD_LIST& CARD_LIST::ac_begin()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).ac_begin();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
CARD_LIST& CARD_LIST::do_ac()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    if (!(**ci).evaluated()) {
      (**ci).do_ac();
    }else{
    }
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
/* ac_load: load list of models to the matrix
 * recursively called to load subcircuits
 */
CARD_LIST& CARD_LIST::ac_load()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    (**ci).ac_load();
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
double CARD_LIST::do_noise() const
{
  double o_power = 0;
  for (const_iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    o_power += (**ci).do_noise();
  }
  return o_power;
}
/*--------------------------------------------------------------------------*/
CARD_LIST& CARD_LIST::do_sens()
{
  for (iterator ci=begin(); ci!=end(); ++ci) {
    trace_func_comp();
    if( (**ci).has_probes() ){
      (**ci).do_sens();
    }
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
// fill _params with params p evaluated in scope.
// this is used to evauluate argument lists
void CARD_LIST::attach_params(PARAM_LIST_BASE* p, const CARD_LIST* scope)
{
  if (p) {
    trace2("CARD_LIST::attach_params", *p, *(scope->params()));
    assert(scope);
    if (_params) {
      // delete _params;
      // _params = NULL;
    }else{
      _params = new PARAM_LIST;
    }
    _params->eval_copy(*p, scope);
    _params->set_try_again(p);
    trace1("CARD_LIST::attach_params done", *_params);
  }else{
  }
}
/*--------------------------------------------------------------------------*/
void CARD_LIST::shallow_copy(const CARD_LIST* p)
{
  assert(p);
  _parent = p;
  for (const_iterator ci = p->begin(); ci != p->end(); ++ci) {
    trace_func_comp();
    if ((**ci).is_device() || dynamic_cast<MODEL_CARD*>(*ci)) {
      trace0("CARD_LIST::shallow_copy cloning " + (*ci)->long_label() );
      CARD* copy = (**ci).clone();
      push_back(copy);
    }else{
      trace0("CARD_LIST::shallow_copy not cloning " + (*ci)->long_label() );
    }
  }
}
/*--------------------------------------------------------------------------*/
// set up the map of external to expanded node numbers
void CARD_LIST::map_subckt_nodes(const CARD* model, const CARD* here)
{
  trace2("CARD_LIST::map_subckt_nodes", here->long_label(), hp(model));
  assert(model);
  assert(model->subckt());
  assert(model->subckt()->nodes());
  assert(here);
  trace1("model nodenames", *(model->subckt()->nodes()));
  trace0("model: "+model->long_label());

  uint_t num_nodes_in_subckt = model->subckt()->nodes()->how_many();
  trace3("CARD_LIST::map_subckt_nodes ", here->long_label(), num_nodes_in_subckt, model->net_nodes());
  uint_t map[num_nodes_in_subckt+1];
  for (uint_t ii = 1;  ii < 1+num_nodes_in_subckt; ++ii) {
    map[ii] = unsigned(-1);
  }
  std::set<unsigned>external; //collect numbers of external nodes
  {
    map[0] = 0;
    // self test: verify that port node numbering is correct
    for (uint_t port = 0; port < model->net_nodes(); ++port) {
      assert(model->n_(port).e_() <= num_nodes_in_subckt);
      //assert(model->n_(port).e_() == port+1);
      trace3("CARD_LIST::map_subckt_nodes ports",
          port, model->n_(port).e_(), here->n_(port).t_());
      assert(model->n_(port).e_()!=INVALID_NODE);
    }
    {
      // take care of the "port" nodes (external connections)
      // map them to what the calling circuit wants
      //
      uint_t i=0;
      for (i=1; i <= model->net_nodes(); ++i) {
        unsigned usernumber = model->n_(i-1).e_();
        external.insert(usernumber);
        assert(usernumber<num_nodes_in_subckt+1);
        map[usernumber] = i; // here->n_(i-1).t_(); // already connected!
        assert(model->n_(i-1).e_() == model->n_(i-1).t_());
      }
      //
      // get new node numbers, and assign them to the remaining
      unsigned labelnumber = 1;
      for (assert(i==model->net_nodes() + 1); i <= num_nodes_in_subckt; ++i) {
	// for each remaining node in card_list
        // these are the internal nodes.

        //the labels are in (*(model->subckt()->nodes())). its hard to tell, which is which.
        //but the order should correspond to i + gap
        while( external.size() && external.find(labelnumber) != external.end() ){
          labelnumber++;
        }
        NODE_BASE* node = (*(model->subckt()->nodes()))[labelnumber];
        IString label = node->short_label();
        if (!node){ untested();
          label = "errornode";
        }else{untested();
	}

        assert(this);
        unsigned k = CKT_BASE::_sim->_total_nodes;
        NODE_MAP* Map = nodes();
        CKT_NODE* n = Map->new_node(label, this); // should increase counter

        if (i==labelnumber){
        }else{ itested();
          // gnucap-geda gets here.
        }
        unsigned newnode = CKT_BASE::_sim->newnode_subckt();
        map[labelnumber] = newnode;
        n->set_user_number(newnode);

        assert (k+1 ==  CKT_BASE::_sim->_total_nodes); USE(k);
        labelnumber++;
      }
      trace1("CARD_LIST::map_subckt_nodes done map", *nodes());
    }
  }
  for (uint_t ii = 1;  ii < 1+num_nodes_in_subckt; ++ii) {
    trace2("map", ii, map[ii]);
  }
  trace0("map done");
  trace1("map done", *(nodes()));
  trace1("map done", *(model->subckt()->nodes()));

  // "map" now contains a translation list,
  // from subckt local numbers to matrix index numbers

  // Device nodes (type node_t) points to the NODE in the parent.
  // Mapping is done in node_t.

  // scan the list, map the nodes
  for (CARD_LIST::iterator ci = begin(); ci != end(); ++ci) {
    if ((**ci).is_device()) { itested();
      const CARD* c = *ci;
      trace2("CARD_LIST::map_subckt_nodes subdevice node ", c->long_label(), here->long_label());
      for (uint_t ii = 0;  ii < (**ci).net_nodes();  ++ii) {
        unsigned n = c->n_(ii).t_();
        // n is the usernumber of the node, where c->n_(ii) connected to.
        // that usernumber is NOT the number in owner->_n[]

        CKT_NODE* nn;
        if( !c->n_(ii).n_() ){ untested();
          error(bDANGER, "somethings wrong " + here->long_label() + " " +
              c->long_label() + ": " + to_string(ii) + "\n");
        } else if( external.find(n) != external.end()) {
          // external
          trace4("CARD_LIST::map_subckt_nodes ext", ii, n, map[n], here->n_( map[n]-1 ).n_());
          assert( map[n] != (uint_t)-1 );
          trace1("CARD_LIST::map_subckt_nodes mod", c->n_(ii).n_());
          c->n_(ii) = here->n_( map[n]-1 );
          trace1("CARD_LIST::map_subckt_nodes now", c->n_(ii).n_());
        } else { // internal node
          IString nodelabel =  c->n_(ii).n_()->short_label() ;
          trace4("CARD_LIST::map_subckt_nodes int", ii, n, nodelabel, here->net_nodes());
        // if a string is the key, a string is the key :|
          nn = prechecked_cast<CKT_NODE*>((*(nodes()))[nodelabel]);
          assert(nn);
          c->n_(ii) = node_t( nn );
        }
        if(!c->n_(ii).n_()) {
          throw Exception(here->long_label() + ": need more nodes");
        }
        assert(c->n_(ii).e_() != INVALID_NODE);
#if 0
      for (int ii = 0;  ii < (**ci).net_nodes();  ++ii) {
	// for each connection node in card
	try{
	  (**ci).n_(ii).map_subckt_node(map, owner);
	}catch(...){
	  delete[] map;
	  throw;
	}
      }
#endif
      }
    }else{ itested();
      assert(dynamic_cast<MODEL_CARD*>(*ci));
    }
  }
}
/*--------------------------------------------------------------------------*/
void CARD_LIST::q_hack(CARD* x)
{ itested();

#ifndef NDEBUG
  for (iterator i=card_list._eq.begin(); i!=card_list._eq.end(); ++i){
    if(*i == x){ unreachable();
      trace1("already toplevel queued", x->long_label());
      exit(1);
    }
  }
#endif

  card_list._eq.push_front(x);
}
/*--------------------------------------------------------------------------*/
///ADP_NODE* CARD_LIST::new_adp_node{
///  assert(d);
///  assert(d->scope());
///
///  NODE_MAP* Map = nodes();
///  assert(Map);
///  ADP_NODE* a Map->new_adp_node(node_name);
///  trace2("CARDLIST::new_adp_node", node_name, _nnn->user_number());
///  assert(_nnn);
///}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
