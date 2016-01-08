/*                                 -*- C++ -*-
 * vim:ts=8:sw=2:et
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
 * aux functions associated with the SIM class
 */
//testing=obsolete
#include "m_wave.h"
#include "e_node.h"
#include "u_nodemap.h"
#include "e_cardlist.h"
#include "u_status.h"
#include "e_subckt.h"
#include "io_misc.h"
using namespace std;
using std::fill_n;
/*--------------------------------------------------------------------------*/
SIM_DATA::SIM_DATA()
  :_time0(0.),
   _Time0(0.),
   _freq(0.),
   _temp_c(0.),
   _damp(0.),
   _dtmin(0.),
   _genout(0.),
   _bypass_ok(true),
   _fulldamp(false),
   _last_time(0.),
   _last_Time(0.),
   _freezetime(false),
   _user_nodes(0),
   _subckt_nodes(0),
   _model_nodes(0),
   _total_nodes(0),
   _jomega(0.,0.),
   _limiting(true),
   _vmax(0.),
   _vmin(0.),
   _age(false),
   _uic(false),
   _cont(false),
   _inc_mode(tsNO),
   //_mode(),
   //_phase(),
   _nm(NULL),
   _i(NULL),
   _v0(NULL),
   _vt1(NULL),
   _ac(NULL),
   _nstat(NULL),
   _aa(),
   _lu(),
   _acx(),
   _eq(),
   _loadq(),
   _acceptq(),
   _evalq1(),
   _evalq2(),
   _late_evalq(),
   _evalq(NULL),
   _evalq_uc(NULL)
{
  _evalq = &_evalq1;
  _evalq_uc = &_evalq2;
  for (unsigned i=0; i<(unsigned)iCOUNT; ++i) {
    _iter[i]=0;
  }
}
/*--------------------------------------------------------------------------*/
SIM_DATA::~SIM_DATA()
{ untested();
  if (_nm) {unreachable();
    delete [] _nm;
    _nm = NULL;
  }else{ untested();
  }
  if (_i) {unreachable();
    delete [] _i;
    _i = NULL;
  }else{ untested();
  }
  if (_v0) {unreachable();
    delete [] _v0;
    _v0 = NULL;
  }else{ untested();
  }
  if (_vt1) {unreachable();
    delete [] _vt1;
    _vt1 = NULL;
  }else{ untested();
  }
  if (_ac) {unreachable();
    delete [] _ac;
    _ac = NULL;
  }else{ untested();
  }
  if (_nstat) {unreachable();
    delete [] _nstat;
    _nstat = NULL;
  }else{ untested();
  }
  if (_vdcstack.size()) {unreachable();
    delete [] _vdcstack.top();
    _vdcstack.pop();
    assert(_vdcstack.empty());
  }else{ untested();
  }
  //assert(_eq.empty()); //not empty means an analysis ended with an unhandled event
			 // could be DC, could be tran with event time past the end
  assert(_loadq.empty());
  assert(_acceptq.empty());
  assert(_evalq1.empty());
  assert(_evalq2.empty());
  assert(_late_evalq.empty());
  assert(_evalq);
  assert(_evalq_uc);
  _evalq = NULL;
  _evalq_uc = NULL;
}
/*--------------------------------------------------------------------------*/
void SIM_DATA::set_limit()
{
  for (uint_t ii = 1;  ii <= _total_nodes;  ++ii) {
    set_limit(_v0[ii]);
  }
}
/*--------------------------------------------------------------------------*/
void SIM_DATA::set_limit(double v)
{
//  trace1("SIM_DATA::set_limit", v);
  if (v+.4 > _vmax) {
    _vmax = v+.5;
    error(bTRACE, "new max = %g, new limit = %g\n", v, _vmax);
  }
  if (v-.4 < _vmin) {
    _vmin = v-.5;
    error(bTRACE, "new min = %g, new limit = %g\n", v, _vmin);
  }
}
/*--------------------------------------------------------------------------*/
void SIM_DATA::clear_limit()
{
  _vmax = OPT::vmax;
  _vmin = OPT::vmin;
}
/*--------------------------------------------------------------------------*/
void SIM_DATA::put_v1_to_v0()
{
  for (uint_t ii = 1;  ii <= _total_nodes;  ++ii) {
    _v0[ii] = _vt1[ii];
  }
}
/*--------------------------------------------------------------------------*/
double SIM_DATA::v0dist()const
{
  double ret = 0;
  double* vdc;
  try{
    vdc = _vdcstack.second();
  }catch (Exception_Cant_Find){ untested();
    return NOT_VALID;
  }
  assert(vdc);

  for (unsigned ii = 1;  ii <= _total_nodes;  ++ii) {
    ret+=(vdc[ii] - _v0[ii]) * (vdc[ii] - _v0[ii]);
  }
  return ret;
}
/*--------------------------------------------------------------------------*/
// compute time that minimizes distance to vdc = _vdcstack.second;
// time is relative to reftime == time1 == time(vt1)
// timescale is "dt0" (after accept!)
double SIM_DATA::v0dist_min(double *what)const
{
  double* vdc;
  try{
    vdc = _vdcstack.second();
  }catch (Exception_Cant_Find){ untested();
    return NOT_VALID;
  }

  double num = 0.;
  double den = 0.;

  for (unsigned ii = 1;  ii <= _total_nodes;  ++ii) {
    // trace4("", ii, _v0[ii], _vt1[ii], vdc[ii]);
    num += (_v0[ii] - _vt1[ii]) * (vdc[ii] - _vt1[ii]);
    den += (_v0[ii] - _vt1[ii]) * (_v0[ii] - _vt1[ii]);
  }
  assert(is_number(num));
  assert(is_number(den));
  double t = num/den;
  trace4("v0dist_min", _time0, num, den, num/den);
  if(den==0){ untested();
    return NOT_VALID;
  }else if(what){
    assert(*what==0);
    for (unsigned ii = 1;  ii <= _total_nodes;  ++ii) {
      double p = t * (_v0[ii] - _vt1[ii]) + _vt1[ii];
      *what += (p - vdc[ii]) * (p - vdc[ii]);
    }
  }else{
  }
  return t;
}
/*--------------------------------------------------------------------------*/
void SIM_DATA::keep_voltages(bool push)
{
  trace2("SIM_DATA::keep_voltages", push, _freezetime);
  assert(_vdcstack.size());
  if(push) {
    double* vdc =  new double[_total_nodes+1+_adp_nodes];
    _vdcstack.push(vdc);
    // _tt = vdc + _total_nodes + 1;
  }else{
  }
  double* vdc = _vdcstack.top();
  if (!_freezetime){
    for (unsigned ii = 1;  ii <= _total_nodes;  ++ii) {
      vdc[ii] = _v0[ii];
    }
    _last_time = (_time0 > 0.) ? _time0 : 0.;
  }else{
    if(push) incomplete();
  }
}
/*--------------------------------------------------------------------------*/
void SIM_DATA::restore_voltages(bool pop)
{
  trace2("SIM_DATA::restore_voltages", _freezetime, pop);
  assert(!_vdcstack.empty());

  double* vdc = _vdcstack.top();
  _nstat[0].set_discont(disNONE);
  for (unsigned ii = 1;  ii <= _total_nodes;  ++ii) {
    _vt1[ii] = _v0[ii] = vdc[ii];
    _nstat[ii].set_discont(disNONE);
    //_nstat[_nm[ii]].set_last_change_time(0);
    //_nstat[_nm[ii]].store_old_last_change_time();
    //_nstat[_nm[ii]].set_final_time(0);
  }
  if(pop){ untested();
    pop_voltages();
  }
}
/*--------------------------------------------------------------------------*/
void SIM_DATA::pop_voltages()
{
  trace1("SIM_DATA::pop_voltages", _vdcstack.size());
  delete[] _vdcstack.top();
  _vdcstack.pop();
  if (_vdcstack.empty()) {
    // unreachable(); reachable if init dc fails. hmmm
    return;
  }
  // _tt = _vdcstack.top() + _total_nodes + 1;
}
/*--------------------------------------------------------------------------*/
// probaly quite stupid thing to do.
void SIM_DATA::zero_currents()
{
  for (uint_t ii = 1;  ii <= _total_nodes;  ++ii) {
    _i[ii] = 0.;
  }
}
/*--------------------------------------------------------------------------*/
void SIM_DATA::zero_some_voltages()
{
  for (uint_t ii = 1;  ii <= _total_nodes;  ++ii) {
    _vt1[ii] = _v0[ii] = _i[ii] = 0.;
  }
}
/*--------------------------------------------------------------------------*/
void SIM_DATA::zero_dc_voltages()
{ untested();
  for (uint_t ii = 1;  ii <= _total_nodes;  ++ii) { untested();
    vdc()[ii] = 0.;
  }
}
/*--------------------------------------------------------------------------*/
void SIM_DATA::zero_voltages()
{
  for (uint_t ii = 1;  ii <= _total_nodes;  ++ii) {
    _vt1[ii] = _v0[ii] = vdc()[ii] = _i[ii] = 0.;
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/* map__nodes: map intermediate node number to internal node number.
 * Ideally, this function would find some near-optimal order
 * and squash out gaps.
 */
void SIM_DATA::map__nodes()
{
  trace2("SIM_DATA::map__nodes", _total_nodes, OPT::order);
  _nm = new unsigned[_total_nodes+1];
  ::status.order.reset().start();
  switch (OPT::order) { //
  default:       unreachable(); error(bWARNING, "invalid order spec: %d\n", OPT::order);
    case oAUTO:		       order_auto();    break;
    case oREVERSE:             order_reverse(); break;
    case oFORWARD:             order_forward(); break;
    case oSINK_F:              order_sink_forward(); break;
    case oSINK_R:              order_sink_reverse(); break;
    case oTREE_DF: untested(); order_tree_df(); break;
    case oTREE_BF:             order_tree_bf(); break;
    case oCOMP:                order_comp(); break;
  }

  _aa.iwant(_nm);
  _lu.iwant(_nm);
  _acx.iwant(_nm);

  ::status.order.stop();
  trace0("map__nodes: done" );

  // delete[_nm]?
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/* order_reverse: force ordering to reverse of user ordering
 *  subcircuits at beginning, results on border at the bottom
 */
void SIM_DATA::order_reverse()
{
  _nm[0] = 0;
  for (uint_t node = 1;  node <= _total_nodes;  ++node) {
    _nm[node] = _total_nodes - node + 1;
  }
}
/*--------------------------------------------------------------------------*/
void SIM_DATA::order_sink_reverse()
{
  _nm[0] = 0;
  // hmm why use _aa?
  _aa.sink_reverse(_nm);
  assert(_nm[0]==0);
}
/*--------------------------------------------------------------------------*/
/* order_forward: use user ordering, with subcircuits added to end
 * results in border at the top (worst possible if lots of subcircuits)
 */
void SIM_DATA::order_forward()
{
  _nm[0] = 0;
  for (uint_t node = 1;  node <= _total_nodes;  ++node) {
    _nm[node] = node;
  }
}
/*--------------------------------------------------------------------------*/
void SIM_DATA::order_sink_forward()
{
  _nm[0] = 0;
  // hmm why use _aa?
  _aa.sink_forward(_nm);
  assert(_nm[0]==0);
}
/*--------------------------------------------------------------------------*/
/* order_auto: full automatic ordering
 * reverse, for now
 */
void SIM_DATA::order_auto()
{
  order_reverse();
}
/*--------------------------------------------------------------------------*/
// component wise node sort. depth-first
void SIM_DATA::order_comp( const CARD_LIST* scope )
{
  static unsigned c;
  static bool* d;
  unsigned t = CKT_BASE::_sim->_total_nodes;
  if (scope==&CARD_LIST::card_list){
    trace2("SIM_DATA::order_comp init", CKT_BASE::_sim->_total_nodes, c );
    assert(!c);
    c = 0;
    d = new bool[t+2]; // +2?
    _nm[0]=0;
    for(unsigned k=1; k<=t ;++k){
      _nm[k]=0;
      d[k]=0;
    }
    d[t+1]=0;
    d[0]=1;
  }

  for (CARD_LIST::const_iterator i = scope->begin(); i != scope->end(); ++i) {
//    trace1("SIM_DATA::order_comp " , (*i)->short_label());
    for (CARD_LIST::const_iterator j = scope->begin(); j != scope->end(); ++j) {
      const BASE_SUBCKT* s = dynamic_cast<const BASE_SUBCKT*>(*j);
      if (s && (*j)->is_device()) {
        order_comp(s->subckt());
      }
    }

    // for(unsigned k=0; k<(*i)->net_nodes();++k){ untested();
    for(int k = int((*i)->net_nodes())-1; k>=0 ;--k) {
      unsigned un = (*i)->n_(unsigned(k)).e_();

      if(!d[un]){
        c++;
        _nm[un] = c;
        d[un] = true;
        trace4("SIM_DATA::order_comp " , (*i)->long_label(), k, un, c);
        assert(un<_total_nodes+1);
      }
    }
  }

  if (scope == &CARD_LIST::card_list){
    assert(c);
    assert(d);
    trace2("SIM_DATA::order_tree", c, t);
    if  (c != CKT_BASE::_sim->_total_nodes ){ untested();
      error(bDANGER, "c=%i, t=%i\n", c, CKT_BASE::_sim->_total_nodes  );
    }
    assert  (c <= CKT_BASE::_sim->_total_nodes );
    delete d;
  }
}
/*--------------------------------------------------------------------------*/
// breadth first tree traversal
// based on NODE_MAP
void SIM_DATA::order_tree_bf( const CARD_LIST* scope)
{
  static unsigned c;
  static bool* d;
  unsigned t = CKT_BASE::_sim->_total_nodes;
  if (scope == &CARD_LIST::card_list){
    assert(c==0);
    assert(d==0);
    for(unsigned k=0; k<=t ;++k){
      _nm[k] = 0;
    }
    d = new bool[t+2]();
  }
  const NODE_MAP * nm = scope->nodes();

  for (NODE_MAP::const_iterator i = nm->begin(); i != nm->end(); ++i) {
    CKT_NODE* s=dynamic_cast<CKT_NODE*>(i->second );
    if (!s) continue;
    if (i->first != "0") {
      unsigned un = s->user_number();
      if(!d[un]){
        c++;
        d[un] = true;
        _nm[un] = c;
      }
      trace3("SIM_DATA::order_tree ", s->long_label(), un, _nm[un]);
    }
  }

  for (CARD_LIST::const_iterator i = scope->begin(); i != scope->end(); ++i) {
    const BASE_SUBCKT* s = dynamic_cast<const BASE_SUBCKT*>(*i);
    if (s && (*i)->is_device()) {
      trace1("SIM_DATA::order_tree child ", s->long_label() );
      order_tree_bf(s->subckt());
    }
  }

  if (scope == &CARD_LIST::card_list){
    trace2("SIM_DATA::order_tree", c,  CKT_BASE::_sim->_total_nodes );
    assert(c <= CKT_BASE::_sim->_total_nodes );
    delete d;
    c=0;
  }
}
/*--------------------------------------------------------------------------*/
// depth-first tree
void SIM_DATA::order_tree_df( const CARD_LIST* scope)
{ untested();
  static unsigned c;
  static bool* d;
  unsigned t = CKT_BASE::_sim->_total_nodes;
  if (scope == &CARD_LIST::card_list){ untested();
    assert(c==0);
    assert(d==0);
    for(unsigned k=0; k<=t ;++k){ untested();
      _nm[k] = 0;
    }
    d = new bool[t+2]();
  }
  const NODE_MAP * nm = scope->nodes();

  for (CARD_LIST::const_iterator i = scope->begin(); i != scope->end(); ++i) { untested();
    const BASE_SUBCKT* s = dynamic_cast<const BASE_SUBCKT*>(*i);
    if (s && (*i)->is_device()) { untested();
      trace1("SIM_DATA::order_tree child ", s->long_label() );
      order_tree_df(s->subckt());
    }
  }

  for (NODE_MAP::const_iterator i = nm->begin(); i != nm->end(); ++i) { untested();
    CKT_NODE* s=dynamic_cast<CKT_NODE*>(i->second );
    if (!s) continue;
    if (i->first != "0") { untested();
      unsigned un = s->user_number();
      if(!d[un]){ untested();
        c++;
        d[un] = true;
        _nm[un] = c;
      }
      trace3("SIM_DATA::order_tree ", s->long_label(), un, _nm[un]);
    }
  }

  if (scope == &CARD_LIST::card_list){ untested();
    trace2("SIM_DATA::order_tree", c,  CKT_BASE::_sim->_total_nodes );
    assert(c <= CKT_BASE::_sim->_total_nodes );
    delete d;
    c=0;
  }
}
/*--------------------------------------------------------------------------
int SIM_DATA::init_node_count( const CARD_LIST* l ) { //int user, int sub, int mod) { untested();
  trace3("SIM_DATA::init_node_count", user, sub, mod);

  _user_nodes=_subckt_nodes=_model_nodes=_adp_nodes=0;

  l->init_node_count( &_user_nodes, &_subckt_nodes, &_model_nodes, &_adp_nodes);

  assert(_subckt_nodes==0);
  assert(_model_nodes==0);

  return (_total_nodes = _user_nodes);
}
--------------------------------------------------------------------------*/
/* init: allocate, set up, etc ... for any type of simulation
 * also called by status and probe for access to internals and subckts
 */
void SIM_DATA::init(bool need_precalc)
{
  if (is_first_expand()) {
    trace2("SIM_DATA::init first", *CARD_LIST::card_list.nodes(), CARD_LIST::card_list.adp_nodes());
    uninit();

    init_node_count(CARD_LIST::card_list.total_nodes(),
                    0, 0,
                    // this is a temporary hack/workaround
                    CARD_LIST::card_list.adp_nodes() );

    trace1("SIM_DATA::init expanding...", _total_nodes);
    CARD_LIST::card_list.expand();
    trace1("SIM_DATA::init expanded", _total_nodes);
    CARD_LIST::card_list.precalc_last();

    _aa.reinit(_total_nodes);
    _lu.reinit(_total_nodes);
    _acx.reinit(_total_nodes);

    { // a hack. iwant_matrix should request the user_number, not m_
      assert(!_nm);
      unsigned id[_total_nodes+1];
      for(unsigned i=0; i<=_total_nodes; i++){
        id[i] = i;
      }
      _nm = &id[0];
      CARD_LIST::card_list.map_nodes(); //tell nodes about their intermediate number
      _nm = 0; 
    }

    CARD_LIST::card_list.tr_iwant_matrix();
    CARD_LIST::card_list.ac_iwant_matrix();

    map__nodes(); //calculate map.
    assert(_nm[0]==0);
    CARD_LIST::card_list.map_nodes(); //tell nodes about their final number
    alloc_hold_vectors();

    _last_time = 0;
    _last_Time = 0;
  }else if(need_precalc){
    CARD_LIST::card_list.precalc_first();
    CARD_LIST::card_list.precalc_last();
  }
  _tt_iter=0;
  _dT0=0;
  _dT1=0;
  _dT2=0;
  _expect_file=0;
  //_waves=0;
  //_waves_tt=0;
}
/*--------------------------------------------------------------------------*/
/* alloc_hold_vectors:
 * allocate space to hold data between commands.
 * for restart, convergence assistance, bias for AC, post-processing, etc.
 * must be done BEFORE deciding what array elements to allocate,
 * but after mapping
 * if they already exist, leave them alone to save data
 */
void SIM_DATA::alloc_hold_vectors()
{
  assert(is_first_expand());

  assert(!_nstat);
  _nstat = new LOGIC_NODE[_total_nodes+1];
  for (unsigned ii=0;  ii <= _total_nodes;  ++ii) {
    _nstat[_nm[ii]].set_user_number(ii);
    assert(!_nstat[ii].discont());
  }

  assert(_vdcstack.empty());
  double* vdc =  new double[_total_nodes+1+_adp_nodes];
  // _tt = vdc + _total_nodes + 1;
  _tt = new double[_adp_nodes];
  _vdcstack.push(vdc);
  std::fill_n(vdc, _total_nodes+1, 0);

  std::fill_n(_tt, _adp_nodes, 0);

  assert(_nstat);
}
/*--------------------------------------------------------------------------*/
/* alloc_vectors:
 * these are new with every run and are discarded after the run.
 */
void SIM_DATA::alloc_vectors()
{
  trace1("SIM_DATA::alloc_vectors",  _total_nodes);
  assert(_evalq1.empty());
  if (! _evalq2.empty() ){ trace1("SIM_DATA::alloc_vectors", _evalq2);}
  assert(_evalq2.empty());
  assert(_evalq != _evalq_uc);

  assert(!_ac);
  assert(!_i);
  assert(!_v0);
  assert(!_vt1);

  _ac = new COMPLEX[_total_nodes+1];
  _sens = new COMPLEX[_total_nodes+1];
  _i   = new double[_total_nodes+1];
  _v0  = new double[2*_total_nodes+2];
  _vt1 = &_v0[_total_nodes+1];
  // _sens = (COMPLEX*) _v0;
	
  trace1("SIM_DATA::alloc_vectors ", hp(_i));
  _tr  = new double[_adp_nodes];

  _tr1 = new double[_adp_nodes];
  _tr2 = new double[_adp_nodes];
  _tr3 = new double[_adp_nodes];
  _tt1 = new double[_adp_nodes];

  std::fill_n(_ac, _total_nodes+1, 0);
  std::fill_n(_sens, _total_nodes+1, 0);
  std::fill_n(_i,  _total_nodes+1, 0);
  std::fill_n(_v0, _total_nodes+1, 0);
  std::fill_n(_vt1,_total_nodes+1, 0);

  invalidate_tt();
}

/*--------------------------------------------------------------------------*/
// debugging flow
void SIM_DATA::invalidate_tt()
{
  std::fill_n(_tr, _adp_nodes, NAN); 
#ifndef NDEBUG
  std::fill_n(_tr1, _adp_nodes, NAN); 
  std::fill_n(_tr2, _adp_nodes, NAN); 
  std::fill_n(_tr3, _adp_nodes, NAN); 
  std::fill_n(_tt1, _adp_nodes, NAN); 
#else
  // BUG. valgrind unhappy in pre-simulation print
#endif
}
/*--------------------------------------------------------------------------*/
void SIM_DATA::unalloc_vectors()
{
  trace1("SIM_DATA::unalloc_vectors",_total_nodes);
  _evalq1.clear();
  _evalq2.clear();
  delete [] _i;
  _i = NULL;
  delete [] _v0;
  _v0 = NULL;
  // delete [] _vt1;
  _vt1 = NULL;
  delete [] _ac;
  _ac = NULL;
  delete [] _sens;
  _sens = NULL;
  delete [] _tr;
  delete [] _tr1;
  delete [] _tr2;
  delete [] _tr3;
  delete [] _tt1;
  _tr=_tr1=_tt1=_tr2=_tr3=NULL;
}
/*--------------------------------------------------------------------------*/
/* uninit: undo all the allocation associated with any simulation
 * called when the circuit changes after a run, so it needs a restart
 * may be called multiple times without damage to make sure it is clean
 */
void SIM_DATA::uninit()
{
  trace1("SIM_DATA::uninit", _vdcstack.size());
  // fixme adp?
  if (_vdcstack.size()) {
    _acx.reinit(0);
    _lu.reinit(0);
    _aa.reinit(0);
    double* vdc = _vdcstack.top();
    delete [] vdc;
    _vdcstack.pop();
    while (_vdcstack.size()) {
      _vdcstack.pop();
    }
    assert(!_vdcstack.size());
    delete [] _tt;
    _tt = NULL;
    delete [] _nstat;
    _nstat = NULL; // triggers first_expand.
    delete [] _nm;
    _nm = NULL;
  }else{
    assert(_acx.size() == 0);
    assert(_lu.size() == 0);
    assert(_aa.size() == 0);
    assert(!_nstat);
    assert(!_nm);
  }
}
/*--------------------------------------------------------------------------*/
void SIM_DATA::update_tt_order()
{
  uint_t new_order=3;
  if(_dT2 == .0 ) new_order=2;
  if(_dT1 == .0 ) new_order=1;
  if(_dT0 == .0 ) new_order=0; // ?
  if( OPT::ttsteporder < new_order ) new_order = OPT::ttsteporder;

  if (_tt_order != new_order ) {
    trace3("SIM_DATA::update_tt_order ", tt_iteration_number(),  new_order, _tt_order );
    _tt_order=min(new_order, _tt_order+1);
  } else {
    trace3("SIM_DATA::update_tt_order unchanged.", _tt_order, _dT0, _dT1);
  }
}
/*--------------------------------------------------------------------------*/
// FIXME: individualize (like trsteporder)
uint_t SIM_DATA::get_tt_order() const {
  trace2("SIM_DATA::get_tt_order", iteration_number(), tt_iteration_number());
  assert (_tt_order <= tt_iteration_number());
  return min(_tt_order, OPT::ttsteporder);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
