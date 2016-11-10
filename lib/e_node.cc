/*                                   -*- C++ -*-
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
 * node probes
 */
#include "u_nodemap.h"
#include "d_logic.h"
#include "e_aux.h"
#include "u_xprobe.h"
#include "e_adp.h"
#include "e_node.h"
#include "io_misc.h"
/*--------------------------------------------------------------------------*/
const _LOGICVAL LOGICVAL::or_truth[lvNUM_STATES][lvNUM_STATES] = {
  {lvSTABLE0, lvRISING,  lvFALLING, lvSTABLE1, lvUNKNOWN},
  {lvRISING,  lvRISING,  lvRISING,  lvSTABLE1, lvRISING},
  {lvFALLING, lvRISING,  lvFALLING, lvSTABLE1, lvUNKNOWN},
  {lvSTABLE1, lvSTABLE1, lvSTABLE1, lvSTABLE1, lvSTABLE1},
  {lvUNKNOWN, lvRISING,  lvUNKNOWN, lvSTABLE1, lvUNKNOWN}
};
/*--------------------------------------------------------------------------*/
const _LOGICVAL LOGICVAL::xor_truth[lvNUM_STATES][lvNUM_STATES] = {
  {lvSTABLE0, lvRISING,  lvFALLING, lvSTABLE1, lvUNKNOWN},
  {lvRISING,  lvFALLING, lvRISING,  lvFALLING, lvUNKNOWN},
  {lvFALLING, lvRISING,  lvFALLING, lvRISING,  lvUNKNOWN},
  {lvSTABLE1, lvFALLING, lvRISING,  lvSTABLE0, lvUNKNOWN},
  {lvUNKNOWN, lvUNKNOWN, lvUNKNOWN, lvUNKNOWN, lvUNKNOWN}
};
/*--------------------------------------------------------------------------*/
const _LOGICVAL LOGICVAL::and_truth[lvNUM_STATES][lvNUM_STATES] = {
  {lvSTABLE0, lvSTABLE0, lvSTABLE0, lvSTABLE0, lvSTABLE0},
  {lvSTABLE0, lvRISING,  lvFALLING, lvRISING,  lvUNKNOWN},
  {lvSTABLE0, lvFALLING, lvFALLING, lvFALLING, lvFALLING},
  {lvSTABLE0, lvRISING,  lvFALLING, lvSTABLE1, lvUNKNOWN},
  {lvSTABLE0, lvUNKNOWN, lvFALLING, lvUNKNOWN, lvUNKNOWN}
};
/*--------------------------------------------------------------------------*/
const _LOGICVAL LOGICVAL::not_truth[lvNUM_STATES] = {
  lvSTABLE1, lvFALLING, lvRISING,  lvSTABLE0, lvUNKNOWN  
};
/*--------------------------------------------------------------------------*/
static _LOGICVAL prop_truth[lvNUM_STATES][lvNUM_STATES] = {
  {lvSTABLE0, lvUNKNOWN, lvUNKNOWN, lvRISING,  lvUNKNOWN},
  {lvUNKNOWN, lvUNKNOWN, lvUNKNOWN, lvUNKNOWN, lvUNKNOWN},
  {lvUNKNOWN, lvUNKNOWN, lvUNKNOWN, lvUNKNOWN, lvUNKNOWN},
  {lvFALLING, lvUNKNOWN, lvUNKNOWN, lvSTABLE1, lvUNKNOWN},
  {lvFALLING, lvUNKNOWN, lvUNKNOWN, lvRISING,  lvUNKNOWN}
};
/*--------------------------------------------------------------------------*/
inline LOGICVAL& LOGICVAL::set_in_transition(LOGICVAL newval)
{
  trace3("LOGICVAL::set_in_transition", _lv ,newval,  prop_truth[_lv][newval]);
     
  LOGICVAL a(_lv);
  _lv = prop_truth[_lv][newval];
  if(_lv == lvUNKNOWN){
    error(bWARNING, "set_in_transition: lv unknown... newval %i from %i\n",
        (int)newval, (int)a );
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
LOGIC_NODE::LOGIC_NODE()
  :NODE(),
   _family(0),
   _d_iter(0),
   _a_iter(0),
   _final_time(0),
   _lastchange(0),
   _old_lastchange(0),
   _mode(moANALOG),
   _lv(),
   _old_lv(),
   _quality(qBAD),
   _failure_mode("initial")
{
}
/*--------------------------------------------------------------------------*/
/* default constructor : unconnected, don't use
 */
NODE_BASE::NODE_BASE() 
  : CKT_BASE(),
  _owner(0),
  _scope(0),
  _next(this),
  _user_number(INVALID_NODE)
{
}
/*--------------------------------------------------------------------------*/
NODE::NODE()
  :NODE_BASE(), _discont(disNONE)
{
}
/*--------------------------------------------------------------------------*/
/* copy constructor : user data only
 */
NODE_BASE::NODE_BASE(const NODE_BASE& p)
  : CKT_BASE(p),
  _owner(p._owner),
  _scope(p._scope),
  _next(this),
  _user_number(p._user_number)
   //_flat_number(p._flat_number)
   //_matrix_number(INVALID_NODE)
{
  assert( p._next == &p );
  unreachable();
}
/*--------------------------------------------------------------------------*/
CKT_NODE::CKT_NODE(const CKT_NODE& p)
  :NODE_BASE(p)
{
  unreachable();
}
/*--------------------------------------------------------------------------*/
NODE_BASE::NODE_BASE(const std::string& s, unsigned n, const CARD_LIST* p)
  :CKT_BASE(s),
   _owner(0),
   _scope(p),
   _next(this),
   _user_number(n)
   //_flat_number(n)
   //_matrix_number(INVALID_NODE)
{
#ifndef NDEBUG
  // has already been caught and reported.
  std::string::size_type first = s.find_first_of(".");
  std::string::size_type last = s.find_last_of(".");
  if(last != std::string::npos && last != first) {
    // unreachable(); coil has no sckt for branch node...
  }
#endif

  trace1("NODE_BASE::NODE_BASE()" + s, n);
}
/*--------------------------------------------------------------------------*/
/* constructor taking a pointer : it must be valid
 * supposedly not used, but used by a required function that is also not used
 */
NODE::NODE(const NODE* p)
  :NODE_BASE(*p)
{
  unreachable();
}
/*--------------------------------------------------------------------------*/
/* usual initializing constructor : name and index
 */
CKT_NODE::CKT_NODE(const string& s, unsigned n, const CARD_LIST*p) : NODE_BASE(s,n,p) { }
/*--------------------------------------------------------------------------*/
node_t::node_t()
  :_nnn(0),
   _ttt(INVALID_NODE),
   _m(INVALID_NODE),
   _disc(d_electrical)
{ }
/*--------------------------------------------------------------------------*/
node_t::node_t(const node_t& p)
  :_nnn(p._nnn),
   _ttt(p._ttt),
   _m(p._m),
   _disc(p._disc)
{
  trace1("node_t::node_t cloning", _ttt);
  //assert(_ttt == _nnn->flat_number());
}
/*--------------------------------------------------------------------------*/
node_t::node_t(NODE* n)
  :_nnn(n),
   _ttt(n->user_number()),
   _m(to_internal(n->user_number())),
   _disc(d_electrical) // incomplete. use n->discipline?
                       // (dont store disc at all (use _nnn instead)?)
{
  trace3("node_t::node_t(NODE*) " , _nnn->long_label(), n->user_number(), to_internal(n->user_number()));
}
/*--------------------------------------------------------------------------*/
node_t& node_t::operator=(const node_t& p)
{
  _disc = p._disc;
  if (p._nnn) {
    //assert(p._ttt == p._nnn->flat_number());
  }else if (is_adp()) {
  }else{
    assert(p._ttt == INVALID_NODE);
    assert(p._m   == INVALID_NODE);
  }
  _nnn   = p._nnn;
  _ttt = p._ttt;
  _m   = p._m;
  return *this;
}
/*--------------------------------------------------------------------------*/
LOGIC_NODE& node_t::data()const
{
  assert(CKT_BASE::_sim->_nstat);
  return CKT_BASE::_sim->_nstat[m_()];
}
/*--------------------------------------------------------------------------*/
double NODE_BASE::tr_probe_num(const std::string& x)const
{
  if (Umatch(x, "v ")) {
    // return v0(); denoised
    return floor(v0()/OPT::vfloor + .5) * OPT::vfloor;
  }else{
    return CKT_BASE::tr_probe_num(x);
  }
}
/*--------------------------------------------------------------------------*/
double NODE::tr_probe_num(const std::string& x)const
{
  if (Umatch(x, "v1 ")) {
    return floor(vt1()/OPT::vfloor + .5) * OPT::vfloor;
  }else if (Umatch(x, "z ")) {
    return port_impedance(node_t(const_cast<NODE*>(this)), node_t(&ground_node), _sim->_lu, 0.);
  }else if (Umatch(x, "l{ogic} |la{stchange} |fi{naltime} |di{ter} |ai{ter} |count ")) {
    assert(_sim->_nstat);
    return _sim->_nstat[matrix_number()].tr_probe_num(x);
  }else if (Umatch(x, "mdy ")) {
    // matrix diagonal admittance
    const BSMATRIX<double>&  aaa = _sim->_aa;
    return aaa.d(m_(),m_());
  }else if (Umatch(x, "mdz ")) {
    // matrix diagonal impedance
    const BSMATRIX<double>&  aaa = _sim->_aa;
    return 1/aaa.d(m_(),m_());
  }else if (Umatch(x, "zero ")) {
    // fake probe: 0.0
    return 0.0;
  }else if (Umatch(x, "pdz ")) {
    // fake probe 1/0 .. positive divide by zero = Infinity
    double z1 = tr_probe_num("zero ");
    return 1.0/z1;
  }else if (Umatch(x, "ndz ")) {
    // fake probe -1/0 .. negative divide by zero = -Infinity
    double z1 = tr_probe_num("zero ");
    return -1.0/z1;
  }else if (Umatch(x, "dv ")) { // differential of v
      return ( vt1() );
  }else if (Umatch(x, "ddv ")) { // divided difference v
    double val = 0.; 
    val = ( vdc() - v0() ) / _sim->_dt0;
    return floor(val/OPT::vfloor + .5) * OPT::vfloor;
  }else if (Umatch(x, "nan ")) {
    // fake probe 0/0 = NaN
    double z1 = tr_probe_num("zero ");
    double z2 = tr_probe_num("zero ");
    return z1/z2;
  }else if (Umatch(x, "dis{cont} ")) {
    // fake probe 0/0 = NaN
//    assert((*this)->discont() == _sim->_nstat[matrix_number()].discont());
    return (unsigned int) _sim->_nstat[matrix_number()].discont();
#ifndef NDEBUG
  }else if (Umatch(x, "n ")) {
    return  matrix_number();
  }else if (Umatch(x, "m ")) {
    return  m_();
#endif
  }else if (Umatch(x, "dis{cont} ")) {
    return _sim->_nstat[matrix_number()]._discont;
  }else{
    return NODE_BASE::tr_probe_num(x);
  }
}
/*--------------------------------------------------------------------------*/
double LOGIC_NODE::tr_probe_num(const std::string& x)const
{
  if (Umatch(x, "l{ogic} ")) {
    return annotated_logic_value();
  }else if (Umatch(x, "la{stchange} ")) {
    return _lastchange;
  }else if (Umatch(x, "fi{naltime} ")) {
    return final_time();
  }else if (Umatch(x, "di{ter} ")) {
    return static_cast<double>(_d_iter);
  }else if (Umatch(x, "ai{ter} ")) {
    return static_cast<double>(_a_iter);
  }else{
    return NODE_BASE::tr_probe_num(x);
  }
}
/*--------------------------------------------------------------------------*/
const std::string NODE_BASE::long_label()const
{
  string ret(short_label());
  if (_scope){
    if( _scope->owner()){
      return (_scope->owner()->long_label() + "." + ret);
    }
    return (ret);
  } else if(_owner) {
    if(dynamic_cast<const ADP_NODE*>(this)){
      return (_owner->long_label() + "." + ret);
    }
  } else {
    trace1("NODE_BASE::long_label, have no parent", short_label());
  }
  return ret;
}
/*--------------------------------------------------------------------------*/
XPROBE NODE::ac_probe_ext(const std::string& x)const
{
  if (Umatch(x, "v ")) {
    return XPROBE(vac());
  }else if (Umatch(x, "n ")) {
    return XPROBE(port_noise(node_t(const_cast<NODE*>(this)),node_t(&ground_node)));
  }else if (Umatch(x, "z ")) {
    return XPROBE(port_impedance(node_t(const_cast<NODE*>(this)),
				 node_t(&ground_node), _sim->_acx, COMPLEX(0.)));
  }else{itested();
    return CKT_BASE::ac_probe_ext(x);
  }
}
/*--------------------------------------------------------------------------*/
/* annotated_logic_value:  a printable value for probe
 * that has secondary info encoded in its fraction part
 */
double LOGIC_NODE::annotated_logic_value()const
{
  return (_lv + (.1 * (OPT::transits - quality())) + (.01 * (2 - _mode)));
}
/*--------------------------------------------------------------------------*/
static bool newly_stable[lvUNKNOWN+1][lvUNKNOWN+1] = { // oldlv, _lv
  /*	   s0	  rise   fall	s1     u */
  /* s0 */{false, false, false, true,  false},
  /*rise*/{false, false, false, true,  false},
  /*fall*/{true,  false, false, false, false},
  /* s1 */{true,  false, false, false, false},
  /* u  */{true,  false, false, true,  false}
};
/*--------------------------------------------------------------------------*/
inline bool LOGIC_NODE::just_reached_stable()const
{
  return newly_stable[old_lv()][lv()];
}
/*--------------------------------------------------------------------------*/
/* to_logic: set up logic data for a node, if needed
 * If the logic data is already up to date, do nothing.
 * else set up: logic value (_lv) and quality.
 * Use and update _d_iter, _lastchange to keep track of what was done.
 */
void LOGIC_NODE::to_logic(const MODEL_LOGIC*f)
{
  assert(f);
  if (process() && process() != f->logic_hash()) {untested();
    set_bad_quality("logic process mismatch");
    error(bWARNING, "node " + long_label() 
	  + " logic process mismatch\nis it %s " + 
	   " or " + f->long_label() + "?\n");
  }
  set_process(f);

  if (is_analog() &&  d_iter() < a_iter()) {
    if (_sim->analysis_is_restore()) {untested();
    }else if (_sim->analysis_is_static()) {
    }else{
    }
    if (_sim->analysis_is_static() || _sim->analysis_is_restore()) {
      set_last_change_time(0);
      store_old_last_change_time();
      set_lv(lvUNKNOWN);
    }else{
    }
    double dt = _sim->_time0 - last_change_time();
    if (dt < 0.) {untested();
      error(bPICKY, "time moving backwards.  was %g, now %g, %g\n",
	    last_change_time(), _sim->_time0, _sim->_Time0);
      dt = _sim->_time0 - old_last_change_time();
      if (dt <= 0.) {untested();
	throw Exception("internal error: time moving backwards, can't recover " + long_label());
      }else{untested();
      }
      assert(dt > 0.);
      restore_lv();			/* skip back one */
    }else{
      store_old_last_change_time();
      store_old_lv();			/* save to see if it changes */
    }
    
    // MIXED_NODE?
    double sv = v0() / f->range;	/* new scaled voltage */
    if (sv >= f->th1) {		/* logic 1 */
      switch (lv()) {
      case lvSTABLE0: dont_set_quality("stable 0 to stable 1");	break;
      case lvRISING:  dont_set_quality("begin stable 1");	break;
      case lvFALLING:untested();set_bad_quality("falling to stable 1"); break;
      case lvSTABLE1: dont_set_quality("continuing stable 1");	break;
      case lvUNKNOWN: set_good_quality("initial 1");		break;
      }
      set_lv(lvSTABLE1);
    }else if (sv <= f->th0) {	/* logic 0 */
      switch (lv()) {
      case lvSTABLE0: dont_set_quality("continuing stable 0");	break;
      case lvRISING: untested();set_bad_quality("rising to stable 0");	break;
      case lvFALLING: dont_set_quality("begin stable 0");	break;
      case lvSTABLE1: dont_set_quality("stable 1 to stable 0");	break;
      case lvUNKNOWN: set_good_quality("initial 0");		break;
      }
      set_lv(lvSTABLE0);
    }else{				/* transition region */
      double oldsv = vt1() / f->range;/* old scaled voltage */
      double diff  = sv - oldsv;
      if (diff > 0) {	/* rising */
	switch (lv()) {
	case lvSTABLE0:
	  dont_set_quality("begin good rise");
	  break;
	case lvRISING:
	  if (diff < dt/(f->mr * f->rise)) {
	    set_bad_quality("slow rise");
	  }else{
	    dont_set_quality("continuing good rise");
	  }
	  break;
	case lvFALLING:
	  untested();
	  set_bad_quality("positive glitch in fall");
	  break;
	case lvSTABLE1:
	  untested();
	  set_bad_quality("negative glitch in 1");
	  break;
	case lvUNKNOWN:
	  set_bad_quality("initial rise");
	  break;
	}
	set_lv(lvRISING);
      }else if (diff < 0) {	/* falling */
	switch (lv()) {
	case lvSTABLE0:
	  untested();
	  set_bad_quality("positive glitch in 0");
	  break;
	case lvRISING:
	  set_bad_quality("negative glitch in rise");
	  break;
	case lvFALLING:
	  if (-diff < dt/(f->mf * f->fall)) {
	    set_bad_quality("slow fall");
	  }else{
	    dont_set_quality("continuing good fall");
	  }
	  break;
	case lvSTABLE1:
	  dont_set_quality("begin good fall");
	  break;
	case lvUNKNOWN:
	  untested();
	  set_bad_quality("initial fall");
	  break;
	}
	set_lv(lvFALLING);
      }else{				/* hanging up in transition */
	untested();
	error(bDANGER, "inflection???\n");
	set_bad_quality("in transition but no change");
	/* state (rise/fall)  unchanged */
      }
    }
    if (sv > 1.+f->over || sv < -f->over) {
      trace2("out of range", sv, f->over);
      // out of range
      set_bad_quality("out of range");
    }
    if (just_reached_stable()) { /* A bad node gets a little better */
      improve_quality();	/* on every good transition.	   */
    }				/* Eventually, it is good enough.  */
				/* A good transition is defined as */
				/* entering a stable state from    */
				/* a transition state.		   */
    set_d_iter();
    set_last_change_time();
    trace3(_failure_mode.c_str(), _lastchange, _quality, _lv);
  }
}
/*--------------------------------------------------------------------------*/
void LOGIC_NODE::set_process(const MODEL_LOGIC* f) {_family = f->logic_hash();}
/*--------------------------------------------------------------------------*/
double LOGIC_NODE::to_analog(const MODEL_LOGIC* f)
{
  assert(f);
  if (process() && process() != f->logic_hash()) {untested();
    error(bWARNING, "node " + long_label() 
	  + " logic process mismatch\nis it %i"
	  + " or " + f->long_label() + "?\n");
  }
  set_process(f);

  double start = NOT_VALID;
  double end = NOT_VALID;
  double del = NOT_VALID; // the analog transition will take this much longer,
                          // until final_time() + del...
  double risefall = NOT_VALID;

  switch (lv()) {
  case lvRISING:
    risefall = f->rise;
    del = risefall * (1 - f->th1);
    set_final_time_a(final_time()+ del);
  case lvSTABLE1:
    risefall = f->rise;
    start = f->vmin;
    end = f->vmax;
    //end = f->vmax;
    break;
  case lvFALLING:
    risefall = f->fall;
    del = risefall *  f->th0;
    set_final_time_a(final_time()+ del);
  case lvSTABLE0:
    risefall = f->fall;
    start = f->vmax;
    end = f->vmin;
    break;
  case lvUNKNOWN:
    set_final_time_a(NEVER);
    return f->unknown;
  }

  if(_sim->_time0 > final_time_a()){
    return end;
  }

  assert(start != NOT_VALID);
  assert(end != NOT_VALID);
  assert(risefall != NOT_VALID);

  if (_sim->_time0 <= (final_time_a()-risefall)) {
    trace1("", final_time_a());
    return start;
  }else if (_sim->_time0 >= final_time_a()) {
    return end;
  }else{
    double share = (final_time_a() - _sim->_time0) / risefall;
    trace4("LOGIC_NODE::to_analog in between", _sim->_time0, final_time(),
        risefall, share );
    double ret = end - (end-start) * share;
    return ret;
  }
}
/*--------------------------------------------------------------------------*/
void LOGIC_NODE::propagate()
{
  assert(in_transit());
  if (lv().is_rising()) {
    set_lv(lvSTABLE1);
  }else if (lv().is_falling()) {
    set_lv(lvSTABLE0);
  }else{
    // lv no change
  }
  set_d_iter();
  set_final_time(NEVER);
  set_last_change_time();
  assert(!(in_transit()) || final_time_a() < NEVER);
}
/*--------------------------------------------------------------------------*/
void LOGIC_NODE::force_initial_value(LOGICVAL v)
{
  trace2("LOGIC_NODE::force_initial_value "+ short_label(), _mode, v);
  if (_sim->analysis_is_restore()) {untested();
  }else if (_sim->analysis_is_static()) {
  }else{untested();
  }
  assert(_sim->analysis_is_static() || _sim->analysis_is_restore());
  assert(_sim->_time0 == 0.);
  assert(is_unknown());
  assert(is_digital()); // fails e.g. if outputs are short.
  set_lv(v); // BUG ??
  set_good_quality("initial dc");
  set_d_iter();
  set_final_time(NEVER);
  set_final_time_a(0); // analog transition has just ended... (good idea?)
  set_last_change_time();
}
/*--------------------------------------------------------------------------*/
bool LOGIC_NODE::in_transit()const
{
  return (final_time() < NEVER) ; // || (final_time_a() < NEVER);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
void LOGIC_NODE::set_event_abs(double time, LOGICVAL v)
{
  trace2("LOGIC_NODE::set_event_abs", time, v);
  _lv.set_in_transition(v);
  if (_sim->analysis_is_tran_dynamic()  &&  in_transit()) {untested();
    set_bad_quality("race");
  }
  set_d_iter();
  set_final_time(time);

  /*
  double del=0;
  switch (lv()) {
  case lvSTABLE1:
  case lvRISING:
    del = f->rise * (1 - f->th1);
    break;
  case lvSTABLE0:
  case lvFALLING:
    del = f->fall * (f->th0);
    break;
  default:
  }
  */


  if (OPT::picky <= bTRACE) {untested();
    error(bTRACE, "%s:%u:%g new event\n",
	  long_label().c_str(), d_iter(), final_time());
  }
  set_last_change_time();
}
/*--------------------------------------------------------------------------*/
void LOGIC_NODE::set_event(double delay, LOGICVAL v)
{
  return set_event_abs(delay+_sim->_time0,v);
}
/*--------------------------------------------------------------------------*/
void node_t::set_to_ground(CARD* d)
{
  assert( is_electrical() );
  //assert(!_nnn); //BUG// fails on MUTUAL_L::expand after clone
  assert(d);

  NODE_MAP* Map = d->scope()->nodes();
  assert(Map);
  _nnn = dynamic_cast<CKT_NODE*>((*Map)["0"]);
  _ttt = 0;
  assert(_nnn);
}
/*--------------------------------------------------------------------------*/
size_t NODE_BASE::how_many()const
{
  unsigned ret=1;
  const NODE_BASE* n = this;
  while( n->_next != this){
    trace1("how_many", n->long_label());
    ret++;
    n = n->_next;
  }
  return ret;
}
/*--------------------------------------------------------------------------*/
void NODE_BASE::collapse(const NODE_BASE& to) const
{
  if(user_number()==to.user_number()){untested();
    return;
  }
  const NODE_BASE* n = this;
  trace4("NODE_BASE::collapse", n, to.user_number(), to.how_many(), how_many());
  set_user_number(to.user_number());
  while( n->_next != this){
    trace2("NODE_BASE::collapse", n, to.user_number());
    n = n->_next;
    n->set_user_number(to.user_number());
  }
  n->_next = to._next;
  to._next = this;
}
/*--------------------------------------------------------------------------*/
void node_t::collapse(CARD* d, const node_t to)
{
  trace3("node_t::collapse", d->long_label(), _ttt, to._ttt);
  // assert(_m == INVALID_NODE); // not true in second init()
  trace2("node_t::collapse", _nnn->user_number(), to.e_());
  assert(d);

  NODE_MAP* Map = d->scope()->nodes();
  assert(Map); USE(Map);
  assert(_nnn);
  // _ttt = to._ttt;
//  const NODE_BASE* n = to._nnn->user_node();
  if(_nnn->user_number()){
    _nnn->collapse( *(to.n_()) );
  } else {
    to.n_()->collapse(*_nnn);
  }
  // _nnn = to._nnn
  trace5("node_t::collapse", d->long_label(), _ttt, to._ttt, _nnn->short_label(), to.n_()->short_label());
  assert(_nnn);
  // assert(e_() == to.e_());
  // assert(t_() == to.t_());
}
/*--------------------------------------------------------------------------*/
/* new_node: a raw new node, as when a netlist is parsed
 */
void node_t::new_node(const std::string& node_name, const CARD* d)
{
  trace2("node_t::new_node", node_name, d->long_label());
  return new_node(node_name, d->scope());
}
/*--------------------------------------------------------------------------*/
// same, but leave choice of scope to user.
void node_t::new_node(const std::string& node_name, const CARD_LIST* scope)
{
  if (is_adp()){
    untested();
    return new_adp_node(node_name, scope);
  }
  assert( is_electrical() );
  //assert(!_nnn); //BUG// fails on MUTUAL_L::expand after clone
  assert(scope);

  NODE_MAP* Map = scope->nodes();
  assert(Map);
  if(!_nnn) {
  }else if(_nnn != (*Map)[node_name]){ untested();
    error(bWARNING, "%s is already a node: %s\n", node_name.c_str(), _nnn->long_label().c_str());
  }

  _nnn = (CKT_NODE*) Map->new_node(node_name, scope);
  _ttt = _nnn->user_number();
  assert(_nnn);
}
/*--------------------------------------------------------------------------*/
void node_t::new_adp_node(const std::string& node_name, const CARD_LIST* scope)
{
  // FIXME: remove, implement more generic new_node
  assert(scope);

  NODE_MAP* Map = scope->nodes();
  assert(Map);
  _nnn = (CKT_NODE*) Map->new_adp_node(node_name, scope);
  assert(_nnn);
  trace2("node_t::new_adp_node", node_name, _nnn->user_number());
  _ttt = _nnn->user_number();
  assert(_nnn);
}
/*--------------------------------------------------------------------------*/
/* new_model_node: a mapped new node, produced through model expansion.
 * Not really a model_node, but a node in the subckt that is made
 * in model expansion.
 * Supposedly equivalent to new_node() then map_subckt_node()
 * but it does it without building a map
 */
void node_t::new_model_node(const std::string& s_in, CARD* d)
{
  string s = s_in;
  std::string::size_type dotplace = s_in.find_last_of(".");
  if(dotplace != std::string::npos && s_in.c_str()[0]!='.'){
    incomplete(); // this is really dangerous...
                  // sometimes dots are there "intentionally":
                  // - hack in spice wrapper
                  // - hack in d_coil.cc
    trace1("node_t::new_model_node too long (BUG).", s_in);
    s = s_in.substr(dotplace+1, std::string::npos);
  }

  if (is_adp()){
    untested();
    return new_model_adp_node(s,d);
  }
  assert( is_electrical() );

  if (d->subckt()){
    new_node(s, d->subckt());
  } else { // happens in spice wrapper and in inductance.
    assert(d->scope());
    new_node(s, d->scope());
  }
  _ttt = CKT_BASE::_sim->newnode_model(); // increase global counter.
  _nnn->set_user_number(_ttt);
  assert(_ttt == _nnn->user_number());
}
/*--------------------------------------------------------------------------*/
void node_t::new_sckt_node(const std::string& node_name, const CARD_LIST* scope)
{
  assert( is_electrical() );
  std::string::size_type dotplace = node_name.find_last_of(".");
  assert(dotplace == std::string::npos); USE(dotplace);

  assert(scope);
  new_node(node_name, scope);
  _ttt = CKT_BASE::_sim->newnode_subckt(); // equal to newnode_model, but uses different counter.
  _nnn->set_user_number(_ttt);
}
/*--------------------------------------------------------------------------*/
void node_t::new_model_adp_node(const std::string& s_in, CARD* d)
{
  string s = s_in;
  std::string::size_type dotplace = s_in.find_last_of(".");
  assert(dotplace == std::string::npos); USE(dotplace);

  // new_node(node_name, d);
  assert (d->subckt());
  new_adp_node(s, d->subckt());
  _ttt = CKT_BASE::_sim->newnode_adp(); // increase global counter.
  _nnn->set_user_number(_ttt);
  trace2("node_t::new_model_adp_node", _ttt, s_in);
}
/*--------------------------------------------------------------------------*/
node_t& node_t::map(){
  if(is_adp()){
    _m = _ttt; // + CKT_BASE::_sim->_adp_nodes;
    // incomplete(); need more generic approach. not now.
    return *this;
  }

  if (t_() != INVALID_NODE) {
    assert(_nnn);
    const NODE* n = prechecked_cast<const NODE*>(_nnn);
    assert(n); USE(n);
    _m = to_internal(t_());
    trace2("node_t::map", t_(), e_());
  }else{
    assert(_m == INVALID_NODE);
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
void node_t::map_subckt_node(uint_t* m, const CARD* owner)
{
  assert( is_electrical() );
  assert(m);
  assert(e_() !=INVALID_NODE);
  trace5("node_t::map_subckt_node", t_(), e_(), m[0], owner->long_label(), n_()->short_label());
  if (node_is_valid(m[e_()])) {
    _ttt = m[e_()];
    assert(_nnn);
  }else{untested();
    throw Exception(owner->long_label() + ": need more nodes");
  }
  //_nnn->set_flat_number(_ttt);
  assert(node_is_valid(_ttt));
}
/*--------------------------------------------------------------------------*/
double	NODE_BASE::tt_probe_num(const std::string& x)const{return tr_probe_num(x);}
XPROBE	NODE_BASE::ac_probe_ext(const std::string&)const{ return XPROBE(0);}
/*--------------------------------------------------------------------------*/

/* FIXME: throw exceptions if device doesnt exist */
NODE_BASE* NODE_BASE::lookup_node(string nodelabel, const CARD_LIST* scope)
{
  if(scope==NULL){untested();
    scope = &CARD_LIST::card_list;
  }else{
  }
  trace1("NODE_BASE::lookup_node", nodelabel);
  std::string::size_type dotplace = nodelabel.find_first_of(".");
  if (dotplace != std::string::npos) {
    string node_tail = nodelabel.substr(dotplace+1, std::string::npos);
    string container = nodelabel.substr(0, dotplace);
    trace2("NODE_BASE::lookup_node has dot ", node_tail, container);
    for (CARD_LIST::const_iterator
        i = scope->begin();  i != scope->end();  ++i) {
      CARD* card = *i;
      trace1("NODE_BASE::lookup_node... ",  card->short_label());

      if (card->is_device()
          && card->subckt()
          && wmatch(card->short_label(), container)) {
        trace1( "NODE_BASE::lookup_node dot cont: " + container + " node_tail " + node_tail ,     card->long_label());
        return lookup_node(node_tail, card->subckt());
      }else{
         // trace1( "NODE_BASE::lookup_node no device", card->short_label());
      }
    }
    trace2("NODE_BASE::lookup_node not found?! ", node_tail, container);
    throw Exception_Cant_Find( "...", container );

  }else{ // no dots, look here
    trace1("PROBELIST::add_branches no dots ", nodelabel );
    NODE_BASE* node = (*scope).node(nodelabel);
    if(!node) throw Exception_Cant_Find( "...", nodelabel );
    return node;

  }
}
/*--------------------------------------------------------------------------*/
void NODE::discont(DISCONT x)
{
  if (this == &ground_node) { untested();
    return;
  }
  _discont |= x;
  for (prop_iterator i=_prop.begin(); i!=_prop.end(); ++i) {
    assert(OPT::disc==dIMM);
    assert(*i!=this);
    (*i)->_discont |= x;
  }
}
/*--------------------------------------------------------------------------*/
void NODE::set_discont(DISCONT x)
{
  _discont = x;
}
/*--------------------------------------------------------------------------*/
void node_t::register_prop(node_t& to)
{
  (*this)->register_prop(to.operator->());
}
/*--------------------------------------------------------------------------*/
void NODE::register_prop(NODE* to)
{
  if (to==this) {
  }else if (OPT::disc != dIMM){
  }else if (to==&ground_node) { untested();
  }else if (this==&ground_node) { untested();
  }else{
    assert(to);
    for (prop_iterator i=to->_prop.begin(); i!=to->_prop.end(); ++i) {
      if (*i!=this) {
	_prop.insert(*i);
	(*i)->_coprop.insert(this);
      }
    }
    _prop.insert(to);

    for (prop_iterator i=_coprop.begin(); i!=_coprop.end(); ++i) {
      if (*i!=to) {
	(*i)->_prop.insert(to);
	to->_coprop.insert(*i);
      }
    }
    to->_coprop.insert(this);
  }

}
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
