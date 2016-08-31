/* $Id: d_rcd.cc,v 1.9 2010-09-07 07:46:21 felix Exp $ -*- C++ -*-
 * vim:ts=8:sw=2:et:
 *
 * RCD device ...
 *
 * (c) 2010 Felix Salfelder
 *
 * GPLv3
 */

#include "e_aux.h"
#include "e_adp.h"
#include "e_storag.h"
#include "globals.h"
#include "e_elemnt.h"
#include "d_rcd.h"
#include "u_nodemap.h"
#include <iomanip>
#ifdef DO_TRACE
# include "io_misc.h"
#endif

using namespace std;
/*--------------------------------------------------------------------------*/
const double _default_value (1); //input scale.
/*--------------------------------------------------------------------------*/
class ADP_NODE;
static bool dummy=false;
const double NA(NOT_INPUT);
const double INF(BIGBIG);
/*--------------------------------------------------------------------------*/
int DEV_BUILT_IN_RCD::_count = -1;
int COMMON_BUILT_IN_RCD::_count = -1;
static COMMON_BUILT_IN_RCD Default_BUILT_IN_RCD(CC_STATIC);
/*--------------------------------------------------------------------------*/
int MODEL_BUILT_IN_RCD::_count = 0;
/*--------------------------------------------------------------------------*/
namespace MODEL_BUILT_IN_RCD_DISPATCHER { //
  static DEV_BUILT_IN_RCD p1d;
  static MODEL_BUILT_IN_RCD_NET p1(&p1d);
  static DISPATCHER<MODEL_CARD>::INSTALL
    d1(&model_dispatcher, "rcdnet|rcdmodel", &p1);
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_RCD::tr_accept()
{
  //trace3("DEV_BUILT_IN_RCD::tr_accept", long_label(), _n[n_u].m_(), _n[n_u].v0());
  //trace3("DEV_BUILT_IN_RCD::tr_accept", _sim->_time0, _sim->_Time0, involts());
  // assert(subckt()); subckt()->tr_accept();
  assert(is_number((double)_tr_fill));
  if (_sim->_dt0){
    tr_stress();
  }else{
  }
  _Ccgfill->tr_lo = min(involts(), _Ccgfill->tr_lo);
  _Ccgfill->tr_hi = max(involts(), _Ccgfill->tr_hi);
  assert(is_number((double)_tr_fill));
//  q_eval();
  //trace5("DEV_BUILT_IN_RCD::tr_accept", _sim->_time0, _sim->_Time0, involts(), _Ccgfill->tr_lo, _Ccgfill->tr_hi);
}
/*--------------------------------------------------------------------------*/
void SDP_BUILT_IN_RCD::init(const COMMON_COMPONENT* cc)
{
  assert(cc);
  SDP_CARD::init(cc);
}
/*--------------------------------------------------------------------------*/
TDP_BUILT_IN_RCD::TDP_BUILT_IN_RCD(const DEV_BUILT_IN_RCD*)
{ untested();
}
/*--------------------------------------------------------------------------*/
MODEL_BUILT_IN_RCD::MODEL_BUILT_IN_RCD(const BASE_SUBCKT* p)
  :MODEL_CARD(p),
   anneal(true),
   Remodel(1e6),
   Re1(1),
   Re0(1e6),
   Rc1(1),
   Rc0(1e6),
   flags(int(USE_OPT)),
   uref(0),
   modelparm(0),
   norm_uin(false)
{
  if (ENV::run_mode != rPRE_MAIN) {
    ++_count;
  }else{
  }
  set_default(&_tnom_c, OPT::tnom_c);
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_RCD::do_tr_stress( const COMPONENT* ) const { untested();
  unreachable(); // stress done  by device.
}
/*--------------------------------------------------------------------------*/
MODEL_BUILT_IN_RCD::MODEL_BUILT_IN_RCD(const MODEL_BUILT_IN_RCD& p)
  :MODEL_CARD(p),
   anneal(p.anneal),
   Remodel(p.Remodel),
   Re1(p.Re1),
   Re0(p.Re0),
   Rc1(p.Rc1),
   Rc0(p.Rc0),
   flags(p.flags),
   uref(p.uref),
   modelparm(p.modelparm),
   norm_uin(p.norm_uin)
{
  if (ENV::run_mode != rPRE_MAIN) {
    ++_count;
  }else{untested();
  }
}
/*--------------------------------------------------------------------------*/
double MODEL_BUILT_IN_RCD::P(const COMPONENT* ) const
{ untested();
  assert(false);
  return ( NAN ); //c->_Ccgfill->get()+ _tr_fill ) * cc->_weight * cc->_wcorr;
}
/*--------------------------------------------------------------------------*/
double MODEL_BUILT_IN_RCD::__Rc(double s, const COMMON_COMPONENT* ccmp)const
{ untested();
  const COMMON_BUILT_IN_RCD* cc = prechecked_cast<const COMMON_BUILT_IN_RCD*>(ccmp);
  double ret = ( cc->_Rc0 + s * cc->_lambda * cc->_Rc1 ); 
  assert (is_number(cc->_lambda));
  assert (is_number(cc->_Rc0));
  assert (is_number(cc->_Rc1));
  assert (is_number(ret));
  return ret;
}
/*--------------------------------------------------------------------------*/
double MODEL_BUILT_IN_RCD::__Re(double , const COMMON_COMPONENT* cc)const
{ untested();
  const COMMON_BUILT_IN_RCD* c = prechecked_cast<const COMMON_BUILT_IN_RCD*>(cc);

  return c->_Re0;
//  return c->__Re(uin);
}
/*--------------------------------------------------------------------------*/
double MODEL_BUILT_IN_RCD::__Ge(double, const COMMON_COMPONENT* )const { untested();
  return NAN;
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_RCD_NET::dev_type()const
{ untested();
  return "rcdnet";
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_RCD::dev_type()const
{ untested();
  unreachable();
  if (dummy == true) { untested();
    return "rcdmodel?";
  }else{untested();//235
    return MODEL_CARD::dev_type();
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_RCD::set_dev_type(const std::string& new_type)
{ untested();
  untested();
  if (Umatch(new_type, "rcdmodel ")) { untested();
    dummy = true;
  }else{ untested();
    MODEL_CARD::set_dev_type(new_type);
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_RCD::precalc_first()
{
    const CARD_LIST* par_scope = scope();
    assert(par_scope);
    MODEL_CARD::precalc_first();
    e_val(&(this->anneal), true, par_scope);
    e_val(&(this->Remodel), 1e6, par_scope);
    e_val(&(this->Re0), 1.0, par_scope);
    e_val(&(this->Rc0), 1.0, par_scope);
    e_val(&(this->Re1), 1.0, par_scope);
    e_val(&(this->Rc1), 1.0, par_scope);
    e_val(&(this->flags), int(USE_OPT), par_scope);
    e_val(&(this->uref), NA, par_scope);
    e_val(&(this->modelparm), 0, par_scope);
    e_val(&(this->norm_uin), false, par_scope);
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_RCD::precalc_last()
{
  MODEL_CARD::precalc_last();
  const CARD_LIST* par_scope = scope();
  e_val(&(this->Re0), 1.0, par_scope);
  e_val(&(this->Rc0), 1.0, par_scope);
  e_val(&(this->Re1), 1.0, par_scope);
  e_val(&(this->Rc1), 1.0, par_scope);
}
/*--------------------------------------------------------------------------*/
SDP_CARD* MODEL_BUILT_IN_RCD::new_sdp(COMMON_COMPONENT* c)const
{
  assert(c);
  if (COMMON_BUILT_IN_RCD* cc = dynamic_cast<COMMON_BUILT_IN_RCD*>(c)) {
    if (cc->_sdp) {
      cc->_sdp->init(cc);
      assert(cc->_sdp);
      return cc->_sdp;
    }else{
      delete cc->_sdp;
      return new SDP_BUILT_IN_RCD(c);
    }
  }else{ untested();
    trace0("MODEL_BUILT_IN_RCD::new_sdp, MODEL_CARD");
    return MODEL_CARD::new_sdp(c);
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_RCD::set_param_by_index(int i, std::string& value, int offset)
{ untested();
  switch (MODEL_BUILT_IN_RCD::param_count() - 1 - i) { untested();
  case 0: untested(); break;
  case 1: _tnom_c = value; break;
  case 2: anneal = value; break;
  case 3: Remodel = value; break;
  case 4: Re0 = value; break; // former Re
  case 5: Rc0 = value; break;
  case 6: flags = value; break;
  case 7: uref = value; break;
  case 8: modelparm = value; break;
  case 9: norm_uin = value; break;
  default: throw Exception_Too_Many(unsigned(i), 7u, offset); break;
  }
}
/*--------------------------------------------------------------------------*/
bool MODEL_BUILT_IN_RCD::param_is_printable(int i)const
{
  switch (MODEL_BUILT_IN_RCD::param_count() - 1 - i) { untested();
  case 0:  return (false);
  case 1:  return (true);
  case 2:  return (true);
  case 3:  return (true);
  case 4:  return (true); // Re
  case 5:  return (true); //Rc
  case 6:  return (!(flags & USE_OPT));
  case 7:  return (uref.has_hard_value());
  case 8:  return (modelparm.has_hard_value());
  case 9:  return (norm_uin.has_hard_value());
  default: return false;
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_RCD::param_name(int i)const
{
  switch (MODEL_BUILT_IN_RCD::param_count() - 1 - i) { untested();
  case 0:  return "=====";
  case 1:  return "tnom";
  case 2:  return "anneal";
  case 3:  return "rem";
  case 4:  return "Re";
  case 5:  return "Rc";
  case 6:  return "flags";
  case 7:  return "uref";
  case 8:  return "modelparm";
  case 9:  return "norm_uin";
  default: return "";
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_RCD::param_name(int i, int j)const
{ untested();
  if (j == 0) { untested();
    return param_name(i);
  }else if (j == 1) { untested();
    switch (MODEL_BUILT_IN_RCD::param_count() - 1 - i) { untested();
    case 4:  return "Re";
    case 5:  return "Rc";
    case 9:  return "pos";
    default: return "";
    }
  }else{ untested();
    return "";
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_RCD::param_value(int i)const
{
  switch (MODEL_BUILT_IN_RCD::param_count() - 1 - i) { untested();
  case 0:  unreachable(); return "";
  case 1:  return _tnom_c.string();
  case 2:  return anneal.string();
  case 3:  return Remodel.string();
  case 4:  return Re0.string(); // former Re
  case 5:  return Rc0.string(); // former Rc
  case 6:  return flags.string();
  case 7:  return uref.string();
  case 8:  return modelparm.string();
  case 10:  return norm_uin.string();
  default: return "";
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_RCD::do_precalc_last(COMMON_COMPONENT* ccc, const CARD_LIST* )const
{ untested();
  COMMON_BUILT_IN_RCD* cc = dynamic_cast<COMMON_BUILT_IN_RCD*>(ccc);
  assert(cc);
  //const MODEL_BUILT_IN_RCD* m = this;
  //
  if (is_number(cc->Uref)) { untested();
    assert( is_number( cc->Rccommon0 ) && (double)cc->Rccommon0 != NA );
    assert( is_number( cc->Rccommon1 ) );
  }

  trace3("MODEL_BUILT_IN_RCD::do_precalc_last", cc->Uref, cc->Recommon0, cc->Rccommon0);
  long double ueff = cc->Uref; // ( exp ( lambda * Uref ) - 1 );

  double up   =  cc->Recommon0;
  double down =  cc->Rccommon0;

  double rad = double(ueff*ueff*up*up + 2.0*(up*up + up*down)*ueff + up*up - 2*up*down + down*down);
  //double s = ueff*up + up - down;
  double up_res = double ( 1.0/2.0*(ueff*up + up - down + sqrt(rad))/ueff );
  double down_res = down;

  cc->_Re0 = up_res;
  cc->_Rc0 = down_res;
  cc->_Rc1 = up_res;
  // double _rr = _rr_.subs(runter=runter, u_gate_=uref)

  // double _rh = _rh_.subs(runter=runter, u_gate_=uref)  
  double Eend_bad = (cc->Uref / (cc->Re(cc->Uref) / cc->Rc(cc->Uref) +1));

  cc->_wcorr = double ( cc->Uref / Eend_bad );
  cc->_weight = cc->weight;
  // sanity check.
  trace3("MODEL_BUILT_IN_RCD::do_precalc_last", cc->__tau_up(cc->Uref), cc->Recommon0, cc->Rccommon0);
  trace3("MODEL_BUILT_IN_RCD::do_precalc_last",cc->_Rc1, cc->_Re0, cc->_Rc0);
  assert( cc->weight != 0 );
  assert( cc->_weight != 0 );
  assert( abs( cc->__tau_up(cc->Uref) - cc->Recommon0)/cc->Recommon0 <1e-6 );
  assert( is_number( cc->_Rc1 ) );
  assert( is_number( cc->_Rc0 ) );
}
/*--------------------------------------------------------------------------*/
bool MODEL_BUILT_IN_RCD::is_valid(const COMPONENT* d)const
{
  assert(d);
  return MODEL_CARD::is_valid(d);
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_RCD::tr_eval(COMPONENT*)const
{untested();//425
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
COMMON_BUILT_IN_RCD::COMMON_BUILT_IN_RCD(int c)
  :COMMON_COMPONENT(c),
   perim(0.0),
   weight(1.0),
   Recommon0(NA),
   Recommon1(NA),
   Rccommon0(NA),
   Rccommon1(NA),
   positive(true),
   Uref(0.0),
   mu(1.0),
   lambda(1.0),
   dummy_capture(false),
   dummy_emit(false),
   _sdp(0),
   cj_adjusted(NA)
{
  trace1("COMMON_BUILT_IN_RCD::COMMON_BUILT_IN_RCD", (double) Uref);
  ++_count;
}
/*--------------------------------------------------------------------------*/
COMMON_BUILT_IN_RCD::COMMON_BUILT_IN_RCD(const COMMON_BUILT_IN_RCD& p)
  :COMMON_COMPONENT(p),
   perim(p.perim),
   weight(p.weight),
   Recommon0(p.Recommon0),
   Recommon1(p.Recommon1),
   Rccommon0(p.Rccommon0),
   Rccommon1(p.Rccommon1),
   positive(p.positive),
   Uref(p.Uref),
   mu(p.mu),
   lambda(p.lambda),
   dummy_capture(p.dummy_capture),
   dummy_emit(p.dummy_emit),
   _sdp(0),
   cj_adjusted(p.cj_adjusted)
{
  trace1("COMMON_BUILT_IN_RCD::COMMON_BUILT_IN_RCD(copy)", (double) Uref);
  ++_count;
}
/*--------------------------------------------------------------------------*/
COMMON_BUILT_IN_RCD::~COMMON_BUILT_IN_RCD()
{
  --_count;
  delete _sdp;
}
/*--------------------------------------------------------------------------*/
bool COMMON_BUILT_IN_RCD::operator==(const COMMON_COMPONENT& x)const
{
  const COMMON_BUILT_IN_RCD* p = dynamic_cast<const COMMON_BUILT_IN_RCD*>(&x);
  return (p
    && perim == p->perim
    && weight == p->weight
    && Rccommon1 == p->Rccommon1
    && Rccommon0 == p->Rccommon0
    && Recommon1 == p->Recommon1
    && Recommon0 == p->Recommon0
    && positive == p->positive
    && Uref == p->Uref
    && mu == p->mu
    && lambda == p->lambda
    && dummy_capture == p->dummy_capture
    && dummy_emit == p->dummy_emit
    && _sdp == p->_sdp
    && COMMON_COMPONENT::operator==(x));
}
/*--------------------------------------------------------------------------*/
int COMMON_BUILT_IN_RCD::param_count()const
{
  return (10 + COMMON_COMPONENT::param_count());
}
/*--------------------------------------------------------------------------*/
void COMMON_BUILT_IN_RCD::set_param_by_index(int I, std::string& Value, int Offset)
{
  trace3("COMMON_BUILT_IN_RCD::set_param_by_index ",I, COMMON_BUILT_IN_RCD::param_count() - 1 - I, Value );
  switch (COMMON_BUILT_IN_RCD::param_count() - 1 - I) { untested();
  case 0:  perim = Value; break;
  case 1:  weight = Value; trace1("wt", Value); break;
//  order          double _hl, double _hc, double _rl, double _rc)
  case 2:  Rccommon1 = Value; break; // swapped
  case 3:  Rccommon0 = Value; break;
  case 4:  Recommon1 = Value; break;
  case 5:  Recommon0 = Value; break;
  case 6:  itested(); positive = Value; break;
  case 7:  untested(); Uref = Value; break;
  case 8:  untested(); mu = Value; break;
  case 9:  untested(); lambda = Value; break;
  case 10: untested(); Uref = Value; break;
  case 11: untested(); dummy_emit = Value; break;
  default: untested(); COMMON_COMPONENT::set_param_by_index(I, Value, Offset);
  }
}
/*--------------------------------------------------------------------------*/
bool COMMON_BUILT_IN_RCD::param_is_printable(int i)const
{
  switch (COMMON_BUILT_IN_RCD::param_count() - 1 - i) { untested();
  case 0:  return (perim != 0.);
  case 1:  return (true);
  case 2:  return (true);
  case 3:  return (true);
  case 4:  return (true);
  case 5:  return (true);
  case 6:  return positive.has_hard_value();
  case 7:  return (true);
  case 8:  return mu.has_hard_value();
  case 9:  return lambda.has_hard_value();
  case 10: return Uref.has_hard_value();
  default: return COMMON_COMPONENT::param_is_printable(i);
  }
}
/*--------------------------------------------------------------------------*/
void COMMON_BUILT_IN_RCD::set_param_by_name(std::string Name, std::string Value)
{
  if (Umatch (Name,"weight ")) {
    weight = Value;
  } else if (Umatch (Name,"re1 ")) {
    Recommon1 = Value;
  } else if (Umatch (Name,"re0 ")) {
    Recommon0 = Value;
  } else if (Umatch (Name,"rc1 ")) {
    Rccommon1 = Value;
  } else if (Umatch (Name,"rc0 ")) {
    Rccommon0 = Value;
  } else if (Umatch (Name,"pos{itive}")) {
    positive = Value;
  } else { untested();
    incomplete();
    trace2("COMMON_BUILT_IN_RCD::set_param_by_name", Name, Value);
    COMMON_COMPONENT::set_param_by_name(Name, Value);
  }
}
/*--------------------------------------------------------------------------*/
std::string COMMON_BUILT_IN_RCD::param_name(int i)const
{
  trace1("COMMON_BUILT_IN_RCD::param_name", i);
  switch (COMMON_BUILT_IN_RCD::param_count() - 1 - i) { untested();
  case 0:  return "perim";
  case 1:  return "weight";
  case 2:  return "rc1";
  case 3:  return "rc0";
  case 4:  return "re1";
  case 5:  return "re0";
  case 6:  return "positive";
  case 7:  return "mu";
  case 8:  return "lam";
  case 9:  return "rcdummy";
  case 10: return "uref";
  //case 10:  return "norm_uin";
  default: return COMMON_COMPONENT::param_name(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string COMMON_BUILT_IN_RCD::param_name(int i, int j)const
{ untested();
  if (j == 0) { untested();
    return param_name(i);
  }else if (j == 1) { untested();
    switch (COMMON_BUILT_IN_RCD::param_count() - 1 - i) { untested();
    case 0:  return "";
    case 1:  return "";
    case 2:  return "";
    case 3:  return "";
    case 4:  return "";
    case 5:  return "";
    case 6:  return "";
    case 7:  return "";
    case 8:  return "";
    default: return "";
    }
  }else{untested();//281
    return COMMON_COMPONENT::param_name(i, j);
  }
}
/*--------------------------------------------------------------------------*/
std::string COMMON_BUILT_IN_RCD::param_value(int i)const
{
  switch (COMMON_BUILT_IN_RCD::param_count() - 1 - i) { untested();
  case 0:  return perim.string();
  case 1:  return weight.string();
  case 2:  return Rccommon1.string();
  case 3:  return Rccommon0.string();
  case 4:  return Recommon1.string();
  case 5:  return Recommon0.string();
  case 6:  return positive.string();
  case 7:  return Uref.string();
  case 8:  return mu.string();
  case 9:  return lambda.string();
  case 10: return Uref.string();
  case 11: return dummy_emit.string();
  default: return COMMON_COMPONENT::param_value(i);
  }
}
/*--------------------------------------------------------------------------*/
void COMMON_BUILT_IN_RCD::expand(const COMPONENT* d)
{
  COMMON_COMPONENT::expand(d);
  attach_model(d);
  COMMON_BUILT_IN_RCD* c = this; USE(c);
  const MODEL_BUILT_IN_RCD* m = dynamic_cast<const MODEL_BUILT_IN_RCD*>(model());
  if (!m) { untested();
    throw Exception_Model_Type_Mismatch(d->long_label(), modelname(), "rcd");
  }else{
  }

  if ((double)Recommon0 == NA) { Recommon0 = m->Re0;}
  if ((double)Rccommon0 == NA) { Rccommon0 = m->Rc0;}

  trace6(("COMMON_BUILT_IN_RCD::expand" + d->short_label()).c_str(), Rccommon0,
      Recommon0, m->Re0, m->Rc0, Uref, m->uref );
  // size dependent
  //delete _sdp;
  _sdp = m->new_sdp(this);
  assert(_sdp);
  const SDP_BUILT_IN_RCD* s = prechecked_cast<const SDP_BUILT_IN_RCD*>(_sdp);
  assert(s); USE(s);

  // subcircuit commons, recursive
  assert(c == this);
}
/*--------------------------------------------------------------------------*/
/* calculate before model? */
/*--------------------------------------------------------------------------*/
void COMMON_BUILT_IN_RCD::precalc_first(const CARD_LIST* par_scope)
{
  trace2("precalc_first", *(par_scope->params()), weight);
  assert(par_scope);
  COMMON_COMPONENT::precalc_first(par_scope);
  e_val(&(this->perim), 0.0, par_scope);
  e_val(&(this->weight), 1.0, par_scope); // fixme: oscale?
  e_val(&(this->Recommon0), 1.0, par_scope);
  e_val(&(this->Recommon1), 1.0, par_scope);
  e_val(&(this->Rccommon0), 1.0, par_scope);
  e_val(&(this->Rccommon1), 1.0, par_scope);
  e_val(&(this->Uref), 0.00001, par_scope);
  trace3("uref...",  Uref, NOT_INPUT, Uref );
  e_val(&(this->mu), 1.0, par_scope);
  e_val(&(this->lambda), 1.0, par_scope);
  e_val(&(this->dummy_capture), false, par_scope);
  e_val(&(this->dummy_emit), false, par_scope);
  e_val(&(this->positive), true, par_scope);
  if(!positive){
  }else{
  }
}
/*--------------------------------------------------------------------------*/
void COMMON_BUILT_IN_RCD::precalc_last(const CARD_LIST* par_scope)
{
  assert(par_scope);
  trace2("precalc_last", *(par_scope->params()), weight);
  COMMON_COMPONENT::precalc_last(par_scope);
  COMMON_BUILT_IN_RCD* cc = this;
  COMMON_BUILT_IN_RCD* c = this;
  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(model());
  assert(m);
  // final adjust: code_pre
  trace2("COMMON_BUILT_IN_RCD::precalc_last", m->v2(), m->uref);
  // final adjust: override
  // final adjust: raw
  e_val(&(this->weight), 1.0, par_scope); // fixme: oscale?
  e_val(&(this->perim), 0.0, par_scope);
  e_val(&(this->weight), 1.0, par_scope);
  e_val(&(this->Recommon0), m->Re0, par_scope);
  e_val(&(this->Recommon1), m->Re1, par_scope);
  e_val(&(this->Rccommon0), m->Rc0, par_scope);
  e_val(&(this->Rccommon1), m->Rc1, par_scope);
  e_val(&(this->Uref), m->uref, par_scope);
  
  trace4("uref...",  m->uref, NOT_INPUT, Uref , double (Uref));
  Uref.e_val( (double)m->uref, par_scope);
  trace3("uref...",  m->uref, NOT_INPUT, Uref );

  e_val(&(this->mu), 1.0, par_scope);
  e_val(&(this->lambda), 1.0, par_scope);
  e_val(&(this->dummy_capture), false, par_scope);
  e_val(&(this->positive), true, par_scope);
  if(!positive){
  }else{
    trace2("", positive, (bool)positive);
  }
  e_val(&(this->dummy_emit), false, par_scope);
  // final adjust: mid
  // final adjust: calculated
  cj_adjusted = 19.0;

  _lambda = 1;
  lambda=1;

  if((double)Uref == NOT_INPUT) { 
    trace2("nin",  m->uref, NOT_INPUT );
    Uref = 0.0;
  }
  if((double)Uref == NA) { 
    trace0("na");
    Uref = 0.0;
  }

  // size dependent
  //delete _sdp;
  _sdp = m->new_sdp(this);
  assert(_sdp);
  const SDP_BUILT_IN_RCD* s = prechecked_cast<const SDP_BUILT_IN_RCD*>(_sdp);
  assert(s); USE(s);

  // subcircuit commons, recursive
  c->_wcorr = 0;
  c->_zero = 0;


  if(m->v2()){
    m->do_precalc_last(cc, par_scope );
  } else if (Uref!=0.0 ){ untested();
    m->do_precalc_last(cc, par_scope );
  } else { // no uref...
    cc->_Re0 = Recommon0;
    cc->_Rc0 = Rccommon0;
    cc->_Rc1 = _Re0;

    if(Rccommon1 != NA && Rccommon1 != NOT_INPUT ) _Rc1 = Rccommon1;

    _weight = weight;
    _wcorr = 1;
    assert (weight != 0);
    trace5("COMMON_BUILT_IN_RCD::precalc_last no uref. simple model", _Re0, _Rc0, _Rc1, Re(1), Rc(0));
    Uref=1;
  }
  trace3("COMMON_BUILT_IN_RCD::precalc_last", Uref*weight, _Re0, _Rc0 );
  trace1("COMMON_BUILT_IN_RCD::precalc_last", cc->_Re0);
  trace1("COMMON_BUILT_IN_RCD::precalc_last done", m->v2());

  if (cc->_Re1 < 0) { untested(); // turnt
    throw Exception_Precalc("nonnegative Re1: " + ::to_string(cc->_Re1) + "\n");
  }
  if (cc->_Rc1 > 0) { untested(); // turnt
    throw Exception_Precalc("nonpositive Rc1: " + ::to_string(cc->_Rc1) + "\n");
  }
}
/*--------------------------------------------------------------------------*/
//double MODEL_BUILT_IN_RCD::__tau_up ( double uin, const COMMON_BUILT_IN_RCD* cc ) const{ untested();
//        return cc->__tau_up(uin);
//}
/*--------------------------------------------------------------------------*/
double COMMON_BUILT_IN_RCD::__tau_up ( double uin ) const{ untested();
  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(model());
  double  rc = m->__Rc(uin, this);
  double  re = m->__Re(uin, this);
  return float( rc / ( 1 + rc/re )  ) ;
} 
/*--------------------------------------------------------------------------*/
namespace DEV_BUILT_IN_RCD_DISPATCHER { 
  static DEV_BUILT_IN_RCD p0;
  static DISPATCHER<CARD>::INSTALL
    d0(&device_dispatcher, "Z|rcd", &p0);
}
/*--------------------------------------------------------------------------*/
double COMMON_BUILT_IN_RCD::Re(double uin) const
{ untested();
  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(model());
  return m->__Re(uin,this);
}
/*--------------------------------------------------------------------------*/
double COMMON_BUILT_IN_RCD::Rc(double uin) const
{ untested();
  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(model());
  return m->__Rc(uin,this);

   //double ret = ( _Rc0 + uin * _lambda * _Rc1 ); 
   //assert (ret==ret);
   //return ret;
}
/*--------------------------------------------------------------------------*/
static EVAL_BUILT_IN_RCD_GRc Eval_GRc(CC_STATIC);
void EVAL_BUILT_IN_RCD_GRc::tr_eval(ELEMENT* d)const
{ untested();
  assert(d);
  DEV_BUILT_IN_RCD* p = prechecked_cast<DEV_BUILT_IN_RCD*>(d->owner());
  assert(p);
  const COMMON_BUILT_IN_RCD* cc = prechecked_cast<const COMMON_BUILT_IN_RCD*>(p->common());
  const SDP_BUILT_IN_RCD* s = prechecked_cast<const SDP_BUILT_IN_RCD*>(cc->sdp());
  assert(s); USE(s);
  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(cc->model());
  assert(m); USE(m);

  // FIXME: merge with __Rc somehow
  double _c[3] = { cc->_Rc0, cc->_Rc1 * cc->_lambda, 0 };
  double x = (d->_y[0].x);
  //    trace1("Rc", x);
  double f0 = 0.;
  double f1 = 0.;
  for (size_t i=1; i>0; --i) { untested();
    f0 += _c[i];
    f0 *= x;
    f1 *= x;
    f1 += _c[i]*int(i);
  }
  f0 += _c[0];
  d->_y[0] = FPOLY1(x, f0, f1);

  //  d->set_converged(d->conv_check());
}
/*--------------------------------------------------------------------------*/
DEV_BUILT_IN_RCD::DEV_BUILT_IN_RCD()
  :BASE_SUBCKT(),
   // input parameters,
   // calculated parameters,
   _lasts(-inf),
   _region(UNKNOWN),
   // netlist,
   _Ccg(0),
   _Ye(0),
   _Re(0),
   _Rc(0),
   _GRc(0),
   _Ccgfill(0),
   _tr_fill(0),
   _tr_dfill(0),
   _iter_newton(0),
   _iter_bisect(0)
{

  _n = _nodes;
  _nodes[n_p].set_adp();
  assert(_n[n_p].is_adp());
  assert(_n[n_u].is_electrical());
  attach_common(&Default_BUILT_IN_RCD);
  // _nodes[n_p].a_()->trhack = 0;

  ++_count;
  // overrides
  set_value(1);
}
/*--------------------------------------------------------------------------*/
DEV_BUILT_IN_RCD::DEV_BUILT_IN_RCD(const DEV_BUILT_IN_RCD& p)
  :BASE_SUBCKT(p),
   // input parameters,
   // calculated parameters,
   _lasts(-inf),
   _region(p._region),
   // netlist,
   _Ccg(0),
   _Ye(0),
   _Re(0),
   _Rc(0),
   _GRc(0),
   _Ccgfill(0),
   _tr_fill(p._tr_fill),
   _tr_dfill(p._tr_dfill),
   _iter_newton(0),
   _iter_bisect(0)
{
  _n = _nodes;
  assert(p._n[n_p].is_adp());
  for (uint_t ii = 0; ii < max_nodes() + int_nodes(); ++ii) {
    _n[ii] = p._n[ii];
    trace2("DEV_BUILT_IN_RCD::DEV_BUILT_IN_RCD(p)", ii, _n[ii].is_adp());
  }
  assert(_n[n_p].is_adp());
  ++_count;
  // overrides
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_RCD::expand()
{
  trace1("DEV_BUILT_IN_RCD::expand", long_label());
  assert(_n[n_p].is_adp());
  BASE_SUBCKT::expand(); // calls common->expand, attaches model
  assert(_n);
  assert(common());
  const COMMON_BUILT_IN_RCD* c = static_cast<const COMMON_BUILT_IN_RCD*>(common());
  assert(c);
  assert(c->model());
  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(c->model());
  assert(m); USE(m);
  assert(c->sdp());
  const SDP_BUILT_IN_RCD* s = prechecked_cast<const SDP_BUILT_IN_RCD*>(c->sdp());
  assert(s); USE(s);

  if (!subckt()) {
    new_subckt();
  }else{ untested();
  }
  if (_sim->is_first_expand()) {
    precalc_first();
    precalc_last();
    trace4("DEV_BUILT_IN_RCD::expand, first", long_label(), hp(this), _n[n_p].n_(), _n[n_p].is_adp());
    assert(_n[n_p].is_adp());
    assert(!(_n[n_p].n_())); // n_ is electrical
    if (!(_n[n_p].a_())){
      trace1("DEV_BUILT_IN_RCD::expand, no 3rd external node connected", hp(_n));
      _n[n_p].new_model_adp_node("c", this);
      trace1("DEV_BUILT_IN_RCD::expand have new adpnode", hp(_n[n_p].a_()));
    } else { itested();
      trace1("DEV_BUILT_IN_RCD::expand, 3rd node present, external", long_label());
    }
  }

  _Ccgfill = _n[n_p].a_(); assert(_Ccgfill);

// idee: _Ccgfill:: tr_value <= udc
//                  tt_value <= E
//
 // _Udc = new ADP_NODE_UDC((const COMPONENT*) common()); //, _Ccgfill);
 // ADP_NODE_LIST::adp_node_list.push_back( _Udc );

  //precalc();
  subckt()->expand();
  //subckt()->precalc();
  assert(!is_constant());
  if ( adp() == NULL ){
  //  attach_adp( m->new_adp( (COMPONENT*) this ) );
  }else{ untested();
    untested(); // rebuild circuit??
  }

  if (m->v2()){
    // incomplete();
    // do_tt?
    // _Ccgfill->tt_set( -c->_wcorr );
  }
}
/*--------------------------------------------------------------------------*/
double DEV_BUILT_IN_RCD::P()const {
  const COMMON_BUILT_IN_RCD* c = static_cast<const COMMON_BUILT_IN_RCD*>(common());
  assert(c);
  assert(c->model());
  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(c->model());
  return m->P(this);
}
/*--------------------------------------------------------------------------*/
double DEV_BUILT_IN_RCD::tr_probe_num(const std::string& x)const
{
  assert(_n);
  const COMMON_BUILT_IN_RCD* c = prechecked_cast<const COMMON_BUILT_IN_RCD*>(common());
  assert(c);
  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(c->model());
  assert(m);
  const SDP_BUILT_IN_RCD* s = prechecked_cast<const SDP_BUILT_IN_RCD*>(c->sdp());
  assert(s); USE(s);
  ADP_NODE* Ccgfill = _n[n_p].a_();
  assert(Ccgfill);

  if (Umatch(x, "region ")) { untested();
    return  static_cast<double>(region());
  }else if (Umatch(x, "ulo |vlo ")) {
    return  ( Ccgfill->tr_lo );
  }else if (Umatch(x, "uin |vin ")) {
//    assert (involts() <= _Ccgfill->tr_hi || _sim->_time0==0);
//    assert (involts() >= _Ccgfill->tr_lo || _sim->_time0==0);
    return involts();
  }else if (Umatch(x, "s{tresslevel} ")) {
//    assert(involts() * value() == _Ccgfill->tr());
//    return Ccgfill->tr(); (hmmm.)
    return involts() * value();
  }else if (Umatch(x, "uhi |vhi ")) {
    return  ( Ccgfill->tr_hi );
  }else if (Umatch(x, "tra ")) { untested();
    return  ( Ccgfill->tr_abs_err() );
  }else if (Umatch(x, "noise ")) { untested();
    assert(Ccgfill);
    return 0; // ( Ccgfill->get_tr_noise() );
  }else if (Umatch(x, "udc |ueff")) {
    return (_Ccgfill->tr());
  }else if (Umatch(x, "dudc ")) {
    return(_Ccgfill->tr1() - _Ccgfill->tr());
  }else if (Umatch(x, "udc1 ")) {
    return(_Ccgfill->tr1());
  }else if (Umatch(x, "trr ")) { untested();
    return  ( _Ccgfill->tr_rel_err() );
  } else if (Umatch(x, "RE0 "    )) { return( c->_Re1 );}
    else if (Umatch(x, "RE1 "    )) { return( c->_Re0 );}
    else if (Umatch(x, "RC0 "    )) { return( c->_Rc1 );}
    else if (Umatch(x, "RC1 "    )) { return( c->_Rc0 );}
    else if (Umatch(x, "E0|zero "    )) {
    return  ( c->_zero );
  }else if (Umatch(x, "dE ")) { untested();
    return (double)_tr_fill; // shifted
  }else if (Umatch(x, "E ")) {
    return (double)_tr_fill; // shifted == unshifted
  }else if (Umatch(x, "E{end} |E_end ")) { itested();
    if(is_number(_Ccgfill->tr()))
      return double(m->E_end((long double)(value()*_Ccgfill->tr()),c));
    return double(m->E_end((long double)(value()*involts()),c));
  }else if (Umatch(x, "te ")) { untested();
      return  ( c->__tau_up(c->Uref) );
  }else if (Umatch(x, "tc ")) { untested();
      return  ( m->__Rc(0., c) );
  }else if (Umatch(x, "re ")) { untested();
      return  ( m->__Re(c->Uref, c) );
  }else if (Umatch(x, "rc ")) { untested();
      return  ( m->__Rc(0. , c) );
  }else if (Umatch(x, "tau ")) {
    return  ( double(m->__tau((long double)(value()*involts()),c)));
  }else if (Umatch(x, "tra ")) { untested();
    return  ( _Ccgfill->tr_abs_err() );
  }else if (Umatch(x, "uref ")) { untested();
    return  ( c->Uref );
  }else if (Umatch(x, "P ")) {
    return  ( P() );
  }else if (Umatch(x, "vw{v} ")) { untested();
    assert (c->_weight != .0);
    return _Ccgfill->get_total() * c->_weight * c->_wcorr;
  }else if (Umatch(x, "tt ")) { itested();
    return (double)_Ccgfill->tt();
  }else if (Umatch(x, "order " )) { return( _Ccgfill->order() );
  }else if (Umatch(x, "fill ")) { untested();
      return (double)_tr_fill;
#ifndef NDEBUG
  }else if (Umatch(x, "chp ")) {
    return hp(_Ccgfill);
  }else if (Umatch(x, "m{atrix} ")) {
    return _Ccgfill->m_();
  }else if (Umatch(x, "v{alue} ")) {
    return value();
#endif
  }else if (Umatch(x, "wt ")) { untested();
    return  c->_weight;
  }else if (Umatch(x, "status ")) { untested();
    return  static_cast<double>(converged() * 2);
  }else if (Umatch(x, "_region ")) { untested();
    return _region;
  }

  return BASE_SUBCKT::tr_probe_num(x);
}
/*--------------------------------------------------------------------------*/
double DEV_BUILT_IN_RCD::tt_probe_num(const std::string& x)const
{
  assert(_n);
  const COMMON_BUILT_IN_RCD* c = prechecked_cast<const COMMON_BUILT_IN_RCD*>(common());
  const COMMON_BUILT_IN_RCD* cc = prechecked_cast<const COMMON_BUILT_IN_RCD*>(common());
  assert(c);
  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(cc->model());
  assert(m);
  const SDP_BUILT_IN_RCD* s = prechecked_cast<const SDP_BUILT_IN_RCD*>(cc->sdp());
  assert(s); USE(s);
  ADP_NODE* Ccgfill = _n[n_p].a_();
  assert(Ccgfill);

  if (Umatch(x, "P "))  {
    return P();
  } else if (Umatch(x, "vc "))  { untested();
    //  return _Ccgfill->tt();
    //  depends on model!
    assert(is_number(  ( P() )));
    return  ( P() );
  }
  else if (Umatch(x, "tr ")) {
    return( _Ccgfill->tr_get() );
  } else if (Umatch(x, "tt ")) { itested();
    return( _Ccgfill->tt() );
  } else if (Umatch(x, "dtt ")) { untested();
    return( _Ccgfill->tt()-_Ccgfill->tt1());
  } else if (Umatch(x, "ttg{ain} ")) { itested();
    return _ttgain;
  } else if (Umatch(x, "ttf{uture} ")) { itested();
    return _ttfuture;
  } else if (Umatch(x, "ttr ")) { untested();
    return( _Ccgfill->tt_rel_err() );
  } else if (Umatch(x, "trr "   )) { untested();
    return( _Ccgfill->tr_rel_err() ); }
  else if (Umatch(x, "v{alue} "  )) { untested();
    return value();
  } else if (Umatch(x, "s{tresslevel} ")) {
    return Ccgfill->tr();
  }
#ifndef NDEBUG
  else if (Umatch(x, "iter "  )) { return( _iter_newton ); }
  else if (Umatch(x, "in "  )) { return( _iter_newton ); }
  else if (Umatch(x, "ib "  )) { return( _iter_bisect ); }
#endif
  else if (Umatch(x, "tra "   )) { return( _Ccgfill->tr_abs_err() ); }
  else if (Umatch(x, "Rc "    )) { return( c->_Rc0 ); }
  else if (Umatch(x, "wt "    )) { return( c->_weight ); }
  else if (Umatch(x, "region ")) { return( m->tt_region( this ) ); }
  else if (Umatch(x, "uref "  )) { return( c->Uref ); }
  else if (Umatch(x, "udc "))    { return( _Ccgfill->tr() ); }
  else if (Umatch(x, "tauinv "))    { untested();
    long double n = _Ccgfill->tr();
    if(!is_number(n)) return NOT_VALID;
    return( double(m->__tau_inv( value()*n,cc))); }
  else if (Umatch(x, "tau "))    {
    long double n = _Ccgfill->tr();
    if(!is_number(n)) return NOT_VALID;
    return( double(m->__tau( value()*n,cc))); }
  else if (Umatch(x, "uend "  )) { return( c->Uref / (c->Re(c->Uref) / c->Rc(c->Uref) +1) * c->_wcorr ) * c->_weight; }
  else if (Umatch(x, "tc "    )) { return( c->Rc(0) ); }
  else if (Umatch(x, "te ")) { untested();
    if (m->v2())
      return( m->__Re(cc->Uref,cc));
    else
    return  ( c->__tau_up( c->Uref ) );
  }
  else if (Umatch(x, "vwtr ")) { untested();
    return  ( _Ccgfill->tr_get() * c->_weight );
  }
  return BASE_SUBCKT::tt_probe_num(x);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// cc_direct

double ADP_BUILT_IN_RCD::tt_probe_num(const std::string& )const
 {untested(); return 888;}
double ADP_BUILT_IN_RCD::tr_probe_num(const std::string& )const
{untested(); return 888;}
void ADP_BUILT_IN_RCD::init(const COMPONENT* )
{ untested();
  untested();
}
/*--------------------------------------------------------------------------*/
region_t MODEL_BUILT_IN_RCD::region( const COMPONENT* )const
{ untested();
  return UNKNOWN;
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
int MODEL_BUILT_IN_RCD::tt_region( const COMPONENT* )const
{ untested();
  return 179;
}
/*--------------------------------------------------------------------------*/
region_t DEV_BUILT_IN_RCD::region(  )const
{ untested();
  return UNKNOWN;
}
/*--------------------------------------------------------------------------*/
int  DEV_BUILT_IN_RCD::tt_region(  )const
{ untested();
  return 14;
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_RCD::tt_eval(COMPONENT* )const
{ untested();
  untested();
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_RCD::do_tt_prepare(COMPONENT*)const
{ untested();
}
///*--------------------------------------------------------------------------*/
ADP_CARD* MODEL_BUILT_IN_RCD::new_adp( COMPONENT* c)const
{ untested();
  trace0("MODEL_BUILT_IN_RCD::new_adp");
  assert(c);
  return MODEL_CARD::new_adp(c);
}
/*--------------------------------------------------------------------------*/
bool DEV_BUILT_IN_RCD::tr_needs_eval()const
{
  return 0; // not needed to call q_accept.
}
/*--------------------------------------------------------------------------*/
long double MODEL_BUILT_IN_RCD::__step(long double s, long double cur, double deltat, const COMMON_COMPONENT* c ) const
{
//  trace3("MODEL_BUILT_IN_RCD::__step", s, cur, deltat);
  assert(deltat>=0.);
  assert(is_number(s));
  assert(is_number(cur));
  assert(cur>-.1 && cur <1.1);
  const COMMON_BUILT_IN_RCD* cc = dynamic_cast<const COMMON_BUILT_IN_RCD*>(c);
  const MODEL_BUILT_IN_RCD* m = static_cast<const MODEL_BUILT_IN_RCD*>(this);
  USE(m);
  if(cc->positive && s<0) {
    s = 0;
  }else if (s<0) {
  }else{
  }
  trace3("", s, E_end(s,cc), E_end(-s,cc));
  assert ( s >= -0.001 || !cc->positive);
  long double Eend = E_end(s,cc);
  assert(Eend<=1.);
  long double tauinv = __tau_inv(s,cc);
  long double ret = (cur-Eend) * expl( -deltat*tauinv ) + Eend;

  assert(ret<=1.0000001);

//  trace5("MODEL_BUILT_IN_RCD::__step", Eend, deltat, s,  ret,  logl(fabsl(cur-Eend ) ) );

  assert(ret>=0.);
  if (ret>1) {
    return 1;
  }else{
    return ret;
  }
}
///*--------------------------------------------------------------------------*/
long double DEV_BUILT_IN_RCD::extrapolate(const double* tr0)const
{
  if(tr0) assert(is_number(*tr0));
  const COMMON_BUILT_IN_RCD* c = static_cast<const COMMON_BUILT_IN_RCD*>(common());
  const COMMON_BUILT_IN_RCD* cc = c;
  assert(c);
  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(c->model());
  assert(m);
  long double E_old = _Ccgfill->tt1() + c->_zero; // E at Time1 or last_Time
  long double fill_new = E_old;
  double toT = _sim->_Time0;
  double fromT = _sim->_Time0 - _sim->_dT0 + _sim->last_time();
  if (!tr0) {
    //pred mode
  }else{
    trace1("", *tr0);
    toT += _sim->last_time();
  }
  double ex_time = toT - fromT;

  trace8("DEV_BUILT_IN_RCD::extrapolate", ex_time, _sim->_dT0, _sim->last_time(), _sim->last_Time(), fromT, toT, E_old, tr0);
  if( ex_time < 1e-18 ){
    assert(is_almost(double(fill_new) , double (E_old)));
  }
  assert(ex_time >= 0);

  long double eff1 = value()*_Ccgfill->tr( fromT + ex_time/3.0, tr0);
  long double eff2 = value()*_Ccgfill->tr( fromT + ex_time*2.0/3.0, tr0);

  if(!is_number(eff1)){
    error(bDANGER, "nan bug %d %f\n", _Ccgfill->order(), tr0);
    assert(0);
  }
  assert(is_number(eff2));

  if(cc->positive) {
    eff1 = max(eff1,0.0L);
    eff2 = max(eff2,0.0L);
  }else{ itested();
  }

  fill_new = m->__step( eff1, fill_new, ex_time*.5, c );
  assert(is_number(fill_new));
  fill_new = m->__step( eff2, fill_new, ex_time*.5, c );
  assert(is_number(fill_new));
  assert(1.01 > fill_new && fill_new > -.01);

  if (fill_new>1) {
    return 1.;
  }else if (fill_new<0) {
    return 0.;
  }else {
    return fill_new;
  }
}
///*--------------------------------------------------------------------------*/
// dashier funktioniert (aus e_subckt.h) 
//   void tr_queue_eval() {assert(subckt()); subckt()->tr_queue_eval();}
// dashier macht der modelgen:
//   void tr_queue_eval()      {if(tr_needs_eval()){q_eval();}}
///*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_RCD::do_tt()
{
  assert(_n);
  assert(common());
  const COMMON_BUILT_IN_RCD* c = static_cast<const COMMON_BUILT_IN_RCD*>(common());
  assert(c);
  assert(c->model());
  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(c->model());
  assert(m); USE(m);
  assert(c->sdp());
  const SDP_BUILT_IN_RCD* s = prechecked_cast<const SDP_BUILT_IN_RCD*>(c->sdp());
  assert(s); USE(s);

  assert(_sim->_time0 == 0 || _sim->_mode==s_TRAN ); //?

  trace5("DEV_BUILT_IN_RCD::do_tt", long_label(), _sim->_dT0, _Ccgfill->tt(), _sim->_Time0, _Ccgfill->tr() );
  if(_sim->phase() == p_PD){
    // assert(!_Ccgfill->tr());
    _Ccgfill->tr() = 0; // check: what is tr()?
                        // why not use tr(0) and leave tr alone?
    assert(_sim->_dT0);
  }

  if (_sim->_dT0==0.) {
    assert(c->_zero); // 0 is very unlikely...
    _tr_fill = _Ccgfill->tt() + c->_zero; // shifted
    assert(is_number(_tr_fill));
    assert(_tr_fill<=1+1e-9);
    if (_tr_fill>1.) { untested();
      _tr_fill = 1.;
    }
    _tr_dfill = 0;
    _nodes[n_p].a_()->trhack = (double) _tr_dfill;
    return;
  } else {
    assert(is_number(_tr_fill));
  }

  double Time1 = _sim->_Time0 - _sim->_dT0;

  trace5("DEV_BUILT_IN_RCD::do_tt", _sim->_Time0, _sim->_dT0, tt_iteration_number(), _sim->_Time0, _sim->_time0);
  trace4("DEV_BUILT_IN_RCD::do_tt ", _Ccgfill->tr() , _Ccgfill->tr(_sim->_Time0 ), _Ccgfill->order(), _sim->_time0);

  if (! ( is_almost( _Ccgfill->tr1() , _Ccgfill->tr(Time1) ))) {
      error(bDANGER, "DEV_BUILT_IN_RCD::tr_stress !almost tr1 %E tr(T1) %E T1 %f. order %d\n",
          double(_Ccgfill->tr1()) , double(_Ccgfill->tr(Time1)), Time1, _Ccgfill->order() );
  }

  if (!is_number(_Ccgfill->tr_rel(_sim->_dT0 ))) { unreachable();
    error(bDANGER, "DEV_BUILT_IN_RCD::do_tt tr_rel nan at %f dT %f, order %d\n",
        _sim->_Time0, _sim->_dT0, _Ccgfill->order());
//    assert(0);
  }

  long double E_old = _Ccgfill->tt1() + c->_zero; // shifted
  assert(E_old<=1+1e-9);
  if (E_old>1.) { untested();
    E_old = 1.;
  }
  _tr_fill = E_old;
  assert(_tr_fill<=1);

  assert (is_number(E_old));

  long double fill_new = extrapolate();
  assert(fill_new <= 1);

  trace3("DEV_BUILT_IN_RCD::do_tt", fill_new, E_old, fill_new-_tr_fill );

  assert(is_number(fill_new));

  trace3("DEV_BUILT_IN_RCD::do_tt", _tr_fill, _sim->tt_iteration_number(), c->_zero );

  assert(is_number(fill_new));
  _Ccgfill->tt() = (double) fill_new - c->_zero;
  assert(_Ccgfill->tt() > -.01);
  assert(_Ccgfill->tt() <= 1.);
  _tr_fill = fill_new;
  assert(is_number(_tr_fill));
  assert(_tr_fill<=1);
  _tr_dfill = 0;
  _nodes[n_p].a_()->trhack = (double) _tr_dfill;
  trace2("DEV_BUILT_IN_RCD::do_tt done", fill_new, _tr_fill );

  // print something useful for "stress" (BUG?)
  _Ccgfill->tr() = _Ccgfill->tr(_sim->_Time0);
}
///*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_RCD_NET::do_stress_apply( COMPONENT*  ) const
{ untested();
}
///*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_RCD::tr_advance()
{
  BASE_SUBCKT::tr_advance(); //?
  const COMMON_BUILT_IN_RCD* c = prechecked_cast<const COMMON_BUILT_IN_RCD*>(common());
  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(c->model());
  assert(c);
  ADP_NODE* a = _n[n_p].a_();
  set_not_converged();
  assert(a);
  double stress = involts() * value();
  _time1 = _time0;
  double _dt = _sim->_time0 - _time1;
  if(_dt < 0){ unreachable();
    // BUG in tr_begin?. tt_advance?. should not be here
    _dt = 0.;
  }
  _time0 = _sim->_time0;
  long double newtrhack;
  assert(_dt>=0.);
  if(_dt){
    trace3("??", _dt, stress, long_label());
    newtrhack = m->__step( stress, a->tt() + _tr_dfill + c->_zero, _dt, c ) - c->_zero;
  }else{
    // BUG?!
    newtrhack = a->tt() + _tr_dfill;
  }
  trace7("tradv",_sim->_time0, a->trhack, newtrhack, _dt, stress, a->tt(), _tr_dfill);

  // predictor...
  if (_dt) {
    a->trhack = (double) newtrhack - a->tt();
  }
}
///*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_RCD::tr_regress()
{
  BASE_SUBCKT::tr_regress(); //??
  const COMMON_BUILT_IN_RCD* c = prechecked_cast<const COMMON_BUILT_IN_RCD*>(common());
  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(c->model());
  assert(c);
  ADP_NODE* a = _n[n_p].a_();
  set_not_converged();
  double stress = involts() * value();
  _time0 = _sim->_time0;
  double _dt = _sim->_time0 - _time1;
  trace8("",_sim->_time0, a->trhack, a->State(), _dt, a->tr(), stress, a->tt(), _tr_dfill);

  if(_dt < 0){ unreachable();
    // BUG in tr_begin?. tt_advance?. should not be here
    _dt = 0.;
  }
  // predictor...
  a->trhack = double(m->__step( stress, a->tt() + _tr_dfill + c->_zero, _dt, c ) - c->_zero);
}
///*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_RCD::tr_stress()
{
#ifndef NDEBUG
  if (_stressiter == (unsigned)(_sim->iteration_tag())){ untested();
    error(bDANGER, "double stress at %E in %s, %i %i\n", _sim->_time0,
                   long_label().c_str(), _stressiter, _sim->iteration_tag());
    assert(0);
  }
  _stressiter = _sim->iteration_tag();
#endif

  // fixme: assert (TTT || trage)
  const DEV_BUILT_IN_RCD* rcd = this;
  double h = _sim->_dt0;
  assert(h || !_sim->_time0);
  const COMMON_BUILT_IN_RCD* c = static_cast<const COMMON_BUILT_IN_RCD*>(common());
  const COMMON_BUILT_IN_RCD* cc = c;
  assert(c);
  assert(c->model());
  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(c->model());
  assert(m);
  assert(c->sdp());
  assert(_Ccgfill);
  ADP_NODE* _c = _Ccgfill;

  long double stress = rcd->involts() * value();

  if( stress < -0.001 && cc->positive ) {
    error(bLOG, "at %f, %f in %s: stress too small: %f, u %f b %f\n",
        _sim->_Time0, _sim->_time0,
        long_label().c_str(), (double)stress, _n[n_u].v0(), _n[n_b].v0());
    stress = 0.;
  } else if ( stress < 0 && cc->positive ) {
    stress = 0.;
  }else{
  }

  assert(value());
  assert (stress==stress);
//  trace7("DEV_BUILT_IN_RCD::tr_stress", long_label(), stress, _sim->_time0, _lasts, _tr_fill, h, value());

  if( _sim->_time0 > _lasts ){
    _lasts = _sim->_time0;
  }else {
    trace1("DEV_BUILT_IN_RCD::tr_stress again?? bug??", _sim->_time0 );
    if (! (_sim->_time0 == _lasts) ){
      error(bDANGER,"DEV_BUILT_IN_RCD::tr_stress critically unequal now: %E lasts: %E at %E\n",
          _sim->_time0, _lasts, _sim->_Time0 );
      unreachable();
      assert(0);
      throw(Exception("time mismatch in %s: time0: %E lasts: %E " , long_label().c_str(), _sim->_time0, _lasts));
    }
    return;
  }

  assert(is_number(involts()));
  assert(fabs(involts()<1e3));

  if( _sim->_time0==0 ){
    assert( _Ccgfill->tr_lo == inf );
    assert( _Ccgfill->tr_hi == -inf );
    // assert !tran_dynamic()?

    _Ccgfill->tr_lo = min(involts(), _Ccgfill->tr_lo);
    _Ccgfill->tr_hi = max(involts(), _Ccgfill->tr_hi);

    // _Ccgfill->set_tr(_Ccgfill->tt());

  } else {
    // do not update tr_lo, tr_hi. could be wrong!

    // assert tran_dynamic()?

    _Ccgfill->tr_lo = min(involts(), _Ccgfill->tr_lo);
    _Ccgfill->tr_hi = max(involts(), _Ccgfill->tr_hi);

    assert( is_number(_Ccgfill->tr_lo) );
    assert( is_number(_Ccgfill->tr_hi) );
    assert(fabs(_Ccgfill->tr_lo) < 1e5);
    assert(fabs(_Ccgfill->tr_hi) < 1e5);

    if( -15<  _Ccgfill->tr_hi &&  _Ccgfill->tr_hi < 100  &&
        -15<  _Ccgfill->tr_lo &&  _Ccgfill->tr_lo < 100 ){
    }else{unreachable();
#ifndef NDEBUG
      error(bDANGER, "%s something wrong with input range %f %f at time %E\n",
          long_label().c_str(), _Ccgfill->tr_lo, _Ccgfill->tr_hi, _sim->_time0 );
      assert(false);
#endif
    }
  }
  if (!h) { untested();
    trace1("not h\n", _tr_fill);
    return;
  }
  if(_sim->phase()==p_PD){
        assert(!_c->tr_lo); USE(_c);
        assert(!_c->tr_hi);
  }

  if (stress >= 0) {
    if (m->E_end(stress,c) <= -1e-16) { // shifted. untested();
        error(bDANGER, "nonmonotony at %f: %f < %f\n", stress, c->_zero, m->E_end(stress,c)  );
    }
    assert(m->E_end(stress,c) > -1e-15); // shifted
  }else if (cc->positive && stress < -1e-14) { itested();
    error(bTRACE, "in %s not positive uin_s=%LE, value %f\n", long_label().c_str(), stress, double(value()));
  }else if(cc->positive) { untested();
    stress = 0;
  }else{
  }
  assert(double(value()));

  assert( is_almost(double(m->E_end_0(c)), c->_zero ));
  assert( is_almost(double(m->E_end(0.l,c)), c->_zero ));

  if( cc->positive) {
    if (_tr_fill < 0) { untested();
      trace1(("DEV_BUILT_IN_RCD::tr_stress fill is negative: " +
            short_label()).c_str() ,  _Ccgfill->get_total() );
    }
    if (value() * involts() < -2e-1) {
      error(bTRACE, "DEV_BUILT_IN_RCD::tr_stress input %s is negative: %f."
          " overshoot?\n", long_label().c_str(), value()*involts() );
    }
  } else {
  }

  /*----------------------------------------------------------------------------*/

  long double fill = _tr_fill;

  fill = _tr_dfill + _Ccgfill->tt() + c->_zero; // shifted, unshifted?

  // trace5("DEV_BUILT_IN_RCD::tr_stress", _tr_dfill, fill, stress, rcd->involts(), iteration_number());
  assert(fill<1.1 && fill>-.1);
  if (fill >= 1.001 ){ untested();
    error(bDANGER, "DEV_BUILT_IN_RCD::tr_stress %s fill %LE big stress=%f\n", long_label().c_str(), _tr_fill, stress );
    fill = 1;
  }

  //double  fill = _n[n_ic].v0();
  assert (fill==fill);

  //tr_stress...
  long double newfill;
  switch(_sim->_stepno){ incomplete();
    case 0:
    case 1:
    default:
      // trace1("DEV_BUILT_IN_RCD::tr_stress calling __step", _sim->_stepno);
      assert(fill<1.1 && fill>-.1);
      newfill = m->__step(stress, fill, h, c);
      assert(newfill<= 1. && newfill>= 0);
  }
  assert( newfill > -0.01 || !cc->positive);
  if (newfill <= 0 && cc->positive) { untested();
    newfill = 0.0;
  } else if(newfill > 1.000001){ untested();
    error(bDANGER, ("* RCD_V3 %f too big\n" + long_label() ).c_str() , newfill );
    newfill=1;
  }
  assert(newfill==newfill);

  _tr_dfill += newfill - _tr_fill;
  _nodes[n_p].a_()->trhack = (double) _tr_dfill;
  // trace7("DEV_BUILT_IN_RCD::tr_stress ", fill, h, (newfill-fill)/h, newfill-fill, _tr_dfill, newfill, _tr_fill);
  _tr_fill = newfill;
  assert(is_number(_tr_fill));
  assert(_tr_dfill<1.1);
  assert(_tr_fill<=1.);
  assert(_tr_dfill>-1.1);

  if (_tr_dfill>1) {
    assert(false);
  } else if (_tr_dfill>.5) {
    error(bWARNING, "at %f: _tr_fill very high in %s: %f\n", _sim->_Time0, long_label().c_str(), (double)_tr_dfill);
  }

  assert(_tr_fill<=1.);
  assert(_tr_fill>=0.);

  assert(h > 0);
}
/*----------------------------------------------------------------------------*/
double DEV_BUILT_IN_RCD::involts() const
{
  double v = _n[n_u].v0() - _n[n_b].v0();
  assert(is_number(v));
  assert(fabs(v) < 1e5);
  return v;
}
/*----------------------------------------------------------------------------*/
// FIXME: move pred/corr to here.
/*----------------------------------------------------------------------------*/
void DEV_BUILT_IN_RCD::tr_stress_last() // "tt_review?".
{
  ADP_NODE* a = _n[n_p].a_();
  trace6("DEV_BUILT_IN_RCD::tr_stress_last", short_label(), _tr_fill,
      _sim->_adp_nodes, _sim->_time0,
      _Ccgfill->tt()+_tr_dfill, _sim->tt_iteration_number());
  assert(_n[n_p].tt() == _n[n_p].tt());
  assert(_n[n_p].tt() == a->tt());
  assert(a);
  const COMMON_BUILT_IN_RCD* c = static_cast<const COMMON_BUILT_IN_RCD*>(common());
  assert(c);
  assert(c->model());
  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(c->model());
  assert(m);
  assert(c->sdp());


  assert(is_number(_tr_dfill));
  _tr_fill = _Ccgfill->tt() + c->_zero + _tr_dfill; // shifted
  assert(is_number(_tr_fill));
  _nodes[n_p].a_()->trhack = (double) _tr_dfill;

  if( _sim->phase()==p_PD){
        assert(!a->tr_lo);
        assert(!a->tr_hi);
  }
  trace2("", a->tr_lo, a->tr_hi );
  // assert( a->tr_lo <= a->tr_hi );
  if(_tr_fill>1){
    error(bTRACE, "at %f: something wrong with _tr_fill in %s\n", _sim->_Time0, long_label().c_str());
    _tr_fill = 1;
  }

  trace5("DEV_BUILT_IN_RCD::tr_stress_last s ", _n[n_p].a_()->tt(),
      _n[n_p].a_()->get_tr(), a->get_total(), _n[n_p].m_(), _n[n_p].t_() );
  // fixme. move to common.
  assert(is_number(_tr_fill));
  assert(is_number(a->tt()));
  double uin_eff;
  try {
    // calculate udc (will be put into _c->tr() ? )
    m->do_tr_stress_last(_tr_fill,_Ccgfill, this);
    uin_eff = a->tr();
  } catch (Exception &e) { untested();
    error(bDANGER, "%s\n", long_label().c_str());
    throw(e);
  }
  assert(is_number(_tr_fill));

  if( _sim->phase()==p_PD ){
    assert(!a->tr());
  }

  if (uin_eff==0) {
  }else if ((uin_eff < a->tr_lo) || (uin_eff > a->tr_hi )) { untested();
    if ( uin_eff - a->tr_lo < -1e-20 || a->tr_hi - uin_eff < -1e-20 ){ untested();
    error(bDANGER, "Time %E: %s order broken, should be %E < %E < %E, is %E %E\n",
        (double)_sim->_Time0,
        long_label().c_str(),
        a->tr_lo, uin_eff, a->tr_hi,
        (uin_eff - a->tr_lo ), (a->tr_hi - uin_eff) );
    }
    a->set_new_order(1); // hmm better 1?
  } else { // tested.
  }

  assert( a->tr_lo <= a->tr_hi );

  if (a->tr_lo > a->tr() ) {
    a->tr() = a->tr_lo ;
    a->set_new_order(1);
  }else if ( a->tr() > a->tr_hi ) { untested();
    a->tr() = a->tr_hi ;
    a->set_new_order(1);
  } else { // tested.
  }

  trace3("DEV_BUILT_IN_RCD::tr_stress_last done",  a->tr_lo ,  a->tr(), a->tr_hi ) ;

  assert(is_number(a->tr()));
  assert(is_number(a->tt()));
  //
  // tt_value not needed for rollback.
  _n[n_p].tt() = double(_tr_fill) - c->_zero; // shifted.
  assert(_n[n_p].tt() > -.01);
  if( _sim->phase()==p_PD ){
    assert(!a->tr());
  }
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_RCD::tt_begin()
{
#ifndef NDEBUG
  _stressiter=0;
#endif
  ADP_NODE* a = _n[n_p].a_();
  assert(a); USE(a);
  const COMMON_BUILT_IN_RCD* c = prechecked_cast<const COMMON_BUILT_IN_RCD*>(common());
  trace1("DEV_BUILT_IN_RCD::tt_begin()", _sim->_tt_uic);

  if (_sim->_tt_uic) {
  } else {
    trace3("DEV_BUILT_IN_RCD::tt_begin()", c->_zero, _n[n_p].m_(), hp(_Ccgfill));
    _n[n_p].tt() = 0; // shifted
  }
  assert(is_number(a->tt()));

  _tr_fill = _n[n_p].tt() + c->_zero;
  assert(_tr_fill<=1.001);
  assert(is_number(_tr_fill));
  if(_tr_fill>1){untested();
    _tr_fill = 1;
  }else{
  }

  _tr_dfill = 0;
  _nodes[n_p].a_()->trhack = (double) _tr_dfill;
  _Ccgfill->set_tr(-inf);

  _Ccgfill->tr_lo = inf;
  _Ccgfill->tr_hi = -inf;

  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(c->model());
  m->do_tt_prepare(this);

  _lasts = -inf;
//  q_eval(); untested();

  trace4("DEV_BUILT_IN_RCD::tt_begin done", long_label(), c->_zero, _sim->_tt_uic, _n[n_p].tt() );
}
///*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_RCD::tr_restore()
{
  ADP_NODE* a = _n[n_p].a_();
  assert(a);
  a->trhack = 0;
}
///*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_RCD::tr_begin()
{
  BASE_SUBCKT::tr_begin();
  trace3("DEV_BUILT_IN_RCD::tr_begin()", _n[n_u].v0(), _n[n_b].v0(), _n[n_u].m_());

  const COMMON_BUILT_IN_RCD* c = prechecked_cast<const COMMON_BUILT_IN_RCD*>(common());
  assert(c->_zero); // 0 is very unlikely...
  _tr_fill = _Ccgfill->tt() + c->_zero; // shifted
  if(_tr_fill>1){untested();
    _tr_fill = 1;
  }else{
  }
  _lasts = -inf; // must be -inf, so time0==0 works out
  _Ccgfill->tr_hi = -inf;
  _Ccgfill->tr_lo = inf;

//  _Ccgfill->tr() = 0; untested();

  // BUG: only do if inside .tw
  //              or trage flag set
  // q_eval(); untested();
  _time0 = _time1 = _sim->_time0;
}
///*--------------------------------------------------------------------------*/
bool DEV_BUILT_IN_RCD::do_tr()
{ untested();
  // is this needed?
  trace3("DEV_BUILT_IN_RCD::do_tr", long_label(), _sim->_time0, _sim->iteration_tag());
  const COMMON_BUILT_IN_RCD* c = static_cast<const COMMON_BUILT_IN_RCD*>(common());
  assert(c); USE(c);
  assert(c->model());
  assert(c->sdp());

  return true;
}
//
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_RCD::precalc_last()
{
  trace1("DEV_BUILT_IN_RCD::precalc_last", value());
  COMPONENT::precalc_last();
  if (!double(value())) { itested();
    set_value(1);
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_RCD::do_tr_stress_last( long double E, ADP_NODE* _c,
    COMPONENT* dd ) const
{
  DEV_BUILT_IN_RCD* d = prechecked_cast<DEV_BUILT_IN_RCD*>(dd);
  const COMMON_BUILT_IN_RCD* cc = prechecked_cast<const COMMON_BUILT_IN_RCD*>(d->common());
  const COMMON_BUILT_IN_RCD* c = cc;
  assert(is_number(_c->tt()));
  double E_old = _c->tt() + cc->_zero;
  double trtime = CKT_BASE::_sim->last_time();
  assert(trtime>=0);

  if (_c->tr() == -inf) _c->tr() = ( _c->tr_lo+_c->tr_hi)/2.0;
  if (!is_number(_c->tr())) _c->tr() = (_c->tr_lo+_c->tr_hi)/2.0;

  long double uin_eff = _c->tr(); // 0 == current estimate

  double v = dd->value(); assert(v);
  trace8("MODEL_BUILT_IN_RCD_EXP::do_tr_stress_last ", _c->label(), E_old,
      tt_iteration_number(), uin_eff, _sim->_time0, _c->tr_lo, _c->tr_hi, v);

  assert(E_old<E_max);
  assert(E<=E_max);

  assert(E<=1);

  long double E_high;
  long double E_test;
  long double E_low;
  long double uin_high_s = (v>0)? _c->tr_hi*v :  _c->tr_lo*v;
  long double uin_low_s = (v>0)? _c->tr_lo*v :  _c->tr_hi*v;
  assert(uin_high_s>=uin_low_s);
  if( _c->tr_hi == _c->tr_lo) {
    // do a sanity check and return
    // incomplete(); possibly ok (constant stress)
  } else {
  }

  if( uin_low_s < -0.001 && cc->positive ){
    // probably numerical problem in transient
    error(bLOG,"%s s too small: %f value %f\n", dd->long_label().c_str(), (double)uin_low_s, double(v));
  }
  if (uin_high_s < -0.001 && cc->positive) {
    error(bLOG,"s too small: %f, value %f\n", (double)uin_high_s, v);
  }

  E_high = __step( uin_high_s, E_old, trtime, c );
  E_low = __step( uin_low_s, E_old, trtime, c );

  // catch numerical problems in __step
  if (E_high<E_low){
    error(bNOERROR,"%s weird. %f>%f\n", dd->long_label().c_str(), (double)E_low, (double)E_high);
    _c->set_new_order(1); // hmmm
    _c->set_tr(double((uin_low_s+uin_high_s)/(2*v))); // assume low ueff
    return;
  } else if (E_low > E) {
    error(bNOERROR,"%s weird. %f>%f\n", dd->long_label().c_str(), (double)E_low, (double)E);
    // low stress ages more than transient.
    _c->set_new_order(1); // hmmm
    _c->set_tr(double((uin_low_s+uin_high_s)/(2*v))); // assume low ueff
    if ( _c->tr() - _c->tr_lo < -1e-20) { untested();
      unreachable();
    }
    if( _c->tr_hi - _c->tr() < -1e-20 ){ untested();
      unreachable();
    }
    return;
  } else if (E_high < E) {
    error(bNOERROR,"%s weird. %f<%f\n", dd->long_label().c_str(), (double)E_high, (double)E);
    _c->set_new_order(1); // hmmm
    _c->set_tr(double((uin_low_s+uin_high_s)/(2*v))); // assume low ueff
    if ( _c->tr() - _c->tr_lo < -1e-20 ){ untested();
      unreachable();
    }
    if ( _c->tr_hi - _c->tr() < -1e-20 ){ untested();
      unreachable();
    }
    return;
  }

  if (_c->tr_lo > uin_eff) {
    error(bNOERROR,"%s do_tr_stress_last uin_eff too low, %f below %f\n", dd->long_label().c_str(), (double)uin_eff, _c->tr_lo);
    uin_eff =  ( _c->tr_lo+_c->tr_hi)/2.0;
  } else if (uin_eff > _c->tr_hi) {
    error(bNOERROR,"%s do_tr_stress_last uin_eff too high, %f above %f\n", dd->long_label().c_str(), (double)uin_eff, _c->tr_hi);
    uin_eff =  ( _c->tr_lo+_c->tr_hi)/2.0;
  }else{
  }

  if (cc->positive ){
    uin_low_s = max(0.L,uin_low_s);
    if(uin_low_s==0){
      uin_eff=0;
    }
    uin_high_s = max(uin_high_s,uin_low_s);
  }else{untested();
  }

  assert(uin_high_s>=uin_low_s);
  assert(v*uin_eff>=uin_low_s);
  assert(uin_high_s>=v*uin_eff);
  assert(uin_low_s>=-0.001 || !cc->positive);
  assert(uin_high_s>=-0.001 || !cc->positive);

  E_high = __step( uin_high_s, E_old, trtime, c );
  E_test = __step( v*uin_eff, E_old, trtime, c );
  E_low = __step( uin_low_s, E_old, trtime, c );

  assert (E_low>=0);

  if((double(E_high-E_low) < - 1e-18)){ untested();
    error(bDANGER,"%s do_tr_stress_last uin_high=%LE uin_low_s=%LE deltaE= %LE; %LE>%LE\n",
        dd->long_label().c_str(),
        uin_high_s, uin_low_s, E_high - E_low, E_high, E_low );

    // this is very bad...
    throw(Exception("monotony violation in __step"));
  }

  bool tight_bounds = false; // effective uin is close to last (within noise margin)

  if ( E_low <= E && E <= E_high ) {
    // invert linearly
    long double dE = (E-E_low) / (E_high-E_low);
    assert(! (dE<0));

    uin_eff = uin_low_s + (uin_high_s - uin_low_s) * dE;
    trace5("MODEL_BUILT_IN_RCD_EXP::do_tr_stress_last linv", 1-E_low, 1-E, E_high-E_low, uin_low_s, uin_high_s );
    if (!(is_number(uin_eff))) {
      uin_eff = (uin_high_s+uin_low_s)/2.; // stupid fallback
    } else if ((uin_eff <= uin_high_s) && (uin_low_s <= uin_eff)) {
      tight_bounds = true;
    } else { untested();
      error(bDANGER, "unreachable? %E %E %E %E %E\n", (double)uin_eff, (double)uin_high_s, (double)uin_low_s, 
          double(uin_high_s - uin_eff), double(uin_eff - uin_low_s));
      assert(false);
    }
    uin_eff /= v;
  }
  //  error(bDANGER,"%s weird uin_high_s=%LE uin_low_s=%LE dus=%LE deltaE= %LE; %LE>%LE %LE %LE uin_eff_s %LE E %LE\n",
  //      dd->long_label().c_str(),
  //      uin_high_s, uin_low_s,
  //      uin_high_s-uin_low_s,
  //      E_high - E_low, E_high, E_low ,E-E_low, E_high-E, uin_eff*v, E);

  long double uin_eff_s = uin_eff * v;
  if(cc->positive && uin_eff_s<0){ untested();
    error(bWARNING, "uin_eff_s too small: %f in [%f, %f]\n", (double) uin_eff_s, _c->tr_lo, _c->tr_hi );
    assert(uin_eff_s >= -.001 || !cc->positive);
  }
  if (tight_bounds){
    uin_eff *= v;
    trace3("MODEL_BUILT_IN_RCD_EXP::do_tr_stress_last tight", uin_low_s, uin_eff_s, uin_high_s);
    assert(double(uin_eff_s)<=double(uin_high_s));
    assert(double(uin_eff_s)>=double(uin_low_s));
//        FIXME: obsolete call? compare
    uin_eff_s = d->seff_bisect( uin_eff_s, (long double)E_old, E, double(uin_low_s), double(uin_high_s));
    assert(uin_low_s <= (double)uin_eff_s);
    assert(uin_high_s >= (double)uin_eff_s);
    trace2("MODEL_BUILT_IN_RCD_EXP::do_tr_stress_last tight", _c->tr_lo, _c->tr_hi);
    assert(uin_eff_s >= -.001 || !cc->positive);
  }else{
    try{
      trace3("MODEL_BUILT_IN_RCD_EXP::do_tr_stress_last", uin_low_s, uin_eff_s, uin_high_s);
      uin_eff_s = d->__uin_iter(uin_eff_s, (long double)E_old, E, double(uin_low_s), double(uin_high_s));
      assert(uin_eff_s >= -.001 || !cc->positive);

      assert(uin_eff_s >= uin_low_s);
      assert(uin_eff_s <= uin_high_s);

    } catch (Exception &e) { untested();
      error(bDANGER, "Exception in %s\n", long_label().c_str());
      throw(e);
    }
    trace4("MODEL_BUILT_IN_RCD_EXP::do_tr_stress_last iteration", E_old, E, uin_eff, E-E_old);
  }
  uin_eff = uin_eff_s/v;

  if(uin_eff > _c->tr_hi){ untested();
    error(bWARNING, "uin_eff too big??  %f in [%f, %f] \n",  (double) uin_eff,  _c->tr_lo, _c->tr_hi );
    itested();
    uin_eff = _c->tr_hi;
  }

  if (uin_eff < _c->tr_lo) {
    uin_eff = _c->tr_lo;
    uin_eff_s = uin_eff*v;
    //assert(v*uin_eff >= -.001 || !cc->positive);
  }

  if (cc->positive && uin_eff*v<0) {
    uin_eff = 0;
  }

  // sanitycheck (monotonicity)
  E_test = this->__step( v*uin_eff, E_old, trtime, c );


  if ( tight_bounds && ((E_low > E_test) || ( E_test > E_high) )){ untested();
    trace3("MODEL_BUILT_IN_RCD_EXP::do_tr_stress_last monotonicity check ",
        uin_low_s, uin_high_s, uin_eff );
    assert( v*uin_eff<=uin_high_s && uin_low_s<=v*uin_eff);
    trace5("MODEL_BUILT_IN_RCD_EXP::do_tr_stress_last ", 1-E, 1-E_test,
        E_test-E_low, E_high-E_test, tight_bounds );

    if (  E_test - E_high <= 1e-30  && E_low - E_test <= 1e-30  ) { untested();

    } else { untested();
      error(bDANGER, "rcd E order problem, %LE, %LE, %LE\n", E_low ,E_test , E_high);
    }
    //assert( E_test - E_high <= 1e-30 );
    //assert( E_low - E_test <= 1e-30);
  }

  uin_high_s = max ( v*uin_eff + OPT::abstol, v*uin_eff * (1 + OPT::reltol) );
  uin_low_s  = min ( v*uin_eff - OPT::abstol, v*uin_eff * (1-OPT::reltol) );

  // FIXME: get rid of __ prefix
  E_high = __step( uin_high_s, E_old, trtime, c );
  E_low  = __step( uin_low_s, E_old, trtime, c );

  if( ( E_old < E_high ) && ( E_low <= E_old )) {
    _c->set_new_order(1);
  }else if (tight_bounds) {
    // _c->set_new_order(1);
  }else{
  }

  if (double(E_high-E_low) < - 1.5e-19){ untested();
    error( bDANGER, "MODEL_BUILT_IN_RCD_EXP:: sanitycheck (delta %LE ) in %s\n",
        E_high - E_low, dd->long_label().c_str());
    error( bDANGER, "MODEL_BUILT_IN_RCD_EXP:: sanitycheck (abs %LE ) in %s\n",
        E_high, dd->long_label().c_str());
    throw(Exception("sanitycheck failed"));
  }

  //assert ( uin_eff == _c->get_tr());
  //
  //
  //
  assert(is_number(E_high-E_low));
  assert(is_number(E));
  //assert( (double) E_high-(double) E_low >= 0);
  ///assert((double) E_high>= (double) E_low);
  if (!(E_low <= E + 1.5e-19 || double(E)==1.0 || double(E_high)==1.0 || !tight_bounds)){ untested();
    error(bWARNING,"someting wrong with bounds in %s, %LE, %LE, %LE\n", dd->long_label().c_str(),
          E_low, E, E_high);
  }

  if(E > E_high && E!=1){
    // this is due to double <=> long double discrepancy
    trace1("MODEL_IN_RCD_EXP::", tight_bounds);
    error( bNOERROR, "%s: Sanitycheck failed ( %LE =E >  E_high=%LE ) del %LE\n",
        dd->long_label().c_str(),
        E, E_high, E_high-E);
  }

  _c->set_tr_noise(0); //(double)E_high-(double)E_low);
  _c->set_tr((double)uin_eff);
  if (CKT_BASE::_sim->tt_iteration_number()>1) {
    if ((_c->tr()-_c->tr1())*(_c->tr1()-_c->tr2())<=0) {
      _c->set_new_order(1);
    }
  }

  trace2("MODEL_BUILT_IN_RCD_EXP::do_tr_stress_last done", _c->get_tr_noise(), uin_eff);
  assert(is_number(uin_eff));

  if ( uin_eff - _c->tr_lo < -1e-20 || _c->tr_hi - uin_eff < -1e-20 ){
    unreachable();
  }
}
/*--------------------------------------------------------------------------*/
long double MODEL_BUILT_IN_RCD::__uin_iter(long double& s, double E_old, long double E_in,
    double bound_lo, double bound_hi, const COMPONENT* dd ) const
{ untested();
  assert(0); //dead code;
  long double E = (long double) E_in;
  const COMMON_BUILT_IN_RCD* c = prechecked_cast<const COMMON_BUILT_IN_RCD*>(dd->common());
  const COMMON_BUILT_IN_RCD* cc = c;
  const MODEL_BUILT_IN_RCD* m = dynamic_cast<const MODEL_BUILT_IN_RCD*>(c->model());
  const DEV_BUILT_IN_RCD* d = prechecked_cast<const DEV_BUILT_IN_RCD*>(dd);
  trace5("COMMON_BUILT_IN_RCD::__uin_iter: ", s, E_old, E, E-E_old,  CKT_BASE::_sim->last_time() );
  assert (E<1.000001);

  if (E>1) { untested();
    trace0("COMMON_BUILT_IN_RCD::__uin_iter aligned E");
    E = 1;
  }
  double h = BASE_SUBCKT::_sim->last_time();

  long double Euin = 0;
  long double Euin_alt = 0;
  if (E < 0.) { untested();
    E = 0.;
  }else{
  }
  //  E = max(E,0.0);
  //if(E<1e-12) return 0;

  assert(( -1 < E ) && (E < 2) ) ;

  long double res=1;   // dx
  long double deltaE =1; // df 
  long double deltaE_alt = 0; // dfold
  long double damp=1;
  long double fu =1; // function we try to find zero of, at u
  long double cres; // cut res
  long double dx_res=1;
  unsigned i=0;
  int hhack=0;
  long double Edu=0;

  double ustart = (double) s;
  bool B=true;

  assert(s>=-0.001 || !cc->positive);
  Euin = m->__step( s, E_old, h, c );
  if(!is_number(Euin)) { untested();
    error( bDANGER, "COMMON_BUILT_IN_RCD::__uin_iter pl cannot evaluate E "
        "at s=%LE (E_old=%E) %i\n", s, E_old, CKT_BASE::_sim->tt_iteration_number());
    assert(false);
  }
  unsigned iterlimit=50;
  if (!_sim->tt_iteration_number()){ untested();
    iterlimit = 500;
  }

  double Esign;
  bool S=false;

// fixme: collect min_{f positive}(s) and max_{f negative}(s)
  while( !S || B ) { untested();
    i++;
    trace7("COMMON_BUILT_IN _RCD::__uin_iter loop", (double)s, (double)res, (double)deltaE, Edu, E, i,Euin);
    if( i >= iterlimit ){ untested();
      d->_iter_newton = i;
      error( bDANGER, "newton does not converge after %u in [%f,%f]: %s in tt %i, s %f", i, dd->long_label().c_str(), bound_lo, bound_hi,_sim->tt_iteration_number(), (double)s );
      throw(Exception("does not converge after %u in [%f,%f]: %s in tt %i", i, dd->long_label().c_str(), bound_lo, bound_hi,_sim->tt_iteration_number() ));
    }
    if(!is_number(s)){ untested();
      error( bDANGER, "COMMON_BUILT_IN_RCD::__uin_iter s wrong %E diff "
          "%E loking for %E \n", Euin, Edu, E );
      assert(false);
      return( inf );
    }
    Edu = m->__dstepds(s, E_old, c); // ??
    if(!is_number(Edu) ){ untested();
      untested();
      error( bDANGER, "COMMON_BUILT_IN_RCD::__uin_iter step %i:%i Edu nan at %LE Euin=%LE C=%LE diff "
          "%LE looking for %E, start %E res %LE\n", CKT_BASE::_sim->tt_iteration_number(),i,
             s, Euin, 1-Euin, Edu, E, ustart, res  );
      // assert(false);
      if (i==1) { untested();
        // first iteration no reliable res
        trace2("Guessing the nan reason",Euin,s);
        if ( Euin >= 0.5 )
        { untested();
          s = s/2; // Guess: s probably much too large
        } else { untested();
          s = s*2; //
        }
      } else {  // next iterations, try half step size 
        res /= 2.0;
        s += res;
        s = std::max(s,0.0L);
        untested();
        hhack++;
      }
      continue;
    }
    assert(Edu>=0);
    if((Edu==0 || fabs(Edu) < 1E-150 ) ){ untested();
      untested();
     // error( bDANGER, "COMMON_BUILT_IN_RCD::__uin_iter step %i:%i Edu 0 at %LE Euin=%LE C=%LE diff "
     //     "%LE looking for %LE, start %E res %LE\n", CKT_BASE::_sim->tt_iteration_number(),i,
     //     s, Euin, 1-Euin, Edu, E, ustart, res  );
      trace0("COMMON_BUILT_IN_RCD::__uin_iter Edu 0");
// assert(false);
      Edu=1;
    }
    assert (is_number (Euin));
    assert (is_number (Edu));
    fu = Euin-E;
    res = damp*fu/Edu;       // dx, Das ist die Differenz in x also: delta x 

    if(!is_number(res) && fabs(Edu) < 1E-150 ) { untested();
      untested();
      trace2("COMMON_BUILT_IN_RCD::__uin_iter",res, Edu);
      Edu=Edu*1E10;
      res = fu/Edu; 
    }

    assert(is_number(res));

    cres= std::min( 1.L,res);
    cres= std::max(-1.L,res);

    dx_res = s;
    double uin1 = double(s);
    s -= cres;
    bool edge = false;
    if(( (double) s > bound_hi )){ untested();
      edge = true;
      s = bound_hi;
    }else if(( (double) s < bound_lo )){ untested();
      edge = true;
      s = bound_lo;
    }
    USE(edge);

    assert(is_number(s));
    if( (s<-0.0000) && cc->positive ) { untested();
      untested();
      trace1( "COMMON_BUILT_IN_RCD::__uin_iter neg s ", s );
      cres = s;
      s = .00;
    }
    dx_res -= s; // Effektives dx da es durch Numerik kaputt gehen kann


    assert(s>=-0.001 || !cc->positive);
    Euin = m->__step( s, E_old, h, c );
    deltaE_alt = deltaE;
    deltaE = Euin - Euin_alt; // Effektive Veraenderung
    if (deltaE > 0) Esign=1; else Esign=-1;
    Euin_alt = Euin;

    // monotony check.
    if(Esign*cres > 0 ){ untested();
      d->_iter_newton = i;
      throw(Exception("monotony violation"));
    }

    if(deltaE * deltaE_alt < 0){ untested();
      damp *= 0.8;
      S = true;
    } else { untested();
      S = 0;
      damp *= 1.2;
    }

    damp=max(damp,0.01L);
    //damp=min(damp,5.L);


    if( !is_number( Euin ) ){ untested();
      error( bDANGER, "COMMON_BUILT_IN_RCD::__uin_iter cannot evaluate E "
          "at s=%LE (E_old=%E) %i:%i\n", s, E_old, CKT_BASE::_sim->tt_iteration_number(), i);
      assert(is_number(Euin)); break;
    }

    // FIXME: use bounds to zero step and break.
    if (Euin > E && s==0 && cc->positive ){ untested();
      break;
    }

  //double reltol = pow(OPT::reltol,1.2);
  //double abstol = OPT::abstol/10.0;

    B = conchk(double(s),uin1);
  }
  d->_iter_newton=i;

  assert(s>=-0.001 || !cc->positive);
  return s;
}
/*-------------------------------------------------------------------------------*/
TIME_PAIR DEV_BUILT_IN_RCD::tt_review()
{
  const COMMON_BUILT_IN_RCD* c = prechecked_cast<const COMMON_BUILT_IN_RCD*>(common());
  assert(c);

  if (_sim->_Time0<=0) {
    return TIME_PAIR(_sim->_dT0, NEVER);
  }

  { // move to tt_accept eventually
    double ueff_guessed = _Ccgfill->tr(_sim->_Time0);
    double ueff_sim = _Ccgfill->tr();

    long double E_from_ueff;

    if (_Ccgfill->order()==1) {
      long double E_sim = _Ccgfill->tt() + c->_zero;
      E_from_ueff = extrapolate(&ueff_sim);
      assert(is_number(E_sim));
      trace5("correct", _sim->_Time0, ueff_guessed, ueff_sim, E_from_ueff, E_sim);
      if(OPT::ttcorr){
        _Ccgfill->tr() = (ueff_guessed + ueff_sim)*.5;
        long double E_new = (E_from_ueff + E_sim)*.5;
        double corr = (double) E_new - c->_zero;
        trace2("correct replacing", _Ccgfill->tt(), corr);
        assert(corr>-.1);
        _Ccgfill->tt() = corr;
      }
    }else{
    }
  }

  double delta = fabs(_Ccgfill->tr1() - _Ccgfill->tr());
  double tol = OPT::tttol * (1./14.);

  if (_Ccgfill->order()==0) {
    /// uuh ooh let simulator decide?
    return TIME_PAIR();
  } else if (delta) {
    error(bNOERROR, "tt step contr" + long_label()
        + " old " + ::to_string(_Ccgfill->tr1()) + " new " + ::to_string(_Ccgfill->tr()) + '\n');
    double timestep = _sim->_dT0 * tol / delta;
    _ttgain = timestep/_sim->_dT0;
    _ttfuture = tt_review_check_and_convert(timestep);
    return TIME_PAIR(_ttfuture,NEVER);
  } else {
    return TIME_PAIR();
  }
}
/*-------------------------------------------------------------------------------*/
TIME_PAIR  DEV_BUILT_IN_RCD::tr_review()
{
  q_accept();
  return TIME_PAIR();
}
/*-------------------------------------------------------------------------------*/
void DEV_BUILT_IN_RCD::tt_regress()
{
  _lasts = -inf;
  _Ccgfill->tr_hi = -inf;
  _Ccgfill->tr_lo = inf;
}
/*-------------------------------------------------------------------------------*/
void DEV_BUILT_IN_RCD::tt_advance()
{
  trace1("DEV_BUILT_IN_RCD::tt_advance", _n[n_p].tt());

  // hmmm
  const COMMON_BUILT_IN_RCD* c = prechecked_cast<const COMMON_BUILT_IN_RCD*>(common());
  USE(c);
  assert(c->_zero); // 0 is very unlikely...

//  _tr_fill = _Ccgfill->tt() + c->_zero; // tt is not valid (should be?!)
//   assert(is_number(_tr_fill));

//  tr_begin();?
  _lasts = -inf;
  _Ccgfill->tr_hi = -inf;
  _Ccgfill->tr_lo = inf;
  // q_eval();
}
/*--------------------------------------------------------------------------*/
#define T long double
template<>
T DEV_BUILT_IN_RCD::__uin_iter<T>(T& uin, T E_old, T E_in, double bound_lo, double bound_hi) const
{
  const DEV_BUILT_IN_RCD* d = this;
  const COMMON_COMPONENT* cc = common();
  const MODEL_BUILT_IN_RCD* m = dynamic_cast<const MODEL_BUILT_IN_RCD*>(cc->model()); USE(m);
  assert(d);

  long double res;
  d->_iter_newton = 0;
  d->_iter_bisect = 0;

  try{
    //res = m->__uin_iter( uin, (double) E_old, E_in, bound_lo, bound_hi, this );
    res = d->seff_bisect( uin, E_old, E_in, bound_lo, bound_hi);

  } catch( Exception e ){ untested();
    assert(uin<=bound_hi);
    assert(uin>=bound_lo);
    res = seff_bisect( uin, E_old, E_in, bound_lo, bound_hi);
  }

  return res;
}
/*--------------------------------------------------------------------------*/
template<>
T DEV_BUILT_IN_RCD::seff_bisect<T>(T& uin, T E_old, T E_in, double bound_lo, double bound_hi) const
{
  assert(bound_lo<=bound_hi);
  assert(bound_hi>=double(uin));
  assert(bound_lo<=double(uin));
  const COMMON_COMPONENT* cc = common();
  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(cc->model());
  bool converged = false;
  const DEV_BUILT_IN_RCD* d = this;

  double dT = BASE_SUBCKT::_sim->last_time();
  trace5("bisect", long_label(), dT, E_old, bound_hi, bound_lo);

  if(bound_lo==bound_hi){
    return uin;
  }

  T hi = bound_hi;
  T lo = bound_lo;
  assert(hi>lo);

  T middle = (hi + lo)/2;
  T E_m = m->__step( middle, (long double)E_old, dT, cc );
#ifndef NDEBUG
#ifdef DO_TRACE
  T Euin_hi = m->__step( hi, (long double)E_old, dT, cc );
  T Euin_lo = m->__step( lo, (long double)E_old, dT, cc );
  USE(Euin_hi);
  USE(Euin_lo);

  if(Euin_hi < E_in){ untested();
  }
  if(Euin_lo > E_in){ untested();
  }
#endif
#endif

  unsigned iter = 0;
  double tolex = 5.;
  while(!converged){

    if (iter++>1000) { untested();
      cerr<<converged<<"\n";
      error(bDANGER, "seff_bisect fail %s, %E, %E, %E, %E, %E\n", d->long_label().c_str(),
            double(lo), double(middle), double(hi), std::abs(double(hi)-double(lo)),
             pow(OPT::abstol, 1.) );
      break;
    }
    assert(hi>=lo);

    if (E_in < E_m) {
      hi = middle;
    } else if (E_m <= E_in) {
      lo = middle;
    }else{ untested();
      unreachable();
    }

    converged = conchk(hi,lo, pow(OPT::abstol, 1.),
                              pow(OPT::reltol, tolex));
    T oldmiddle = middle;
    middle = (hi + lo)/2.;
    T E_m1 = E_m;
    E_m = m->__step( middle, E_old, dT, cc);
    if((oldmiddle-middle)*(E_m1-E_m) < 0){
      error(bWARNING,"%s: bisect monotony violation at %f, delta %f\n",
          d->long_label().c_str(),
          (double)middle, double(oldmiddle-middle));
      break;
    }
  }
  trace5("bisect done ", long_label(), E_old, E_m, dT, middle);
  d->_iter_bisect = iter;
  return middle;
}
#undef T
/*-------------------------------------------------------------------------------*/
MODEL_BUILT_IN_RCD_NET::MODEL_BUILT_IN_RCD_NET(const BASE_SUBCKT* p)
  : MODEL_BUILT_IN_RCD(p){ }
/*-------------------------------------------------------------------------------*/
MODEL_BUILT_IN_RCD_NET::MODEL_BUILT_IN_RCD_NET(const MODEL_BUILT_IN_RCD_NET& p)
  : MODEL_BUILT_IN_RCD(p){ }
/*-------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------*/
