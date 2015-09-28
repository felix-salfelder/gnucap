/* vim:ts=8:sw=2:et:
 *
 * (c)2010 felix salfelder
 * nonGPL ?
 */


#include "e_aux.h"
#include "e_storag.h"

#define ADD_VERSION
#include "globals.h"
#include "e_elemnt.h"
#include "e_adp.h"
#include "d_rcd_sym.h"

/*--------------------------------------------------------------------------*/
double MODEL_BUILT_IN_RCD_SYM::P( const COMPONENT* brh) const
{
  const DEV_BUILT_IN_RCD* c = prechecked_cast<const DEV_BUILT_IN_RCD*>(brh);
  const COMMON_BUILT_IN_RCD* cc = prechecked_cast<const COMMON_BUILT_IN_RCD*>(c->common());
  if (_sim->analysis_is_tt()) { untested();
    return (cc->_zero + c->_Ccgfill->tt()) * cc->_weight * cc->_wcorr; // shifted
  } else { untested();
    assert(is_number( c->_Ccgfill->tt() * cc->_weight * cc->_wcorr));
    // return c->_Ccgfill->tt() * cc->_weight * cc->_wcorr;
    //
    return double((c->_tr_fill + c->_Ccgfill->tt() + cc->_zero ) * cc->_weight * cc->_wcorr); // shifted
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_RCD_SYM::do_stress_apply( COMPONENT* ) const
{
  if (!_sim->analysis_is_tt()){
//        _Ccgfill->
  } 
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_RCD_SYM::tr_stress()
{
#if 0
  const COMMON_BUILT_IN_RCD* c = static_cast<const COMMON_BUILT_IN_RCD*>(common());
  assert(c);
  assert(c->model());
  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(c->model());
  assert(m);
  assert(c->sdp());
#endif

  unreachable(); //obsolete.
  assert(false);
}
///*--------------------------------------------------------------------------*/
namespace MODEL_BUILT_IN_RCD_DISPATCHER { 
  static DEV_BUILT_IN_RCD_SYM p2d;
  static MODEL_BUILT_IN_RCD_SYM p2(&p2d);
  static DISPATCHER<MODEL_CARD>::INSTALL
    d2(&model_dispatcher, "rcdsym_base", &p2);
}
///*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_RCD_SYM::do_expand( COMPONENT* ) const
{
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_RCD_SYM::dev_type()const
{
  return "rcdsym";
}
/*--------------------------------------------------------------------------*/
MODEL_BUILT_IN_RCD_SYM::MODEL_BUILT_IN_RCD_SYM(const BASE_SUBCKT* p) 
  : MODEL_BUILT_IN_RCD(p){ }
/*--------------------------------------------------------------------------*/
MODEL_BUILT_IN_RCD_SYM::MODEL_BUILT_IN_RCD_SYM(const MODEL_BUILT_IN_RCD_SYM& p)  
  : MODEL_BUILT_IN_RCD(p){ }
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_RCD_SYM::expand() {
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

  trace0("DEV_BUILT_IN_RCD_SYM::expand()");

  if (_sim->is_first_expand()) {
    precalc_first();
    precalc_last();
    // local nodes
    //assert(!(_n[n_ic].n_()));
    //BUG// this assert fails on a repeat elaboration after a change.
    //not sure of consequences when new_model_node called twice.
   #if 0
    if (!(_n[n_ic].n_())) {
      if (false) {
        _n[n_ic] = _n[n_b];
      }else{
        _n[n_ic].new_model_node("." + long_label() + ".ic", this);
      }
    }else{
      if (false) {
        assert(_n[n_ic] == _n[n_b]);
      }else{
        //_n[n_ic].new_model_node("ic." + long_label(), this);
      }
    }
#endif
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_RCD_SYM::do_tr_stress( const COMPONENT*) const
{
  assert(false); // use DEV::tr_stress
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_RCD_SYM::do_tt_prepare( COMPONENT* brh) const
{
  const DEV_BUILT_IN_RCD* c = prechecked_cast<const DEV_BUILT_IN_RCD*>(brh);
  const COMMON_BUILT_IN_RCD* cc = prechecked_cast<const COMMON_BUILT_IN_RCD*>(c->common());
  assert(is_number(cc->_zero));

  trace7( "MODEL_BUILT_IN_RCD_SYM_V2::do_tt_prepare", brh->short_label(),
      -cc->_wcorr, cc->_zero, c->_Ccgfill->tt(), hp(c->_Ccgfill), c->_Ccgfill->m_(), _sim->_adp_nodes);

  if(!_sim->_tt_uic) {
    assert(is_number(cc->_zero));
    assert(is_number(c->_Ccgfill->tt()));
    assert(c->_Ccgfill->tt() == 0); //shifted
  } else {
  }
}
/*--------------------------------------------------------------------------*/
DEV_BUILT_IN_RCD_SYM::DEV_BUILT_IN_RCD_SYM()
  :DEV_BUILT_IN_RCD()
{
  //  _n = _nodes;
  //  attach_common(&Default_BUILT_IN_RCD);

  //  ++_count;
  // overrides
}
/*--------------------------------------------------------------------------*/
int  MODEL_BUILT_IN_RCD_SYM::tt_region(const COMPONENT* brh) const
{ untested();
  const DEV_BUILT_IN_RCD* c = (const DEV_BUILT_IN_RCD*) brh;

  assert(c);
  return ( (c->_Ccgfill)->region() );
}
/*--------------------------------------------------------------------------*/
//double MODEL_BUILT_IN_RCD_SYM::__Edu(double s, const COMMON_COMPONENT* c ) const 
//{
//  const COMMON_BUILT_IN_RCD* cc = dynamic_cast<const COMMON_BUILT_IN_RCD*>(c) ;
//  double te = __Re( s, cc);
//  double tc = __Rc( s, cc);
//  double dte = __dRe( s, cc);
//  double dtc = __dRc( s, cc);
//  double ret =  (dtc * ( te + tc ) - ( dte + dtc ) * tc) / ( dte + dtc ) / ( dte + dtc )  ;
//  assert(is_number(ret));
//  return ret;
//}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
//ADP_NODE_RCD* MODEL_BUILT_IN_RCD_SYM::new_adp_node(const COMPONENT* c) const{
//  return new ADP_NODE_RCD(c);
//}
