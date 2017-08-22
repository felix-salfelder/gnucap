/*
 *
 * (c)2010 felix salfelder
 */

// generic RCD cell (exponential taus)

#include "e_aux.h"
#include "e_storag.h"
// #include "d_mos_base.h"

#include "globals.h"
#include "e_elemnt.h"
#include "e_adp.h"

#include "d_rcd_sym.h"

using std::max;
using std::min;

class DEV_BUILT_IN_RCD;

class MODEL_BUILT_IN_RCD_EXP : public MODEL_BUILT_IN_RCD_SYM { //
  protected:
    explicit MODEL_BUILT_IN_RCD_EXP(const MODEL_BUILT_IN_RCD_EXP& p);
  public:
    explicit MODEL_BUILT_IN_RCD_EXP(const BASE_SUBCKT* p);
    //virtual void do_tt_prepare(COMPONENT*)const;
    virtual void do_precalc_last(COMMON_COMPONENT* , const CARD_LIST* par_scope)const;
    virtual bool v2() const{return true;}
    // ~MODEL_BUILT_IN_RCD_EXP() : ~MODEL_BUILT_IN_RCD {}
    void do_stress_apply( COMPONENT* d ) const;
    // void do_tr_stress( const COMPONENT*) const;
    std::string dev_type()const ;
    void      set_dev_type(const std::string& nt )
    {
     trace0(("MODEL_BUILT_IN_RCD_EXP::set_dev_type() " + nt).c_str()); 
    };
    CARD* clone()const {return new MODEL_BUILT_IN_RCD_EXP(*this);}
    void do_expand( COMPONENT*) const;
//    ADP_NODE_RCD* new_adp_node(const COMPONENT*) const;
//    region_t region(const COMPONENT*) const;
    double P( const COMPONENT* brh) const;
  private:
    double __dRe(double uin, const COMMON_COMPONENT* cc)const; // unneeded?
    double __dRc(double uin, const COMMON_COMPONENT* cc)const;
    double __Re(double uin, const COMMON_COMPONENT* cc)const;
    long double __Re(long double uin, const COMMON_COMPONENT* cc)const;
    double __Rc(double uin, const COMMON_COMPONENT* cc)const;
    long double __Rc(long double uin, const COMMON_COMPONENT* cc)const;
    double __Ge(double uin, const COMMON_COMPONENT* cc)const;
    double __tau(double uin, const COMMON_COMPONENT* cc)const;
  public:
//    virtual void do_tr_stress_last(long double fill, ADP_NODE*, COMPONENT*  ) const; //base
  private:
    double __E(double uin, const COMMON_COMPONENT* cc)const;

    long double __E(double uin, long double cur, const COMMON_COMPONENT* cc)const;
    long double __dstepds(long double uin, long double cur, const COMMON_COMPONENT* cc)const;
//    template<class T>
//    BUG: _step in MODEL_RCD not generic.
//  protected:
//    long double __step(long double s, long double cur,  double dt, const COMMON_COMPONENT* ) const ; // MODEL_RCD
};
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_RCD_EXP::do_stress_apply( COMPONENT*  ) const
{ untested();
  unreachable();
}
/*--------------------------------------------------------------------------*/
//void DEV_BUILT_IN_RCD_EXP::tr_stress()
//{ untested();
//  unreachable(); // obsolet....
//
//  const COMMON_BUILT_IN_RCD* c = static_cast<const COMMON_BUILT_IN_RCD*>(common());
//  assert(c);
//  assert(c->model());
//  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(c->model());
//  assert(m);
//  assert(c->sdp());
//
//  assert(false);
//
//}
///*--------------------------------------------------------------------------*/
namespace MODEL_BUILT_IN_RCD_DISPATCHER
{ 
  static DEV_BUILT_IN_RCD p4d;
  static MODEL_BUILT_IN_RCD_EXP p4(&p4d);
  static DISPATCHER<MODEL_CARD>::INSTALL
    d3(&model_dispatcher, "rcd_exp", &p4);
}
///*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_RCD_EXP::do_expand(  COMPONENT* c) const
{ untested();
  DEV_BUILT_IN_RCD_SYM* d = dynamic_cast<DEV_BUILT_IN_RCD_SYM*>(c);
  assert(d); USE(d);
  // d->expand();
  
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_RCD_EXP::dev_type()const
{
  return "rcd";
}
/*--------------------------------------------------------------------------*/
MODEL_BUILT_IN_RCD_EXP::MODEL_BUILT_IN_RCD_EXP(const BASE_SUBCKT* p) 
  : MODEL_BUILT_IN_RCD_SYM(p)
{ 
  // uref=1; 
}
/*--------------------------------------------------------------------------*/
MODEL_BUILT_IN_RCD_EXP::MODEL_BUILT_IN_RCD_EXP(const MODEL_BUILT_IN_RCD_EXP& p)  
  : MODEL_BUILT_IN_RCD_SYM(p)
{ 
  
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
#if 0 // cleanup
double MODEL_BUILT_IN_RCD_EXP::__tau(double uin, const COMMON_COMPONENT* cc)const
{ untested();
  double tau_e = __Re( uin, cc);
  double tau_c = __Rc( uin, cc);

  if (tau_c < tau_e){ untested();
    return tau_c  / ( tau_c/tau_e +1 );
  } else{ untested();
    return tau_e  / ( tau_e/tau_c +1 );
  }
}
#endif
/*--------------------------------------------------------------------------*/
double MODEL_BUILT_IN_RCD_EXP::__Re(double s, const COMMON_COMPONENT* c) const
{ untested();
  assert(0); // long..
  const COMMON_BUILT_IN_RCD* cc = dynamic_cast<const COMMON_BUILT_IN_RCD*>(c) ;
  assert(is_number(s));
  assert(is_number(cc->_Re1));
  assert(is_number(cc->_Re0));
  double r = exp( cc->_Re1 * s + cc->_Re0 );
  assert(is_number(r));
  return r;
}
/*--------------------------------------------------------------------------*/
long double MODEL_BUILT_IN_RCD_EXP::__Re(long double s, const COMMON_COMPONENT* c) const
{
  const COMMON_BUILT_IN_RCD* cc = dynamic_cast<const COMMON_BUILT_IN_RCD*>(c) ;
  assert(is_number(s));
  assert(is_number(cc->_Re1));
  assert(is_number(cc->_Re0));
  long double r = expl( cc->_Re1 * s + cc->_Re0 );
  assert(is_number(r));
  return r;
}
/*--------------------------------------------------------------------------*/
double MODEL_BUILT_IN_RCD_EXP::__dRe(double s, const COMMON_COMPONENT* c) const
{ untested();
  const COMMON_BUILT_IN_RCD* cc = dynamic_cast<const COMMON_BUILT_IN_RCD*>(c) ;
  return cc->_Re1 * exp( cc->_Re1 * s + cc->_Re0 );
}
/*--------------------------------------------------------------------------*/
double MODEL_BUILT_IN_RCD_EXP::__Rc(double s, const COMMON_COMPONENT* c ) const
{ untested();
  assert(0); // long
  const COMMON_BUILT_IN_RCD* cc = dynamic_cast<const COMMON_BUILT_IN_RCD*>(c) ;
  return exp( cc->_Rc1 * s + cc->_Rc0 );
}
/*--------------------------------------------------------------------------*/
long double MODEL_BUILT_IN_RCD_EXP::__Rc(long double s, const COMMON_COMPONENT* c ) const
{
  const COMMON_BUILT_IN_RCD* cc = dynamic_cast<const COMMON_BUILT_IN_RCD*>(c) ;
  return expl( cc->_Rc1 * s + cc->_Rc0 );
}
/*--------------------------------------------------------------------------*/
double MODEL_BUILT_IN_RCD_EXP::__dRc(double s, const COMMON_COMPONENT* c ) const
{ untested();
  const COMMON_BUILT_IN_RCD* cc = dynamic_cast<const COMMON_BUILT_IN_RCD*>(c) ;
  return cc->_Rc1 * exp( cc->_Rc1 * s  + cc->_Rc0 );
}
/*--------------------------------------------------------------------------*/
double MODEL_BUILT_IN_RCD_EXP::__Ge(double s, const COMMON_COMPONENT* c ) const
{ untested();
  assert(0);
  const COMMON_BUILT_IN_RCD* cc = dynamic_cast<const COMMON_BUILT_IN_RCD*>(c) ;
  return exp( - cc->_Re1 * s - cc->_Re0 );
}
/*--------------------------------------------------------------------------*/
long double MODEL_BUILT_IN_RCD_EXP::__dstepds(long double uin, long double cur, const COMMON_COMPONENT* cc ) const
{ untested();
  const COMMON_BUILT_IN_RCD* c = dynamic_cast<const COMMON_BUILT_IN_RCD*>(cc);
  long double Re0=c->_Rc0; // swapped
  long double Rc0=c->_Re0; // swapped
  long double Re1=c->_Rc1; // swapped
  long double Rc1=c->_Re1; // swapped
  long double t=_sim->last_time();

  unreachable();

  return 1;

  long double ret =
    -((cur - 1.L)*Rc0*Rc0*Rc0*Re1*t*expl(3.L*(Rc1 + Re1)*uin) - Rc1*Re0*Re0*Re0*cur*t - ((cur - 1.L)*Rc0*Rc0*Rc1*Re0 - 
          (2.L*cur - 1.L)*Rc0*Rc0*Re0*Re1)*t*expl(2.L*(Rc1 + Re1)*uin) - (((2.L*cur - 1.L)*Rc0*Rc1*Re0*Re0
                - Rc0*Re0*Re0*Re1*cur)*t*expl(Rc1*uin) - (Rc0*Rc0*Rc1*Re0*Re0 +
                  Rc0*Rc0*Re0*Re0*Re1)*expl(2.L*Rc1*uin))*expl(Re1*uin) - (Rc0*Rc0*Rc1*Re0*Re0 +
                  Rc0*Rc0*Re0*Re0*Re1)*expl(((2*Rc0*Rc1*Re0 + Rc0*Re0*Re1)*uin*expl(Rc1*uin) +
                      Rc0*t*expl((Rc1 + Re1)*uin) +
                      Re0*t)*expl(-Rc1*uin)/(Rc0*Re0)))
    *expl(-(Rc0*t*expl((Rc1 + Re1)*uin) + Re0*t)*expl(-Rc1*uin)/(Rc0*Re0))
    /(Rc0*Rc0*Rc0*Re0*expl((3.L*Rc1 + 2.L*Re1)*uin) +
                      2.L*Rc0*Rc0*Re0*Re0*expl((2.L*Rc1 + Re1)*uin) + Rc0*Re0*Re0*Re0*expl(Rc1*uin));

  
  assert(ret>=0);
  return(ret);
  /*sage
  ret =
    -((cur - 1)*Rc0^3*Re1*t*e^(3*(Rc1 +
            Re1)*uin) - Rc1*Re0^3*cur*t - ((cur - 1)*Rc0^2*Rc1*Re0 - (2*cur -
              1)*Rc0^2*Re0*Re1)*t*e^(2*(Rc1 + Re1)*uin) - (((2*cur - 1)*Rc0*Rc1*Re0^2
                - Rc0*Re0^2*Re1*cur)*t*e^(Rc1*uin) - (Rc0^2*Rc1*Re0^2 +
                  Rc0^2*Re0^2*Re1)*e^(2*Rc1*uin))*e^(Re1*uin) - (Rc0^2*Rc1*Re0^2 +
                  Rc0^2*Re0^2*Re1)*e^(((2*Rc0*Rc1*Re0 + Rc0*Re0*Re1)*uin*e^(Rc1*uin) +
                      Rc0*t*e^((Rc1 + Re1)*uin) +
                      Re0*t)*e^(-Rc1*uin)/(Rc0*Re0)))*e^(-(Rc0*t*e^((Rc1 + Re1)*uin) +
                      Re0*t)*e^(-Rc1*uin)/(Rc0*Re0))/(Rc0^3*Re0*e^((3*Rc1 + 2*Re1)*uin) +
                      2*Rc0^2*Re0^2*e^((2*Rc1 + Re1)*uin) + Rc0*Re0^3*e^(Rc1*uin))
*/

}
/*--------------------------------------------------------------------------*/
// solve E(uin)-E==0
//long double MODEL_BUILT_IN_RCD_EXP::__uin_iter(long double& uin, double
//    E_old, double E, COMPONENT* dd ) const { const
//  COMMON_BUILT_IN_RCD* c = dynamic_cast<const COMMON_BUILT_IN_RCD*>(dd->common());
//  assert(false);
//  return c->__uin_iter( uin, E_old ,E,0,0,dd ); 
//}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_RCD_EXP::do_precalc_last(COMMON_COMPONENT* ccc, const CARD_LIST*)const
{
  COMMON_BUILT_IN_RCD* cc = prechecked_cast<COMMON_BUILT_IN_RCD*>(ccc);
  COMMON_BUILT_IN_RCD* c=cc;
  assert(cc);
  //const MODEL_BUILT_IN_RCD_EXP* m=this;

  cc->Uref = 0.;

  trace5("MODEL_BUILT_IN_RCD_EXP::do_precalc_last", cc->Uref,
      c->Recommon1,
      c->Recommon0,
      c->Rccommon1,
      c->Rccommon0);

  //double up   =  cc->Recommon0;

  c->_Re1 = cc->Recommon1;
  c->_Re0 = cc->Recommon0;
  c->_Rc1 = cc->Rccommon1;
  c->_Rc0 = cc->Rccommon0;

  c->_zero = (double)E_end_0(c);

  cerr.precision(150);
  trace1("MODEL_BUILT_IN_RCD_EXP::do_precalc_last", cc->_zero);
  cerr.precision(16);

  cc->_wcorr = 1;
  cc->_weight = cc->weight;

  // assert( cc->weight != 0 );
  // assert( cc->_weight != 0 );
  assert( is_number( cc->_Rc1 ) );
  assert( is_number( cc->_Rc0 ) );
  assert( is_number( cc->_Re1 ) );
  assert( is_number( cc->_Re0 ) );
}
/*-------------------------------------------------*/
double MODEL_BUILT_IN_RCD_EXP::P( const COMPONENT* brh) const
{
  const DEV_BUILT_IN_RCD* d = prechecked_cast<const DEV_BUILT_IN_RCD*>(brh);
  const COMMON_BUILT_IN_RCD* cc = prechecked_cast<const COMMON_BUILT_IN_RCD*>(brh->common());
  const COMMON_BUILT_IN_RCD* c = cc;
  assert(cc);
  if( _sim->analysis_is_tt() ){
    if(cc->positive && (d->_Ccgfill->tt() < 0)) {
      error(bTRACE,"RCD_EXP: %d not positive %s. tt: %f zero: %f\n",
          d->_Ccgfill->m_(),
          brh->long_label().c_str(), double(d->_Ccgfill->tt()), double(c->_zero));
    }
    return double (d->_tr_fill - (long double) c->_zero) * c->_weight;
    // return (d->_Ccgfill->tt() - c->_zero) * c->_weight;
  }else{
    assert(is_number( d->_Ccgfill->tt() * c->_weight));
    // return c->_Ccgfill->tt() * c->_weight * c->_wcorr;
    //
    //FIXME. _tr_fill must be part of an ADP_NODE
    trace2("MODEL_BUILT_IN_RCD_EXP::P",  d->_tr_fill,  c->_weight  );
    assert( d->_Ccgfill->tt() <=1 );
    assert( d->_tr_fill <=1 );
    return double((d->_tr_fill - c->_zero) * c->_weight);
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

// vim:ts=8:sw=2:et:
