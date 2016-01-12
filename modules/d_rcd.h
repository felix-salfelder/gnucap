/* $Id: d_rcd.h,v 1.7 2010-09-07 07:46:22 felix Exp $ -*- C++ -*-
 *
 * RCD device header.
 *
 * (c) 2010 Felix Salfelder
 *
 * GPLv3
 */

#ifndef D_RCD_H_INCLUDED
#define D_RCD_H_INCLUDED

#include "d_diode.h"
#include "e_adp.h"
#include "e_aux.h"
#include "d_rcd.h"
#include "d_bti.h"
#include "u_limit.h"
#include "u_sdp.h"
#include "e_node.h"
#include "e_subckt.h"
#include "e_model.h"
//static bool dummy=false;
 enum {USE_OPT = 0x8000};
/*--------------------------------------------------------------------------*/
using std::cerr;
/*--------------------------------------------------------------------------*/
class SDP_BUILT_IN_RCD
  :public SDP_CARD{
public:
  explicit SDP_BUILT_IN_RCD(const COMMON_COMPONENT* c) : SDP_CARD(c) {init(c);}
  void init(const COMMON_COMPONENT*);
public:
};
/*--------------------------------------------------------------------------*/
class DEV_BUILT_IN_RCD;
class COMMON_BUILT_IN_RCD;
class TDP_BUILT_IN_RCD{
public:
  explicit TDP_BUILT_IN_RCD(const DEV_BUILT_IN_RCD*);
public:
};
/*--------------------------------------------------------------------------*/
class MODEL_BUILT_IN_RCD : public MODEL_CARD {
  friend class COMMON_BUILT_IN_RCD;
  friend class DEV_BUILT_IN_RCD;
  protected:
    explicit MODEL_BUILT_IN_RCD(const MODEL_BUILT_IN_RCD& p);
  public:
    explicit MODEL_BUILT_IN_RCD(const BASE_SUBCKT*);
    ~MODEL_BUILT_IN_RCD() {--_count;}
    /// virtual void     do_expand( COMPONENT*)const {}; // doesnt work, not in use
    //
//    virtual void     do_precalc_last(COMPONENT*, const );
    virtual void     do_precalc_last(COMMON_COMPONENT*, const CARD_LIST*)const;
    virtual void     do_tt_prepare(COMPONENT*)const;
//    virtual ADP_NODE_RCD* new_adp_node(const COMPONENT*) const;
  public: // override virtual
    virtual bool v2() const{return false;}
    virtual double P(const COMPONENT*) const = 0;
    virtual std::string dev_type()const;
    virtual void      set_dev_type(const std::string& nt);
    virtual CARD*     clone()const = 0;//{return new MODEL_BUILT_IN_RCD(*this);}
    void      precalc_first();
    void      precalc_last();
    SDP_CARD* new_sdp(COMMON_COMPONENT* c)const;
    ADP_CARD* new_adp(COMPONENT* c)const;
    void      set_param_by_index(int, std::string&, int);
    bool      param_is_printable(int)const;
    std::string param_name(int)const;
    std::string param_name(int,int)const;
    std::string param_value(int)const;
    int param_count()const {return (10);}
    void tt_eval(COMPONENT*)const;
    bool      is_valid(const COMPONENT*)const;
    void      tr_eval(COMPONENT*)const;
    virtual void      do_stress_apply( COMPONENT*)const {unreachable();}
    virtual void      do_tr_stress( const COMPONENT*) const;
  public: // not virtual
    static int count() {return _count;}
  private: // strictly internal
    static int _count;
  public: // input parameters
    PARAMETER<bool> anneal;	// flag: anneal
    PARAMETER<double> Remodel;	// emit Resistance
    PARAMETER<double> Re1;	//
    PARAMETER<double> Re0;	// emit res
    PARAMETER<double> Rc1;	//.
    PARAMETER<double> Rc0;	// capt res.
    PARAMETER<int> flags;	// 
    PARAMETER<double> uref;	// 
    PARAMETER<int> modelparm;	// just a test
    PARAMETER<bool> norm_uin;	// normalized uin
    PARAMETER<double> iscale;   // scale stress before doing anything.
  public: // calculated parameters
    virtual bool use_net() const {unreachable(); return 1;}
  public: // probes
    virtual region_t region(const COMPONENT* ) const;
    virtual int tt_region(const COMPONENT* ) const;

  protected:
    virtual      double __Re(double uin, const COMMON_COMPONENT* cc)const;
    virtual long double __Re(long double uin, const COMMON_COMPONENT* cc)const
    {unreachable(); return __Re(double(uin), cc);}
    virtual      double __Rc(double uin, const COMMON_COMPONENT* cc)const;
    virtual long double __Rc(long double uin, const COMMON_COMPONENT* cc)const
    {unreachable(); return __Rc(double(uin), cc);}
    virtual double __Ge(double uin, const COMMON_COMPONENT* )const;
    virtual double __dRe(double, const COMMON_COMPONENT* )const { unreachable(); return 0; }
    virtual double __dRc(double, const COMMON_COMPONENT* )const { unreachable(); return 0; }

    virtual long double __step(long double uin, long double cur, double deltat, const COMMON_COMPONENT* ) const;
  
  protected: // functions using the virtual __Re __Rc
    template<class T>
      T __tau(T, const COMMON_COMPONENT* )const; 
    template<class T>
      T __tau_inv(T, const COMMON_COMPONENT* )const ;
//    template<class T> // final state at fixed stress
    long double E_end_0( const COMMON_COMPONENT* c ) const;
    template<class T> // final state at fixed stress
      T E_end(T s, const COMMON_COMPONENT* c ) const;
    
  protected:
    virtual long double __dstepds( long double, long double , const COMMON_COMPONENT*)const { return 0; }
    //virtual double __dstepds( double, double , const COMMON_COMPONENT*)const { return 0; }
    //  dstep/du
    //
  private:
    long double __uin_iter(long double& uin, double E_old, long double E_in, double bound_lo, double bound_hi, const COMPONENT* dd ) const;

  public:
    virtual void do_tr_stress_last( long double , ADP_NODE* , COMPONENT* ) const;
};
/*--------------------------------------------------------------------------*/
template <class T>
inline T MODEL_BUILT_IN_RCD::__tau_inv(T s, const COMMON_COMPONENT* cc)const
{
  const MODEL_BUILT_IN_RCD* m = this;
  assert(m);
  T tau_e = m->__Re( s, cc );
  T tau_c = m->__Rc( s, cc );
  return  ( tau_e + tau_c ) /  tau_c / tau_e;
  if (tau_c < tau_e){
    return ( tau_c/tau_e + 1 ) / tau_c;
  } else if(tau_c>=tau_e){
    return ( tau_e/tau_c + 1 ) / tau_e;
  }
}
///*--------------------------------------------------------------------------*/
template <class T>
inline T MODEL_BUILT_IN_RCD::__tau(T s, const COMMON_COMPONENT* cc)const
{
  const MODEL_BUILT_IN_RCD* m = this;
  assert(m);
  T tau_e = m->__Re( s, cc );
  T tau_c = m->__Rc( s, cc );
  return  ( tau_e * tau_c ) / ( tau_c + tau_e );
  if (tau_c < tau_e){
    return tau_c  / ( tau_c/tau_e + 1 );
  } else{
    return tau_e  / ( tau_e/tau_c + 1 );
  }
//  return tau_c*tau_e / ( tau_e + tau_c );
}
///*--------------------------------------------------------------------------*/
class MODEL_BUILT_IN_RCD_NET : public MODEL_BUILT_IN_RCD {
  protected:
    explicit MODEL_BUILT_IN_RCD_NET(const MODEL_BUILT_IN_RCD_NET& p);
  public:
    explicit MODEL_BUILT_IN_RCD_NET(const BASE_SUBCKT* p);
    // ~MODEL_BUILT_IN_RCD_NET() {  ~MODEL_BUILT_IN_RCD() };
    bool use_net() const { return(1); }
    virtual  void do_stress_apply( COMPONENT* ) const ;
    std::string dev_type()const;
    CARD*     clone()const {return new MODEL_BUILT_IN_RCD_NET(*this);}
    // void     do_expand(const COMPONENT*) const;
    double P(const COMPONENT* )const{ return 17; }
};
/*--------------------------------------------------------------------------*/
class COMMON_BUILT_IN_RCD :public COMMON_COMPONENT{
  public:
    explicit COMMON_BUILT_IN_RCD(const COMMON_BUILT_IN_RCD& p);
    explicit COMMON_BUILT_IN_RCD(int c=0);
    ~COMMON_BUILT_IN_RCD();
    bool     operator==(const COMMON_COMPONENT&)const;
    COMMON_COMPONENT* clone()const {return new COMMON_BUILT_IN_RCD(*this);}
    void     set_param_by_index(int, std::string&, int);
    bool     param_is_printable(int)const;
    std::string param_name(int)const;
    std::string param_name(int,int)const;
    void set_param_by_name(std::string Name, std::string Value);
    std::string param_value(int)const;
    int param_count() const;
    void     precalc_first(const CARD_LIST*);
    void     expand(const COMPONENT*);
    void     precalc_last(const CARD_LIST*);
    std::string name()const {itested();return "rcd";}
    const SDP_CARD* sdp()const {return _sdp;}
    bool     has_sdp()const {untested();return _sdp;}
    static int  count() {return _count;}
  private: // strictly internal
    static int _count;
  public: // input parameters
    PARAMETER<double> perim;	// perimeter factor
    PARAMETER<double> weight;	// cap weight
    PARAMETER<double> Recommon0;	// emit resistance
    PARAMETER<double> Recommon1;	// emit resistance
    PARAMETER<double> Rccommon0;	// capt resistance abs
    PARAMETER<double> Rccommon1;	// capt resistance steepness
    PARAMETER<bool>   positive;	// treat negative stress as zero
    PARAMETER<double> Uref;	// reference voltage
    PARAMETER<double> mu;	// mu
    PARAMETER<double> lambda;	// lambda
    PARAMETER<bool> dummy_capture;// use dummy as capture resistor
    PARAMETER<bool> dummy_emit;	// use dummy as emit resistor
  public: // calculated parameters
    double X;
    double _Re1;
    double _Re0;
    double _Rc0;
    double _Rc1;
    double _weight;
    double _wcorr; // correction for uref fit. FIXME: remove
    double _zero;  // this much fill is dvth=0 ("E0")
    double _lambda;
    SDP_CARD* _sdp;
    double cj_adjusted;	// 
  public: // functions... (protected?)  should be in MODEL??
    //double _Ueff( double ug);
    double Re(double s)const; // just forwarding to model
    double Rc(double s)const; // just forwarding to model
//    virtual double __Re(double s)const;
//    virtual double __Rc(double s)const;
    double tau ( double s) const; 
    double __tau_up ( double ueff ) const;
    long double __uin_iter(long double& uin, double, double, double bhi, double blo,COMPONENT*)const;
    long double seff_bisect(long double& uin, double, double, double bhi, double blo,COMPONENT*)const;

  public: // attached commons
};
/*--------------------------------------------------------------------------*/
class EVAL_BUILT_IN_RCD_GRc : public COMMON_COMPONENT {
private:
  explicit EVAL_BUILT_IN_RCD_GRc(const EVAL_BUILT_IN_RCD_GRc& p)
    :COMMON_COMPONENT(p) {}
public:
  explicit EVAL_BUILT_IN_RCD_GRc(int c=0) :COMMON_COMPONENT(c) {}
  bool operator==(const COMMON_COMPONENT& x)const {return COMMON_COMPONENT::operator==(x);}
  COMMON_COMPONENT* clone()const {return new EVAL_BUILT_IN_RCD_GRc(*this);}
  std::string name()const {untested(); return "EVAL_BUILT_IN_RCD_GRc";}
  void tr_eval(ELEMENT*d)const;
  bool has_tr_eval()const {return true;}
  bool has_ac_eval()const {return false;}
};
/*--------------------------------------------------------------------------*/
class EVAL_BUILT_IN_RCD_Ye : public COMMON_COMPONENT {
private:
  explicit EVAL_BUILT_IN_RCD_Ye(const EVAL_BUILT_IN_RCD_Ye& p)
    :COMMON_COMPONENT(p) {}
public:
  explicit EVAL_BUILT_IN_RCD_Ye(int c=0) :COMMON_COMPONENT(c) {}
  bool operator==(const COMMON_COMPONENT& x)const {return COMMON_COMPONENT::operator==(x);}
  COMMON_COMPONENT* clone()const {return new EVAL_BUILT_IN_RCD_Ye(*this);}
  std::string name()const {untested(); return "EVAL_BUILT_IN_RCD_Ye";}
  void tr_eval(ELEMENT*d)const;
  bool has_tr_eval()const {return true;}
  bool has_ac_eval()const {return false;}
};
/*--------------------------------------------------------------------------*/
class DEV_BUILT_IN_RCD : public BASE_SUBCKT { // BUG: COMPONENT (fix later)
  friend class COMMON_BUILT_IN_RCD;
  friend class MODEL_BUILT_IN_RCD;
  private:
    explicit DEV_BUILT_IN_RCD(const DEV_BUILT_IN_RCD& p);
    double _lasts; // last time tr_stress was called
  public:
    explicit DEV_BUILT_IN_RCD();
    ~DEV_BUILT_IN_RCD() {
      --_count;
      if( _Ccgfill ) {
        //   ADP_NODE_LIST::adp_node_list.erase( _Ccgfill );
      }
      if( _Udc )  {
        //   ADP_NODE_LIST::adp_node_list.erase( _Udc );
      }
    }
  protected: // override virtual
    char id_letter()const {return '\0';}
    bool print_type_in_spice()const {return false;}
    std::string value_name()const  {return "sensitivity";}
    //std::string dev_type()const;   //BASE_SUBCKT
    uint_t       max_nodes()const     {return 3;}
    uint_t       min_nodes()const     {return 2;}
    //int     matrix_nodes()const; //BASE_SUBCKT
    uint_t       net_nodes()const     {return _net_nodes;}
    uint_t       int_nodes()const     {return 0;}
    CARD*     clone()const         {return new DEV_BUILT_IN_RCD(*this);}
    void      precalc_first() {COMPONENT::precalc_first(); if(subckt()) subckt()->precalc_first();}
    virtual   void      expand(); // virtual??
    void      precalc_last();
    //void    map_nodes();         //BASE_SUBCKT
    void    tr_begin();    // BASE_SUBCKT
    void    tr_stress(); //BASE_SUBCKT
    void    tr_stress_();
    void    tr_restore();
    void do_tt();
    void tt_begin();
    void tt_prepare();
    void tt_advance();
    void tt_regress();
    void dc_advance() {set_not_converged(); BASE_SUBCKT::dc_advance();}
    void tr_advance();
    void tr_regress();
    bool tr_needs_eval()const;
    // ????? what is this good for?
    void      tr_queue_eval()      {if(tr_needs_eval()){q_eval();}}
    bool      do_tr();
    //void    tr_load();           //BASE_SUBCKT
    TIME_PAIR tr_review();
    void      tr_accept(); // 	{assert(subckt()); subckt()->tr_accept();}
    //void    tr_unload();         //BASE_SUBCKT
    double    tr_probe_num(const std::string&)const;
    double    tt_probe_num(const std::string&)const;
    //void    ac_begin();          //BASE_SUBCKT
    //void    do_ac();             //BASE_SUBCKT
    //void    ac_load();           //BASE_SUBCKT
    //XPROBE  ac_probe_ext(CS&)const;//CKT_BASE/nothing
    TIME_PAIR tt_review();
    double _uin_eff;
  public:
    double    involts()const;
    virtual double E()const {return 3;}
    virtual double P() const;
    static int  count() {return _count;}
  public: // may be used by models
  private: // not available even to models
    void     expand_net();
    void     expand_sym();
    static int _count;
  public: // input parameters
    region_t _region;	// fwd, reverse, unknown
  private: // calculated parameters
    int _tt_region;	// fwd, reverse, unknown
  public: // netlist
    COMPONENT* _Ccg;
    COMPONENT* _Ye;
    COMPONENT* _Re;
    COMPONENT* _Rc;
    COMPONENT* _GRc;
    ADP_NODE* _Ccgfill;
    ADP_NODE_UDC* _Udc;
    long double _tr_fill;
    long double _tr_dfill;
  protected: // node list
    enum {n_u, n_b, n_p};
    node_t _nodes[4];
    std::string port_name(uint_t i)const {
      assert(i != INVALID_NODE);
      assert(i < 3);
      static std::string names[] = {"u", "b", "p", ""};
      return names[i];
    }
    double _time0;
    double _time1;
  public:
    virtual region_t region() const ;
    virtual int tt_region() const ;
    virtual void tr_stress_last();
    mutable unsigned _iter_newton;
    mutable unsigned _iter_bisect;
  private:
    template<class T>
      T __uin_iter(T& uin, T E_old, T E_in, double bound_lo, double bound_hi) const;
    template<class T>
      T seff_bisect(T& uin, T E_old, T E_in, double bound_lo, double bound_hi) const;

    long double extrapolate(const double* tr0=0) const;

#ifndef NDEBUG
    unsigned _stressiter;
#endif
    double _ttgain;
    double _ttfuture;
};
/*--------------------------------------------------------------------------*/
#define T long double
template<>
T DEV_BUILT_IN_RCD::__uin_iter(T& uin, T E_old, T E_in, double bound_lo, double bound_hi) const;
#undef T
#define T long double
template<>
T DEV_BUILT_IN_RCD::seff_bisect<T>(T& uin, T E_old, T E_in, double bound_lo, double bound_hi) const;
#undef T
/*--------------------------------------------------------------------------*/
class ADP_BUILT_IN_RCD
  :public ADP_CARD{
public:
  explicit ADP_BUILT_IN_RCD( COMPONENT* c) : ADP_CARD(c)
    {init(c);}

  void init(const COMPONENT*);
public:
  ADP_NODE* ids_stress;
  ADP_NODE* igd_stress;
public:

  ADP_NODE* bti_stress; // FIXME, BTI_?
  double tr_probe_num(const std::string& x)const;
  double tt_probe_num(const std::string& x)const;

private:

public:
 // virtual void tt_accept();
  double vthscale_bti ; //  exp ( 10000. * a->hci_stress->get() / c->w_in );
  double vthdelta_bti ;
  double eff(){return 0.0; }
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
template<class T>
inline T MODEL_BUILT_IN_RCD::E_end(T s, const COMMON_COMPONENT* c ) const
{
  assert(is_number(s));
  T tau_e = __Re(s, c);
  assert(is_number(tau_e));
  T tau_c = __Rc(s, c);
  if (!is_number(tau_c)){
    error(bDANGER, "nan in tau_c at %f, %E\n", s, s);
  }
  assert(is_number(tau_c));
  T ret =  tau_e / ( tau_e + tau_c ); // swapped
  // T ret = T( 1.l / ( 1.l +  tau_e / tau_c ));
  assert(is_number(ret));
  return ret;
}
/*--------------------------------------------------------------------------*/
//template<class T>
#define T long double
inline T MODEL_BUILT_IN_RCD::E_end_0( const COMMON_COMPONENT* cc)const
{
  return E_end(0.L,cc);
}
#undef T
/*--------------------------------------------------------------------------*/
#endif

#if 0
double COMMON_BUILT_IN_RCD::__tau(double uin)const
{
  const COMMON_COMPONENT* cc=this;
  const MODEL_BUILT_IN_RCD* m = prechecked_cast<const MODEL_BUILT_IN_RCD*>(cc->model());
  assert(m);
  double tau_e = m->__Re( uin, cc );
  double tau_c = m->__Rc( uin, cc );
  if (tau_c < tau_e){
    return tau_c  / ( tau_c/tau_e + 1 );
  } else{
    return tau_e  / ( tau_e/tau_c + 1 );
  }
}
///*--------------------------------------------------------------------------*/
#endif
//* vim:ts=8:sw=2:et:
