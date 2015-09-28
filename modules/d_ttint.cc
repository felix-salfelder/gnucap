/*                                 -*- C++ -*-
 *
 * tt int device ...
 *
 * (c) 2014 Felix Salfelder
 *
 * GPLv3
 */

#include "e_aux.h"
#include "e_adp.h"
#include "e_storag.h"
#include "globals.h"
#include "e_elemnt.h"
#include "u_nodemap.h"
#include <iomanip>
#ifdef DO_TRACE
# include "io_misc.h"
#endif

/*--------------------------------------------------------------------------*/
class ADP_NODE;
/*--------------------------------------------------------------------------*/
namespace { //
using namespace std;
/*--------------------------------------------------------------------------*/
const double _default_value (1); //input scale.
/*--------------------------------------------------------------------------*/
const double INF(BIGBIG);
/*--------------------------------------------------------------------------*/
class DEV_TTINT : public ELEMENT { //
  private:
    explicit DEV_TTINT(const DEV_TTINT& p);
    double _lasts; // last time tr_stress was called
  public:
    explicit DEV_TTINT();
    ~DEV_TTINT() { itested();
      --_count;
    }
  protected: // override virtual
    void tr_iwant_matrix(){}
    void ac_iwant_matrix(){}
    double  tr_involts()const{return 0;}
    double    involts()const {return 0;}
    double  tr_involts_limited()const{return 0;}
    COMPLEX ac_involts()const{return 0;}
    char id_letter()const {return '\0';}
    bool print_type_in_spice()const {return false;}
    std::string value_name()const  {return "sensitivity";}
    uint_t       max_nodes()const     {return 1;}
    uint_t       min_nodes()const     {return 0;}
    //int     matrix_nodes()const; //BASE_SUBCKT
    uint_t       net_nodes()const     {return _net_nodes;}
    uint_t       int_nodes()const     {return 0;}
    CARD*     clone()const         {return new DEV_TTINT(*this);}
    void      precalc_first() {COMPONENT::precalc_first(); if(subckt()) subckt()->precalc_first();}
    virtual   void      expand(); // virtual??
    void      precalc_last();
    //void    map_nodes();         //
    void tr_begin();
    //void    tr_restore();
    void do_tt();
    void tt_begin();
    void tt_advance();
    void tt_regress();
    bool do_tr();
    TIME_PAIR tr_review();
    void      tr_accept(); // 	{assert(subckt()); subckt()->tr_accept();}
    double    tr_probe_num(const std::string&)const;
    TIME_PAIR tt_review();
  public:
  private: // not available even to models
    static int _count;
  protected: // node list
    enum {n_p};
    node_t _nodes[1];
    std::string port_name(uint_t i)const { untested();
      assert(i != INVALID_NODE);
      assert(i < 1);
      static string names[] = {"p"};
      return names[i];
    }
  public:
    void tr_stress_last();
  private:
    long double extrapolate(const double* tr0=0) const;
    double _tt0;
};
/*--------------------------------------------------------------------------*/
int DEV_TTINT::_count = -1;
/*--------------------------------------------------------------------------*/
void DEV_TTINT::tr_accept()
{
  double h = _sim->_dt0;
  _n[n_p].a_()->tt() += _n[n_p].a_()->tr() * h;
}
/*--------------------------------------------------------------------------*/
namespace DEV_TTINT_DISPATCHER { //
  static DEV_TTINT p0;
  static DISPATCHER<CARD>::INSTALL
    d0(&device_dispatcher, "ttint", &p0);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
DEV_TTINT::DEV_TTINT()
  :ELEMENT(),
   _lasts(-inf)
{ itested();
  _n = _nodes;
  _nodes[n_p].set_adp();
  assert(_n[n_p].is_adp());
  // _nodes[n_p].a_()->trhack = 0;

  ++_count;
  // overrides
  set_value(1);
}
/*--------------------------------------------------------------------------*/
DEV_TTINT::DEV_TTINT(const DEV_TTINT& p)
  :ELEMENT(p),
   _lasts(-inf)
{ itested();
  _n = _nodes;
  assert(p._n[n_p].is_adp());
  for (uint_t ii = 0; ii < max_nodes() + int_nodes(); ++ii) { itested();
    _n[ii] = p._n[ii];
    trace2("DEV_TTINT::DEV_TTINT(p)", ii, _n[ii].is_adp());
  }
  assert(_n[n_p].is_adp());
  ++_count;
  // overrides
}
/*--------------------------------------------------------------------------*/
void DEV_TTINT::expand()
{ itested();
  trace1("DEV_TTINT::expand", long_label());
  assert(_n[n_p].is_adp());
  assert(_n);

  if (_sim->is_first_expand()) { itested();
    precalc_first();
    precalc_last();
    assert(_n[n_p].is_adp());
    assert(!(_n[n_p].n_())); // n_ is electrical
    if (!(_n[n_p].a_())){ untested();
      _n[n_p].new_model_adp_node("c", this);
    } else { itested();
    }
  }
}
/*--------------------------------------------------------------------------*/
double DEV_TTINT::tr_probe_num(const std::string& x)const
{ itested();
  ADP_NODE* Ccgfill = _n[n_p].a_();
  assert(Ccgfill);

  if (Umatch(x, "tr ")) { itested();
    return Ccgfill->tr();
  }else if (Umatch(x, "tt ")) { itested();
    return Ccgfill->tt();
  }else if (Umatch(x, "s{tresslevel} ")) { untested();
    return Ccgfill->tr();
  }else if (Umatch(x, "uhi |vhi ")) { untested();
    return  ( Ccgfill->tr_hi );
  }

  return ELEMENT::tr_probe_num(x);
}
/*--------------------------------------------------------------------------*/
long double DEV_TTINT::extrapolate(const double* tr0)const
{ itested();
  trace2("DEV_TTINT::extrapolate", hp(tr0), _sim->_Time0);
  if(tr0){ untested();
    assert(is_number(*tr0));
  }else{ itested();
  }

  double E_old = _n[n_p].a_()->tt1();
  assert(is_number(E_old));

  double ex_time = _sim->_dT0 - _sim->last_time();
  assert(ex_time >= 0);

  double Time1 = _sim->_Time0 - _sim->_dT0;
  double eff1 = value()*_n[n_p].a_()->tr( Time1 + ex_time/3.0, tr0);
  double eff2 = value()*_n[n_p].a_()->tr( Time1 + ex_time*2.0/3.0, tr0);

  if(!is_number(eff1)){ untested();
    error(bDANGER, "nan bug %d\n", _n[n_p].a_()->order());
    assert(0);
  }
  assert(is_number(eff2));

  double fill_new = E_old + (eff1 + eff2 )*.5*ex_time;
  assert(is_number(fill_new));

  return fill_new;
}
/*--------------------------------------------------------------------------*/
bool DEV_TTINT::do_tr()
{ itested();
  return true;
}
/*--------------------------------------------------------------------------*/
void DEV_TTINT::do_tt()
{ itested();
  assert(_n);

  assert(_sim->_time0 == 0 || _sim->_mode==s_TRAN ); //?

  if(_sim->phase() == p_PD){ untested();
    // assert(!_n[n_p].a_()->tr());
    _n[n_p].a_()->tr() = 0; // check: what is tr()?
                        // why not use tr(0) and leave tr alone?
    assert(_sim->_dT0);
  }

  assert(_nodes[n_p].a_()->trhack == 0);

  double Time1 = _sim->_Time0 - _sim->_dT0;

  if (! ( is_almost( _n[n_p].a_()->tr1() , _n[n_p].a_()->tr(Time1) ))) { untested();
  }

  if (!is_number(_n[n_p].a_()->tr_rel(_sim->_dT0 ))) { unreachable();
  }

  double fill_new = 0;
  if (_sim->_dT0==0.) { itested();
    return;
  }else{ itested();
    fill_new = (double) extrapolate();
  }

  assert(is_number(fill_new));
  _n[n_p].a_()->tt() = (double) fill_new;
  _tt0 = fill_new;
  // print something useful for "stress" (good idea?)
  _n[n_p].a_()->tr() = _n[n_p].a_()->tr(_sim->_Time0);
}
/*----------------------------------------------------------------------------*/
void DEV_TTINT::tr_stress_last()
{ itested();
  ADP_NODE* a = _n[n_p].a_();
  assert(_n[n_p].tt() == _n[n_p].tt());
  assert(_n[n_p].tt() == a->tt());
  assert(a);

  if( _sim->phase()==p_PD){ untested();
        assert(!a->tr_lo);
        assert(!a->tr_hi);
  }

  trace5("DEV_BUILT_IN_RCD::tr_stress_last s ", _n[n_p].a_()->tt(),
      _n[n_p].a_()->get_tr(), a->get_total(), _n[n_p].m_(), _n[n_p].t_() );
  // fixme. move to common.
  assert(is_number(a->tt()));

  a->tr() = (a->tt()-_tt0) / _sim->last_time();

  if( _sim->phase()==p_PD ){ untested();
  }
}
/*--------------------------------------------------------------------------*/
void DEV_TTINT::tt_begin()
{ itested();
  ADP_NODE* a = _n[n_p].a_();
  assert(a); USE(a);

  if (_sim->_tt_uic) { untested();
  } else { itested();
    _n[n_p].tt() = 0;
  }
  assert(is_number(a->tt()));
  _n[n_p].a_()->tr() = 0;
  _tt0 = a->tt();

  _n[n_p].a_()->tr_lo = inf;
  _n[n_p].a_()->tr_hi = -inf;

  _lasts = -inf;
}
/*--------------------------------------------------------------------------*/
void DEV_TTINT::tr_begin()
{ itested();
  ELEMENT::tr_begin();

  _lasts = -inf; // must be -inf, so time0==0 works out
}
/*--------------------------------------------------------------------------*/
void DEV_TTINT::precalc_last()
{ itested();
  trace1("DEV_BUILT_IN_RCD::precalc_last", value());
  COMPONENT::precalc_last();
  if (!double(value())) { itested();
    set_value(1);
  }
}
/*--------------------------------------------------------------------------*/
TIME_PAIR DEV_TTINT::tt_review()
{ incomplete();
  return TIME_PAIR(NEVER,NEVER);
}
/*--------------------------------------------------------------------------*/
TIME_PAIR  DEV_TTINT::tr_review()
{ itested();
  q_accept();
  return TIME_PAIR();
}
/*--------------------------------------------------------------------------*/
void DEV_TTINT::tt_regress()
{ untested();
  _lasts = -inf;
  _n[n_p].a_()->tr_hi = -inf;
  _n[n_p].a_()->tr_lo = inf;
}
/*--------------------------------------------------------------------------*/
void DEV_TTINT::tt_advance()
{ itested();
  trace1("DEV_TTINT::tt_advance", _n[n_p].tt());

  _lasts = -inf;
  _n[n_p].a_()->tr_hi = -inf;
  _n[n_p].a_()->tr_lo = inf;
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
}
// vim:ts=8:sw=2:noet
