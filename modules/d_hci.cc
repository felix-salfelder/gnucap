/*                     -*- C++ -*-
 * Copyright (C) 2013 Felix Salfelder
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
 * simple hci device
 */
#include "e_storag.h"
#include "d_mos.h"
#include "d_mos8.h" // hack
#include "e_adp.h"
namespace { //
using std::string;
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// FIXME: this is not a storage ELEMENT
class DEV_HCI : public ELEMENT { //
  enum {n_hci=0};
  double _hci_tr;
  double _L;
protected:
  explicit DEV_HCI(const DEV_HCI& p) :ELEMENT(p) {
    _hci_node.set_adp();
    _n = &_hci_node;
  }
public:
  explicit DEV_HCI()	:ELEMENT() {
    _hci_node.set_adp();
    _n = &_hci_node;
  }
 // void      keep_ic();
protected: // override virtual
  node_t _hci_node;
  std::string value_name()const {return "dvth";} // quick hack
                                                 // implement interface for adp::apply
  // std::string dev_type()const	{return "capacitor";}
  uint_t	   max_nodes()const	{return 0;}
  uint_t	   min_nodes()const	{return 0;}
  uint_t	   matrix_nodes()const	{return 0;}
  uint_t	   net_nodes()const	{return 0;}
  // bool	   has_iv_probe()const  {return true;}
  // bool	   use_obsolete_callback_parse()const {return true;}
  CARD*	   clone()const		{return new DEV_HCI(*this);}
  // void	   tr_iwant_matrix()	{tr_iwant_matrix_passive();}
  bool	   do_tr();
  void	   tr_accept(); // uic. possibly a hack
  void	   tr_stress_last();
  // void	   tr_load()		{tr_load_passive();}
  // void	   tr_unload()		{tr_unload_passive();}
  TIME_PAIR tr_review();
/*--------------------------------------------------------------------------*/
public:
  double    tr_involts()const	{ return 0;}
  void dc_advance();
  void tr_begin();
  void tr_advance();
  void tr_regress();
public: // tt
  void tt_begin();
  TIME_PAIR tt_review();
  void tt_accept();
  void do_tt();
  void tt_advance();
  void tt_regress();

  int order()const		{return min(ELEMENT::order(), 1);}
  double error_factor()const	{return err*OPT::trstepcoef[3];}
protected:
  hp_float_t   tr_involts_limited()const {return tr_outvolts_limited();}
  double   tr_probe_num(const std::string&)const;
  double   tt_probe_num(const std::string&)const;
  void	   tr_iwant_matrix(){}
  void	   ac_iwant_matrix(){}
//  void	   ac_begin()		{_ev = _y[0].f1;}
//  void	   do_ac();
//  void	   ac_load()		{ac_load_passive();}
  COMPLEX  ac_involts()const	{return 0;}

  string port_name(uint_t)const { untested();
    return "unknown";
  }
private:
  double err;
  double _ttfuture;
  FPOLY1   _i[OPT::_keep_time_steps];
  void precalc_first();
  void precalc_last();
  void expand();
  double vthdelta_hci;
  double tr_stress_() const;

  double _c1;
  double _c2;

  double _expected_tr;

};
/*--------------------------------------------------------------------------*/
TIME_PAIR DEV_HCI::tr_review()
{
  const DEV_BUILT_IN_MOS* d = prechecked_cast<const DEV_BUILT_IN_MOS*>(owner());
  assert(d);
  _time_by.reset();

  if(!_sim->_age) { // BUG
   // return _time_by;
  }

  _y[0].x = tr_stress_();

  _y[0].f0 = 5e-8*d->ids;
  _y[0].f0 = 1e-8*d->ids;

  err = 30;
  _c1 = tr_review_trunc_error(_y);

  if (_c1 < _sim->_dtmin) {
    error(bTRACE, "hci zero time step\n");
    error(bTRACE, "c1=%e  rejectedtime=%e  time1=%e\n",
        _c1, _sim->_time0, _time[1]);
    _c1 = _sim->_dtmin*1.000001;
  }

  if (_c1 < _dt * 0.75) { // .8 is hardcoded limit in s_tr_swp
    error(bTRACE, "step rejected:" + long_label() + '\n');
    error(bTRACE, "new=%g  old=%g  required=%g, ord=%d\n",
                  _c1, _dt, _dt * OPT::trreject, order());
    error(bTRACE, "%g %g\n",  _y[0].f0,  _y[1].f0);
    _c1 = _time[1] + _c1;
  }else{
    _c1 = _time[0] + _c1;
  }

  _time_by.min_error_estimate(_c1);
  _c1 -= _time[0];

  double Isub;
  if (d->reversed) {
    Isub = d->isb;
  } else {
    Isub = d->idb;
  }
  _i[0].f0 = Isub;

  err = 200;
  _c2 = tr_review_trunc_error(_i);

  if (_c2 < _sim->_dtmin) {
    error(bTRACE, "hci zero time step\n");
    error(bTRACE, "c2=%e  rejectedtime=%e  time1=%e\n",
        _c2, _sim->_time0, _time[1]);
    _c2 = _sim->_dtmin*1.000001;
    int error_deriv = order()+1;
    trace4("..", error_deriv, _i[0].f0, _i[1].f0, _i[2].f0);
  }

  if (_c2 < _dt * .79) {
    error(bTRACE, "step rejected2:" + long_label() + '\n');
    error(bTRACE, "new=%g  old=%g  required=%g\n",
                  _dt, _sim->_dt0, _dt * OPT::trreject);
    _c2 = _time[1] + _c2;
  }else{
    _c2 = _time[0] + _c2;
  }
  _time_by.min_error_estimate(_c2);
  _c2 -= _time[0];

  return _time_by;
}
/*--------------------------------------------------------------------------*/
void DEV_HCI::precalc_first()
{
  trace1("DEV_HCI::precalc_first", long_label());
//  ELEMENT::precalc_first();
}
/*--------------------------------------------------------------------------*/
void DEV_HCI::precalc_last()
{
  trace1("DEV_HCI::precalc_last", long_label());
//  ELEMENT::precalc_last();
}
/*--------------------------------------------------------------------------*/
void DEV_HCI::expand()
{
  if (!subckt()) {
    new_subckt();
  }else{ untested();
  }
  if (_sim->is_first_expand()) {
    // precalc_first();
    // precalc_last();
    trace4("DEV_HCI::expand, first", long_label(), hp(this), _n[n_hci].n_(), _n[n_hci].is_adp());
    assert(_n[n_hci].is_adp());
    assert(!(_n[n_hci].n_())); // n_ is electrical
    if (!(_n[n_hci].a_())){
      _n[n_hci].new_model_adp_node("raw_hci", this);
    } else { untested();
      unreachable();
    }
  }
}
/*--------------------------------------------------------------------------*/
bool DEV_HCI::do_tr()
{
  q_accept();
  return true;
}
/*--------------------------------------------------------------------------*/
void DEV_HCI::tt_begin()
{
  ADP_NODE* _raw_hci_node = _hci_node.a_();
  assert(_raw_hci_node);
  _i[0].f0 = 0;
  _i[0].f0 = 0;
  _y[1].f0 = 0;
  _y[1].f0 = 0;
  assert(_raw_hci_node);
  if (_sim->_tt_uic) {
    _raw_hci_node->set_tr(0);
  } else {
    _raw_hci_node->set_tr(0);
    _raw_hci_node->set_tt(0);
  }
  _hci_tr = 0;
  vthdelta_hci = 0;
}
/*--------------------------------------------------------------------------*/
void DEV_HCI::dc_advance()
{
  ELEMENT::dc_advance();

  for (int i = 1;  i < OPT::_keep_time_steps;  ++i) {
    _i[i] = _i[0];
  }
}
/*--------------------------------------------------------------------------*/
void DEV_HCI::tr_advance()
{
  ELEMENT::tr_advance();

  for (int i=OPT::_keep_time_steps-1; i>0; --i) {
    _i[i] = _i[i-1];
  }
}
/*--------------------------------------------------------------------------*/
void DEV_HCI::tr_regress()
{
  ELEMENT::tr_regress();
}
/*--------------------------------------------------------------------------*/
void DEV_HCI::tt_advance()
{ 
  ELEMENT::tt_advance();
  ADP_NODE* _raw_hci_node = _hci_node.a_();
  assert(_raw_hci_node);
  _hci_tr = 0;
  _raw_hci_node->set_tr(0);
}
/*--------------------------------------------------------------------------*/
void DEV_HCI::tt_regress()
{ untested();
  ELEMENT::tt_regress();
  ADP_NODE* _raw_hci_node = _hci_node.a_();
  assert(_raw_hci_node);
  _hci_tr = 0;
  _raw_hci_node->set_tr(0);
}
/*--------------------------------------------------------------------------*/
void DEV_HCI::tr_begin()
{
  ELEMENT::tr_begin();
  if(_sim->_cont){ // FIXME. move to tr_accept
    for (int i = 1;  i < OPT::_keep_time_steps;  ++i) {
      _i[i] = FPOLY1(0., 0., 0.); // state
    }
    trace7("STORAGE::tr_begin", _i[0].x, tr_involts(), _i[0].f0, _m0.x, _i[0].f1, _m0.c0, _m0.c1);
     // m.x = volts, m.c0 = amps, acg = m.c1 = mhos
  } else {
    for (int i = 0;  i < OPT::_keep_time_steps;  ++i) {
      _i[i] = FPOLY1(0., 0., 0.); // state
    }
  }
}
/*--------------------------------------------------------------------------*/
void DEV_HCI::do_tt()
{
  const COMMON_BUILT_IN_MOS* c = prechecked_cast<const COMMON_BUILT_IN_MOS*>(common());
  assert(c);
  const MODEL_BUILT_IN_MOS8* m = prechecked_cast<const MODEL_BUILT_IN_MOS8*>(c->model());
  assert(m);

  ADP_NODE* a = _hci_node.a_();
  assert(a);

  double eff_now = a->tr( _sim->_Time0 ); // fixme: faster?
  if(eff_now<0){ unreachable();
    trace2("do_tt", long_label(), eff_now);
    eff_now = 0;
  }
  double eff_last_timeframe = a->tr(_sim->last_Time());
  if(eff_last_timeframe<0){ unreachable();
    eff_last_timeframe = 0;
  }

  double ex_time = _sim->_dT0 - _sim->last_time(); // stress that long
  if(fabs(ex_time)<=1e-18){
    ex_time=0;
  }else{
  }
  if(_sim->phase() == p_PD){
    ex_time = 0;
  }
  trace6("ADP_BUILT_IN_MOS8::do_tt", eff_last_timeframe, eff_now, a->tt1(), a->order(), a->tr1(), a->tr());
  trace6("ADP_BUILT_IN_MOS8::do_tt", a->tr(), _sim->last_Time(), ex_time, _sim->_dT0, _sim->last_time(), _sim->_Time0);

  if (ex_time<0) {
    error(bDANGER, "extime %f?! %f last_time %f last_Time %f\n",ex_time, _sim->_dT0, _sim->last_time(), _sim->last_Time());
    unreachable(); // assert(false);
    ex_time = 0.;
  }

  double hci_new;
  if (!a->order()) {
    hci_new = a->tt();
    trace2("ADP_BUILT_IN_MOS8::do_tt", a->tt(), a->tt1());
  } else {
    trace2("ADP_BUILT_IN_MOS8::do_tt", _sim->last_time(), _sim->_dTmin);
    hci_new = a->tt1(); // tt @ last_Time (?)
  }

  if (!is_number(hci_new) && _sim->phase() != p_PD){
    error(bDANGER, "ADP_BUILT_IN_MOS8::do_tt hci hist bug Time0 %E last %E ordr %i \n", _sim->_Time0, _sim->_last_Time, a->order() );
  }

  if(_sim->phase() == p_PD) {
    //hack ?
    hci_new = a->tt1();
    assert(is_number(hci_new));
  } else {
    
    switch (a->order()) {
      case 4:
      case 3:
      case 2:
	// quadratic extrapolation
	assert(is_number(eff_last_timeframe));
	hci_new = a->tt_integrate_2_poly(_sim->last_time());
	trace6("pre integrate", _sim->_Time0, hci_new, a->tt(), a->tt1(), a->tr1(), a->tr());
	_expected_tr = eff_now;
	break;
      case 1:
	assert(is_number(eff_last_timeframe));
	hci_new += ex_time * (  eff_last_timeframe + eff_now ) / 2.0 ;
	_expected_tr = eff_now; // hack??
	break;
      case 0:
	assert(is_number(hci_new));
	hci_new += ex_time * eff_now;
	break;
      default:
	assert(0);
    }
  }
  assert(is_number(hci_new));

  if (!is_number(hci_new)) { untested();
    error(bDANGER, "mos8 do_tt Time0 %E last %E ordr %i, %f, "
	"ex_time %f eff_now %f \n", _sim->_Time0, _sim->_last_Time,
	a->order(), a->tt1(), ex_time, eff_now
	);
    assert(is_number(hci_new));
  }

  a->set_tt(hci_new);
  a->set_tr(_expected_tr); // need default in case we don't intend to recompute
  assert(is_number(a->tt()));

  double H = m->h0;
  const SDP_BUILT_IN_MOS_BASE* s = prechecked_cast<const SDP_BUILT_IN_MOS_BASE*>(c->sdp());
  double W = s->w_eff;

  vthdelta_hci = pow(a->tt()/(H*W), m->hci_n);
  trace4("", a->tt(), H, W, m->hci_n);
  assert(is_number(vthdelta_hci));
  set_value(vthdelta_hci);
  assert(!_hci_tr);
}
/*--------------------------------------------------------------------------*/
void DEV_HCI::tr_stress_last()
{
  ADP_NODE* _raw_hci_node = _hci_node.a_();
  assert(_raw_hci_node);
  _raw_hci_node->set_tr(_hci_tr/_sim->last_time());
  q_tt_accept();
  _raw_hci_node->set_tr_noise(0); // incomplete;
}
/*--------------------------------------------------------------------------*/
double DEV_HCI::tr_stress_() const
{
  double exponent = 3; // which is m in [4]
  double hcis = 0;
  double Wg = 0.8;
  const DEV_BUILT_IN_MOS* d = prechecked_cast<const DEV_BUILT_IN_MOS*>(owner());
  const COMMON_BUILT_IN_MOS* c = asserted_cast<const COMMON_BUILT_IN_MOS*>(d->common());
  const MODEL_BUILT_IN_MOS8* m = asserted_cast<const MODEL_BUILT_IN_MOS8*>(c->model());
  const SDP_BUILT_IN_MOS_BASE* s = prechecked_cast<const SDP_BUILT_IN_MOS_BASE*>(c->sdp());
  assert(d);
  double Ids = fabs(d->ids); // fabs?

  //  std::cerr << "DEV_BUILT_IN_MOS8::h0 of "<<  d->short_label() << " " <<   m->h0 << "\n";
  double H = m->h0;
  double W = s->w_eff;
  double Hg = m->h0;
  //  double m0 = m->m0;

  if( Ids < 1e-40) {
    trace1("MODEL_BUILT_IN_MOS8::tr_stress ids too small: ", d->ids );
  } else {
    double Isub;
    if (d->reversed) {
      Isub = d->isb;
    } else {
      Isub = d->idb;
    }

    assert(Isub >= 0);
    //  assert(d->ids >= 0); isnich
    //  assert( dt >= 0 )

    switch (m->polarity) { untested();
      case dunno:
	unreachable();
      case pN:
	hcis = Ids * pow( Isub / Ids, m->hci_m);
	assert(is_number(hcis));
	break;
      case pP:
	double mg = 3.0;
	double ig = d->probe_num("ig"); // BUG
	hcis = ( Wg/Hg * pow( fabs(ig)/W, mg ) * H * W
	    + (1-Wg)*Ids * pow(Isub/fabs(Ids), exponent));
	trace8( "ADP_BUILT_IN_MOS8::tr_accept pmos", d->long_label(),
	    d->isb, d->idb, _sim->_time0, Isub, Ids, H, hcis);
	assert(is_number(hcis));
	break;
    }
    if (hcis > 1e-10)
    {
    }
  }
  return hcis;
}
/*--------------------------------------------------------------------------*/
void DEV_HCI::tr_accept()
{
  const DEV_BUILT_IN_MOS* d = prechecked_cast<const DEV_BUILT_IN_MOS*>(owner());
  assert(d);
  const COMMON_BUILT_IN_MOS* c = asserted_cast<const COMMON_BUILT_IN_MOS*>(d->common());
  const MODEL_BUILT_IN_MOS8* m = asserted_cast<const MODEL_BUILT_IN_MOS8*>(c->model());
  const SDP_BUILT_IN_MOS_BASE* s = prechecked_cast<const SDP_BUILT_IN_MOS_BASE*>(c->sdp());
  assert(d);
  double H = m->h0;
  double W = s->w_eff;

  { // must not accept twice!
    // effectively add up on tt, and take difference in the end (numerically unstable...)
    double hcis = tr_stress_();
    _y[0].x = hcis;

    if (_sim->_dt0) {
      _hci_tr += (_y[1].x + _y[0].x) * _sim->_dt0;
    }else{
    }

    ADP_NODE* _raw_hci_node = _n[n_hci].a_(); assert(_raw_hci_node);
    vthdelta_hci = pow((_raw_hci_node->tt()+_hci_tr) / (H*W) ,m->hci_n);
    //delta_vth=vthdelta_hci;
  } // end hci
  q_eval();
}
/*--------------------------------------------------------------------------*/
TIME_PAIR DEV_HCI::tt_review()
{
  ADP_NODE* a = _hci_node.a_();
  assert(a);
  double c[OPT::_keep_time_steps];
  double timestep = NEVER;

  a->derivatives(c);
  double tol = OPT::tttol;

  double delta = fabs(_expected_tr - a->tr());
  trace7("DEV_HCI::tt_review", a->order(), c[0], c[1], c[2], a->tr(), _expected_tr, delta);

  switch(a->order()){
    case 0: incomplete();
      break;
    case 1:
      timestep = _sim->last_time();
      assert(timestep);
      break;
    case 2:
      timestep = tol / fabs(c[2]) * 1e-4;
      break;
    case 3: untested();
      if(OPT::ttcorr) {
	assert(a->tr() == _hci_tr/_sim->last_time());
	trace2("hci::review corr", a->tt(), _sim->_dT0 * (a->tr() - _expected_tr));
	double fac = .2;
	a->tt() += _sim->_dT0 * (a->tr() - _expected_tr) * fac;
	fac = .0;
	a->tr() = ((1.-fac) * a->tr() + fac * _expected_tr);
      }
      timestep = 1e4 * sqrt( tol / delta );
      break;
    default: incomplete();
  }

  _ttfuture = tt_review_check_and_convert(timestep);

  return TIME_PAIR(_ttfuture,NEVER);
}
/*--------------------------------------------------------------------------*/
void DEV_HCI::tt_accept()
{
  DEV_BUILT_IN_MOS* d = asserted_cast<DEV_BUILT_IN_MOS*>(owner());
  const COMMON_BUILT_IN_MOS* c = asserted_cast<const COMMON_BUILT_IN_MOS*>(d->common());
  const MODEL_BUILT_IN_MOS8* m = asserted_cast<const MODEL_BUILT_IN_MOS8*>(c->model());
  ELEMENT::tt_accept();
  ADP_NODE* _raw_hci_node = _hci_node.a_();
  assert(_raw_hci_node);
  const SDP_BUILT_IN_MOS_BASE* s = prechecked_cast<const SDP_BUILT_IN_MOS_BASE*>(c->sdp());
  assert(d);
  double H = m->h0;
  double W = s->w_eff;

  _raw_hci_node->tt() += _hci_tr;
  if (m->hci_n!=.3) {untested();}
  vthdelta_hci = pow(_raw_hci_node->tt()/(H*W), m->hci_n);
  assert(is_number(vthdelta_hci));
  assert(fabs(vthdelta_hci)<10);
  set_value(vthdelta_hci);

  _L = _hci_tr / _sim->last_time();
  _hci_tr = 0.;
}
/*--------------------------------------------------------------------------*/
double DEV_HCI::tr_probe_num(const std::string& x)const
{
  const DEV_BUILT_IN_MOS* d = prechecked_cast<const DEV_BUILT_IN_MOS*>(owner());
  assert(d);
  const COMMON_BUILT_IN_MOS* c = asserted_cast<const COMMON_BUILT_IN_MOS*>(d->common());
  const MODEL_BUILT_IN_MOS8* m = asserted_cast<const MODEL_BUILT_IN_MOS8*>(c->model());
  const SDP_BUILT_IN_MOS_BASE* s = prechecked_cast<const SDP_BUILT_IN_MOS_BASE*>(c->sdp());
  ADP_NODE* a = _n[n_hci].a_(); assert(a);
  double H = m->h0; USE(H);
  double W = s->w_eff; USE(W);
  const ADP_NODE* _raw_hci_node = _hci_node.a_();
  assert(_raw_hci_node);
  if (Umatch(x, "hci")) {
    return _hci_tr + _raw_hci_node->tt();
  } else if (Umatch(x, "dhci")) {
    return _hci_tr;
  } else if (Umatch(x, "dvth")) {
  //  double vthdelta_hci = pow(_raw_hci_node->tt()/(H*W), m->hci_n);
    return vthdelta_hci;
  } else if (Umatch(x, "stress{level}")) {
    return _y[0].x; // (H*W);
  } else if (Umatch(x, "tr ")) { untested();
    return _hci_tr;
  } else if (Umatch(x, "tt ")) {
    return a->tt();
  } else if (Umatch(x, "t0 ")) { untested();
    return _time_by._error_estimate;
  } else if (Umatch(x, "t1 ")) { untested();
    return _time[1];
  } else if (Umatch(x, "m ")) { itested();
    return m->hci_m;
  } else if (Umatch(x, "n ")) { itested();
    return m->hci_n;
  } else if (Umatch(x, "c1 ")) { itested();
    return _c1;
  } else if (Umatch(x, "c2 ")) { itested();
    return _c2;
  } else if (Umatch(x, "ttord{er}")) {
    return _raw_hci_node->order();
  }
  return ELEMENT::tr_probe_num(x);
}
/*--------------------------------------------------------------------------*/
double DEV_HCI::tt_probe_num(const std::string& x)const
{
  if (Umatch(x, "ttf{uture}")) { untested();
    return _ttfuture;
  } else if (Umatch(x, "stress{level} ")) {
    return _L;
  }
  return ELEMENT::tt_probe_num(x);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
DEV_HCI p1;
DISPATCHER<CARD>::INSTALL
  d1(&device_dispatcher, "hci",	    &p1);
}
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet
