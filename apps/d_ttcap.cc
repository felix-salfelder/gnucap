/*                              -*- C++ -*-
 * Copyright (C) 2014 Felix Salfelder
 * Author: <felix@salfelder.org>
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
 * tt capacitance device. experimental
 */
#include "e_storag.h"
#include "m_divdiff.h"
#include "config.h"
#ifdef HAVE_GSL_FIT_H
# include <gsl/gsl_multifit.h>
#endif
namespace { //
#define TTCAP_HACK
#include "d_cap.h"
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
using SOME_CAP_HACK::DEV_CAPACITANCE;

bool DEV_CAPACITANCE::has_ic() const
{
  if ( const EVAL_BM_ACTION_BASE* x = dynamic_cast<const EVAL_BM_ACTION_BASE*>(common()) ){
    if (x->_ic != NOT_INPUT) return 1;
  }
  return 0;
}
/*--------------------------------------------------------------------------*/
bool DEV_CAPACITANCE::do_tr()
{
  // FPOLY1* q=_y;
  //trace0(("DEV_CAPACITANCE::do_tr " + long_label()));
  //trace3(("DEV_CAPACITANCE::do_tr " + long_label()).c_str(), _y[0].f1, value(), tr_input() );
  if (using_tr_eval()) {
    _y[0].x = tr_input_limited();
    tr_eval();
  }else{ itested();
    _y[0].x = tr_input(); // tr_involts();
    assert(_y[0].f1 == value());
    _y[0].f0 = _y[0].x * _y[0].f1;
    assert(converged());
  }
  store_values();
  q_load();

  if (_sim->uic_now() && has_ic()) {
    double G = 1./OPT::shortckt;
    assert(_time[0] == 0.);

    // imitate voltage source... (d_vs.cc)
    _i[0] = FPOLY1( CPOLY1( 0., -_y[0].x * G,         G  ) ); 
//    trace2("2 quotienten", ( -_y[0].f0 / (OPT::shortckt)  ) /-_y[0].x  / (OPT::shortckt) ,
//			         (_y[0].f1 / (OPT::shortckt))/    1/(OPT::shortckt)  );

    // picked up by tr_load_*
    _loss0 = G;

  } else {
    _i[0] = differentiate(_y, _i, _time, _method_a);
  }
  _m0 = CPOLY1(_i[0]);

  if (_sim->analysis_is_restore()) {
//    _ttstate[1] = _y[0].f0;
  }
  return converged();
}
/*--------------------------------------------------------------------------*/
void DEV_CAPACITANCE::tr_restore()
{
//  _ttstate[1] = _y[0].f0;
  q_accept();
}
/*--------------------------------------------------------------------------*/
void DEV_CAPACITANCE::tr_accept()
{
  ELEMENT::tr_accept(); // commmon->accept..

  if (_loss0) {
    // _loss0 should be irrelevant
    // sources might have been abused to enforce initial conditions
    // during dc analysis (otherwise this should be unreachable)

    if (0 && !_sim->uic_now()) { // s_ddc
      _m0.c0 = 0;
    } else {
      // vera hack. something about qdot
      // also used in dc uic + tr cont ...
      FPOLY1 m(_m0);

      m.x = tr_input();
      m.f0 = - _loss0 * tr_outvolts() - _m0.c0;
      m.f1 = 0;
      _m0 = CPOLY1(m);
      _i[0].f0=0;
      _i[0].f1=0;
    }

    _loss0 = 0.;
    set_converged(false);
  } else {
  }

  if (_sim->is_advance_iteration()) {
  }
}
/*--------------------------------------------------------------------------*/
void DEV_CAPACITANCE::tr_stress_last()
{
  q_tt_accept();
}
/*--------------------------------------------------------------------------*/
void DEV_CAPACITANCE::tt_regress()
{
//  _ttsteporder = _new_ttsteporder;

  // for now
//  _ttmethod = _ttsteporder;
  _ttstate[0] = NAN;
  _dttstate[0] = NAN;

  _b[0] = _vy[0] = _vt[0] = NAN;
}
/*--------------------------------------------------------------------------*/
// necessary? move to adv.
void DEV_CAPACITANCE::tt_accept()
{
  ELEMENT::tt_accept();
  assert(_y[0].f1>=0);

  _dT2 = _dT1;
  _dT1 = _sim->_dT0;

// set in review.
//  _ttstate_ = _ttstate[0] = _y[0].f0; // extrapolate from here...
}
/*--------------------------------------------------------------------------*/
void DEV_CAPACITANCE::tt_advance()
{
  _ttsteporder = _new_ttsteporder;

  for (unsigned i=3; i>0; --i) {
    _dttstate[i] = _dttstate[i-1];
    _ttstate[i] = _ttstate[i-1];
  }

  _b[1] = _b[0];
  _vy[1] = _vy[0];
  _vt[1] = _vt[0];

  _b[0] = _vy[0] = _vt[0] = NAN;

  // for now
  _ttmethod = _ttsteporder;
}
/*--------------------------------------------------------------------------*/
void DEV_CAPACITANCE::tt_begin()
{
  _ttsteporder = 0;
  _new_ttsteporder = 0;
  _ttstate[0] = _y[0].f0;

  _ttmethod = 0; // hmm advance bug?
  _dT1 = 0;
  _dT2 = 0;

  _ttstate[1] = NAN;
  _ttstate[2] = NAN;
  _ttstate[3] = NAN;
  _ttstate_ = NAN;
  _dttstate[0] = NAN;
  _dttstate[1] = NAN;
  _dttstate[2] = NAN;
  _dttstate[3] = NAN;
  _dttstate_ = NAN;
  _pred=NAN;
  _corr=NAN;
  _preddt=NAN;
  _ttstep=NAN;
  _ttstep0=NAN;
  _ttstep1=NAN;
  _ttstep2=NAN;
  _ttstep3=NAN;
  _ttfuture=NAN;

  _vt[0] = _vy[0] = _b[0] = 0;
  _vt[1] = _vy[1] = _b[1] = NAN;
  _chisq = 0;
  _dv = 0;
}
/*--------------------------------------------------------------------------*/
#ifdef HAVE_GSL_FIT_H
// find _vy[0], _vt[0], _b[0]
double DEV_CAPACITANCE::gsl_fit_3()
{
  static gsl_multifit_linear_workspace* work;
  static gsl_matrix *X, *cov;
  static gsl_vector *c, *y, *r;
  if (!work) {
    work = gsl_multifit_linear_alloc (3, 3);
    cov = gsl_matrix_alloc(3,3);
    X = gsl_matrix_alloc(3,3);
    gsl_matrix_set(X,0,0,0);
    gsl_matrix_set(X,0,2,1);
    gsl_matrix_set(X,1,2,1);
    gsl_matrix_set(X,2,2,1);

    c = gsl_vector_alloc(3);
    y = gsl_vector_alloc(3);
    r = gsl_vector_alloc(3);
  }
  double time[3];
  time[0] = 0;
  time[1] = - _sim->_dT0;
  time[2] = - _sim->_dT0 - _dT1;

  double tscale = 1e8;

  gsl_matrix_set(X,1,0,time[1]*tscale);
  gsl_matrix_set(X,2,0,time[2]*tscale);

  gsl_matrix_set(X,0,1,_ttstate[0]);
  gsl_matrix_set(X,1,1,_ttstate[1]);
  gsl_matrix_set(X,2,1,_ttstate[2]);

  gsl_vector_set(y,0,_dttstate[0]);
  gsl_vector_set(y,1,_dttstate[1]);
  gsl_vector_set(y,2,_dttstate[2]);

  int fit;
  double m0, m1, m2;
  double d1;

  fit = gsl_multifit_linear (X, y, c, cov, &_chisq, work);
  USE(fit);
  assert(!fit);
  _vt[0] = gsl_vector_get(c,0) * tscale;
  _vy[0] = gsl_vector_get(c,1);
  _b[0] = gsl_vector_get(c,2);
  d1 = time[1] * _vt[0] + _ttstate[1] * _vy[0] + _b[0];
  _dv = - _vt[0]/_vy[0];
  m0 = - _b[0] / _vy[0];
  m1 = ( - _b[0] - time[1]*_vt[0] ) / _vy[0];
  m2 = ( - _b[0] - time[2]*_vt[0] ) / _vy[0];

  trace3("corr", time[0], time[1], time[2]);
  trace3("corr", _dttstate[0], _dttstate[1], _dttstate[2]);
  trace3("corr", _ttstate[0], _ttstate[1], _ttstate[2]);
  trace3("fit", _vt, _vy, _b);
  trace3("corr", m0, m1, m2);
  trace4("corr", d1, _dv, _vt, _vy);
  trace1("tau", _dv-d1);

  gsl_multifit_linear_residuals (X, y, c, r);
  trace3("res", gsl_vector_get(r,0), gsl_vector_get(r,1), gsl_vector_get(r,2));

  double e = exp(-_sim->_dT0*fabs((_dv-d1)/(_ttstate[1]-m1)));
  double pred = (_ttstate[1] - m1) * e + m0;
  trace4("ex", e, pred, m0, m1);
  if(0){
    gsl_multifit_linear_free(work);
    gsl_matrix_free(cov);
    gsl_matrix_free(X);

    gsl_vector_free(c);
    gsl_vector_free(y);
    gsl_vector_free(r);
  }

  assert(is_number(_vt[0]));
  assert(is_number(_b[0]));
  assert(is_number(_vy[0]));
  return pred;
}
/*--------------------------------------------------------------------------*/
// find _vy[0], _vt[0], _b[0]
double DEV_CAPACITANCE::gsl_fit_4()
{
  static gsl_multifit_linear_workspace* work;
  static gsl_matrix *X, *cov;
  static gsl_vector *c, *y, *r;
  if (!work){ // memory leak!
    work = gsl_multifit_linear_alloc (4, 3);
    cov = gsl_matrix_alloc(3,3);
    X = gsl_matrix_alloc(4,3);
    gsl_matrix_set(X,0,0,0);

    gsl_matrix_set(X,0,2,1);
    gsl_matrix_set(X,1,2,1);
    gsl_matrix_set(X,2,2,1);
    gsl_matrix_set(X,3,2,1);

    c = gsl_vector_alloc(3);
    y = gsl_vector_alloc(4);
    r = gsl_vector_alloc(4);
    // incomplete(); yes.
  }
  double time[4];
  time[0] = 0;
  time[1] = - _sim->_dT0;
  time[2] = - _sim->_dT0 - _dT1;
  time[3] = - _sim->_dT0 - _dT2 - _dT1;

  double tscale=1e8;

  for (unsigned i=1; i<4; ++i) {
    gsl_matrix_set(X,i,0,time[i]*tscale);
  }

  for (unsigned i=0; i<4; ++i) {
    gsl_matrix_set(X,i,1,_ttstate[i]);
    gsl_vector_set(y,i,_dttstate[i]);
  }

  int fit;
  double m0, m1, m2;
  double d1;

  fit = gsl_multifit_linear (X, y, c, cov, &_chisq, work);
  USE(fit);
  assert(!fit);
  _vt[0] = gsl_vector_get(c,0) * tscale;
  _vy[0] = gsl_vector_get(c,1);
  _b[0] = gsl_vector_get(c,2);
  d1 = time[1] * _vt[0] + _ttstate[1] * _vy[0] + _b[0];
  _dv = - _vt[0]/_vy[0];
//  m0 = - _b / _vy;
  m0 = ( - _b[0] - _sim->last_time()*_vt[0] ) / _vy[0];
  m0 = ( - _b[0] ) / _vy[0];
  m1 = ( - _b[0] - time[1]*_vt[0] ) / _vy[0];
  m2 = ( - _b[0] - time[2]*_vt[0] ) / _vy[0];

  trace3("corr4", time[0], time[1], time[2]);
  trace3("corr4", _dttstate[0], _dttstate[1], _dttstate[2]);
  trace3("corr4", _ttstate[0], _ttstate[1], _ttstate[2]);
  trace3("fit4", _vt, _vy, _b);
  trace3("corr4", m0, m1, m2);
  trace4("corr4", d1, _dv, _vt[0], _vy[0]);
  trace1("tau4", _dv-d1);

  gsl_multifit_linear_residuals (X, y, c, r);
  trace4("res4", gsl_vector_get(r,0), gsl_vector_get(r,1), gsl_vector_get(r,2), gsl_vector_get(r,3));

  double extime = _sim->_dT0;

  double e = exp(-extime*fabs((_dv-_dttstate[1])/(_ttstate[1]-m1)));
  e = exp(-extime*fabs((_dv-d1)/(_ttstate[1]-m1)));
  double corr = (_ttstate[1] - m1) * e + m0;
  trace4("ex4", e, corr, m0, m1);
  return corr;
}
#endif
/*--------------------------------------------------------------------------*/
TIME_PAIR DEV_CAPACITANCE::tt_review()
{
  _new_ttsteporder = _ttsteporder + 1;;
  _pred = _ttstate[0]; // predicted state
  _preddt = _dttstate[0];
  _dttstate[0] = (_y[0].f0 - _ttstate[0]) / _sim->_time0;
  _dttstate_ = _y[0].f0 - _ttstate[0];
  _ttstate[0] = _y[0].f0;
  assert(is_number(_y[0].f0));
  double ds[4];
  ds[0] = _dttstate[0];
  ds[1] = _dttstate[1];
  ds[2] = _dttstate[2];
  ds[3] = _dttstate[3];
  double time[4];
  time[0] = 0;
  time[1] = - _sim->_dT0;
  time[2] = - _sim->_dT0 - _dT1;
  time[3] = - _sim->_dT0 - _dT1 - _dT2;
  double t, d1;
  USE(t); USE(d1);
  _corr = NAN; // nan("");

  //  - _ttstate[0]/ _y[0].f1;
  switch (_ttmethod) {
    case 0:
      break;
    case 1:

//       _ttstate[0] += (_ttstate[1] + _sim->_dT0 * (_dttstate[0] + _dttstate[1])/2);
//       _ttstate[0] /= 2;
//       _dttstate[0] = (_ttstate[0] - _ttstate[1])/_sim->_dT0;
      break;
#ifdef HAVE_GSL_FIT_H
    case 2:
      _corr = gsl_fit_3();
    break;
    case 3:
      _corr = gsl_fit_3();
      break;
    case 4:
      _corr = gsl_fit_4();
      break;
#else
    case 2:
    case 3:
    case 4:
      incomplete();
      break;
#endif
    default: incomplete();
  }

  assert(is_number(_ttstate[0]));
  assert(is_number(_dttstate[0]));

//  _dttstate[0] *= 1/1.5;


  trace4("rev", long_label(), ds[0], ds[1], ds[2]);
  trace4("rev", long_label(), time[0], time[1], time[2]);

  derivatives(ds, _ttsteporder, time);

  _ttfuture = NEVER;

  double errf[5];
  errf[0] = 1e-2;
  errf[1] = 1e-9;
  errf[2] = 1e-10; // gsl3
  errf[3] = 1e-8; // gsl4
  errf[4] = 1e-12;

  double tol = errf[_ttsteporder] * OPT::tttol * (1./7.);
  double denom = 0;
  if (_ttsteporder>0) {
    denom = fabs(ds[_ttsteporder-1]);
  }
  double chargetol;

  _ttstep = NEVER;

  switch (_ttmethod) {
    case 0:
      _ttstep = _sim->_dT0;
      _ttfuture = tt_review_check_and_convert(_ttstep);
      break;
    case 1:
      assert(denom != 0);
      _ttstep = sqrt( tol / denom );
      _ttfuture = tt_review_check_and_convert(_ttstep);
      assert(_ttfuture >= _sim->_dTmin);
      break;
    case 2:
      assert(denom != 0);
      _ttstep = OPT::tttol / 7. * _sim->_dTmin; // nonsense.
      _ttfuture = tt_review_check_and_convert(_ttstep);
      break;
    case 3:
      assert(is_number(denom));
      _ttstep = OPT::tttol / 7. *_sim->_dT0 * 1.3; // hack. necessary?
      if(denom != 0){
	_ttstep1 = pow( tol / denom, 1/3. );
	_ttstep = min(_ttstep, _ttstep1);
	assert(is_number(_ttstep));
      }else{
	assert(is_number(_ttstep));
      }

      _ttstep3 = (fabs(OPT::tttol/(7*_dv)))*_sim->_dTmin;

      denom = fabs(ds[2]);
      _ttstep2 = pow( OPT::tttol / (7*denom), 1/2. );

      denom = fabs(_pred - _corr);
      chargetol = std::max(OPT::chgtol, OPT::reltol * std::max((double)std::abs(_corr), (double) std::abs(_pred)));
      _ttstep1 = pow( OPT::tttol * chargetol / denom, 1/2. );

      _ttstep = min(_ttstep, _ttstep1);
      _ttstep = min(_ttstep, _ttstep2);
      _ttstep = min(_ttstep, _ttstep3);

      _ttfuture = tt_review_check_and_convert(_ttstep);

      assert(_ttfuture >= _sim->_dTmin);
      break;
    case 4:
      denom = fabs(_dttstate[0] - _preddt);
      _ttstep0 = pow( 1e-6 / denom, .5);

      chargetol = std::max(OPT::chgtol, OPT::reltol * std::max((double)std::abs(_corr), (double) std::abs(_pred)));
      assert(is_number(_corr));
      denom = fabs(_pred - _corr);
      _ttstep1 = pow( OPT::tttol * chargetol / denom, 1/2. );

      denom = fabs(ds[3]);
      _ttstep2 = pow( 1e6 / denom, 1/2. );

      _ttstep3 = (fabs(1/_dv))*_sim->_dTmin;

//      _ttstep = _sim->_dT0 * 1.3;
      _ttstep = min(_ttstep, _ttstep3);
      _ttstep = min(_ttstep, _ttstep2);
      _ttstep = min(_ttstep, _ttstep1);
      _ttstep = min(_ttstep, _ttstep0);
      _ttfuture = tt_review_check_and_convert(_ttstep);
      break;
    default:
      incomplete();
  }

  if (_new_ttsteporder>int(OPT::ttsteporder)) {
    _new_ttsteporder = OPT::ttsteporder;
  }else if (_new_ttsteporder>4) { untested();
    _new_ttsteporder = 4;
  }
 // return TIME_PAIR(NEVER,NEVER);
  return TIME_PAIR(_ttfuture,NEVER);
}
/*--------------------------------------------------------------------------*/
void DEV_CAPACITANCE::do_tt()
{
  // set _ttstate[0] (temporarily) and y0.f0.
  // and _dttstate[0]
  //
  // use {_vt _vy _b} [1]
  //
  // _ttstate[1] is state at end of last timeframe. i.e.
  // at Time0 - dT0 + last_time()

  STORAGE::do_tt();

  double m1, m0, e;
  assert(_sim->last_time() || !_ttmethod);
  _extime = _sim->_dT0 - _sim->last_time();
  assert(_extime>=0);
  switch (_ttmethod) {
    case 0:
      _y[0].f0 = _ttstate[0];
      break;
    case 1:
      _y[0].f0 = _ttstate[1] + _extime * _dttstate[1];
      break;
    case 2:
      _y[0].f0 = _ttstate[1] + _extime * .5 * (_dttstate[1] + _dttstate[2]);
      break;
    case 3:
    //  _y[0].f0 = _ttstate[1] + _extime * .5 * (_dttstate[1] + _dttstate[2]);
    //  assert(is_number(_y[0].f0));
    //  break;
    case 4:
//      m1 = - _b / _vy;
      m1 = ( - _b[1] - _sim->last_time()*_vt[1] ) / _vy[1];
      m0 = ( - _b[1] - _sim->_dT0*_vt[1] ) / _vy[1];
      _dv = - _vt[1]/_vy[1];
      assert(is_number(_vt[1]));
      assert(is_number(_b[1]));
      assert(is_number(_vy[1]));

      trace4("do_tt 3", _sim->_Time0, _sim->_dT0, _ttstate[1], _ttstate[2]);
      trace2("do_tt 3", _dttstate[1], _dttstate[2]);
      trace3("do_tt 3", _vy, _vt, _b);
      trace5("do_tt 3", m1, m0, _dv, _sim->_dT0, _extime);
      assert(_extime>=0);
      e = exp(-_extime*fabs((_dv-_dttstate[1])/(_ttstate[1]-m1)));

      trace2("do_tt 3", e, _ttstate[1] - m1);

      _y[0].f0 = (_ttstate[1] - m1) * e + m0;

      assert(is_number(_y[0].f0));

      _preddt =
      _dttstate[0] = _sim->_dT0*_vt[1]+_y[0].f0*_vy[1]+_b[1];
      break;
    default: unreachable();
  }
  _ttstate[0] = _y[0].f0;
  trace3("do_tt done", _ttstate[1], _ttstate[0], _ttmethod);
  store_values();
}
/*--------------------------------------------------------------------------*/
void DEV_CAPACITANCE::do_ac()
{ untested();
  if (using_ac_eval()) { untested();
    ac_eval();
  }else{ untested();
    assert(_ev == _y[0].f1);
    assert(has_tr_eval() || _ev == hp_float_t(value()));
  }
  _acg =  (COMPLEX)_ev * _sim->_jomega;
}
/*--------------------------------------------------------------------------*/
double DEV_CAPACITANCE::tr_probe_num(const std::string& x)const
{
  if (Umatch(x, "q{cap} |ch{arge} ")) {
    return _y[0].f0;
  }else if (Umatch(x, "c{apacitance} ")) { untested();
    return _y[0].f1;
  }else if (Umatch(x, "dcdt ")) {untested();
    return (_y[0].f1 - _y[1].f1) / _dt;
  }else if (Umatch(x, "dc ")) {untested();
    return (_y[0].f1 - _y[1].f1);
  }else if (Umatch(x, "dqdt ")) { untested();
    return (_y[0].f0 - _y[1].f0) / _dt;
  }else if (Umatch(x, "dq ")) { untested();
    return (_y[0].f0 - _y[1].f0);
  }else if (Umatch(x, "ttstep ")) {
    return _ttstep;
  }else if (Umatch(x, "ttgain ")) {
    return _ttstep/_sim->_dTmin;
  }else if (Umatch(x, "ttstep0 ")) {
    return _ttstep0/_sim->_dTmin;
  }else if (Umatch(x, "ttstep1 ")) {
    return _ttstep1/_sim->_dTmin;
  }else if (Umatch(x, "ttstep2 ")) {
    return _ttstep2/_sim->_dTmin;
  }else if (Umatch(x, "ttstep3 ")) {
    return _ttstep3/_sim->_dTmin;
  }else if (Umatch(x, "tt{future} ")) {
    return _ttfuture;
  }else if (Umatch(x, "ttstate")) {
    return (_ttstate[0]);
  }else if (Umatch(x, "ttord{er} ")) {
    return (_ttsteporder);
  }else if (Umatch(x, "ttdelta")) {
    return (_dttstate_);
  }else if (Umatch(x, "vt")) {
    return (_vt[0]);
  }else if (Umatch(x, "vy")) {
    return (_vy[0]);
  }else if (Umatch(x, "dv")) {
    return (_dv);
  }else if (Umatch(x, "extime")) {
    // bug hunting...
    return (_extime);
  }else if (Umatch(x, "ttstdt")) {
    return (_dttstate[0]);
  }else if (Umatch(x, "preddt ")) {
    return (_preddt);
  }else if (Umatch(x, "corr ")) {
    return (_corr);
  }else if (Umatch(x, "pred ")) {
    return (_pred);
  }else if (Umatch(x, "chisq ")) {
    return (_chisq);
  }else if (Umatch(x, "ttmeth{od} ")) {
    return (_ttmethod);
  }else{
    return STORAGE::tr_probe_num(x);
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
DEV_CAPACITANCE p1;
DISPATCHER<CARD>::INSTALL
  d1(&device_dispatcher, "ttcap",	    &p1);
}
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet
