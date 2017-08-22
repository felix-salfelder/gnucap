/*                     -*- C++ -*-
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
 * mos model equations: spice level 1 equivalent
 */
/* This file is not automatically generated. */

#include "globals.h"
#include "e_elemnt.h"
#include "d_mos1.h"
namespace UF{
/*--------------------------------------------------------------------------*/
const double NA(NOT_INPUT);
const double INF(BIGBIG);
/*--------------------------------------------------------------------------*/
int MODEL_BUILT_IN_MOS1::_count = 0;
/*--------------------------------------------------------------------------*/
const int LEVEL(1);
/*--------------------------------------------------------------------------*/
namespace MODEL_BUILT_IN_MOS1_DISPATCHER { 
  static DEV_BUILT_IN_MOS p1d;
  static MODEL_BUILT_IN_MOS1 p1(&p1d);
  static DISPATCHER<MODEL_CARD>::INSTALL
    d1(&model_dispatcher, "nmos1|pmos1|nmos|pmos", &p1);
}
/*--------------------------------------------------------------------------*/
void SDP_BUILT_IN_MOS1::init(const COMMON_COMPONENT* cc)
{
  assert(cc);
  SDP_BUILT_IN_MOS123::init(cc);
}
/*--------------------------------------------------------------------------*/
TDP_BUILT_IN_MOS1::TDP_BUILT_IN_MOS1(const DEV_BUILT_IN_MOS* d)
  :TDP_BUILT_IN_MOS123(d)
{
  assert(d);
  const COMMON_BUILT_IN_MOS* c = prechecked_cast<const COMMON_BUILT_IN_MOS*>(d->common());
  assert(c);
  const SDP_BUILT_IN_MOS1* s = prechecked_cast<const SDP_BUILT_IN_MOS1*>(c->sdp());
  assert(s);
  const MODEL_BUILT_IN_MOS1* m = prechecked_cast<const MODEL_BUILT_IN_MOS1*>(c->model());
  assert(m);
  const CARD_LIST* par_scope = d->scope();
  assert(par_scope);
  USE(par_scope);
    // final adjust: code_pre

      double temp = d->_sim->_temp_c + P_CELSIUS0;
      double tempratio  = temp / m->tnom_k;
      double tempratio4 = tempratio * sqrt(tempratio);
      double kt = temp * P_K;
      double vt = temp * P_K_Q;
      double egap_ = 1.16 - (7.02e-4*temp*temp) / (temp+1108.);
      double arg = (m->egap*tempratio - egap_) / (2*kt);
    // final adjust: override
    // final adjust: raw
    // final adjust: mid
    // final adjust: calculated
    phi = m->phi*tempratio + (-2*vt*(1.5*log(tempratio)+P_Q*(arg)));
    beta = (m->kp / tempratio4) * s->w_eff / s->l_eff;
    dbetadw_eff = (m->kp / tempratio4) / s->l_eff;
    dbetadl_eff = - (m->kp / tempratio4) * s->w_eff / s->l_eff / s->l_eff;
    sqrt_phi = sqrt(phi);
    egap = egap_;
    // final adjust: post
    // final adjust: done
}
/*--------------------------------------------------------------------------*/
MODEL_BUILT_IN_MOS1::MODEL_BUILT_IN_MOS1(const BASE_SUBCKT* p)
  :MODEL_BUILT_IN_MOS123(p),
   kp(NA),
   calc_kp(false)
{
  if (ENV::run_mode != rPRE_MAIN) {
    ++_count;
  }else{
  }
  set_default(&mjsw, .5);
  set_default(&cox, NA);
  set_default(&vto, NA);
  set_default(&gamma, NA);
  set_default(&phi, NA);
  set_default(&mos_level, LEVEL);
}
/*--------------------------------------------------------------------------*/
MODEL_BUILT_IN_MOS1::MODEL_BUILT_IN_MOS1(const MODEL_BUILT_IN_MOS1& p)
  :MODEL_BUILT_IN_MOS123(p),
   kp(p.kp),
   calc_kp(p.calc_kp)
{
  if (ENV::run_mode != rPRE_MAIN) {
    ++_count;
  }else{untested();//194
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_MOS1::dev_type()const
{
  if (polarity == pN) {
    return "nmos1";
  }else if (polarity == pP) {
    return "pmos1";
  }else if (polarity == pN) {
    return "nmos";
  }else if (polarity == pP) {
    return "pmos";
  }else{untested();//235
    return MODEL_BUILT_IN_MOS123::dev_type();
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_MOS1::set_dev_type(const std::string& new_type)
{
  if (Umatch(new_type, "nmos1 ")) {
    polarity = pN;
  }else if (Umatch(new_type, "pmos1 ")) {
    polarity = pP;
  }else if (Umatch(new_type, "nmos ")) {
    polarity = pN;
  }else if (Umatch(new_type, "pmos ")) {
    polarity = pP;
  }else{
    MODEL_BUILT_IN_MOS123::set_dev_type(new_type);
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_MOS1::precalc_first()
{
    const CARD_LIST* par_scope = scope();
    assert(par_scope);
    MODEL_BUILT_IN_MOS123::precalc_first();
    e_val(&(this->kp), NA, par_scope);
    // final adjust: code_pre

      if (tox != NA) {
	cox = P_EPS_OX / tox;
	if (kp == NA) {
	  kp = uo * cox;
	  calc_kp = true;
	}
	if (nsub != NA) {
	  if (phi == NA) {
	    phi = (2. * P_K_Q) * tnom_k * log(nsub/NI);
	    if (phi < .1) {
	      untested();
	      error(((!_sim->is_first_expand()) ? (bDEBUG) : (bWARNING)),
		    long_label() + ": calculated phi too small, using .1\n");
	      phi = .1;
	    }
	    calc_phi = true;
	  }
	  if (gamma == NA) {
	    gamma = sqrt(2. * P_EPS_SI * P_Q * nsub) / cox;
	    calc_gamma = true;
	  }
	  if (vto == NA) {
	    double phi_ms = (tpg == gtMETAL)
	      ? polarity * (-.05 - (egap + polarity * phi) / 2.)
	      : -(tpg * egap + phi) / 2.;
	    double vfb = phi_ms - polarity * P_Q * nss / cox;
	    vto = vfb + phi + gamma * sqrt(phi);
	    calc_vto = true;
	  }
	}else{
	  // tox is input, nsub isn't
	}
      }
    // final adjust: override
    if (cox == NA) {
      cox = 0.;
    }else{
    }
    if (vto == NA) {
      vto = 0.;
    }else{
    }
    if (gamma == NA) {
      gamma = 0.;
    }else{
    }
    if (phi == NA) {
      phi = .6;
    }else{
    }
    // final adjust: raw
    e_val(&(this->kp), 2e-5, par_scope);
    // final adjust: mid
    // final adjust: calculated
    // final adjust: post
    // final adjust: done
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_MOS1::precalc_last()
{
    MODEL_BUILT_IN_MOS123::precalc_last();
}
/*--------------------------------------------------------------------------*/
SDP_CARD* MODEL_BUILT_IN_MOS1::new_sdp(COMMON_COMPONENT* c)const
{
  assert(c);
  if (COMMON_BUILT_IN_MOS* cc = dynamic_cast<COMMON_BUILT_IN_MOS*>(c)) {
    if (cc->_sdp) {
      cc->_sdp->init(cc);
      return cc->_sdp;
    }else{
      delete cc->_sdp;
      return new SDP_BUILT_IN_MOS1(c);
    }
  }else{
    return MODEL_BUILT_IN_MOS123::new_sdp(c);
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_MOS1::set_param_by_index(int i, std::string& value, int offset)
{
  switch (MODEL_BUILT_IN_MOS1::param_count() - 1 - i) {
  case 0: level = value; break; //1
  case 1: unreachable(); break;
  case 2: unreachable(); break;
  case 3: unreachable(); break;
  case 4: unreachable(); break;
  case 5: unreachable(); break;
  case 6: mos_level = value; break;
  case 7: kp = value; break;
  default: MODEL_BUILT_IN_MOS123::set_param_by_index(i, value, offset); break;
  }
}
/*--------------------------------------------------------------------------*/
bool MODEL_BUILT_IN_MOS1::param_is_printable(int i)const
{
  switch (MODEL_BUILT_IN_MOS1::param_count() - 1 - i) {
  case 0:  return (true);
  case 1:  return (false);
  case 2:  return (false);
  case 3:  return (false);
  case 4:  return (false);
  case 5:  return (false);
  case 6:  return (mos_level != LEVEL);
  case 7:  return (!calc_kp);
  default: return MODEL_BUILT_IN_MOS123::param_is_printable(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_MOS1::param_name(int i)const
{
  switch (MODEL_BUILT_IN_MOS1::param_count() - 1 - i) {
  case 0:  return "level";
  case 1:  return "=====";
  case 2:  return "=====";
  case 3:  return "=====";
  case 4:  return "=====";
  case 5:  return "=====";
  case 6:  return "diodelevel";
  case 7:  return "kp";
  default: return MODEL_BUILT_IN_MOS123::param_name(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_MOS1::param_name(int i, int j)const
{
  if (j == 0) {
    return param_name(i);
  }else if (j == 1) {
    switch (MODEL_BUILT_IN_MOS1::param_count() - 1 - i) {
    case 0:  return "";
    case 1:  return "";
    case 2:  return "";
    case 3:  return "";
    case 4:  return "";
    case 5:  return "";
    case 6:  return "";
    case 7:  return "";
    default: return MODEL_BUILT_IN_MOS123::param_name(i, j);
    }
  }else if (i < 8) {
    return "";
  }else{
    return MODEL_BUILT_IN_MOS123::param_name(i, j);
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_MOS1::param_value(int i)const
{
  switch (MODEL_BUILT_IN_MOS1::param_count() - 1 - i) {
  case 0:  return "1";
  case 1:  unreachable(); return "";
  case 2:  unreachable(); return "";
  case 3:  unreachable(); return "";
  case 4:  unreachable(); return "";
  case 5:  unreachable(); return "";
  case 6:  return mos_level.string();
  case 7:  return kp.string();
  default: return MODEL_BUILT_IN_MOS123::param_value(i);
  }
}
/*--------------------------------------------------------------------------*/
bool MODEL_BUILT_IN_MOS1::is_valid(const COMPONENT* d)const
{
  assert(d);
  return MODEL_BUILT_IN_MOS123::is_valid(d);
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_MOS1::tr_eval(COMPONENT* brh)const
{
  DEV_BUILT_IN_MOS* d = prechecked_cast<DEV_BUILT_IN_MOS*>(brh);
  assert(d);
  const COMMON_BUILT_IN_MOS* c = prechecked_cast<const COMMON_BUILT_IN_MOS*>(d->common());
  assert(c);
  const ADP_BUILT_IN_MOS* a = prechecked_cast<const ADP_BUILT_IN_MOS*>(d->adp());
  assert(a);
  const SDP_BUILT_IN_MOS1* s = prechecked_cast<const SDP_BUILT_IN_MOS1*>(c->sdp());
  assert(s); USE(s);
  const MODEL_BUILT_IN_MOS1* m = this;
  const TDP_BUILT_IN_MOS1 T(d);
  const TDP_BUILT_IN_MOS1* t = &T;

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
    trace0(d->long_label().c_str());
    trace3("", d->vds, d->vgs, d->vbs);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
    d->reverse_if_needed();
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
    double sarg, dsarg_dvbs;
    {
      if (d->vbs <= 0.) {
	sarg = sqrt(t->phi - d->vbs);
	dsarg_dvbs = -.5 / sarg;
	d->sbfwd = false;
	trace2("sb-ok", sarg, dsarg_dvbs);
      }else{
	sarg = t->sqrt_phi / (1. + .5 * d->vbs / t->phi);
	dsarg_dvbs = -.5 * sarg * sarg / t->phi*t->sqrt_phi; /* is wrong!! */
	d->sbfwd = true;
	trace2("***sb-reversed***", sarg, dsarg_dvbs);
      }
    }
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
    //d->vto += ((const DEV_BUILT_IN_BTI*)(d->_BTI))->dvth();

    d->von = m->vto + m->gamma * (sarg - sqrt(m->phi))
      + .5 * (m->egap - t->egap) + .5 * (t->phi - m->phi);

    d->von += a->delta_vth;   

    d->vgst = d->vdsat = d->vgs - d->von;
    if (d->vdsat < 0.) {
      d->vdsat = 0.;
    }
    d->cutoff = (d->vgst < 0.);
    d->saturated = (d->vds > d->vdsat);
    trace3("MODEL_BUILT_IN_MOS1", d->von, d->vgst, d->vdsat);
    double Lambda = (m->lambda != NA) ? m->lambda : 0.;
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
    if (d->cutoff) {
      d->gds = d->gmf = d->ids = d->gmbf = 0.;
      trace4("cut", d->ids, d->gmf, d->gds, d->gmbf);
    }else if (d->saturated) {
      d->gmf  = t->beta * d->vgst * (1. + Lambda * d->vds);
      d->ids = d->gmf * (.5 * d->vgst);
      d->gds = .5 * t->beta * Lambda * d->vgst * d->vgst;
      d->gmbf = - d->gmf * m->gamma * dsarg_dvbs;
      trace4("sat", d->ids, d->gmf, d->gds, d->gmbf);
    }else{ /* triode */
      d->gmf  = t->beta * d->vds * (1. + Lambda * d->vds);
      d->ids = d->gmf * (d->vgst - .5*d->vds);
      d->gds = t->beta * ((d->vgst - d->vds) 
			 + Lambda * d->vds * (2.*d->vgst - 1.5*d->vds));
      d->gmbf = -d->gmf * m->gamma * dsarg_dvbs;
      trace4("lin", d->ids, d->gmf, d->gds, d->gmbf);
    }
    if (d->reversed) {
      d->ids *= -1;
      d->gmr = d->gmf;
      d->gmbr = d->gmbf;
      d->gmf = d->gmbf = 0;
    }else{
      d->gmr = d->gmbr = 0.;
    }
}
/*--------------------------------------------------------------------------*/
double MODEL_BUILT_IN_MOS1::dgmdl_eff(const DEV_BUILT_IN_MOS* brh) const
{
  const DEV_BUILT_IN_MOS* d = prechecked_cast<const DEV_BUILT_IN_MOS*>(brh);
  assert(d);
  const MODEL_BUILT_IN_MOS1* m = this;
  const TDP_BUILT_IN_MOS1 T(d);
  const TDP_BUILT_IN_MOS1* t = &T;
  double Lambda = (m->lambda != NA) ? m->lambda : 0.;
  if (d->cutoff) {
    return 0;
  }else if (d->saturated) {
    return t->dbetadl_eff * d->vgst * (1. + Lambda * d->vds);
  }else{ /* triode */
    return t->dbetadl_eff * d->vds * (1. + Lambda * d->vds);
  }
}
/*--------------------------------------------------------------------------*/
double MODEL_BUILT_IN_MOS1::dgmdw_eff(const DEV_BUILT_IN_MOS* brh) const
{
  const DEV_BUILT_IN_MOS* d = prechecked_cast<const DEV_BUILT_IN_MOS*>(brh);
  assert(d);
  const MODEL_BUILT_IN_MOS1* m = this;
  const TDP_BUILT_IN_MOS1 T(d);
  const TDP_BUILT_IN_MOS1* t = &T;
  double Lambda = (m->lambda != NA) ? m->lambda : 0.;
  if (d->cutoff) {
    return 0;
  }else if (d->saturated) {
    return t->dbetadw_eff * d->vgst * (1. + Lambda * d->vds);
  }else{ /* triode */
    return t->dbetadw_eff * d->vds * (1. + Lambda * d->vds);
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
}
// vim:sw=2:ts=8:
