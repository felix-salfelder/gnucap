/*                        -*- C++ -*-
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
 * data structures and defaults for mos model.
 * internal units all mks (meters)
 * but some user input parameters are in cm.
 *
 * netlist syntax:
 * device:  mxxxx d g s b mname <device args> <model card args>
 * model:   .model mname NMOS <args>
 *	or  .model mname PMOS <args>
 */
/* This file is not automatically generated */

#include "u_limit.h"
#include "io_trace.h"
#include "e_storag.h"
#include "d_mos_base.h"
#include "e_adp.h"
#include "e_adp_mos.h"
#include "e_adplist.h"
#include "globals.h"
#include "e_elemnt.h"
#include "d_mos.h"
#include "u_lang.h"
#include "u_xprobe.h"

// doesnt work without yet.
// ( maybe better, probes will not work otherwise... )
#define BTI_IN_SUBCKT
#define HCI_IN_SUBCKT

#define BTI_LATE_EVAL
namespace UF{
/*--------------------------------------------------------------------------*/
const double NA(NOT_INPUT);
const double INF(BIGBIG);
/*--------------------------------------------------------------------------*/
int DEV_BUILT_IN_MOS::_count = -1;
int COMMON_BUILT_IN_MOS::_count = -1;
static COMMON_BUILT_IN_MOS Default_BUILT_IN_MOS(CC_STATIC);
/*--------------------------------------------------------------------------*/
#ifndef BTI_IN_SUBCKT
#define BTI_LATE_EVAL
#endif
/*--------------------------------------------------------------------------*/
COMMON_BUILT_IN_MOS::COMMON_BUILT_IN_MOS(int c)
  :COMMON_COMPONENT(c),
   l_in(OPT::defl),
   w_in(OPT::defw),
   ad_in(OPT::defad),
   as_in(OPT::defas),
   pd(0.0),
   ps(0.0),
   nrd(1.0),
   nrs(1.0),
   _sdp(0),
   _db(0),
   _sb(0)
{
  ++_count;
}
/*--------------------------------------------------------------------------*/
COMMON_BUILT_IN_MOS::COMMON_BUILT_IN_MOS(const COMMON_BUILT_IN_MOS& p)
  :COMMON_COMPONENT(p),
   l_in(p.l_in),
   w_in(p.w_in),
   ad_in(p.ad_in),
   as_in(p.as_in),
   pd(p.pd),
   ps(p.ps),
   nrd(p.nrd),
   nrs(p.nrs),
   _sdp(0),
   _db(0),
   _sb(0)
{
//  trace0("some mos copy");
  ++_count;
}
/*--------------------------------------------------------------------------*/
COMMON_BUILT_IN_MOS::~COMMON_BUILT_IN_MOS()
{
  detach_common(&_db);
  detach_common(&_sb);
  --_count;
  delete _sdp;
}
/*--------------------------------------------------------------------------*/
bool COMMON_BUILT_IN_MOS::operator==(const COMMON_COMPONENT& x)const
{
  const COMMON_BUILT_IN_MOS* p = dynamic_cast<const COMMON_BUILT_IN_MOS*>(&x);
  return (p
    && l_in == p->l_in
    && w_in == p->w_in
    && ad_in == p->ad_in
    && as_in == p->as_in
    && pd == p->pd
    && ps == p->ps
    && nrd == p->nrd
    && nrs == p->nrs
    && _sdp == p->_sdp
    && COMMON_COMPONENT::operator==(x));
}
/*--------------------------------------------------------------------------*/
void COMMON_BUILT_IN_MOS::set_param_by_index(int I, std::string& Value, int Offset)
{ unreachable(); // obsolete.
  switch (COMMON_BUILT_IN_MOS::param_count() - 1 - I) { untested();
  case 0:  l_in = Value; break;
  case 1:  w_in = Value; break;
  case 2:  ad_in = Value; break;
  case 3:  as_in = Value; break;
  case 4:  pd = Value; break;
  case 5:  ps = Value; break;
  case 6:  nrd = Value; break;
  case 7:  nrs = Value; break;
  default: COMMON_COMPONENT::set_param_by_index(I, Value, Offset);
  }
}
/*--------------------------------------------------------------------------*/
// FIXME: use map
void COMMON_BUILT_IN_MOS::set_param_by_name(string Name, string Value)
{
  trace2("COMMON_BUILT_IN_MOS::set_param_by_name", Name, Value);
  if (Umatch (Name,"l")) { l_in = Value; }
  else if (Umatch (Name,"w"))  { w_in  = Value; }
  else if (Umatch (Name,"ad")) { ad_in = Value; }
  else if (Umatch (Name,"as")) { as_in = Value; }
  else if (Umatch (Name,"pd")) { pd    = Value; }
  else if (Umatch (Name,"ps")) { ps    = Value; }
  else if (Umatch (Name,"nrd")){ nrd   = Value; }
  else if (Umatch (Name,"nrs")){ nrs   = Value; }
  else { untested();
    return COMMON_COMPONENT::set_param_by_name(Name, Value);
  }
}
/*--------------------------------------------------------------------------*/
bool COMMON_BUILT_IN_MOS::param_is_printable(int i)const
{
  switch (COMMON_BUILT_IN_MOS::param_count() - 1 - i) { untested();
  case 0:  return (true);
  case 1:  return (true);
  case 2:  return (has_hard_value(ad_in));
  case 3:  return (has_hard_value(as_in));
  case 4:  return (has_hard_value(pd));
  case 5:  return (has_hard_value(ps));
  case 6:  return (has_hard_value(nrd));
  case 7:  return (has_hard_value(nrs));
  default: return COMMON_COMPONENT::param_is_printable(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string COMMON_BUILT_IN_MOS::param_name(int i)const
{
  switch (COMMON_BUILT_IN_MOS::param_count() - 1 - i) { untested();
  case 0:  return "l";
  case 1:  return "w";
  case 2:  return "ad";
  case 3:  return "as";
  case 4:  return "pd";
  case 5:  return "ps";
  case 6:  return "nrd";
  case 7:  return "nrs";
  default: return COMMON_COMPONENT::param_name(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string COMMON_BUILT_IN_MOS::param_name(int i, int j)const
{ untested();
  if (j == 0) { untested();
    return param_name(i);
  }else if (j == 1) { untested();
    switch (COMMON_BUILT_IN_MOS::param_count() - 1 - i) { untested();
    case 0:  return "";
    case 1:  return "";
    case 2:  return "";
    case 3:  return "";
    case 4:  return "";
    case 5:  return "";
    case 6:  return "";
    case 7:  return "";
    default: return "";
    }
  }else{untested();//281
    return COMMON_COMPONENT::param_name(i, j);
  }
}
/*--------------------------------------------------------------------------*/
std::string COMMON_BUILT_IN_MOS::param_value(int i)const
{
  switch (COMMON_BUILT_IN_MOS::param_count() - 1 - i) { untested();
  case 0:  return l_in.string();
  case 1:  return w_in.string();
  case 2:  return ad_in.string();
  case 3:  return as_in.string();
  case 4:  return pd.string();
  case 5:  return ps.string();
  case 6:  return nrd.string();
  case 7:  return nrs.string();
  default: return COMMON_COMPONENT::param_value(i);
  }
}
/*--------------------------------------------------------------------------*/
void COMMON_BUILT_IN_MOS::expand(const COMPONENT* d)
{
  COMMON_COMPONENT::expand(d);
  attach_model(d);
  COMMON_BUILT_IN_MOS* c = this;
  const MODEL_BUILT_IN_MOS_BASE* m = dynamic_cast<const MODEL_BUILT_IN_MOS_BASE*>(model());
  if (!m) {
    throw Exception_Model_Type_Mismatch(d->long_label(), modelname(), "mosfet");
  }else{
  }
  // size dependent
  //delete _sdp;
  _sdp = m->new_sdp(this);
  assert(_sdp);
  const SDP_BUILT_IN_MOS_BASE* s = prechecked_cast<const SDP_BUILT_IN_MOS_BASE*>(_sdp);
  assert(s);

  // subcircuit commons, recursive
  COMMON_BUILT_IN_DIODE* db = new COMMON_BUILT_IN_DIODE;
  db->area = double(s->ad);
  db->perim = double(c->pd);
  db->is_raw = double(s->idsat);
  db->cj_raw = double(m->cbd);
  db->cjsw_raw = NA;
  db->off = true;
  db->set_modelname(modelname());
  db->attach(model());
  attach_common(db, &_db);

  COMMON_BUILT_IN_DIODE* sb = new COMMON_BUILT_IN_DIODE;
  sb->area = double(s->as);
  sb->perim = double(c->ps);
  sb->is_raw = double(s->issat);
  sb->cj_raw = double(m->cbs);
  sb->cjsw_raw = NA;
  sb->off = true;
  sb->set_modelname(modelname());
  sb->attach(model());
  attach_common(sb, &_sb);
}
/*--------------------------------------------------------------------------*/
void COMMON_BUILT_IN_MOS::precalc_first(const CARD_LIST* par_scope)
{
  assert(par_scope);
  COMMON_COMPONENT::precalc_first(par_scope);
    e_val(&(this->l_in), OPT::defl, par_scope);
    e_val(&(this->w_in), OPT::defw, par_scope);
    e_val(&(this->ad_in), OPT::defad, par_scope);
    e_val(&(this->as_in), OPT::defas, par_scope);
    e_val(&(this->pd), 0.0, par_scope);
    e_val(&(this->ps), 0.0, par_scope);
    e_val(&(this->nrd), 1.0, par_scope);
    e_val(&(this->nrs), 1.0, par_scope);
}
/*--------------------------------------------------------------------------*/
void COMMON_BUILT_IN_MOS::precalc_last(const CARD_LIST* par_scope)
{
  assert(par_scope);
  COMMON_COMPONENT::precalc_last(par_scope);
  COMMON_BUILT_IN_MOS* c = this;
  const MODEL_BUILT_IN_MOS_BASE* m = prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(model());
    // final adjust: code_pre
    // final adjust: override
    // final adjust: raw
    e_val(&(this->l_in), OPT::defl, par_scope);
    e_val(&(this->w_in), OPT::defw, par_scope);
    e_val(&(this->ad_in), OPT::defad, par_scope);
    e_val(&(this->as_in), OPT::defas, par_scope);
    e_val(&(this->pd), 0.0, par_scope);
    e_val(&(this->ps), 0.0, par_scope);
    e_val(&(this->nrd), 1.0, par_scope);
    e_val(&(this->nrs), 1.0, par_scope);
    // final adjust: mid
    // final adjust: calculated
    // final adjust: post
    // final adjust: done

  // size dependent
  //delete _sdp;
  _sdp = m->new_sdp(this);
  assert(_sdp);
  const SDP_BUILT_IN_MOS_BASE* s = prechecked_cast<const SDP_BUILT_IN_MOS_BASE*>(_sdp);
  assert(s);

  // subcircuit commons, recursive
  COMMON_BUILT_IN_DIODE* db = new COMMON_BUILT_IN_DIODE;
  db->area = double(s->ad);
  db->perim = double(c->pd);
  db->is_raw = double(s->idsat);
  db->cj_raw = double(m->cbd);
  db->cjsw_raw = NA;
  db->off = true;
  db->set_modelname(modelname());
  db->attach(model());
  attach_common(db, &_db);

  COMMON_BUILT_IN_DIODE* sb = new COMMON_BUILT_IN_DIODE;
  sb->area = double(s->as);
  sb->perim = double(c->ps);
  sb->is_raw = double(s->issat);
  sb->cj_raw = double(m->cbs);
  sb->cjsw_raw = NA;
  sb->off = true;
  sb->set_modelname(modelname());
  sb->attach(model());
  attach_common(sb, &_sb);

  assert(c == this);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
namespace DEV_BUILT_IN_MOS_DISPATCHER { 
  static DEV_BUILT_IN_MOS p0;
  static DISPATCHER<CARD>::INSTALL
    d0(&device_dispatcher, "M|mosfet", &p0);
}
/*--------------------------------------------------------------------------*/
static EVAL_BUILT_IN_MOS_Cgb Eval_Cgb(CC_STATIC);
void EVAL_BUILT_IN_MOS_Cgb::tr_eval(ELEMENT* d)const
{
  assert(d);
  DEV_BUILT_IN_MOS* p = prechecked_cast<DEV_BUILT_IN_MOS*>(d->owner());
  assert(p);
  const COMMON_BUILT_IN_MOS* c = prechecked_cast<const COMMON_BUILT_IN_MOS*>(p->common());
  assert(c);
  const SDP_BUILT_IN_MOS_BASE* s = prechecked_cast<const SDP_BUILT_IN_MOS_BASE*>(c->sdp());
  assert(s);
  const MODEL_BUILT_IN_MOS_BASE* m = prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(c->model());
  assert(m);

    STORAGE* brh = prechecked_cast<STORAGE*>(d);
    assert(brh);

    double cap = brh->value();

    if (m->cmodel != 0) {
      if (p->vgst < - s->phi) { 		/* accumulation */
	cap += s->cgate;
      }else if (p->vgst < 0.) {			/* depletion */
	cap += s->cgate * (-p->vgst) / s->phi;
      }else{					/* active, overlap only */
      }
    }
    brh->_y[0].f1 = cap;
    if (d->_sim->analysis_is_tran_dynamic()) {
      cap = (brh->_y[0].f1 + brh->_y[1].f1) / 2;
      brh->_y[0].f0 = (brh->_y[0].x - brh->_y[1].x) * cap + brh->_y[1].f0;
    }else{
      assert(d->_sim->analysis_is_static() || d->_sim->analysis_is_restore());
      brh->_y[0].f0 = brh->_y[0].x * brh->_y[0].f1;
    }
    trace3(brh->long_label().c_str(), brh->_y[0].x, brh->_y[0].f0, brh->_y[0].f1);
}
/*--------------------------------------------------------------------------*/
static EVAL_BUILT_IN_MOS_Cgd Eval_Cgd(CC_STATIC);
void EVAL_BUILT_IN_MOS_Cgd::tr_eval(ELEMENT* d)const
{
  assert(d);
  DEV_BUILT_IN_MOS* p = prechecked_cast<DEV_BUILT_IN_MOS*>(d->owner());
  assert(p);
  const COMMON_BUILT_IN_MOS* c = prechecked_cast<const COMMON_BUILT_IN_MOS*>(p->common());
  assert(c);
  const SDP_BUILT_IN_MOS_BASE* s = prechecked_cast<const SDP_BUILT_IN_MOS_BASE*>(c->sdp());
  assert(s);
  const MODEL_BUILT_IN_MOS_BASE* m = prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(c->model());
  assert(m);

    STORAGE* brh = prechecked_cast<STORAGE*>(d);
    assert(brh);

    double cap = 0;
    if (m->cmodel != 0) {
      assert(p->vdsat >= 0.);
      assert(p->vds >= 0.);
      double vbs    = (m->cmodel == 3) ? 0. : p->vbs;
      double vdbsat = p->vdsat - vbs;
      double vdb    = p->vds   - vbs;
      double ddif   = 2. * vdbsat - vdb;
      
      if (!p->reversed) { // treat as Cgs
	if (p->vgst >= 0.) {
	  if (p->vdsat > p->vds) {		/* linear */
	    cap = (2./3.) * s->cgate * (1. - (vdbsat*vdbsat)/(ddif*ddif));
	    if (p->vgst <= .1) {
	      cap *= 10. * p->vgst;	// smooth discontinuity
	    }
	  }
	}
      }else{ // treat as Cgs
	if (p->vgst >= -s->phi/2.) {		/* depletion  or active */
	  cap = (2./3.) * s->cgate;
	  if (p->vdsat > p->vds) {			/* linear */
	    double ndif   = p->vdsat - p->vds;
	    cap *= 1. - (ndif*ndif)/(ddif*ddif);
	  }
	  if (p->vgst <= 0) {
	    cap *= 1. + p->vgst / (s->phi);
	    cap *= 1. + p->vgst / (s->phi);
	  }
	}
      }
    }
    cap += brh->value();		/* else overlap only */
    
    brh->_y[0].f1 = cap;
    if (d->_sim->analysis_is_tran_dynamic()) {
      cap = (brh->_y[0].f1 + brh->_y[1].f1) / 2;
      brh->_y[0].f0 = (brh->_y[0].x - brh->_y[1].x) * cap + brh->_y[1].f0;
    }else{
      assert(d->_sim->analysis_is_static() || d->_sim->analysis_is_restore());
      brh->_y[0].f0 = brh->_y[0].x * brh->_y[0].f1;
    }
    trace3(brh->long_label().c_str(), brh->_y[0].x, brh->_y[0].f0, brh->_y[0].f1);
}
/*--------------------------------------------------------------------------*/
static EVAL_BUILT_IN_MOS_Cgs Eval_Cgs(CC_STATIC);
void EVAL_BUILT_IN_MOS_Cgs::tr_eval(ELEMENT* d)const
{
  assert(d);
  DEV_BUILT_IN_MOS* p = prechecked_cast<DEV_BUILT_IN_MOS*>(d->owner());
  assert(p);
  const COMMON_BUILT_IN_MOS* c = prechecked_cast<const COMMON_BUILT_IN_MOS*>(p->common());
  assert(c);
  const SDP_BUILT_IN_MOS_BASE* s = prechecked_cast<const SDP_BUILT_IN_MOS_BASE*>(c->sdp());
  assert(s);
  const MODEL_BUILT_IN_MOS_BASE* m = prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(c->model());
  assert(m);

    STORAGE* brh = prechecked_cast<STORAGE*>(d);
    assert(brh);

    double cap = 0;
    if (m->cmodel != 0) {
      assert(p->vdsat >= 0.);
      assert(p->vds >= 0.);
      double vbs    = (m->cmodel == 3) ? 0. : p->vbs;
      double vdbsat = p->vdsat - vbs;
      double vdb    = p->vds   - vbs;
      double ddif   = 2. * vdbsat - vdb;
      
      if (p->reversed) { // treat as Cgd
	if (p->vgst >= 0.) {
	  if (p->vdsat > p->vds) {		/* linear */
	    cap = (2./3.) * s->cgate * (1. - (vdbsat*vdbsat)/(ddif*ddif));
	    if (p->vgst <= .1) {
	      cap *= 10. * p->vgst;	// smooth discontinuity
	    }
	  }
	}
      }else{ // treat as Cgs
	if (p->vgst >= -s->phi/2.) {		/* depletion  or active */
	  cap = (2./3.) * s->cgate;
	  if (p->vdsat > p->vds) {			/* linear */
	    double ndif   = p->vdsat - p->vds;
	    cap *= 1. - (ndif*ndif)/(ddif*ddif);
	  }
	  if (p->vgst <= 0) {
	    cap *= 1. + p->vgst / (s->phi);
	    cap *= 1. + p->vgst / (s->phi);
	  }
	}
      }
    }
    cap += brh->value();		/* else overlap only */
    
    brh->_y[0].f1 = cap;
    if (d->_sim->analysis_is_tran_dynamic()) {
      cap = (brh->_y[0].f1 + brh->_y[1].f1) / 2;
      brh->_y[0].f0 = (brh->_y[0].x - brh->_y[1].x) * cap + brh->_y[1].f0;
    }else{
      assert(d->_sim->analysis_is_static() || d->_sim->analysis_is_restore());
      brh->_y[0].f0 = brh->_y[0].x * brh->_y[0].f1;
    }
    trace3(brh->long_label().c_str(), brh->_y[0].x, brh->_y[0].f0, brh->_y[0].f1);
}
/*--------------------------------------------------------------------------*/
DEV_BUILT_IN_MOS::DEV_BUILT_IN_MOS()
  :BASE_SUBCKT(),
   // input parameters,
   // calculated parameters,
   ids(0.),
   idsxxx(NA),
   gds(0.),
   gmf(0.),
   gmr(0.),
   gmbf(0.),
   gmbr(0.),
   idb(0.),
   idbxxx(0.),
   gdbdb(0.),
   gdbds(0.),
   gdbgs(0.),
   gdbbs(0.),
   isb(0.),
   isbxxx(0.),
   gsbsb(0.),
   gsbsd(0.),
   gsbgd(0.),
   gsbbd(0.),
   qgate(0.),
   cgs(0.),
   cggb(0.),
   cgsb(0.),
   cgdb(0.),
   qgs(0.),
   cgsgs(0.),
   cgsgb(0.),
   cgssb(0.),
   cgsdb(0.),
   qgd(0.),
   cgdgd(0.),
   cgdgb(0.),
   cgdsb(0.),
   cgddb(0.),
   qdrn(0.),
   cdsds(0.),
   cdgb(0.),
   cdsb(0.),
   cddb(0.),
   qbulk(0.),
   cbs(0.),
   cbgb(0.),
   cbsb(0.),
   cbdb(0.),
   qbs(0.),
   cbsbs(0.),
   cbsgb(0.),
   cbssb(0.),
   cbsdb(0.),
   qbd(0.),
   cbdbd(0.),
   cbdgb(0.),
   cbdsb(0.),
   cbddb(0.),
   gtau(0.),
   cqgb(0.),
   cqsb(0.),
   cqdb(0.),
   cqbb(0.),
   vgs(0.),
   vds(0.),
   vbs(0.),
   vdsat(0.),
   vgst(0.),
   von(0.),
   reversed(false),
   cutoff(false),
   subthreshold(false),
   saturated(false),
   sbfwd(false),
   punchthru(false),
   // netlist,
   _Rs(0),
   _Rd(0),
   _Ddb(0),
   _Dsb(0),
   _Cgs(0),
   _Cgd(0),
   _Cgb(0),
   _Cqgs(0),
   _Cqgd(0),
   _Cqds(0),
   _Cqbs(0),
   _Cqbd(0),
   _BTI(0),
   _HCI(0),
   _Ids(0),
   _Idb(0),
   _Isb(0)
{
  _n = _nodes;
  attach_common(&Default_BUILT_IN_MOS);
  ++_count;
  // overrides
}
/*--------------------------------------------------------------------------*/
DEV_BUILT_IN_MOS::DEV_BUILT_IN_MOS(const DEV_BUILT_IN_MOS& p)
  :BASE_SUBCKT(p),
   // input parameters,
   // calculated parameters,
   ids(p.ids),
   idsxxx(p.idsxxx),
   gds(p.gds),
   gmf(p.gmf),
   gmr(p.gmr),
   gmbf(p.gmbf),
   gmbr(p.gmbr),
   idb(p.idb),
   idbxxx(p.idbxxx),
   gdbdb(p.gdbdb),
   gdbds(p.gdbds),
   gdbgs(p.gdbgs),
   gdbbs(p.gdbbs),
   isb(p.isb),
   isbxxx(p.isbxxx),
   gsbsb(p.gsbsb),
   gsbsd(p.gsbsd),
   gsbgd(p.gsbgd),
   gsbbd(p.gsbbd),
   qgate(p.qgate),
   cgs(p.cgs),
   cggb(p.cggb),
   cgsb(p.cgsb),
   cgdb(p.cgdb),
   qgs(p.qgs),
   cgsgs(p.cgsgs),
   cgsgb(p.cgsgb),
   cgssb(p.cgssb),
   cgsdb(p.cgsdb),
   qgd(p.qgd),
   cgdgd(p.cgdgd),
   cgdgb(p.cgdgb),
   cgdsb(p.cgdsb),
   cgddb(p.cgddb),
   qdrn(p.qdrn),
   cdsds(p.cdsds),
   cdgb(p.cdgb),
   cdsb(p.cdsb),
   cddb(p.cddb),
   qbulk(p.qbulk),
   cbs(p.cbs),
   cbgb(p.cbgb),
   cbsb(p.cbsb),
   cbdb(p.cbdb),
   qbs(p.qbs),
   cbsbs(p.cbsbs),
   cbsgb(p.cbsgb),
   cbssb(p.cbssb),
   cbsdb(p.cbsdb),
   qbd(p.qbd),
   cbdbd(p.cbdbd),
   cbdgb(p.cbdgb),
   cbdsb(p.cbdsb),
   cbddb(p.cbddb),
   gtau(p.gtau),
   cqgb(p.cqgb),
   cqsb(p.cqsb),
   cqdb(p.cqdb),
   cqbb(p.cqbb),
   vgs(p.vgs),
   vds(p.vds),
   vbs(p.vbs),
   vdsat(p.vdsat),
   vgst(p.vgst),
   von(p.von),
   reversed(p.reversed),
   cutoff(p.cutoff),
   subthreshold(p.subthreshold),
   saturated(p.saturated),
   sbfwd(p.sbfwd),
   punchthru(p.punchthru),
   // netlist,
   _Rs(0),
   _Rd(0),
   _Ddb(0),
   _Dsb(0),
   _Cgs(0),
   _Cgd(0),
   _Cgb(0),
   _Cqgs(0),
   _Cqgd(0),
   _Cqds(0),
   _Cqbs(0),
   _Cqbd(0),
   _BTI(0),
   _HCI(0),
   _Ids(0),
   _Idb(0),
   _Isb(0)
{
  _n = _nodes;
  for (uint_t ii = 0; ii < max_nodes() + int_nodes(); ++ii) {
    _n[ii] = p._n[ii];
  }
  ++_count;
  // overrides
}
/*--------------------------------------------------------------------------*/
DEV_BUILT_IN_MOS::~DEV_BUILT_IN_MOS()
{
  --_count;
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_MOS::expand()
{
  BASE_SUBCKT::expand(); // COMPON::expand()
  assert(_n);
  assert(common());
  const COMMON_BUILT_IN_MOS* c = static_cast<const COMMON_BUILT_IN_MOS*>(common());
  assert(c);
  assert(c->model());
  const MODEL_BUILT_IN_MOS_BASE* m = prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(c->model());
  assert(m);
  assert(c->sdp());
  const SDP_BUILT_IN_MOS_BASE* s = prechecked_cast<const SDP_BUILT_IN_MOS_BASE*>(c->sdp());
  assert(s);
  if (!subckt()) {
    new_subckt(scope()->params());
  }else{
  }

  if (_sim->is_first_expand()) {
    precalc_first();
    precalc_last();
    // local nodes
    //assert(!(_n[n_id].n_()));
    //BUG// this assert fails on a repeat elaboration after a change.
    //not sure of consequences when new_model_node called twice.
    if (!(_n[n_id].n_())) {
      if (!OPT::rstray || s->rd == 0.) {
        _n[n_id] = _n[n_d];
      }else{
        trace0("DEV_BUILT_IN_MOS::expand new_model_node id");
        _n[n_id].new_model_node("id", this); // increase global counter.
      }
    }else{
      if (!OPT::rstray || s->rd == 0.) {
        assert(_n[n_id] == _n[n_d]);
      }else{
        //_n[n_id].new_model_node("id." + long_label(), this);
      }
    }
    //assert(!(_n[n_is].n_()));
    //BUG// this assert fails on a repeat elaboration after a change.
    //not sure of consequences when new_model_node called twice.
    if (!(_n[n_is].n_())) {
      if (!OPT::rstray || s->rs == 0.) {
        _n[n_is] = _n[n_s];
      }else{
        _n[n_is].new_model_node("is", this);
      }
    }else{
      if (!OPT::rstray || s->rs == 0.) {
        assert(_n[n_is] == _n[n_s]);
      }else{ untested();
        //_n[n_is].new_model_node("is." + long_label(), this);
      }
    }

    // clone subckt elements
    if (!OPT::rstray || s->rs == 0.) {
      if (_Rs) { untested();
        subckt()->erase(_Rs);
        _Rs = NULL;
      }else{
      }
    }else{
      if (!_Rs) {
        const CARD* p = device_dispatcher["resistor"];
        assert(p);
        _Rs = dynamic_cast<COMPONENT*>(p->clone());
        assert(_Rs);
        subckt()->push_front(_Rs);
      }else{ untested();
      }
      {
        node_t nodes[] = {_n[n_s], _n[n_is]};
        _Rs->set_parameters("Rs", this, NULL, s->rs, 0, NULL, 2, nodes);
      }
    }
    if (!OPT::rstray || s->rd == 0.) {
      if (_Rd) { untested();
        subckt()->erase(_Rd);
        _Rd = NULL;
      }else{
      }
    }else{
      if (!_Rd) {
        const CARD* p = device_dispatcher["resistor"];
        assert(p);
        _Rd = dynamic_cast<COMPONENT*>(p->clone());
        assert(_Rd);
        subckt()->push_front(_Rd);
      }else{
      }
      {
        node_t nodes[] = {_n[n_d], _n[n_id]};
        _Rd->set_parameters("Rd", this, NULL, s->rd, 0, NULL, 2, nodes);
      }
    }
    if (_n[n_b].n_() == _n[n_d].n_() || s->idsat == 0.) {
      if (_Ddb) { untested();
        subckt()->erase(_Ddb);
        _Ddb = NULL;
      }else{
      }
    }else{
      if (!_Ddb) {
        const CARD* p = device_dispatcher["diode"];
        assert(p);
        _Ddb = dynamic_cast<COMPONENT*>(p->clone());
        assert(_Ddb);
        subckt()->push_front(_Ddb);
      }else{
      }
      if (m->polarity==pP) {
        node_t nodes[] = {_n[n_id], _n[n_b]  };
        _Ddb->set_parameters("Ddb", this, c->_db, 0., 0, NULL, 2, nodes);
      }else{
        node_t nodes[] = {_n[n_b], _n[n_id]};
        _Ddb->set_parameters("Ddb", this, c->_db, 0., 0, NULL, 2, nodes);
      }
    }
    if (_n[n_b].n_() == _n[n_s].n_() || s->issat == 0.) {
      if (_Dsb) { untested();
        subckt()->erase(_Dsb);
        _Dsb = NULL;
      }else{
      }
    }else{
      if (!_Dsb) {
        const CARD* p = device_dispatcher["diode"];
        assert(p);
        _Dsb = dynamic_cast<COMPONENT*>(p->clone());
        assert(_Dsb);
        subckt()->push_front(_Dsb);
      }else{
      }
      if (m->polarity==pP) {
        node_t nodes[] = {_n[n_is], _n[n_b]  };
        _Dsb->set_parameters("Dsb", this, c->_sb, 0., 0, NULL, 2, nodes);
      }else{
        node_t nodes[] = {_n[n_b], _n[n_is]};
        _Dsb->set_parameters("Dsb", this, c->_sb, 0., 0, NULL, 2, nodes);
      }
    }
    if (m->cmodel != 0 || !OPT::cstray 
		|| _n[n_g].n_() == _n[n_is].n_()) {
      if (_Cqgs) {
        subckt()->erase(_Cqgs);
        _Cqgs = NULL;
      }else{
      }
    }else{
      if (!_Cqgs) {
        const CARD* p = device_dispatcher["fpoly_cap"];
        assert(p);
        _Cqgs = dynamic_cast<COMPONENT*>(p->clone());
        assert(_Cqgs);
        subckt()->push_front(_Cqgs);
      }else{
      }
      {
        node_t nodes[] = {_n[n_g], _n[n_is], _n[n_g], _n[n_b], _n[n_is], _n[n_b], _n[n_id], _n[n_b]};
      _Cqgs->set_parameters("Cqgs", this, NULL, 0., 5, &qgs, 8, nodes);
      }
    }
    if (m->cmodel != 0 || !OPT::cstray 
		|| _n[n_g].n_() == _n[n_id].n_()) {
      if (_Cqgd) {
        subckt()->erase(_Cqgd);
        _Cqgd = NULL;
      }else{
      }
    }else{
      if (!_Cqgd) {
        const CARD* p = device_dispatcher["fpoly_cap"];
        assert(p);
        _Cqgd = dynamic_cast<COMPONENT*>(p->clone());
        assert(_Cqgd);
        subckt()->push_front(_Cqgd);
      }else{
      }
      {
        node_t nodes[] = {_n[n_g], _n[n_id], _n[n_g], _n[n_b], _n[n_id], _n[n_b], _n[n_is], _n[n_b]};
      _Cqgd->set_parameters("Cqgd", this, NULL, 0., 5, &qgd, 8, nodes);
      }
    }
    if (m->cmodel != 0 || !OPT::cstray 
             || _n[n_id].n_() == _n[n_is].n_()) {
      if (_Cqds) {
        subckt()->erase(_Cqds);
        _Cqds = NULL;
      }else{
      }
    }else{
      if (!_Cqds) {
        const CARD* p = device_dispatcher["fpoly_cap"];
        assert(p);
        _Cqds = dynamic_cast<COMPONENT*>(p->clone());
        assert(_Cqds);
        subckt()->push_front(_Cqds);
      }else{
      }
      {
        node_t nodes[] = {_n[n_id], _n[n_is], _n[n_g], _n[n_b], _n[n_is], _n[n_b], _n[n_id], _n[n_b]};
      _Cqds->set_parameters("Cqds", this, NULL, 0., 5, &qdrn, 8, nodes);
      }
    }
    if (m->cmodel != 0 || !OPT::cstray 
             || _n[n_b].n_() == _n[n_is].n_()) {
      if (_Cqbs) {
        subckt()->erase(_Cqbs);
        _Cqbs = NULL;
      }else{
      }
    }else{
      if (!_Cqbs) {
        const CARD* p = device_dispatcher["fpoly_cap"];
        assert(p);
        _Cqbs = dynamic_cast<COMPONENT*>(p->clone());
        assert(_Cqbs);
        subckt()->push_front(_Cqbs);
      }else{
      }
      {
        node_t nodes[] = {_n[n_b], _n[n_is], _n[n_g], _n[n_b], _n[n_is], _n[n_b], _n[n_id], _n[n_b]};
      _Cqbs->set_parameters("Cqbs", this, NULL, 0., 5, &qbs, 8, nodes);
      }
    }
    if (m->cmodel != 0 || !OPT::cstray 
             || _n[n_b].n_() == _n[n_id].n_()) {
      if (_Cqbd) {
        subckt()->erase(_Cqbd);
        _Cqbd = NULL;
      }else{
      }
    }else{
      if (!_Cqbd) {
        const CARD* p = device_dispatcher["fpoly_cap"];
        assert(p);
        _Cqbd = dynamic_cast<COMPONENT*>(p->clone());
        assert(_Cqbd);
        subckt()->push_front(_Cqbd);
      }else{
      }
      {
        node_t nodes[] = {_n[n_b], _n[n_id], _n[n_g], _n[n_b], _n[n_id], _n[n_b], _n[n_is], _n[n_b]};
      _Cqbd->set_parameters("Cqbd", this, NULL, 0., 5, &qbd, 8, nodes);
      }
    }
    if (!OPT::cstray || _n[n_g].n_() == _n[n_s].n_()) {
      if (_Cgs) {
        subckt()->erase(_Cgs);
        _Cgs = NULL;
      }else{
      }
    }else{
      if (!_Cgs) {
        const CARD* p = device_dispatcher["capacitor"];
        assert(p);
        _Cgs = dynamic_cast<COMPONENT*>(p->clone());
        assert(_Cgs);
        subckt()->push_front(_Cgs);
      }else{
      }
      {
        node_t nodes[] = {_n[n_g], _n[n_is]};
        _Cgs->set_parameters("Cgs", this, &Eval_Cgs, s->cgso, 0, NULL, 2, nodes);
      }
    }
    if (!OPT::cstray || _n[n_g].n_() == _n[n_d].n_()) {
      if (_Cgd) {
        subckt()->erase(_Cgd);
        _Cgd = NULL;
      }else{
      }
    }else{
      if (!_Cgd) {
        const CARD* p = device_dispatcher["capacitor"];
        assert(p);
        _Cgd = dynamic_cast<COMPONENT*>(p->clone());
        assert(_Cgd);
        subckt()->push_front(_Cgd);
      }else{
      }
      {
        node_t nodes[] = {_n[n_g], _n[n_id]};
        _Cgd->set_parameters("Cgd", this, &Eval_Cgd, s->cgdo, 0, NULL, 2, nodes);
      }
    }
    if (!OPT::cstray || _n[n_b].n_() == _n[n_g].n_()) {
      if (_Cgb) {
        subckt()->erase(_Cgb);
        _Cgb = NULL;
      }else{
      }
    }else{
      if (!_Cgb) {
        const CARD* p = device_dispatcher["capacitor"];
        assert(p);
        _Cgb = dynamic_cast<COMPONENT*>(p->clone());
        assert(_Cgb);
        subckt()->push_front(_Cgb);
      }else{
      }
      {
        node_t nodes[] = {_n[n_g], _n[n_b]};
        _Cgb->set_parameters("Cgb", this, &Eval_Cgb, s->cgbo, 0, NULL, 2, nodes);
      }
    }

    {
      if (!_Ids) {
        const CARD* p = device_dispatcher["cpoly_g"];
        assert(p);
        _Ids = dynamic_cast<COMPONENT*>(p->clone());
        assert(_Ids);
        subckt()->push_front(_Ids);
      }else{
      }
      {
        node_t nodes[] = {_n[n_id], _n[n_is], _n[n_g], _n[n_is], _n[n_id],
          _n[n_g], _n[n_b], _n[n_is], _n[n_id], _n[n_b]};
        _Ids->set_parameters("Ids", this, NULL, 0., 6, &idsxxx, 10, nodes);
      }
    }
    if (!(m->needs_isub) || _n[n_d].n_() == _n[n_b].n_()) {
      if (_Idb) { untested();
        subckt()->erase(_Idb);
        _Idb = NULL;
      }else{
      }
    }else{ untested();
      if (!_Idb) { untested();
        const CARD* p = device_dispatcher["cpoly_g"];
        assert(p);
        _Idb = dynamic_cast<COMPONENT*>(p->clone());
        assert(_Idb);
        subckt()->push_front(_Idb);
      }else{ untested();
      }
      { untested();
        node_t nodes[] = {_n[n_id], _n[n_b], _n[n_id], _n[n_is], _n[n_g], _n[n_is], _n[n_b], _n[n_is]};
        _Idb->set_parameters("Idb", this, NULL, 0., 5, &idbxxx, 8, nodes);
      }
    }
    if (!(m->needs_isub) || _n[n_s].n_() == _n[n_b].n_()) {
      if (_Isb) { untested();
        subckt()->erase(_Isb);
        _Isb = NULL;
      }else{
      }
    }else{ untested();
      if (!_Isb) { untested();
        const CARD* p = device_dispatcher["cpoly_g"];
        assert(p);
        _Isb = dynamic_cast<COMPONENT*>(p->clone());
        assert(_Isb);
        subckt()->push_front(_Isb);
      }else{ untested();
      }
      { untested();
        node_t nodes[] = {_n[n_is], _n[n_b], _n[n_is], _n[n_id], _n[n_g], _n[n_id], _n[n_b], _n[n_id]};
        _Isb->set_parameters("Isb", this, NULL, 0., 5, &isbxxx, 8, nodes);
      }
    }

    trace2("DEV_BUILT_IN_MOS::expand bti expand", m->polarity, m->use_bti());
    if (!m->use_bti()) {
      if (_BTI) { untested();
#ifdef BTI_IN_SUBCKT
        subckt()->erase(_BTI);
        _BTI = NULL;
#else
        delete _BTI;
#endif
      }
    } else {
      if (!_BTI) {
	string btimodelname = m->bti_model_name.value();
        const CARD* p = LANGUAGE::find_proto(btimodelname, this);
	if(!p) error(bDANGER, "mos: cannot find btimodel %s\n", btimodelname.c_str());
	if (!p) { itested();
	  throw(Exception_Cant_Find(long_label(), btimodelname));
	}
        _BTI = dynamic_cast<COMPONENT*>(p->clone_instance());
	_BTI->set_dev_type(btimodelname); // check: could be obsolete...
        assert(_BTI);
#ifdef BTI_IN_SUBCKT
        subckt()->push_front(_BTI);
#endif
      }
      {
	node_t nodes[] = {_n[n_id], _n[n_g], _n[n_is]};
	double Lcoeff = 1;
	_BTI->set_parameters("BTI", this, _BTI->mutable_common(), m->polarity * Lcoeff,
	    0 /*states*/, &vgs, 3, nodes);
      }
    }
    if (!m->use_hci()) {
      if (_HCI) { untested();
#ifdef HCI_IN_SUBCKT
        subckt()->erase(_HCI);
#else
	delete dynamic_cast<CARD*>(_HCI);
#endif
        _HCI = NULL;
      }
    } else {
      if (!_HCI) {
        const CARD* p = m->hci_device();
	if (!p) { itested();
	  throw(Exception_Cant_Find(long_label(), "hcimodel"));
	}
        _HCI = dynamic_cast<COMPONENT*>(p->clone_instance());
        assert(_HCI);
//	_HCI->attach_common(mutable_common());
#ifdef HCI_IN_SUBCKT
        subckt()->push_front(_HCI);
#endif
      }
      {
	// attaching same common. make sure hci::precalc_* is shortened....
	_HCI->set_parameters("HCI", this, mutable_common(), 0,
	    0 /*states*/, NULL, 0 /*numnodes*/, NULL /*nodes*/); // BUG!!!
      }
    }
  }else{ untested();
    //precalc();
  }

  assert(!is_constant());
  subckt()->set_slave();
  trace0("mos: expanding sckt");
  subckt()->expand();

#ifndef HCI_IN_SUBCKT
  if(_HCI){ untested();
    _HCI->expand();
  }
#endif

#ifndef BTI_IN_SUBCKT
  if( m->use_bti() ){ untested();
    _BTI->expand();
    _BTI->set_slave();
  }
#endif


  trace0(("DEV_BUILT_IN_MOS::expand, ADP things " + long_label()).c_str());
  if (! adp()) {
    ADP_CARD* A = m->new_adp( (COMPONENT*) this );
    if(A) {
      attach_adp(A);
    }else{
    }
  }

  if (adp()) {
    // adp()->expand();
  }else{
  }

  // adp()->expand();
//  adp()->q_accept();
//  adp()->q_eval();
  trace0(("DEV_BUILT_IN_MOS::expanded, ADP things " + long_label()).c_str());
}
/*--------------------------------------------------------------------------*/
double DEV_BUILT_IN_MOS::tr_probe_num(const std::string& x)const
{
  const DEV_BUILT_IN_MOS* d = this;
  assert(_n);
  USE(d);
  const COMMON_BUILT_IN_MOS* c = prechecked_cast<const COMMON_BUILT_IN_MOS*>(common());
  assert(c);
  const MODEL_BUILT_IN_MOS_BASE* m = prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(c->model());
  assert(m);
  const SDP_BUILT_IN_MOS_BASE* s = prechecked_cast<const SDP_BUILT_IN_MOS_BASE*>(c->sdp());
  assert(s);
  const ADP_BUILT_IN_MOS* a = prechecked_cast<const ADP_BUILT_IN_MOS*>(adp());

  if (Umatch(x, "v ")) {
    return  _n[n_d].v0() - _n[n_s].v0();
  }else if (Umatch(x, "vds ")) {
    return  _n[n_d].v0() - _n[n_s].v0();
  }else if (Umatch(x, "dvth ")) {
    assert(a);
    return a->delta_vth;
  }else if (Umatch(x, "bti_stress ")) { untested();
    return  a->bti_stress->tr_get();
  }else if (Umatch(x, "vgs ")) {
    return  _n[n_g].v0() - _n[n_s].v0();
  }else if (Umatch(x, "vbs ")) {
    return  _n[n_b].v0() - _n[n_s].v0();
  }else if (Umatch(x, "vdsi{nt} ")) {
    return  vds;
  }else if (Umatch(x, "vgsi{nt} ")) { untested();
    return  vgs;
  }else if (Umatch(x, "vbsi{nt} ")) { untested();
    return  vbs;
  }else if (Umatch(x, "use_bti ")) { untested();
    return  m->use_bti();
  }else if (Umatch(x, "use_hci ")) {
    return  m->use_hci();
  }else if (Umatch(x, "vgd ")) { untested();
    return  _n[n_g].v0() - _n[n_d].v0();
  }else if (Umatch(x, "vbd ")) {
    return  _n[n_b].v0() - _n[n_d].v0();
  }else if (Umatch(x, "vsd ")) { untested();
    return  _n[n_s].v0() - _n[n_d].v0();
  }else if (Umatch(x, "vdm ")) {
    return  ( _n[n_d].v0() - _n[n_s].v0() + _n[n_d].v0() - _n[n_d].v0() ) / 2.;
  }else if (Umatch(x, "vgm ")) {
    return  ( _n[n_g].v0() - _n[n_s].v0() + _n[n_g].v0() - _n[n_d].v0() ) / 2.;
  }else if (Umatch(x, "vbm ")) {
    return  ( _n[n_b].v0() - _n[n_s].v0() + _n[n_b].v0() - _n[n_d].v0() ) / 2.;
  }else if (Umatch(x, "vsm ")) {
    return  ( _n[n_s].v0() - _n[n_s].v0() + _n[n_s].v0() - _n[n_d].v0() ) / 2.;
  }else if (Umatch(x, "vdg ")) { untested();
    return  _n[n_d].v0() - _n[n_g].v0();
  }else if (Umatch(x, "vbg ")) { untested();
    return  _n[n_b].v0() - _n[n_g].v0();
  }else if (Umatch(x, "vsg ")) { untested();
    return  _n[n_s].v0() - _n[n_g].v0();
  }else if (Umatch(x, "vdb ")) { untested();
    return  _n[n_d].v0() - _n[n_b].v0();
  }else if (Umatch(x, "vgb ")) { untested();
    return  _n[n_g].v0() - _n[n_b].v0();
  }else if (Umatch(x, "vsb ")) { untested();
    return  _n[n_s].v0() - _n[n_b].v0();
  }else if (Umatch(x, "vd ")) {
    return  _n[n_d].v0();
  }else if (Umatch(x, "vid ")) {
    return  _n[n_id].v0();
  }else if (Umatch(x, "vg ")) {
    return  _n[n_g].v0();
  }else if (Umatch(x, "vb ")) {
    return  _n[n_b].v0();
  }else if (Umatch(x, "vs ")) {
    return  _n[n_s].v0();
  }else if (Umatch(x, "vis ")) {
    return  _n[n_is].v0();
  }else if (Umatch(x, "i{d} ")) {
    return  (_Rd) ? CARD::probe(_Rd,"I") : ( CARD::probe(_Ids,"I") - CARD::probe(_Cgd,"I") - CARD::probe(_Ddb,"I") * m->polarity);
  }else if (Umatch(x, "is ")) {
    return  (_Rs) ? CARD::probe(_Rs,"I") : (- CARD::probe(_Ids,"I") - CARD::probe(_Cgs,"I") - CARD::probe(_Dsb,"I") * m->polarity);
  }else if (Umatch(x, "ig ")) {
    return  CARD::probe(_Cgs,"I") + CARD::probe(_Cgd,"I") + CARD::probe(_Cgb,"I");
  }else if (Umatch(x, "ib ")) {
    return  - CARD::probe(_Ddb,"I") * m->polarity - CARD::probe(_Dsb,"I") * m->polarity - CARD::probe(_Cgb,"I");
  }else if (Umatch(x, "ibd ")) {
    return  CARD::probe(_Ddb,"I");
  }else if (Umatch(x, "ibs ")) {
    return  CARD::probe(_Dsb,"I");
  }else if (Umatch(x, "cgso{vl} ")) {
    return  CARD::probe(_Cgs,"NV");
  }else if (Umatch(x, "cgdo{vl} ")) {
    return  CARD::probe(_Cgd,"NV");
  }else if (Umatch(x, "cgbo{vl} ")) {
    return  CARD::probe(_Cgb,"NV");
  }else if (Umatch(x, "cgst ")) {
    return  CARD::probe(_Cgs,"EV");
  }else if (Umatch(x, "cgdt ")) {
    return  CARD::probe(_Cgd,"EV");
  }else if (Umatch(x, "cgbt ")) {
    return  CARD::probe(_Cgb,"EV");
  }else if (Umatch(x, "cgs{m} ")) {
    return  CARD::probe(_Cgs,"EV") - CARD::probe(_Cgs,"NV");
  }else if (Umatch(x, "cgd{m} ")) {
    return  CARD::probe(_Cgd,"EV") - CARD::probe(_Cgd,"NV");
  }else if (Umatch(x, "cgb{m} ")) {
    return  CARD::probe(_Cgb,"EV") - CARD::probe(_Cgb,"NV");
  }else if (Umatch(x, "cbd ")) {
    return  CARD::probe(_Ddb,"Cap");
  }else if (Umatch(x, "cbs ")) {
    return  CARD::probe(_Dsb,"Cap");
  }else if (Umatch(x, "cgate ")) {
    return  s->cgate;
  }else if (Umatch(x, "gm ")) {
    return  (reversed) ? gmr : gmf;
  }else if (Umatch(x, "gmb{s} ")) {
    return  (reversed) ? gmbr : gmbf;
  }else if (Umatch(x, "gbd ")) {
    return  CARD::probe(_Ddb,"G");
  }else if (Umatch(x, "gbs ")) {
    return  CARD::probe(_Dsb,"G");
  }else if (Umatch(x, "vth ")) {
    return  von * m->polarity;
  }else if (Umatch(x, "ids ")) {
    return  m->polarity * ((reversed) ? -ids : ids);
  }else if (Umatch(x, "idst{ray} ")) {
    return  - CARD::probe(_Cgd,"I") + CARD::probe(_Ddb,"I") * m->polarity;
  }else if (Umatch(x, "p ")) {
    return  CARD::probe(_Rs,"P") + CARD::probe(_Rd,"P") + CARD::probe(_Ddb,"P") + CARD::probe(_Dsb,"P") + CARD::probe(_Cgs,"P") + CARD::probe(_Cgd,"P") + CARD::probe(_Cgb,"P") + CARD::probe(_Ids,"P");
  }else if (Umatch(x, "pd ")) {
    return  CARD::probe(_Rs,"PD") + CARD::probe(_Rd,"PD") + CARD::probe(_Ddb,"PD") + CARD::probe(_Dsb,"PD") + CARD::probe(_Cgs,"PD") + CARD::probe(_Cgd,"PD") + CARD::probe(_Cgb,"PD") + CARD::probe(_Ids,"PD");
  }else if (Umatch(x, "ps ")) {
    return  CARD::probe(_Rs,"PS") + CARD::probe(_Rd,"PS") + CARD::probe(_Ddb,"PS") + CARD::probe(_Dsb,"PS") + CARD::probe(_Cgs,"PS") + CARD::probe(_Cgd,"PS") + CARD::probe(_Cgb,"PS") + CARD::probe(_Ids,"PS");
  }else if (Umatch(x, "REgion ")) {
    return  static_cast<double>((!cutoff) + (!subthreshold * 2) + (saturated * 4) + (sbfwd * 10) + ((vbs > vds) * 20) + (punchthru * 40)) * ((reversed)? -1 : 1);
  }else if (Umatch(x, "SUBthreshold ")) { untested();
    return  static_cast<double>(subthreshold);
  }else if (Umatch(x, "CUToff ")) { untested();
    return  static_cast<double>(cutoff);
  }else if (Umatch(x, "SATurated ")) { untested();
    return  static_cast<double>(saturated);
  }else if (Umatch(x, "TRIode ")) { untested();
    return  static_cast<double>(!saturated && !subthreshold);
  }else if (Umatch(x, "SBFwd ")) { untested();
    return  static_cast<double>(sbfwd);
  }else if (Umatch(x, "DBFwd ")) { untested();
    return  static_cast<double>(vbs > vds);
  }else if (Umatch(x, "REVersed ")) {
    return  static_cast<double>(reversed);
  }else if (Umatch(x, "status ")) { untested();
    return  static_cast<double>(converged() * 2);
  }else if (Umatch(x, "ids ")) { untested();
    return ids;
  }else if (Umatch(x, "idsxxx ")) {
    return idsxxx;
  }else if (Umatch(x, "gds ")) {
    return gds;
  }else if (Umatch(x, "gmf ")) { untested();
    return gmf;
  }else if (Umatch(x, "gmr ")) { untested();
    return gmr;
  }else if (Umatch(x, "gmbf ")) { untested();
    return gmbf;
  }else if (Umatch(x, "gmbr ")) { untested();
    return gmbr;
  }else if (Umatch(x, "idb ")) {
    return idb;
  }else if (Umatch(x, "idbxxx ")) { untested();
    return idbxxx;
  }else if (Umatch(x, "gdbdb ")) { untested();
    return gdbdb;
  }else if (Umatch(x, "gdbds ")) {
    return gdbds;
  }else if (Umatch(x, "gdbgs ")) {
    return gdbgs;
  }else if (Umatch(x, "gdbbs ")) {
    return gdbbs;
  }else if (Umatch(x, "isb ")) {
    return isb;
  }else if (Umatch(x, "isbxxx ")) { untested();
    return isbxxx;
  }else if (Umatch(x, "gsbsb ")) { untested();
    return gsbsb;
  }else if (Umatch(x, "gsbsd ")) {
    return gsbsd;
  }else if (Umatch(x, "gsbgd ")) {
    return gsbgd;
  }else if (Umatch(x, "gsbbd ")) {
    return gsbbd;
  }else if (Umatch(x, "qgate ")) { untested();
    return qgate;
  }else if (Umatch(x, "cgs ")) { untested();
    return cgs;
  }else if (Umatch(x, "cggb ")) { untested();
    return cggb;
  }else if (Umatch(x, "cgsb ")) { untested();
    return cgsb;
  }else if (Umatch(x, "cgdb ")) { untested();
    return cgdb;
  }else if (Umatch(x, "qgs ")) {
    return qgs;
  }else if (Umatch(x, "cgsgs ")) { untested();
    return cgsgs;
  }else if (Umatch(x, "cgsgb ")) { untested();
    return cgsgb;
  }else if (Umatch(x, "cgssb ")) { untested();
    return cgssb;
  }else if (Umatch(x, "cgsdb ")) { untested();
    return cgsdb;
  }else if (Umatch(x, "qgd ")) {
    return qgd;
  }else if (Umatch(x, "cgdgd ")) { untested();
    return cgdgd;
  }else if (Umatch(x, "cgdgb ")) { untested();
    return cgdgb;
  }else if (Umatch(x, "cgdsb ")) { untested();
    return cgdsb;
  }else if (Umatch(x, "cgddb ")) { untested();
    return cgddb;
  }else if (Umatch(x, "qdrn ")) { untested();
    return qdrn;
  }else if (Umatch(x, "cdsds ")) { untested();
    return cdsds;
  }else if (Umatch(x, "cdgb ")) { untested();
    return cdgb;
  }else if (Umatch(x, "cdsb ")) { untested();
    return cdsb;
  }else if (Umatch(x, "cddb ")) { untested();
    return cddb;
  }else if (Umatch(x, "qbulk ")) { untested();
    return qbulk;
  }else if (Umatch(x, "cbs ")) { untested();
    return cbs;
  }else if (Umatch(x, "cbgb ")) { untested();
    return cbgb;
  }else if (Umatch(x, "cbsb ")) { untested();
    return cbsb;
  }else if (Umatch(x, "cbdb ")) { untested();
    return cbdb;
  }else if (Umatch(x, "qbs ")) {
    return qbs;
  }else if (Umatch(x, "cbsbs ")) { untested();
    return cbsbs;
  }else if (Umatch(x, "cbsgb ")) { untested();
    return cbsgb;
  }else if (Umatch(x, "cbssb ")) { untested();
    return cbssb;
  }else if (Umatch(x, "cbsdb ")) { untested();
    return cbsdb;
  }else if (Umatch(x, "qbd ")) {
    return qbd;
  }else if (Umatch(x, "cbdbd ")) { untested();
    return cbdbd;
  }else if (Umatch(x, "cbdgb ")) { untested();
    return cbdgb;
  }else if (Umatch(x, "cbdsb ")) { untested();
    return cbdsb;
  }else if (Umatch(x, "cbddb ")) { untested();
    return cbddb;
  }else if (Umatch(x, "gtau ")) { untested();
    return gtau;
  }else if (Umatch(x, "cqgb ")) {
    return cqgb;
  }else if (Umatch(x, "cqsb ")) { untested();
    return cqsb;
  }else if (Umatch(x, "cqdb ")) { untested();
    return cqdb;
  }else if (Umatch(x, "cqbb ")) { untested();
    return cqbb;
  }else if (Umatch(x, "vgs ")) { untested();
    return vgs;
  }else if (Umatch(x, "vds ")) { untested();
    return vds;
  }else if (Umatch(x, "vbs ")) { untested();
    return vbs;
  }else if (Umatch(x, "vdsat ")) {
    return vdsat;
  }else if (Umatch(x, "vgst ")) {
    return vgst;
  }else if (Umatch(x, "von ")) {
    return von;
  }else if (Umatch(x, "reversed ")) { untested();
    return reversed;
  }else if (Umatch(x, "cutoff ")) { untested();
    return cutoff;
  }else if (Umatch(x, "subthreshold ")) { untested();
    return subthreshold;
  }else if (Umatch(x, "saturated ")) { untested();
    return saturated;
  }else if (Umatch(x, "sbfwd ")) { untested();
    return sbfwd;
  }else if (Umatch(x, "punchthru ")) { untested();
    return punchthru;
  }else if (Umatch(x, "w ")) {
    return c->w_in;
  }else if (Umatch(x, "l ")) {
    return c->l_in;
  }else {
    return BASE_SUBCKT::tr_probe_num(x);
  }
}
/*--------------------------------------------------------------------------*/
XPROBE DEV_BUILT_IN_MOS::ac_probe_ext(const std::string& x)const
{ itested();
  const ADP_BUILT_IN_MOS* a = prechecked_cast<const ADP_BUILT_IN_MOS*>(adp());
  assert(a);
  if (Umatch(x, "dvth ")) { untested();
    assert(a);
    return XPROBE(a->delta_vth);
  }
  return BASE_SUBCKT::ac_probe_ext(x);
}
/*--------------------------------------------------------------------------*/
double DEV_BUILT_IN_MOS::tt_probe_num(const std::string& x)const
{
  assert(_n);
  const COMMON_BUILT_IN_MOS* c = prechecked_cast<const COMMON_BUILT_IN_MOS*>(common());
  assert(c);
  const MODEL_BUILT_IN_MOS_BASE* m = prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(c->model());
  assert(m);
  const SDP_BUILT_IN_MOS_BASE* s = prechecked_cast<const SDP_BUILT_IN_MOS_BASE*>(c->sdp());
  assert(s); USE(s);

  if (Umatch(x, "use_bti ")) { itested();
    return  m->use_bti();
  }else if (Umatch(x, "stress ")) { untested();
    return 19999;
  }
  return BASE_SUBCKT::tt_probe_num(x);
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_MOS::reverse_if_needed()
{
    if (vds < 0) {
      error(bTRACE, long_label() + ": reversing\n");
      error(bTRACE, "before: vds=%g vgs=%g vbs=%g\n", vds, vgs, vbs);
      reversed = !reversed;
      vgs -= vds;
      vbs -= vds;
      vds = -vds;
      error(bTRACE, "after: vds=%g vgs=%g vbs=%g\n", vds, vgs, vbs);
      if (OPT::dampstrategy & dsREVERSE) { untested();
	_sim->_fulldamp = true;
	untested();
	error(bTRACE, long_label() + ":reverse damp\n");
      }
      if (!(OPT::mosflags & 0040)) {
	vbs = std::min(vbs,(hp_float_t)0.);
      }else{ untested();
	untested();
      }
    }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
bool DEV_BUILT_IN_MOS::tr_needs_eval()const
{
  if (is_q_for_eval()) { untested();
    return false;
  }else if (!converged()) {
    return true;
  }else{
    const COMMON_BUILT_IN_MOS* c = prechecked_cast<const COMMON_BUILT_IN_MOS*>(common());
    assert(c);
    const MODEL_BUILT_IN_MOS_BASE* m=prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(c->model());
    assert(m);
    polarity_t polarity = m->polarity;
    node_t& eff_s((reversed) ? _n[n_id] : _n[n_is]);
    node_t& eff_d((reversed) ? _n[n_is] : _n[n_id]);
    if (_BTI) {
      if ( _BTI->tr_needs_eval()) { untested();
	return true;
      }
    }
    // HCI?
    return !(conchk(vds,polarity*(eff_d.v0()-eff_s.v0()),OPT::vntol)
	     && conchk(vgs, polarity*(_n[n_g].v0()-eff_s.v0()),
		       OPT::vntol)
	     && conchk(vbs, polarity*(_n[n_b].v0()-eff_s.v0()),
		       OPT::vntol));
  }
}
/*--------------------------------------------------------------------------*/
bool DEV_BUILT_IN_MOS::do_tr()
{
  const COMMON_BUILT_IN_MOS* c = prechecked_cast<const COMMON_BUILT_IN_MOS*>(common());
  assert(c);
  const MODEL_BUILT_IN_MOS_BASE* m = prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(c->model());
  assert(m);

  bool was_cutoff = cutoff;
  bool was_subthreshold = subthreshold;
  bool was_saturated = saturated;
  bool was_reversed = reversed;
  bool was_sbfwd = sbfwd;
  polarity_t polarity = m->polarity;

  if (adp()) {
    if (_sim->_age) {
      adp()->apply(this);
    }else if (_sim->analysis_is_dcop()) {
      adp()->apply(this);
    }else{
    }
  }else{
  }

  if (_sim->is_initial_step()) {
    reversed = false;
    vds = vgs = vbs = 0.;
  }else{
    double Vds, Vgs, Vbs;
    if (reversed) {
      Vds = polarity * volts_limited(_n[n_is],_n[n_id]);
      Vgs = polarity * volts_limited(_n[n_g],_n[n_id]);
      Vbs = polarity * volts_limited(_n[n_b],_n[n_id]);
    }else{
      Vds = polarity * volts_limited(_n[n_id],_n[n_is]);
      Vgs = polarity * volts_limited(_n[n_g],_n[n_is]);
      Vbs = polarity * volts_limited(_n[n_b],_n[n_is]);
    }
    vgs = fet_limit_vgs(Vgs, vgs, von);
    if (_n[n_d].n_() == _n[n_g].n_()) {
      vds = Vds + (vgs - Vgs);
    }else{
      // Spice hacks Vds here, but my tests show that it often makes
      // convergence worse, and never improves it.
      // I am guessing that it does help when drain and gate are connected,
      // and Spice does it here in case they are and cannot be determined
      // whether they are or not.
      // The hack maintains Vdg after Vgs limiting.
      //Vds = Vds + (vgs - Vgs);
      vds = fet_limit_vds(Vds, vds);
    }
    vbs = std::min(Vbs, 0.);
    //vbs = pnj_limit(double Vbs, double vbs, double vt, double vcrit);
    //vds = Vds;
    //vgs = Vgs;
    //vbs = Vbs;
  }

  assert(qgate == qgate);
  assert(qgs == qgs);
  assert(qgd == qgd);
  assert(qdrn == qdrn);
  assert(qbulk == qbulk);
  assert(qbs == qbs);
  assert(qbd == qbd);

  m->tr_eval(this);

  assert(qgate == qgate);
  assert(qgs == qgs);
  assert(qgd == qgd);
  assert(qdrn == qdrn);
  assert(qbulk == qbulk);
  assert(qbs == qbs);
  assert(qbd == qbd);

  if (reversed) {
    idsxxx = ids + vds*gds + vgs*gmr + vbs*gmbr;
    isbxxx = isb - vds*gsbsd - vgs*gsbgd - vbs*gsbbd;
    idbxxx = 0.;
  }else{
    idsxxx = ids - vds*gds - vgs*gmf - vbs*gmbf;
    idbxxx = idb - vds*gdbds - vgs*gdbgs - vbs*gdbbs;
    isbxxx = 0.;
  }
  ids *= polarity;
  idsxxx *= polarity;
  assert(subckt());

#ifdef BTI_LATE_EVAL
  // ~ SUBCKT
  bool isconverged = true;
  if (OPT::bypass) {
    for (std::list<CARD*>::iterator ci=subckt()->begin(); ci!=subckt()->end(); ++ci) {
      if ((**ci).tr_needs_eval()) {
	if (*ci !=_BTI ) isconverged &= (**ci).do_tr();
      }else{
      }
    }
  }else{
    for (std::list<CARD*>::iterator ci=subckt()->begin(); ci!=subckt()->end(); ++ci) {
      if (*ci !=_BTI )isconverged &= (**ci).do_tr();
    }
  }
//=============
  if( m->use_bti() ){
    if(  isconverged ){
      //std::cout << "* btieval " << _sim->iteration_number() << " \n";
//      _BTI->do_tr();
    } else {
      //std::cout << "* not btieval " << _sim->iteration_number() << " \n";
    }
  }
  set_converged(isconverged);
#else
  set_converged(subckt()->do_tr());
#endif

  // BTI_HACK
#ifndef BTI_IN_SUBCKT
 // if( m->use_bti &&  converged() )
   // set_converged(_BTI->do_tr());
#endif
  
  trace3(long_label().c_str(), vds, vgs, vbs);
  trace4("", ids, gmf, gds, gmbf);
  trace4("", ids, gmr, gds, gmbr);
  if (was_cutoff != cutoff  ||  was_subthreshold != subthreshold  
  	||  was_saturated != saturated  ||  was_reversed != reversed  
	||  was_sbfwd != sbfwd) {
    if (OPT::dampstrategy & dsDEVREGION) {
      _sim->_fulldamp = true;
    }else{
    }
    trace1("region change", long_label());
  }else{
  }

#ifdef BTI_IN_SUBCKT
  return converged();
#else
  return (converged() && _BTI->converged());
  _BTI->q_accept();
#endif
//  q_accept();
//
// necessary? adp should care for itself.
//  adp()->q_accept();
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_MOS::do_sens()
{
  const COMMON_BUILT_IN_MOS* c = prechecked_cast<const COMMON_BUILT_IN_MOS*>(common());
  assert(c);
  const MODEL_BUILT_IN_MOS_BASE* m = prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(c->model());
  assert(m);

  m->sens_eval(this);
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_MOS::tt_begin()
{
  BASE_SUBCKT::tt_begin();

  const COMMON_BUILT_IN_MOS* c = (const COMMON_BUILT_IN_MOS*)(common());
  assert(c);
  assert(c->model());

  const MODEL_BUILT_IN_MOS_BASE* m = (const MODEL_BUILT_IN_MOS_BASE*)(c->model());
  assert(m);

  adp()->tt_begin();

  m->do_tt_prepare(this);
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_MOS::tr_stress_last()
{
  q_tt_accept(); // acceptq ist lifo
//  BASE_SUBCKT::tr_stress_last();
  if (_HCI) {
    _HCI->tr_stress_last(); // calls q_tt_accept
  }
  if (_BTI) {
    _BTI->tr_stress_last(); // calls q_tt_accept
  }
  if (adp()) {
//    adp()->apply(this);
  }
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_MOS::tt_accept()
{
  if (_HCI) {
//    _HCI->tt_accept();
  }
  if (_BTI) {
//    _BTI->tt_accept();
  }
  if (adp()) {
    adp()->apply(this);
  }
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_MOS::tt_advance()
{
  if (adp()) {
    adp()->tt_advance();
  } else { untested();
  }

  // tt_advance also resets time[] in elements.
  BASE_SUBCKT::tt_advance();
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_MOS::tr_save_amps(int trstep)
{ untested();
  trace1(("DEV_BUILT_IN_MOS::tr_save_amps()" + long_label()).c_str(),trstep);
  const COMMON_COMPONENT* cc = common();
  const MODEL_BUILT_IN_MOS_BASE* m = prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(cc->model());
  assert(m); USE(m);
  int j = 3;
  int k = j;

  incomplete();

  hp_float_t tramps0= tr_probe_num("id");
  trace2(("DEV_BUILT_IN_MOS::tr_save_amps " + long_label()).c_str(),trstep, tramps0);
  hp_float_t tramps1=  ((ELEMENT*)_Cgb)->tr_amps() 
                     +((_Cgd) ? ((ELEMENT*)_Cgd)->tr_amps() : 0 )
                     +((_Cgs) ? ((ELEMENT*)_Cgs)->tr_amps() : 0 )                   ;  
  trace2(("DEV_BUILT_IN_MOS::tr_save_amps " + long_label()).c_str(),trstep, tramps1);

  // bloeder bug
  hp_float_t tramps2 = 0;
  //hp_float_t tramps2=  (_Rs) ? ((ELEMENT*)_Rs )->tr_amps() :
  //                           - ((ELEMENT*)_Ids)->tr_amps()
  //                           - ((ELEMENT*)_Cgs)->tr_amps()
  //                           - ((ELEMENT*)_Dsb)->tr_amps() * m->polarity;

  hp_float_t tramps3= tramps0 + tramps1 ;
  trace2(("DEV_BUILT_IN_MOS::tr_save_amps" + long_label()).c_str(),trstep, tramps3);
  //
  // double tramps0= tr_probe_num("id");
  // double tramps1=0; // tr_probe_num("ig");
  // double tramps2=0; // tr_probe_num("is");
  // double tramps3= tramps0 + tramps1 + tramps2;


  //std::cerr << "saving _amps[ " << n << " ]. have " << net_nodes() << " nodes\n";


  if (_amps==0) _tr_amps_diff_cur = 0;
  _tr_amps_diff_max = 0; // maximum der delta i in diesem zeitschritt.
  _tr_amps_scale_max = 0; 


  if(_amps!=0) { untested();
    double _amps3 = _amps[trstep*k + 0] + _amps[trstep*k + 1] + _amps[trstep*k + 2] ;

    double diff0= _amps[trstep*k + 0] - tramps0;
    double diff1= _amps[trstep*k + 1] - tramps1;
    double diff2= _amps[trstep*k + 2] - tramps2;
    double diff3= _amps3              - tramps3;

    double sum0= fabs(_amps[trstep*k + 0]) + fabs(tramps0);
    double sum1= fabs(_amps[trstep*k + 1]) + fabs(tramps1);
    double sum2= fabs(_amps[trstep*k + 2]) + fabs(tramps2);
    double sum3= fabs(_amps3)              + fabs(tramps3);

    _tr_amps_diff_max =  fabs(diff0);
    _tr_amps_diff_max =  fmax( _tr_amps_diff_max, fabs(diff1));
    _tr_amps_diff_max =  fmax( _tr_amps_diff_max, fabs(diff2));
    _tr_amps_diff_max =  fmax( _tr_amps_diff_max, fabs(diff0 + diff1 + diff2));

    _tr_amps_scale_max =                           fabs(diff0/sum0) ;
    _tr_amps_scale_max =  fmax( _tr_amps_diff_max, fabs(diff1/sum1) ) ;
    _tr_amps_scale_max =  fmax( _tr_amps_diff_max, fabs(diff2/sum2) ) ;
    _tr_amps_scale_max =  fmax( _tr_amps_diff_max, fabs(diff3/sum3) ) ;
  }


 // _amps_new[ trstep*k+0] = tramps0;
 // _amps_new[ trstep*k+1] = tramps1;
 // _amps_new[ trstep*k+2] = tramps2;

  tr_behaviour_del = _tr_amps_diff_max;
  tr_behaviour_rel = _tr_amps_scale_max;

  assert( tr_behaviour_del >= 0 );

  CKT_BASE::tr_behaviour_del = fmax( CKT_BASE::tr_behaviour_del, tr_behaviour_del);
  CKT_BASE::tr_behaviour_rel = fmax( CKT_BASE::tr_behaviour_rel, tr_behaviour_rel);
//  trace10 << "DEV_BUILT_IN_MOS::tr_save_amps " << trstep << short_label() << 
//    " " << tr_behaviour_del<< "rel: " << tr_behaviour_rel << "\n";

  trace1(("DEV_BUILT_IN_MOS::tr_save_amps tr_behaviour" + long_label()).c_str(),trstep);

  tr_behaviour();

  trace1(("DEV_BUILT_IN_MOS::tr_save_amps done " + long_label()).c_str(),trstep);
}
/*--------------------------------------------------------------------------*/
void      DEV_BUILT_IN_MOS::dc_advance() {set_not_converged(); BASE_SUBCKT::dc_advance();
#ifndef BTI_IN_SUBCKT
  const COMMON_COMPONENT* cc = common();
  MODEL_BUILT_IN_MOS_BASE* m = ( MODEL_BUILT_IN_MOS_BASE*)(cc->model());
  if(_BTI) { untested();
    _BTI->dc_advance();
  }
#endif
  
  }
  void      DEV_BUILT_IN_MOS::tr_advance() {set_not_converged(); BASE_SUBCKT::tr_advance();}
  void      DEV_BUILT_IN_MOS::tr_regress() {set_not_converged(); BASE_SUBCKT::tr_regress();}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_MOS::precalc_last()
{
  COMPONENT::precalc_last();
  assert(subckt());
  if(subckt()) {
    subckt()->precalc_last();
  }
#ifndef BTI_IN_SUBCKT
  trace0("DEV_BUILT_IN_MOS::precalc_last bti notin sckt");
  const COMMON_COMPONENT* cc = common();
  MODEL_BUILT_IN_MOS_BASE* m = ( MODEL_BUILT_IN_MOS_BASE*)(cc->model());
  //if(m->use_bti())
  //  _BTI->precalc_last();
#endif
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_MOS::precalc_first()
{
  COMPONENT::precalc_first();
  if(subckt()) {
    subckt()->precalc_first();
  }
#ifndef BTI_IN_SUBCKT
  trace0("DEV_BUILT_IN_MOS::precalc_first bti notin sckt");
  // const COMMON_COMPONENT* cc = common();
  //MODEL_BUILT_IN_MOS_BASE* m = ( MODEL_BUILT_IN_MOS_BASE*)(cc->model());
  //if(m->use_bti())
  //  _BTI->precalc_first();
#endif
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_MOS::map_nodes()
{
  trace0("DEV_BUILT_IN_MOS::map_nodes");
  BASE_SUBCKT::map_nodes();
  const COMMON_BUILT_IN_MOS* c = static_cast<const COMMON_BUILT_IN_MOS*>(common());
  assert(c); USE(c);
  assert(c->model());
#ifndef BTI_IN_SUBCKT
  trace0("DEV_BUILT_IN_MOS::map_nodes, handle BTI...");
  if (_BTI) { untested();
    _BTI->map_nodes();
  }
#endif
  if(_HCI) {
  //  adp()->map_nodes();
  } else {
  }
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_MOS::tr_begin()
{
  BASE_SUBCKT::tr_begin();
#ifndef BTI_IN_SUBCKT
  const COMMON_BUILT_IN_MOS* c = static_cast<const COMMON_BUILT_IN_MOS*>(common());
  assert(c);
  assert(c->model());
  const MODEL_BUILT_IN_MOS_BASE* m = prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(c->model());
  assert(m);
  if (_BTI) { untested();
    _BTI->tr_begin();
  }
#endif
#ifndef HCI_IN_SUBCKT
  if (_HCI) {
    _HCI->tr_begin();
  }
#endif
}
/*--------------------------------------------------------------------------*/
void    DEV_BUILT_IN_MOS::tr_restore(){
  BASE_SUBCKT::tr_restore();
#ifndef BTI_IN_SUBCKT
  const COMMON_BUILT_IN_MOS* c = static_cast<const COMMON_BUILT_IN_MOS*>(common());
  assert(c);
  assert(c->model());
  const MODEL_BUILT_IN_MOS_BASE* m = prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(c->model());
  assert(m);
  if(_BTI){ untested();
    _BTI->tr_restore();
  }
#endif
  if(_HCI){
    _HCI->tr_restore();
  }
}
/*-------------------------------------------------------*/
void DEV_BUILT_IN_MOS::tr_load()
{
  BASE_SUBCKT::tr_load();
#ifndef BTI_IN_SUBCKT
  const COMMON_BUILT_IN_MOS* c = static_cast<const COMMON_BUILT_IN_MOS*>(common());
  assert(c);
  assert(c->model());
  const MODEL_BUILT_IN_MOS_BASE* m = prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(c->model());
  assert(m);
  if (_BTI) { untested();
    _BTI->tr_load();
  }
#endif
}
/*-------------------------------------------------------*/
TIME_PAIR  DEV_BUILT_IN_MOS::tr_review()
{
  const COMMON_BUILT_IN_MOS* c = static_cast<const COMMON_BUILT_IN_MOS*>(common());
  assert(c);
  assert(c->model());
  const MODEL_BUILT_IN_MOS_BASE* m = prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(c->model());
  assert(m); USE(m);

  TIME_PAIR time_by;
  if (adp()) {
//    time_by = adp()->tr_review(); // not part of subckt...
  } else {
  }
  time_by.min(BASE_SUBCKT::tr_review());
#ifndef HCI_IN_SUBCKT
  if(_HCI){
    time_by.min(_HCI->tr_review());
  }
#endif
#ifndef BTI_IN_SUBCKT
  if (_BTI) { untested();
    time_by.min(_BTI->tr_review());
  }
#endif
  if (adp()) {
    // incomplete();
    if (_sim->analysis_is_dcop()) { untested();
      adp()->apply(this);
    }else{
    }
  }else{
  }
  return time_by;
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_MOS::tr_unload()
{
  BASE_SUBCKT::tr_unload();
#ifndef BTI_IN_SUBCKT
  const COMMON_BUILT_IN_MOS* c = static_cast<const COMMON_BUILT_IN_MOS*>(common());
  assert(c);
  assert(c->model());
  const MODEL_BUILT_IN_MOS_BASE* m = prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(c->model());
  assert(m);
  if (_BTI) { untested();
    _BTI->tr_unload();
  }
#endif
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_MOS::do_tt()
{
//  BASE_SUBCKT::do_tt();

  if (_BTI) {
    _BTI->do_tt();
  } else {
  }
  if (_HCI) {
    _HCI->do_tt();
  } else {
  }
  if (adp()) {
    adp()->apply(this);
  } else { untested();
  }
}
/*-------------------------------------------------------*/
void DEV_BUILT_IN_MOS::tr_accept()
{
  BASE_SUBCKT::tr_accept();

  const COMMON_COMPONENT* cc = common();
  const MODEL_BUILT_IN_MOS_BASE* m = asserted_cast<const MODEL_BUILT_IN_MOS_BASE*>(cc->model());
#ifndef BTI_IN_SUBCKT
  if (_BTI) { untested();
    assert(reversed);
    _BTI->tr_accept();
  } else { untested();
  }
#endif
  if (_BTI) {
    if (reversed && m->polarity<0) { itested();

    }
  }

#ifndef HCI_IN_SUBCKT
  if (_HCI) { untested();
    _HCI->tr_accept();
  }
#endif
}
/*-------------------------------------------------------*/
/*-------------------------------------------------------*/
XPROBE DEV_BUILT_IN_MOS::sens_probe_ext(const std::string& x)const
{
  // BUG delegate to subdevices!
  const COMMON_BUILT_IN_MOS* c = prechecked_cast<const COMMON_BUILT_IN_MOS*>(common());
  assert(c);
  const MODEL_BUILT_IN_MOS_BASE* m = asserted_cast<const MODEL_BUILT_IN_MOS_BASE*>(c->model());
  trace1("DEV_BUILT_IN_MOS::sens_probe_ext", x);
  unsigned n1 = _n[n_is].m_();
  unsigned n2 = _n[n_id].m_();
  COMPLEX a = CKT_BASE::_sim->_sens[n1];
  COMPLEX b = CKT_BASE::_sim->_sens[n2];
  double vn = CKT_BASE::_sim->vdc()[n1];
  double vp = CKT_BASE::_sim->vdc()[n2];
  double dv = vp - vn;

  COMPLEX gmsens = (a-b) * dv * vgs;
  if (Umatch(x, "gm ")) {
    return XPROBE( gmsens );
  } else if (Umatch(x, "w{idth} ")) {
    return XPROBE( gmsens * m->dgmdw_in(this) );
  } else if (Umatch(x, "l{ength} ")) {
    return XPROBE( gmsens * m->dgmdl_in(this) );
  } else { untested();
  }

  if (Umatch(x, "d1 ")) {	 // debug...
    return XPROBE( dv );
  }
  return BASE_SUBCKT::sens_probe_ext(x);
}
/*-------------------------------------------------------*/
/*-------------------------------------------------------*/
}
// * vim:ts=8:sw=2:noet:
