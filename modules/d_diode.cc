/* $Id: d_diode.cc,v 1.8 2010-07-16 14:49:39 felix Exp $ -*- C++ -*-
 * vim:et:sw=2:ts=8:
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
 * diode model.
 * netlist syntax:
 * device:  dxxxx n+ n- mname <area> <off> <ic=vd> <model-card-args>
 * model:   .model mname D <args>
 *
 * The section "eval Yj" is a big mess.
 * It will be redone using multiple files, like the MOS models.
 */
/* This file is automatically generated. DO NOT EDIT */

#include "e_aux.h"
#include "e_storag.h"
  static bool dummy=false;
  enum {USE_OPT = 0x8000};
#include "globals.h"
#include "e_elemnt.h"
#include "d_diode.h"
namespace UF{
/*--------------------------------------------------------------------------*/
using notstd::to_lower;
/*--------------------------------------------------------------------------*/
int DEV_BUILT_IN_DIODE::_count = -1;
int COMMON_BUILT_IN_DIODE::_count = -1;
static COMMON_BUILT_IN_DIODE Default_BUILT_IN_DIODE(CC_STATIC);
/*--------------------------------------------------------------------------*/
const double NA(NOT_INPUT);
const double INF(BIGBIG);
/*--------------------------------------------------------------------------*/
int MODEL_BUILT_IN_DIODE::_count = 0;
/*--------------------------------------------------------------------------*/
namespace MODEL_BUILT_IN_DIODE_DISPATCHER { 
  static DEV_BUILT_IN_DIODE p1d;
  static MODEL_BUILT_IN_DIODE p1(&p1d);
  static DISPATCHER<MODEL_CARD>::INSTALL
    d1(&model_dispatcher, "d", &p1);
}
/*--------------------------------------------------------------------------*/
void SDP_BUILT_IN_DIODE::init(const COMMON_COMPONENT* cc)
{
  assert(cc);
  SDP_CARD::init(cc);
}
/*--------------------------------------------------------------------------*/
TDP_BUILT_IN_DIODE::TDP_BUILT_IN_DIODE(const DEV_BUILT_IN_DIODE*)
{
}
/*--------------------------------------------------------------------------*/
MODEL_BUILT_IN_DIODE::MODEL_BUILT_IN_DIODE(const BASE_SUBCKT* p)
  :MODEL_CARD(p),
   js(1e-14),
   rs(0.0),
   n_factor(1.0),
   tt(0.0),
   cjo(NA),
   pb(NA),
   mj(0.5),
   eg(1.11),
   xti(3.0),
   kf(NA),
   af(NA),
   fc(0.5),
   bv(NA),
   ibv(1e-3),
   cjsw(0.0),
   pbsw(NA),
   mjsw(NA),
   gparallel(0.0),
   flags(int(USE_OPT)),
   mos_level(0)
{
  if (ENV::run_mode != rPRE_MAIN) {
    ++_count;
  }else{
  }
  set_default(&_tnom_c, OPT::tnom_c);
}
/*--------------------------------------------------------------------------*/
MODEL_BUILT_IN_DIODE::MODEL_BUILT_IN_DIODE(const MODEL_BUILT_IN_DIODE& p)
  :MODEL_CARD(p),
   js(p.js),
   rs(p.rs),
   n_factor(p.n_factor),
   tt(p.tt),
   cjo(p.cjo),
   pb(p.pb),
   mj(p.mj),
   eg(p.eg),
   xti(p.xti),
   kf(p.kf),
   af(p.af),
   fc(p.fc),
   bv(p.bv),
   ibv(p.ibv),
   cjsw(p.cjsw),
   pbsw(p.pbsw),
   mjsw(p.mjsw),
   gparallel(p.gparallel),
   flags(p.flags),
   mos_level(p.mos_level)
{
  if (ENV::run_mode != rPRE_MAIN) {
    ++_count;
  }else{untested();//194
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_DIODE::dev_type()const
{
  if (dummy == true) {
    return "d";
  }else{untested();//235
    return MODEL_CARD::dev_type();
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_DIODE::set_dev_type(const std::string& new_type)
{
  if (Umatch(new_type, "d ")) {
    dummy = true;
  }else{
    MODEL_CARD::set_dev_type(new_type);
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_DIODE::precalc_first()
{
    const CARD_LIST* par_scope = scope();
    assert(par_scope);
    MODEL_CARD::precalc_first();
    e_val(&(this->js), 1e-14, par_scope);
    e_val(&(this->rs), 0.0, par_scope);
    e_val(&(this->n_factor), 1.0, par_scope);
    e_val(&(this->tt), 0.0, par_scope);
    e_val(&(this->cjo), NA, par_scope);
    e_val(&(this->pb), NA, par_scope);
    e_val(&(this->mj), 0.5, par_scope);
    e_val(&(this->eg), 1.11, par_scope);
    e_val(&(this->xti), 3.0, par_scope);
    e_val(&(this->kf), NA, par_scope);
    e_val(&(this->af), NA, par_scope);
    e_val(&(this->fc), 0.5, par_scope);
    e_val(&(this->bv), NA, par_scope);
    e_val(&(this->ibv), 1e-3, par_scope);
    e_val(&(this->cjsw), 0.0, par_scope);
    e_val(&(this->pbsw), NA, par_scope);
    e_val(&(this->mjsw), NA, par_scope);
    e_val(&(this->gparallel), 0.0, par_scope);
    e_val(&(this->flags), int(USE_OPT), par_scope);
    e_val(&(this->mos_level), 0, par_scope);
    // final adjust: code_pre
    // final adjust: override
    // final adjust: raw
    e_val(&(this->js), 1e-14, par_scope);
    e_val(&(this->rs), 0.0, par_scope);
    e_val(&(this->n_factor), 1.0, par_scope);
    e_val(&(this->tt), 0.0, par_scope);
    e_val(&(this->cjo), 0.0, par_scope);
    e_val(&(this->pb), 1.0, par_scope);
    e_val(&(this->mj), 0.5, par_scope);
    e_val(&(this->eg), 1.11, par_scope);
    e_val(&(this->xti), 3.0, par_scope);
    e_val(&(this->kf), NA, par_scope);
    e_val(&(this->af), NA, par_scope);
    e_val(&(this->fc), 0.5, par_scope);
    e_val(&(this->bv), NA, par_scope);
    e_val(&(this->ibv), 1e-3, par_scope);
    e_val(&(this->cjsw), 0.0, par_scope);
    e_val(&(this->pbsw), pb, par_scope);
    e_val(&(this->mjsw), 0.33, par_scope);
    e_val(&(this->gparallel), 0.0, par_scope);
    e_val(&(this->flags), int(USE_OPT), par_scope);
    e_val(&(this->mos_level), 0, par_scope);
    // final adjust: mid
    // final adjust: calculated
    // final adjust: post

      if (bv == 0.) {
	bv = NA;
      }
    // final adjust: done
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_DIODE::precalc_last()
{
    MODEL_CARD::precalc_last();
}
/*--------------------------------------------------------------------------*/
SDP_CARD* MODEL_BUILT_IN_DIODE::new_sdp(COMMON_COMPONENT* c)const
{
  assert(c);
  if (COMMON_BUILT_IN_DIODE* cc = dynamic_cast<COMMON_BUILT_IN_DIODE*>(c)) {
    if (cc->_sdp) {
      cc->_sdp->init(cc);
      return cc->_sdp;
    }else{
      delete cc->_sdp;
      return new SDP_BUILT_IN_DIODE(c);
    }
  }else{
    return MODEL_CARD::new_sdp(c);
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_DIODE::set_param_by_index(int i, std::string& value, int offset)
{
  switch (MODEL_BUILT_IN_DIODE::param_count() - 1 - i) {
  case 0: untested(); break;
  case 1: _tnom_c = value; break;
  case 2: js = value; break;
  case 3: rs = value; break;
  case 4: n_factor = value; break;
  case 5: tt = value; break;
  case 6: cjo = value; break;
  case 7: pb = value; break;
  case 8: mj = value; break;
  case 9: eg = value; break;
  case 10: xti = value; break;
  case 11: kf = value; break;
  case 12: af = value; break;
  case 13: fc = value; break;
  case 14: bv = value; break;
  case 15: ibv = value; break;
  case 16: cjsw = value; break;
  case 17: pbsw = value; break;
  case 18: mjsw = value; break;
  case 19: gparallel = value; break;
  case 20: flags = value; break;
  case 21: mos_level = value; break;
  default: throw Exception_Too_Many(unsigned(i), 21u, offset); break;
  }
}
/*--------------------------------------------------------------------------*/
bool MODEL_BUILT_IN_DIODE::param_is_printable(int i)const
{
  switch (MODEL_BUILT_IN_DIODE::param_count() - 1 - i) {
  case 0:  return (false);
  case 1:  return (true);
  case 2:  return (true);
  case 3:  return (true);
  case 4:  return (true);
  case 5:  return (true);
  case 6:  return (true);
  case 7:  return (true);
  case 8:  return (true);
  case 9:  return (true);
  case 10:  return (true);
  case 11:  return (kf.has_hard_value());
  case 12:  return (af.has_hard_value());
  case 13:  return (true);
  case 14:  return (bv.has_hard_value());
  case 15:  return (has_good_value(bv));
  case 16:  return (cjsw != 0.);
  case 17:  return (cjsw != 0.);
  case 18:  return (cjsw != 0.);
  case 19:  return (gparallel != 0.);
  case 20:  return (!(flags & USE_OPT));
  case 21:  return (mos_level.has_hard_value());
  default: return false;
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_DIODE::param_name(int i)const
{
  switch (MODEL_BUILT_IN_DIODE::param_count() - 1 - i) {
  case 0:  return "=====";
  case 1:  return "tnom";
  case 2:  return "is";
  case 3:  return "rs";
  case 4:  return "n";
  case 5:  return "tt";
  case 6:  return "cjo";
  case 7:  return "pb";
  case 8:  return "mj";
  case 9:  return "egap";
  case 10:  return "xti";
  case 11:  return "kf";
  case 12:  return "af";
  case 13:  return "fc";
  case 14:  return "bv";
  case 15:  return "ibv";
  case 16:  return "cjsw";
  case 17:  return "pbsw";
  case 18:  return "mjsw";
  case 19:  return "gparallel";
  case 20:  return "flags";
  case 21:  return "mos_level";
  default: return "";
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_DIODE::param_name(int i, int j)const
{
  if (j == 0) {
    return param_name(i);
  }else if (j == 1) {
    switch (MODEL_BUILT_IN_DIODE::param_count() - 1 - i) {
    case 0:  return "";
    case 1:  return "";
    case 2:  return "";
    case 3:  return "";
    case 4:  return "";
    case 5:  return "";
    case 6:  return "";
    case 7:  return "vj";
    case 8:  return "m";
    case 9:  return "eg";
    case 10:  return "";
    case 11:  return "";
    case 12:  return "";
    case 13:  return "";
    case 14:  return "";
    case 15:  return "";
    case 16:  return "cjs";
    case 17:  return "pbs";
    case 18:  return "mjs";
    case 19:  return "gp";
    case 20:  return "";
    case 21:  return "";
    default: return "";
    }
  }else{
    return "";
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_DIODE::param_value(int i)const
{
  switch (MODEL_BUILT_IN_DIODE::param_count() - 1 - i) {
  case 0:  unreachable(); return "";
  case 1:  return _tnom_c.string();
  case 2:  return js.string();
  case 3:  return rs.string();
  case 4:  return n_factor.string();
  case 5:  return tt.string();
  case 6:  return cjo.string();
  case 7:  return pb.string();
  case 8:  return mj.string();
  case 9:  return eg.string();
  case 10:  return xti.string();
  case 11:  return kf.string();
  case 12:  return af.string();
  case 13:  return fc.string();
  case 14:  return bv.string();
  case 15:  return ibv.string();
  case 16:  return cjsw.string();
  case 17:  return pbsw.string();
  case 18:  return mjsw.string();
  case 19:  return gparallel.string();
  case 20:  return flags.string();
  case 21:  return mos_level.string();
  default: return "";
  }
}
/*--------------------------------------------------------------------------*/
bool MODEL_BUILT_IN_DIODE::is_valid(const COMPONENT* d)const
{
  assert(d);
  return MODEL_CARD::is_valid(d);
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_DIODE::tr_eval(COMPONENT*)const
{untested();//425
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
COMMON_BUILT_IN_DIODE::COMMON_BUILT_IN_DIODE(int c)
  :COMMON_COMPONENT(c),
   area(1.0),
   perim(0.0),
   off(false),
   ic(NA),
   is_raw(NA),
   rs_raw(NA),
   cj_raw(NA),
   cjsw_raw(NA),
   gparallel_raw(NA),
   _sdp(0),
   is_adjusted(NA),
   rs_adjusted(NA),
   cj_adjusted(NA),
   cjsw_adjusted(NA),
   gparallel_adjusted(NA)
{
  ++_count;
}
/*--------------------------------------------------------------------------*/
COMMON_BUILT_IN_DIODE::COMMON_BUILT_IN_DIODE(const COMMON_BUILT_IN_DIODE& p)
  :COMMON_COMPONENT(p),
   area(p.area),
   perim(p.perim),
   off(p.off),
   ic(p.ic),
   is_raw(p.is_raw),
   rs_raw(p.rs_raw),
   cj_raw(p.cj_raw),
   cjsw_raw(p.cjsw_raw),
   gparallel_raw(p.gparallel_raw),
   _sdp(0),
   is_adjusted(p.is_adjusted),
   rs_adjusted(p.rs_adjusted),
   cj_adjusted(p.cj_adjusted),
   cjsw_adjusted(p.cjsw_adjusted),
   gparallel_adjusted(p.gparallel_adjusted)
{
  ++_count;
}
/*--------------------------------------------------------------------------*/
COMMON_BUILT_IN_DIODE::~COMMON_BUILT_IN_DIODE()
{
  --_count;
  delete _sdp;
}
/*--------------------------------------------------------------------------*/
bool COMMON_BUILT_IN_DIODE::operator==(const COMMON_COMPONENT& x)const
{
  const COMMON_BUILT_IN_DIODE* p = dynamic_cast<const COMMON_BUILT_IN_DIODE*>(&x);
  return (p
    && area == p->area
    && perim == p->perim
    && off == p->off
    && ic == p->ic
    && is_raw == p->is_raw
    && rs_raw == p->rs_raw
    && cj_raw == p->cj_raw
    && cjsw_raw == p->cjsw_raw
    && gparallel_raw == p->gparallel_raw
    && _sdp == p->_sdp
    && COMMON_COMPONENT::operator==(x));
}
/*--------------------------------------------------------------------------*/
map<string, PARA_BASE COMMON_BUILT_IN_DIODE::*> COMMON_BUILT_IN_DIODE::param_dict
  = boost::assign::map_list_of
("Area", (PARA_BASE COMMON_BUILT_IN_DIODE::*)  (&COMMON_BUILT_IN_DIODE::area))
("Perim", (PARA_BASE COMMON_BUILT_IN_DIODE::*)  (&COMMON_BUILT_IN_DIODE::perim))
("OFF", (PARA_BASE COMMON_BUILT_IN_DIODE::*)  (&COMMON_BUILT_IN_DIODE::off))
("IC", (PARA_BASE COMMON_BUILT_IN_DIODE::*)  (&COMMON_BUILT_IN_DIODE::ic))
("IS", (PARA_BASE COMMON_BUILT_IN_DIODE::*)  (&COMMON_BUILT_IN_DIODE::is_raw))
("Rs", (PARA_BASE COMMON_BUILT_IN_DIODE::*)  (&COMMON_BUILT_IN_DIODE::rs_raw))
("Cjo", (PARA_BASE COMMON_BUILT_IN_DIODE::*)  (&COMMON_BUILT_IN_DIODE::cj_raw))
("CJSW", (PARA_BASE COMMON_BUILT_IN_DIODE::*)  (&COMMON_BUILT_IN_DIODE::cjsw_raw))
("GParallel", (PARA_BASE COMMON_BUILT_IN_DIODE::*)  (&COMMON_BUILT_IN_DIODE::gparallel_raw))
;map<string, PARA_BASE COMMON_BUILT_IN_DIODE::*> COMMON_BUILT_IN_DIODE::param_dict_low
  = boost::assign::map_list_of
("area", (PARA_BASE COMMON_BUILT_IN_DIODE::*)  (&COMMON_BUILT_IN_DIODE::area))
("perim", (PARA_BASE COMMON_BUILT_IN_DIODE::*)  (&COMMON_BUILT_IN_DIODE::perim))
("off", (PARA_BASE COMMON_BUILT_IN_DIODE::*)  (&COMMON_BUILT_IN_DIODE::off))
("ic", (PARA_BASE COMMON_BUILT_IN_DIODE::*)  (&COMMON_BUILT_IN_DIODE::ic))
("is", (PARA_BASE COMMON_BUILT_IN_DIODE::*)  (&COMMON_BUILT_IN_DIODE::is_raw))
("rs", (PARA_BASE COMMON_BUILT_IN_DIODE::*)  (&COMMON_BUILT_IN_DIODE::rs_raw))
("cjo", (PARA_BASE COMMON_BUILT_IN_DIODE::*)  (&COMMON_BUILT_IN_DIODE::cj_raw))
("cjsw", (PARA_BASE COMMON_BUILT_IN_DIODE::*)  (&COMMON_BUILT_IN_DIODE::cjsw_raw))
("gparallel", (PARA_BASE COMMON_BUILT_IN_DIODE::*)  (&COMMON_BUILT_IN_DIODE::gparallel_raw))
;
/*--------------------------------------------------------------------------*/
void COMMON_BUILT_IN_DIODE::set_param_by_name(string Name, string Value)
{
  PARA_BASE COMMON_BUILT_IN_DIODE::* x = (OPT::case_insensitive)?
     (param_dict_low[to_lower(Name)]) : (param_dict[Name]);
  if(x) { PARA_BASE* p = &(this->*x); *p = Value; return; }
  COMMON_COMPONENT::set_param_by_name(Name, Value);
}
/*--------------------------------------------------------------------------*/
void COMMON_BUILT_IN_DIODE::set_param_by_index(int i, std::string& Value, int Offset)
{
  switch (COMMON_BUILT_IN_DIODE::param_count() - 1 - i) {
  case 0:  area = Value; break;
  case 1:  perim = Value; break;
  case 2:  off = Value; break;
  case 3:  ic = Value; break;
  case 4:  is_raw = Value; break;
  case 5:  rs_raw = Value; break;
  case 6:  cj_raw = Value; break;
  case 7:  cjsw_raw = Value; break;
  case 8:  gparallel_raw = Value; break;
  default: COMMON_COMPONENT::set_param_by_index(i, Value, Offset);
  }
}
/*--------------------------------------------------------------------------*/
bool COMMON_BUILT_IN_DIODE::param_is_printable(int i)const
{
  switch (COMMON_BUILT_IN_DIODE::param_count() - 1 - i) {
  case 0:  return (true);
  case 1:  return (perim != 0.);
  case 2:  return (off);
  case 3:  return (has_good_value(ic));
  case 4:  return (has_good_value(is_raw));
  case 5:  return (has_good_value(rs_raw));
  case 6:  return (has_good_value(cj_raw));
  case 7:  return (has_good_value(cjsw_raw));
  case 8:  return (has_good_value(gparallel_raw));
  default: return COMMON_COMPONENT::param_is_printable(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string COMMON_BUILT_IN_DIODE::param_name(int i)const
{
  switch (COMMON_BUILT_IN_DIODE::param_count() - 1 - i) {
  case 0:  return "area";
  case 1:  return "perim";
  case 2:  return "off";
  case 3:  return "ic";
  case 4:  return "is";
  case 5:  return "rs";
  case 6:  return "cjo";
  case 7:  return "cjsw";
  case 8:  return "gparallel";
  default: return COMMON_COMPONENT::param_name(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string COMMON_BUILT_IN_DIODE::param_name(int i, int j)const
{
  if (j == 0) {
    return param_name(i);
  }else if (j == 1) {
    switch (COMMON_BUILT_IN_DIODE::param_count() - 1 - i) {
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
std::string COMMON_BUILT_IN_DIODE::param_value(int i)const
{
  switch (COMMON_BUILT_IN_DIODE::param_count() - 1 - i) {
  case 0:  return area.string();
  case 1:  return perim.string();
  case 2:  return off.string();
  case 3:  return ic.string();
  case 4:  return is_raw.string();
  case 5:  return rs_raw.string();
  case 6:  return cj_raw.string();
  case 7:  return cjsw_raw.string();
  case 8:  return gparallel_raw.string();
  default: return COMMON_COMPONENT::param_value(i);
  }
}
/*--------------------------------------------------------------------------*/
void COMMON_BUILT_IN_DIODE::expand(const COMPONENT* d)
{
  COMMON_COMPONENT::expand(d);
  attach_model(d);
  COMMON_BUILT_IN_DIODE* c = this; USE(c);
  const MODEL_BUILT_IN_DIODE* m = dynamic_cast<const MODEL_BUILT_IN_DIODE*>(model());
  if (!m) {
    throw Exception_Model_Type_Mismatch(d->long_label(), modelname(), "diode");
  }else{
  }
  // size dependent
  //delete _sdp;
  _sdp = m->new_sdp(this);
  assert(_sdp);
  const SDP_BUILT_IN_DIODE* s = prechecked_cast<const SDP_BUILT_IN_DIODE*>(_sdp);
  assert(s); USE(s);

  // subcircuit commons, recursive
  assert(c == this);
}
/*--------------------------------------------------------------------------*/
void COMMON_BUILT_IN_DIODE::precalc_first(const CARD_LIST* par_scope)
{
  assert(par_scope);
  COMMON_COMPONENT::precalc_first(par_scope);
    e_val(&(this->area), 1.0, par_scope);
    e_val(&(this->perim), 0.0, par_scope);
    e_val(&(this->off), false, par_scope);
    e_val(&(this->ic), NA, par_scope);
    e_val(&(this->is_raw), NA, par_scope);
    e_val(&(this->rs_raw), NA, par_scope);
    e_val(&(this->cj_raw), NA, par_scope);
    e_val(&(this->cjsw_raw), NA, par_scope);
    e_val(&(this->gparallel_raw), NA, par_scope);
}
/*--------------------------------------------------------------------------*/
void COMMON_BUILT_IN_DIODE::precalc_last(const CARD_LIST* par_scope)
{
  assert(par_scope);
  COMMON_COMPONENT::precalc_last(par_scope);
  COMMON_BUILT_IN_DIODE* c = this;
  const MODEL_BUILT_IN_DIODE* m = prechecked_cast<const MODEL_BUILT_IN_DIODE*>(model());
    // final adjust: code_pre
    // final adjust: override
    // final adjust: raw
    e_val(&(this->area), 1.0, par_scope);
    e_val(&(this->perim), 0.0, par_scope);
    e_val(&(this->off), false, par_scope);
    e_val(&(this->ic), NA, par_scope);
    e_val(&(this->is_raw), NA, par_scope);
    e_val(&(this->rs_raw), NA, par_scope);
    e_val(&(this->cj_raw), NA, par_scope);
    e_val(&(this->cjsw_raw), NA, par_scope);
    e_val(&(this->gparallel_raw), NA, par_scope);
    // final adjust: mid
    // final adjust: calculated
    is_adjusted = ((!has_good_value(c->is_raw))?(m->js*c->area):(c->is_raw));
    rs_adjusted = ((!has_good_value(c->rs_raw))
		? (m->rs / (c->area+1e-20)) : (c->rs_raw));
    cj_adjusted = ((!has_good_value(c->cj_raw))?(m->cjo*c->area):(c->cj_raw));
    cjsw_adjusted = ((!has_good_value(c->cjsw_raw))
		? (m->cjsw * c->perim) : (c->cjsw_raw));
    gparallel_adjusted = ((!has_good_value(c->gparallel_raw))
		? (m->gparallel*c->area) : (c->gparallel_raw));
    // final adjust: post
    // final adjust: done

  // size dependent
  //delete _sdp;
  _sdp = m->new_sdp(this);
  assert(_sdp);
  const SDP_BUILT_IN_DIODE* s = prechecked_cast<const SDP_BUILT_IN_DIODE*>(_sdp);
  assert(s); USE(s);

  // subcircuit commons, recursive
  assert(c == this);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
namespace DEV_BUILT_IN_DIODE_DISPATCHER { 
  static DEV_BUILT_IN_DIODE p0;
  static DISPATCHER<CARD>::INSTALL
    d0(&device_dispatcher, "D|diode", &p0);
}
/*--------------------------------------------------------------------------*/
static EVAL_BUILT_IN_DIODE_Cj Eval_Cj(CC_STATIC);
void EVAL_BUILT_IN_DIODE_Cj::tr_eval(ELEMENT* d)const
{
  assert(d);
  DEV_BUILT_IN_DIODE* p = prechecked_cast<DEV_BUILT_IN_DIODE*>(d->owner());
  assert(p);
  const COMMON_BUILT_IN_DIODE* c = prechecked_cast<const COMMON_BUILT_IN_DIODE*>(p->common());
  assert(c);
  const SDP_BUILT_IN_DIODE* s = prechecked_cast<const SDP_BUILT_IN_DIODE*>(c->sdp());
  assert(s); USE(s);
  const MODEL_BUILT_IN_DIODE* m = prechecked_cast<const MODEL_BUILT_IN_DIODE*>(c->model());
  assert(m);

    hp_float_t& volts = d->_y[0].x;
    trace1(d->long_label().c_str(), volts);
    
    double cb;
    if (c->cj_adjusted != 0.) {
      if (volts < m->fc * m->pb) {
	cb = c->cj_adjusted / pow(1. - (volts / m->pb),  m->mj);
      }else{
	cb = (c->cj_adjusted / pow(1. - m->fc, 1. + m->mj))
	  * (1. - m->fc*(1.+m->mj) + (volts/m->pb)*m->mj);
      }
    }else{
      cb = 0.;
    }
    assert(cb >= 0.);
    
    double csw;
    if (c->cjsw_adjusted != 0.) {
      if (volts < m->fc * m->pbsw) {
	csw = c->cjsw_adjusted / pow(1. - (volts / m->pbsw),  m->mjsw);
      }else{
	csw = (c->cjsw_adjusted / pow(1. - m->fc, 1. + m->mjsw))
	  * (1. - m->fc*(1.+m->mjsw) + (volts/m->pbsw)*m->mjsw);
      }
    }else{
      csw = 0.;
    }
    assert(csw >= 0.);
    
    double ctt;
    if (m->tt != 0.) {
      ctt = p->_gd * m->tt;
    }else{
      ctt = 0.;
    }
    assert(ctt >= 0.);
    
    trace4("", cb, csw, ctt, cb+csw+ctt);
    d->_y[0].f1 = cb + csw + ctt;
    if (d->_sim->analysis_is_tran_dynamic()) {
      const STORAGE* dd = prechecked_cast<const STORAGE*>(d);
      assert(dd);
      double cap = (d->_y[0].f1 + dd->_y[1].f1) / 2;
      d->_y[0].f0 = (d->_y[0].x - dd->_y[1].x) * cap + dd->_y[1].f0;
    }else{
      assert(d->_sim->analysis_is_static() || d->_sim->analysis_is_restore());
      d->_y[0].f0 = d->_y[0].x * d->_y[0].f1;
    }
    trace3(d->long_label().c_str(), d->_y[0].x, d->_y[0].f0, d->_y[0].f1);
}
/*--------------------------------------------------------------------------*/
static EVAL_BUILT_IN_DIODE_Yj Eval_Yj(CC_STATIC);
void EVAL_BUILT_IN_DIODE_Yj::tr_eval(ELEMENT* d)const
{
  assert(d);
  DEV_BUILT_IN_DIODE* p = prechecked_cast<DEV_BUILT_IN_DIODE*>(d->owner());
  assert(p);
  const COMMON_BUILT_IN_DIODE* c = prechecked_cast<const COMMON_BUILT_IN_DIODE*>(p->common());
  assert(c);
  const SDP_BUILT_IN_DIODE* s = prechecked_cast<const SDP_BUILT_IN_DIODE*>(c->sdp());
  assert(s); USE(s);
  const MODEL_BUILT_IN_DIODE* m = prechecked_cast<const MODEL_BUILT_IN_DIODE*>(c->model());
  assert(m);

    FPOLY1& y = d->_y[0];
    double volts = y.x;
    double amps  = y.f0;
    trace2(d->long_label().c_str(), volts, amps);
    
    int flags = (m->flags & USE_OPT) ? OPT::diodeflags : m->flags;
    double tempratio = (d->_sim->_temp_c+P_CELSIUS0) / (m->_tnom_c+P_CELSIUS0);
    double vt = P_K_Q * (d->_sim->_temp_c+P_CELSIUS0) * m->n_factor;
    region_t oldregion = p->_region;
    p->_isat = c->is_adjusted * pow(tempratio, m->xti)
      * exp((m->eg/vt) *(tempratio-1));
    trace4("", tempratio, vt, oldregion, p->_isat);
    
    if (m->mos_level > 0 || flags & 0040) { // Spice style limiting
      double vcrit = vt * log(vt / (M_SQRT2 * p->_isat));
      double vold = d->_y1.f0;
      if((volts > vcrit) && (std::abs(volts - vold) > (vt + vt))) {
	if(vold > 0) {
	  double arg = 1 + (volts - vold) / vt;
	  if(arg > 0) {
	    volts = vold + vt * log(arg);
	  }else{
	    volts = vcrit;
	  }
	}else{
	  volts = vt *log(volts/vt);
	}
      }else{
	// leave volts as is
      }
    }
    
    if (m->mos_level > 0) {
      switch (m->mos_level) {
      case 1:
      case 2:
      case 3:
      case 6:
      case 4:
      case 5:
	if (volts <= 0.) {
	  p->_region = REVERSE;
	  y.f1 = p->_isat / vt + OPT::gmin;
	  y.f0 = y.f1 * volts;
	}else{
	  p->_region = FORWARD;
	  double ev = exp(volts/vt);
	  y.f1 = p->_isat * ev / vt + OPT::gmin;
	  y.f0 = p->_isat * (ev - 1) + OPT::gmin * volts;
	}
	break;
      case 7:
      case 8:
	if (volts < .5) {
	  p->_region = REVERSE;
	  double ev = exp(volts/vt);
	  y.f1 = p->_isat * ev / vt + OPT::gmin;
	  y.f0 = p->_isat * (ev - 1) + OPT::gmin * volts;
	}else{
	  p->_region = FORWARD;
	  double ev = exp(.5/vt);
	  double t0 = p->_isat * ev / vt;
	  y.f1 = t0 + OPT::gmin;
	  y.f0 = p->_isat * (ev - 1) + t0 * (volts - .5) + OPT::gmin * volts;
	}
	break;
      default:
	unreachable();
	y.f1 = OPT::gmin;
	y.f0 = volts * y.f1;
      }
    }else if (flags & 0040) { // exact Spice model
      if (volts >= -3*vt) { // forward and weak reversed
	double evd = exp(volts/vt);
	y.f0 = p->_isat * (evd-1);
	y.f1 = p->_isat * evd/vt;
      }else if (has_good_value(m->bv) || volts >= m->bv) {
	double arg = 3 * vt / (volts * M_E); // strong reversed
	arg = arg * arg * arg;
	y.f0 = -p->_isat * (1+arg);
	y.f1 = p->_isat * 3 * arg / volts;
      }else{
	incomplete();
	double evrev = exp(-(m->bv+volts)/vt);
	y.f0 = -p->_isat * evrev;
	y.f1 = p->_isat * evrev / vt;
      }
      y.f0 += OPT::gmin * volts;
      y.f1 += OPT::gmin;
    }else{
       if (c->off  &&  d->_sim->is_initial_step()) { /*initially guess off*/
	p->_region = INITOFF;
	y.f1 = 0.;
	y.f0 = 0.;
	if (flags & 0020) {
	  untested();
	  y.f1 = OPT::gmin;
	}
	trace2("initoff", y.f0, y.f1);
      }else if (volts <= 0. /* &&  amps < 0.*/) {    	  /* reverse biased */
	p->_region = REVERSE;	    		  /* x = volts, f(x) = amps */
	if (flags & 0010) {
	  untested();
	  y.f1 = y.f0 = 0.;
	}else{
	  double expterm = p->_isat * exp(volts/vt);	
	  y.f0 = expterm - p->_isat;/* i = f(x) = _isat * (exp(volts/vt)-1) */
	  y.f1 = expterm / vt;	    /* f'(x) = (_isat/vt) * exp(volts/vt)   */
	}
	
	if (flags & 0002) {	// g = gmin, maintain actual current
	  y.f1 += OPT::gmin;	// 3 is a resistor, R=1/gmin
	  y.f0 += OPT::gmin * volts;
	}
	if (flags & 0004) {	// 5 is a resistor, R=vt/_isat
	  double x = p->_isat / vt;
	  y.f1 += x;
	  y.f0 += x * volts;
	}
	if (flags & 0001) {
	  //y.f0 = y.f1 * volts;	// a resistor, R=1/f1
	}
	
	trace2("reverse", y.f0, y.f1);
      }else if (volts >= 0.  &&  amps >= 0.) {		  /* forward biased */
				    /* x = amps, f(x) = volts */
	/* derivation: */	    /* if f(x) = log(u): f'(x)=(1/u)(du/dx) */
	/* poly1 r; */
	/* r.f0 = vt * log(amps/p->_isat +1.); */
	/* r.f1 = vt / (_isat + amps); */
	/* y.f1 = 1. / r.f1; */
	/* y.f0 = amps - r.f0*y.f1 + volts*y.f1; */
	
	p->_region = FORWARD;
	y.f1 = (p->_isat + amps) / vt;
	y.f0 = amps - log(amps/p->_isat +1.)*(p->_isat + amps) + volts*y.f1;
	trace2("forward", y.f0, y.f1);
      }else{			    /* non-converged, inconsistent	    */
	p->_region = UNKNOWN;	    /* volts and amps have different signs  */
	y.f1 = p->_isat/vt;	    /* guess that the voltage should be 0   */
	y.f0 = 0.;		    /* (it usually is very close)	    */
	if (flags & 0001) {	    /* use the correct value there	    */
	  y.f0 = volts * y.f1;
	}
	trace2("unknown", y.f0, y.f1);
      }
      y.f1 += c->gparallel_adjusted;
      y.f0 += c->gparallel_adjusted * volts;
      
      if (oldregion != p->_region  &&  OPT::dampstrategy & dsDEVLIMIT) {
	d->_sim->_fulldamp = true;
	error(bTRACE, p->long_label() + ":device limit damp\n");
      }
      if (flags & 0100) {		// twist g to guarantee g >= gmin
	if (y.f1 < OPT::gmin) {	// without changing i
	  y.f1 = OPT::gmin;
	  untested();
	}else{
	  untested();
	}
      }
      if (flags & 0200) {		// add a gmin in parallel
	y.f1 += OPT::gmin;
	y.f0 += OPT::gmin * volts;
	untested();
      }
      if (flags & 0400) {		// linearize .. shift I to pass thru 0
	untested();
	y.f0 = y.f1 * volts;
      }
    }
    trace3(d->long_label().c_str(), y.x, y.f0, y.f1);
    p->_gd = y.f1;
}
/*--------------------------------------------------------------------------*/
DEV_BUILT_IN_DIODE::DEV_BUILT_IN_DIODE()
  :BASE_SUBCKT(),
   // input parameters,
   // calculated parameters,
   _region(UNKNOWN),
   _gd(NA),
   _isat(NA),
   // netlist,
   _Cj(0),
   _Yj(0),
   _Rs(0)
{
  _n = _nodes;
  attach_common(&Default_BUILT_IN_DIODE);
  ++_count;
  // overrides
}
/*--------------------------------------------------------------------------*/
DEV_BUILT_IN_DIODE::DEV_BUILT_IN_DIODE(const DEV_BUILT_IN_DIODE& p)
  :BASE_SUBCKT(p),
   // input parameters,
   // calculated parameters,
   _region(p._region),
   _gd(p._gd),
   _isat(p._isat),
   // netlist,
   _Cj(0),
   _Yj(0),
   _Rs(0)
{
  _n = _nodes;
  for (uint_t ii = 0; ii < max_nodes() + int_nodes(); ++ii) {
    _n[ii] = p._n[ii];
  }
  ++_count;
  // overrides
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_DIODE::expand()
{
  BASE_SUBCKT::expand(); // calls common->expand, attached model
  assert(_n);
  assert(common());
  const COMMON_BUILT_IN_DIODE* c = static_cast<const COMMON_BUILT_IN_DIODE*>(common());
  assert(c);
  assert(c->model());
  const MODEL_BUILT_IN_DIODE* m = prechecked_cast<const MODEL_BUILT_IN_DIODE*>(c->model());
  assert(m);
  assert(c->sdp());
  const SDP_BUILT_IN_DIODE* s = prechecked_cast<const SDP_BUILT_IN_DIODE*>(c->sdp());
  assert(s); USE(s);
  if (!subckt()) {
    new_subckt();
  }else{
  }

  if (_sim->is_first_expand()) {
    precalc_first();
    precalc_last();
    // local nodes
    //assert(!(_n[n_ia].n_()));
    //BUG// this assert fails on a repeat elaboration after a change.
    //not sure of consequences when new_model_node called twice.
    if (!(_n[n_ia].n_())) {
      if (!OPT::rstray || c->rs_adjusted==0.) {
        _n[n_ia] = _n[n_a];
      }else{
        _n[n_ia].new_model_node("." + long_label() + ".ia", this);
      }
    }else{
      if (!OPT::rstray || c->rs_adjusted==0.) {
        assert(_n[n_ia] == _n[n_a]);
      }else{
        //_n[n_ia].new_model_node("ia." + long_label(), this);
      }
    }

    // clone subckt elements
    if (c->cj_adjusted == 0. && c->cjsw_adjusted == 0. && m->tt == 0.) {
      if (_Cj) {
        subckt()->erase(_Cj);
        _Cj = NULL;
      }else{
      }
    }else{
      if (!_Cj) {
        const CARD* p = device_dispatcher["capacitor"];
        assert(p);
        _Cj = dynamic_cast<COMPONENT*>(p->clone());
        assert(_Cj);
        subckt()->push_front(_Cj);
      }else{
      }
      {
        node_t nodes[] = {_n[n_ia], _n[n_c]};
        _Cj->set_parameters("Cj", this, &Eval_Cj, 0., 0, NULL, 2, nodes);
      }
    }
    {
      if (!_Yj) {
        const CARD* p = device_dispatcher["admittance"];
        assert(p);
        _Yj = dynamic_cast<COMPONENT*>(p->clone());
        assert(_Yj);
        subckt()->push_front(_Yj);
      }else{
      }
      {
        node_t nodes[] = {_n[n_ia], _n[n_c]};
        _Yj->set_parameters("Yj", this, &Eval_Yj, 0., 0, NULL, 2, nodes);
      }
    }
    if (!OPT::rstray || c->rs_adjusted==0.) {
      if (_Rs) {
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
      }else{
      }
      {
        node_t nodes[] = {_n[n_a], _n[n_ia]};
        _Rs->set_parameters("Rs", this, NULL, c->rs_adjusted, 0, NULL, 2, nodes);
      }
    }
  }else{
    //precalc();
  }
  //precalc();
  subckt()->expand();
  //subckt()->precalc();
  assert(!is_constant());
  trace0(("DEV_BUILT_IN_DIODE::expand " + long_label()).c_str() );
  if ( adp() == NULL ){
    assert(common());
    assert(c);
    //attach_adp( m->new_adp( (const COMPONENT*) this ) );
  }else{
    assert(false);
  }
  trace0(("DEV_BUILT_IN_DIODE::expand done " + long_label()).c_str() );
}
/*--------------------------------------------------------------------------*/
double DEV_BUILT_IN_DIODE::tr_probe_num(const std::string& x)const
{
  trace2("tr_probe_num" + short_label(), (intptr_t) common() %1000,
                                         (intptr_t)( common()->model()) %1000 );
  assert(_n);
  const COMMON_BUILT_IN_DIODE* c = prechecked_cast<const COMMON_BUILT_IN_DIODE*>(common());
  assert(c);
  const MODEL_BUILT_IN_DIODE* m = prechecked_cast<const MODEL_BUILT_IN_DIODE*>(c->model());
  assert(m); USE(m);
  const SDP_BUILT_IN_DIODE* s = prechecked_cast<const SDP_BUILT_IN_DIODE*>(c->sdp());
  assert(s); USE(s);
//  const ADP_BUILT_IN_DIODE* a = prechecked_cast<const ADP_BUILT_IN_DIODE*>(adp());

  if (Umatch(x, "v{d} ")) {
    return  _n[n_a].v0() - _n[n_c].v0();
  }else if (Umatch(x, "i{d} ")) {
    return  CARD::probe(_Yj,"I") + CARD::probe(_Cj,"I");
  }else if (Umatch(x, "vj ")) {
    return  _n[n_ia].v0() - _n[n_c].v0();
  }else if (Umatch(x, "vsr ")) {
    return  _n[n_a].v0() - _n[n_ia].v0();
  }else if (Umatch(x, "vrs ")) {
    return  _n[n_a].v0() - _n[n_ia].v0();
  }else if (Umatch(x, "ij ")) {
    return  CARD::probe(_Yj,"I");
  }else if (Umatch(x, "ic ")) {
    return  CARD::probe(_Cj,"I");
  }else if (Umatch(x, "capcur ")) {
    return  CARD::probe(_Cj,"I");
  }else if (Umatch(x, "p ")) {
    return  CARD::probe(_Yj,"P") + CARD::probe(_Cj,"P") + CARD::probe(_Rs,"P");
  }else if (Umatch(x, "pd ")) {
    return  CARD::probe(_Yj,"PD") + CARD::probe(_Cj,"PD") + CARD::probe(_Rs,"PD");
  }else if (Umatch(x, "ps ")) {
    return  CARD::probe(_Yj,"PS") + CARD::probe(_Cj,"PS") + CARD::probe(_Rs,"PS");
  }else if (Umatch(x, "pj ")) {
    return  CARD::probe(_Yj,"P");
  }else if (Umatch(x, "pc ")) {
    return  CARD::probe(_Cj,"P");
  }else if (Umatch(x, "c{apacitance} ")) {
    return  CARD::probe(_Cj,"Capacitance");
  }else if (Umatch(x, "cd ")) {
    return  CARD::probe(_Cj,"Capacitance");
  }else if (Umatch(x, "charge ")) {
    return  CARD::probe(_Cj,"Charge");
  }else if (Umatch(x, "r{eq} ")) {
    return  CARD::probe(_Yj,"R") + CARD::probe(_Rs,"R");
  }else if (Umatch(x, "g{eq} ")) {
    return  (( CARD::probe(_Yj,"R") + CARD::probe(_Rs,"R") ) != 0) ? (1./( CARD::probe(_Yj,"R") + CARD::probe(_Rs,"R") )) : CARD::probe(_Yj,"Y");
  }else if (Umatch(x, "gd ")) {
    return  CARD::probe(_Yj,"Y");
  }else if (Umatch(x, "y ")) {
    return  ( CARD::probe(_Rs,"R") != 0. && ( CARD::probe(_Yj,"Y") + CARD::probe(_Cj,"Y") ) != 0) ? 1./((1./( CARD::probe(_Yj,"Y") + CARD::probe(_Cj,"Y") )) + CARD::probe(_Rs,"R") ) : CARD::probe(_Yj,"Y") + CARD::probe(_Cj,"Y");
  }else if (Umatch(x, "z ")) {
    return  port_impedance( _n[n_a] , _n[n_c] ,_sim->_lu,mfactor()*tr_probe_num("Y"));
  }else if (Umatch(x, "zraw ")) {
    return  port_impedance( _n[n_a] , _n[n_c] , _sim->_lu, 0.);
  }else if (Umatch(x, "region ")) {
    return  static_cast<double>(_region);
  }else if (Umatch(x, "_region ")) {
    return _region;
  }else if (Umatch(x, "_gd ")) {
    return _gd;
  }else if (Umatch(x, "_isat ")) {
    return _isat;
  }else {
    return BASE_SUBCKT::tr_probe_num(x);
  }
}
/*--------------------------------------------------------------------------*/
double DEV_BUILT_IN_DIODE::tt_probe_num(const std::string& x)const
{
  assert(_n);
  const COMMON_BUILT_IN_DIODE* c = prechecked_cast<const COMMON_BUILT_IN_DIODE*>(common());
  assert(c);
  const MODEL_BUILT_IN_DIODE* m = prechecked_cast<const MODEL_BUILT_IN_DIODE*>(c->model());
  assert(m); USE(m);
  const SDP_BUILT_IN_DIODE* s = prechecked_cast<const SDP_BUILT_IN_DIODE*>(c->sdp());
  assert(s); USE(s);
  // const ADP_BUILT_IN_DIODE* a = prechecked_cast<const ADP_BUILT_IN_DIODE*>(adp());
  //  if(!a)untested0("no a");

  {
    return BASE_SUBCKT::tt_probe_num(x);
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// cc_direct



void ADP_BUILT_IN_DIODE::init(const COMPONENT* cc)
{
  trace0(( "ADP_BUILT_IN_DIODE::init " + cc->short_label() ).c_str() );
  assert(cc); USE(cc);
}

/*--------------------------------------------------------------------------*/
ADP_CARD* MODEL_BUILT_IN_DIODE::new_adp( COMPONENT* c)const
{
  assert(c);
  return MODEL_CARD::new_adp(c);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
}
