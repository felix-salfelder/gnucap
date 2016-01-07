/* $Id: d_diode.h,v 1.7 2010-07-16 14:49:39 felix Exp $ -*- C++ -*-
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
#ifndef D_DIODE_H_INCLUDED
#define D_DIODE_H_INCLUDED

  enum region_t {INITOFF=-2, REVERSE=-1, UNKNOWN=0, FORWARD=1};
  //enum polarity_t {pP = -1, dunno=0, pN = 1};
#include "e_adp.h"

#include "u_sdp.h"
#include "e_node.h"
#include "e_subckt.h"
#include "e_model.h"
#include <map>
/*--------------------------------------------------------------------------*/
class SDP_BUILT_IN_DIODE
  :public SDP_CARD{
public:
  explicit SDP_BUILT_IN_DIODE(const COMMON_COMPONENT* c) : SDP_CARD(c) {init(c);}
  void init(const COMMON_COMPONENT*);
public:
};
/*--------------------------------------------------------------------------*/
class DEV_BUILT_IN_DIODE;
class TDP_BUILT_IN_DIODE{
public:
  explicit TDP_BUILT_IN_DIODE(const DEV_BUILT_IN_DIODE*);
public:
};
/*--------------------------------------------------------------------------*/
class MODEL_BUILT_IN_DIODE
  :public MODEL_CARD{
protected:
  explicit MODEL_BUILT_IN_DIODE(const MODEL_BUILT_IN_DIODE& p);
public:
  explicit MODEL_BUILT_IN_DIODE(const BASE_SUBCKT*);
  ~MODEL_BUILT_IN_DIODE() {--_count;}
public: // override virtual
  std::string dev_type()const;
  void      set_dev_type(const std::string& nt);
  CARD*     clone()const {return new MODEL_BUILT_IN_DIODE(*this);}
  void      precalc_first();
  void      precalc_last();
  SDP_CARD* new_sdp(COMMON_COMPONENT* c)const;
  ADP_CARD* new_adp(COMPONENT* c)const;
  void      set_param_by_index(int, std::string&, int);
  bool      param_is_printable(int)const;
  std::string param_name(int)const;
  std::string param_name(int,int)const;
  std::string param_value(int)const;
  int param_count()const {return (22);}
  bool      is_valid(const COMPONENT*)const;
  void      tr_eval(COMPONENT*)const;
  virtual void      do_stress_apply(COMPONENT*)const{ std::cerr<<"virtual stress apply(C)\n" ;}
public: // not virtual
  static int count() {return _count;}
private: // strictly internal
  static int _count;
public: // input parameters
  PARAMETER<double> js;	// = is, saturation current (per area)
  PARAMETER<double> rs;	// ohmic resistance (per area)
  PARAMETER<double> n_factor;	// emission coefficient
  PARAMETER<double> tt;	// transit time
  PARAMETER<double> cjo;	// cj, zero-bias jct capacitance (per area)
  PARAMETER<double> pb;	// vj, junction potential
  PARAMETER<double> mj;	// m, grading coefficient
  PARAMETER<double> eg;	// activation energy
  PARAMETER<double> xti;	// saturation-current temp. exp.
  PARAMETER<double> kf;	// flicker noise coefficient
  PARAMETER<double> af;	// flicker noise exponent
  PARAMETER<double> fc;	// coef for fwd bias depl cap formula
  PARAMETER<double> bv;	// reverse breakdown voltage
  PARAMETER<double> ibv;	// current at reverse breakdown
  PARAMETER<double> cjsw;	// zero bias sidewall cap (per perim.)
  PARAMETER<double> pbsw;	// sidewall junction potential
  PARAMETER<double> mjsw;	// sidewall grading coefficient
  PARAMETER<double> gparallel;	// parallel conductance
  PARAMETER<int> flags;	// 
  PARAMETER<int> mos_level;	// 
public: // calculated parameters
};
/*--------------------------------------------------------------------------*/
class COMMON_BUILT_IN_DIODE
  :public COMMON_COMPONENT{
public:
  explicit COMMON_BUILT_IN_DIODE(const COMMON_BUILT_IN_DIODE& p);
  explicit COMMON_BUILT_IN_DIODE(int c=0);
           ~COMMON_BUILT_IN_DIODE();
  bool     operator==(const COMMON_COMPONENT&)const;
  COMMON_COMPONENT* clone()const {return new COMMON_BUILT_IN_DIODE(*this);}
  void     set_param_by_index(int, std::string&, int);
  void     set_param_by_name(std::string, std::string);
  bool     param_is_printable(int)const;
  std::string param_name(int)const;
  std::string param_name(int,int)const;
  std::string param_value(int)const;
  int param_count()const {return (9 + COMMON_COMPONENT::param_count());}
  void     precalc_first(const CARD_LIST*);
  void     expand(const COMPONENT*);
  void     precalc_last(const CARD_LIST*);
  std::string name()const {return "diode";}
  const SDP_CARD* sdp()const {return _sdp;}
  bool     has_sdp()const {untested();return _sdp;}
  static int  count() {return _count;}
private: // strictly internal
  static std::map<std::string, PARA_BASE COMMON_BUILT_IN_DIODE::*> param_dict;
  static std::map<std::string, PARA_BASE COMMON_BUILT_IN_DIODE::*> param_dict_low;
  static int _count;
public: // input parameters
  PARAMETER<double> area;	// area factor
  PARAMETER<double> perim;	// perimeter factor
  PARAMETER<bool> off;	// flag: assume reverse biased
  PARAMETER<double> ic;	// initial voltage
  PARAMETER<double> is_raw;	// saturation current
  PARAMETER<double> rs_raw;	// series resistance
  PARAMETER<double> cj_raw;	// zero bias jct capacitance
  PARAMETER<double> cjsw_raw;	// zero bias sidewall capacitance
  PARAMETER<double> gparallel_raw;	// parallel conductance
public: // calculated parameters
  SDP_CARD* _sdp;
  double is_adjusted;	// 
  double rs_adjusted;	// 
  double cj_adjusted;	// 
  double cjsw_adjusted;	// 
  double gparallel_adjusted;	// 
public: // attached commons
};
/*--------------------------------------------------------------------------*/
class EVAL_BUILT_IN_DIODE_Cj : public COMMON_COMPONENT {
private:
  explicit EVAL_BUILT_IN_DIODE_Cj(const EVAL_BUILT_IN_DIODE_Cj& p)
    :COMMON_COMPONENT(p) {}
public:
  explicit EVAL_BUILT_IN_DIODE_Cj(int c=0) :COMMON_COMPONENT(c) {}
  bool operator==(const COMMON_COMPONENT& x)const {return COMMON_COMPONENT::operator==(x);}
  COMMON_COMPONENT* clone()const {return new EVAL_BUILT_IN_DIODE_Cj(*this);}
  std::string name()const {untested(); return "EVAL_BUILT_IN_DIODE_Cj";}
  void tr_eval(ELEMENT*d)const;
  bool has_tr_eval()const {return true;}
  bool has_ac_eval()const {return false;}
};
/*--------------------------------------------------------------------------*/
class EVAL_BUILT_IN_DIODE_Yj : public COMMON_COMPONENT {
private:
  explicit EVAL_BUILT_IN_DIODE_Yj(const EVAL_BUILT_IN_DIODE_Yj& p)
    :COMMON_COMPONENT(p) {}
public:
  explicit EVAL_BUILT_IN_DIODE_Yj(int c=0) :COMMON_COMPONENT(c) {}
  bool operator==(const COMMON_COMPONENT& x)const {return COMMON_COMPONENT::operator==(x);}
  COMMON_COMPONENT* clone()const {return new EVAL_BUILT_IN_DIODE_Yj(*this);}
  std::string name()const {untested(); return "EVAL_BUILT_IN_DIODE_Yj";}
  void tr_eval(ELEMENT*d)const;
  bool has_tr_eval()const {return true;}
  bool has_ac_eval()const {return false;}
};
/*--------------------------------------------------------------------------*/
class DEV_BUILT_IN_DIODE : public BASE_SUBCKT {
private:
  explicit DEV_BUILT_IN_DIODE(const DEV_BUILT_IN_DIODE& p);
public:
  explicit DEV_BUILT_IN_DIODE();
           ~DEV_BUILT_IN_DIODE() {--_count;}
private: // override virtual
  char      id_letter()const     {untested();return 'D';}
  bool      print_type_in_spice()const {return true;}
  std::string value_name()const  {return "area";}
  //std::string dev_type()const;   //BASE_SUBCKT
  uint_t       max_nodes()const     {return 2;}
  uint_t       min_nodes()const     {return 2;}
  //int     matrix_nodes()const; //BASE_SUBCKT
  uint_t       net_nodes()const     {return 2;}
  uint_t       int_nodes()const     {return 1;}
  CARD*     clone()const         {return new DEV_BUILT_IN_DIODE(*this);}
  void      precalc_first() {COMPONENT::precalc_first(); if(subckt()) subckt()->precalc_first();}
  void      expand();
  void      precalc_last()  {COMPONENT::precalc_last(); assert(subckt()); subckt()->precalc_last();}
  //void    map_nodes();         //BASE_SUBCKT
  //void    tr_begin();          //BASE_SUBCKT
  //void    tr_restore();        //BASE_SUBCKT
  //void    tt_commit();         //BASE_SUBCKT
  //void    tt_prepare();         //BASE_SUBCKT
  //void    dc_advance();        //BASE_SUBCKT
  //void    tr_advance();        //BASE_SUBCKT
  //void    tr_regress();        //BASE_SUBCKT
  //bool    tr_needs_eval()const;//BASE_SUBCKT
  //void    tr_queue_eval();     //BASE_SUBCKT
  //bool    do_tr();             //BASE_SUBCKT
  //void    tr_load();           //BASE_SUBCKT
  //double  tr_review();         //BASE_SUBCKT
  //void    tr_accept();         //BASE_SUBCKT
  //void    tr_unload();         //BASE_SUBCKT
  double    tr_probe_num(const std::string&)const;
  double    tt_probe_num(const std::string&)const;
  //void    ac_begin();          //BASE_SUBCKT
  //void    do_ac();             //BASE_SUBCKT
  //void    ac_load();           //BASE_SUBCKT
  //XPROBE  ac_probe_ext(CS&)const;//CKT_BASE/nothing
public:
  static int  count() {return _count;}
public: // may be used by models
private: // not available even to models
  static int _count;
public: // input parameters
public: // calculated parameters
  region_t _region;	// fwd, reverse, unknown
  double _gd;	// conductance to pass to capacitor
  double _isat;	// is adjusted for temp, etc.
public: // netlist
  COMPONENT* _Cj;
  COMPONENT* _Yj;
  COMPONENT* _Rs;
private: // node list
  enum {n_a, n_c, n_ia};
  node_t _nodes[3];
  std::string port_name(uint_t i)const {
    assert(i < 2);
    static std::string names[] = {"a", "c", ""};
    return names[i];
  }
};
/*--------------------------------------------------------------------------*/
// h_direct


class ADP_BUILT_IN_DIODE :public ADP_CARD{
public:
  explicit ADP_BUILT_IN_DIODE( COMPONENT* c, const std::string n) :
    ADP_CARD(c,n) 
    {init(c);}
protected:
  void init(const COMPONENT*);
};

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
#endif
