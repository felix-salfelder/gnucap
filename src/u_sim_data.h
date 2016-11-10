/*                             -*- C++ -*-
 * Copyright (C) 2009 Albert Davis
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
 * real base for anything to do with a circuit
 */
//testing=obsolete
#ifndef U_SIM_DATA_H
#define U_SIM_DATA_H
#include "constant.h"
#include "l_compar.h"
#include "u_opt.h"
#include "m_matrix.h"
// #include "e_cardlist.h"
#include "e_base.h"
#include <stack>
/*--------------------------------------------------------------------------*/
#include "io_error.h"
/*--------------------------------------------------------------------------*/
class CARD_LIST;
/*--------------------------------------------------------------------------*/
template<class T>
class MYSTACK : public std::stack<T>
{
public:
    const T& second() const throw(Exception_Cant_Find);
    // T& second();
};
/*--------------------------------------------------------------------------*/
template<class T>
const T& MYSTACK<T>::second() const throw(Exception_Cant_Find)
    {
      if(std::stack<T>::size()<2){ untested();
        throw Exception_Cant_Find("","");
      }else{
        typename std::stack<T>::container_type::const_reverse_iterator a = std::stack<T>::c.rbegin();
        return *(++a);
      }
    }
/*--------------------------------------------------------------------------*/
typedef MYSTACK<double*> VOLTAGE_STACK;
/*--------------------------------------------------------------------------*/
// external
class WAVE;
class CARD;
class LOGIC_NODE;
/*--------------------------------------------------------------------------*/
enum TRI_STATE {tsNO=0, tsYES=1, tsBAD=-1};
/*--------------------------------------------------------------------------*/
#define HAVE_SIM_DT0
struct INTERFACE SIM_DATA {
  double _dt0; // hack. remove (sooner or later..)
  double _time0;	/* time now */
  double _Time0;	/* Time now */
  double _dT0;	        /* last step length (including transient time) */
  double _dT1;	
  double _dT2;	/* HACK? */
  double _dT3;	/* HACK? */
  //union{
    double _freq;		/* AC frequency to analyze at (Hertz) */
    double _dxm;
  //};
  double _temp_c;	/* ambient temperature, actual */
  double _damp;		/* Newton-Raphson damping coefficient actual */
  double _dtmin;	/* min internal step size */
  double _dTmin;	/* min internal step size */
  double _genout;	/* tr dc input to circuit (generator) */
  bool   _bypass_ok;	/* flag: ok to bypass model evaluation */
  bool   _fulldamp; 	/* flag: big iter. jump. use full (min) damp */
  double _last_time;	/* time at which "volts" is valid */
public:
  double last_time()const {return _last_time;}
  double _last_Time;	/* end of last accepted timeframe */
  double last_Time()const {return _last_Time;} // recheck: reuse _last_time?
                                               // need some keep_voltages kludge
  void tr_reset(){_last_time=0;}
public:
  bool _freezetime;	/* flag: don't advance stored time */
  int _iter[iCOUNT];
  uint_t _user_nodes;
  uint_t _subckt_nodes;
  uint_t _model_nodes;
  uint_t _total_nodes;
  uint_t _adp_nodes;
  COMPLEX _jomega;	/* AC frequency to analyze at (radians) */
  bool _limiting;	/* flag: node limiting */
  double _vmax;
  double _vmin;
  bool _age;		/* flag: enable aging capabilities */
  bool _uic;		/* flag: use initial conditions (spice-like) */
  bool _tt_uic;         /* flag: use saved state in tt */
  bool _cont;           /* continue from dc point */
  TRI_STATE _inc_mode;	/* flag: make incremental changes (3 state) */
  SIM_MODE _mode;	/* simulation type (AC, DC, ...) */
  SIM_PHASE _phase;	/* phase of simulation (iter, init-dc,) */
  unsigned *_nm;	/* node map (external to internal)	*/
  double *_i;		/* dc-tran current (i) vector		*/
  double *_v0;		/* dc-tran voltage, new			*/
  double *_vt1;		/* dc-tran voltage, 1 time ago		*/
			/*  used to restore after rejected step	*/
  double *_vt2;		/* dc-tran voltage, 1 time step after dc    necessary??*/ 
  COMPLEX *_ac;		/* ac right side			*/
  COMPLEX *_sens;       /* sensitivity (noise etc.) */
  double *_tr;
  double *_tr1;
  double *_tr2;
  double *_tr3;

  double *_tt1;
  //double (* (_tr))[3];
  //double (* (_tt))[3];
  LOGIC_NODE* _nstat;	/* digital data				*/
  double *_tt;          /* aging node data */
  VOLTAGE_STACK _vdcstack;
  double* vdc() const{ assert(_vdcstack.size()); return _vdcstack.top(); }
  double v0dist() const;
  double v0dist_min(double* what=NULL) const;
  BSMATRIX<double> _aa;	/* raw matrix for DC & tran */
  BSMATRIX<double> _lu;	/* decomposed matrix for DC & tran */
  BSMATRIX<COMPLEX> _acx;/* raw & decomposed matrix for AC */
  std::priority_queue<double, std::vector<double> > _eq; /*event queue*/
  std::vector<CARD*> _loadq;
  std::vector<CARD*> _acceptq;
  std::vector<CARD*> _tt_acceptq;
  std::deque<CARD*>  _pre_evalq; /* candidates for evaluation */
  std::deque<CARD*>  _evalq1; /* evaluate queues -- alternate between */
  std::deque<CARD*>  _evalq2; /* build one while other is processed */
  std::deque<CARD*>  _late_evalq; /* eval after everything else */
  std::deque<CARD*>* _evalq;   /* pointer to evalq to process */
  std::deque<CARD*>* _evalq_uc;/* pointer to evalq under construction */
  std::map<std::string,std::map<std::string,WAVE> > _waves; /* wave storage, "store" command */
  std::string _label;
  SIM_MODE _has_op;
  SIM_DATA();
  ~SIM_DATA();
  bool is_first_expand() {return !_nstat;}
  void alloc_hold_vectors(); /* s__init.cc */
  void alloc_vectors();
  void unalloc_vectors();
  void init(bool need_precalc=1);
  void uninit();
  void set_limit();  /* s__aux.cc */
  void set_limit(double v);
  void clear_limit();
  void put_v1_to_v0(); /// temporary hack to buffer voltages
  void keep_voltages(bool push=false);
  void restore_voltages(bool pop=false);
  void pop_voltages();
  void zero_voltages();
  void zero_dc_voltages();
  void zero_some_voltages();
  void zero_currents();
  void map__nodes();		/* s__map.cc */
  void order_reverse();
  void order_sink_reverse();
  void order_forward();
  void order_sink_forward();
  void order_auto();
  void order_tree_bf(const CARD_LIST* scope=NULL);
  void order_tree_df(const CARD_LIST* scope=NULL);
  void order_comp(const CARD_LIST* scope=NULL);

  //int init_node_count( const CARD_LIST* l=&CARD_LIST::card_list); 
  unsigned init_node_count(unsigned user, unsigned sub, unsigned mod, unsigned adp) {
    _user_nodes=user; _subckt_nodes=sub; _model_nodes=mod; 
    _adp_nodes=adp;
    return (_total_nodes=user+sub+mod);
  }
  unsigned newnode_subckt() {++_subckt_nodes; return ++_total_nodes;}
  unsigned newnode_model()  {++_model_nodes;  return ++_total_nodes;}
  unsigned newnode_adp()  {return _adp_nodes++;}
  bool is_inc_mode()	 {return _inc_mode;}
  bool inc_mode_is_no()	 {return _inc_mode == tsNO;}
  bool inc_mode_is_bad() {return _inc_mode == tsBAD;}
  void set_inc_mode_bad() {_inc_mode = tsBAD;}
  void set_inc_mode_yes() {_inc_mode = tsYES;}
  void set_inc_mode_no()  {_inc_mode = tsNO;}
  void mark_inc_mode_bad() {
    switch (_inc_mode) {
    case tsYES: _inc_mode = tsBAD; break;
    case tsBAD: break;
    case tsNO:  break;
    }
  }
  void new_event(double etime) {
    if (etime <= BIGBIG) {
      trace2("SIM_DATA queueing event", etime, _time0);
      _eq.push(etime);
    }else{
    }
  }
  void set_command_none() {_mode = s_NONE; _label="";}
  void set_command_ac()	  {_mode = s_AC; _label="ac";}
  void set_command_dc()	  {_mode = s_DC; _label="dc";}
  void set_command_ddc()  {_mode = s_DDC; _label="ddc";}
  void set_command_op()	  {_mode = s_OP; _label="op";}
  void set_command_tran() {_mode = s_TRAN;}
  void set_command_fourier() {_mode = s_FOURIER; _label="fourier";}
  void set_command_sens() {_mode = s_SENS; _label="sens";}
  void set_command_tt() {_mode = s_TTT;}
  void set_label(std::string x) { _label=x; }
  SIM_MODE sim_mode()	   {return _mode;}
  SIM_PHASE phase()	   {return _phase;}
  bool command_is_ac()	   {return _mode == s_AC;}
  bool command_is_dc()	   {return _mode == s_DC;}
  bool command_is_op()	   {return _mode == s_OP;}
  //bool command_is_tran()    {return _mode == s_TRAN;}
  //bool command_is_fourier() {return _mode == s_FOURIER;}
  bool analysis_is_ac()      {return _mode == s_AC || 
                                    (_phase == p_AC && _mode == s_SOCK);}
  bool analysis_is_dcop()    {return _mode == s_DC || _mode == s_OP || _mode==s_DDC ;}
  bool analysis_is_static()const  {return _phase == p_INIT_DC || _phase == p_DC_SWEEP;}
  bool analysis_is_restore() {return _phase == p_RESTORE;}
  bool analysis_is_tran()    {return _mode == s_TRAN || _mode == s_FOURIER || _mode == s_TTT;}
  bool analysis_is_tt()      {return _mode == s_TTT;}
  bool analysis_is_tran_static()  {return analysis_is_tran() && _phase == p_INIT_DC;}
  bool analysis_is_tran_restore() {return analysis_is_tran() && _phase == p_RESTORE;}
  bool analysis_is_tran_dynamic() {return analysis_is_tran() && _phase == p_TRAN;}
  bool analysis_is_ddc() {return _mode == s_DDC; }
  bool analysis_is_ddc_op() {return analysis_is_ddc() && _phase == p_INIT_DC;}
  bool analysis_is_sens() {return _mode == s_SENS; }

  void reset_iteration_counter(int i) {assert(up_order(0,i,iCOUNT-1)); _iter[i] = 0;}
  void count_iterations(int i)	{assert(up_order(0,i,iCOUNT-1)); ++_iter[i];}
  int iteration_tag()const      {return _iter[iTOTAL];}
  int iteration_number()const   {return _iter[iSTEP];}
  bool is_initial_step()	{return (_iter[_mode] <= 1  && analysis_is_static());}
  bool is_advance_iteration()const   {return (_iter[iSTEP] == 0);}
  bool is_advance_or_first_iteration()const {assert(_iter[iSTEP]>=0); return (_iter[iSTEP]<=1);}
  bool is_first_iteration()const     {assert(_iter[iSTEP] > 0); return (_iter[iSTEP] == 1);}
  bool is_second_iteration()const    {assert(_iter[iSTEP] > 0); return (_iter[iSTEP] == 2);}
  bool is_iteration_number(int n)const    {return (_iter[iSTEP] == n);}
  bool exceeds_iteration_limit(OPT::ITL itlnum)const {return(_iter[iSTEP] > OPT::itl[itlnum]);}
  bool uic_now() const {return _uic && analysis_is_static() && _time0==0.;}
  SIM_MODE has_op()const {return _has_op;}

private: //  TT stuff and hacks.
  uint_t _tt_order;
  uint_t last_order_tt;
  CS* _expect_file;

  public: // tt hacks/extensions partly lacking clean implementation
  unsigned _tt_iter;
  unsigned _tt_done;
  unsigned _tt_accepted;
  unsigned _tt_rejects;
  unsigned _tt_rejects_total;

  unsigned tt_iteration_number()const {return _tt_iter;} // accepted steps
  unsigned tt_iteration_tries()const {return _tt_done;} // accepted steps
  bool is_initial_step_tt()const {return (_tt_iter < 1 );}

  public:
  int _stepno; // number of transient steps accepted.
  void update_tt_order();
  uint_t get_tt_order() const;
  void invalidate_tt();
  void force_tt_order(uint_t i){_tt_order = i;}
  unsigned total_outsteps()const;
  std::vector<double>_expect_raw;

public: // expect (dead code...)
  void expect(CS* x){ assert (!_expect_file); _expect_file = x; }
  CS* expect_file()const { return _expect_file; }
};
/*--------------------------------------------------------------------------*/
inline int CKT_BASE::iteration_number()const
{ untested();
  return _sim->_iter[iSTEP];
}
/*--------------------------------------------------------------------------*/
inline unsigned CKT_BASE::tt_iteration_number()const
{
  // accepted steps...
  return _sim->_tt_iter;
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
#endif
// vim:ts=8:sw=2:noet:
