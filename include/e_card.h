/*                           -*- C++ -*-
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
 * base class for anything in a netlist or circuit file
 */
#ifndef E_CARD_H
#define E_CARD_H
#define HAVE_TT 2
#include "e_base.h"
#include "u_time_pair.h"
/*--------------------------------------------------------------------------*/
// this file
class CARD;
#define PARAM_LIST PARAM_LIST_MAP
/*--------------------------------------------------------------------------*/
// external
class node_t;
class CARD_LIST;
class PARAM_LIST_BASE;
class PARAM_LIST;
class LANGUAGE;
class COMPONENT;
/*--------------------------------------------------------------------------*/
class INTERFACE CARD : public CKT_BASE {
private:
  mutable int	_evaliter;	// model eval iteration number
  CARD_LIST*	_subckt;
  CARD* 	_owner;
  bool		_constant;	// eval stays the same every iteration
  std::string	_comment;
protected:
  node_t*	_n;
public:
  uint_t 	_net_nodes;	// actual number of "nodes" in the netlist
  //--------------------------------------------------------------------
public:   				// traversal functions
  CARD* find_in_my_scope(const IString& name);
  const CARD* find_in_my_scope(const IString& name)const;
  const CARD* find_in_parent_scope(const IString& name)const;
  const CARD* find_looking_out(const IString& name)const;
  //--------------------------------------------------------------------
protected: // create and destroy.
  explicit CARD();
  explicit CARD(const CARD&);
public:
  virtual  ~CARD();
  virtual CARD*	 clone()const = 0;
  virtual CARD*	 clone_instance()const  {return clone();}
  //--------------------------------------------------------------------
public:	// "elaborate"
  virtual void	 precalc_first()	{}
  virtual void	 expand_first()		{}
  virtual void	 expand()		{}
  virtual void	 expand_last()		{}
  virtual void	 precalc_last()		{}
  virtual void	 map_nodes()		{}
  //--------------------------------------------------------------------
public:	// dc-tran
  virtual void	 tr_iwant_matrix()	{}
  virtual void	 tr_begin()		{}
  virtual void	 tr_restore()		{}
  virtual void	 dc_advance()		{}
  virtual void	 tr_advance()		{}
  virtual void	 tr_regress()		{}
  virtual void	 keep_ic()		{}

  virtual bool	 tr_needs_eval()const	{return false;}
  virtual void	 tr_queue_eval()	{} // not const, would need mutable iteration_tag
  virtual bool	 do_tr()		{return true;}
  virtual bool	 do_tr_last()		{return true;}
  virtual void	 tr_load()		{}
  virtual TIME_PAIR tr_review();	//{return TIME_PAIR(NEVER,NEVER);}
  virtual void	 tr_accept()		{}
  virtual void	 tr_unload()		{untested(); assert(false);}
  //--------------------------------------------------------------------
public:	// ac
  virtual void	 ac_iwant_matrix()	{}
  virtual void	 ac_begin()		{}
  virtual void	 do_ac()		{}
  virtual void	 ac_load()		{}
  //--------------------------------------------------------------------
public:	// noise
  virtual double do_noise() const {return 0; incomplete(); trace1("not implemented", typeid(*this).name()); }
  //--------------------------------------------------------------------
public:	// sens
  //  virtual void sens_load() {} // ckt_base
  virtual void do_sens() {}
  //--------------------------------------------------------------------
public:	// state, aux data
  // not unreachable. some devices are simply not spice.
  virtual char id_letter()const	{return '\0';}
  virtual uint_t  net_nodes()const {return 0;}
  virtual bool is_device()const	{return false;}
  virtual void set_slave()	{untested(); assert(!subckt());}
	  bool evaluated()const;

  void	set_constant(bool c)	{_constant = c;}
  bool	is_constant()const	{return _constant;}
  //--------------------------------------------------------------------
public: // owner, scope
  virtual CARD_LIST*	   scope();
  virtual const CARD_LIST* scope()const;
  virtual bool		   makes_own_scope()const  {return false;}
  CARD*		owner()		   {return _owner;}
  const CARD*	owner()const	   {return _owner;}
  void		set_owner(CARD* o) {assert(!_owner||_owner==o); _owner=o;}
  //--------------------------------------------------------------------
public: // subckt
  CARD_LIST*	     subckt()		{return _subckt;}
  const CARD_LIST*   subckt()const	{return _subckt;}
  void	  new_subckt(PARAM_LIST* p=NULL);
  void	  new_subckt(const CARD* model, PARAM_LIST* p);
  void	  renew_subckt(const CARD* model, PARAM_LIST* p);
  void    new_subckt(const CARD* model, CARD* Owner, const CARD_LIST* Scope, PARAM_LIST* p){
    USE(Scope); assert(Scope==scope());
    USE(Owner); assert(Owner==this);
    new_subckt(model, p);
  }
  void    renew_subckt(const CARD* model, CARD* Owner, const CARD_LIST* Scope, PARAM_LIST* p){
    USE(Scope); assert(Scope==scope());
    USE(Owner); assert(Owner==this);
    renew_subckt(model, p);
  }

  //--------------------------------------------------------------------
public:	// type
  virtual std::string dev_type()const	{unreachable(); return "";}
  virtual void set_dev_type(const std::string&);
  //--------------------------------------------------------------------
public:	// label -- in CKT_BASE
  // non-virtual void set_label(const std::string& s) //BASE
  // non-virtual const std::string& short_label()const //BASE
  /*virtual*/ const std::string long_label()const; // no further override
  //--------------------------------------------------------------------
public:	// ports -- mostly defer to COMPONENT
  node_t& n_(unsigned i)const;
  int     connects_to(const node_t& node)const;
  //--------------------------------------------------------------------
public: // parameters
  virtual void set_param_by_name(std::string, std::string);
  virtual void set_param_by_index(int i, std::string&, int offset)
				{untested(); throw Exception_Too_Many(unsigned(i), 0u, offset);}
  virtual int  param_count_dont_print()const	   {return 0;}
  virtual int  param_count()const		   {return 0;}
  virtual bool param_is_printable(int)const	   {untested(); return false;}
  virtual std::string param_name(int)const	   {return "";}
  virtual std::string param_name(int i,int j)const {return (j==0) ? param_name(i) : "";}
  virtual std::string param_value(int)const	   {untested(); return "";}
  virtual std::string value_name()const = 0;
  //--------------------------------------------------------------------
public:	// obsolete -- do not use in new code
  virtual bool use_obsolete_callback_parse()const {return false;}
  virtual bool use_obsolete_callback_print()const {return false;}
  virtual void print_args_obsolete_callback(OMSTREAM&,LANGUAGE*)const {unreachable();}
  //--------------------------------------------------------------------
public:	// tt
  virtual void	 tt_begin(){}
  virtual void	 tt_restore(){}
  virtual void   tt_advance(){}  // prepare for next sweep
  virtual void	 tt_regress(){}  // same but backwards in time
  virtual void   tt_accept(){}
  virtual void  tt_init_i(){}       // save unstressed parameters
//  virtual void	 tr_stress(){}      // calculate stress during tr. tr_accept?
  virtual void	 do_tt(){};         // called before tr_begin.
  virtual void	 tr_stress_last(){} // calculate stress during tr. tt_review? do_tt_last?
  virtual TIME_PAIR tt_review()		{return TIME_PAIR(NEVER,NEVER);}
public: /// experimental & cruft
  std::string comment() const {return _comment;}
  void set_comment(std::string s) {_comment = s;}
  virtual void	 tr_save_amps( int ){ } // behaviouir??
  hp_float_t tr_behaviour_del; // behaviour now.
  hp_float_t tt_behaviour_del;
  hp_float_t tr_behaviour_rel; 
  hp_float_t tt_behaviour_rel;
  void tt_behaviour_reset() { tt_behaviour_del=0; tt_behaviour_rel=0; }
  void tt_behaviour_commit(){ tt_behaviour_reset(); }
};
/*--------------------------------------------------------------------------*/

template <class S>
inline S& operator<<( S& o, const  std::deque<CARD*> &d){
  for(std::deque<CARD*>::const_iterator i=d.begin(); i!=d.end(); ++i){
       o << "\n" << (*i)->long_label() << " " << (*i)->short_label() ; 
  }
  return o<<"\n";

}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
#endif
// vim:ts=8:sw=2:noet:
