/*                              -*- C++ -*-
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
 * current controlled source base
 */
#ifndef E_CCSRC_H
#define E_CCSRC_H
#include "e_elemnt.h"
#include <string>
/*--------------------------------------------------------------------------*/
class INTERFACE CCSRC_BASE : public ELEMENT {
private:
protected:
  explicit	CCSRC_BASE()
    :ELEMENT(), _input_label(), _input(0) {}
  explicit	CCSRC_BASE(const CCSRC_BASE& p)
    :ELEMENT(p), _input_label(p._input_label), _input(p._input) {}
  ~CCSRC_BASE() {}
protected: // override virtual
  virtual uint_t	   max_nodes()const	{assert(input_order()); return 2 + 1*input_order();}
  virtual uint_t	   net_nodes()const	{return 2;} // before:2
  virtual uint_t	   ext_nodes()const     {return 2*(1+input_order());} // for matrix allocation
  virtual uint_t	   num_current_ports()const {return input_order();}
  virtual uint_t	   ports_per_input()const {return 0;}
  const std::string current_port_value(uint_t i)const {
    if(i<_input_label.size())
      return _input_label[i];
    return "";
  }
  // bool     do_tr()		{_sim->_late_evalq.push_back(this); return true;}
  bool do_tr_hack() {return do_tr_last();}
  //void   precalc_first();	//ELEMENT
  void	   expand_last();
  //void   precalc_last();	//ELEMENT
  bool	   tr_needs_eval()const	{assert(!is_q_for_eval()); return true;}
  //void   tr_queue_eval()	//ELEMENT
  void	   tr_unload()		{untested(); tr_unload_active();}
  hp_float_t   tr_involts()const	{untested();return dn_diff(_n[IN1].v0(), _n[IN2].v0());}
  hp_float_t   tr_input()const	{untested(); return _input->tr_amps();}
  hp_float_t   tr_involts_limited()const {return volts_limited(_n[IN1],_n[IN2]);}
  hp_float_t   tr_input_limited()const {return _input->tr_amps();}
  COMPLEX  ac_involts()const	{untested();return _n[IN1].vac()-_n[IN2].vac();}
  void	   set_port_by_index(uint_t index, std::string& value);
  bool	   node_is_connected(uint_t i)const;
  bool _master;
public:
  void set_inputs(unsigned node_count, node_t* Nodes);
  void	set_parameters_cc(const std::string& Label, CARD* Parent,
		       COMMON_COMPONENT* Common, double Value,
		       const node_t& N0, const node_t& N1,
		       ELEMENT* Input);
public: //BUG// for language plugin
  std::vector<std::string>	 _input_label;
  const ELEMENT* _input;
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
#endif
// vim:ts=8:sw=2:noet:
