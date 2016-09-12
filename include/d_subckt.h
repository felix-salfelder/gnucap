/*$Id: d_subckt.h,v 1.3 2009-12-13 17:55:01 felix Exp $ -*- C++ -*-
 * vim:ts=8:sw=2:et:
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
 * data structures for subcircuits
 */
//testing=script,sparse 2006.07.17
#ifndef D_SUBCKT_H
#define D_SUBCKT_H
#include "e_node.h"
#include "e_subckt.h"
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
class INTERFACE COMMON_SUBCKT : public COMMON_COMPONENT {
private:
  explicit COMMON_SUBCKT(const COMMON_SUBCKT& p)
    :COMMON_COMPONENT(p), _params(p._params) {++_count;}
public:
  explicit COMMON_SUBCKT(int c=0)	:COMMON_COMPONENT(c) {++_count;}
	   ~COMMON_SUBCKT()		{--_count;}
  bool operator==(const COMMON_COMPONENT&)const;
  COMMON_COMPONENT* clone()const	{return new COMMON_SUBCKT(*this);}
  std::string	name()const		{itested(); return "subckt";}
  static int	count()			{return _count;}

  void set_param_by_name(std::string Name, std::string Value) {_params.set(Name, Value);}
  bool		param_is_printable(int)const;
  std::string	param_name(int)const;
  std::string	param_name(int,int)const;
  std::string	param_value(int)const;
  int param_count()const
	{return (static_cast<int>(_params.size()) + COMMON_COMPONENT::param_count());}

  void		precalc_first(const CARD_LIST*);
  void		precalc_last(const CARD_LIST*);
private:
  static int	_count;
public:
  PARAM_LIST	_params;
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
#endif
// vim:ts=8:sw=2:noet:
