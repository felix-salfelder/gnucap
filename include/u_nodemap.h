/*$Id: u_nodemap.h,v 1.1 2009-10-23 12:01:45 felix Exp $ -*- C++ -*-
 * Copyright (C) 2002 Albert Davis
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
 * node name to number mapping -- for named nodes
 */
//testing=script,complete 2006.07.14
#ifndef U_NODEMAP_H
#define U_NODEMAP_H
#include "md.h"
#include "l_istring.h"
/*--------------------------------------------------------------------------*/
class NODE_BASE;
class COMPONENT;
class ADP_NODE; // fixme
class CKT_NODE;
class CARD_LIST;
/*--------------------------------------------------------------------------*/
class NODE_MAP {
private:
  std::map<IString, NODE_BASE*> _node_map;
  explicit  NODE_MAP(const NODE_MAP&);
  unsigned ckt;
  unsigned adp;
public:
  //  NODE_MAP( const NODE_MAP& p) : _node_map(p._node_map) {}
  explicit  NODE_MAP();
	   ~NODE_MAP();
  NODE_BASE*     operator[](IString);
  NODE_BASE*     operator[](unsigned) const;
  CKT_NODE*     new_node(IString, const CARD_LIST* p=0);
  ADP_NODE*     new_adp_node(IString, const COMPONENT* p);
  ADP_NODE*     new_adp_node(IString, const CARD_LIST* p=0);

  typedef std::map<IString, NODE_BASE*>::iterator iterator;
  typedef std::map<IString, NODE_BASE*>::const_iterator const_iterator;

  const_iterator begin()const		{return _node_map.begin();}
  const_iterator end()const		{return _node_map.end();}
  uint_t how_many()const	{assert(_node_map.size()>0); return static_cast<uint_t>(_node_map.size()-1);}
  uint_t how_many_ckt()const	{return ckt;}
  uint_t how_many_adp()const	{return adp;}

};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
#endif
// vim:ts=8:sw=2:noet:
