/*$Id: e_aux.h,v 1.1 2009-10-23 12:01:45 felix Exp $ -*- C++ -*-
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
 * helper functions, etc., that sort of belong to circuit elements
 */
//testing=script 2007.07.13
#ifndef E_AUX_H
#define E_AUX_H
#include "e_node.h"
#include "u_status.h"
/*--------------------------------------------------------------------------*/
template <class T>
T port_impedance(const node_t& n1, const node_t& n2,
		 BSMATRIX<T>& mat, const T& parallel)
{
  T* zapit = new T[mat.size()+2];

  for (unsigned ii = 0;  ii < mat.size()+2;  ++ii) {
    zapit[ii] = 0.;
  }
  if (n1.m_() != 0) {
    zapit[n1.m_()] =  1.;
  }else{untested();
  }
  if (n2.m_() != 0) {untested();
    zapit[n2.m_()] = -1.;
  }else{
  }

  mat.fbsub(zapit);
  T raw_z = zapit[n1.m_()] - zapit[n2.m_()];
  delete [] zapit;
  return (parallel != 0.) 
    ? 1. / ((1./raw_z)-parallel)
    : raw_z;
}
/*--------------------------------------------------------------------------*/
inline void set_sens_port(const node_t& n1, const node_t& n2)
{
  std::fill_n(CKT_BASE::_sim->_sens, 1*CKT_BASE::_sim->_total_nodes+1, 0);
  CKT_BASE::_sim->_sens[n1.m_()] = 1;
  if(n2.m_()){
    CKT_BASE::_sim->_sens[n2.m_()] = -1;
  } else {
    // weird, _sens[0] shouldnt do anything
  }
}
/*--------------------------------------------------------------------------*/
inline COMPLEX port_noise(const node_t& n1, const node_t& n2){
  set_sens_port(n1,n2);
  ::status.back.start();
  CKT_BASE::_sim->_acx.fbsubt(CKT_BASE::_sim->_sens);
  ::status.back.stop();
  double a = CARD_LIST::card_list.do_noise();
  return a;
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
#endif
// vim:ts=8:sw=2:noet:
