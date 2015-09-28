/*                               -*- C++ -*-
 * Copyright (C) 2001 Albert Davis
 *           (C) 2014 Felix Salfelder
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
 * inductor base class used by various inductor implementations
 */
#include "e_storag.h"

namespace INDUCTANCE {

class DEV_INDUCTANCE : public STORAGE {
protected:
  explicit DEV_INDUCTANCE(const DEV_INDUCTANCE& p) 
    :STORAGE(p), _c_model(p._c_model) {}
public:
  explicit DEV_INDUCTANCE()
    :STORAGE(), _c_model(false) {}
public: // override virtual
  char	   id_letter()const	{return 'L';}
  std::string value_name()const {return "l";}
  std::string dev_type()const	{return "inductor";}
  uint_t	   max_nodes()const	{return 2;}
  uint_t	   min_nodes()const	{return 2;}
  uint_t	   net_nodes()const	{return 2;}
  uint_t	   int_nodes()const     {return (!_c_model) ? 0 : 1;}
  uint_t	   matrix_nodes()const	{return net_nodes() + int_nodes();}

  bool	   has_inode()const	{return _c_model;}
  bool	   has_iv_probe()const  {return true;}
  bool	   use_obsolete_callback_parse()const {return true;}
  CARD*	   clone()const		{return new DEV_INDUCTANCE(*this);}
  void     expand();
  void	   tr_iwant_matrix();
  void     tr_begin();
  void tr_accept();
  bool	   do_tr();
  void	   tr_load();
  void	   tr_unload();
  double   tr_involts()const	{return tr_outvolts();}
  double   tr_input()const;
  double   tr_involts_limited()const {return tr_outvolts_limited();}
  double   tr_input_limited()const;
  double   tr_amps()const;
  double   tr_probe_num(const std::string&)const;
  void	   ac_iwant_matrix();
  void	   ac_begin()		{_loss1 = _loss0 = ((!_c_model) ? 0. : 1.); _ev = _y[0].f1;}
  void	   do_ac();
  void	   ac_load();
  COMPLEX  ac_involts()const	{return ac_outvolts();}
  COMPLEX  ac_amps()const;

  void set_param_by_name(string Name, string Value);

  std::string port_name(uint_t i)const {itested();
    assert(i != INVALID_NODE);
    assert(i < 2);
    static std::string names[] = {"p", "n"};
    return names[i];
  }
  bool _c_model;
private:
  bool has_ic() const;
};

}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
