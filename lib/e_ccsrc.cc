/*                         -*- C++ -*-
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
#include "e_ccsrc.h"
/*--------------------------------------------------------------------------*/
void CCSRC_BASE::set_inputs(unsigned node_count, node_t* Nodes)
{
  notstd::copy_n(Nodes, node_count, _n+2);
}
/*--------------------------------------------------------------------------*/
void CCSRC_BASE::expand_last()
{
  ELEMENT::expand_last();
  unsigned bm_ord = 1;

  if (const EVAL_BM_ACTION_BASE* e = dynamic_cast<const EVAL_BM_ACTION_BASE*>(common())) {
    bm_ord = e->input_order();
  }
  _master = input_order() == bm_ord;
  
  trace5("CCSRC_BASE::expand_last", long_label(), input_order(), bm_ord, _input_label.size(), _input_label[0]);

  CARD* master = this;
  for ( unsigned i=0; i<bm_ord-input_order(); i++) {
    master = master->owner();
  }

  if (_input_label.size() && _input_label[0]!="") {
    trace3("looking out", long_label(), _input_label.size(), _input_label[0]);
    _input = dynamic_cast<const ELEMENT*>(master->find_in_my_scope(_input_label[0]));
  }else{untested();
    // _input already set, an internal element.  example: mutual L.
  }

  assert(bm_ord);	
  // if(bm_ord == input_order()) {untested();
  if (_master) {
    for(unsigned i=bm_ord; i>0;) {
      i--;

      assert(master==this || bm_ord!=input_order());
      if (_input_label.size()>i && _input_label[i]!="") {
	_input = dynamic_cast<const ELEMENT*>(master->find_in_my_scope(_input_label[i]));
      }else{ untested();
	trace4("no input label, strange", long_label(), i, _input_label.size(), bm_ord);
	// _input already set, an internal element.  example: mutual L.
      }

      if (!_input) {untested();
	throw Exception(long_label() + ": " + _input_label[i] + " cannot be used as current probe");
      }else if (_input->subckt()) {untested();
	throw Exception(long_label() + ": " + _input_label[i]
	    + " has a subckt, cannot be used as current probe");
      }else if (_input->has_inode()) {untested();
	_n[2*i+IN1] = _input->n_(IN1);
	_n[2*i+IN2].set_to_ground(this);
      }else if (_input->has_iv_probe()) {
	_n[2*i+IN1] = _input->n_(OUT1);
	_n[2*i+IN2] = _input->n_(OUT2);
	trace4("expand", long_label(), i, _n[2*i+IN1].short_label(),  _n[2*i+IN2].short_label());
      }else{
	throw Exception(long_label() + ": " + _input_label[i] + " cannot be used as current probe");
      }
      trace4("CCSRC_BASE::expand_last", long_label(), _n[2*i+IN1].e_(), _n[2*i+IN2].e_(), _net_nodes);

    }
  }

  trace3("CCSRC_BASE::expand_last branch", long_label(), hp(common()), input_order());
  if (input_order()>1) {
    assert(subckt());
    CARD* c = *(subckt()->begin());
    CCSRC_BASE* more = prechecked_cast<CCSRC_BASE*>(c);
    assert(more);

    node_t* nodes = new node_t[input_order()*2];
    nodes[0] = _n[0];
    nodes[1] = _n[1];
    notstd::copy_n(_n+4, input_order()*2-2, nodes+2);
    assert(mutable_common()==common());
    assert(more->input_order() == input_order()-1);
    more->set_parameters("branch", this, mutable_common(), 0., 0, NULL, 2, _n);
    more->set_inputs((input_order()-1)*2, _n+4); // doesnt really help, need device names (?)

    assert(more->common() == common());
    trace2("DEV_CCSRC::expand branch", hp(more->common()), hp(common()));
    subckt()->expand();
    // assert(more->common() == common()); maybe not
    delete[] nodes;
  }else{
  }
}
/*--------------------------------------------------------------------------*/
void CCSRC_BASE::set_port_by_index(uint_t num, std::string& Value)
{
  if (num > 1) {
    _input_label.resize(num-1);
    _input_label[num-2] = Value;
  }else{
    ELEMENT::set_port_by_index(num, Value);
  }
}
/*--------------------------------------------------------------------------*/
bool CCSRC_BASE::node_is_connected(uint_t i)const
{
  if (i > 1) {
    if (_input_label.size() > i-2)
      return _input_label[i-2] != "";
    else
      return false;
  }else{
    return ELEMENT::node_is_connected(i);
  }
}
/*--------------------------------------------------------------------------*/
void CCSRC_BASE::set_parameters_cc(const std::string& Label, CARD *Owner,
			       COMMON_COMPONENT *Common, double Value,
			       const node_t& N0, const node_t& N1,
			       ELEMENT* Input)
{ untested();
  node_t nodes[] = {N0, N1};
  COMPONENT::set_parameters(Label, Owner, Common, Value, 0, 0, 2, nodes);
  _input = Input;
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
