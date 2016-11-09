/*                                -*- C++ -*-
 * Copyright (C) 2011 Felix Salfelder
 * Author: Felix Salfelder
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
 */

#include "u_sim_data.h"
#include "e_base.h"

template<>
void BSMATRIX<double>::sink_forward(unsigned* nm){
	assert(nm);
	assert(nm[0]==0);
	unsigned i=0;
	for (uint_t node = 1;  node <= CKT_BASE::_sim->_total_nodes;  ++node) {
		nm[node]=0;
		if (_adj[node].size()){
			nm[node]=++i;
			trace2("BSMATRIX::sink_forward", i, node);
		}

	}
}
template<>
void BSMATRIX<double>::sink_reverse(unsigned* nm){
	assert(nm);
	assert(nm[0]==0);
	unsigned i = 0;
	unsigned t = CKT_BASE::_sim->_total_nodes;
	for (uint_t node = 1;  node <= CKT_BASE::_sim->_total_nodes;  ++node) {
		nm[CKT_BASE::_sim->_total_nodes - node + 1 ] = 0;
	}
	for (uint_t node = 1;  node <= CKT_BASE::_sim->_total_nodes;  ++node) {
		if (_adj[node].size()){
			nm[node] = t- ++i + 1;
			trace2("BSMATRIX::sink_reverse", t-i+1, node);
		}
	}
}


#if 0
#include "s_rcm.h"
void rcm(){
	unsigned n = size() - 1; // numer of nodes. 0 doesnt count

	int* xadj = new int[n];
	signed char* mask = new signed char[n];
	for (unsigned x = 0; x<n; mask[x++]=0);

	int edges=0;

	
	// FIXME: count edges during iwant.
	for(unsigned i=0; i < n; ++i ){
		edges += _adj[i].size();
	}

	int* adj = new int[edges];

	unsigned here = 0;
	for(unsigned i=0; i<n; ++i ){
		unsigned s = _adj[j].size();

		for(unsigned j=0; i<s; ++j ){
			adjx[i]=here; here+=s;
		}
	}
	assert(i == size());
	adjx[i] = here;

	_rcm = new int[n+1];
	_rcm[0] = 0;
	int* deg = new int[n];
	genrcmi(n, 0, xadj, adj, _rcm+1, mask, deg);

	delete deg;
	delete adjx;


}
#endif
