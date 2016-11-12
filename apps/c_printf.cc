/*$Id: c_exp.cc,v 1.4 2010-09-17 12:25:56 felix Exp $ -*- C++ -*-
 * Copyright (C) 2007 Albert Davis
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
 */
//testing=none
#include "globals.h"
#include "m_expression.h"
#include "c_comand.h"
/*--------------------------------------------------------------------------*/
namespace {
/*--------------------------------------------------------------------------*/
class CMD_ : public CMD {
public:
	void do_it(CS& cmd, CARD_LIST* Scope)
	{
		char s1[2048];
		char* oseek=s1;
		string s;
		cmd >> s;
		OMSTREAM _out = IO::mstdout;
		char* s0 = strdup(s.c_str());
		char* iseek = s0;
		char* percend = strchrnul(iseek, '%');
		strncpy(s1, s0, unsigned(percend-iseek));
		oseek+= percend-iseek;
		iseek = percend;

		while(cmd.more()){
			cmd >> ',';
			char c = cmd.peek();
			if (c=='>') {
				_out.outset(cmd);
				break;
			} else {
			}
			Expression e(cmd);
			Expression r(e, Scope);
			double d=0.;
			try{
				d = r.eval();
			} catch (Exception& e) { untested();
				incomplete();
				throw e;
			}

			percend = strchrnul(percend+1, '%');
			*percend = 0;
			oseek+= sprintf(oseek, iseek, d );
			*percend = '%';
			iseek = percend;
		}
		*oseek = 0;
		cmd.check(bDANGER, "syntax error parsing expression");
		trace1("done", s1);
		_out << string(s1);
		_out.outreset();
		free(s0);
	}
} p0;
DISPATCHER<CMD>::INSTALL d0(&command_dispatcher, "printf", &p0);
/*--------------------------------------------------------------------------*/
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
