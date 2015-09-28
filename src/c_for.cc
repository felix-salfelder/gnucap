/*                       -*- C++ -*-
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
 * Step a parameter and repeat a group of commands
 */
#include "c_comand.h"
#include "globals.h"
#include "e_cardlist.h"
#include "u_parameter.h"
#include "u_lang.h"
/*--------------------------------------------------------------------------*/
const int swp_nest_max = 10;
extern int swp_count[], swp_steps[];
extern int swp_type[];
extern int swp_nest;
/*--------------------------------------------------------------------------*/
namespace {
  using std::string;
  static string tempfile = STEPFILE;
  char my_tempfile[swp_nest_max][128];
  vector<string> body[swp_nest_max];
  std::string para_name[swp_nest_max];
  double start[swp_nest_max];
  double last[swp_nest_max];
/*--------------------------------------------------------------------------*/
/* sweep_fix: fix the value for sweep command.
 * (find value by interpolation)
 * if not sweeping, return "start" (the first arg).
 */

  // copy from c_modify
double sweep_fix(CS&, double start, double last)
{
  double value = start;
  if (swp_steps[swp_nest] != 0) {
    double offset = static_cast<double>(swp_count[swp_nest]) 
      / static_cast<double>(swp_steps[swp_nest]);
    if (swp_type[swp_nest]=='L') {
      untested();
      if (start == 0.) {
	untested();
	throw Exception("log sweep can't pass zero");
	value = 0;
      }else{
	untested();
	value = start * pow( (last/start), offset );
      }
    }else{
      value = start + (last-start) * offset;
    }
  }
  return value;
}
/*--------------------------------------------------------------------------*/
static void setup(CS& cmd, CARD_LIST* scope)
{
  for (;;) {
    if (cmd.umatch("li{near} ")) {
      swp_type[swp_nest] = 0;
    }else if (cmd.umatch("lo{g} ")) {
      swp_type[swp_nest] = 'L';
    }else{
      break;
    }
  }
  PARAMETER<double> s, l;
  PARAMETER<uint_t> c;
  cmd >> c;
  cmd >> para_name[swp_nest];
  cmd >> s;
  cmd >> l;

  s.e_val(0., scope);
  l.e_val(1., scope);
  c.e_val(2, scope);

  start[swp_nest] = s;
  last[swp_nest] = l;
  assert (c>0);
  swp_steps[swp_nest] = int(c)-1;


  trace3("got", para_name[swp_nest], start[swp_nest], last[swp_nest] );
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
class CMD_SWEEP : public CMD {
public:
  void do_it(CS& cmd, CARD_LIST* Scope)
  {
    if(swp_nest>=swp_nest_max) { untested();
      error(bDANGER, "nesting too deep\n");
    } else if (cmd.more()) {
      buildfile(cmd, Scope);
      doit(Scope);
      unlink(my_tempfile[swp_nest]);
    }else{ incomplete();
    }
  }
/*--------------------------------------------------------------------------*/
  void doit(CARD_LIST* scope)
  {
    trace2("sweep doit", head, swp_nest);
    FILE *fptr = 0;
    double para_value;

    for (swp_count[swp_nest]=0; swp_count[swp_nest]<=swp_steps[swp_nest];
        ++swp_count[swp_nest]) {
      trace1("for loop", head);
      if (fptr) {
        fclose(fptr);
      }else{
      }
      fptr = fopen(my_tempfile[swp_nest], "r");
      if (!fptr) {
        throw Exception_File_Open("can't open " + string(my_tempfile[swp_nest]));
      }else{
      }
      char buffer[BUFLEN];
      fgets(buffer,BUFLEN,fptr);

      trace1("header check", head);
      CS cmd(CS::_STRING, buffer); //fgets from local file, obsolete
      if (cmd.umatch("sw{eep} ") || cmd.umatch(".sw{eep} ") || (std::string(cmd).c_str())[0] == '*') {
        // setup(cmd);
        trace2("found sweep", cmd, head);
      }else{
        throw Exception("bad file format: " + string(my_tempfile[swp_nest]) + " (" + std::string(cmd) + ")" );
      }

      para_value = sweep_fix( cmd, start[swp_nest], last[swp_nest] );

      trace2("setting", para_name, para_value);
      scope->params()->set( para_name[swp_nest], para_value );

      trace3("doit excuting ", string(my_tempfile[swp_nest]), (long int)ENV::run_mode, head);

      // FIXME? makes interactive mode impossble
      SET_RUN_MODE xx(rBATCH);
      ++swp_nest;
      trace4("running", swp_nest, my_tempfile[swp_nest-1], body[swp_nest-1].size(), head);
      CMD::command(std::string("< ") + string(my_tempfile[swp_nest-1]) , scope);
      --swp_nest;
      fclose(fptr);
      fptr = NULL;
    }
    if(fptr) fclose(fptr);
    fptr = NULL;
    swp_count[swp_nest] = 0;
  }
/*--------------------------------------------------------------------------*/
  void buildfile(CS& cmd, CARD_LIST* scope)
  {
    static FILE *fptr;

    setup(cmd, scope);
    if (fptr) {
      fclose(fptr);
    }else{
    }
    sprintf(my_tempfile[swp_nest], "%s", tempfile.c_str());
    int fd = mkstemp(my_tempfile[swp_nest]);
    fptr = fdopen( fd, "w+");
    if (!fptr) {
      throw Exception_File_Open("can't open temporary file:" + std::string(my_tempfile[swp_nest]));
    }else{
    }

    fprintf(fptr, "* %s\n", cmd.fullstring().c_str());
    unsigned nest = 0;

    for (;;) {
      char buffer[BUFLEN];
      std::string sbuffer;
      trace1("getting things", head);

      switch (ENV::run_mode) {
        case rSCRIPT:
          trace0("in script mode");
          sbuffer = std::string( cmd.get_line("") );
          break;
         case rBATCH:
           trace0("in batch mode");
           sbuffer = std::string( cmd.get_line("") );
           trace1("in batch mode", sbuffer);
 	  break;
        default:
          trace1("not in script mode", ENV::run_mode);
          getcmd(">>>", buffer, BUFLEN);
          sbuffer=buffer;
      }

      if (Umatch(sbuffer.c_str(),".sweep ")) {
	nest++;
      } else if (Umatch(sbuffer.c_str(),"go ")) {
        trace0(("got go " + std::string(my_tempfile[swp_nest])).c_str());
	if (nest) {
	  nest--;
	}else {
	  break;
	}
      }else{
      }
      fprintf(fptr, "%s\n", sbuffer.c_str());
      if( sbuffer.size() && sbuffer[0] == '.') {
	body[swp_nest].push_back(sbuffer);
	//body[swp_nest].push_back(sbuffer.substr(1));
      } else {
	body[swp_nest].push_back(sbuffer);
      }
    }
    fclose(fptr);
    trace2("closed tmp ", my_tempfile[swp_nest], head);
    fptr = NULL;
  }
} p;
DISPATCHER<CMD>::INSTALL d(&command_dispatcher, "sweep| for ", &p);
/*--------------------------------------------------------------------------*/
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet
