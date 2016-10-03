/*                            -*- C++ -*-
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
 * tr,dc analysis output functions (and some ac)
 */
//testing=obsolete,script 2005.09.17
#include "u_sim_data.h"
#include "u_status.h"
#include "m_wave.h"
#include "u_prblst.h"
#include "declare.h"	/* plottr, plopen */
#include "s__.h"
/*--------------------------------------------------------------------------*/
/* SIM::____list: access probe lists
 */
const PROBELIST& SIM::alarmlist()const
{
  return _probe_lists->alarm[_sim->_mode];
}
const PROBELIST& SIM::plotlist()const
{
  return _probe_lists->plot[_sim->_mode];
}
const PROBELIST& SIM::printlist()const
{
  trace2("SIM::printlist", _sim->_mode, _probe_lists->print[_sim->_mode].size() );
  return _probe_lists->print[_sim->_mode];
}
const PROBELIST& SIM::storelist()const
{
  return _probe_lists->store[_sim->_mode];
}
/*--------------------------------------------------------------------------*/
void SIM::expect_results(double t){
  CS* c = _sim->expect_file();

  std::vector<double> last_raw(_sim->_expect_raw);

  if (!c) return;
  trace2("SIM::expect_results",t,  c->fullstring() );


  double d;

  double et = -1;

  while(  et-t < -1e-20 ){ // FIXME rounding bug.
    trace2("SIM::expect_results loop", et, t);
    c->get_line(I_PROMPT);
    _sim->_expect_raw.clear();

    while (c->umatch("'|*|#|//|\"")){
      c->get_line(I_PROMPT);
      c->eat_lines();
    }
    c->eat_lines();

    if(!c->more()){
      incomplete();
      assert(false);
    }

    *c >> et;
    trace1("SIM::expect_results time", et);
    assert(( et==0) == (t==0));
    _sim->_expect_raw.push_back(et);

    while( c->more()){
      *c >> d;
      trace1("SIM::expect_results found", d);
      _sim->_expect_raw.push_back(d);

      if (_sim->_expect_raw.size() > 100) assert(false); // safety net.
    }

    trace0("need more input");
    last_raw = _sim->_expect_raw;

    {

      //interpolate
      //
    }


  }


  trace1("SIM::expect_results size", _sim->_expect_raw.size());

}
/*--------------------------------------------------------------------------*/
/* SIM::out: output the data, "keep" for ac reference
 * x = the x coordinate
 * print = selected points, "print" to screen, files, etc.
 * store = all points, for internal postprocessing, measure
 * keep = after the command is done, dcop for ac
 */
void SIM::outdata(double x, int outflags)
{

  trace1("SIM::outdata ", x);

  ::status.output.start();
  if (outflags & ofKEEP) {
    _sim->keep_voltages();
  }else{
  }
  if (outflags & ofPRINT) {
    plottr(x, plotlist());
    print_results(x);
    _sim->reset_iteration_counter(iPRINTSTEP);
    ::status.hidden_steps = 0;
  }else{
    ++::status.hidden_steps;
  }
  if (outflags & ofSTORE) {
    alarm();
    store_results(x);
  }else{
  }
  ::status.output.stop();
}
/*--------------------------------------------------------------------------*/
/* SIM::head: print column headings and draw plot borders
 */
void SIM::head(double start, double stop, const std::string& col1)
{
  _sim->_waves[_sim->_label].clear();

  if(_wavep){
    delete[] _wavep;
  }else{
  }
  _wavep = new WAVE*[storelist().size()];

  unsigned ii = 0;
  for (PROBELIST::const_iterator
	 p=storelist().begin();  p!=storelist().end();  ++p) {
    string l = (*p)->label();
    if (OPT::case_insensitive) {
      notstd::to_upper(&l);
    }
    WAVE* newwave = &(_sim->_waves[_sim->_label][l]);

    // UGLY: linear.
    for(unsigned j = 0; j<ii; ++j){
      if(_wavep[j] == newwave){
	newwave = NULL;
	error(bPICKY,"cowardly refusing to register twice: " + (*p)->label() + " at index " + to_string(ii) + "\n");
	break;
      }
    }
    _wavep[ii++] = newwave;
  }

  if (!plopen(start, stop, plotlist())) {
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int width = std::min(OPT::numdgt+5, BIGBUFLEN-10);
    char format[20];
    //sprintf(format, "%%c%%-%u.%us", width, width);
    sprintf(format, "%%c%%-%us", width);

    _out.form(format, '#', col1.c_str());

    for (PROBELIST::const_iterator
	   p=printlist().begin();  p!=printlist().end();  ++p) {
      _out.form(format, ' ', (*p)->label().c_str());
    }
    //    _out <<  printlist().size() << '\n';
    _out << '\n';
  }else{
  }
}
/*--------------------------------------------------------------------------*/
/* SIM::print_results: print the list of results (text form) to _out
 * The argument is the first column (independent variable, aka "x")
 */
void SIM::print_results(double x)
{
  trace2("SIM::print_results", printlist().size(), x);
  unsigned i=0;

  //if (_sim->expect_file()){
  //  assert( fabs (_sim->_expect_raw[i] - x) < 1e-20);
  //}

  if (!IO::plotout.any()) {
    _out.setfloatwidth(OPT::numdgt, OPT::numdgt+6);
    assert(x != NOT_VALID);
    _out << x;
    for (PROBELIST::const_iterator
	   p=printlist().begin();  p!=printlist().end();  ++p) {
      trace1("SIM::print_results", (*p)->label());
      double v= (*p)->value(); 
      if(_sim->expect_file()){ 
        i++;
        trace3("SIM::print_results", v,_sim->_expect_raw[i],i);
        trace1("SIM::print_results", fabs(v-_sim->_expect_raw[i]));
          
        assert( fabs(v-_sim->_expect_raw[i]) < pow(10., 1-OPT::numdgt) );
      }
      _out << v;
    }
    _out << '\n';
  }else{
  }
  trace0("SIM::print_results done");

}
/*--------------------------------------------------------------------------*/
/* SIM::alarm: print a message when a probe is out of range
 */
void SIM::alarm(void)
{
  _out.setfloatwidth(OPT::numdgt, OPT::numdgt+6);
  for (PROBELIST::const_iterator
	 p=alarmlist().begin();  p!=alarmlist().end();  ++p) {
    if (!(*p)->in_range()) {
      _out << (*p)->label() << '=' << (*p)->value() << '\n';
    }else{
    }
  }
}
/*--------------------------------------------------------------------------*/
/* SIM::store: store data in preparation for post processing
 */
void SIM::store_results(double x)
{
  int ii = 0;
  for (PROBELIST::const_iterator
	 p=storelist().begin();  p!=storelist().end();  ++p) {
    if(_wavep[ii]){
      _wavep[ii]->push(x, (*p)->value());
    }else{
      // unneeded.
    }
    ++ii;
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
