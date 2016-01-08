/*$Id: e_base.cc 2015/02/05 al $ -*- C++ -*-
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
 * Base class for "cards" in the circuit description file
 */
//testing=script 2014.07.04
#include "ap.h"
#include "u_sim_data.h"
#include "m_wave.h"
#include "u_prblst.h"
#include "u_xprobe.h"
#include "e_base.h"
#include <typeinfo>
#include "ap.h"
/*--------------------------------------------------------------------------*/
static char fix_case(char c)
{
  return ((OPT::case_insensitive) ? (static_cast<char>(tolower(c))) : (c));
}
/*--------------------------------------------------------------------------*/
double CKT_BASE::tr_probe_num(const std::string&)const {return NOT_VALID;}
double CKT_BASE::tt_probe_num(const std::string&)const {return NOT_VALID;}
XPROBE CKT_BASE::ac_probe_ext(const std::string&)const {return XPROBE(NOT_VALID, mtNONE);}
XPROBE CKT_BASE::sens_probe_ext(const std::string&)const {
  trace1("CKT_BASE::sens_probe_ext", typeid(*this).name());
  return XPROBE(NOT_VALID, mtNONE);
}
/*--------------------------------------------------------------------------*/
SIM_DATA* CKT_BASE::_sim = NULL;
PROBE_LISTS* CKT_BASE::_probe_lists = NULL;
/*--------------------------------------------------------------------------*/
double	CKT_BASE::tt_behaviour =  0;
double	CKT_BASE::tr_behaviour_del =  0;
double	CKT_BASE::tr_behaviour_rel =  0;
double	CKT_BASE::tt_behaviour_del =  0;
double	CKT_BASE::tt_behaviour_rel =  0;
/*--------------------------------------------------------------------------*/
CKT_BASE::~CKT_BASE()
{
  trace1("~CKT_BASE", _probes);
  if (_probes == 0) {
  }else if (!_probe_lists) {untested();
  }else if (!_sim) {untested();
  }else{
    _probe_lists->purge(this);
  }
}
/*--------------------------------------------------------------------------*/
const std::string CKT_BASE::long_label()const
{
  trace0("CKT_BASE::long_label");
  //incomplete();
  std::string buffer(short_label());
  //for (const CKT_BASE* brh = owner(); exists(brh); brh = brh->owner()) {untested();
  //  buffer += '.' + brh->short_label();
  //}
  return buffer;
}
/*--------------------------------------------------------------------------*/
bool CKT_BASE::help(CS& Cmd, OMSTREAM& Out)const
{
  if (help_text() != "") {
    unsigned here = Cmd.cursor();
    std::string keyword;
    Cmd >> keyword;
    CS ht(CS::_STRING, help_text());
    if (keyword == "") {
      Out << ht.get_to("@@");
    }else if (ht.scan("@@" + keyword + ' ')) {
      Out << ht.get_to("@@");
    }else if (keyword == "?") {
      while (ht.scan("@@")) {
	Out << "  " << ht.get_to("\n") << '\n';
      }
    }else{
      Cmd.warn(bWARNING, here, "no help on subtopic " + Cmd.substr(here));
    }
    return true;
  }else{
    return false;
  }
}
/*--------------------------------------------------------------------------*/
double CKT_BASE::probe_num(const std::string& what)const
{
  trace2("CKT_BASE::probe_num", what, long_label());
  double x;
  if (_sim->analysis_is_tt()){
    x = tt_probe_num(what) ;
  }else  if (_sim->analysis_is_ac()) {
    x = ac_probe_num(what);
  }else  if (_sim->analysis_is_sens()) {
    x = ac_probe_num(what);
  }else{
    x = tr_probe_num(what);
  }
  // FIXME, HACK
  return x; // (std::abs(x)>=1) ? x : floor(x/OPT::floor + .5) * OPT::floor;
}
/*--------------------------------------------------------------------------*/
double CKT_BASE::ac_probe_num(const std::string& what)const
{
  trace1("CKT_BASE::ac_probe_num", what);
  size_t length = what.length();
  mod_t modifier = mtNONE;
  bool want_db = false;
  char parameter[BUFLEN+1];
  strcpy(parameter, what.c_str());

  if (length > 2  &&  Umatch(&parameter[length-2], "db ")) {
    want_db = true;
    length -= 2;
  }
  if (length > 1) { // selects modifier based on last letter of parameter
    switch (fix_case(parameter[length-1])) {
      case 'm': modifier = mtMAG;   length--;	break;
      case 'p': modifier = mtPHASE; length--;	break;
      case 'r': modifier = mtREAL;  length--;	break;
      case 'i': modifier = mtIMAG;  length--;	break;
      default:  modifier = mtNONE;		break;
    }
  }
  parameter[length] = '\0'; // chop
  
  // "p" is "what" with the modifier chopped off.
  // Try that first.
  XPROBE xp(0);
  if (_sim->analysis_is_ac()) {
    xp = XPROBE(ac_probe_ext(parameter));
  } else {
    xp = XPROBE(sens_probe_ext(parameter));
  }

  // If we don't find it, try again with the full string.
  if (!xp.exists()) {
    xp = ac_probe_ext(what);
    if (!xp.exists()) {
      // Still didn't find anything.  Print "??".
    }else{untested();
      // The second attempt worked.
    }
  }
  return xp(modifier, want_db);
}
/*--------------------------------------------------------------------------*/
/*static*/ double CKT_BASE::probe(const CKT_BASE *This, const std::string& what)
{
  if (This) {
    return This->probe_num(what);
  }else{				/* return 0 if doesn't exist */
    return 0.0;				/* happens when optimized models */
  }					/* don't have all parts */
}
/*--------------------------------------------------------------------------*/
/*static*/ WAVE_LIST& CKT_BASE::create_waves(const std::string& coll_name)
{
  assert(_sim);
  _sim->_label = coll_name;
  return _sim->_waves[coll_name];
}
/*--------------------------------------------------------------------------*/
/*static*/ WAVE_LIST* CKT_BASE::find_waves(const std::string& coll_name)
{
  std::map<std::string,WAVE_LIST>::iterator i;
  string n = coll_name;
  if (0 && OPT::case_insensitive) { untested();
    notstd::to_upper(&n);
  }else{
  }
  i = _sim->_waves.find(n);
  if(i!=_sim->_waves.end()){
    return &i->second;
  }else{
    return NULL;
  }
}
/*--------------------------------------------------------------------------*/
/*static*/ WAVE& CKT_BASE::create_wave(const std::string& wave_name, std::string coll_name)
{
  assert(_sim);
  if(coll_name==""){untested();
    coll_name = _sim->_label;
  }else{
  }
  std::string n = wave_name;
  if (OPT::case_insensitive) {
    notstd::to_upper(&n);
  }else{ untested();
  }
  return _sim->_waves[coll_name][n];
}
/*--------------------------------------------------------------------------*/
/*static*/ WAVE* CKT_BASE::find_wave(const std::string& probe_name)
{
  trace2("find_wave", probe_name, _sim->_label);
  std::map<std::string,WAVE>* w = &_sim->_waves[_sim->_label];
  std::map<std::string,WAVE>::iterator i;
  string n = probe_name;
  if (OPT::case_insensitive) {
    notstd::to_upper(&n);
  }else{ untested();
  }
  i = w->find(n);
  if(i!=_sim->_waves[_sim->_label].end()){
    return &i->second;
  }else{
  }
  std::string prefix;
  std::string suffix;
  CS cmd(CS::_STRING, probe_name);
  prefix = cmd.ctos(":");

  if(_sim->_waves.find(prefix) == _sim->_waves.end()) {
    return NULL;
  }else{
  }

  cmd >> ":";
  suffix = cmd.ctos("");
  if (OPT::case_insensitive) {
    notstd::to_upper(&suffix);
  }else{untested();
  }
  trace0(("\"" + prefix + "\":\"" + suffix + "\"").c_str());
  i = _sim->_waves[prefix].find(suffix);
  if(i!=_sim->_waves[prefix].end()){
    return &i->second;
  }

  return NULL;
}
/*--------------------------------------------------------------------------*/
bool CKT_BASE::operator!=(const std::string& n)const {return strcasecmp(_label.c_str(),n.c_str())!=0;}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
