/*$Id: u_lang.cc 2016/09/22 $ -*- C++ -*-
 * Copyright (C) 2006 Albert Davis
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
 */
// testing=script 2015.01.27
#include "globals.h"
#include "c_comand.h"
#include "d_dot.h"
#include "d_coment.h"
#include "e_subckt.h"
#include "e_model.h"
#include "u_lang.h"
/*--------------------------------------------------------------------------*/
LANGUAGE::~LANGUAGE()
{
  if (OPT::language == this) {
    OPT::language = NULL;
  }else{
  }
}
/*--------------------------------------------------------------------------*/
std::string LANGUAGE::getlines(FILE *fileptr) const
{
  // probably missing getlines in your current lang.
  assert(fileptr);
  const int buffer_size = BIGBUFLEN;

  char buffer[buffer_size+1];
  buffer[buffer_size] = '\n';

  char* line = fgets(buffer, buffer_size, fileptr);

  if(line==NULL){
    throw Exception_End_Of_Input("");
  }else if(line[buffer_size]!='\n'){ incomplete();
    // join lines? just truncate for now.
    line[buffer_size] = '\0'; // might be already.
  }else{
    assert(strlen(line)>0);
    line[strlen(line)-1] = '\0';
  }
  return std::string(line);
}
/*--------------------------------------------------------------------------*/
void LANGUAGE::parse_top_item(CS& cmd, CARD_LIST* Scope)
{
  cmd.get_line(I_PROMPT);
  CMD::cmdproc(cmd, Scope);
  trace0("done top item");
}
/*--------------------------------------------------------------------------*/
const CARD* LANGUAGE::find_card(string name, CARD_LIST* Scope, bool nondevice) {
  if (!Scope) Scope = &CARD_LIST::card_list;
  CARD_LIST::const_iterator i = Scope->find_(name);
  if(nondevice){
    while (i!=Scope->end()) {
      if((*i)->is_device()){
        i = Scope->find_again(name, ++i); // skip
      } else {
        break;
      }
    }
  }
  if (i == Scope->end()) {
    throw Exception_Cant_Find(name, "scope");
  }
  return *i;
}
/*--------------------------------------------------------------------------*/
const CARD* LANGUAGE::find_proto(const IString& Name, const CARD* Scope)
{ untested();
  trace2("LANGUAGE::find_proto", Name, Scope);
  const CARD* p = NULL;
  if (Scope) { untested();
    try { untested();
      p = Scope->find_looking_out(Name);
    }catch (Exception_Cant_Find& e) { untested();
      assert(!p);
    }
  }else{ untested();
    CARD_LIST::const_iterator i = CARD_LIST::card_list.find_(Name);
    if (i != CARD_LIST::card_list.end()) { untested();
      p = *i;
    }else{ untested();
      assert(!p);
    }
  }
  
  trace2("lookout", Name, model_dispatcher.size());
  if (p) { untested();
    trace1("LANGUAGE::find_proto found something", prechecked_cast<const COMPONENT*>(p));
    trace1("LANGUAGE::find_proto found something", prechecked_cast<const MODEL_CARD*>(p));
    return p;
  }else if ((command_dispatcher[Name])) { untested();
    trace1("command_dispatcher", Name);
    return new DEV_DOT;	//BUG// memory leak
  }else if ((p = device_dispatcher[Name])) { untested();
    trace0("LANGUAGE::find_proto found device " +Name);
    return p;
  }else if ((p = model_dispatcher[Name])) { untested();
    trace0("LANGUAGE::find_proto found model " +Name);
    return p;
  }else{ untested();
    assert(!p);
    IString s;
    /* */if (Umatch(Name, "b{uild} "))      {untested(); s = "build";}
    else if (Umatch(Name, "del{ete} "))     {untested(); s = "delete";}
    else if (Umatch(Name, "fo{urier} "))    {untested(); s = "fourier";}
    else if (Umatch(Name, "gen{erator} "))  {		 s = "generator";}
    else if (Umatch(Name, "inc{lude} "))    {itested();  s = "include";}
    else if (Umatch(Name, "l{ist} "))       {            s = "list";}
    else if (Umatch(Name, "m{odify} "))     {            s = "modify";}
    else if (Umatch(Name, "opt{ions} "))    {            s = "options";}
    else if (Umatch(Name, "par{ameter} "))  {            s = "param";}
    else if (Umatch(Name, "pr{int} "))      {            s = "print";}
    else if (Umatch(Name, "q{uit} "))       {		 s = "quit";}
    else if (Umatch(Name, "st{atus} "))     {            s = "status";}
    else if (Umatch(Name, "te{mperature} ")){untested(); s = "temperature";}
    else if (Umatch(Name, "tr{ansient} "))  {            s = "transient";}
    else if (Umatch(Name, "tw{otimetran} ")){            s = "twotimetran";}
    else if (Umatch(Name, "ttr "))          {            s = "twotimetran";}
    else if (Umatch(Name, "!"))		    {untested(); s = "system";}
    else if (Umatch(Name, "<"))		    {untested(); s = "<";}
    else if (Umatch(Name, ">"))		    {untested(); s = ">";}
    else{ /* no shortcut available */
      s = Name;
    }
    if ((command_dispatcher[s])) { untested();
      trace1("command_dispatcher", Name);
      return new DEV_DOT; //BUG// we will look it up twice, //BUG// memory leak
    }else{
#ifndef NDEBUG
      trace1("notmodel ", model_dispatcher.size());
      for(DISPATCHER<CARD*>::const_iterator i=model_dispatcher.begin();
	  i!=model_dispatcher.end(); ++i)
      {
	trace2("notmodel ", Name, i->first);
      }
#endif
      return NULL;
    }
  }
}
/*--------------------------------------------------------------------------*/
void LANGUAGE::new__instance(CS& cmd, BASE_SUBCKT* owner, CARD_LIST* Scope)
{
  trace4("LANGUAGE::new__instance", cmd.fullstring(), name(), hp(Scope), owner);

  if (cmd.is_end()) {
    // nothing
  }else{
    IString type = IString(find_type_in_string(cmd));
    if (const CARD* proto = find_proto(type, owner)) {
      trace3("LANGUAGE::new__instance", type, name(), prechecked_cast<const DEV_DOT*>(proto));
      CARD* new_instance = proto->clone_instance();
      trace2("LANGUAGE::new__instance", hp(owner), (owner?int(hp(owner->scope())):0));
      assert(new_instance);
      new_instance->set_owner(owner);
      CARD* x = parse_item(cmd, new_instance);
      if (x) {
	assert(Scope);
	Scope->push_back(x);
        trace3("LANGUAGE::new__instance pushback", x->long_label(), hp(Scope), hp(x->scope()));
      }else{
      }
    }else{
      cmd.warn(bDANGER, type + ": no match");
    }
  }
}
/*--------------------------------------------------------------------------*/
CARD* LANGUAGE::parse_item(CS& cmd, CARD* c)
{
  trace1("LANGUAGE::parse_item", cmd.fullstring());
  // See Stroustrup 15.4.5
  // If you can think of a better way, tell me.
  // It must be in the LANGUAGE class, not CARD.

#if 0 // recheck against upstream
  assert(c);
  if (c->is_device()){
    assert(dynamic_cast<COMPONENT*>(c));
    return parse_instance(cmd, prechecked_cast<COMPONENT*>(c));
  }else if (BASE_SUBCKT* s = dynamic_cast<BASE_SUBCKT*>(c)) {
    return parse_module(cmd, s);
#endif
  if (BASE_SUBCKT* s=dynamic_cast<BASE_SUBCKT*>(c)) {
    if (c->is_device()){
      return parse_instance(cmd, s);
    }else{ untested();
      return parse_module(cmd, s);
    }
  }else if (dynamic_cast<COMPONENT*>(c)) {
    trace0("LANGUAGE::parse_item: COMPONENT");
    return parse_instance(cmd, prechecked_cast<COMPONENT*>(c));
  }else if (MODEL_CARD* m = dynamic_cast<MODEL_CARD*>(c)) {untested();
    return parse_paramset(cmd, m);
  }else if (DEV_COMMENT* com = dynamic_cast<DEV_COMMENT*>(c)) {
    return parse_comment(cmd, com);
  }else if (DEV_DOT* d = dynamic_cast<DEV_DOT*>(c)) {
    return parse_command(cmd, d);
  }else{untested();
    incomplete();
    unreachable();
    return NULL;
  }
}
/*--------------------------------------------------------------------------*/
void LANGUAGE::print_item(OMSTREAM& o, const CARD* c)
{
  // See Stroustrup 15.4.5
  // If you can think of a better way, tell me.
  // It must be in the LANGUAGE class, not CARD.

  assert(c);
  if (c->is_device()){
    assert(dynamic_cast<const COMPONENT*>(c));
    print_instance(o, prechecked_cast<const COMPONENT*>(c));
  }else if (const BASE_SUBCKT* s = dynamic_cast<const BASE_SUBCKT*>(c)) {
    print_module(o, s);
  }else if (const MODEL_CARD* m = dynamic_cast<const MODEL_CARD*>(c)) {
    print_paramset(o, m);
  }else if (const DEV_COMMENT* com = dynamic_cast<const DEV_COMMENT*>(c)) {
    print_comment(o, com);
  }else if (const DEV_DOT* d = dynamic_cast<const DEV_DOT*>(c)) {untested();
    print_command(o, d);
  }else{untested();
    incomplete();
    unreachable();
  }
}
/*--------------------------------------------------------------------------*/
OMSTREAM& operator<<(OMSTREAM& o, LANGUAGE* x)
{
  if (x) {
    return (o << x->name());
  }else{untested();
    return (o << "none");
  }
}
/*--------------------------------------------------------------------------*/
bool Get(CS& cmd, const IString& key, LANGUAGE** val)
{
  if (cmd.umatch(key + " {=}")) {
    LANGUAGE* lang = language_dispatcher[cmd];
    if (lang) {
      *val = lang;
    }else{untested();
      std::string choices;
      for(DISPATCHER<LANGUAGE>::const_iterator
	  i = language_dispatcher.begin(); i != language_dispatcher.end(); ++i) {untested();
	if (i->second) {untested();
	  choices += i->first + ' ';
	}else{untested();
	}
      }
      cmd.warn(bWARNING, "need a language (" + choices + ")");
    }
    return true;
  }else{
    return false;
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:sw=2:ts=8:noet:
