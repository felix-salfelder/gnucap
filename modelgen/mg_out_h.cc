/*$Id: mg_out_h.cc,v 1.5 2010-07-14 15:17:30 felix Exp $ -*- C++ -*-
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
 */
//testing=script 2006.11.01
#include "mg_out.h"
/*--------------------------------------------------------------------------*/
static void make_header(std::ofstream& out, const File& in,
			const std::string& dump_name)
{
  out << in.head()
      << "/* This file is automatically generated. DO NOT EDIT */\n"
    "#ifndef " << to_upper(dump_name) << "_H_INCLUDED\n"
    "#define " << to_upper(dump_name) << "_H_INCLUDED\n"
      << in.h_headers() <<
    "#include \"u_sdp.h\"\n"
    "#include \"e_node.h\"\n"
    "#include \"e_subckt.h\"\n"
    "#include \"e_model.h\"\n"
    "#include <boost/assign.hpp>\n"
    "#include <boost/algorithm/string.hpp>\n"
    "/*--------------------------------------"
    "------------------------------------*/\n";
}
/*--------------------------------------------------------------------------*/
static void make_sdp(std::ofstream& out, const Model& m)
{
  out << "class SDP_" << m.name()
      << "\n  :public SDP_" << m.inherit()
      << "{\n"
    "public:\n"
    "  explicit SDP_" << m.name() << "(const COMMON_COMPONENT* c) : SDP_" 
      << m.inherit() << "(c) {init(c);}\n"
    "  void init(const COMMON_COMPONENT*);\n"
    "public:\n";
  Parameter_List::const_iterator p = m.size_dependent().raw().begin();
  for (;;) {
    if (p == m.size_dependent().raw().end()) {
      p = m.size_dependent().calculated().begin();
    }else{
    }
    if (p == m.size_dependent().calculated().end()) {
      break;
    }else{
    }
    out << "  " << (**p).type() << " " << (**p).code_name() 
	<< ";\t// " << (**p).comment() << '\n';
    ++p;
  }
  out << "};\n"
    "/*--------------------------------------"
    "------------------------------------*/\n";
}
/*--------------------------------------------------------------------------*/
static void make_tdp(std::ofstream& out, const Model& m)
{
  out << "class DEV_" << m.dev_type() << ";\n";
  out << "class TDP_" << m.name();
  if (!m.hide_base()) {
    out << "\n  :public TDP_" << m.inherit();
  }else{
  }
  out <<
    "{\n"
    "public:\n"
    "  explicit TDP_"<< m.name() <<"(const DEV_" << m.dev_type() << "*);\n"
    "public:\n";
  for (Parameter_List::const_iterator
       p = m.temperature().calculated().begin();
       p != m.temperature().calculated().end();
       ++p) {
    out << "  " << (**p).type() << " " << (**p).code_name() 
	<< ";\t// " << (**p).comment() << '\n';
  }
  out << "};\n"
    "/*--------------------------------------"
    "------------------------------------*/\n";
}
/*--------------------------------------------------------------------------*/
static void make_model(std::ofstream& out, const Model& m)
{
  std::string class_name = "MODEL_" + m.name().to_string();
  out <<
    "class " << class_name << "\n"
    "  :public MODEL_" << m.inherit() << "{\n"
    "protected:\n"
    "  explicit " << class_name << "(const " << class_name << "& p);\n"
    "public:\n"
    "  explicit " << class_name << "(const BASE_SUBCKT*);\n"
    "  ~" << class_name << "() {--_count;}\n"
    "public: // override virtual\n"
    "  std::string dev_type()const;\n"
    "  void      set_dev_type(const std::string& nt);\n"
    "  CARD*     clone()const {return new " << class_name << "(*this);}\n"
    "  void      precalc_first();\n"
    "  void      precalc_last();\n"
    "  SDP_CARD* new_sdp(COMMON_COMPONENT* c)const;\n"
    "  void      set_param_by_index(int, std::string&, int);\n"
    "  bool      param_is_printable(int)const;\n"
    "  std::string param_name(int)const;\n"
    "  std::string param_name(int,int)const;\n"
    "  std::string param_value(int)const;\n"
    "  int param_count()const {return (" << 1 + m.independent().override().size()
		+ 4 * m.size_dependent().raw().size() + m.independent().raw().size();
  if (!m.hide_base()) {
    out << " + MODEL_" << m.inherit() << "::param_count());}\n";
  }else{
    out << ");}\n";
  }
  if (!m.tt_eval().is_empty()) {
    out << "void tt_eval(COMPONENT*)const;\n";
  }else{
    out << "";
  }
  out <<
    "  bool      is_valid(const COMPONENT*)const;\n"
    "  void      tr_eval(COMPONENT*)const;\n"
    "  virtual void      stress_apply(COMPONENT*)const{ std::cerr<<\"virtual stress apply(C)\\n\" ;}\n"
    "public: // not virtual\n"
    "  static int count() {return _count;}\n"
    "private: // strictly internal\n";
  out <<
    "  static int _count;\n"
    "public: // input parameters\n";
  for (Parameter_List::const_iterator
       p = m.size_dependent().raw().begin();
       p != m.size_dependent().raw().end();
       ++p) {
    out << "  " << "SDP" << " " << (**p).code_name()
	<< ";\t// " << (**p).comment() << '\n';
  }
  for (Parameter_List::const_iterator
       p = m.independent().raw().begin();
       p != m.independent().raw().end();
       ++p) {
    out << "  PARAMETER<" << (**p).type() << "> " << (**p).code_name()
	<< ";\t// " << (**p).comment() << '\n';
  }
  out << "public: // calculated parameters\n";
  for (Parameter_List::const_iterator
       p = m.independent().calculated().begin();
       p != m.independent().calculated().end();
       ++p) {
    out << "  " << (**p).type() << " " << (**p).code_name()
	<< ";\t// " << (**p).comment() << '\n';
  }
  out << "};\n"
    "/*--------------------------------------"
    "------------------------------------*/\n";
}
/*--------------------------------------------------------------------------*/
static void make_common(std::ofstream& out, const Device& d)
{
  std::string class_name = "COMMON_" + d.name().to_string();
  out <<
    "class " << class_name << "\n"
    "  :public COMMON_COMPONENT{\n"
    "public:\n"
    "  explicit " << class_name << "(const " << class_name << "& p);\n"
    "  explicit " << class_name << "(int c=0);\n"
    "           ~" << class_name << "();\n"
    "  bool     operator==(const COMMON_COMPONENT&)const;\n"
    "  COMMON_COMPONENT* clone()const {return new "<<class_name<<"(*this);}\n"
    "  void     set_param_by_index(int, std::string&, int);\n"
    "  void     set_param_by_name(std::string, std::string);\n"
    "  bool     param_is_printable(int)const;\n"
    "  std::string param_name(int)const;\n"
    "  std::string param_name(int,int)const;\n"
    "  std::string param_value(int)const;\n"
    "  int param_count()const {return (" 
	     << d.common().override().size() + d.common().raw().size()
	     << " + COMMON_COMPONENT::param_count());}\n"
    "  void     precalc_first(const CARD_LIST*);\n"
    "  void     expand(const COMPONENT*);\n"
    "  void     precalc_last(const CARD_LIST*);\n"
    "  std::string name()const {itested();return \"" << d.parse_name() << "\";}\n"
    "  const SDP_CARD* sdp()const {return _sdp;}\n"
    "  bool     has_sdp()const {untested();return _sdp;}\n"
    "  static int  count() {return _count;}\n"
    "private: // strictly internal\n"
	 "  static map<std::string, PARA_BASE " << class_name << "::*> param_dict;\n"
    "  static map<std::string, PARA_BASE " << class_name << "::*> param_dict_low;\n"
    "  static int _count;\n"
    "public: // input parameters\n";
  for (Parameter_List::const_iterator
       p = d.common().raw().begin();
       p != d.common().raw().end();
       ++p) {
    out << "  PARAMETER<" << (**p).type() << "> " << (**p).code_name()
	<< ";\t// " << (**p).comment() << '\n';
  }
  out <<
    "public: // calculated parameters\n"
    "  SDP_CARD* _sdp;\n";
  for (Parameter_List::const_iterator
       p = d.common().calculated().begin();
       p != d.common().calculated().end();
       ++p) {
    out << "  " << (**p).type() << " " << (**p).code_name()
	<< ";\t// " << (**p).comment() << '\n';
  }
  out << "public: // attached commons\n";
  for (Args_List::const_iterator
       p = d.circuit().args_list().begin();
       p != d.circuit().args_list().end();
       ++p) {
    out << "  COMMON_COMPONENT* _" << (**p).name() << ";\n";
  }
  out << "};\n"
    "/*--------------------------------------"
    "------------------------------------*/\n";
}
/*--------------------------------------------------------------------------*/
static void make_device(std::ofstream& out, const Device& d)
{
  std::string class_name = "DEV_" + d.name().to_string();
  out <<
    "class " << class_name << " : public BASE_SUBCKT {\n"
    "private:\n"
    "  explicit " << class_name << "(const " << class_name << "& p);\n"
    "public:\n"
    "  explicit " << class_name << "();\n"
    "           ~" << class_name << "() {--_count;}\n"
    "private: // override virtual\n"
    "  char      id_letter()const     {untested();return '" << d.id_letter() << "';}\n"
    "  bool      print_type_in_spice()const {return true;}\n"
    "  std::string value_name()const  {return \"area\";}\n"
    "  //std::string dev_type()const;   //BASE_SUBCKT\n"
    "  uint_t       max_nodes()const     {return " << d.max_nodes() << ";}\n"
    "  uint_t       min_nodes()const     {return " << d.min_nodes() << ";}\n";
  if (d.max_nodes() != d.min_nodes()) {
    out <<
      "  //int     matrix_nodes()const; //BASE_SUBCKT\n"
      "  //int     net_nodes()const;    //BASE_SUBCKT\n";
  }else{
    out <<
      "  //int     matrix_nodes()const; //BASE_SUBCKT\n"
      "  uint_t       net_nodes()const     {return " << d.max_nodes() << ";}\n";
  }
  out << 
    "  uint_t       int_nodes()const     {return " 
      << d.circuit().local_nodes().size() << ";}\n"
    "  CARD*     clone()const         {return new "
      << class_name << "(*this);}\n"
    "  void      precalc_first() {COMPONENT::precalc_first(); if(subckt()) subckt()->precalc_first();}\n"
    "  void      expand();\n"
    "  void      precalc_last()  {COMPONENT::precalc_last(); assert(subckt()); subckt()->precalc_last();}\n"
    "  //void    map_nodes();         //BASE_SUBCKT\n"
    "  //void    tr_begin();          //BASE_SUBCKT\n"
    "  //void    tr_restore();        //BASE_SUBCKT\n";
  if (d.tt_eval().is_empty()) {
	 trace0( "tt_eval isempty" );
    out <<
      "  //void    tt_commit();         //BASE_SUBCKT\n"
      "  //void    tt_prepare();         //BASE_SUBCKT\n";
  }else{
    out <<
      "  void    stress_apply();         //BASE_SUBCKT\n"
      "  void    tt_commit();         //BASE_SUBCKT\n"
      "  void    tt_prepare();         //BASE_SUBCKT\n";
  }
  if (d.tr_eval().is_empty()) {
	 trace0( "tr_eval isempty" );
    out <<
      "  //void    dc_advance();        //BASE_SUBCKT\n"
      "  //void    tr_advance();        //BASE_SUBCKT\n"
      "  //void    tr_regress();        //BASE_SUBCKT\n"
      "  //bool    tr_needs_eval()const;//BASE_SUBCKT\n"
      "  //void    tr_queue_eval();     //BASE_SUBCKT\n"
      "  //bool    do_tr();             //BASE_SUBCKT\n";
  }else{
    out <<
      "  void      dc_advance() {set_not_converged(); BASE_SUBCKT::dc_advance();}\n"
      "  void      tr_advance() {set_not_converged(); BASE_SUBCKT::tr_advance();}\n"
      "  void      tr_regress() {set_not_converged(); BASE_SUBCKT::tr_regress();}\n"
      "  bool      tr_needs_eval()const;\n"
		"  // ????? what is this good for?\n"
      "  void      tr_queue_eval()      {if(tr_needs_eval()){q_eval();}}\n"
      "  bool      do_tr();\n";
  }
  out <<
    "  //void    tr_load();           //BASE_SUBCKT\n"
    "  //double  tr_review();         //BASE_SUBCKT\n"
    "  //void    tr_accept();         //BASE_SUBCKT\n"
    "  //void    tr_unload();         //BASE_SUBCKT\n"
    "  double    tr_probe_num(const std::string&)const;\n"
    "  double    tt_probe_num(const std::string&)const;\n"
    "  //void    ac_begin();          //BASE_SUBCKT\n"
    "  //void    do_ac();             //BASE_SUBCKT\n"
    "  //void    ac_load();           //BASE_SUBCKT\n"
    "  //XPROBE  ac_probe_ext(CS&)const;//CKT_BASE/nothing\n"
    "public:\n"
    "  static int  count() {return _count;}\n"
    "public: // may be used by models\n";
  for (Function_List::const_iterator
       p = d.function_list().begin();
       p != d.function_list().end();
       ++p) {
    out << "  void " << (**p).name() << ";\n";
  }
  out << 
    "private: // not available even to models\n"
    "  static int _count;\n";
  out <<  "public: // input parameters\n";
  for (Parameter_List::const_iterator
       p = d.device().raw().begin();
       p != d.device().raw().end();
       ++p) {untested();
    untested();
    out << "  PARAMETER<" << (**p).type() << "> " << (**p).code_name()
	<< ";\t// " << (**p).comment() << '\n';
  }
  out << "public: // calculated parameters\n";
  for (Parameter_List::const_iterator
       p = d.device().calculated().begin();
       p != d.device().calculated().end();
       ++p) {
    out << "  " << (**p).type() << " " << (**p).code_name()
	<< ";\t// " << (**p).comment() << '\n';
  }
  out << "public: // netlist\n";
  for (Element_List::const_iterator
       p = d.circuit().elements().begin();
       p != d.circuit().elements().end();
       ++p) {
    out << "  COMPONENT* _" << (**p).name() << ";\n";
  }
  out << "private: // node list\n"
    "  enum {";
  for (Port_List::const_iterator
       p = d.circuit().req_nodes().begin();
       p != d.circuit().req_nodes().end();
       ++p) {
    if (p != d.circuit().req_nodes().begin()) {
      out << ", ";
    }else{
    }
    out << "n_" << (**p).name();
  }
  for (Port_List::const_iterator
       p = d.circuit().opt_nodes().begin();
       p != d.circuit().opt_nodes().end();
       ++p) {
    out << ", ";
    out << "n_" << (**p).name();
  }
  for (Port_List::const_iterator
       p = d.circuit().local_nodes().begin();
       p != d.circuit().local_nodes().end();
       ++p) {
    out << ", n_" << (**p).name();
  }
  size_t total_nodes = d.circuit().req_nodes().size() + d.circuit().opt_nodes().size()
    + d.circuit().local_nodes().size();
  out << "};\n"
    "  node_t _nodes[" << total_nodes << "];\n"
    "  std::string port_name(uint_t i)const {\n"
    "    assert(i < " << d.circuit().req_nodes().size() + d.circuit().opt_nodes().size() << ");\n"
    "    static std::string names[] = {";
  for (Port_List::const_iterator
	 p = d.circuit().req_nodes().begin();
       p != d.circuit().req_nodes().end();
       ++p) {
    out << '"' << (**p).name() << "\", ";
  }
  for (Port_List::const_iterator
       p = d.circuit().opt_nodes().begin();
       p != d.circuit().opt_nodes().end();
       ++p) {
    out << '"' << (**p).name() << "\", ";
  }
  out << "\"\"};\n"
    "    return names[i];\n"
    "  }\n"
    "};\n"
    "/*--------------------------------------"
    "------------------------------------*/\n";
}
/*--------------------------------------------------------------------------*/
static void make_eval(std::ofstream& out, const Eval& e,
		      const String_Arg& dev_name)
{
  std::string class_name = "EVAL_" + dev_name.to_string() + '_' 
    + e.name().to_string();
  out <<
    "class " << class_name << " : public COMMON_COMPONENT {\n"
    "private:\n"
    "  explicit "<< class_name << "(const "<< class_name << "& p)\n"
    "    :COMMON_COMPONENT(p) {}\n"
    "public:\n"
    "  explicit "<< class_name << "(int c=0) :COMMON_COMPONENT(c) {}\n"
    "  bool operator==(const COMMON_COMPONENT& x)const "
		"{return COMMON_COMPONENT::operator==(x);}\n"
    "  COMMON_COMPONENT* clone()const {return new "<<class_name<<"(*this);}\n"
    "  std::string name()const {untested(); return \""<< class_name << "\";}\n"
    "  void tr_eval(ELEMENT*d)const;\n"
    "  bool has_tr_eval()const {return true;}\n"
    "  bool has_ac_eval()const {return false;}\n"
    "};\n"
    "/*--------------------------------------"
    "------------------------------------*/\n";
}
/*--------------------------------------------------------------------------*/
static void make_evals(std::ofstream& out, const Device& d)
{
  for (Eval_List::const_iterator
       e = d.eval_list().begin();
       e != d.eval_list().end();
       ++e) {
    make_eval(out, **e, d.name());
  }
}
/*--------------------------------------------------------------------------*/
static void make_tail(std::ofstream& out, const File& in)
{
  out << "// h_direct\n" << in.h_direct() <<
    "/*--------------------------------------"
    "------------------------------------*/\n"
    "/*--------------------------------------"
    "------------------------------------*/\n"
    "#endif\n";
}
/*--------------------------------------------------------------------------*/
void make_h_file(const File& in)
{
  std::string dump_name = in.name();
  { // chop suffix .model
    std::string::size_type loc = dump_name.rfind(".model");
    if (loc == std::string::npos) {
      untested();
      loc = dump_name.rfind(".vams");
    }else{
    }
    if (loc != std::string::npos) {
      dump_name.erase(loc);
    }else{untested();
    }
  }
  { // chop prefix path
    std::string::size_type loc = dump_name.find_last_of(ENDDIR);
    if (loc != std::string::npos) {
      dump_name.erase(0, loc+1);
    }else{itested();
    }
  }

  // open file
  std::ofstream out((dump_name+".h").c_str());
  if (!out) {untested();
    os_error(dump_name);
  }else{
  }

  make_header(out, in, dump_name);

  for (Model_List::const_iterator
       m = in.models().begin();
       m != in.models().end();
       ++m) {
    make_sdp(out, **m);
    make_tdp(out, **m);
    make_model(out, **m);
  }
  for (Device_List::const_iterator
       m = in.devices().begin();
       m != in.devices().end();
       ++m) {
    make_common(out, **m);
    make_evals(out, **m);
    make_device(out, **m);
  }
  make_tail(out, in);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
