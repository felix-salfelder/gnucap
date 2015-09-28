/*                                -*- C++ -*-
 * Copyright (C) 2005 Albert Davis
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
 * A class for parameterized values
 * Used for spice compatible .param statements
 * and passing arguments to models and subcircuits
 */
#ifndef U_PARAMETER_H
#define U_PARAMETER_H
#define PARAM_LIST PARAM_LIST_MAP
#include "md.h"
#include "globals.h"
#include "m_expression.h"
#include "u_opt.h"
#include "e_cardlist.h"
/*--------------------------------------------------------------------------*/
class LANGUAGE;
class PARAM_LIST; //unnnec?
//class MODEL_BUILT_IN_RCD;
//class MODEL_BUILT_IN_BTI;
/*--------------------------------------------------------------------------*/

using ::to_string;
namespace PARM{
  enum NA_ {NA};
}

/*--------------------------------------------------------------------------*/
#define HAVE_PARA_BASE
#define HAVE_PARA_VEC
class PARA_BASE {
protected:
  std::string _s;

public:
  PARA_BASE( ): _s(){}
  PARA_BASE(const PARA_BASE& p): _s(p._s) {}
  PARA_BASE(const std::string s): _s(s){}
  virtual ~PARA_BASE(){}

  bool	has_hard_value()const {return (_s != "");}
  virtual bool	has_good_value()const = 0;

  virtual void	parse(CS& cmd) = 0;
  virtual void	operator=(const std::string& s) = 0;

  void	print(OMSTREAM& o)const		{o << string(*this);}
  void	print(ostream& o)const		{o << string(*this);}
  virtual operator std::string()const = 0;
};
/*--------------------------------------------------------------------------*/
template <class T>
class PARAMETER : public PARA_BASE {
private:
  mutable T _v;
  virtual T my_infty() const; // HACK?
public:
  T _NOT_INPUT() const;

  explicit PARAMETER() : PARA_BASE() {_v=_NOT_INPUT();}
  PARAMETER(const PARAMETER<double>& p): PARA_BASE(p), _v(p._v) {}
  explicit PARAMETER(T v) :PARA_BASE(), _v(v) {}
  //explicit PARAMETER(T v, const std::string& s) :_v(v), _s(s) {untested();}
  ~PARAMETER() {}

  bool	has_good_value()const {return (_v != NOT_INPUT);}
  //bool has_soft_value()const {untested(); return (has_good_value() && !has_hard_value());}

  operator T()const {return _v;}
  // not a good idea:
  // const T* operator&() const {return &_v;}

  const T*	pointer()const	 {return &_v;}
  T	e_val(const T& def, const CARD_LIST* scope, bool try_=false)const;
  void	parse(CS& cmd);

  std::string debugstring()const {
    return(_s + " -> " + to_string(_v));
  }
  virtual std::string string( )const {return *this;}
  operator std::string()const;
  void	print(OMSTREAM& o)const		{o << string();}
  void	print(ostream& o)const		{o << string();}
  void	set_default(const T& v)		{_v = v; _s = "";}
  void	operator=(const PARAMETER& p)	{_v = p._v; _s = p._s;}
  void	operator=(const T& v)		{_v = v; _s = "#";}
  //void	operator=(const std::string& s)	{untested();_s = s;}

  void	operator=(const std::string& s) {
    if (strchr("'\"{", s[0])) {
      CS cmd(CS::_STRING, s);
      _s = cmd.ctos("", "'\"{", "'\"}");
    }else if (s == "NA") {
      _s = "";
    }else{
      _s = s;
    }
  }
  bool  operator==(const PARAMETER& p)const {
    return (_v == p._v  &&  _s == p._s);
  }
  bool  operator==(const T& v)const {
    if (v != NOT_INPUT) {
      return _v == v;
    }else{
      return (_v == NOT_INPUT || !has_hard_value());
    }
  }
  //bool	operator!=(const PARAMETER& p)const {untested();
  //  return !(*this == p);
  //}
  //bool	operator!=(const T& v)const {untested();
  //  return !(*this == v);
  //}
  T*	pointer_hack()	 {return &_v;}
private:
  T lookup_solve(const T& def, const CARD_LIST* scope, bool try_=false)const;
};
/*=========================*/
/*=========================*/
template <>
PARAMETER<double>::PARAMETER(const PARAMETER<double>& p);
//inline std::string PARAMETER<T>::string()const {
//  if (_s == "#") {
//    return to_string(_v);
//  }else if (_s == "") {
//    return "NA(" + to_string(_v) + ")";
//  }else{
//    return _s;
//  }
//}
/*=========================*/
template <class T>
inline  PARAMETER<T>::operator std::string()const {
  if (_s == "#") {
    return to_string(_v);
  }else if (_s == "") {
    return "NA(" + to_string(_v) + ")";
  }else{
    return _s;
  }
}
/*=========================*/
/*--------------------------------------------------------------------------*/
template<>
inline int64_t PARAMETER<int64_t>::_NOT_INPUT() const { return 0 ;} // BUG.
/*--------------------------------------------------------------------------*/
template <>
inline bool PARAMETER<bool>::_NOT_INPUT() const { return true; untested();} // BUG
/*--------------------------------------------------------------------------*/
template <class T>
inline T PARAMETER<T>::_NOT_INPUT() const { return NOT_INPUT; untested();}
/*--------------------------------------------------------------------------*/
template <>
inline int PARAMETER<int>::_NOT_INPUT()const { return NOT_INPUT_INT;} //BUG. magic number?
/*--------------------------------------------------------------------------*/
template <>
inline unsigned int PARAMETER<unsigned int>::_NOT_INPUT() const{
  // stupid();
  return NOT_INPUT_INT;
} //BUG. magic number?
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
template <>
inline std::list<double> PARAMETER< std::list<double> >::_NOT_INPUT() const {
  untested();
  return std::list<double>(0);
}
/*--------------------------------------------------------------------------*/

// ugly hack, probably
template <>
class PARAMETER<string> : public PARA_BASE {
private:
  mutable std::string _v; //value
  std::string my_infty() const;
private:
  std::string lookup_solve(const std::string& def, const CARD_LIST* scope, bool try_=false)const;
public:
  explicit PARAMETER() :_v("") {}
  PARAMETER(const PARAMETER<std::string>& p) : PARA_BASE(p), _v(p._v) {}
  explicit PARAMETER(std::string v) :PARA_BASE(v) {}
  //explicit PARAMETER(T v, const std::string& s) :_v(v), _s(s) {untested();}
  ~PARAMETER() {}

  bool	has_good_value()const {return (true);}

  operator std::string()const {
    trace0(("PARAMETER::string " + _s + " -> " + _v ).c_str());
    return _v;
  }

  std::string e_val(const std::string& def, const CARD_LIST* scope, bool try_=false)const;
  std::string e_val_strange(const std::string& def, const CARD_LIST* scope)const;
  void	parse(CS& cmd);

  std::string value()const;
  std::string string()const;
  void	print(OMSTREAM& o)const	{o << string();}
  void	set_default(const std::string& v)		{_v = v;}
  void	operator=(const PARAMETER& p)	{ _v = p._v; _s = p._s; }
 // argh
  void	operator=(const std::string& s)	{
    if (strchr("'\"{", s[0])) {
      CS cmd(CS::_STRING, s);
      _s = cmd.ctos("", "'\"{", "'\"}");
    }else if (s == "NA") {
      _s = "";
    }else{
      _s = s;
    }
  }
  //void	operator=(const std::string& s)	{untested();_s = s;}

  bool  operator==(const PARAMETER& p)const;
//  {
//    return (_v == p._v && _s == p._s );
//  }
  bool  operator==(const std::string& v)const {
      return _v == v;
  }
  //bool	operator!=(const PARAMETER& p)const {untested();
  //  return !(*this == p);
  //}
  //bool	operator!=(const T& v)const {untested();
  //  return !(*this == v);
  //}
  std::string*	pointer_hack()	 {return &_v;}
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
//template <>
//std::string PARAMETER<std::string>::e_val(const std::string& def, const CARD_LIST* scope)const
/*--------------------------------------------------------------------------*/
/* non-class interface, so non-paramaters can have same syntax */
/* It is needed by the model compiler */
#if 0
inline bool operator==(const PARAMETER<int>& p, double v)
{untested();
  if (v != NOT_INPUT) {untested();
    return p.operator==(static_cast<int>(v));
  }else{untested();
    return (!(p.has_value()));
  }
}
#endif

inline bool has_hard_value(const PARA_BASE& p)
{
  return p.has_hard_value();
}

inline bool has_good_value(const PARA_BASE& p)
{
  return p.has_good_value();
}

#if 0
template <class T>
bool has_soft_value(const PARA_BASE& p)
{untested();
  return p.has_soft_value();
}
#endif

template <class T>
bool has_nz_value(const T& p)
{
  return (has_good_value(p) && p != 0);
}

template <class T>
void set_default(PARAMETER<T>* p, const T& v)
{
  assert(p);
  p->set_default(v);
}

template <class T>
void set_default(T* p, const T& v)
{
  assert(p);
  *p = v;
}

// making readonly pointer transparent.
template <class T>
const T* get_pointer(const PARAMETER<T>& p)
{
  assert(p);
  return p.pointer();
}

template <class T>
const T* get_pointer(const T& p)
{
  assert(p);
  return &p;
}

template <class T>
void e_val(PARAMETER<T>* p, const PARAMETER<T>& def, const CARD_LIST* scope)
{
  assert(p);
  p->e_val(def, scope);
}

template <class T>
void e_val(PARAMETER<T>* p, const T& def, const CARD_LIST* scope)
{
  assert(p);
  p->e_val(def, scope);
}

template <class T>
void e_val(T* p, const T& def, const CARD_LIST* /*scope*/)
{
  assert(p);
  *p=def;
}

template <class T>
void e_val(T* p, const PARAMETER<T>& def, const CARD_LIST*)
{
  assert(p);
  *p=def;
}

//e_val(PARAMETER<MODEL_BUILT_IN_BTI>*, NULL, const CARD_LIST*&)’

#if 0
template <class T>
void e_val(T* p, const T& def, const CARD_LIST*)
{untested();
  assert(p);
  if (*p == NOT_INPUT) {untested();
    *p = def;
  }else{untested();
  }
}
#endif
/*--------------------------------------------------------------------------*/
class PARAM_LIST_COPY;
/*--------------------------------------------------------------------------*/
// base class for param lists
class PARAM_LIST_BASE {
public:
  typedef std::map<const std::string, PARAMETER<double> >::const_iterator
		const_iterator;
  typedef std::map<const std::string, PARAMETER<double> >::iterator
		iterator;
protected:
  PARAM_LIST_BASE(): _try_again(NULL) {}
  PARAM_LIST_BASE* _try_again; // if you don't find it, also look here
protected:
  virtual std::map<const std::string, PARAMETER<double> >& pl() const = 0;
  virtual std::map<const std::string, PARAMETER<double> >& pl() = 0;
public:
  virtual PARAM_LIST_BASE* try_again()const {return _try_again;}
  virtual bool operator==(const PARAM_LIST_BASE& )const {return 0;}
  void set_try_again(PARAM_LIST_BASE* t) {_try_again = t;}
  const PARAMETER<double>& deep_lookup(std::string)const;
  const PARAMETER<double>& operator[](std::string i)const {return deep_lookup(i);}
  PARAMETER<double>& find(std::string); // hack?

  //ddc? PARAMETER<double>& deep_lookup(std::string);
  //ddc? PARAMETER<double>& operator[](std::string i) {return deep_lookup(i);}

  iterator begin() {return pl().begin();}
  const_iterator begin() const {return pl().begin();}
  iterator end() {return pl().end();}
  const_iterator end() const {return pl().end();}

  size_t size()const {return pl().size();}
  bool	 is_empty()const {return pl().empty();}
  bool	 is_printable(int)const;
  std::string name(int)const;
  std::string value(int)const;
};
/*--------------------------------------------------------------------------*/
// actual map
class PARAM_LIST_MAP : public PARAM_LIST_BASE {
private:
  mutable std::map<const std::string, PARAMETER<double> > _pl;
protected:
  std::map<const std::string, PARAMETER<double> >& pl() {return _pl;}
public:
  std::map<const std::string, PARAMETER<double> >& pl() const {return _pl;}
  explicit PARAM_LIST() : PARAM_LIST_BASE() {}
  explicit PARAM_LIST(const PARAM_LIST& p)
				:_pl(p._pl){ _try_again=(p._try_again);};
  //explicit PARAM_LIST(PARAM_LIST* ta) :_try_again(ta) {untested();}
  virtual ~PARAM_LIST() {}
  void	parse(CS& cmd);
  void	print(OMSTREAM&, LANGUAGE*)const;

  size_t size()const {return _pl.size();}
  bool	 is_empty()const {return _pl.empty();}
  bool	 is_printable(int)const;
  std::string name(int)const;
  std::string value(int)const;

  void	eval_copy(PARAM_LIST&, const CARD_LIST*);
  bool  operator==(const PARAM_LIST& p)const {return _pl == p._pl;}
  void	eval_copy(PARAM_LIST_BASE&, const CARD_LIST*);
  bool  operator==(const PARAM_LIST_BASE& p)const {
    const PARAM_LIST_MAP* q = dynamic_cast<const PARAM_LIST_MAP*>(&p);
    return q && _pl == q->_pl;
  }

//  const PARAMETER<std::string>& operator[](std::string i)const {return string_lookup(i);}

  void set(std::string, const std::string&);
  void set(std::string, const double);

public:
  // return a lined up copy of *this
  PARAM_LIST_COPY* copy(PARAM_LIST_BASE* try_again)const;
};
/*--------------------------------------------------------------------------*/
class PARAM_LIST_COPY : public PARAM_LIST_BASE {
public:
  PARAM_LIST_COPY( ) : PARAM_LIST_BASE(), _parent(&_emptyparent) {}
  PARAM_LIST_COPY( const PARAM_LIST_COPY& x ) : PARAM_LIST_BASE(x),
  _parent(x._parent){}
  PARAM_LIST_COPY( const PARAM_LIST* x );
  ~PARAM_LIST_COPY() {}
  virtual std::map<const std::string, PARAMETER<double> >& pl() const
  {
    assert(_parent);
    return _parent->pl();
  }
  PARAM_LIST_COPY& operator=(const PARAM_LIST* c) {_parent=c; return *this;}
protected:
  virtual std::map<const std::string, PARAMETER<double> >& pl()
  {
    assert(_parent);
    return _parent->pl();
  }
private:
  const PARAM_LIST* _parent;
  static PARAM_LIST _emptyparent;
};
/*--------------------------------------------------------------------------*/
inline
PARAM_LIST_COPY::PARAM_LIST_COPY( const PARAM_LIST* x ) : PARAM_LIST_BASE(), _parent(x)
{
  if(!x){
    _parent = &_emptyparent;
  }else{untested();
  }
}
/*--------------------------------------------------------------------------*/
inline string PARAMETER<string>::lookup_solve(const std::string& def,
    const CARD_LIST* scope, bool )const
{
  CS cmd(CS::_STRING, _s);
  trace1("lookup_solve ", _s);
  const PARAM_LIST_MAP* pl = scope->params();
  PARAMETER<double> x =  (pl->deep_lookup(_s));
  return ::string(x);

  return def;
}
/*--------------------------------------------------------------------------*/
// hmm why not lookup bools?!
template <>
inline bool PARAMETER<bool>::lookup_solve(const bool&, const CARD_LIST*, bool)const
{
  CS cmd(CS::_STRING, _s);
  return cmd.ctob();
}
/*--------------------------------------------------------------------------*/
// fixme throw exception in NOT_INPUT case. no magic numbers!!1
template <class T>
inline T PARAMETER<T>::lookup_solve(const T& def, const CARD_LIST* scope, bool try_)const
{
  trace1("PARAMETER<T>::lookup_solve", def);
  CS cmd(CS::_STRING, _s);
  Expression e(cmd);
  Expression reduced(e, scope);
  T v = T(reduced.eval());
  if (v != NOT_INPUT) {
    return v;
  }else{
    trace1("PARAMETER<T>::lookup_solve no luck", _s);
    const PARAM_LIST* pl = scope->params();
    return T(pl->deep_lookup(_s).e_val(def, scope, try_));
  }
}
/*--------------------------------------------------------------------------*/
#if 0
template <class T>
inline T PARAMETER<T>::lookup_solve(const T& def, const CARD_LIST* scope)const
{
  const PARAM_LIST* pl = scope->params();
  return T(pl->deep_lookup(_s).e_val(def, scope));
}
#endif
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
template <>
inline double PARAMETER<double>::my_infty()const{ return inf; }
/*--------------------------------------------------------------------------*/
  enum polarity_t {pP = -1, dunno=0, pN = 1};
template <>
inline polarity_t PARAMETER<polarity_t>::my_infty()const{ return dunno; }
/*--------------------------------------------------------------------------*/
template <class T>
inline T PARAMETER<T>::my_infty()const{ untested(); return 0; }
/*--------------------------------------------------------------------------*/
template <>
inline std::list<double> PARAMETER<std::list<double> >::e_val(const
    std::list<double>& , const CARD_LIST*, bool )const
{
  trace1("PARAMETER<std::list<double> >::e_val", _s);
  double d;

  CS c(CS::_STRING,_s);
  std::list<double>::iterator a;

  if (_v.size()!=0){
    a=_v.begin();
    _v.erase(a);
  }

  while ( c.more() ){
    incomplete();
    d=c.ctof();
    trace1("PARAMETER dl", d);

    _v.push_back( d );
  }
  return _v;
}
/*--------------------------------------------------------------------------*/
// fallback e_val
template <class T>
T PARAMETER<T>::e_val(const T& def, const CARD_LIST* scope, bool try_)const
{
  trace2("PARAMETER<T>::e_val", _s, *this);
  assert(scope);

  static int recursion=0;
  static const std::string* first_name = NULL;
  if (recursion == 0) {
    first_name = &_s;
  }else{
  }
  assert(first_name);

  if (_v == inf) {
    return my_infty();
  }
  if (_s == "inf") {
    return my_infty();
  }
  ++recursion;
  if (_s == "") {
    // blank string means to use default value
    _v = def;
    if (recursion > 1) {
      // temporarily removed...
      if(!try_)
        error(bWARNING, "parameter " + *first_name + " has no value\n");
    }else{
    }
  }else if (_s != "#") {
    // anything else means look up the value
    if (recursion <= OPT::recursion) {
      _v = lookup_solve(def, scope, try_); // FIXME: use try/catch (how?)
      if (_v == _NOT_INPUT()) {
        if(!try_) {
	  error(bDANGER, "parameter " + *first_name + " value is \"NOT_INPUT\"\n");
	}
	//BUG// needs to show scope
	//BUG// it is likely to have a numeric overflow resulting from the bad value
      }else{
      }
    }else{untested();
      _v = def;
      error(bDANGER, "parameter " + *first_name + " recursion too deep\n");
    }
  }else{
    // start with # means we have a final value
  }
  --recursion;
  // trace1("PARAMETER<T>::e_val ", _v);
  return _v;
}
/*--------------------------------------------------------------------------*/
template <>
bool PARAMETER<bool>::e_val(const bool& def, const CARD_LIST* scope, bool try_)const;
/*--------------------------------------------------------------------------*/
template <>
inline void PARAMETER<bool>::parse(CS& cmd)
{
  bool new_val;
  cmd >> new_val;
  if (cmd) {
    _v = new_val;
    _s = "#";
  }else{untested();
    std::string name;
    //cmd >> name;
    name = cmd.ctos(",=();", "'{\"", "'}\"");
    if (cmd) {untested();
      if (name == "NA") {untested();
	_s = "";
      }else{untested();
	_s = name;
      }
    }else{untested();
    }
  }
  if (!cmd) {untested();
    _v = true;
    _s = "#";
  }else{
  }
}
/*--------------------------------------------------------------------------*/
//template <>
//void PARAMETER<vector<PARAMETER<PARAMETER< double > > > >::parse(CS& cmd) ;
/*--------------------------------------------------------------------------*/
template <class T>
inline void PARAMETER<T>::parse(CS& cmd)
{
  trace0(("PARAMETER<T>::parse " + cmd.tail()).c_str());
  T new_val;
  //try
  cmd >> new_val;
  if (cmd) {
    _v = new_val;
    _s = "#";
  //except
  }else{
    std::string name;
    //cmd >> name;
    name = cmd.ctos(",=();", "'{\"", "'}\"");
    if (cmd) {
      if (cmd.match1('(')) {
	_s = name + '(' + cmd.ctos("", "(", ")") + ')';
      }else{
	_s = name;
      }
      if (name == "NA") {
        _s = "";
      }else{
      }
    }else{
    }
  }
}
/*--------------------------------------------------------------------------*/
// hacked e_val. does strange things
inline std::string PARAMETER<std::string>::e_val_strange(const std::string& /*def*/, const CARD_LIST* scope)const
{
  trace0("PARAMETER<std::string>::e_val_strange");
  assert(scope);

  static int recursion=0;
  static const std::string* first_name = NULL;
  if (recursion == 0) {
    first_name = &_s;
  }else{
  }
  assert(first_name);

  ++recursion;
  if (_s == "") {
    // blank string means to use default value
    _v = _s; // where does it come from?
    if (recursion > 1) {
      error(bWARNING, "parameter " + *first_name + " not specified, using default\n");
    }else{
    }
  }else if (_s != "#") {
    // anything else means look up the value
    trace0(("looking up value for "+_s).c_str());
    if (recursion <= OPT::recursion) {
      _v = lookup_solve(_s, scope);
      if (_v == "") {untested();itested();
        error(bDANGER, "parameter " + *first_name + " has no value\n");
      }else{
      }
    }else{untested();
      _v = _s;
      error(bDANGER, "parameter " + *first_name + " recursion too deep\n");
    }
  }else{
    // start with # means we have a final value
  }
  --recursion;

  // ARGH?
  // // cleanup needed
  if (_v=="NA") _v=_s;
  if (_v=="NA( NA)") _v=_s;

  if (_v=="empty") _v="";
  if (_v=="none") _v="";

  trace0(("Evaluated " +_s + " to " + _v).c_str());
  return _v;
}
/*--------------------------------------------------------------------------*/
INTERFACE bool Get(CS& cmd, const std::string& key, PARAMETER<bool>* val);
INTERFACE bool Get(CS& cmd, const std::string& key, PARAMETER<int>* val);
/*--------------------------------------------------------------------------*/
template <class T>
inline ostream& operator<<(ostream& o, const PARAMETER<T> p)
{
  p.print(o);
  return o;
}
/*--------------------------------------------------------------------------*/
template <class T>
inline OMSTREAM& operator<<(OMSTREAM& o, const PARAMETER<T> p)
{
  p.print(o);
  return o;
}
/*--------------------------------------------------------------------------*/
typedef struct{
  double tr_sum; // sum up stress during tr
  double tt_now; // total stress
  double tt_old; // old total
} stress;
/*--------------------------------------------------------------------------*/
template <class T>
inline std::string to_string( PARAMETER<T> n) {
  return string(n);
}

/*--------------------------------------------------------------------------*/

#endif
// a specialization for vectors
// parameters are assigned like
// PARAMETER<vector<PARAMETER<double> > > x;
// x = "1,2,3,4" // without brackets

#include "e_cardlist.h"
#ifndef PV_H
#define PV_H

#include <vector>
#include "md.h"

typedef vector<PARAMETER<double> > dpv;
typedef vector<PARAMETER<dpv> > dpvv;

template <class T>
class PARAMETER<vector<PARAMETER<T> > > : public PARA_BASE{
  private:
    mutable vector<PARAMETER<T> > _v;
    vector<PARAMETER<T> > _NOT_INPUT() const;
  public:
    operator vector<PARAMETER<T> >()const { return _v;}
    explicit PARAMETER(T v) : PARA_BASE("#"), _v(v) {}
    PARAMETER() : PARA_BASE(), _v(vector<PARAMETER<T> >()) {}
    PARAMETER(const PARAMETER<vector<PARAMETER<T> > >& p) :
      PARA_BASE(p), _v(p._v){ }

    //		void	print(OMSTREAM& o)const		{o << string();}
    //		void	print(ostream& o)const		{o << string();}

    std::string string()const;
    //vector<PARAMETER<T> >  _NOT_INPUT() const;
    void	operator=(const std::string& s);
    void	operator=(const PARAMETER<vector<PARAMETER<T> > >& p)	{_v = p._v; _s = p._s;}
    void	operator=(const vector<PARAMETER<T> >& v)		{_v = v; _s = "#";}
    vector<PARAMETER<T> >	e_val(const vector<PARAMETER<T> >& def,
        const CARD_LIST* scope)const;
    std::string to_string(vector< PARAMETER<T> > n) const;

    operator std::string()const;
    size_t size()const{return _v.size();}

    bool has_good_value()const {incomplete(); return false;}
    void parse(CS&) {incomplete();}
};
/*--------------------------------------------------------------------------*/
template <class T>
PARAMETER<vector<PARAMETER<T> > >::operator std::string()const
{
  return string();
}
/*--------------------------------------------------------------------------*/
template <class T>
inline std::string PARAMETER<vector<PARAMETER<T> > >::string()const{
  std::string ret("");
  if (PARAMETER<vector<PARAMETER<T> > >::_s == "#") {
    ret+= "(";
  }else if (_s == "") {
    ret+= "NA(";
  }else{
    return _s;
  }
  for(unsigned  i=0; i<_v.size(); i++){
    ret+= (i)?",":"";
    ret+= std::string(_v[i]);
  }
  ret+=")";
  return ret;
}
/*===============================*/
//typedef PARAMETER<double> dp;
//============================================================
template <class T>
inline vector<PARAMETER<T> > PARAMETER<vector<PARAMETER<T> > >::_NOT_INPUT() const {
  return vector<PARAMETER< T> > ();
}
/*--------------------------------------------------------------------------*/
template<class T>
void PARAMETER<vector<PARAMETER<T> > >::operator=(const std::string& s){
  trace1("PARAMETER dv::operator=" , s);

  CS cmd(CS::_STRING, s);
  _v.clear();
  std::string compon;

  cmd.skipbl();

  for(;;){
    //if (! cmd.umatch("(")){
    //   untested();
    //   break;
    //}
    if(!cmd.more()) break;
    compon = cmd.ctos(",","({",")}",""); // More sorts of braces?

    PARAMETER<T> d;

    //d =  '(' + compon + ')';
    d = compon;
    _v.push_back( d );

    cmd.skip1b(')');
    trace1("PARAMETER vector loop pushed back ", d);
  }
  _s = "#";
  trace2("PARAMETER done vector loop", cmd.tail(), *this);
}
// FIXME: templatify
// template<class T>
//string to_string(vector<PARAMETER<double> > n);
//string to_string(vector< PARAMETER< vector< PARAMETER<double> > > > n);
/*--------------------------------------------------------------------------*/


//#include "u_parameter.h"
template<>
inline CS&     CS::operator>>(vector<PARAMETER<vector<PARAMETER<double> > > >& x);
//----------------------------------------------------------------
/*--------------------------------------------------------------------------*/
  template<class T>
inline string to_string(vector<PARAMETER<T> > n)
{

  return n; // use operator.
  //  string buf("");
  //  // FIXME: remove one ,
  //  if (n.size()==0){return "( )";}
  //
  //  vector<PARAMETER<T> >::iterator i=n.begin();
  //  buf += string("(")+ftos((double)*i, 0, 7, 0);
  //  ++i;
  //
  //  while (i!=n.end()){
  //    buf += std::string(", ") + ftos((double)*i, 0, 7, 0);
  //    ++i;
  //  }
  //  return buf + std::string(" )");;
}
//========================================================================
//template<class T>
//inline vector<PARAMETER<T> > PARAMETER<vector<PARAMETER<T> > >::e_val(
//    const vector<PARAMETER<T> >& def, const CARD_LIST* scope)const{
//
//  trace2("dp::e_val", _s, _v.size());
//  assert(def.size() >= _v.size());
//  assert(_s=="" ||  _v.size());
//
//  for(unsigned  i=0; i<_v.size()  ; i++){
//    _v[i].e_val(def[i], scope);
//  }
//
//  return _v;
//
//}
/*-----------------------------------*/
//template<class T>
//CS&     CS::operator>>(vector<PARAMETER<vector<PARAMETER<double> > > >& ){
//
//  trace0("CS::operator");
//  incomplete();
//  return *this;
//}
/*-----------------------------------*/
template <class T>
inline vector<PARAMETER<T> >
PARAMETER<vector<PARAMETER<T> > >::e_val(const vector<PARAMETER<T> >& def, const CARD_LIST* scope)const
{
  trace2("PARAMETER dv::e_val", _s, _v.size());
  trace1("PARAMETER dv::e_val", (std::string)(*this));
  trace1("PARAMETER dv::e_val", def.size() );

  //  CS c(CS::_STRING,_s);
  // FIXME: accept strings and parse...
  for(unsigned  i=0; i<_v.size()  ; i++){
    PARAMETER<T> D;
    if (i < def.size()){
      D = def[i];
    }
    trace3("PARAMETER vector eval", i, _v[i], D);
    _v[i].e_val(D, scope);
    trace2("PARAMETER vector eval", i, _v[i]);
  }
  return _v;
}
/*--------------------------------------------------------------------------*/
template <class S>
inline S& operator<<( S& o, const vector<PARAMETER<double> >  &m)
{
  o << "(";

  for ( vector<PARAMETER<double> >::const_iterator ci=m.begin();
      ci!=m.end();)
  {
    o << " " << *(ci) << " ";
    ++ci;
  }
  o << ")";
  return o;
}
/*-----------------------------------*/
#endif
// vim:ts=8:sw=2:noet:
