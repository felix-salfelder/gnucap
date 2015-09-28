/*                          -*- C++ -*-
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
 * "wave" class, for transmission lines and delays
 */
//testing=script 2006.07.13
#include "l_denoise.h"
#include "m_interp.h"
#include "io_misc.h"
#include <boost/range/iterator_range.hpp>
/*--------------------------------------------------------------------------*/
class WAVE {
private:
  std::deque<DPAIR> _w;
  double _delay;
public:
  typedef std::deque<DPAIR>::iterator iterator;
  typedef std::deque<DPAIR>::const_iterator const_iterator;

  explicit WAVE(double d=0);
  explicit WAVE(const WAVE&);
	  ~WAVE() {}
  WAVE&	   set_delay(double d);
  WAVE&	   initialize();
  WAVE&	   push(double t, double v);
  FPOLY1   v_out(double t)const;
  FPOLY1   operator()(double xx)const {return v_out(xx);}
  double   v_reflect(double t, double v_total)const;
  WAVE&	   operator+=(const WAVE& x);
  WAVE&	   operator+=(double x);
  WAVE&	   operator*=(const WAVE& x);
  WAVE&	   operator*=(double x);
  WAVE&	   warp(double x);
  double   dhd_linear(const WAVE& x, std::pair<DPAIR, DPAIR>* where=NULL) const;
  double   dhd_discrete(const WAVE& x, std::pair<DPAIR, DPAIR>* where=NULL) const;
  const_iterator begin()const {return _w.begin();}
  const_iterator end()const {return _w.end();}
  size_t size()const {return _w.size();}

private:
  double segment2point(const_iterator l, const DPAIR& p, DPAIR* where=NULL) const;
  double norm2(const DPAIR&, const DPAIR&) const;
  double norm2(const DPAIR&) const;
  double point2chain(DPAIR p, const_iterator B, DPAIR* w=NULL) const; // FIXME: pass right end
  double edge2chain(DPAIR l, DPAIR r, const_iterator& L, const_iterator& R, std::pair<DPAIR, DPAIR>* where=NULL) const;
  FPOLY1 bisector_dual(const const_iterator& a, const const_iterator& b) const;
  double bisector_distance(const const_iterator& a, const const_iterator& b, DPAIR l, DPAIR r) const;
  double endpoints2chain(DPAIR l, DPAIR r, const_iterator& L, const_iterator& R, std::pair<DPAIR, DPAIR>* where=NULL) const;
private:
  class RANGE {
  private:
    const_iterator _a;
    const_iterator _b;
  public:
    RANGE(const const_iterator& a, const const_iterator& b) : _a(a), _b(b){ untested();}
    const_iterator begin() const {return _a;}
    const_iterator end() const {return _b;}
  };
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// push: insert a signal on the "input" end.
// args: t = the time now
//       v = the value to push
//
inline WAVE& WAVE::push(double t, double v)
{
  _w.push_back(DPAIR(t+_delay, v));
  return *this;
}
/*--------------------------------------------------------------------------*/
// initialize: remove all info, fill it with all 0.
//
inline WAVE& WAVE::initialize()
{
  _w.clear();
  return *this;
}
/*--------------------------------------------------------------------------*/
inline WAVE::WAVE(const WAVE& w)
  :_w(w._w),
   _delay(w._delay)
{
}
/*--------------------------------------------------------------------------*/
// constructor -- argument is the delay
//
inline WAVE::WAVE(double d)
  :_w(),
   _delay(d)
{
  initialize();
}
/*--------------------------------------------------------------------------*/
inline WAVE& WAVE::set_delay(double d) 
{
  _delay = d; 
  return *this;
}
/*--------------------------------------------------------------------------*/
// v_out: return the value at the "output" end
// args: t = the time now
//
inline FPOLY1 WAVE::v_out(double t)const
{
  return interpolate(_w.begin(), _w.end(), t, 0., 0.);
}
/*--------------------------------------------------------------------------*/
// reflect: calculate a reflection
// args: t = the time now
//       v_total = actual voltage across the termination
// returns: the value (voltage) to send back as the reflection
//
inline double WAVE::v_reflect(double t, double v_total)const
{
  // return (v_total*2 - v_out(t)); // de-noised
  return dn_diff(v_total*2, v_out(t).f0);
}
/*--------------------------------------------------------------------------*/
inline WAVE& WAVE::operator+=(const WAVE& x)
{
  untested();
  for (std::deque<DPAIR>::iterator
	 i = _w.begin(); i != _w.end(); ++i) {
    untested();
    (*i).second += x.v_out((*i).first).f0;
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
inline WAVE& WAVE::operator+=(double x)
{
  untested();
  for (std::deque<DPAIR>::iterator
	 i = _w.begin(); i != _w.end(); ++i) {
    untested();
    (*i).second += x;
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
inline WAVE& WAVE::operator*=(const WAVE& x)
{
  untested();
  for (std::deque<DPAIR>::iterator
	 i = _w.begin(); i != _w.end(); ++i) {
    untested();
    (*i).second *= x.v_out((*i).first).f0;
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
inline WAVE& WAVE::operator*=(double x)
{
  for (std::deque<DPAIR>::iterator
	 i = _w.begin(); i != _w.end(); ++i) {
    (*i).second *= x;
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
inline WAVE& WAVE::warp(double t)
{
  assert(t>0);
  for (std::deque<DPAIR>::iterator
	 i = _w.begin(); i != _w.end(); ++i) {
    (*i).first *= t;
  }
  return *this;
}
/*--------------------------------------------------------------------------*/
inline double WAVE::norm2(const DPAIR& a) const
{
  double sum = (a.first) * (a.first);
  assert(is_number(sum));
  sum += (a.second) * (a.second);
  assert(is_number(sum));
  return sqrt(sum);
}
/*--------------------------------------------------------------------------*/
inline double WAVE::norm2(const DPAIR& a, const DPAIR& b) const
{
  double sum = (a.first - b.first) * (a.first - b.first);
  assert(is_number(sum));
  sum += (a.second - b.second) * (a.second - b.second);
  assert(is_number(sum));
  return sqrt(sum);
}
/*--------------------------------------------------------------------------*/
inline DPAIR operator+(const DPAIR& a, const DPAIR& b)
{
  DPAIR sum(a);
  sum.first += b.first;
  sum.second += b.second;
  return sum;
}
/*--------------------------------------------------------------------------*/
inline DPAIR operator-(const DPAIR& a, const DPAIR& b)
{
  DPAIR sub(a);
  sub.first -= b.first;
  sub.second -= b.second;
  return sub;
}
/*--------------------------------------------------------------------------*/
inline DPAIR operator*(const double& t, const DPAIR& a)
{
  DPAIR s(a);
  s.first *= t;
  s.second *= t;
  return s;
}
#if __cplusplus > 199711L
// compute segment to point distance
// take left endpoint of segment
inline double WAVE::segment2point(const_iterator li, const DPAIR& p, DPAIR* where) const
{
  WAVE::const_iterator ri = next(li);
  const DPAIR l = *li;
  const DPAIR r = *ri;
//  trace3("norm2", l, r, p);
  double ret = 0;

  if(li == end()){ unreachable();
  }else if(ri == end()){
    ret = norm2(l,p);
    if(where) {
      *where = *li;
    }
  }else{
    assert(l != r);
    if(l.first >= r.first){
      error(bDANGER, "something wrong %E %E\n", l.first, r.first);
    }
    assert(l.first < r.first);
    // solve l + t*(r-l) closest to p
    double t;
    double p1 = p.first;
    double p2 = p.second;
    double l1 = l.first;
    double l2 = l.second;
    double r1 = r.first;
    double r2 = r.second;
    t = - (p1 - l1) * (l1 - r1) - (p2 - l2) * (l2 - r2);
    assert(is_number(t));
    t /= (l1 - r1) * (l1 - r1) + (l2 - r2) * (l2 - r2);
    assert(is_number(t));

    DPAIR op;
    if(t<=0) {
      op = l;
    } else if(t>=1) {
      op = r;
    }else{
      op = l+(t*(r-l));
    }
    ret = norm2(op,p);
    if(where){
      *where = op;
    }
//    trace2("norm2 found", t, ret);
  }
  return ret;
}
/*--------------------------------------------------------------------------*/
// compute directed distance from point to chain
// distance is min_{x \in chain} |b-x|
// chain is iterator of this starting at B
// edge has vertices l and r
// FIXME: pass right end.
inline double WAVE::point2chain(DPAIR p, const_iterator L, DPAIR* w) const
{
  assert(L!=end());

  double ret = norm2(*L,p);
  if(w){
    *w = *L;
    trace2("point2chain", ret, *w);
  }

  while(next(L) != end() && L->first < p.first+ret){
    trace2("point2chain", *L, *next(L));
    assert(next(L)!=end());
    DPAIR w2;
    DPAIR* pw2 = (w)? &w2 : NULL;
    double distance = segment2point(L,p,pw2);
    trace4("point2chain", p, distance, *L, w2);
    if(distance<ret){
      ret = distance;
      if (w){
	*w = w2;
      }
    }
    ++L;
  }

  trace4("point2chain", p, *L, *next(L), ret);
  if(w) trace1("point2chain", *w);
  return(ret);
}
/*--------------------------------------------------------------------------*/
// compute bisector of 2 edges
inline FPOLY1 WAVE::bisector_dual(const const_iterator& a, const const_iterator& b) const
{

  FPOLY1 fa(a->first, a->second, (next(a)->second - a->second)/(next(a)->first - a->first));
  FPOLY1 fb(b->first, b->second, (next(b)->second - b->second)/(next(b)->first - b->first));

//  trace2("bisector_dual", *a, *next(a));
//  trace2("bisector_dual", *b, *next(b));
//  trace2("bisector_dual", fa, fb);
//  trace2("bisector_dual", CPOLY1(fa), CPOLY1(fb));

  DPAIR t = CPOLY1(fa).intersect(CPOLY1(fb));
//   trace1("bisector_dual", t);

  double norma = 1./norm2(*next(a) - *a);
  DPAIR A = norma*(*next(a) - *a);
  double normb = 1./norm2(*next(b) - *b);
  DPAIR B = normb*(*next(b) - *b);
  DPAIR sum = A+B;

  return FPOLY1(t.first, t.second, sum.second/sum.first );
}
/*--------------------------------------------------------------------------*/
// find distance from endpoints to chain
// update chain iterators to relevant region for inner points
inline double WAVE::endpoints2chain(DPAIR l, DPAIR r, const_iterator& L, const_iterator& R, std::pair<DPAIR, DPAIR>* where) const
{
  trace4("endpoints2chain start....", l, r, *L, *R);
  if (L == R){
    // empty range
  }
  assert(R != end()); // bad right end.
  auto Lsave = L; // mark start, so we don't sweep twice
  auto Rsave = R;
  auto Lmin = L; // found minimal distance here (points to left vertex of edge)
  DPAIR hmmmL, hmmmR; // FIXME: merge. later.
  DPAIR* hmmmLp = NULL;
  DPAIR* hmmmRp = NULL;
  if(where){
    hmmmL = *L;
    hmmmR = *R;
    hmmmLp = &hmmmL;
    hmmmRp = &hmmmR;
  }
  double distL = norm2(*L, l);
  hmmmL = *L;

  // shift L to left if necessary
  while(L!=begin()){
    double tmp = segment2point(prev(L), l, hmmmLp);
    if(tmp < distL){ untested();
      distL = tmp;
      if(where){ untested();
	where->second = hmmmL;
      }
      Lmin = L;
      --L;
    } else if( prev(L)->first > l.first - distL) {
      --L;
    }else{
      break;
    }
  }

  trace4("endpoints2chain first L pass", *Lsave, *Lmin, distL, hmmmL);
  // shift L to the right as much as possible
  while(next(Lsave)!=end()){
    trace3("L ->", *Lsave, *next(Lsave), *R);
    double tmp = segment2point(Lsave, l, hmmmLp);
    if(tmp < distL){
      distL = tmp;
      if(where){
	where->second = hmmmL;
      }
      Lmin = Lsave;
      ++Lsave;
    } else if( Lsave->first < l.first + distL) {
      ++Lsave;
    }else{
      break;
    }
  }

  trace4("endpoints2chain second L pass", *Lsave, *Lmin, distL, hmmmL);
  L = Lmin;

  double distR = norm2(*R, r);
  if (distR>distL) {
    if(where){
      where->second = *R;
    }
  }
  auto Rmin = R; // found minimal distance here (points to right vertex of edge)

  trace2("guess", distR, *R);
  // shift R to right if neccessary
  while(next(R)!=end()){
    trace3("R ->", *R, r, *prev(R));
    double tmp = segment2point(R, r, hmmmRp); // check edge right to range
    if(tmp < distR){
      distR = tmp;
      if(where && (distR > distL)){
	trace1("endpoints2chain R->", hmmmR);
	where->second = hmmmR;
      }else if(where){
	trace3("endpoints2chain R-> nosecond", hmmmR, distL, distR);
      }
      ++R;
      Rmin = R;
    } else if( R->first < r.first + distR) {
      ++R;
    }else{
      break;
    }
  }

  trace6("endpoints2chain done first R pass", *Rsave, *Rmin, *next(Rmin), distR, *R, hmmmR);


  // shift R to the left as much as possible
  while(Rsave != Lmin){
    assert(Rsave!=begin());
    double tmp = segment2point(prev(Rsave), r, hmmmRp);
    if(tmp < distR){
      Rmin = Rsave;
      --Rsave;
      distR = tmp;
      if(where && (distR > distL)){
	trace1("endpoints2chain <-R", hmmmR);
	where->second = hmmmR;
      }else if(where){
	trace3("endpoints2chain <-R nosecond", hmmmR, distL, distR);
      }
    } else if( prev(Rsave)->first < r.first + distR) {
      --Rsave;
    }else{ untested();
      break;
    }
  }
  trace6("endpoints2chain done second R pass", *Rsave, *Rmin, *next(Rmin), distR, *R, hmmmR);

  R = Rmin;

  trace4("endpoints2chain check", *L, *R, distL, distR);
  assert(L->first <= R->first);
  assert(L != end());

//  auto x = L;
//  for(; x!=R; ++x){
//    assert(x!=end());
//  }
//  assert(x==R);

  trace8("endpoints2chain done", l, r, *Lmin, *Rmin, distL, distR, *L, *R);
  if (distL > distR){
    if(where){
      where->first = l;
      //where->second = hmmmL;
      trace1("endpoints2chain done l", where->second);
    }
    return distL;
  }else{
    if(where){
      where->first = r;
      trace1("endpoints2chain done r", where->second);
      //where->second = hmmmR;
    }
    return distR;
  }
}
/*--------------------------------------------------------------------------*/
// compute worst case distance from edge to chain, leaving out right (?) endpoint
// chain is subsegment of this incluging L .. R
// edge has vertices l and r
// B will be incremented as needed
inline double WAVE::edge2chain(DPAIR l, DPAIR r, const_iterator& L, const_iterator& R, std::pair<DPAIR, DPAIR>* where) const
{
  double ret = 0;
  trace4("edge2chain...", l, r, *L, *R);

  // need to check (left) endpoint of edge against range
  if(l == r){ incomplete();
    return point2chain(l,L);
  }

  ret = endpoints2chain(l, r, L, R, where);
  trace5("edge2chain new endpoints", l, r, *L, *R, ret);
  if (where) trace2("", where->first, where->second);

  std::set<DPAIR> intersections;

  for(auto a = L; a != R; ++a){
    trace2("iterating left...", *a, *next(a));
    if (next(a)==end()){
    }else {
      for(auto b = next(a); b != R; ++b){
	assert (next(b)!=end());
	trace2("iterating right...", *b, *next(b));

	if (b==end()){ incomplete();

	} else if (next(b)==end()){ incomplete();
	  // we have only one edge left

	} else {
	  // need to check bisector

	  FPOLY1 bis = bisector_dual(a,b); // cache?
	  trace5("bisector...", *a, *next(a), *b, *next(b), bis);
	  FPOLY1 edge(l,r);

	  double t = CPOLY1(edge).zero(bis);

	  if(t>=l.first && t<=r.first){
	    double ft = edge(t);
	    trace6("found intersection...", t, ft, *a, *b, l, r);
	    DPAIR ww;
	    DPAIR* w2 = (where)? &ww : NULL;
	    double d = point2chain(DPAIR(t,ft), a, w2); // FIXME. pass right end
	    if(w2) trace1("point2chain done", *w2);
	    if(d<ret){
	      trace1("intersection is closer", d);
	    }else{
	      if(where){
		trace2("further intersection...", *w2, d);
		where->first = DPAIR(t,ft);
		where->second = ww;
	      }
	      ret = d;
	    }
	  }else{
	    trace3("no intersection", t, bis, edge);
	  }
	}
      }
    }
  }
  return ret;
}
/*--------------------------------------------------------------------------*/
// compute directed hausdorff distance from this to that using linear
// interpolations
// == max_{i \in lin(this)} min_{j \in lin(that)} d(i,j)
inline double WAVE::dhd_linear(const WAVE& that, std::pair<DPAIR,DPAIR>* wp) const
{
  auto L = that.begin();
  auto R = L;
  assert(L!=that.end()); // for now.

  double ret = 0.; // fabs(that.v_out(begin()->first).f0 - begin()->second);

  for(auto ii = begin(); ; ++ii){
    auto other = next(ii);

//    ret = std::max(ret, fabs(other->second - that.v_out(other->first).f0)); // hmmm
    trace3("dhd outer loop", ret, *ii, *other);

    if (other==end()){
      // necessary?!
      trace3("dhd, end", *ii, ret, *L);
      DPAIR ww;
      auto wp2 = (wp)? &ww : NULL;
      double r = that.point2chain(*ii, L, wp2);
      if(ret<r){
	ret = r;
	if(wp){
	  wp->first = *ii;
	  wp->second = ww;
	}
      }else{
      }
      trace3("dhd end", ii->first, ret, *L);
      if(wp) trace3("dhd end", ii->second, wp->first, wp->second);
      return ret;
    } else if (R==that.end()) { untested();
      DPAIR ww;
      auto wp2 = (wp)? &ww : NULL;
      double r = that.segment2point(ii, *L, wp2);
      if(ret<r){untested();
	ret = r;
	if(wp){untested();
	  wp->first = *ii;
	  wp->second = ww;
	}
      }
      trace3("dhd singleton chain", ret, *L, *R);
    } else {
      // perform windowed sweep.
      std::pair<DPAIR,DPAIR> ww;
      auto wp2 = (wp)? &ww : NULL;
      double r = that.edge2chain(*ii, *other, L, R, wp2);

      if(r>ret){
	ret = r;
	if(wp){
	  *wp = ww;
	}
	trace5("dhd windowed sweep, found", *ii, *other, ret, *L, *R);
	if(wp) trace3("dhd windowed sweep", ret, wp->first, wp->second);
      }else{
      }
    }

  }

  return ret;
}
/*--------------------------------------------------------------------------*/
// compute directed hausdorff distance from points of this to points of that
// this is useless, so lets not try to be efficient.
inline double WAVE::dhd_discrete(const WAVE& that, std::pair<DPAIR, DPAIR>* where) const
{
  double ret = 0;
  DPAIR here;
  for(auto ii : _w) {
    double dist = inf;
    for(auto j : that._w) {
      double norm = norm2(ii,j);
      if(norm < dist){
	dist = norm;
	if(where) {
	  here = j;
	}
      }
    }
    if(ret<dist){
      ret = dist;
      if(where){
	where->first = ii;
	where->second = here;
      }
    }
  }
  return ret;
}
#endif
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
