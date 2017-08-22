// some more output stream helpers
// these are all inline templates, because ostream is not a base of OMSTREAM
//
#ifndef IO_MISC_H
#define IO_MISC_H

#include "e_node.h"
#include "u_nodemap.h"
#include "e_cardlist.h"

class NODE;
class NODE_MAP;
/*----------------------------------------------*/
template <class S>
inline S& operator<<( S& o, const NODE* m){
	if(m) o << m->long_label();
	return o;
}
template <class S>
inline S& operator<<( S& o, const NODE_BASE* m){
	if(m) o << m->long_label();
	return o;
}
/*----------------------------------------------*/
template <class S>
inline S& operator<<( S& o, const CARD_LIST &m)
{
	for (CARD_LIST::const_iterator ci=m.begin(); ci!=m.end(); ++ci) {
		o << (**ci).long_label();
		o << "(" << (*ci) << ")";
		o << " ";
	}
	return o;
}
/*----------------------------------------------*/
template<class S>
inline S& operator<<( S& o, const NODE_MAP& nm){
	o<<"NODE_MAP\n";
	for (NODE_MAP::const_iterator ni = nm.begin(); ni != nm.end(); ++ni) {
		NODE_BASE* n = (*ni).second;
		assert(n);
		string label = (*ni).first;
		o << "(" << n->user_number() << ":" << label << ":" << n->long_label() << ")";
	}
	return o;
}
/*----------------------------------------------*/
template <class S>
inline S& operator<<( S& o, const PARAM_LIST_BASE &m)
{
	o << "[";
	for (PARAM_LIST::const_iterator ci=m.begin(); ci!=m.end(); ++ci) {
		o << ci->first;
		o << "(" << ci->second << ")";
		o << " ";
	}
	o << "]";
	if(m.try_again())
		o << "->" << *(m.try_again());
	return o;
}
/*----------------------------------------------*/
template <class S>
inline S& operator<<( S& o, const DPAIR &p)
{
	o<<"("<<p.first << ", " << p.second<<")";
	return o;
}

#endif
