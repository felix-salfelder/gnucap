/*                  -*- C++ -*-
 *
 * (c) 2010-13 Felix Salfelder
 *
 * GPLv3
 *
 */

#include "e_aux.h"
#include "e_storag.h"
#include "globals.h"
#include "e_elemnt.h"
#include "d_bti.h"
#include "io_trace.h"
#include "u_parameter.h"
#include "d_mos.h"
#include "d_mos8.h"
/*--------------------------------------------------------------------------*/
template <>
inline polarity_t PARAMETER<polarity_t>::_NOT_INPUT() const { return pP; untested();} // BUG
template <>
inline void PARAMETER<polarity_t>::parse(CS&)
{ incomplete(); }
/*--------------------------------------------------------------------------*/
static bool dummy=false;
/*--------------------------------------------------------------------------*/
using std::string;
/*--------------------------------------------------------------------------*/
class COMMON_COMPONENT;
/*--------------------------------------------------------------------------*/
//int  MODEL_BUILT_IN_BTI_MATRIX::foo(){return 1;}
//int  MODEL_BUILT_IN_BTI_SINGLE::foo(){return 1;}
/*--------------------------------------------------------------------------*/
int DEV_BUILT_IN_BTI::_count = -1;
int COMMON_BUILT_IN_BTI::_count = -1;
static COMMON_BUILT_IN_BTI Default_BUILT_IN_BTI(CC_STATIC);
/*--------------------------------------------------------------------------*/

//const double NA(NOT_INPUT);
const double INF(BIGBIG);
/*--------------------------------------------------------------------------*/
int MODEL_BUILT_IN_BTI::_count = 0;
/*--------------------------------------------------------------------------*/
namespace { //
  static DEV_BUILT_IN_BTI p3d;
  static MODEL_BUILT_IN_BTI p3(&p3d);
  static DISPATCHER<MODEL_CARD>::INSTALL
    d3(&model_dispatcher, "bti_default", &p3);
}
/*--------------------------------------------------------------------------*/
void SDP_BUILT_IN_BTI::init(const COMMON_COMPONENT* cc)
{
  assert(cc);
  SDP_CARD::init(cc);
}
/*--------------------------------------------------------------------------*/
TDP_BUILT_IN_BTI::TDP_BUILT_IN_BTI(const DEV_BUILT_IN_BTI*)
{ untested();
}
/*--------------------------------------------------------------------------*/
MODEL_BUILT_IN_BTI::MODEL_BUILT_IN_BTI(const BASE_SUBCKT* p)
  :MODEL_CARD(p),
   gparallel(0.0),
   flags(int(USE_OPT)),
   mos_level(0),
   rcd_number(1),
   bti_type(0),
   bti_base(10.0),
   anneal(true),
   rcd_model_name( (std::string) "rcdmodel"),
   weight(0.199),
   symmetric(false)
{
  if (ENV::run_mode != rPRE_MAIN) {
    ++_count;
  }else{
  }
  set_default(&_tnom_c, OPT::tnom_c);
}
/*--------------------------------------------------------------------------*/
MODEL_BUILT_IN_BTI::MODEL_BUILT_IN_BTI(const MODEL_BUILT_IN_BTI& p)
  :MODEL_CARD(p),
   gparallel(p.gparallel),
   flags(p.flags),
   mos_level(p.mos_level),
   rcd_number(p.rcd_number),
   bti_type(p.bti_type),
   bti_base(p.bti_base),
   anneal(p.anneal),
   rcd_model_name(p.rcd_model_name),
   weight(p.weight),
   symmetric(p.symmetric)
{
  if (ENV::run_mode != rPRE_MAIN) {
    ++_count;
  }else{ untested();
    untested();//194
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_BTI_SINGLE::attach_rcds(COMMON_BUILT_IN_RCD** _RCD) const
{ untested();
  trace0("MODEL_BUILT_IN_BTI_SINGLE::attach_rcds()");
  int i=0;
  assert((int)rcd_number==1);
  COMMON_BUILT_IN_RCD* RCD1 = new COMMON_BUILT_IN_RCD;
  RCD1->set_modelname( rcd_model_name ); // <= !
  RCD1->Uref = uref ;
  trace0(("MODEL_BUILT_IN_BTI_SINGLE::attach_rcds set_modelname " + (std::string) rcd_model_name).c_str());
  RCD1->attach(this); // ?
  COMMON_COMPONENT::attach_common(RCD1, (COMMON_COMPONENT**)&(_RCD[i]));
  // assert (RCD1 == _RCD[i]);
}
/*--------------------------------------------------------------------------*/
 MODEL_BUILT_IN_BTI_SINGLE::MODEL_BUILT_IN_BTI_SINGLE(const BASE_SUBCKT* p)
  :MODEL_BUILT_IN_BTI(p),
  fooo(0)
{

}
/*--------------------------------------------------------------------------*/
MODEL_BUILT_IN_BTI_SINGLE::MODEL_BUILT_IN_BTI_SINGLE(const
    MODEL_BUILT_IN_BTI_SINGLE& p)
  :MODEL_BUILT_IN_BTI(p),
  fooo(0)
{ untested();
}
/*--------------------------------------------------------------------------*/
ADP_CARD* MODEL_BUILT_IN_BTI::new_adp( COMPONENT* c)const
{ untested();
  assert(c);
  return MODEL_CARD::new_adp(c);
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_BTI::dev_type()const
{ untested();
  return "btimodel ";
  if (dummy == true) { untested();
    return "btimodel";
  }else{untested();//235
    return MODEL_CARD::dev_type();
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_BTI::set_dev_type(const std::string& new_type)
{
  if (Umatch(new_type, "btimodel ")) { untested();
    dummy = true;
  }else{
    MODEL_CARD::set_dev_type(new_type);
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_BTI::param_value(int i)const
{
  switch (MODEL_BUILT_IN_BTI::param_count() - 1 - i) { untested();
  case 0:  unreachable(); return "";
  case 1:  return _tnom_c.string();
  case 2:  return gparallel.string();
  case 3:  return flags.string();
  case 4:  return mos_level.string();
  case 5:  return rcd_number.string();
  case 6:  return anneal.string();
  case 7:  return rcd_model_name.string();
  case 8:  return weight.string();
  case 11:  return symmetric.string();
//  case 8:  return "some string"; // rcd_model.string();
  default: return "unknown";
  }
}
/*--------------------------------------------------------------------------*/
bool MODEL_BUILT_IN_BTI::is_valid(const COMPONENT* d)const
{
  assert(d);
  return MODEL_CARD::is_valid(d);
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_BTI::tr_eval(COMPONENT*)const
{untested();//425
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_BTI::precalc_first()
{
  trace0(("MODEL_BUILT_IN_BTI::precalc_first" + rcd_model_name.string()).c_str());
  const CARD_LIST* par_scope = scope();
  assert(par_scope);
  MODEL_CARD::precalc_first();
  e_val(&(this->gparallel), 0.0, par_scope);
  e_val(&(this->mos_level), 0, par_scope);
  e_val(&(this->anneal), true, par_scope);
  e_val(&(this->flags), int(USE_OPT), par_scope);
  e_val(&(this->mos_level), 0, par_scope);
  rcd_number.e_val(1, par_scope);
  rcd_model_name.e_val_strange( std::string("rcd_model_hc"), par_scope);
  e_val(&(this->weight), 0.123, par_scope);
  e_val(&(this->uref), NOT_INPUT, par_scope);
  e_val(&(this->v2), false, par_scope);
  symmetric.e_val(false, par_scope);
}
/*--------------------------------------------------------------------------*/
MODEL_BUILT_IN_BTI_MATRIX::MODEL_BUILT_IN_BTI_MATRIX(const BASE_SUBCKT* p)
  :MODEL_BUILT_IN_BTI_SUM(p),
  cols(0),
  rows(0),
  base(10)
{
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_BTI_SUM::precalc_first()
{
  trace0("MODEL_BUILT_IN_BTI_SUM::precalc_first");
  const CARD_LIST* par_scope = scope();
  dpvv x; // empty params.

  rcdparm.e_val( x, par_scope);
  MODEL_BUILT_IN_BTI::precalc_first();
  if(symmetric){
    rcd_number = 2*(uint_t(rcdparm.size()));
  }else{
    rcd_number = uint_t(rcdparm.size());
  }
  e_val(&(this->rcdparm), x, par_scope);
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_BTI_SUM::set_dev_type(const std::string& new_type)
{
  if (Umatch(new_type, "matrix ")) { untested();
  }
  MODEL_BUILT_IN_BTI::set_dev_type(new_type);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_BTI_MATRIX::set_dev_type(const std::string& new_type)
{ untested();
  trace0(("MODEL_BUILT_IN_BTI_MATRIX::set_dev_type " + new_type).c_str());
  if (Umatch(new_type, "matrix ")) { untested();
  }
  MODEL_BUILT_IN_BTI_SUM::set_dev_type(new_type);
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_BTI_MATRIX::precalc_first()
{ untested();
  trace0(("MODEL_BUILT_IN_BTI_MATRIX::precalc_first" + rcd_model_name.string()).c_str());
  const CARD_LIST* par_scope = scope();
  assert(par_scope);
  rows.e_val( 2, par_scope);
  cols.e_val( 3, par_scope);
  e_val(&(this->base), 11.0, par_scope);
  MODEL_BUILT_IN_BTI_SUM::precalc_first();
  rcd_number=rows*cols;
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_BTI::precalc_last()
{
    MODEL_CARD::precalc_last();
}
/*--------------------------------------------------------------------------*/
SDP_CARD* MODEL_BUILT_IN_BTI::new_sdp(COMMON_COMPONENT* c)const
{
  assert(c);
  if (COMMON_BUILT_IN_BTI* cc = dynamic_cast<COMMON_BUILT_IN_BTI*>(c)) {
    if (cc->_sdp) {
      cc->_sdp->init(cc);
      return cc->_sdp;
    }else{
      delete cc->_sdp;
      return new SDP_BUILT_IN_BTI(c);
    }
  }else{ untested();
    return MODEL_CARD::new_sdp(c);
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_BTI::set_param_by_index(int i, std::string& value, int offset)
{
  trace2("MODEL_BUILT_IN_BTI::set_param_by_index", i, offset );
  switch (MODEL_BUILT_IN_BTI::param_count() - 1 - i) { untested();
  case 0: untested(); break;
  case 1: _tnom_c = value; break;
  case 2: gparallel = value; break;
  case 3: flags = value; break;
  case 4: mos_level = value; break;
  case 5: rcd_number = value; break;
  case 6: anneal = value; break;
  case 7: rcd_model_name = value; break;
  case 8: weight = value; break;
  case 9: uref = value; break;
  case 10: v2 = value; break;
  case 11: symmetric = value; break;
  default: throw Exception_Too_Many(unsigned(i), unsigned(MODEL_BUILT_IN_BTI::param_count()), offset);
           break;
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_BTI::RCD_name(uint_t i) const
{
  stringstream a;
  a << "RCD" << i;
  return a.str();
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_BTI::param_name(int i)const
{
  switch (MODEL_BUILT_IN_BTI::param_count() - 1 - i) { untested();
  case 0:  return "=====";
  case 1:  return "tnom";
  case 2:  return "gparallel";
  case 3:  return "flags";
  case 4:  return "mos_level";
  case 5:  return "rcd_number";
  case 6:  return "anneal";
  case 7:  return "rcd_model_name";
  case 8:  return "weight";
  case 9:  return "uref";
  case 10: return "v2";
  case 11: return "symmetric";
  default: return "";
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_BTI::param_name(int i, int j)const
{
  trace2("MODEL_BUILT_IN_BTI::param_name", i, j );
  if (j == 0) { untested();
    return param_name(i);
  }else if (j == 1) {
    switch (MODEL_BUILT_IN_BTI::param_count() - 1 - i) { untested();
    case 0:  return "=====";
    case 1:  return "";
    case 2:  return "gp";
    case 3:  return "";
    case 4:  return "";
    case 5:  return "";
    case 6:  return "";
    case 7:  return "rcd_model";
    case 8:  return "";
    case 9:  return "uref";
    case 10: return "v2";
    case 11: untested(); return "";
    default: return "";
    }
  }else{ untested();
    return "";
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_BTI_SUM::param_name(int i)const
{
  switch (MODEL_BUILT_IN_BTI_SUM::param_count() - 1 -i ){ untested();
        case 0: return "=====";
        case 1: return "params";
  }
  return MODEL_BUILT_IN_BTI::param_name(i);
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_BTI_SUM::param_value(int i)const
{
  switch (MODEL_BUILT_IN_BTI_SUM::param_count() - 1 - i) { untested();
    case 0:  unreachable(); return "";
    case 1:  return string(rcdparm);
  }
  return MODEL_BUILT_IN_BTI::param_value(i);
}
/*--------------------------------------------------------------------------*/
MODEL_BUILT_IN_BTI_SUM::MODEL_BUILT_IN_BTI_SUM(const BASE_SUBCKT* p)
  :MODEL_BUILT_IN_BTI(p),
  rcdparm( PARAMETER< dpvv>() )
{
}
/*--------------------------------------------------------------------------*/
MODEL_BUILT_IN_BTI_SUM::MODEL_BUILT_IN_BTI_SUM(const MODEL_BUILT_IN_BTI_SUM& p)
  :MODEL_BUILT_IN_BTI(p),
   rcdparm(p.rcdparm)
{
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_BTI_SUM::param_name(int i, int j)const
{
  if (j == 0) {
    return param_name(i);
  }else if (j == 1) {
    switch (MODEL_BUILT_IN_BTI_SUM::param_count() - 1 - i) { untested();
      case 0:  return "=====";
      case 1:  return "parms";
    }
    return MODEL_BUILT_IN_BTI::param_name( i, j  );
  }else{
    return "";
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_BTI_SUM::set_param_by_index(int i, std::string& value, int offset)
{
  switch (MODEL_BUILT_IN_BTI_SUM::param_count() - 1 - i) { untested();
    case 0: untested(); break;
    case 1: trace1("MODEL_BUILT_IN_BTI_SUM::set_param_by_index", value);
            rcdparm = value; 
            break;
    default: MODEL_BUILT_IN_BTI::set_param_by_index(i,value,offset);
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_BTI_SUM::attach_rcds(COMMON_BUILT_IN_RCD** _RCD) const
{
  trace0("MODEL_BUILT_IN_BTI_SUM()");
  trace0(rcd_model_name.string().c_str());

  size_t Nr = rcdparm.size();
  if(!Nr) { untested();
    untested();
  }

  trace1("MODEL_BUILT_IN_BTI_SUM::attach_rcds", Nr);

  for (size_t n=0; n<Nr; n++) {
    dpvv P = dpvv(rcdparm);
    PARAMETER<dpv> Q = P[n];
    assert(Q.size()==5 || Q.size()==6); // "positive" is optional

    COMMON_BUILT_IN_RCD* RCD1 = new COMMON_BUILT_IN_RCD;
    RCD1->set_modelname( rcd_model_name ); // <= !
    RCD1->attach(this); // ?

    size_t Np = dpv(Q).size();
    for (uint_t m=0; m<Np ; m++ ) {
      trace1("MODEL_BUILT_IN_BTI_SUM::attach_rcds ", m);
      PARAMETER<double> pk = dpv(Q)[m];
      string s = pk.string();
      uint_t index =  uint_t(RCD1->param_count() - int(m) - 2);
      trace3("MODEL_BUILT_IN_BTI_SUM::attach_rcds ", n, index, s);
      RCD1->set_param_by_index( int(index), s, 0 );
    }
    COMMON_COMPONENT::attach_common(RCD1, (COMMON_COMPONENT**)&(_RCD[n]));
//    _RCD[n]->precalc_first(scope());
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
MODEL_BUILT_IN_BTI_MATRIX::MODEL_BUILT_IN_BTI_MATRIX(const MODEL_BUILT_IN_BTI_MATRIX& p)
  :MODEL_BUILT_IN_BTI_SUM(p),
   cols(1),
   rows(1),
   base(1)
{ untested();
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_BTI_MATRIX::attach_rcds(COMMON_BUILT_IN_RCD** _RCD) const
{ untested();
  trace0("MODEL_BUILT_IN_BTI_MATRIX::attach_rcds()");
  trace0(rcd_model_name.string().c_str());
  uint_t row, col , k;
  assert ((uint_t)rcd_number==(rows*cols));
  trace3("attaching", rows, cols, rcd_number);

  long double up = 1;
  long double down = 1;
  //double lambda=1;
  //double mu=1;

  // k
  // 1 2 3
  // 4 5 6
  // 7 8 9

  for(col=0; col < cols; col++ ){ untested();
    up *= base;
    for(row=0; row < rows; row++ ){ untested();
      down *=base;
      k=rows*col+row;
      COMMON_BUILT_IN_RCD* RCD1 = new COMMON_BUILT_IN_RCD;
      RCD1->set_modelname( rcd_model_name ); // <= !
      RCD1->attach(this); // ?
//  assert ((double)uref != NOT_INPUT );
      RCD1->Uref = uref;
      RCD1->Recommon0 = double(up);
      RCD1->Rccommon0 = double(down);

      trace5("MODEL_BUILT_IN_BTI_MATRIX::attach_rcds ", row, col, k, up, down); 

      COMMON_COMPONENT::attach_common(RCD1, (COMMON_COMPONENT**)&(_RCD[k]));
      // assert (RCD1 == _RCD[i]);
    }
    down=1;

  }

}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_BTI_MATRIX::set_param_by_index(int i, std::string& value, int offset)
{ untested();
  trace2(("MATRIX::set_param_by_index" + value).c_str(), i, offset );
  switch (MODEL_BUILT_IN_BTI_MATRIX::param_count() - 1 - i) { untested();
    case 0: untested(); break;
    case 1: rows = value; break;
    case 2: cols = value; break;
    case 3: base = value; break;
    default: MODEL_BUILT_IN_BTI_SUM::set_param_by_index(i,value,offset);
  }
}
/*--------------------------------------------------------------------------*/
int MODEL_BUILT_IN_BTI_MATRIX::param_count( )const
  {return (3 + MODEL_BUILT_IN_BTI_SUM::param_count());}
/*--------------------------------------------------------------------------*/
bool MODEL_BUILT_IN_BTI_SUM::param_is_printable(int i)const
{
  switch (MODEL_BUILT_IN_BTI_SUM::param_count() - 1 - i) { untested();
  case 0: return (false);
  case 1: return (true);
  default: return MODEL_BUILT_IN_BTI::param_is_printable( i );
  }
}
/*--------------------------------------------------------------------------*/
namespace { //
  static DEV_BUILT_IN_BTI p1d;
  static MODEL_BUILT_IN_BTI_SUM p1(&p1d);
  static DISPATCHER<MODEL_CARD>::INSTALL
    d1(&model_dispatcher, "bti_sum|btisum", &p1);
}
/*--------------------------------------------------------------------------*/
bool MODEL_BUILT_IN_BTI::param_is_printable(int i)const
{
  trace1("MODEL_BUILT_IN_BTI::param_is_printable", i );
  switch (MODEL_BUILT_IN_BTI::param_count() - 1 - i) { untested();
  case 0:  return (false);
  case 1:  return (true);
  case 2:  return (gparallel != 0.);
  case 3:  return (!(flags & USE_OPT));
  case 4:  return (mos_level.has_hard_value());
  case 5:  return (true);
  case 6:  return (true);
  case 7:  return (true);
  case 8:  return (true);
  default: return false;
  }
}
/*--------------------------------------------------------------------------*/
bool MODEL_BUILT_IN_BTI_MATRIX::param_is_printable(int i)const
{ untested();
  trace1("MODEL_BUILT_IN_BTI_MATRIX::param_is_printable", i );
  switch (MODEL_BUILT_IN_BTI_MATRIX::param_count() - 1 - i) { untested();
  case 0:  return (false);
  case 1:
  case 2:
  case 3:
           return (true);
  default: return MODEL_BUILT_IN_BTI_SUM::param_is_printable( i );
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_BTI_MATRIX::param_name(int i)const
{ untested();
  trace1("MODEL_BUILT_IN_BTI_MATRIX::param_name",i);
  switch (MODEL_BUILT_IN_BTI_MATRIX::param_count() - 1 -i ){ untested();
        case 0: return "=====";
        case 1: return "rows";
        case 2: return "cols";
        case 3: return "base";
  }
  return MODEL_BUILT_IN_BTI::param_name( i  );

}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_BTI_MATRIX::param_name(int i, int j)const
{ untested();
  if (j == 0) { untested();
    return param_name(i);
  }else if (j == 1) { untested();
    switch (MODEL_BUILT_IN_BTI_MATRIX::param_count() - 1 - i) { untested();
    case 0:  return "=====";
    default: return "";
    }
  }else{ untested();
    return "";
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_BTI_MATRIX::param_value(int i)const
{ untested();
  switch (MODEL_BUILT_IN_BTI_MATRIX::param_count() - 1 - i) { untested();
    case 0:  return "???";
    case 1:  return rows.string();
    case 2:  return cols.string();
    case 3:  return base.string();
  }
  return MODEL_BUILT_IN_BTI_SUM::param_value(i);
}
/*--------------------------------------------------------------------------*/
namespace MODEL_BUILT_IN_BTI_MATRIX_DISPATCHER { //
  static DEV_BUILT_IN_BTI p1d;
  static MODEL_BUILT_IN_BTI_MATRIX p1(&p1d);
  static DISPATCHER<MODEL_CARD>::INSTALL
    d1(&model_dispatcher, "bti_matrix_old", &p1);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_BTI_SINGLE::param_value(int i)const
{ untested();
  switch (MODEL_BUILT_IN_BTI_SINGLE::param_count() - 1 - i) { untested();
  case 0:  return fooo.string();
  }
  return MODEL_BUILT_IN_BTI::param_value(i);
}
/*--------------------------------------------------------------------------*/
int MODEL_BUILT_IN_BTI_SINGLE::param_count( )const
  {return (1 + MODEL_BUILT_IN_BTI::param_count());}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_BTI_SINGLE::set_param_by_index(int i, std::string& value, int offset)
{ untested();
  switch (MODEL_BUILT_IN_BTI_SINGLE::param_count() - 1 - i) { untested();
    case 0: untested(); break;
    case 1: fooo = value; break;
    default: MODEL_BUILT_IN_BTI::set_param_by_index(i,value,offset);
  }
}
/*--------------------------------------------------------------------------*/
std::string MODEL_BUILT_IN_BTI_SINGLE::param_name(int i)const
{ untested();
  switch (MODEL_BUILT_IN_BTI_SINGLE::param_count() - 1 -i ){ untested();
        case 0: return "======";
        case 1: return "fooo";
  }
  return MODEL_BUILT_IN_BTI::param_name(i );

}
/*--------------------------------------------------------------------------*/
bool MODEL_BUILT_IN_BTI_SINGLE::param_is_printable(int i)const
{ untested();
  switch (MODEL_BUILT_IN_BTI_SINGLE::param_count() - 1 - i) { untested();
  case 0:  return (false);
  case 1:  return (true);
  default: return MODEL_BUILT_IN_BTI::param_is_printable( i );
  }
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_BTI_SINGLE::precalc_first()
{ untested();
    const CARD_LIST* par_scope = scope();
    assert(par_scope); USE(par_scope);
    MODEL_BUILT_IN_BTI::precalc_first();
}
/*--------------------------------------------------------------------------*/
namespace MODEL_BUILT_IN_BTI_SINGLE_DISPATCHER { 
  static DEV_BUILT_IN_BTI p2d;
  static MODEL_BUILT_IN_BTI_SINGLE p2(&p2d);
  static DISPATCHER<MODEL_CARD>::INSTALL
    d2(&model_dispatcher, "bti_single", &p2);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
COMMON_BUILT_IN_BTI::COMMON_BUILT_IN_BTI(int c)
  :COMMON_COMPONENT(c),
   lambda(1.0),
   number(0),
   weight(1.0),
   _sdp(0),
   _RCD(0)
{
  polarity = pN;
  trace1("COMMON_BUILT_IN_BTI::COMMON_BUILT_IN_BTI new bti", c);
  ++_count;
}
/*--------------------------------------------------------------------------*/
COMMON_BUILT_IN_BTI::COMMON_BUILT_IN_BTI(const COMMON_BUILT_IN_BTI& p)
  :COMMON_COMPONENT(p),
   lambda(p.lambda),
   number(p.number),
   weight(p.weight),
   polarity(p.polarity),
   _sdp(0),
   _RCD(0)
{
   trace0("some bti copy");      
  ++_count;
}
/*--------------------------------------------------------------------------*/
COMMON_BUILT_IN_BTI::~COMMON_BUILT_IN_BTI()
{
// assert (_RCD1 == _RCD[0]);
//  assert(_RCD1);
if (_RCD) {  detach_common(&(_RCD[0]));
  //detach_common(&_RCD2);
  //detach_common(&_RCD3);
  delete []  _RCD;
  _RCD=NULL;

}
  --_count;
  delete _sdp;
}
/*--------------------------------------------------------------------------*/
bool COMMON_BUILT_IN_BTI::operator==(const COMMON_COMPONENT& x)const
{
  const COMMON_BUILT_IN_BTI* p = dynamic_cast<const COMMON_BUILT_IN_BTI*>(&x);
  trace0("COMMON_BUILT_IN_BTI::operator==");

  return (p
    && lambda == p->lambda
    && number == p->number
    && 0
    && weight == p->weight
    && _sdp == p->_sdp
    && COMMON_COMPONENT::operator==(x));
}
/*--------------------------------------------------------------------------*/
void COMMON_BUILT_IN_BTI::set_param_by_index(int I, std::string& Value, int Offset)
{ untested();
  trace1("COMMON_BUILT_IN_BTI::set_param_by_index " + Value , I);
  switch (COMMON_BUILT_IN_BTI::param_count() - 1 - I) { untested();
    case 0:  lambda = Value; break;
    case 1:  number = Value; break;
    case 2:  weight = Value; break;
    default: COMMON_COMPONENT::set_param_by_index(I, Value, Offset);
  }
}
/*--------------------------------------------------------------------------*/
bool COMMON_BUILT_IN_BTI::param_is_printable(int i)const
{ untested();
  switch (COMMON_BUILT_IN_BTI::param_count() - 1 - i) { untested();
    case 0:  return (true);
    case 1:  return (true);
    case 2:  return (true);
    default: return COMMON_COMPONENT::param_is_printable(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string COMMON_BUILT_IN_BTI::param_name(int i)const
{ untested();
  switch (COMMON_BUILT_IN_BTI::param_count() - 1 - i) { untested();
  case 0:  return "lambda";
  case 1:  return "num";
  case 2:  return "weight";
  default: return COMMON_COMPONENT::param_name(i);
  }
}
/*--------------------------------------------------------------------------*/
std::string COMMON_BUILT_IN_BTI::param_name(int i, int j)const
{ untested();
  if (j == 0) { untested();
    return param_name(i);
  }else if (j == 1) { untested();
    switch (COMMON_BUILT_IN_BTI::param_count() - 1 - i) { untested();
    case 0:  return "";
    case 1:  return "";
    case 2:  return "";
    default: return "";
    }
  }else{untested();//281
    return COMMON_COMPONENT::param_name(i, j);
  }
}
/*--------------------------------------------------------------------------*/
std::string COMMON_BUILT_IN_BTI::param_value(int i)const
{ untested();
  switch (COMMON_BUILT_IN_BTI::param_count() - 1 - i) { untested();
  case 0:  return lambda.string();
  case 1:  return number.string();
  case 2:  return weight.string();
  default: return COMMON_COMPONENT::param_value(i);
  }
}
/*--------------------------------------------------------------------------*/
void COMMON_BUILT_IN_BTI::precalc_first(const CARD_LIST* par_scope)
{
//  const MODEL_BUILT_IN_BTI* m = prechecked_cast<const MODEL_BUILT_IN_BTI*>(model());
//  assert(m); not yet?
  assert(par_scope);
  COMMON_COMPONENT::precalc_first(par_scope);

  e_val(&(this->lambda), 1.0, par_scope);
  e_val(&(this->number), uint_t(0), par_scope);
  e_val(&(this->weight), 17.0, par_scope);

}
/*--------------------------------------------------------------------------*/
void COMMON_BUILT_IN_BTI::expand(const COMPONENT* d)
{
  trace0("COMMON_BUILT_IN_BTI::expand");
  COMMON_COMPONENT::expand(d);
  trace0(("COMMON_BUILT_IN_BTI::expand attaching " + modelname() ).c_str());
  attach_model(d);
  COMMON_BUILT_IN_BTI* c = this; USE(c);
  const MODEL_BUILT_IN_BTI* m = dynamic_cast<const MODEL_BUILT_IN_BTI*>(model());
  if (!m) { untested();
    throw Exception_Model_Type_Mismatch(d->long_label(), modelname(), "bti");
  }
  trace1("COMMON_BUILT_IN_BTI ", m->rcd_number);
  
  if (!_RCD) {
    trace1("alloc", m->rcd_number);
    assert(!_RCD);
    _RCD = new COMMON_COMPONENT*[m->rcd_number];
    uint_t i;
    for(i=0; i<(uint_t)m->rcd_number; i++ ) {
      _RCD[i]=NULL;
    }
  }else{ untested();
    trace1("allocd", m->rcd_number);
  }
  // size dependent
  //delete _sdp;
  _sdp = m->new_sdp(this);
  assert(_sdp);
  const SDP_BUILT_IN_BTI* s = prechecked_cast<const SDP_BUILT_IN_BTI*>(_sdp);
  assert(s); USE(s);

  // subcircuit commons, recursive
  //

  m->attach_rcds( (COMMON_BUILT_IN_RCD**) _RCD);
  trace0("COMMON_BUILT_IN_BTI::expand attached rcds");

  //weight=m->weight;

  assert(c == this);
}
/*--------------------------------------------------------------------------*/
void COMMON_BUILT_IN_BTI::precalc_last(const CARD_LIST* par_scope)
{
  assert(par_scope);
  COMMON_COMPONENT::precalc_last(par_scope);
  COMMON_BUILT_IN_BTI* c = this; USE(c);
  const MODEL_BUILT_IN_BTI* m = prechecked_cast<const MODEL_BUILT_IN_BTI*>(model());
  assert(m);
  e_val(&(this->lambda), 1.0, par_scope);
  e_val(&(this->number), uint_t(0), par_scope);
  e_val(&weight, m->weight, par_scope);
  _sdp = m->new_sdp(this);
  assert(_sdp);
  const SDP_BUILT_IN_BTI* s = prechecked_cast<const SDP_BUILT_IN_BTI*>(_sdp);
  assert(s); USE(s);

  //
  m->attach_rcds( (COMMON_BUILT_IN_RCD**) _RCD);

  assert(c == this);
  return;
}
/*--------------------------------------------------------------------------*/
void MODEL_BUILT_IN_BTI::attach_rcds(COMMON_BUILT_IN_RCD** _RCDc) const
{ untested();
  trace1("MODEL_BUILT_IN_BTI::attach_rcds()", (long int)_RCDc[0]);
  uint_t i;
  trace1(("attach_rcds "+ std::string(rcd_model_name)).c_str() , rcd_number);
  for(i=0; i < rcd_number; i++ ){ untested();
    COMMON_BUILT_IN_RCD* RCD1 = new COMMON_BUILT_IN_RCD;
    trace0(("MODEL_BUILT_IN_BTI::attach_rcds rcd modelname set to " 
          + std::string(rcd_model_name) ).c_str());
    RCD1->set_modelname( rcd_model_name ); // <= !
    RCD1->attach(this); // ?
    RCD1->weight = 1;
    RCD1->Recommon0 = hoch(int(i));
    RCD1->Rccommon0 = runter(int(i));
    RCD1->Rccommon1 = hoch(int(i));
    COMMON_COMPONENT::attach_common(RCD1, (COMMON_COMPONENT**)&(_RCDc[i]));
    trace1("", i );
    // assert (RCD1 == _RCD[i]);
  }

}
/*--------------------------------------------------------------------------*/
double MODEL_BUILT_IN_BTI::hoch(int i) const { untested();
  uint_t n = rcd_number;
  //int type = m->bti_type;
  int a = int(sqrt(n));
  int col=i % a;
  // int row=i / a;

  return 100 * pow(10,col);
}
/*--------------------------------------------------------------------------*/
double MODEL_BUILT_IN_BTI::runter(int i) const{ untested();
  uint_t n = rcd_number;
  //int type = m->bti_type;
  int a = int(sqrt(n));
  // int col=i % a;
  int row=i / a;

  return 100 * pow(10,row);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
namespace DEV_BUILT_IN_BTI_DISPATCHER { 
  static DEV_BUILT_IN_BTI p0;
  static DISPATCHER<CARD>::INSTALL
    dd0(&device_dispatcher, "bti", &p0);
}
/*--------------------------------------------------------------------------*/
DEV_BUILT_IN_BTI::DEV_BUILT_IN_BTI()
  :BASE_SUBCKT(),
   // input parameters,
   // calculated parameters,
   cutoff(false),
   // netlist,
   _RCD(0)
//   _RCD2(0),
//   _RCD3(0)
{
  _n = _nodes;
  trace1("DEV_BUILT_IN_BTI::DEV_BUILT_IN_BTI calling attach", hp(&Default_BUILT_IN_BTI));
  trace1("DEV_BUILT_IN_BTI::DEV_BUILT_IN_BTI" , Default_BUILT_IN_BTI.attach_count());
  trace0("DEV_BUILT_IN_BTI::DEV_BUILT_IN_BTI" + Default_BUILT_IN_BTI.modelname());
  attach_common(&Default_BUILT_IN_BTI);
  ++_count;
  // overrides
}
/*--------------------------------------------------------------------------*/
DEV_BUILT_IN_BTI::DEV_BUILT_IN_BTI(const DEV_BUILT_IN_BTI& p)
  :BASE_SUBCKT(p),
   // input parameters,
   // calculated parameters,
   cutoff(p.cutoff),
   // netlist,
   _RCD(0)
{
  _n = _nodes;
  for (uint_t ii = 0; ii < max_nodes() + int_nodes(); ++ii) {
    _n[ii] = p._n[ii];
  }
  ++_count;
  // overrides
}
/*--------------------------------------------------------------------------*/
bool DEV_BUILT_IN_BTI::do_tr()
{
  assert(subckt());set_converged(subckt()->do_tr());

  q_accept();
  return converged();
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_BTI::set_parameters(const std::string& Label, CARD *Owner,
			       COMMON_COMPONENT *Common, double Value,
			       uint_t /*ns*/, double s[],
			       uint_t node_count, const node_t Nodes[])
{
  COMPONENT::set_parameters(Label, Owner, Common, Value, 0, NULL, node_count, Nodes);
  _voltages_p = s;
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_BTI::precalc_first()
{
  COMPONENT::precalc_first();
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_BTI::expand()
{
  trace1("DEV_BUILT_IN_BTI::expand", long_label());
  BASE_SUBCKT::expand(); // calls common->expand, attaches model
  assert(_n);
  assert(common());
  const COMMON_BUILT_IN_BTI* c = static_cast<const COMMON_BUILT_IN_BTI*>(common());
  assert(c);
  assert(c->model());
  const MODEL_BUILT_IN_BTI* m = prechecked_cast<const MODEL_BUILT_IN_BTI*>(c->model());
  assert(m);
  assert(c->sdp());
  const SDP_BUILT_IN_BTI* s = prechecked_cast<const SDP_BUILT_IN_BTI*>(c->sdp());
  assert(s); USE(s);

  if (!subckt()) {
    new_subckt();
  }else{ untested();
  }

  if (!_RCD) {
    _RCD = new COMPONENT*[m->rcd_number];
    for(uint_t i=0; i<m->rcd_number; i++){
      _RCD[i]=NULL;
    }
  }else{ untested();
    trace1("allocd", m->rcd_number);
  }

  if (_sim->is_first_expand()) {
    precalc_first();
    precalc_last();

    unsigned n = m->rcd_number;
    unsigned K = 1;
    if(m->symmetric){
      n/=2;
      K=2;
    }
    unsigned offset=0;

    node_t dsnode = _n[n_s];

    const CARD* p = device_dispatcher["rcd"];
    for(unsigned k=0; k<K; ++k){
      for(unsigned i=0; i<n ; i++ ) {
        if (!_RCD[i+n*k]) {
          assert(p);
          _RCD[i+offset] = dynamic_cast<COMPONENT*>(p->clone());
          assert(_RCD[i+offset]);
          subckt()->push_front(_RCD[i+offset]);
        } else { untested();
        }
        {
          node_t nodes[] = {_n[n_g], dsnode };
          string a(m->RCD_name(i+offset));
          _RCD[i+offset]->set_owner(this);
          _RCD[i+offset]->set_parameters(a , this, c->_RCD[i], value(), 0, NULL, 2, nodes);
          trace2("BTI::expand setpar ", a, i);

          if(m->symmetric){

          }
        }
      }
      offset+=n;
      dsnode = _n[n_d];
      assert(k<2); // incomplete;
    }
  }else{ untested();
    //precalc();
  }
  //precalc();
  subckt()->params()->set_try_again(scope()->params()); // here?
  subckt()->expand();
  //subckt()->precalc();
  assert(!is_constant());
  //if ( adp() == NULL ){ untested();
  //  attach_adp( m->new_adp( (COMPONENT*) this ) );
  //}else{ untested();
  //}
  //assert(adp());
  assert(subckt()->size() == size_t (m->rcd_number));
  // necessary?
  subckt()->set_owner(this);
  trace0("DEV_BUILT_IN_BTI::expand");
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_BTI::precalc_last()
{
  COMPONENT::precalc_last();
  assert(subckt());
  // incomplete?
}
/*--------------------------------------------------------------------------*/
// resultung dvth.
double DEV_BUILT_IN_BTI::vw()const
{ untested();
  assert(_n);
  const COMMON_BUILT_IN_BTI* c = prechecked_cast<const COMMON_BUILT_IN_BTI*>(common());
  assert(c);
  const MODEL_BUILT_IN_BTI* m = prechecked_cast<const MODEL_BUILT_IN_BTI*>(c->model());
  assert(m);
  const SDP_BUILT_IN_BTI* s = prechecked_cast<const SDP_BUILT_IN_BTI*>(c->sdp());
  assert(s); USE(s);
//  const ADP_BUILT_IN_BTI* a = prechecked_cast<const ADP_BUILT_IN_BTI*>(adp());
  double buf = 0;
  uint_t i = m->rcd_number;
  for(; i-->0; ) { untested();
    buf += CARD::probe(_RCD[i],"vw");
  }
  return buf*c->weight;
}
/*--------------------------------------------------------------------------*/
double DEV_BUILT_IN_BTI::vth()const
{
  assert(owner());
  const DEV_BUILT_IN_MOS* md = prechecked_cast<const DEV_BUILT_IN_MOS*>(owner());
  assert(md);
  const COMMON_BUILT_IN_MOS* mc = prechecked_cast<const COMMON_BUILT_IN_MOS*>(md->common());
  assert(mc);
//  const MODEL_BUILT_IN_MOS8* mm = prechecked_cast<const MODEL_BUILT_IN_MOS8*>(mc->model());
//  assert(mm);
  const SDP_BUILT_IN_MOS8* ms = prechecked_cast<const SDP_BUILT_IN_MOS8*>(mc->sdp());
  if(!ms){incomplete();
  }else{
  }

  assert(ms);
  return(ms->vth0);
}
/*--------------------------------------------------------------------------*/
static double triangle_intersection_area(double vgseff, double vgdeff, double btis, double btid, double* dads=NULL, double* dadd=NULL)
{
  trace2("triangle_intersection_area", vgseff, vgdeff);
  assert(vgseff>=vgdeff);

  if(btid < 0){
    // incomplete(); so what. use gnucap-adms.
    btid = 0;
  }

  assert(btis>=0);
  double area = -1;
  if(vgseff<=0){
    area = 0.;
    if(dads){
      *dads = *dadd = 0;
    }
  }else if (vgdeff<=0){
    double x = vgseff/(vgseff-vgdeff);

    double all = x*vgseff;
    double small = 0;
    if(btis<vgseff){
      double y = (btis - vgseff) / ( vgdeff - vgseff - btid + btis);
      small = (vgseff - btis)*y;
      assert(small <= all);
    }else{
      trace5("triangle_intersection_area", vgseff, vgdeff, x, small, all);
    }
    area = (all-small)*.5;
  }else{
    // both positive...
    if((btid - vgdeff) * (btis - vgseff) > 0){
      area = (min(btis,vgseff)+min(btid,vgdeff))*.5;
    }else{
      double t = (btis-vgseff) /  ( vgdeff - vgseff - btid + btis);
      double big = (max(btis,vgseff) + min(btid,vgdeff))*.5;
      double small = fabs(vgseff - btis)*t;
      area = big-small;
    }
  }
  assert(area>=0);
  assert(2*area<=btis+btid);
  return area;
}
/*--------------------------------------------------------------------------*/
static double dvth_from_area(double vgs, double vgd, double vth, double area)
{
  assert(vth>0);
  assert(vgs>=vgd);
  assert(area>=0);

  double ret;
  if(vgs<=vth && vgd<=vth){
    ret = 0;
  }else if (vgd<=vth){
    double x = (vgs-vth)/(vgs-vgd); USE(x);
    double rad = 2*area*(vgd-vgs) + (vgs-vth)*(vgs-vth); // vgs - 2*vgs*vth + vth*vth;
    assert(rad>-1e-15);
    ret = vgs - vth - sqrt(max(0.,rad));
    trace6("dvth_from_area", area, vgs, vgd, vth, ret, rad);
    assert(2*area <= x*(vgs-vth));
    assert(ret>=-1e-15);
    assert(ret*.9999<vgs-vth);
  }else if(vgd - vth > area){
    ret = area;
  }else{
    double done = (vgd - vth);
    ret = dvth_from_area( vgs-done, vth, vth, area-done ) + done;
    assert(done<=area);
    trace6("dvth_from_area", area, vgs, vgd, vth, ret, done);
    assert(ret*.9999<vgs-vth);
  }
//  assert (ret < area); wrong!
  assert (ret >= -1e-15);
  return max(ret,0.);
}
/*--------------------------------------------------------------------------*/
double DEV_BUILT_IN_BTI::dvths()const
{
  assert(_n);
  const COMMON_BUILT_IN_BTI* c = prechecked_cast<const COMMON_BUILT_IN_BTI*>(common());
  assert(c);
  const MODEL_BUILT_IN_BTI* m = prechecked_cast<const MODEL_BUILT_IN_BTI*>(c->model());
  assert(m);
  double btis = 0;
  uint_t i = m->rcd_number;
  if(m->symmetric){
    assert(!(i&1));
    i/=2;
  }else{
    return NOT_VALID;
  }
  for(; i-->0; ){
    const DEV_BUILT_IN_RCD* r = asserted_cast< const DEV_BUILT_IN_RCD* > ( _RCD[i] );
    double cont = r->P();
    if (cont<-1e-5) {
      error(bWARNING, "strange, %s contributes too negative %e\n", _RCD[i]->long_label().c_str(), cont);
      cont = 0;
    } else if (cont<-1e-10) {
      error(bTRACE, "strange, %s contributes (slightly) negative %e\n", _RCD[i]->long_label().c_str(), cont);
    }
    btis += cont;
  }
  btis *= c->weight;
  return btis;
}
/*--------------------------------------------------------------------------*/
double DEV_BUILT_IN_BTI::dvthd()const
{
  assert(_n);
  const COMMON_BUILT_IN_BTI* c = prechecked_cast<const COMMON_BUILT_IN_BTI*>(common());
  assert(c);
  const MODEL_BUILT_IN_BTI* m = prechecked_cast<const MODEL_BUILT_IN_BTI*>(c->model());
  assert(m);
  double bti = 0;
  uint_t i = m->rcd_number;
  if(m->symmetric){
    assert(!(i&1));
  }else{
    return NOT_VALID;
  }
  for(i=m->rcd_number/2; i<m->rcd_number; ++i ){
    const DEV_BUILT_IN_RCD* r = asserted_cast< const DEV_BUILT_IN_RCD* > ( _RCD[i] );
    double cont = r->P();
    if (cont<-1e-5) {
      error(bWARNING, "strange, %s contributes too negative %e\n", _RCD[i]->long_label().c_str(), cont);
      cont = 0;
    } else if (cont<-1e-10) {
      error(bTRACE, "strange, %s contributes (slightly) negative %e\n", _RCD[i]->long_label().c_str(), cont);
    }
    bti += cont;
  }
  bti *= c->weight;
  return bti;
}
/*--------------------------------------------------------------------------*/
double DEV_BUILT_IN_BTI::dvth()const
{
  assert(_n);
  const COMMON_BUILT_IN_BTI* c = prechecked_cast<const COMMON_BUILT_IN_BTI*>(common());
  assert(c);
  const MODEL_BUILT_IN_BTI* m = prechecked_cast<const MODEL_BUILT_IN_BTI*>(c->model());
  assert(m);
  double btis = 0;
  uint_t i = m->rcd_number;
  if(m->symmetric){
    assert(!(i&1));
    i/=2;
  }
  for(; i-->0; ){
    const DEV_BUILT_IN_RCD* R = asserted_cast< const DEV_BUILT_IN_RCD* > ( _RCD[i] );
    double cont = R->P();
    if (cont<-1e-5) {
      error(bWARNING, "Strange, %s contributes too negative %E\n", _RCD[i]->long_label().c_str(), cont);
      cont = 0;
    } else if (cont<-1e-10) {
      error(bTRACE, "Strange, %s contributes (slightly) negative %E\n", _RCD[i]->long_label().c_str(), cont);
    }
    btis += cont;
  }
  btis *= c->weight;

  double btid = 0.; // drain bti
  double area = 0.;
  if(m->symmetric){

    for(i=m->rcd_number/2; i<m->rcd_number; i++ ){
      const DEV_BUILT_IN_RCD* R = asserted_cast< const DEV_BUILT_IN_RCD* > ( _RCD[i] );
      double cont = R->P();
      if (cont<-1e-5) {
        error(bWARNING, "Strange, %s contributes too negative %E\n", _RCD[i]->long_label().c_str(), cont);
        cont = 0;
      } else if (cont<-1e-10) {
        error(bTRACE, "Strange, %s contributes (slightly) negative %E\n", _RCD[i]->long_label().c_str(), cont);
      }
      btid += cont;
    }
    
    btid *= c->weight;
    double Vgs = vgs(); // this is real vgs, not polarity dependent
    double Vgd = vgd();
    double Vth = vth();
    if(Vth<0){
      // p-type fet
      Vgs = -Vgs;
      Vgd = -Vgd;
      Vth = -Vth;
    }else{untested();
      // n-type fet
    }

    double ret;
    if(Vgs>=Vgd){ // *this* is polarity
      area = triangle_intersection_area(Vgs-Vth, Vgd-Vth, btis, btid);
      ret = dvth_from_area(Vgs, Vgd, Vth, area);
    }else{ untested();
      area = triangle_intersection_area(Vgd-Vth, Vgs-Vth, btid, btis);
      ret = dvth_from_area(Vgd, Vgs, Vth, area);
    }
    _area = area;
    _max = (btis + btid) *.5;
    double err = .5*(btis+btid) - ret;
    trace8("symm", Vgd, Vgs, Vth, area, btid, btis, ret, err);
    assert(ret<=max(btis,btid)+1e-13);
    if (ret>max(btis,btid)){
      ret = max(btis,btid);
    }
    return ret;
  }
  
  return btis;
}
/*--------------------------------------------------------------------------*/
double DEV_BUILT_IN_BTI::tr_probe_num(const std::string& x)const
{
  assert(_n);
  const COMMON_BUILT_IN_BTI* c = prechecked_cast<const COMMON_BUILT_IN_BTI*>(common());
  assert(c);
  const MODEL_BUILT_IN_BTI* m = prechecked_cast<const MODEL_BUILT_IN_BTI*>(c->model());
  assert(m);
//  const SDP_BUILT_IN_BTI* s = prechecked_cast<const SDP_BUILT_IN_BTI*>(c->sdp());
//  assert(s);
//  const ADP_BUILT_IN_BTI* a = prechecked_cast<const ADP_BUILT_IN_BTI*>(adp());

  if (Umatch(x, "vc{v} |fill ")) { 
    double buf = 0;
    unsigned i = m->rcd_number;
    for(; i-->0; ){ untested();
      buf += CARD::probe(_RCD[i],"vc");
    }
    return buf;
  }else if (Umatch(x, "rcdnumber ")) {
    return m->rcd_number;
  }else if (Umatch(x, "dvth |vth ")) {
    double buf = 0;
    unsigned i = m->rcd_number;
    while ( i-->0 ){
      buf += CARD::probe(_RCD[i],"P");
    }
    assert(is_number(buf));
    return buf * c->weight;
  }else if (Umatch(x, "vw ")) { untested();
    return vw();
  }else if (Umatch(x, "wt ")) { untested();
    return  c->weight;
  }else if (Umatch(x, "uin|vin ")) { untested();
    return  _n[n_g].v0() - (_n[n_d].v0()+_n[n_s].v0())/2.0;
  }else if (Umatch(x, "s{tress} ")) {
    return  (_n[n_g].v0() - _n[n_s].v0())*value();
    return  (_n[n_g].v0() - (_n[n_d].v0()+_n[n_s].v0())/2.0)*value();
  }else if (Umatch(x, "iscale ")) {
    return value();
  }else if (Umatch(x, "pol{arity} ")) { untested();
    return c->polarity;
  }else if (Umatch(x, "dvthd ")) {
    return dvthd();
  }else if (Umatch(x, "dvths ")) {
    return dvths();
  }else if (Umatch(x, "max ")) {
    return _max;
  }else if (Umatch(x, "area ")) { itested();
    return _area;
  }else {
    return BASE_SUBCKT::tr_probe_num(x);
  }
}
/*--------------------------------------------------------------------------*/
double DEV_BUILT_IN_BTI::tt_probe_num(const std::string& x)const
{
  assert(_n);
  const COMMON_BUILT_IN_BTI* c = prechecked_cast<const COMMON_BUILT_IN_BTI*>(common());
  assert(c);
  const MODEL_BUILT_IN_BTI* m = prechecked_cast<const MODEL_BUILT_IN_BTI*>(c->model());
  assert(m);
  const SDP_BUILT_IN_BTI* s = prechecked_cast<const SDP_BUILT_IN_BTI*>(c->sdp());
  assert(s); USE(s);
//  const ADP_BUILT_IN_BTI* a = prechecked_cast<const ADP_BUILT_IN_BTI*>(adp());
//  if(!a)untested0(("no a"+long_label()).c_str());

  if (Umatch(x, "vc{v} |fill ")) { incomplete();
    double buf = 0;
    unsigned i = m->rcd_number;
    while ( i-->0 ) { untested();
      buf += CARD::probe(_RCD[i],"vc");
    }
    return buf * c->weight;
  }else if (Umatch(x, "dvth |vth ")) {
    return dvth(); // cache?
  }else if (Umatch(x, "pol{arity} ")) { untested();
    return  c->polarity;
  }else if (Umatch(x, "num ")) { untested();
    return  m->rcd_number;
  }else {
    return BASE_SUBCKT::tt_probe_num(x);
  }
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_BTI::tr_stress_last()
{
  trace0("DEV_BUILT_IN_BTI::tr_stress_last");
  assert(subckt()); 
  subckt()->do_forall( &CARD::tr_stress_last );
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_BTI::tr_stress()
{ untested();
  assert(_n[n_s].v0() >= _n[n_d].v0());
  
  trace0("DEV_BUILT_IN_BTI::tr_stress()");
  unreachable();
//   subckt()->do_forall( &CARD::tr_stress );
}
/*--------------------------------------------------------------------------*/
void DEV_BUILT_IN_BTI::do_tt()
{
  //FIXME, subckt default
  //        RCD reicht!
  subckt()->do_tt();
}
/*--------------------------------------------------------------------------*/
// cc_direct

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
void ADP_BUILT_IN_BTI::init(const COMPONENT* ) { }
void ADP_BUILT_IN_BTI::tt_accept() {}
double ADP_BUILT_IN_BTI::tt_probe_num(const std::string& )const {untested(); return 888;}
double ADP_BUILT_IN_BTI::tr_probe_num(const std::string& )const {untested(); return 888;}

// vim:ts=8:sw=2:et:
