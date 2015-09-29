/* Copyright (C) 2011-2013 Felix Salfelder, Lars Hedrich
 * Author: Felix Salfelder
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
 * a remote control socket. used by verification tools
 */
#include "globals.h"
#include "u_status.h"
#include <unistd.h>
#include "u_prblst.h"
#include "u_cardst.h"
#include "u_nodemap.h"
#include "e_elemnt.h"
#include "e_storag.h"
#include "e_subckt.h"
#include "d_subckt.h"
#include "u_sock.h"
#include "s_ddc.h"
#include "io_error.h"
#include "resolv.h"
#include "s__.h"
#include "io_matrix.h"
#include "m_matrix_extra.h"
#include <iomanip>
#include "d_cap.h"
#include <fcntl.h>
/*--------------------------------------------------------------------------*/
#define userinfo( a,b,c,d,e,f )
/*--------------------------------------------------------------------------*/
using namespace std;
using namespace SOME_CAP_HACK; // FIXME. (maybe use STORAGE device interface?)
/*--------------------------------------------------------------------------*/
namespace {
/*--------------------------------------------------------------------------*/
typedef union {
  double double_val;
  int32_t   int_val;
} di_union_t;

/*--------------------------------------------------------------------------*/
class SOCK : public DDC_BASE {
public:
  explicit SOCK();
  ~SOCK();
  string status()const;
protected:
  void	options(CS&, int x=0);
private:
  void	sweep();
  void	sweep_recursive(int);
  void	first(int);
  bool	next(int);
  void send_matrix();
  void undo_time_step();
  void tr_reject();
  static unsigned tr_steps_rejected_;
protected:
  enum {DCNEST = 4};
  int _n_sweeps;
//  typedef void (*p)(double);
  //ELEMENT* (_pushel[DCNEST]);	/* pointer to thing to sweep, dc command */
  CARDSTASH _stash[DCNEST];	/* store std values of elements being swept */
  bool _loop[DCNEST];		/* flag: do it again backwards */
  bool _reverse_in[DCNEST];	/* flag: sweep backwards, input */
  bool _reverse[DCNEST];	/* flag: sweep backwards, working */
  bool _cont;			/* flag: continue from previous run */
  TRACE _trace;			/* enum: show extended diagnostics */
  enum {ONE_PT, LIN_STEP, LIN_PTS, TIMES, OCTAVE, DECADE} _stepmode[DCNEST];
  static double temp_c_in;	/* ambient temperature, input and sweep variable */
  bool _dump_matrix;
  unsigned port;
  double* U;
  double* CU;
  double* CUTCU;
public:
  void	do_it(CS&, CARD_LIST*);
private:
  void	setup(CS&);
  void fillnames( const CARD_LIST* scope);
  void findcaps( CARD_LIST* scope);
  void cap_prepare();
  void cap_reset();
  vector<string> var_namen_arr;
  vector<COMPONENT*> _caplist; // FIXME: use cardlist
  CARDSTASH* _capstash;
  uint16_t var_namen_total_size;

private: //vera stuff.
  void main_loop();
  void verainit(unsigned, unsigned, unsigned);
  void verakons();
  void veraop();
  void verainit_reply();
  void verakons_reply();
  void veraop_reply();
  void ev_iter();
  void set_param();
  unsigned transtep(unsigned init, double dt);
  void transtep_reply(unsigned, bool eol=true);
  void transtep_gc_reply(unsigned);

  char* var_names_buf;

  unsigned _verbose;
  size_t total;
  size_t n_inputs() const{return _input_names.size();}
  vector<string> _input_names;
  vector<ELEMENT*> _input_devs;
  unsigned n_vars;
  unsigned n_vars_square;
  uint16_t n_eingaenge;
  uint16_t length;

  SocketStream stream;
  //unsigned BUFSIZE;
  unsigned n_bytes;

  uint16_t _error; // transport between method and method_reply
  double _dthack; // same

  double *dc_werteA,*dc_loesungA,*kons_loesungA,*kons_residuumA;

  int channel;
  int frame_number;
  Socket* socket;
  short unsigned _port_range;
  string _port;   // kommt ueber die Uebergabeparamter port
                  // globale Variable daher (default: port=1400)
  bool _client_mode;
  bool _server_mode;
  unsigned _bufsize;
  bool _bigarg;
  string _host;
  int reuseaddr;
  struct sockaddr_in sin;

  double *matrixg, *matrixc,*vectorq;

  static const int printlevel=0;

};
/*--------------------------------------------------------------------------*/
unsigned SOCK::tr_steps_rejected_ = 0;
double	SOCK::temp_c_in = 0.;
/*--------------------------------------------------------------------------*/
void SOCK::do_it(CS& Cmd, CARD_LIST* Scope)
{
  trace0("SOCK::do_it");
  IO::error.detach(stdout);
  IO::error.attach(stderr);
  _scope = Scope;
  _sim->_time0 = 0.;
  //_sim->set_command_ddc();
  _sim->set_command_dc();
  _sim->_phase = p_INIT_DC;
  //_sim->_ddc = true;
  ::status.ddc.reset().start();
  _sim->_temp_c = temp_c_in;
  _dump_matrix = 0;
  reuseaddr = 0;
  _port = "1400";
  _port_range = 1;
  _client_mode = false;
  _host = "localhost";

  command_base(Cmd);

  //cleanup
  cap_reset();
}
/*--------------------------------------------------------------------------*/
SOCK::SOCK() :DDC_BASE(),
   _n_sweeps(1),
   _cont(false),
   _trace(tNONE)
{

  for (int ii = 0; ii < DCNEST; ++ii) {
    _loop[ii] = false;
    _reverse_in[ii] = false;
    _reverse[ii] = false;
    _step[ii] = 0.;
    _linswp[ii] = true;
    _zap[ii]=NULL;
    _stepmode[ii] = ONE_PT;
  }
  temp_c_in=OPT::temp_c;
  _out=IO::mstdout;
}
/*--------------------------------------------------------------------------*/
void SOCK::setup(CS& Cmd)
{
  _sim->_uic = false;
  trace0("SOCK::setup");
  _cont = false;
  _trace = tNONE;
  _out = IO::mstdout;
  _out.reset(); //BUG// don't know why this is needed */
  bool ploton = IO::plotset  &&  plotlist().size() > 0;

  options(Cmd);

  if (Cmd.more()) {
  }else{
  }
  Cmd.check(bWARNING, "what's this?");

  _sim->_uic = true;
  CKT_BASE::_sim->init();

  IO::plotout = (ploton) ? IO::mstdout : OMSTREAM();
  initio(_out);

  _error = 0; /* verainit(v_flag, n_inputs, &n_vars, charbuf, &length); */
  n_vars = static_cast<uint16_t>( _sim->_total_nodes) ; // _sim->total_nodes doesnt include gnd
  // do this also in sweep to initialize the socket streaM
  var_namen_arr.resize( n_vars, string("unset"));
//  var_namen_arr[0]="0";
  fillnames( &CARD_LIST::card_list );
  n_vars_square = n_vars * n_vars;

  findcaps(&CARD_LIST::card_list);

  // BUG! this kills cap bm's
  cap_prepare(); // attach value common if !has_common. stash old common

#ifndef NDEBUG
    for (unsigned i=0; i < n_vars; i++)
      trace0("name: " + var_namen_arr[i]);
#endif

  assert(_n_sweeps > 0);
  _sim->_freq = 0;
}
/*--------------------------------------------------------------------------*/
void SOCK::options(CS& Cmd, int Nest)
{
  trace0("SOCK::options... ");

  _sim->_uic = _loop[Nest] = _reverse_in[Nest] = false;
  _port = "1400";
  _bufsize = BUFSIZE;
  _bigarg = true;
  unsigned here = Cmd.cursor();
  do{
    ONE_OF
      || (Cmd.match1("'\"({")	&& ((Cmd >> _step_in[Nest]), (_stepmode[Nest] = LIN_STEP)))
      || (Cmd.is_float()	&& ((Cmd >> _step_in[Nest]), (_stepmode[Nest] = LIN_STEP)))
      || (Get(Cmd, "*",		  &_step_in[Nest]) && (_stepmode[Nest] = TIMES))
      || (Get(Cmd, "+",		  &_step_in[Nest]) && (_stepmode[Nest] = LIN_STEP))
      || Get(Cmd, "b{ufsize}",    &_bufsize)
      || Get(Cmd, "c{ontinue}",   &_cont)
      || Get(Cmd, "port" ,        &_port)
      || Get(Cmd, "listen{port}", &_port)
      || Get(Cmd, "bigarg",       &_bigarg)
      || Get(Cmd, "host" ,        &_host)
      || Get(Cmd, "tr{s}",        &_do_tran_step)
      || Get(Cmd, "dm",           &_dump_matrix)
      || Get(Cmd, "client",       &_client_mode)
      || Get(Cmd, "server",       &_server_mode)
      || Get(Cmd, "dt{emp}",	  &temp_c_in,   mOFFSET, OPT::temp_c)
      || Get(Cmd, "lo{op}", 	  &_loop[Nest])
      || Get(Cmd, "re{verse}",	  &_reverse_in[Nest])
      || Get(Cmd, "te{mperature}",&temp_c_in)
      || (Cmd.umatch("tr{ace} {=}") &&
	  (ONE_OF
	   || Set(Cmd, "n{one}",      &_trace, tNONE)
	   || Set(Cmd, "o{ff}",       &_trace, tNONE)
	   || Set(Cmd, "w{arnings}",  &_trace, tUNDER)
	   || Set(Cmd, "i{terations}",&_trace, tITERATION)
	   || Set(Cmd, "v{erbose}",   &_trace, tVERBOSE)
	   || Cmd.warn(bWARNING,
		       "need none, off, warnings, iterations, verbose")
	   )
	  )
      || _out.outset(Cmd);
  }while (Cmd.more() && !Cmd.stuck(&here));
}
/*--------------------------------------------------------------------------*/
static void register_status();
/*--------------------------------------------------------------------------*/
void SOCK::sweep()
{
  register_status();
  // later... FIXME
  //
  frame_number = 0;
  n_bytes = 0;
  trace1("SOCK::do_it", _port);

  // Calculate buffersize
  n_vars = _sim->_total_nodes;
  // _sim->total_nodes doesnt include gnd
  n_vars_square = n_vars * n_vars;
  _bufsize = n_vars_square * 3 * (unsigned)sizeof(double) + 32;
  // 3 quadratic matrices * doubles + some safety margin (BUG)

  if (_client_mode){
    trace1("bufsize Client ", _bufsize);
    socket = new ClientSocket(Socket::TCP, _port, _host, _bufsize);
    trace1("connected to "+ _host, _port );
    stream = *socket;
  } else if (_server_mode) {
    stringstream p(_port);
    uint16_t _port_;
    p >> _port_;
    trace1("bufsize Server ", _bufsize);
    socket = new ServerSocket(Socket::TCP, _port_, 1u, _bufsize);
    ServerSocket* sock=prechecked_cast<ServerSocket*>(socket);

    trace0("SOCK::do_it waiting");
    stream = sock->listen();
  } else {
    fflush( stdout );
    fflush( stdin );
    trace0("SOCK::sweep simple i/o");
    socket=0;
    trace1("bufsize Stdin ", _bufsize);
    stream = SocketStream( STDOUT_FILENO, STDIN_FILENO, _bufsize);
    stream << "gnucap sock ready";
  }

  main_loop();
  delete socket;
  return;
}
/*--------------------------------------------------------------------------*/
string SOCK::status()const
{
  // incomplete();
  return "sock: tr rejected=" + to_string(tr_steps_rejected_) + "\n";
}
/*--------------------------------------------------------------------------*/
void SOCK::sweep_recursive(int )
{
  assert(false);
}
/*--------------------------------------------------------------------------*/
// fetch names from circuit recursively. fill into local vector.
void SOCK::fillnames( const CARD_LIST* scope){
  trace0("SOCK::fillnames");

  const NODE_MAP * nm = scope->nodes();
  for (NODE_MAP::const_iterator i = nm->begin(); i != nm->end(); ++i) {
    if (i->first != "0") {
      if (const NODE* a= dynamic_cast<const NODE*>(i->second)){
        stringstream s;
        string myname(a->long_label());
        var_namen_arr[a->matrix_number()-1] = myname;
        var_namen_total_size = static_cast<uint16_t>( var_namen_total_size + static_cast<uint16_t>(myname.length()) + 1 );
      }
    }else{
    }
  }

  for (CARD_LIST::const_iterator i = scope->begin(); i != scope->end(); ++i) {
    if(const COMPONENT* s = dynamic_cast<const COMPONENT*>(*i)){
      if (!s->is_device()){ untested();
      }else if ( s->subckt() ) {
        fillnames( s->subckt() );
      }
    }else{ untested();
    }
  }
}
/*--------------------------------------------------------------------------*/
void SOCK::findcaps( CARD_LIST* scope){
  for (CARD_LIST::iterator i = scope->begin(); i != scope->end(); ++i) {
    if ( COMPONENT* c = dynamic_cast< COMPONENT*>(*i) ) { untested();
      if (c->is_device() && c->has_memory()){
        trace1("found cap", c->long_label());
        _caplist.push_back( c );
      }
    }
    if (!(*i)->is_device()){ untested();
    } else if ( BASE_SUBCKT* s = dynamic_cast< BASE_SUBCKT*>(*i) ) { untested();
      trace1("going down", s->long_label());
      findcaps( s->subckt() );
    }
  }
}
/*--------------------------------------------------------------------------*/
static SOCK p2;
static DISPATCHER<CMD>::INSTALL d2(&command_dispatcher, "sock", &p2);
/*--------------------------------------------------------------------------*/
static void register_status()
{
  static bool done;

  if(!done) {
    new DISPATCHER<CKT_BASE>::INSTALL(&status_dispatcher, "sock", &p2);
    done = true;
  }else{ untested();
  }
}
/*--------------------------------------------------------------------------*/
}
/*--------------------------------------------------------------------------*/
static unsigned argc(unsigned opcode)
{
  switch(opcode){
    case 51: untested();
      return 3;
    case 52: untested();
      return 0;
    case 53: untested();
      return 0;
    case 104: untested();
    case 102:
      return 1;
    default:
      return 0;
  }
}
/*--------------------------------------------------------------------------*/
void SOCK::main_loop()
{
  trace0("SOCK::main_loop");

  uint16_t opcode=-1;
  uint16_t arg[3];
  unsigned char tmp;
  arg[0] = -1;
  bool init_done=false;

  while(true) {
    trace0("SOCK::main_loop waiting for opcode");

    if(_bigarg){
      stream >> tmp >> 7; // sic!
      opcode = tmp;
      for(unsigned i=0; i<argc(opcode); ++i){
        stream >> arg[i] >> 6;
      }
    }else{
      stream >> opcode;
      stream >> arg[0];
      stream >> arg[1];
      stream >> arg[2];
    }

    _sim->_mode = s_SOCK; // nonsense.
                          // use respective built-in mode!
                          //
    if(_bigarg) {
      ::error(bDEBUG, "sock opcode %d\n", opcode);
    }else{
      ::error(bDEBUG, "sock opcode %d %d %d %d\n", opcode, arg[0], arg[1], arg[2]);
    }

    double dt;
    unsigned status;
    switch (opcode) {
      case '\0': // 0
        return;
        break;
      case '3': // 51
        if(init_done) throw Exception("init twice??");
        verainit(arg[0], arg[1], arg[2]);
        verainit_reply();
        init_done=true;
        break;
      case '4':  // 52
        veraop();
        veraop_reply();
        break;
      case '5':  // 53
        verakons();
        verakons_reply();
        break;
      case 101: untested();
        ev_iter();
        break;
      case 102:
        stream >> dt;
        status = transtep(arg[0], dt);
        transtep_reply(status);
        break;
      case 103: untested();
        set_param();
        break;
      case 104: itested();
        stream >> dt;
        status = transtep(arg[0], dt);
        transtep_gc_reply(status);
	break;
      default:
        ::error(bDANGER, "unknown opcode %i\n", opcode);
        throw Exception("unknown opcode");
        break;
    }
    _sim->reset_iteration_counter(iPRINTSTEP); // used by solve_with_homotopy
                                               // reset in outdata (not called)
  }
}
/*--------------------------------------------------------------------------*/
// very clever way to transfer strings.
static void putstring8(SocketStream* s, const string x)
{
  const char* A = x.c_str();

  while(*A){
    *s<<*A<<*A<<*A<<*A<<*A<<*A<<*A<<*A;
    A++;
  }
}
/*--------------------------------------------------------------------------*/
void SOCK::verainit(unsigned verbose, unsigned n_in, unsigned length)
{ untested();
  _verbose = verbose;
  char input_namen[length+1];
  unsigned here =0;
  unsigned n=0;
  _input_names.resize(n_in);
  _input_devs.resize(n_in);
  trace3("verainit: ", verbose,n_inputs(),length);
  assert(n_inputs()==n_in);

  for (unsigned i=0; i < length; i++) {
    stream >> input_namen[i] >> 7;
    if(input_namen[i] == '\t'){
      input_namen[i] = 0;
      trace1("input_namen", string(input_namen+here));
      trace5("input_namen",input_namen[i-2],input_namen[i-1], input_namen[i-0],here,i );
      _input_names[n] = string(input_namen+here);
      here = i+1;

      trace1("looking out for", _input_names[n]);
      CARD_LIST::fat_iterator ci = findbranch(_input_names[n], &CARD_LIST::card_list);
      if (ci.is_end()){
        throw Exception("cannot find voltage source \"" + _input_names[n] +"\", giving up");
      }
      trace0("found input device");
      ELEMENT* d = prechecked_cast<ELEMENT*>(*ci);
      if (!d){
        throw Exception("not something we can use as source, \"" +  _input_names[n] +"\", giving up");
      }
      _input_devs[n] = d;
      assert(_input_devs[n]);

      unsigned ii=n;
      // see s_cd. FIXME: revert.
      _stash[ii] = _input_devs[ii];
      _input_devs[ii]->inc_probes();
      _input_devs[ii]->set_value(_input_devs[ii]->value(),0);
      _input_devs[ii]->set_constant(false);

      ++n;
    }
  }

  //trace0("input_namen " + string(input_namen) );
  total = (unsigned) (length+4);
  assert(stream.bufsize() >= total);

  if (!stream.at_end())
  {
    printf("Error in Verainit! no of bytes received %i <> expected %i\n",
        n_bytes, (int)(total*sizeof(di_union_t)));
    throw Exception("bloed\n");
  }

  assert(!var_names_buf);
  var_names_buf = (char*) malloc( n_vars * 128 * sizeof(char));
  strcpy(var_names_buf,"");
  for (unsigned i=0; i < n_vars; i++)
  {
    trace1("SOCK::verainit ", var_namen_arr[i]);
    strcat(var_names_buf, var_namen_arr[i].c_str());
    strcat(var_names_buf, "\t");
  }
  //length = static_cast<uint16_t>( strlen(var_names_buf) );
  // userinfo(1,"vera_titan_ak","Variablennamen %s\n",var_names_buf);
}
/*--------------------------------------------------------------------------*/
void SOCK::veraop()
{
  trace1("SOCK::veraop",n_vars);
  total = n_vars;
  _sim->set_command_dc();
  assert(stream.bufsize() >= total);

  dc_werteA = (double*) malloc(sizeof(double)*n_vars);
  trace1("fetching ",n_vars);
  assert(_sim->vdc()[0] == 0 );
  for (unsigned i=0; i < n_inputs(); i++) {
    double d;
    stream >> d;
    trace3("SOCK::veraop setting input", _input_devs[i]->long_label(), i, d);
    _input_devs[i]->set_value(d);
  }

  _error = 0; /* veraop(sweep_val, x_new, G, C); */

  // ================do_dc========
  _sim->_uic = false;
  _sim->_bypass_ok = false;
  _sim->set_inc_mode_bad();
  OPT::ITL itl = OPT::DCBIAS;
  _sim->_phase = p_INIT_DC;

  trace2("SOCK::veraop homotopy", _sim->_phase, _sim->_mode);
  if(_dump_matrix) {
    _trace=tVERBOSE;
  }
  CARD_LIST::card_list.tr_begin();  // hier muesste eigentlich eine dc hin.
  if(printlist().size()) {
    head(0,0," ");
  }
  bool converged = false;
  try{
    converged = solve_with_homotopy(itl,_trace);
  }catch( Exception e) {
    ::error(bDANGER, "hot failed\n");
    throw e;
  }

  if(!converged){ itested();
    _error = 1;
  }else{

    ::status.accept.start();
    _sim->set_limit();
    trace0("SOCK::veraop, accepting");
    CARD_LIST::card_list.tr_accept();
    ::status.accept.stop();

    _sim->keep_voltages();
  }
  //========================

  // Die Variablenwerte stehen schon durch solve_system in var_werte
  // in dc_sysA und muessen noch nach A kopiert werden:
  //
  //A->var_werte = x_neu;
  //free(dc_loesungA);
  //
  //        A->var_werte = dc_loesungA;
  // diff_werte, par_werte muessen auch gesetzt sein!
  // in diesem Falle sind sie per default = 0
  //
  //  solve  somthing
  //
  //  A->eval_lin_gl(dc_werteA,&matrixg,&matrixc,&vectorq);

  trace1("matrix", 1);
}
/*--------------------------------------------------------------------------*/
void SOCK::verakons()
{
  //        n_eingaenge == #caps?
  _sim->set_command_dc();
  trace1("SOCK::kons", _caplist.size());
  _sim->_time0 = 0; // just assert?
  _sim->_dt0 = 0; // just assert?
  advance_time();
  _sim->_phase = p_INIT_DC;
  total =  n_eingaenge + n_vars + 1;

  for (unsigned i=0; i < n_inputs(); i++) {
    double d;
    stream >> d;
    trace3("SOCK::verakons setting input", _input_devs[i]->long_label(), i, d);
    _input_devs[i]->set_value(d);
  }

  assert(_sim->_v0);
  for (unsigned i=1; i <= n_vars; i++) {
    stream >> _sim->_v0[i];
    //    trace2("SOCK::kons start ", i,  _sim->_v0[i] );
  }
  if (printlist().size()) { untested();
    outdata(0.);
  }
  _sim->keep_voltages(); // v0->vdc

  _error = 0; /* verakons(Dwork, x_new, q_dot, G, C); */
  //	n_vars = A->n_var;

  for( unsigned i = 0; i < _caplist.size(); i++) {
    trace1("SOCK::kons",_caplist[i]->long_label());
    _caplist[i]->keep_ic(); // latch voltage applied to _v0
    _caplist[i]->set_constant(true);
    _caplist[i]->q_eval();		// so it will be updated
  }
  //
  //================do_dc========
  _sim->_uic = true;
  _sim->_bypass_ok = false;
  _sim->set_inc_mode_bad();

  OPT::ITL itl = OPT::DCBIAS;

  if(printlist().size()) { untested();
    head(0,0," ");
  }
  CARD_LIST::card_list.tr_begin();
  bool converged = false;
  try{
    converged = solve(itl,_trace);
    //    solve(OPT::TRHIGH,_trace);
    //    solve_with_homotopy(itl,_trace);
    // homotopy is to much effort calling
    // procedures must catch problems with convergence.
  }catch( Exception e) {
    ::error(bDANGER, "hot failed\n");
    throw e;
  }
  if (!converged) {
    ::error(bWARNING, "s_sock::verakons: solve did not converge\n");
    converged = solve_with_homotopy(itl,_trace);
    if (!converged) {
      ::error(bWARNING, "s_sock::verakons: solve did not converge even with homotopy\n");
      _error = 1;
    }
  }
  ::status.accept.start();
  assert(_sim->_uic);
  CARD_LIST::card_list.tr_accept();
  ::status.accept.stop();

  //  assert(_sim->_mode==s_SOCK);
  if (printlist().size()) { untested();
    outdata(_sim->_time0);
  }
  _sim->_mode = s_SOCK; // for now.

  _sim->keep_voltages();
  _sim->zero_currents();
  _sim->set_inc_mode_no();

  // vera wants just cap stamps
  for (unsigned i = 0; i < _caplist.size(); i++) { itested();
    CARD* c = prechecked_cast<CARD*>(_caplist[i]);
    assert(c);
    trace1("verakons, loading cap", c->long_label());
    _sim->_damp = 1.; // need raw stamps
    c->tr_load();
  }
}
/*--------------------------------------------------------------------------*/
#undef rescale
#define normalize
void SOCK::ev_iter()
{
  trace1("SOCK::ev_iter", n_vars);
  static double* aug;
  double lambda;
  if(!aug){ untested(); // quick hack!
    assert(_sim->_aa.size()==n_vars);
    aug = new double[n_vars*2+1];
  }
#ifdef DO_TRACE
  for (unsigned i=0; i < 1+2*n_vars; i++){ aug[i]=i; }
#endif
  _sim->_aa = -(_sim->_acx.real());
  assert(_sim->_aa.size()==n_vars);
  double tmp[n_vars+1];
  double xiin[n_vars+1];
  double norm = 0;
  for (unsigned i=0; i < n_vars; i++) {
    stream >> tmp[i];
    trace1("SOCK::ev_iter got", tmp[i]);
    norm+=tmp[i]*tmp[i];
  }
#ifdef normalize
  norm=1./sqrt(norm);
#endif
  for (unsigned i=0; i < n_vars; i++) {
#ifdef normalize
    tmp[i]*=norm;
#endif
    aug[n_vars*2-i] = tmp[i];
    xiin[i]=tmp[i];
  }
#ifdef normalize
  aug[n_vars] = 0;
#else
  aug[n_vars] = norm-1;
#endif
  trace1("SOCK filled xi aa", _sim->_aa);
  untested();
  const BSMATRIX<double> C = _sim->_acx.imag();
  C.rmul(aug-1,xiin-1);
  stream >> lambda;

#ifdef rescale
  for (unsigned i=0; i < n_vars; i++) {aug[i]*=lambda;}
#endif

  _sim->_aa+= lambda * _sim->_acx.imag();

  trace0("SOCK::ev_iter rmul");
  _sim->_aa.rmul(tmp-1,xiin-1);
  untested();
  tmp[n_vars] = 0;
  for (unsigned i=0; i <= n_vars; i++){ untested();
    trace2("rhs",i,tmp[i]);
  }

  _sim->_aa.augment(aug);
  assert(_sim->_aa.size()==n_vars+1);
  trace1("SOCK augmented aa", _sim->_aa);
  _sim->_aa.lu_decomp();
  trace1("SOCK decomposed aa", _sim->_aa);
  untested();

  _sim->_aa.fbsub(tmp-1);
  for (unsigned i=0; i <= n_vars; i++){
    trace2("res",i,tmp[i]);
  }
  untested();
  trace0("SOCK::ev_iter done");
  _sim->_aa.deaugment();
  assert(_sim->_aa.size()==n_vars);
  stream << ((uint16_t)_error); stream.pad(6);
  double resnormsq = 0;
  for (unsigned i=0; i < n_vars; i++) { untested();
    resnormsq+=tmp[i]*tmp[i];
    stream << xiin[i]-tmp[i];
  }
#ifdef rescale
  resnormsq+=lambda*tmp[n_vars]*lambda*tmp[n_vars];
  stream << lambda*(1-tmp[n_vars]);
#else
  resnormsq+=tmp[n_vars]*tmp[n_vars];
  stream << lambda-tmp[n_vars];
#endif
  stream << resnormsq;
  stream << SocketStream::eol;
  trace0("SOCK::ev_iter sent");
}
/*--------------------------------------------------------------------------*/
void SOCK::set_param()
{ untested();
  string tmp[2];
  stream >> tmp[0];
  stream >> tmp[1];
  PARAM_LIST* pl = CARD_LIST::card_list.params();
  string paramcommand = tmp[0]+"={"+tmp[1]+"}";
  trace3("set_param", tmp[0], tmp[1], paramcommand);
  CS cmd(CS::_STRING, paramcommand);
  pl->parse(cmd);
  CARD_LIST::card_list.precalc_first(); // BUG? params don't propagate into sckts otherwise
  CARD_LIST::card_list.precalc_last();
}
/*--------------------------------------------------------------------------*/
enum{
  sOK=0,          // ok, returned time step indicates next prefered
  sNOCONV=1,      // no convergence for second or 3. newton iteration: hope
                  // for success with shorter time step
  sTRUNC=2,       // Truncation error (time step not accepted (to long!))
                  // restart with shorter returned time step
  sFIRSTNOCONV=3, // no convergence for first newton iteration after cons
                  // time step has no meaning, but reduction of time step should
                  // help
  sTROUBLE=99     // exception in solve.
};
/*--------------------------------------------------------------------------*/
// init=0 continue with ruse of all old data (do 1 step)
//      1 restore op and do ~5 steps (could be used directly after a not accepted
//        time step)
//      2 do cons_op followed by do ~5 step (depending on Euler,Gear,
//        Trapezoid, options)
//
// return:
// ret (see above)
// dt. positive timestep, use instead in next call.
//
//
unsigned SOCK::transtep(unsigned init, double dt)
{
  unsigned ret = sOK;
  _error = 0;
  trace3("SOCK::transtep", n_vars, init, dt);
  _sim->set_command_tran();
  _sim->restore_voltages(); //     _vt1[ii] = _v0[ii] = vdc[ii];
  _sim->_uic = false;
  _sim->_phase = p_RESTORE;

  unsigned stepno=-1;
  double reftime = _sim->_time0;

  if (dt==0.) {
    dt = OPT::dtmin;
  }else{
  }

  assert(dt);

  if (init==0) { itested();         // continue
    stepno = 1;
    ::error(bDEBUG, "transient step, at %f, dt %f\n", _sim->_time0, dt);
    _sim->_dt0 = dt;
    assert(_sim->_time0);
    _sim->_time0 += _sim->_dt0;
    // CARD_LIST::card_list.tr_restore();
  } else if (init==1) { itested(); // restore
    assert(_sim->_time0 > 0);
    //    CARD_LIST::card_list.tr_restore();
    stepno = 1;
    _sim->_dt0 = dt/stepno;
    ::error(bDEBUG, "restore step, dt %e\n", dt);
    _sim->_time0 += _sim->_dt0;
  } else if (init==2) {
    _sim->_time0 = reftime = 0;
    verakons(); // do we need it ?
    stepno = 5;
    _sim->_dt0 = dt/stepno;
    ::error(bDEBUG, "initial step, dt %e\n", dt);
    _sim->_time0 += _sim->_dt0;
  } else { unreachable();
  }


  trace2("SOCK::transtep ", _sim->_phase, _sim->_mode);

  advance_time(); // sink _time0 into devices
  _sim->set_command_tran();

  if (printlist().size()) {
    outdata(reftime);
  }

  _sim->_phase = p_TRAN;

  //_sim->_bypass_ok = false;
  assert(_sim->_dt0 == dt/stepno);
  //_sim->_genout = gen();

  assert(_sim->analysis_is_tran());


  bool tr_converged;
  for (unsigned i = stepno; i>0; --i) {
    tr_converged = false;
    try {
      tr_converged = solve(OPT::TRHIGH, _trace);
    }catch (Exception e) { incomplete();
      ret = sTROUBLE;
    }
    if (!tr_converged) { itested();
      ::error(bDEBUG, "s_sock::transtep did not converge\n");
      if (i==stepno) { untested();
        ret = sFIRSTNOCONV;
        break;
      } else { untested();
        ret = sNOCONV;
	break;
      }
    }else{
    }

    if(stepno>1) {
      /// uuh ooh HACK.
      assert(ret!=sFIRSTNOCONV);
      ::status.accept.start(); // accounting
      _sim->set_limit();
      CARD_LIST::card_list.tr_accept(); // accept solution in all components
      ::status.accept.stop();
    }
    if (printlist().size()) {
      outdata(_sim->_time0);
    }
    if(i>1){
      _sim->_time0 += _sim->_dt0;
    }
  }

  TIME_PAIR time_by = CARD_LIST::card_list.tr_review();
  double time_by_error_estimate = time_by._error_estimate;
  assert(time_by_error_estimate>=0);

  // if(!tr_converged?) { ... } else
  if (time_by_error_estimate > _sim->_time0) {
    dt = time_by_error_estimate - _sim->_time0;
    ::status.accept.start();
    _sim->set_limit();
    CARD_LIST::card_list.tr_accept();
    ::status.accept.stop();
    _sim->keep_voltages(); //  vdc  = v0
    assert(dt>0);
  }else{ untested();
    double t1 = _sim->_time0 - _sim->_dt0;
    dt = time_by_error_estimate - t1;
    assert(dt);
    ::error(bTRACE,"sock step not accepted at %f, need %f\n", reftime, dt);
    ret = sTRUNC;
    assert (dt>0); // hopefully.
    assert (dt < _sim->_dt0); // hopefully.
    dt *= stepno; // stepno is valid for next step.
    _sim->_time0 = reftime;
    tr_reject();
  }

#ifndef NDEBUG
  for (unsigned i=1; i <= n_vars; i++) {
    if(isnan(_sim->_i[i])||isnan(_sim->_vdcstack.top()[i]) ) {
      _error = 1;
      assert(!tr_converged);
    }
  }
#endif

  _dthack = dt;
  return ret; // FIXME: proper status.
}
/*--------------------------------------------------------------------------*/
void SOCK::tr_reject()
{
  ::status.accept.start();
  _sim->_acceptq.clear();
  ++tr_steps_rejected_;
  ::status.accept.stop();
}
/*--------------------------------------------------------------------------*/
void SOCK::verakons_reply()
{
  stream << ((uint16_t)_error); stream.pad(6);
  for (unsigned i=1; i <= n_vars; i++) {
    trace1("SOCK::kons_reply", _sim->vdc()[i]);
    stream << _sim->vdc()[i];
  }

  //
  // now i contains the negative sum of the cap value.
  // does this break anything?

  for (unsigned i=1; i <= n_vars; i++)         /* Q-Punkt == I == RS */
  {
    trace1("SOCK::verakons_reply  dqdt", _sim->_i[i]);
    stream << _sim->_i[i];
  }

  if(_dump_matrix){
    _out << "u,dq=i: \n";
    for (unsigned i=0; i <= n_vars; i++)         /* Variablen-Werte */
    {
      _out << _sim->vdc()[i] << "," << _sim->_i[i] <<  "\n" ;
    }
  }

  trace0("SOCK::verakons_reply  ac_snap");
  ac_snapshot(); // FIXME: to verakons()

  send_matrix();

  assert(stream.bufsize() >= total);
  stream << SocketStream::eol;

//  for( unsigned i = 0; i < _caplist.size(); i++)
//  {
//    trace1("unlatching", i);
//
//    dashier umgekehrt.
//    commonstash[i] = _caplist[i]->common();
//    detach(_caplist[i]);
//
//    COMMON_COMPONENT* c = bm_dispatcher.clone("eval_bm_value");
//    COMMON_COMPONENT* dc = c->deflate();
//
//    _caplist[i]->attach_common(dc);
//
//    _caplist[i]->keep_ic(); // latch voltage applied to _v0
//
//
//  }

}

void SOCK::transtep_reply(unsigned ret, bool eol)
{
  _error = 0;
  stream << _error; stream.pad(6);
  stream << (uint64_t)ret;
  stream << _dthack;

  for (unsigned i=1; i <= n_vars; i++) {
    stream << _sim->_vdcstack.top()[i];
  }

  // OUCH, move to main loop.
  if(eol) {
    stream << SocketStream::eol;
  }
}

void SOCK::transtep_gc_reply(unsigned ret)
{ itested();
  transtep_reply(ret, false);

  for (unsigned i=1; i <= n_vars; i++)
  {
    trace1("SOCK::transtep_gc_reply  dqdt", _sim->_i[i]);
    stream << _sim->_i[i];
  }

  trace0("SOCK::transtep_gc_reply  ac_snap");
  ac_snapshot(); // BUG: calls tr_begin
  send_matrix();
  assert(stream.bufsize() >= total);
  stream << SocketStream::eol;
}
/*--------------------------------------------------------------------------*/
void SOCK::verainit_reply()
{
  trace4("SOCK::verainit_reply ", _error, n_vars, length, var_namen_total_size);
  stream << _error;   stream.pad(6);
  stream << int32_t (n_vars);  stream.pad(4);
  stream << int32_t (var_namen_total_size);  stream.pad(4);

  assert(stream.tcur() == 24);

  for (unsigned i = 0; i < n_vars; i++) {
    trace2("SOCK::verainit_reply",  i, var_namen_arr[i]);
    putstring8( &stream, var_namen_arr[i]);
    stream << '\t'; stream.pad(7);
  }

//good idea?? not yet. vera insists on branch node
//
//  for (unsigned i=0; i < n_inputs; i++)
//  {
//    trace1("putting name " +  input_names[i], i);
//    putstring8( &stream, input_names[i]);
//    stream << '\t'; stream.pad(7);
//  }
  stream.flush();

  trace0("done verainit_reply");
}
/*--------------------------------------------------------------------------*/
void SOCK::veraop_reply()
{
  assert(n_vars==_sim->_total_nodes);
  trace1("SOCK::veraop_reply ", n_vars);
  stream << _error; stream.pad(6);
  for (unsigned i=1; i <= n_vars; i++)         /* Variablen-Werte */
  {
    stream << _sim->vdc()[i];
  }
  if(_dump_matrix){
    _out << "i,u: \n";
    for (unsigned i=0; i <= n_vars; i++)         /* Variablen-Werte */
    {
      _out << _sim->_i[i] << "," <<  _sim->vdc()[i] << "\n" ;
    }
  }

  trace0("SOCK::veraop_reply taking ac snapshot (matrix only?).");
  ac_snapshot();
  send_matrix();
  stream << SocketStream::eol;

  trace0("SOCK::veraop_reply sent");

  total = 2*n_vars_square+n_vars+1;
  assert(stream.bufsize() >= total);
  if (printlevel >= 1)
  {
    userinfo(1,"vera_titan_ak","Sende: Error %i Framenumber %i, Laenge %i\n",
        _error,frame_number,total);
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
void SOCK::send_matrix()
{
  const BSMATRIX<double> G = _sim->_acx.real();
  const BSMATRIX<double> C = _sim->_acx.imag();
  if(_dump_matrix){
     _out << "G\n" << G << "\n";
     _out << "C\n" << C << "\n";
  }
  trace0("SOCK::send_matrix G");
  for (unsigned i=1; i <= n_vars; i++){
    for (unsigned j=1; j <= n_vars; j++) {
      stream << G.s(j,i);
    }
  }
  trace0("SOCK::send_matrix C");
  for (unsigned i=1; i <= n_vars; i++){
    for (unsigned j=1; j <= n_vars; j++) {
      stream << C.s(j,i);
    }
  }
}
/*--------------------------------------------------------------------------*/
void SOCK::cap_prepare(void)
{ untested();
  trace1("SOCK::cap_prepare", _caplist.size() );
  assert(!_capstash);
  _capstash = new CARDSTASH[_caplist.size()];

  for (unsigned ii = 0;  ii < _caplist.size();  ++ii) { untested();
    _caplist[ii]->inc_probes();			// we need to keep track of it
    _capstash[ii] = _caplist[ii];			// stash the std value

    if(_caplist[ii]->has_common()){ untested();
      trace1("SOCK::cap_prepare. have common", _caplist[ii]->value());
//      _caplist[ii]->set_value(_caplist[ii]->value(), 0); // zap out extensions
//      _caplist[ii]->set_constant(false);		 // update HACK?
//
//    all devices need individual commons.
      _caplist[ii]->attach_common(_caplist[ii]->common()->clone());
      assert(_caplist[ii]->has_common());
    }else{ untested();
      //      _sweepval[ii] = _zap[ii]->set__value();	// point to value to patch
      COMMON_COMPONENT* c = bm_dispatcher.clone("eval_bm_value");
      double capval = _caplist[ii]->value();
      trace2("SOCK::cap_prepare, attaching common", _caplist[ii]->long_label(), capval);
      c->set_value(capval);
      COMMON_COMPONENT* dc = c->deflate();
      assert(dc);
      //
      _caplist[ii]->set_value(_caplist[ii]->value(), dc);	// zap out extensions
      _caplist[ii]->set_constant(false);		// so it will be updated
      trace1("SOCK::cap_prepare calling precalc_first", _caplist[ii]->long_label());
      _caplist[ii]->precalc_first();
      trace0("SOCK::cap_prepare calling precalc_last");
      _caplist[ii]->precalc_last();
      trace0("SOCK::cap_prepare calling trbegin");
      _caplist[ii]->tr_begin();
    }
  }
}
/*--------------------------------------------------------------------------*/
void SOCK::cap_reset(void)
{
  trace0("SOCK::cap_reset");
  for (unsigned ii = 0;  ii < _caplist.size();  ++ii) {
      _capstash[ii].restore();
      _caplist[ii]->dec_probes();
//      _caplist[ii]->precalc_first();
//      _caplist[ii]->precalc_last();
  }
  delete[] _capstash;
  _capstash = NULL;
  trace0("SOCK::cap_reset done");
}
/*--------------------------------------------------------------------------*/
SOCK::~SOCK()
{
  trace0("SOCK::~SOCK()");
}
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
