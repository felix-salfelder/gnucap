

#include "s_tt.h"

namespace TT {

void TTT::options(CS& Cmd)
{
	trace0(( "TTT::options rest ||| " +Cmd.tail() ).c_str());
	_out = TTT::mstdout;
	_out.reset();
	bool tr = false;
	bool cont = false;
	_sim->_age = true;
	bool noage = false;
	_no_act = false;
	_evt = true;
	// _stepmode = LIN;
	// cout << "stepmode " << _stepmode << "\n";

	_sim->_temp_c = OPT::temp_c;
	_trace = tNONE;
	_alarm = aREPORT;
	_agemode = amALWAYS;
	bool depcont = false;

	_power_down = _quiet = false;
	unsigned here = Cmd.cursor();
	do{
		ONE_OF
			|| common_options(Cmd)
			|| Get(Cmd, "noage", &noage)
			|| Get(Cmd, "p{owerdown}", &_power_down)
			|| Get(Cmd, "pd",          &_power_down)
			|| Get(Cmd, "evt",         &_evt)
			|| Get(Cmd, "ev{ents}",    &_evt)
			|| Get(Cmd, "new",         &_new)
			|| Get(Cmd, "n{o-act}",    &_no_act)
			|| Get(Cmd, "q{uiet}",     &_quiet)
			|| Get(Cmd, "c{ont}",      &cont) // continue from adp_node values.
			|| Get(Cmd, "cont_dc",     &depcont) // don't recompute dc point
			|| Get(Cmd, "uic",         &_sim->_uic) // wrong/incomplete
//			|| (Get(Cmd, "+",           &_Tstep) && 
//				 ( Cmd.is_float()	&& (Cmd >> _Tstep ) && (_stepmode = LIN)))
			|| (Get(Cmd, "*",           &_Tstep) &&  (_stepmode = tts_MUL) &&
				 ( Cmd.is_float()	&& (Cmd >> _Tstep ) ))
			|| (Cmd.umatch("agemode {=}") &&
					(ONE_OF
					 || Set(Cmd, "n{one}",     &_agemode, amNONE)
					 || Set(Cmd, "o{nce}",     &_agemode, amONCE)
					 || Set(Cmd, "b{nce}",     &_agemode, amBYPASS)
					 || Set(Cmd, "a{lways}",   &_agemode, amALWAYS)
					)
				)
			|| (Cmd.umatch("alarm {=}") &&
					(ONE_OF
					 || Set(Cmd, "r{eport}",   &_alarm, aREPORT)
					 || Set(Cmd, "n{one}",     &_alarm, aNONE)
					 || Set(Cmd, "red{irect}", &_alarm, aREDIR)
					 || Set(Cmd, "a{bort}",    &_alarm, aABORT)
					)
				)
			|| (Cmd.umatch("tr{ace} {=}") &&
					(ONE_OF
					 || Set(Cmd, "n{one}",      &_trace, tNONE)
					 || Set(Cmd, "o{ff}",       &_trace, tNONE)
					 || Set(Cmd, "w{arnings}",  &_trace, tUNDER)
					 || Set(Cmd, "g{uess}",     &_trace, tGUESS)
					 || Set(Cmd, "a{lltime}",   &_trace, tALLTIME)    //print also if stepcause!=user
					 || Set(Cmd, "r{ejected}",  &_trace, tREJECTED)
					 || Set(Cmd, "i{terations}",&_trace, tITERATION)
					 || Set(Cmd, "v{erbose}",   &_trace, tVERBOSE)
					 || Set(Cmd, "d{ebug}",     &_trace, tDEBUG)
					 || Cmd.warn(bWARNING, "need none, off, warnings, alltime, "
						 "rejected, iterations, verbose")
					)
				)
			;
		if (!( Get(Cmd , "tran",	   &tr) ) ) {
			// tr means 'trace'
			_out.outset(Cmd);
		}
	}while (Cmd.more() && !Cmd.stuck(&here) && !tr);

	initio(_out);
	TRANSIENT::options(Cmd);
	if (noage) {
		_sim->_age = false;
	}
	if (_stepmode == tts_MUL && _Tstep<=1 ) { untested();
		throw Exception("Multiplier too small: %f", double(_Tstep));
	}

	if (_new){
		_cont_tt = false;
		if(cont) {
			incomplete();
		}
		cont = false;
	}
	if(cont){
		trace0("continuing with aging state...");
		_cont_tt = true;
	}

	_dtmax_in.e_val(BIGBIG, _scope);
	// _dTmin_in.e_val(OPT::dTmin, _scope);
	_dtratio_in.e_val(OPT::dtratio, _scope);
	_skip_in.e_val(1, _scope);

	trace3("TTT::options done", _stepmode, _Tstep, _cont_dc);

	if(depcont){
		error(bWARNING, "cont_dc deprecated. use tran cont to continue tran\n");
		_cont_dc = true;
	}
}
/*--------------------------------------------------------------------------*/
void TTT::setup(CS& Cmd)
{
	trace3("TTT::setup", _cont_tt, Cmd.tail(), Cmd.fullstring());
	Cmd.reset(0);
	Cmd.skip1('.');

	if (Cmd.umatch("tw")) {
		return setup_tw(Cmd); // old setup, weird, inconsistent syntax
	} else if (Cmd.umatch("ttr{an} ")) {
	} else if (Cmd.umatch("tt")) { untested();
		return setup_tw(Cmd); // old setup, weird, inconsistent syntax
	} else { unreachable();
		trace2("problem", Cmd.cursor(), Cmd.fullstring());
	}

	_Tstart.e_val(NOT_INPUT, _scope);
	_Tstop.e_val(NOT_INPUT, _scope);
	_Tstep.e_val(NOT_INPUT, _scope);
	_tstart.e_val(NOT_INPUT, _scope);
	_tstop.e_val(NOT_INPUT, _scope);
	_tstep.e_val(NOT_INPUT, _scope);

	_new = false;
	_cont = true;
	_cont_tt = true;

	if (Cmd.match1("'\"({") || Cmd.is_pfloat()) {
		PARAMETER<double> arg1, arg2, arg3, arg4, arg5, arg6;
		Cmd >> arg1;
		arg1.e_val(0.0,_scope);
		if (Cmd.match1("'\"({") || Cmd.is_float()) {
			Cmd >> arg2;
			arg2.e_val(0.0,_scope);
		}else{
		}
		if (Cmd.match1("'\"({") || Cmd.is_float()) {
			Cmd >> arg3;
			arg3.e_val(0.0,_scope);
		}else{
		}
		if (Cmd.match1("'\"({") || Cmd.is_float()) {
			Cmd >> arg4;
			arg4.e_val(0.0,_scope);
		}else{
		}
		if (Cmd.match1("'\"({") || Cmd.is_float()) {
			Cmd >> arg5;
			arg5.e_val(0.0,_scope);
		}else{
		}
		if (Cmd.match1("'\"({") || Cmd.is_float()) {
			Cmd >> arg6;
			arg6.e_val(0.0,_scope);
		}else{
		}

		trace5(("TTT::setup args " + std::string(Cmd)).c_str(), arg1, arg2 , arg3 , arg4, arg5);

		if (arg5.has_hard_value()) {
			trace0("5 args");
			_stepmode = tts_LIN;
			// 5 args  TTstart TTstop tstep tstop TTstep
			// 5 args  TTstart TTstop tstop tstep TTstep
			assert(arg5.has_hard_value());
			assert(arg4.has_hard_value());
			assert(arg3.has_hard_value());
			assert(arg2.has_hard_value());
			assert(arg1.has_hard_value());
			arg1.e_val(0.,_scope);
			arg2.e_val(0.,_scope);
			arg3.e_val(0.,_scope);
			arg4.e_val(0.,_scope);
			arg5.e_val(0.,_scope);
			_sim->_last_Time = .0; // ?


			if( arg3 < arg4 ) {
				_tstep  = arg3;
				_tstop  = arg4;
			} else {
				_tstep  = arg4;
				_tstop  = arg3;
			}

			_Tstart = arg1;
			_Tstep  = arg5;
			_Tstop  = arg2;

			if (double(_Tstart) == 0) {
				_cont_tt = false;
			} else {
			}
		} else if (arg4.has_hard_value()) {
			trace0("4 args");
			assert(arg3.has_hard_value());
			assert(arg2.has_hard_value());
			assert(arg1.has_hard_value());
			arg1.e_val(0.,_scope);
			arg2.e_val(0.,_scope);
			arg3.e_val(0.,_scope);
			arg4.e_val(0.,_scope);

//			_sim->_last_Time = _arg1;
			if( arg1 < arg2 ) {
				// 4 args: TTstart TTstop tstep tstop
				// 4 args: TTstart TTstop tstop tstep
				if(_Tstep==NOT_INPUT){
					_Tstart = 0.;
					_Tstep = arg1;
				}else if(_Tstep==0.){ untested();
					_Tstart = 0.;
					_Tstep = arg1;
				}else{
					assert(_Tstep>0.);
					_Tstart = arg1;
				}
				_Tstop = arg2;
				if( arg3 < arg4 ) {
					_tstep = arg3;
					_tstop = arg4;
				} else {
					_tstep = arg4;
					_tstop = arg3;
				}
			} else { untested();
				incomplete();
				// 4 args: TTstop tstep tstop TTstep
			}
			/*
			if (arg4.has_hard_value()) {
				arg4.e_val(0.,_scope);
				_Tstep = arg4;
				_stepmode = LIN;
			} else {
				_Tstep = OPT::ttstepgrow;
				_stepmode = MUL;
			}
			*/
			if (!_stepmode){
				_stepmode = tts_LIN;
			} else {
			}
		} else if (arg3.has_hard_value()) {
			trace0("3 args");
			assert(arg3.has_hard_value());
			assert(arg2.has_hard_value());
			assert(arg1.has_hard_value());
			arg1.e_val(0.,_scope);
			arg2.e_val(0.,_scope);
			arg3.e_val(0.,_scope);

			_Tstop = arg1;
			_tstep = arg2;
			_tstop = arg3;

			if( arg1 < arg2 ) { untested();
				// 3 args: TTstart TTstop tstop (tstep from previous run or =tstop)
				_Tstart = arg1;
				_Tstop = arg2;
				_tstop = arg3;
				_tstep = _tstop; incomplete(); // bug, take from last if exists
			} else if ( arg2 < arg3 ) {
				// 3 args: TTstop tstep tstop   (TTstart=0 or previous)
				_Tstart = _sim->_last_Time;
				_tstep = arg2;
				_tstop = arg3;
			} else {
				// 3 args: TTstop tstop tstep   (TTstart=0 or previous)
				_Tstart = _sim->_last_Time;
				_tstep = arg3;
				_tstop = arg2;
			}

			if (!_Tstep.has_hard_value()) {
				_Tstep = _tstop;
				_stepmode = tts_LIN;
			} else {
			}
			assert((double)_tstep!=0 || !_tstop);
			if (!_tstop) {
			} else {
			}

			_cont_tt = false;
			trace6("TTT::setup ", _tstep, _tstop, _Tstep, _Tstop, _sim->last_time(), _sim->last_Time());

			if(_sim->last_time() && !_sim->last_Time()){
				// a transient has been run, but no tt yet
				trace1("TTT::options setting cont_tt", _sim->last_time());
				_cont_dc = true;
				_new = true;
				_sim->_last_Time = _sim->last_time();
			}
			_sim->tr_reset();

		} else if (arg2.has_hard_value() ) {
			// Tstart Tstop, previous tstop/tstep
			trace0("TTT::setup have 2");
			_Tstart = _sim->_last_Time;
			_sim->_time0 = 0;
			_sim->tr_reset();
			if ((double)_Tstart == 0){
				trace0("TTT::setup latching tr times");
				_tstep = arg1;
				_tstop = arg2;
				_Tstop = 0;

			}else{
				trace1("TTT::setup ran already", (double)_Tstart );
				if((double)arg1==0){
					_Tstop  = arg2; 
					_Tstart =0;
				}else if(arg1<arg2){
					_Tstop  = arg2; 
					_Tstep = arg1;
				}else{
					_Tstop  = arg1; 
					_Tstep = arg2;
				}
			}

			_cont_tt = true;
			_cont = true;
			trace4("TTT::setup 2 args ", _tstep, _tstop , _Tstep , _Tstop);

		} else if (arg1.has_hard_value() ) {
			trace1("TTT::setup same tr, new Tend", _sim->_last_Time);
			_Tstart = _sim->_last_Time;
			_Tstop  = arg1; // as tran
			_sim->_time0 = 0;
			_sim->tr_reset();

			// to trigger prints... (hack?)
			if(double(_Tstop) == 0) _Tstop = double( _Tstart );

			if((!_Tstep.has_hard_value() )|| ((double)_Tstep == 0)) {
				trace1("set Tstep ", _Tstop);
				_Tstep=_Tstop;
			} else {
				trace1("set Tstep ", _Tstep);
			}
			if (_Tstart!=0) { // obsolete?
				_cont_tt = true;
				_cont = true;
			}

		} else {
			unreachable();
			assert (!arg1.has_hard_value());   // for now...

			double oldrange = _Tstop - _Tstart;
			_Tstart = _sim->_last_Time;
			_Tstop  = _sim->_last_Time + oldrange;
		}
	} else {

	  	/* no args */
		// std::cerr << "setup ttt -- no args\n";
		double oldrange = _Tstop - _Tstart;
		_Tstart = _sim->_last_Time;
		_Tstop  = _sim->_last_Time + oldrange;

		if(_sim->_last_Time==0){
			trace1("TTT::setup no args at beginning", _cont_tt);
			_Tstop=0;
			_tstep=0;
			_tstop=0;

		}
		if (!_stepmode){
			_stepmode = tts_LIN;
		}else{ untested();
		}
	}

	if (_tstep>_tstop){
		error(bWARNING, "_tstep > _tstop. really?\n");
	}

	options(Cmd);

	_tstart.e_val(0., _scope);
	_tstop.e_val(NOT_INPUT, _scope);
	_tstep.e_val(NOT_INPUT, _scope);
	_Tstop.e_val(NOT_INPUT, _scope);
	_Tstep.e_val(NOT_INPUT, _scope);

	//Time1 = 
	_sim->_Time0 = _Tstart;
	_Time1 = _Tstart;


	if (!_tstep.has_good_value()) {
		throw Exception("transient: Time step is required");
	}else if (_tstep==0. && _tstop ) {
		untested();
		throw Exception("Time step == 0 while tend");
	}else{
	}

//		if (_dtmax_in.has_hard_value()) {
//			_dtmax = _dtmax_in;
//		}else if (_skip_in.has_hard_value()) {
//			_dtmax = _tstep / double(_skip_in);
//		}else{
//			_dtmax = std::min(_dtmax_in, _tstep);
//		}

	_dTmin= _tstop;
	_sim->_dTmin= _tstop * .5;

	if (_dTmin_in.has_hard_value()) { untested();
		_dTmin = _dTmin_in;
	}else if (_dtratio_in.has_hard_value()) { untested();
		_dTmin = _dTmax / _dTratio_in;
	}else{
		// use larger of soft values
		// _dTmin = std::max(double(_dTmin_in), _dTmax/_dTratio_in);
		// _dTmin=0.5; // HACK
	}


	if  (_Tstart < _sim->_last_Time) {
		_cont_tt = false;
		_Time1 = _sim->_Time0 = _Tstart;
	}else if (_Tstart < _sim->_last_Time  ||  _sim->_last_Time <= 0.) {
		//    _out << "* last_Time " << _sim->_last_Time << "\n";
		trace3("TTT::setup no cont ", _Tstart, _sim->_last_Time, _cont_tt );
		//_cont_tt = false;
		_Time1 = _sim->_Time0 = _Tstart;
	}else{
		_cont_tt = true;
		_Time1 = _sim->_Time0 = _sim->_last_Time;
	}

	// TRANSIENT setup

	_tstart.e_val(0., _scope);
	_tstop.e_val(NOT_INPUT, _scope);
	_tstep.e_val(NOT_INPUT, _scope);

	if  (_cont_tt){
		_cont = false;

	}else{
		_cont = false;
	}

	_time1 = _sim->_time0 = 0.;
	//{}else{
	//  untested();
	// _cont = true;
	//  time1 = _sim->_time0 = _sim->_last_time;
	// }
	_sim->_freq = ((_tstop > _tstart) ? (1 / (_tstop - _tstart)) : (0.));

	if (!_tstep.has_good_value()) {
		throw Exception("transient: time step is required");
	}else if (_tstep==0. && _tstop ) { itested();
		throw Exception("time step == 0 and tstop");
	}else{
	}

	if (_dtmax_in.has_hard_value()) {
		_dtmax = _dtmax_in;
	}else if (_skip_in.has_hard_value()) {
		_dtmax = _tstep / double(_skip_in);
	}else{
		_dtmax = std::min(_dtmax_in, _tstep);
	}

	if (_dtmin_in.has_hard_value()) {
		_sim->_dtmin = _dtmin_in;
	}else if (_dtratio_in.has_hard_value()) {
		_sim->_dtmin = _dtmax / _dtratio_in;
	}else{
		// use larger of soft values
		_sim->_dtmin = std::max(double(_dtmin_in), _dtmax/_dtratio_in);
	}

	assert(_stepmode);

	steps_total_out_ = (uint_t) (1 + ceil( ( (_tstop - _tstart ) / _tstep ) ));
	trace7( "TTT::setup done ",  _sim->_Time0, _Tstart, _Tstep ,_stepmode, _cont, _cont_tt, _Tstop );
	trace7( "TTT::setup done ",  steps_total_out_ , (double)_tstep , _tstop ,_tstart, _cont, _cont_tt, (double)_Tstop );
	allocate();
}
/*--------------------------------------------------------------------------*/
void TTT::setup_tw(CS& Cmd)
{
	trace3("TTT::setup_old", _cont_tt, Cmd.tail(), Cmd.fullstring());

	_tstart.e_val(NOT_INPUT, _scope);
	_tstop.e_val(NOT_INPUT, _scope);
	_tstep.e_val(NOT_INPUT, _scope);
	_Tstop.e_val(NOT_INPUT, _scope);
	_Tstep.e_val(NOT_INPUT, _scope);

	_new = false;
	_cont = true;
	_cont_tt = true;
	_Tstep = 2.0;

	if (Cmd.match1("'\"({") || Cmd.is_pfloat()) {
		PARAMETER<double> arg1, arg2, arg3, arg4, arg5, arg6;
		Cmd >> arg1;
		arg1.e_val(0.0,_scope);
		if (Cmd.match1("'\"({") || Cmd.is_float()) {
			Cmd >> arg2;
			arg2.e_val(0.0,_scope);
		}else{
		}
		if (Cmd.match1("'\"({") || Cmd.is_float()) {
			Cmd >> arg3;
			arg3.e_val(0.0,_scope);
		}else{
		}
		if (Cmd.match1("'\"({") || Cmd.is_float()) {
			Cmd >> arg4;
			arg4.e_val(0.0,_scope);
		}else{
		}
		if (Cmd.match1("'\"({") || Cmd.is_float()) {
			Cmd >> arg5;
			arg5.e_val(0.0,_scope);
		}else{
		}
		if (Cmd.match1("'\"({") || Cmd.is_float()) {
			Cmd >> arg6;
			arg6.e_val(0.0,_scope);
		}else{
		}

		trace5(("TTT::setup_tw args " + std::string(Cmd)).c_str(), arg1, arg2 , arg3 , arg4, arg5);

		if (arg5.has_hard_value()) {	    // 5 args  tt tt TTstart TTstop TTstep
			trace0("TTT::setup_tw have 5");
			untested();
			assert(arg5.has_hard_value());
			assert(arg4.has_hard_value());
			assert(arg3.has_hard_value());
			assert(arg2.has_hard_value());
			assert(arg1.has_hard_value());
			arg1.e_val(0.,_scope);
			arg2.e_val(0.,_scope);
			arg3.e_val(0.,_scope);
			arg4.e_val(0.,_scope);
			_sim->_last_Time = .0; // ?
			_sim->tr_reset();

			_tstep  = arg1;
			_tstop  = arg2;

			_Tstart = arg3;
			_Tstep  = arg4;
			_Tstop  = arg5;
			assert((double)_tstep!=0 || !_tstop);
			if (!_tstop) { untested(); }

			if(double( _Tstart) == 0) 
				_cont_tt = false;

			trace4("TTT::setup_tw ", _tstep, _tstop , _Tstep , _Tstop);

		} else if (arg3.has_hard_value()) {	    // 4 args: tt tt TT [TTstep]
			trace0("TTT::setup_tw have 3");
			assert(arg3.has_hard_value());
			assert(arg2.has_hard_value());
			assert(arg1.has_hard_value());
			arg1.e_val(0.,_scope);
			arg2.e_val(0.,_scope);
			arg3.e_val(0.,_scope);
			_sim->_last_Time = .0;

			_tstep  = arg1;
			_tstop  = arg2;

			if ( _tstep>_tstop ) {
				double x=_tstep;
				_tstep=_tstop;
				_tstop=x;
			} else {
			}
			_Tstop = arg3;
			if (arg4.has_hard_value()) {
				arg4.e_val(0.,_scope);
				_Tstep = arg4;
				_stepmode = tts_LIN;
			} else {
				_Tstep = OPT::ttstepgrow;
				_stepmode = tts_MUL;
			}
			_Tstart = 0; //HACK?
			assert((double)_tstep!=0 || !_tstop);
			if (!_tstop) { untested();
			} else {
			}

			_cont_tt = false;
			trace6("TTT::setup_tw ", _tstep, _tstop, _Tstep, _Tstop, _sim->last_time(), _sim->last_Time());

			if(_sim->last_time() && !_sim->last_Time()){
				// a transient has been run, but no tt yet
				trace1("TTT::options setting cont_tt", _sim->last_time());
				_cont_dc = true;
				_new = true;
				_sim->_last_Time = _sim->last_time();
			}
			_sim->tr_reset();

		} else if (arg2.has_hard_value() ) {
			trace0("TTT::setup_tw have 2");
			_Tstart = _sim->_last_Time;
			_sim->_time0 = 0;
			_sim->tr_reset();
			if ((double)_Tstart == 0){
				trace0("TTT::setup_tw latching tr times");
				_tstep = arg1;
				_tstop = arg2;
				_Tstop = 0;

			}else{
				trace1("TTT::setup_tw ran already", (double)_Tstart );
				if((double)arg1==0){
					_Tstop  = arg2; 
					_Tstart =0;
				}else if(arg1<arg2){
					_Tstop  = arg2; 
					_Tstep = arg1;
				}else{
					_Tstop  = arg1; 
					_Tstep = arg2;
				}
			}

			_cont_tt = true;
			_cont = true; untested();
			trace4("TTT::setup_tw 2 args ", _tstep, _tstop , _Tstep , _Tstop);

		} else if (arg1.has_hard_value() ) {
			trace1("TTT::setup_tw same tr, new Tend", _sim->_last_Time);
			_Tstart = _sim->_last_Time;
			_Tstop  = arg1; // as tran
			_sim->_time0 = 0;
			_sim->tr_reset();

			// to trigger prints... (hack?)
			if(double(_Tstop) == 0) _Tstop = double( _Tstart );

			if ((!_Tstep.has_hard_value() )|| ((double)_Tstep == 0)) {
				trace1("set Tstep ", _Tstop);
				_Tstep=_Tstop;
			} else {
				trace1("set Tstep ", _Tstep);
			}
			if (_Tstart!=0) {
				_cont_tt = true;
				_cont = true;
			}

		} else {
			unreachable();
			assert (!arg1.has_hard_value());   // for now...

			double oldrange = _Tstop - _Tstart;
			_Tstart = _sim->_last_Time;
			_Tstop  = _sim->_last_Time + oldrange;
		}
	}else{ /* no args */
		// std::cerr << "setup ttt -- no args\n";
		double oldrange = _Tstop - _Tstart;
		_Tstart = _sim->_last_Time;
		_Tstop  = _sim->_last_Time + oldrange;

		if(_sim->_last_Time==0){
			trace1("TTT::setup_tw no args at beginning", _cont_tt);
			_Tstop=0;
			_tstep=0;
			_tstop=0;

		}
	}

	if (_tstep>_tstop){
		error(bWARNING, "_tstep > _tstop. really?\n");
	}

	options(Cmd);

	_tstart.e_val(0., _scope);
	_tstop.e_val(NOT_INPUT, _scope);
	_tstep.e_val(NOT_INPUT, _scope);
	_Tstop.e_val(NOT_INPUT, _scope);
	_Tstep.e_val(NOT_INPUT, _scope);

	//Time1 = 
	_sim->_Time0 = _Tstart;
	_Time1 = _Tstart;


	if (!_tstep.has_good_value()) {
		throw Exception("transient: Time step is required");
	}else if (_tstep==0. && _tstop ) {
		untested();
		throw Exception("Time step == 0 while tend");
	}else{
	}

//		if (_dtmax_in.has_hard_value()) {
//			_dtmax = _dtmax_in;
//		}else if (_skip_in.has_hard_value()) {
//			_dtmax = _tstep / double(_skip_in);
//		}else{
//			_dtmax = std::min(_dtmax_in, _tstep);
//		}

	_dTmin= _tstop * .5;
	_sim->_dTmin= _tstop; // FIXME: don't use

	if (_dTmin_in.has_hard_value()) {
		_dTmin = _dTmin_in;
	}else if (_dtratio_in.has_hard_value()) {
		_dTmin = _dTmax / _dTratio_in;
	}else{
		// use larger of soft values
		// _dTmin = std::max(double(_dTmin_in), _dTmax/_dTratio_in);
		// _dTmin=0.5; // HACK
	}


	if  ( _Tstart < _sim->_last_Time  ||  _sim->_last_Time <= 0.) {
		//    _out << "* last_Time " << _sim->_last_Time << "\n";
		trace3("TTT::setup_tw no cont ", _Tstart, _sim->_last_Time, _cont_tt );
		//_cont_tt = false;
		_Time1 = _sim->_Time0 = 0.;
	}else{
		_cont_tt = true;
		_Time1 = _sim->_Time0 = _sim->_last_Time;
	}

	// TRANSIENT setup

	_tstart.e_val(0., _scope);
	_tstop.e_val(NOT_INPUT, _scope);
	_tstep.e_val(NOT_INPUT, _scope);

	if  (_cont_tt){
		_cont = false;

	}else{
		_cont = false;
	}

	_time1 = _sim->_time0 = 0.;
	//{}else{
	//  untested();
	// _cont = true;
	//  time1 = _sim->_time0 = _sim->_last_time;
	// }
	_sim->_freq = ((_tstop > _tstart) ? (1 / (_tstop - _tstart)) : (0.));

	if (!_tstep.has_good_value()) {
		throw Exception("transient: time step is required");
	}else if (_tstep==0. && _tstop ) { itested();
		throw Exception("time step == 0 and tstop");
	}else{
	}

	if (_dtmax_in.has_hard_value()) {
		_dtmax = _dtmax_in;
	}else if (_skip_in.has_hard_value()) {
		_dtmax = _tstep / double(_skip_in);
	}else{
		_dtmax = std::min(_dtmax_in, _tstep);
	}

	if (_dtmin_in.has_hard_value()) {
		_sim->_dtmin = _dtmin_in;
	}else if (_dtratio_in.has_hard_value()) {
		_sim->_dtmin = _dtmax / _dtratio_in;
	}else{
		// use larger of soft values
		_sim->_dtmin = std::max(double(_dtmin_in), _dtmax/_dtratio_in);
	}

	steps_total_out_ = (uint_t) (1 + ceil( ( (_tstop - _tstart ) / _tstep ) ));
	trace7( "TTT::setup_tw done ",  steps_total_out_ , _tstep , _Tstep ,_stepmode, _cont, _cont_tt, _Tstop );
	trace7( "TTT::setup_tw done ",  steps_total_out_ , (double)_tstep , _tstop ,_tstart, _cont, _cont_tt, (double)_Tstop );
	allocate();
}
	/*--------------------------------------------------------------------------*/
}
