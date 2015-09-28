
#ifndef S_TT_H
#define S_TT_H

#include "s_tr.h"
#include "u_probe.h"
#include "u_prblst.h"

namespace TT {

enum stepmode_t {tts_UNKNOWN=0, tts_LIN, tts_MUL};
enum agemode_t {amNONE=0, amONCE, amBYPASS, amALWAYS};

class TTT : public TRANSIENT {
  public:
    void	do_it(CS&, CARD_LIST*);
    std::string status()const;

    explicit TTT():
      TRANSIENT(),
      _trace(tNONE),
      _alarm(aREPORT),
      _evt(true),
      _Tstart(0.), // ?? unused?
      _Tstop(0.),
      _Tstep(0.),
      _timesteps(0),
		_new(false),
		_cont_tt(false),
      _fdata_tt(NULL),
      _tt_store(0),
		_no_act(false)
  { }
    ~TTT() {}
    void tt_advance_time();
  public:
    int step_cause()const;
    void	advance_Time();

  private:
    static OMSTREAM mstdout;
    OMSTREAM _out;
    explicit TTT(const TTT&): 
      TRANSIENT(),
      _trace(tNONE),
      _Tstart(0.),
      _Tstop(0.),
      _Tstep(0.),
      _timesteps(0),
      _fdata_tt(NULL)
    {incomplete();} // reachable?
    void	setup(CS&);
    void	setup_tw(CS&);
    void	allocate();
    void  rescale_behaviour();
    void	finish();
    bool	next();
    void	out_tt();
    void	outdata_tt(double);
    void	outdata_b4(double);
    bool	review();
    void tt_begin(); // from BASE
    void tt_cont();
    void	options(CS&);
    void	sweep_tr();
	 void do_initial_dc();
    void	sweep();
    double get_new_dT();
    void	accept();
    void	print_head_tr();
    void	print_foot_tr();

	 // FIXME: private could inhibit call from within TRANSIENT....
    void	head(double,double,const std::string&);

    void	head_tt(double,double,const std::string&);
    void	set_step_cause(STEP_CAUSE);
    void	first();
    void after_interruption_prep();
    void first_after_interruption();
    void	fohead(const PROBE&);
    void	store_results(double); // override virtual
    void	store_results_tt(double); 
    void	print_stored_results_tt(double); 
    void	power_down(double ); 
  private:
	 bool conchk() const;
    TRACE _trace;		// enum: show extended diagnostics
    ALARM _alarm;
    bool _power_down;
	 bool _evt;
    PARAMETER<double> _Tstart;	// unused?
    PARAMETER<double> _Tstop;	/* user stop Time */
    PARAMETER<double> _Tstep;	/* user Tstep */
    PARAMETER<double> _dTmin_in;
    PARAMETER<double> _dTratio_in;
    int    _timesteps;	/* number of time steps in tran analysis, incl 0 */
    int    _Timesteps;	/* number of time steps in tt analysis, incl 0 */
    int print_tr_probe_number;
    double _Time1;
    uint_t steps_total_tt;
	 double behaviour_time();
	 bool _new;	 // reset to fresh devices. (and do initial dc)
	 bool _cont_tt; // continue with adps
    double time_by_voltages();
    double _time_by_adp;
    double _dT_by_adp;
    double _dTmin;
    double _dTmax;
    void   outdata(double);
    void   print_results(double); // stupid?
    void   print_results_tr(double);
    void   print_results_tt(double);
    double _Time_by_user_request;
	 double _Time_by_error_estimate;
	 stepmode_t _stepmode;
	 agemode_t _agemode;

	 unsigned _printed_steps;
    bool _accepted_tt;
	 bool _quiet;
    COMPLEX** _fdata_tt;	/* storage to allow postprocessing */
    double*   _tt_store;
    bool _no_act;
	 PROBELIST oldstore; //save tr_store (which is abused for caching/measurements)
	 void probeexpand();
  private: // sim override (hack)
    void tt_alarm(ALARM a=aREPORT, OMSTREAM*o=NULL);
  public: //status
	 unsigned steps_accepted()const{return _sim->_tt_accepted;}
	 unsigned steps_rejected()const{return _sim->_tt_rejects;}
	 unsigned steps_total()const{return steps_total_tt;}
  protected:
	 WAVE** _wavep_tt; // FIXME: use SIM::_wavep
}; // TTT : TRANSIENT
/*--------------------------------------------------------------------------*/
struct Exception_Alarm :public Exception{
  Exception_Alarm(const std::string& Message)
    :Exception(Message) {
  }
};
/*--------------------------------------------------------------------------*/
template<class T>
T& operator<<(T& o, const stepmode_t& x){
	switch(x){
		case tts_UNKNOWN: return (o << "unknown");
		case tts_LIN: return (o << "lin");
		case tts_MUL: return (o << "mul");
	}
	return o;
}

}
#endif

// vim:ts=3:sw=3:
