
#ifndef DDC_H__
#define DDC_H__
class DDC_BASE : public SIM { // public DCOP??
public:
  void	finish();
protected:
  void	fix_args(int);
  void	options(CS&, int);
  bool _old_solver;
private:
  void	sweep();
  void	sweep_recursive(int);
  void	first(int);
  bool	next(int);
  void do_tran_step();
  void undo_time_step();
  explicit DDC_BASE(const DDC_BASE&): SIM() {unreachable(); incomplete();}
protected:
  explicit DDC_BASE();
  ~DDC_BASE() {}

  void block_solver(); // solve block matrix ( C 0; G C )
  void old_solver(); // solve using G^-1
  
protected:
  enum {DCNEST = 4};
  int _n_sweeps;
  PARAMETER<double> _start[DCNEST];
  PARAMETER<double> _stop[DCNEST];
  PARAMETER<double> _step_in[DCNEST];
  PARAMETER<double> _para[DCNEST];
  double _step[DCNEST];
  bool _linswp[DCNEST];
  double* _sweepval[DCNEST];	/* pointer to thing to sweep, dc command */
//  typedef void (*p)(double);
  ELEMENT* _pushel[DCNEST]; /* pointer to thing to sweep, dc command */
  ELEMENT* _zap[DCNEST]; /* to branch to zap, for re-expand */
  std::string _para_name[DCNEST];
  PARAMETER<double>* _param[DCNEST];

  std::vector<STORAGE*> _uic_caplist;
  void set_uic_caps_constant(bool x=true);

  CARDSTASH _stash[DCNEST];	/* store std values of elements being swept */
  bool _loop[DCNEST];		/* flag: do it again backwards */
  bool _reverse_in[DCNEST];	/* flag: sweep backwards, input */
  bool _reverse[DCNEST];	/* flag: sweep backwards, working */
  bool _cont;			/* flag: continue from previous run */
  TRACE _trace;			/* enum: show extended diagnostics */
  enum {ONE_PT, LIN_STEP, LIN_PTS, TIMES, OCTAVE, DECADE} _stepmode[DCNEST];
  static double temp_c_in;	/* ambient temperature, input and sweep variable */
  bool _do_tran_step;
  double _tran_step;
  bool _dump_matrix;
  double* U;
  double* CU;
  double* CUTCU;
  double* A;
  double* y;
  void ac_snapshot();

private: // cleanup later
  double* _Gu;
  double* _dv;
};
#endif
