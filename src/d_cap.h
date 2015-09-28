
// this is a workaround, to have a cast to DEV_CAPACITANCE in s_sock.cc
//
#include "e_storag.h"
namespace SOME_CAP_HACK 
{
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
class DEV_CAPACITANCE : public STORAGE {
protected:
  explicit DEV_CAPACITANCE(const DEV_CAPACITANCE& p) :STORAGE(p) {}
public:
  explicit DEV_CAPACITANCE()	:STORAGE() {}
 // void      keep_ic();
protected: // override virtual
  char	   id_letter()const	{return 'C';}
  std::string value_name()const {return "c";}
  std::string dev_type()const	{return "capacitor";}
  uint_t	   max_nodes()const	{return 2;}
  uint_t	   min_nodes()const	{return 2;}
  uint_t	   matrix_nodes()const	{return 2;}
  uint_t	   net_nodes()const	{return 2;}
  bool	   has_iv_probe()const  {return true;}
  bool	   use_obsolete_callback_parse()const {return true;}
  CARD*	   clone()const		{return new DEV_CAPACITANCE(*this);}
  void	   tr_iwant_matrix()	{tr_iwant_matrix_passive();}
  bool	   do_tr();
  void	   tr_accept(); // uic. possibly a hack
  void	   tr_load()		{tr_load_passive();}
  void	   tr_unload()		{tr_unload_passive();}
public:
  double    tr_involts()const	{ return tr_outvolts(); }
#ifdef TTCAP_HACK
public: //tt
  void tt_begin();
  TIME_PAIR tt_review();
  void tt_accept();
  void tt_advance();
  void tt_regress();
  void tr_stress_last();
  void tr_restore();
  void do_tt();
private: // another hack -> storag
  double _ttstate[4];
  double _ttstate_;
  double _dttstate[4];
  double _dttstate_;
  double _pred;
  double _corr;
  double _preddt;
  double _ttstep;
  double _ttstep0;
  double _ttstep1;
  double _ttstep2;
  double _ttstep3;
  double _ttfuture;
  int _ttsteporder;
  int _new_ttsteporder;
  int _ttmethod;
  double _dT1;
  double _dT2;

  double _vt[2];
  double _vy[2];
  double _b[2];
  double _chisq;
  double _dv;
  double gsl_fit_3();
  double gsl_fit_4();

  double _extime;
#endif
protected:
  hp_float_t   tr_involts_limited()const {return tr_outvolts_limited();}
  double   tr_probe_num(const std::string&)const;
  void	   ac_iwant_matrix()	{ac_iwant_matrix_passive();}
  void	   ac_begin()		{_ev = _y[0].f1;}
  void	   do_ac();
  void	   ac_load()		{ac_load_passive();}
  COMPLEX  ac_involts()const	{itested();return ac_outvolts();}

  std::string port_name(uint_t i)const {
    assert(i != INVALID_NODE);
    assert(i < 2);
    static std::string names[] = {"p", "n"};
    return names[i];
  }
private: // hack
  bool has_ic()const;

};

// inline ostream& operator<<( ostream& o, const DEV_CAPACITANCE &c){
// 	o << "cap " << c.long_label() << " " << c.value();
// 	return o;
// }


} // namespace
