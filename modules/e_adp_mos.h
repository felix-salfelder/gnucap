#ifndef ADPMOS_H
#define ADPMOS_H

class ADP_BUILT_IN_MOS :public ADP_BUILT_IN_DIODE{
public:
  explicit ADP_BUILT_IN_MOS( COMPONENT* c, const std::string n);
  virtual ADP_CARD* clone()const{ return new ADP_BUILT_IN_MOS(*this); }
protected:
  ADP_BUILT_IN_MOS(const ADP_BUILT_IN_MOS& a): ADP_BUILT_IN_DIODE(a),
	vthscale_bti(a.vthscale_bti),
	vthdelta_bti(a.vthdelta_bti),
	delta_vth(a.delta_vth)
	{
		_n = _nodes;
	};

  void init(const COMPONENT*);
public:
//  ADP_NODE* ids_stress; // obsolete?
//  ADP_NODE* igd_stress;
public:

  ADP_NODE* bti_stress; // FIXME, BTI_?
  double tr_probe_num(const std::string& x)const;
  double tt_probe_num(const std::string& x)const;

private:
//   double btistress_taken; // unused?
//   double bti_eff_voltage; // unused?
  double _tr_last_acc; // unused?

protected:
  enum {n_bti};
  node_t _nodes[1];
public:
  unsigned max_nodes()const{return 1;}
  unsigned min_nodes()const{return 1;}

public:
  void expand();
  virtual void tt_accept();
  virtual void tr_stress_last();
  virtual void apply(const COMPONENT*);
  double vthscale_bti ; //  exp ( 10000. * a->hci_stress->get() / c->w_in );
  double vthdelta_bti ;
//  virtual double wdT()const;
  double delta_vth ;

  virtual void tr_accept();
};
/*--------------------------------------------------------------------------*/
#endif
