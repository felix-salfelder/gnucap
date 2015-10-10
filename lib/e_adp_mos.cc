
//#include "d_mos.h"
#include "d_mos8.h"
#include "e_adp_mos.h"
#include "u_nodemap.h" // fixme?
#include "d_bti.h"

// aging helpers.
// alpha version
/*--------------------------------------------------------------------------*/
// double ADP_BUILT_IN_MOS::wdT() const{ untested();
// 	return ids_stress->wdT();
// }
/*--------------------------------------------------------------------------*/
void ADP_BUILT_IN_MOS::tt_accept()
{ untested();
	trace0("ADP_BUILT_IN_MOS::tt_accept");
	//FIXME: move c to ADP_CARD. merge ADP_card with DEV?
	// const DEV_BUILT_IN_MOS* c = (const DEV_BUILT_IN_MOS*) (bti_stress->c());
	//SIM_DATA* sim = c->_sim;
	//  std::cerr << "ADP_BUILT_IN_MOS::tt_accept " << c->long_label() << "\n";
	//  std::cerr << "ADP_BUILT_IN_MOS::tt_accept time " << sim->_Time0 << "\n";
	//  std::cerr << "ADP_BUILT_IN_MOS::tt_accept stress " << bti_stress->get() << "\n";
	const COMMON_COMPONENT* cc = prechecked_cast<const COMMON_COMPONENT*>(_c->common());
	const MODEL_BUILT_IN_MOS_BASE* m = prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(cc->model());
	assert(m);
	assert(cc);

	if (m->use_bti()) { untested();
		// btistress_taken  = bti_stress->get();
	} else { untested();
	}
}
/*--------------------------------------------------------------------------*/
double ADP_BUILT_IN_MOS::tt_probe_num(const std::string& x)const
{ untested();
	unreachable(); //?
	if( Umatch("bti ", x) ){ untested();
		return bti_stress->get();
	} else if( Umatch("dvth_bti ", x) ) { untested();
		return bti_stress->get();
	} else if( Umatch("dvth_bti ", x) ) { untested();
		return vthdelta_bti;
	}else{ untested();
		return 999;    //    return ADP_BUILT_IN_MOS::tr_probe_num(x); diode?
	}
}
/*--------------------------------------------------------------------------*/
void ADP_BUILT_IN_MOS::apply(const COMPONENT* x)
{
	const DEV_BUILT_IN_MOS* d = prechecked_cast<const DEV_BUILT_IN_MOS*>(x);
	assert(d);

	delta_vth = 0; // BUG?

	if (d->bti_device()) {
		const DEV_BUILT_IN_BTI* bti = prechecked_cast<const DEV_BUILT_IN_BTI*>(d->bti_device());
		assert(bti);
		delta_vth += bti->dvth();
		assert(is_number(delta_vth));
	} else {
	}
	trace2("ADP_BUILT_IN_MOS::apply done", long_label(), delta_vth);
}
/*--------------------------------------------------------------------------*/
// obsolete. cleanup later.
void ADP_BUILT_IN_MOS::tr_accept()
{
	const DEV_BUILT_IN_MOS* d = prechecked_cast<const DEV_BUILT_IN_MOS*>(owner());
	assert(d); USE(d);
//	const COMMON_BUILT_IN_MOS* c = asserted_cast<const COMMON_BUILT_IN_MOS*>(d->common());
//	const MODEL_BUILT_IN_MOS_BASE* m = asserted_cast<const MODEL_BUILT_IN_MOS_BASE*>(c->model());
#if 0
	if (m->use_bti()) { untested();
		assert(bti_stress);
		bti_stress->tt() = 0; // not in use (yet?)
		bti_stress->set_tr(0); // not in use (yet?)
	}
#endif

	if(_sim->_time0 <= _tr_last_acc) {
		if(_sim->_time0){
			cerr << "acc twice\n";
			error(bTRACE, "%s, acc twice %f, %f\n", long_label().c_str(), _sim->_time0, _tr_last_acc);
			trace2("tr ADP accepting twice %s time0 %E\n", d->long_label().c_str(), _sim->_time0);
		}else{
		}
	}

	_tr_last_acc = _sim->_time0;
}
/*--------------------------------------------------------------------------*/
double ADP_BUILT_IN_MOS::tr_probe_num(const std::string& x)const
{ unreachable();
	if( Umatch("bti ", x) ){ untested();
		if(bti_stress) return bti_stress->tr_get();
		return 9999;
	} else if( Umatch("dvth_bti ", x) ) { untested();
		return vthdelta_bti;
	}else{ untested();
		untested();
		return 888;   //    return ADP_BUILT_IN_MOS::tr_probe_num(x); diode?
	}

}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
ADP_BUILT_IN_MOS::ADP_BUILT_IN_MOS( COMPONENT* c, const std::string n) :
	ADP_BUILT_IN_DIODE(c,n),
	bti_stress(0),
  	_tr_last_acc(-inf)
{
	_n = _nodes;
	init(c);
}
/*--------------------------------------------------------------------------*/
void ADP_BUILT_IN_MOS::init(const COMPONENT* c)
{
	trace1("ADP_BUILT_IN_MOS::init", c->long_label());
	const COMMON_COMPONENT* cc = prechecked_cast<const COMMON_COMPONENT*>(c->common());
	const MODEL_BUILT_IN_MOS_BASE* m = prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(cc->model());
	assert(m); USE(m);
	assert(cc);

	vthscale_bti = 1;
	vthdelta_bti = 0;
	delta_vth = 0;
}
/*--------------------------------------------------------------------------*/
void ADP_BUILT_IN_MOS::expand()
{ untested();
	unreachable();
	return;
	assert(0);
	DEV_BUILT_IN_MOS* d = asserted_cast<DEV_BUILT_IN_MOS*>(owner());
	const COMMON_COMPONENT* cc = prechecked_cast<const COMMON_COMPONENT*>(d->common());
	const MODEL_BUILT_IN_MOS_BASE* m = prechecked_cast<const MODEL_BUILT_IN_MOS_BASE*>(cc->model());
	assert(m);
	assert(cc);

	if (!subckt()) { untested();
		new_subckt();
	}else{ untested();
	}
	trace2("ADP_BUILT_IN_MOS::expand", _sim->is_first_expand(), m->use_bti());

	if (_sim->is_first_expand()) { untested();
		if( m->use_bti()){ untested();

			if(bti_stress){ untested();
				error(bDANGER, "%s, btistress: %f\n", long_label().c_str(), bti_stress);
				bti_stress = 0;
			}

			// precalc_first();
			// precalc_last();
			trace4("ADP_BUILT_IN_MOS::expand, first", long_label(), hp(this), _n[n_bti].n_(), _n[n_bti].is_adp());
			assert(_n[n_bti].is_adp());
			assert(!(_n[n_bti].n_())); // n_ is electrical
			if (!(_n[n_bti].a_())){ untested();
				trace0("ADP_BUILT_IN_MOS::expand, no 3rd external node connected");
				_n[n_bti].new_model_adp_node("bti", d);
				trace1("ADP_BUILT_IN_MOS::expand have new adpnode", hp(_n[n_bti].a_()));
			} else { untested();
				trace1("ADP_BUILT_IN_MOS::expand, 3rd node present, external", hp(_n[n_bti].a_()));
				unreachable();
			}
		}
	}

	if (m->use_bti()) { untested();
		bti_stress = _n[n_bti].a_(); assert(bti_stress);
	} else { untested();
	}
}
/*--------------------------------------------------------------------------*/
void ADP_BUILT_IN_MOS::tr_stress_last()
{ untested();
	trace1("ADP_BUILT_IN_MOS::tr_stress_last", long_label());
	// HIER
	const DEV_BUILT_IN_MOS* d = asserted_cast<const DEV_BUILT_IN_MOS*>(owner());
	const COMMON_BUILT_IN_MOS* c = asserted_cast<const COMMON_BUILT_IN_MOS*>(d->common());
	const MODEL_BUILT_IN_MOS_BASE* m = asserted_cast<const MODEL_BUILT_IN_MOS_BASE*>(c->model());
	assert(m);

	//delta_vth = 0;

	if (m->use_bti()){ untested();
		const DEV_BUILT_IN_BTI* bti = prechecked_cast<const DEV_BUILT_IN_BTI*>(d->BTI()); // hack, in a way
		if (!bti){ untested();
			error(bDANGER,"no bti here? %i %i %f\n", _sim->tt_iteration_number(), _sim->iteration_number(),
					_sim->_Time0);
			assert(bti);
		}

		bti_stress->tt() = bti->dvth();
		bti_stress->set_tr(0); // not in use (yet?)
		incomplete(); // call apply or do similarly
	}
}

/*--------------------------------------------------------------------------*/
// MOS8
/*--------------------------------------------------------------------------*/
double ADP_BUILT_IN_MOS8::tt_probe_num(const std::string& x)const
{ untested();
	trace1("ADP_BUILT_IN_MOS8::tt_probe_num", x);
	if (Umatch("hci|dvth_hci ", x)) { untested();
		assert(fabs(vthdelta_hci)<10);
		return vthdelta_hci;
	} else { untested();
		return ADP_BUILT_IN_MOS::tt_probe_num(x);
	}
	return -1;
}
/*--------------------------------------------------------------------------*/
void ADP_BUILT_IN_MOS8::tr_accept() // tr_stress() ...
{
	const DEV_BUILT_IN_MOS* d = prechecked_cast<const DEV_BUILT_IN_MOS*>(owner());
	assert(d);
	const COMMON_BUILT_IN_MOS* c = asserted_cast<const COMMON_BUILT_IN_MOS*>(d->common());
	const MODEL_BUILT_IN_MOS8* m = asserted_cast<const MODEL_BUILT_IN_MOS8*>(c->model());
	//ADP_BUILT_IN_MOS8* a = (ADP_BUILT_IN_MOS8*) this;
	ADP_BUILT_IN_MOS::tr_accept();
	const SDP_BUILT_IN_MOS8* s = prechecked_cast<const SDP_BUILT_IN_MOS8*>(c->sdp());

	if (0 && m->use_hci()) { // must not accept twice!
		double exponent = 3; // which is m in [4]
		double hcis = 0;
		double Wg = 0.8;
		double Ids = fabs(d->ids); // fabs?

		//  std::cerr << "DEV_BUILT_IN_MOS8::h0 of "<<  d->short_label() << " " <<   m->h0 << "\n";
		double H = m->h0;
		double W = s->w_eff;
		double Hg = m->h0;
		//  double m0 = m->m0;

		if( Ids < 1e-40) { untested();
			trace1("MODEL_BUILT_IN_MOS8::tr_stress ids too small: ", d->ids );
		} else { untested();
			hp_float_t Isub;
			if (d->reversed) { untested();
				Isub = d->isb;
			} else { untested();
				Isub = d->idb;
			}

			assert(Isub >= 0);
			//  assert(d->ids >= 0); isnich
			//  assert( dt >= 0 )

			switch(m->polarity){ untested();
				case dunno:
					unreachable();
				case pN:
					hcis = Ids * pow( Isub / Ids, exponent)/H/W;
					assert(is_number(hcis));
					trace1( "ADP_BUILT_IN_MOS8::tr_accept nmos", hcis );
					break;
				case pP:
					double mg = 3.0;
					double ig = d->probe_num("ig"); // BUG
					hcis = ( Wg/Hg * pow( fabs(ig)/W, mg )
							+ (1-Wg)*Ids/H/W * pow(Isub/fabs(Ids), exponent));
					trace8( "ADP_BUILT_IN_MOS8::tr_accept pmos", d->long_label(),
							d->isb, d->idb, _sim->_time0, Isub, Ids, H, hcis);
					assert(is_number(hcis));
					break;
			}
			if (hcis > 1e-10)
			{ untested();
			}
		}

		// a->_raw_hci_node->add_tr(hcis * _sim->_dt0); doesnt work, a->tr must be effective mean
		//
		// effectively add up on tt, and take difference in the end (numerically unstable...)
		long double hcis1 = _hci_stress;
		_hci_stress = hcis;
		//cout << _hci_tr << "+=" << hcis * _sim->_dt0 << " -> ";
		//long double hold = _hci_tr;
// 	   cout << _sim->_dt0 << " ";
// 		cout << _sim->_time0 << "  --  " << hcis << "\n";
		if (_sim->_dt0) { untested();
			_hci_tr += (hcis+hcis1) * _sim->_dt0;
		}else{untested();
		}
		//cout << _hci_tr - hold << "\n";

// debug
//		vthdelta_hci = pow(_raw_hci_node->tt()+_hci_tr,0.3);
		//delta_vth=vthdelta_hci;
	} // end hci
	q_eval();
}
/*--------------------------------------------------------------------------*/
double ADP_BUILT_IN_MOS8::tr_probe_num(const std::string&)const
{ untested();
	return 1e99;
#if 0
	double ret=771;
	if( Umatch(x, "hci_raw ") ){ untested();
		if(_raw_hci_node) return _raw_hci_node->tt()+_hci_tr;
		if(_raw_hci_node) return _hci_tr;
	}else if( Umatch(x, "hcistress ") ){ untested();
		if(_raw_hci_node) return _hci_stress;
	}else	if( Umatch(x, "hci_tt ") ){ untested();
		if(_raw_hci_node) return _raw_hci_node->tt();
	}else	if( Umatch(x, "hci ") ){ untested();
		// assert(fabs(vthdelta_hci)<10);
		return vthdelta_hci;
		//if(!_raw_hci_node) return -1;
		//ret= _raw_hci_node->tr();
	} else { untested();
		trace0("ADP_BUILT_IN_MOS8::tr_probe_num, fallback");
		ret= ADP_BUILT_IN_MOS::tr_probe_num(x);
	}

	// maybe too small to output?
	//  std::cerr << "ADP_BUILT_IN_MOS8::tr_probe_num " << x << ": " << ret << "\n";

	return ret;

#endif
}
/*--------------------------------------------------------------------------*/
void ADP_BUILT_IN_MOS8::init(const COMPONENT*)
{
#if 0
  const COMMON_BUILT_IN_MOS* c = asserted_cast<const COMMON_BUILT_IN_MOS*>(d->common());
  const MODEL_BUILT_IN_MOS8* m = asserted_cast<const MODEL_BUILT_IN_MOS8*>(c->model());

  assert(c);
  trace1( "ADP_BUILT_IN_MOS8::init", m->use_hci() );
  
  // constructor does that. (init is intentionally non-virtual)
  // ADP_BUILT_IN_MOS::init(d);

  if (m->use_hci()) { untested();
    assert(!_raw_hci_node);
	 _nodes[n_bti].set_adp();
//    _raw_hci_node = scope()->nodes()->new_adp_node( "raw_hci", d );
//    _nodes[0] = _raw_hci_node;

    vthdelta_hci = 0;
    vthscale_hci = 1;

  }else
  { untested();
    vthdelta_hci = vthscale_hci = NAN;
    _raw_hci_node = NULL;
  }
#endif
  // _n[n_hci].set_adp();
}
/*--------------------------------------------------------------------------*/
void ADP_BUILT_IN_MOS8::expand()
{ untested();
#if 0
	DEV_BUILT_IN_MOS* d = asserted_cast<DEV_BUILT_IN_MOS*>(owner());
	const COMMON_BUILT_IN_MOS* c = asserted_cast<const COMMON_BUILT_IN_MOS*>(d->common());
	const MODEL_BUILT_IN_MOS8* m = asserted_cast<const MODEL_BUILT_IN_MOS8*>(c->model());

	trace1("DEV_BUILT_IN_RCD::expand", long_label());
	// assert(_n[n_hci].is_adp());
	ADP_BUILT_IN_MOS::expand();
	assert(_n);

	if (!subckt()) { untested();
		new_subckt();
	}else{ untested();
	}

	if (_sim->is_first_expand()) { untested();
		// precalc_first();
		// precalc_last();
		if(m->use_hci()){ untested();
			trace4("DEV_BUILT_IN_RCD::expand, first", long_label(), hp(this), _n[n_hci].n_(), _n[n_hci].is_adp());
			assert(_n[n_hci].is_adp());
			assert(!(_n[n_hci].n_())); // n_ is electrical
			if (!(_n[n_hci].a_())){ untested();
				trace0("DEV_BUILT_IN_RCD::expand, no 3rd external node connected");
				_n[n_hci].new_model_adp_node("raw_hci", d);
				trace1("DEV_BUILT_IN_RCD::expand have new adpnode", hp(_n[n_hci].a_()));
			} else { untested();
				trace0("DEV_BUILT_IN_RCD::expand, 3rd node present, external");
				unreachable();
			}
		}
	}

	if(m->use_hci()){ untested();
		_raw_hci_node = _n[n_hci].a_(); assert(_raw_hci_node);
	}
#endif
}
/*--------------------------------------------------------------------------*/
void ADP_BUILT_IN_MOS8::tt_begin()
{
	return; // BUG. has no tt_begin.
#if 0
	const DEV_BUILT_IN_MOS* d = asserted_cast<const DEV_BUILT_IN_MOS*>(owner());
	const COMMON_BUILT_IN_MOS* c = asserted_cast<const COMMON_BUILT_IN_MOS*>(d->common());
	const MODEL_BUILT_IN_MOS8* m = asserted_cast<const MODEL_BUILT_IN_MOS8*>(c->model());
	trace3("ADP_BUILT_IN_MOS8::tt_begin", long_label(), m->use_hci(), _sim->tt_iteration_number());


	if (_sim->_tt_uic) { untested();
		delta_vth = 0;
	}else{ untested();
		delta_vth = 0;
	}
	return;

	if (m->use_hci()){ untested();
		_hci_stress = 0;
		assert(_raw_hci_node);
		if (_sim->_tt_uic) { untested();
			_raw_hci_node->set_tr(0);
		} else { untested();
			_raw_hci_node->set_tr(0);
			_raw_hci_node->set_tt(0);
		}
		_hci_tr = 0;
		_hci_stress = 0;
		vthdelta_hci = 0;
		vthscale_hci = 0;
	} else { untested();
	}
	q_eval();
#endif
}
/*------------------------------------------------------------------*/
void ADP_BUILT_IN_MOS8::tt_advance()
{
#if 0
	const DEV_BUILT_IN_MOS* d = asserted_cast<const DEV_BUILT_IN_MOS*>(owner());
	const COMMON_BUILT_IN_MOS* c = asserted_cast<const COMMON_BUILT_IN_MOS*>(d->common());
	const MODEL_BUILT_IN_MOS8* m = asserted_cast<const MODEL_BUILT_IN_MOS8*>(c->model());
	if (m->use_hci()) { untested();
		assert(_raw_hci_node);
		_hci_tr = 0;
		_raw_hci_node->set_tr(0);
	}
#endif
}
/*--------------------------------------------------------------------------*/
void ADP_BUILT_IN_MOS8::tr_stress_last()
{ untested();
	ADP_BUILT_IN_MOS::tr_stress_last();

	assert(_sim->last_time());

#if 0
	if(m->use_hci()) { untested();
		_raw_hci_node->set_tr(_hci_tr/_sim->last_time());
		q_tt_accept();
		_raw_hci_node->set_tr_noise(0); // incomplete;
	} else if(0) { // old hci

		trace2("ADP_BUILT_IN_MOS8::tr_stress_last hci", _sim->last_time(), _hci_tr/_sim->last_time() );

		// ??!
		_raw_hci_node->tt() += _hci_tr; // tt at last_Time.

		if (!_sim->last_time()) { untested();
			assert(_sim->phase() == p_PD);
			_raw_hci_node->set_tr(0);
		} else { untested();
			_raw_hci_node->set_tr(double(_hci_tr/_sim->last_time()));  // tr= d(tt)/dt
		}

		assert(is_number(_raw_hci_node->tr()));
		trace2("ADP_BUILT_IN_MOS8::tr_stress_last", _raw_hci_node->tt(), _raw_hci_node->tr());
		_hci_tr = 0; // tt_advance!
		// obsolete?
		vthdelta_hci = pow(_raw_hci_node->tt(), .3); // this is wrong :/
		delta_vth += vthdelta_hci;
	} else { untested();
	}
#endif
}
/*------------------------------------------------------------------*/
void ADP_BUILT_IN_MOS8::apply(const COMPONENT* x)
{
	ADP_BUILT_IN_MOS::apply(x);
	const DEV_BUILT_IN_MOS* d = prechecked_cast<const DEV_BUILT_IN_MOS*>(x);
	assert(d);

	if (const COMPONENT* hci = d->hci_device()) {
		trace3("ADP_BUILT_IN_MOS8::apply", x->long_label(), hci->value(), delta_vth);
		assert(is_number(hci->value()));
		delta_vth += hci->value();
	} else { untested();
	}
	trace2("ADP_BUILT_IN_MOS8::apply done", long_label(), delta_vth);
// 	q_eval(); // ?
}
/*------------------------------------------------------------------*/
// before accept: tt() = extrapolated value before transient run
//                tt1() accepted value at last_Time
// during tran:  tt() whatever.
// after accept: tt() = accepted value after,
void ADP_BUILT_IN_MOS8::tt_accept()
{ untested();
	unreachable();
#if 0
	DEV_BUILT_IN_MOS* d = asserted_cast<DEV_BUILT_IN_MOS*>(owner());
	const COMMON_BUILT_IN_MOS* c = asserted_cast<const COMMON_BUILT_IN_MOS*>(d->common());
	const MODEL_BUILT_IN_MOS8* m = asserted_cast<const MODEL_BUILT_IN_MOS8*>(c->model());
	ADP_BUILT_IN_MOS::tt_accept();
	if (0 && m->use_hci()) { untested();
		_raw_hci_node->tt() += _hci_tr;
		vthdelta_hci = pow(_raw_hci_node->tt(), .3);
	}
#endif
}
/*------------------------------------------------------------------*/
// vim works with any ts=sw
