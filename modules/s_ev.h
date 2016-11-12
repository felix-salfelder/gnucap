#ifndef EV_H__
#define EV_H__

#include "s__.h"
#include "u_parameter.h"
#include "e_elemnt.h"
#include "e_storag.h"
#include "u_status.h"
#include "u_prblst.h"
#include "u_cardst.h"
#include "io_matrix.h"
#include "e_aux.h" // set_sens_port...

#include <gsl/gsl_multifit.h>

const gsl_complex one = gsl_complex_rect(1.,0.);
const gsl_complex zero = gsl_complex_rect(0.,0.);

class EV_BASE : public SIM { // public DDC_BASE?
	public:
		void	finish();
		string status()const;
		TRACE trace() const{return _trace;}
	protected:
		void	fix_args(int);
		virtual void options(CS&, int);
	private:
		void	sweep();
	protected:
		void	first(int);
		bool	next(int);
		virtual void sweep_recursive(int){};
	private:
		void do_tran_step();
		void undo_time_step();
		explicit EV_BASE(const EV_BASE&): SIM() {unreachable(); incomplete();}
	protected:
		explicit EV_BASE();
		~EV_BASE() {}

	protected:
		enum {DCNEST = 5};
		int _n_sweeps; // all sweeps
		unsigned _n_inputs; // just input devices
		unsigned _inputno[DCNEST]; // map sweeps to inputs
		unsigned _n_states;
		PARAMETER<double> _start[DCNEST];
		PARAMETER<double> _stop[DCNEST];
		PARAMETER<double> _step_in[DCNEST];
		PARAMETER<double> _para[DCNEST];
		double _step[DCNEST];
		bool _linswp[DCNEST];
		double* _sweepval[DCNEST];	/* pointer to thing to sweep, dc command */
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
		bool _uic;
		bool _quantize_states;
		double _slew_rate;
		double* U;
		double* CU;
		double* CUTCU;
		double* A;
		double* y;
		void ac_snapshot();

		SSP_STATE& ev_solver(SSP_CHART*& c);
		unsigned sens_inf_vectors(gsl_matrix_complex_view EV, gsl_vector_view grid, gsl_matrix_complex* dc_sens);
		int process_tail(gsl_matrix_complex_view& relevant_EV);
		void output_sens_grid(gsl_matrix_complex_view&, gsl_vector*, const gsl_vector_complex* lambda);
		void do_ddc(gsl_vector_complex_view z0dot_);
		void compute_zdot(gsl_vector_complex_const_view z0_left, gsl_vector_complex_view z0dot_);
		void compute_canonical_op(gsl_vector_complex_view op, gsl_matrix_complex_const_view dcsens);
		void quantize(gsl_vector_complex_view z0, gsl_vector_complex_const_view lambda, gsl_vector_const_view tol, index_t* index);
		//void compute_op(OPT::ITL);
		void get_op(OPT::ITL); // compute_op.
		void get_op_();

		struct output_t{
			string label;
			CKT_BASE* brh[2];
		};
		std::vector<output_t> _output;
		std::vector<output_t>::iterator _output_iter;

	protected: //sim_data?
		static void set_sens(std::vector<output_t>::iterator _output_iter);
		static void clear_sens();

	protected: // cleanup later
		double* _Gu;
		double* _dv;


		gsl_vector_complex *ztmp_, *sens_;
		gsl_vector *beta_, *betabar_;
		gsl_vector *rhs_;
		gsl_matrix *Q_, *Z_; // probably unneeded.
		gsl_matrix *G_, *C_;
		gsl_matrix_complex *EV_;
		gsl_matrix_complex *EVLU_;
		gsl_permutation *LU_perm;

		gsl_vector_complex *z0_, *z0dot_; // EV * z0 = v0
		gsl_vector_complex *alpha_, *lambda_; //  *mu_;
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// sim_data?
inline void EV_BASE::clear_sens()
{
	for(unsigned i=0; i<_sim->_total_nodes+1;i++){
		_sim->_sens[i] = 0.;
	}
}
/*--------------------------------------------------------------------------*/
inline void EV_BASE::set_sens(std::vector<output_t>::iterator _output_iter)
{
  CKT_NODE* np = dynamic_cast<CKT_NODE*>((*_output_iter).brh[0]);
  assert(np);
  CKT_NODE* nn = dynamic_cast<CKT_NODE*>((*_output_iter).brh[1]);
  if(!nn) nn = &ground_node;
  clear_sens();

  trace2("set_sens_port", node_t(np).m_(), node_t(nn).m_());
  set_sens_port(node_t(np), node_t(nn));
}
/*--------------------------------------------------------------------------*/
inline void EV_BASE::do_ddc(gsl_vector_complex_view z0dot)
{
	trace1("do_ddc G^-1 I", *rhs_);
	size_t fin_dim = z0dot.vector.size;
	gsl_vector_complex_const_view lambda__ = gsl_vector_complex_const_subvector(lambda_, 0, fin_dim);
	gsl_vector_complex_view ztmp = gsl_vector_complex_subvector(ztmp_, 0, fin_dim);
	gsl_vector_view r = gsl_vector_complex_real(ztmp_);
	gsl_vector_view i = gsl_vector_complex_imag(ztmp_);
	gsl_vector_set_zero(&i.vector);
	gsl_vector_memcpy(&r.vector, rhs_);
	trace1("untransformed", *ztmp_);
	gsl_linalg_complex_LU_svx(EVLU_, LU_perm, ztmp_);
	trace1("transformed", *ztmp_);
	//err+= gsl_blas_zgemv( CblasNoTrans, one, EV_, rhs_, zero, ztmp_ );
	//trace1("backtest", *ztmp_);
	gsl_vector_memcpy(rhs_, &r.vector);// for now.
	trace2("newgrid", *rhs_, *lambda_);
	trace1("", z0dot.vector);
	gsl_vector_sub(z0dot, ztmp);
	gsl_vector_mul(z0dot, lambda__);
//	gsl_vector_mul(ztmp, lambda__);
	trace1("limitprime", ztmp);
	trace1("done", z0dot);
}
/*--------------------------------------------------------------------------*/
// compute z0dot from z0 and G^-1 rhs
// need EVLU
// zdot - lambda z = \tilde I
// \tilde I = diag(\lambda) EV^-1 G^-1 I
inline void EV_BASE::compute_zdot(gsl_vector_complex_const_view z0_left, gsl_vector_complex_view z0dot)
{
	trace1("do_ddc G^-1 I", *rhs_);
	size_t fin_dim = z0dot.vector.size;
	gsl_vector_complex_const_view lambda__ = gsl_vector_complex_const_subvector(lambda_, 0, fin_dim);
	gsl_vector_complex_view ztmp = gsl_vector_complex_subvector(ztmp_, 0, fin_dim);
	gsl_vector_view r = gsl_vector_complex_real(ztmp_);
	gsl_vector_view i = gsl_vector_complex_imag(ztmp_);
	gsl_vector_set_zero(&i.vector);
	gsl_vector_memcpy(&r.vector, rhs_); // ztmp = G^-1 I
	trace1("untransformed", *ztmp_);
	gsl_linalg_complex_LU_svx(EVLU_, LU_perm, ztmp_); // ztmp = EV^-1 G^-1 I
	trace1("transformed", *ztmp_);
	gsl_vector_memcpy(rhs_, &r.vector);// for now.
	trace2("newgrid", *rhs_, *lambda_);
	gsl_vector_complex_memcpy(z0dot, z0_left);
	trace1("", z0dot.vector);
	gsl_vector_sub(z0dot, ztmp);
	gsl_vector_mul(z0dot, lambda__);
//	gsl_vector_mul(ztmp, lambda__);
	trace1("limitprime", ztmp);
	trace1("done", z0dot);
}
/*--------------------------------------------------------------------------*/
#if 0
EV_BASE::compute_canonical_op()
{
	if(_n_inputs){
		// minimize | DCS z - z0limit |
		size_t numEV = min(n,_n_inputs + fin_dim + 1); // FIXME: _n_inputs + fin_dim == n should be simpler.
		untested();
/*--------------------------------------------------------------------------*/
// compute z0dot from z0 and G^-1 rhs
// need EVLU
// zdot - lambda z = \tilde I
// \tilde I = diag(\lambda) EV^-1 G^-1 I
		gsl_vector_view l = gsl_vector_complex_real(z0limit_);
		gsl_matrix* s = gsl_matrix_alloc(numEV, _n_inputs);
		gsl_matrix_complex_const_view dcs_left = gsl_matrix_complex_const_submatrix(dc_sens_, 0, 0, numEV, _n_inputs);
		gsl_matrix_memcpy_real(s, dcs_left);
		gsl_multifit_linear_workspace *w = gsl_multifit_linear_alloc(numEV, _n_inputs);

		gsl_vector_view t = gsl_vector_subvector(tmp_,0,_n_inputs);
		trace1("mf", *s);
		trace1("mf", l);
		trace1("mf", t);
		gsl_matrix* cov = gsl_matrix_alloc(_n_inputs, _n_inputs);
		double chisq;
		gsl_multifit_linear(s, &l.vector, &t.vector, cov, &chisq, w);

		gsl_multifit_linear_free(w);
		gsl_matrix_free(cov);
		trace1("coords", *tmp_);

		gsl_vector *tmp2 = gsl_vector_alloc(numEV);
		gsl_multifit_linear_residuals ( s, &l.vector, &t.vector, tmp2);
		trace1("coords", *tmp2);
		gsl_vector_free(tmp2);
	}else{ incomplete();
	}
}
#elif 0
// orthogonal approach
inline void EV_BASE::compute_canonical_op(gsl_vector_complex_view op, gsl_matrix_complex_const_view dcsens)
{
	trace1("compute_canonical_op", *z0_);
//	size_t n = _sim->_total_nodes;

	gsl_complex* op_k = gsl_vector_complex_ptr(z0_, _n_states + _n_inputs);
	// fixme: use Q and Z somehow?
	gsl_complex sp;
	for(unsigned i=0; i<_n_inputs; ++i) {
		gsl_vector_complex_const_view
		s = gsl_matrix_complex_const_column(&dcsens.matrix,i);
		gsl_blas_zdotc(&op.vector, &s.vector, &sp);
	   gsl_vector_complex_scale(&op.vector, gsl_complex_sub(one, sp) );

		gsl_complex* k = gsl_vector_complex_ptr(z0_, _n_states + i);
		*k = gsl_complex_add(*k, gsl_complex_mul(sp, *op_k));
	}

	// now blast is minimal. (and orthogonal to dcsens).
	
	for(unsigned i=0; i<_n_inputs; ++i) {
		gsl_complex* op_k = gsl_vector_complex_ptr(z0_, _n_states + i);
		trace2("inputs vs inputs", *op_k, *_sweepval[_inputno[i]]);
	}

// 	gsl_vector_complex_const_view b = gsl_vector_complex_const_subvector(z0_, );
}
#else
#endif
/*--------------------------------------------------------------------------*/
// gsl_complex& operator[] (gsl_vector_complex_const_view v, const int i)
// {
// 	return *gsl_vector_complex_ptr(v,i);
// }
/*--------------------------------------------------------------------------*/

#if 1
#else
// quantize z0, stash index into index
// z0 = \pm (2**(\pm index)) * tol
// index = floor ( ln (z0 + reltol/ reltol) / ln(2) )
inline void EV_BASE::quantize_state(gsl_vector_complex_view z0, gsl_vector_complex_const_view /*lambda*/, gsl_vector_const_view tol, int* index)
{ incomplete();
	for(unsigned i=0; i<_n_states; ++i) {
		gsl_complex z = gsl_vector_complex_get(z0,i);
		double t_ = gsl_vector_get(tol,i);
		gsl_complex t = gsl_complex_rect(t_,0); // hmmm...
		int sign = 1-2* int(bool(signbit(GSL_REAL(z))));
		trace2("log discretized", z, sign);
		z = gsl_complex_mul_real(z, double(sign));
		trace2("log discretized", z, sign);

		z = gsl_complex_div(z, t);
		trace2("scale", z, t);
//		z = gsl_complex_add(z,gsl_complex_rect(1,0));
		z = gsl_complex_log(z);

		double ln2 = logf(2.);
		z = gsl_complex_add( gsl_complex_rect(ln2, 0. ),z) ;
		trace2("hmm", z, ln2);
		z = gsl_complex_div( z, gsl_complex_rect(ln2 ,0.));

//		GSL_SET_REAL(z, max(0.,GSL_REAL(*z))); // negative real means too close to zero.
		if (GSL_REAL(z) < 0){
			sign = 0;
		}

		index[i] = sign*int(floor(GSL_REAL(z)));
		GSL_SET_REAL(&z, sign * GSL_REAL(z));
		trace4("log quantized", z, index[i], t_, sign);

		if (index[i]){
			z = gsl_complex_rect( sign*t_ * (pow(2,sign*index[i])), 0 );
		} else {
			z = zero;
		}
		trace3("log back", z, index[i], gsl_vector_complex_get(z0,i));
// 		assert(sign * GSL_REAL(z) >= 0.);
		trace3("log back", z, index[i], gsl_vector_complex_get(z0,i));
		gsl_vector_complex_set(&z0.vector,i,z);
	}
	// if err; throw...
}
#endif
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
#endif
