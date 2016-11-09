# device models hand coded in C++
D_SRCS = \
d_admit.cc d_cap.cc d_cccs.cc d_ccvs.cc d_coil.cc \
d_meas.cc \
d_meter.cc \
d_ttcap.cc \
d_cs.cc d_poly_g.cc d_poly_cap.cc d_res.cc d_switch.cc \
d_trln.cc d_vcg.cc d_vcr.cc d_vcvs.cc d_vs.cc  \
d_branchvs.cc

# matrix_example_LDADD = libgnucap.la

# behavioral modeling functions
BM_SRCS = \
bm_complex.cc bm_exp.cc bm_fit.cc bm_generator.cc \
bm_model.cc bm_poly.cc bm_posy.cc bm_pulse.cc bm_pwl.cc bm_sffm.cc \
bm_sin.cc bm_tanh.cc bmm_table.cc bmm_semi.cc
# parameter functions
FUNC_SRCS = \
func_core.cc measure_eval.cc \
measure_max.cc measure_min.cc measure_slewrate.cc measure_cross.cc \
measure_integral.cc measure_average.cc measure_rms.cc measure_at.cc

# apps
plugindir = ${libdir}/$(shell echo gnucap | sed -e "$(transform)")
plugin_LTLIBRARIES = gnucap-default-plugins.la

# simulation commands, usually multiple files for one command
S_SRCS = \
	s_ac.cc s_dc.cc s_fo.cc s_sens.cc \
	s_tr.cc s_tr_set.cc s_tr_swp.cc \
	s_ttt.cc s_tt_set.cc


gnucap_default_plugins_la_LDFLAGS = -shared -module -avoid-version
gnucap_default_plugins_la_SOURCES = \
	$(BM_SRCS) $(D_SRCS) $(FUNC_SRCS) \
	$(S_SRCS) \
	patchlev.cc \
	signon.cc
gnucap_default_plugins_la_CPPFLAGS = \
	-I$(top_builddir)/include -I$(top_srcdir)/include \
	-I$(top_builddir)/src -I$(top_srcdir)/src

p:
	echo $(gnucap_default_plugins_la_SOURCES)
