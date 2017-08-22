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

# utility commands, usually one file per command
# some have multiple commands in a file
C_SRCS = \
c_clear.cc c_comand.cc c_delete.cc c_exp.cc  c_genrat.cc \
c_help.cc \
c_list.cc c_measure.cc c_modify.cc c_param.cc c_prbcmd.cc \
c_status.cc c_sweep.cc c_sim.cc c_system.cc c_for.cc \
c_printf.cc

# apps
plugindir = ${libdir}/$(shell echo gnucap | sed -e "$(transform)")
plugin_LTLIBRARIES = gnucap-default-plugins.la

# simulation commands, usually multiple files for one command
S_SRCS = \
	s_ac.cc s_dc.cc s_fo.cc s_sens.cc \
	s_tr.cc s_tr_set.cc s_tr_swp.cc \
	s_ttt.cc s_tt_set.cc

# netlist languages
LANG_SRCS = \
lang_spice.cc lang_spectre.cc lang_verilog.cc

## The modelgen generated files
MODELSRCS= ${MODELS:.model=.cc} ${MODELS:.model=.h}

gnucap_default_plugins_la_LDFLAGS = -shared -module -avoid-version
gnucap_default_plugins_la_SOURCES = \
	$(BM_SRCS) $(D_SRCS) $(FUNC_SRCS) \
	$(C_SRCS) \
	$(S_SRCS) $(LANG_SRCS) \
	${MODELS} \
	patchlev.cc \
	signon.cc
gnucap_default_plugins_la_CPPFLAGS = \
	-I$(top_builddir)/include -I$(top_srcdir)/include \
	-I$(top_builddir)/src -I$(top_srcdir)/src

EXTRA_DIST = $(MODELS)

### Clean out the modelgen generated files
CLEANFILES = ${MODELSRCS}
# make all of the model object files depend on all of the model sources
# this is a bit overkill, but for example the bjt model makes use of the
# diode model and we want to be sure and capture that dependency
#
${MODELOBJS}: ${MODELSRCS}

## Suffix rules for modelgen.
##
%.h : %.model ${MODELGEN}${EXEEXT}
	${MODELGEN} -h $<

MODELS = \
d_mos1.model d_mos2.model d_mos3.model \
d_mos4.model d_mos5.model d_mos6.model d_mos7.model d_mos8.model \
d_mos123.model d_mos_base.model d_mos.model \
d_bjt.model d_diode.model 
#d_rcd.model d_bti.model 

.model.cc:

BUILT_SOURCES = ${MODELS:%.model=%.h} ${MODELS:%.model=%.cc}


SUFFIXES= .model
MODELGEN= @MODELGEN@
%.cc : %.model %.h ${MODELGEN}${EXEEXT}
	${MODELGEN} -cc $<

MAKEDEPEND = gcc -M $(CPPFLAGS) -o $*.d $<

p:
	echo $(gnucap_default_plugins_la_SOURCES)
