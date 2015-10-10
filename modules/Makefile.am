# this file is part of gnucap-uf
# (c) 2009-2013 Felix Salfelder
# license: GPLv3+
#
#------------------------------------------------------------------------
plugindir = ${libdir}/$(shell echo gnucap | sed -e "$(transform)")
plugin_LTLIBRARIES =

AM_LDFLAGS = -shared -module -avoid-version

AM_CPPFLAGS = \
	-I$(top_builddir)/include -I$(top_srcdir)/include \
	-I$(top_builddir)/src -I$(top_srcdir)/src

AM_DEFAULT_SOURCE_EXT = .cc

#------------------------------------------------------------------------
plugin_LTLIBRARIES += d_gpolyk_wrap.la
d_gpolyk_wrap_la_CPPFLAGS = -DADD_VERSION $(AM_CPPFLAGS)
#------------------------------------------------------------------------
plugin_LTLIBRARIES += rcd_sym.la
rcd_sym_la_SOURCES = d_rcd_sym.cc d_rcd_exp.cc
#------------------------------------------------------------------------
plugin_LTLIBRARIES += bti_inf.la
bti_inf_la_SOURCES = d_bti_inf.cc
bti_inf_la_CPPFLAGS = -DADD_VERSION $(AM_CPPFLAGS)
#------------------------------------------------------------------------
plugin_LTLIBRARIES += d_brl.la
d_brl_la_SOURCES = d_brl.cc
d_brl_la_CPPFLAGS = -DADD_VERSION $(AM_CPPFLAGS)
#------------------------------------------------------------------------
plugin_LTLIBRARIES += d_ttint.la
d_ttint_la_SOURCES = d_ttint.cc
d_ttint_la_CPPFLAGS = -DADD_VERSION $(AM_CPPFLAGS)
#------------------------------------------------------------------------
plugin_LTLIBRARIES += d_branchvs.la
d_branchvs_la_SOURCES = d_branchvs.cc
d_branchvs_la_CPPFLAGS = -DADD_VERSION $(AM_CPPFLAGS)
#------------------------------------------------------------------------
plugin_LTLIBRARIES += s_sock.la
s_sock_la_SOURCES = s_sock.cc
s_sock_la_CPPFLAGS = -DADD_VERSION $(AM_CPPFLAGS)
#------------------------------------------------------------------------
if USE_GSL
SEVLA = s_ev.la
s_ev_la_SOURCES = s_ev.cc s_reach.cc u_gsl.cc u_ssp.cc
# s_ev_la_CPPFLAGS = -DADD_VERSION $(AM_CPPFLAGS)
else
SEVLA =
endif
plugin_LTLIBRARIES += $(SEVLA)
noinst_HEADERS = u_ssp.h u_gsl.h s_ev.h
#------------------------------------------------------------------------
plugin_LTLIBRARIES += c_wave.la
c_wave_la_SOURCES = c_wave.cc
c_wave_la_CPPFLAGS = -DADD_VERSION $(AM_CPPFLAGS)
c_wave_la_CXXFLAGS = -std=c++11 $(AM_CXXFLAGS)
#------------------------------------------------------------------------
plugin_LTLIBRARIES += c_debug.la
c_debug_la_LDFLAGS = -module -avoid-version
c_debug_la_CPPFLAGS = -DADD_VERSION $(AM_CPPFLAGS)
#------------------------------------------------------------------------
plugin_LTLIBRARIES += c_nodeset.la
c_nodeset_la_LDFLAGS = -module -avoid-version
c_nodeset_la_CPPFLAGS = -DADD_VERSION $(AM_CPPFLAGS)
#------------------------------------------------------------------------
plugin_LTLIBRARIES += models-uf.la
models_uf_la_SOURCES = \
	patchlev.cc \
	e_adp_mos.cc \
	d_diode.cc \
	d_mos_base.cc \
	d_mos.cc \
	d_mos123.cc \
	d_mos1.cc \
	d_mos8.cc \
	d_bti.cc \
	d_hci.cc \
	d_rcd.cc
noinst_HEADERS += \
	d_bjt.h \
	d_bti.h \
	d_diode.h \
	d_rcd.h \
	e_adp_mos.h \
	d_mos_base.h \
	d_mos.h \
	d_mos1.h \
	d_mos123.h \
	d_mos8.h
#------------------------------------------------------------------------
plugin_LTLIBRARIES += lars_bjt.la
lars_bjt_la_LDFLAGS = -module -avoid-version
lars_bjt_la_SOURCES = d_bjt.cc
lars_bjt_la_CPPFLAGS = -DLARS_BJT_HACK -DADD_VERSION $(AM_CPPFLAGS)

#needed?
includedir_transformed = ${includedir}/$(shell echo gnucap | sed -e "$(transform)")
mincludedir = ${includedir}
minclude_HEADERS = d_rcd_sym.h d_bti_inf.h d_rcd_sym_v2.h

# cruft, cleanup later
EXTRA_DIST = \
d_rcd_sym_v2.cc \
d_rcd_sym_v3.cc \
d_rcd_sym_v3.h \
d_rcd_sym_v4.cc

