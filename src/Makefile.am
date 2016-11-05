## -*- Makefile -*-
##
## $Id: Makefile.am,v 1.8 2010-07-29 12:37:18 felix Exp $
##
##                            COPYRIGHT
##
## This file is part of "Gnucap", the Gnu Circuit Analysis Package
##
##  Copyright (C) 2005, 2006 Dan McMahill
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
##  02110-1301, USA.


#------------------------------------------------------------------------
#------------------------------------------------------------------------
# These source files can be omitted, and supplied as plugins.
# device models defined in high level .model files
MODELS = \
d_mos1.model d_mos2.model d_mos3.model \
d_mos4.model d_mos5.model d_mos6.model d_mos7.model d_mos8.model \
d_mos123.model d_mos_base.model d_mos.model \
d_bjt.model d_diode.model 
#d_rcd.model d_bti.model 

.model.cc:

BUILT_SOURCES = ${MODELS:%.model=%.h} ${MODELS:%.model=%.cc}


# utility commands, usually one file per command
# some have multiple commands in a file
C_SRCS = \
c_clear.cc c_comand.cc c_delete.cc c_exp.cc c_list.cc c_save.cc \
c_measure.cc c_modify.cc c_param.cc c_prbcmd.cc c_status.cc \
c_sweep.cc c_sim.cc c_system.cc c_for.cc c_printf.cc \
c_help.cc

# simulation commands, usually multiple files for one command
S_SRCS = \
s_ac.cc s_dc.cc s_fo.cc s_sens.cc \
s_ddc.cc s_ttt.cc s_tt_set.cc \
s_tr.cc s_tr_set.cc s_tr_swp.cc

# netlist languages
LANG_SRCS = \
lang_spice.cc lang_spectre.cc lang_verilog.cc

#------------------------------------------------------------------------
RAW_SRCS = \
$(C_SRCS) \
$(S_SRCS) \
$(LANG_SRCS)

#------------------------------------------------------------------------
RAW_HDRS = \
startup.h md.h ap.h mode.h constant.h declare.h \
m_cpoly.h m_divdiff.h m_interp.h m_matrix.h m_spline.h m_wave.h \
m_base.h m_expression.h m_phase.h \
io_.h io_error.h io_matrix.h io_misc.h \
u_prblst.h u_probe.h u_sim_data.h u_sdp.h u_status.h u_time_pair.h u_xprobe.h \
s__.h s_tr.h e_adp.h s_tt.h\
m_matrix_extra.h \
e_adplist.h \
d_cap.h \
extlib.h \
u_sock.h s_ddc.h u_atlas.h

noinst_HEADERS = d_coil.h
#------------------------------------------------------------------------ 
## The modelgen generated files
MODELSRCS= ${MODELS:.model=.cc} ${MODELS:.model=.h}

# extlib_o_SOURCES = extlib.h 
#
plugindir_transformed = ${libdir}/$(shell echo gnucap | sed -e "$(transform)")
plugindir = $(plugindir_transformed)

# henjo:
MODELCODE= ${MODELS:.model=.cc}
MODELHDR= ${MODELS:.model=.h}

MODELOBJS= ${MODELS:.model=.${OBJEXT}}

EXTRA_PROGRAMS = m_matrix_example
m_matrix_example_SOURCES = m_matrix_example.cc # io.cc io_error.cc
m_matrix_example_CPPFLAGS = $(AM_CPPFLAGS)


SRCS1= 
SRCS2= ${RAW_SRCS} 

# make all of the model object files depend on all of the model sources
# this is a bit overkill, but for example the bjt model makes use of the
# diode model and we want to be sure and capture that dependency
#
${MODELOBJS}: ${MODELSRCS}

## Include the .model files in the distfile as well as the rest of the
## non-generated sources

EXTRA_DIST= ${MODELS} ${SRCS1} ${SRCS2} \
test_readline.cc spice-wrapper.cc \
Makefile.am

GNUCAP_MAKE=$(shell echo gnucap.mk | sed -e '$(transform)')

noinst_SCRIPTS = gnucap
bin_PROGRAMS = gnucap.real
transform=s&\.real$$&&;$(program_transform_name)
gnucap_real_SOURCES = ${SRCS2} startup.cc main.cc ${SRCS1} ${MODELS}

gnucap_real_CPPFLAGS = -I$(top_builddir)/include -I$(top_srcdir)/include
gnucap_real_CPPFLAGS+= $(AM_CPPFLAGS) -DGNUCAP_MAKE=\"$(GNUCAP_MAKE)\"

#for now (ddc). fix later
gnucap_real_LDFLAGS = $(LAPACK_LIBS)

gnucap_real_LDADD = ../lib/libgnucap-uf.la
#libgnucap_la_CFLAGS = $(AM_CFLAGS)


includedir_transformed = ${includedir}/$(shell echo gnucap | sed -e "$(transform)")
libgnucap_includedir = $(includedir_transformed)
libgnucap_include_HEADERS = ${RAW_HDRS} spice-wrapper.cc
 
#gnucap_include_HEADERS = ${RAW_HDRS}

gnucap_includedir = $(pkgincludedir)


### Clean out the modelgen generated files
CLEANFILES = ${MODELSRCS}

## Suffix rules for modelgen.
##

SUFFIXES= .model
MODELGEN= @MODELGEN@

MAKEDEPEND = gcc -M $(CPPFLAGS) -o $*.d $<

%.cc : %.model %.h ${MODELGEN}${EXEEXT}
	${MODELGEN} -cc $<

# #############################################

%.h : %.model ${MODELGEN}${EXEEXT}
	${MODELGEN} -h $<

$(D_IVL_O): CPPFLAGS+= IVL_CPPFLAGS

# FIXME: -I iverilog only for .v stuff

.DELETE_ON_ERROR: *.cc

# include $(srcdir)/depend
#
# %.o: .deps/%.

#tags:
#	ctags -R --sort=yes --c++-kinds=+p --fields=+iaS --extra=+q --language-force=C++ -f tags .
#
.PRECIOUS: %.cc %.h %.hxx
#------------------------------------------------------------------------
#------------------------------------------------------------------------
EXTRA_DIST += c_expect.cc


# embers of vvp attempts. sort out later
EXTRA_DIST+= \
 array.h \
 d_clock.v \
 d_digital.cc \
 d_digital.h \
 d_invert.v \
 d_ivl_ports.cc \
 d_ivl_ports.h \
 d_latch.v \
 d_lfsr.v \
 d_vvp.cc \
 d_vvp.h \
 delay.h \
 e_ivl_compile.cc \
 e_ivl_compile.h \
 extlib.cc \
 pointers.h


AM_DEFAULT_SOURCE_EXT = .cc
