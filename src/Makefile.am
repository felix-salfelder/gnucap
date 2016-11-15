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
u_prblst.h u_sim_data.h u_sdp.h u_status.h u_time_pair.h u_xprobe.h \
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


## Include the .model files in the distfile as well as the rest of the
## non-generated sources

EXTRA_DIST = ${SRCS1} ${SRCS2} \
test_readline.cc spice-wrapper.cc \
Makefile.am

GNUCAP_MAKE=$(shell echo gnucap.mk | sed -e '$(transform)')

noinst_SCRIPTS = gnucap
bin_PROGRAMS = gnucap.real
transform=s&\.real$$&&;$(program_transform_name)
gnucap_real_SOURCES = ${SRCS2} startup.cc main.cc ${SRCS1}

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



# #############################################

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
