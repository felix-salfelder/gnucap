
dist_bin_SCRIPTS = gnucap-conf gnucap.mk
CLEANFILES = gnucap-conf gnucap.mk

bin_PROGRAMS = gnucap.real
transform=s&\.real$$&&;$(program_transform_name)
gnucap_real_SOURCES = startup.cc main.cc

gnucap_real_CPPFLAGS = -I$(top_builddir)/include -I$(top_srcdir)/include
gnucap_real_CPPFLAGS+= $(AM_CPPFLAGS) -DGNUCAP_MAKE=\"$(GNUCAP_MAKE)\"

#for now (ddc). fix later
gnucap_real_LDFLAGS = $(LAPACK_LIBS)

gnucap_real_LDADD = ../lib/libgnucap-uf.la
noinst_SCRIPTS = gnucap

# FIXME: obsolete
EXTRA_DIST = startup.h
