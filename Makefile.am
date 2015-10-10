## $Id: Makefile.am,v 1.1 2009-10-23 12:01:23 felix Exp $
##
## Top level automake file for gnucap

ACLOCAL_AMFLAGS=-I m4

dist_bin_SCRIPTS = gnucap-conf gnucap.mk
doc_DATA = COPYING
CLEANFILES = gnucap-conf gnucap.mk

if MODULES
  MODULES_DIR = modules
endif

if DOC
  DOC_DIR = doc
endif

#later..
MODELGEN =
DIST_SUBDIRS = modelgen $(SUBDIRS) doc

SUBDIRS = \
    include \
    lib \
    $(DOC_DIR) \
    examples \
    $(MODELGEN) \
    src \
    $(MODULES_DIR) \
    plugins \
    tests

EXTRA_DIST = \
    extern/README \
    Doxyfile \
    bootstrap \
    gnuplot.mk \
    BUGS \
    CHANGES.uf
