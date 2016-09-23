## $Id: Makefile.am,v 1.1 2009-10-23 12:01:23 felix Exp $
##
## Top level automake file for gnucap

ACLOCAL_AMFLAGS=-I m4

doc_DATA = COPYING

if MODULES
  MODULES_DIR = modules
endif

if DOC
  DOC_DIR = doc
endif

MODELGEN = modelgen
DIST_SUBDIRS = $(SUBDIRS) doc

SUBDIRS = \
    include \
    lib \
    main \
    $(DOC_DIR) \
    examples \
    $(MODELGEN) \
    src \
    $(MODULES_DIR) \
    apps \
    plugins \
    tests

EXTRA_DIST = \
    extern/README \
    Doxyfile \
    bootstrap \
    gnuplot.mk \
    BUGS \
    CHANGES.uf
