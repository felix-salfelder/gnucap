## -*- Makefile -*-
##
##                            COPYRIGHT
##
## This file is part of "Gnucap", the Gnu Circuit Analysis Package
##
##  Copyright (C) 2015 Felix Salfelder
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

includedir_transformed = ${includedir}/$(shell echo gnucap | sed -e "$(transform)")
libgnucap_includedir = $(includedir_transformed)
libgnucap_include_HEADERS = \
	patchlev.h m_mvpoly.h io_trace.h \
	d_coment.h d_dot.h d_logic.h d_subckt.h \
	u_cardst.h u_function.h u_lang.h u_limit.h u_nodemap.h u_opt.h u_parameter.h

EXTRA_DIST = patchlev.h.in
BUILT_SOURCES = patchlev.h
if HAVE_GIT_REPO
CLEANFILES = patchlev.h
# the values for GIT_BRANCH and GIT_COMMIT are unknown to configure, so this is
# done manually (for now). see config.status for a more portable way.
#
# don't touch patchlev.h if it's not necessary, to prevent rebuilds.
#
patchlev.h: $(srcdir)/patchlev.h.in $(top_srcdir)/.git
	@echo "// this file has been created from patchlev.h.in, don't edit" > patchlev.tmp
	@echo "#line 1 \"$<\"" >> patchlev.tmp
	@sed -e "s\\@""GIT_BRANCH@\\`git rev-parse --abbrev-ref HEAD`\\" \
	     -e "s/@""INTERFACE_VERSION@/0x$(GIT_COMMIT_SHORT)/" \
	     -e "s/@""PACKAGE_VERSION@/$(PACKAGE_VERSION)/" \
	     -e "s\\@""INTERFACE_NAME@\\\"`git rev-parse --abbrev-ref HEAD`\"\\" \
	     -e "s/@""GIT_COMMIT@/`git rev-parse HEAD`/" $< >> patchlev.tmp
	@diff patchlev.tmp patchlev.h 2>/dev/null >/dev/null; \
	ret=$$?; \
	if [ $$ret -eq 0 ]; then \
	  echo git revision unchanged; \
	  rm patchlev.tmp; \
	else \
	  echo updating $@; \
	  mv patchlev.tmp $@; \
	fi;
# else patchlev.h is already there
endif
