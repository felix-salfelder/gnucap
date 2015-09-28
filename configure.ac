#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
#                            COPYRIGHT
#
#  Copyright (C) 2005, 2006 Dan McMahill
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
#  02110-1301, USA.

AC_PREREQ(2.67)
AC_INIT([gnucap-uf],[0.0.2-rc2],[felix@salfelder.org])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([-Wno-portability silent-rules])


AC_PROG_CXX
AC_PROG_MAKE_SET

AC_CONFIG_SRCDIR([src/main.cc])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([src/config.h])
AC_LANG(C++)

# this will break Makefile dependencies (d'0h)
# AM_MAINTAINER_MODE

AC_PREFIX_DEFAULT("/usr/local")

# libdir_const="${prefix}/lib"

#sysconfdir='/bar/etc'

#AC_SUBST([PREFIX])
#AC_SUBST([LIBDIR])

dnl ----------------- doc ------------------
AC_MSG_CHECKING([doc])
AC_ARG_ENABLE([doc],
[  --enable-doc            Enable building of modules. [[default: disabled]]],
[
if test "X$enable_doc" = "Xno" ; then
	AC_MSG_RESULT([no])
else
	AC_MSG_RESULT([yes])
	enable_doc=yes
fi
],
[
	AC_MSG_RESULT([no])
	enable_doc=no
])
AC_OUTPUT([doc/Makefile])
AM_CONDITIONAL([DOC], [test x$enable_doc != xno] )
dnl ----------------- man ------------------
AC_MSG_CHECKING([man])
AC_ARG_ENABLE([man],
[  --enable-man            Enable building of modules. [[default: disabled]]],
[
if test "X$enable_man" = "Xno" ; then
	AC_MSG_RESULT([no])
else
	AC_MSG_RESULT([yes])
	enable_man=yes
fi
],
[
	AC_MSG_RESULT([no])
	enable_man=no
])

dnl ------------- modules -------------------
dnl MODULES=
AC_MSG_CHECKING([if shared modules should be compiled and installed])
AC_ARG_ENABLE([modules],
[  --enable-modules            Enable building of modules. [[default: disabled]]],
[
if test "X$enable_modules" = "Xno" ; then
	AC_MSG_RESULT([no])
else
	AC_MSG_RESULT([yes])
	enable_modules=yes
fi
],
[
	AC_MSG_RESULT([yes])
	enable_modules=yes
])
AC_OUTPUT([modules/Makefile])
AM_CONDITIONAL([MODULES], [test x$enable_modules != xno] )
dnl --------------- lib -----------------
dnl SHARED_LIB=
AC_MSG_CHECKING([if library shall be installed])
AC_ARG_ENABLE([lib],
[  --enable-lib            Build and install shared library. [[default: disabled]]],
[
if test "X$enable_lib" = "Xno" ; then
	AC_MSG_RESULT([no])
else
	AC_MSG_RESULT([yes])
	enable_lib=yes
fi
],
[
	AC_MSG_RESULT([no])
	enable_lib=no
])
AM_CONDITIONAL([SHARED_LIB], [test x$enable_lib != xno] )

dnl ------------------- debugging -------------------
AC_MSG_CHECKING([if debug code should be compiled in])
AC_ARG_ENABLE([debug],
[  --enable-debug            Enable building of debug code. [[default: disabled]]],
[
if test "X$enable_debug" = "Xno" ; then
	AC_MSG_RESULT([no])
else
	AC_MSG_RESULT([yes])
	enable_debug=yes
fi
],
[
	AC_MSG_RESULT([no])
	enable_debug=no
])

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_MKDIR_P
AC_CHECK_PROG(GIT,git,yes)
AC_MSG_CHECKING([for git repo])
AS_IF([test x"$GIT" = "xyes" && test -d $srcdir/.git],
      [have_git_repo=yes],
      [have_git_repo=no])
AC_MSG_RESULT([$have_git_repo])

WHICH_MAKE=`which make`
AC_SUBST([WHICH_MAKE])

AH_TEMPLATE(HAVE_GIT_REPO, "whether this is built from git")
AH_TEMPLATE(GIT_BRANCH, "the current branch")
AH_TEMPLATE(GIT_TAG, "the description")
AH_TEMPLATE(GIT_COMMIT, "the commit id")
AH_TEMPLATE(GIT_COMMIT_SHORT, "the commit id short")

AM_CONDITIONAL([HAVE_GIT_REPO], [test yes = $have_git_repo])
AS_IF([test yes = $have_git_repo],
      [AC_MSG_NOTICE([building from git])
       git_commit=`git rev-parse HEAD`
       git_commit_short=`git rev-parse --short HEAD`
       git_tag=`git describe`
       git_branch=`git rev-parse --abbrev-ref HEAD`
       AC_DEFINE([HAVE_GIT_REPO],[1])
      ],[
       AC_DEFINE([HAVE_GIT_REPO],[0])
      ])


# hmmm...
# AC_SUBST(GIT_COMMIT, $git_commit)
AC_SUBST(GIT_COMMIT_SHORT, $git_commit_short)
# AC_SUBST(GIT_BRANCH, $git_branch)
# AC_SUBST(GIT_TAG, $git_tag)

AM_CONDITIONAL([VPATH_BUILD], [test ! x$srcdir = "x."])

ac_script='s/\$\$/$/g;s/\\\\/\\/'
reverse=`$as_echo "$program_transform_name" | sed "$ac_script"`
GNUCAP_CPPFLAGS=-I$includedir/`echo gnucap | sed -e "$reverse"`
AC_SUBST([GNUCAP_CPPFLAGS])

# required by attach command
AH_TEMPLATE([GNUCAP_CPPFLAGS])
AC_DEFINE_UNQUOTED([GNUCAP_CPPFLAGS],["$GNUCAP_CPPFLAGS"])

# If we are cross compiling, then we need to search for a
# gnucap-modelgen program to use for our build.  This can
# either be an installed modelgen or it can be specified
# like:
#   env MODELGEN=/build/i686--linux/modelgen/gnucap-modelgen /srcs/gnucap/configure
#           --host=alpha--netbsd --build=i686--linux
#
if test "$cross_compiling" = yes; then
	AC_PATH_PROG(MODELGEN, gnucap-modelgen)
else
	MODELGEN=../modelgen/gnucap-modelgen
	AC_SUBST([MODELGEN])
fi

# for building the documentation
#AC_PATH_PROG(DVIPDFM, dvipdfm, notfound)
#AM_CONDITIONAL(MISSING_DVIPDFM, test x$DVIPDFM = xnotfound)
#AC_PATH_PROG(HACHA, hacha, notfound)
#AC_PATH_PROG(HEVEA, hevea, notfound)
#AM_CONDITIONAL(MISSING_HEVEA, test x$HEVEA = xnotfound -o x$HACHA = xnotfound)
#AC_PATH_PROG(LATEX, latex, notfound)
#AC_PATH_PROG(MAKEINDEX, makeindex, notfound)
#AM_CONDITIONAL(MISSING_LATEX, test x$LATEX = xnotfound -o x$MAKEINDEX = xnotfound)

# Checks for libraries.

AC_LANG([C++])
LT_INIT()
AC_CHECK_LIB([m], [cos])

# probably unneeded.
AC_CHECK_LIB([gslcblas], [cblas_dgemv])
AC_CHECK_LIB([gsl], [gsl_blas_dgemm])
AC_CHECK_LIB([gsl], [gsl_fit_linear])
AC_CHECK_FUNCS([sqrt])


dnl ----------------------------------------------------
have_fftw=no
AC_CHECK_HEADERS([rfftw.h], [AC_DEFINE([HAVE_FFTW_H], [1],
                           [Define to 1 if you have fftw])
                            have_fftw=yes
                            ] )

# fixme: complain if HAVE_FFTW
AC_CHECK_LIB([fftw], [fftw_create_plan], LIBS="$LIBS -lfftw")
AC_CHECK_LIB([rfftw], [rfftw_create_plan], LIBS="-lrfftw $LIBS")

dnl ----------------------------------------------------
# some more checks. some should depend on switches. but the switches arent
# there yet
HAVE_CLAPACK_H=0
AC_CHECK_HEADER([atlas/clapack.h],
                [AC_DEFINE([HAVE_CLAPACK_H], [1],
                           [Define to 1 if you have atlas/clapack.h])
                 HAVE_CLAPACK_H=1],
                [AC_MSG_NOTICE([please install package libatlas-dev])])
AM_CONDITIONAL([USE_ATLAS], [test $HAVE_CLAPACK_H -eq 1])

LAPACK_LIBS=
AC_CHECK_LIB([lapack], [dgelss_],
             [LAPACK_LIBS="-llapack"
              AC_DEFINE([HAVE_LIBLAPACK], [1], [Define to 1 if there is liblapack])])
AC_SUBST([LAPACK_LIBS])

AC_CHECK_HEADER([boost/algorithm/string.hpp],
                [AC_DEFINE([HAVE_STRING_HPP], [1],
                           [Define to 1 if there is boost/algorithm/string.hpp])],
                [AC_MSG_NOTICE([please install the boost package])])

# fixme: refine a bit?
have_boost=no
AC_CHECK_HEADERS([boost/assign.hpp],
                 [have_boost=yes
                 AC_DEFINE([HAVE_BOOST], [1], [have boost])
                 ],
                 [AC_MSG_NOTICE(You need the Boost libraries.)
                 ])

AC_CHECK_HEADER([boost/algorithm/string/split.hpp],
                [AC_DEFINE([HAVE_SPLIT_HPP], [1],
                           [Define to 1 if there is boost/algorithm/string/split.hpp])],
                [AC_MSG_NOTICE([please install the boost package])])

AC_CHECK_HEADER([readline/readline.h],
                [AC_DEFINE([HAVE_READLINE_H], [1],
                           [Define to 1 if you have readline/readline.h])],
                [AC_MSG_ERROR([please install the libreadline-dev package])])

AC_CHECK_HEADER([readline/history.h],
                [AC_DEFINE([HAVE_HISTORY_H], [1],
                           [Define to 1 if you have readline/history.h])],
                [AC_MSG_ERROR([please install the libreadline-dev package])])


AC_ARG_WITH([gsl],
[AS_HELP_STRING([--with-gsl],
  [support gls @<:@default=yes@:>@])],
[],
[with_gsl=yes])

HAVE_GSL_H=0
if test "x$with_gsl" != xno ; then
	AC_CHECK_HEADER([gsl/gsl.h],
	[AC_DEFINE([HAVE_GSL_FIT_H], [1],
		[Define to 1 if you have gsl/gsl.h])
		HAVE_GSL_H=1
	],
	[AC_MSG_NOTICE([consider installing the libgsl0-dev package])])

	AC_CHECK_HEADER([gsl/gsl_fit.h],
	[AC_DEFINE([HAVE_GSL_FIT_H], [1],
		[Define to 1 if you have gsl/gsl_fit.h])],
	[AC_MSG_NOTICE([consider installing the libgsl0-dev package])])
fi
AM_CONDITIONAL([USE_GSL], [test $HAVE_GSL_H -eq 1])

AC_ARG_WITH([readline],
[AS_HELP_STRING([--with-readline],
  [support command line editing @<:@default=yes@:>@])],
[],
[with_readline=yes])

if test "x$with_readline" != xno ; then
	AC_CHECK_LIB([termcap], [main])	
	AC_CHECK_LIB([readline], [main])
fi

AC_CHECK_LIB([dl], [dlopen])

# LT_INIT([static shared])

# shared.
#AC_DISABLE_STATIC
#AC_PROG_LIBTOOL

# Checks for header files.
#AC_CHECK_HEADERS([fcntl.h unistd.h])

AC_CHECK_HEADERS([boost/assign.hpp], [],
                 [AC_MSG_ERROR(You need the Boost libraries.)])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

if test "$enable_debug" = "yes" ; then
	CPPFLAGS="$CPPFLAGS -DTRACE_UNTESTED"
else
	CPPFLAGS="$CPPFLAGS -DNDEBUG"
fi

CXXFLAGS+=" -Wall \
   -W -Wno-sign-compare \
   -Wpointer-arith -Wcast-qual \
   -Wwrite-strings -Wconversion \
   -Woverloaded-virtual -O2 -Wlong-long \
   -Wsign-compare -Wcast-align"

# if we have gcc and we've asked for debugging then add lots of -W
if test "x$GCC" = "xyes" -a "$enable_debug" = "yes"; then
	for flag in -DTRACE_UNTESTED -Wall -W -Wno-sign-compare \
		-Wpointer-arith -Wcast-qual \
		-Wwrite-strings -Wconversion \
		-Woverloaded-virtual -O2 -Wlong-long \
		-Wsign-compare -Wcast-align ; do
		case " ${CFLAGS} " in
			*\ ${flag}\ *)
				# flag is already present
				;;
			*)
				CFLAGS="$CFLAGS ${flag}"
				;;
		esac
		case " ${CXXFLAGS} " in
			*\ ${flag}\ *)
				# flag is already present
				;;
			*)
				CXXFLAGS="$CXXFLAGS ${flag}"
				;;
		esac
	done

# remove additional added  -OX
    FOUNDOX="no"
    NEWCXXFLAGS=""
	for flag in ${CXXFLAGS} ; do
        case " -O0 -O1 -O2 -O3 " in
			*\ ${flag}\ *)
              if test "${FOUNDOX}" = "no"; then
                  NEWCXXFLAGS="${NEWCXXFLAGS} ${flag}"
              fi
              FOUNDOX="yes"
              ;;
            *)
              NEWCXXFLAGS="${NEWCXXFLAGS} ${flag}"
              ;;
        esac
    done
    CXXFLAGS=${NEWCXXFLAGS}
#    echo "CXXFLAGS: ${CXXFLAGS}"
    FOUNDOX="no"
    NEWCFLAGS=""
	for flag in ${CFLAGS} ; do
        case " -O0 -O1 -O2 -O3 " in
			*\ ${flag}\ *)
              if test "${FOUNDOX}" = "no"; then
                  NEWCFLAGS="${NEWCFLAGS} ${flag}"
              fi
              FOUNDOX="yes"
              ;;
            *)
              NEWCFLAGS="${NEWCFLAGS} ${flag}"
              ;;
        esac
    done
    CFLAGS=${NEWCFLAGS}
#    echo "CFLAGS: ${CFLAGS}"

fi

# exports symbols to plugins
LDFLAGS="$LDFLAGS -rdynamic"

# stupid hack (needed?)
if test "x${prefix}" = "xNONE"; then
		AC_DEFINE_UNQUOTED( PREFIX, "${ac_default_prefix}", "Prefix")
		AC_DEFINE_UNQUOTED( LIBDIR, "${ac_default_prefix}/lib", "Libdir")
		AC_DEFINE_UNQUOTED( SODIR, "${ac_default_prefix}/lib", "Sodir")
else
		AC_DEFINE_UNQUOTED( PREFIX, "${prefix}", "Prefix")
		AC_DEFINE_UNQUOTED( LIBDIR, "${prefix}/lib", "Libdir")
		AC_DEFINE_UNQUOTED( SODIR, "${prefix}/lib", "Sodir")
fi

# todo
#if test "X$enable_man" = "Xyes"; then
	#AC_OUTPUT([
		#man/Makefile
		#man/Addmodel/Makefile
		#man/Behavior/Makefile
		#man/Circuit/Makefile
		#man/Commands/Makefile
		#man/Tech/Makefile
	#])
#fi

# workaround, as we don't want to know the names literally
AC_CONFIG_SUBDIRS()
for i in `shopt -s nullglob; echo extern/*`; do
	AS_IF([test -d $i -a -x $i/configure], [subdirs="$subdirs $i"])
done

#AC_OUTPUT([
#	Makefile
#	doc/Makefile
#	examples/Makefile
#	man/Makefile
#	man/Addmodel/Makefile
#	man/Behavior/Makefile
#	man/Circuit/Makefile
#	man/Commands/Makefile
#	man/Tech/Makefile
#	modelgen/Makefile
#	src/Makefile
#	])

AC_CONFIG_FILES([gnucap-conf], [chmod +x gnucap-conf])
AC_CONFIG_FILES([gnucap.mk], [chmod +x gnucap.mk])
AC_CONFIG_FILES([src/gnucap], [chmod +x src/gnucap])

AC_OUTPUT([
	Makefile
	modelgen/Makefile
	plugins/Makefile
	src/Makefile
	tests/Makefile
	examples/Makefile])

plugindir=$libdir/$package

AC_MSG_RESULT([
** Configuration summary for $PACKAGE $VERSION:

   prefix:     $prefix
   suffix:     $program_suffix
   exec_prefix:$exec_prefix
   plugindir:  $plugindir
   libdir:     $libdir
   CPPFLAGS:   $CPPFLAGS
   CFLAGS:     $CFLAGS
   CXXFLAGS:   $CXXFLAGS
   LDFLAGS:    $LDFLAGS
   LIBS:       $LIBS
   debug:      $enable_debug
   modules:    $enable_modules
   shared_lib: $enable_lib
   doc:        $enable_doc
   fftw:       $have_fftw
   boost:      $have_boost
   externals:   $subdirs
   git repo:   $have_git_repo

])
