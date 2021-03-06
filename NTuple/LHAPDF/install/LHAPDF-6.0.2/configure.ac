## Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([LHAPDF],[6.0.2],[lhapdf@projects.hepforge.org],[LHAPDF])
if test "$prefix" = "$PWD"; then
  AC_MSG_ERROR([Installation into the build directory is not supported: use a different --prefix argument])
fi
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE()
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])


## Library version flags (at least update the version comment with each new release version)
# TODO: add -avoid-version to the LDFLAGS to get rid of the versioned libs (we don't use them anyway)
VERSIONFLAGS="-version-info 2:0:0" ## for 6.0.2
AC_SUBST(VERSIONFLAGS)
AC_DEFINE_UNQUOTED(LHAPDF_VERSION, "$PACKAGE_VERSION", "LHAPDF version string")
# TODO: improve this sed'ing so that it'll work with version digits > 9
PACKAGE_VERSION_CODE=[`echo $PACKAGE_VERSION | sed -e 's/\./0/g' -e 's/\([0-9]\+\).*/\1/'`]
AC_DEFINE_UNQUOTED(LHAPDF_VERSION_CODE, $PACKAGE_VERSION_CODE, "LHAPDF version as an int")


## Mac OS X compatibility
AC_CHECK_TOOL(SWVERS, sw_vers)
if test x$SWVERS != x; then
  PROD_NAME=$($SWVERS -productName | cut -f 2 -d:)
fi
AM_CONDITIONAL(WITH_OSX, [test "$PROD_NAME" == "Mac OS X"])
if test "$PROD_NAME" == "Mac OS X"; then
  MACOSX_DEPLOYMENT_TARGET=$($SWVERS -productVersion | cut -f 1,2 -d.)
  AC_MSG_NOTICE([MACOSX_DEPLOYMENT_TARGET = $MACOSX_DEPLOYMENT_TARGET])
  AM_CXXFLAGS="$AM_CXXFLAGS -Dunix"
  AM_CXXFLAGS="$AM_CXXFLAGS -flat_namespace"
fi


## Checks and config headers
AC_LANG(C++)
AC_CONFIG_SRCDIR([src/GridPDF.cc])
AC_CONFIG_HEADERS([config/config.h include/LHAPDF/Version.h])

## Checks for programs
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_LIBTOOL


## Enable debugging symbols.
AC_ARG_ENABLE([debug],
  [AC_HELP_STRING(--enable-debug, [build with debugging symbols (default=no)])],
  [enable_debug=yes], [enable_debug=no])
if test x$enable_debug == xyes; then
  [AM_CXXFLAGS="$AM_CXXFLAGS -g"]
fi


## Enable LHAGLUE compatibility fns for Fortran and the old C++ interface
AC_ARG_ENABLE([lhaglue],
  [AC_HELP_STRING(--disable-lhaglue, [build without LHAPDF5 compatibility routines])],
  [], [enable_lhaglue=yes])
if test x$enable_lhaglue == xyes; then
   AC_MSG_NOTICE(Building LHAGLUE wrappers)
else
   AC_MSG_NOTICE(Not building LHAGLUE wrappers)
fi
AM_CONDITIONAL(ENABLE_LHAGLUE, [test x$enable_lhaglue == xyes])


## Find pdflatex if possible
AC_PATH_PROG(PDFLATEX, pdflatex)
AM_CONDITIONAL(WITH_PDFLATEX, [test x$PDFLATEX != x])


## Find Doxygen if possible
AC_PATH_PROG(DOXYGEN, doxygen)
AM_CONDITIONAL(WITH_DOXYGEN, [test x$DOXYGEN != x])


## Find Python and Cython if possible
AC_ARG_ENABLE([python],
  [AC_HELP_STRING(--disable-python, [disable the build of the Python interface (default=enabled)])],
  [], [enable_python=yes])
if test x$enable_python == xyes; then
  AC_PATH_PROG(PYTHON, python)
  if test x$PYTHON != x; then
    PYTHON_VERSION=`$PYTHON -c "import sys; print '.'.join(map(str, sys.version_info@<:@:2@:>@));"`
    AC_SUBST(PYTHON_VERSION)
    PYTHON_PATH=`$PYTHON -c "import distutils.sysconfig; print distutils.sysconfig.get_python_lib(prefix='$prefix', plat_specific=True);"`
    AC_SUBST(PYTHON_PATH)
  fi
  AM_CHECK_CYTHON([0.18], [:], [:])
  if test x$CYTHON_FOUND = xyes; then
    AC_MSG_NOTICE([Cython >= 0.18 found: Python extension source can be rebuilt (for developers)])
  fi
fi
AM_CONDITIONAL(WITH_PYTHON, [test x$PYTHON != x])
AM_CONDITIONAL(WITH_CYTHON, [test x$CYTHON_FOUND = xyes])


dnl ## Build Octave interface if possible
dnl AC_ARG_ENABLE([octave],
dnl   [AC_HELP_STRING(--disable-octave, [don't try to make Octave interface])],
dnl   [], [enable_octave=yes])
dnl if test x$enable_octave == xyes; then
dnl    AC_PATH_PROG(OCTAVE, octave)
dnl    AC_PATH_PROG(MKOCTFILE, mkoctfile)
dnl fi
dnl AM_CONDITIONAL(WITH_OCTAVE, [test x$OCTAVE != x && test x$MKOCTFILE != x])


## Set base compiler flags
AC_CEDAR_CHECKCXXFLAG([-std=c++98], [AM_CXXFLAGS="$AM_CXXFLAGS -std=c++98 "])
AC_CEDAR_CHECKCXXFLAG([-pedantic], [AM_CXXFLAGS="$AM_CXXFLAGS -pedantic "])
AC_CEDAR_CHECKCXXFLAG([-Wall], [AM_CXXFLAGS="$AM_CXXFLAGS -Wall "])
AC_CEDAR_CHECKCXXFLAG([-Wno-long-long], [AM_CXXFLAGS="$AM_CXXFLAGS -Wno-long-long "])


## Include $prefix in the compiler flags for the rest of the configure run
if test x$prefix != xNONE; then
  CPPFLAGS="$CPPFLAGS -I$prefix/include"
  LDFLAGS="$LDFLAGS -L$prefix/lib"
fi

## Boost utility library
BOOST_REQUIRE([1.49.0])
BOOST_FOREACH
BOOST_SMART_PTR
BOOST_STRING_ALGO
BOOST_FIND_HEADER([boost/lexical_cast.hpp])
BOOST_FIND_HEADER([boost/range.hpp])
dnl BOOST_FIND_HEADER([boost/range/algorithm/sort.hpp])
BOOST_FILESYSTEM
BOOST_LDFLAGS="$BOOST_FILESYSTEM_LDFLAGS $BOOST_SYSTEM_LDFLAGS"
BOOST_LIBS="$BOOST_FILESYSTEM_LIBS $BOOST_SYSTEM_LIBS"
AC_SUBST(BOOST_LDFLAGS)
AC_SUBST(BOOST_LIBS)


# yaml-cpp parser library
AC_ARG_WITH([yaml-cpp], AS_HELP_STRING([--with-yaml-cpp], [indicate the position of the yaml-cpp installation (should contain both lib and include dirs]))
AC_ARG_WITH([yaml-cpp-lib], AS_HELP_STRING([--with-yaml-cpp-lib], [indicate the position of the yaml-cpp lib directory (the one containing libyaml-cpp.so)]))
AC_ARG_WITH([yaml-cpp-inc], AS_HELP_STRING([--with-yaml-cpp-inc], [indicate the position of the yaml-cpp include directory (the one *containing* the 'yaml-cpp' directory)]))
YAMLCPP_CPPFLAGS=
test "$with_yaml_cpp" && YAMLCPP_CPPFLAGS="-I$with_yaml_cpp/include"
test "$with_yaml_cpp_inc" && YAMLCPP_CPPFLAGS="-I$with_yaml_cpp_inc"
YAMLCPP_LDFLAGS=
test "$with_yaml_cpp" && test -d $with_yaml_cpp/lib && YAMLCPP_LDFLAGS="-L$with_yaml_cpp/lib"
test "$with_yaml_cpp" && test -d $with_yaml_cpp/lib64 && YAMLCPP_LDFLAGS="-L$with_yaml_cpp/lib64"
test "$with_yaml_cpp_lib" && YAMLCPP_LDFLAGS="-L$with_yaml_cpp_lib"
YAMLCPP_LIBS="-lyaml-cpp"
oldCPPFLAGS=$CPPFLAGS
oldLDFLAGS=$LDFLAGS
CPPFLAGS="$AM_CPPFLAGS $CPPFLAGS $YAMLCPP_CPPFLAGS $BOOST_CPPFLAGS"
LDFLAGS="$AM_LDFLAGS $LDFLAGS $YAMLCPP_LDFLAGS"
AC_CHECK_HEADER([yaml-cpp/yaml.h], [], [AC_MSG_ERROR([yaml-cpp/yaml.h header not found.])])
AC_CHECK_HEADER([yaml-cpp/node.h], [YAMLCPP_MAJOR_VERSION=3])
AC_CHECK_HEADER([yaml-cpp/node/node.h], [YAMLCPP_MAJOR_VERSION=5])
if test "x$YAMLCPP_MAJOR_VERSION" = x; then
  AC_MSG_ERROR([yaml-cpp API version could not be determined])
fi
AC_MSG_NOTICE([yaml-cpp API version = $YAMLCPP_MAJOR_VERSION])
AC_SUBST(YAMLCPP_MAJOR_VERSION)
CPPFLAGS=$oldCPPFLAGS
LDFLAGS=$oldLDFLAGS
AC_SUBST(YAMLCPP_CPPFLAGS)
AC_SUBST(YAMLCPP_LDFLAGS)
AC_SUBST(YAMLCPP_LIBS)


## Export base compiler flags
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_CXXFLAGS)
AC_SUBST(AM_LDFLAGS)


## Output
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([include/Makefile include/LHAPDF/Makefile])
AC_CONFIG_FILES([examples/Makefile])
AC_CONFIG_FILES([tests/Makefile])
AC_CONFIG_FILES([doc/Makefile doc/Doxyfile])
AC_CONFIG_FILES([lhapdf-config])
AC_CONFIG_FILES([wrappers/Makefile])
AC_CONFIG_FILES([wrappers/python/Makefile wrappers/python/setup.py])
dnl AC_CONFIG_FILES([wrappers/octave/Makefile])
AC_OUTPUT

## Message about getting PDF sets
pkgdatadir=`eval "echo $datarootdir/$PACKAGE_TARNAME"`
echo
echo "****************************************************************"
echo "IMPORTANT INFORMATION ABOUT PDF SETS"
echo
echo "LHAPDF no longer bundles PDF set data in the package tarball."
echo "The sets are instead all stored online at"
echo "  http://www.hepforge.org/archive/lhapdf/pdfsets/"
echo "and you should install those that you wish to use into"
echo "  $pkgdatadir"
echo
echo "In LHAPDFv6, downloadable sets are packaged as tarballs, which"
echo "must be expanded to be used. Here is an example of how to"
echo "retrieve and install a PDF set from the command line:"
echo
echo " wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.0/CT10nlo.tar.gz -O- | tar xz -C $pkgdatadir"
echo
echo "****************************************************************"
echo
