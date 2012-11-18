
# Keep the version in sync with src/fftw3.mli
AC_INIT(ocaml-fftw3, 0.5.2, chris_77@users.sourceforge.net)

AC_MSG_RESULT([configuring $PACKAGE_STRING])

OCAMLFIND_LDCONF=""
AC_ARG_ENABLE([ldconf],
  AC_HELP_STRING([--disable-ldconf],
    [don't modify the dynamic loader configuration file (default is enable)]),
  [ac_enable_ldconf=$enableval],
  [ac_enable_ldconf=$enableval],
  [ac_enable_ldconf=yes])
if test "$ac_enable_ldconf" = no ; then
	AC_MSG_RESULT([disabling modification of ld.conf])
	OCAMLFIND_LDCONF=dummy
fi

# Check for Ocaml compilers

# we first look for ocamlc in the path; if not present, we fail
AC_PATH_PROG(OCAMLC,ocamlc,no)
if test "$OCAMLC" = no ; then
	AC_MSG_ERROR(Cannot find ocamlc.)
fi

# Directory of ocamlc in $OCAMLC
OCAMLBIN=`dirname $OCAMLC`

# Ocaml version number and library path
OCAML_VERSION=`$OCAMLC -version`
echo "ocaml version is $OCAML_VERSION"
OCAMLLIBPATH=`$OCAMLC -where`
echo "ocaml library path is $OCAMLLIBPATH"

# We look for ocamlopt; if not present, we issue a warning
# if the version is not the same, we also discard it
# we set OCAMLBEST to "opt" or "byte", whether ocamlopt is available or not
AC_PATH_PROG(OCAMLOPT,ocamlopt,no)
OCAMLBEST=byte
if test "$OCAMLOPT" = no ; then
	AC_MSG_WARN(Cannot find ocamlopt; bytecode compilation only.)
else
	AC_MSG_CHECKING(ocamlopt version)
	TMPVERSION=`$OCAMLOPT -v | sed -n -e 's|.*version* *\(.*\)$|\1|p' `
	if test "$TMPVERSION" != "$OCAML_VERSION" ; then
	    AC_MSG_RESULT(differs from ocamlc; ocamlopt discarded.)
	    OCAMLOPT=no
	else
	    AC_MSG_RESULT(ok)
	    OCAMLBEST=opt
	fi
fi

# Checking for ocamlc.opt
AC_PATH_PROG(OCAMLCDOTOPT,ocamlc.opt,no)
if test "$OCAMLCDOTOPT" != no ; then
	AC_MSG_CHECKING(ocamlc.opt version)
	TMPVERSION=`$OCAMLCDOTOPT -version`
	if test "$TMPVERSION" != "$OCAML_VERSION" ; then
	    AC_MSG_RESULT(differs from ocamlc; ocamlc.opt discarded.)
	else
	    AC_MSG_RESULT(ok)
	    OCAMLC=$OCAMLCDOTOPT
	fi
fi

# Checking for ocamlopt.opt
if test "$OCAMLOPT" != no ; then
    AC_PATH_PROG(OCAMLOPTDOTOPT,ocamlopt.opt,no)
    if test "$OCAMLOPTDOTOPT" != no ; then
	AC_MSG_CHECKING(ocamlc.opt version)
	TMPVER=`$OCAMLOPTDOTOPT -version`
	if test "$TMPVER" != "$OCAML_VERSION" ; then
	    AC_MSG_RESULT(differs from ocamlc; ocamlopt.opt discarded.)
	else
	    AC_MSG_RESULT(ok)
	    OCAMLOPT=$OCAMLOPTDOTOPT
	fi
    fi
fi


AC_PATH_PROG(OCAMLDEP,ocamldep,no)
if test "$OCAMLDEP" = no ; then
	AC_MSG_ERROR(Cannot find ocamldep.)
fi

AC_PATH_PROG(OCAMLDOC,ocamldoc,no)

AC_PATH_PROG(OCAMLFIND,ocamlfind,no)
if test "$OCAMLFIND" = no ; then
	AC_MSG_ERROR(Cannot find ocamlfind.)
fi

AC_PROG_CC()

# Add the math lib
CLIBS="$CLIBS m"

# libfftw3 detection (put in CLIBS)
AC_MSG_CHECKING([for libfftw3 (double precision)])
LIBS="$LIBS -lfftw3"
AC_TRY_LINK([#include <fftw3.h>],
            [const int n = 10;
	    fftw_complex in[10], out[10];
	    fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);],
            [AC_MSG_RESULT([ok]); CLIBS="$CLIBS fftw3"],
            [AC_MSG_ERROR([Check for fftw3 (double precision) failed.])])

# Keep the order of "type r2r_kind" in fftw3SD.ml in sync with fftw3.h values
AC_TRY_RUN([#include <fftw3.h>
            int main(void) {
              if(FFTW_R2HC==0 && FFTW_HC2R==1 && FFTW_DHT==2 &&
                 FFTW_REDFT00==3 && FFTW_REDFT01==4 && FFTW_REDFT10==5 &&
                 FFTW_REDFT11==6 && FFTW_RODFT00==7 && FFTW_RODFT01==8 &&
                 FFTW_RODFT10==9 && FFTW_RODFT11==10)
                return(0);
              else
                return(1);
            }],
	[],
	[AC_MSG_FAILURE([The values of the fields FFTW_R2HC,... in fftw3.h has
           changed. Please use a newer version of the OCaml-fftw3 library
           or contact AC_PACKAGE_BUGREPORT to report the issue.])])


AC_MSG_CHECKING([for libfftw3f (single precision)])
LIBS="$LIBS -lfftw3f"
AC_TRY_LINK([#include <fftw3.h>],
            [const int n = 10;
	    fftwf_complex in[10], out[10];
	    fftwf_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);],
            [AC_MSG_RESULT([ok]); CLIBS="$CLIBS fftw3f";
	     FFTW3F_EXISTS="-DFFTW3F_EXISTS"],
            [AC_MSG_RESULT(
	      [failed, the single precision module will be disabled])])

# Check for OCaml libraries

# OCAML_LIB(lib, action-if-found, action-if-not-found)
# action-if-found can use the variable ocaml_lib_path
AC_DEFUN([OCAML_LIB],[
  AC_ARG_WITH([$1],
    AC_HELP_STRING([--with-$1=PATH], [set the path of the OCaml library $1]),
    [], [])

  AC_MSG_CHECKING(for OCaml library $1)
  # $1 substituted when the script configure is written
  if test "$with_$1" != "" ; then
    if test "$with_$1" = "yes" ; then
      AC_MSG_ERROR(Option --with-$1=PATH with PATH not provided)
    else
      if test "$with_$1" = "no" ; then
        AC_MSG_RESULT(disabled)
	$3
      else
        ocaml_lib_path="$with_$1"
	AC_MSG_RESULT($ocaml_lib_path)
        $2
      fi
    fi
  else
    ocaml_lib_path=`$OCAMLFIND query $1 2> /dev/null || echo no`
    if test "$ocaml_lib_path" != "no" ; then
      AC_MSG_RESULT($ocaml_lib_path)
      $2
    else
      # Look in the standard OCaml directory
      ocaml_lib=`find $OCAMLLIBPATH -name $1.cma | head -n 1`
      if test "$ocaml_lib" != "" ; then
        ocaml_lib_path=`dirname $ocaml_lib`
	AC_MSG_RESULT($ocaml_lib_path)
        $2
      else
	AC_MSG_RESULT(not found; use --with-$1=PATH to indicate its location $3)
      fi
    fi
  fi
])

OCAML_LIB(lacaml, [LACAML_PATH=$ocaml_lib_path],
		  [(Some examples will not be compiled.)])
OCAML_LIB(gnuplot, [GNUPLOT_PATH=$ocaml_lib_path],
		   [(Examples will not produce graphics.)])



if test "$OCAMLOPT" = no ; then
    BEST=byte
    ALL=byte
else
    BEST=opt
    ALL="byte opt"
fi

# substitutions to perform
AC_SUBST(PACKAGE_VERSION)
AC_SUBST(OCAMLC)
AC_SUBST(OCAMLOPT)
AC_SUBST(OCAMLDEP)
AC_SUBST(OCAMLBEST)
AC_SUBST(OCAML_VERSION)
AC_SUBST(OCAMLLIBPATH)
AC_SUBST(OCAMLBIN)
AC_SUBST(OCAMLDOC)
AC_SUBST(OCAMLFIND)
AC_SUBST(OCAMLFIND_LDCONF)
AC_SUBST(OCAMLCP)
AC_SUBST(BEST)
AC_SUBST(ALL)
AC_SUBST(INC)
AC_SUBST(CLIBS)
AC_SUBST(CFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(FFTW3F_EXISTS)
AC_SUBST(LACAML_PATH)
AC_SUBST(GNUPLOT_PATH)

# Finally create the Makefile and samples (read only to avoid changing
# them by mistake)
AC_CONFIG_FILES([Makefile],[chmod a-w Makefile])
AC_CONFIG_FILES([src/META],[chmod a-w src/META])
AC_CONFIG_FILES([src/Makefile],[chmod a-w src/Makefile])
AC_CONFIG_FILES([examples/Makefile],[chmod a-w examples/Makefile])
AC_OUTPUT
