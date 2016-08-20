AC_PREREQ(2.63)
AC_INIT(mpi-examples, 0.0.1+git, me@nowhere.com)

LDFLAGS=""

AC_PROG_CXX([CC mpiCC])
AC_PROG_CC([cc mpicc])

AC_SEARCH_LIBS([hbw_check_available], [memkind])
dnl ac_cv_search_hbw_check_available=-lmemkind
AH_TEMPLATE([HAVE_MEMKIND], [whether the memkind library is available])
if test "x${ac_cv_search_hbw_check_available}" != "xno"; then
   AC_DEFINE(HAVE_MEMKIND)
   LDFLAGS="${LDFLAGS} ac_cv_search_hbw_check_available"
fi

dnl check if a cray prgenv wrapper compiler is being used
cray_prgenv_wrapper=""
if test "x${PE_ENV}" != "x" ; then
  case "${CC}" in
    CC|cc)
	cray_prgenv_wrapper="yes"
	;;
    *)
       ;;
  esac
fi

AC_CACHE_CHECK([whether compiler is cray],
              [me_cv_cray],
              [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
                                                  [
#ifndef _CRAYC
  int fail[-1];
#endif
])],
                            [me_cv_cray=yes],
                            [me_cv_cray=no])])

if test "x${me_cv_cray}" = "xyes" ; then
  CXXFLAGS="${CXXFLAGS} -h std=c++11"
else
  CXXFLAGS="${CXXFLAGS} -fopenmp"
  LDFLAGS="${LDFLAGS} -fopenmp -Wl,--whole-archive,-ldmapp,--no-whole-archive"
fi

AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

AC_SUBST(CXX)
AC_SUBST(CXXFLAGS)
AC_SUBST([LDFLAGS])

AC_CONFIG_HEADERS(config.h)
AC_CONFIG_FILES(Makefile)

AC_OUTPUT
