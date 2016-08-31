AC_PREREQ(2.63)
AC_INIT(mpi-examples, 0.0.1+git, me@nowhere.com)

LDFLAGS=""

AX_PROG_CXX_MPI(true, [], AC_MSG_FAILURE([no mpi compiler found]))

AC_C_RESTRICT


AC_SEARCH_LIBS([hbw_check_available], [memkind])
dnl ac_cv_search_hbw_check_available=-lmemkind
AH_TEMPLATE([HAVE_MEMKIND], [whether the memkind library is available])
if test "x${ac_cv_search_hbw_check_available}" != "xno"; then
   AC_DEFINE(HAVE_MEMKIND)
fi


PASSIVE_QTHREADS=''
AC_SEARCH_LIBS([qthread_initialize], [qthread])
dnl ac_cv_search_qthread_initialize=-lqthread
AH_TEMPLATE([HAVE_QTHREAD], [whether the qthread library is available])
if test "x${ac_cv_search_qthread_initialize}" != "xno"; then
   AC_DEFINE(HAVE_QTHREAD)
   PASSIVE_QTHREADS='cart-passive-qthreads.x'
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

AC_CACHE_CHECK([whether compiler is intel],
              [me_cv_intel],
              [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
                                                  [
#ifndef __INTEL_COMPILER
  int fail[-1];
#endif
])],
                            [me_cv_intel=yes],
                            [me_cv_intel=no])])


AC_CACHE_CHECK([whether compiler is gnu],
              [me_cv_gnu],
              [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
                                                  [
#ifndef __GNUC__
  int fail[-1];
#endif
])],
                            [me_cv_gnu=yes],
                            [me_cv_gnu=no])])


AC_CACHE_CHECK([whether compiler is clang],
              [me_cv_clang],
              [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
                                                  [
#ifndef __clang__
  int fail[-1];
#endif
])],
                            [me_cv_clang=yes],
                            [me_cv_clang=no])])


if test "x${me_cv_cray}" = "xyes" ; then
  CXXFLAGS="${CXXFLAGS} -h std=c++11"
elif test "x${me_cv_intel}" = "xyes" ; then
  CXXFLAGS="${CXXFLAGS} -qopenmp -restrict -qopt-report=5"
  LDFLAGS="${LDFLAGS} -qopenmp"
else
  CXXFLAGS="${CXXFLAGS} -fopenmp"
  LDFLAGS="${LDFLAGS} -fopenmp"
fi

if test "x${cray_prgenv_wrapper}" = "xyes" ; then
  LDFLAGS="${LDFLAGS} -Wl,--whole-archive,-ldmapp,--no-whole-archive"
  MTLDFLAGS="${LDFLAGS} -Wl,--whole-archive,-ldmapp,--no-whole-archive -craympich-mt"
  else
  MTLDFLAGS="${LDFLAGS}"
fi

AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])


LDFLAGS="${LDFLAGS} ${LIBS}"
MTLDFLAGS="${MTLDFLAGS} ${LIBS}"

AC_SUBST(CXX)
AC_SUBST(CXXFLAGS)
AC_SUBST([LDFLAGS])
AC_SUBST([MTLDFLAGS])
AC_SUBST([PASSIVE_QTHREADS])

AC_CONFIG_HEADERS(config.h)
AC_CONFIG_FILES(Makefile)

AC_OUTPUT
