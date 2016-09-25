AC_PREREQ(2.63)
AC_INIT(mpi-examples, 0.0.1+git, me@nowhere.com)

LDFLAGS=""

CXXSHELL=$SHELL
if test -f /usr/local/bin/bash; then
  CXXSHELL=/usr/local/bin/bash
elif test -f /bin/bash; then
  CXXSHELL=/bin/bash
fi

AX_PROG_CC_MPI(true, [], AC_MSG_FAILURE([no c mpi compiler found]))
AX_PROG_CXX_MPI(true, [], AC_MSG_FAILURE([no c++ mpi compiler found]))

AC_C_RESTRICT

AX_OPENMP
AC_CHECK_HEADER(mpp/shmem.h)

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
AH_TEMPLATE([HAVE_CRAY], [whether the cray prgenv wrapper is available])
cray_prgenv_wrapper=""
if test "x${PE_ENV}" != "x" ; then
  case "${CC}" in
    CC|cc)
    AC_DEFINE(HAVE_CRAY)
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
#ifdef _CRAYC
  int fail[-1];
#endif
#ifdef __INTEL_COMPILER
  int fail[-1];
  #endif
#ifdef __clang__
  int fail[-1];
#endif
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
elif test "x${me_cv_gnu}" = "xyes" ; then
  CXXFLAGS="${CXXFLAGS}"
  LDFLAGS="${LDFLAGS}"
elif test "x${me_cv_clang}" = "xyes" ; then
  CXXFLAGS="${CXXFLAGS} -std=c++11"
  LDFLAGS="${LDFLAGS}"
else
  CXXFLAGS="${CXXFLAGS} -fopenmp"
  LDFLAGS="${LDFLAGS} -fopenmp"
fi

if test "x${cray_prgenv_wrapper}" = "xyes" ; then
  LDFLAGS="${LDFLAGS} -Wl,--whole-archive,-ldmapp,--no-whole-archive"
fi

AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])


if test "x${ac_cv_header_mpp_shmem_h}" = "xyes"; then
  SHMEM="$(SHMEM)"
else
  SHMEM=""
fi

if test "x${ax_cv_c_openmp}" = "xunknown"; then
  OMP=""
else
  OMP="$(OMP)"
fi


CXXFLAGS="${CXXFLAGS} ${OPENMP_CXXFLAGS}"
LDFLAGS="${LDFLAGS} ${LIBS} ${OPENMP_CXXFLAGS}"

AC_SUBST(CXX)
AC_SUBST(CXXFLAGS)
AC_SUBST([LDFLAGS])
AC_SUBST([PASSIVE_QTHREADS])
AC_SUBST(SHMEM)
AC_SUBST(OMP)
AC_SUBST(CXXSHELL)

AC_CONFIG_HEADERS(config.h)
AC_CONFIG_FILES(Makefile)

AC_OUTPUT
