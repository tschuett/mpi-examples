#pragma once

#define string1(x) #x
#define stringize(x) string1(x)

#define PRAGMA(...) _Pragma(stringize(__VAR_ARGS))

#ifdef WITH_OPENMP

#define PARALLEL(...) PRAGMA(omp parallel __VAR_ARGS)
#define FOR(...) PRAGMA(omp for __VAR_ARGS)
#define MASTER PRAGMA(omp master)
#define BARRIER PRAGMA(omp barrier)

#else

#define PARALLEL(...)
#define FOR(...)
#define MASTER
#define BARRIER

#endif
