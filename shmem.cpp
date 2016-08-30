#include "shmem.hpp"
#include "config.h"
#include "openmp.hpp"

#include <mpi.h>

#include <cassert>

#define IN(X, Y, Z) in[(X) + (N + 2) * (Y) + (N + 2) * (N + 2) * (Z)]
#define OTHER(X, Y, Z) other[(X) + (N + 2) * (Y) + (N + 2) * (N + 2) * (Z)]

void copy_in_neighbors_data_shmem(double* restrict in,
                                  const double* restrict shmem_in,
                                  int rank_shmem, int i, int N) {
  const double* restrict other = shmem_in + rank_shmem * (N + 2) * (N + 2) * (N + 2);
  
  switch (i) {
  case 0: // xplus
    FOR(collapse(2) simd)
    for (int y = 1; y < N + 1; y++)
      for (int z = 1; z < N + 1; z++) {
        IN(N + 1, y, z) = OTHER(1, y, z);
      }
    break;
  case 1: // xminus
    FOR(collapse(2) simd)
    for (int y = 1; y < N + 1; y++)
      for (int z = 1; z < N + 1; z++) {
        IN(0, y, z) = OTHER(N, y, z);
      }
    break;
  case 2: // yplus
    FOR(collapse(2) simd)
    for (int x = 1; x < N + 1; x++)
      for (int z = 1; z < N + 1; z++) {
        IN(x, N + 1, z) = OTHER(x, 1, z);
      }
    break;
  case 3: // yminus
    FOR(collapse(2) simd)
    for (int x = 1; x < N + 1; x++)
      for (int z = 1; z < N + 1; z++) {
        IN(x, 0, z) = OTHER(x, N, z);
      }
    break;
  case 4: // zplus
    FOR(collapse(2) simd)
    for (int x = 1; x < N + 1; x++)
      for (int y = 1; y < N + 1; y++) {
        IN(x, y, N + 1) = OTHER(x, y, 1);
      }
    break;
  case 5: // zminus
    FOR(collapse(2) simd)
    for (int x = 1; x < N + 1; x++)
      for (int y = 1; y < N + 1; y++) {
        IN(x, y, 0) = OTHER(x, y, N);
      }
    break;
  }
}

void parse_argv_shmem(int argc, char** argv, int rank, int& N, int& iterations,
                      int& rankspershmem) {
  if (rank == 0) {
    if (argc != 4) {
      printf("usage: %s <dimension> <iterations> <rankspershmem>\n", *argv);
      exit(1);
    }
    N = atoi(*++argv);
    iterations = atoi(*++argv);
    rankspershmem = atoi(*++argv);
  }

  int res = MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
  assert(res == 0);
  res = MPI_Bcast(&iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);
  assert(res == 0);
  res = MPI_Bcast(&rankspershmem, 1, MPI_INT, 0, MPI_COMM_WORLD);
  assert(res == 0);
}
