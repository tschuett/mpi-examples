#include "shmem.hpp"
#include "config.h"
#include "openmp.hpp"

#include <mpi.h>

#include <cassert>

#define IN(X, Y, Z)                                                            \
  in[(X) + dims[1] * (N + 2) * (Y) +                                           \
     dims[1] * dims[2] * (N + 2) * (N + 2) * (Z)]
#define OUT(X, Y, Z)                                                           \
  out[(X) + dims[1] * (N + 2) * (Y) +                                          \
      dims[1] * dims[2] * (N + 2) * (N + 2) * (Z)]

#define OTHER(X, Y, Z)                                                         \
  other[(X) + dims[1] * (N + 2) * (Y) +                                        \
        dims[1] * dims[2] * (N + 2) * (N + 2) * (Z)]

/*
void copy_in_neighbors_data_shmem(double* restrict in,
                                  const double* restrict shmem_in,
                                  int rank_shmem, int i, int N) {
  const double* restrict other =
      shmem_in + rank_shmem * (N + 2) * (N + 2) * (N + 2);

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
*/
void copy_in_neighbors_data_shmem(double* in, const double* surface_data_in,
                                  int i, int N, const int* coords,
                                  const int* dims) {}

void memset_shmem(double* in, int N, int* coords, int* dims) {
  /* @todo somebody has to memset the halo */
  for (int x = 1 + coords[0] * N; x < N + coords[0] * N + 1; x++)
    for (int y = 1 + coords[1] * N; y < N + coords[1] * N + 1; y++)
      for (int z = 1 + coords[2] * N; z < N + coords[2] * N + 1; z++)
        IN(x, y, z) = 0.0;
}

void pack_surface_shmem(const double* in, double* surface_data_out, int i,
                        int N, const int* coords, const int* dims) {}

void update_surface_shmem(double* shmem_out, const double* shmem_in, int i,
                          int N, const int* coords_smp, const int* dims_smp) {
  what ? ? switch (i) {
  case 0: // xplus
    FOR(collapse(2) simd)
    for (int y = 1; y < N + 1; y++)
      for (int z = 1; z < N + 1; z++) {
        int x = 1;
        OUT(x, y, z) = 1 +
                       (IN(x - 1, y, z) + IN(x + 1, y, z) + IN(x, y - 1, z) +
                        IN(x, y + 1, z) + IN(x, y, z - 1) + IN(x, y, z + 1)) /
                           6.0;
        x = N;
        OUT(x, y, z) = 1 +
                       (IN(x - 1, y, z) + IN(x + 1, y, z) + IN(x, y - 1, z) +
                        IN(x, y + 1, z) + IN(x, y, z - 1) + IN(x, y, z + 1)) /
                           6.0;
      }
    break;
  case 1: // xminus
    FOR(collapse(2) simd)
    for (int y = 1; y < N + 1; y++)
      for (int z = 1; z < N + 1; z++) {
        int x = 1;
        OUT(x, y, z) = 1 +
                       (IN(x - 1, y, z) + IN(x + 1, y, z) + IN(x, y - 1, z) +
                        IN(x, y + 1, z) + IN(x, y, z - 1) + IN(x, y, z + 1)) /
                           6.0;
        x = N;
        OUT(x, y, z) = 1 +
                       (IN(x - 1, y, z) + IN(x + 1, y, z) + IN(x, y - 1, z) +
                        IN(x, y + 1, z) + IN(x, y, z - 1) + IN(x, y, z + 1)) /
                           6.0;
      }
    break;
  case 2: // yplus
    FOR(collapse(2) simd)
    for (int x = 1; x < N + 1; x++)
      for (int z = 1; z < N + 1; z++) {
        int y = 1;
        OUT(x, y, z) = 1 +
                       (IN(x - 1, y, z) + IN(x + 1, y, z) + IN(x, y - 1, z) +
                        IN(x, y + 1, z) + IN(x, y, z - 1) + IN(x, y, z + 1)) /
                           6.0;
        y = N;
        OUT(x, y, z) = 1 +
                       (IN(x - 1, y, z) + IN(x + 1, y, z) + IN(x, y - 1, z) +
                        IN(x, y + 1, z) + IN(x, y, z - 1) + IN(x, y, z + 1)) /
                           6.0;
      }
    break;
  case 3: // yminus
    FOR(collapse(2) simd)
    for (int x = 1; x < N + 1; x++)
      for (int z = 1; z < N + 1; z++) {
        int y = 1;
        OUT(x, y, z) = 1 +
                       (IN(x - 1, y, z) + IN(x + 1, y, z) + IN(x, y - 1, z) +
                        IN(x, y + 1, z) + IN(x, y, z - 1) + IN(x, y, z + 1)) /
                           6.0;
        y = N;
        OUT(x, y, z) = 1 +
                       (IN(x - 1, y, z) + IN(x + 1, y, z) + IN(x, y - 1, z) +
                        IN(x, y + 1, z) + IN(x, y, z - 1) + IN(x, y, z + 1)) /
                           6.0;
      }
    break;
  case 4: // zplus
    FOR(collapse(2) simd)
    for (int x = 1; x < N + 1; x++)
      for (int y = 1; y < N + 1; y++) {
        int z = 1;
        OUT(x, y, z) = 1 +
                       (IN(x - 1, y, z) + IN(x + 1, y, z) + IN(x, y - 1, z) +
                        IN(x, y + 1, z) + IN(x, y, z - 1) + IN(x, y, z + 1)) /
                           6.0;
        z = N;
        OUT(x, y, z) = 1 +
                       (IN(x - 1, y, z) + IN(x + 1, y, z) + IN(x, y - 1, z) +
                        IN(x, y + 1, z) + IN(x, y, z - 1) + IN(x, y, z + 1)) /
                           6.0;
      }
    break;
  case 5: // zminus
    FOR(collapse(2) simd)
    for (int x = 1; x < N + 1; x++)
      for (int y = 1; y < N + 1; y++) {
        int z = 1;
        OUT(x, y, z) = 1 +
                       (IN(x - 1, y, z) + IN(x + 1, y, z) + IN(x, y - 1, z) +
                        IN(x, y + 1, z) + IN(x, y, z - 1) + IN(x, y, z + 1)) /
                           6.0;
        z = N;
        OUT(x, y, z) = 1 +
                       (IN(x - 1, y, z) + IN(x + 1, y, z) + IN(x, y - 1, z) +
                        IN(x, y + 1, z) + IN(x, y, z - 1) + IN(x, y, z + 1)) /
                           6.0;
      }
    break;
  }
}

void update_local_grid_shmem(double* out, const double* in, int N,
                             const int* coords, const int* dims,
                             const int* neighbors) {
  for (int x = 2 + coords[0] * N; x < N + coords[0] * N; x++)
    for (int y = 2 + coords[1] * N; y < N + coords[1] * N; y++)
      for (int z = 2 + coords[2] * N; z < N + coords[2] * N; z++)
        OUT(x, y, z) = 1 +
                       (IN(x - 1, y, z) + IN(x + 1, y, z) + IN(x, y - 1, z) +
                        IN(x, y + 1, z) + IN(x, y, z - 1) + IN(x, y, z + 1)) /
                           6.0;
}

void verify_result_shmem(double* in, int iterations, int rank, int N,
                         const int* coords, const int* dims) {
  double sum = 0;
  double global_sum = 0.0;
  size_t errors = 0;

  for (int x = 1 + coords[0] * N; x < N + coords[0] * N + 1; x++)
    for (int y = 1 + coords[1] * N; y < N + coords[1] * N + 1; y++)
      for (int z = 1 + coords[2] * N; z < N + coords[2] * N + 1; z++) {
        sum += IN(x, y, z) - iterations;
        if (rank == 0)
          if (IN(x, y, z) != iterations) {
            errors++;
            printf("%d %d %d %f %f\n", x, y, z, IN(x, y, z),
                   (double)iterations);
          }
      }

  int res =
      MPI_Reduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  assert(res == 0);

  if (rank == 0)
    printf("error: %f %f %ld %d\n", sum, sum / (N * N * N), errors, iterations);
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
