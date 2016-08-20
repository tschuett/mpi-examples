#include "simple.hpp"

#include <mpi.h>

#include <cassert>
#include <cstdio>
#include <cstdlib>

#define IN(X,Y,Z) in[(X) + (N+2)*(Y) + (N+2)*(N+2)*(Z)]
#define OUT(X,Y,Z) out[(X) + (N+2)*(Y) + (N+2)*(N+2)*(Z)]

#define SURFACE(S, X,Y) (S)[(X) + (N)*(Y)]


void pack_surface_simple(const double *in, double *surface_data_out[6], int N) {
  //xplus, xminus
  for(int y = 1; y < N+1; y++)
  for(int z = 1; z < N+1; z++) {
    SURFACE(surface_data_out[0], y-1, z-1) = IN(N,y,z);
    SURFACE(surface_data_out[1], y-1, z-1) = IN(1,y,z);
  }

  //yplus, yminus
  for(int x = 1; x < N+1; x++)
  for(int z = 1; z < N+1; z++) {
    SURFACE(surface_data_out[2], x-1, z-1) = IN(x,N,z);
    SURFACE(surface_data_out[3], x-1, z-1) = IN(x,1,z);
  }

  //zplus, zminus
  for(int x = 1; x < N+1; x++)
  for(int y = 1; y < N+1; y++) {
    SURFACE(surface_data_out[4], x-1, y-1) = IN(x,y,N);
    SURFACE(surface_data_out[5], x-1, y-1) = IN(x,y,1);
  }
}

void update_local_grid_simple(double *out, const double *in, int N) {
  for(int x = 2; x < N; x++)
  for(int y = 2; y < N; y++)
  for(int z = 2; z < N; z++)
    OUT(x,y,z) = 1 + (IN(x-1,y,z) + IN(x+1,y,z)
                    + IN(x,y-1,z) + IN(x,y+1,z)
                    + IN(x,y,z-1) + IN(x,y,z+1)) / 6.0;
}

void copy_in_neighbors_data_simple(double *in, double *surface_data_in[6], int N) {
  //xplus, xminus
  for(int y = 1; y < N+1; y++)
  for(int z = 1; z < N+1; z++) {
    IN(N+1,y,z) = SURFACE(surface_data_in[0], y-1, z-1);
    IN(0  ,y,z) = SURFACE(surface_data_in[1], y-1, z-1);
  }

  //yplus, yminus
  for(int x = 1; x < N+1; x++)
  for(int z = 1; z < N+1; z++) {
    IN(x,N+1,z) = SURFACE(surface_data_in[2], x-1, z-1);
    IN(x,  0,z) = SURFACE(surface_data_in[3], x-1, z-1);
  }

  //zplus, zminus
  for(int x = 1; x < N+1; x++)
  for(int y = 1; y < N+1; y++) {
    IN(x,y,N+1) = SURFACE(surface_data_in[4], x-1, y-1);
    IN(x,y,  0) = SURFACE(surface_data_in[5], x-1, y-1);
  }
}

void update_surface_simple(double *out, const double *in, int N) {
  // xplus,xminus
  for(int y = 1; y < N+1; y++)
    for(int z = 1; z < N+1; z++) {
      int x = 1;
      OUT(x,y,z) = 1 + (IN(x-1,y,z) + IN(x+1,y,z)
                      + IN(x,y-1,z) + IN(x,y+1,z)
                      + IN(x,y,z-1) + IN(x,y,z+1)) / 6.0;
      x = N;
      OUT(x,y,z) = 1 + (IN(x-1,y,z) + IN(x+1,y,z)
                      + IN(x,y-1,z) + IN(x,y+1,z)
                      + IN(x,y,z-1) + IN(x,y,z+1)) / 6.0;
    }

  // yplus,yminus
  for(int x = 1; x < N+1; x++)
    for(int z = 1; z < N+1; z++) {
      int y = 1;
      OUT(x,y,z) = 1 + (IN(x-1,y,z) + IN(x+1,y,z)
                      + IN(x,y-1,z) + IN(x,y+1,z)
                      + IN(x,y,z-1) + IN(x,y,z+1)) / 6.0;
      y = N;
      OUT(x,y,z) = 1 + (IN(x-1,y,z) + IN(x+1,y,z)
                      + IN(x,y-1,z) + IN(x,y+1,z)
                      + IN(x,y,z-1) + IN(x,y,z+1)) / 6.0;
    }

  // zplus,zminus
  for(int x = 1; x < N+1; x++)
    for(int y = 1; y < N+1; y++) {
      int z = 1;
      OUT(x,y,z) = 1 + (IN(x-1,y,z) + IN(x+1,y,z)
                      + IN(x,y-1,z) + IN(x,y+1,z)
                      + IN(x,y,z-1) + IN(x,y,z+1)) / 6.0;
      z = N;
      OUT(x,y,z) = 1 + (IN(x-1,y,z) + IN(x+1,y,z)
                      + IN(x,y-1,z) + IN(x,y+1,z)
                      + IN(x,y,z-1) + IN(x,y,z+1)) / 6.0;
    }
}

void verify_result_simple(double *in, int iterations, int rank, int N) {
  double sum = 0;
  double global_sum = 0.0;
  size_t errors = 0;

  for(int x = 1; x < N+1; x++)
  for(int y = 1; y < N+1; y++)
  for(int z = 1; z < N+1; z++) {
    sum += IN(x,y,z) - iterations;
    if(rank == 0)
      if(IN(x,y,z) != iterations) {
        errors++;
        printf("%d %d %d %f %f\n", x, y, z, IN(x,y,z), (double)iterations);
      }
  }

  int res = MPI_Reduce(&sum, &global_sum, 1,
                       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  assert(res == 0);

  if(rank == 0)
    printf("error: %f %f %ld %d\n", sum, sum / (N*N*N), errors, iterations);
}

void parse_argv_simple(int argc, char **argv, int rank, int& N) {
  if(rank == 0) {
    if(argc != 2) {
      printf("usage: %s <dimension>\n", *argv);
      exit(1);
    }
    N  = atoi(*++argv);
  }

  int res = MPI_Bcast(&N, 1, MPI_INT, 0,
                      MPI_COMM_WORLD);
  assert(res == 0);
}
