#include <mpi.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <cassert>
#include <utility>

#include "my-malloc.hpp"
#include "simple.hpp"

int iterations = 1000;

using namespace std;

int main(int argc, char** argv) {
  int rank, size;              // rank and size in MPI_COMM_WORLD
  int dims[3] = {0, 0, 0};     // dimensions for MPI_Dims_create
  int periods[3] = {1, 1, 1};  // periods for the cartesion grid
  MPI_Comm comm_cart;          // cartesian communicator
  double* surface_data_out[6]; // surface pointer
  double* surface_data_in[6];  // surface pointer
  int xminus, xplus, yminus,   // neighbors in the cartesian grid
      yplus, zminus, zplus, rank_source;
  int N = 0;                 // size of the local grid N^3
  double* in = nullptr;      // input grid
  double* out = nullptr;     // output grid
  MPI_Win data_win;          // window for data exchange
  double* baseptr = nullptr; // base pointer of the data window

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  parse_argv_simple(argc, argv, rank, N, iterations);

  // create topology
  int res = MPI_Dims_create(size, 3, dims);
  assert(res == 0);

  if (rank == 0)
    printf("dims : %dx%dx%d\n", dims[0], dims[1], dims[2]);

  res = MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &comm_cart);
  assert(comm_cart != MPI_COMM_NULL);
  assert(res == 0);

  // find neighbors (x,y,z)
  res = MPI_Cart_shift(comm_cart, 0, -1, &rank_source, &xminus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_cart, 0, +1, &rank_source, &xplus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_cart, 1, -1, &rank_source, &yminus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_cart, 1, +1, &rank_source, &yplus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_cart, 2, -1, &rank_source, &zminus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_cart, 2, +1, &rank_source, &zplus);
  assert(res == 0);

  int neighbors[] = {xplus, xminus, yplus, yminus, zplus, zminus};

  int remote_offset[] = {1, 0, 3, 2, 5, 4};

  // init data
  for (int i = 0; i < 6; i++) {
    surface_data_out[i] = (double*)my_malloc(N * N * sizeof(double));
    memset(surface_data_out[i], 0, N * N * sizeof(double));
  }

  in = (double*)my_malloc((N + 2) * (N + 2) * (N + 2) * sizeof(double));
  out = (double*)my_malloc((N + 2) * (N + 2) * (N + 2) * sizeof(double));

  memset(in, 0, (N + 2) * (N + 2) * (N + 2) * sizeof(double));
  memset(out, 0, (N + 2) * (N + 2) * (N + 2) * sizeof(double));

  res = MPI_Win_allocate(6 * N * N * sizeof(double), sizeof(double),
                         MPI_INFO_NULL, comm_cart, &baseptr, &data_win);
  assert(res == 0);

  memset(baseptr, 0, 6 * N * N * sizeof(double));

  // data input comes from the window
  for (int i = 0; i < 6; i++)
    surface_data_in[i] = baseptr + i * N * N;

  double start = omp_get_wtime();
  for (int epoch = 1; epoch < iterations + 1; epoch++) {
    // 1. open exposure and access epoch
    res = MPI_Win_fence(0 /*assert*/, data_win); //@todo
    assert(res == 0);

    // 2. pack surface
    pack_surface_simple(in, surface_data_out, N); // local -> surface_data_out[]

    // 3. puts into neigboring nodes
    for (int i = 0; i < 6; i++) {
      res = MPI_Put(surface_data_out[i], N * N, MPI_DOUBLE, neighbors[i],
                    remote_offset[i] * (N * N), N * N, MPI_DOUBLE, data_win);
      assert(res == 0);
    }

    // 4. update local grid
    update_local_grid_simple(out, in, N);

    // 5. close exposure and access epoch
    res = MPI_Win_fence(0 /*assert*/, data_win); //@todo
    assert(res == 0);

    // 7. update surface
    copy_in_neighbors_data_simple(in, surface_data_in, N);
    update_surface_simple(out, in, N);

    // 8. swap in and out
    swap(in, out);
  }
  double stop = omp_get_wtime();

  verify_result_simple(in, iterations, rank, N);

  // cleanup
  my_free(in);
  my_free(out);
  for (int i = 0; i < 6; i++)
    my_free(surface_data_out[i]);
  MPI_Win_free(&data_win);

  double local_duration = stop - start;
  double max_duration;
  res = MPI_Reduce(&local_duration, &max_duration, 1, MPI_DOUBLE, MPI_MAX, 0,
                   MPI_COMM_WORLD);
  assert(res == 0);

  if (rank == 0)
    printf("time : %fs\n", max_duration);
  MPI_Finalize();
  return 0;
}
