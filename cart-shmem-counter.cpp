#include <mpi.h>
#include <mpp/shmem.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <utility>

#include "my-malloc.hpp"
#include "simple.hpp"

int iterations = 1000;

const uint64_t one = 1;

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
  double* baseptr = nullptr; // base pointer of the data window
  volatile long long* counter_baseptr =
      nullptr; // base pointer of the counter window
  MPI_Group world_group;
  MPI_Group cart_group;
  int neighbors[6];
  int remote_offset[] = {1, 0, 3, 2, 5, 4};

  MPI_Init(&argc, &argv);
  shmem_init();

  rank = shmem_my_pe();
  size = shmem_n_pes();

  parse_argv_simple(argc, argv, rank, N, iterations);

  // create topology
  int res = MPI_Dims_create(size, 3, dims);
  assert(res == 0);

  if (rank == 0)
    printf("dims : %dx%dx%d\n", dims[0], dims[1], dims[2]);

  res = MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, true, &comm_cart);
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

  int cart_neighbors[] = {xplus, xminus, yplus, yminus, zplus, zminus};

  res = MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  assert(res == 0);
  res = MPI_Comm_group(comm_cart, &cart_group);
  assert(res == 0);
  res = MPI_Group_translate_ranks(cart_group, 6, cart_neighbors, world_group,
                                  neighbors);
  assert(res == 0);

  // init data
  for (int i = 0; i < 6; i++) {
    surface_data_in[i] = (double*)my_malloc(N * N * sizeof(double));
    memset(surface_data_in[i], 0, N * N * sizeof(double));
  }

  in = (double*)my_malloc((N + 2) * (N + 2) * (N + 2) * sizeof(double));
  out = (double*)my_malloc((N + 2) * (N + 2) * (N + 2) * sizeof(double));

  memset(in, 0, (N + 2) * (N + 2) * (N + 2) * sizeof(double));
  memset(out, 0, (N + 2) * (N + 2) * (N + 2) * sizeof(double));

  // data window
  baseptr = (double*)shmem_malloc(6 * N * N * sizeof(double));
  assert(baseptr != nullptr);
  memset(baseptr, 0, 6 * N * N * sizeof(double));

  // create counter window
  counter_baseptr = (long long*)shmem_malloc(12 * sizeof(long long));
  assert(counter_baseptr != nullptr);
  memset((void*)counter_baseptr, 0, 12 * sizeof(uint64_t));

  // data output goes into the window
  for (int i = 0; i < 6; i++)
    surface_data_out[i] = baseptr + i * N * N;

  double start = omp_get_wtime();

  shmem_barrier_all();

  for (uint64_t epoch = 1; epoch < iterations + 1; epoch++) {
    // 1. pack surface
    pack_surface_simple(in, surface_data_out, N); // local -> surface_data_out[]

    // 2. signal data availability
    for (int i = 0; i < 6; i++)
      shmem_longlong_add((long long*)counter_baseptr + remote_offset[i], 1,
                         neighbors[i]);

    shmem_quiet(); //?

    // 3. wait for data availability from neighbors
    // is that legal? atomic read? ordering for atomicity in MPI?
    bool have_data = false;
    while (!have_data) {
      have_data =
          (counter_baseptr[0] == epoch) && (counter_baseptr[1] == epoch) &&
          (counter_baseptr[2] == epoch) && (counter_baseptr[3] == epoch) &&
          (counter_baseptr[4] == epoch) && (counter_baseptr[5] == epoch);
    }

    // 4. get data from neighbors
    for (int i = 0; i < 6; i++) {
      shmem_double_get_nbi(surface_data_in[i],
                           baseptr + remote_offset[i] * N * N, N * N,
                           neighbors[i]);
      assert(res == 0);
    }

    // 5. update local grid
    update_local_grid_simple(out, in, N);

    // 6. wait for neigbor's data
    shmem_quiet();

    // 7. signal completion
    for (int i = 0; i < 6; i++)
      shmem_longlong_add((long long*)counter_baseptr + 6 + remote_offset[i], 1,
                         neighbors[i]);

    shmem_quiet();

    // 8. update surface
    copy_in_neighbors_data_simple(in, surface_data_in, N);
    update_surface_simple(out, in, N);

    // 9. wait for completion
    have_data = false;
    while (!have_data) {
      have_data =
          (counter_baseptr[6] == epoch) && (counter_baseptr[7] == epoch) &&
          (counter_baseptr[8] == epoch) && (counter_baseptr[9] == epoch) &&
          (counter_baseptr[10] == epoch) && (counter_baseptr[11] == epoch);
    }

    // 8. swap in and out
    swap(in, out);
  }

  double stop = omp_get_wtime();

  verify_result_simple(in, iterations, rank, N);

  // cleanup
  my_free(in);
  my_free(out);
  for (int i = 0; i < 6; i++)
    my_free(surface_data_in[i]);
  // shmem_free(&baseptr);
  // shmem_free(&counter_baseptr);

  double local_duration = stop - start;
  double max_duration;
  res = MPI_Reduce(&local_duration, &max_duration, 1, MPI_DOUBLE, MPI_MAX, 0,
                   MPI_COMM_WORLD);
  assert(res == 0);

  if (rank == 0)
    printf("time : %fs\n", max_duration);
  shmem_finalize();
  MPI_Finalize();
  return 0;
}
