#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <cassert>
#include <chrono>
#include <cmath>
#include <utility>

#include "my-malloc.hpp"
#include "shmem.hpp"
#include "simple.hpp"

int iterations = 1000;
int rankspershmem = 2;

const int ndims = 3;

using namespace std;

int main(int argc, char** argv) {
  int rank, size; // rank and size in MPI_COMM_WORLD
  int rank_smp, size_smp, rank_nodes, size_nodes, rank_global;
  int size_smp_min, size_smp_max;
  MPI_Comm comm_nodes_flat, comm_smp_flat, comm_nodes_cart, comm_smp_cart,
      comm_global_flat, comm_global_cart;
  double* surface_data_out[6]; // surface pointer
  double* surface_data_in[6];  // surface pointer
  int xminus, xplus, yminus,   // neighbors in the shmem cartesian grid
      yplus, zminus, zplus, rank_source;
  int global_xminus, global_xplus,
      global_yminus, // neighbors in the global cartesian grid
      global_yplus, global_zminus, global_zplus;
  int N = 0;                   // size of the local grid N^3
  double* shmem_in = nullptr;  // input grid
  double* shmem_out = nullptr; // output grid
  MPI_Win data_win;            // window for global data exchange
  MPI_Win shared_win;          // shared memory window for data exchange
  double* baseptr = nullptr;   // base pointer of the data window
  double* shared_baseptr = nullptr;
  double* local_shared_baseptr = nullptr;

  int dims_nodes[] = {0, 0, 0};
  int periods_nodes[] = {1, 1, 1};

  int dims_smp[] = {0, 0, 0};
  int periods_smp[] = {0, 0, 0}; // non-periodic

  int dims_global[] = {0, 0, 0};    // ?
  int periods_global[] = {1, 1, 1}; // periodic

  int coords_smp[] = {0, 0, 0};
  int coords_nodes[] = {0, 0, 0};
  int coords_global[] = {0, 0, 0};

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  parse_argv_shmem(argc, argv, rank, N, iterations, rankspershmem);

  // create topology
  /* find shared memory nodes */
  int res = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                                MPI_INFO_NULL, &comm_smp_flat);
  assert(res == 0);

  res = MPI_Comm_size(comm_smp_flat, &size_smp);
  assert(res == 0);

  /* cartesian topology per SMP node */
  res = MPI_Dims_create(size_smp, ndims, dims_smp);
  res = MPI_Cart_create(comm_smp_flat, ndims, dims_smp, periods_smp,
                        /*reorder*/ 1, &comm_smp_cart);
  res = MPI_Comm_free(&comm_smp_flat);
  /* local rank and cartesian coordinate*/
  res = MPI_Comm_rank(comm_smp_cart, &rank_smp);
  res = MPI_Cart_coords(comm_smp_cart, rank_smp, ndims, coords_smp);

  /* check that all shared memory nodes have the same size */
  res = MPI_Allreduce(&size_smp, &size_smp_min, 1, MPI_INT, MPI_MIN,
                      MPI_COMM_WORLD);
  res = MPI_Allreduce(&size_smp, &size_smp_max, 1, MPI_INT, MPI_MAX,
                      MPI_COMM_WORLD);
  assert(size_smp_min == size_smp_max);

  /* find the number of shared memory nodes (size_nodes) */
  res = MPI_Comm_split(MPI_COMM_WORLD, rank_smp, 0, &comm_nodes_flat);
  res = MPI_Comm_size(comm_nodes_flat, &size_nodes);

  /* create dimensions of the cartesian grid *over* the nodes */
  MPI_Dims_create(size_nodes, ndims, dims_nodes);
  if (rank_smp == 0) {
    /* find the rank of each shared memory node in the cartesian grid over the
     * nodes */
    res = MPI_Cart_create(comm_nodes_flat, ndims, dims_nodes, periods_nodes, 1,
                          &comm_nodes_cart);
    res = MPI_Comm_rank(comm_nodes_cart, &rank_nodes);
    res = MPI_Comm_free(&comm_nodes_cart);
  }
  res = MPI_Comm_free(&comm_nodes_flat);
  /* broadcast the rank within each shared memory node */
  res = MPI_Bcast(&rank_nodes, 1, MPI_INT, 0, comm_smp_cart);

  /* */
  res = MPI_Comm_split(MPI_COMM_WORLD, rank_smp, rank_nodes, &comm_nodes_flat);
  res = MPI_Cart_create(comm_nodes_flat, ndims, dims_nodes, periods_nodes,
                        0 /*!*/, &comm_nodes_cart);
  res = MPI_Cart_coords(comm_nodes_cart, rank_nodes, ndims, coords_nodes);
  res = MPI_Comm_free(&comm_nodes_flat);

  /* establish the global cart communicator */
  for (int i = 0; i < ndims; i++) {
    dims_global[i] = dims_smp[i] * dims_nodes[i];
    coords_global[i] = coords_nodes[i] * dims_smp[i] + coords_smp[i];
  }
  rank_global = coords_global[0];
  for (int i = 1; i < ndims; i++)
    rank_global = rank_global * dims_global[i] + coords_global[i];

  /* create communicator with ranks sorted by rank_global (color=0) */
  MPI_Comm_split(MPI_COMM_WORLD, 0 /*color*/, rank_global, &comm_global_flat);
  /* create global cartesian communicator (the goal) */
  MPI_Cart_create(comm_global_flat, ndims, dims_global, periods_global, 0,
                  &comm_global_cart);
  MPI_Comm_free(&comm_global_flat);

  // find neighbors in the shared memory node cartesian communicator (x,y,z)
  // neighbors can be MPI_PROC_NULL if the neighbor is in another shared memory
  // node
  res = MPI_Cart_shift(comm_smp_cart, 0, -1, &rank_source, &xminus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_smp_cart, 0, +1, &rank_source, &xplus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_smp_cart, 1, -1, &rank_source, &yminus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_smp_cart, 1, +1, &rank_source, &yplus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_smp_cart, 2, -1, &rank_source, &zminus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_smp_cart, 2, +1, &rank_source, &zplus);
  assert(res == 0);

  // look for neighbors in comm_global_cart
  res = MPI_Cart_shift(comm_global_cart, 0, -1, &rank_source, &global_xminus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_global_cart, 0, +1, &rank_source, &global_xplus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_global_cart, 1, -1, &rank_source, &global_yminus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_global_cart, 1, +1, &rank_source, &global_yplus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_global_cart, 2, -1, &rank_source, &global_zminus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_global_cart, 2, +1, &rank_source, &global_zplus);
  assert(res == 0);

  int smp_neighbors[] = {xplus, xminus, yplus, yminus, zplus, zminus};
  int global_neighbors[] = {global_xplus,  global_xminus, global_yplus,
                            global_yminus, global_zplus,  global_zminus};

  int remote_offset[] = {1, 0, 3, 2, 5, 4};

  // init data
  for (int i = 0; i < 6; i++) {
    surface_data_out[i] = (double*)my_malloc(N * N * sizeof(double));
    memset(surface_data_out[i], 0, N * N * sizeof(double));
  }

  // each rank can have up to 6 non-shmem neighbors
  res = MPI_Win_allocate(6 * N * N * sizeof(double), sizeof(double),
                         MPI_INFO_NULL, comm_global_cart, &baseptr, &data_win);
  assert(res == 0);

  memset(baseptr, 0, 6 * N * N * sizeof(double));

  // init shmem window
  // we need two volumes: in and out
  // the size of one volume is (2 + dims_smp[0]*N)*(2 +*N)*(2 +x*N)
  size_t size_of_volume =
      (2 + dims_smp[0] * N) * (2 + dims_smp[1] * N) * (2 + dims_smp[0] * N);
  // all ranks have to provide the same value to allocate_shared
  size_t local_shmem_size = ceil(size_of_volume / (double)size_smp);
  res = MPI_Win_allocate_shared(local_shmem_size * sizeof(double),
                                sizeof(double), MPI_INFO_NULL, comm_smp_cart,
                                &local_shared_baseptr, &shared_win);
  assert(res == 0);

  int disp_unit;
  MPI_Aint wsize;
  res =
      MPI_Win_shared_query(shared_win, 0, &wsize, &disp_unit, &shared_baseptr);
  assert(res == 0);
  shmem_in = shared_baseptr;
  shmem_out = shared_baseptr + size_smp * (N + 2) * (N + 2) * (N + 2);

  memset_shmem(shmem_in, N, coords_smp, dims_smp);
  memset_shmem(shmem_out, N, coords_smp, dims_smp);

  // data input comes from the window
  for (int i = 0; i < 6; i++)
    surface_data_in[i] = baseptr + i * N * N;

  auto start = std::chrono::high_resolution_clock::now();
  for (int epoch = 1; epoch < iterations + 1; epoch++) {
    // 1. open exposure and access epoch
    res = MPI_Win_fence(0 /*assert*/, data_win); //@todo
    assert(res == 0);

    // 2. pack surface
    for (int i = 0; i < 6; i++)
      if (smp_neighbors[i] == MPI_PROC_NULL)
        pack_surface_shmem(shmem_in, surface_data_out[i], i, N, coords_smp,
                           dims_smp); // local -> surface_data_out[i]

    // 3. puts into neigboring nodes
    for (int i = 0; i < 6; i++) {
      if (smp_neighbors[i] == MPI_PROC_NULL) {
        res =
            MPI_Put(surface_data_out[i], N * N, MPI_DOUBLE, global_neighbors[i],
                    remote_offset[i] * (N * N), N * N, MPI_DOUBLE, data_win);
        assert(res == 0);
      }
    }

    // 4. update local grid
    update_local_grid_shmem(shmem_out, shmem_in, N, coords_smp, dims_smp,
                            smp_neighbors);
    // @todo move computation here
    // for (int i = 0; i < 6; i++)
    //  if (smp_neighbors[i] != MPI_PROC_NULL)
    //    update_local_grid_shmem(out, in, shmem_in, smp_neighbors[i], N);

    // 5. close exposure and access epoch
    res = MPI_Win_fence(0 /*assert*/, data_win); //@todo
    assert(res == 0);

    // 7. update surface
    for (int i = 0; i < 6; i++) {
      if (smp_neighbors[i] != MPI_PROC_NULL)
        copy_in_surface_data_shmem(shmem_in, surface_data_in[i], i, N,
                                   coords_smp, dims_smp);
      update_surface_shmem(shmem_out, shmem_in, i, N, coords_smp, dims_smp);
    }

    // 8. swap in and out
    swap(shmem_in, shmem_out);
  }
  auto stop = std::chrono::high_resolution_clock::now();

  verify_result_shmem(shmem_in, iterations, rank, N, coords_smp, dims_smp);

  // cleanup
  for (int i = 0; i < 6; i++)
    my_free(surface_data_out[i]);
  MPI_Win_free(&data_win);
  MPI_Win_free(&shared_win);

  const std::chrono::duration<double> diff = stop - start;
  const double local_duration = diff.count();
  double max_duration;
  res = MPI_Reduce(&local_duration, &max_duration, 1, MPI_DOUBLE, MPI_MAX, 0,
                   MPI_COMM_WORLD);
  assert(res == 0);

  if (rank == 0)
    printf("time : %fs\n", max_duration);
  MPI_Finalize();
  return 0;
}
