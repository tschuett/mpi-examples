#include <pmi2.h>
#include <mpi.h>
#include <string.h>

#include <rdma/fi_cm.h>
#include <rdma/fi_rma.h>
#include <rdma/fi_atomic.h>
#include <rdma/fi_errno.h>

#include <sys/uio.h>

#include <cassert>
#include <chrono>
#include <utility>
#include <iostream>

#include "grid.hpp"
#include "simple.hpp"
#include "my-malloc.hpp"
#include "pmi.hpp"
#include "ofi.hpp"

int iterations = 1000;

struct fid_mr *mr_surface; // memory registration key of the surface window
struct fid_mr *mr_counters;

double *surface = nullptr;
double *counters = nullptr;

fi_addr_t neighbors[6];
uintptr_t counter_baseptrs[6];
uint64_t counter_keys[6];
uintptr_t surface_baseptrs[6];
uint64_t surface_keys[6];
size_t gni_address_size = 0;

int N = 0;             // size of the local grid N^3
const uint64_t one = 1;

using namespace std;

void setup_pmi(int rank) {
  uint8_t gni_address[1000];
  char buf[1000];
  size_t value_size = 1000;

  // publish network address
  int res = fi_getname((fid*)ep, gni_address, &value_size);
  assert(res == 0);
  gni_address_size = value_size;
  pmi_put("endpoint", gni_address, gni_address_size, rank);

  value_size = 1000;

  // mr info for surface
  size_t surface_size = 6 * N * N * sizeof(double);
  surface = (double*)my_malloc(surface_size);
  res = fi_mr_reg(domain, surface, surface_size, FI_WRITE | FI_REMOTE_WRITE, 0, 0, 0, &mr_surface, NULL);
  assert(res == 0);

  uint64_t mr_key = fi_mr_key(mr_surface);
  pmi_put("mr-key-surface", &mr_key, sizeof(mr_key), rank);
  pmi_put("mr-baseptr-surface", &surface, sizeof(uintptr_t), rank);

  // mr info for counters
  size_t counters_size = 6 * sizeof(double);
  counters = (double*)my_malloc(counters_size);
  res = fi_mr_reg(domain, counters, counters_size, FI_WRITE | FI_REMOTE_WRITE, 0, 1, 0, &mr_counters, NULL);
  assert(res == 0);

  mr_key = fi_mr_key(mr_counters);
  pmi_put("mr-key-counters", &mr_key, sizeof(uint64_t), rank);
  pmi_put("mr-baseptr-counters", &counters, sizeof(uintptr_t), rank);

  res = PMI2_KVS_Fence();
  assert(res == MPI_SUCCESS);
}

void setup_ofi_post_pmi(int rank, int size) {
  uint8_t value[128];
  uint32_t valuelen = 128;

  Grid g = get_3d_grid(rank, size);

  // x+
  int other = g.xplus();
  pmi_get(other, "endpoint", value, gni_address_size);
  int res = fi_av_insert(av, value, 1, &neighbors[0], 0, NULL);
  assert(res == 1);
  pmi_get(other, "mr-key-surface", &surface_keys[0], sizeof(uint64_t));
  pmi_get(other, "mr-key-counters", &counter_keys[0], sizeof(uint64_t));

  pmi_get(other, "mr-baseptr-surface", &surface_baseptrs[0], sizeof(uintptr_t));
  pmi_get(other, "mr-baseptr-counters", &counter_baseptrs[0], sizeof(uintptr_t));

  // x-
  other = g.xminus();
  pmi_get(other, "endpoint", value, gni_address_size);
  res = fi_av_insert(av, value, 1, &neighbors[1], 0, NULL);
  assert(res == 1);
  pmi_get(other, "mr-key-surface", &surface_keys[1], sizeof(uint64_t));
  pmi_get(other, "mr-key-counters", &counter_keys[1], sizeof(uint64_t));

  pmi_get(other, "mr-baseptr-surface", &surface_baseptrs[1], sizeof(uintptr_t));
  pmi_get(other, "mr-baseptr-counters", &counter_baseptrs[1], sizeof(uintptr_t));

  // y+
  other = g.yplus();
  pmi_get(other, "endpoint", value, gni_address_size);
  res = fi_av_insert(av, value, 1, &neighbors[2], 0, NULL);
  assert(res == 1);
  pmi_get(other, "mr-key-surface", &surface_keys[2], sizeof(uint64_t));
  pmi_get(other, "mr-key-counters", &counter_keys[2], sizeof(uint64_t));

  pmi_get(other, "mr-baseptr-surface", &surface_baseptrs[2], sizeof(uintptr_t));
  pmi_get(other, "mr-baseptr-counters", &counter_baseptrs[2], sizeof(uintptr_t));

  // y-
  other = g.yminus();
  pmi_get(other, "endpoint", value, gni_address_size);
  res = fi_av_insert(av, value, 1, &neighbors[3], 0, NULL);
  assert(res == 1);
  pmi_get(other, "mr-key-surface", &surface_keys[3], sizeof(uint64_t));
  pmi_get(other, "mr-key-counters", &counter_keys[3], sizeof(uint64_t));

  pmi_get(other, "mr-baseptr-surface", &surface_baseptrs[3], sizeof(uintptr_t));
  pmi_get(other, "mr-baseptr-counters", &counter_baseptrs[3], sizeof(uintptr_t));

  // z+
  other = g.zplus();
  pmi_get(other, "endpoint", value, gni_address_size);
  res = fi_av_insert(av, value, 1, &neighbors[4], 0, NULL);
  assert(res == 1);
  pmi_get(other, "mr-key-surface", &surface_keys[4], sizeof(uint64_t));
  pmi_get(other, "mr-key-counters", &counter_keys[4], sizeof(uint64_t));

  pmi_get(other, "mr-baseptr-surface", &surface_baseptrs[4], sizeof(uintptr_t));
  pmi_get(other, "mr-baseptr-counters", &counter_baseptrs[4], sizeof(uintptr_t));

  // z-
  other = g.zminus();
  pmi_get(other, "endpoint", value, gni_address_size);
  res = fi_av_insert(av, value, 1, &neighbors[5], 0, NULL);
  assert(res == 1);
  pmi_get(other, "mr-key-surface", &surface_keys[5], sizeof(uint64_t));
  pmi_get(other, "mr-key-counters", &counter_keys[5], sizeof(uint64_t));

  pmi_get(other, "mr-baseptr-surface", &surface_baseptrs[5], sizeof(uintptr_t));
  pmi_get(other, "mr-baseptr-counters", &counter_baseptrs[5], sizeof(uintptr_t));
}

int main(int argc, char **argv) {
  int spawned;
  int size;
  int rank;
  int appnum;
  struct iovec msg_iov[6];
  struct fi_msg_rma msg[6];
  struct fi_rma_iov rma_iov[6];
  double* surface_data_out[6]; // surface pointer
  double* surface_data_in[6];  // surface pointer
  double* in = nullptr;      // input grid
  double* out = nullptr;     // output grid

  parse_argv_simple(argc, argv, rank, N, iterations);

  init_pmi(&spawned, &size, &rank, &appnum);

  setup_ofi_with_cntr(false);

  setup_pmi(rank);

  setup_ofi_post_pmi(rank, size);

  int remote_offset[] = {1, 0, 3, 2, 5, 4};

    // init data
  for (int i = 0; i < 6; i++) {
    surface_data_in[i] = (double*)my_malloc(N * N * sizeof(double));
    memset(surface_data_in[i], 0, N * N * sizeof(double));
  }

  in = (double*)my_malloc((N + 2) * (N + 2) * (N + 2) * sizeof(double));
  out = (double*)my_malloc((N + 2) * (N + 2) * (N + 2) * sizeof(double));

  memset(in, 0, (N + 2) * (N + 2) * (N + 2) * sizeof(double));
  memset(out, 0, (N + 2) * (N + 2) * (N + 2) * sizeof(double));

  // data input is in the surface
  for (int i = 0; i < 6; i++)
    surface_data_out[i] = surface + i * N * N;

  auto start = std::chrono::high_resolution_clock::now();

  uint64_t txrqs = 0;

  for (uint64_t epoch = 1; epoch < iterations + 1; epoch++) {
    // 1. pack surface
    pack_surface_simple(in, surface_data_out, N); // local -> surface_data_out[]

    // 2. signal data availability
    for (int i = 0; i < 6; i++) {
      ssize_t res = fi_inject_atomic(ep, &one,
                                     1, neighbors[i],
                                     counter_baseptrs[i] + remote_offset[i]*sizeof(uint64_t),
                                     counter_keys[i],
                                     FI_UINT64, FI_SUM);
      assert(res == 0);
    }

    // 3. wait for data availability from neighbors
    // is that legal? atomic read? ordering for atomicity in MPI?
    bool have_data = false;
    while (!have_data) {
      have_data =
          (counters[0] == epoch) && (counters[1] == epoch) &&
          (counters[2] == epoch) && (counters[3] == epoch) &&
          (counters[4] == epoch) && (counters[5] == epoch);
    }


    // 4. get data from neighbors
    for (int i = 0; i < 6; i++) {
      int res = fi_write(ep, surface_data_in[i], N * N * sizeof(double),
                         fi_mr_desc(mr_surface), neighbors[i],
                         surface_baseptrs[i] + remote_offset[i] * N * N * sizeof(double),
                         surface_keys[i], NULL);
      assert(res == 0);
      txrqs++;
    }

    // 5. update local grid
    update_local_grid_simple(out, in, N);

    // 6. wait for neigbor's data
    int res = fi_cntr_wait(txcntr, txrqs, -1);
    assert(res == 0);

    // 7. signal completion
    for (int i = 0; i < 6; i++) {
      ssize_t res = fi_inject_atomic(ep, &one,
                                     1, neighbors[i],
                                     counter_baseptrs[i] + 6*sizeof(uint64_t) + remote_offset[i]*sizeof(uint64_t),
                                     counter_keys[i],
                                     FI_UINT64, FI_SUM);
      assert(res == 0);
    }

    // 8. update surface
    copy_in_neighbors_data_simple(in, surface_data_in, N);
    update_surface_simple(out, in, N);

    // 9. wait for completion
    have_data = false;
    while (!have_data) {
      have_data =
          (counters[6] == epoch) && (counters[7] == epoch) &&
          (counters[8] == epoch) && (counters[9] == epoch) &&
          (counters[10] == epoch) && (counters[11] == epoch);
    }

    // 10. swap in and out
    swap(in, out);
  }
  auto stop = std::chrono::high_resolution_clock::now();

  const std::chrono::duration<double> diff = stop - start;
  const double local_duration = diff.count();
  if (rank == 0)
    printf("time : %fs\n", local_duration);

  int res = fi_close((fid*)rxcntr);
  res = fi_close((fid*)txcntr);
  res = fi_close((fid*)av);
  res = fi_close((fid*)ep);
  res = fi_close((fid*)domain);
  res = fi_close((fid*)fabric);
  res = PMI2_Finalize();
}
