#pragma once

#include <array>

extern void copy_in_neighbors_data_shmem(double* in,
                                         const double* surface_data_in, int i,
                                         int N, const int* coords,
                                         const int* dims);

extern void parse_argv_shmem(int argc, char** argv, int rank, int& N,
                             int& iterations, int& rankspershmem);

extern void memset_shmem(double* in, int N, int* coords, int* dims);

extern void pack_surface_shmem(const double* in, double* surface_data_out,
                               int i, int N, const int* coords,
                               const int* dims);

extern void update_surface_shmem(double* out, const double* in, int i, int N,
                                 const int* coords, const int* dims);

extern void update_local_grid_shmem(double* shmem_out, const double* shmem_in,
                                    int N, const int* coords,
                                    const int* dims_smp, const int* neighbors);

extern void verify_result_shmem(double* in, int iterations, int rank, int N,
                                const int* coords, const int* dims);
