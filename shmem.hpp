#pragma once

extern void copy_in_neighbors_data_shmem(double* in, const double* shmem_in,
                                         int rank_shmem, int i, int N);

extern void parse_argv_shmem(int argc, char** argv, int rank, int& N,
                             int& iterations, int& rankspershmem);
