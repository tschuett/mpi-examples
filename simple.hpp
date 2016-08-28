#pragma once

extern void pack_surface_simple(const double* in, double* surface_data_out[6],
                                int N);

extern void pack_surface_simple(const double* in, double* surface_data_out,
                                int i, int N);

extern void update_local_grid_simple(double* out, const double* in, int N);

extern void copy_in_neighbors_data_simple(double* in,
                                          double* surface_data_in[6], int N);
extern void copy_in_neighbors_data_simple(double* in, double* surface_data_in,
                                          int i, int N);

extern void update_surface_simple(double* out, const double* in, int N);
extern void update_surface_simple(double* out, const double* in, int i, int N);

extern void verify_result_simple(double* in, int iterations, int rank, int N);

extern void parse_argv_simple(int argc, char** argv, int rank, int& N,
                              int& iterations);
