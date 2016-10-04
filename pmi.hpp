#pragma once

extern void pmi_get(int pe, char *key, void *value, size_t valuelen);
extern void pmi_put(char *key, void *value, size_t valuelen, int rank);

extern void init_pmi(int *spawned, int *size, int *rank, int *appnum);
