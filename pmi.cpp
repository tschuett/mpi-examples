#include <pmi2.h>
#include <mpi.h>

#include <string.h>

#include <cassert>
#include <cstdio>
#include <cstdlib>

#ifndef PMI2_SUCCESS
#define PMI2_SUCCESS 0
#endif

char *kvs_name, *kvs_key, *kvs_value;
int max_name_len, max_key_len, max_val_len;

int encode(const void *inval, int invallen, char *outval, int outvallen) {
    static unsigned char encodings[] = {
        '0','1','2','3','4','5','6','7', \
        '8','9','a','b','c','d','e','f' };
    int i;

    if (invallen * 2 + 1 > outvallen) {
        return 1;
    }

    for (i = 0; i < invallen; i++) {
        outval[2 * i] = encodings[((unsigned char *)inval)[i] & 0xf];
        outval[2 * i + 1] = encodings[((unsigned char *)inval)[i] >> 4];
    }

    outval[invallen * 2] = '\0';

    return 0;
}

int decode(const char *inval, void *outval, int outvallen) {
    int i;
    char *ret = (char*) outval;

    if (outvallen != strlen(inval) / 2) {
        return 1;
    }

    for (i = 0 ; i < outvallen ; ++i) {
        if (*inval >= '0' && *inval <= '9') {
            ret[i] = *inval - '0';
        } else {
            ret[i] = *inval - 'a' + 10;
        }
        inval++;
        if (*inval >= '0' && *inval <= '9') {
            ret[i] |= ((*inval - '0') << 4);
        } else {
            ret[i] |= ((*inval - 'a' + 10) << 4);
        }
        inval++;
    }

    return 0;
}

void pmi_put(char *key, void *value, size_t valuelen, int rank) {
  int res = snprintf(kvs_key, max_key_len, "ofi-%lu-%s", (long unsigned) rank, key);
  assert(res > 0);
  res = encode(value, valuelen, kvs_value, max_val_len);
  assert(res == 0);
  res = PMI2_KVS_Put(kvs_key, kvs_value);
  assert(res == PMI2_SUCCESS);
}

void pmi_get(int pe, char *key, void *value, size_t valuelen) {
    int len;

    snprintf(kvs_key, max_key_len, "ofi-%lu-%s", (long unsigned) pe, key);
    int res = PMI2_KVS_Get(kvs_name, PMI2_ID_NULL,
                           kvs_key, kvs_value, max_val_len, &len);
    assert(res == PMI2_SUCCESS);
    res = decode(kvs_value, value, valuelen);
    assert(res == 0);
}


void init_pmi(int *spawned, int *size, int *rank, int *appnum) {
  int res = PMI2_Init(spawned, size, rank, appnum);
  assert(res == MPI_SUCCESS);

  max_name_len = PMI2_MAX_VALLEN;
  kvs_name = (char*) malloc(max_name_len);
  assert(kvs_name != NULL);

  max_key_len = PMI2_MAX_KEYLEN;
  kvs_key = (char*) malloc(max_key_len);
  assert(kvs_key != NULL);

  max_val_len = PMI2_MAX_VALLEN;
  kvs_value = (char*) malloc(max_val_len);
  assert(kvs_value != NULL);

  res = PMI2_Job_GetId(kvs_name, max_name_len);
  assert(res == MPI_SUCCESS);
}
