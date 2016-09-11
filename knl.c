//#define _GNU_SOURCE
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <unistd.h>

int get_node(char* node) {
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    int found_node = 0;

    fp = fopen("/proc/cray_xt/cname", "r");
    if (fp == NULL) {
      printf("could not open /proc/cray_xt/cname\n");
      return -1;
    }

    while ((read = getline(&line, &len, fp)) != -1) {
      if(strlen(line) >= 2) {
        *node = line[strlen(line)-2];
        found_node = 1;
      }
    }

    fclose(fp);
    if (line)
        free(line);

    if(found_node)
      return 0;
    else
      return -1;
}

int knl_mode(char *numa, char *mcdram, size_t n) {
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char mcdram_pattern[16];
    char numa_pattern[16];
    char node = 0;
    int found_mcdram = 0;
    int found_numa = 0;

    int res = get_node(&node);
    if(res != 0)
      return res;

    res = snprintf(mcdram_pattern, sizeof(mcdram_pattern), "mcdram_cfg[%c", node);
    if(res < 0)
      return res;
    res = snprintf(numa_pattern, sizeof(numa_pattern), "numa_cfg[%c", node);
    if(res < 0)
      return res;

    fp = fopen("/.hwinfo.cray", "r");
    if (fp == NULL) {
      printf("could not open /.hwinfo.cray\n");
      return -1;
    }

    while ((read = getline(&line, &len, fp)) != -1) {
      if(strstr(line, mcdram_pattern) == line) {
        char *s = strchr(line, '=');
        if(s != NULL) {
          strncpy(mcdram, s + 1, MIN(n, strlen(s+1) - 1));
          found_mcdram = 1;
        }
      }
      if(strstr(line, numa_pattern) == line) {
        char *s = strchr(line, '=');
        if(s != NULL) {
          strncpy(numa, s + 1, MIN(n, strlen(s+1) - 1));
          found_numa = 1;
        }
      }
    }

    fclose(fp);
    if (line)
        free(line);

    if(found_numa && found_mcdram)
      return 0;
    else
      return -1;
}

#ifdef HAVE_CRAY
void print_mode() {
  char numa[32];
  char mcdram[32];
  char hostname[32];

  memset(numa, 0, sizeof(numa));
  memset(mcdram, 0, sizeof(mcdram));

  if(access("/proc/cray_xt/cname", F_OK) == 0 && access("/.hwinfo.cray", F_OK) == 0) {  
    int res0 = knl_mode(numa, mcdram, sizeof(numa));
    if(res != 0)
      return 1;

    int res1 = gethostname(hostname, sizeof(hostname));
    if(res0 == 0 && res1 == 0)
      printf("%s: NUMA: %s     MCDRAM: %s\n", hostname, numa, mcdram);
  }
}

void check_mode() {
}
#else
void print_mode() {
}

void check_mode() {
}
#endif
