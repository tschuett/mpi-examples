#include "my-malloc.hpp"

#include "config.h"

#ifdef HAVE_MEMKIND
#include <hbwmalloc.h>
#endif

#include <cassert>

#ifdef HAVE_MEMKIND
#else
void *my_malloc(size_t bytes) {
  void *buffer = nullptr;
  int res = posix_memalign(&buffer, 4096, bytes);
  assert(res == 0);

  return buffer;
}
#endif

#ifdef HAVE_MEMKIND
#else
void my_free(void *buffer) {
  free(buffer);
}
#endif
