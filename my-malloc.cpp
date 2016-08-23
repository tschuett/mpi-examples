#include "my-malloc.hpp"

#include "config.h"

#ifdef HAVE_MEMKIND
#include <hbwmalloc.h>
#include <errno.h>
#endif

#include <cassert>

#ifdef HAVE_MEMKIND
void *my_malloc(size_t bytes) {
  void *buffer = nullptr;
  if(hbw_check_available()) {
    int res = hbw_posix_memalign_psize(&buffer, 2*1024*1024ULL, bytes,
				       HBW_PAGESIZE_2MB);
    if(res == 0)
      return buffer;
    else if(res == ENOMEM) {
      res = hbw_posix_memalign_psize(&buffer, 2*1024*1024ULL, bytes,
				     HBW_PAGESIZE_4KB);
      assert(res == 0);
      return buffer;
    }
  } else {
    int res = posix_memalign(&buffer, 4096, bytes);
    assert(res == 0);

    return buffer;
  } 
}
#else
void *my_malloc(size_t bytes) {
  void *buffer = nullptr;
  int res = posix_memalign(&buffer, 4096, bytes);
  assert(res == 0);

  return buffer;
}
#endif

#ifdef HAVE_MEMKIND
void my_free(void *buffer) {
  hbw_free(buffer);
}
#else
void my_free(void *buffer) {
  free(buffer);
}
#endif
