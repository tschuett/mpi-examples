/*
 use DMAPP-based version of MPI RMA
export MPICH_RMA_OVER_DMAPP=1

aprun -N24 -n 96 ./cart-pscw.x
 */

#include <assert.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>

const int iterations = 100;
const int message_size = 12*2*3*8*8; // MPI_COMPLEX
const int window_size = 16 * 6 * (12*2*3*8*8); //bytes
const int work_time = 500;
const int longer_time = 1000;
const int longer_time2 = 0;

int main(int argc, char **argv) {
  int rank, size;
  MPI_Win win;
  MPI_Group world_group, group;
  void *baseptr;
  int dims[3] = {0, 0, 0};
  int periods[3] = {1,1,1};
  int group_members[6];
  void *data[6];
  MPI_Comm comm_cart;
  int xminus, xplus, yminus, yplus, zminus, zplus, rank_source;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // init data
  for(int i = 0; i < 6; i++) {
    posix_memalign(&data[i], 4096, message_size*16);
  }

  // create topology
  int res = MPI_Dims_create(size, 3, dims);
  assert(res == 0);

  res = MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &comm_cart);
  assert(comm_cart != MPI_COMM_NULL);
  assert(res == 0);

  // find neighbors (x,y,z)
  res = MPI_Cart_shift(comm_cart, 0, -1, &rank_source, &xminus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_cart, 0, +1, &rank_source, &xplus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_cart, 1, -1, &rank_source, &yminus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_cart, 1, +1, &rank_source, &yplus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_cart, 2, -1, &rank_source, &zminus);
  assert(res == 0);
  res = MPI_Cart_shift(comm_cart, 2, +1, &rank_source, &zplus);
  assert(res == 0);

  // create window
  res = MPI_Win_allocate(window_size, 2*sizeof(double), MPI_INFO_NULL,
                         comm_cart, &baseptr, &win);
  assert(res == 0);


  // create neighbors group
  res = MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  assert(res == 0);

  group_members[0] = xminus;
  group_members[1] = xplus;
  group_members[2] = yminus;
  group_members[3] = yplus;
  group_members[4] = zminus;
  group_members[5] = zplus;
  res = MPI_Group_incl(world_group, 6, group_members, &group);
  assert(res == 0);


  for(int i = 0; i < iterations; i++) {
    // post
    res = MPI_Win_post(group, 0, win);
    assert(res == 0);

    usleep(work_time); // work

    // start
    res = MPI_Win_start(group, 0, win);
    assert(res == 0);

    // put xminus
    res = MPI_Put(data[0], message_size, MPI_DOUBLE_COMPLEX, xminus, 0, message_size,
                  MPI_DOUBLE_COMPLEX, win);
    assert(res == 0);
    // put xplus
    res = MPI_Put(data[1], message_size, MPI_DOUBLE_COMPLEX, xplus, 0, message_size,
                  MPI_DOUBLE_COMPLEX, win);
    assert(res == 0);

    usleep(work_time); // work

    // put yminus
    res = MPI_Put(data[2], message_size, MPI_DOUBLE_COMPLEX, yminus, 0, message_size,
                  MPI_DOUBLE_COMPLEX, win);
    assert(res == 0);
    // put yplus
    res = MPI_Put(data[3], message_size, MPI_DOUBLE_COMPLEX, yplus, 0, message_size,
                  MPI_DOUBLE_COMPLEX, win);
    assert(res == 0);

    usleep(work_time); // work

    // put zminus
    res = MPI_Put(data[4], message_size, MPI_DOUBLE_COMPLEX, zminus, 0, message_size,
                  MPI_DOUBLE_COMPLEX, win);
    assert(res == 0);
    // put zplus
    res = MPI_Put(data[5], message_size, MPI_DOUBLE_COMPLEX, zplus, 0, message_size,
                  MPI_DOUBLE_COMPLEX, win);
    assert(res == 0);

    usleep(longer_time); // work

    // complete
    res = MPI_Win_complete(win);
    assert(res == 0);

    usleep(longer_time2); // work

    // wait
    res = MPI_Win_wait(win);
    assert(res == 0);
  }

  // cleanup
  res = MPI_Win_free(&win);
  assert(res == 0);

  MPI_Finalize();
  return 0;
}
