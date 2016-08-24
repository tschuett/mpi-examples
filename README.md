6-point stencil on a 3 dimensional grid. Each rank has a local grid of N^3.

module load gcc/6.1.0 ?
./configure

make

msub bench.pbs

aprun -n24 ./cart-immediate-msg.x 100
