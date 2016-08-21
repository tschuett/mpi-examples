6-point stencil on a 3 dimensional grid. Each rank has a local grid of N^3.

./configure

make

msub bench.pbs

aprun -n24 ./cart-immediate-msg.x 100|200
