cd ..
make -j &&
cd test/cha
mpirun -n 8 cans &&
cd ..
