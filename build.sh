rm lib/libgromacs.so.3.0.0
ninja
mv lib/libgromacs.so.3.0.0 lib/libgromacs_SSE2.so.3.0.0
/usr/local/gcc-7.2.0/bin/g++ -shared ../src/dispatch.cc -o lib/libgromacs.so.3.0.0 -fPIC -I ../src/ ../src/gromacs/hardware/cpuinfo.cpp -Isrc

