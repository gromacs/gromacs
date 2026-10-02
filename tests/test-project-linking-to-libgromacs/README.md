Simple CMake project to test that linking to an installed
GROMACS version works.

Usage
-----

```
# configure (to check GROMACS is found)
export gmx_install=<path-to-installed-gromacs>
cmake -S . -B build-dir -DCMAKE_PREFIX_PATH=$gmx_install

# build (to check compilation succeeds)
cmake --build build-dir

# run (to check dynamic linking succeeds)
build-dir/test-gromacs
```
