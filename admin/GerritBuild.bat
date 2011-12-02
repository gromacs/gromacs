cmake -G "Visual Studio 10 Win64" -D GMX_DOUBLE=%GMX_DOUBLE% -D GMX_MPI=%GMX_MPI% -D GMX_OPENMP=%GMX_OPENMP% -DGMX_DEFAULT_SUFFIX=off -DCMAKE_BUILD_TYPE=Debug . && ^
msbuild && ^
ctest -D ExperimentalTest
