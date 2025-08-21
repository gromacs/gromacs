# Building GROMACS on Apple Silicon (OpenCL/Metal)

This HowTo documents a proven configuration for macOS ARM (M1/M2/M3) using OpenCL (Metal backend) with VkFFT.

## Prerequisites
- Xcode CLT
- ARM Homebrew at /opt/homebrew
- Packages: `opencl-headers`, `libomp`, `fftw`, `hwloc`, `cmake`

## Key pitfalls
- Avoid mixing `/usr/local` (x86_64) and `/opt/homebrew` (arm64)
- Use `CL_TARGET_OPENCL_VERSION=120` and silence deprecation warnings
- For OpenCL build, use `-bonded cpu -update cpu` in `mdrun`

## Toolchain usage
```
cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE=cmake/toolchains/apple-arm64-opencl.cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/opt/gromacs \
  -DGMX_MPI=OFF -DGMX_OPENMP=ON -DGMX_GPU=OpenCL -DGMX_HWLOC=ON \
  -DGMX_BUILD_OWN_FFTW=ON
cmake --build build -j12
./build/bin/gmx --version
```

## GPU smoke
```
gmx grompp -f water.mdp -c water.gro -p water.top -o water.tpr -maxwarn 2
gmx mdrun -s water.tpr -nb gpu -pme gpu -pmefft gpu -bonded cpu -update cpu -gpu_id 0
```
