SZ3: A Modular Error-bounded Lossy Compression Framework for Scientific Datasets
=====

SZ3 compresses floating-point and integer arrays from simulations and instruments, and guarantees that every
decompressed value differs from the original by no more than the error bound you set. It is a header-only C++17
library. It can also be used from C, Python, the `sz3` command line, and HDF5 through a filter.

## Installation

With Homebrew (macOS and Linux): `brew install szcompressor/tap/sz3`.

To build from source you need a C++17 compiler and CMake 3.19 or newer. Zstd and OpenMP are used when found.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
cmake --build build -j
cmake --install build
```
The tools go to `<INSTALL_DIR>/bin` and the headers to `<INSTALL_DIR>/include`.

## Quick start

### Command line

```bash
# Compress 8x8x128 float data (dimensions fastest-varying first) with a relative error bound of 1e-3
sz3 -f -i tools/sz3/testfloat_8_8_128.dat -z test.sz -3 8 8 128 -M REL 1e-3
# Decompress to test.out; -a compares with the original (-i) and prints the maximum error and compression ratio
sz3 -f -z test.sz -o test.out -i tools/sz3/testfloat_8_8_128.dat -a
```
Run `sz3 -h` for every option; the common ones:

| Option | Meaning |
|---|---|
| `-f`, `-d`, `-I 32`, `-I 64` | data type: float, double, int32, int64 |
| `-i <file>` | original data, raw binary |
| `-z <file>` | compressed file: written when compressing, read when decompressing |
| `-o <file>` | decompressed data, raw binary (`-t` writes text) |
| `-1 nx`, `-2 nx ny`, `-3 nx ny nz`, `-4 nx ny nz np` | dimensions, fastest-varying first: `-3 nx ny nz` is `data[nz][ny][nx]` |
| `-M <mode> <bound>` | error-bound mode and bound: `ABS`, `REL`, `PSNR`, `NORM`; `ABS_AND_REL` and `ABS_OR_REL` take `-A <abs> -R <rel>` |
| `-c <file>` | configuration file, such as [tools/sz3/sz3.config](tools/sz3/sz3.config); choose the algorithm here |
| `-a` | after decompression, print error statistics (needs `-i`) |
| `-p` | print the configuration of the compressed data |
| `-v` | print the SZ3 version and the data-format version |

### C++

```cpp
#include <SZ3/api/sz.hpp>
#include <vector>

int main() {
    std::vector<float> data(100 * 200 * 300, 1.0f);
    SZ3::Config conf(100, 200, 300);  // 300 is the fastest-varying dimension
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = 1e-3;

    size_t cmpSize;
    char *cmpData = SZ_compress(conf, data.data(), cmpSize);

    SZ3::Config decConf;  // filled from the compressed data
    float *decData = SZ_decompress<float>(decConf, cmpData, cmpSize);

    delete[] cmpData;
    delete[] decData;
}
```

In a CMake project, add `<INSTALL_DIR>` to `CMAKE_PREFIX_PATH`, then `find_package(SZ3)` and link `SZ3::SZ3`.

## Interfaces

| Interface | How to use / where | Maintained by |
|---|---|---|
| C++ | `#include <SZ3/api/sz.hpp>`; see [above](#c) | SZ3 |
| C | [tools/sz3c/include/sz3c.h](tools/sz3c/include/sz3c.h), library `SZ3c`; SZ2-compatible functions | SZ3 |
| Python | `pip install pysz`; [tools/pysz](tools/pysz/README.md) | SZ3 |
| Command line | `sz3`; see [above](#command-line) | SZ3 |
| HDF5 filter | H5Z-SZ3, filter ID 32024; [tools/H5Z-SZ3](tools/H5Z-SZ3/README.md) | SZ3 |
| ParaView | SZ3Reader plugin; [tools/paraview](tools/paraview/README.md) | SZ3 |
| Fortran | [ofmla/sz3_simple_example](https://github.com/ofmla/sz3_simple_example) | [Oscar Mojica](https://github.com/ofmla) |
| Rust | [sz3-rs](https://github.com/apertus-open-source-cinema/sz3-rs) | [Juniper Tyree](https://github.com/juntyr) and [Robin Heinemann](https://github.com/rroohhh) |
| Python numcodecs | [numcodecs-rs codecs/sz3](https://github.com/juntyr/numcodecs-rs/blob/main/codecs/sz3/) | [Juniper Tyree](https://github.com/juntyr) |

## Further documentation

* [docs/advanced-guide.md](docs/advanced-guide.md): algorithms and error-bound modes, build options, CMake and
  packaging, the C++ API in depth, OpenMP, floating point, and the rules for contributors.
* [docs/molecular-dynamics.md](docs/molecular-dynamics.md): molecular-dynamics data and GROMACS.
* [CHANGELOG.md](CHANGELOG.md): version history and which data formats each version reads.

## Citing SZ3

[//]: # (**Kindly note**: If you mention SZ3 in your paper, the most appropriate citation is to include these three references &#40;**TBD22, ICDE21, Bigdata18**&#41; because they cover the design and implementation of the latest version of SZ.)
* QOZv2 (the enhanced interpolation-based algorithm): [High-performance Effective Scientific Error-bounded Lossy Compression with Auto-tuned Multi-component Interpolation](https://dl.acm.org/doi/10.1145/3639259).
* SZ3's interpolation-based algorithm: [Optimizing Error-Bounded Lossy Compression for Scientiﬁc Data by Dynamic Spline Interpolation](https://ieeexplore.ieee.org/document/9458791).
* The software engineering design of SZ3: [SZ3: A modular framework for composing prediction-based error-bounded lossy compressors](https://ieeexplore.ieee.org/abstract/document/9866018).

## License and contact

SZ3 is released under a BSD license; see [copyright-and-BSD-license.txt](copyright-and-BSD-license.txt).
The vendored Zstd in `tools/zstd` keeps its own license.

* Lead developer and maintainer: Kai Zhao
* Contributors: Robert Underwood, Xin Liang, Jinyang Liu, Sheng Di, and [everyone else on GitHub](https://github.com/szcompressor/SZ3/graphs/contributors)
* SZ project lead: Franck Cappello

Report bugs and ask questions in [GitHub issues](https://github.com/szcompressor/SZ3/issues).
