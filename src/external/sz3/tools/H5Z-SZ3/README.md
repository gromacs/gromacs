# H5Z-SZ3 Filter Integration

The H5Z-SZ3 filter integrates the SZ3 compression library with HDF5, providing an efficient way to compress and decompress data within HDF5 files.

Use the filter from SZ3 3.4.0 or later.

## Table of Contents
- [Installation](#installation)
- [H5Z-SZ3 cd_values](#h5z-sz3-cd_values)
- [Usage](#usage)
  - [HDF5 Executables](#hdf5-executables)
  - [Python (h5py)](#python-h5py)
  - [C/C++](#cc)
- [Using the filter in an application](#using-the-filter-in-an-application)
- [Rules for the filter's code](#rules-for-the-filters-code)

## Installation

### Step 1: Build and Install SZ3
Compile and install SZ3 with the H5Z-SZ3 filter enabled:
```bash
cmake -S . -B build -DBUILD_H5Z_FILTER=ON -DCMAKE_INSTALL_PREFIX=<PREFIX>
cmake --build build -j && cmake --install build
```
A filter built against HDF5 1.14.5 or newer needs HDF5 1.14.5 or newer at run time.

### Step 2: Configure Environment
Installing puts a copy of the filter in `<prefix>/lib/plugin` (`H5Z_SZ3_PLUGIN_INSTALL_DIR`; empty for no copy), a
directory that holds nothing else. Point `HDF5_PLUGIN_PATH` at that, not at `<prefix>/lib` or
`<prefix>/bin`: HDF5 dlopens every `lib*.so` in each directory on the path, and every `*.dll` on
Windows. The variable also replaces HDF5's own compiled-in plugin directory rather than adding to
it, so list that too if other filters are in use.
```bash
export HDF5_PLUGIN_PATH=<PREFIX>/lib/plugin
```
An application that ships the plugin with itself can call `H5PLprepend("<its plugin directory>")` at
startup instead and set nothing; one that links `SZ3::hdf5sz3` can call
`H5Zregister(H5PLget_plugin_info())` (before any `H5Zfilter_avail`, which may load another SZ3 plugin).
On Windows it still has to find the library: the install puts `hdf5sz3.dll` (`libhdf5sz3.dll` under
MinGW) in `<prefix>/bin` and only the import library in `<prefix>/lib`, and there is no RPATH to
record either, so `<prefix>/bin` has to be on `PATH` before the application runs.

## H5Z-SZ3 cd_values
* HDF5 restricts the parameters that can be passed to filters through an integers array called `cd_values`.
* H5Z-SZ3 uses `cd_values` to pass the desired compression settings (e.g., algorithm, error bounds) to the compression process. `cd_values[0]` is the SZ3 data version, `(major << 24) | (minor << 16) | (patch << 8)` (`0x03040000` for 3.4.0), followed by the `Config` object serialized with `save()`.
* `cd_values` are read only when compressing, including appends; decompression reads the configuration from the compressed data.

## Usage

### HDF5 Executables

Note: if the HDF5 in your system contains both `h5repack-shared` and `h5repack`, you need to use `h5repack-shared` for filters, as outlined in [H5Z-ZFP Issue #137](https://github.com/LLNL/H5Z-ZFP/issues/137).

#### Compression

**Method 1: Compression with default settings (see SZ3/utils/Config.hpp for defaults)**
```bash
h5repack-shared -f UD=32024,0 data.h5 data.sz3.h5
```

**Method 2: Compression with customized settings in configuration file**

Step 1: Generate H5Z-SZ3 cd_values from SZ3 configuration file:
```bash
cdvalueHelper -c sz3.config
```
Step 2: Compress by h5repack-shared with generated cd_values:
```bash
h5repack-shared -f UD=[put cd_values here] data.h5 data.sz3.h5
```

#### Decompression
```bash
h5repack-shared -f NONE data.sz3.h5 data.sz3_decompressed.h5
```

### Python (h5py)
To use H5Z-SZ3 with h5py in Python:

1. Ensure h5py is linked against the same HDF5 library used to build the SZ3 filter (typically system HDF5).
2. Set `HDF5_PLUGIN_PATH` to the directory containing `libhdf5sz3.so` or `libhdf5sz3.dylib` **before importing h5py**.
3. Use h5py's compression interface with SZ3 filter ID 32024 and appropriate compression options.

```python
import os
os.environ["HDF5_PLUGIN_PATH"] = "/path/to/prefix/lib/plugin"
import h5py
import numpy as np

from cdvalueHelper import SZ3
config = SZ3('ALGO_INTERP_LORENZO', absolute=1e-3)

with h5py.File('data.h5', 'w') as f:
    f.create_dataset('dataset', data=np.random.rand(100, 100), 
                     compression=32024, compression_opts=tuple(config.cd_values))
```

**Note on HDF5 Runtime Compatibility**: If h5py was installed with its bundled HDF5 (common on macOS), it may not load plugins built against system HDF5 due to library conflicts. To resolve:
- Rebuild h5py against system HDF5: `export HDF5_DIR=/path/to/system/hdf5 && pip install --force-reinstall --no-binary h5py h5py`
- Ensure HDF5 versions match between h5py and the plugin.

### C/C++
A C program sets the filter on a dataset creation property list with `H5Pset_sz3`, declared in `H5Z_SZ3.hpp`. A CMake project links `SZ3::hdf5sz3` and requires it with `find_package(SZ3 COMPONENTS hdf5sz3)`. See the examples `sz3ToHDF5.cpp` and `dsz3FromHDF5.cpp`.

## Using the filter in an application

### Setting and registering it
- Set the filter with `H5Pset_sz3(plist, cmprAlgo, errorBoundMode, abs, rel, psnr, l2norm)`, not by writing `Config::save()` into `cd_values` yourself.
- An application that links the filter registers it with `H5Zregister(H5PLget_plugin_info())`, before any `H5Zfilter_avail` (which may load another SZ3 plugin), and then needs no `HDF5_PLUGIN_PATH`. One that ships the plugin calls `H5PLprepend(dir)`.
- Register when the application opens any file, not only when it creates an SZ3 dataset, so that a program that only reads can read SZ3 data. A failed registration should not stop files without SZ3 data from opening.
- Registering replaces, for the whole process, an SZ3 filter loaded from `HDF5_PLUGIN_PATH`, also a newer one. State that choice where it is made.
- `h5repack` exits 0 and writes an uncompressed copy when it cannot load the filter. Check the output, for example with `h5dump -pH`, which shows the filter and its version.

### Datatypes, chunks and other filters
- The dataset's type must be the host's float, double or 1, 2, 4 or 8-byte integer, in host byte order; other types are refused when the dataset is created. Chunks with more than 4 dimensions longer than 1 are refused too.
- For molecular-dynamics coordinates use `H5Z_SZ3_ALGO_BIOMD` (better with several frames per chunk). A `H5Z_SZ3_ALGO_BIOMD` chunk must have 3 as its last dimension, `(frames, atoms, 3)`; creating the dataset fails otherwise, and the error names the chunk shape. Set the chunk shape yourself (in h5py `chunks=(k, atoms, 3)`): h5py's automatic chunking splits the atoms and xyz. See [docs/molecular-dynamics.md](../../docs/molecular-dynamics.md).
- Do not put a filter that rearranges bytes, such as shuffle, before SZ3: SZ3 then compresses the rearranged bytes as values. The filter refuses a chunk whose size an earlier filter changed.
- SZ3 has no checksum; add `H5Pset_fletcher32` after `H5Pset_sz3` to detect damaged chunks.
- The filter keeps the flags it was set with. With `H5Z_FLAG_OPTIONAL` (h5py sets it), HDF5 stores a chunk raw when SZ3 fails on it.
- The filter does not use OpenMP.

### Writing a chunk again and appending
- A chunk that is filtered and then written again is decompressed and compressed once more, and its error can then exceed the bound, except with `H5Z_SZ3_ALGO_NOPRED` with an absolute bound, and with `H5Z_SZ3_ALGO_BIOMD` with an absolute bound for float coordinates within B, 2^21 to 2^22 times the bound (2048 nm at 5e-4 nm, 256 nm at 1e-4 nm, 128 nm at 5e-5 nm, 32 nm at 1e-5 nm). A chunk with a coordinate past B, or with NaN or Inf, is stored losslessly, at a ratio of about 1.1 to 1.5 instead of 3 to 6; that also keeps a rewritten chunk within the bound.
- When appending frames, `H5D_CHUNK_DONT_FILTER_PARTIAL_CHUNKS` alone does not stop recompression: once the dataset is extended past a chunk, HDF5 filters it, and the next write decompresses and recompresses it unless the chunk cache holds the whole chunk. Set `H5D_CHUNK_DONT_FILTER_PARTIAL_CHUNKS` and give the dataset a chunk cache that holds a whole chunk (`H5Pset_chunk_cache`), also when the file is reopened. The last, partial chunk is then stored uncompressed at full chunk size.
- A test that writes and reads in the same open file reads HDF5's chunk cache and never decompresses. Close and reopen the file, and compare within the bound. To test that reading registers the filter, unregister it (`H5Zunregister`) first.

## Rules for the filter's code
- `cd_values[0]` is `versionInt(SZ3_DATA_VER)`, followed by the `Config` bytes, little-endian. Tools that write `cd_values` (`cdvalueHelper.py`) follow the same layout.
- `H5Z_SZ3.hpp` is the only public header; everything outside its `__cplusplus` block stays plain C (CI compiles it as C99 with `-pedantic-errors -Werror`). The `H5Z_SZ3_*` constants are `static_assert`ed against the C++ enums.
- To read `cd_values` back, ask HDF5 for the count, allocate, then retrieve.
- Do not use `H5Zfilter_avail()` to decide anything about a property list: it answers for the library. Look the filter up on the list.
- `set_local` refuses `cd_values` it cannot read, rather than continuing with defaults.
- No exception may leave `H5Z_filter_sz3` or `set_local`: catch, push the message with `H5E_ERR_CLS` (`e.what()` through `"%s"`), and return 0 (filter) or a negative value (`set_local`). Print to stderr only when HDF5 has paused the error stack, never to stdout.
- The buffer returned to HDF5 comes from `malloc` and is freed if compression throws.
- The filter links `SZ3core` and does not use OpenMP. With `SZ3_DEBUG_TIMINGS` it still compiles every datatype `set_local` accepts. Test against HDF5 1.10, 1.14 and 2.x.
