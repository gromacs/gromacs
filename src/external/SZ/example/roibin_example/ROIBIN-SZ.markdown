# Using ROIBIN-SZ

## Installing ROIBIN-SZ

ROIBIN-SZ is included in SZ after release 2.1.11.1.

The easiest and recommended way to use ROIBIN-SZ is to use libpressio which provides additional error checking and ease of use above the SZ C interface.

They can both can be installed easily via [spack](https://github.com/spack/spack) on Linux and MacOS:

```bash
git clone https://github.com/spack/spack
git clone https://github.com/robertu94/spack_packages robertu94_packages # only for libpressio
source ./spack/share/spack/setup-env.sh
spack compiler find
spack repo add ./robertu94_packages

spack install libpressio^sz
```

## Using ROIBIN-SZ


Examples for both can be found in `example_lp.c` and `example_sz.c` respectively.
These files will need to be customized to load data for your specific detector.

Here are a few of the important configuration parameters:

+ `bin_size <uint32>` for ROIBIN-SZ, the size of the binning applied
+ `calib_panel <uint8_t[]>` for ROIBIN-SZ the size of the calibration panel.  Expects a array of the same dimensions as the input with a 1 to indicate that an entry should be ignored.
+ `num_peaks <uint32>` for ROIBIN-SZ the number of peaks.
+ `peak_size <uint32>` for ROIBIN-SZ the size of the region of interest around the peaks.  Must be odd.
+ `peaks_cols <uint16[]>` for ROIBIN-SZ the list of columns peaks appear in
+ `peaks_rows <uint16[]>` for ROIBIN-SZ the list of rows peaks appear in
+ `peaks_segs <uint16[]>` for ROIBIN-SZ the segments peaks appear in
+ `sz_dim <uint32>` for ROIBIN-SZ the SZ dimensionality prefered. 3 is recommended
+ `tolerance <double>` for ROIBIN-SZ the tolerance used after binning

