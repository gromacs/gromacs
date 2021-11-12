The NxM atom-cluster non-bonded algorithm {#page_nbnxm}
=========================================

The algorithm
=============

Computing non-bonded pair interactions is the most time consuming part
of most molecular dynamics simulations. It is therefore necessary to
(highly) optimize this to achieve good simulation performance.
The standard atom pair lists are not a good match for modern SIMD
(single-instruction multiple-data) CPU architectures, nor for GPGPU
accelerators. To achieve higher (cache) data reuse and instruction
parallelism, we cluster atoms in groups of size N and M, where N and M
are either the same or differ by a factor of 2. This is done by spatial
gridding and binning. We then construct a pair list between these
clusters instead of single atoms. This not only leads to a smaller list,
but also regularizes the data access. The non-bonded pair-interaction
kernels can then compute interactions between all atom-pairs in
a cluster-pair simultaneously. For GPUs there is another layer:
superclusters consisting of multiple clusters to increase data reuse.

Architecture support
====================

Currently the module supports 5 different kernel architectures:
* Plain C++: slow, only for reference.
* SIMD 4xM: targets CPUs using SIMD intrinsics with N=4 and M=2, 4 or 8, SIMD width 2, 4 or 8.
* SIMD 2xMM: targets CPUs using SIMD intrinsics with N=4 and M=4 or 8, SIMD width 8 or 16.
* GPU: targets GPUs with N=M=8 or N=M=4, depending on
  `GMX_GPU_NB_CLUSTER_SIZE` compilation option value.
