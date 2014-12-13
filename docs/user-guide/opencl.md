TABLE OF CONTENTS
1. SUPPORTED OPENCL DEVICES
2. OPENCL SETUP
3. KNOWN LIMITATIONS
4. TESTED CONFIGURATIONS

1. SUPPORTED OPENCL DEVICES
   ========================
The current version works with NVIDIA GPUs and GCN based AMD GPUs.
Make sure that you have the latest drivers installed.
The minimum OpenCL version required is 1.1
Also check "Known Limitations" chapter.

2. OPENCL SETUP
   ============
Build Gromacs with OpenCL support enabled
-----------------------------------------
To build Gromacs with OpenCL support enabled, an OpenCL SDK must be installed
and the following cmake flags must be set:
	GMX_GPU
	GMX_USE_OPENCL
After setting these flags, if an OpenCL implementation has been detected on
your system, the following cmake entries will be defined:
	OPENCL_INCLUDE_DIR - the OpenCL include directory
	OPENCL_LIBRARY - the OpenCL library directory
The two paths will be automatically used to configure the project. 

Run Gromacs with OpenCL accelerations enabled
---------------------------------------------
Gromacs loads and builds at runtime the OpenCL kernels. To do so, it needs to
know the location of the OpenCL source files.
If you want to run the installed version, the path to the OpenCL files is
automatically defined.
If you do not wish to install Gromacs, but run the version built from sources,
you need to provide the path to the source tree with the OpenCL kernels like
below:
	export GMX_OCL_FILE_PATH=/path-to-gromacs/src/

Caching options for building the kernels
----------------------------------------
Building an OpenCL program can take a significant amount of time. NVIDIA
implements a mechanism to cache the result of the build. As a consequence,
only the first build of the OpenCL kernels will take longer, the following
builds will be very fast. AMD drivers, on the other hand, implement no
caching and building a program can be very slow. That's why we have started
implementing our own caching. Caching for OpenCL kernel builds is by default
enabled. To disable it, set GMX_OCL_NOGENCACHE environment variable.

 If you plan to modify the OpenCL kernels, you should disable any caching:
 * add GMX_OCL_NOGENCACHE environment variable and set it to 1
 * for NVIDIA cards: add CUDA_CACHE_DISABLE environment variable and set it to 1
 
OpenCL Device Selection
-----------------------
The same option used to select CUDA devices or to define a mapping of GPUs to
PP ranks can also be used for OpenCL devices: -gpu_id

Environment Variables For OpenCL
--------------------------------
Currently, several environment variables exist that help customize some aspects
of the OpenCL version of Gromacs. They are mostly related to the runtime
compilation of OpenCL kernels, but they are also used on the device selection.

   GXM_OCL_FILE_PATH: Is the full path to Gromacs src folder. Useful when gmx
   is called from a folder other than the installation/bin folder.
   
   GMX_OCL_NOGENCACHE: Disable caching for OpenCL kernel builds.
   
   GMX_OCL_NOFASTGEN: Generates and compiles all algorithm flavors, otherwise
   only the flavor required for the simulation is generated and compiled.
   
   GMX_OCL_FASTMATH: Adds the option cl-fast-relaxed-math to the compiler
   options (in the CUDA version this is enabled by default, it is likely that
   the same will happen with the OpenCL version soon)
   
   GMX_OCL_DUMP_LOG: If defined, the OpenCL build log is always written to file.
   The file is saved in the current directory with the name
   OpenCL_kernel_file_name.build_status where OpenCL_kernel_file_name is the name
   of the file containing the OpenCL source code (usually nbnxn_ocl_kernels.cl)
   and build_status can be either SUCCEEDED or FAILED. If this environment
   variable is not defined, the default behavior is the following:
      - Debug build: build log is always written to file
	  - Release build: build log is written to file only in case of errors.
   
   GMX_OCL_VERBOSE: If defined, it enables verbose mode for OpenCL kernel build.
   Currently available only for NVIDIA GPUs. See GMX_OCL_DUMP_LOG for details
   about how to obtain the OpenCL build log.
   
   GMX_OCL_DUMP_INTERM_FILES: If defined, intermediate language code corresponding
   to the OpenCL build process is saved to file. Caching has to be turned off in
   order for this option to take effect (see GMX_OCL_NOGENCACHE).
      - NVIDIA GPUs: PTX code is saved in the current directory with the name
	  device_name.ptx
	  - AMD GPUs: .IL/.ISA files will be created for each OpenCL kernel built.
	  For details about where these files are created check AMD documentation
	  for -save-temps compiler option.
   
   GMX_OCL_DEBUG: Use in conjunction with OCL_FORCE_CPU or with an AMD device.
   It adds the debug flag to the compiler options (-g).
   
   GMX_OCL_NOOPT: Disable optimisations. Adds the option cl-opt-disable to the
   compiler options.
   
   GMX_OCL_FORCE_CPU: Force the selection of a CPU device instead of a GPU.
   This exists only for debugging purposes. Do not expect Gromacs to function
   properly with this option on, it is solely for the simplicity of stepping
   in a kernel and see what is happening.
   
   GMX_OCL_NB_ANA_EWALD: Forces the use of analytical Ewald kernels.
   Equivalent of CUDA env var GMX_CUDA_NB_ANA_EWALD
   
   GMX_OCL_NB_TAB_EWALD: Forces the use of tabulated Ewald kernel. Equivalent
   of CUDA env var GMX_OCL_NB_TAB_EWALD
   
   GMX_OCL_NB_EWALD_TWINCUT: Forces the use of twin-range cutoff kernel.
   Equivalent of CUDA env var GMX_CUDA_NB_EWALD_TWINCUT
   
   GMX_DISABLE_OCL_TIMING: Disables timing for OpenCL operations

3. KNOWN LIMITATIONS
   =================
- Intel CPUs are not supported
- Intel GPUs are not supported
- The current implementation is not compatible with OpenCL devices that are
  not using warp/wavefronts or for which the warp/wavefront size is not a
  multiple of 32
- The following kernels are known to produce incorrect results:
	nbnxn_kernel_ElecEwQSTab_VdwLJ_VF_prune_opencl
	nbnxn_kernel_ElecEwQSTab_VdwLJ_F_prune_opencl
- Due to blocking behavior of clEnqueue functions in NVIDIA driver, there is
  almost no performance gain when using NVIDIA GPUs. A bug report has already
  been filled on about this issue. A possible workaround would be to have a
  separate thread for issuing GPU commands. However this hasn't been implemented
  yet.	

4. TESTED CONFIGURATIONS
   =====================
Tested devices:
	NVIDIA GPUs: GeForce GTX 660M, GeForce GTX 750Ti, GeForce GTX 780
	AMD GPUs: FirePro W5100, HD 7950, FirePro W9100, Radeon R7 M260, R9 290
	
Tested kernels:
Kernel                                          |Benchmark test                                 |Remarks
--------------------------------------------------------------------------------------------------------
nbnxn_kernel_ElecCut_VdwLJ_VF_prune_opencl      |d.poly-ch2                                     |
nbnxn_kernel_ElecCut_VdwLJ_F_opencl             |d.poly-ch2                                     |
nbnxn_kernel_ElecCut_VdwLJ_F_prune_opencl       |d.poly-ch2                                     |
nbnxn_kernel_ElecCut_VdwLJ_VF_opencl            |d.poly-ch2                                     |
nbnxn_kernel_ElecRF_VdwLJ_VF_prune_opencl       |adh_cubic with rf_verlet.mdp                   |
nbnxn_kernel_ElecRF_VdwLJ_F_opencl              |adh_cubic with rf_verlet.mdp                   |
nbnxn_kernel_ElecRF_VdwLJ_F_prune_opencl        |adh_cubic with rf_verlet.mdp                   |
nbnxn_kernel_ElecEwQSTab_VdwLJ_VF_prune_opencl  |adh_cubic_vsites with pme_verlet_vsites.mdp    |Failed
nbnxn_kernel_ElecEwQSTab_VdwLJ_F_prune_opencl   |adh_cubic_vsites with pme_verlet_vsites.mdp    |Failed
nbnxn_kernel_ElecEw_VdwLJ_VF_prune_opencl       |adh_cubic_vsites with pme_verlet_vsites.mdp	|
nbnxn_kernel_ElecEw_VdwLJ_F_opencl              |adh_cubic_vsites with pme_verlet_vsites.mdp	|
nbnxn_kernel_ElecEw_VdwLJ_F_prune_opencl        |adh_cubic_vsites with pme_verlet_vsites.mdp	|
nbnxn_kernel_ElecEwTwinCut_VdwLJ_F_prune_opencl	|adh_cubic_vsites with pme_verlet_vsites.mdp	|
nbnxn_kernel_ElecEwTwinCut_VdwLJ_F_opencl       |adh_cubic_vsites with pme_verlet_vsites.mdp    |

Input data used for testing - Benchmark data sets available here:
ftp://ftp.gromacs.org/pub/benchmarks