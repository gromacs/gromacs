Performance improvements
^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

Up to a factor 2.5 speed-up of the non-bonded free-energy kernel
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The non-bonded free-energy kernel is a factor 2.5 faster with non-zero A and B
states and a factor 1.5 with one zero state. This especially improves the run
performance when non-perturbed non-bondeds are offloaded to a GPU. In that case
the PME-mesh calculation now always takes the most CPU time.


Proper dihedrals of Fourier type and improper dihedrals of preriodic type are SIMD accelerated
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Avoid configuring the own-FFTW with AVX512 enabled when |Gromacs| does not use AVX512
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Previously if |Gromacs| was configured to use any AVX flavor, the internally built FFTW
would be configured to also contain AVX512 kernels. This could cause performance loss
if the (often noisy) FFTW auto-tuner picks an AVX512 kernel in a run that otherwise 
only uses AVX/AVX2 which could run at higher CPU clocks without AVX512 clock speed limitation.
Now AVX512 is only used for the internal FFTW if |Gromacs| is also configured with
the same SIMD flavor.

Update and constraints can run on a GPU
"""""""""""""""""""""""""""""""""""""""

For standard simulations (see the user guide for more details),
update and constraints can be offloaded to a GPU with CUDA. Thus all compute
intensive parts of a simulation can be offloaded, which provides
better performance when using a fast GPU combined with a slow CPU.
By default, update will run on the CPU, to use GPU in single rank simulations,
one can use new '-update gpu' command line option.
For use with domain decomposition, please see below.

GPU Direct Communications
"""""""""""""""""""""""""

When running on multiple GPUs with CUDA, communication operations can
now be performed directly between GPU memory spaces (automatically
routed, including via NVLink where available). This behaviour is not
yet enabled by default: the new codepaths have been verified by the
standard |Gromacs| regression tests, but (at the time of release) still
lack substantial "real-world" testing. They can be enabled by setting
the following environment variables to any non-NULL value in your
shell: GMX_GPU_DD_COMMS (for halo exchange communications between PP
tasks); GMX_GPU_PME_PP_COMMS (for communications between PME and PP
tasks); GMX_FORCE_UPDATE_DEFAULT_GPU can also be set in
order to combine with the new GPU update feature (above). The
combination of these will (for many common simulations) keep data
resident on the GPU across most timesteps, avoiding expensive data
transfers. Note that these currently require |Gromacs| to be built
with its internal thread-MPI library rather than any external MPI
library, and are limited to a single compute node. We stress that
users should carefully verify results against the default path, and
any reported issues will be gratefully received to help us mature the
software.


Bonded kernels on GPU have been fused
"""""""""""""""""""""""""""""""""""""

Instead of launching one GPU kernel for each listed interaction type there is now one
GPU kernel that handles all listed interactions. This improves the performance when
running bonded calculations on a GPU.

Delay for ramp-up added to PP-PME tuning
""""""""""""""""""""""""""""""""""""""""

Modern CPUs and GPUs can take a few seconds to ramp up their clock speeds.
Therefore the PP-PME load balancing now starts after 5 seconds instead
of after a few MD steps. This avoids sub-optimal performance settings.
