Performance improvements
^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on redmine, without the
   a space between the colon and number!

Up to a factor 2.5 speed-up of the non-bonded free-energy kernel
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The non-bonded free-energy kernel is a factor 2.5 faster with non-zero A and B
states and a factor 1.5 with one zero state. This especially improves the run
performance when non-perturbed non-bondeds are offloaded to a GPU. In that case
the PME-mesh calculation now always takes the most CPU time.


Proper dihedrals of Fourier type and improper dihedrals of preriodic type are SIMD accelerated
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Avoid configuring the own-FFTW with AVX512 enabled when GROMACS does not use AVX512
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Previously if GROMACS was configured to use any AVX flavor, the internally built FFTW
would be configured to also contain AVX512 kernels. This could cause performance loss
if the (often noisy) FFTW auto-tuner picks an AVX512 kernel in a run that otherwise 
only uses AVX/AVX2 which could run at higher CPU clocks without AVX512 clock speed limitation.
Now AVX512 is only used for the internal FFTW if GROMACS is also configured with
the same SIMD flavor.

Update and constraints can run on a (single) GPU
""""""""""""""""""""""""""""""""""""""""""""""""

For standard constant-NVT runs (see the user guide for more details),
update and constraints can be offloaded to a GPU with CUDA. Thus all compute
intensive parts of a simulation can be offloaded, which provides
better performance when using a fast GPU combined with a slow CPU.
Note that this does not work with domain decomposition yet.

Bonded kernels on GPU have been fused
"""""""""""""""""""""""""""""""""""""

Instead of launching one GPU kernel for each listed interaction type there is now one
GPU kernel that handles all listed interactions. This improves the performance when
running bonded calculations on a GPU.
