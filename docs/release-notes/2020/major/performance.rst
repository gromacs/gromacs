Performance improvements
^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on redmine, without the
   a space between the colon and number!

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

Bonded kernels on GPU have been fused
"""""""""""""""""""""""""""""""""""""

Instead of launching one GPU kernel for each listed interaction type there is now one
GPU kernel that handles all listed interactions. This improves the performance when
running bonded calculations on a GPU.
