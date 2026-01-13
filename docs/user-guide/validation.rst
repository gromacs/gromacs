Validation
==========

As with any code, both scientific or commercial, results produced by |Gromacs|
can be incorrect. We provide no guarantees on the correctness of any output.
With that said, because of extensive test coverage and its large user base,
the results produced by |Gromacs| are very reliable. However, new code and features
have not had as extensive testing as established parts of code and there is
therefore a higher risk for bugs, which the user needs to be aware of.
We mark such features with two different tags: "experimental" and "validation pending".
:ref:`gmx mdrun` will notify users in the log file and on the standard output
when such features are in use.

Note that the responsibility of the |Gromacs| developers is only the correctness
of unmodified |Gromacs| code distributed as part of an official release, not that of
modified versions, external packages or plug-ins. An example is
QM/MM simulations. If the results are incorrect because of an issue in the QM-package
used, that is not something the |Gromacs| developers are responsible for.
What the |Gromacs| developers are responsible for in this case is the correctness
of the QM/MM interface (indeed there are tests that cover the QM/MM interface alone).

Experimental features
---------------------

Features labeled "experimental" have not received sufficient testing. Such features
are not suitable for production simulations. We try to keep code belonging to
this category to a minimum. Such features are included in the main codebase
primarily because it simplifies the development process.
Experimental features always need to be activated using an environment variable
and/or by activating a ``cmake`` configuration option.
Users can try such features and are encouraged to provide feedback, in particular
in case of issues.

The current features with experimental status are:

* CUDA graph code, activated by the ``GMX_CUDA_GRAPH`` environment variable
* NVSHMEM communication, enabled at build with ``-DGMX_NVSHMEM=ON`` CMake variable
* OneAPI graph support, enabled at build with ``-DGMX_SYCL_ENABLE_GRAPHS=ON`` CMake variable
* Parallel PME over multiple GPUs with the HIP GPU backend targeting AMD GPUs, enabled at build with ``-DGMX_GPU=HIP`` and ``-DGMX_USE_HEFFTE=ON`` CMake variables
* Direct halo communication, activated by the ``GMX_FILLERS_IN_LOCAL_STATE`` environment variable
* NBNxM 1x1 non-bonded kernels, activated by the ``GMX_NBNXN_PLAINC_1X1`` environment variable
* H5MD trajectory output, activated by selecting the ``h5md`` type for trajectory file output in mdrun
* Using wave64 execution on wave32 (RDNA) devices with HIP, activated by passing ``-mwavefrontsize64`` as an additional compilation flag.
* VkFFT support for evaluating 3D Fast Fourier Transforms (on non-AMD and non-Apple platforms), enabled at build time with ``-DGMX_GPU_FFT_LIBRARY=vkFFT`` CMake variable
* oneMath support for evaluating 3D Fast Fourier Transforms (on any platform), enabled at build time with ``-DGMX_GPU_FFT_LIBRARY=oneMath`` CMake variable
* Double-batched FFT library support for evaluating 3D Fast Fourier Transforms, enabled at build time with ``-DGMX_GPU_FFT_LIBRARY=BBFFT`` CMake variable

Features with validation pending
--------------------------------

Ideally, we would have a validation test suite that covers all combinations of features
of |Gromacs| as well as all possible CPU architectures and GPU backends. Currently we
do not have this, but we are working on a collection of validation systems. Still it
will be challenging to cover the combinatorial explosion of features and options, as well as
all supported hardware. When a new feature or acceleration backend has successfully
passed our internal testing and validation, it will in most cases enter in a |Gromacs|
release in the category "validation pending". For the user this means that we expect
the feature to give correct results, but there is a possibility of incorrect results
with certain combinations of features and/or backends. Such features will never be
active by default; the user will have to actively set an :ref:`gmx mdrun` command line
option or an environment variable, and in some cases an ``mdp`` option.
We expect features to be fully validated over the course of one or two years.
Users are encouraged to test features with validation pending to aid the validation process.
Carefully check your results and please report any issues on the |Gromacs|
`user forum <http://forums.gromacs.org/>`_ or by opening an issue on
`GitLab issues <https://gitlab.com/gromacs/gromacs/-/issues>`_.

The current features with status validation pending are:

* The modular simulator with an integrator different from Velocity Verlet
* The Colvars interface, activated by the ``colvars-active`` ``mdp`` option
* The PLUMED interface, activated by the ``-plumed`` option of :ref:`gmx mdrun`
* The neural network potential interface, activated by the ``nnpot-active`` ``mdp`` option and configuring with LibTorch
* The SYCL GPU backend for non-AMD and non-Intel GPU platforms, activated by choosing the ``SYCL`` option for ``GMX_GPU`` in ``cmake``
* HIP GPU backend targeting AMD GPUs, activated by choosing the ``HIP`` option for ``GMX_GPU`` in ``cmake``
* The Fast Multipole Method interface, enabled at build time by the ``-DGMX_USE_EXT_FMM`` CMake variable
* AMBER LEaP-compatible dihedral reordering by grompp, activated by the ``_FF_AMBER_LEAP_ATOM_REORDERING`` preprocessor define
* Non-bonded free-energy calculations on a GPU, activated by the ``-nbfe gpu`` mdrun option

Feature lifecycle stages
------------------------

1. Development/Experimental: A new feature that may be incomplete, under active development, or not sufficiently tested.

* Risk: High. May contain bugs, the API or behavior may change, and the feature could be removed without notice.
* Usage recommendation: Not for production use. Use only for testing and feedback.

2. Validation pending: A feature that is code-complete and has passed initial tests, but has not been validated across the full matrix of scientific cases, hardware, and parallelization options.

* Risk: Moderate. The feature is expected to be correct, but carries a risk of incorrect results exists in specific, untested scenarios.
* Usage recommendation: Use with caution. Users are encouraged to test and aid in validation by carefully checking results and reporting any issues.

3. Stable: The default state for a feature, validated,  and considered reliable for production simulations.

* Risk: Low (within the bounds of the standard GROMACS disclaimer).
* Usage recommendation: Recommended for all users.

4. Legacy: A feature that is functional but has not been extensively used or tested in recent releases, or it has  been replaced by a newer alternative.

* Risk: Low to moderate (for correctness), but it will not receive updates and may be less performant or flexible than alternatives.
* Usage recommendation: Discouraged for new work. Users should migrate to the recommended modern alternative (if applicable).

5. Deprecated: A formal, end-of-life warning that the feature is no longer recommended and is scheduled for removal in a future major release.

* Risk: Low to moderate (for correctness), high (for future-proofing/continuity). It will receive no bug fixes or maintenance.
* Usage recommendation: Do not use.

6. Removed: The feature's code has been permanently deleted from the GROMACS codebase, documentation still may refer to it as a past feature.

* Risk: N/A
* Usage recommendation: If required, use an older release which still contains it.
