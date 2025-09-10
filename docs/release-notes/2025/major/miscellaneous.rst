Miscellaneous
^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Internal build of FFTW now uses version 3.3.10 
""""""""""""""""""""""""""""""""""""""""""""""

When the ``-DGMX_BUILD_OWN_FFTW=ON`` option is enabled, we now compile FFTW 3.3.10 instead
of the previous 3.3.8 version. For more details, please refer to the
`FFTW release notes <https://www.fftw.org/release-notes.html>`_.

Increased AWH parameter 'awh-nsamples-update' default value from 10 to 100
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This is to decrease the overhead of updating the bias, in particular with multiple walkers.

Support for continuing expanded ensemble equilibration across simulation parts
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The mdp options `init-lambda-counts` and `init-wl-histogram-counts`
can now initialize the number of counts at each sampled lambda state
and the Wang-Landau histograms used to determine simulation
equilibration. These are most useful when running short simulation
parts, so that the information about how the system is equilibrating
can be carried over between simulations.  Otherwise, chains of short
simulation parts would never converge when using expanded-ensemble methods.
The information is now also carried over through the checkpoint file.

CTest can now run tests in parallel
"""""""""""""""""""""""""""""""""""

Properly express the resource requirements of test binaries and
their dependencies so that ``ctest --parallel`` can be used safely.
