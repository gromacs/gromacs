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

