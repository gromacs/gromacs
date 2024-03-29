.. _gmx-reference-manual-rst:

****************
Reference Manual
****************

.. highlight:: bash

.. todo:: this needs to be carefully checked that I didn't mess anything up too bad

.. only:: gmx_image_convert_possible

    This part of the documentation covers implementation details of |Gromacs|.

.. only:: gmx_image_convert_impossible

    The manual could not be properly build because we could not
    convert the images for proper display without Imagemagick.
    The rest of the documentation is still there, but this part here will be missing.

For quick simulation set-up and short explanations,
please refer to the :ref:`User guide <user guide>`.

Help with the installation of |Gromacs| can be found in the
:ref:`Install guide <install guide>`.

If you want to help with developing |Gromacs|, your are most welcome
to read up on the :ref:`Developer Guide <dev guide>` and continue
right away with coding for |Gromacs|.

.. toctree::
    :maxdepth: 2
    :glob:

    preface
    introduction
    definitions
    algorithms*/algorithms
    functions*/functions
    topologies/topologies
    file-formats
    special*/special
    run-parameters
    analysis*/analysis
    details
    averages
    references
