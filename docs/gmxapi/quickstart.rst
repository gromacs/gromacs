===========
Quick start
===========

This document provides simple recipes for some specific installation scenarios.
Refer to :doc:`userguide/install` for complete documentation about prerequisites and
building from source code.

.. contents:: Quick start recipes
    :local:
    :depth: 2

Please report bugs or updates at https://redmine.gromacs.org

The following instructions assume the ``bash`` shell and

::

    SOURCE=$PWD/gromacs-gmxapi

Get GROMACS
^^^^^^^^^^^

::

    git clone https://gerrit.gromacs.org/gromacs.git $SOURCE

Build and install GROMACS
^^^^^^^^^^^^^^^^^^^^^^^^^

Refer to standard GROMACS instructions, but be sure to
enable the gmxapi library with the ``-DGMXAPI=ON`` argument.

::

    mkdir build
    cd build
    cmake .. -DGMXAPI=ON
    make install
    cd ..

After installing GROMACS, be sure to "source" the GMXRC. E.g. if you used
``-DCMAKE_INSTALL_PREFIX=/usr/local/gromacs`` as a CMake argument to configure
the install location, in a ``bash`` shell run ``source /usr/local/gromacs/bin/GMXRC``
before proceeding.

Build and install the gmxapi Python package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Assuming ``python`` and ``pip`` are from a Python 3 installation (Python 3.5 or higher)::

    cd $SOURCE/python_packaging
    pip install -r src/requirements.txt
    (cd src && pip install .)

For more detailed instructions, refer to the ``README.md`` file in the ``python_packaging``
directory.

..  todo::

    Document sample_restraint package. Reference issue
    `2893 <https://redmine.gromacs.org/issues/2893>`_ and change
    `11483 <https://gerrit.gromacs.org/c/gromacs/+/11483>`_
