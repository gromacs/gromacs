.. _user guide:

.. This file is needed because in some causes we might not be able to build the 
   reference manual (e.g. LaTeX or ImageMagick are missing). Then the file formats
   need to be still available in there previous place in the user guide, instead of
   in the reference manual.

**********
User guide
**********

.. highlight:: bash

This guide provides

* material introducing |Gromacs|
* practical advice for making effective use of |Gromacs|.

For getting, building and installing |Gromacs|, see the
:doc:`/install-guide/index`.
For background on algorithms and implementations, see the
:ref:`reference manual part <gmx-reference-manual-rst>` of the documentation.

|GMX_MANUAL_DOI_STRING|

|GMX_SOURCE_DOI_STRING|

.. TODO This is going to require more organization now that
   we are getting more content available.

.. toctree::
   :maxdepth: 2

   getting-started
   system-preparation
   managing-simulations
   faq
   force-fields
   cutoff-schemes
   cmdline
   mdp-options
   mdrun-features
   mdrun-performance
   run-time-errors
   file-formats
   cmdline
   terminology
   environment-variables
   floating-point
   deprecation-policy
