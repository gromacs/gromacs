***************
Developer guide
***************

This set of pages contains guidelines, instructions, and explanations related
to |Gromacs| development.  The actual code is documented in Doxygen
documentation linked from the front page of the documentation.

The focus is (at least for now) on things that are tightly tied to the code
itself, such as helper scripts that reside in the source repository and
organization of the code itself, and may require the documentation to be
updated in sync.

The guide is currently split into a few main parts:

* Overview of the |Gromacs| codebase.
* Collection of overview pages that describe some important implementation
  aspects.
* Generic guidelines to follow when developing |Gromacs|.
  For some of the guidelines, scripts exist (see below) to automatically
  reformat the code and/or enforce the guidelines for each commit.
* Instructions on what tools are used, and how to use them.

.. toctree::
   :maxdepth: 2

   overview
   build-system
   relocatable-binaries
   style
   tools
