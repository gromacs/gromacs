***************
Developer guide
***************

This set of pages contains guidelines, instructions, and explanations related
to |Gromacs| development.  The actual code is documented in Doxygen
documentation linked from the front page of the documentation.

The focus is (at least for now) on things that are tightly tied to the code
itself, such as helper scripts that reside in the source repository, and may
require the documentation to be updated in sync.  Wiki pages at
http://www.gromacs.org/Developer_Zone contain additional information (much of
it outdated, though), and can be linked from relevant locations in this guide.

The guide is currently split into two main parts:

 * Generic guidelines to follow when developing |Gromacs|.
   For some of the guidelines, scripts exist (see below) to automatically
   reformat the code and/or enforce the guidelines for each commit.
 * Instructions on what tools are used, and how to use them.

.. toctree::
   :maxdepth: 2

   style
   tools
