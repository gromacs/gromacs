.. _dev guide:

***************
Developer Guide
***************

.. highlight:: bash

This set of pages contains guidelines, instructions, and explanations related
to |Gromacs| development.  The actual code is documented in Doxygen
documentation linked below.

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

.. only:: html

        In addition to this, Doxygen documentation extracted from the comments
        in the C/C++ code is available to document the actual existing code.

.. only:: latex

        **The full code documentation generated from Doxygen can be found in the online
        documentation. It is not included here in order to save the trees.**

Some overview documentation that is closely related to the actual C/C++ code
appears in the Doxygen documentation, while some other overview content is in
the developer guide.  The reasons are partially technical, but crosslinks
between the developer guide and the Doxygen documentation are provided whenever
related content appears split between the two sources.

The documentation does not yet cover all areas, but more content is being
(slowly) added.

.. toctree::
   :maxdepth: 2

   contribute
   overview
   build-system
   change-management
   relocatable-binaries
   documentation-generation
   style
   tools
   known-issues

*********************
Doxygen documentation
*********************

.. only:: html

  * `Public API documentation <../doxygen/html-user/index.xhtml>`_
       This contains documentation for code that appears in installed headers,
       as well as some overview documentation to understand those parts of the
       code.
       Please note that the definition of the public API is very preliminary
       and subject to change, in particular for parts that have not been
       documented.
  * `Code documentation <../doxygen/html-lib/index.xhtml>`_
       This contains the public API documentation as a subset, but also has more
       details on the internal implementation of |Gromacs|.  This is a good
       place to start to understand some specific area of the code in general
       to, e.g., contribute.
  * `Full documentation <../doxygen/html-full/index.xhtml>`_
       This contains every single documented function from the codebase,
       including internal  There can be an overwhelming amount of detail, but
       this can be useful if trying to find some specific function within the
       codebase.

.. only:: latex

    The doxygen code documentation is available on the |Gromacs| webpage.


