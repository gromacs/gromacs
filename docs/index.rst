Welcome to the |Gromacs| documentation!
=======================================

..  TODO : consolidate at least some of the material in the
    Documentation links below into the new user guide, along with all
    of http://www.gromacs.org/Documentation/Cut-off_schemes,
    http://www.gromacs.org/Documentation/Acceleration_and_parallelization
    and http://www.gromacs.org/Documentation/Performance_checklist)

Contents:

.. toctree::
   :maxdepth: 2

   download
   install-guide/index

   user-guide/index

* `Reference Manual`_ (PDF format)

* Answers to `Frequently Asked Questions <http://www.gromacs.org/Documentation/FAQs>`_

* Coping with `errors while using GROMACS <http://www.gromacs.org/Documentation/Errors>`_

* Links to `tutorial material <http://www.gromacs.org/Documentation/Tutorials>`_

Documentation for developers
----------------------------

The developer documentation currently consists of two parts:

* A developer guide that provides an overview of the |Gromacs| codebase, and
  includes more detailed resouces such as guidelines and information on tools
  used during development.
* Doxygen documentation extracted from comments in C/C++ code, documenting the
  actual C/C++ code.

Some overview documentation that is closely related to the actual C/C++ code
appears in the Doxygen documentation, while some other overview content is in
the developer guide.  The reasons are partially technical, but crosslinks
between the developer guide and the Doxygen documentation are provided whenever
related content appears split between the two sources.

The documentation does not yet cover all areas, but more content is being
(slowly) added.
Wiki pages at <http://www.gromacs.org/Developer_Zone> may contain additional
information (much of it outdated, though), and can be linked from relevant
locations in the developer guide.

Contents:

.. toctree::
   :maxdepth: 3

   dev-manual/index

* Doxygen documentation

  * `Public API documentation <doxygen/html-user/index.xhtml>`_
       This contains documentation for code that appears in installed headers,
       as well as some overview documentation to understand those parts of the
       code.
       Please note that the definition of the public API is very preliminary
       and subject to change, in particular for parts that have not been
       documented.
  * `Code documentation <doxygen/html-lib/index.xhtml>`_
       This contains the public API documentation as a subset, but also has more
       details on the internal implementation of |Gromacs|.  This is a good
       place to start to understand some specific area of the code in general
       to, e.g., contribute.
  * `Full documentation <doxygen/html-full/index.xhtml>`_
       This contains every single documented function from the codebase,
       including internal  There can be an overwhelming amount of detail, but
       this can be useful if trying to find some specific function within the
       codebase.


Indices and tables
==================

* :ref:`genindex`

.. _Reference Manual: `gmx-manual`_
