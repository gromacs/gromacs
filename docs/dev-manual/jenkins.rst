Understanding Jenkins builds
============================

This page documents what different Jenkins builds actually run from the
|Gromacs| source tree.  The purpose is two-fold:

* Provide information on how to interpret Jenkins failures and how to run the
  same tasks locally to diagnose issues.
* Provide information on what changes in the build system (or other parts of
  the repository) need special care to not break Jenkins builds.

.. TODO: Add a link to a wiki page about general Jenkins documentation, once
   there is more of that.

cppcheck build
--------------

This build runs the :command:`cppcheck` static analysis tool.  Any issues found
mark the build unstable, and can be browsed in Jenkins.

It runs :command:`cmake` to generate the build system, and then builds the
``cppcheck`` target.  Nothing is compiled by this target, it only runs
:command:`cppcheck` for the designated source files.  The CMake configuration
options do not affect the set of files checked, but they do affect the checked
code through :file:`config.h` and such.

Documentation build
-------------------

This build builds various types of documentation:

* PDF reference manual using LaTeX
* Doxygen documentation extracted from the source code
* Set of HTML pages containing an installation guide, a user guide, and a
  developer guide, as well as links to the above.  This set of HTML pages can
  be browsed from Jenkins.
* Man pages
* INSTALL text file

The last three require building the :file:`gmx` binary and running it, so
compilation failures will also show in this build.
All log files that contain warnings are archived as artifacts in the build, and
presence of any warnings marks the build unstable.  Brief description of which
part failed is reported back to Gerrit.

Additionally, the build runs some source code checks that rely on the Doxygen
documentation.  See the description of the ``check-source`` target in
:doc:`gmxtree`.

:doc:`doxygen` provides general guidelines for Doxygen usage, which can be
helpful in understanding and solving Doxygen warnings and some of the
``check-source`` issues.
:doc:`includestyle` provides guidelines for #include order and style, which is
another part of ``check-source`` checks.

Technically, the documentation build runs :file:`admin/build-docs.sh`.
See that file for details of what it exactly builds and how.  Most changes in the
documentation build system will require changes in this script, but Jenkins
configuration should be more static.

uncrustify build
----------------

This build checks the source code for formatting such as consistent indentation
and use of braces, as well as for copyright headers.  See :doc:`formatting` for
the guidelines that are enforced.

Technically, the build simply runs :file:`admin/run-uncrustify.sh` to check the
formatting in the changes for the commit.
If the any changes are required, the build is marked unstable.
If the script completely fails (should be rare), the build fails.
A file with issues found by the script is archived as an artifact in the build,
and a summary is reported back to Gerrit (or the actual issues if there are
only a few).
See :doc:`uncrustify` for more details on uncrustify and on scripts to run it.

.. TODO: Provide links to the build system page, once there are on the git
   commit chain...

.. TODO: Document all the rest.
