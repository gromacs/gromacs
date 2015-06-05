Understaning Jenkins builds
===========================

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

The cppcheck build builds the ``cppcheck`` target.

Documentation build
-------------------

The documentation build runs :file:`admin/build-docs.sh`.
See that file for details of which targets it builds.  Most changes in the
documentation build system will require changes in this script, but Jenkins
configuration should be more static.
All log files that contain warnings are archived as artifacts in the build.

uncrustify build
----------------

The uncrustify build runs :file:`admin/uncrustify.sh` as ::

  admin/uncrustify.sh check --rev=HEAD^

to check the formatting in the changes for the commit.
If the return value is non-zero, the build is marked unstable.
Currently, only the console log of the build shows what actually was the issue.
See :doc:`uncrustify` for more details on the script and uncrustify.

.. TODO: Provide links to the build system page, once there are on the git
   commit chain...

.. TODO: Document all the rest.
