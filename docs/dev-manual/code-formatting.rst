.. _gmx-codeformatting:

Automatic source code formatting
================================

.. highlight:: bash

The source code can be automatically formatted using clang-format
(GROMACS 2020 and later)
or uncrustify (GROMACS 2019 and earlier).
Both are formatting tools that apply the guidelines in :doc:`formatting`.
Additionally, other Python scripts are used for a few other automatic
formatting/checking tasks.  The overview tools page contains a list of these
tools: :ref:`dev-formatting-tools`.
This page provides more details for clang-format, uncrustify and copyright scripts.

Jenkins uses these same scripts (in particular, ``clang-format.sh``,
``copyright.sh`` and the ``check-source`` target) to enforce that
the code stays invariant under such formatting.

.. _gmx-uncrustify:

Setting up uncrustify
---------------------

A patched version of uncrustify is used for |Gromacs|.  To set this up, you need
to do these (once):

1. Change to a directory under which you want to build uncrustify and run::

     git clone -b gromacs git://github.com/rolandschulz/uncrustify.git
     cd uncrustify
     ./configure
     make

2. Copy the binary ``src/uncrustify`` into a directory of your choice
   (``/path/to/uncrustify`` below).

Alternatively, if you are running Linux, you can try whether the binary from
http://redmine.gromacs.org/issues/845 works for you.

In order to use the binary for ``uncrustify.sh`` and for the pre-commit hook, you
also need to run this in each of your |Gromacs| repositories::

  git config hooks.uncrustifypath /path/to/uncrustify

Alternatively, if you just want to use ``uncrustify.sh``, you can set the
``UNCRUSTIFY`` environment variable to ``/path/to/uncrustify``.

Using the pre-commit hook or git filters needs additional setup; see the
respective sections below.

Note that Jenkins now only allows formatting using ``clang-format``.

.. _gmx-clang-format:

Setting up clang-format
-----------------------

|Gromacs| formatting is enforced with clang-format 7.0.1.
:command:`clang-format` is one of the core *clang* tools.
It may be included in a *clang* or *llvm* package from your favorite packaging
system or you may find a standalone *clang-format* package,
but you should confirm that the provided command is version 7.0.1 or
7.1.0. Example::

    $ clang-format --version
    clang-format version 7.1.0 (tags/RELEASE_710/final)

If you use a different version of clang-format,
you will likely get different formatting results than
the |Gromacs| continuous integration testing system,
and the commits that you push will fail the automated tests.

.. note::

    Refer to `LLVM <http://releases.llvm.org/download.html#7.1.0>`__ for
    source and binary downloads.
    If downloading sources, note that you will need to download both the
    *LLVM source code* and the *Clang source code*.
    As per the clang
    `INSTALL.txt <https://github.com/llvm/llvm-project/blob/release/7.x/clang/INSTALL.txt>`__,
    place the expanded clang source into a :file:`tools/clang` subdirectory within
    the expanded llvm archive, then run CMake against the llvm source directory.

.. todo::

    Consider referencing or providing binary packages and/or checking/managing
    the executable from an :file:`admin/` script.
    Reference: https://github.com/mongodb/mongo/blob/master/buildscripts/clang_format.py

In order to use the installed version of clang-format for ``clang-format.sh``
and for the pre-commit hook, you also need to run this in each of your |Gromacs| repositories::

  git config hooks.clangformatpath /path/to/clang-format

Alternatively, if you just want to use ``clang-format.sh``, you can set the
``CLANG_FORMAT`` environment variable to ``/path/to/clang-format``.

As above, see the sections below for using the pre-commit hook or git filters.

clang-format discovers which formatting rules to apply from the
:file:`.clang-format` configuration file(s) in project directories,
which will be automatically updated (if necessary) when you :command:`git pull`
from the |Gromacs| repository.
For more about the tool and the :file:`.clang-format` configuration file,
visit https://releases.llvm.org/7.0.1/tools/clang/docs/ClangFormat.html

What is automatically formatted?
--------------------------------

To identify which files are subject to automatic formatting, the scripts use
git filters, specified in ``.gitattributes`` files.  Only files that have the
attribute ``filter`` set to one of the below values are processed:

- ``filter=complete_formatting``: Performs all formatting. Uses clang-format for code formatting.
- ``filter=uncrustify``: uncrustify is run. Deprecated and here for historical reasons.
- ``filter=clangformat``: clang-format is run.
- ``filter=includesort``: include order is enforced and copyright headers are checked.
- ``filter=copyright``: only copyright headers are checked.

Other files are ignored by ``uncrustify.sh``, ``clang-format.sh``,
``copyright.sh`` and ``reformat_all.sh`` scripts (see below).


Scripts
-------

``copyright.py``
^^^^^^^^^^^^^^^^

This script provides low-level functionality to check and update copyright
headers in C/C++ source files, as well as in several other types of files like
CMake and Python scripts.

This file is also used as a loadable Python module for kernel generators, and
provides the functionality to generate conformant copyright headers for such
scripts.

You should rarely need to run this
directly, but instead the bash scripts below use it internally.  You can run
the script with ``--help`` option if you want to see what all options it provides
if you need to do some maintenance on the copyright headers themselves.

``uncrustify.sh``
^^^^^^^^^^^^^^^^^

The information for ``uncrustify`` is mainly provided for historical reasons,
as the actual code formatting is now done using ``clang-format``.

This script runs ``uncrustify`` on modified files and reports/applies the results.
By default, the current HEAD commit is compared to the work tree,
and files that

1. are different between these two trees and
2. change under uncrustify

are reported.  This behavior can be changed by

1. Specifying an ``--rev=REV`` argument, which uses ``REV`` instead of HEAD as
   the base of the comparison.  A typical use case is to specify ``--rev=HEAD^``
   to check the HEAD commit.
2. Specifying an action:

   - ``check-*``:   reports the files that uncrustify changes
   - ``diff-*``:    prints the actual diff of what would change
   - ``update-*``:  applies the changes to the repository
   - ``*-workdir``: operates on the working directory (files on disk)
   - ``*-index``:   operates on the index of the repository

   For convenience, if you omit the workdir/index suffix, workdir is assumed
   (i.e., ``diff`` equals ``diff-workdir``).
3. Specifying ``--uncrustify=off``, which does not run uncrustify.

By default, ``update-*`` refuses to update dirty files (i.e., that differ
between the disk and the index) to make it easy to revert the changes.
This can be overridden by adding a ``-f``/``--force`` option.

``copyright.sh``
^^^^^^^^^^^^^^^^

This script runs ``copyright.py`` on modified files and reports/applies the results.
By default, the current HEAD commit is compared to the work tree,
and files that

1. are different between these two trees and
2. change under have outdated copyright header

are reported.  This behavior can be changed by

1. Specifying an ``--rev=REV`` argument, which uses ``REV`` instead of HEAD as
   the base of the comparison.  A typical use case is to specify ``--rev=HEAD^``
   to check the HEAD commit.
2. Specifying ``--copyright=<mode>``, which alters the level of copyright
   checking is done:

   ``off``
     does not check copyright headers at all
   ``year``
     only update copyright year in new-format copyright headers
   ``add``
     in addition to ``year``, add copyright headers to files that do not
     have any
   ``update``
     in addition to ``year`` and ``add``, also update new-format copyright
     headers if they are broken or outdated
   ``replace``
     replace any copyright header with a new-format copyright header
   ``full``
     do all of the above

By default, ``update-*`` refuses to update dirty files (i.e., that differ
between the disk and the index) to make it easy to revert the changes.
This can be overridden by adding a ``-f``/``--force`` option.

``clang-format.sh``
^^^^^^^^^^^^^^^^^^^

This script runs ``clang-format`` on modified files and reports/applies the results.
By default, the current HEAD commit is compared to the work tree,
and files that

1. are different between these two trees and
2. change under clang-format

are reported.  This behavior can be changed by

1. Specifying an ``--rev=REV`` argument, which uses ``REV`` instead of HEAD as
   the base of the comparison.  A typical use case is to specify ``--rev=HEAD^``
   to check the HEAD commit.
2. Specifying an action:

   - ``check-*``:   reports the files that clang-format changes
   - ``diff-*``:    prints the actual diff of what would change
   - ``update-*``:  applies the changes to the repository
   - ``*-workdir``: operates on the working directory (files on disk)
   - ``*-index``:   operates on the index of the repository

   For convenience, if you omit the workdir/index suffix, workdir is assumed
   (i.e., ``diff`` equals ``diff-workdir``).
3. Specifying ``--format=off``, which does not run clang-format.

By default, ``update-*`` refuses to update dirty files (i.e., that differ
between the disk and the index) to make it easy to revert the changes.
This can be overridden by adding a ``-f``/``--force`` option.

git pre-commit hook
^^^^^^^^^^^^^^^^^^^

If you want to run ``uncrustify.sh``, ``copyright.sh`` and/or
``clang-format.sh`` automatically for changes you make, you can
configure a pre-commit hook using ``admin/git-pre-commit``:

1. Copy the ``git-pre-commit`` script to .git/hooks/pre-commit.

2. Specify the paths to ``uncrustify`` and ``clang-format`` for the hook if you have not already done
   so::

     git config hooks.uncrustifypath /path/to/uncrustify
     git config hooks.clangformatpath /path/to/clang-format

3. Set the operation modes for the hook::

     git config hooks.uncrustifymode check
     git config hooks.clangformatmode check
     git config hooks.copyrightmode  update

With this configuration, all source files modified in the commit are run
through the respective code formatting tool and checked for correct copyright headers.
If any file would be changed by ``uncrustify.sh``, ``clang-format.sh`` or ``copyright.sh``,
the names of those files are reported and the commit is prevented.
The issues can be fixed by running the scripts manually.

To disable the hook without removing the ``pre-commit`` file, you can set ::

  git config hooks.uncrustifymode off
  git config hooks.copyrightmode off
  git config hooks.clangformatmode off

To disable it temporarily for a commit, set NO_FORMAT_CHECK environment
variable.  For example, ::

    NO_FORMAT_CHECK=1 git commit -a

You can also run ``git commit --no-verify``, but that also disables other hooks,
such as the Change-Id ``commit-msg`` hook used by Gerrit.

Note that when you run ``git commit --amend``, the hook is only run for the
changes that are getting amended, not for the whole commit.  During a rebase,
the hook is not run.

The actual work is done by the ``admin/uncrustify.sh``, ``admin/clang-format.sh``
and ``admin/copyright.sh`` scripts, which get run with the ``check-index`` action,
and with ``--uncrustify``, ``--copyright`` and ``--format`` getting set according
to the ``git config`` settings.

``reformat_all.sh``
^^^^^^^^^^^^^^^^^^^

This script runs uncrustify, clang-format, ``copyright.py``, or the include sorter for all
applicable files in the source tree.  See ``reformat_all.sh -h`` for the
invocation.

The script can also produce the list of files for which these commands would be
run.  To do this, specify ``list-files`` on the command line and use
``--filter=<type>`` to specify which command to get the file list for.  This can
be used together with, e.g., ``xargs`` to run other scripts on the same set of
files.

For all the operations, it is also possible to apply patters (of the same style
that various git commands accept, i.e., ``src/*.cpp`` matches all ``.cpp`` files
recursively under ``src/``).  The patterns can be specified with
``--pattern=<pattern>``, and multiple ``--pattern`` arguments can be given.

``-f``/``--force`` is necessary if the working tree and
the git index do not match.


Using git filters
-----------------

An alternative to using a pre-commit hook to automatically apply uncrustify or
clang-format on changes is to use a git filter (does not require either of the scripts,
only the ``.gitattributes`` file).  You can run ::

  git config filter.complete_formatting.clean \
      "/path/to/uncrustify -c admin/uncrustify.cfg -q -l cpp"
  git config filter.clangformat.clean \
      "/path/to/clang-format -i"

To configure a filter for all files that specify ``filter=complete_formatting`` attribute
that indicates that all formatting steps should be performed.

The pre-commit hook + manually running the scripts gives better/more
intuitive control (with the filter, it is possible to have a work tree that is
different from HEAD and still have an empty ``git diff``) and provides better
performance for changes that modify many files.  It is the only way that
currently also checks the copyright headers.

The filter allows one to transparently merge branches that have not been run
through the source checkers, and is applied more consistently (the pre-commit hook is
not run for every commit, e.g., during a rebase).
