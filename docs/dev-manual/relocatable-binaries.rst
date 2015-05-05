Relocatable binaries
====================

|Gromacs| (mostly) implements the concept of relocatable binaries, i.e., that
after initial installation to ``CMAKE_INSTALL_PREFIX`` (or binary packaging
with CPack), the whole installation tree can be moved to a different folder and
|Gromacs| continues to work without further changes to the installation tree.
This page explains how this is implemented, and the known limitations in the
implementation.  This information is mainly of interest to developers who need
to understand this or change the code, but it can also be useful for people
installing or packaging |Gromacs|.

A related feature that needs to be considered in all the code related to this
is that the executables should work directly when executed from the build tree,
before installation.  In such a case, the data files should also be looked up
from the source tree to make development easy.

Finding shared libraries
------------------------

If |Gromacs| is built with dynamic linking, the first part of making the
binaries relocatable is to make it possible for the executable to find
``libgromacs``, no matter how it is executed.  On platforms that support
a relative RPATH, this is used to make the |Gromacs| executables find the
``libgromacs`` from the same installation prefix.  This makes the executables
fully relocatable when it comes to linking, as long as the relative folder
structure between the executables and the library is kept the same.

If the RPATH mechanism does not work, ``GMXRC`` also adds the absolute path to
the ``libgromacs`` installed with it to ``LD_LIBRARY_PATH``.  On platforms that
support this, this makes the linker search for the library here, but it is less
robust, e.g., when mixing calls to different versions of |Gromacs|.  Note that
``GMXRC`` is currently not relocatable, but hardcodes the absolute path.

On native Windows, DLLs are not fully supported; it is currently only possible
to compile a DLL with MinGW, not with Visual Studio or with Intel compilers.
In this case, the DLLs are placed in the ``bin/`` directory instead of
``lib/`` (automatically by CMake, based on the generic binary type assignment
in ``CMakeLists.txt``).  Windows automatically searches DLLs from the
executable directory, so the correct DLL should always be found.

For external libraries, standard CMake linking mechanisms are used and RPATH
for the external dependencies is included in the executable; on Windows,
dynamic linking may require extra effort to make the loader locate the correct
external libraries.

To support executing the built binaries from the build tree without
installation (critical for executing tests during development), standard CMake
mechanism is used: when the binaries are built, the RPATH is set to the build
tree, and during installation, the RPATH in the binaries is rewritten by CMake
to the final (relative) value.  As an extra optimization, if the installation
tree has the same relative folder structure as the build tree, the final
relative RPATH is used already during the initial build.

The RPATH settings are in the root ``CMakeLists.txt``.  It is possible to
disable the use of RPATH during installation with standard CMake variables,
such as setting ``CMAKE_SKIP_INSTALL_RPATH=ON``.

Finding data files
------------------

The other, |Gromacs|-specific part, of making the binaries relocatable is
to make them able to find data files from the installation tree.  Such data
files are used for multiple purposes, including showing the quotes at the end
of program execution.  If the quote database is not found, the quotes are
simply not printed, but other files (mostly used by system preparation tools
like :ref:`gmx pdb2gmx` and :ref:`gmx grompp`, and by various analysis tools
for static data) will cause fatal errors if not found.

There are several considerations here:

* For relocation to work, finding the data files cannot rely on any hard-coded
  absolute path, but it must find out the location of the executing code by
  inspecting the system.  As a fallback, environment variables or such set by
  ``GMXRC`` or similar could be used (but currently are not).
* When running executables from the build tree, it is desirable that they will
  automatically use the data files from the matching source tree to facilitate
  easy testing.  The data files are not copied into the build tree, and the
  user is free to choose any relative locations for the source and build trees.
  Also, the data files are not in the same relative path in the source tree and
  in the installation tree (the source tree has ``share/top/``, the
  installation tree ``share/gromacs/top/``; the latter is customizable during
  CMake configuration).
* In addition to |Gromacs| executables, programs that link against
  ``libgromacs`` need to be able to find the data files if they call certain
  functions in the library.  In this case, the executable may not be in the
  same directory where |Gromacs| is.  In case of static linking, no part of the
  code is actually loaded from the |Gromacs| installation prefix, which makes
  it impossible to find the data files without external information.
* The user can always use the ``GMXLIB`` environment variable to provide
  alternative locations for the data files, but ideally this should never be
  necessary for using the data files from the installation.

Not all the above considerations are fully addressed by the current
implementation, which works like this:

1. It finds the path to the current executable based on ``argv[0]``.  If the
   value contains a directory, this is interpreted as absolute or as relative
   to the current working directory.  If there is no directory, then a file by
   that name is searched from the directories listed in ``PATH``.  On Windows,
   the current directory is also searched before ``PATH``.  If a file with a
   matching name is found, this is used without further checking.
2. If the executable is found and is a symbolic link, the symbolic links are
   traversed until a real file is found.  Note that links in the directory name
   are not resolved, and if some of the links contain relative paths, the end
   result may contain ``..`` components and such.
3. If an absolute path to the executable was found, the code checks whether the
   executable is located in the build output directory (using ``stat()`` or
   similar to account for possible symbolic links in the directory components).
   If it is, then the hard-coded source tree location is returned.
4. If an absolute path to the executable was found and it was not in the build
   tree, then all parent directories are checked.  If a parent directory
   contains :file:`share/{gromacs}/top/gurgle.dat`, this directory is returned
   as the installation prefix.  The file name ``gurgle.dat`` and the location
   are considered unique enough to ensure that the correct directory has been
   found.  The installation directory for the data files can be customized
   during CMake configuration by setting ``GMX_DATA_INSTALL_DIR``, which
   affects the *gromacs* part in the above path (both in the installation
   structure and in this search logic).

   Note that this search does not resolve symbolic links or normalize the input
   path beforehand: if there are ``..`` components *and* symbolic links in the
   path, the search may proceed to unexpected directories, but this should not
   be an issue as the correct installation prefix should be found before
   encountering such symbolic links (as long as the ``bin/`` directory is not a
   symbolic link).
5. If the data files have not been found yet, try a few hard-coded guesses
   (like the original installation ``CMAKE_INSTALL_PREFIX`` and
   ``/usr/local/``).  The first guess that contains suitable files
   (``gurgle.dat``) is returned.
6. If still nothing is found, return ``CMAKE_INSTALL_PREFIX`` and let the
   subsequent data file opening fail.

The above logic to find the installation prefix is in
``src/gromacs/commandline/cmdlineprogramcontext.cpp``.  Note that code that
links to ``libgromacs`` can provide an alternative implementation for
``gmx::ProgramContextInterface`` for locating the data files, and is then fully
responsible of the above considerations.

Information about the used data directories is printed into the console output
(unless run with ``-quiet``), as well as to (some) error messages when locating
data files, to help diagnosing issues.

There is no mechanism to disable this probing search or affect the process
during compilation time, except for the CMake variables mentioned above.

Known issues
------------

* ``GMXRC`` is not relocatable: it hardcodes the absolute installation path in
  one assignment within the script, which no longer works after relocation.
  Contributions to get rid of this on all the shells the ``GMXRC`` currently
  supports are welcome.
* There is no version checking in the search for the data files; in case of
  issues with the search, it may happen that the installation prefix from some
  other installation of |Gromacs| is returned instead, and only cryptic errors
  about missing or invalid files may reveal this.
* If the searching for the installation prefix is not successful, hard-coded
  absolute guesses are used, and one of those returned.  These guesses include
  the absolute path in ``CMAKE_INSTALL_PREFIX`` used during compilation of
  ``libgromacs``, which will be incorrect after relocation.
* The search for the installation prefix is based on the locating the
  executable.  This does not work for programs that link against
  ``libgromacs``, but are not installed in the same prefix.  For such cases,
  the hard-coded guesses will be used, so the search will not find the correct
  data files after relocation.  The calling code can, however, programmatically
  provide the |Gromacs| installation prefix, but ideally this would work
  without offloading work to the calling code.
* One option to (partially) solve the two above issues would be to use the
  ``GMXDATA`` environment variable set by ``GMXRC`` as the fallback (currently
  this environment variable is set, but very rarely used).
* Installed ``pkg-config`` files are not relocatable: they hardcode the
  absolute installation path.
