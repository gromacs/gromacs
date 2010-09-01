                    Welcome to GROMACS!

*Note*: Detailed, step-by-step installation instructions
are available on the website
http://www.gromacs.org/Downloads/Installation_Instructions.

*Note*: If you want to use automake for building look at INSTALL.
        However, automake will be deprecated in releases after 4.5



Cmake  (cross-platform make) is a relatively new build system that
is gaining in popularity. One of the biggest selling points is
compatibility with MS Windows. Starting with the 4.5 release,
it is possible to configure and compile the source code with it.

GETTING CMAKE

Packages for various platforms can be found on the project's download page.
Most of the Linux distributions come with packages available through the
corresponding package manage. Make sure the installed version is 2.6 or later.
Using CMake

Please read carefully the documentation on the CMake website. Developers
may look at some of the online tutorials.

CONFIGURING WITH DEFAULTS SETTINGS

It is advisable that the configuration and the build of the binaries are done
outside of the source tree. On Linux/Mac, the following will configure
the build with the default settings:

$ tar xvfz gromacs-4.5.tar.gz
$ ls
gromacs-4.5
$ mkdir build
$ cd build
$ cmake ../gromacs-4.5
$ make

On multi-core CPU systems (very likely nowadays), the parallel make will do the job much faster:

$ make -j 4
$ make install

Substitute 4 with the number of available cores.

CONFIGURING WITH CUSTOM OPTIONS

Custom options can be set in a few different ways.A list of the more commonly
used ones can be found at http://www.gromacs.org/Developer_Zone/Cmake/Custom_options.

 *command line flag

The values of custom options are supplied with the -D flag. Note that the source path should
be the last argument (otherwise remove the space between -D and the option)!

$ cmake -D GMX_DOUBLE=ON ../gromacs-4.5

 *interactive CMake session

$ cmake -i ../gromacs-4.5

 *curses cmake interface (ccmake)

$ ccmake ../gromacs-4.5

 *CMake GUI

$ cmake-gui ../gromacs-4.5

Explanation about the different options will be presented when using any of the
interactive, curses or gui methods.

All configure options are saved in the CMakeCache.txt file in the build directory.
The file can be edited using a text editor, but after that cmake should be run again.

$ vim CMakeCache.txt
$ cmake ../gromacs-4.5
