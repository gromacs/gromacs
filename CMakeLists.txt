cmake_minimum_required(VERSION 2.6)
# Keep CMake suitably quiet on Cygwin
set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required

PROJECT(gromacs-manual NONE)
INCLUDE(UseLATEX.cmake)
LATEX_GET_OUTPUT_PATH(output_dir)

# Some variables for manual version numbers. These get transferred
# into the LaTeX.
set(MANUAL_MAJOR_VERSION "4")
set(MANUAL_VERSION "4.5.4")

# Some cache variables that the user needs to set
set(
  GMXSRC ""
  CACHE PATH "Path to GROMACS source distribution"
  )
set(
  GMXBIN "/usr/local/gromacs/bin"
  CACHE PATH "Path to GROMACS installed binaries (i.e. including the final '/bin')"
  )

# This helps make the index look decent.
set(MAKEINDEX_COMPILER_FLAGS -s ${CMAKE_CURRENT_SOURCE_DIR}/hfill.ist)
# UseLATEX.cmake is hacked to call Anton's subindex script to improve
# the index generation. Don't replace UseLATEX.cmake blindly!

# We need the name of the program that will print the GROMACS standard
# and hidden options.
set(OPTIONS_PROGRAM_NAME g_options)

if(NOT EXISTS ${GMXSRC})
  MESSAGE(FATAL_ERROR "Please point the GMXSRC CMake variable to the path of the GROMACS source directory.")
endif(NOT EXISTS ${GMXSRC})

if(NOT EXISTS ${GMXBIN} OR NOT EXISTS "${GMXBIN}/GMXRC")
  MESSAGE(FATAL_ERROR "Please point the GMXBIN CMake variable at the path where GROMACS binaries are installed. There should not be multiple GROMACS programs installed there (e.g. with _mpi, _d suffixes).")
endif(NOT EXISTS ${GMXBIN} OR NOT EXISTS "${GMXBIN}/GMXRC")

set(REMOVE_OLD_FILES "")
if(EXISTS ${CMAKE_SOURCE_DIR}/g_options.tex
    OR EXISTS ${CMAKE_SOURCE_DIR}/proglist.tex
    OR EXISTS ${CMAKE_SOURCE_DIR}/progman.tex
    OR EXISTS ${CMAKE_SOURCE_DIR}/mdp_opt.tex)
  message(FATAL_ERROR "You have at least one of g_options.tex, proglist.tex, progman.tex or mdp_opt.tex in your manual directory, probably from a pre-CMake build of the manual. These must be removed in order for CMake to build the dependencies correctly. One way to do this is to switch to the manual directory and use \"make clean\".")
endif(EXISTS ${CMAKE_SOURCE_DIR}/g_options.tex
  OR EXISTS ${CMAKE_SOURCE_DIR}/proglist.tex
  OR EXISTS ${CMAKE_SOURCE_DIR}/progman.tex
  OR EXISTS ${CMAKE_SOURCE_DIR}/mdp_opt.tex)

# Here follows various custom CMake "commands" that call our various
# GROMACS helper scripts to make bits of the manual from the installed
# binaries and/or source distribution, as appropriate.
#
# Only files that can be built should be listed in DEPENDS. Makefile
# rules are generated for those files.
#
# The .mdp options section is now intended to be built directly
# from a GROMACS source directory. This facilitates drafting
# updates in that HTML file without copying files by hand.

ADD_CUSTOM_COMMAND(OUTPUT ${output_dir}/mdp_opt.tex
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/mkmdp ARGS ${GMXSRC}/share/html
  DEPENDS mkmdp ${GMXSRC}/share/html/online/mdp_opt.html
  )

ADD_CUSTOM_COMMAND(OUTPUT ${output_dir}/proglist.tex
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/mkproglist ARGS ${GMXSRC}/admin/programs.txt
  DEPENDS mkproglist ${GMXSRC}/admin/programs.txt
  )

ADD_CUSTOM_COMMAND(OUTPUT ${output_dir}/progman.tex
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/mkman ARGS ${GMXBIN}
  DEPENDS mkman
  )
# In theory, each executable should be a dependency for the above command...
# I can't be bothered about it

# Finally, the command to build the manual.

ADD_LATEX_DOCUMENT(gromacs.tex
  # Normal LaTeX \included files
  INPUTS algorithms.tex  defunits.tex    implement.tex  macros.tex     special.tex
         analyse.tex     files.tex       gmxpar.tex     install.tex    topology.tex
         averages.tex    forcefield.tex  gromacs.tex    intro.tex      programs.tex

 # CMake does variable replacement in these files
  CONFIGURE macros.tex

  # These files we're responsible for creating. They should not be in
  # the git repository, or its directory, or the dependencies will not
  # work properly. Using "make clean" from the manual directory will
  # clean these up if you have them lying around from old-style
  # manual builds.
  DEPENDS proglist.tex mdp_opt.tex progman.tex

  BIBFILES monster.bib unpubl.bib
  IMAGE_DIRS plots
  DEFAULT_PDF
  USE_INDEX
  )