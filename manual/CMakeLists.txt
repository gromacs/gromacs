INCLUDE(UseLATEX.cmake)
if("${PDFLATEX_COMPILER}" MATCHES "NOTFOUND" OR
   "${IMAGEMAGICK_CONVERT}" MATCHES "NOTFOUND")

  # No pdflatex was found, so don't build the manual.
  set(BUILD_MANUAL OFF PARENT_SCOPE)

  # TODO Identify other dependencies like bibtex, make_index, some
  # graphics conversion program, etc. Perhaps patch UseLATEX.cmake and
  # contribute upstream.
else()
  set(BUILD_MANUAL ON PARENT_SCOPE)

  LATEX_GET_OUTPUT_PATH(output_dir)

  # This helps make the index look decent.
  set(MAKEINDEX_COMPILER_FLAGS -s ${CMAKE_CURRENT_SOURCE_DIR}/hfill.ist)
  # UseLATEX.cmake is hacked to call Anton's subindex script to improve
  # the index generation. Don't replace UseLATEX.cmake blindly!

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
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/mkmdp ARGS ${CMAKE_SOURCE_DIR}/share/html
    DEPENDS mkmdp ${CMAKE_SOURCE_DIR}/share/html/online/mdp_opt.html
    )

  ADD_CUSTOM_COMMAND(OUTPUT ${output_dir}/proglist.tex
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/mkproglist ARGS ${CMAKE_SOURCE_DIR}/admin/programs.txt
    DEPENDS mkproglist ${CMAKE_SOURCE_DIR}/admin/programs.txt
    )

  # Helper macro for making a unique temporary filename
  macro(temp_name fname)
    if(${ARGC} GREATER 1) # Have to escape ARGC to correctly compare
      set(_base ${ARGV1})
    else()
      set(_base ".cmake-tmp")
    endif()
    set(_counter 0)
    while(EXISTS "${_base}${_counter}")
      math(EXPR _counter "${_counter} + 1")
    endwhile()
    set(${fname} "${_base}${_counter}")
  endmacro(temp_name)

  # Now get all the binaries that have been built to produce their .tex
  # output and orchestrate that into the appendix.
  #
  # There's no way for CMake to produce a list of targets from which
  # we'd want to build the manual, even if we set a custom
  # label/property on the right ones. So as we set up executable targets
  # in the rest of the repo, we need to add the names of the targets to
  # a variable that is accessible here. The only ways to get such a
  # global variable is to keep pushing it to parent scope or have ugly
  # hacks with set(... CACHE INTERNAL "" FORCE). We are doing the former
  # with the (list) variable TARGETS_NEEDED_BY_MANUAL_BUILD.

  # do alphabetical sort of targets
  list(SORT TARGETS_NEEDED_BY_MANUAL_BUILD)

  file(MAKE_DIRECTORY ${output_dir}/progman)

  # Later, we append commands that append actual content
  add_custom_command(OUTPUT ${output_dir}/progman/all.tex
    COMMAND rm ARGS -f progman/all.tex
    )

  foreach(program ${TARGETS_NEEDED_BY_MANUAL_BUILD})
    # Build useful .tex file for this program
    temp_name(tempfn ${program}.tex)
    set(program_tex_filename ${program}.tex)
    # This is not portable to a system that lacks grep (e.g. Windows)
    add_custom_command(OUTPUT ${output_dir}/progman/${program_tex_filename}
      COMMAND ${program} ARGS -quiet -man tex
      # Get rid of junk not needed in the manual
      COMMAND grep ARGS -v "Use dialog box GUI to edit command line options" ${program_tex_filename} > ${tempfn}
      COMMAND grep ARGS -v "non blocking priority" ${tempfn} > ${program_tex_filename}
      COMMAND rm ARGS -f ${tempfn}
      DEPENDS ${program}
      WORKING_DIRECTORY ${output_dir}/progman
      )

    # Append the file we just created to the file that will control
    # \inputting it into the LaTeX build. Can't use WORKING_DIRECTORY
    # with APPEND as of CMake 2.8.9.
    add_custom_command(OUTPUT ${output_dir}/progman/all.tex
      COMMAND echo ARGS "\"\\\\input{progman/${program}}\"" >> ${output_dir}/progman/all.tex
      DEPENDS ${output_dir}/progman/${program_tex_filename}
      APPEND
      )
  endforeach()

  temp_name(tempfn g_options.tex)
  add_custom_command(OUTPUT g_options.tex
    COMMAND g_options ARGS -hidden -man tex &> /dev/null
    # Remove header line not suitable for appearing in the manual.
    # This is not portable to a system that lacks grep (e.g. Windows)
    COMMAND grep ARGS -v "normindex" g_options.tex > ${tempfn}
    COMMAND mv ARGS ${tempfn} g_options.tex
    DEPENDS g_options
    WORKING_DIRECTORY ${output_dir}/progman
    )

  # Finally, the command to build the manual.

  ADD_LATEX_DOCUMENT(gromacs.tex
    # Normal LaTeX \included files
    INPUTS algorithms.tex  defunits.tex    implement.tex  macros.tex     special.tex
    analyse.tex     files.tex       gmxpar.tex     install.tex    topology.tex
    averages.tex    forcefield.tex  gromacs.tex    intro.tex      programs.tex
    progman.tex

    # CMake does variable replacement in these files
    CONFIGURE macros.tex

    # These files we're responsible for creating in the
    # add_custom_targets() above. They should not be in the git
    # repository, or its directory, or the dependencies will not work
    # properly.
    DEPENDS proglist.tex mdp_opt.tex ${output_dir}/progman/all.tex g_options.tex

    BIBFILES monster.bib unpubl.bib
    IMAGE_DIRS plots
    DEFAULT_PDF
    USE_INDEX
    )

endif()
