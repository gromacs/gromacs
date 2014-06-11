if(GMX_USE_TNG)
    # This code decides whether the user wanted to use an external TNG
    # library installed on the system, and if not builds and installs
    # an internal version. Upon exit, TNG_IO_LIBRARIES and
    # TNG_IO_INCLUDE_DIR are always set in the cache, and should be
    # used by GROMACS CMake code to access on TNG functionality

    option(GMX_EXTERNAL_TNG "Use external TNG instead of compiling the version shipped with GROMACS."
           OFF)
    if(GMX_EXTERNAL_TNG)
        # Detect TNG if GMX_EXTERNAL_TNG is explicitly ON
        find_package(TNG_IO 1.6.0)
        if(NOT TNG_IO_FOUND)
            message(FATAL_ERROR
                "TNG >= 1.6.0 not found. "
                "You can set GMX_EXTERNAL_TNG=OFF to compile TNG.")
        endif()
    else()
        # Do the internal build of TNG
        include(ExternalProject)
        include(gmxTestZLib)

        # If TNG finds this, then GROMACS will need to link to it also
        # (and it is useful for the containers that use gzip compression)
        find_package(ZLIB QUIET)
        gmx_test_zlib(ZLIB_WORKS)

        # Build the TNG project the way GROMACS wants it, in particular
        # so that it shares build type and shared-library type
        ExternalProject_Add(
            tng_io
            PREFIX gromacs-tng_io
            SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tng_io
            CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
            -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
            -DTNG_BUILD_WITH_ZLIB=${ZLIB_WORKS}
            -DTNG_BUILD_FORTRAN=off
            -DTNG_BUILD_EXAMPLES=off
            -DTNG_BUILD_COMPRESSION_TESTS=off
            -DTNG_BUILD_DOCUMENTATION=off
            -DTNG_BUILD_TEST=off
            )

        ExternalProject_Get_Property(tng_io install_dir)

        if (BUILD_SHARED_LIBS)
            set(type SHARED)
        else()
            set(type STATIC)
        endif()
        # TODO Later, we can work out how to query the soversion so
        # that this does not require a manual bump alongside a TNG
        # repo change that bumps these. We can't let GROMACS CMake
        # "just handle it," because it doesn't formally know about the
        # libraries that are built by the ExternalProject.
        set(tng_io_so_major_version "1")
        set(tng_io_so_full_version "1.6.0")
        set(tng_compress_so_major_version "1")
        set(tng_compress_so_full_version "1.0.0")

        # Create some fake targets so that GROMACS targets
        # have a place to hang their dependencies. This
        # is standard practice with ExternalProject, and
        # necessary unless the whole build is structured around
        # the use of ExternalProject for everything.
        foreach(component tng_io tng_compress)
            add_library(gromacs_${component} ${type} IMPORTED GLOBAL)
            set(library_file ${install_dir}/${CMAKE_INSTALL_LIBDIR}/lib${component}${CMAKE_${type}_LIBRARY_SUFFIX})
            set_property(TARGET gromacs_${component} PROPERTY
                IMPORTED_LOCATION ${library_file})
            add_dependencies(gromacs_${component} tng_io)

            # Make sure the internal TNG build installs the libraries
            # behind the fake targets, plus the shared-library special
            # files (if appropriate).
            #
            # TODO This install will over-write an existing system TNG
            # at that location. That's nasty, but unlikely to cause a
            # problem in practice during the lifetime of 5.0. If there
            # was a decent way to rename the libraries that
            # ExternalProject (or via CMake in TNG) installs into the
            # GROMACS build tree, then we could eventually install
            # libraries that won't clash.
            install(FILES ${library_file} DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries)
            if(${type} STREQUAL "SHARED")
                # Fortunately, shared libraries are not supported by
                # GROMACS on Windows, because this will not work
                # there.
                install(FILES ${library_file}.${${component}_so_major_version} ${library_file}.${${component}_so_full_version} DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries)
            endif()
        endforeach()

        # Set up include directory
        set(TNG_IO_INCLUDE_DIR ${install_dir}/include CACHE INTERNAL "Path to TNG include files")

        # Set up the list of library targets that GROMACS targets
        # (e.g. libgromacs) can use
        set(TNG_IO_LIBRARIES gromacs_tng_io gromacs_tng_compress)
        if (ZLIB_WORKS)
            list(APPEND TNG_IO_LIBRARIES ${ZLIB_LIBRARIES})
        endif()
        set(TNG_IO_LIBRARIES ${TNG_IO_LIBRARIES} CACHE INTERNAL "TNG libraries")

        # TODO Perhaps add logic to run the TNG tests if we're on a system that can do so
    endif()
endif()
