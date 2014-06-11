if(GMX_USE_TNG AND NOT GMX_EXTERNAL_TNG)
    include(ExternalProject)

    # If TNG finds this, then GROMACS will need to link to it also
    # (and it is useful for the containers that use gzip compression)
    find_package(ZLIB QUIET)

    ExternalProject_Add(
        tng_io
        PREFIX gromacs-tng_io
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tng_io
        CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DBUILD_SHARED_LIBS=off
            -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
            -DTNG_BUILD_WITH_ZLIB=${ZLIB_FOUND}
            -DTNG_BUILD_FORTRAN=off
            -DTNG_BUILD_EXAMPLES=off
            -DTNG_BUILD_COMPRESSION_TESTS=off
            -DTNG_BUILD_DOCUMENTATION=off
            -DTNG_BUILD_TEST=off
        )

    ExternalProject_Get_Property(tng_io install_dir)

    foreach(component tng_io tng_compress)
        add_library(gromacs_${component} STATIC IMPORTED GLOBAL)
        set_property(TARGET gromacs_${component} PROPERTY 
            IMPORTED_LOCATION ${install_dir}/${CMAKE_INSTALL_LIBDIR}/lib${component}${CMAKE_STATIC_LIBRARY_SUFFIX})
        add_dependencies(gromacs_${component} tng_io)
    endforeach()

    set(TNG_IO_LIBRARIES gromacs_tng_io gromacs_tng_compress ${ZLIB_LIBRARIES} PARENT_SCOPE)
    # TODO later, add logic to run the TNG tests if we're on a system that can do so

endif()