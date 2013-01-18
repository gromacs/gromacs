#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.
#

set(GMXPREPROCESS_SOURCES 
    add_par.c       
    calc_verletbuf.c
    compute_io.c    
    convparm.c      
    gen_ad.c        
    gen_vsite.c     
    genhydro.c   
    gpp_atomtype.c  
    gpp_bond_atomtype.c     
    h_db.c          
    hackblock.c             
    hizzie.c        
    nm2type.c
    pdb2top.c       
    pgutil.c        
    readir.c        
    readpull.c      
    readadress.c      
    readrot.c
    resall.c        
    sorting.c       
    specbond.c      
    ter_db.c        
    tomorse.c       
    topdirs.c       
    topexcl.c       
    topio.c         
    toppush.c       
    topshake.c      
    toputil.c       
    tpbcmp.c        
    vsite_parm.c    
    fflibutil.c
    xlate.c)

set(MDRUN_SOURCES 
    gctio.c    ionize.c runner.c
    do_gct.c     repl_ex.c  xutils.c pme_loadbal.c
    md.c         mdrun.c    genalg.c membed.c)

add_library(gmxpreprocess ${GMXPREPROCESS_SOURCES})
target_link_libraries(gmxpreprocess md)
set_target_properties(gmxpreprocess PROPERTIES OUTPUT_NAME "gmxpreprocess${GMX_LIBS_SUFFIX}" SOVERSION ${SOVERSION} INSTALL_NAME_DIR "${LIB_INSTALL_DIR}"
    COMPILE_FLAGS "${OpenMP_C_FLAGS}")


if(GMX_GPU)
    include_directories(${CMAKE_SOURCE_DIR}/src/gmxlib/gpu_utils)
endif()

if(GMX_OPENMM)
    # Even though the OpenMM build has "moved to contrib", many things
    # have be be done from within the scope of the CMakeLists.txt that
    # builds its mdrun, and that is here
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/src/contrib)
    find_package(OpenMM)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})
    include(${CMAKE_SOURCE_DIR}/src/contrib/BuildMdrunOpenMM.cmake)
endif(GMX_OPENMM)

if(GMX_GPU OR GMX_FORCE_CXX)
    set_source_files_properties(main.c PROPERTIES LANGUAGE CXX)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set_source_files_properties(main.c PROPERTIES COMPILE_FLAGS "-x c++")
    endif()
endif()

if(GMX_FAHCORE)
  add_library(fahcore ${MDRUN_SOURCES})
else(GMX_FAHCORE)

set(GMX_KERNEL_PROGRAMS
    grompp tpbconv pdb2gmx g_protonate gmxdump g_x2top gmxcheck)
if (NOT GMX_NO_QUOTES)
  set(GMX_KERNEL_PROGRAMS ${GMX_KERNEL_PROGRAMS} g_luck)
endif (NOT GMX_NO_QUOTES)


foreach(PROGRAM ${GMX_KERNEL_PROGRAMS})
    add_executable(${PROGRAM} ${PROGRAM}.c main.c)
    if (NOT ${PROGRAM} STREQUAL "g_luck")
        gmx_add_man_page(${PROGRAM})
    endif()
    target_link_libraries(${PROGRAM} gmxpreprocess md gmx ${OpenMP_LINKER_FLAGS})
    set_target_properties(${PROGRAM} PROPERTIES OUTPUT_NAME "${PROGRAM}${GMX_BINARY_SUFFIX}")
endforeach()

add_executable(mdrun ${MDRUN_SOURCES} main.c)
gmx_add_man_page(mdrun)
target_link_libraries(mdrun gmxpreprocess md gmx ${OpenMP_LINKER_FLAGS})
set_target_properties(mdrun PROPERTIES OUTPUT_NAME "mdrun${GMX_BINARY_SUFFIX}" COMPILE_FLAGS "${OpenMP_C_FLAGS}")

if(GMX_OPENMM)
    target_link_libraries(mdrun openmm_api_wrapper)
endif()

# Construct component groups for installation; note that a component may
# belong to only one group
foreach(PROGRAM ${GMX_KERNEL_PROGRAMS})
    set(CPACK_COMPONENT_${PROGRAM}_GROUP tools)
endforeach()
set(CPACK_COMPONENT_MDRUN_GROUP mdrun)

foreach(PROGRAM ${GMX_KERNEL_PROGRAMS} mdrun)
    # Manage CPack component dependencies
    set(CPACK_COMPONENT_${PROGRAM}_DEPENDS libraries libraries-gmxpreprocess)

    # Create custom install-xxxx target
    if (BUILD_SHARED_LIBS)
    # If shared libraries are used, we need to install the libraries in
    # addition to the mdrun binary.
       add_custom_target(install-${PROGRAM}
           COMMAND ${CMAKE_COMMAND} -DCOMPONENT=libraries
                   -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
           COMMAND ${CMAKE_COMMAND} -DCOMPONENT=libraries-gmxpreprocess
                   -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
           COMMAND ${CMAKE_COMMAND} -DCOMPONENT=${PROGRAM}
                   -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
           COMMENT "Installing ${PROGRAM}")
    else()
       add_custom_target(install-${PROGRAM}
           COMMAND ${CMAKE_COMMAND} -DCOMPONENT=${PROGRAM}
                   -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
           COMMENT "Installing ${PROGRAM}")
    endif()
    add_dependencies(install-${PROGRAM} ${PROGRAM})

    # Finally, trigger installation
    install(
      TARGETS ${PROGRAM}
      COMPONENT ${PROGRAM}
      DESTINATION ${BIN_INSTALL_DIR}
      )
endforeach()

install(TARGETS gmxpreprocess DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries-gmxpreprocess)

if (INSTALL_CUDART_LIB) #can be set manual by user
    if (GMX_GPU)
        foreach(CUDA_LIB ${CUDA_LIBRARIES})
            string(REGEX MATCH "cudart" IS_CUDART ${CUDA_LIB})
            if(IS_CUDART) #libcuda should not be installed
                #install also name-links (linker uses those)
                file(GLOB CUDA_LIBS ${CUDA_LIB}*)
                install(FILES ${CUDA_LIBS} DESTINATION
                    ${LIB_INSTALL_DIR} COMPONENT libraries)
            endif()
        endforeach()
    else()
        message(WARNING "INSTALL_CUDART_LIB only makes sense with GMX_GPU")
    endif()
endif ()
endif(GMX_FAHCORE)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libgmxpreprocess.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libgmxpreprocess.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libgmxpreprocess.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "libgmxpreprocess${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)
