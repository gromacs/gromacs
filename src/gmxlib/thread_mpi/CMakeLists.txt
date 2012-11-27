#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012, by the GROMACS development team, led by
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

# Note that not all .c files are compiled directly: some of them 
# are #included (some multiple times) from other source files.
set(THREAD_MPI_LIB_SOURCE 
    alltoall.c      reduce.c        scan.c
    barrier.c       list.c          reduce_fast.c
    bcast.c         lock.c          scatter.c
    collective.c    once.c          tmpi_init.c
    comm.c          p2p_protocol.c  topology.c
    errhandler.c    p2p_send_recv.c type.c
    event.c         p2p_wait.c      tmpi_malloc.c
    gather.c        profile.c
    group.c         numa_malloc.c   )


if (THREAD_PTHREADS)
    set(THREAD_MPI_LIB_SOURCE pthreads.c ${THREAD_MPI_LIB_SOURCE})
else (THREAD_PTHREADS)
    if (THREAD_WINDOWS)
        set(THREAD_MPI_LIB_SOURCE winthreads.c ${THREAD_MPI_LIB_SOURCE})
    endif (THREAD_WINDOWS)
endif(THREAD_PTHREADS)


add_library(thread_mpi STATIC ${THREAD_MPI_LIB_SOURCE})

# make it link to the threads library (-lpthreads, for example)
target_link_libraries(thread_mpi ${THREAD_LIB})

if (TMPI_CXX_LIB)
    set(THREAD_MPI_CXX_LIB_SOURCE 
        system_error.cpp )
    add_library(thread_mpi_cxx STATIC ${THREAD_MPI_CXX_LIB_SOURCE})
    target_link_libraries(thread_mpi_cxx thread_mpi)
endif (TMPI_CXX_LIB)

#configure_file(tmpi_config.h.cmakein  tmpi_config.h)
#add_definitions(-DHAVE_TMPI_CONFIG_H)
#include_directories( . )

# this configure file only works here. If you include the source
# files into your project, insert the contest of the tmpi_config.h.cmakein
# into your config.h.cmakein file
configure_file(tmpi_config.h.cmakein  ${CMAKE_CURRENT_BINARY_DIR}/tmpi_config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_definitions(-DHAVE_TMPI_CONFIG_H)

install(TARGETS thread_mpi DESTINATION lib)
