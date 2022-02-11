#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2012- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

if(GMX_GPU_CUDA)
    if(NOT GMX_CUDA_NB_SINGLE_COMPILATION_UNIT)
        set(NBNXM_CUDA_KERNEL_SOURCES
                nbnxm_cuda_kernel_F_noprune.cu
                nbnxm_cuda_kernel_F_prune.cu
                nbnxm_cuda_kernel_VF_noprune.cu
                nbnxm_cuda_kernel_VF_prune.cu
                nbnxm_cuda_kernel_pruneonly.cu)
    endif()

    file(GLOB NBNXM_CUDA_SOURCES
         nbnxm_cuda.cu
         nbnxm_cuda_data_mgmt.cu
         nbnxm_cuda_jit_support.cu
         nbnxm_gpu_buffer_ops_internal.cu
         ${NBNXM_CUDA_KERNEL_SOURCES})
    set(NBNXM_SOURCES ${NBNXM_SOURCES} ${NBNXM_CUDA_SOURCES} PARENT_SCOPE)
endif()
