#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
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

# If the user did not set GMX_GPU we'll consider this option to be
# in "auto" mode meaning that we will:
# - search for CUDA and set GMX_GPU=ON we find it
# - check whether GPUs are present
# - if CUDA is not found but GPUs were detected issue a warning
if (NOT DEFINED GMX_GPU)
    set(GMX_GPU_AUTO TRUE CACHE INTERNAL "GPU acceleration will be selected automatically")
else()
    set(GMX_GPU_AUTO FALSE CACHE INTERNAL "GPU acceleration will be selected automatically")
endif()
option(GMX_GPU "Enable GPU acceleration" OFF)

option(GMX_GPU_USE_AMD "GPU from AMD Radeon" OFF) 
option(GMX_GPU_USE_NVIDIA "GPU from NVIDIA" OFF) 

option(GMX_CLANG_CUDA "Use clang for CUDA" OFF)

if(GMX_GPU AND GMX_DOUBLE)
    message(FATAL_ERROR "GPU acceleration is not available in double precision!")
endif()
if(GMX_GPU_AUTO AND GMX_DOUBLE)
    message(WARNING "GPU acceleration is not available in double precision, disabled!")
    set_property(CACHE GMX_GPU PROPERTY VALUE OFF)
    set_property(CACHE GMX_GPU_AUTO PROPERTY VALUE OFF)
endif()

if(GMX_GPU OR GMX_GPU_AUTO )
    if(NOT GMX_GPU_USE_AMD AND NOT GMX_GPU_USE_NVIDIA)
       message(FATAL_ERROR "Either AMD Radeon or Nvidia GPU should be used") 
    endif() 
    if(GMX_GPU_USE_AMD AND GMX_GPU_USE_NVIDIA) 
       message(FATAL_ERROR "AMD Radeon or Nvidia GPU, only one can be choosen, but not both") 
    endif()

    if(GMX_GPU_USE_AMD) 
        include(gmxManageGPUAmd) 
    else()
        include(gmxManageGPUNvidia) 
    endif() 
endif() 

