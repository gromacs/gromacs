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
# Check GCC version and if any of the 4.1.x family compiler suites is found
# quit the build system generating process. 
#
# The GCC 4.1.x compilers contain an optimization related bug which might 
# results in code that exhibits incorrect behaviour and often leads to 
# exploding systems or crashes. 
#
# For further details see e.g. 
# https://bugs.launchpad.net/ubuntu/+source/gcc-4.1/+bug/158799
#
# Szilard Pall (pszilard@cbr.su.se)
#

if(NOT GMX_DISABLE_GCC41_CHECK)

if(CMAKE_COMPILER_IS_GNUCC)
    # if we have -dumpversion flag use that, otherwise try the --version
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
        RESULT_VARIABLE _gcc_dumpversion_res
        OUTPUT_VARIABLE _gcc_dumpversion_out
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    # if gcc returned with error the -dumpversion is not available 
    if(${_gcc_dumpversion_res} EQUAL 0)
        if(${_gcc_dumpversion_out} MATCHES ".*4\\.1\\.[0-9]+.*")
            message(FATAL_ERROR " The GCC compiler in use seems to belong to the 4.1.x 
                family (detected version: ${_gcc_dumpversion_out}). These compilers 
                contain an optimization related bug which might results in code that 
                exhibits incorrect behaviour and often leads to exploding systems or 
                crashes. To disable this check set GMX_DISABLE_GCC41_CHECK=YES.")
        endif()
    else()    
        message(WARNING " The GCC compiler in use does not support the -dumpversion flag. 
            Will attempt parsing the version from the \"gcc --version\" output.")        
        execute_process(COMMAND ${CMAKE_C_COMPILER} --version
            OUTPUT_VARIABLE _gcc_version_out
            OUTPUT_STRIP_TRAILING_WHITESPACE)            
        if("${_gcc_version_out}" MATCHES ".*4\\.1\\.[0-9]+.*")
            message(FATAL_ERROR " The GCC compiler in use seems to belong to the 4.1.x 
                family. These compiler  compilers contain an optimization related bug 
                which might results in code that exhibits incorrect behaviour and 
                often leads to exploding systems or crashes. To disable this check set 
                GMX_DISABLE_GCC41_CHECK=YES.")
        endif()
    endif()
endif()

endif()
