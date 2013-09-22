
#=============================================================================
# Copyright 2010 Kitware, Inc.
# Copyright 2010 Todd Gamblin <tgamblin@llnl.gov>
# Copyright 2012 Julien Bigot <julien.bigot@cea.fr>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

include(BlueGeneQ-static)
__BlueGeneQ_set_static_flags(XL C)
__BlueGeneQ_set_static_flags(XL CXX)

# This suppression stops the following information message from
# almost every source file at -O3:
#   1500-036: (I) The NOSTRICT option (default at OPT(3)) has the potential to alter the semantics of a program.  Please refer to documentation on the STRICT/NOSTRICT option for more information.
set(COMPILER_SUPPRESSION "-qsuppress=1500-036")

set(CMAKE_SYSTEM_NAME BlueGeneQ-static CACHE STRING "Cross-compiling for BlueGene/Q" FORCE)
# xl.ndebug is appropriate for production calculations. For debugging,
# use xl to add back error checks and assertions
set(CMAKE_C_COMPILER /bgsys/drivers/ppcfloor/comm/xl.ndebug/bin/mpixlc_r)
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG ${COMPILER_SUPPRESSION}" CACHE STRING "Compiler optimization flags")
set(CMAKE_CXX_COMPILER /bgsys/drivers/ppcfloor/comm/xl.ndebug/bin/mpixlcxx_r)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG ${COMPILER_SUPPRESSION}" CACHE STRING "Compiler optimization flags")

mark_as_advanced(CMAKE_XL_CreateExportList) # No idea what spams this
