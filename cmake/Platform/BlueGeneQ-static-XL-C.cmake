#=============================================================================
# CMake - Cross Platform Makefile Generator
# Copyright 2000-2011 Kitware, Inc., Insight Software Consortium
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#
# * Neither the names of Kitware, Inc., the Insight Software Consortium,
#   nor the names of their contributors may be used to endorse or promote
#   products derived from this software without specific prior written
#   permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#=============================================================================

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
