#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2014, by the GROMACS development team, led by
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

# - Define macro to check if all of the following work:
# sched_getaffinity()
# sched_setaffinity()
# CPU_ZERO()
# CPU_SET()
# CPU_ISSET()
# CPU_CLR()
# CPU_COUNT()

#  test_sched_affinity(VARIABLE)
#
#  VARIABLE will be set to true if all of the functions link fine.

MACRO(test_sched_affinity VARIABLE)

  if(NOT DEFINED sched_affinity_compile)
    MESSAGE(STATUS "Checking for sched.h GNU affinity API")

    check_c_source_compiles(
      "#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <sched.h>
int main(void) {
  int i;
  cpu_set_t mask;
  CPU_ZERO(&mask);
  sched_getaffinity(0, sizeof(cpu_set_t), &mask);
  if(CPU_ISSET(0,&mask))
  {
    CPU_CLR(0,&mask);
    CPU_SET(0,&mask);
  }
  sched_setaffinity(0, sizeof(cpu_set_t), &mask);
  return CPU_COUNT(&mask);
}" sched_affinity_compile)
  endif(NOT DEFINED sched_affinity_compile)

  if(sched_affinity_compile)
    set(${VARIABLE} 1 CACHE INTERNAL "Result of test for sched.h GNU affinity API" FORCE)
  else()
    set(${VARIABLE} 0 CACHE INTERNAL "Result of test for sched.h GNU affinity API" FORCE)
  endif()
ENDMACRO(test_sched_affinity VARIABLE)
