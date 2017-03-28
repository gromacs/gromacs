#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014,2016, by the GROMACS development team, led by
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

# "tests" target builds all the separate test binaries.
add_custom_target(tests)

# "run-ctest" is an internal target that actually runs the tests.
# This is necessary to be able to add separate targets that execute as part
# of 'make check', but are ensured to be executed after the actual tests.
add_custom_target(run-ctest
                  COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
                  COMMENT "Running all tests"
                  USES_TERMINAL VERBATIM)
add_dependencies(run-ctest tests)
# "check-all" target builds and runs all tests.
add_custom_target(check-all DEPENDS run-ctest)

# "run-ctest-nophys" is an internal target that actually runs the tests.
# It replicates the behavior of "run-ctest" and 'make check'
# before physical validation tests were added.
add_custom_target(run-ctest-nophys
                  COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -E physvaltests
                  COMMENT "Running all tests except physical validation"
                  USES_TERMINAL VERBATIM)
add_dependencies(run-ctest-nophys tests)
# "check" target builds and runs all tests except physical validation.
add_custom_target(check DEPENDS run-ctest-nophys)

# "run-ctest-phys" is an internal target that actually runs the tests.
# It only runs the physical validation tests.
add_custom_target(run-ctest-phys
                  COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -R physvaltests
                  COMMENT "Running physical validation tests"
                  USES_TERMINAL VERBATIM)
# "check-phys" target runs only physical validation tests
add_custom_target(check-phys DEPENDS run-ctest-phys)

# Global property for collecting notices to show at the end of the targets.
# Expanded to avoid to show messages about physical validation when only
# the other tests are ran, and vice-versa.
set_property(GLOBAL PROPERTY GMX_TESTS_NOTICE)
set_property(GLOBAL PROPERTY GMX_TESTS_NOTICE_NOPHYS)
set_property(GLOBAL PROPERTY GMX_TESTS_NOTICE_PHYS)

# Original function separated into two adding to the different properties.
function (gmx_add_missing_tests_notice_nophys TEXT)
    set_property(GLOBAL APPEND PROPERTY GMX_TESTS_NOTICE ${TEXT})
    set_property(GLOBAL APPEND PROPERTY GMX_TESTS_NOTICE_NOPHYS ${TEXT})
endfunction()
function (gmx_add_missing_tests_notice_phys TEXT)
    set_property(GLOBAL APPEND PROPERTY GMX_TESTS_NOTICE ${TEXT})
    set_property(GLOBAL APPEND PROPERTY GMX_TESTS_NOTICE_PHYS ${TEXT})
endfunction()

# Original function expanded to handle different targets / properties,
# and added a check to avoid the target being added if no notices are
# to be displayed.
function (gmx_create_missing_tests_notice_target)
    get_property(_text GLOBAL PROPERTY GMX_TESTS_NOTICE)
    set(_cmds)
    foreach (_line ${_text})
        list(APPEND _cmds COMMAND ${CMAKE_COMMAND} -E echo "NOTE: ${_line}")
    endforeach()
    list(LENGTH _cmds n)
    # checking whether any messages should be displayed
    if(${n})
        add_custom_target(missing-tests-notice
                ${_cmds}
                DEPENDS run-ctest
                COMMENT "Some tests not available" VERBATIM)
        add_dependencies(check-all missing-tests-notice)
    endif()
    get_property(_text GLOBAL PROPERTY GMX_TESTS_NOTICE_NOPHYS)
    set(_cmds)
    foreach (_line ${_text})
        list(APPEND _cmds COMMAND ${CMAKE_COMMAND} -E echo "NOTE: ${_line}")
    endforeach()
    list(LENGTH _cmds n)
    # checking whether any messages should be displayed
    if(${n})
        add_custom_target(missing-tests-notice-nophys
                ${_cmds}
                DEPENDS run-ctest-nophys
                COMMENT "Some tests not available" VERBATIM)
        add_dependencies(check missing-tests-notice-nophys)
    endif()
    get_property(_text GLOBAL PROPERTY GMX_TESTS_NOTICE_PHYS)
    set(_cmds)
    foreach (_line ${_text})
        list(APPEND _cmds COMMAND ${CMAKE_COMMAND} -E echo "NOTE: ${_line}")
    endforeach()
    list(LENGTH _cmds n)
    # checking whether any messages should be displayed
    if(${n})
        add_custom_target(missing-tests-notice-phys
                ${_cmds}
                DEPENDS run-ctest-phys
                COMMENT "Some tests not available" VERBATIM)
        add_dependencies(check-phys missing-tests-notice-phys)
    endif()
endfunction()
