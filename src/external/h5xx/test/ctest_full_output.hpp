/*
 * Copyright © 2013  Felix Höfling
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#ifndef CTEST_FULL_OUTPUT_HPP
#define CTEST_FULL_OUTPUT_HPP

#include <boost/test/unit_test_log.hpp>

/**
 * Print CTEST_FULL_OUTPUT to avoid truncation of output by CTest.
 *
 * Requires --log_level=message or --log_level=test_suite.
 *
 * Adding BOOST_GLOBAL_FIXTURE( ctest_full_output ) to the test source triggers
 * the output before any test run.
 */
struct ctest_full_output
{
    ctest_full_output()
    {
        BOOST_TEST_MESSAGE( "Avoid truncation of output by CTest: CTEST_FULL_OUTPUT" );
    }
};

#endif // ! CTEST_FULL_OUTPUT_HPP
