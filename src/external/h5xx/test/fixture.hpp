/*
 * Copyright © 2013-2014 Felix Höfling
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#ifndef TEST_FIXTURE_HPP
#define TEST_FIXTURE_HPP

#include <boost/test/unit_test_log.hpp>

#include <h5xx/file.hpp>

/**
 * fixture that provides and cleans up an HDF5 file
 */
template <char const* filename>
struct h5file
{
    h5file()
      : file(filename, h5xx::file::trunc)
    {
//      The following syntax is preferred, but seems to be broken with Boost 1.65.0
//      BOOST_TEST_MESSAGE("HDF5 file created: " << filename);
        BOOST_TEST_MESSAGE(std::string("HDF5 file created: ") + filename);
    }

    ~h5file()
    {
        file.close(true);
#ifdef NDEBUG
        remove(filename);
        BOOST_TEST_MESSAGE(std::string("HDF5 file removed: ") + filename);
#endif
    }

    h5xx::file file;
};

#endif // ! TEST_FIXTURE_HPP
