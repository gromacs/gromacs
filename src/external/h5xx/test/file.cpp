/*
 * Copyright © 2013      Manuel Dibak
 * Copyright © 2013-2014 Felix Höfling
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#define BOOST_TEST_MODULE h5xx_file
#include <boost/test/unit_test.hpp>

#include <h5xx/file.hpp>

#include <unistd.h>
#include <test/ctest_full_output.hpp>
#include <test/catch_boost_no_throw.hpp>

char const* name = "test_h5xx_file.h5";

using namespace h5xx;

BOOST_GLOBAL_FIXTURE( ctest_full_output );

// test standard use cases
BOOST_AUTO_TEST_CASE( use_cases )
{
    file h5file;
    BOOST_CHECK_NO_THROW(h5file.open(name));
    BOOST_CHECK_THROW(h5file.open(name), h5xx::error);         // open again
    BOOST_CHECK_NO_THROW(h5file.flush());
    BOOST_CHECK_EQUAL(h5file.name(), std::string(name));
    BOOST_CHECK_NO_THROW(h5file.close());
    BOOST_CHECK_NO_THROW(h5file.close());                     // close is silent for closed files
    BOOST_CHECK_NO_THROW(h5file.flush());                     // flush is silent for closed files
    BOOST_CHECK(h5file.hid() == -1);
    BOOST_CHECK(is_hdf5_file(name) > 0);
    unlink(name);
}

// test opening modes
BOOST_AUTO_TEST_CASE( open )
{
    BOOST_CHECK_NO_THROW(file());                              // default trivial constructor and destructor
//    BOOST_CHECK_NO_THROW(file(name));                        // warning: this means a local declaration of "name"
    BOOST_CHECK_NO_THROW(file f(name));                        // default: in | out, create file
    BOOST_CHECK(is_hdf5_file(name) > 0);
    BOOST_CHECK_NO_THROW(file f(name));                        // re-open existing file
    BOOST_CHECK_NO_THROW(file(name, file::in));                // read-only
    BOOST_CHECK_NO_THROW(file(name, file::in | file::out));    // append
    BOOST_CHECK_NO_THROW(file(name, file::out));               // append
    BOOST_CHECK_NO_THROW(file(name, file::trunc));             // write and truncate
    BOOST_CHECK_NO_THROW(file(name, file::out | file::trunc)); // write and truncate

    // remove and recreate file
    unlink(name);
    BOOST_CHECK_THROW(file(name, file::in), error);                // read non-existing file
    BOOST_CHECK_NO_THROW(file(name, file::out | file::excl));      // create new file
    unlink(name);                                                  // remove file ...
    BOOST_CHECK_NO_THROW(file(name, file::excl));                  // ... and create again
    BOOST_CHECK_THROW(file(name, file::out | file::excl), error);  // don't overwrite existing file
    unlink(name);

    // check conflicting modes, should not create a file
    BOOST_CHECK_THROW(file(name, file::trunc | file::excl), error);             // "trunc" and "excl" are conflicting
    BOOST_CHECK_THROW(file(name, file::out | file::trunc | file::excl), error); // "trunc" and "excl" are conflicting

    // file should not exist here
    BOOST_CHECK(is_hdf5_file(name) < 0);
}

// test copying and moving
BOOST_AUTO_TEST_CASE( copy_move )
{
    file foo(name);

    BOOST_CHECK(foo.valid());

    BOOST_CHECK_THROW(file f(foo), h5xx::error);   // copying is not allowed
    BOOST_CHECK(foo.valid());

    BOOST_CHECK_NO_THROW(file f(move(foo)));       // copying from temporary is allowed (with move semantics)
    BOOST_CHECK(!foo.valid());

    foo.open(name);
    hid_t hid = foo.hid();
    file bar;
    BOOST_CHECK_THROW(bar = foo, h5xx::error);     // assignment is not allowed (it makes a copy)
    BOOST_CHECK_NO_THROW(bar = move(foo));         // move assignment
    BOOST_CHECK_EQUAL(bar.hid(), hid);
    BOOST_CHECK(!foo.valid());

    bar.close();
    unlink(name);
}

// test hdf5_id
BOOST_AUTO_TEST_CASE( hdf5_id )
{
    // create file and test HDF5 object ID
    BOOST_CHECK_EQUAL(file().hid(), -1);
    BOOST_CHECK(file(name).hid() >= 0);
    BOOST_CHECK(is_hdf5_file(name) > 0);
    unlink(name);
}
