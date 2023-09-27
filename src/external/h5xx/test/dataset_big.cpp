/*
 * Copyright © 2010-2016 Felix Höfling
 * Copyright © 2014-2016 Klaus Reuter
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#define BOOST_TEST_MODULE h5xx_dataset
#include <boost/test/unit_test.hpp>
#include <boost/shared_ptr.hpp>

#include <h5xx/h5xx.hpp>
#include <test/ctest_full_output.hpp>
#include <test/catch_boost_no_throw.hpp>
#include <test/fixture.hpp>

#include <unistd.h>
#include <cstdio>
#include <cmath>
#include <string>

using namespace h5xx;

BOOST_GLOBAL_FIXTURE( ctest_full_output );

namespace fixture { // preferred over BOOST_FIXTURE_TEST_SUITE

template <typename T>
void zero_multi_array(T &array) {
    for (unsigned int i = 0; i < array.num_elements(); i++)
        array.data()[i] = 0;
}

char filename[] = "test_h5xx_dataset_big.h5";
typedef h5file<filename> BOOST_AUTO_TEST_CASE_FIXTURE;

typedef boost::multi_array<int, 1> array_1d_t;
typedef boost::multi_array<int, 2> array_2d_t;


// big test case, chunked dataset
BOOST_AUTO_TEST_CASE( boost_multi_array_chunked_big )
{
    const int NI=10000;
    const int NJ=NI;
    boost::array<size_t, 2> chunkDims = {{32,32}};
    std::string name;
    array_2d_t arrayRead(boost::extents[NJ][NI]);
    array_2d_t arrayWrite(boost::extents[NJ][NI]);
    {
        const int nelem = NI*NJ;
        for (int i = 0; i < nelem; i++)
            arrayWrite.data()[i] = i;
    }

    {
        name = "boost multi array, big, int, chunked";
        h5xx::policy::storage::chunked storagePolicy(chunkDims);
        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
        BOOST_CHECK(arrayRead == arrayWrite);
        zero_multi_array(arrayRead);
    }

    {
        name = "boost multi array, big, int, chunked, deflate";
        h5xx::policy::storage::chunked storagePolicy(chunkDims);
        storagePolicy.add(h5xx::policy::filter::deflate());
        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
        BOOST_CHECK(arrayRead == arrayWrite);
        zero_multi_array(arrayRead);
    }

//    {
//        name = "boost multi array, int, chunked, szip";
//        h5xx::policy::storage::chunked storagePolicy(chunkDims);
//        //storagePolicy.add(h5xx::policy::filter::szip());  // most (?) HDF5 builds do not support SZIP
//        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
//        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
//        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
//        BOOST_CHECK(arrayRead == arrayWrite);
//        zero_multi_array(arrayRead);
//    }
//
//    {
//        name = "boost multi array, int, chunked, shuffle";
//        h5xx::policy::storage::chunked storagePolicy(chunkDims);
//        storagePolicy.add(h5xx::policy::filter::shuffle());
//        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
//        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
//        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
//        BOOST_CHECK(arrayRead == arrayWrite);
//        zero_multi_array(arrayRead);
//    }
//
//    {
//        name = "boost multi array, int, chunked, fletcher32";
//        h5xx::policy::storage::chunked storagePolicy(chunkDims);
//        storagePolicy.add(h5xx::policy::filter::fletcher32());
//        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
//        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
//        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
//        BOOST_CHECK(arrayRead == arrayWrite);
//        zero_multi_array(arrayRead);
//    }
//
//    {
//        name = "boost multi array, int, chunked, scaleoffset";
//        h5xx::policy::storage::chunked storagePolicy(chunkDims);
//        storagePolicy.add(h5xx::policy::filter::scaleoffset<int>(0));
//        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
//        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
//        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
//        BOOST_CHECK(arrayRead == arrayWrite);
//        zero_multi_array(arrayRead);
//    }
//
//    {
//        name = "boost multi array, int, chunked, nbit";
//        h5xx::policy::storage::chunked storagePolicy(chunkDims);
//        storagePolicy.add(h5xx::policy::filter::nbit());
//        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
//        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
//        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
//        BOOST_CHECK(arrayRead == arrayWrite);
//        zero_multi_array(arrayRead);
//    }
}


// TODO : add more slicing tests here

} //namespace fixture
