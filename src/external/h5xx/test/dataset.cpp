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

char filename[] = "test_h5xx_dataset.h5";
typedef h5file<filename> BOOST_AUTO_TEST_CASE_FIXTURE;

typedef boost::multi_array<int, 1> array_1d_t;
typedef boost::multi_array<int, 2> array_2d_t;


BOOST_AUTO_TEST_CASE( construction )
{
    BOOST_CHECK_NO_THROW(dataset()); // default constructor
    BOOST_CHECK_NO_THROW(create_dataset<int>(file, "foo"));
    BOOST_CHECK_NO_THROW(write_dataset(file, "foo", 1)); // create dataset in a file's root group
    BOOST_CHECK_NO_THROW(dataset(file, "foo")); // open existing attribute on-the-fly

    dataset foo(file, "foo");
    BOOST_CHECK_EQUAL(get_name(foo), "/foo"); // full path of the dataset
    BOOST_CHECK(foo.valid());

    hid_t hid = foo.hid();
    dataset bar;
    BOOST_CHECK_THROW(bar = foo, h5xx::error); // assignment is not allowed (it makes a copy)
    BOOST_CHECK_NO_THROW(bar = move(foo)); // move assignment
    BOOST_CHECK_EQUAL(bar.hid(), hid);
    BOOST_CHECK(!foo.valid());

    BOOST_CHECK_THROW(dataset g(bar), h5xx::error); // copying is not allowed
    BOOST_CHECK_NO_THROW(dataset g(move(bar))); // copying from temporary is allowed (with move semantics)
    BOOST_CHECK(!bar.valid());
}

BOOST_AUTO_TEST_CASE( scalar_fundamental )
{
    bool bool_value = true;
    std::string bool_name = "bool, scalar";
    BOOST_CHECK_NO_THROW(create_dataset<bool>(file, bool_name));
    BOOST_CHECK_NO_THROW(write_dataset(file, bool_name, bool_value));
    BOOST_CHECK_NO_THROW(read_dataset<bool>(file, bool_name));
    BOOST_CHECK(read_dataset<bool>(file, bool_name) == bool_value);
    BOOST_CHECK(exists_dataset(file, bool_name) == true);
    H5E_BEGIN_TRY{
        BOOST_CHECK_THROW(read_attribute<bool>(file, "X"+bool_name), h5xx::error);
    } H5E_END_TRY

    double double_value = std::sqrt(2.L);
    std::string double_name = "double, scalar";
//    BOOST_CHECK_NO_THROW(
//            write_dataset(file, double_name, bool_value)  // cannot change datatype of dataset
//    );
    BOOST_CHECK_NO_THROW(create_dataset<double>(file, double_name));
    BOOST_CHECK_NO_THROW(write_dataset(file, double_name, double_value));
    BOOST_CHECK_NO_THROW(write_dataset(file, double_name, 0.5*double_value));
    BOOST_CHECK_NO_THROW(read_dataset<double>(file, double_name));
    BOOST_CHECK(exists_dataset(file, double_name) == true);

    uint64_t uint64_value = 9223372036854775783LLU;  // largest prime below 2^63
    std::string uint64_name = "uint64, scalar";
    BOOST_CHECK_NO_THROW(create_dataset<uint64_t> (file, uint64_name));
    BOOST_CHECK_NO_THROW(write_dataset(file, uint64_name, uint64_value));
    BOOST_CHECK_NO_THROW(read_dataset<uint64_t>(file, uint64_name));
    BOOST_CHECK(exists_dataset(file, uint64_name) == true);
}

// test default creation of a dataset (no storage policy provided)
BOOST_AUTO_TEST_CASE( boost_multi_array_simple )
{
    typedef boost::multi_array<int, 3> multi_array3;
    int data3[] = {99,98,97,96,95,94,93,92,91,90,89,88,87,86,85,84,83,82,81,80,79,78,77,76};
    multi_array3 multi_array_value(boost::extents[2][3][4]);
    multi_array_value.assign(data3, data3 + 2 * 3 * 4);
    multi_array3 arrayRead(boost::extents[2][3][4]);
    const std::string name = "boost multi array, int, default";
    BOOST_CHECK_NO_THROW(create_dataset(file, name, multi_array_value));
    BOOST_CHECK_NO_THROW(write_dataset(file, name, multi_array_value));
    BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
    BOOST_CHECK(arrayRead == multi_array_value);
}

// test chunked dataset and the filters (compression, etc) it can use
BOOST_AUTO_TEST_CASE( boost_multi_array_chunked )
{
    const int NI=10;
    const int NJ=NI;
    boost::array<size_t, 2> chunkDims = {{2,2}};
    std::string name;
    array_2d_t arrayRead(boost::extents[NJ][NI]);
    array_2d_t arrayWrite(boost::extents[NJ][NI]);
    {
        const int nelem = NI*NJ;
        int data[nelem];
        for (int i = 0; i < nelem; i++) data[i] = i;
        arrayWrite.assign(data, data + nelem);
    }

    {
        name = "boost multi array, int, chunked";
        h5xx::policy::storage::chunked storagePolicy(chunkDims);
        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
        BOOST_CHECK(arrayRead == arrayWrite);
        zero_multi_array(arrayRead);
    }

    {
        name = "boost multi array, int, chunked, deflate";
        h5xx::policy::storage::chunked storagePolicy(chunkDims);
        storagePolicy.add(h5xx::policy::filter::deflate());
        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
        BOOST_CHECK(arrayRead == arrayWrite);
        zero_multi_array(arrayRead);
    }

    {
        name = "boost multi array, int, chunked, szip";
        h5xx::policy::storage::chunked storagePolicy(chunkDims);
        //storagePolicy.add(h5xx::policy::filter::szip());  // most (?) HDF5 builds do not support SZIP
        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
        BOOST_CHECK(arrayRead == arrayWrite);
        zero_multi_array(arrayRead);
    }

    {
        name = "boost multi array, int, chunked, shuffle";
        h5xx::policy::storage::chunked storagePolicy(chunkDims);
        storagePolicy.add(h5xx::policy::filter::shuffle());
        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
        BOOST_CHECK(arrayRead == arrayWrite);
        zero_multi_array(arrayRead);
    }

    {
        name = "boost multi array, int, chunked, fletcher32";
        h5xx::policy::storage::chunked storagePolicy(chunkDims);
        storagePolicy.add(h5xx::policy::filter::fletcher32());
        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
        BOOST_CHECK(arrayRead == arrayWrite);
        zero_multi_array(arrayRead);
    }

    {
        name = "boost multi array, int, chunked, scaleoffset";
        h5xx::policy::storage::chunked storagePolicy(chunkDims);
        storagePolicy.add(h5xx::policy::filter::scaleoffset<int>(0));
        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
        BOOST_CHECK(arrayRead == arrayWrite);
        zero_multi_array(arrayRead);
    }

    {
        name = "boost multi array, int, chunked, nbit";
        h5xx::policy::storage::chunked storagePolicy(chunkDims);
        storagePolicy.add(h5xx::policy::filter::nbit());
        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
        BOOST_CHECK(arrayRead == arrayWrite);
        zero_multi_array(arrayRead);
    }
}


// test contiguous dataset and some storage policies it can use
BOOST_AUTO_TEST_CASE( boost_multi_array_contiguous )
{
    const int NI=10;
    const int NJ=NI;
    std::string name;
    array_2d_t arrayRead(boost::extents[NJ][NI]);
    array_2d_t arrayWrite(boost::extents[NJ][NI]);
    {
        const int nelem = NI*NJ;
        int data[nelem];
        for (int i = 0; i < nelem; i++) data[i] = i;
        arrayWrite.assign(data, data + nelem);
    }

    {
        name = "boost multi array, int, contiguous";
        h5xx::policy::storage::contiguous storagePolicy;
        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
        BOOST_CHECK(arrayRead == arrayWrite);
        zero_multi_array(arrayRead);
    }

    {
        name = "boost multi array, int, contiguous, fill_value";
        const int fillValue = 4711;
        h5xx::policy::storage::contiguous storagePolicy;
        storagePolicy.set(h5xx::policy::storage::fill_value(fillValue));
        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
//        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
        array_2d_t array4711(boost::extents[NJ][NI]);
        for (unsigned int i = 0; i < array4711.num_elements(); i++)
            array4711.data()[i] = fillValue;
        BOOST_CHECK(arrayRead == array4711);
        zero_multi_array(arrayRead);
    }

    {
        name = "boost multi array, int, contiguous, track_times";
        h5xx::policy::storage::contiguous storagePolicy;
        storagePolicy.set(h5xx::policy::storage::track_times());
        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
        BOOST_CHECK(arrayRead == arrayWrite);
        zero_multi_array(arrayRead);
    }
}


// test compact dataset and some storage policies it can use
BOOST_AUTO_TEST_CASE( boost_multi_array_compact )
{
    const int NI=10;
    const int NJ=NI;
    std::string name;
    array_2d_t arrayRead(boost::extents[NJ][NI]);
    array_2d_t arrayWrite(boost::extents[NJ][NI]);
    {
        const int nelem = NI*NJ;
        int data[nelem];
        for (int i = 0; i < nelem; i++) data[i] = i;
        arrayWrite.assign(data, data + nelem);
    }

    {
        name = "boost multi array, int, compact";
        h5xx::policy::storage::compact storagePolicy;
        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
        BOOST_CHECK(arrayRead == arrayWrite);
        zero_multi_array(arrayRead);
    }

    {
        name = "boost multi array, int, compact, fill_value";
        const int fillValue = 4711;
        h5xx::policy::storage::compact storagePolicy;
        storagePolicy.set(h5xx::policy::storage::fill_value(fillValue));
        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
//        BOOST_CHECK_NO_THROW(
//                write_dataset(file, name, arrayWrite)
//        );
        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
        array_2d_t array4711(boost::extents[NJ][NI]);
        for (unsigned int i = 0; i < array4711.num_elements(); i++)
            array4711.data()[i] = fillValue;
        BOOST_CHECK(arrayRead == array4711);
        zero_multi_array(arrayRead);
    }

    {
        name = "boost multi array, int, compact, track_times";
        h5xx::policy::storage::compact storagePolicy;
        storagePolicy.set(h5xx::policy::storage::track_times());
        BOOST_CHECK_NO_THROW(create_dataset(file, name, arrayWrite, storagePolicy));
        BOOST_CHECK_NO_THROW(write_dataset(file, name, arrayWrite));
        BOOST_CHECK_NO_THROW(read_dataset(file, name, arrayRead));
        BOOST_CHECK(arrayRead == arrayWrite);
        zero_multi_array(arrayRead);
    }
}



BOOST_AUTO_TEST_CASE( h5xx_dataset_hyperslab )
{
    const int NI=10;
    const int NJ=NI;
    std::string name;
    array_2d_t arrayRead(boost::extents[NJ][NI]);
    array_2d_t arrayWrite(boost::extents[NJ][NI]);
    {
        const int nelem = NI*NJ;
        int data[nelem];
        for (int i = 0; i < nelem; i++) data[i] = i;
        arrayWrite.assign(data, data + nelem);
    }

    name = "integer array";
    h5xx::create_dataset(file, name, arrayWrite);
    h5xx::write_dataset(file, name, arrayWrite);

    std::vector<int> offset; { int offset_raw[2] = {4,4}; offset.assign(offset_raw, offset_raw + 2); }
    std::vector<int> count;  { int count_raw[2] = {2,2};  count.assign(count_raw, count_raw + 2); }
    h5xx::slice slice(offset, count);

    array_1d_t slice_data(boost::extents[4]);
    {
        int data_raw[4] = {-1,-2,-3,-4};
        slice_data.assign(data_raw, data_raw+4);
    }

    BOOST_CHECK_NO_THROW(h5xx::write_dataset(file, name, slice_data, slice));

    BOOST_CHECK_NO_THROW(h5xx::read_dataset(file, name, arrayRead));

    BOOST_CHECK(arrayRead != arrayWrite);

    // manipulate arrayWrite such that it should be equal to arrayRead at this point
    {
        int neg = -1;
        for (int j = 0; j < count[0]; j++) {
            for (int i = 0; i < count[1]; i++) {
                arrayWrite[ j+offset[0] ][ i+offset[1] ] = neg;
                neg--;
            }
        }
    }

    BOOST_CHECK(arrayRead == arrayWrite);
}

// TODO : add more slicing tests here

} //namespace fixture
