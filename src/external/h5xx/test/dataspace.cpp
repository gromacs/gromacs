/*
 * Copyright © 2010-2014 Felix Höfling
 * Copyright © 2014-2016 Klaus Reuter
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */


#define BOOST_TEST_MODULE h5xx_dataspace

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
    BOOST_CHECK_NO_THROW(dataspace());

    {
        dataspace ds;
        BOOST_CHECK(ds.valid() == false);
    }

    // --- creation of dataspace from a std::vector
    {
        std::vector<int> vec;
        vec.push_back(2);
        vec.push_back(4);
        vec.push_back(6);
        vec.push_back(8);
        dataspace ds;
        BOOST_CHECK_NO_THROW(ds = create_dataspace(vec));
        BOOST_CHECK(ds.valid() == true);
        std::vector<std::size_t> xts;
        BOOST_CHECK_NO_THROW(xts = to_size_t(ds.extents()));
        BOOST_CHECK_EQUAL(xts.size(), std::size_t(1));
        BOOST_CHECK_EQUAL(xts[0], vec.size());
    }

    // --- creation of dataspace from a boost::array
    {
        boost::array<int, 4> vec = {{2,4,6,8}};
        dataspace ds;
        BOOST_CHECK_NO_THROW(ds = create_dataspace(vec));
        BOOST_CHECK(ds.valid() == true);
        std::vector<std::size_t> xts;
        BOOST_CHECK_NO_THROW(xts = to_size_t(ds.extents()));
        BOOST_CHECK_EQUAL(xts.size(), std::size_t(1));
        BOOST_CHECK_EQUAL(xts[0], vec.size());
    }


    // --- creation of dataspace from a boost::multiarray
    {
        const int NI=10;
        const int NJ=NI;
        array_2d_t arr(boost::extents[NJ][NI]);
        {
            const int nelem = NI*NJ;
            int data[nelem];
            for (int i = 0; i < nelem; i++) data[i] = i;
            arr.assign(data, data + nelem);
        }
        dataspace ds;
        BOOST_CHECK_NO_THROW(ds = create_dataspace(arr));
        BOOST_CHECK(ds.valid() == true);
        std::vector<std::size_t> xts;
        BOOST_CHECK_NO_THROW(xts = to_size_t(ds.extents()));
        BOOST_CHECK_EQUAL(xts.size(), std::size_t(2));
        BOOST_CHECK_EQUAL(xts[0], NI);
        BOOST_CHECK_EQUAL(xts[1], NJ);
    }


    {
        std::vector<int> vec;
        vec.push_back(2);
        vec.push_back(4);
        vec.push_back(6);
        dataspace ds;
        BOOST_CHECK_NO_THROW(ds = create_dataspace(vec));
        BOOST_CHECK(ds.valid() == true);
        BOOST_CHECK_EQUAL(ds.rank(), 1);
        BOOST_CHECK(ds.is_scalar() == false);
        BOOST_CHECK(ds.is_simple() == true);
    }


    // --- test basic select facilities
    {
        const int NI=10;
        const int NJ=NI;
        array_2d_t arr(boost::extents[NJ][NI]);
        {
            const int nelem = NI*NJ;
            int data[nelem];
            for (int i = 0; i < nelem; i++) data[i] = i;
            arr.assign(data, data + nelem);
        }
        dataspace ds;
        BOOST_CHECK_NO_THROW(ds = create_dataspace(arr));
        BOOST_CHECK(ds.valid() == true);
        // --- everything is selected by default
        BOOST_CHECK(ds.get_select_npoints() == NI*NJ);
        // --- define a range by using the numpy-arange-like specification
        {
            slice slc("1:3,3:5");  // 2x2 patch
            BOOST_CHECK_NO_THROW(ds.select(slc));
        }
        BOOST_CHECK(ds.get_select_npoints() == 2*2);
        {
            slice slc("7,7");  // separate single item
            BOOST_CHECK_NO_THROW(ds.select(slc, dataspace::OR));  // add slice to selection, "OR" works like a binary OR
        }
        //std::cout << ds.get_select_npoints() << std::endl;
        BOOST_CHECK(ds.get_select_npoints() == 2*2 + 1);
        {
            slice slc("2:4,4:6");
            BOOST_CHECK_NO_THROW(ds.select(slc, dataspace::AND));  // add slice to selection, "AND" works like a binary AND
        }
        // std::cout << ds.get_select_npoints() << std::endl;
        BOOST_CHECK(ds.get_select_npoints() == 1);
    }

}


} //namespace fixture
