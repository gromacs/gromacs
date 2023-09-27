/**
 * Copyright © 2013-2014 Felix Höfling
 * Copyright © 2014      Klaus Reuter
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#include <h5xx/h5xx.hpp>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <iostream>
#include <vector>
#include <cstdio>

const int NI=10;
const int NJ=NI;

typedef boost::array<int, NI> array_t;
typedef boost::multi_array<int, 2> array_2d_t;

void print_array(array_2d_t const& array)
{
    for (unsigned int j = 0; j < array.shape()[1]; j++)
    {
        for (unsigned int i = 0; i < array.shape()[0]; i++)
        {
            printf("%2d ", array[j][i]);
        }
        printf("\n");
    }
}

template <typename ArrayT>
void print_array(ArrayT const& array)
{
    for (unsigned int i = 0; i < array.size(); i++)
    {
        printf("%2d ", array[i]);
    }
    printf("\n");
}


// --- run some tests of the string slicing notation
void write_int_data(std::string const& filename, array_t const& array)
{
    h5xx::file file(filename, h5xx::file::out);
    std::string name;

    {
        // --- create dataset and fill it with the default array data (positive values)
        name = "integer array";
        h5xx::create_dataset(file, name, array);
        h5xx::write_dataset(file, name, array);

        // --- create a slice object (aka hyperslab) in the conventional way
        boost::array<int,1> offset = {{4}};
        boost::array<int,1> count = {{2}};
        h5xx::slice slice(offset, count);

        // --- data to be written to the slice (negative values)
        boost::array<int,2> data = {{-1,-2}};

        // --- overwrite part of the dataset as specified by slice
        h5xx::write_dataset(file, name, data, slice);
    }


    // --- in the code blocks below, we test the Python/NumPy-like slicing notation ---
    {
        name = "integer array 2";
        h5xx::create_dataset(file, name, array);
        h5xx::write_dataset(file, name, array);

        h5xx::slice slice("4:5");

        boost::array<int,1> data = {{-1}};

        h5xx::write_dataset(file, name, data, slice);
    }

    {
        name = "integer array 3";
        h5xx::create_dataset(file, name, array);
        h5xx::write_dataset(file, name, array);

        h5xx::slice slice("2:6:2");

        boost::array<int,2> data = {{-1,-2}};

        h5xx::write_dataset(file, name, data, slice);
    }

    {
        name = "integer array 4";
        h5xx::create_dataset(file, name, array);
        h5xx::write_dataset(file, name, array);

        h5xx::slice slice(":");

        boost::array<int,NI> data = {{-1,-2,-3,-4,-5,-6,-7,-8,-9,0}};

        h5xx::write_dataset(file, name, data, slice);
    }

    {
        name = "integer array 5";
        h5xx::create_dataset(file, name, array);
        h5xx::write_dataset(file, name, array);

        h5xx::slice slice(":5"); // select all elements up to index 5, ie 0,1,2,3,4

        boost::array<int,5> data = {{-1,-2,-3,-4,-5}};

        h5xx::write_dataset(file, name, data, slice);
    }

    {
        name = "integer array 6";
        h5xx::create_dataset(file, name, array);
        h5xx::write_dataset(file, name, array);

        h5xx::slice slice("5:"); // select all elements starting from index 5, ie 5,6,7,8,9

        boost::array<int,5> data = {{-1,-2,-3,-4,-5}};

        h5xx::write_dataset(file, name, data, slice);
    }

    {
        name = "integer array 7";
        h5xx::create_dataset(file, name, array);
        h5xx::write_dataset(file, name, array);

        h5xx::slice slice("::2"); // select every second element, ie 0,2,4,6,8

        boost::array<int,5> data = {{-1,-2,-3,-4,-5}};

        h5xx::write_dataset(file, name, data, slice);
    }

    {
        name = "integer array 8";
        h5xx::create_dataset(file, name, array);
        h5xx::write_dataset(file, name, array);

        h5xx::slice slice("::3"); // select every second element, ie 0,3,6,9

        boost::array<int,4> data = {{-1,-2,-3,-4}};

        h5xx::write_dataset(file, name, data, slice);
    }
}



void write_2D_int_data(std::string const& filename, array_2d_t const& array)
{
    h5xx::file file(filename, h5xx::file::out);
    std::string name;

    // --- write a 2x2 patch in the 2D dataset using a slice/hyperslab
    {
        name = "2D integer array";
        h5xx::create_dataset(file, name, array);
        h5xx::write_dataset(file, name, array);

//        std::vector<int> offset; int offset_raw[2] = {4,4}; offset.assign(offset_raw, offset_raw + 2);
//        std::vector<int> count;  int count_raw[2] = {2,2}; count.assign(count_raw, count_raw + 2);
//        h5xx::slice slice(offset, count);
        h5xx::slice slice("4:6,4:6");

        boost::array<int,4> data = {{-1,-2,-3,-4}};

        h5xx::write_dataset(file, name, data, slice);
    }

    // --- overwrite part of the third row of the 2D dataset using a slice/hyperslab
    {
        name = "2D integer array 1";
        h5xx::create_dataset(file, name, array);
        h5xx::write_dataset(file, name, array);

        h5xx::slice slice("2,4:8");

        boost::array<int,4> data = {{-1,-2,-3,-4}};

        h5xx::write_dataset(file, name, data, slice);
    }

    // --- overwrite part of the third column of the 2D dataset using a slice/hyperslab
    {
        name = "2D integer array 2";
        h5xx::create_dataset(file, name, array);
        h5xx::write_dataset(file, name, array);

        h5xx::slice slice("4:8,2");

        boost::array<int,4> data = {{-1,-2,-3,-4}};

        h5xx::write_dataset(file, name, data, slice);
    }

    // --- overwrite the fifth row of the 2D dataset using a slice/hyperslab
    {
        name = "2D integer array 3";
        h5xx::create_dataset(file, name, array);
        h5xx::write_dataset(file, name, array);

        h5xx::slice slice("4,:");

        boost::array<int,10> data = {{-1,-2,-3,-4,-5,-6,-7,-8,-9,0}};

        h5xx::write_dataset(file, name, data, slice);
    }

    // --- overwrite the fifth column of the 2D dataset using a slice/hyperslab
    {
        name = "2D integer array 4";
        h5xx::create_dataset(file, name, array);
        h5xx::write_dataset(file, name, array);

        h5xx::slice slice(":,4");

        boost::array<int,10> data = {{-1,-2,-3,-4,-5,-6,-7,-8,-9,0}};

        h5xx::write_dataset(file, name, data, slice);
    }

    // --- overwrite a strided pattern of the 2D dataset using a slice/hyperslab
    {
        name = "2D integer array 5";
        h5xx::create_dataset(file, name, array);
        h5xx::write_dataset(file, name, array);

        h5xx::slice slice("4:6,2:7:2");

        boost::array<int,6> data = {{-1,-2,-3,-4,-5,-6}};

        h5xx::write_dataset(file, name, data, slice);
    }
}



void read_int_data(std::string const& filename)
{
    h5xx::file file(filename, h5xx::file::in);
    std::string name = "integer array";

    // read and print the full dataset
    {
        array_t data;
        // --- read the complete dataset into data
        h5xx::read_dataset(file, name, data);
        printf("original integer array read from file, negative number patch was written using a slice\n");
        print_array(data);
        printf("\n");
    }

    // read and print a subset of the dataset
    {
        boost::array<int,6> data;

        // --- create a slice object (aka hyperslab) to specify the patch to be read from the dataset
//        boost::array<int,1> offset; offset[0] = 2;
//        boost::array<int,1> count; count[0] = 6;
//        h5xx::slice slice(offset, count);
        h5xx::slice slice("2:8");

        h5xx::read_dataset(file, name, data, slice);

        printf("1D slice of the integer array, zoom on the negative number patch\n");
        print_array(data);
        printf("\n");
    }
}


int main(int argc, char** argv)
{
    std::string filename = argv[0];
    filename.append(".h5");

    {
        h5xx::file file(filename, h5xx::file::trunc);
    }

    // --- do a few 1D tests using integers
    {
        array_t array;
        for (int i = 0; i < NI; i++)
            array[i] = i;

        write_int_data(filename, array);

        read_int_data(filename);
    }

    // --- 2D tests
    {
        array_2d_t array(boost::extents[NJ][NI]);
        {
            const int nelem = NI*NJ;
            int data[nelem];
            for (int i = 0; i < nelem; i++)
                data[i] = i;
            array.assign(data, data + nelem);
        }

        write_2D_int_data(filename, array);
    }

    return 0;
}
