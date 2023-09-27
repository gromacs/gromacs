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
#include <iostream>
#include <cstdio>

const int NI=10;

void print_array(std::vector<int> const& array)
{
    for (unsigned int i = 0; i < array.size(); i++)
    {
        printf("%2d ", array[i]);
    }
    printf("\n");
}


void write_int_data(std::string const& filename, std::vector<int> const& array)
{
    h5xx::file file(filename, h5xx::file::trunc);
    std::string name;

    {
        // --- create dataset and fill it with the default array data (positive values)
        name = "integer array";
        h5xx::create_dataset(file, name, array);
        h5xx::write_dataset(file, name, array);

        // --- create a slice object (aka hyperslab) to specify the location in the dataset to be overwritten
        std::vector<int> offset; int offset_raw[2] = {4}; offset.assign(offset_raw, offset_raw + 1);
        std::vector<int> count;  int count_raw[2] = {2}; count.assign(count_raw, count_raw + 1);
        h5xx::slice slice(offset, count);

        // --- data to be written to the slice (negative values)
        std::vector<int> data;
        data.push_back(-1);
        data.push_back(-2);

        // --- overwrite part of the dataset as specified by slice
        h5xx::write_dataset(file, name, data, slice);
    }
}


void read_int_data(std::string const& filename)
{
    h5xx::file file(filename, h5xx::file::in);
    std::string name = "integer array";

    // read and print the full dataset
    {
        std::vector<int> data;
        // --- read the complete dataset into data, the vector is resized and overwritten internally
        h5xx::read_dataset(file, name, data);
        printf("original integer array read from file, negative number patch was written using a slice\n");
        print_array(data);
        printf("\n");
    }

    // read and print a subset of the dataset
    {
        std::vector<int> data;

        // --- create a slice object (aka hyperslab) to specify the patch to be read from the dataset
        std::vector<int> offset; int offset_raw[2] = {2}; offset.assign(offset_raw, offset_raw + 1);
        std::vector<int> count;  int count_raw[2] = {6}; count.assign(count_raw, count_raw + 1);
        h5xx::slice slice(offset, count);

        // --- allocate memory for the slice (user's responsibility, is not done by read_dataset if slicing is used)
        data.resize(count[0]);

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

    // --- do a few demos/tests using integers
    {
        std::vector<int> array;
        for (int i = 0; i < NI; i++)
            array.push_back(i);

        write_int_data(filename, array);

        read_int_data(filename);
    }

    return 0;
}
