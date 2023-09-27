/*
 * Copyright © 2013-2014 Felix Höfling
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#include <boost/multi_array.hpp>
#include <h5xx/h5xx.hpp>
#include <iostream>

typedef boost::multi_array<int, 3> array_t;

void write_attribute(std::string const& filename, array_t const& array)
{
    // open group within HDF5 file
    h5xx::file file(filename, h5xx::file::trunc); // truncate existing file
    h5xx::group group(file, "group");

    // attach string attribute to file root
    h5xx::write_attribute(file, "location", "Here is the file root.");

    // attach array data as attribute to group
    h5xx::write_attribute(group, "integer array", array);
}

void read_attribute(std::string const& filename)
{
    // open HDF5 file read-only
    h5xx::file file(filename, h5xx::file::in);

    // read and print string attribute
    std::cout << h5xx::read_attribute<std::string>(file, "location") << std::endl;

    // read array data from attribute
    array_t array = h5xx::read_attribute<array_t>(h5xx::group(file, "group"), "integer array");

    // close file explicitly
    file.close();

    // print array along second dimension
    for (unsigned int i = 0; i < array.shape()[1]; ++i) {
        std::cout << "/group/integer array[0, " << i << ", 0] = " << array[0][i][0] << std::endl;
    }

    // TODO read dataset of rank 3 by iterating over the first index (see master branch)
}


int main(int argc, char** argv)
{
    std::string filename = argv[0];
    filename.append(".h5");

    // set-up data as Boost.MultiArray
    array_t array(boost::extents[2][3][2]);
    int data[] = {99, 98, 97, 96, 95, 94, 93, 92, 91, 90, 89, 88};
    array.assign(data, data + sizeof(data) / sizeof(int));

    // write to HDF5 file
    write_attribute(filename, array);

    // read from HDF5 file
    read_attribute(filename);

    return 0;
}
