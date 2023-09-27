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

#include <mpi.h>

const int NI=10;
const int NJ=NI;

typedef boost::array<int, NI> array_t;

// --- global MPI variables, for simple use in read,write functions
int mpi_size, mpi_rank;
MPI_Comm comm;
MPI_Info info;
// ---


template <typename ArrayT>
void print_array(ArrayT const& array)
{
    for (unsigned int i = 0; i < array.size(); i++)
    {
        printf("%2d ", array[i]);
    }
    printf("\n");
}


// --- write MPI-rank individual data to the dataset
void write_int_data(std::string const& filename, array_t const& array)
{
    h5xx::file file(filename, comm, info, h5xx::file::out);
    std::string name;

    {
        name = "integer array";
        boost::array<size_t, 1> chunk_dims = {{2}};
        h5xx::create_dataset(file, name, array, h5xx::policy::storage::chunked(chunk_dims));
        h5xx::write_dataset(file, name, array);

        boost::array<size_t,1> offset;
        offset[0] = mpi_rank;
        boost::array<size_t,1> count = {{1}};
        h5xx::slice slice(offset, count);

        boost::array<size_t,1> data;
        data[0] = mpi_rank;

        h5xx::write_dataset(file, name, data, slice);
    }

}




void read_int_data(std::string const& filename)
{
    h5xx::file file(filename, h5xx::file::in);
    std::string name = "integer array";

    {
        array_t data;
        h5xx::read_dataset(file, name, data);
        printf("integer array read from file, numbers <10 were written by separate MPI ranks\n");
        print_array(data);
        printf("\n");
    }
}


int main(int argc, char** argv)
{
    std::string filename = argv[0];
    filename.append(".h5");

    MPI_Init(&argc, &argv);
    // --- initialize MPI variables (global variables, for convenience)
    comm = MPI_COMM_WORLD;
    info = MPI_INFO_NULL;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    if (mpi_size > 10)
    {
        MPI_Finalize();
        return 1;
    } else {
        {
            h5xx::file file(filename, comm, info, h5xx::file::trunc);
        }
        {
            array_t array;
            for (int i = 0; i < NI; i++)
                array[i] = 10+i;

            write_int_data(filename, array);

            if (mpi_rank == 0)
                read_int_data(filename);
        }
        MPI_Finalize();
        return 0;
    }
}
