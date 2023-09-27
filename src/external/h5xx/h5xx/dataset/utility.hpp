/*
 * Copyright © 2014-2015 Klaus Reuter
 * Copyright © 2014      Manuel Dibak
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#ifndef H5XX_DATASET_UTILITY_HPP
#define H5XX_DATASET_UTILITY_HPP

#include <h5xx/error.hpp>
#include <h5xx/utility.hpp>

namespace h5xx {

/**
 * Check whether a dataset of the given name is attached to the h5xx object.
 */
template <typename h5xxObject>
inline bool exists_dataset(h5xxObject const& object, std::string const& name)
{
    hid_t hid;
    H5E_BEGIN_TRY {
        hid = H5Dopen(object.hid(), name.c_str(), H5P_DEFAULT);
        if (hid > 0) {
            H5Dclose(hid);
        }
    } H5E_END_TRY
    return (hid > 0);
}

} // namespace h5xx

#endif /* ! H5XX_DATASET_UTILITY_HPP */
