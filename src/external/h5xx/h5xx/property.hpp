/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#ifndef H5XX_PROPERTY_HPP
#define H5XX_PROPERTY_HPP

#include <h5xx/error.hpp>
#include <h5xx/hdf5_compat.hpp>

namespace h5xx {

///**
// * Create link creation property list and set create intermediate group property.
// */
//inline H5::PropList create_intermediate_group_property()
//{
//    hid_t pl = H5Pcreate(H5P_LINK_CREATE);
//    if (pl < 0) {
//        throw error("failed to create link creation property list");
//    }
//    herr_t err = H5Pset_create_intermediate_group(pl, 1);
//    if (err < 0) {
//        throw error("failed to set group intermediate creation property");
//    }
//    return H5::PropList(pl);
//}

} // namespace h5xx

#endif /* ! H5XX_PROPERTY_HPP */
