/*
 * Copyright © 2010  Peter Colberg
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#ifndef H5XX_ERROR_HPP
#define H5XX_ERROR_HPP

#include <stdexcept>

namespace h5xx {

/**
 * h5xx wrapper error
 */
class error
  : virtual public std::exception
{
public:
    error(std::string const& desc)
        : desc_(desc) {}

    virtual ~error() throw() {}

    char const* what() const throw()
    {
        return desc_.c_str();
    }

private:
    std::string desc_;
};

} // namespace h5xx

#endif /* ! H5XX_ERROR_HPP */
