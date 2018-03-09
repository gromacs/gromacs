/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
 
 
#ifndef POLDATA_XML_H
#define POLDATA_XML_H

#include <string>

#include "gromacs/topology/atomprop.h"

namespace alexandria
{
    class Poldata;

    /*! \brief Store the Poldata class to an XML file
     *
     * \param[in] fileName The filename to save to
     * \param[in] pd       Pointer to a Poldata class instance
     * \param[in] compress Whether or not to write a compressed file
     */
    void writePoldata(const std::string &fileName,
                      const Poldata &pd,
                      bool compress = true);

    /*! \brief Read a Poldata class from an XML file
     *
     * \param[in]  fileName The filename to read from
     * \param[out] pd       The Poldata class instance
     * \param[in]  aps      Atom properties
     */
    void readPoldata(const std::string &fileName,
                     Poldata &pd,
                     const gmx_atomprop_t aps);

} // namespace alexandria

#endif
