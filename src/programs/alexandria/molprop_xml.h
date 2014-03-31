/*
 * This source file is part of the Alexandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef MOLPROP_XML_H
#define MOLPROP_XML_H

#include "molprop.h"

/*! \brief
 * Write a vector of molprops to an XML file 
 *
 * Function that uses the libxml2 library to store an array of molprops
 * to file.
 *
 * \todo Implement using a serialized protocol rather than reading the
 * whole file into memory. Since we have a vector of molprops this should
 * be straightforward. 
 * \param[in] fn         The file name to write to
 * \param[in] mpt        The vector of MolProp
 * \param[in] bCompress  Determines whether zlib compression is used when writing
 * \ingroup module_alexandria
 */
void MolPropWrite(const char *fn,
                  const std::vector<alexandria::MolProp> mpt,
                  gmx_bool bCompress);
  
/*! \brief
 * Reads a vector of molprops from an XML file 
 *
 * Function that uses the libxml2 library to store an array of molprops
 * to file.
 *
 * \todo Implement using a serialized protocol rather than reading the
 * whole file into memory. Since we have a vector of molprops this should
 * be straightforward. 
 * \param[in]  fn         The file name to read from
 * \param[out] mpt        The vector of MolProp
 * \ingroup module_alexandria
 */
void MolPropRead(const char *fn,
                 std::vector<alexandria::MolProp>& mpt);

#endif
