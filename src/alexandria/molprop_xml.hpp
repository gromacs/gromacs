/*
 * $Id: molprop_xml.hpp,v 1.3 2009/04/05 11:46:57 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */

#ifndef _molprop_xml_hpp
#define _molprop_xml_hpp

#include "molprop.hpp"

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
