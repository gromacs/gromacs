/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef GAUSS_IO_H
#define GAUSS_IO_H

#include "molprop.h"

/*! \brief
 * Read a Gaussian log file either using home grown methods or using OpenBabel
 *
 *
 * \param[in] g98        The gaussian log file, or in case OpenBabel is used anything
 *                       that can be read by OpenBabel
 * \param[out] mpt       The MolProp
 * \param[in] molnm      Molecule name to override the one from the filename [ maybe NULL ]
 * \param[in] iupac      IUPAC name to override the one from the filename [ maybe NULL ]
 * \param[in] conformation  Conformation the molecule is in [ maybe NULL ]
 * \param[in] basisset   Basis set used for the calculation [ maybe NULL ]
 * \param[in] maxpot     Maximum number of electrostatic potential data points to store
 * \param[in] nsymm      Symmetry number for this molecule. If zero it will be detected from
 *                       the input.
 * \param[in] forcefield One of the force fields supported by OpenBabel used for atomtypes
 * \ingroup module_alexandria
 */
void ReadGauss(const char          *g98,
               alexandria::MolProp &mp,
               const char          *molnm,
               const char          *iupac,
               const char          *conf,
               const char          *basis,
               int                  maxpot,
               int                  nsymm,
               const char          *forcefield,
               const char          *jobtype);

#endif
