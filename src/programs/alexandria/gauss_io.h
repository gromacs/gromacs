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
 
 
#ifndef GAUSS_IO_H
#define GAUSS_IO_H

#include "molprop.h"


enum BabelFileType {
    ebfPDB  = 0,
    ebfXYZ  = 1,
    ebfSDF  = 2,
    ebfMOL  = 3,
    ebfMOL2 = 4,
    ebfG09  = 5,
    ebfNR   = 6
};

/*! \brief
 * Read a Gaussian log file either using home grown methods or using OpenBabel
 *
 *
 * \param[in] g98        The gaussian log file, or in case OpenBabel is used anything
 *                       that can be read by OpenBabel
 * \param[out] mpt       The MolProp
 * \param[in] molnm      Molecule name to override the one from the filename [ maybe nullptr ]
 * \param[in] iupac      IUPAC name to override the one from the filename [ maybe nullptr ]
 * \param[in] conformation  Conformation the molecule is in [ maybe nullptr ]
 * \param[in] basisset   Basis set used for the calculation [ maybe nullptr ]
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
