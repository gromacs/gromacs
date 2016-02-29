/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
/*! \file
 * \brief
 * Provides routines for computing Coulomb integrals analytically.
 * The Slater based code depends on the optional CLN library for arbitrary
 * precision arithmetic.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inpublicapi
 * \ingroup module_alexandria
 */
#ifndef _COULOMBINTEGRALS_H
#define _COULOMBINTEGRALS_H

#define SLATER_MAX 3

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief 
 * Compute the Coulomb overlap integral between two gaussian distributed charges.
 * The widths xi and xj may both be zero.
 *
 * \param[in] r  distance in nm
 * \param[in] xi gaussian width in 1/nm 
 * \param[in] xj gaussian width in 1/nm 
 * \return    Integral value
 */    
double Coulomb_GG(double r,double xi,double xj);

/*! \brief 
 * Compute the Coulomb overlap integral between a gaussian distributed charge
 * and a point charge. The width xi may be zero.
 *
 * \param[in] r  distance in nm
 * \param[in] xi gaussian width in 1/nm 
 * \return    Integral value
 */    
double Nuclear_GG(double r,double xi);

/*! \brief 
 * Compute the derivative of the Coulomb overlap integral between two gaussian 
 * distributed charges with respect to the distance.
 * The widths xi and xj may both be zero.
 *
 * \param[in] r  distance in nm
 * \param[in] xi gaussian width in 1/nm 
 * \param[in] xj gaussian width in 1/nm 
 * \return    Integral value
 */    
double DCoulomb_GG(double r,double xi,double xj);

/*! \brief 
 * Compute the derivative of the Coulomb overlap integral between a gaussian 
 * distributed charge and a point charge. The width xi may be zero.
 *
 * \param[in] r  distance in nm
 * \param[in] xi gaussian width in 1/nm 
 * \return    Integral value
 */    
double DNuclear_GG(double r,double xi);

/*! \brief 
 * Compute the Slater overlap integral between two Slater distributed charges.
 * The widths xi and xj may both be zero. The Slater function types i and j
 * should be in the range 1 .. SLATER_MAX
 *
 * \param[in] r  distance in nm
 * \param[in] i  Slater function type
 * \param[in] j  Slater function type
 * \param[in] xi gaussian width in 1/nm 
 * \param[in] xj gaussian width in 1/nm 
 * \return    Integral value
 */    
double Coulomb_SS(double r,int i,int j,double xi,double xj);

/*! \brief 
 * Compute the Slater overlap integral between a Slater distributed charge and a
 * point charge. The width xi may be zero. The Slater function type i should be in the
 * range 1 .. SLATER_MAX
 *
 * \param[in] r  distance in nm
 * \param[in] i  Slater function type
 * \param[in] xi gaussian width in 1/nm 
 * \return    Integral value
 */    
double Nuclear_SS(double r,int i,double xi);

/*! \brief 
 * Compute the derivative of the Slater overlap integral between two Slater 
 * distributed charges with respect to r.
 * The widths xi and xj may both be zero. The Slater function types i and j
 * should be in the range 1 .. SLATER_MAX
 *
 * \param[in] r  distance in nm
 * \param[in] i  Slater function type
 * \param[in] j  Slater function type
 * \param[in] xi gaussian width in 1/nm 
 * \param[in] xj gaussian width in 1/nm 
 * \return    Integral value
 */    
double DCoulomb_SS(double r,int i,int j,double xi,double xj);

/*! \brief 
 * Compute the derivative of the Slater overlap integral between a Slater 
 * distributed charge and a point charge with respect to r.
 * The width xi may be zero. The Slater function type i should be in the
 * range 1 .. SLATER_MAX
 *
 * \param[in] r  distance in nm
 * \param[in] i  Slater function type
 * \param[in] xi gaussian width in 1/nm 
 * \return    Integral value
 */    
double DNuclear_SS(double r,int i,double xi);

#ifdef __cplusplus
}
#endif

#endif
