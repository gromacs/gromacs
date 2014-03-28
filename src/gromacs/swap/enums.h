/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Enumerated types used in the "Computational Electrophysiology" module.
 *
 * The following enums are mainly used for indexing arrays and when
 * looping over the available ions, channels, or compartments. This hopefully
 * adds to the code's readability because it makes clear which object is dealt
 * with in a block of code.
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 * \inlibraryapi
 * \ingroup module_swap
 */
#ifndef GMX_SWAP_ENUMS_H_
#define GMX_SWAP_ENUMS_H_

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief The two compartments for CompEL setups. */
enum eCompartment {
    eCompA, eCompB, eCompNR
};

/*! \brief The positive and negative ions CompEL setups.
 *
 * Future versions of the protocol might consider more than two types of ions.
 */
enum eIontype {
    eIonNEG, eIonPOS, eIonNR
};

/*! \brief The channels that define with their COM the compartment boundaries in CompEL setups.
 *
 * In principle one could also use modified setups with more than two channels.
 */
enum eChannel {
    eChan0, eChan1, eChanNR
};

#ifdef __cplusplus
}
#endif

#endif
