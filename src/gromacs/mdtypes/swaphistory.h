/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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

/*
 * This file contains data types containing ioin/water position exchange
 * data to be stored in the checkpoint file.
 */

#ifndef GMX_MDLIB_SWAPHISTORY_H
#define GMX_MDLIB_SWAPHISTORY_H

#include "gromacs/mdtypes/md_enums.h"

/* History of an ion type used in position swapping
 */
typedef struct swapstateIons_t
{
    int  nMolReq[eCompNR];                    // Requested # of molecules per compartment
    int *nMolReq_p[eCompNR];                  // Pointer to this data (for checkpoint writing)
    int  inflow_net[eCompNR];                 // Flux determined from the # of swaps
    int *inflow_net_p[eCompNR];               // Pointer to this data
    int *nMolPast[eCompNR];                   // Array with nAverage entries for history
    int *nMolPast_p[eCompNR];                 // Pointer points to the first entry only

    // Channel flux detection, this is counting only and has no influence on whether swaps are performed or not:                                                                 */
    int            fluxfromAtoB[eCompNR];     // Flux determined from the split cylinders
    int           *fluxfromAtoB_p[eCompNR];   // Pointer to this data
    int            nMol;                      // Number of molecules, size of the following arrays
    unsigned char *comp_from;                 // Ion came from which compartment?
    unsigned char *channel_label;             // Through which channel did this ion pass?
} swapstateIons_t;

/* Position swapping state
 *
 * To also make multimeric channel proteins whole, we save the last whole configuration
 * of the channels in the checkpoint file. If we have no checkpoint file, we assume
 * that the starting configuration has the correct PBC representation after making the
 * individual molecules whole
 *
 * \todo move out of this file to ObservablesHistory
 *
 */
typedef struct swaphistory_t
{
    int              eSwapCoords;             // Swapping along x, y, or z-direction?
    int              nIonTypes;               // Number of ion types, this is the size of the following arrays
    int              nAverage;                // Use average over this many swap attempt steps when determining the ion counts
    int              fluxleak;                // Ions not going through any channel (bad!)
    int             *fluxleak_p;              // Pointer to this data
    gmx_bool         bFromCpt;                // Did we start from a checkpoint file?
    int              nat[eChanNR];            // Size of xc_old_whole, i.e. the number of atoms in each channel
    rvec            *xc_old_whole[eChanNR];   // Last known whole positions of the two channels (important for multimeric ch.!)
    rvec           **xc_old_whole_p[eChanNR]; // Pointer to these positions
    swapstateIons_t *ionType;                 // History information for one ion type
}
swaphistory_t;

#endif
