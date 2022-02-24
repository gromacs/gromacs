/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/*
 * This file contains data types containing ioin/water position exchange
 * data to be stored in the checkpoint file.
 */

#ifndef GMX_MDLIB_SWAPHISTORY_H
#define GMX_MDLIB_SWAPHISTORY_H

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/enumerationhelpers.h"

enum class Domain : int;
enum class ChannelHistory : int;
/* History of an ion type used in position swapping
 */
struct swapstateIons_t
{
    gmx::EnumerationArray<Compartment, int> nMolReq; // Requested # of molecules per compartment
    gmx::EnumerationArray<Compartment, int*> nMolReq_p; // Pointer to this data (for checkpoint writing)
    gmx::EnumerationArray<Compartment, int>  inflow_net;   // Flux determined from the # of swaps
    gmx::EnumerationArray<Compartment, int*> inflow_net_p; // Pointer to this data
    gmx::EnumerationArray<Compartment, int*> nMolPast;   // Array with nAverage entries for history
    gmx::EnumerationArray<Compartment, int*> nMolPast_p; // Pointer points to the first entry only

    // Channel flux detection, this is counting only and has no influence on whether swaps are performed or not:                                                                 */
    gmx::EnumerationArray<Channel, int>  fluxfromAtoB;   // Flux determined from the split cylinders
    gmx::EnumerationArray<Channel, int*> fluxfromAtoB_p; // Pointer to this data
    int                                  nMol; // Number of molecules, size of the following arrays
    Domain*                              comp_from;     // Ion came from which compartment?
    ChannelHistory*                      channel_label; // Through which channel did this ion pass?
};

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
    SwapType eSwapCoords; // Swapping along x, y, or z-direction?
    int      nIonTypes;   // Number of ion types, this is the size of the following arrays
    int  nAverage; // Use average over this many swap attempt steps when determining the ion counts
    int  fluxleak; // Ions not going through any channel (bad!)
    int* fluxleak_p;                         // Pointer to this data
    bool bFromCpt;                           // Did we start from a checkpoint file?
    gmx::EnumerationArray<Channel, int> nat; // Size of xc_old_whole, i.e. the number of atoms in each channel
    gmx::EnumerationArray<Channel, rvec*> xc_old_whole; // Last known whole positions of the two channels (important for multimeric ch.!)
    gmx::EnumerationArray<Channel, rvec**> xc_old_whole_p; // Pointer to these positions
    swapstateIons_t*                       ionType;        // History information for one ion type
} swaphistory_t;

#endif
