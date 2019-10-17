/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
#ifndef GMX_FILEIO_TPXIO_H
#define GMX_FILEIO_TPXIO_H

#include <cstdio>

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct t_atoms;
struct t_block;
struct t_inputrec;
class t_state;
struct t_topology;

/*! \libinternal
 * \brief
 * First part of the TPR file structure containing information about
 * the general aspect of the system.
 */
struct TpxFileHeader
{
    //! Non zero if input_rec is present.
    bool bIr = false;
    //! Non zero if a box is present.
    bool bBox = false;
    //! Non zero if a topology is present.
    bool bTop = false;
    //! Non zero if coordinates are present.
    bool bX = false;
    //! Non zero if velocities are present.
    bool bV = false;
    //! Non zero if forces are present (no longer supported, but retained so old .tpr can be read)
    bool bF = false;
    //! The total number of atoms.
    int natoms = 0;
    //! The number of temperature coupling groups.
    int ngtc = 0;
    //! Current value of lambda.
    real lambda = 0;
    //! Current value of the alchemical state - not yet printed out.
    int fep_state = 0;
    /*a better decision will eventually (5.0 or later) need to be made
       on how to treat the alchemical state of the system, which can now
       vary through a simulation, and cannot be completely described
       though a single lambda variable, or even a single state
       index. Eventually, should probably be a vector. MRS*/
    //! Size of the TPR body in chars (equal to number of bytes) during I/O.
    int64_t sizeOfTprBody = 0;
    //! File version.
    int fileVersion = 0;
    //! File generation.
    int fileGeneration = 0;
    //! If the tpr file was written in double precision.
    bool isDouble = false;
};

/*! \brief
 * Contains the partly deserialized contents of a TPR file.
 *
 * Convenience struct that holds a fully deserialized TPR file header,
 * and the body of the TPR file as char buffer that can be deserialized
 * independently from the header.
 */
struct PartialDeserializedTprFile
{
    //! The file header.
    TpxFileHeader header;
    //! The file body.
    std::vector<char> body;
    //! Flag for PBC needed by legacy implementation.
    int ePBC = -1;
};

/*
 * These routines handle reading and writing of preprocessed
 * topology files in any of the following formats:
 * TPR : topology in XDR format, portable accross platforms
 *
 * Files are written in the precision with which the source are compiled,
 * but double and single precision can be read by either.
 */

/*! \brief
 * Read the header from a tpx file and then close it again.
 *
 * By setting \p canReadTopologyOnly to true, it is possible to read future
 * versions too (we skip the changed inputrec), provided we havent
 * changed the topology description. If it is possible to read
 * the inputrec it will still be done even if canReadTopologyOnly is true.
 *
 * \param[in] fileName The name of the input file.
 * \param[in] canReadTopologyOnly If reading the inputrec can be skipped or not.
 * \returns An initialized and populated TPX File header object.
 */
TpxFileHeader readTpxHeader(const char* fileName, bool canReadTopologyOnly);

void write_tpx_state(const char* fn, const t_inputrec* ir, const t_state* state, const gmx_mtop_t* mtop);
/* Write a file, and close it again.
 */

/*! \brief
 * Complete deserialization of TPR file into the individual data structures.
 *
 * If \p state is nullptr, only populates ir and mtop.
 *
 * \param[in] partialDeserializedTpr Struct with header and char buffer needed to populate system.
 * \param[out] ir Input rec to populate.
 * \param[out] state System state variables to populate.
 * \param[out] x Separate vector for coordinates, deprecated.
 * \param[out] v Separate vector for velocities, deprecated.
 * \param[out] mtop Global topology to populate.
 *
 * \returns PBC flag.
 */
int completeTprDeserialization(PartialDeserializedTprFile* partialDeserializedTpr,
                               t_inputrec*                 ir,
                               t_state*                    state,
                               rvec*                       x,
                               rvec*                       v,
                               gmx_mtop_t*                 mtop);

//! Overload for final TPR deserialization when not using state vectors.
int completeTprDeserialization(PartialDeserializedTprFile* partialDeserializedTpr,
                               t_inputrec*                 ir,
                               gmx_mtop_t*                 mtop);

/*! \brief
 * Read a file to set up a simulation and close it after reading.
 *
 * Main function used to initialize simulations. Reads the input \p fn
 * to populate the \p state, \p ir and \p mtop needed to run a simulations.
 *
 * This function returns the partial deserialized TPR file
 * that can then be communicated to set up non-master nodes to run simulations.
 *
 * \param[in] fn Input file name.
 * \param[out] ir Input parameters to be set, or nullptr.
 * \param[out] state State variables for the simulation.
 * \param[out] mtop Global simulation topolgy.
 * \returns Struct with header and body in char vector.
 */
PartialDeserializedTprFile read_tpx_state(const char* fn, t_inputrec* ir, t_state* state, gmx_mtop_t* mtop);

/*! \brief
 * Read a file and close it again.
 *
 * Reads a topology input file and populates the fields if the passed
 * variables are valid. It is possible to pass \p ir, \p natoms,
 * \p x, \p v or \p mtop as nullptr to the function. In those cases,
 * the variables will not be populated from the input file. Passing both
 * \p x and \p v as nullptr is not supported. If both \p natoms and
 * \p mtop are passed as valid objects to the function, the total atom
 * number from \p mtop will be set in \p natoms. Otherwise \p natoms
 * will not be changed. If \p box is valid, the box will be set from
 * the information read in from the file.
 *
 * \param[in] fn Input file name.
 * \param[out] ir Input parameters to be set, or nullptr.
 * \param[out] box Box matrix.
 * \param[out] natoms Total atom numbers to be set, or nullptr.
 * \param[out] x Positions to be filled from file, or nullptr.
 * \param[out] v Velocities to be filled from file, or nullptr.
 * \param[out] mtop Topology to be populated, or nullptr.
 * \returns ir->ePBC if it was read from the file.
 */
int read_tpx(const char* fn, t_inputrec* ir, matrix box, int* natoms, rvec* x, rvec* v, gmx_mtop_t* mtop);

int read_tpx_top(const char* fn, t_inputrec* ir, matrix box, int* natoms, rvec* x, rvec* v, t_topology* top);
/* As read_tpx, but for the old t_topology struct */

gmx_bool fn2bTPX(const char* file);
/* return if *file is one of the TPX file types */

void pr_tpxheader(FILE* fp, int indent, const char* title, const TpxFileHeader* sh);

#endif
