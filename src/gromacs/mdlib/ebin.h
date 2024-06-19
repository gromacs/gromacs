/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
/* \internal \file
 * \brief
 * Declares structures and interfaces to store, compute and print current and average values for thermodynamics properties.
 *
 * The word 'energy' is used here in wide scope and refer to any thermodynamic quantity that can benefit from
 * averaging (e.g. temperature, pressure).
 *
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_EBIN_H
#define GMX_MDLIB_EBIN_H

#include <cstdint>
#include <cstdio>

#include "gromacs/fileio/enxio.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

namespace gmx
{
template<typename>
class ArrayRef;
}

/* \brief Running averaging structure ('energy bin') to store thermodynamic values.
 *
 * Collects the data on thermodynamic parameters (energy terms, temperature,
 * pressure etc.) during the run, including their current and average values.
 *
 * \todo Clean this structure from unused values.
 */
struct t_ebin
{
    //! Number of thermodynamic terms
    int nener;
    //! Name and units for each term
    gmx_enxnm_t* enm;
    //! Number of steps used for sum (for energy history)
    int64_t nsteps;
    //! Number of values added to the sum so far
    int64_t nsum;
    //! Term values: each structure stores current, running average and sum.
    t_energy* e;
    //! Total number of steps saved (for energy history)
    int64_t nsteps_sim;
    //! Total number of values added to sum (used when printing average values at the end of the run)
    int64_t nsum_sim;
    //! Energy values throughout the entire simulation: structure stores current, average and sum, but only sum value is used to compute averages
    t_energy* e_sim;
};

/* \brief Type of expected output: normal or average.
 */
enum
{
    eprNORMAL,
    eprAVER,
    eprNR
};

/*! \brief Create the structure to store thermodynamic properties*/
t_ebin* mk_ebin();

/*! \brief Destroy the \c eb structure.
 *
 * \param[in,out] eb  Pointer to the structure to destroy.
 */
void done_ebin(t_ebin* eb);

/*! \brief Create space for the extra thermodynamic term(s) and register its(their) name(s).
 *
 * The enm array must be static, because the contents are not copied, only the pointers.
 *
 * \param[in] eb     Srtucture in which the space for the termodynamic terms shall be created..
 * \param[in] nener  Number of thermodyamic terms to allocate memory for.
 * \param[in] enm    Names of the terms.
 * \param[in] unit   Units.
 *
 * \returns          A serial number (index) for the newly allocated terms.
 */
int get_ebin_space(t_ebin* eb, int nener, const char* const enm[], const char* unit);

/*! \brief Add current value of the thermodynamic term(s) to the bin(s).
 *
 * Add nener reals (eg. energies, box-lengths, pressures) to the at position specified by \c
 * entryIndex. If bSum is TRUE then the reals are also added to the sum and sum of squares.
 *
 * \param[in] eb          Structure that stores the thermodynamic values.
 * \param[in] entryIndex  Internal index of the term(s) to add.
 * \param[in] nener       Number of the terms to add.
 * \param[in] ener        Value(s) of thermodynamic term(s) (nener ptc.)
 * \param[in] bSum        If the average value should be accumulated for this term(s).
 */
void add_ebin(t_ebin* eb, int entryIndex, int nener, const real ener[], gmx_bool bSum);


/*! \brief Add values from array to the bins if the matching entry in \c shouldUse is true.
 *
 * Caller must ensure that \c shouldUse and \c ener to have the same
 * size, and that \c eb has enough room for the number of true
 * entries in \c shouldUse.
 *
 * \param[in] eb          Structure that stores the thermodynamic values.
 * \param[in] entryIndex  Internal index of the term(s).
 * \param[in] shouldUse   Array of booleans that indicate which terms should be used.
 * \param[in] ener        Values of thermodinamic terms to add.
 * \param[in] bSum        If the average value should be accumulated for these terms.
 */
void add_ebin_indexed(t_ebin*                   eb,
                      int                       entryIndex,
                      gmx::ArrayRef<bool>       shouldUse,
                      gmx::ArrayRef<const real> ener,
                      gmx_bool                  bSum);

/*! \brief Increase the counters for the sums.
 *
 * This routine should be called after all add_ebin calls for this step.
 *
 * \param[in] increment   How much counts should be increased
 * \param[in] eb          Structure that stores the thermodynamic values.
 * \param[in] bSum        If the sums counters should be increased as well.
 */
void ebin_increase_count(int increment, t_ebin* eb, gmx_bool bSum);


/*! \brief Reset the average and fluctuation sums.
 *
 * \param[in] eb          Structure that stores the thermodynamic values.
 */
void reset_ebin_sums(t_ebin* eb);


/*! \brief Print the contents of some energy bins.
 *
 * We will print \c nperline entries on a text line (advisory <=
 * 5). \c prmode may be any of the above listed enum values. \c tsteps is
 * used only when \c eprAVER is set. If \c bPrHead than the
 * header is printed.
 *
 * \c entryIndex and \c nener must be in [0,\c eb->nener), except that \c
 * nener -1 is interpreted as \c eb->nener.
 *
 * \todo Callers should be refactored to pass \c eb->nener, rather than
 *       us implement and rely on this special behavior of -1.
 *
 * \param[in] fp          I/O pointer to print to.
 * \param[in] eb          Structure that stores the thermodynamic values.
 * \param[in] entryIndex  Internal index of the term(s).
 * \param[in] nener       Number of the terms to print.
 * \param[in] nperline    Number of values per line.
 * \param[in] prmode      Print current (eprNORMAL) or average (eprAVER) values.
 * \param[in] bPrHead     If the header should be printed.
 */
void pr_ebin(FILE* fp, t_ebin* eb, int entryIndex, int nener, int nperline, int prmode, gmx_bool bPrHead);

#endif
