/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2018,2019,2020, by the GROMACS development team, led by
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

#include "gmxpre.h"

#include "seed.h"

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"

namespace gmx
{


/*! \brief Check if the RDRAND random device functioning correctly
 *
 * Due to a bug in AMD Ryzen microcode, RDRAND may always return -1 (0xFFFFFFFF).
 * To avoid that, fall back to using PRNG instead of RDRAND if this happens.
 *
 * \param[in] rd Random device to check.
 *
 * \returns The result of the checks.
 */
static bool checkIfRandomDeviceIsFunctional(std::random_device* rd)
{
    uint64_t randomNumber1 = static_cast<uint64_t>((*rd)());
    uint64_t randomNumber2 = static_cast<uint64_t>((*rd)());

    // Due to a bug in AMD Ryzen microcode, RDRAND may always return -1 (0xFFFFFFFF).
    // To avoid that, fall back to using PRNG instead of RDRAND if this happens.
    if (randomNumber1 == 0xFFFFFFFF && randomNumber2 == 0xFFFFFFFF)
    {
        if (debug)
        {
            fprintf(debug,
                    "Hardware random number generator (RDRAND) returned -1 (0xFFFFFFFF) twice in\n"
                    "a row. This may be due to a known bug in AMD Ryzen microcode.");
        }
        return false;
    }
    return true;
}

/*! \brief Get the next pure or pseudo-random number
 *
 * Returns the next random number taken from the hardware generator or from standard PRNG.
 *
 * \param[in] useRdrand  Whether the hardware source of random numbers should be used.
 * \param[in] rd         Pointer to the random device to use.
 *
 * \return Random or pseudo-random number.
 */
static inline uint64_t getNextRandomNumber(bool useRdrand, std::random_device* rd)
{
    if (debug && !useRdrand)
    {
        fprintf(debug,
                "Will use pseudo-random number generator (PRNG) rather than hardware device.");
    }
    return useRdrand ? static_cast<uint64_t>((*rd)()) : static_cast<uint64_t>(rand());
}

uint64_t makeRandomSeed()
{
    std::random_device rd;

    bool useRdrand = checkIfRandomDeviceIsFunctional(&rd);

    const std::size_t resultBits = std::numeric_limits<uint64_t>::digits;
    const std::size_t numBitsInRandomNumber =
            useRdrand ? std::numeric_limits<std::random_device::result_type>::digits
                      : std::numeric_limits<int>::digits;

    uint64_t result = getNextRandomNumber(useRdrand, &rd);
    for (std::size_t bits = numBitsInRandomNumber; bits < resultBits; bits += numBitsInRandomNumber)
    {
        result = (result << numBitsInRandomNumber) | getNextRandomNumber(useRdrand, &rd);
    }
    return result;
}

} // namespace gmx
