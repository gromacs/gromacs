/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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
/*! \defgroup module_random Random engines and distributions (random)
 * \ingroup group_utilitymodules
 * \brief
 * Provides efficient and portable random generators and distributions
 *
 * <H3>Basic Use</H3>
 *
 * \Gromacs relies on random numbers in several different modules, and in
 * particular for methods that influence the integration we both require the
 * generation to be very fast and the resulting numbers of high quality.
 * In addition, it is highly desirable that we generate the same trajectories
 * in parallel as for a single-core run.
 *
 * To realize this, we have implemented the ThreeFry2x64 counter-based random
 * engine. In contrast to a normal random engine that is seeded and then keeps
 * an internal state, ThreeFry2x64 is derived from cryptographic applications
 * where we use a key to turn a highly regular counter int a stream of random
 * numbers. This makes it possible to quickly set the counter in the random
 * engine based e.g. on the timestep and atom index, and get the same random
 * numbers regardless of parallelization.
 *
 * The TreeFry2x64 engine has been implemented to be fully compatible with
 * standard C++11 random engines. There is a gmx::ThreeFry2x64General class that
 * allows full control over the accuracy (more iterations means higher quality),
 * and gmx::ThreeFry2x64 and gmx::ThreeFry2x64Fast that are specialized to 20
 * and 13 iterations, respectively. With 20 iterations this engine passes all
 * tests in the standard BigCrush test, and with 13 iterations only a single
 * test fails (in comparision, Mersenne Twister fails two).
 *
 * All these engines take a template parameter that specifies the number of
 * bits to reserve for an internal counter. This is based on an idea of
 * John Salmon, and it makes it possible to set your external counter based
 * on two simple values (usually timestep and particle index), but then it is
 * still possible to draw more than one value for this external counter since
 * the internal counter increments. If you run out of internal counter space
 * the class will throw an exception to make sure you don't silently end up
 * with corrupted/overlapping random data.
 *
 * <H3>But what if I just want a vanilla random number generator?</H3>
 *
 * We've thought about that. Just use the gmx::DefaultRandomEngine class and
 * forget everything about counters. Initialize the class with a single value
 * for the seed (up to 64 bits), and you are good to go.
 *
 * <H3>Random number distributions</H3>
 *
 * The ThreeFry random engine is fully compatible with all distributions from
 * the C++11 standard library, but unfortunately implementation differences
 * (and bugs) mean you will typically not get the same sequence of numbers from
 * two different library implementations. Since this causes problems for our
 * unit tests, we prefer to use our own implementations - they should work
 * exactly like the corresponding C++11 versions.
 *
 * The normal distribution is frequently used in integration, and it can be
 * a performance bottleneck. To avoid this, we use a special tabulated
 * distribution gmx::TabulatedNormalDistribution that provides very high
 * performance at the cost of slightly discretized values; the default 14-bit
 * table gives us 16,384 unique values, but this has been thoroughly tested to
 * be sufficient for all integration usage.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 */
/*! \file
 * \brief
 * Public API convenience header for random engines and distributions.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \inpublicapi
 * \ingroup module_random
 */
#ifndef GMX_RANDOM_H
#define GMX_RANDOM_H

#include "gromacs/random/exponentialdistribution.h"
#include "gromacs/random/gammadistribution.h"
#include "gromacs/random/normaldistribution.h"
#include "gromacs/random/seed.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformintdistribution.h"
#include "gromacs/random/uniformrealdistribution.h"

#endif
