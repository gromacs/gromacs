/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2010,2014,2015, by the GROMACS development team, led by
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

#ifndef GMX_RANDOM_RANDOM_H
#define GMX_RANDOM_RANDOM_H

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! Fixed random number seeds for different types of RNG */
#define RND_SEED_UPDATE    1 /**< For coordinate update (sd, bd, ..) */
#define RND_SEED_REPLEX    2 /**< For replica exchange */
#define RND_SEED_VRESCALE  3 /**< For V-rescale thermostat */
#define RND_SEED_ANDERSEN  4 /**< For Andersen thermostat */
#define RND_SEED_TPI       5 /**< For test particle insertion */
#define RND_SEED_EXPANDED  6 /**< For expanded emseble methods */

/*! \brief Abstract datatype for a random number generator
 *
 * This is a handle to the full state of a random number generator.
 * You can not access anything inside the gmx_rng structure outside this
 * file.
 */
typedef struct gmx_rng *
    gmx_rng_t;


/*! \brief Returns the size of the RNG integer data structure
 *
 * Returns the size of the RNG integer data structure.
 * \threadsafe Yes.
 */
int
gmx_rng_n(void);


/*! \brief Create a new RNG, seeded from a single integer.
 *
 * If you dont want to pick a seed, just call it as
 * rng=gmx_rng_init(gmx_rng_make_seed()) to seed it from
 * the system time or a random device.
 *
 * \param seed Random seed, unsigned 32-bit integer.
 *
 * \return Reference to a random number generator, or NULL if there was an
 *         error.
 *
 * \threadsafe Yes.
 */
gmx_rng_t
gmx_rng_init(unsigned int seed);


/*! \brief Generate a 'random' RNG seed.
 *
 * This routine tries to get a seed from /dev/random if present,
 * and if not it uses time-of-day and process id to generate one.
 *
 * \return 32-bit unsigned integer random seed.
 *
 * Tip: If you use this in your code, it is a good idea to write the
 * returned random seed to a logfile, so you can recreate the exact sequence
 * of random number if you need to reproduce your run later for one reason
 * or another.
 *
 * \threadsafe Yes.
 */
unsigned int
gmx_rng_make_seed(void);


/*! \brief Initialize a RNG with 624 integers (>32 bits of entropy).
 *
 *  The Mersenne twister RNG used in Gromacs has an extremely long period,
 *  but when you only initialize it with a 32-bit integer there are only
 *  2^32 different possible sequences of number - much less than the generator
 *  is capable of.
 *
 *  If you really need the full entropy, this routine makes it possible to
 *  initialize the RNG with up to 624 32-bit integers, which will give you
 *  up to 2^19968 bits of entropy.
 *
 *  \param seed Array of unsigned integers to form a seed
 *  \param seed_length Number of integers in the array, up to 624 are used.
 *
 * \return Reference to a random number generator, or NULL if there was an
 *         error.
 *
 * \threadsafe Yes.
 */
gmx_rng_t
gmx_rng_init_array(unsigned int    seed[],
                   int             seed_length);


/*! \brief Release resources of a RNG
 *
 *  This routine destroys a random number generator and releases all
 *  resources allocated by it.
 *
 *  \param rng Handle to random number generator previously returned by
 *		       gmx_rng_init() or gmx_rng_init_array().
 *
 * \threadsafe Function itself is threadsafe, but you should only destroy a
 *             certain RNG once (i.e. from one thread).
 */
void
gmx_rng_destroy(gmx_rng_t rng);


/*! \brief Get the state of a RNG
 *
 * This routine stores the random state in \p mt and \p mti.
 *
 * \param[in]  rng Handle to random number generator previously returned by
 *     gmx_rng_init() or gmx_rng_init_array().
 * \param[out] mt  Array of at least 624 integers to receive state.
 * \param[out] mti Pointer to an integer to receive state.
 */
void
gmx_rng_get_state(gmx_rng_t rng, unsigned int *mt, int *mti);


/*! \brief Set the state of a RNG
 *
 * This routine sets the random state from \p mt and \p mti.
 *
 * \param rng Handle to random number generator previously returned by
 *     gmx_rng_init() or gmx_rng_init_array().
 * \param[in]  mt  Array of at least 624 integers.
 * \param[in]  mti Additional integer.
 */
void
gmx_rng_set_state(gmx_rng_t rng, unsigned int *mt, int mti);


/*! \brief Random 32-bit integer from a uniform distribution
 *
 *  This routine returns a random integer from the random number generator
 *  provided, and updates the state of that RNG.
 *
 *  \param rng Handle to random number generator previously returned by
 *		       gmx_rng_init() or gmx_rng_init_array().
 *
 *  \return 32-bit unsigned integer from a uniform distribution.
 *
 *  \threadsafe Function yes, input data no. You should not call this function
 *	        from two different threads using the same RNG handle at the
 *              same time. For performance reasons we cannot lock the handle
 *              with a mutex every time we need a random number - that would
 *              slow the routine down a factor 2-5. There are two simple
 *		solutions: either use a mutex and lock it before calling
 *              the function, or use a separate RNG handle for each thread.
 */
unsigned int
gmx_rng_uniform_uint32(gmx_rng_t rng);


/*! \brief Random gmx_real_t 0<=x<1 from a uniform distribution
 *
 *  This routine returns a random floating-point number from the
 *  random number generator provided, and updates the state of that RNG.
 *
 *  \param rng Handle to random number generator previously returned by
 *		       gmx_rng_init() or gmx_rng_init_array().
 *
 *  \return floating-point number 0<=x<1 from a uniform distribution.
 *
 *  \threadsafe Function yes, input data no. You should not call this function
 *		from two different threads using the same RNG handle at the
 *              same time. For performance reasons we cannot lock the handle
 *              with a mutex every time we need a random number - that would
 *              slow the routine down a factor 2-5. There are two simple
 *		solutions: either use a mutex and lock it before calling
 *              the function, or use a separate RNG handle for each thread.
 */
real
gmx_rng_uniform_real(gmx_rng_t rng);


/*! \brief Random gmx_real_t from a gaussian distribution
 *
 *  This routine returns a random floating-point number from the
 *  random number generator provided, and updates the state of that RNG.
 *
 *  The Box-Muller algorithm is used to provide gaussian random numbers. This
 *  is not the fastest known algorithm for gaussian numbers, but in contrast
 *  to the alternatives it is very well studied and you can trust the returned
 *  random numbers to have good properties and no correlations.
 *
 *  \param rng Handle to random number generator previously returned by
 *			  gmx_rng_init() or gmx_rng_init_array().
 *
 *  \return Gaussian random floating-point number with average 0.0 and
 *	    standard deviation 1.0. You can get any average/mean you want
 *          by first multiplying with the desired average and then adding
 *          the average you want.
 *
 *  \threadsafe Function yes, input data no. You should not call this function
 *		from two different threads using the same RNG handle at the
 *              same time. For performance reasons we cannot lock the handle
 *              with a mutex every time we need a random number - that would
 *              slow the routine down a factor 2-5. There are two simple
 *		solutions: either use a mutex and lock it before calling
 *              the function, or use a separate RNG handle for each thread.
 *
 *  It works perfectly to mix calls to get uniform and gaussian random numbers
 *  from the same generator, but since it will affect the sequence of returned
 *  numbers it is probably better to use separate random number generator
 *  structures.
 */
real
gmx_rng_gaussian_real(gmx_rng_t rng);



/* Return a new gaussian random number with expectation value
 * 0.0 and standard deviation 1.0. This routine uses a table
 * lookup for maximum speed.
 *
 * WARNING: The lookup table is 16k by default, which means
 *          the granularity of the random numbers is coarser
 *	    than what you get from gmx_rng_gauss_real().
 *          In most cases this is no problem whatsoever,
 *          and it is particularly true for BD/SD integration.
 *	    Note that you will NEVER get any really extreme
 *          numbers: the maximum absolute value returned is
 *          4.0255485.
 *
 * threadsafe: yes
 */
real
gmx_rng_gaussian_table(gmx_rng_t rng);


/* The stateless cycle based random number generators below,
 * which all use threefry2x64, take the following arguments:
 *
 * ctr1: In mdrun the step counter, in tools the frame(-step)
 *       counter, so we can ensure reproducible results, even
 *       we starting at different steps/frames. Might need to be
 *       multiplied by a constant if we need more random numbers.
 * ctr2: A local counter, in mdrun often a global atom index.
 *       If any algorithm needs a variable number of random numbers,
 *       the second counter is usually a function of the local
 *       counter.
 * key1: A user provided random seed.
 * key2: A fixed seed which is particular for the algorithm,
 *       as defined at the top of this file, to ensure different
 *       random sequences when the same user seed is used for
 *       different algorithms.
 */

/* Return two uniform random numbers with 2^53 equally
 * probable values between 0 and 1 - 2^-53.
 * It uses a stateless counter based random number generator
 * (threefry2x64).
 */
void
gmx_rng_cycle_2uniform(gmx_uint64_t ctr1, gmx_uint64_t ctr2,
                       gmx_uint64_t key1, gmx_uint64_t key2,
                       double* rnd);

/* Return three Gaussian random numbers with expectation value
 * 0.0 and standard deviation 1.0. This routine uses a table
 * lookup for maximum speed. It uses a stateless counter
 * based random number generator (threefry2x64). See
 * gmx_rng_gaussian_table for a warning about accuracy of the table.
 *
 * threadsafe: yes
 */
void
gmx_rng_cycle_3gaussian_table(gmx_uint64_t ctr1, gmx_uint64_t ctr2,
                              gmx_uint64_t key1, gmx_uint64_t key2,
                              real* rnd);

/* As gmx_rng_3gaussian_table, but returns 6 Gaussian numbers. */
void
gmx_rng_cycle_6gaussian_table(gmx_uint64_t ctr1, gmx_uint64_t ctr2,
                              gmx_uint64_t key1, gmx_uint64_t key2,
                              real* rnd);

#ifdef __cplusplus
}
#endif

#endif
