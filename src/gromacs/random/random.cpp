/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2012,2014,2015, by the GROMACS development team, led by
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

#include "random.h"

#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "external/Random123-1.08/include/Random123/threefry.h"

#include "gromacs/math/utilities.h"
#include "gromacs/random/random_gausstable.h"
#include "gromacs/utility/sysinfo.h"

#define RNG_N 624
#define RNG_M 397
#define RNG_MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define RNG_UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define RNG_LOWER_MASK 0x7fffffffUL /* least significant r bits */

/* Note that if you change the size of the Gaussian table you will
 * also have to generate new initialization data for the table in
 * gmx_random_gausstable.h
 *
 * A routine print_gaussian_table() is in contrib/random.c
 * for convenience - use it if you need a different size of the table.
 */
#define GAUSS_TABLE 14 /* the size of the gauss table is 2^GAUSS_TABLE */
#define GAUSS_MASK  ((1 << GAUSS_TABLE) - 1)


struct gmx_rng {
    unsigned int  mt[RNG_N];
    int           mti;
    int           has_saved;
    double        gauss_saved;
};



int
gmx_rng_n(void)
{
    return RNG_N;
}


gmx_rng_t
gmx_rng_init(unsigned int seed)
{
    struct gmx_rng *rng;

    if ((rng = (struct gmx_rng *)malloc(sizeof(struct gmx_rng))) == NULL)
    {
        return NULL;
    }

    rng->has_saved = 0; /* no saved gaussian number yet */

    rng->mt[0] = seed & 0xffffffffUL;
    for (rng->mti = 1; rng->mti < RNG_N; rng->mti++)
    {
        rng->mt[rng->mti] =
            (1812433253UL * (rng->mt[rng->mti-1] ^
                             (rng->mt[rng->mti-1] >> 30)) + rng->mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        rng->mt[rng->mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
    return rng;
}

gmx_rng_t
gmx_rng_init_array(unsigned int seed[], int seed_length)
{
    int       i, j, k;
    gmx_rng_t rng;

    if ((rng = gmx_rng_init(19650218UL)) == NULL)
    {
        return NULL;
    }

    i = 1; j = 0;
    k = (RNG_N > seed_length ? RNG_N : seed_length);
    for (; k; k--)
    {
        rng->mt[i] = (rng->mt[i] ^ ((rng->mt[i-1] ^
                                     (rng->mt[i-1] >> 30)) * 1664525UL))
            + seed[j] + j;          /* non linear */
        rng->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i >= RNG_N)
        {
            rng->mt[0] = rng->mt[RNG_N-1]; i = 1;
        }
        if (j >= seed_length)
        {
            j = 0;
        }
    }
    for (k = RNG_N-1; k; k--)
    {
        rng->mt[i] = (rng->mt[i] ^ ((rng->mt[i-1] ^
                                     (rng->mt[i-1] >> 30)) *
                                    1566083941UL))
            - i;                    /* non linear */
        rng->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i >= RNG_N)
        {
            rng->mt[0] = rng->mt[RNG_N-1]; i = 1;
        }
    }

    rng->mt[0] = 0x80000000UL;
    /* MSB is 1; assuring non-zero initial array */
    return rng;
}


void
gmx_rng_destroy(gmx_rng_t rng)
{
    if (rng)
    {
        free(rng);
    }
}


void
gmx_rng_get_state(gmx_rng_t rng, unsigned int *mt, int *mti)
{
    int i;

    for (i = 0; i < RNG_N; i++)
    {
        mt[i] = rng->mt[i];
    }
    *mti = rng->mti;
}


void
gmx_rng_set_state(gmx_rng_t rng,  unsigned int *mt, int mti)
{
    int i;

    for (i = 0; i < RNG_N; i++)
    {
        rng->mt[i] = mt[i];
    }
    rng->mti = mti;
}


unsigned int
gmx_rng_make_seed(void)
{
    FILE              *fp;
    unsigned int       data;
    int                ret;
    long               my_pid;

#ifdef HAVE_UNISTD_H
    /* We never want Gromacs execution to halt 10-20 seconds while
     * waiting for enough entropy in the random number generator.
     * For this reason we should NOT use /dev/random, which will
     * block in cases like that. That will cause all sorts of
     * Gromacs programs to block ~20 seconds while waiting for a
     * super-random-number to generate cool quotes. Apart from the
     * minor irritation, it is really bad behavior of a program
     * to abuse the system random numbers like that - other programs
     * need them too.
     * For this reason, we ONLY try to get random numbers from
     * the pseudo-random stream /dev/urandom, and use other means
     * if it is not present (in this case fopen() returns NULL).
     */
    fp = fopen("/dev/urandom", "rb");
#else
    fp = NULL;
#endif
    if (fp != NULL)
    {
        ret = fread(&data, sizeof(unsigned int), 1, fp);
        GMX_IGNORE_RETURN_VALUE(ret);
        fclose(fp);
    }
    else
    {
        /* No random device available, use time-of-day and process id */
        my_pid = gmx_getpid();
        data   = (unsigned int)(((long)time(NULL)+my_pid) % (long)1000000);
    }
    return data;
}


/* The random number state contains RNG_N entries that are returned one by
 * one as random numbers. When we run out of them, this routine is called to
 * regenerate RNG_N new entries.
 */
static void
gmx_rng_update(gmx_rng_t rng)
{
    unsigned int       x1, x2, y, *mt;
    int                k;
    const unsigned int mag01[2] = {0x0UL, RNG_MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    /* update random numbers */
    mt  = rng->mt;  /* pointer to array - avoid repeated dereferencing */

    x1        = mt[0];
    for (k = 0; k < RNG_N-RNG_M-3; k += 4)
    {
        x2      = mt[k+1];
        y       = (x1 & RNG_UPPER_MASK) | (x2 & RNG_LOWER_MASK);
        mt[k]   = mt[k+RNG_M]   ^ (y >> 1) ^ mag01[y & 0x1UL];
        x1      = mt[k+2];
        y       = (x2 & RNG_UPPER_MASK) | (x1 & RNG_LOWER_MASK);
        mt[k+1] = mt[k+RNG_M+1] ^ (y >> 1) ^ mag01[y & 0x1UL];
        x2      = mt[k+3];
        y       = (x1 & RNG_UPPER_MASK) | (x2 & RNG_LOWER_MASK);
        mt[k+2] = mt[k+RNG_M+2] ^ (y >> 1) ^ mag01[y & 0x1UL];
        x1      = mt[k+4];
        y       = (x2 & RNG_UPPER_MASK) | (x1 & RNG_LOWER_MASK);
        mt[k+3] = mt[k+RNG_M+3] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    x2        = mt[k+1];
    y         = (x1 & RNG_UPPER_MASK) | (x2 & RNG_LOWER_MASK);
    mt[k]     = mt[k+RNG_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    k++;
    x1        = mt[k+1];
    y         = (x2 & RNG_UPPER_MASK) | (x1 & RNG_LOWER_MASK);
    mt[k]     = mt[k+RNG_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    k++;
    x2        = mt[k+1];
    y         = (x1 & RNG_UPPER_MASK) | (x2 & RNG_LOWER_MASK);
    mt[k]     = mt[k+RNG_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    k++;
    for (; k < RNG_N-1; k += 4)
    {
        x1      = mt[k+1];
        y       = (x2 & RNG_UPPER_MASK) | (x1 & RNG_LOWER_MASK);
        mt[k]   = mt[k+(RNG_M-RNG_N)]   ^ (y >> 1) ^ mag01[y & 0x1UL];
        x2      = mt[k+2];
        y       = (x1 & RNG_UPPER_MASK) | (x2 & RNG_LOWER_MASK);
        mt[k+1] = mt[k+(RNG_M-RNG_N)+1] ^ (y >> 1) ^ mag01[y & 0x1UL];
        x1      = mt[k+3];
        y       = (x2 & RNG_UPPER_MASK) | (x1 & RNG_LOWER_MASK);
        mt[k+2] = mt[k+(RNG_M-RNG_N)+2] ^ (y >> 1) ^ mag01[y & 0x1UL];
        x2      = mt[k+4];
        y       = (x1 & RNG_UPPER_MASK) | (x2 & RNG_LOWER_MASK);
        mt[k+3] = mt[k+(RNG_M-RNG_N)+3] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    y           = (x2 & RNG_UPPER_MASK) | (mt[0] & RNG_LOWER_MASK);
    mt[RNG_N-1] = mt[RNG_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

    rng->mti = 0;
}


real
gmx_rng_gaussian_real(gmx_rng_t rng)
{
    real x, y, r;

    if (rng->has_saved)
    {
        rng->has_saved = 0;
        return rng->gauss_saved;
    }
    else
    {
        do
        {
            x = 2.0*gmx_rng_uniform_real(rng)-1.0;
            y = 2.0*gmx_rng_uniform_real(rng)-1.0;
            r = x*x+y*y;
        }
        while (r > 1.0 || r == 0.0);

        r                = sqrt(-2.0*log(r)/r);
        rng->gauss_saved = y*r; /* save second random number */
        rng->has_saved   = 1;
        return x*r;             /* return first random number */
    }
}




/* Return a random unsigned integer, i.e. 0..4294967295
 * Provided in header file for performace reasons.
 * Unfortunately this function cannot be inlined, since
 * it needs to refer the internal-linkage gmx_rng_update().
 */
unsigned int
gmx_rng_uniform_uint32(gmx_rng_t rng)
{
    unsigned int y;

    if (rng->mti == RNG_N)
    {
        gmx_rng_update(rng);
    }
    y = rng->mt[rng->mti++];

    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}





/* Return a uniform floating point number on the interval 0<=x<1 */
real
gmx_rng_uniform_real(gmx_rng_t rng)
{
    if (sizeof(real) == sizeof(double))
    {
        return ((double)gmx_rng_uniform_uint32(rng))*(1.0/4294967296.0);
    }
    else
    {
        return ((float)gmx_rng_uniform_uint32(rng))*(1.0/4294967423.0);
    }
    /* divided by the smallest number that will generate a
     * single precision real number on 0<=x<1.
     * This needs to be slightly larger than MAX_UNIT since
     * we are limited to an accuracy of 1e-7.
     */
}
real

gmx_rng_gaussian_table(gmx_rng_t rng)
{
    unsigned int i;

    i = gmx_rng_uniform_uint32(rng);

    /* The Gaussian table is a static constant in this file */
    return gaussian_table[i >> (32 - GAUSS_TABLE)];
}

void
gmx_rng_cycle_2uniform(gmx_uint64_t ctr1, gmx_uint64_t ctr2,
                       gmx_uint64_t key1, gmx_uint64_t key2,
                       double* rnd)
{
    const gmx_int64_t  mask_53bits     = 0x1FFFFFFFFFFFFFULL;
    const double       two_power_min53 = 1.0/9007199254740992.0;

    threefry2x64_ctr_t ctr  = {{ctr1, ctr2}};
    threefry2x64_key_t key  = {{key1, key2}};
    threefry2x64_ctr_t rand = threefry2x64(ctr, key);

    rnd[0] = (rand.v[0] & mask_53bits)*two_power_min53;
    rnd[1] = (rand.v[1] & mask_53bits)*two_power_min53;
}

void
gmx_rng_cycle_3gaussian_table(gmx_uint64_t ctr1, gmx_uint64_t ctr2,
                              gmx_uint64_t key1, gmx_uint64_t key2,
                              real* rnd)
{
    threefry2x64_ctr_t ctr  = {{ctr1, ctr2}};
    threefry2x64_key_t key  = {{key1, key2}};
    threefry2x64_ctr_t rand = threefry2x64(ctr, key);

    rnd[0] = gaussian_table[(rand.v[0] >> 48) & GAUSS_MASK];
    rnd[1] = gaussian_table[(rand.v[0] >> 32) & GAUSS_MASK];
    rnd[2] = gaussian_table[(rand.v[0] >> 16) & GAUSS_MASK];
}

void
gmx_rng_cycle_6gaussian_table(gmx_uint64_t ctr1, gmx_uint64_t ctr2,
                              gmx_uint64_t key1, gmx_uint64_t key2,
                              real* rnd)
{
    threefry2x64_ctr_t ctr  = {{ctr1, ctr2}};
    threefry2x64_key_t key  = {{key1, key2}};
    threefry2x64_ctr_t rand = threefry2x64(ctr, key);

    rnd[0] = gaussian_table[(rand.v[0] >> 48) & GAUSS_MASK];
    rnd[1] = gaussian_table[(rand.v[0] >> 32) & GAUSS_MASK];
    rnd[2] = gaussian_table[(rand.v[0] >> 16) & GAUSS_MASK];
    rnd[3] = gaussian_table[(rand.v[1] >> 48) & GAUSS_MASK];
    rnd[4] = gaussian_table[(rand.v[1] >> 32) & GAUSS_MASK];
    rnd[5] = gaussian_table[(rand.v[1] >> 16) & GAUSS_MASK];
}
