/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
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
#ifndef _LAPACK_LIMITS_H_
#define _LAPACK_LIMITS_H_

#define DSTEBZ_BLOCKSIZE  1

#define DORGBR_BLOCKSIZE    32
#define DORGBR_MINBLOCKSIZE 2
#define DORGBR_CROSSOVER    128

#define DGEQRF_BLOCKSIZE    32
#define DGEQRF_MINBLOCKSIZE 2  
#define DGEQRF_CROSSOVER    128

#define DORGQR_BLOCKSIZE    32
#define DORGQR_MINBLOCKSIZE 2
#define DORGQR_CROSSOVER    128

#define DORMLQ_BLOCKSIZE    32
#define DORMLQ_MINBLOCKSIZE 2
#define DORMLQ_CROSSOVER    128

#define DORMQL_BLOCKSIZE    32
#define DORMQL_MINBLOCKSIZE 2
#define DORMQL_CROSSOVER    128

#define DSYTRD_BLOCKSIZE    32
#define DSYTRD_MINBLOCKSIZE 2
#define DSYTRD_CROSSOVER    128

#define DGEBRD_BLOCKSIZE    32
#define DGEBRD_MINBLOCKSIZE 2
#define DGEBRD_CROSSOVER    128

#define DORMQR_BLOCKSIZE    32
#define DORMQR_MINBLOCKSIZE 2
#define DORMQR_CROSSOVER    128

#define DGELQF_BLOCKSIZE    32
#define DGELQF_MINBLOCKSIZE 2  
#define DGELQF_CROSSOVER    128

#define DGETRF_BLOCKSIZE    64
#define DGETRF_MINBLOCKSIZE 2

#define DGETRI_BLOCKSIZE    64
#define DGETRI_MINBLOCKSIZE 2

#define DTRTRI_BLOCKSIZE    64

#define DBDSDC_SMALLSIZE 25

#define DBDSQR_MAXITR 6



#endif /* _LAPACK_LIMITS_H_ */
