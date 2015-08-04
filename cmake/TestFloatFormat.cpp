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
#include <stdio.h>

volatile const double abc [10] = {
    /* Zero-terminated strings encoded as floating-point numbers */
    /* "GROMACSX" in ascii    */
    (double)  3.80279098314984902657e+35, (double) 0.0,
    /* "GROMACSX" in ebcdic   */
    (double) -1.37384666579378297437e+38, (double) 0.0,
    /* "D__float" (vax)       */
    (double)  3.53802595280598432000e+18, (double) 0.0,
    /* "IBMHEXFP" s390/ascii  */
    (double)  1.77977764695171661377e+10, (double) 0.0,
    /* "IBMHEXFP" s390/ebcdic */
    (double) -5.22995989424860458374e+10, (double) 0.0,
};

/* Check that a double is 8 bytes - compilation dies if it isnt */
extern char xyz [sizeof(double) == 8 ? 1 : -1];

int
main()
{
    int         i;
    double      d;

    /* Make sure some compilers do not optimize away the entire structure
     * with floating-point data by using it to produce a return value.
     */
    for (i = 0, d = 0; i < 10; i++)
    {
        d += abc[i];
    }
    return (d == 12345.0);
}
