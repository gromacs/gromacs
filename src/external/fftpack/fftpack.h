/*

 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation

 ************************************************************/

#ifndef _fftpack_h
#define _fftpack_h

#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

#define Treal real

    void fftpack_cffti1(int n, Treal wa[], int ifac[]);
    void fftpack_cfftf1(int n, Treal c[], Treal ch[], const Treal wa[], const int ifac[], int isign);
    void fftpack_rffti1(int n, Treal wa[], int ifac[]);
    void fftpack_rfftf1(int n, Treal c[], Treal ch[], const Treal wa[], const int ifac[]);
    void fftpack_rfftb1(int n, Treal c[], Treal ch[], const Treal wa[], const int ifac[]);

#ifdef __cplusplus
}
#endif

#endif
