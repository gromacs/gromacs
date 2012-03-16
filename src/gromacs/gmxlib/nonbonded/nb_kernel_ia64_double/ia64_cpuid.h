/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

/* IA64 provides 4 different CPUID registers,
 * CPUID[x], where is from 0 to at least 4.
 * The argument to this function is the index of the desired
 * register, whose contents will be returned as a 64-bit integer.
 *
 * Normally we are interested in the CPU generation numbers in
 * CPUID[3]. The contents of this register is:
 *
 * Bits 63-40: Reserved
 * Bits 39-32: Architecture revision.
 * Bits 31-24: Processor family number.
 * Bits 23-16: Processor model number.
 * Bits 15-8:  Processor revision number.
 * Bits 7-0:   Number of available CPUID registers.
 *
 * As an example, my 1.3GHz Madison has
 */
unsigned long long
ia64_cpuid(int reg);

