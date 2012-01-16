/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * GRoups of Organic Molecules in ACtion for Science
 */
#ifndef _globsig_h
#define _globsig_h

#ifdef __cplusplus
extern "C" {
#endif

#if 0
}
/* Hack to make automatic indenting work */
#endif

/*! \brief Enumeration of the types of inter-process signal that can
 *         occur during a simulation.
 *
 * The esignalXXX_NR elements allow for looping over only some types
 * of signal.
 */
enum {
    /* Intra-simulation signal types that are local to this
     * simulation. These used to be called "global" signals. */
    esignalNABNSB,
    esignalRESETCOUNTERS,
    esignalINTRA_NR,
    /* Inter-simulation signal types that might require coordinating
     * multiple simulations. */
    esignalCHKPT = esignalINTRA_NR,
    esignalSTOPCOND,
    esignalINTER_NR,
    esignalNR = esignalINTER_NR
};

/*! \brief Under various conditions, processes may need to signal
 *         others. This can occur in an inter-simulation or
 *         intra-simulation sense.
 *
 * Keep in mind that they are transmitted to other processes through
 * an MPI_Reduce after casting them to a real (so the signals can be
 * sent together with other data). This means that the only meaningful
 * input values are positive, negative or zero, because the output
 * values are sums over the local values of these.
 */

typedef struct {
    int nst_signal_intra;
    /* The number of steps between intra-simulation communication. */
    int nst_signal_inter;
    /* The number of steps between inter-simulation communication. */

    /* The next two data members are arrays of signal values. These
     * are 0 when not set. */

    int init[esignalNR];
    /* The signal initiated by this process, which will be transferred
     * to the set array during inter-simulation communication, and
     * then reset to zero here. */
    int set[esignalNR];
    /* The communicated signal, always equal for all
     * processes. Updated only during inter-simulation
     * communication. */
} gmx_signal;

#ifdef __cplusplus
}
#endif

#endif
