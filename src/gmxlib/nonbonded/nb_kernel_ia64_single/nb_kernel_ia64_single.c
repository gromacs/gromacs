/* -*- mode: c; tab-width: 4; indent-tabs-mode: n; c-basic-offset: 4 -*- 
 *
 * $Id$
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2004
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */

#include <stdio.h>

#include <types/nrnb.h>

#include "../nb_kerneltype.h"

/* Include kernel headers in local directory.
* We can only have one routine in each file due to a bug
* in the intel assembler program...
*/
#include "nb_kernel010_ia64_single.h"
#include "nb_kernel010nf_ia64_single.h"
#include "nb_kernel030_ia64_single.h"
#include "nb_kernel030nf_ia64_single.h"
#include "nb_kernel100_ia64_single.h"
#include "nb_kernel100nf_ia64_single.h"
#include "nb_kernel110_ia64_single.h"
#include "nb_kernel110nf_ia64_single.h"
#include "nb_kernel130_ia64_single.h"
#include "nb_kernel130nf_ia64_single.h"
#include "nb_kernel200_ia64_single.h"
#include "nb_kernel200nf_ia64_single.h"
#include "nb_kernel210_ia64_single.h"
#include "nb_kernel210nf_ia64_single.h"
#include "nb_kernel230_ia64_single.h"
#include "nb_kernel230nf_ia64_single.h"
#include "nb_kernel300_ia64_single.h"
#include "nb_kernel300nf_ia64_single.h"
#include "nb_kernel310_ia64_single.h"
#include "nb_kernel310nf_ia64_single.h"
#include "nb_kernel330_ia64_single.h"
#include "nb_kernel330nf_ia64_single.h"
#include "nb_kernel400_ia64_single.h"
#include "nb_kernel400nf_ia64_single.h"
#include "nb_kernel410_ia64_single.h"
#include "nb_kernel410nf_ia64_single.h"
#include "nb_kernel430_ia64_single.h"
#include "nb_kernel430nf_ia64_single.h"



static nb_kernel_t *
kernellist_ia64_single[eNR_NBKERNEL_NR] = 
{
    nb_kernel010_ia64_single,
    NULL,
    nb_kernel030_ia64_single,
    nb_kernel100_ia64_single,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel110_ia64_single,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel130_ia64_single,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel200_ia64_single,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel210_ia64_single,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel230_ia64_single,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel300_ia64_single,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel310_ia64_single,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel330_ia64_single,
    NULL,
    NULL,
    NULL,
    NULL,
    nb_kernel400_ia64_single,
    nb_kernel410_ia64_single,
    nb_kernel430_ia64_single    
};




void
nb_kernel_setup_ia64_single(FILE *log,nb_kernel_t **list)
{
    int i;
    nb_kernel_t *p;
        
	if(log)
        fprintf(log,"Using single precision ia64 assembly kernels.");
	
    for(i=0;i<eNR_NBKERNEL_NR;i++)
    {
        p = kernellist_ia64_single[i];
        if(p!=NULL)
            list[i] = p; 
    }
}    









