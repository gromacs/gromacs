/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
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
 * Groningen Machine for Chemical Simulation
 */
#include <stdio.h>

#include <types/nrnb.h>

#include "../nb_kerneltype.h"

#include "nb_kernel.h"

/* Include standard kernel headers in local directory */
#include "nb_kernel010.h"
#include "nb_kernel020.h"
#include "nb_kernel030.h"
#include "nb_kernel100.h"
#include "nb_kernel101.h"
#include "nb_kernel102.h"
#include "nb_kernel103.h"
#include "nb_kernel104.h"
#include "nb_kernel110.h"
#include "nb_kernel111.h"
#include "nb_kernel112.h"
#include "nb_kernel113.h"
#include "nb_kernel114.h"
#include "nb_kernel120.h"
#include "nb_kernel121.h"
#include "nb_kernel122.h"
#include "nb_kernel123.h"
#include "nb_kernel124.h"
#include "nb_kernel130.h"
#include "nb_kernel131.h"
#include "nb_kernel132.h"
#include "nb_kernel133.h"
#include "nb_kernel134.h"
#include "nb_kernel200.h"
#include "nb_kernel201.h"
#include "nb_kernel202.h"
#include "nb_kernel203.h"
#include "nb_kernel204.h"
#include "nb_kernel210.h"
#include "nb_kernel211.h"
#include "nb_kernel212.h"
#include "nb_kernel213.h"
#include "nb_kernel214.h"
#include "nb_kernel220.h"
#include "nb_kernel221.h"
#include "nb_kernel222.h"
#include "nb_kernel223.h"
#include "nb_kernel224.h"
#include "nb_kernel230.h"
#include "nb_kernel231.h"
#include "nb_kernel232.h"
#include "nb_kernel233.h"
#include "nb_kernel234.h"
#include "nb_kernel300.h"
#include "nb_kernel301.h"
#include "nb_kernel302.h"
#include "nb_kernel303.h"
#include "nb_kernel304.h"
#include "nb_kernel310.h"
#include "nb_kernel311.h"
#include "nb_kernel312.h"
#include "nb_kernel313.h"
#include "nb_kernel314.h"
#include "nb_kernel320.h"
#include "nb_kernel321.h"
#include "nb_kernel322.h"
#include "nb_kernel323.h"
#include "nb_kernel324.h"
#include "nb_kernel330.h"
#include "nb_kernel331.h"
#include "nb_kernel332.h"
#include "nb_kernel333.h"
#include "nb_kernel334.h"
#include "nb_kernel400.h"
#include "nb_kernel410.h"
#include "nb_kernel420.h"
#include "nb_kernel430.h"


static nb_kernel_t *
kernellist[eNR_NBKERNEL_NR] = 
{
    F77_OR_C_FUNC_(nb_kernel010,NB_KERNEL010),
    F77_OR_C_FUNC_(nb_kernel020,NB_KERNEL020),
    F77_OR_C_FUNC_(nb_kernel030,NB_KERNEL030),
    F77_OR_C_FUNC_(nb_kernel100,NB_KERNEL100),
    F77_OR_C_FUNC_(nb_kernel101,NB_KERNEL101),
    F77_OR_C_FUNC_(nb_kernel102,NB_KERNEL102),
    F77_OR_C_FUNC_(nb_kernel103,NB_KERNEL103),
    F77_OR_C_FUNC_(nb_kernel104,NB_KERNEL104),
    F77_OR_C_FUNC_(nb_kernel110,NB_KERNEL110),
    F77_OR_C_FUNC_(nb_kernel111,NB_KERNEL111),
    F77_OR_C_FUNC_(nb_kernel112,NB_KERNEL112),
    F77_OR_C_FUNC_(nb_kernel113,NB_KERNEL113),
    F77_OR_C_FUNC_(nb_kernel114,NB_KERNEL114),
    F77_OR_C_FUNC_(nb_kernel120,NB_KERNEL120),
    F77_OR_C_FUNC_(nb_kernel121,NB_KERNEL121),
    F77_OR_C_FUNC_(nb_kernel122,NB_KERNEL122),
    F77_OR_C_FUNC_(nb_kernel123,NB_KERNEL123),
    F77_OR_C_FUNC_(nb_kernel124,NB_KERNEL124),
    F77_OR_C_FUNC_(nb_kernel130,NB_KERNEL130),
    F77_OR_C_FUNC_(nb_kernel131,NB_KERNEL131),
    F77_OR_C_FUNC_(nb_kernel132,NB_KERNEL132),
    F77_OR_C_FUNC_(nb_kernel133,NB_KERNEL133),
    F77_OR_C_FUNC_(nb_kernel134,NB_KERNEL134),
    F77_OR_C_FUNC_(nb_kernel200,NB_KERNEL200),
    F77_OR_C_FUNC_(nb_kernel201,NB_KERNEL201),
    F77_OR_C_FUNC_(nb_kernel202,NB_KERNEL202),
    F77_OR_C_FUNC_(nb_kernel203,NB_KERNEL203),
    F77_OR_C_FUNC_(nb_kernel204,NB_KERNEL204),
    F77_OR_C_FUNC_(nb_kernel210,NB_KERNEL210),
    F77_OR_C_FUNC_(nb_kernel211,NB_KERNEL211),
    F77_OR_C_FUNC_(nb_kernel212,NB_KERNEL212),
    F77_OR_C_FUNC_(nb_kernel213,NB_KERNEL213),
    F77_OR_C_FUNC_(nb_kernel214,NB_KERNEL214),
    F77_OR_C_FUNC_(nb_kernel220,NB_KERNEL220),
    F77_OR_C_FUNC_(nb_kernel221,NB_KERNEL221),
    F77_OR_C_FUNC_(nb_kernel222,NB_KERNEL222),
    F77_OR_C_FUNC_(nb_kernel223,NB_KERNEL223),
    F77_OR_C_FUNC_(nb_kernel224,NB_KERNEL224),
    F77_OR_C_FUNC_(nb_kernel230,NB_KERNEL230),
    F77_OR_C_FUNC_(nb_kernel231,NB_KERNEL231),
    F77_OR_C_FUNC_(nb_kernel232,NB_KERNEL232),
    F77_OR_C_FUNC_(nb_kernel233,NB_KERNEL233),
    F77_OR_C_FUNC_(nb_kernel234,NB_KERNEL234),
    F77_OR_C_FUNC_(nb_kernel300,NB_KERNEL300),
    F77_OR_C_FUNC_(nb_kernel301,NB_KERNEL301),
    F77_OR_C_FUNC_(nb_kernel302,NB_KERNEL302),
    F77_OR_C_FUNC_(nb_kernel303,NB_KERNEL303),
    F77_OR_C_FUNC_(nb_kernel304,NB_KERNEL304),
    F77_OR_C_FUNC_(nb_kernel310,NB_KERNEL310),
    F77_OR_C_FUNC_(nb_kernel311,NB_KERNEL311),
    F77_OR_C_FUNC_(nb_kernel312,NB_KERNEL312),
    F77_OR_C_FUNC_(nb_kernel313,NB_KERNEL313),
    F77_OR_C_FUNC_(nb_kernel314,NB_KERNEL314),
    F77_OR_C_FUNC_(nb_kernel320,NB_KERNEL320),
    F77_OR_C_FUNC_(nb_kernel321,NB_KERNEL321),
    F77_OR_C_FUNC_(nb_kernel322,NB_KERNEL322),
    F77_OR_C_FUNC_(nb_kernel323,NB_KERNEL323),
    F77_OR_C_FUNC_(nb_kernel324,NB_KERNEL324),
    F77_OR_C_FUNC_(nb_kernel330,NB_KERNEL330),
    F77_OR_C_FUNC_(nb_kernel331,NB_KERNEL331),
    F77_OR_C_FUNC_(nb_kernel332,NB_KERNEL332),
    F77_OR_C_FUNC_(nb_kernel333,NB_KERNEL333),
    F77_OR_C_FUNC_(nb_kernel334,NB_KERNEL334),
    F77_OR_C_FUNC_(nb_kernel400,NB_KERNEL400),
    F77_OR_C_FUNC_(nb_kernel410,NB_KERNEL410),
    F77_OR_C_FUNC_(nb_kernel430,NB_KERNEL430),
	F77_OR_C_FUNC_(nb_kernel010nf,NB_KERNEL010NF),
    F77_OR_C_FUNC_(nb_kernel020nf,NB_KERNEL020NF),
    F77_OR_C_FUNC_(nb_kernel030nf,NB_KERNEL030NF),
    F77_OR_C_FUNC_(nb_kernel100nf,NB_KERNEL100NF),
    F77_OR_C_FUNC_(nb_kernel101nf,NB_KERNEL101NF),
    F77_OR_C_FUNC_(nb_kernel102nf,NB_KERNEL102NF),
    F77_OR_C_FUNC_(nb_kernel103nf,NB_KERNEL103NF),
    F77_OR_C_FUNC_(nb_kernel104nf,NB_KERNEL104NF),
    F77_OR_C_FUNC_(nb_kernel110nf,NB_KERNEL110NF),
    F77_OR_C_FUNC_(nb_kernel111nf,NB_KERNEL111NF),
    F77_OR_C_FUNC_(nb_kernel112nf,NB_KERNEL112NF),
    F77_OR_C_FUNC_(nb_kernel113nf,NB_KERNEL113NF),
    F77_OR_C_FUNC_(nb_kernel114nf,NB_KERNEL114NF),
    F77_OR_C_FUNC_(nb_kernel120nf,NB_KERNEL120NF),
    F77_OR_C_FUNC_(nb_kernel121nf,NB_KERNEL121NF),
    F77_OR_C_FUNC_(nb_kernel122nf,NB_KERNEL122NF),
    F77_OR_C_FUNC_(nb_kernel123nf,NB_KERNEL123NF),
    F77_OR_C_FUNC_(nb_kernel124nf,NB_KERNEL124NF),
    F77_OR_C_FUNC_(nb_kernel130nf,NB_KERNEL130NF),
    F77_OR_C_FUNC_(nb_kernel131nf,NB_KERNEL131NF),
    F77_OR_C_FUNC_(nb_kernel132nf,NB_KERNEL132NF),
    F77_OR_C_FUNC_(nb_kernel133nf,NB_KERNEL133NF),
    F77_OR_C_FUNC_(nb_kernel134nf,NB_KERNEL134NF),
    F77_OR_C_FUNC_(nb_kernel200nf,NB_KERNEL200NF),
    F77_OR_C_FUNC_(nb_kernel201nf,NB_KERNEL201NF),
    F77_OR_C_FUNC_(nb_kernel202nf,NB_KERNEL202NF),
    F77_OR_C_FUNC_(nb_kernel203nf,NB_KERNEL203NF),
    F77_OR_C_FUNC_(nb_kernel204nf,NB_KERNEL204NF),
    F77_OR_C_FUNC_(nb_kernel210nf,NB_KERNEL210NF),
    F77_OR_C_FUNC_(nb_kernel211nf,NB_KERNEL211NF),
    F77_OR_C_FUNC_(nb_kernel212nf,NB_KERNEL212NF),
    F77_OR_C_FUNC_(nb_kernel213nf,NB_KERNEL213NF),
    F77_OR_C_FUNC_(nb_kernel214nf,NB_KERNEL214NF),
    F77_OR_C_FUNC_(nb_kernel220nf,NB_KERNEL220NF),
    F77_OR_C_FUNC_(nb_kernel221nf,NB_KERNEL221NF),
    F77_OR_C_FUNC_(nb_kernel222nf,NB_KERNEL222NF),
    F77_OR_C_FUNC_(nb_kernel223nf,NB_KERNEL223NF),
    F77_OR_C_FUNC_(nb_kernel224nf,NB_KERNEL224NF),
    F77_OR_C_FUNC_(nb_kernel230nf,NB_KERNEL230NF),
    F77_OR_C_FUNC_(nb_kernel231nf,NB_KERNEL231NF),
    F77_OR_C_FUNC_(nb_kernel232nf,NB_KERNEL232NF),
    F77_OR_C_FUNC_(nb_kernel233nf,NB_KERNEL233NF),
    F77_OR_C_FUNC_(nb_kernel234nf,NB_KERNEL234NF),
    F77_OR_C_FUNC_(nb_kernel300nf,NB_KERNEL300NF),
    F77_OR_C_FUNC_(nb_kernel301nf,NB_KERNEL301NF),
    F77_OR_C_FUNC_(nb_kernel302nf,NB_KERNEL302NF),
    F77_OR_C_FUNC_(nb_kernel303nf,NB_KERNEL303NF),
    F77_OR_C_FUNC_(nb_kernel304nf,NB_KERNEL304NF),
    F77_OR_C_FUNC_(nb_kernel310nf,NB_KERNEL310NF),
    F77_OR_C_FUNC_(nb_kernel311nf,NB_KERNEL311NF),
    F77_OR_C_FUNC_(nb_kernel312nf,NB_KERNEL312NF),
    F77_OR_C_FUNC_(nb_kernel313nf,NB_KERNEL313NF),
    F77_OR_C_FUNC_(nb_kernel314nf,NB_KERNEL314NF),
    F77_OR_C_FUNC_(nb_kernel320nf,NB_KERNEL320NF),
    F77_OR_C_FUNC_(nb_kernel321nf,NB_KERNEL321NF),
    F77_OR_C_FUNC_(nb_kernel322nf,NB_KERNEL322NF),
    F77_OR_C_FUNC_(nb_kernel323nf,NB_KERNEL323NF),
    F77_OR_C_FUNC_(nb_kernel324nf,NB_KERNEL324NF),
    F77_OR_C_FUNC_(nb_kernel330nf,NB_KERNEL330NF),
    F77_OR_C_FUNC_(nb_kernel331nf,NB_KERNEL331NF),
    F77_OR_C_FUNC_(nb_kernel332nf,NB_KERNEL332NF),
    F77_OR_C_FUNC_(nb_kernel333nf,NB_KERNEL333NF),
    F77_OR_C_FUNC_(nb_kernel334nf,NB_KERNEL334NF),
    F77_OR_C_FUNC_(nb_kernel400nf,NB_KERNEL400NF),
    F77_OR_C_FUNC_(nb_kernel410nf,NB_KERNEL410NF),
    F77_OR_C_FUNC_(nb_kernel430nf,NB_KERNEL430NF)
};


void
nb_kernel_setup(FILE *log,nb_kernel_t **list)
{
	int i;
    nb_kernel_t *p;

	for(i=0;i<eNR_NBKERNEL_NR;i++)
    {
        p = kernellist[i];
        if(p!=NULL)
            list[i] = p; 
	}
}    
