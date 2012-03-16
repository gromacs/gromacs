/*
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 *
 * Gromacs is a library for molecular simulation and trajectory analysis,
 * written by Erik Lindahl, David van der Spoel, Berk Hess, and others - for
 * a full list of developers and information, check out http://www.gromacs.org
 *
 * This program is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the Free 
 * Software Foundation; either version 2 of the License, or (at your option) any 
 * later version.
 * As a special exception, you may use this file as part of a free software
 * library without restriction.  Specifically, if other files instantiate
 * templates or use macros or inline functions from this file, or you compile
 * this file and link it with other files to produce an executable, this
 * file does not by itself cause the resulting executable to be covered by
 * the GNU Lesser General Public License.  
 *
 * In plain-speak: do not worry about classes/macros/templates either - only
 * changes to the library have to be LGPL, not an application linking with it.
 *
 * To help fund GROMACS development, we humbly ask that you cite
 * the papers people have written on it - you can find them on the website!
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#include "types/nrnb.h"
#include "nb_kernel_f77_single.h"
#include "../nb_kerneltype.h"


/* Include standard kernel headers in local directory */
#include "nb_kernel010_f77_single.h"
#include "nb_kernel020_f77_single.h"
#include "nb_kernel030_f77_single.h"
#include "nb_kernel100_f77_single.h"
#include "nb_kernel101_f77_single.h"
#include "nb_kernel102_f77_single.h"
#include "nb_kernel103_f77_single.h"
#include "nb_kernel104_f77_single.h"
#include "nb_kernel110_f77_single.h"
#include "nb_kernel111_f77_single.h"
#include "nb_kernel112_f77_single.h"
#include "nb_kernel113_f77_single.h"
#include "nb_kernel114_f77_single.h"
#include "nb_kernel120_f77_single.h"
#include "nb_kernel121_f77_single.h"
#include "nb_kernel122_f77_single.h"
#include "nb_kernel123_f77_single.h"
#include "nb_kernel124_f77_single.h"
#include "nb_kernel130_f77_single.h"
#include "nb_kernel131_f77_single.h"
#include "nb_kernel132_f77_single.h"
#include "nb_kernel133_f77_single.h"
#include "nb_kernel134_f77_single.h"
#include "nb_kernel200_f77_single.h"
#include "nb_kernel201_f77_single.h"
#include "nb_kernel202_f77_single.h"
#include "nb_kernel203_f77_single.h"
#include "nb_kernel204_f77_single.h"
#include "nb_kernel210_f77_single.h"
#include "nb_kernel211_f77_single.h"
#include "nb_kernel212_f77_single.h"
#include "nb_kernel213_f77_single.h"
#include "nb_kernel214_f77_single.h"
#include "nb_kernel220_f77_single.h"
#include "nb_kernel221_f77_single.h"
#include "nb_kernel222_f77_single.h"
#include "nb_kernel223_f77_single.h"
#include "nb_kernel224_f77_single.h"
#include "nb_kernel230_f77_single.h"
#include "nb_kernel231_f77_single.h"
#include "nb_kernel232_f77_single.h"
#include "nb_kernel233_f77_single.h"
#include "nb_kernel234_f77_single.h"
#include "nb_kernel300_f77_single.h"
#include "nb_kernel301_f77_single.h"
#include "nb_kernel302_f77_single.h"
#include "nb_kernel303_f77_single.h"
#include "nb_kernel304_f77_single.h"
#include "nb_kernel310_f77_single.h"
#include "nb_kernel311_f77_single.h"
#include "nb_kernel312_f77_single.h"
#include "nb_kernel313_f77_single.h"
#include "nb_kernel314_f77_single.h"
#include "nb_kernel320_f77_single.h"
#include "nb_kernel321_f77_single.h"
#include "nb_kernel322_f77_single.h"
#include "nb_kernel323_f77_single.h"
#include "nb_kernel324_f77_single.h"
#include "nb_kernel330_f77_single.h"
#include "nb_kernel331_f77_single.h"
#include "nb_kernel332_f77_single.h"
#include "nb_kernel333_f77_single.h"
#include "nb_kernel334_f77_single.h"
#include "nb_kernel400_f77_single.h"
#include "nb_kernel410_f77_single.h"
#include "nb_kernel420_f77_single.h"
#include "nb_kernel430_f77_single.h"


static nb_kernel_t *
kernellist[eNR_NBKERNEL_NR] = 
{
    nb_kernel010_f77_single,
    nb_kernel020_f77_single,
    nb_kernel030_f77_single,
    nb_kernel100_f77_single,
    nb_kernel101_f77_single,
    nb_kernel102_f77_single,
    nb_kernel103_f77_single,
    nb_kernel104_f77_single,
    nb_kernel110_f77_single,
    nb_kernel111_f77_single,
    nb_kernel112_f77_single,
    nb_kernel113_f77_single,
    nb_kernel114_f77_single,
    nb_kernel120_f77_single,
    nb_kernel121_f77_single,
    nb_kernel122_f77_single,
    nb_kernel123_f77_single,
    nb_kernel124_f77_single,
    nb_kernel130_f77_single,
    nb_kernel131_f77_single,
    nb_kernel132_f77_single,
    nb_kernel133_f77_single,
    nb_kernel134_f77_single,
    nb_kernel200_f77_single,
    nb_kernel201_f77_single,
    nb_kernel202_f77_single,
    nb_kernel203_f77_single,
    nb_kernel204_f77_single,
    nb_kernel210_f77_single,
    nb_kernel211_f77_single,
    nb_kernel212_f77_single,
    nb_kernel213_f77_single,
    nb_kernel214_f77_single,
    nb_kernel220_f77_single,
    nb_kernel221_f77_single,
    nb_kernel222_f77_single,
    nb_kernel223_f77_single,
    nb_kernel224_f77_single,
    nb_kernel230_f77_single,
    nb_kernel231_f77_single,
    nb_kernel232_f77_single,
    nb_kernel233_f77_single,
    nb_kernel234_f77_single,
    nb_kernel300_f77_single,
    nb_kernel301_f77_single,
    nb_kernel302_f77_single,
    nb_kernel303_f77_single,
    nb_kernel304_f77_single,
    nb_kernel310_f77_single,
    nb_kernel311_f77_single,
    nb_kernel312_f77_single,
    nb_kernel313_f77_single,
    nb_kernel314_f77_single,
    nb_kernel320_f77_single,
    nb_kernel321_f77_single,
    nb_kernel322_f77_single,
    nb_kernel323_f77_single,
    nb_kernel324_f77_single,
    nb_kernel330_f77_single,
    nb_kernel331_f77_single,
    nb_kernel332_f77_single,
    nb_kernel333_f77_single,
    nb_kernel334_f77_single,
    nb_kernel400_f77_single,
    nb_kernel410_f77_single,
    nb_kernel430_f77_single,
	nb_kernel010nf_f77_single,
    nb_kernel020nf_f77_single,
    nb_kernel030nf_f77_single,
    nb_kernel100nf_f77_single,
    nb_kernel101nf_f77_single,
    nb_kernel102nf_f77_single,
    nb_kernel103nf_f77_single,
    nb_kernel104nf_f77_single,
    nb_kernel110nf_f77_single,
    nb_kernel111nf_f77_single,
    nb_kernel112nf_f77_single,
    nb_kernel113nf_f77_single,
    nb_kernel114nf_f77_single,
    nb_kernel120nf_f77_single,
    nb_kernel121nf_f77_single,
    nb_kernel122nf_f77_single,
    nb_kernel123nf_f77_single,
    nb_kernel124nf_f77_single,
    nb_kernel130nf_f77_single,
    nb_kernel131nf_f77_single,
    nb_kernel132nf_f77_single,
    nb_kernel133nf_f77_single,
    nb_kernel134nf_f77_single,
    nb_kernel200nf_f77_single,
    nb_kernel201nf_f77_single,
    nb_kernel202nf_f77_single,
    nb_kernel203nf_f77_single,
    nb_kernel204nf_f77_single,
    nb_kernel210nf_f77_single,
    nb_kernel211nf_f77_single,
    nb_kernel212nf_f77_single,
    nb_kernel213nf_f77_single,
    nb_kernel214nf_f77_single,
    nb_kernel220nf_f77_single,
    nb_kernel221nf_f77_single,
    nb_kernel222nf_f77_single,
    nb_kernel223nf_f77_single,
    nb_kernel224nf_f77_single,
    nb_kernel230nf_f77_single,
    nb_kernel231nf_f77_single,
    nb_kernel232nf_f77_single,
    nb_kernel233nf_f77_single,
    nb_kernel234nf_f77_single,
    nb_kernel300nf_f77_single,
    nb_kernel301nf_f77_single,
    nb_kernel302nf_f77_single,
    nb_kernel303nf_f77_single,
    nb_kernel304nf_f77_single,
    nb_kernel310nf_f77_single,
    nb_kernel311nf_f77_single,
    nb_kernel312nf_f77_single,
    nb_kernel313nf_f77_single,
    nb_kernel314nf_f77_single,
    nb_kernel320nf_f77_single,
    nb_kernel321nf_f77_single,
    nb_kernel322nf_f77_single,
    nb_kernel323nf_f77_single,
    nb_kernel324nf_f77_single,
    nb_kernel330nf_f77_single,
    nb_kernel331nf_f77_single,
    nb_kernel332nf_f77_single,
    nb_kernel333nf_f77_single,
    nb_kernel334nf_f77_single,
    nb_kernel400nf_f77_single,
    nb_kernel410nf_f77_single,
    nb_kernel430nf_f77_single,
};


void
nb_kernel_setup_f77_single(FILE *log,nb_kernel_t **list)
{
  int i;
  nb_kernel_t *p;

  if(log)
    fprintf(log,"Configuring single precision Fortran kernels...\n");

  for(i=0;i<eNR_NBKERNEL_NR;i++)
  {
    p = kernellist[i];
    if(p!=NULL)
      list[i] = p; 
  }
}    
