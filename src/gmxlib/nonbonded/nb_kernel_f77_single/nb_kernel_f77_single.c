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
#include "nbkernel010_f77_single.h"
#include "nbkernel020_f77_single.h"
#include "nbkernel030_f77_single.h"
#include "nbkernel100_f77_single.h"
#include "nbkernel101_f77_single.h"
#include "nbkernel102_f77_single.h"
#include "nbkernel103_f77_single.h"
#include "nbkernel104_f77_single.h"
#include "nbkernel110_f77_single.h"
#include "nbkernel111_f77_single.h"
#include "nbkernel112_f77_single.h"
#include "nbkernel113_f77_single.h"
#include "nbkernel114_f77_single.h"
#include "nbkernel120_f77_single.h"
#include "nbkernel121_f77_single.h"
#include "nbkernel122_f77_single.h"
#include "nbkernel123_f77_single.h"
#include "nbkernel124_f77_single.h"
#include "nbkernel130_f77_single.h"
#include "nbkernel131_f77_single.h"
#include "nbkernel132_f77_single.h"
#include "nbkernel133_f77_single.h"
#include "nbkernel134_f77_single.h"
#include "nbkernel200_f77_single.h"
#include "nbkernel201_f77_single.h"
#include "nbkernel202_f77_single.h"
#include "nbkernel203_f77_single.h"
#include "nbkernel204_f77_single.h"
#include "nbkernel210_f77_single.h"
#include "nbkernel211_f77_single.h"
#include "nbkernel212_f77_single.h"
#include "nbkernel213_f77_single.h"
#include "nbkernel214_f77_single.h"
#include "nbkernel220_f77_single.h"
#include "nbkernel221_f77_single.h"
#include "nbkernel222_f77_single.h"
#include "nbkernel223_f77_single.h"
#include "nbkernel224_f77_single.h"
#include "nbkernel230_f77_single.h"
#include "nbkernel231_f77_single.h"
#include "nbkernel232_f77_single.h"
#include "nbkernel233_f77_single.h"
#include "nbkernel234_f77_single.h"
#include "nbkernel300_f77_single.h"
#include "nbkernel301_f77_single.h"
#include "nbkernel302_f77_single.h"
#include "nbkernel303_f77_single.h"
#include "nbkernel304_f77_single.h"
#include "nbkernel310_f77_single.h"
#include "nbkernel311_f77_single.h"
#include "nbkernel312_f77_single.h"
#include "nbkernel313_f77_single.h"
#include "nbkernel314_f77_single.h"
#include "nbkernel320_f77_single.h"
#include "nbkernel321_f77_single.h"
#include "nbkernel322_f77_single.h"
#include "nbkernel323_f77_single.h"
#include "nbkernel324_f77_single.h"
#include "nbkernel330_f77_single.h"
#include "nbkernel331_f77_single.h"
#include "nbkernel332_f77_single.h"
#include "nbkernel333_f77_single.h"
#include "nbkernel334_f77_single.h"
#include "nbkernel400_f77_single.h"
#include "nbkernel410_f77_single.h"
#include "nbkernel420_f77_single.h"
#include "nbkernel430_f77_single.h"


static nb_kernel_t *
kernellist[eNR_NBKERNEL_NR] = 
{
    nbkernel010_f77_single,
    nbkernel020_f77_single,
    nbkernel030_f77_single,
    nbkernel100_f77_single,
    nbkernel101_f77_single,
    nbkernel102_f77_single,
    nbkernel103_f77_single,
    nbkernel104_f77_single,
    nbkernel110_f77_single,
    nbkernel111_f77_single,
    nbkernel112_f77_single,
    nbkernel113_f77_single,
    nbkernel114_f77_single,
    nbkernel120_f77_single,
    nbkernel121_f77_single,
    nbkernel122_f77_single,
    nbkernel123_f77_single,
    nbkernel124_f77_single,
    nbkernel130_f77_single,
    nbkernel131_f77_single,
    nbkernel132_f77_single,
    nbkernel133_f77_single,
    nbkernel134_f77_single,
    nbkernel200_f77_single,
    nbkernel201_f77_single,
    nbkernel202_f77_single,
    nbkernel203_f77_single,
    nbkernel204_f77_single,
    nbkernel210_f77_single,
    nbkernel211_f77_single,
    nbkernel212_f77_single,
    nbkernel213_f77_single,
    nbkernel214_f77_single,
    nbkernel220_f77_single,
    nbkernel221_f77_single,
    nbkernel222_f77_single,
    nbkernel223_f77_single,
    nbkernel224_f77_single,
    nbkernel230_f77_single,
    nbkernel231_f77_single,
    nbkernel232_f77_single,
    nbkernel233_f77_single,
    nbkernel234_f77_single,
    nbkernel300_f77_single,
    nbkernel301_f77_single,
    nbkernel302_f77_single,
    nbkernel303_f77_single,
    nbkernel304_f77_single,
    nbkernel310_f77_single,
    nbkernel311_f77_single,
    nbkernel312_f77_single,
    nbkernel313_f77_single,
    nbkernel314_f77_single,
    nbkernel320_f77_single,
    nbkernel321_f77_single,
    nbkernel322_f77_single,
    nbkernel323_f77_single,
    nbkernel324_f77_single,
    nbkernel330_f77_single,
    nbkernel331_f77_single,
    nbkernel332_f77_single,
    nbkernel333_f77_single,
    nbkernel334_f77_single,
    nbkernel400_f77_single,
    nbkernel410_f77_single,
    nbkernel430_f77_single,
	nbkernel010nf_f77_single,
    nbkernel020nf_f77_single,
    nbkernel030nf_f77_single,
    nbkernel100nf_f77_single,
    nbkernel101nf_f77_single,
    nbkernel102nf_f77_single,
    nbkernel103nf_f77_single,
    nbkernel104nf_f77_single,
    nbkernel110nf_f77_single,
    nbkernel111nf_f77_single,
    nbkernel112nf_f77_single,
    nbkernel113nf_f77_single,
    nbkernel114nf_f77_single,
    nbkernel120nf_f77_single,
    nbkernel121nf_f77_single,
    nbkernel122nf_f77_single,
    nbkernel123nf_f77_single,
    nbkernel124nf_f77_single,
    nbkernel130nf_f77_single,
    nbkernel131nf_f77_single,
    nbkernel132nf_f77_single,
    nbkernel133nf_f77_single,
    nbkernel134nf_f77_single,
    nbkernel200nf_f77_single,
    nbkernel201nf_f77_single,
    nbkernel202nf_f77_single,
    nbkernel203nf_f77_single,
    nbkernel204nf_f77_single,
    nbkernel210nf_f77_single,
    nbkernel211nf_f77_single,
    nbkernel212nf_f77_single,
    nbkernel213nf_f77_single,
    nbkernel214nf_f77_single,
    nbkernel220nf_f77_single,
    nbkernel221nf_f77_single,
    nbkernel222nf_f77_single,
    nbkernel223nf_f77_single,
    nbkernel224nf_f77_single,
    nbkernel230nf_f77_single,
    nbkernel231nf_f77_single,
    nbkernel232nf_f77_single,
    nbkernel233nf_f77_single,
    nbkernel234nf_f77_single,
    nbkernel300nf_f77_single,
    nbkernel301nf_f77_single,
    nbkernel302nf_f77_single,
    nbkernel303nf_f77_single,
    nbkernel304nf_f77_single,
    nbkernel310nf_f77_single,
    nbkernel311nf_f77_single,
    nbkernel312nf_f77_single,
    nbkernel313nf_f77_single,
    nbkernel314nf_f77_single,
    nbkernel320nf_f77_single,
    nbkernel321nf_f77_single,
    nbkernel322nf_f77_single,
    nbkernel323nf_f77_single,
    nbkernel324nf_f77_single,
    nbkernel330nf_f77_single,
    nbkernel331nf_f77_single,
    nbkernel332nf_f77_single,
    nbkernel333nf_f77_single,
    nbkernel334nf_f77_single,
    nbkernel400nf_f77_single,
    nbkernel410nf_f77_single,
    nbkernel430nf_f77_single,
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
