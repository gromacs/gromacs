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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _writeps_h
#define _writeps_h


#include <stdio.h>
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum { 
  eXCenter, eXLeft, eXRight
} eXPos;

typedef enum { 
  eYCenter, eYTop,  eYBottom
} eYPos; 

enum { 
  efontTIMES, efontTIMESITALIC, efontTIMESBOLD, efontTIMESBOLDITALIC,
  efontHELV,  efontHELVITALIC,  efontHELVBOLD,  efontHELVBOLDITALIC,
  efontCOUR,  efontCOURITALIC,  efontCOURBOLD,  efontCOURBOLDITALIC,
  efontNR };


typedef struct t_int_psdata *t_psdata;
/* Only use t_psdata - it is a pointer to an abstract datatype
 * that maintains the state of the postscript currently written.
 */

extern const char *fontnm[efontNR];

t_psdata ps_open(const char *fn,real x1,real y1,real x2,real y2);

void ps_linewidth(t_psdata ps,int lw);
void ps_color(t_psdata ps,real r,real g,real b);
void ps_rgb(t_psdata ps,t_rgb *rgb);

void ps_rgb_box(t_psdata ps,t_rgb *rgb);
void ps_rgb_nbox(t_psdata ps,t_rgb *rgb,real n);
void ps_init_rgb_box(t_psdata ps,real xbox, real ybox);
void ps_init_rgb_nbox(t_psdata ps,real xbox, real ybox);

void ps_lineto(t_psdata ps,real x,real y);
void ps_linerel(t_psdata ps,real dx,real dy);

void ps_moveto(t_psdata ps,real x,real y);
void ps_moverel(t_psdata ps,real dx,real dy);

void ps_line(t_psdata ps,real x1,real y1,real x2,real y2);

void ps_box(t_psdata ps,real x1,real y1,real x2,real y2);
void ps_fillbox(t_psdata ps,real x1,real y1,real x2,real y2);

void ps_arc(t_psdata ps,real x1,real y1,real rad,real a0,real a1);
void ps_fillarc(t_psdata ps,real x1,real y1,real rad,real a0,real a1);
void ps_arcslice(t_psdata ps,real xc,real yc,
			real rad1,real rad2,real a0,real a1);
void ps_fillarcslice(t_psdata ps,real xc,real yc,
			    real rad1,real rad2,real a0,real a1);

void ps_circle(t_psdata ps,real x1,real y1,real rad);

void ps_font(t_psdata ps,int font,real size);
void ps_strfont(t_psdata ps,char *font,real size);

void ps_text(t_psdata ps,real x1,real y1,const char *str);
void ps_ctext(t_psdata ps,real x1,real y1,const char *str,int expos);

void ps_close(t_psdata ps);

void ps_flip(t_psdata ps,gmx_bool bPlus);
/* Rotate over 90 (bPlus) or -90 (!bPlus) degrees */

void ps_rotate(t_psdata ps,real angle);

void ps_translate(t_psdata ps,real x,real y);

void ps_setorigin(t_psdata ps);
void ps_unsetorigin(t_psdata ps);

void viewps(char *fn);

void ps_comment(t_psdata ps,const char *s);

#ifdef __cplusplus
}
#endif

#endif	/* _writeps_h */

