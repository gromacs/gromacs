/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
 */

#ifndef _corrutil_h
#define _corrutil_h

static char *SRCID_corrutil_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) corrutil.h 1.22 9/30/97"
#endif /* HAVE_IDENT */
#include "sysstuff.h"
#include "statutil.h"

#define FACTOR  1000.0	/* Convert nm^2/ps to 10e-5 cm^2/s */

extern float dim_factor;

class c_corr {
public:
  int  natoms,nrestart;
  real t0,delta_t;
  int  nframes,maxframes,nlast,ngrp;
  int  *n_offs;
  int  **ndata;
  real **data,*time,*data_x,*data_y,*data_z,*data_xy;
  rvec **x0;
public:
  /* Constructor, destructor */
  c_corr(int nrgrp,int nrframes);
  virtual ~c_corr();
  
  /* Normal methods */
  void normalise();
  void init_restart();
  void print(char *fn,char *title,char *yaxis,bool bXvgr);
  void subtitle(FILE *out);
  void calc(int nr,int nx,atom_id index[],rvec xc[]);
  real thistime() { return time[nframes]; };
  void loop(char *fn,int gnx[],atom_id *index[],
	    t_first_x *fx,t_next_x *nx);
  int  in_data(int nx00) { return nframes-n_offs[nx00]; };
  
  /* Virtual methods, may/must be overridden by children */
  virtual real calc1(int nx,atom_id index[],int nx0,rvec xc[]) = 0;
  virtual void prep_data(int gnx,atom_id index[],rvec xcur[],rvec xprev[],
			 matrix box) = 0;
};


#endif	/* _corrutil_h */







