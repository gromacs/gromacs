/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * GRoups of Organic Molecules in ACtion for Science
 */
#ifndef _G_HBOND_H_
#define _G_HBOND_H_


/* TODO list */
/*
 1. check generating search list 
 group 9 geeft vreemde fouten
 2. selected zowiezo alle hbonds meenemen
 */



typedef enum {
  INTER_GRP,INTRA_GRP,SELECTED,INSERT,COUNT,userNR
} t_mode;

typedef enum {
  INTER,INTRA,ALL,typeNR
} t_type;

typedef enum {
  NN0,NN3,NN4,NN5,eNNr
} t_helix_hb;


/* global variables */
extern matrix box;       
extern rvec   *x;        
extern real   rcut;      
extern real   rcut2;
extern real   alfcut;    
extern t_topology *top;  
extern real   this_time; 
extern int    nr_frames; 
extern int    this_frame;
extern t_mode mode;
extern int    nr_hbonds;
extern int    natoms;

#endif





