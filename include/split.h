/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _split_h
#define _split_h

static char *SRCID_split_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) split.h 1.20 12/16/92"
#endif /* HAVE_IDENT */

/*
 * Determine on which node a particle should reside and on which
 * node is also should be available. The distribution algorithm
 * should account for the actual ring architecture and how nodes
 * are numbered. The typedef t_splitd has two separate structures that
 * describe the distribution:
 *
 * The nodeinfo part describes which node containst which particles, 
 * while the nodeids part describes on which node(s) a particle can be 
 * found and what local particle number is assigned to it.
 *
 */

#include <stdio.h>
#include "typedefs.h"

typedef enum {SPLIT_NONE,SPLIT_SORTX,SPLIT_REDUCE,SPLIT_NR} t_splitalg;

typedef struct
{
  int hid;
  atom_id *nodeid;
} t_nodeids;

typedef struct
{
  int nr;		/* Length of the long list.                         */
  int *lst;		/* The actual list.                                 */
} t_nlist;

typedef struct
{
  t_nlist home;		/* List of home particles.                          */
} t_nodeinfo;

typedef struct
{
  int nnodes;		/* Number of nodes this splitinfo is for.      */
  t_nodeinfo *nodeinfo;	/* Home and available particles for each node. */
  int nnodeids;		/* Number of particles this splitinfo is for.       */
  t_nodeids *nodeids;	/* List of node id's for every particle,       */
  			/* entry[nodeid] gives the local atom id (NO_ATID if*/
			/* not available). Entry[MAXNODES] contains home    */
                        /* node's id.                                  */
} t_splitd;

extern void init_splitd(t_splitd *splitd,int nnodes,int nnodeids);
     /*
      * Initialises the splitd data structure for the specified number of
      * nodes (nnodes) and number of atoms (nnodeids).
      */
 
extern void make_splitd(t_splitalg algorithm,int nnodes,t_topology *top,
                        rvec *x,t_splitd *splitd,char *loadfile);
     /*
      * Initialises the splitd data structure for the specified number of
      * nodes (nnodes) and number of atoms (top) and fills it using
      * the specified algorithm (algorithm):
      *
      *    SPLIT_NONE   : Generate partial systems by dividing it into nnodes
      *                   consecutive, equal, parts without any intelligence.
      *    SPLIT_SORTX  : Like SPLIT_NONE but sort the coordinates before 
      *                   dividing the system into nnodes consecutive, equal, 
      *                   parts.
      *    SPLIT_REDUCE : Like SPLIT_NONE but minimise the bond lengths, i.e
      *                   invoke the reduce algorithm before dividing the 
      *                   system into nnodes consecutive, equal, parts.
      *
      * The topology (top) and the coordinates (x) are not modified. The 
      * calculations of bonded forces are assigned to the node with
      * the highest id that has one of the needed particles as home particle.
      */
                  
extern long wr_split(FILE *fp,t_splitd *splitd);
     /*
      * Writes the split descriptor (splitd) to the file specified by fp.
      */

extern long rd_split(FILE *fp,t_splitd *splitd);
     /*
      * Reads the split descriptor (splitd) from the file specified by fp.
      */

extern void rm_splitd(t_splitd *splitd);
     /*
      * Frees all allocated space for the splitd data structure.
      */

extern void pr_splitd(FILE *fp,int indent,char *title,t_splitd *splitd);
     /*
      * This routine prints out a (human) readable representation of 
      * the split descriptor to the file fp. Ident specifies the
      * number of spaces the text should be indented. Title is used
      * to print a header text.
      */

extern void split_topology(t_splitalg algorithm,int nnodes,t_topology *top,
                           rvec x[],char *loadfile);
     /*
      * Distributes the non-bonded forces defined in top over nnodes nodes
      * using the algoritm specified by algorithm. The distribution is made
      * by creating a split descriptor and then putting a bonded force on the 
      * highest home node number of the paricles involved.
      */
      

#endif	/* _split_h */
