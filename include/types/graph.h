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

#include "fatal.h"

typedef struct {
  int      maxbond;     /* Max number of bonds per atom                 */
  int      nnodes;	/* The number of nodes				*/
  int      nbound;	/* The number of nodes with edges		*/
  int      start;	/* The first atom in this graph			*/
  int      end;		/* The last atom in this graph			*/
  int      *nedge;	/* For each node the number of edges		*/
  atom_id  **edge;	/* For each node, the actual edges (bidirect.)	*/
  int      *ishift;	/* The actual shift codes			*/
} t_graph;


#ifdef HALLO
#define SHIFT_INDEX(g,i) ((g)->ishift((i)-(g)->start))
#else
static int shift_index(t_graph *g,int i)
{
  if ((i < g->start) || (i > g->end))
    fatal_error(0,"Shift index out of range: i=%d, start=%d,end=%d",
		i,g->start,g->end);
  return g->ishift[i-g->start];
}
#define SHIFT_INDEX(g,i) shift_index(g,i)
#endif
