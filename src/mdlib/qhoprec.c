/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: poldata.c,v 1.20 2009/05/17 13:56:55 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
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
#include <string.h>
#include "smalloc.h"
#include "txtdump.h"
#include "futil.h"
#include "qhoprec.h"
#include "qhop_toputil.h"
#include "types/gmx_qhop_types.h"
#include "types/mdatom.h"
#include "types/topology.h"

t_qhoprec *mk_qhoprec(void)
{
  t_qhoprec *qr;

  snew(qr,1);

  return qr;
}  /* mk_qhoprec */

/* Sets the bqhopdonor[] and bqhopacceptor[] arrays in a t_mdatoms. */
static void qhop_atoms2md(t_mdatoms *md, const t_qhoprec *qr)
{
  int i, j;
  t_qhop_atom *a;
  qhop_db *db;

  db = qr->db;

  /* Should probably set the massT array too. */

  for (i=0; i < qr->nr_qhop_atoms; i++)
    {
      a = &(qr->qhop_atoms[i]);

      md->bqhopacceptor[a->atom_id] = (a->state & eQACC) != 0;
      md->bqhopdonor[a->atom_id]    = (a->state & eQDON) != 0;

      for (j=0; j < a->nr_protons; j++)
	{
	  if (db->H_map.H[db->H_map.atomid2H[a->protons[j]]] == 0)
	    {
	      md->massT[a->protons[j]] = 0;
	    }
	}
    }
}

/* Sets the interactions according to the hydrogen existence map.
 * This requires a finalized t_mdatoms. */
void finalize_qhoprec(t_qhoprec *qhoprec, gmx_localtop_t *top, 
		      t_mdatoms *md, t_commrec *cr)
{
  int i, j, nb, ft;
  qhop_db *db;
  t_qhop_residue *qres;
  /* #define DUMPTOP 1 */
#ifdef DUMPTOP 
  FILE *fp;
  if (NULL != (fp=ffopen("before.txt","w"))) 
    {
      pr_ltop(fp,0,"before",top,TRUE);
      fclose(fp);
    }
#endif
  db = qhoprec->db;

  for (i=0; i < qhoprec->nr_qhop_residues; i++)
    {
      qres = &(qhoprec->qhop_residues[i]);

      /* What flavour of the residue shall we set it to? */
      qres->subres = which_subRes(qhoprec, db, i);
      
      /* Alocate memory for the bonded index. */
      /* This should probably go to qhop_toputil.c */
      snew(qres->bindex.ilist_pos, F_NRE);
      snew(qres->bindex.indexed,   F_NRE);
      for (j=0; j < ebtsNR; j++)
	{
	  ft = db->rb.btype[j];
	  
	  if (ft >= 0)
	    {
	      nb = db->rtp[db->rb.qrt[qres->rtype].irtp].rb[j].nb;
	      qres->bindex.nr[ft] = nb;

	      snew((qres->bindex.ilist_pos[ft]), nb);

	      snew((qres->bindex.indexed[ft]), nb);
	    }
	}
      set_interactions(qhoprec, db, top, md, qres, cr);
    }
#ifdef DUMPTOP 
  if (NULL != (fp=ffopen("after.txt","w"))) 
    {
      pr_ltop(fp,0,"after",top,TRUE);
      fclose(fp);
    }
#endif

  qhop_atoms2md(md, qhoprec);
}
