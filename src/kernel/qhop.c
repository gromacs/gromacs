#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "smalloc.h"
#include "assert.h"
#include "physics.h"
#include "macros.h"
#include "vec.h"
#include "force.h"
#include "invblock.h"
#include "confio.h"
#include "names.h"
#include "network.h"
#include "pbc.h"
#include "ns.h"
#include "nrnb.h"
#include "bondf.h"
#include "mshift.h"
#include "txtdump.h"
#include "copyrite.h"
#include <stdio.h>
#include <string.h>
#include "gmx_fatal.h"
#include "typedefs.h"
#include <stdlib.h>
#include "mtop_util.h"
#include "qhop.h"
#include "random.h"
#include "gmx_random.h"
#include "constr.h"
#include "types/idef.h"
#include "hackblock.h"
#include "gmx_qhop_db.h"
#include "types/qhoprec.h"
#include "types/gmx_qhop_types.h"
#include "lambertw.h"

#define VERBOSE_QHOP
/* THIS IS JUST FOR QUICK AND DIRTY TEST WORK!!! */
/* typedef struct{ */
  
/*   real alpha, beta, gamma; */
/*   real  k_1, k_2, k_3, m_1, m_2, m_3; */
/*   real  s_A, t_A, v_A, s_B, s_C, t_C, v_C; */
/*   real  f, g, h; */
/*   real  p_1, q_1, q_2, q_3, r_1, r_2, r_3; */


/* } t_qhop_parameters; */

static const char *waterName = "HXO";

qhop_parameters *get_qhop_params(char *donor_name,char *acceptor_name, qhop_db_t db){
  /* parameters of the current donor acceptor combination are based on 
   * the residue name.
   */
  qhop_parameters
    *q;
  snew(q,1);

#ifdef VERBOSE_QHOP
  fprintf(stderr, "--- PARAMETERS FOR DONOR %s - ACCEPTOR %s ---\n",
	  donor_name, acceptor_name);
#endif /* VERBOSE_QHOP */
  
  if (db == NULL)
    gmx_fatal(FARGS, "Qhop database not initialized.");
  else
    if (qhop_db_get_parameters(db,donor_name,acceptor_name,q) != 1)
      gmx_fatal(FARGS, "Parameters not found in qhop database.");
#ifdef VERBOSE_QHOP
  qhop_db_print(q);
#endif /* VERBOSE_QHOP */
  return (q);
}

/* Make sure that the protonation state can be dealt with correctly */
static void qhop_check_validity(int nqhopatoms, unsigned int protonation_state)
{
    if (nqhopatoms < 0)
        gmx_fatal(FARGS, "Cannot accept a negative number of qhop atoms.");
    else if (nqhopatoms > 8*sizeof(protonation_state))
        gmx_fatal(FARGS, "%d bits are insufficient for %d possible protonation states.", 
                8*sizeof(protonation_state), nqhopatoms);
}


/* Translates a protonation state into an integer number */
static unsigned int qhop_get_protonation_state(int nqhopatoms, gmx_bool *bHyd)
{
    int i;
    unsigned int protonation_state=0;
    

    qhop_check_validity(nqhopatoms, protonation_state);
    
    for (i=0; i<nqhopatoms; i++)
        protonation_state |=  bHyd[i] << i; 
    
    
    return protonation_state;
}


/* OBSOLETE!
 * Reuse the name for similar purpose later. */

/* Translates an integer number into a protonation state */
/* static int qhop_set_protonation_state(int nqhopatoms, unsigned int protonation_state, gmx_bool *bHyd) */
/* { */
/*     int i, nchanges=0; */
/*     unsigned int new_state; */
    
    
/*     qhop_check_validity(nqhopatoms, protonation_state); */
    
/*     for (i=0; i<nqhopatoms; i++) */
/*     { */
/*         new_state = (protonation_state & (1<<i)); */
/*         if (bHyd[i]<<i != new_state) */
/*             nchanges++; */
/*         bHyd[i] = new_state; */
/*     } */
    
/*     return nchanges; */
/* } */

static void set_proton_presence(qhop_H_exist *Hext, int atomid, gmx_bool present)
{
  Hext->H[Hext->atomid2H[atomid]] == present ? (char)1: (char)0;
}


static gmx_bool get_proton_presence(qhop_H_exist *Hext, int atomid)
{
  return Hext->H[Hext->atomid2H[atomid]] == (char)1;
}


/* This function now makes get_qhop_residue obsolete. */
/* remove commrec? */
static int get_qhop_atoms(t_commrec *cr, 
			  gmx_mtop_t *mtop,
			  t_qhoprec *qhoprec,
			  t_inputrec *ir,
			  qhop_db *qdb){

  /* reads in the qhop donor and acceptor atoms, finds their protons
     (the ones that can move), identifies their residue numbers and
     stores the information in the qhop_atoms struct.
  */

  /* Uses the qhop_resblocks in qdb as a basis for finding acceptors and donors. */
  int rb, r, a, reac,
    **H, nH, Ha, *nHref, nHext;
  char *Hname;
  gmx_bool matchRB, match, same, exists;

  qhop_reactant *qreac[2];
  int nreac[2], react;
  gmx_groups_t *groups;
  gmx_mtop_atomloop_all_t aloop;
  t_atom    *atom;
  char *resname, *atomname;
  int
    q_atoms_nr=0,q_atoms_max=1000,
    q_residue_nr=0,q_residue_max=1000,
    i,j,jmax,resnr;
  t_qhop_atom
    *q_atoms;
  t_qhop_residue
    *q_residue;
  //  snew(resname,6);
  // snew(atomname,6);
  snew(q_atoms,q_atoms_max);
  snew(q_residue,q_residue_max);
    
/*   groups = &mtop->groups; */

/*   jmax = ir->opts.ngqhopdonors; */
  fprintf(stderr,
	  "ir->opts.ngqhopdonors = %d\nir->opts.ngqhopacceptors = %d\n",
	  ir->opts.ngqhopdonors,ir->opts.ngqhopacceptors);

  /*for(j=0;j<jmax;j++){ We don't loop over all atoms in a group anymore. We scan the resblocks instead.*/
  aloop = gmx_mtop_atomloop_all_init(mtop);
  while (gmx_mtop_atomloop_all_next(aloop,&i,&atom)) {
    //      fprintf(stderr, "don, i = %d\n",i);

    /*       if (ggrpnr(groups,egcqhopdonors ,i) == j) { */
    /* 	/\* we can complete the rest of the qhop_atoms struct now: */
    /* 	 *\/ */
    /* 	if(qhop_atoms_nr >= qhop_atoms_max){ */
    /* 	  qhop_atoms_max += 1000; */
    /* 	  srenew(qhop_atoms,qhop_atoms_max); */
    /* 	} */
    /* 	gmx_mtop_atomloop_all_names(aloop,&atomname,&resnr,&resname); */
    /* 	qhop_atoms[qhop_atoms_nr].atom_id = i; */
    /* 	qhop_atoms[qhop_atoms_nr].res_id = resnr; */
    /* 	qhop_atoms[qhop_atoms_nr].state = 1;/\* active *\/ */
    /* 	/\*qhop_atoms[qhop_atoms_nr].bdonor = TRUE;*\/ */
    /* 	strncpy(qhop_atoms[qhop_atoms_nr].resname,resname,3); */
    /* 	strncpy(qhop_atoms[qhop_atoms_nr].atomname,atomname,3); */
    /* 	if(!strncmp(resname,"SOL",3)||!strncmp(resname,"HOH",3)){ */
    /* 	  qhop_atoms[qhop_atoms_nr].bWater=TRUE; */
    /* 	} else { */
    /* 	  qhop_atoms[qhop_atoms_nr].bWater=FALSE; */
    /* 	} */
    /* 	qhop_atoms_nr++;	 */
    /*       } */

    /* Do the resname and atomname show up in the qhop_resblocks?
     * The loops below replace the loops above (and the donor loop further down). */

    if(q_atoms_nr >= q_atoms_max)
      {
	q_atoms_max += 1000;
	srenew(q_atoms,q_atoms_max);
      }

    if(q_residue_nr >= q_residue_max)
      {
	q_residue_max += 1000;
	srenew(q_residue,q_residue_max);
      }

    gmx_mtop_atomloop_all_names(aloop,&atomname,&resnr,&resname);
    matchRB = FALSE;
    match = FALSE;
    qhoprec->global_atom_to_qhop_atom[i] = NOTSET;

    for (rb=0; rb < qdb->rb.nrestypes && !matchRB; rb++)
      {
	if (strcmp(qdb->rb.restype[rb], resname) == 0)
	  {
	    /* We have a matching res(block) */
	    matchRB = TRUE;

	    for (r=0; r < qdb->rb.nres[rb] && !match; r++)
	      {
		/* Is the atom found among the donors/acceptors ?*/
		nreac[0] = qdb->rb.res[rb][r].na;
		nreac[1] = qdb->rb.res[rb][r].nd;
		qreac[0] = qdb->rb.res[rb][r].acc;
		qreac[1] = qdb->rb.res[rb][r].don;

		/* react==0 => acceptors, 1 => donors. This enables looping over both. */
		for (react=0;
		     !match && react<2;
		     react++)		  
		  {
	      
		    for (reac=0;
			 !match && reac<nreac[react];
			 reac++)
		      {
			for (a=0;
			     !match && a<qreac[react][reac].nname;
			     a++)
			  {
			    if (strcmp(qreac[react][reac].name[a], atomname) == 0)
			      {
				match = TRUE;
				q_atoms[q_atoms_nr].resname      = qdb->rb.restype[rb];
				q_atoms[q_atoms_nr].atomname     = qreac[react][reac].name[a];
				q_atoms[q_atoms_nr].res_id       = resnr;
				q_atoms[q_atoms_nr].atom_id      = i;
				q_atoms[q_atoms_nr].nr_protons   = NOTSET; /* To be decided */
				q_atoms[q_atoms_nr].protons      = NULL;
				q_atoms[q_atoms_nr].nr_acceptors = NOTSET;
				q_atoms[q_atoms_nr].acceptors    = NULL;
				qhoprec->global_atom_to_qhop_atom[i]   = q_atoms_nr;

				/* Short circuiting makes this conditional valid */
				if ((q_atoms_nr==0)
				    || (q_atoms[q_atoms_nr].res_id !=
					q_atoms[q_atoms_nr-1].res_id))
				  {
				    /* New residue. */
				    q_atoms[q_atoms_nr].qres_id                = q_residue_nr;
				    q_residue[q_residue_nr].rtype              = rb;
				    q_residue[q_residue_nr].res                = NOTSET;     /* To be decided. */
				    q_residue[q_residue_nr].atoms              = NULL;     /* To be set. */
				    q_residue[q_residue_nr].nr_atoms           = 0;     /* To be decided. */
				    q_residue[q_residue_nr].nr_titrating_sites = 0; /* To be decided. */
				    q_residue[q_residue_nr].res_nr             = NOTSET;  /* To be set */
				    q_residue[q_residue_nr].res_nr             = resnr;

				    q_residue_nr++;
				  }

				if (q_residue_nr > 0)
				  {
				    srenew(q_residue[q_residue_nr-1].titrating_sites,
					   q_residue[q_residue_nr-1].nr_titrating_sites + 1);

				    q_residue[q_residue_nr-1].titrating_sites[q_residue[q_residue_nr-1].nr_titrating_sites] = q_atoms_nr;
				    q_residue[q_residue_nr-1].nr_titrating_sites++;
				  }
				else
				  {
				    gmx_fatal(FARGS, "q_residue_nr = %i, which shouldn't be able to happen!\n", q_residue_nr);
				  }

				q_atoms_nr++;		    
			      }
			  }
		      }
		  }
	      }
	  }
      }
  }
  fprintf(stderr,"jmax = %d, nr qhop donors = %d\n",
	  /*jmax*/ -1, q_atoms_nr);
  /* now we have the donors stored, we now check out the acceptors
   */
/*   jmax = ir->opts.ngqhopacceptors; */
/*   fprintf(stderr,"jmax acc = %d\n", */
/* 	  jmax);  */
/*   for(j=0;j<jmax;j++){ */
/*     aloop = gmx_mtop_atomloop_all_init(mtop); */
/*     while (gmx_mtop_atomloop_all_next(aloop,&i,&atom)) { */
    
/*       if (ggrpnr(groups,egcqhopacceptors ,i) == j) { */
/* 	if(qhop_atoms_nr >= qhop_atoms_max){ */
/* 	  qhop_atoms_max += 1000; */
/* 	  srenew(qhop_atoms,qhop_atoms_max); */
/* 	} */
/* 	/\* we can complete the rest of the qhop_atoms struct now: */
/* 	 *\/ */
/* 	gmx_mtop_atomloop_all_names(aloop,&atomname,&resnr,&resname); */
/* 	qhop_atoms[qhop_atoms_nr].atom_id = i; */
/* 	qhop_atoms[qhop_atoms_nr].res_id = resnr; */
/* 	qhop_atoms[qhop_atoms_nr].state = 0;/\* inactive *\/ */
/* 	/\* qhop_atoms[qhop_atoms_nr].bdonor = FALSE;/\\* inactive *\\/ *\/ */
/* 	strncpy(qhop_atoms[qhop_atoms_nr].resname,resname,3); */
/* 	strncpy(qhop_atoms[qhop_atoms_nr].atomname,atomname,3); */
/* /\* 	if(!strncmp(resname,"SOL",3)||!strncmp(resname,"HOH",3)){ *\/ */
/* /\* 	  qhop_atoms[qhop_atoms_nr].bWater=TRUE; *\/ */
/* /\* 	} else { *\/ */
/* /\* 	  qhop_atoms[qhop_atoms_nr].bWater=FALSE; *\/ */
/* /\* 	} *\/ */
/* 	qhop_atoms_nr++; */
/*       } */
/*     } */
/*   }  */


  snew(H, q_residue_nr);
  snew(nHref, q_residue_nr);

  for (r=0; r < q_residue_nr; r++)
    {
      nHref[r] = 0;
    }
  
  /* Now find the atoms belonging to each res. */
  aloop = gmx_mtop_atomloop_all_init(mtop);

  while (gmx_mtop_atomloop_all_next(aloop, &i, &atom))
    {
      gmx_mtop_atomloop_all_names(aloop, &atomname, &resnr, &resname);

      for (r=0; r < q_residue_nr; r++)
	{
	  if (resnr == q_residue[r].res_nr)
	    {
	      srenew(q_residue[r].atoms, q_residue[r].nr_atoms + 1);

	      q_residue[r].atoms[q_residue[r].nr_atoms] = i;
	      q_residue[r].nr_atoms++;

	      /* Is it a hydrogen? */
	      if (atom->atomnumber == 1)
		{
		  srenew(H[r], nHref[r]+1);
		  nHref[r]++;

		  /* Does it exist? H = -1 if not.*/
		  H[r][nHref[r]] = i;
		}
	    }
	}
    }

  /* Figure out the protonation state and which version of the residue to use
   * by looking at the names of the existing/non-existing protons and match them
   * with the rtp-entries. */
  for (rb=0; rb < q_residue_nr; rb++)
    {
      match = FALSE;

      /* Loop over possible protonation states (residue subtypes)*/
      for (r=0;
	   r < qdb->rb.nres[q_residue[rb].rtype] && !match;
	   r++)
	{
	  /* Loop over possible hydrogens */
	  for (Ha=0; Ha < nHref[rb]; Ha++)
	    {
	      /* get hydrogen NAMES
	       * look for them in residue r (Need the rtp data)
	       * if it's a match, then we've found our restype.
	       * Maybe this needs to be adjusted for tautomerism later,
	       * but let's keep it for now.
	       *
	       * This is how it can be done later on:
	       * Examine the titrating sites and count the number of protons on them.
	       * This informations is available through the mtop, qhop_residue, qhop_atoms,
	       * H_map, rtp and rb.res. Or so I think... */
	      
	      /* Find the res in qdb->rtp[] */
	      for (i=0; i < qdb->nrtp; i++)
		{
		  if (strcmp(qdb->rtp[i].resname,
			     qdb->rb.res[rb][r].name) == 0)
		    {
		      break;
		    }
		}

	      Hname = q_atoms[H[rb][Ha]].atomname;
	      /* Does it exist? */
	      exists = get_proton_presence(&(qdb->H_map), H[rb][Ha]);

	      /* can we find it in the rtp entry? */
	      match = TRUE; /* Assume that we can, and see if that holds. */
	      nH = 0;
	      j = 0;

	      for (a=0; qdb->rtp[i].natom; a++)
		{
		  /* double check if it's really a hydrogen. Just to be sure. */
		  if (qdb->rtp[i].atom[a].atomnumber == 1)
		    { /* yes */
		      same = strcmp(*(qdb->rtp[i].atomname[a]), Hname) == 0;

		      if (same)
			{
			  break;
			}
		    }
		  else
		    {
		      gmx_fatal(FARGS, "Atom %d in res %s is not a hydrogen: %s",
				a, qdb->rtp[i].resname, Hname);
		    }
		}

	      if ((exists && !same) || (!exists && same))
		{
		  /* Oops. This H was found but doesn't exist, or the other way arond.
		   * This can't be the residue subtype we're looking for. Stop scanning this
		   * rtp entry and look at the next one.*/
		  match = FALSE;
		  break;
		}
	      /* If the previous conditional was false, then this may still be the
	       * residue subtype we're looking for, and the matchhypothesis remains TRUE.
	       * Keep scanning the rest of the hydrogens in the "superresidue" */
	    }
	}

      if (!match)
	{
	  gmx_fatal(FARGS, "Could not find a residue subtype that matches the given proton configuraion:\n"
		    " - %s %i", qdb->rb.restype[q_residue[rb].rtype], q_residue[rb].res_nr);
	}

      q_residue[rb].res = i;
    }

  /* Clean up */
  for (r=0; r<q_residue_nr; r++)
    {
      sfree(H[r]);
    }
  sfree(H);
  sfree(nHref);
  
  srenew(q_atoms, q_atoms_nr);     /* size it down, as it is likely oversized. */
  srenew(q_residue, q_residue_nr); /* size it down, as it is likely oversized. */
  qhoprec->qhop_atoms       = q_atoms;
  qhoprec->qhop_residues    = q_residue;
  qhoprec->nr_qhop_atoms    = q_atoms_nr;
  qhoprec->nr_qhop_residues = q_residue_nr;
  //  free(atomname);
  //  free(resname);
  return(q_atoms_nr);
} /* get_qhop_atoms */


/* Reads the info in the qhop_H_exist and finds the corresponding
   residue subtype in the qhop_db. */
static int H_exist_2_subRes(const gmx_mtop_t *top, const t_qhoprec *qr,
			    const qhop_db *db, const int resnr)
{
  int
    h, i, j, k, nH, *nHres,
    atomid, nres, natoms, qrtp;
  char **Hlist=NULL, *Hname;
  t_atom atom;
  t_qhop_residue *qres;
 /*  qhop_res *res = db->rb.res[qres->rtype]; */
  gmx_bool bNewH, bSameRes;

  qres = &(qr->qhop_residues[resnr]);
  nres = db->rb.nres[qres->rtype];

  /* Make a list of all hydrogens for this restype */
  qrtp = db->rb.rtp[qres->rtype]; /* The rtp data */
  natoms = db->rtp[qrtp].natom;

  for (i=0, nH=0; i<natoms; i++)
    {
      if (db->rtp[qrtp].atom[i].atomnumber == 1)
	{	  
	  srenew(Hlist, nH);
	  Hlist[nH] = *(db->rtp[qrtp].atomname[i]);
	  nH++;
	}
    }

  /* Loop over the residue subtypes and see which one matches */
  for (i=0; i<nres; i++)
    {
      /* scan the rtp */
      for (j=0; j < db->nrtp; j++)
	{
	  if (gmx_strcasecmp(db->rtp[j].resname,
			     db->rb.res[qres->rtype][i]) == 0)
	    { /* Mathcing resnames! Now look at the protons. */
	      bSameRes = TRUE;
	    }
	}
    }
}

/* Qr    = qhop residue
 * to    = to what residue subtype  */
static void set_interactions(t_qhoprec *qhoprec, qhop_db *qdb, int Qr, int to)
{
  int i, nri, bt, b, Gr, RB, R;
  t_qhop_residue *qres;
  qhop_res *reac;
  t_restp *res;
  
  /* point pointers at the right residue */
  qres = &(qhoprec->qhop_residues[Qr]);
  Gr   = qres->res_nr; /* global res id. */
  RB   = qres->rtype;
  R = qres->res;
  reac = &qdb->rb.res[RB][R];
  
  /* find the res in the rtp */
  for (i=0; i < qdb->nrtp; i++)
    {
      if (strcmp(reac->name,
		 qdb->rtp[i].resname) == 0)
	{
	  res = &(qdb->rtp[i]);
	  break;
	}
    }
  
  /* --------- Bonded interactions --------- */

  /* loop over bonded types */
  for (bt = 0; bt<ebtsNR; bt++)
    {
      /* loop over interactions of type bt */
      for (b=0; b<res->rb[bt].nb; b++)
	{
	  
	}
    }


    /* Non-bonded */
}

static int int_comp(const void *a, const void *b){
  return (*(int *)a) - (*(int *)b);
} 

static int qhop_atom_comp(const void *a, const void *b){
  
  return (int)(((t_qhop_atom *)a)->atom_id)-(int)(((t_qhop_atom *)b)->atom_id);
  
}

/* OBSOLETE! */
/* remove commrec? */
/* static void get_qhop_residue(t_commrec *cr, gmx_mtop_t *mtop,  */
/* 			     t_qhop_atom    *qhop_atom, */
/* 			     t_qhop_residue *qhop_residue){ */
  
/*   /\* from the qhopatoms struct, we now retreive the corresponding */
/*      residues. Of these residues we store the charges for each of the */
/*      states, the residue can be in. The latter are picked up from a */
/*      file somewhere. */
/*   *\/   */

/*  /\* We'll do this in an entirely new way. We no longer have a t_qhop_residue per t_qhop_atom, */
/*   * but one per qhopable residue. As some residues may have several titrating sites the numbers may differ. *\/ */

/* /\*   gmx_mtop_atomloop_all_t aloop; *\/ */
/* /\*   t_atom *\/ */
/* /\*     *atom; *\/ */
/*   int */
/*     maxatom = 100, nr_atoms=0,j,resnr; */
/*   char *resname, *atomname; */
  //  snew(resname,6);
  //snew(atomname,6);
/*   snew(qhop_residue->atoms,maxatom); */
/*   qhop_residue->res_nr    = qhop_atom->res_id; */
/*   qhop_residue->atoms[0] = qhop_atom->atom_id; /\* we quicksort later... *\/ */

/*   aloop = gmx_mtop_atomloop_all_init(mtop); */
/*   while (gmx_mtop_atomloop_all_next(aloop,&j,&atom)) { */
/*     gmx_mtop_atomloop_all_names(aloop,&atomname,&resnr,&resname); */
/*     if(resnr == qhop_atom->res_id){ */
/*       if(nr_atoms >= maxatom){ */
/* 	maxatom += 100; */
/* 	srenew(qhop_residue->atoms,maxatom); */
/*       } */
/*       qhop_residue->atoms[nr_atoms++] = j; */
/*     } */
/*   }     */
/*   qhop_residue->nr_atoms = nr_atoms; */
/*   srenew(qhop_residue->atoms,nr_atoms); */
/*   qsort(qhop_residue->atoms,nr_atoms, */
/* 	(size_t)sizeof(qhop_residue->atoms[0]),int_comp); */
  //  free(atomname);
  //free(resname);
/* } /\* get_qhop_residue *\/ */

/* Call this one to decide wether e.g. LYS or LYSH should be used,
 * based upon the hydrogen existence function.
 * This should only be required at t=0, because at other times we know
 * from and to what residue subtype a transfer is going. */

/* OBSOLETE! */
/* static void determine_protonation_state(t_qhoprec *qhoprec, qhop_db *qdb, gmx_mtop_t *top,  int resnr) */
/* { */
/*   int a, r=0, nres=0, i=0, nH, nHref; */
/*   int *H; */
/*   gmx_bool match; */
/*   t_qhop_residue *res=NULL; */
/*   qhop_res *qres; */
/*   gmx_mtop_atomloop_all_t aloop; */
/*   t_atom    *atom; */

/*   for (r=0; r<qhoprec->nr_qhop_residues; r++) */
/*     if (resnr = qhoprec->qhop_residues[r].res_nr) */
/*       { */
/* 	res = &(qhoprec->qhop_residues[r]); */
/* 	break; */
/*       } */
  
/*   if (res == NULL) */
/*     gmx_fatal(FARGS, "Couldn't find residue with global resnr %i in the qhoprec."); */

/*   /\* extract the hydrogens in res to a short list *\/ */

/*   snew(H, 10); */
/*   aloop = gmx_mtop_atomloop_all_init(mtop); */
/*   while (gmx_mtop_atomloop_all_next(aloop,&i,&atom)) */
/*     { */
/*       gmx_mtop_atomloop_all_names(aloop,&atomname,&resnr,&resname); */
/*       if (nH % 10 == 0) */
/* 	{ */
/* 	  nH+=10; */
/* 	  srenew(H, nH); */
/* 	} */
      
/*       for (a=0; a<res->nr_atoms; a++) */
/* 	{ */
	  
/* 	} */
/*     } */
/*   /\* Which restypes are available? *\/ */
/*   nres = qdb->rb.nres[res->rtype]; /\* How many to choose from? *\/ */
/*   match = FALSE; */
/*   for (r=0; r<nres && !match; r++) */
/*     { */
/*       /\* compare the names of the existing and non-existing hydrogens */
/*        * with those in the residue subtypes *\/ */
      
/*     } */
  

/*   /\* specify residue type from the H-existence array. *\/ */
/*   sfree(H); */
/* } */

static void get_protons(t_commrec *cr, gmx_mtop_t *mtop, 
			t_qhop_atom    *qhop_atom,
			t_qhop_residue *qhop_residue,
			rvec *x, t_pbc pbc){
  /* searching through the residue atoms to find the protons. By using
     a distance criteriumm we avoid the graph (don't know if there is
     such thing with v-sites, which we eventually want to use for
     qhop).
  */
  int
    i,nr_protons,max_protons;
  t_atom
    *atom;
  rvec
    vec;
  real
    dist;

  nr_protons = 0;
  max_protons = 10;

  snew(qhop_atom->protons,max_protons);

  for(i=0; i < qhop_residue->nr_atoms; i++)
    {
      gmx_mtop_atomnr_to_atom(mtop,qhop_residue->atoms[i],&atom);

      if (atom->atomnumber <= 1)
	{
	  /* we need  v-sites as well */
	  pbc_dx(&pbc,
		 x[qhop_atom->atom_id],
		 x[qhop_residue->atoms[i]],
		 vec);

	  dist = sqrt(iprod(vec,vec));

	  if(dist < BOND )
	    {
	      if(nr_protons >= max_protons)
		{
		  max_protons += 10;
		  srenew(qhop_atom->protons,max_protons);
		}
	      qhop_atom->protons[nr_protons++] = qhop_residue->atoms[i];
	    }
	}
    }

  qhop_atom->nr_protons = nr_protons;
  srenew(qhop_atom->protons, nr_protons);
} /* get_protons */ 

/* This must go! */
static void create_res_links(t_commrec *cr, t_qhoprec *qhoprec){
  /* if a residue has more than one titratable atoms (stored in the
     qhop_atom, we need to connect these residues, such that we can
     modifiy the state and charges simultaneously if one of the atoms
     changes protonation state. 
  */
  int
    i,j;
  
  /* we also quicksort the qhop atoms, which is important later on
     when we want to find out what chargeset we are supposed to use
     for the residue */
  
  qsort(qhoprec->qhop_atoms,qhoprec->nr_qhop_atoms,
	(size_t)sizeof(qhoprec->qhop_atoms[0]),qhop_atom_comp);
  

  for(i=0;i<qhoprec->nr_qhop_atoms;i++){
    qhoprec->global_atom_to_qhop_atom[qhoprec->qhop_atoms[i].atom_id]=i;
/*     qhoprec->qhop_atoms[i].nr_links=0; */
/*     snew(qhoprec->qhop_atoms[i].links,qhoprec->nr_qhop_atoms); */
  }
  /* double loop, there might be more efficient ways, but I do not
     care... ;-)
   */
  for(i=0;i<qhoprec->nr_qhop_atoms;i++){
    for(j=i+1;j<qhoprec->nr_qhop_atoms;j++){
      if (qhoprec->qhop_atoms[i].res_id==qhoprec->qhop_atoms[j].res_id){
/* 	qhoprec->qhop_atoms[i].links[qhoprec->qhop_atoms[i].nr_links++] */
/* 	  = j; */
/* 	qhoprec->qhop_atoms[j].links[qhoprec->qhop_atoms[j].nr_links++] */
/* 	  = i; */
      }
    }
  }  
  /* clean up some memory
   */
  for(i=0;i<qhoprec->nr_qhop_atoms;i++){
/*     srenew(qhoprec->qhop_atoms[i].links,qhoprec->qhop_atoms[i].nr_links); */
/*     qhoprec->qhop_residues[i].nr_titrating_sites =  */
/*       qhoprec->qhop_atoms[i].nr_links+1; */
  }
} /* create_res_links */

static void get_chargesets(t_commrec *cr, t_qhoprec *qhoprec){

  int 
    i,j,k,l,states;
  FILE
    *in;
  char
    filename[4],file[8],line[300],atomname[3];

  for (i=0; i < qhoprec->nr_qhop_atoms; i++)
    {
      states = 1;
/*     for(j=0;j<qhoprec->qhop_atoms[i].nr_links+1;j++){ */
/*       states*=2; */
/*     } */
/*     qhoprec->qhop_residues[i].max_state = states; */
/*   for(i=0;i<qhoprec->nr_qhop_atoms;i+=(qhoprec->qhop_atoms[i].nr_links+1)){ */
    /* each unique residue is visited once in this loop
     */
    /* open db file, GG style 
     */

    
    strncpy(filename, qhoprec->qhop_atoms[i].resname, 3);
    sprintf(file, "%3s.dat",qhoprec->qhop_atoms[i].resname);
    //  fprintf(stderr,"trying to open %s for reading\n",file);

    in = fopen (file,"r");
    
    /* allocate the chargesets, one for each protonation state 
     */ 
/*     for(j=0;j<qhoprec->qhop_atoms[i].nr_links+1;j++){ */
/*       snew(qhoprec->qhop_residues[i+j].charge_set, */
/* 	   qhoprec->qhop_residues[i+j].max_state); */
/*       for(k=0;k<qhoprec->qhop_residues[i+j].max_state;k++){ */
/*       	snew(qhoprec->qhop_residues[i+j].charge_set[k], */
/* 	     qhoprec->qhop_residues[i+j].nr_atoms); */
/* 	/\* read in the chargesets *\/ */
/* 	for(l=0;l<qhoprec->qhop_residues[i+j].nr_atoms;l++){ */
/* 	  /\* copy the charges from the qhop database SOMEHOW *\/ */
/* 	  fgets(line,200,in); */
/* 	  sscanf(line,"%s%f",atomname, */
/* 		 &qhoprec->qhop_residues[i+j].charge_set[k][l]); */
/* 	} */
/* 	fgets(line,200,in);/\* read in a white line after every charges */
/* 			      section *\/ */
/*       } */
/*     } */
    fclose(in);
  }
} /* get_chargesets */


static void get_protonation_state(t_commrec *cr, t_qhoprec *qhoprec){
  /* Here we can find out the state of the residue: We make use
     of the order to decide what the state is.... We use a the bit
     definition to decide on the protonation state.  In case of three
     sites, there are 8 possibilites, and therefore 8 chargesets to be
     retreived from file.

     0 1 2 3 4 5 6 7 
     ---------------
     0 1 0 1 0 1 0 1
     0 0 1 1 0 0 1 1
     0 0 0 0 1 1 1 1
     ---------------
     
     Thus, if the second titrating atom (in global atom numbering) is
     protonated and the third one too, we take the seventh (6, C
     counting) row of charges from the external qhop data base.
     Anyway, something like this should work I hope. We need to make
     sure there is a good manual chapter explaining this if someone
     wants to add a new residue.
  */


  int 
    i,j,k,states;
  int
    nr_protons;
  gmx_bool
    *protonated;
  
  /* now we need to decide the protonation state of the residue and
     retreive the charges. The qhop_atoms are sorted already.
  */
/*   for(i=0;i<qhoprec->nr_qhop_atoms;i+=(qhoprec->qhop_atoms[i].nr_links+1)){ */
/*     snew(protonated,qhoprec->qhop_atoms[i].nr_links+1); */
/*     nr_protons = 0; */
/*     /\* for(j=0;j<qhoprec->qhop_atoms[i].nr_links+1;j++){ *\/ */
/* /\*       protonated[j]=qhoprec->qhop_atoms[i+j].bdonor; *\/ */
/* /\*       nr_protons+=qhoprec->qhop_atoms[i+j].bdonor; *\/ */
/* /\*     } *\/ */
/*     for(j=0;j<qhoprec->qhop_atoms[i].nr_links+1;j++){ */
/*       srenew(qhoprec->qhop_residues[j+i].protonated, */
/* 	   qhoprec->qhop_atoms[i].nr_links+1); */
/*       for(k=0;k<qhoprec->qhop_atoms[i].nr_links+1;k++){ */
/* 	qhoprec->qhop_residues[j+i].protonated[k]=protonated[k]; */
/*       } */
/*       /\* compute the integer corresponding to the proton byte, or pryte  */
/*        *\/ */
/*       qhoprec->qhop_residues[i+j].pryte =  */
/* 	qhop_get_protonation_state(qhoprec->qhop_atoms[i].nr_links+1,  */
/* 				       protonated); */
/*     } */
/*     /\* with the pryte, we can decide what chargeset to use during the */
/*        simulation: charge[pryte][.] */
/*      *\/  */
/*     free(protonated); */
/*   } */
} /* get_protonation_state */
 
t_qhoprec *mk_qhoprec(void)
{
  t_qhoprec *qr;

  snew(qr,1);
  memset((void *)qr, 0, sizeof(t_qhoprec)); /* Just in case */
  return (qr);
}  /* mk_qhoprec */

static void set_charges(t_commrec *cr, t_qhoprec *qhoprec, t_mdatoms *md)
{
  /* check and correct the charges of the donors etc, using the pryte
   */
  t_qhop_residue
    *res;
  int
    i,j;
  
  for(i=0; i < qhoprec->nr_qhop_residues; i++)
    {
      res = &qhoprec->qhop_residues[i];
      for(j=0; j < res->nr_atoms; j++){
	/*       md->chargeA[res->atoms[j]] = res->charge_set[res->pryte][j]; */
      }
    }
} /* set_charges */

int init_qhop(t_commrec *cr, gmx_mtop_t *mtop, t_inputrec *ir, 
	      t_forcerec *fr, rvec *x,matrix box,t_mdatoms *md,
	      qhop_db_t *db){
  /* initialize qhoprec, picks out the atoms that are
     donors/acceptors, creates a bqhop donor array in mdatoms to be
     used by nbsearch, completes the qhop residues array and reads in
     the parameters and charges of the residues. THere is a seperate
     function to handle the io with external files.
  */
  int 
    nr_qhop_atoms, nr_qhop_residues, i, target;
  t_qhoprec
    *qhoprec; 
  t_pbc
    pbc;
  
  qhoprec = NULL;
  nr_qhop_atoms    = 0;
  nr_qhop_residues = 0;

  if (ir->qhopmode == etQhopNONE)
    {
      return 0;
    }

  if ((*db = qhop_db_read("ffoplsaa", mtop, md)) == NULL)
    {
      gmx_fatal(FARGS,"Can not read qhop database information");
    }

  set_pbc_dd(&pbc,fr->ePBC,DOMAINDECOMP(cr) ? cr->dd : NULL, FALSE, box);
  
  qhoprec = fr->qhoprec; /* the mk_qhoprec is called in init_forcerec.
			  * No, it's done from md.c.
			  */
  qhoprec->qhopfreq = ir->qhopfreq;
  qhoprec->qhopmode = ir->qhopmode;

  snew(qhoprec->global_atom_to_qhop_atom, mtop->natoms);
  nr_qhop_atoms = get_qhop_atoms(cr, mtop, qhoprec ,ir, *db);
  /* now trace down the atoms belonging to these atoms to create an
     array with nr_qhop_atoms qhop_residue elements. We also set a
     corresponding residue ID, so that later we can change the states
     of atoms on the same residue instantly.
  */
  /* qhoprec->nr_qhop_residues = nr_qhop_atoms; */
/*   qhoprec->nr_qhop_atoms    = nr_qhop_atoms; */
/*   srenew(qhoprec->qhop_atoms, nr_qhop_atoms); */
/*   snew(qhoprec->qhop_residues,nr_qhop_atoms); */
  /* for(i=0;i<nr_qhop_atoms;i++){ */
/*     get_qhop_residue(cr,mtop,&qhoprec->qhop_atoms[i], */
/* 		     &qhoprec->qhop_residues[i]); */
/*     /\* now we have to find the protons that can hop. The in-active */
/*        protons that reside on the donor side are also included. For */
/*        water we need to be careful, as we will be using David's */
/*        settled/v-site water/hydronium hybrid model */
/*     *\/ */
/*     get_protons(cr,mtop,&(qhoprec->qhop_atoms[i]), */
/* 		&(qhoprec->qhop_residues[i]),x,pbc); */
/*   } */
/*   create_res_links(cr, qhoprec); */
/*   fprintf(stderr,"links created\n"); */
/*   get_chargesets(cr, qhoprec); */
/*   get_protonation_state(cr, qhoprec); */
/*   set_charges(cr, qhoprec, md); */
  for (i=0; i < qhoprec->nr_qhop_residues; i++)
    {
      /* What flavour of the residue shall we set it to? */
      target = H_exist_2_subRes(mtop, qhoprec, *db, i);
      set_interactions(qhoprec, *db, i, target);
    }

  qhoprec->hop = NULL;
  
  return (nr_qhop_atoms);
} /* init_hop */

static t_hop *find_acceptors(t_commrec *cr, t_forcerec *fr, rvec *x, t_pbc pbc,
			     int *nr, t_mdatoms *md){
  /* using the qhop nblists with i particles the donor atoms and
     j-particles the acceptor atoms, we search which of the acceptor
     atoms are available to accept a proton from the donor atoms. We
     make use of simple geometric criteria, such as distance and DHA
     angles. The acceptor atoms are stored in acceptor arrays of each
     of the donor atoms.
  */
  int
    i,j,jmax=100,k,acc,don;
  t_nblist
    qhoplist;
  t_qhoprec
    *qhoprec;
  t_qhop_atom
    donor,acceptor;
  t_hop
    *hop;
  rvec
    dx,veca,vecb,vecc;
  int
    nr_hops=0,max_hops=10, closest;
  gmx_bool 
    found_proton;
  real 
    ang, r_close;
  qhoprec = fr->qhoprec;
  qhoplist = fr->qhopnblist;
  snew(hop, max_hops);
  for(i=0; i < qhoplist.nri; i++)
    {
      /* every qhop atom should be a unique i
	 particle I think  */

      if(!md->bqhopdonor[qhoplist.iinr[i]])
	{
	  continue;
	}

      don      = qhoprec->global_atom_to_qhop_atom[qhoplist.iinr[i]];
      donor    = qhoprec->qhop_atoms[don];

      for(j = qhoplist.jindex[i];
	  j < qhoplist.jindex[i+1];
	  j++)
	{
	  if(!md->bqhopacceptor[qhoplist.jjnr[j]])
	    {
	      continue;
	    }

	  acc =  qhoprec->global_atom_to_qhop_atom[qhoplist.jjnr[j]];
	  acceptor = qhoprec->qhop_atoms[acc];

	  /* check whether the acceptor can take one of the donor's protons
	   */
	  pbc_dx(&pbc, x[acceptor.atom_id], x[donor.atom_id], dx);
	  if (norm(dx) <= DOO)
	    {
	      /* find the proton that can go. I assume that there can only be one!
	       */
	      found_proton = FALSE;
	      r_close = 100000;
	      for (k=0; k < donor.nr_protons; k++)
		{
		  /* TRICKY.... Only if proton carries a charge we allow the hop */
		  if( md->chargeA[donor.protons[k]] > 0.000001 )
		    {
		      pbc_dx(&pbc,x[donor.protons[k]],x[acceptor.atom_id],dx);
		      /* if(norm(dx)< HBOND) */
		      {
			found_proton = TRUE;
			if (norm(dx) < r_close)
			  {
			    r_close =  norm(dx);
			    closest = k;
			  }
			/* break; */
		      }
		    }
		}

	      k = closest;

	      if(found_proton)
		{
		  /* compute DHA angels */
		  pbc_dx(&pbc, x[donor.atom_id],    x[donor.protons[k]], veca);
		  pbc_dx(&pbc, x[donor.protons[k]], x[acceptor.atom_id], vecb);
		  pbc_dx(&pbc, x[acceptor.atom_id], x[donor.atom_id],    vecc);
		
		  ang = acos( (norm2(vecb) + norm2(veca) - norm2(vecc) )/
			      (2 * norm(veca) * norm(vecb)));

		  if (ang*180.0/(M_PI) >= (HOPANG))
		    {
		      /* the geometric conditions are right for a hop to take place
		       */
		      if(nr_hops>=max_hops)
			{
			  max_hops += 10;
			  srenew(hop, max_hops);
			}
		      hop[nr_hops].donor_id    = don;
		      hop[nr_hops].acceptor_id = acc;
		      hop[nr_hops].proton_id   = k;
		      hop[nr_hops].rda         = norm(vecc);
		      hop[nr_hops].prob        = 0;
		      nr_hops++;
		    }
		}
	    }
	}
    }

  /* clean up some memory */

  srenew(hop,nr_hops);
  *nr = nr_hops;
  return (hop);
} /* find_acceptors */

static void rotate_water(rvec *x,t_qhop_atom *res,real angle,t_mdatoms *md){
  int
    i,j;
  matrix
    r;
  rvec
    axis,xnew,temp;
  real
    n,theta;
  clear_rvec(axis);
  clear_rvec(xnew);
  clear_rvec(temp);

  /* find rotation axis */
  for(j=2; j < res->nr_protons; j++)
    {
      for(i=0; i<DIM; i++)
	{
	  axis[i] += x[res->protons[j]][i] / (3.0);
	}
    }

  for(i=0; i<DIM; i++)
    {
      axis[i] -= x[res->atom_id][i];
    }

  n=norm(axis);

  for(i=0; i<DIM; i++)
    {
      axis[i] = axis[i]/n;
    }

  theta = angle*M_PI/180.0;
  /* construct rotation matrix. I hope there is no mistake here.... 
   */
  r[0][0] = axis[0]*axis[0]+(1.0-axis[0]*axis[0])*cos(theta);
  r[0][1] = axis[0]*axis[1]*(1.0-cos(theta))-axis[2]*sin(theta);
  r[0][2] = axis[0]*axis[2]*(1.0-cos(theta))+axis[1]*sin(theta);
  
  r[1][0] = axis[0]*axis[1]*(1.0-cos(theta))+axis[2]*sin(theta);
  r[1][1] = axis[1]*axis[1]+(1.0-axis[1]*axis[1])*cos(theta);
  r[1][2] = axis[1]*axis[2]*(1.0-cos(theta))-axis[0]*sin(theta);
  
  r[2][0] = axis[0]*axis[2]*(1.0-cos(theta))-axis[1]*sin(theta);
  r[2][1] = axis[1]*axis[2]*(1.0-cos(theta))+axis[0]*sin(theta);
  r[2][2] = axis[2]*axis[2]+(1.0-axis[2]*axis[2])*cos(theta);

  for(j=0; j < res->nr_protons; j++)
    {
      for(i=0; i<DIM; i++)
	{
	  temp[i] = x[res->protons[j]][i] - x[res->atom_id][i];
	}

      mvmul(r,temp,xnew);     

      for(i=0; i<DIM; i++)
	{
	  x[res->protons[j]][i] = xnew[i] + x[res->atom_id][i];
	} 
    }
} /* rotate_water */

static void qhop_swap_bondeds(t_qhop_residue* swapres, qhop_res* prod)
{
  /* placeholder. */
  int i;
  i=0;
  
}

static void qhop_swap_nonbondeds(t_qhop_residue* swapres, qhop_res* prod)
{
  /* placeholder. */
  int i;
  i=0;
  
}


/* Supreseeded by flip-routines in change_protonation() */
/* static void invert_hydronium(){ */

/*   /\* if the accepting water would actually result in a inverted */
/*      hydronium, we need to re-invert it */
/*    *\/ */
/* } */

/* Arguments like qhop_(de)protonate(),
 * with the addition of bDonor which,
 * if TRUE, makes a deprotonation and,
 * if FALSE, a protonation. */
static void qhop_titrate(qhop_db *db, t_qhoprec *qr, t_qhop_atom *qatom, gmx_bool bWater, gmx_bool bSwapBondeds, gmx_bool bDonor)
{
  int
    i, rt, r, reactant, t;
  qhop_res
    *product_res;

  product_res = NULL;
  
  rt = qr->qhop_residues[qatom->qres_id].rtype;
  r = qr->qhop_residues[qatom->qres_id].res;

  /* Find out what the product is. */

  /* Loop over acceptors or donors */
  for (reactant = 0;
       (reactant < bDonor ? (db->rb.res[rt][r].nd) : (db->rb.res[rt][r].na))
	 && product_res==NULL;
       reactant++)
    {
      /* Loop over tautomers */
      for (t = 0;
	   t < bDonor ? (db->rb.res[rt][r].don[reactant].nname) : (db->rb.res[rt][r].acc[reactant].nname);
	   t++)
	{
	  /* Is this the reactant? */
	  if (!strcasecmp(db->rb.res[rt][r].acc[reactant].name[t],
			  qatom->atomname));
	    {
	      /* This is the reactant! */
	      product_res = db->rb.res[rt][r].acc[reactant].productdata;
	      break;
	    }
	}
    }

  if (product_res == NULL)
    {
      gmx_fatal(FARGS, "Could not find the product for this reaction.");
    }

  /* Change interactions */
  qhop_swap_bondeds   (&(qr->qhop_residues[qatom->qres_id]), product_res);
  qhop_swap_nonbondeds(&(qr->qhop_residues[qatom->qres_id]), product_res);

}


/* Wrappers for qhop_titrate() */
extern void qhop_protonate(qhop_db *db, t_qhoprec *qr, t_qhop_atom *qatom,
			   gmx_bool bWater, gmx_bool bSwapBondeds)
{
  qhop_titrate(db, qr, qatom, bWater, bSwapBondeds, FALSE);
}

extern void qhop_deprotonate(qhop_db *db, t_qhoprec *qr, t_qhop_atom *qatom,
			     gmx_bool bWater, gmx_bool bSwapBondeds)
{
  qhop_titrate(db, qr, qatom, bWater, bSwapBondeds, TRUE);
}


static gmx_bool change_protonation(t_commrec *cr, t_qhoprec *qhoprec, 
				   t_mdatoms *md, t_hop *hop, rvec *x,
				   gmx_bool bUndo,gmx_mtop_t *mtop,
				   rvec x_old, t_pbc *pbc,
				   qhop_db *db){
  /* alters the topology in order to evaluate the new energy. In case
     of hydronium donating a proton, we might have to change the
     position of the resulting HW1, and HW2 in order to create a
     correct water confiuration after the hop.  This is easiest done
     by a reflection in the plane of Hw1, OW, HW2 in the line
     OW-HW. We do not need to put this all back, if the proposal is
     not accepted. Water-Hydronium hops will be topology
     swaps. hydronium-residue hops will lead to annihilation of
     hydronium, and its index. residue-water hops lead to creation of
     a new hydronium.
  */

  /* Note that we still have special treatment of water and its symmetry.
   * This should probably be generalized, since similar things happen
   * with carboxylic acids, ammonia, etc. */

  int
    i;
  t_qhop_residue
    *donor,*acceptor;
  t_qhop_atom
    *donor_atom,*acceptor_atom;
  gmx_bool
    bAccWater, bDonWater;
  double
    ang;
  rvec
    x_tmp,
    w_normal,
    w_don,
    OH1,
    OH2,
    w_third;
  real
    planedist, d;
  const int
    /* 0 and 1 are ordinary water protons, 3 is the "out of plane" one
     * in hydronium, while 3 and 4 are in the water plane */
    out_of_plane_proton = 2;

  
  donor_atom    = &qhoprec->qhop_atoms[hop->donor_id];
  acceptor_atom = &qhoprec->qhop_atoms[hop->acceptor_id];
  donor         = &qhoprec->qhop_residues[hop->donor_id];
  acceptor      = &qhoprec->qhop_residues[hop->acceptor_id];

  bAccWater = db->rb.bWater[acceptor->rtype];
  bDonWater = db->rb.bWater[donor->rtype];

  /* in case of our special hydronium/waters, where the protons are
     simply a dummy charge, we might have to alter the positions of
     the real atoms.
  */
  
  /* check if we need to rotation */ 

  if(bUndo)
    {
      ang = -120;
    }
  else
    {
      ang = 120;
    }

  if(bDonWater){
    /* assumuming the user uses the correct water topology..... */

    switch (hop->proton_id)
      {
      case 3:
	rotate_water(x,donor_atom,-ang,md);
	break;
	
      case 4:
	rotate_water(x,donor_atom, ang,md);
	break;
	
      default:
	break;
      }
  }


  /* Check if the protonends up on the wrong side of the acceptor */
  if(bAccWater)
    {
      if (!bUndo)
	/* Do we need to flip? */
	{

	  pbc_dx(pbc,
		 x[donor_atom->atom_id],
		 x[acceptor_atom->atom_id],
		 w_don);
	  pbc_dx(pbc,
		 x[acceptor_atom->protons[0]],
		 x[acceptor_atom->atom_id],
		 OH1);
	  pbc_dx(pbc,
		 x[acceptor_atom->protons[1]],
		 x[acceptor_atom->atom_id],
		 OH2);
	  pbc_dx(pbc,
		 x[acceptor_atom->protons[out_of_plane_proton]],
		 x[acceptor_atom->atom_id],
		 w_third);

	  cprod(OH1, OH2, w_normal);
	  
	  hop->bFlip = (iprod(w_don, w_normal) * iprod(w_third, w_normal) < 0);
	  
	  /* So, if <AD, (DH1 x DH2)> is negative,
	     then the vsite construction needs flipping. */

	  /* ERRATA: above line should be correct, but the c-parameter
	   * is *negative* in the water/H3O model I'm using.
	   * Therefore I check if the projection of the out-of-plane-proton
	   * on the water normal points in the same
	   * direction as the projection of the AD-vector.
	   * <AD, w_normal> * <w_third, w_normal> > 0  <==>  no flipping needed.*/
	}
      
    /* If bUndo, use the existing value for hop->bFlip */
      
      if (hop->bFlip)
	{
	  /* Flip the water molecule */
	  
	  /* =============== Strategy #1 ================
	   * Swap hydrogen coordinates to have the third
	   * hydrogen appear on the other side.
	   * Problem is that velocities are determined
	   * from (x(t)-x(t-1)) / delta_t*, so velocities
	   * for those protons will be wrong for that
	   * timestep. This can be compensated/circumvented.
	   */
	
	  copy_rvec(x[acceptor_atom->protons[0]], x_tmp);
	  copy_rvec(x[acceptor_atom->protons[1]], x[acceptor_atom->protons[0]]);
	  copy_rvec(x_tmp, x[acceptor_atom->protons[1]]);

	  /* Update vsite position by calling construct_vsite() or by mirroring in the water plane */
	  /* Try manually, since we have a vector normal to the water plane half of the time.
	   * It will be done automatically in all concurrent time steps. */
	

	  if (bUndo)
	    {
	      /* We don't have the water normal, nor dowe have the other
	       * vectors that we need. Let's calculate them. */

	      pbc_dx(pbc,
		     x[acceptor_atom->protons[0]],
		     x[acceptor_atom->atom_id],
		     OH1);
	      pbc_dx(pbc,
		     x[acceptor_atom->protons[1]],
		     x[acceptor_atom->atom_id],
		     OH2);
	      pbc_dx(pbc,
		     x[acceptor_atom->protons[out_of_plane_proton]],
		     x[acceptor_atom->atom_id],
		     w_third);

	      cprod(OH1, OH2, w_normal);
	    }

	  unitv(w_normal, w_normal); /* Make water normal unit vector */

	  planedist = iprod(w_normal, w_third);  /* project out-of-plane proton onto unit water normal*/
	  svmul(-2*planedist, w_normal, w_normal);
	  rvec_inc(w_third, w_normal);
	  copy_rvec(w_third, x[acceptor_atom->protons[out_of_plane_proton]]); /* Mirror the position. */
	}
    }

  /* alter the protonation state of the donor and acceptor atoms, but
     only if donor is still a donor, and acceptor still an acceptor 
  */

  /* Change the following part according to the new machinery. */

/*   if( (donor_atom->bdonor && !acceptor_atom->bdonor) ||  */
/*      (acceptor_atom->bdonor && !donor_atom->bdonor) ){ */
/*     b                    = donor_atom->bdonor; */
/*     donor_atom->bdonor    = acceptor_atom->bdonor; */
/*     acceptor_atom->bdonor = b; */
    
/*     b = md->bqhopdonor[donor_atom->atom_id]; */
/*     md->bqhopdonor[donor_atom->atom_id]= */
/*       md->bqhopdonor[acceptor_atom->atom_id]; */
/*     md->bqhopdonor[acceptor_atom->atom_id] = b; */
    
/*     b = md->bqhopacceptor[donor_atom->atom_id]; */
/*     md->bqhopacceptor[donor_atom->atom_id]= */
/*       md->bqhopacceptor[acceptor_atom->atom_id]; */
/*     md->bqhopacceptor[acceptor_atom->atom_id] = b;  */
    
/*     get_protonation_state(cr, qhoprec); */
/*     /\* alter the charges in the mdatoms struct  */
/*      *\/ */
/*     for(i=0;i<donor->nr_atoms;i++){ */
/*       md->chargeA[donor->atoms[i]] = donor->charge_set[donor->pryte][i];  */
/*     } */
/*     for(i=0;i<acceptor->nr_atoms;i++){ */
/*       md->chargeA[acceptor->atoms[i]]=acceptor->charge_set[acceptor->pryte][i]; */
/*     } */
/*     return (TRUE); */
/*   } */
/*   else{ */
/*     return(FALSE); */
/*   } */

  if( ((donor_atom->state & eQDON) && (acceptor_atom->state & eQACC))
      || ((acceptor_atom->state & eQDON) && (donor_atom->state & eQACC)) )
    {
      if (bUndo)
	{
	  /* These functions aren't really implemented yet. */
	  qhop_protonate  (db, qhoprec, donor_atom,    bDonWater, FALSE);
	  qhop_deprotonate(db, qhoprec, acceptor_atom, bAccWater, FALSE);
	}
      else
	{
	  /* These functions aren't really implemented yet. */
	  qhop_protonate  (db, qhoprec, acceptor_atom, bAccWater, FALSE);
	  qhop_deprotonate(db, qhoprec, donor_atom,    bDonWater, FALSE);
	}
    }
  
  if (!bAccWater)
    {
      if (bUndo)
	{
	  /* Put proton back */
	  copy_rvec(x_old, x[acceptor_atom->protons[hop->proton_id]]);
	}
      else
	{
	  /* We need to position the hydrogen correctly */
	  
	  /* Which proton?
	   * For now, assume it's only one. */
	  

	  if (x_old != NULL)
	    {
	      copy_rvec(x[acceptor_atom->protons[hop->proton_id]], x_old); /* Save old position. Needed for undo */
	    }
	  
	  /* Get old AH-distance (assume it's correct) */
	  pbc_dx(pbc,
		 x[donor_atom->atom_id],
		 x[acceptor_atom->atom_id],
		 x_tmp);
	  d = norm(x_tmp);

	  /* Reposition the proton on the DA-line */
	  pbc_dx(pbc,
		 x[donor_atom->atom_id],
		 x[acceptor_atom->atom_id],
		 x_tmp);
	  unitv(x_tmp, x_tmp);
	  svmul(d, x_tmp, x_tmp);
	  rvec_add(x[acceptor_atom->atom_id], x_tmp, x[acceptor_atom->protons[0]]);
	}
    }

  /* Change md->chargeA[], md->invmass[] upon (de)protonation. */

} /* change_protonation */

static real evaluate_energy(t_commrec *cr,t_inputrec *ir, t_nrnb *nrnb,
			    gmx_wallcycle_t wcycle, 
			    gmx_localtop_t *top,gmx_mtop_t *mtop, 
			    gmx_groups_t *groups,t_state *state,
			    t_mdatoms *md,t_fcdata *fcd,
			    t_graph *graph,
			    t_forcerec *fr,gmx_vsite_t *vsite,rvec mu_tot,
			    /*gmx_genborn_t *born,*/
			    /* gmx_bool bBornRadii, */int step){
  real 
    t,etot;
  real 
    terminate;
  gmx_enerdata_t 
    *enerd;
  tensor 
    force_vir;

  terminate = 0;
  snew  (enerd,1);
  init_enerdata(groups->grps[egcENER].nr, ir->n_flambda, enerd);

  if (vsite)
    {
      construct_vsites(NULL,vsite,state->x,nrnb,1,NULL,
		       top->idef.iparams,top->idef.il,
		       fr->ePBC,fr->bMolPBC,graph,cr,state->box);
    }
  
  do_force(stderr,cr,ir,step,nrnb,wcycle,top,mtop,groups,
	   state->box,state->x,&state->hist,
	   /*NULL,*/NULL,force_vir,md,enerd,fcd,
	   state->lambda,graph,
	   fr,vsite,mu_tot,step*ir->delta_t,NULL,NULL,/*born*//* bBornRadii */ FALSE,
	   GMX_FORCE_STATECHANGED|GMX_FORCE_NS );
  //	   GMX_FORCE_ALLFORCES | (GMX_FORCE_VIRIAL ));

  /* Communicate stuff when parallel */
#ifdef GMX_MPI

  if (PAR(cr))
    {
    //    wallcycle_start(wcycle,ewcMoveE);

      global_stat(NULL,cr,enerd,NULL,NULL,mu_tot,
		  ir,NULL,FALSE,NULL,NULL,NULL,NULL,&terminate);

    //    wallcycle_stop(wcycle,ewcMoveE);
  }
#endif

  etot = enerd->term[F_EPOT];
  free(enerd);
  return(etot);
} /* evaluate_energy */

static void print_qhop_info(t_qhoprec *qhoprec){
  int
    i,j,k;
  /* print qhop atoms */
  for(i=0; i < qhoprec->nr_qhop_atoms; i++)
    {
      fprintf(stderr,"qhop atom [%d], atom_id [%d], res [%d] bDonor = %d, %d protons\n",i,
	      qhoprec->qhop_atoms[i].atom_id,
	      qhoprec->qhop_atoms[i].res_id,
	      /*qhoprec->qhop_atoms[i].bdonor*/0,
	      qhoprec->qhop_atoms[i].nr_protons);
    }
  for(i=0; i < qhoprec->nr_qhop_residues; i++)
    {
/*     fprintf(stderr,"\nqhop residue %d [%d]: charges, set used: %d\n\n", */
/* 	    i,qhoprec->qhop_residues[i].res_nr,qhoprec->qhop_residues[i].pryte); */
      for(k=0; k < qhoprec->qhop_residues[i].nr_atoms; k++)
	{
	  fprintf(stderr,"%3d ",qhoprec->qhop_residues[i].atoms[k]);
/*       for(j=0;j<qhoprec->qhop_residues[i].max_state;j++){ */
/* 	fprintf(stderr,"%6.3f ",qhoprec->qhop_residues[i].charge_set[j][k]); */
/*       } */
	  fprintf(stderr,"\n");
	}
    }
} /* print_qhop_info */

void print_hop(t_hop hop)
{
  fprintf(stderr,"donor atom %d, acceptor atom %d, proton %d\n",
	  hop.donor_id,hop.acceptor_id,hop.proton_id);
} /* print_hop */

static real get_self_energy(t_commrec *cr, t_qhop_residue donor, 
			    t_qhop_residue acceptor, t_mdatoms *md, 
			    rvec *x, t_pbc pbc, t_forcerec *fr){
  /* the internal energy of the hoppers together Is difficult because
     the internal energy of the individal residues, eg. if part of
     protein, will also change upon hopping. For water does not
     matter.
   */
  int
    i,j;
  real
    qq, Eself, ek, ec, eps, drsq;
  rvec
    vec;

  Eself = 0.0;
  ek    = 0.0;
  ec    = 0.0;
  
  if (EEL_RF(fr->eeltype))
    {
      ek = fr->k_rf;
      ec = fr->c_rf;
    }

  eps = fr->epsfac;

  for(i=0; i < donor.nr_atoms; i++)
    {
      for(j=0; j < acceptor.nr_atoms; j++)
	{
	  qq = md->chargeA[donor.atoms[i]] * md->chargeA[acceptor.atoms[j]];
	  pbc_dx(&pbc,
		 x[donor.atoms[i]],
		 x[acceptor.atoms[j]],
		 vec);
	  drsq = norm2(vec);
	  Eself +=  qq*(gmx_invsqrt(drsq)+ek*drsq-ec);
	}
    }

  Eself *= eps;

  /* reaction field correction STILL needed!!!! nu even geen zin in*/
  return(Eself);
} /* evaluate energy */

static real calc_S(const qhop_parameters *p,
		   const t_hop *hop)
{
  real d = hop->rda - (p->t_A);

  return (p->s_A) * (d * d) + (p->v_A);
}

static real calc_V(const qhop_parameters *p,
		   const t_hop *hop)
{
  return (p->s_C) * exp( -(p->t_C) * (hop->rda - 0.20)) + (p->v_C); 
}


/**
 * Return the transfer probability for an arbitrary time window
 *
 * Models the transfer as a poisson process,
 * p(t) =  1-exp(-k*t) = 1 - ( 1-p(0.01 ps) )^t
 * 
 */
static real poisson_prob(real p,            ///< Hopping probability over a 0.01 ps window
			 real qhop_interval ///< The interval between qhop attempts in picoseconds
			 )
{
  return 1-pow((1-p), (0.01 / qhop_interval));
}


/**
 * Calculates the SE-regime validity limit and stores it in hop->El.
 */
static void compute_E12_left(qhop_parameters *p, ///< Pointer to hopping parameters			      
		      t_hop *hop          ///< The hop
		      )
{
  
  real
    K, M;

  /* compute the value of E12 that is required for 
   * p_SE(rda,10fs) = 0.1. 
   */
  K = p->k_1 * exp(-p->k_2 * (hop->rda-0.23)) + p->k_3;
  M = p->m_1 * exp(-p->m_2 * (hop->rda-0.23)) + p->m_3;

  hop->El = -(atanh( 2*(0.1-0.5) )-M) / K;
}


/**
 * Calculates the TST-regime validity limit and stores it in hop->Er.
 */
static void compute_E12_right(qhop_parameters *p, ///< Pointer to hopping parameters
		       t_hop *hop,         ///< The hop
		       real Temp	   ///< Temperature
		       )
{
  
  real
    S,T,V,Eb,d;

  /* compute the value of E12 at which Exp(Eb/kT)=100
   * Eb = S +T*E12 + V*E12*E12
   *
   *      -T + Sqrt(T*T - 4*V*(S-Eb))
   * E12 = --------------------------
   *                 2*V
   *
   * where Eb = 1/beta*ln(100)
   */
  
  /*  Eb = (BOLTZ*Temp)*log(100);
   *  No, we must use RGAS, since the units are in kJ/mol.
   */
  Eb = RGAS * Temp * log(100);
  d = hop->rda - (p->t_A);
  S = calc_S(p, hop);
  T = p->s_B;
  V = calc_V(p, hop);

  d = T*T - 4*V*(S-Eb); /* I want to make sure this is non-negative. */

  hop->Er = (-T + sqrt(d)) / (2*V);

  /*   E12_right = ((-99.50-0.0760*Temp)*(rda*10.-(2.3450+0.000410*Temp))*(rda*10.-(2.3450+0.000410*Temp))+10.3)*CAL2JOULE; */
} /* compute_E12_right */


/**
 * Calculates the barrier height.
 */
static real compute_Eb(qhop_parameters *p, ///< Pointer to qhop_parameters			   
		       t_hop *hop
		       )
{
  real
    Eb, S, T, V, temp;
  
  /* temp = rda - p->t_A; */
  S = calc_S(p, hop); /* p->s_A * (temp*temp) + p->v_A; */
  T = p->s_B;
  V = calc_V(p, hop); /* p->s_C*exp(-p->t_C*(rda-0.20))+p->v_C; */
  Eb = S + T*hop->E12 + V*hop->E12*hop->E12;
  return(Eb);
}

/**
 * Calculates the TST transfer probability for a 10 fs time window.
 */
real compute_rate_TST(qhop_parameters *p, ///< Pointer to qhop_parameters
		      t_hop *hop,         ///< The hop
		      real T              ///< Temperature
		      )
{
  real
    Q, R, kappa, half_hbar_omega, ETST, pTST, Emax, E_M, Eb, E12, rda;

  rda = hop->rda;
  E12 = hop->Eb;
  Eb = hop->Eb;

  if(E12 > 0.0)
    {
      Emax = Eb;
      E_M  = Emax - E12;
    }
  else
    {
      Emax = Eb - E12;
      E_M = Eb;
    }

  Q =
    p->q_1
    + p->q_2 * T
    + p->q_3 * T*T;

  R =
    p->r_1
    + p->r_2 * T
    + p->r_3 * T*T;

  kappa =
    exp(p->p_1 + Q*E_M + R*E_M*E_M);

  half_hbar_omega = p->f * exp(-p->g * Eb) + p->h;
  ETST =-(Eb-half_hbar_omega)/(RGAS * T);
  pTST = (kappa*BOLTZ*T/PLANCK) * exp(ETST) * 0.01;

  hop->hbo = half_hbar_omega;

  return (pTST);
} /* compute_prob_TST */

/**
 * Calculates the barrierless transfer probability for a 10 fs time window.
 */
static real compute_rate_SE(qhop_parameters *p, ///< Pointer to qhop_parameters
			    t_hop *hop          ///< The hop
			    )
{
  real
    K, M, pSE;
  
  K = p->k_1 * exp(-p->k_2 * (hop->rda - 0.23)) + p->k_3;
  M = p->m_1 * exp(-p->m_2 * (hop->rda - 0.23)) + p->m_3;
  pSE = (0.5 * tanh(-K * hop->E12 + M) + 0.5);
  return(pSE);
} /* compute_prob_SE */


/**
 * Calculates the transfer probability for the intermediate regime over a 10 fs time window.
 */
real compute_rate_log(qhop_parameters *p, ///< Pointer to qhop_parameters
		      t_hop *hop,         ///< The hop 
		      real T              ///< Temperature
		      )
{
  
   /* we now compute the probability over a 10fs interval. We return the
   * probability per timestep.
   */
  real 
    rSE,rTST,rate;
  
  rSE  = compute_rate_SE(p, hop);
  rTST = compute_rate_TST(p, hop, T);
  rate = rSE * pow((rTST/rSE), (hop->E12 - hop->El)/(hop->Er - hop->El));
  /* See, the log-space interpolation translates to the expression in the first argument above.
   * No need to log/exp back and forth. */

  return(rate);
} /* compute_prob_log */

static real compute_E12_0(t_qhop_params *p, t_hop *hop)
{
  return
    p->alpha
    + p->beta * hop->rda
    + p->gamma * hop->rda * hop->rda;
}

static real get_hop_prob(t_commrec *cr,t_inputrec *ir, t_nrnb *nrnb,
			 gmx_wallcycle_t wcycle, 
			 gmx_localtop_t *top,gmx_mtop_t *mtop, 
			 gmx_groups_t *groups,t_state *state,
			 t_mdatoms *md,t_fcdata *fcd,
			 t_graph *graph,
			 t_forcerec *fr,gmx_vsite_t *vsite,rvec mu_tot,
			 /* gmx_genborn_t *born, gmx_bool bBornRadii, */
			 t_hop *hop, real T,
			 t_qhoprec *qhoprec,t_pbc pbc,int step, qhop_db *db){
  
  /* compute the hopping probability based on the Q-hop criteria
   */
  real
    Ebefore_tot,Ebefore_self,Ebefore,
    Eafter_tot, Eafter_self, Eafter,
    E12,E12_right,E12_left,Eb,
    r_TST,r_SE,r_log;
  qhop_parameters 
    *p;
  rvec
    x_old;

  /* liever lui dan moe */
  p = get_qhop_params(qhoprec->qhop_atoms[hop->donor_id].resname,
		      qhoprec->qhop_atoms[hop->acceptor_id].resname, db);
  
  Ebefore_tot = evaluate_energy(cr, ir, nrnb, wcycle,top, mtop,
				groups, state, md, fcd, graph,
				fr,vsite,mu_tot, /* born, bBornRadii, */step);
  Ebefore_self = get_self_energy(cr, fr->qhoprec->qhop_residues[hop->donor_id],
				 fr->qhoprec->qhop_residues[hop->acceptor_id],
				 md, state->x, pbc, fr);
  Ebefore = Ebefore_tot-Ebefore_self;

  if(change_protonation(cr, fr->qhoprec, md, hop, state->x,FALSE,mtop, x_old, &pbc, db))
    {
      Eafter_tot = evaluate_energy(cr,ir,nrnb,wcycle,top,mtop,groups,state,md,
				   fcd,graph,fr,vsite,mu_tot,/*  born, bBornRadii, */
				   step);
      Eafter_self = get_self_energy(cr, fr->qhoprec->qhop_residues[hop->donor_id],
				    fr->qhoprec->qhop_residues[hop->acceptor_id],
				    md, state->x, pbc, fr);
      Eafter = Eafter_tot-Eafter_self;
    
      E12 = compute_E12_0(p, hop) + Eafter-Ebefore;

      compute_E12_right(p, hop, T);
    
      compute_E12_left (p, hop);
      Eb        = compute_Eb(p,hop->rda,E12);  
      fprintf(stderr,"E12 = %f\n",E12);

      

      r_TST = compute_rate_TST(p, hop, T);
      r_SE  = compute_rate_SE (p, hop);

      if(E12 > hop->Er)
	{
	  /* Classical TST regime */
	  hop->prob = poisson_prob(r_TST, (ir->delta_t * qhoprec->qhopfreq)); /* (1-pow((1-r_TST),fr->qhoprec->qhopfreq)); */
	  hop->regime = etQhopTST;
	}w
      else 
	{
	  if (E12 < hop->El)
	    {
	      /* Schroedinger regime */
	      /* using volkhards probality propagation:
	       */
	      hop->prob = poisson_prob(r_SE, (ir->delta_t * qhoprec->qhopfreq));
	      hop->regime = etQhopSE;
	      fprintf(stderr, "TEMPORARY:       r_SE = %e\n", r_SE);
	    }
	  else
	    {
	      r_log = compute_rate_log(p, hop, T);
	      /* intermediate regime */
	      /* using volkhards probality propagation */
	      hop->prob = poisson_prob(r_log, (ir->delta_t * qhoprec->qhopfreq));
	      hop->regime = etQhopI;
	    }
	}

      /* now undo the move */
      if(!change_protonation(cr, fr->qhoprec, md, hop, state->x, TRUE, mtop, x_old, &pbc, db))
	{
	  gmx_fatal(FARGS,"Oops, cannot undo the change in protonation.");
	}
    }
  else
    {
      fprintf(stderr, "change_protonation() returned FALSE!");
    }
  
  /* return the FF energy difference, as this is what we need to
     compensate by scaling velocities
  */
  hop->DE_MM = Eafter_tot-Ebefore_tot;
  hop->E12 = E12;
  hop->Eb  = Eb;
  /* free(p); */

  /* return(Eafter_tot-Ebefore_tot); */
}


static gmx_bool swap_waters(t_commrec *cr, t_qhoprec *qhoprec, 
			t_mdatoms *md, t_hop *hop, rvec *x,
			gmx_mtop_t *mtop){
  int
    i;
  rvec
    temp_vec;
  
  /* swap the donor and acceptor atom
   */  
  /* assumuming the user uses the correct water topology..... */
  switch (hop->proton_id)
    {
    case 3:
      rotate_water(x,&qhoprec->qhop_atoms[hop->donor_id],-120,md);
      break;
    case 4:
      rotate_water(x,&qhoprec->qhop_atoms[hop->donor_id], 120,md);
      break;
    default:
      break;
    }
  
  copy_rvec(x[qhoprec->qhop_atoms[hop->donor_id].atom_id],
	    temp_vec);
  copy_rvec(x[qhoprec->qhop_atoms[hop->acceptor_id].atom_id],
	    x[qhoprec->qhop_atoms[hop->donor_id].atom_id]);
  copy_rvec(temp_vec,
	    x[qhoprec->qhop_atoms[hop->acceptor_id].atom_id]);
  
  if(qhoprec->qhop_atoms[hop->donor_id].nr_protons ==
     qhoprec->qhop_atoms[hop->acceptor_id].nr_protons)
    {
      for(i=0; i < qhoprec->qhop_atoms[hop->donor_id].nr_protons; i++)
	{
	  copy_rvec(x[qhoprec->qhop_atoms[hop->donor_id].protons[i]],
		    temp_vec);
	  copy_rvec(x[qhoprec->qhop_atoms[hop->acceptor_id].protons[i]],
		    x[qhoprec->qhop_atoms[hop->donor_id].protons[i]]);
	  copy_rvec(temp_vec,
		    x[qhoprec->qhop_atoms[hop->acceptor_id].protons[i]]);
	}

      return(TRUE);
	
    }
  else
    {
      gmx_fatal(FARGS,"Oops, donor and acceptor do not have same number of protons!\n");
      return(FALSE);
    }
}
	    
static gmx_bool do_hop(t_commrec *cr, t_qhoprec *qhoprec, 
		   t_mdatoms *md, t_hop *hop, rvec *x,
		   gmx_mtop_t *mtop){
  /* change the state of the system, such that it corresponds to the
     situation after a proton transfer between hop.donor_id and
     hop.acceptor_id. For hops from hydronium to water we swap donor
     and acceptor. All other hops are carried out by changes in the
     charges in mdatoms.
  */
  gmx_bool
    bOk, b;
  
/*   if(qhoprec->qhop_atoms[hop->donor_id].bWater &&  */
/*      qhoprec->qhop_atoms[hop->acceptor_id].bWater ){ */

/*     if(qhoprec->qhop_atoms[hop->donor_id].bdonor && */
/*        !qhoprec->qhop_atoms[hop->acceptor_id].bdonor ){ */
/*           fprintf(stderr,"swappen maar!!!\n");   */
/* 	  //      b = qhoprec->qhop_atoms[hop->donor_id].bdonor; */
/* 	  //    qhoprec->qhop_atoms[hop->donor_id].bdonor    =  */
/* 	  //	qhoprec->qhop_atoms[hop->acceptor_id].bdonor; */
/* 	  // qhoprec->qhop_atoms[hop->acceptor_id].bdonor = b; */
/* 	  bOk = swap_waters(cr, qhoprec, md, hop, x,mtop); */
/*     } */
/*     else bOk=FALSE; */
/*   } */
/*   else{ */
/*     bOk = change_protonation(cr, qhoprec, md, hop,  */
/* 			     x,FALSE,mtop); */
/*   } */
  return(bOk);
}

real distance_dependence(rvec x, rvec c, real rc, t_pbc pbc){
  real
    r,f;
  rvec
    dx;
  
  f = 0;

  pbc_dx(&pbc, x, c, dx);
  r = norm(dx);

  if(r < rc)
    {
      f = 1 - exp(-r/rc);
    }
  return(f);
}  /* distance_dependence */




static real check_ekin(rvec *v, t_mdatoms *md){
  int
    i,j;
  real
    vsq,    ekin;

  for(i=0;i<md->nalloc;i++){
    vsq = 0;
    for(j=0;j<DIM;j++){
      vsq += (v[i][j]*v[i][j]);
    }
    ekin += 0.5*md->massT[i]*vsq;
  }
  
  return(ekin);
} /* check_ekin */

static void scale_v(rvec *x, rvec *v, t_mdatoms *md, 
		    int donor_atom, int acceptor_atom,
		    real rc, real DE,t_pbc pbc){
  /* computes new velocities, to get ridf of the DE. Returns the new
     kinetic energy;
  */
  real
    Cinv, C, a, b, c, *f, f2, mass, lambda;
  rvec
    center, p, q, dv;
  int
    i, j;

  /* make sure everythins is set to zero 
   */
  Cinv = 0;
  a = 0;
  b = 0;
  c = 0;
  clear_rvec(p);
  clear_rvec(q);
  clear_rvec(dv);
  snew(f, md->nalloc);

  /* compute the center of the reactants. 
   */
  for(j=0; j<DIM ; j++)
    {
      center[j] = 0.5*(x[donor_atom][j] + x[acceptor_atom][j]);
    }

  /* first compute C, p and q */  
  for(i=0; i < md->nalloc; i++)
    {
      mass = md->massT[i];
      f[i] = distance_dependence(x[i], center, rc, pbc);
      Cinv += mass*f[i];
      for(j=0; j < DIM; j++)
	{
	  p[j] += mass * v[i][j];
	  q[j] += mass * v[i][j] * f[i];
	}	
    }

  C = 1/Cinv;

  /* with C, p and q, we can compute a, b and c that we need to find
     the solution for the quadratic equation of David and Erik M.
  */

  a = 0;
  b = 0;
  c = 2*DE;
  for(i=0; i < md->nalloc; i++)
    {
      mass = md->massT[i];
      f2 = f[i]*f[i];
      a += mass * f2 * ( iprod(v[i], v[i])
			 - 2 * C * iprod(v[i], q)
			 + C * C * iprod(q, q));
    
      b -= 2 * mass * f2 * ( C * C * iprod(p, q)
			     - C * iprod(v[i], p));
    
      c += mass * ( f2 * C * C * iprod(p, p)
		    - iprod(v[i], v[i]));
  }
  
  /* with a,b,c we can find lambda. We need only the positive solution : */
  lambda = (-b + sqrt(b*b - 4*a*c))/(2*a);

  if(lambda < 0.0000000 )
    {
      free(f);
      return;
    }

  /* with lambda, we compute dv:
   */
  
  for(j=0; j<DIM; j++)
    {
      dv[j] = C * (p[j] - lambda*q[j]) / lambda;
    }

  /* with lambda and dv we scale the velocities of all particles
   */
  for(i=0; i < md->nalloc; i++)
    {
      for(j=0; j < DIM; j++)
	{
	  v[i][j] = lambda * f[i] * (v[i][j] + dv[j]);
	}
    }

  free(f);
} /* scale_v */

static gmx_bool scale_velocities(t_commrec *cr,t_inputrec *ir, t_nrnb *nrnb,
				 gmx_wallcycle_t wcycle,  gmx_localtop_t *top,
				 gmx_mtop_t *mtop, gmx_groups_t *groups,
				 t_state *state,  t_mdatoms *md, 
				 t_qhoprec *qhoprec,int step,
				 gmx_constr_t constr,real DE,t_pbc pbc,  
				 t_hop *hop,gmx_ekindata_t *ekindata,
				 gmx_bool bPscal,real veta,real vetanew)
{
  /* takes as input the total MM potential energy change for the hop
     and alters the kinetic energy so that the total energy remains
     constant.
  */
  gmx_bool
    bConverged, bConstrain;
  int
    donor_atom, acceptor_atom, iter;
  real
    ekin_old, ekin_new, dvdl;
  real
    ekin_before, DEpot;

  iter = 0;
  bConverged = FALSE;
  bConstrain = FALSE;

  /* iterate until the velocities satisfy both constraints and
   * energy/momentum conservation
   */
  if (constr != NULL)
    {
      bConstrain=TRUE;
    }

  donor_atom    = qhoprec->qhop_atoms[hop->donor_id].atom_id;
  acceptor_atom = qhoprec->qhop_atoms[hop->acceptor_id].atom_id;

  ekin_before = check_ekin(state->v,md);
  ekin_old = ekin_before;

  /* copy the old velocities 
   */
  DEpot = DE;
  do
    {

      scale_v(state->x, state->v,
	      md, donor_atom, acceptor_atom,
	      qhoprec->qhop_rc, DE, pbc);
    
      /* constrain the velocities 
       */ 
      fprintf(stderr,"before constr. ekin_new = %f\n",
	      check_ekin(state->v ,md));

      if(bConstrain)
	{
	  constrain(NULL,FALSE,FALSE,constr,&top->idef,ir,ekindata,cr,
		    step,1,md,state->x,state->v,state->v,state->box,
		    state->lambda,&dvdl,
		    NULL,NULL,nrnb,
		    econqVeloc,bPscal,veta,vetanew);
	  /* and send the veloities around again, bit expensive, 
	     so there must be an easier way
	     of doing this...
	  */
	  iter++;
	}

      ekin_new = check_ekin(state->v ,md);   
      fprintf(stderr,"iteration %d, ekin_new = %f, ekin_old = %f\n",
	      iter,ekin_new,ekin_old);


      DE = DE - (ekin_new - ekin_old);
      ekin_old = ekin_new;
      bConverged = (DE*DE < 0.001);/* hard coded treshold.... */
    } while( !bConverged && iter < 100 );

  fprintf(stderr,"totat energy correction: %f, DE_MM: %f\n",
	  ekin_new-ekin_before, DEpot);

  return (bConverged);
} /* scale_velocities */
  
/**
 * Changes the order of the hops
 */
static t_hop* scramble_hops(t_hop *hop,    ///< List of hops
			    int nr_hops,   ///< Number of hops in the list
			    int mode,      ///< etQhopMode????
			    gmx_rng_t rng  ///< Random number generator
			    )
{
  int       i, j, r, hops_left;
  gmx_bool  bSorted;
  t_hop     *new_hop, tmp_hop;

  switch (mode)
    {
    case etQhopModeOne:
    case etQhopModeGillespie:
      new_hop = hop;
      break;

    case etQhopModeRandlist:

      /* make random order */
      {
	snew(new_hop, nr_hops);
	for (i=0; i<nr_hops; i++)
	  {
	    r = gmx_rng_uniform_uint32(rng) % (nr_hops - i);

	    new_hop[i] = hop[r];

	    /* Remove that hop from the list */
	    hop[r] = hop[nr_hops-i-1];
	  }

	sfree(hop);

	break;
      }

    case etQhopModeList:

      /* Sort according to probability
       * There's not gonna be a lot of hops, so just bubble sort them. */
      {
	bSorted = FALSE;
	while (!bSorted)
	  {
	    bSorted = TRUE;
	    for (i=0; i<nr_hops-1; i++)
	      {
		if (hop[i].prob < hop[i+1].prob)
		  {
		    tmp_hop  = hop[i];
		    hop[i]   = hop[i+1];
		    hop[i+1] = tmp_hop;
		    bSorted = FALSE;
		  }
	      }
	  }

	new_hop = hop;

	break;
      }

    default:
      gmx_fatal(FARGS, "That qhopMode is unsupported: %i", mode);
    }

  return new_hop;   
}

void do_qhop(FILE *fplog, t_commrec *cr,t_inputrec *ir, t_nrnb *nrnb,
	     gmx_wallcycle_t wcycle, 
	     gmx_localtop_t *top,gmx_mtop_t *mtop, gmx_groups_t *groups,
	     t_state *state,
	     t_mdatoms *md, t_fcdata *fcd,t_graph *graph, t_forcerec *fr,
	     gmx_vsite_t *vsite,rvec mu_tot,/*gmx_genborn_t *born, */
	     gmx_bool bBornRadii,real T, int step,
	     tensor force_vir, qhop_db_t db){

  int
    nr_hops,i,j,a;
  t_pbc
    pbc;
  static gmx_bool 
    bFirst=TRUE;
  gmx_bool
    bAgain,bHop;
  real 
    rnr;
  static int 
    start_seed=0;
  static gmx_rng_t 
    rng,rng_int;
  t_qhoprec
    *qr;

  set_pbc_dd(&pbc,fr->ePBC,DOMAINDECOMP(cr) ? cr->dd : NULL,FALSE,state->box);
  
  if(bFirst)
    {
      if(MASTER(cr))
	{
	  start_seed = gmx_rng_make_seed();
	}
#ifdef GMX_MPI
      if(PAR(cr))
	{
	  MPI_Bcast(&start_seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
#endif
      rng     = gmx_rng_init(12/*start_seed*/); 
      rng_int = gmx_rng_init(12/*start_seed*/); 
      bFirst  = FALSE;
    } 

  qr->hop = find_acceptors(cr, fr, state->x, pbc, &nr_hops, md);

  if(qr->nr_hops > 0)
    {
      for(i=0; i<nr_hops; i++)
	{
	  if (qr->qhopmode == etQhopModeOne)
	    {
	      i = (int) floor (nr_hops * (gmx_rng_uniform_uint32(rng_int))/4294967296.0);
	    }
	
	  qr->hop[i].DE_MM = get_hop_prob(cr,ir,nrnb,wcycle,top,mtop,groups,state,md,
					  fcd,graph,fr,vsite,mu_tot,/*born,bBornRadii,*/
					  &(qr->hop[i]),T,&(qr->hop[i].E12), &(qr->hop[i].Eb),
					  fr->qhoprec,pbc,step, db);
	  
	  if (qr->qhopmode == etQhopModeOne)
	    {
	      break;
	    }

	}

      /* now we have for all hops the energy difference and the
	 probability. For the moment I loop over them in order given and
	 keep looping until I am done */
    
      /* INstead, to perserve deaitled ballance, we need to randomly select one
       * 
       * Sorry, but that doesn't preserve detailed balance.
       */
    
      qr->hop = scramble_hops(qr->hop, qr->nr_hops, qr->qhopmode, rng);

      for(i=0; i < nr_hops; i++)
	{
	  if (qr->hop[i].regime != etQhopNONE)
	    {
	      rnr = gmx_rng_uniform_real(rng); 
	      if(MASTER(cr)){
		/* some printing to the log file */
		fprintf(fplog,
			"\n%d. don %d acc %d. E12 = %18.12f, DE_MM = %18.12f, Eb = %18.12f, prob. = %f, ran. = %f (%s)", i,
			fr->qhoprec->qhop_atoms[qr->hop[i].donor_id].res_id,
			fr->qhoprec->qhop_atoms[qr->hop[i].acceptor_id].res_id,
			qr->hop[i].E12, qr->hop[i].DE_MM, qr->hop[i].Eb,
			qr->hop[i].prob, rnr, qhopregimes[qr->hop[i].regime]);
	      }
 
	      /* Attempt hop */
	      if(qr->hop[i].prob > rnr)
		{
		  /* hoppenmaar! */
      
		  /* Move the actual hopping to update() */
		  bHop = do_hop(cr, fr->qhoprec, md, &(qr->hop[i]),state->x, mtop);
		  fprintf(stderr,"hopping!\n");
		  /*      scale_velocities();*/

		  if (qr->qhopmode != etQhopModeOne)
		    {
		      /* Zap all hops whose reactants just was consumed. */
		      for (j = i+1; j < nr_hops; j++)
			{
			  if (qr->hop[j].donor_id    == qr->hop[i].donor_id    ||
			      qr->hop[j].acceptor_id == qr->hop[i].acceptor_id ||
			      qr->hop[j].donor_id    == qr->hop[i].acceptor_id ||
			      qr->hop[j].acceptor_id == qr->hop[i].donor_id)
			    {
			      qr->hop[j].regime = etQhopNONE;
			      fprintf(stderr, " * Zapping hop number %i *\n", j);
			    }
			}
		    }
		}
	      else
		{
		  bHop = FALSE;
		}

	      if(MASTER(cr) && bHop)
		{
		  fprintf(fplog,"\n\nQ-hop is TRUE at step %d!\nE12 = %f, hopper: %d don: %d (%s) acc: %d (%s)\n",
			  step,qr->hop[i].E12, i,
			  fr->qhoprec->qhop_atoms[qr->hop[i].donor_id].res_id,
			  fr->qhoprec->qhop_atoms[qr->hop[i].donor_id].resname,
			  fr->qhoprec->qhop_atoms[qr->hop[i].acceptor_id].res_id,
			  fr->qhoprec->qhop_atoms[qr->hop[i].acceptor_id].resname);
		}
	    }
	}
      fprintf(fplog,"\n");
    }
}

/* patamStr is the string form the rtp file.
 * bts = ebtsBONDS, ..., ebtsCMAP
 * bt is the tyep found in t_rbonded.type
 */

static void str2bonded(const char *paramStr,
		       const int bt, const int bts,
		       t_iparams *params, int *iatoms)
{

#ifndef _ALL_BPARAMS
#define _ALL_BPARAMS &p[0],&p[1],&p[2],&p[3],&p[4],&p[5],&p[6],&p[7],&p[8],&p[9],&p[10],&p[11] /* p[0:MAXFORCEPARAM] */
#endif

  int ft, i, niatoms=btsNiatoms[bts];
  real p[MAXFORCEPARAM];

  int np, strOffset;
  char format[1+2+MAXATOMLIST*2+MAXFORCEPARAM*2];
  /* 2 for the functype,
   * MAXATOMS*2 for the %i belonging to each atom,
   * MAXFORCEPARAMS*2 for the %f belonging to each forceparameter,
   * 1 for the '\0'
   */
  memset(p, 0, sizeof(real)*MAXFORCEPARAM);
  memset(&params, 0, sizeof(t_iparams)); /* zap the params to zero */
  strOffset = btsNiatoms[bts]*2+1; /* +1 for the initial whitespace, *2 for the "%i" */
  sprintf(format," %s","%i%i%i%i%i%i%i"); /* the atoms */
  for (i=0; i<MAXFORCEPARAM; i++)
    sprintf(&(format[i*2+2+strOffset]),"%%f");
  
  switch(btsNiatoms[bts])
    {
    case 2:
      np = sscanf(paramStr, format,
		  iatoms[0],
		  iatoms[1],
		  _ALL_BPARAMS);
      break;
    case 3:
      np = sscanf(paramStr, format,
		  iatoms[0],
		  iatoms[1],
		  iatoms[2],
		  _ALL_BPARAMS);
      break;
    case 4:
      np = sscanf(paramStr, format,
		  iatoms[0],
		  iatoms[1],
		  iatoms[2],
		  iatoms[3],
		  _ALL_BPARAMS);
      break;
    case 5:
      np = sscanf(paramStr, format,
		  iatoms[0],
		  iatoms[1],
		  iatoms[2],
		  iatoms[3],
		  iatoms[4],
		  _ALL_BPARAMS);
      break;
    case 6:
      np = sscanf(paramStr, format,
		  iatoms[0],
		  iatoms[1],
		  iatoms[2],
		  iatoms[3],
		  iatoms[4],
		  iatoms[5],
		  _ALL_BPARAMS);
      break;
    default:
      /* This reeeeally shouldn't happen, unless a new bonded type is added with more than 5 atoms. */
      gmx_fatal(FARGS, "Wrong number of atoms requested: %d", btsNiatoms[bts]);
    }

  /*   ====  Now set the params  ====
   * Note that tabulated potentials, vsiten, orires and cmap
   * are not supported at this stage, simply because such information
   * is not available in the rtp file. Other bonded types are striktly
   * speaking also not available, e.g. vsite, but since their parameters
   * all reals it's trivial to just pass the parameters to the t_iparams.
   */
  if (bts==ebtsPDIHS && bt==1) /* Ordinary proper dihedral. */
    {
      params->pdihs.phiA = p[0];
      params->pdihs.cpA  = p[1];
      params->pdihs.mult = (int)p[2]; /* <- This is why this one needs special treatment */
      params->pdihs.phiA = p[3];
      params->pdihs.cpA  = p[4];
    }
  else
    {
      for (i=0; i<MAXFORCEPARAM; i++)
	(*params).generic.buf[i] = p[i];
    }
}

extern void qhop_stash_bonded(qhop_db_t db, gmx_mtop_t *mtop)
{
  int rt, res, nrt, nres, bt, bi,rtpIndex,
    natoms, nqatoms, i, j, ft,
    iatoms[MAXATOMLIST],
    iatomsMapped[MAXATOMLIST],
    *atomMap=NULL,
    ilistEntry[MAXATOMLIST+1];

  /* WE NEED TO DETERMINE THE FUNCTYPE ft !!!!*/

  t_iparams params;

  nrt = db->rb.nrestypes;
  for (rt=0; rt<nrt; rt++)
    {
      nres = db->rb.nres[rt];
      nqatoms = db->rtp[db->rb.rtp[rt]].natom;
      for (res=0; res<nres; res++)
	{
	  rtpIndex = db->rb.res[rt][res].rtp;

	  /* Now map the atoms in the residue subtype (XXX)
	   * to atoms in the residue type (qXXX) */
	  natoms = db->rtp[rtpIndex].natom;
	  srenew(atomMap, natoms);
	  for (i=0; i<natoms; i++)
	    {
	      atomMap[i]=-1;
	      for (j=0; j<nqatoms; j++)
		{
		  /* Match them on a name basis */
		  if (gmx_strcasecmp(*(db->rtp[rtpIndex].atomname[i]),
				     *(db->rtp[db->rb.rtp[rt]].atomname[j])) == 0)
		    {
		      /* It's a match! */
		      atomMap[i] = j;
		      break;
		    }
		  if (atomMap[i] == -1)
		    gmx_fatal(FARGS, "Could not map atom %s in %s onto the atoms in %s.",
			      *(db->rtp[rtpIndex].atomname[i]),
			      db->rtp[rtpIndex].resname,
			      db->rtp[db->rb.rtp[rt]].resname);
		}
	    }

	  /* Loop over all bonded interactions for this residue */
	  for (bt=0; bt<ebtsNR; bt++)
	    {
	      for (bi=0;
		   bi<db->rtp[rtpIndex].rb[bt].nb;
		   bi++)
		{
		  /* Parse the params and atoms. */
		  str2bonded(db->rtp[rtpIndex].rb[bt].b[bi].s,
			     bt,
			     bt<=ebtsIDIHS ? db->bts[bt] : 0,
			     &params,
			     iatoms);

		  /* Store the functype/iparams index */
		  ilistEntry[0] = gmx_mtop_append_itype(mtop, ft, params);

		  /* find corresponding atoms in qXXX via iatomMap.
		   * note that the atoms in ilistEntry[1,...] are relative
		   * to the start of the residue, so they must be offset when
		   * actually adding this stuff to the ilist. */
		  for (i=0; i<btsNiatoms[bt]; i++)
		    ilistEntry[i+1] = db->rb.res[rt][res].iatomMap[iatoms[i]];
		}
	    }
	}
    }
}

/* /\* Goes through the t_ilist and finds the bonded interactions */
/*  * that can be changed *\/ */
/* extern void qhop_index_bondeds(t_ilist *ilist, qhop_db_t db, */
/* 			       t_qhoprec *qr, gmx_bool bGlobal) */
/* { */
/*   gmx_bool bQhop, bInAtomList; */
/*   int i, j, k, qres, res, rt, bt, bi, ftype, btype, niatoms, */
/*     iatoms[MAXATOMLIST]; /\* The atoms involved in an interaction *\/ */
  
/*   t_restp *rtp; */

/*   /\* Loop over all qhopable residues *\/ */
/*   for (qres=0; qres<qr->nr_qhop_residues; qres++) */
/*     { */
/*       rt = qr->qhop_residues[qres].rtype; */
/*       rtp = &(db->rtp[db->rb.rtp[rt]]); */
/*       /\* Loop over bonded types *\/ */
/*       for (bt=0; bt<ebtsNR; bt++) */
/* 	{ */
/* 	  /\* Loop over bonded interactions of type bt *\/ */
/* 	  for (bi=0; bi<rtp->rb[bt].nb; bi++) */
/* 	    { */
/* 	      /\* Scan the ilist *\/ */
/* 	      for (i=0, bQhop=FALSE; */
/* 		   bQhop==FALSE && i<ilist[bt].nr; */
/* 		   i += btsNiatoms[bt]) */
/* 		{ */
/* 		  bQhop = TRUE; */
/* 		  ftype = ilist[bt].iatoms[i++]; */
/* 		  /\* j : Loop over the atoms in this interaction *\/ */
/* 		  for (j=0; j<btsNiatoms[bt]; j++) */
/* 		    { */
/* 		      /\* k : Loop over the atoms in this restype *\/ */
/* 		      for (k=0, bInAtomList == FALSE; */
/* 			   k<btsNiatoms[bt] && bInAtomList==FALSE && bQhop==TRUE; */
/* 			   k++) */
/* 			{ */
/* 			  if (db->rb.iatoms[rt][bt][bi][k] < 0) */
/* 			    { */
/* 			      /\* This means that the interaction goes beyond the */
/* 			       * residue boundry. Don't bother changing it. */
/* 			       * This may have to change eventually. *\/ */
/* 			      bQhop = FALSE; */
/* 			    } */
/* 			  /\* This asumes that the atoms[] is sorted in ascending order! */
/* 			   * Make sure to order it! *\/ */
/* 			  bInAtomList = */
/* 			    ilist[bt].iatoms[i] == */
/* 			    (db->rb.iatoms[rt][bt][bi][k] + qr->qhop_residues[qres].atoms[0]); */
/* 			} /\* k *\/ */
/* 		      bQhop &= bInAtomList; */
/* 		    } /\* j *\/ */
/* 		} /\* bi *\/ */
/* 	    } /\* i *\/ */
/* 	} /\* bt *\/ */
/*     } /\* qres *\/ */
/* } */
