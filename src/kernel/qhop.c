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
/* #include "qhop.h" */
#include "random.h"
#include "gmx_random.h"
#include "constr.h"
#include "types/idef.h"
#include "types/ifunc.h"
#include "hackblock.h"
#include "gmx_qhop_db.h"
#include "types/qhoprec.h"
#include "types/gmx_qhop_types.h"
#include "lambertw.h"

/* #define VERBOSE_QHOP */


#define DOO   0.35
#define BOND  0.15
#define HBOND  0.2
#define HOPANG 120


const char *waterName = "HXO";

static qhop_parameters *get_qhop_params (t_qhoprec *qr,t_hop *hop, qhop_db_t db)
{/* (char *donor_name,char *acceptor_name, qhop_db_t db){ */

  /* parameters of the current donor acceptor combination are based on 
   * the residue name.
   */

  char *donor_name, *acceptor_name;

  qhop_parameters
    *q;

  snew(q,1);

  donor_name    = db->rb.res
    [qr->qhop_residues[qr->qhop_atoms[hop->donor_id].qres_id].rtype]
    [qr->qhop_residues[qr->qhop_atoms[hop->donor_id].qres_id].res].name;
  acceptor_name = db->rb.res
    [qr->qhop_residues[qr->qhop_atoms[hop->acceptor_id].qres_id].rtype]
    [qr->qhop_residues[qr->qhop_atoms[hop->acceptor_id].qres_id].res].name;

/* 		      [(qhoprec->qhop_atoms[hop->donor_id].res)].name

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

static void set_proton_presence(qhop_H_exist *Hext, const int atomid, const gmx_bool present)
{
  Hext->H[Hext->atomid2H[atomid]] = present ? (char)1: (char)0;
}

static gmx_bool get_proton_presence(const qhop_H_exist *Hext, const int atomid)
{
  return Hext->H[Hext->atomid2H[atomid]] == (char)1;
}

static int ilist_iterator_start(const t_ilist *ilist, const int itype, int* index, t_iatom *atomlist)
{
  if (ilist->nr == 0)
    {
      *index = 0;
      atomlist = NULL;
      return 0;
    }

  *index = 0;
  atomlist = &(ilist->iatoms[0]);
  return interaction_function[itype].nratoms;
}

/* Returns the length of the interaction in atoms, including the type. */
static int ilist_iterator(const t_ilist *ilist, const int itype, int *index, t_iatom *atomlist)
{
  int ilen;
  if (*index < ilist->nr)
    {
      return 0;
    }
  
  /* Go to next interaction */
  ilen = interaction_function[itype].nratoms;
  *index += ilen + 1;         /* +1 for the type. */
  atomlist = &(ilist->iatoms[*index]);

  return ilen;
}

/* Is atom id among the qhop_atoms?
 * Returns the matching qhop_atom,
 * or -1 if not found.*/
static int is_qhopatom(const t_qhoprec *qr, const int id)
{
  int i;
  for (i=0; i<qr->nr_qhop_atoms; i++)
    {
      if (qr->qhop_atoms[i].atom_id == id)
	{
	  return i;
	}
    }

  return -1;
}

static int is_restype(const qhop_db *db, const t_qhoprec *qr, const char *rname)
{
  int i;

  if (rname == NULL || rname[0]=='\0')
    {
      return -1;
    }

  for (i=0; i<qr->nr_qhop_residues; i++)
    {
      if (strcmp(rname, db->rb.restype[qr->qhop_residues[i].rtype]) == 0)
	{
	  return i;
	}
    }

  return -2;
}

static void attach_proton(t_qhop_atom *qatom, const int id)
{
  int i;

  for (i=0; i<qatom->nr_protons; i++)
    {
      if (qatom->protons[i] == id)
	{
	  return;
	}
    }

  srenew(qatom->protons, qatom->nr_protons+1);
  qatom->protons[qatom->nr_protons++] = id;
}

static void scan_for_H(t_qhoprec *qr, const qhop_db *db)
{

#ifdef COMPLICATED_SCAN
  /*REDUNDANT (I hope)*/
  int i, j, a,  natoms, tsite, a1, a2, a3, qa;
  t_iatom *atom;

  natoms = ilist_iterator_start(ilist, itype, &i, atom);
   
  if (natoms == 0)
    {
      return;
    }

  /* Uses global atom numbers. Make this work with node-local numbering! */

  do
    {
      /* Remember that the interaction's first element is the interaction type,
       * so the first atom is ilist[1] */

      switch(itype)	
	{
	case F_BONDS:
	case F_G96BONDS:
	case F_MORSE:
	case F_CUBICBONDS:
	case F_CONSTR:
	case F_HARMONIC:
	case F_FENEBONDS:
	case F_TABBONDS:
	  for(a=1; a<=2; a++)
	    {
	      a1 = atom[a];
	      if (mtop->moltype[mt].atoms.atom[a1].atomnumber = 1) /* Make sure atomnumber is set! */
		{
		  /* is it connected to a t_qhop_atom? */
		  for (i=0; i<qr->nr_qhop_residues; i++)
		    {
		      if (qr->qhop_residues[i].rtype = rt)
			{
			  for (tsite = 0;
			       tsite<qr->qhop_residues[i].nr_titrating_sites;
			       tsite++)
			    {
			      a2 = atom[a%2 + 1];
			      qa = qr->qhop_residues[i].titrating_sites[tsite];
			      if (strcmp(qr->qhop_atoms[qa].atomname,
					 *(mtop->moltype[mt].atoms.atomname[a2])) == 0)
				{
				  attach_proton(&(qr->qhop_atoms[qa]), a2);
				  break;
				}
			    }
			}
		    }
		}
	    }
	  break;
	  
	case F_SETTLE:
	  /* This implies an oxygen and two attached hydrogens */
	  a1 = atom[1]; /* O */
	  a2 = atom[2]; /* H */
	  a3 = atom[3]; /* H */

	  for (i=0; i<qr->nr_qhop_residues; i++)
	    {
	      if (qr->qhop_residues[i].rtype = rt)
		{
		  for (tsite = 0;
		       tsite<qr->qhop_residues[i].nr_titrating_sites;
		       tsite++)
		    {
		      qa = qr->qhop_residues[i].titrating_sites[tsite];
		      if (strcmp(qr->qhop_atoms[qa].atomname,
				 *(mtop->moltype[mt].atoms.atomname[a1])) == 0)
			{
			  attach_proton(&(qr->qhop_atoms[qa]), a2);
			  attach_proton(&(qr->qhop_atoms[qa]), a3);
			}
		    }
		}
	    }
	  break;

	case F_VSITE2:
	  /* Can't see when this vsite type would be applicable,
	   * but let's implement it anyway. */
	  a1 = atom[1]; /* H */
	  if (mtop->moltype[mt].atoms.atom[a1].atomnumber = 1)
	    {
	      /* is it connected to a t_qhop_atom? */
	      for (i=0; i<qr->nr_qhop_residues; i++)
		{
		  if (qr->qhop_residues[i].rtype = rt)
		    {
		      for (tsite = 0;
			   tsite<qr->qhop_residues[i].nr_titrating_sites;
			   tsite++)
			{
			  for (a=2; a<=3; a++)
			    {
			      a2 = atom[a];
			      qa = qr->qhop_residues[i].titrating_sites[tsite];
			      if (strcmp(qr->qhop_atoms[qa].atomname,
					 *(mtop->moltype[mt].atoms.atomname[a2])) == 0)
				{
				  attach_proton(&(qr->qhop_atoms[qa]), a2);
				  break;
				}
			    }
			}
		    }
		}
	    }
	  break;

	case F_VSITE3:
	case F_VSITE3FD:
	case F_VSITE3FAD:
	case F_VSITE3OUT: /* <- Not entirely sure this works with an amine */
	case F_VSITE4FD:
	  a1 = atom[1]; /* H */
	  if (mtop->moltype[mt].atoms.atom[a1].atomnumber = 1)
	    {
	      /* is it connected to a t_qhop_atom? */
	      for (i=0; i<qr->nr_qhop_residues; i++)
		{
		  if (qr->qhop_residues[i].rtype = rt)
		    {
		      for (tsite = 0;
			   tsite<qr->qhop_residues[i].nr_titrating_sites;
			   tsite++)
			{
			  a2 = atom[2];
			  qa = qr->qhop_residues[i].titrating_sites[tsite];
			}
		    }
		}
	    }

	  if (db->H_map.atomid2H[atom[1]] >= 0)
	    {
	      /* Second atom is attached to the H. */
	      tsite = is_qhopatom(qr, atom[2]);
	      if (tsite != -1)
		{
		  attach_proton(&(qr->qhop_atoms[tsite]), atom[1]);
		  break;
		}
	    }
	  break;

	default:
	  fprintf(stderr, "Unsupported connection type for qhop.\n");
	}

    }
  while(natoms = ilist_iterator(ilist, itype, &i, atom) != 0);

#else
  /* Here's the simple scan */

  gmx_bool bNext;
  int  a, aa, r, qa;
  t_qhop_atom *qatoms;
  t_qhop_residue *qresidues;

  qatoms = qr->qhop_atoms;
  qresidues = qr->qhop_residues;

  for (qa=0; qa < qr->nr_qhop_atoms; qa++)
    {
      bNext = FALSE;
      r = qatoms[qa].qres_id;

      for (a=0; a < qresidues[r].nr_atoms && !bNext; a++)
	{
	  if (qresidues[r].atoms[a] == qatoms[qa].atom_id)
	    {
	      /* Atart scanning for hydrogens. Attach all H found
	       * to the qhop_atom qa and break as soon as a heavy
	       * atom is encountered.
	       * This assumes that all hydrogens, and the hydrogens only,
	       * have names starting with 'H', and that hydrogens come
	       * immediately after the heavy atoms they're attached to. */

	      for (aa = a+1; aa < qresidues[r].nr_atoms; aa++)
		{
		  if (qresidues[r].atomnames[aa] != NULL &&
		      qresidues[r].atomnames[aa][0] == 'H')
		    {
		      attach_proton(&(qatoms[qa]), qresidues[r].atoms[aa]);
		    }
		  else
		    {
		      bNext = TRUE;
		      break;
		    }
		}
	    }
	}
    }
    
#endif 

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
  int
    rb, r, a,
    reac, nreac[2], AD,
    **H, nH,
    q_atoms_nr,   q_atoms_max,
    q_residue_nr, q_residue_max,
    i, j, k, o, ngrps, resnr, mt;

  gmx_bool
    bMatchRB, bMatch;

  char
    *resname,
    *atomname;

  t_qhop_atom
    *q_atoms,
    *qa;

  t_qhop_residue          *q_residue;
  t_restp                 *rtp;
  qhop_reactant           *qreac[2];
  gmx_groups_t            *groups;
  gmx_mtop_atomloop_all_t aloop;
  t_atom                  *atom;
  t_ilist                 *ilist;
  q_atoms_nr = 0;
  q_residue_nr = 0;
  q_atoms_max = 1000;
  q_residue_max = 1000;

  snew(q_atoms, q_atoms_max);
  snew(q_residue, q_residue_max);
    
  groups = &mtop->groups;

  ngrps = ir->opts.ngqhopH;

  fprintf(stderr, "ir->opts.ngqhopH = %d", ir->opts.ngqhopH);

  aloop = gmx_mtop_atomloop_all_init(mtop);

  while (gmx_mtop_atomloop_all_next(aloop, &i, &atom)) {


    /* ***************************
       Set the H-existence array */

    /* group loop */
    for(j=0; j<ngrps; j++)
      {
	if (ggrpnr(groups, egcqhopH ,i) == j)
	  {
	    qdb->H_map.H[qdb->H_map.atomid2H[i]]        = 1;
	  }
      }

    /* *********************
       Find the qhop atoms */

    /* Do the resname and atomname show up in the qhop_resblocks? */

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

    gmx_mtop_atomloop_all_names(aloop, &atomname, &resnr, &resname);

    bMatchRB = FALSE;
    bMatch = FALSE;
    qhoprec->global_atom_to_qhop_atom[i] = NOTSET;

    for (rb=0; rb < qdb->rb.nrestypes && !bMatchRB; rb++)
      {
	if (strcmp(qdb->rb.restype[rb], resname) == 0)
	  {
	    /* We have a matching res(block) */
	    bMatchRB = TRUE;

	    for (r=0; r < qdb->rb.nres[rb] && !bMatch; r++)
	      {
		/* Is the atom found among the donors/acceptors ?*/
		nreac[0] = qdb->rb.res[rb][r].na;
		nreac[1] = qdb->rb.res[rb][r].nd;
		qreac[0] = qdb->rb.res[rb][r].acc;
		qreac[1] = qdb->rb.res[rb][r].don;

		/* AD==0 => acceptors, 1 => donors. This enables looping over both. */
		for (AD=0;
		     !bMatch && AD<2;
		     AD++)		  
		  {
	      
		    for (reac = 0;
			 !bMatch && (reac < nreac[AD]);
			 reac++)
		      {
			for (a=0;
			     !bMatch && a<qreac[AD][reac].nname;
			     a++)
			  {
			    if (strcmp(qreac[AD][reac].name[a], atomname) == 0)
			      {
				bMatch = TRUE;
				q_atoms[q_atoms_nr].resname      = qdb->rb.restype[rb];
				q_atoms[q_atoms_nr].atomname     = qreac[AD][reac].name[a];
				q_atoms[q_atoms_nr].res_id       = resnr-1;
				q_atoms[q_atoms_nr].atom_id      = i;
				q_atoms[q_atoms_nr].nr_protons   = 0; /* To be decided */
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
				    q_residue[q_residue_nr].atomnames          = NULL;     /* To be set. */
				    q_residue[q_residue_nr].nr_atoms           = 0;     /* To be decided. */
				    q_residue[q_residue_nr].nr_titrating_sites = 0; /* To be decided. */
				    q_residue[q_residue_nr].res_nr             = NOTSET;  /* To be set */
				    q_residue[q_residue_nr].res_nr             = resnr-1;

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

  fprintf(stderr,"nr qhop donors = %d\n", q_atoms_nr);

  /* Assign atoms to each qhop_residue */
  aloop = gmx_mtop_atomloop_all_init(mtop);

  while (gmx_mtop_atomloop_all_next(aloop, &i, &atom))
    {
      gmx_mtop_atomloop_all_names(aloop, &atomname, &resnr, &resname);

      for (r=0; r < q_residue_nr; r++)
	{
	  if (resnr-1 == q_residue[r].res_nr)
	    {
	      rtp = &(qdb->rtp[q_residue[r].rtype]);
	      
	      srenew(q_residue[r].atoms,     q_residue[r].nr_atoms + 1);
	      srenew(q_residue[r].atomnames, q_residue[r].nr_atoms + 1);

	      q_residue[r].atoms    [q_residue[r].nr_atoms] = i;
	      q_residue[r].atomnames[q_residue[r].nr_atoms] = strdup(atomname);
	      q_residue[r].nr_atoms++;

	      /* The following block was made redundant by the call to scan_for_H in qhop_init() */

/* 	      /\* Is it a hydrogen? *\/ */
/* 	      if (atom->atomnumber == 1) */
/* 		{ */
/* 		  /\* Find the right titrating site and add the H to the proton array. *\/ */
/* 		  bMatch = FALSE; */

/* 		  for (j=0; j < q_residue[r].nr_titrating_sites && !bMatch; j++) */
/* 		    { */
/* 		      qa = &(q_atoms[q_residue[r].titrating_sites[j]]); */

/* 		      /\* Loop over the t_restp to see which donor/acceptor it is bound to *\/ */
/* 		      for (k=0; k < rtp->rb[ebtsBONDS].nb && !bMatch; k++) */
/* 			{ */
/* 			  for (a=0; a<2; a++) */
/* 			    { */
/* 			      if (strcmp(rtp->rb[ebtsBONDS].b[k].a[a%2], qa->atomname) == 0 && */
/* 				  strcmp(rtp->rb[ebtsBONDS].b[k].a[(a+1)%2], atomname) == 0) */
/* 				{ */
/* 				  srenew(qa->protons, qa->nr_protons+1); */
/* 				  qa->protons[qa->nr_protons] = i; */
/* 				  qa->nr_protons++; */

/* 				  bMatch = TRUE; */
/* 				  break; */
/* 				} */
/* 			    } */
/* 			} */
/* 		    } */
/* 		} */
	    }
	}
    }  

  /* Clean up */

  srenew(q_atoms, q_atoms_nr);     /* size it down, as it is most likely oversized. */
  srenew(q_residue, q_residue_nr); /* size it down, as it is most likely oversized. */
  qhoprec->qhop_atoms       = q_atoms;
  qhoprec->qhop_residues    = q_residue;
  qhoprec->nr_qhop_atoms    = q_atoms_nr;
  qhoprec->nr_qhop_residues = q_residue_nr;


  /* Scan the idef for the hydrogens we may have missed. */

  /* Tried gmx_mtop_ilistloop_*, but I couldn't get it to work.*/


  /* loop over moltypes */
/*   for (mt=0; mt < mtop->nmoltype; mt++) */
/*     { */
/*       /\* loop over residues in molecule *\/ */
/*       for (r=0; r < mtop->moltype[mt].atoms.nres; r++) */
/* 	{ */
/* 	  rt = is_restype(qr, *(mtop->moltype[mt].atoms.resinfo[r].name)); */

/* 	  if (rt >= 0) */
/* 	    { */
/* 	      for(a=0; a < mtop->moltype[mt].atoms.nr; a++) */
/* 		{ */
/* 		  if (mtop->moltype[mt].atoms.atom[a].res == r) */
/* 		    { */
		     
/* 		    } */
/* 		} */
	      /* Try the g_hbond way. */

/* 	      for (i=0; i<F_NRE; i++) */
/* 		{ */
/* 		  if (i == F_POSRES) */
/* 		    { */
/* 		      continue; */
/* 		    } */
		  
/* 		  for (j=0; */
/* 		       j < ilist[i]->nr; */
/* 		       j+= interaction_function[i].nratoms+1) */
/* 		    { */
/* 		      if (i == F_SETTLE) */
/* 			{ */
/* 			  a1 = ilist[i]->iatoms[j+1]; */
/* 			  a2 = ilist[i]->iatoms[j+2]; */
/* 			  a3 = ilist[i]->iatoms[j+3]; */
/* 			  for (k=0; k < qhoprec->nr_qhop_residues; k++) */
/* 			    { */
/* 			      if (qhoprec->qhop_residues[k] == rt) */
/* 				{ */
/* 				  for (a=0; a < qhoprec->qhop_residues[k].nr_titratable_sites; a++) */
/* 				    { */
/* 				      qa = qhoprec->qhop_residues[k].titratable_sites[a]; */
/* 				      if (strcmp(qr->qhop_atoms[qa], */
/* 						 *(mtop->moltype[mt].t_atoms.atomname[a1])) == 0) */
/* 					{ */
/* 					  attach_proton(&(qhoprec->qhop_atoms[qa]), a2); */
/* 					  attach_proton(&(qhoprec->qhop_atoms[qa]), a3); */
/* 					} */
/* 				    } */
/* 				} */
/* 			    } */
/* 			  continue; */
/* 			} */
		      
/* 		      if (IS_CHEMBOND(i)) */
/* 			{ */
/* 			  for (o=0; o<2; o++) */
/* 			    { */
/* 			      a1 = ilist[i]->iatoms[j+1+o]; */
/* 			      a2 = ilist[i]->iatoms[j+1+((o+1)%2)]; */
/* 			      if (mtop->moltype[mt].t_atoms.atom[a1].atomnumber = 1) */
/* 				{ */
/* 				  for (k=0; k < qhoprec->nr_qhop_residues; k++) */
/* 				    { */
/* 				      if (qhoprec->qhop_residues[k] == rt) */
/* 					{ */
/* 					  for (a=0; a < qhoprec->qhop_residues[k].nr_titratable_sites; a++) */
/* 					    { */
/* 					      qa = qhoprec->qhop_residues[k].titratable_sites[a]; */
/* 					      if (strcmp(qr->qhop_atoms[qa], */
/* 							 *(mtop->moltype[mt].t_atoms.atomname[a2])) == 0) */
/* 						{ */
/* 						  attach_proton(&(qhoprec->qhop_atoms[qa]), a1); */
/* 						} */
/* 					    } */
/* 					} */
/* 				    } */
/* 				  break; */
/* 				} */
/* 			    } */
/* 			  continue; */
/* 			} */

/* 		      if (IS_VSITE(i)) */
/* 			{ */
/* 			  a1 = ilist[i]->iatoms[j+1]; */
/* 			  /\* Now step backwards until a non-hydrogen is found.*\/ */
/* 			  a2 = a1-1; */
			  
/* 			} */
/* 		    } */
/* 		} */

/* 	      for (i=0; i<F_NRE; i++) */
/* 		{ */
/* 		  if (IS_CHEMBOND(i) */
/* 		      || i == F_SETTLE */
/* 		      || (i >= F_VSITE2 && i <= F_VSITE4FD)) */
/* 		    { */
/* 		      scan_ilist_for_hydrogens(mtop, qhoprec, &(qdb->H_map), rt, ilist[i], i); */
/* 		    } */
/* 		} */
/* 	    } */
/* 	} */
/*     } */
  return(q_atoms_nr);
} /* get_qhop_atoms */


/* Reads the info in the db->qhop_H_exist and finds the corresponding
 * residue subtype in the qhop_db.
 * Returns the index in db->res[rt][]
 * This function is used to set the correct protonation states
 * at the start of the simulation.
 * Right now matching is based on counting the protons on
 * titrating sites in the qhop_res and in the rtp. Thus,
 * the existence map may need slight reshuffling to really
 * represent the global protonation state. */
static int which_subRes(const gmx_mtop_t *top, const t_qhoprec *qr,
			qhop_db *db, const int resnr)
{
  int
    h, i, j, k, t, nH, *nHres, a, r, b, aa,
    atomid, nres, natoms, qrtp;
  char
    *Hname, *DAname;
  t_atom
    atom;
  t_qhop_residue
    *qres;
  gmx_bool
    bSameRes, bNew, bMatch;
  char
    **DAlist, **DAlist2;
  int
    *n_H, *n_H2, nDA, nDA2, DA;
  qhop_reactant
    *reac;
  t_restp
    *rtp;

  qres = &(qr->qhop_residues[resnr]);
  nres = db->rb.nres[qres->rtype];


  /* To better handle the tautomerism, maybe one should
     count hydrogens on the donors/acceptors instead
     of matching the names. Here we match the names.
  */

  nDA = qres->nr_titrating_sites;
  snew(n_H, nDA);
  snew(DAlist, nDA);

  /* ************************************************
     Count the protons present on each titrating site
     (qhop_atom) for this qhop_residue.
  */

  /* loop over donors/acceptors in the residue */
  for (i=0; i < nDA; i++)
    {
      k = qres->titrating_sites[i];

      /* Stash the donor/acceptor name */
      DAlist[i] = qr->qhop_atoms[k].atomname;

      /* Proton roll call */
      for (j=0; j < qr->qhop_atoms[k].nr_protons; j++)
	{
	  /* get the atom id of the protons */
	  h = qr->qhop_atoms[k].protons[j];

	  /* is it present? */
	  if (get_proton_presence(&(db->H_map), h))
	    {
	      n_H[i]++;
	    }
	}
    }

  /* Ok. n_H now stores the number of active protons on all donors/acceptors.
     Donor/acceptor names are stored in DAlist. */


  /* ********************************************************
     Count the protons on each titrating site for different
     residue subtypes and see if any matches the qhop_residue
  */

  /* Loop over residue subtypes */
  bSameRes = FALSE;
  for (r=0; r<nres && !bSameRes; r++)
    {
      n_H2 = NULL;
      DAlist2 = NULL;
      nDA2 = 0;

      bMatch = TRUE;

      /* DA==0 => acceptors, DA==1 => donors */
      for (DA=0; DA<2 && bMatch; DA++)
	{
	  /* Go through the donors and acceptors and make a list */
	  for (j=0;
	       (j < (DA==0 ?
		     db->rb.res[qres->rtype][r].na :
		     db->rb.res[qres->rtype][r].nd))
		 && bMatch;
	       j++)
	    {
	      reac = (DA==0 ?
		      &(db->rb.res[qres->rtype][r].acc[j]) :
		      &(db->rb.res[qres->rtype][r].don[j]));

	      /* Tautomer loop. There may be chemically equivalent
	       * sites, e.g. the oxygens in a carboxylate. */
	      for (t=0; t<reac->nname && bMatch; t++)
		{
		  /* Have we already seen this donor/acceptor? */
		  bNew = TRUE;
		  for (k=0; (k < nDA2) && bNew; k++)
		    {
		      if (strcmp(DAlist2[k], reac->name[t]) == 0)
			{
			  /* We already had it in the list */
			  bNew = FALSE;
			}
		    }
		  
		  if (bNew)
		    {
		      srenew(DAlist2, nDA2+1);
		      srenew(n_H2, nDA2+1);

		      n_H2[nDA2] = 0;
		      DAlist2[nDA2] = reac->name[t];
		      
		      rtp = &(db->rtp[db->rb.res[qres->rtype][r].rtp]);
		      for (a=0; a<rtp->natom; a++)
			{
			  if (strcmp(reac->name[t], *(rtp->atomname[a])) == 0)
			    {
			      /* Here's the titrating site. Find hydrogens.
			       * Assume a contiguous strtch of hydrogens followed
			       * by the next heavy atom. */
			      for (aa=a+1; aa<rtp->natom; aa++)
				{
				  if (*(rtp->atomname[aa]) && *(rtp->atomname[aa][0]) == 'H')
				    {
				      /* A(nother) hydrogen. */
				      n_H2[nDA2]++;
				    }
				  else
				    {
				      /* Next heavy atom. */
				      break;
				    }
				}
			      break;
			    }
			}


		      /* How do we know the number of protons? We'll have to scan the rtp. */
/* 		      for (b=0; b < db->rtp[reac->rtp].rb[ebtsBONDS].nb; b++) */
/* 			{ */
/* 			  rtp = &(db->rtp[reac->rtp]); */
/* 			  n_H2[nDA2] = 0; */
/* 			  for (a=0; a<2; a++) */
/* 			    { */
/* 			      /\* If the bond is between the donor and a H* then we ++ the hydrogen count *\/ */
/* 			      DAname = rtp->rb[ebtsBONDS].b[b].a[a%2]; */
/* 			      Hname  = rtp->rb[ebtsBONDS].b[b].a[(a+1)%2]; */

/* 			      if ((strcmp(DAname, DAlist2[nDA2]) == 0) && */
/* 				  (strncmp(Hname, "H", 1) == 0)) */
/* 				{ */
/* 				  /\* Yes, this H was connected to this donor/acceptor *\/ */
/* 				  n_H2[nDA2]++; */
/* 				} */
/* 			    } */
/* 			} */


/* 		      /\* Is the proton count the same for the titrating_sites[] in qhop_residue? *\/ */
/* 		      bMatch = FALSE; */
/* 		      for (i=0; i<nDA; i++) */
/* 			{ */
/* 			  if (strcmp(DAlist2[nDA2], DAlist[i]) == 0) */
/* 			    { */
/* 			      /\* Same donor/acceptor *\/ */
/* 			      if (n_H2[nDA2] == n_H[i]) */
/* 				{ */
/* 				  /\* Same proton count *\/ */
/* 				  bMatch = TRUE; */
/* 				} */
/* 			    } */
/* 			} */

/* 		      if (!bMatch) */
/* 			{ */
/* 			  break; */
/* 			} */

		      nDA2++;
		    }
		}
	    }
	}

      /* Now compare the hydrogen counts to see if this is the flavor we want. */
      /* First check if the total number of hydrogens are the same for both lists */
      nH = 0;
      for (j=0; j<nDA; j++)
	{
	  nH += n_H[j];
	}

      for (j=0; j<nDA2; j++)
	{
	  nH -= n_H2[j];
	}
      
      bSameRes = (nH==0); /* Same total number of hydrogens? */

      /* Sum up hydrogens on tautomeric sites and compare */

      /* DA==0 => acceptors, DA==1 => donors */
      for (DA=0; DA<2 && bSameRes; DA++)
	{
	  for (j=0;
	       (j < (DA==0 ?
		     db->rb.res[qres->rtype][r].na :
		     db->rb.res[qres->rtype][r].nd))
		 && bSameRes;
	       j++)
	    {
	      reac = (DA==0 ?
		      &(db->rb.res[qres->rtype][r].acc[j]) :
		      &(db->rb.res[qres->rtype][r].don[j]));

	      nH = 0;
	      
	      for (t=0; t<reac->nname && bSameRes; t++)
		{
		  for (i=0; i<nDA; i++)
		    {
		      if (strcmp(DAlist[i], reac->name[t]) == 0)
			{
			  nH += n_H[i];
			  break;
			}
		    }

		  for (i=0; i<nDA2; i++)
		    {
		      if (strcmp(DAlist2[i], reac->name[t]) == 0)
			{
			  nH -= n_H2[i];
			  break;
			}
		    }

		  /* If they are the same, then nH should be zero here. */
		  bSameRes = (nH==0);
		}
	    }
	}

      if (bSameRes)
	{
	  /* Adjust the existence map */
	  
	  /* Zero all hydrogens */
	  for (i=0; i < qres->nr_titrating_sites; i++)
	    {
	      for (j=0; j < qr->qhop_atoms[qres->titrating_sites[i]].nr_protons; j++)
		{
		  set_proton_presence(&(db->H_map),
				      qr->qhop_atoms[qres->titrating_sites[i]].protons[j],
					 FALSE);
		}
	    }
	  
	  /* Go through the rtp */
	  rtp = &(db->rtp[db->rb.res[qres->rtype][r].rtp]);
	  
	  for (i=0; i < rtp->natom; i++)
	    {
	      /* Is it a hydrogen atom? */
	      if (rtp->atom[i].atomnumber == 1)
		{
		  /* Match with atoms in qres */
		  for (j=0; j < qres->nr_atoms; j++)
		    {
		      if (strcmp(*(rtp->atomname[i]), qres->atomnames[j]) == 0)
			{
			  set_proton_presence(&(db->H_map),
					      qres->atoms[j],
					      TRUE);
			}
		    }
		}
	    }
	}


      /* Clean up */
      if (n_H2 != NULL)
	{
	  sfree(n_H2);
	}
      
      if (DAlist2 != NULL)
	{
	  sfree(DAlist2);
	}
      
      if (bSameRes)
	{
	  break;
	}
    }

  if (r >= nres)
    {
      gmx_fatal(FARGS, "Didn't find the subres.");
    }

  return r;

 
  /* Here's the other method, matching H-names exactly to the
     residue subtypes. Makes gromacs a bit too strict perhaps. */

  /* Make a list of all hydrogens for this restype */
/*   qrtp = db->rb.rtp[qres->rtype]; /\* The rtp data *\/ */
/*   natoms = db->rtp[qrtp].natom; */

/*   for (i=0, nH=0; */
/*        i < natoms; */
/*        i++) */
/*     { */
/*       if (db->rtp[qrtp].atom[i].atomnumber == 1) */
/* 	{ */
/* 	  srenew(Hlist, nH+1); */
/* 	  Hlist[nH] = *(db->rtp[qrtp].atomname[i]); */
/* 	  nH++; */
/* 	} */
/*     } */

  /* Loop over the residue subtypes and see which one matches */
  /* for (r=0; r<nres; r++) */
/*     { */
/*       /\* scan the rtp *\/ */
/*       for (j=0; j < db->nrtp; j++) */
/* 	{ */
/* 	  if (gmx_strcasecmp(db->rtp[j].resname, */
/* 			     db->rb.res[qres->rtype][r]) == 0) */
/* 	    { /\* Mathcing resnames! Now look at the protons. *\/ */
/* 	      bSameRes = TRUE; */
/* 	    } */
/* 	} */
/*     } */
}


static void qhop_swap_bondeds(t_qhop_residue *swapres, qhop_res *prod)
{
  /* placeholder. */
  int i;
  i=0;
  
}

static void qhop_swap_vdws(t_qhop_residue *swapres, qhop_res *prod)
{
  /* placeholder. */
  int i;
  i=0;
  
}

static void low_level_swap_m_and_q(t_mdatoms *md, t_atom *atom, const int atomid)
{
  real m, q;
  
  if (atom == NULL)
    {
      m = 1;
      q = 0;
    }
  else
    {
      m = atom->m;
      q = atom->q;
    }

  md->chargeA[atomid] = q;

/*   if (md->massA) */
/*     { */
/*       md->massA[atomid] = m; */
/*     } */

/*   md->invmass[atomid] = (m==0 ? 0 : 1/m); */
}

/* Hacks mdatoms such that the masses and
 * charges of swapres match those of prod.
 * Works on an atomname basis
 */
static void qhop_swap_m_and_q(const t_qhop_residue *swapres,
			      const qhop_res *prod,
			      t_mdatoms *md,
			      const qhop_db *db)
{
  int i, j;
  t_restp *rtp;
  real m;

  rtp = &(db->rtp[prod->rtp]);

  /* Compare atomnames with atomnames in rtp.
   * Change */

  for (i=0; i < swapres->nr_atoms; i++)
    {
      /* Zero m and q for all atoms
       * Otherwise any atoms that was
       * just turned off will retain
       * their m and q. */

      low_level_swap_m_and_q(md, NULL, swapres->atoms[i]);

      for (j=0; j < rtp->natom; j++)
	{
	  if (strcmp(*(rtp->atomname[j]), swapres->atomnames[i]) == 0)
	    {
	      /* Now set the q and m for real */
	      low_level_swap_m_and_q(md, &(rtp->atom[j]), swapres->atoms[i]);
	      break;
	    }
	}
    }
}

static void set_interactions(qhop_db *qdb, t_mdatoms *md, t_qhop_atom *QA, t_qhop_residue *qres)
{
  int i, j, k, nri, bt, b, Gr, RB, R;
  qhop_res *reac;
  t_restp *res;
  
  /* point pointers at the right residue */
  Gr   = qres->res_nr; /* global res id. */
  RB   = qres->rtype;
  R    = qres->res;
  reac = &(qdb->rb.res[RB][R]);

  for (i=0; i<qres->nr_titrating_sites; i++)
    {
      /* find among qdb->res[][].acc/.don */
      for (j=0; j<reac->na; j++)
	{
	  /* loop over tautomeric sites */
	  for (k=0; k<reac->acc[j].nname; k++)
	    {
	      if (strcmp(reac->acc[j].name[k], QA[qres->titrating_sites[i]].atomname) == 0)
		{
		  QA[qres->titrating_sites[i]].state |= eQACC;
		  break;
		}
	    }
	}

      /* find among qdb->res[][].acc/.don */
      for (j=0; j<reac->nd; j++)
	{
	  /* loop over tautomeric sites */
	  for (k=0; k<reac->don[j].nname; k++)
	    {
	      if (strcmp(reac->don[j].name[k], QA[qres->titrating_sites[i]].atomname) == 0)
		{
		  QA[qres->titrating_sites[i]].state |= eQDON;
		  break;
		}
	    }
	}
    }

  qhop_swap_bondeds(qres, reac);
  qhop_swap_vdws(qres, reac);
  qhop_swap_m_and_q(qres, reac, md, qdb);
  
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


/* Obsolete */
/* static void get_chargesets(t_commrec *cr, t_qhoprec *qhoprec){ */

/*   int  */
/*     i,j,k,l,states; */
/*   FILE */
/*     *in; */
/*   char */
/*     filename[4],file[8],line[300],atomname[3]; */

/*   for (i=0; i < qhoprec->nr_qhop_atoms; i++) */
/*     { */
/*       states = 1; */
/* /\*     for(j=0;j<qhoprec->qhop_atoms[i].nr_links+1;j++){ *\/ */
/* /\*       states*=2; *\/ */
/* /\*     } *\/ */
/* /\*     qhoprec->qhop_residues[i].max_state = states; *\/ */
/* /\*   for(i=0;i<qhoprec->nr_qhop_atoms;i+=(qhoprec->qhop_atoms[i].nr_links+1)){ *\/ */
/*     /\* each unique residue is visited once in this loop */
/*      *\/ */
/*     /\* open db file, GG style  */
/*      *\/ */

    
/*     strncpy(filename, qhoprec->qhop_atoms[i].resname, 3); */
/*     sprintf(file, "%3s.dat",qhoprec->qhop_atoms[i].resname); */
/*     //  fprintf(stderr,"trying to open %s for reading\n",file); */

/*     in = fopen (file,"r"); */
    
/*     /\* allocate the chargesets, one for each protonation state  */
/*      *\/  */
/* /\*     for(j=0;j<qhoprec->qhop_atoms[i].nr_links+1;j++){ *\/ */
/* /\*       snew(qhoprec->qhop_residues[i+j].charge_set, *\/ */
/* /\* 	   qhoprec->qhop_residues[i+j].max_state); *\/ */
/* /\*       for(k=0;k<qhoprec->qhop_residues[i+j].max_state;k++){ *\/ */
/* /\*       	snew(qhoprec->qhop_residues[i+j].charge_set[k], *\/ */
/* /\* 	     qhoprec->qhop_residues[i+j].nr_atoms); *\/ */
/* /\* 	/\\* read in the chargesets *\\/ *\/ */
/* /\* 	for(l=0;l<qhoprec->qhop_residues[i+j].nr_atoms;l++){ *\/ */
/* /\* 	  /\\* copy the charges from the qhop database SOMEHOW *\\/ *\/ */
/* /\* 	  fgets(line,200,in); *\/ */
/* /\* 	  sscanf(line,"%s%f",atomname, *\/ */
/* /\* 		 &qhoprec->qhop_residues[i+j].charge_set[k][l]); *\/ */
/* /\* 	} *\/ */
/* /\* 	fgets(line,200,in);/\\* read in a white line after every charges *\/ */
/* /\* 			      section *\\/ *\/ */
/* /\*       } *\/ */
/* /\*     } *\/ */
/*     fclose(in); */
/*   } */
/* } /\* get_chargesets *\/ */

/* Obsolete */
/* static void get_protonation_state(t_commrec *cr, t_qhoprec *qhoprec){ */
/*   /\* Here we can find out the state of the residue: We make use */
/*      of the order to decide what the state is.... We use a the bit */
/*      definition to decide on the protonation state.  In case of three */
/*      sites, there are 8 possibilites, and therefore 8 chargesets to be */
/*      retreived from file. */

/*      0 1 2 3 4 5 6 7  */
/*      --------------- */
/*      0 1 0 1 0 1 0 1 */
/*      0 0 1 1 0 0 1 1 */
/*      0 0 0 0 1 1 1 1 */
/*      --------------- */
     
/*      Thus, if the second titrating atom (in global atom numbering) is */
/*      protonated and the third one too, we take the seventh (6, C */
/*      counting) row of charges from the external qhop data base. */
/*      Anyway, something like this should work I hope. We need to make */
/*      sure there is a good manual chapter explaining this if someone */
/*      wants to add a new residue. */
/*   *\/ */


/*   int  */
/*     i,j,k,states; */
/*   int */
/*     nr_protons; */
/*   gmx_bool */
/*     *protonated; */
  
/*   /\* now we need to decide the protonation state of the residue and */
/*      retreive the charges. The qhop_atoms are sorted already. */
/*   *\/ */
/* /\*   for(i=0;i<qhoprec->nr_qhop_atoms;i+=(qhoprec->qhop_atoms[i].nr_links+1)){ *\/ */
/* /\*     snew(protonated,qhoprec->qhop_atoms[i].nr_links+1); *\/ */
/* /\*     nr_protons = 0; *\/ */
/* /\*     /\\* for(j=0;j<qhoprec->qhop_atoms[i].nr_links+1;j++){ *\\/ *\/ */
/* /\* /\\*       protonated[j]=qhoprec->qhop_atoms[i+j].bdonor; *\\/ *\/ */
/* /\* /\\*       nr_protons+=qhoprec->qhop_atoms[i+j].bdonor; *\\/ *\/ */
/* /\* /\\*     } *\\/ *\/ */
/* /\*     for(j=0;j<qhoprec->qhop_atoms[i].nr_links+1;j++){ *\/ */
/* /\*       srenew(qhoprec->qhop_residues[j+i].protonated, *\/ */
/* /\* 	   qhoprec->qhop_atoms[i].nr_links+1); *\/ */
/* /\*       for(k=0;k<qhoprec->qhop_atoms[i].nr_links+1;k++){ *\/ */
/* /\* 	qhoprec->qhop_residues[j+i].protonated[k]=protonated[k]; *\/ */
/* /\*       } *\/ */
/* /\*       /\\* compute the integer corresponding to the proton byte, or pryte  *\/ */
/* /\*        *\\/ *\/ */
/* /\*       qhoprec->qhop_residues[i+j].pryte =  *\/ */
/* /\* 	qhop_get_protonation_state(qhoprec->qhop_atoms[i].nr_links+1,  *\/ */
/* /\* 				       protonated); *\/ */
/* /\*     } *\/ */
/* /\*     /\\* with the pryte, we can decide what chargeset to use during the *\/ */
/* /\*        simulation: charge[pryte][.] *\/ */
/* /\*      *\\/  *\/ */
/* /\*     free(protonated); *\/ */
/* /\*   } *\/ */
/* } /\* get_protonation_state *\/ */
 

/* THIS ONE MOVED TO src/gmxlib/qhoprec.c */
/* t_qhoprec *mk_qhoprec(void) */
/* { */
/*   t_qhoprec *qr; */

/*   snew(qr,1); */
/*   memset((void *)qr, 0, sizeof(t_qhoprec)); /\* Just in case *\/ */
/*   return (qr); */
/* }  /\* mk_qhoprec *\/ */

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

/* Sets the bqhopdonor[] and bqhopacceptor[] arrays in a t_mdatoms. */
void qhop_atoms2md(t_mdatoms *md, t_qhoprec *qr)
{
  int i;
  t_qhop_atom *a;

  for (i=0; i < qr->nr_qhop_atoms; i++)
    {
      a = &(qr->qhop_atoms[i]);

      md->bqhopacceptor[a->atom_id] = (a->state & eQACC) != 0;
      md->bqhopdonor[a->atom_id]    = (a->state & eQDON) != 0;
    }
}

int init_qhop(t_commrec *cr, gmx_mtop_t *mtop, t_inputrec *ir, 
	      t_forcerec *fr, rvec *x,matrix box,t_mdatoms *md,
	      qhop_db **db){

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

  if ((*db = qhop_db_read("qamber99sb.ff", mtop, md)) == NULL)
    {
      gmx_fatal(FARGS,"Can not read qhop database information");
    }

  set_pbc_dd(&pbc,fr->ePBC,DOMAINDECOMP(cr) ? cr->dd : NULL, FALSE, box);
  
  qhoprec = fr->qhoprec;

  qhoprec->qhopfreq = ir->qhopfreq;
  qhoprec->qhopmode = ir->qhopmode;

  snew(qhoprec->global_atom_to_qhop_atom, mtop->natoms);

  nr_qhop_atoms = get_qhop_atoms(cr, mtop, qhoprec ,ir, *db);

  /* Find hydrogens and fill out the qr->qhop_atoms[].protons */
  scan_for_H(qhoprec, *db);

  for (i=0; i < qhoprec->nr_qhop_residues; i++)
    {
      /* What flavour of the residue shall we set it to? */
      target = which_subRes(mtop, qhoprec, *db, i);
      qhoprec->qhop_residues[i].res = target;
      set_interactions(*db, md, qhoprec->qhop_atoms, &(qhoprec->qhop_residues[i]));
    }

  /* Complete the t_mdatoms */
  qhop_atoms2md(md, qhoprec);

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

static gmx_bool is_DA(qhop_db *db, t_qhop_residue *qr, t_qhop_atom *qa, int prod, int DA)
{
  int i, j, rt;
  qhop_reactant *qres;

  if (DA != eQDON && DA != eQACC)
    {
      gmx_fatal(FARGS, "Called is_DA() with DA != eQDON || eQACC");
    }

  rt = qr[qa->qres_id].rtype;

  if ((DA == eQACC ?
       db->rb.res[rt][prod].na :
       db->rb.res[rt][prod].nd)
      > 1)
    {
      /* Loop over acceptors on this residue */
      for (i=0;
	   i < (DA == eQACC ?
		db->rb.res[rt][prod].na :
		db->rb.res[rt][prod].nd);
	   i++)
	{
	  /* Loop over tautomers */
	  qres =  (DA == eQACC) ?
	    &(db->rb.res[rt][prod].acc[i]) :
	    &(db->rb.res[rt][prod].don[i]);

	  /* See if an atom matches */
	  for (j=0; j < qres->nname; j++)
	    {
	      if (strcmp(qres->name[j], qa->atomname) == 0)
		{
		  return TRUE;
		}
	    }
	}
    }

  return FALSE;
}

/* is_donor() and is_acceptor() are wrappers for is_DA() */
static gmx_bool is_donor(qhop_db *db, t_qhop_residue *qr, t_qhop_atom *qa, int prod)
{return is_DA(db, qr, qa, prod, eQDON);}
static gmx_bool is_acceptor(qhop_db *db, t_qhop_residue *qr, t_qhop_atom *qa, int prod)
{return is_DA(db, qr, qa, prod, eQACC);}

/* Don't use qhop_titrate directly. Use qhop_protonate()
 * or qhop_deprotonate().
 *
 * Arguments are like qhop_(de)protonate(),
 * with the addition of DA which,
 * if DA==eQDON, makes a deprotonation and,
 * if DA==eQACC, a protonation.
 * Returns the subres being the product. */
static int qhop_titrate(qhop_db *db, t_qhoprec *qr,
			t_qhop_atom *qatom, t_mdatoms *md,
			gmx_bool bWater, gmx_bool bSwapBondeds, int DA)
{
  int
    i, rt, r, reactant, t;
  qhop_res
    *product_res;
  gmx_bool bDonor;

  if (DA != eQDON && DA != eQACC)
    {
      gmx_fatal(FARGS, "Called qhop_titrate() with DA != eQDON || eQACC");
    }
  
  bDonor = (DA & eQDON);

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
	  if (!strcasecmp(bDonor ?
			  db->rb.res[rt][r].don[reactant].name[t] :
			  db->rb.res[rt][r].acc[reactant].name[t],
			  qatom->atomname));
	    {
	      /* This is the reactant! */
	      product_res = bDonor ?
		db->rb.res[rt][r].don[reactant].productdata :
		db->rb.res[rt][r].acc[reactant].productdata;
	      break;
	    }
	}
    }

  if (product_res == NULL)
    {
      gmx_fatal(FARGS, "Could not find the product for this reaction.");
    }

  /* Find the index in db->rb.res[rt][] */
  for (i=0; i<db->rb.nres[rt]; i++)
    {
      if (strcmp(product_res->name,
		 db->rb.res[rt][i].name) == 0)
	{
	  break;
	}
    }

  if (i == db->rb.nres[rt])
    {
      gmx_fatal(FARGS, "Could not find the product for this reaction.");
    }

  /* Change interactions */
  /* For water we don't change any bonded or VdW. For now. */
  if (!bWater)
    {
      if (bSwapBondeds)
	{
	  qhop_swap_bondeds(&(qr->qhop_residues[qatom->qres_id]), product_res);
	}
      
      qhop_swap_vdws   (&(qr->qhop_residues[qatom->qres_id]), product_res);
    }

  qhop_swap_m_and_q(&(qr->qhop_residues[qatom->qres_id]), product_res, md, db);

  /* Set the residue subtype */
  qr->qhop_residues[qatom->qres_id].res = i;

  return i;
}

/* Wrappers for qhop_titrate() */
extern void qhop_protonate(qhop_db *db, t_qhoprec *qr, t_qhop_atom *qatom,
			   t_mdatoms *md, gmx_bool bWater, gmx_bool bSwapBondeds)
{
  int prod;
  gmx_bool isacc;

  prod = qhop_titrate(db, qr, qatom, md, bWater, bSwapBondeds, eQACC);
  
  /* Since it just got a proton it must become a donor. */
  md->bqhopdonor[qatom->atom_id] = TRUE;
  qatom->state |= eQDON;

  /* Is is still an acceptor? */
  isacc =  is_acceptor(db, qr->qhop_residues, qatom, prod);
  md->bqhopacceptor[qatom->atom_id] = isacc;

  if (isacc)
    {
      qatom->state |= eQACC;
    }
  else
    {
      qatom->state &= ~eQACC;
    }

}

extern void qhop_deprotonate(qhop_db *db, t_qhoprec *qr, t_qhop_atom *qatom,
			     t_mdatoms *md, gmx_bool bWater, gmx_bool bSwapBondeds)
{
  int prod;
  gmx_bool isdon;

  prod = qhop_titrate(db, qr, qatom, md, bWater, bSwapBondeds, eQDON);

  /* Since it just lost a proton it must become an acceptor. */
  md->bqhopacceptor[qatom->atom_id] = TRUE;
  qatom->state |= eQACC;

  /* Is it still a donor? */
  isdon = is_donor(db, qr->qhop_residues, qatom, prod);
  md->bqhopdonor[qatom->atom_id] = isdon;

  if (isdon)
    {
      qatom->state |= eQDON;
    }
  else
    {
      qatom->state &= ~eQDON;
    }

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
    /* 0 and 1 are ordinary water protons, 2 is the "out of plane" one
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
	   * timestep. This can be compensated/circumvented
	   * by making the coordinate change in the update
	   * routines.
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

  /* (de)protonate */
  if( ((donor_atom->state & eQDON) && (acceptor_atom->state & eQACC))
      || ((acceptor_atom->state & eQDON) && (donor_atom->state & eQACC)) )
    {
      if (bUndo)
	{
	  /* These functions aren't fully implemented yet. */
	  qhop_protonate  (db, qhoprec, donor_atom, md,    bDonWater, FALSE);
	  qhop_deprotonate(db, qhoprec, acceptor_atom, md, bAccWater, FALSE);
	}
      else
	{
	  /* These functions aren't fully implemented yet. */
	  qhop_protonate  (db, qhoprec, acceptor_atom, md, bAccWater, FALSE);
	  qhop_deprotonate(db, qhoprec, donor_atom,    md, bDonWater, FALSE);
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
  sfree(enerd);
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
      for(k=0; k < qhoprec->qhop_residues[i].nr_atoms; k++)
	{
	  fprintf(stderr,"%3d ",qhoprec->qhop_residues[i].atoms[k]);
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
  return 1-pow((1-p), (qhop_interval/0.01));
}

static real calc_K(qhop_parameters *p, t_hop *hop)
{return p->k_1 * exp(-p->k_2 * (hop->rda-0.23)) + p->k_3;}
static real calc_M(qhop_parameters *p, t_hop *hop)
{return p->m_1 * exp(-p->m_2 * (hop->rda-0.23)) + p->m_3;}
/**
 * \brief Calculates the SE-regime validity limit and stores it in hop->El.
 *
 * \param p    Pointer to hopping parameters
 * \param hop  Pointer to a t_hop */
static void compute_E12_left(qhop_parameters *p,
		      t_hop *hop
		      )
{
  
  real
    K, M;

  /* compute the value of E12 that is required for 
   * p_SE(rda,10fs) = 0.1. 
   */
  K = calc_K(p, hop);
  M = calc_M(p, hop);

  hop->El = -(atanh( 2*(0.1-0.5) )-M) / K;
}


/**
 * Calculates the TST-regime validity limit and stores it in hop->Er.
 */
static void compute_E12_right(qhop_parameters *p, ///< Pointer to hopping parameters
		       t_hop *hop         ///< The hop
		       )
{
  
  real
    S, T, V, Eb, f, g, A;

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


  /* It gets worse. The quantity which is supposed to be
   * a factor log(100) larger than the thermal energy is
   * not Eb, but Eb-hbo. This leads to the following equation:
   *  Eb = f * exp(-g *Eb) + A  (1),
   * where A := h+log(100). Basically one needs to take the
   * ground state vibration into account. Eq. (1) can be
   * solved with the Labert W function (a.k.a. the product
   * logarithm). Using that Eb, we can solve the quadratic
   * equation described at the top of this function definition.
   */

  A = p->h + (log(100) * BOLTZ * hop->T);
  g = p->g;
  f = p->f;
  Eb = (A * g + LambertW(exp(-A * g) * f * g)) / g;

    /* Eb = RGAS * Temp * log(100); */
  /* d = hop->rda - (p->t_A); */
  S = calc_S(p, hop);
  T = p->s_B;
  V = calc_V(p, hop);

  hop->Er = (-T + sqrt(T*T - 4*V*(S-Eb))) / (2*V);
  
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
 * \brief Calculates the TST transfer probability for a 10 fs time window and hbar omega/2.
 */
static real compute_rate_TST(qhop_parameters *p, ///< Pointer to qhop_parameters
		      t_hop *hop         ///< The hop
		      )
{
  real
    Q, R, kappa, half_hbar_omega, ETST, pTST, E_M, Eb, E12, rda, T;

  rda = hop->rda;
  E12 = hop->Eb;
  Eb = hop->Eb;
  T = hop->T;

  if(E12 > 0.0)
    {
      E_M  = Eb - E12;
    }
  else
    {
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
  ETST =-(Eb-half_hbar_omega)/(BOLTZ * T);
  pTST = (kappa*BOLTZMANN*T/(PLANCK1 * 1e12)) * exp(ETST) * 0.01;

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
  
  K = calc_K(p, hop);
  M = calc_M(p, hop);
  pSE = (0.5 * tanh(-K * hop->E12 + M) + 0.5);
  return(pSE);
} /* compute_prob_SE */


/**
 * Calculates the transfer probability for the intermediate regime over a 10 fs time window.
 */
real compute_rate_log(qhop_parameters *p, ///< Pointer to qhop_parameters
		      t_hop *hop         ///< The hop 
		      )
{
  
   /* we now compute the probability over a 10fs interval. We return the
   * probability per timestep.
   */
  real 
    rSE,rTST,rate;
  
  rSE  = compute_rate_SE(p, hop);
  rTST = compute_rate_TST(p, hop);
  rate = rSE * pow((rTST/rSE), (hop->E12 - hop->El)/(hop->Er - hop->El));
  /* See, the log-space interpolation translates to the expression in the first argument above.
   * No need to log/exp back and forth. */

  return(rate);
} /* compute_prob_log */

static real compute_E12_0(qhop_parameters *p, t_hop *hop)
{
  return
    p->alpha
    + p->beta  * hop->rda
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
			 t_hop *hop, t_qhoprec *qhoprec,t_pbc pbc,int step, qhop_db *db){
  
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
/*   p = get_qhop_params(db->rb.res[qhoprec->qhop_atoms[hop->donor_id].rtype] */
/* 		      [(qhoprec->qhop_atoms[hop->donor_id].res)].name, */
/* 		      db->rb.res[qhoprec->qhop_atoms[hop->acceptor_id].rtype] */
/* 		      [(qhoprec->qhop_atoms[hop->acceptor_id].res)].name, */
/* 		      db); */

  p = get_qhop_params(qhoprec, hop, db);
    
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
    
      hop->E12 = compute_E12_0(p, hop) + Eafter-Ebefore;

      compute_E12_right(p, hop);
    
      compute_E12_left (p, hop);
      hop->Eb        = compute_Eb(p,hop);
      /* fprintf(stderr,"E12 = %f\n",hop->E12); */

      

      r_TST = compute_rate_TST(p, hop);
      r_SE  = compute_rate_SE (p, hop);

      if(E12 > hop->Er)
	{
	  /* Classical TST regime */
	  hop->prob = poisson_prob(r_TST, (ir->delta_t * qhoprec->qhopfreq)); /* (1-pow((1-r_TST),fr->qhoprec->qhopfreq)); */
	  hop->regime = etQhopTST;
	}
      else 
	{
	  if (E12 < hop->El)
	    {
	      /* Schroedinger regime */
	      /* using volkhards probality propagation:
	       */
	      hop->prob = poisson_prob(r_SE, (ir->delta_t * qhoprec->qhopfreq));
	      hop->regime = etQhopSE;
	      /* fprintf(stderr, "TEMPORARY:       r_SE = %e\n", r_SE); */
	    }
	  else
	    {
	      r_log = compute_rate_log(p, hop);
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
  /* hop->E12 = E12; */
/*   hop->Eb  = Eb; */
  free(p);

  return(Eafter_tot-Ebefore_tot);
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
  gmx_bool
    bAgain,bHop;
  real 
    rnr;
  int 
    start_seed=0;
  t_qhoprec
    *qr;

  qr = fr->qhoprec;

  set_pbc_dd(&pbc,fr->ePBC,DOMAINDECOMP(cr) ? cr->dd : NULL,FALSE,state->box);
  
  if(qr->rng == NULL)
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
      qr->rng     = gmx_rng_init(12/*start_seed*/);
      qr->rng_int = gmx_rng_init(12/*start_seed*/);
    } 

  qr->hop = find_acceptors(cr, fr, state->x, pbc, &(qr->nr_hops), md);

  if(qr->nr_hops > 0)
    {
      for(i=0; i<qr->nr_hops; i++)
	{
	  qr->hop[i].T = T;
	  if (qr->qhopmode == etQhopModeOne)
	    {
	      i = (int) floor (nr_hops * (gmx_rng_uniform_uint32(qr->rng_int))/4294967296.0);
	    }
	
	  qr->hop[i].DE_MM = get_hop_prob(cr,ir,nrnb,wcycle,top,mtop,groups,state,md,
					  fcd,graph,fr,vsite,mu_tot,/*born,bBornRadii,*/
					  &(qr->hop[i]),
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
    
      qr->hop = scramble_hops(qr->hop, qr->nr_hops, qr->qhopmode, qr->rng);

      for(i=0; i <qr-> nr_hops; i++)
	{
	  if (qr->hop[i].regime != etQhopNONE)
	    {
	      rnr = gmx_rng_uniform_real(qr->rng); 
	      if(MASTER(cr)){
		/* some printing to the log file */
		fprintf(fplog,
			"\n%d. don %d acc %d. E12 = %18.12f, DE_MM = %18.12f, Eb = %18.12f, rda = %3.4f, prob. = %f, ran. = %f (%s)", i,
			fr->qhoprec->qhop_atoms[qr->hop[i].donor_id].res_id,
			fr->qhoprec->qhop_atoms[qr->hop[i].acceptor_id].res_id,
			qr->hop[i].E12, qr->hop[i].DE_MM, qr->hop[i].Eb, qr->hop[i].rda,
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
		      for (j = i+1; j < qr->nr_hops; j++)
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

/* paramStr is the string from the rtp file.
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


