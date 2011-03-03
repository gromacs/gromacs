#include <stdlib.h>
#include <string.h>
#include "smalloc.h"
#include "string2.h"
#include "gmx_fatal.h"
#include "futil.h"
#include "gmxfio.h"
#include "hackblock.h"
#include "resall.h"
#include "gpp_atomtype.h"
#include "types/gmx_qhop_types.h"
#include "types/idef.h"
#include "gmx_qhop_parm.h"

#define assign_str(dst,src)  if (NULL != src) { if (NULL != dst) *dst = strdup(src); } else { *dst = NULL; }
#define assign_scal(dst,src) if (NULL != dst) *dst = src

/* Return a new qhop_resblocks structure */
extern qhop_resblocks_t qhop_resblock_init()
{
  qhop_resblocks *rb;
  snew(rb,1);
  rb->nrestypes = 0;
  return rb;
}

/* Adds parameters for an interaction to an array of parameters, */
/* or finds the parameter in case it's already in the array. */
/* Returns the index in params where it was found or inserted. */

/* btype is ebtsBONDS, ebtsANGLES, ...
 * bflavour specifies the functype. See the npbon[] in gmx_qhop_parm.h */
static int qhop_add_iparam(t_iparams **params, t_functype **ft, int *nparams, t_iparams *interaction, int btype, int bflavor, int ftype)
{
  int p, np, i;
  t_iparams par;
  gmx_bool bMatch;

  np = *nparams;

  par = *interaction;
  
  if (np == 0)
    {
      snew((*params), 1);
      snew((*ft), 1);
      *nparams = 1;
      (*params)[0] = par; /* copy params */
      (*ft)[0] = ftype;
       return 0;
    }

  for (p=0; p < np; p++)
    {
      bMatch = TRUE;
      /* Assume that the params are type real.
       * This means that the following types of params are disallowed:
       *   tab, vsiten, disres, dihres, orires, and cmap */
      for (i=0; i < npbon[btype][bflavor-1][1]; i++)
	{
	  if (par.generic.buf[i] != (*params)[p].generic.buf[i])
	    {
	      bMatch = FALSE;

	      break;
	    }
	}

      if (bMatch)
	{
	  /* It was already in the list */
	  return p;
	}
    }

  /* It was not found, let's add it. */
  srenew((*params), np+1);      /* expand the paramarray*/
  srenew((*ft), np+1);          /* expand the functypearray */
  (*nparams)++;
  (*params)[np] = *interaction; /* copy params */
  (*ft)[np]     = ftype;
  return np;
}

/* Reads the next entry in the string of parameters s. */
static real read_next_param(const char *s, int *pos)
{
  int n, i, p;
  double d;
  real r;

  
  n = sscanf(&(s[*pos]), "%lg", &d);
  r = (real)d;
    
  if (n <= 0)
    {
      /* No parameters read */
      *pos = -2;
      return r;
    }

  /* Find next param in s */
  p = *pos;

  /* find next space */
  while (s[p] != '\0')
    {
      if (s[p] == ' ')
	{
	  break;
	}
      p++;
    }
  
  /* find next param */
  while (s[p] != '\0')
    {
      if (s[p] != ' ')
	{
	  /* Here starts the next param. */
	  *pos = p;
	  return r;
	}
      p++;
    }
  
  /* We've reached the end of the string */
  *pos = -1;
  return r;
}

static t_iparams str2iparams(const char *s, int bt, int ft, int *nread)
{
  int i, pos, np;
  t_iparams params;
  real p[MAXFORCEPARAM+1];

  if (s==NULL || s[0]=='\0')
    {
      gmx_fatal(FARGS, "Interaction parameters were empty.");
    }

  memset(&params, 0, sizeof(t_iparams));

  pos = 0;
  np = 0;

  /* Keep reading until end of string or to last parameter in string */
  while(pos >= 0 && np < MAXFORCEPARAM)
    {
      p[np++] = read_next_param(s, &pos);

      if (pos==-2)
	{
	  /* No parameters were read. */
	  p[--np]=0;
	}
    }

/*   if (bt==ebtsBONDS) */
/*     { */
/*       /\* make it a constraint. *\/ */
/*       np = 2; */
/*     } */

  /* copy to params */
  for (i=0; i<np; i++)
    {
      params.generic.buf[i] = p[i];
    }

  if (ft == 8)
    {
      /* It's a table */
      params.tab.table  = (int)p[0];
      params.tab.kA     = (int)p[1];
      params.tab.kB     = (np==3) ? (int)p[2] : (int)p[1];
    }

  if ((bt == ebtsPDIHS && (ft==1 || ft==9)) ||
      (bt == ebtsIDIHS && ft==4))
    {
      params.pdihs.phiA  = p[0];
      params.pdihs.cpA   = p[1];
      params.pdihs.mult  = (int)p[2];
      params.pdihs.phiB  = (np==5) ? p[3] : p[1];
      params.pdihs.cpB   = (np==5) ? p[4] : p[2];
    }

  *nread = np;
  return params;
}

/* Copies A to B parameters (when needed) */
static void postprocess_params(t_iparams *params, int nread, int bt, int ft)
{
  t_iparams p, op;
  int i;

  op = *params;
  p = op;

  switch (bt)
    {
    case ebtsBONDS:
      switch (ft)
	{
	case 1: /* bond */
	case 2: /* G96 bond */
	case 6: /* harmonic */
	  if (nread == 2)
	    {
	      p.harmonic.rB   = op.harmonic.rA;
	      p.harmonic.krB  = op.harmonic.krA;
	    }
	  break;

	case 10: /* restraint */
	  if (nread == 4)
	    {
	      p.restraint.lowB  = op.restraint.lowA;
	      p.restraint.up1B  = op.restraint.up1A;
	      p.restraint.up2B  = op.restraint.up2A;
	      p.restraint.kB    = op.restraint.kA;
	    }
	  break;
	  
	default:
	  break;
	}

    case ebtsANGLES:
      switch (ft)
	{
	case 1: /* angle */
	case 2: /* G96 angle */
	  if (nread == 2)
	    {
	      p.harmonic.rB   = op.harmonic.rA;
	      p.harmonic.krB  = op.harmonic.krA;
	    }
	  break;

	default:
	  break;
	}

    case ebtsPDIHS:
    case ebtsIDIHS:
      switch (ft)
	{
	case 2: /* improper */
	  if (nread == 2)
	    {
	      p.harmonic.rB   = op.harmonic.rA;
	      p.harmonic.krB  = op.harmonic.krA;
	    }
	  break;

	case 3: /* RB */
	  if (nread == 6)
	    {
	      for (i=0; i<NR_RBDIHS; i++)
		{
		  p.rbdihs.rbcB[i] = p.rbdihs.rbcA[i];
		}
	    }
	  break;

	case 5: /* Fourier */
	  if (nread == 4)
	    {
	      for (i=0; i<4; i++)
		{
		  p.rbdihs.rbcB[i] = p.rbdihs.rbcA[i];
		}
	    }
	  break;
	  
	default:
	  break;
	}

    }

  *params = p;
}

/* Makes a library of bonded interaction parameters from rtp-data. */
extern void make_ilib(qhop_db *db)
{
  int rt, r, bt, b, bf, np[ebtsNR], nptot, fi, i, nread;
  int btype[ebtsNR] = {
    ebtypeBOND,
    ebtypeANGLE,
    ebtypeDIHEDRAL,
    ebtypeDIHEDRAL,
    -1,
    -1
  }; /* must reflect the enums in hackblock.h */

  for (bt=0; bt<ebtsNR; bt++)
    {
      np[bt] = 0;
    }

  t_restp *rtp;
  t_iparams *ilib, **il, params;
  t_functype *ftlib, **ft;

  /* Start with an empty ilib */
  ilib  = NULL;
  ftlib = NULL;
  il = NULL;
  ft = NULL;

  snew(il, ebtsNR);
  snew(ft, ebtsNR);

  for (rt=0; rt < db->rb.nrestypes; rt++)
    {
      if ((db->rb.bWater[rt] && db->constrain!=0) || !db->rb.bInTop[rt])
	{
	  continue;
	}

      for (r=0; r < db->rb.nres[rt]; r++)
	{
	  rtp = &(db->rtp[db->rb.res[rt][r].irtp]);

	  for (bt=0; bt < ebtsNR; bt++)
	    {
	      for (b=0; b < rtp->rb[bt].nb; b++)
		{
		  if (btype[bt] == -1)
		    {
		      gmx_fatal(FARGS, "Bonded type not supported. bt = %i", bt);
		    }

		  params = str2iparams(rtp->rb[bt].b[b].s, bt, db->rb.ftype[bt], &nread);

		  /* Sometimes not all parameters are perturbable,
		   * so we may have to shuffle around the params
		   * a bit */
		  postprocess_params(&params, nread, bt, db->rb.ftype[bt]);

		  fi =  qhop_add_iparam(&(il[bt]), &(ft[bt]), &(np[bt]), &params, btype[bt], db->rb.ftype[bt], db->rb.btype[bt]);
		  db->rb.res[rt][r].findex[bt][b] = fi;
		}
	    }
	}
    }

  /* /\* Add ebtsNR null interactions *\/ */
/*   for (bt=0; bt < ebtsNR; bt++) */
/*     { */
/*       switch (bt) */
/* 	{ */
/* 	case ebtsBONDS: */
/* 	  btype = ebtypeBOND; */
	  
/* 	  break; */
/* 	case ebtsANGLES: */
/* 	  btype = ebtypeANGLE; */
/* 	  break; */
/* 	case ebtsPDIHS: */
/* 	case ebtsIDIHS: */
/* 	  btype = ebtypeDIHEDRAL; */
/* 	  break; */
/* 	default: */
/* 	  gmx_fatal(FARGS, "Bonded type not supported. bt = %i", bt); */
/* 	} */

/*       if (bt == ebtsBONDS) */
/* 	{ */
/* 	  /\* Since we're constraining bonds, skip this one. *\/ */
/* 	  continue; */
/* 	} */

/*       /\* All zeroes. *\/ */
/*       memset(&params, 0, sizeof(t_iparams)); */

/*       db->rb.ft_null qhop_add_iparam(&ilib, &np, &params, btype, db->rb.btype[bt]); */
/*     } */

/*   We really just need one null interaction. */
  memset(&params, 0, sizeof(t_iparams));

  /* Stitch together the ilib and ftlib while adding dummy interactions. */
  for (bt=0, nptot=0; bt < ebtsNR; bt++)
    {
      db->rb.inull[bt] = NOTSET;

      if (bt <= ebtsIDIHS)
	{
	  db->rb.inull[bt] = qhop_add_iparam(&(il[bt]), &(ft[bt]), &(np[bt]), &params, btype[bt], db->rb.ftype[bt], db->rb.btype[bt]) + nptot;

	  /* Must shift the findex, now that we know how many params we have for each btype */
	  for (rt=0; rt < db->rb.nrestypes; rt++)
	    {
	      if ((db->rb.bWater[rt] && db->constrain!=0) || !db->rb.bInTop[rt])
		{
		  continue;
		}
	      
	      for (r=0; r < db->rb.nres[rt]; r++)
		{
		  rtp = &(db->rtp[db->rb.res[rt][r].irtp]);

		  for (b=0; b < rtp->rb[bt].nb; b++)
		    {
		      db->rb.res[rt][r].findex[bt][b] += nptot;
		    }
		}
	    }


	  nptot += np[bt];

	  srenew(ilib,  nptot);
	  srenew(ftlib, nptot);

	  for (i=0; i<np[bt]; i++)
	    {
	      /* memcpy(&(ilib[nptot-np[bt]]),  &(il[bt]), np[bt]*sizeof(t_iparams)); */
	      /* memcpy(&(ftlib[nptot-np[bt]]), &(ft[bt]), np[bt]*sizeof(t_functype)); */
	      ilib[nptot - np[bt] + i]  = il[bt][i];
	      ftlib[nptot - np[bt] + i] = ft[bt][i];
	    }

	  sfree(il[bt]);
	  sfree(ft[bt]);
	}
    }


  db->rb.ilib  = ilib;
  db->rb.ftlib = ftlib;
  db->rb.ni = nptot;
}

#if 0
t_idef* qhop_build_interaction_lib(char *ff, qhop_db *qdb, gpp_atomtype_t atype)
{
  FILE *fp;
  /* Chop up qdb a bit.*/
  t_restp *rtp = qdb->rtp;
  int nfiles = qdb->rb.nf;
/*   char **files = qdb->rb.files; */
  qhop_res *qres;
  gmx_bool match;
  int rr, bt, bi, i, m, rb, r,
    *a[5], ft,  ni, np, read, nbparams=0;
  size_t lineLen=0;
  char *dot, format[128];
  char const nt[]="%s";
  char pt[MAXFORCEPARAM*2+1];
  char line[1024], **bparams, *s, bondedname[128], p[MAXFORCEPARAM];
#define ALL_BPARAMS &ft, &p[0],&p[1],&p[2],&p[3],&p[4],&p[5],&p[6],&p[7],&p[8],&p[9],&p[10],&p[11] /* p[0:MAXFORCEPARAM] */
  t_idef *idef;

  sprintf(pt,"%%i"); /* ftype */
  for (i=0; i<MAXFORCEPARAM; i++)
    sprintf(&(pt[i*2+2]),"%%f");

  snew(idef,1);
  if (strcmp(ff, "ffqamber99sb") != 0)
    gmx_fatal(FARGS, "Qhop will only work with ffqamber99sb for now.");

  /* We already have the rtp data that we need */
  
  /* rtp[].atom                 : Non-bonded stuff
   *                              ptype, atomnumber, and maybe a few others are not
   *                              set up as they should at this point.
   *
   * rtp[].rb[]                 : The bonded stuff
   * rtp[].rb[].b.a[]           : The atoms for one bonded interaction
   * rtp[].rb[].b.s             : If an interaction parameter set is defined in the rtp, 
   *                              this is where to find it. This is typically a string
   *                              #defined in ffXXXXbon.itp.
   *
   * For now, let's stick to the following:
   * - All bonded interactions for qhop-related residues have to be set in the rtp-file
   * - All such definitions are to be found in ffXXXXbon.itp.
   * This may change later on, so that all force field files are scanned for parmeters,
   * or this may all go into grompp, but this allows me to go on with other parts of
   * qhop with less effort. Maybe that allows us to ignore the itp-scanning below for the time being. */


  /* How do we store stuff? In a t_idef I guess,
   * either in the topology itself (t_topology),
   * or in a separate libry. Let's go for a t_idef
   * in the qhop_resblocks structure for now. */

  /* Find strings from the rtp file defined in the bonded parameter file. */

  fp = libopen("ffqamber99sbbon.itp");
  while (get_a_line(fp, &line, 1023)) /* kills comments and leading spaces. */
    {
#define BONDED_REC_MIN_LEN 7 /* This is the length of "#define" */
      read = strlen(line);
      if (read >= BONDED_REC_MIN_LEN) {
	if (strncmp(line, "#define", BONDED_REC_MIN_LEN) == 0) {
	  s = strpbrk(line, " \t"); /* go past the #define token*/
	  if (s != NULL) {
	    srenew(bparams, nbparams+1);
	    bparams[nbparams++] = strdup(s);
	    /* This stores the whole remainder of the line as a string */
	  }
	}
      }
#undef BONDED_REC_MIN_LEN
    }

  /* now, go through the bonded definitions */
  for (rr=0; rr<qdb->nrtp; rr++) { /* Different residues. */
    for (bt=0; bt<ebtsNR; bt++) { /* The different bonded types */
      for (bi=0; bi<rtp[rr].rb[bt].nb; bi++) { /* All bonded interactions of that type */
	if (rtp[rr].rb[bt].b[bi].s != NULL) {
	  for (m=0; m<nbparams; m++) { /* Loop over the macro definitions found in ffXXXbon.itp */
	    s = bparams[m];

	    /* Just match the name to begin with */
	    if (sscanf(s, " %s", bondedname) > 0) {
	      if (0 == strcmp(rtp[rr].rb[bt].b[bi].s, bondedname)) { /* We have a match */
		sprintf(format," %s","%s%s%s%s");
		memcpy((void *)(format+(+btsNiatoms[bt]*2+1)), (void *)pt, 17);
	  
		switch(btsNiatoms[bt])
		  {
		  case 2:
		    np = sscanf(s, format, bondedname,
				a[0], a[1], ALL_BPARAMS);
		    break;
		  case 3:
		    np = sscanf(s, format, bondedname,
				a[0], a[1], a[2], ALL_BPARAMS);
		    break;
		  case 4:
		    np = sscanf(s, format, bondedname,
				a[0], a[1], a[2], a[3], ALL_BPARAMS);
		    break;
		  case 5:
		    np = sscanf(s, format, bondedname,
				a[0], a[1], a[2], a[3], a[4], ALL_BPARAMS);
		    break;
		  default:
		    /* This reeeeally shouldn't happen, unless a new bonded type is added with more than 5 atoms. */
		    gmx_fatal(FARGS, "Wrong number of atoms requested: %d", btsNiatoms[bt]);
		  }
		/* Build the interaction */
		if (np > 0)
		  {
		    /* We're gonna make constraints outta every bond, at least for now. */
		    /* What about vsites? Not yet. */
		    srenew(qdb->rb.ilib, ++qdb->rb.ni);
		    switch(bt)
		      {
		      case ebtsBONDS:			/* bonds */
			qdb->rb.ilib[ni].constr.dA = p[0];
			break;
		      case ebtsANGLES:			/* angles */
			switch (ft)
			  {
			  case 1: /* harmonic */
			  case 2: /* G96 */
			    qdb->rb.ilib[ni].harmonic.rA = p[0];
			    qdb->rb.ilib[ni].harmonic.krA = p[1];
			    break;
			  case 3: /* Cross bond-bond */
			    qdb->rb.ilib[ni].cross_bb.r1e = p[0];
			    qdb->rb.ilib[ni].cross_bb.r2e = p[1];
			    qdb->rb.ilib[ni].cross_bb.krr = p[2];
			    break;
			  case 4: /* Cross bond-angle */
			    qdb->rb.ilib[ni].cross_ba.r1e = p[0];
			    qdb->rb.ilib[ni].cross_ba.r2e = p[1];
			    qdb->rb.ilib[ni].cross_ba.r3e = p[2];
			    qdb->rb.ilib[ni].cross_ba.krt = p[3];
			    break;
			  case 5: /* UB */
			    qdb->rb.ilib[ni].u_b.theta  = p[0];
			    qdb->rb.ilib[ni].u_b.ktheta = p[1];
			    qdb->rb.ilib[ni].u_b.r13    = p[2];
			    qdb->rb.ilib[ni].u_b.kUB    = p[3];
			    break;
			  case 6: /* qangle */
			    qdb->rb.ilib[ni].qangle.theta  = p[0];
			    for (i=0; i<6;i++)
			      qdb->rb.ilib[ni].qangle.c[i]   = p[i+1];
			    break;
			  case 8: /* table */
			    qdb->rb.ilib[ni].tab.table = (int)p[0];
			    qdb->rb.ilib[ni].tab.kA    = p[1];
			    qdb->rb.ilib[ni].tab.kB    = p[2];
			    break;
			  default:
			    gmx_fatal(FARGS, "Unknown angle type.");
			  }
			break;
		      case ebtsPDIHS:			/* proper dihedrals */
		      case ebtsIDIHS:			/* improper dihedrals */
			switch (ft)
			  {
			  case 1:
			    qdb->rb.ilib[ni].pdihs.phiA = p[0];
			    qdb->rb.ilib[ni].pdihs.cpA  = p[1];
			    qdb->rb.ilib[ni].pdihs.mult = (int)p[2];
			    break;
			  case 2: /* improper */
			    qdb->rb.ilib[ni].harmonic.rA  = p[0];
			    qdb->rb.ilib[ni].harmonic.krA = p[1];
			    break;
			  case 3: /* RB */
			    for (i=0; i<NR_RBDIHS;i++)
			      qdb->rb.ilib[ni].rbdihs.rbcA[i] = p[i];
			    break;
			  case 5: /* Fourier. Dunno what to do here yet. */
			    break;
			  case 8:
			    qdb->rb.ilib[ni].tab.table = (int)p[0];
			    qdb->rb.ilib[ni].tab.kA    = p[1];
			    qdb->rb.ilib[ni].tab.kB    = p[2];
			    break;
			  default:
			    gmx_fatal(FARGS, "Unknown %sdihedral.", bt==ebtsIDIHS ? "improper ":"");
			  }
			break;
		      default:
			gmx_fatal(FARGS, "In the current implementation of qhop %s are unsupported.", btsNames[bt]);
		      }
		    /* find the res in qdb->rb*/
		    match = FALSE;
		    for (rb=0; rb<qdb->rb.nrestypes && !match; rb++)
		      for (r=0; r<qdb->rb.nres[rb] && !match; r++)
			{
			  qres = &(qdb->rb.res[rb][r]);
			  match = (strcmp(rtp[rr].resname, qres->name)==0);
			  if (match)
			    {
			      srenew(qres->ft, qres->nft + 1);
			      qres->ft[qres->nft] = ni;
			    }
			}
		    if (!match)
		      gmx_fatal(FARGS, "Residue %s not found in the qhop database.", rtp[rr].resname);
		    ni++;
		  }
		else
		  {
		    gmx_fatal(FARGS, "No parameters read!");
		  }
	      }
	    }
	  }
	}
      }
    }
  }

  ffclose(fp);

  /* Add dummy interactions.
   * We need only two: one for bonds, with rA=0.1nm and one with zeroes all over. */

  qdb->rb.ni+=2;
  srenew(qdb->rb.ilib, qdb->rb.ni);
  memset(&(qdb->rb.ilib[qdb->rb.ni-2]), 0, sizeof(t_iparams)*2); /* zeroes all over */
  qdb->rb.ilib[qdb->rb.ni-2].harmonic.rA = 0.1; /* Dunno if this is needed, but it can sure be used for dummy constraints.
						 * question is if we need such? */


  /* -----------------------------------------.
   * Now scan extra itp-files                  \
   */
#if 0
  /* Let's ignore this for now, as everything is in the rtp right now. */
  for (i=0; i<nfiles; i++){
    /* determine filetype */
    if (dot = rindex(files[i],'.'))
      if (strcmp(".itp", dot) == 0)
	{
	  /* It's an itp file. */
	}
      else
	dot = NULL;

    if (!dot)
      gmx_fatal(FARGS, "QHOP: %s appears not to be an itp-file, by judging from its name.", files[i]);
  }
  /* itp-files scanned.                        /
   * -----------------------------------------'
   */
#endif
}
#endif

/* Add a res to a restype. Requires a set of residues collected in a qhop_res_t. */
extern void qhop_add_res(qhop_resblocks_t rb, int resblocknr, qhop_res_t res, int nres)
{
  int i,n;
  if (!rb->res)
    snew(rb->res, 1);
  else
    /* Let's not assume that we're adding stuff to a new resblock. */
    if (resblocknr > rb->nrestypes+1)
      srenew(rb->res,rb->nrestypes+1);
  
  /* The following allows for addition of residues to a non-empty resblock. */
  n = rb->nres[resblocknr];
  if (n > 0)
    n--;

  for (i=0; i<nres; i++) {
    rb->res[resblocknr][i + n] = res[i];
    rb->nres[resblocknr]++;
  }
}

/* Add a restype. Requires a set of residues collected in a qhop_res_t. */
extern void qhop_add_restype(qhop_resblocks_t rb, char *name, int nres, qhop_res_t res)
{
  int i, j;
  if (!rb)
    gmx_fatal(FARGS,"The resblock is not initialized!\n");

  /* (re)allocate restype */
  if (!rb->restype) {
    snew(rb->restype,1);
    i = 1;
  } else {
    i = rb->nrestypes+1;
    srenew(rb->restype, i);
  }

  /* (re)allocate nres */
  if (!rb->nres)
    snew(rb->nres, 1);
  else
    srenew(rb->nres, i);

  /*rb->nres[rb->nrestypes] = nres;*/

  qhop_add_res(rb, rb->nrestypes, res, nres);

  rb->nrestypes++;
}

static void set_H_exist(const qhop_H_exist *H_map, const atom_id H, const gmx_bool bON)
{
  H_map->H[H_map->atomid2H[H]] = (char) bON ? 1:0;
}


/* Redirects the interaction parameters to the ones dictated by the qres */
static void qhop_change_interactions(t_ilist *ilist, qhop_db *db, t_qhop_residue *qres)
{
  printf("Not implemented yet!\n");
}

/* H is global atom number */
extern void qhop_set_protonation(const qhop_db *db, t_qhop_residue *qres,
				 const atom_id H)
{
  const gmx_bool bON = db->H_map.H[db->H_map.atomid2H[H]] == 0; /* if it's zero, then the proton will appear */
  gmx_bool bNotFound = TRUE;
  int rt, Hloc, i, j, reactant;
  char *Hname;
  qhop_res *res, *product;
  qhop_reactant *qreac;
  const qhop_resblocks *rb = &(db->rb);
  
  qreac = NULL;
  j  =-1;
  rt = qres->rtype;
  
  if (db->H_map.atomid2H < 0)
    gmx_fatal(FARGS, "In qhop_set_protonation(): Atom %i is not a hydrogen!", H);

  /* Hloc is the local atom number of H */
  Hloc = H-qres->atoms[0];
  Hname = *(db->rtp[rt].atomname[Hloc]);
 
  /* What residue subtype will this result in? */
  res = &(rb->res[qres->rtype][qres->res]);
  for (reactant=0;
       reactant<bON ? res->na : res->nd;
       reactant++)
    {
      qreac = bON ? &(res->acc[reactant]) : &(res->don[reactant]);
      for (j=0; j<qreac->nH; j++)
	{
	  if (gmx_strcasecmp(Hname, qreac->H[j]) == 0)
	    {
	      bNotFound = FALSE;
	      break;
	    }
	}
    }
  if (bNotFound)
    gmx_fatal(FARGS, "In qhop_set_protonation(): "
	      "Hydrogen %s not found in residue %s!",
	      Hname, res->name);
  if (NULL == qreac)
    gmx_fatal(FARGS,"qreac was never initialized!");
  product = &(qreac->productdata[j]);

  /* Must find the index to this residue subtype
   * Reuse bNotFound */

  bNotFound = TRUE;
  for (j=0; j<rb->nres[rt]; j++)
    {
      if (gmx_strcasecmp(product->name, rb->res[rt][j].name) == 0)
	{
	  bNotFound = FALSE;
	  break;
	}
    }
     
  if (bNotFound)
    gmx_fatal(FARGS, "In qhop_set_protonation(): "
	      "Product %s not found among the residue subtypes of %s!",
	      product->name, rb->restype[rt]);

  set_H_exist(&db->H_map, H, bON);
  /* NOW SWAP INTERACTIONS */
}

extern int qhop_get_residue_subtype(const t_qhop_residue *qres)
{
  return qres->res;
}

/* Return a new qhop structure */
extern qhop_t qhop_init()
{
  struct qhop *qht;
  
  snew(qht,1);

  return qht;
}

extern void qhop_set_donor(qhop_t gqh, const char *donor)
{
  if (donor != NULL) {
    printf(" donor %s,", donor);
    gqh->donor = strdup(donor);
  }
}

extern void qhop_set_acceptor(qhop_t gqh, const char *acceptor)
{
  if (acceptor != NULL) {
    printf(" acceptor %s\n", acceptor);
    gqh->acceptor = strdup(acceptor);
  }
}

extern void qhop_set_don_atom(qhop_t gqh, const char *don_atom)
{
  if (don_atom != NULL)
    {
      printf("Setting don_atom %s\n", don_atom);
      gqh->don_atom = strdup(don_atom);
    }
}

extern void qhop_set_acc_atom(qhop_t gqh, const char *acc_atom)
{
  if (acc_atom != NULL)
    {
      printf("Setting acc_atom %s\n", acc_atom);
      gqh->acc_atom = strdup(acc_atom);
    }
}

extern char *qhop_get_don_atom(const qhop_t gqh)
{
  return gqh->don_atom;
}

extern char *qhop_get_acc_atom(const qhop_t gqh)
{
  return gqh->acc_atom;
}

extern char *qhop_get_donor(const qhop_t gqh)
{
  return gqh->donor;
}

extern char *qhop_get_acceptor(const qhop_t gqh)
{
  return gqh->acceptor;
}

/* Add parameter to gqh, return 1 if OK, 0 if not OK */
extern int qhop_add_param(qhop_t gqh,char *name,char *value,char *unit)
{
  srenew(gqh->name,gqh->nparam+1);
  srenew(gqh->value,gqh->nparam+1);
  srenew(gqh->unit,gqh->nparam+1);
  gqh->name[gqh->nparam]  = strdup(name);
  gqh->value[gqh->nparam] = strdup(value);
  gqh->unit[gqh->nparam]  = strdup(unit);
  gqh->nparam++;
  
  return 1;
}

/* Lists the parameters, one by one on repeatedly calling the
   function. Returns 1 if OK, 0 if not OK */
extern int qhop_get_param(qhop_t gqh,char **name,char **value,char **unit)
{
  if (gqh->nparam_c < gqh->nparam) {
    assign_str(name,gqh->name[gqh->nparam_c]);
    assign_str(value,gqh->value[gqh->nparam_c]);
    assign_str(unit,gqh->unit[gqh->nparam_c]);
    gqh->nparam_c++;
    
    return 1;
  }
  else
    gqh->nparam_c = 0;
    
  return 0;
}

/* Return a value corresponding to name */
extern int qhop_get_value(qhop_t gqh,char *name,double *x)
{
  int i;
  
  for(i=0; (i<gqh->nparam); i++) 
    if (gmx_strcasecmp(gqh->name[i],name) == 0) {
      *x = strtod(gqh->value[i],NULL);
      return 1;
    }
    
  return 0;
}

/* Liberate memory */
extern void qhop_done(qhop_t gqh)
{
  int i;
  
  for(i=0; (i<gqh->nparam); i++) {
    sfree(gqh->name[i]);
    sfree(gqh->value[i]);
    sfree(gqh->unit[i]);
  }
  if (gqh->nparam > 0) {
    sfree(gqh->name);
    sfree(gqh->value);
    sfree(gqh->unit);
  }
  if (gqh->donor)
    sfree(gqh->donor);
  if (gqh->acceptor)
    sfree(gqh->acceptor);

}
