/*
 * 
 *                This source code is NOT part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "typedefs.h"
#include "macros.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "copyrite.h"
#include "gromacs/utility/futil.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/smalloc.h"
#include "nrnb.h"
#include "gstat.h"
#include "gromacs/utility/fatalerror.h"


#define G_REF1      0
#define G_REF2      1
#define G_REF3      2
#define G_SDF       3
#define NDX_REF1    4
#define NDX_REF2    5
#define NDX_REF3    6
#define G_REFMOL    7

static void i_write(FILE *output, int value)
{
  if(fwrite(&value,sizeof(int),1,output) != 1)
  {
    gmx_fatal(FARGS,"Error writing to output file");
  }
}


static void f_write(FILE *output,float value)
{
  if(fwrite(&value,sizeof(float),1,output) != 1)
  {
    gmx_fatal(FARGS,"Error writing to output file");
  }
}


static void do_sdf(const char *fnNDX,const char *fnTPS,const char *fnTRX, 
                   const char *fnSDF, const char *fnREF, gmx_bool bRef, 
                   rvec cutoff, real binwidth, int mode, rvec triangle, 
                   rvec dtri, const output_env_t oenv)
{
  FILE       *fp;
  t_trxstatus *status;
  int        ng,natoms,i,j,k,l,X,Y,Z,lc,dest;
  ivec       nbin;
  int        ***count;
  /* real       ***sdf; */
  real       sdf,min_sdf=1e10,max_sdf=0;
  char       **grpname;
  int        *isize;
  int        isize_cg=0;
  int        isize_ref=3;
  int        ref_resind[3]={0};
  int        nrefmol=0,refc=0;
  atom_id    **index;
  atom_id    *index_cg=NULL;
  atom_id    *index_ref=NULL;
  real       t,boxmin,hbox,normfac;
  real       invbinw;
  rvec       tri_upper,tri_lower;
  rvec       *x,xcog,dx,*x_i1,xi,*x_refmol;
  matrix     box;
  matrix     rot; /* rotation matrix := unit vectors for the molecule frame */
  rvec       k_mol,i1_mol,i2_mol,dx_mol;
  real       delta;
  atom_id    ix,jx;
  t_topology top;
  gmx_rmpbc_t  gpbc=NULL;
  int        ePBC=-1;
  t_pbc      pbc;
  gmx_bool       bTop=FALSE,bRefDone=FALSE,bInGroup=FALSE;
  char       title[STRLEN];


  /* Read Topology */
  if (fnTPS) {
    bTop=read_tps_conf(fnTPS,title,&top,&ePBC,&x,NULL,box,TRUE);
  }
  


  if ( !bTop ) {
    fprintf(stderr,"\nNeed tpr-file to make a reference structure.\n");
    fprintf(stderr,"Option -r will be ignored!\n");
    bRef = FALSE;
  }


  /* Allocate memory for 4 groups, 3 dummy groups and a group for the ref 
structure if needed */
  ng = 4;
  snew(grpname,ng);
  /* the dummy groups are used to dynamically store triples of atoms */
  /* for molecular coordinate systems */
  if ( bRef )
    {
      snew(isize,ng+4);
      snew(index,ng+4);
    }
  else 
    {
      snew(isize,ng+3);
      snew(index,ng+3);
    }


  /* Read the index groups */
  fprintf(stderr,"\nSelect the 3 reference groups and the SDF group:\n");
  if (fnTPS)
    get_index(&top.atoms,fnNDX,ng,isize,index,grpname);
  else
    rd_index(fnNDX,ng,isize,index,grpname);


  isize[NDX_REF1]=isize[G_REF1];
  for (i=NDX_REF1; i<=NDX_REF3; i++)
    snew(index[i],isize[NDX_REF1]);


  /* Read first frame and check it */
  natoms=read_first_x(oenv,&status,fnTRX,&t,&x,box);
  if ( !natoms )
    gmx_fatal(FARGS,"Could not read coordinates from statusfile!\n");


  /* check with topology */
  if (fnTPS)
    if ( natoms > top.atoms.nr )
      gmx_fatal(FARGS,
                "Trajectory (%d atoms) does not match topology (%d atoms)!\n",
                natoms,top.atoms.nr);


  /* check with index groups */
  for (i=0; i<ng; i++)
    for (j=0; j<isize[i]; j++)
      if ( index[i][j] >= natoms )
        gmx_fatal(FARGS,"Atom index (%d) in index group %s (%d atoms) larger "
                    "than number of atoms in trajectory (%d atoms)!\n",
                    index[i][j],grpname[i],isize[i],natoms);


  /* check reference groups */
  if ( mode == 1 )
    {
      if ( isize[G_REF1] != isize[G_REF2] || isize[G_REF1] != isize[G_REF3] || 
           isize[G_REF2] != isize[G_REF3] )
        gmx_fatal(FARGS,"For single particle SDF, all reference groups"
                    "must have the same size.\n");


      /* for single particle SDF dynamic triples are not needed */
      /* so we build them right here */


      /* copy all triples from G_REFx to NDX_REFx */    
      for (i=0; i<isize[G_REF1]; i++)
        {
          /* check if all three atoms come from the same molecule */
          for (j=G_REF1; j<=G_REF3; j++)
            ref_resind[j] = top.atoms.atom[index[j][i]].resind;


          if ( ref_resind[G_REF1] != ref_resind[G_REF2] ||
                 ref_resind[G_REF2] != ref_resind[G_REF3] ||
                 ref_resind[G_REF3] != ref_resind[G_REF1] )
              {
                fprintf(stderr,"\nWarning: reference triple (%d) will be skipped.\n",i);
                fprintf(stderr,  "         resnr[1]: %d, resnr[2]: %d, resnr[3]: %d\n",
                        ref_resind[G_REF1],ref_resind[G_REF2], ref_resind[G_REF3]);
                isize[NDX_REF1]--;
                for (j=NDX_REF1; j<=NDX_REF3; j++)
                  srenew(index[j],isize[NDX_REF1]);
                continue;
              }
          else
            /* check if all entries are unique*/
            if ( index[G_REF1][i] == index[G_REF2][i] ||
                 index[G_REF2][i] == index[G_REF3][i] ||
                 index[G_REF3][i] == index[G_REF1][i] )
              {
                fprintf(stderr,"Warning: reference triple (%d) will be skipped.\n",i);
                fprintf(stderr,  "         index[1]: %d, index[2]: %d, index[3]: %d\n",
                        index[G_REF1][i],index[G_REF2][i],index[G_REF3][i]);    
                isize[NDX_REF1]--;
                for (j=NDX_REF1; j<=NDX_REF3; j++)
                  srenew(index[j],isize[NDX_REF1]);
                continue;
              }
            else /* everythings fine, copy that one */
              for (j=G_REF1; j<=G_REF3; j++)
                index[j+4][i] = index[j][i];
        }
    }
  else if ( mode == 2 )
    {
      if ( isize[G_REF1] != isize[G_REF2] )
        gmx_fatal(FARGS,"For two particle SDF, reference groups 1 and 2"
                    "must have the same size.\n");


      for (i=0; i<isize[G_REF1]; i++)
        {
          /* check consistency for atoms 1 and 2 */
          for (j=G_REF1; j<=G_REF2; j++)
            ref_resind[j] = top.atoms.atom[index[j][i]].resind;


          if ( ref_resind[G_REF1] != ref_resind[G_REF2] ||
               index[G_REF1][i] == index[G_REF2][i] )
            {
              if ( ref_resind[G_REF1] != ref_resind[G_REF2] )
                {
                  fprintf(stderr,"\nWarning: bond (%d) not from one molecule."
                          "Will not be used for SDF.\n",i);
                  fprintf(stderr,  "         resnr[1]: %d, resnr[2]: %d\n",
                          ref_resind[G_REF1],ref_resind[G_REF2]);
                }
              else
                {
                  fprintf(stderr,"\nWarning: atom1 and atom2 are identical."
                          "Bond (%d) will not be used for SDF.\n",i);
                  fprintf(stderr,  "         index[1]: %d, index[2]: %d\n",
                          index[G_REF1][i],index[G_REF2][i]);
                }
              for (j=NDX_REF1; j<=NDX_REF2; j++)
                {
                  for (k=i; k<isize[G_REF1]-1; k++)
                    index[j][k]=index[j][k+1];
                  isize[j]--;
                  srenew(index[j],isize[j]);
                }
            }
        }
    }


  /* Read Atoms for refmol group */
  if ( bRef )
    {
      snew(index[G_REFMOL],1);


      for (i=G_REF1; i<=G_REF3; i++)
        ref_resind[i] = top.atoms.atom[index[i][0]].resind;


      for (i=0; i<natoms; i++)
        {
          if (  ref_resind[G_REF1] == top.atoms.atom[i].resind ||
                ref_resind[G_REF2] == top.atoms.atom[i].resind ||
                ref_resind[G_REF3] == top.atoms.atom[i].resind )
            nrefmol++;
        }
      srenew(index[G_REFMOL],nrefmol);
      isize[G_REFMOL] = nrefmol;
      nrefmol = 0;


      for (i=0; i<natoms; i++)
        {
          if (  ref_resind[G_REF1] == top.atoms.atom[i].resind ||
                ref_resind[G_REF2] == top.atoms.atom[i].resind ||
                ref_resind[G_REF3] == top.atoms.atom[i].resind )
            {
              index[G_REFMOL][nrefmol] = i;
              nrefmol++;
            }
        }
    }


  /* initialize some stuff */
  boxmin = min( norm(box[XX]), min( norm(box[YY]), norm(box[ZZ]) ) );
  hbox   = boxmin / 2.0;


  for (i=0; i<DIM; i++)
    {
      cutoff[i] = cutoff[i] / 2;
      nbin[i]   = (int)(2 * cutoff[i] / binwidth) + 1;
      invbinw = 1.0 / binwidth;
      tri_upper[i] = triangle[i] + dtri[i];
      tri_upper[i] = sqr(tri_upper[i]);
      tri_lower[i] = triangle[i] - dtri[i];
      tri_lower[i] = sqr(tri_lower[i]);
    }


  /* Allocate the array's for sdf */
  snew(count,nbin[0]+1);
  for(i=0; i<nbin[0]+1; i++) 
    {
      snew(count[i],nbin[1]+1);
      for (j=0; j<nbin[1]+1; j++) 
        snew(count[i][j],nbin[2]+1);
    }


  /* Allocate space for the coordinates */
  snew(x_i1,isize[G_SDF]);
  snew(x_refmol,isize[G_REFMOL]);
  for (i=0; i<isize[G_REFMOL]; i++)
    for (j=XX; j<=ZZ; j++)
      x_refmol[i][j] = 0;


  normfac = 0;

  gpbc = gmx_rmpbc_init(&top.idef,ePBC,natoms,box);
  
  do {
    /* Must init pbc every step because of pressure coupling */
    set_pbc(&pbc,ePBC,box);
    gmx_rmpbc(gpbc,natoms,box,x);
  
    /* Dynamically build the ref triples */
    if ( mode == 2 )
      {
        isize[NDX_REF1]=0;
        for (j=NDX_REF1; j<=NDX_REF3; j++)
          srenew(index[j],isize[NDX_REF1]+1);


        /* consistancy of G_REF[1,2] has already been check */
        /* hence we can look for the third atom right away */


        for (i=0; i<isize[G_REF1]; i++)
          {
            for (j=0; j<isize[G_REF3]; j++)
              {
                /* Avoid expensive stuff if possible */
                if ( top.atoms.atom[index[G_REF1][i]].resind != 
                     top.atoms.atom[index[G_REF3][j]].resind &&
                     index[G_REF1][i] != index[G_REF3][j] &&
                     index[G_REF2][i] != index[G_REF3][j] )
                  {
                    pbc_dx(&pbc,x[index[G_REF1][i]],x[index[G_REF3][j]],dx);
                    delta = norm2(dx);
                    if ( delta < tri_upper[G_REF1] &&
                         delta > tri_lower[G_REF1] )
                      {
                        pbc_dx(&pbc,x[index[G_REF2][i]],x[index[G_REF3][j]],dx);
                        delta = norm2(dx);
                        if ( delta < tri_upper[G_REF2] &&
                             delta > tri_lower[G_REF2] )
                          {
                            /* found triple */
                            index[NDX_REF1][isize[NDX_REF1]]=index[G_REF1][i];
                            index[NDX_REF2][isize[NDX_REF1]]=index[G_REF2][i];
                            index[NDX_REF3][isize[NDX_REF1]]=index[G_REF3][j];


                            /* resize groups */
                            isize[NDX_REF1]++;
                            for (k=NDX_REF1; k<=NDX_REF3; k++)
                              srenew(index[k],isize[NDX_REF1]+1);
                          }
                      }
                  }
              }
          }
      }
    else if ( mode ==3 )
      {
        isize[NDX_REF1]=0;
        for (j=NDX_REF1; j<=NDX_REF3; j++)
          srenew(index[j],isize[NDX_REF1]+1);

        /* consistancy will be checked while searching */


        for (i=0; i<isize[G_REF1]; i++)
          {
            for (j=0; j<isize[G_REF2]; j++)
              {
                /* Avoid expensive stuff if possible */
                if ( top.atoms.atom[index[G_REF1][i]].resind != 
                     top.atoms.atom[index[G_REF2][j]].resind &&
                     index[G_REF1][i] != index[G_REF2][j] )
                  {
                    pbc_dx(&pbc,x[index[G_REF1][i]],x[index[G_REF2][j]],dx);
                    delta = norm2(dx);
                    if ( delta < tri_upper[G_REF3] &&
                         delta > tri_lower[G_REF3] )
                      {
                        for (k=0; k<isize[G_REF3]; k++)
                          {
                            if ( top.atoms.atom[index[G_REF1][i]].resind != 
                                 top.atoms.atom[index[G_REF3][k]].resind &&
                                 top.atoms.atom[index[G_REF2][j]].resind != 
                                 top.atoms.atom[index[G_REF3][k]].resind &&
                                 index[G_REF1][i] != index[G_REF3][k] &&
                                 index[G_REF2][j] != index[G_REF3][k])
                              {
                                pbc_dx(&pbc,x[index[G_REF1][i]],x[index[G_REF3][k]],dx);
                                delta = norm2(dx);
                                if ( delta < tri_upper[G_REF1] &&
                                     delta > tri_lower[G_REF1] )
                                  {
                                    pbc_dx(&pbc,x[index[G_REF2][j]],x[index[G_REF3][k]],dx);
                                    delta = norm2(dx);
                                    if ( delta < tri_upper[G_REF2] &&
                                         delta > tri_lower[G_REF2] )
                                      {
                                        /* found triple */
                                        index[NDX_REF1][isize[NDX_REF1]]=index[G_REF1][i];
                                        index[NDX_REF2][isize[NDX_REF1]]=index[G_REF2][j];
                                        index[NDX_REF3][isize[NDX_REF1]]=index[G_REF3][k];
                                    
                                        /* resize groups */
                                        isize[NDX_REF1]++;
                                        for (l=NDX_REF1; l<=NDX_REF3; l++)
                                          srenew(index[l],isize[NDX_REF1]+1);
                                      }
                                  }
                              }
                          }
                      }
                  }
              }
          }
      }
 
    for (i=0; i<isize[NDX_REF1]; i++)
      {
        /* setup the molecular coordinate system (i',j',k') */
        /* because the coodinate system of the box forms a unit matrix */
        /* (i',j',k') is identical with the rotation matrix */
        clear_mat(rot);


        /* k' = unitv(r(atom0) - r(atom1)) */
        pbc_dx(&pbc,x[index[NDX_REF1][i]],x[index[NDX_REF2][i]],k_mol);
        unitv(k_mol,rot[2]);
        
        /* i' = unitv(k' x (r(atom2) - r(atom1))) */
        pbc_dx(&pbc,x[index[NDX_REF3][i]],x[index[NDX_REF2][i]],i1_mol);
        cprod(i1_mol,rot[2],i2_mol);
        unitv(i2_mol,rot[0]);
      
        /* j' = k' x i' */
        cprod(rot[2],rot[0],rot[1]);


        /* set the point of reference */
        if ( mode == 2 )
          copy_rvec(x[index[NDX_REF3][i]],xi);
        else
          copy_rvec(x[index[NDX_REF1][i]],xi);


        /* make the reference */
        if ( bRef )
          {
            for (j=0; j<isize[G_REFMOL]; j++)
              {
                pbc_dx(&pbc,xi,x[index[G_REFMOL][j]],dx);
                mvmul(rot,dx,dx_mol);
                rvec_inc(x_refmol[j],dx_mol);
                for(k=XX; k<=ZZ; k++)
                   x_refmol[j][k] = x_refmol[j][k] / 2;
              }
          }


        /* Copy the indexed coordinates to a continuous array */
        for(j=0; j<isize[G_SDF]; j++)
          copy_rvec(x[index[G_SDF][j]],x_i1[j]);
        
        /* count the SDF */
        for(j=0; j<isize[G_SDF]; j++) 
          {
            /* Exclude peaks from the reference set */
            bInGroup=FALSE;
            for (k=NDX_REF1; k<=NDX_REF3; k++)
              if ( index[G_SDF][j] == index[k][i] )
                bInGroup=TRUE;


            if ( !bInGroup )
              {
                pbc_dx(&pbc,xi,x_i1[j],dx);
            
                /* transfer dx to the molecular coordinate system */
                mvmul(rot,dx,dx_mol);


                /* check cutoff's and count */
                if ( dx_mol[XX] > -cutoff[XX] && dx_mol[XX] < cutoff[XX] )
                  if ( dx_mol[YY] > -cutoff[YY] && dx_mol[YY] < cutoff[YY] )
                    if ( dx_mol[ZZ] > -cutoff[ZZ] && dx_mol[ZZ] < cutoff[ZZ] )
                      {
                        X = (int)(floor(dx_mol[XX]*invbinw)) + (nbin[XX]-1)/2 
+1;
                        Y = (int)(floor(dx_mol[YY]*invbinw)) + (nbin[YY]-1)/2 
+1;
                        Z = (int)(floor(dx_mol[ZZ]*invbinw)) + (nbin[ZZ]-1)/2 
+1;
                        count[X][Y][Z]++;
                        normfac++;
                      }
              }
          }
      }
  } while (read_next_x(oenv,status,&t,natoms,x,box));
  fprintf(stderr,"\n");
  
  gmx_rmpbc_done(gpbc);


  close_trj(status);
  
  sfree(x);


  /* write the reference strcture*/
  if ( bRef )
    {
      fp=gmx_ffopen(fnREF,"w"); 
      fprintf(fp,"%s\n",title);
      fprintf(fp,"  %d\n",isize[G_REFMOL]);


      for (i=0; i<isize[G_REFMOL]; i++)
        fprintf(fp,"%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",
                top.atoms.resinfo[top.atoms.atom[index[G_REFMOL][i]].resind].nr,
                *(top.atoms.resinfo[top.atoms.atom[index[G_REFMOL][i]].resind].name),
                *(top.atoms.atomname[index[G_REFMOL][i]]),i+1,
                -1*x_refmol[i][XX],-1*x_refmol[i][YY],-1*x_refmol[i][ZZ]);
      /* Inserted -1* on the line above three times */
      fprintf(fp,"   10.00000   10.00000   10.00000\n");
      gmx_ffclose(fp);
      fprintf(stderr,"\nWrote reference structure. (%d Atoms)\n",isize[G_REFMOL]);
    }


  /* Calculate the mean probability density */
  fprintf(stderr,"\nNumber of configurations used for SDF: %d\n",(int)normfac);


  normfac = nbin[0]*nbin[1]*nbin[2] / normfac;


  fprintf(stderr,"\nMean probability density: %f\n",1/normfac);


  /* normalize the SDF and write output */
  /* see http://www.csc.fi/gopenmol/index.phtml for documentation */
  fp=gmx_ffopen(fnSDF,"wb"); 


  /* rank */
  i_write(fp,3);


  /* Type of surface */
  i_write(fp,42);


  /* Zdim, Ydim, Xdim */
  for (i=ZZ; i>=XX; i--)
    i_write(fp,nbin[i]);


  /* [Z,Y,X][min,max] (box corners in Angstroem)*/
  for (i=ZZ; i>=XX; i--)
    {
      f_write(fp,-cutoff[i]*10);
      f_write(fp,cutoff[i]*10);
    }


/* Original Code
  for (i=1; i<nbin[2]+1; i++)
    for (j=1; j<nbin[1]+1; j++)
      for (k=1; k<nbin[0]+1; k++) 
        {
          sdf = normfac * count[k][j][i];
          if ( sdf < min_sdf ) min_sdf = sdf;
          if ( sdf > max_sdf ) max_sdf = sdf;
          f_write(fp,sdf);
        }*/
/* Changed Code to Mirror SDF to correct coordinates */
  for (i=nbin[2]; i>0; i--)
    for (j=nbin[1]; j>0; j--)
      for (k=nbin[0]; k>0; k--)
        {
          sdf = normfac * count[k][j][i];
          if ( sdf < min_sdf ) min_sdf = sdf;
          if ( sdf > max_sdf ) max_sdf = sdf;
          f_write(fp,sdf);
        }

  fprintf(stderr,"\nMin: %f Max: %f\n",min_sdf,max_sdf);


  gmx_ffclose(fp); 


  /* Give back the mem */
  for(i=0; i<nbin[0]+1; i++)
    {
      for (j=0; j<nbin[1]+1; j++)
        {
          sfree(count[i][j]);
        }
      sfree(count[i]);
    }
  sfree(count);
}

int gmx_sdf(int argc,char *argv[])
{
  const char *desc[] = {
    "[TT]g_sdf[tt] calculates the spatial distribution function (SDF) of a set of atoms",
    "within a coordinate system defined by three atoms. There is single body, ",
    "two body and three body SDF implemented (select with option [TT]-mode[tt]). ",
    "In the single body case the local coordinate system is defined by using",
    "a triple of atoms from one single molecule, for the two and three body case",
    "the configurations are dynamically searched complexes of two or three",
    "molecules (or residues) meeting certain distance consitions (see below).[PAR]",
    "The program needs a trajectory, a GROMACS run input file and an index ",
    "file to work. ",
    "You have to setup 4 groups in the index file before using g_sdf: [PAR]",
    "The first three groups are used to define the SDF coordinate system.",
    "The program will dynamically generate the atom triples according to ",
    "the selected [TT]-mode[tt]: ", 
    "In [TT]-mode[tt] 1 the triples will be just the 1st, 2nd, 3rd, ... atoms from ",
    "groups 1, 2 and 3. Hence the nth entries in groups 1, 2 and 3 must be from the",
    "same residue. In [TT]-mode[tt] 2 the triples will be 1st, 2nd, 3rd, ... atoms from",
    "groups 1 and 2 (with the nth entries in groups 1 and 2 having the same res-id).",
    "For each pair from groups 1 and 2  group 3 is searched for an atom meeting the",
    "distance conditions set with [TT]-triangle[tt] and [TT]-dtri[tt] relative to atoms 1 and 2. In",
    "[TT]-mode[tt] 3 for each atom in group 1 group 2 is searched for an atom meeting the",
    "distance condition and if a pair is found group 3 is searched for an atom",
    "meeting the further conditions. The triple will only be used if all three atoms",
    "have different res-id's.[PAR]",
    "The local coordinate system is always defined using the following scheme:",
    "Atom 1 will be used as the point of origin for the SDF. ",
    "Atom 1 and 2 will define the principle axis (Z) of the coordinate system.",
    "The other two axis will be defined inplane (Y) and normal (X) to the plane through",
    "Atoms 1, 2 and 3. ",
    "The fourth group",
    "contains the atoms for which the SDF will be evaluated.[PAR]",
    "For [TT]-mode[tt] 2 and 3 you have to define the distance conditions for the ",
    "2 resp. 3 molecule complexes to be searched for using [TT]-triangle[tt] and [TT]-dtri[tt].[PAR]",
    "The SDF will be sampled in cartesian coordinates.",
    "Use [TT]-grid x y z[tt] to define the size of the SDF grid around the ",
    "reference molecule. ",
    "The Volume of the SDF grid will be V=x*y*z (nm^3). ",
    "Use [TT]-bin[tt] to set the binwidth for grid.[PAR]",
    "The output will be a binary 3D-grid file ([TT]gom_plt.dat[tt]) in the [REF].plt[ref] format that can be be",
    "read directly by gOpenMol. ",
    "The option [TT]-r[tt] will generate a [REF].gro[ref] file with the reference molecule(s) transferred to",
    "the SDF coordinate system. Load this file into gOpenMol and display the",
    "SDF as a contour plot (see http://www.csc.fi/gopenmol/index.phtml for ",
    "further documentation). [PAR]",
    "For further information about SDF's have a look at: A. Vishnyakov, JPC A, 105,",
    "2001, 1702 and the references cited within."
  };
  output_env_t oenv;
  static gmx_bool bRef=FALSE;
  static int mode=1;
  static rvec triangle={0.0,0.0,0.0};
  static rvec dtri={0.03,0.03,0.03};
  static rvec cutoff={1,1,1};
  static real binwidth=0.05;
  t_pargs pa[] = {
    { "-mode",     FALSE, etINT, {&mode},
      "SDF in [1,2,3] particle mode"},
    { "-triangle", FALSE, etRVEC, {&triangle},
      "r(1,3), r(2,3), r(1,2)"},
    { "-dtri",     FALSE, etRVEC, {&dtri},
      "dr(1,3), dr(2,3), dr(1,2)"},
    { "-bin",      FALSE, etREAL, {&binwidth},
      "Binwidth for the 3D-grid (nm)" },
    { "-grid",      FALSE, etRVEC, {&cutoff},
      "Size of the 3D-grid (nm,nm,nm)"}
  };
#define NPA asize(pa)
  const char       *fnTPS,*fnNDX,*fnREF;
  
  t_filenm   fnm[] = {
    { efTRX, "-f",  NULL,     ffREAD },
    { efNDX, NULL,  NULL,     ffREAD },
    { efTPS, NULL,  NULL,     ffOPTRD },
    { efDAT, "-o",  "gom_plt",     ffWRITE },
    { efSTO, "-r",  "refmol",  ffOPTWR },
  };
#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME,
                    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL,&oenv);


  fnTPS = ftp2fn_null(efTPS,NFILE,fnm);
  fnNDX = ftp2fn_null(efNDX,NFILE,fnm);
  fnREF = opt2fn_null("-r",NFILE,fnm);
  bRef  = opt2bSet("-r",NFILE,fnm);


  
  if (!fnNDX)
    gmx_fatal(FARGS,"No index file specified\n"
                "             Nothing to do!");


  if (!fnTPS)
    gmx_fatal(FARGS,"No topology file specified\n"
                "             Nothing to do!");


  if ( bRef && (fn2ftp(fnREF) != efGRO))
    {
      fprintf(stderr,"\nOnly GROMACS format is supported for reference structures.\n");
      fprintf(stderr,"Option -r will be ignored!\n");
      bRef = FALSE;
    }


  if ((mode < 1) || (mode > 3))
    gmx_fatal(FARGS,"Wrong -mode selection. Chose 1-, 2- oder 3-particle mode.\n");


  do_sdf(fnNDX,fnTPS,ftp2fn(efTRX,NFILE,fnm),opt2fn("-o",NFILE,fnm),
         fnREF,bRef,cutoff,binwidth,mode,triangle,dtri,oenv);


  gmx_thanx(stderr);
  
  return 0;
}

int
main(int argc, char *argv[])
{
  gmx_sdf(argc,argv);
  return 0;
}
