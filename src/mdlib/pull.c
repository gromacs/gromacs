#include "fatal.h"
#include "futil.h"
#include "rdgroup.h"
#include "statutil.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "vec.h" 
#include "typedefs.h"
#include "network.h"
#include "filenm.h"
#include "string.h"
#include "smalloc.h"
#include "pull.h"

/* check if we're close enough to the target position */
static bool check_convergence(t_pull *pull) {
  bool bTest = TRUE;
  real tol;
  int i,m;
  rvec tmp1;
  rvec dr; /* distance between target and current position */

  tol = pull->tolerance;

  for (i=0;i<pull->pull.n;i++) {
    rvec_sub(pull->ref.x_unc[0],pull->pull.x_unc[i],tmp1);
    rvec_add(pull->pull.xtarget[i],tmp1,dr);
    /* multiply by pull->dims to pick the elements we want */
    for(m=0;m<DIM;m++) 
      dr[m] = pull->dims[m]*dr[m];
    bTest = bTest && (norm(dr) < tol);
  }
  
  return bTest;
}

static void do_start(t_pull *pull, rvec *x, matrix box, t_mdatoms *md, 
	      real dt, int step, t_topology *top) 
{
  int i,j,ii,m;
  rvec dr,dx,tmp;     
  bool bThereYet,bDump;
  static int nout = 0;
  rvec ds;

  /* pull.xtarget[i] is the desired position of group i with respect
     to the reference group */

  /* bThereYet is true if all groups have reached their desired position */
  bThereYet = check_convergence(pull);
  
  if (pull->bVerbose) {
    for (i=0;i<pull->pull.n;i++) {
      copy_rvec(pull->pull.xtarget[i],tmp);
      fprintf(stderr,"Group %d goal:%8.3f%8.3f%8.3f from reference\n",
	      i,tmp[0],tmp[1],tmp[2]);
    }
  }

  /* some groups are not there yet! */
  if (!bThereYet) { 
    for (i=0;i<pull->pull.n;i++) {
      /* move this group closer to its target position */
      rvec_sub(pull->ref.x_unc[0],pull->pull.x_unc[i],tmp);
      rvec_add(pull->pull.xtarget[i],tmp,dr);

      /* multiply by pull->dims to pick the elements we want */
      for(m=0;m<DIM;m++) 
	dr[m] = pull->dims[m]*dr[m];

      /* move over pull->xlt_rate in the direction of dr */
      svmul(pull->xlt_rate/norm(dr),dr,dx);

      if (pull->bVerbose) 
	fprintf(stderr,"To go:%10.2e %10.2e %10.2e\n",
		dr[XX],dr[YY],dr[ZZ]);

      /* update the actual positions */
      for (j=0;j<pull->pull.ngx[i];j++) {
	ii = pull->pull.idx[i][j];
	rvec_add(x[ii],dx,x[ii]);
      }

      /* get the new center of mass in pull.x_unc[i] */
      (void)calc_com(x,pull->pull.ngx[i],pull->pull.idx[i],
		     md,pull->pull.x_unc[i],box);
    }
  }

  bDump = check_convergence(pull); 
  
  if (bDump) {
    /* we have to write a coordinate file and reset the desired position */
    for (i=0;i<pull->pull.n;i++) {
      for (m=0;m<DIM;m++) {
	pull->pull.xtarget[i][m] += pull->pull.dir[i][m]*pull->xlt_incr/
	  norm(pull->pull.dir[i]);
      }
      
      if (pull->bVerbose) 
	fprintf(stderr,"New target position: %8.3f%8.3f%8.3f\n",
		pull->pull.xtarget[i][0],pull->pull.xtarget[i][1],
		pull->pull.xtarget[i][2]);
    }
    dump_conf(pull,x,box,top,nout,step*dt/1000); 
    nout++;
  }
}

static void do_umbrella(t_pull *pull, rvec *x,rvec *f,matrix box, 
			t_mdatoms *md) 
{
  int i,j,m,g;
  real mi;
  rvec cm;    /* center of mass displacement of reference */
  rvec dr;    /* distance from com to potential center */

  /* loop over the groups that are being umbrella sampled */
  for (i=0;i<pull->pull.n;i++) {
    
    /* pull->dims selects which components (x,y,z) we want */
    for (m=0;m<DIM;m++) {
      dr[m]=pull->dims[m]*(pull->pull.x_ref[i][m]-pull->pull.x_unc[i][m]); 
      if (dr[m] >  box[m][m]/2) dr[m]-=box[m][m];
      if (dr[m] < -box[m][m]/2) dr[m]+=box[m][m];
    }
    
    /* f = um_width*x */
    svmul(pull->um_width,dr,pull->pull.f[i]);

    /* distribute the force over all atoms in the group */
    for (j=0;j<pull->pull.ngx[i];j++) 
      for (m=0;m<DIM;m++) {
	mi = md->massT[pull->pull.idx[i][j]];
	f[pull->pull.idx[i][j]][m] += 
	  mi*(pull->pull.f[i][m])/(pull->pull.tmass[i]);
      }
  }
  
  /* reset center of mass of the reference group to its starting value */
  rvec_sub(pull->ref.x_unc[0],pull->ref.x_con[0],cm);
  for (i=0;i<pull->ref.ngx[0];i++)
    rvec_sub(x[pull->ref.idx[0][i]],cm,x[pull->ref.idx[0][i]]);
}

/* this implements a constraint run like SHAKE does. */
static void do_constraint(t_pull *pull, rvec *x, matrix box, t_mdatoms *md, 
			  real dt) 
{
  
  rvec r_ij, /* x_con[i] com of i in prev. step. Obeys constr. -> r_ij   */
    unc_ij,  /* x_unc[i] com of i this step, before constr.    -> unc_ij */
    ref_ij;  /* x_ref[i] com of i at t0, not updated           -> ref_ij */

  rvec *rinew;           /* current 'new' position of group i */
  rvec *rjnew;           /* current 'new' position of group j */
  real *direction;       /* direction of dr relative to r_ij */
  double lambda, rm, mass, tol = 0.00001;
  bool bConverged = FALSE;
  int n=0,i,ii,j,m,max_iter=1000;
  int ref;
  double x1,x2,q,a,b,c;  /* for solving the quadratic equation, 
                            see Num. Recipes in C ed 2 p. 184 */
  rvec *dr;              /* correction for group i */
  rvec *ref_dr;          /* correction for group j */
  rvec tmp,tmp2,tmp3,sum;

  if (pull->bCyl) {
    snew(ref_dr,pull->pull.n);
    snew(rjnew,pull->pull.n);
  } else {
    snew(ref_dr,1);
    snew(rjnew,1);
  }
  snew(dr,pull->pull.n);
  snew(rinew,pull->pull.n);
  snew(direction,pull->pull.n);

  /* copy the current unconstraint positions for use in iterations. We 
     iterate until rinew[i] and rjnew[j] obey the constraints. Then
     rinew - pull.x_unc[i] is the correction dr to group i */
  for (i=0;i<pull->pull.n;i++) 
    copy_rvec(pull->pull.x_unc[i],rinew[i]);
  if (pull->bCyl)
    for (i=0;i<pull->pull.n;i++)
      copy_rvec(pull->dyna.x_unc[i],rjnew[i]);
  else
    copy_rvec(pull->ref.x_unc[0],rjnew[0]);

  while (!bConverged && n<max_iter) {
    /* loop over all constraints */
    for (i=0;(i<pull->pull.n);i++) {
      
      if (pull->bVerbose)
	fprintf(stderr,"group %d, iteration %d\n",i,n);

      if (pull->bCyl) {
	rvec_sub(pull->dyna.x_con[i],pull->pull.x_con[i],r_ij);
	rvec_sub(rjnew[i],rinew[i],unc_ij);
	rvec_sub(pull->dyna.x_ref[i],pull->pull.x_ref[i],ref_ij);
      } else {
	rvec_sub(pull->ref.x_con[0],pull->pull.x_con[i],r_ij);
	rvec_sub(rjnew[0],rinew[i],unc_ij);
	rvec_sub(pull->ref.x_ref[0],pull->pull.x_ref[i],ref_ij);
      }

      /* select components we want */
      for (m=0;m<DIM;m++) {
	r_ij[m]     *= pull->dims[m];
	unc_ij[m]   *= pull->dims[m];
	ref_ij[m]   *= pull->dims[m];

	/* correct for PBC. Is this correct? */
	if (r_ij[m]   < -0.5*box[m][m]) r_ij[m]   += box[m][m];
	if (r_ij[m]   >  0.5*box[m][m]) r_ij[m]   -= box[m][m];
	if (unc_ij[m] < -0.5*box[m][m]) unc_ij[m] += box[m][m];
	if (unc_ij[m] >  0.5*box[m][m]) unc_ij[m] -= box[m][m];
	if (ref_ij[m] < -0.5*box[m][m]) ref_ij[m] += box[m][m];
	if (ref_ij[m] >  0.5*box[m][m]) ref_ij[m] -= box[m][m];
      }

      
      if (pull->bCyl) 
	rm = 1/pull->pull.tmass[i] + 1/pull->dyna.tmass[i];
      else
	rm = 1/pull->pull.tmass[i] + 1/pull->ref.tmass[0];

      a = iprod(r_ij,r_ij)*dt*dt*dt*dt*rm*rm; 
      b = iprod(unc_ij,r_ij)*2*dt*dt*rm;
      c = iprod(unc_ij,unc_ij) - norm2(ref_ij);

      if (b<0) 
	q = -0.5*(b-sqrt(b*b-4*a*c));
      else
	q = -0.5*(b+sqrt(b*b-4*a*c));
      x1 = q/a; x2 = c/q;
      lambda = x1 > 0 ? x1 : x2;

      if (pull->bVerbose) 
	fprintf(stderr,"\nax^2+bx+c=0: a=%e b=%e c=%e\n"
		"x1=%e x2=%e sum:%e,%e, lambda:%e\n",a,b,c,x1,x2,
		a*x1*x1+b*x1+c,a*x2*x2+b*x2+c,lambda);

      /* the position corrections dr due to the constraint are: */
      if (pull->bCyl) {
	svmul(-dt*dt*lambda/pull->pull.tmass[i],r_ij,dr[i]);
	svmul(dt*dt*lambda/pull->dyna.tmass[i],r_ij,ref_dr[i]);
      } else {
	svmul(-dt*dt*lambda/pull->pull.tmass[i],r_ij,dr[i]);
	svmul(dt*dt*lambda/pull->ref.tmass[0],r_ij,ref_dr[0]);
      }

      /* and the direction of the constraint force: */
      direction[i] = cos_angle(r_ij,dr[i]);

      /* DEBUG */
      if (pull->bVerbose) { 
	fprintf(stderr,"Direction: %f\n",direction[i]);
	if (pull->bCyl) {
	  rvec_sub(rinew[i],rjnew[i],tmp);
	  rvec_sub(pull->pull.x_ref[i],pull->dyna.x_ref[i],tmp2);
	} else {
	  rvec_sub(pull->pull.x_ref[i],pull->ref.x_ref[0],tmp2);
	  rvec_sub(rinew[i],rjnew[0],tmp);
	}
	rvec_sub(dr[i],ref_dr[0],tmp3);
	for (m=0;m<DIM;m++) {
	  tmp[m]  *= pull->dims[m];
	  tmp2[m] *= pull->dims[m];
	  tmp3[m] *= pull->dims[m];
	  if (tmp[m]  < -0.5*box[m][m]) tmp[m]  += box[m][m];
	  if (tmp[m]  >  0.5*box[m][m]) tmp[m]  -= box[m][m];
	  if (tmp2[m] < -0.5*box[m][m]) tmp2[m] += box[m][m];
	  if (tmp2[m] >  0.5*box[m][m]) tmp2[m] -= box[m][m];
	  if (tmp3[m] < -0.5*box[m][m]) tmp3[m] += box[m][m];
	  if (tmp3[m] >  0.5*box[m][m]) tmp3[m] -= box[m][m];
	}
	
	if (pull->bCyl) 
	  fprintf(stderr,
		  "cur. i:%f %f %f j:%f %f %f d: %f\n"
		  "ref. i:%f %f %f j:%f %f %f d: %f\n"
		  "cor. i:%f %f %f j:%f %f %f d: %f\n",
		  rinew[i][0],rinew[i][1],rinew[i][2], 
		  rjnew[i][0],rjnew[i][1],rjnew[i][2], norm(tmp),
		  pull->pull.x_ref[i][0],pull->pull.x_ref[i][1],
		  pull->pull.x_ref[i][2],pull->dyna.x_ref[i][0],
		  pull->dyna.x_ref[i][1],pull->dyna.x_ref[i][2],
		  norm(tmp2),
		  dr[i][0],dr[i][1],dr[i][2],
		  ref_dr[0][0],ref_dr[0][1],ref_dr[0][2],
		  norm(tmp3));
	else 
	  fprintf(stderr,
		  "cur. i:%f %f %f j:%f %f %f d: %f\n"
		  "ref. i:%f %f %f j:%f %f %f d: %f\n"
		  "cor. i:%f %f %f j:%f %f %f d: %f\n",
		  rinew[i][0],rinew[i][1],rinew[i][2], 
		  rjnew[0][0],rjnew[0][1],rjnew[0][2], norm(tmp),
		  pull->pull.x_ref[i][0],pull->pull.x_ref[i][1],
		  pull->pull.x_ref[i][2],pull->ref.x_ref[0][0],
		  pull->ref.x_ref[0][1],pull->ref.x_ref[0][2],
		  norm(tmp2),
		  dr[i][0],dr[i][1],dr[i][2],
		  ref_dr[0][0],ref_dr[0][1],ref_dr[0][2],
		  norm(tmp3));
      } /* END DEBUG */

      /* update the positions with dr */
      rvec_add(rinew[i],dr[i],rinew[i]);

      if (pull->bCyl) {
	rvec_add(rjnew[i],ref_dr[i],rjnew[i]);
	
	/* calculate new distance between the two groups */
	rvec_sub(rjnew[i],rinew[i],unc_ij);
	
	/* select components and check PBC again */
	for (m=0;m<DIM;m++) {
	  unc_ij[m] *= pull->dims[m];
	  if (unc_ij[m] < -0.5*box[m][m]) unc_ij[m] += box[m][m];
	  if (unc_ij[m] >  0.5*box[m][m]) unc_ij[m] -= box[m][m];
	}
      } else {
	rvec_add(rjnew[0],ref_dr[0],rjnew[0]);

	/* calculate new distance between the two groups */
	rvec_sub(rjnew[0],rinew[i],unc_ij);

	/* select components again and check PBC again */
	for (m=0;m<DIM;m++) {
	  unc_ij[m] *= pull->dims[m];
	  if (unc_ij[m] < -0.5*box[m][m]) unc_ij[m] += box[m][m];
	  if (unc_ij[m] >  0.5*box[m][m]) unc_ij[m] -= box[m][m];
	}
      }
    }

    /* check if all constraints are fullfilled now */
    bConverged = TRUE;
    for (i=0;i<pull->pull.n;i++) {
      if (pull->bCyl) {
	rvec_sub(rjnew[i],rinew[i],unc_ij);
	rvec_sub(pull->dyna.x_ref[i],pull->pull.x_ref[i],ref_ij);
      } else {
	rvec_sub(rjnew[0],rinew[i],unc_ij);
	rvec_sub(pull->ref.x_ref[0],pull->pull.x_ref[i],ref_ij);
      }

      for (m=0;m<DIM;m++) { 
	ref_ij[m] *= pull->dims[m];
	unc_ij[m] *= pull->dims[m];
	if (unc_ij[m] < -0.5*box[m][m]) unc_ij[m] += box[m][m];
	if (unc_ij[m] >  0.5*box[m][m]) unc_ij[m] -= box[m][m];
	if (ref_ij[m] < -0.5*box[m][m]) ref_ij[m] += box[m][m];
	if (ref_ij[m] >  0.5*box[m][m]) ref_ij[m] -= box[m][m];
      }

      bConverged = bConverged && (fabs(norm(unc_ij)-norm(ref_ij)) < tol);
    }
    

    /* DEBUG */
    if (pull->bVerbose) {
      if (!bConverged)
	fprintf(stderr,"NOT CONVERGED YET: Group %d (%s):"
		"d_ref = %f, current d = %f\n",
		i,pull->pull.grps[i], norm(ref_ij),norm(unc_ij));
    } /* END DEBUG */
    

    n++;
    /* if after all constraints are dealt with and bConverged is still TRUE
       we're finished, if not we do another iteration */
  }
  if (n>max_iter) 
    fatal_error(0,"Too many iterations for constraint run");

  /* DONE ITERATING, NOW UPDATE COORDINATES AND CALC. CONSTRAINT FORCES */

  /* update the normal groups */
  for (i=0;i<pull->pull.n;i++) {
    /* get the final dr and constraint force for group i */
    rvec_sub(rinew[i],pull->pull.x_unc[i],dr[i]);
    /* select components of dr */
    for (m=0;m<DIM;m++) dr[i][m] *= pull->dims[m];
    svmul(pull->pull.tmass[i]/(dt*dt),dr[i],tmp);
    /* get the direction of dr */
    pull->pull.f[i][ZZ] = -norm(tmp)*direction[i];

    /* copy the new x_unc to x_con */
    copy_rvec(rinew[i],pull->pull.x_con[i]);
    
    /* update the atom positions */
    clear_rvec(sum);
    for (j=0;j<pull->pull.ngx[i];j++) {
      ii = pull->pull.idx[i][j];
      rvec_add(x[ii],dr[i],x[ii]);
      svmul(md->massT[ii],dr[i],tmp);
      rvec_add(tmp,sum,sum);
    }
    if (pull->bVerbose) 
      fprintf(stderr,"Group %i: correction %e %e %e\n",
	      i,sum[0],sum[1],sum[2]);
  }
  
  /* update the reference groups */
  if (pull->bCyl) {
    /* update the dynamic reference groups */
    for (i=0;i<pull->pull.n;i++) {
      rvec_sub(rjnew[i],pull->dyna.x_unc[i],ref_dr[i]);
      /* copy the new x_unc to x_con */
      copy_rvec(rjnew[i],pull->dyna.x_con[i]);
      /* select components of ref_dr */
      for (m=0;m<DIM;m++) ref_dr[i][m] *= pull->dims[m];

      clear_rvec(sum);
      for (j=0;j<pull->dyna.ngx[i];j++) {
	/* reset the atoms with dr, weighted by w_i */
	svmul(pull->dyna.weights[i][j],ref_dr[i],tmp); 
	ii = pull->dyna.idx[i][j];
	rvec_add(x[ii],tmp,x[ii]);
	svmul(md->massT[ii],tmp,tmp2);
	rvec_add(tmp2,sum,sum);
      }
      if (pull->bVerbose) 
	fprintf(stderr,"Dyna grp %i: correction %e %e %e\n",
		i,sum[0],sum[1],sum[2]);
    }
  } else { 
    /* update the reference group */
    rvec_sub(rjnew[0],pull->ref.x_unc[0], ref_dr[0]); 
    /* copy the new x_unc to x_con */
    copy_rvec(rjnew[0],pull->ref.x_con[0]);
    /* select components of ref_dr */
    for (m=0;m<DIM;m++) ref_dr[0][m] *= pull->dims[m];

    clear_rvec(sum);
    for (j=0;j<pull->ref.ngx[0];j++) {
      ii = pull->ref.idx[0][j];
      rvec_add(x[ii],ref_dr[0],x[ii]);
      svmul(md->massT[ii],ref_dr[0],tmp);
      rvec_add(tmp,sum,sum);
    }
    if (pull->bVerbose) 
      fprintf(stderr,"Reference: correction %e %e %e\n",
	      sum[0],sum[1],sum[2]);
    
  }

  /* finished! I hope. Give back some memory */
  sfree(ref_dr);
  sfree(rinew);
  sfree(rjnew);
  sfree(dr);
  sfree(direction);
}

/* mimicks an AFM experiment, groups are pulled via a spring */
static void do_afm(t_pull *pull,rvec *f,matrix box,t_mdatoms *md)
{
  int i,ii,j,m,g;
  real mi; 
  rvec cm;     /* center of mass displacement of reference */
  rvec dr;     /* extension of the springs */

  /* loop over the groups that are being pulled */
  for (i=0;i<pull->pull.n;i++) {
    /* compute how far the springs are stretched */
    for (m=0;m<DIM;m++) {
      dr[m]=pull->dims[m]*(pull->pull.spring[i][m]-pull->pull.x_unc[i][m]); 
      while (dr[m] >  box[m][m]/2) dr[m]-=box[m][m];
      while (dr[m] < -box[m][m]/2) dr[m]+=box[m][m];
    }

    /* calculate force from the spring on the pull group: f = - k*dr */
    for (m=0;m<DIM;m++)
      pull->pull.f[i][m] = pull->k*dr[m];
    
    /* distribute force on com over individual atoms in the group */
    for (j=0;j<pull->pull.ngx[i];j++) {
      ii = pull->pull.idx[i][j];
      for (m=0;m<DIM;m++) {
	mi = md->massT[ii];
	f[ii][m] +=  mi*(pull->pull.f[i][m])/(pull->pull.tmass[i]);
      }
    }
    
    /* move pulling spring along dir, over pull->rate  */
    for (m=0;m<DIM;m++) 
      pull->pull.spring[i][m]+=pull->pull.dir[i][m]*pull->rate;
  }
  /* done */
}

void pull(t_pull *pull,rvec *x,rvec *f,matrix box, t_topology *top, 
	  real dt, int step, int natoms, t_mdatoms *md) 
{
  int i;
  static rvec *x_s = NULL;

  if (!x_s)
    snew(x_s,md->nr); /* can't rely on natoms */

  /* copy x to temp array x_s. We assume all molecules are whole already */
  for (i=0;i<md->nr;i++) 
    copy_rvec(x[i],x_s[i]);  
  
  switch (pull->runtype) {
  case eAfm:
    /* calculate center of mass of the pull groups */
    for (i=0;i<pull->pull.n;i++) 
      (void)calc_com(x_s,pull->pull.ngx[i],pull->pull.idx[i],md,
		     pull->pull.x_unc[i],box);
    do_afm(pull,f,box,md);
    print_afm(pull,step);
    break;
    
  case eStart:
    for (i=0;i<pull->pull.n;i++) 
      (void)calc_com(x_s,pull->pull.ngx[i],pull->pull.idx[i],md,
		     pull->pull.x_unc[i],box);
    if (pull->bCyl)
      make_refgrps(pull,x_s,md);
    else 
      (void)calc_com(x_s,pull->ref.ngx[0],pull->ref.idx[0],md,
		     pull->ref.x_unc[0],box);
    do_start(pull,x,box,md,dt,step,top);
    print_start(pull,step);
    break; 
    
  case eConstraint:
    /* if necessary, correct for particles jumping across the box 
       this makes sure pull->ref.x0 has the pbc-corrected coordinates
       Else just copy the normal coordinates to ref.x0
     */
    if (pull->reftype == eComT0 || pull->reftype == eDynT0) {
      if (pull->bVerbose)
	fprintf(stderr,"\nCalling correct_t0_pbc\n");
      correct_t0_pbc(pull,x_s,md,box);
    } else {
      for (i=0;i<pull->ref.ngx[0];i++) 
	copy_rvec(x_s[pull->ref.idx[0][i]],pull->ref.x0[0][i]);
    }

    /* get centers of mass for the pull groups. Does this work correctly
       with pbc? */
    for (i=0;i<pull->pull.n;i++) {
      (void)calc_com(x_s,pull->pull.ngx[i],pull->pull.idx[i],md,
		     pull->pull.x_unc[i],box);
    }
    
    /* get new centers of mass for reference groups, from the, possibly 
       corrected, pull->ref.x0 */
    /* dynamic case */
    if (pull->bCyl) {
      if (step % pull->update == 0)    /* make new ones? */
	make_refgrps(pull,x_s,md);
      else {
	for (i=0;i<pull->pull.n;i++) {
	  (void)calc_com2(pull->ref.x0[0],pull->dyna.ngx[i],pull->dyna.idx[i],
			  md,pull->dyna.x_unc[i],box);
	  if (pull->bVerbose) 
	    fprintf(stderr,"dynacom: %8.3f%8.3f%8.3f\n",pull->dyna.x_unc[i][0],
		    pull->dyna.x_unc[i][1],pull->dyna.x_unc[i][2]);
	}
      }
    } 

    /* normal case */
    if (!pull->bCyl)
      (void)calc_com2(pull->ref.x0[0],pull->ref.ngx[0],pull->ref.idx[0],
		      md,pull->ref.x_unc[0],box);
    
    /* if necessary, do a running average over the last reflag steps for com */
    if (pull->reflag > 1) {
      if (pull->bVerbose) 
	fprintf(stderr,"Calling calc_running_com\n");
      calc_running_com(pull);
    }

    /* print some debug info if necessary */
    if (pull->bVerbose) {
      if (pull->bCyl) {
	for (i=0;i<pull->pull.n;i++) {
	  fprintf(stderr,"I      :%9.6f %9.6f %9.6f\n",pull->pull.x_unc[i][0],
		  pull->pull.x_unc[i][1],pull->pull.x_unc[i][2]);
	  fprintf(stderr,"dyna xref: unconstr. com:%9.6f %9.6f %9.6f\n",
		  pull->dyna.x_unc[i][0],
		  pull->dyna.x_unc[i][1],pull->dyna.x_unc[i][2]);
	}
      } else {
	fprintf(stderr,"xref: unconstr. com:%9.6f %9.6f %9.6f\n",
		pull->ref.x_unc[0][0], pull->ref.x_unc[0][1],
		pull->ref.x_unc[0][2]);
      }
    }
    
    /* do the actual constraint calculation */
    do_constraint(pull,x,box,md,dt);
    print_constraint(pull,f,step,box); 
    break;
    
  case eUmbrella:
    do_umbrella(pull,x,f,box,md);
    print_umbrella(pull,step);
    break;
    
  case eTest:
    /* code to test reference groups, without actually doing anything 
       else 
    */
    (void)calc_com(x,pull->ref.ngx[0],pull->ref.idx[0],
		   md,pull->ref.x_unc[0],box);
    fprintf(stderr,"ref: %8.3f %8.3f %8.3f\n",pull->ref.x_unc[0][XX],
	    pull->ref.x_unc[0][YY],pull->ref.x_unc[0][ZZ]);
    correct_t0_pbc(pull,x,md,box);
    (void)calc_com2(pull->ref.x0[0],pull->ref.ngx[0],pull->ref.idx[0],
		    md,pull->ref.x_unc[0],box);
    fprintf(stderr,"ref_t0: %8.3f %8.3f %8.3f\n",pull->ref.x_unc[0][XX],
	    pull->ref.x_unc[0][YY],pull->ref.x_unc[0][ZZ]);
    break;

  default:
    fatal_error(0,"undetermined runtype");
  }
}









