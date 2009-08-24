/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2004
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */

/* This file is NOT threadsafe, but it is only used to create
 * the nonbonded functions during the build process, so it will never be
 * executed by multiple threads.
 */

#include <string.h>
#include <mknb_common.h>
#include <mknb_innerloop.h>
#include <mknb_metacode.h>


int
mknb_load_shift_vector()
{
	mknb_comment("Load shift vector for this list");

	/* index in the shift vector for this list */
	mknb_assign("is3","3*%s%s",mknb_array("shift","n"), 
				(mknb_fortran) ? "+1" : "");

	/* load vector corresponding to this index */
	mknb_assign("shX",mknb_array("shiftvec","is3"));
	mknb_assign("shY",mknb_array("shiftvec","is3+1"));
	mknb_assign("shZ",mknb_array("shiftvec","is3+2"));

	return 0; /* no flops in this function */
}

int
mknb_zero_outer_potential()
{
	mknb_comment("Zero the potential energy for this list");
	if(mknb_func.coul)
		mknb_assign("vctot","0"); /* zero local potentials */
	if(mknb_func.vdw)
		mknb_assign("Vvdwtot","0");
	if(mknb_func.coul==MKNB_COUL_GB && mknb_func.do_force)
		mknb_assign("dvdasum","0");
  
	return 0; /* no flops in this function */
}

int
mknb_load_outer_coordinates()
{
	int i,firsti,nflops;
	char tmp[255];

	/* TIP4P water doesnt have any coulomb interaction
	 * on atom 1, so we skip it if we dont do LJ
	 */
	firsti = ((mknb_func.vdw==MKNB_VDW_NO) &&
			  (mknb_func.water==MKNB_WATER_TIP4P_SINGLE || 
			   mknb_func.water==MKNB_WATER_TIP4P_PAIR)) ? 2 : 1;
  
	mknb_comment("Load i atom data, add shift vector");
	nflops = 0;
  
	for(i=firsti;i<=mknb_func.ni;i++) {
		sprintf(tmp,"ii3+%d",3*(i-1));
		mknb_assign("ix%d","shX + %s",i,mknb_array("pos",tmp));
		sprintf(tmp,"ii3+%d",3*(i-1)+1);
		mknb_assign("iy%d","shY + %s",i,mknb_array("pos",tmp));
		sprintf(tmp,"ii3+%d",3*(i-1)+2);
		mknb_assign("iz%d","shZ + %s",i,mknb_array("pos",tmp));
		nflops += 3; /* three additions per iteration */
	}
	return nflops;
}


int
mknb_zero_outer_forces()
{
	int i,firsti;

	/* TIP4P water doesnt have any coulomb interaction
	 * on atom 1, so we skip it if we dont do LJ
	 */
	firsti = ((mknb_func.vdw==MKNB_VDW_NO) &&
			  (mknb_func.water==MKNB_WATER_TIP4P_SINGLE || 
			   mknb_func.water==MKNB_WATER_TIP4P_PAIR)) ? 2 : 1;

	mknb_comment("Clear i atom forces");
	if(mknb_func.do_force)
		for(i=firsti;i<=mknb_func.ni;i++) {
			mknb_assign("fix%d","0",i);
			mknb_assign("fiy%d","0",i);
			mknb_assign("fiz%d","0",i);
		}
	return 0; /* no flops in this function */
}


int
mknb_load_outer_parameters()
{
	/* If we do water or water-water interactions, the
	 * parameters for the outer atoms are constant and
	 * already calculated outside the loops.
	 *
	 * For normal interactions we have to load e.g. charge
	 * if we do coulomb, and the atomtype for VdW.
	 * The actual VdW parameters for the pair of atoms
	 * is looked up in the core of the inner loop, but the
	 * constant 2*ntype*type[ii] can be calculated now.
	 *
	 * For generalized-born loops we also load 1/sqrt(a)
	 */
	int nflops = 0;
  
	if(mknb_func.water==MKNB_WATER_NO) {
		mknb_comment("Load parameters for i atom");

		/* Coulomb parameters */
		if(mknb_func.coul) {
			mknb_assign("iq","facel*%s",mknb_array("charge","ii"));
			nflops++;
		}
    
		/* GB parameters - inverse sqrt of born radius */
		if(mknb_func.coul==MKNB_COUL_GB)
			mknb_assign("isai",mknb_array("invsqrta","ii"));

		/* VdW parameters */
		if(mknb_func.vdw)
			mknb_assign("nti","%d*ntype*%s",mknb_func.nvdw_parameters,
						mknb_array("type","ii"));

	}
	return nflops;
}




int
mknb_update_outer_forces()
{
	int i,firsti,nflops=0;
	char tmp[255],fsumx[255],fsumy[255],fsumz[255];

	mknb_comment("Add i forces to mem and shifted force list");

	/* TIP4P water doesnt have any coulomb interaction
	 * on atom 1, so we skip it if we dont do LJ
	 */
	firsti = ((mknb_func.vdw==MKNB_VDW_NO) &&
			  (mknb_func.water==MKNB_WATER_TIP4P_SINGLE || 
			   mknb_func.water==MKNB_WATER_TIP4P_PAIR)) ? 2 : 1;
  
	if(mknb_func.do_force) {
		fsumx[0]=fsumy[0]=fsumz[0]=0;
		for(i=firsti;i<=mknb_func.ni;i++) {
			sprintf(tmp,"ii3+%d",3*(i-1));

			mknb_assign(mknb_array("faction",tmp),"%s + fix%d",
						mknb_array("faction",tmp),i);

			sprintf(tmp,"ii3+%d",3*(i-1)+1);

			mknb_assign(mknb_array("faction",tmp),"%s + fiy%d",
						mknb_array("faction",tmp),i);

			sprintf(tmp,"ii3+%d",3*(i-1)+2);

			mknb_assign(mknb_array("faction",tmp),"%s + fiz%d",
						mknb_array("faction",tmp),i);

			sprintf(tmp,"+fix%d",i); strcat(fsumx,tmp);
			sprintf(tmp,"+fiy%d",i); strcat(fsumy,tmp);
			sprintf(tmp,"+fiz%d",i); strcat(fsumz,tmp);
			nflops += 6; /*  2*3 additions per iteration */
		}

		mknb_assign(mknb_array("fshift","is3"),"%s%s",
					mknb_array("fshift","is3"),fsumx);

		mknb_assign(mknb_array("fshift","is3+1"),"%s%s",
					mknb_array("fshift","is3+1"),fsumy);

		mknb_assign(mknb_array("fshift","is3+2"),"%s%s",
					mknb_array("fshift","is3+2"),fsumz);

	}
	return nflops;
}

int
mknb_update_outer_potential()
{
	char tmp[32];
	int nflops = 0;
  
	mknb_comment("Add potential energies to the group for this list");

	mknb_assign("ggid", "%s%s", mknb_array("gid","n"), 
				(mknb_fortran) ? "+1" : "");

	if(mknb_func.coul) {
		sprintf(tmp,"%s",mknb_array("Vc","ggid"));
		mknb_assign(tmp,"%s + vctot",tmp);
		nflops++;
	}
	if(mknb_func.vdw) {
		sprintf(tmp,"%s",mknb_array("Vvdw","ggid"));
		mknb_assign(tmp,"%s + Vvdwtot",tmp);
		nflops++;
	}
	/* Update dVda=dVda+0.5*dvdasum*(1/a) for Generalized-born.
	 * To save a couple of flops in the inner loop, each element
	 * in this list should be divided by the born radius after
	 * calling the nonbonded routine.
	 */
	if(mknb_func.coul==MKNB_COUL_GB && mknb_func.do_force) {
		mknb_assign(mknb_array("dvda","ii"),
					"%s + dvdasum*isai*isai",
					mknb_array("dvda","ii"));
		nflops++;
	}

	return nflops;
}


void
mknb_outerloop(void) {
	int i,nflops = 0;
	char tmp[255];
    int indent;
	
	
	if(mknb_options.threads) {
		mknb_comment("Loop over thread workunits");
		if(mknb_fortran) {
			char space[25];
			indent = MKNB_FORTRAN_INDENT_STEP*mknb_indent_level;
			for(i=0 ; i<indent ; i++)
				space[i]=' ';
			space[i]=0;
			fprintf(mknb_output,
					"   10 %scall f77kernelsync(mtx,count,nri,nthreads,nn0,nn1)\n",space);
			/* since f77 use call-by-reference we can send the pointer of the
			 * count variable and the mutex to c without fortran knowing about it!
			 */
			mknb_indent_level++;
			mknb_code("if(nn1.gt.nri) nn1=nri");
			mknb_comment("Start outer loop over neighborlists");
			mknb_start_loop("n", "nn0+1", "nn1");
		} else {
			/* C */
  		    mknb_code("");
			mknb_code("do");
			mknb_code("{");
			mknb_indent_level++;
			mknb_code("tMPI_Thread_mutex_lock((tMPI_Thread_mutex_t *)mtx);");
			mknb_assign("nn0","*count");
			mknb_comment("Take successively smaller chunks (at least 10 lists)");
			mknb_assign("nn1","nn0+(nri-nn0)/(2*nthreads)+10");
			/* take sucessively smaller chunks */
			mknb_assign("*count","nn1");
			mknb_code("tMPI_Thread_mutex_unlock((tMPI_Thread_mutex_t *)mtx);");
			mknb_code("if(nn1>nri) nn1=nri;");
			mknb_comment("Start outer loop over neighborlists");
			mknb_start_loop("n", "nn0", "nn1");
		}

	} else {
		mknb_comment("Start outer loop over neighborlists");
		mknb_start_loop("n", (mknb_fortran) ? "1" : "0", "nri");
	}
  
	/* load shift index and then shift vector for this list */
	nflops += mknb_load_shift_vector();

	mknb_comment("Load limits for loop over neighbors");

	mknb_assign("nj0", "%s%s", mknb_array("jindex","n"), 
				(mknb_fortran) ? "+1" : "");

	mknb_assign("nj1", mknb_array("jindex","n+1"));

	mknb_comment("Get outer coordinate index");

	mknb_assign("ii", "%s%s", mknb_array("iinr","n"), 
				(mknb_fortran) ? "+1" : "");

	mknb_assign("ii3", "3*ii%s", (mknb_fortran) ? "-2" : "");
	nflops += mknb_load_outer_coordinates();
	nflops += mknb_load_outer_parameters();
  
	nflops += mknb_zero_outer_potential();
	nflops += mknb_zero_outer_forces();

	/* do the inner loop ( separate nflops, so no return value) */
	mknb_innerloop();
	nflops += mknb_update_outer_forces();
	nflops += mknb_update_outer_potential();

	/* Update innerloop counter (for accounting) if we did threads */
	mknb_comment("Increment number of inner iterations");
	if(mknb_options.threads)
		mknb_assign("ninner","ninner + nj1 - nj0");
  
	/* close outer loop */
	sprintf(tmp,"Outer loop uses %d flops/iteration",nflops);
	mknb_comment(tmp);
	mknb_end_loop();
  
	if(mknb_options.threads) {
		mknb_comment("Increment number of outer iterations");
		mknb_assign("nouter","nouter + nn1 - nn0");
		mknb_indent_level--;
		if(mknb_fortran)
			mknb_code("if(nn1.lt.nri) goto 10");
		else
		{
		    mknb_code("}");
			mknb_code("while (nn1<nri);");
			mknb_code("");
        }
	} 
}
 
