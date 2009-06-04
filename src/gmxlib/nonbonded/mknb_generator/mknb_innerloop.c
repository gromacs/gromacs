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
#include <mknb_common.h>
#include <mknb_interactions.h>
#include <mknb_metacode.h>

int
ppc_invsqrt(char *rsq, char *rinv)
{
	int nflops = 1;

	mknb_comment("PowerPC intrinsics 1/sqrt lookup table");

	/* fsqrte step */

	if (mknb_double) {
		if (mknb_fortran) {
			mknb_assign (rinv,"frsqrte(%s)",rsq);
		} else {
			mknb_assign (rinv,"__frsqrte(%s)",rsq);
		}
	} else {
		/* Frsqrtes is only supported on Power5 and higher, but no use in optimizing old HW! */
		if (mknb_fortran) {
			mknb_assign (rinv,"frsqrtes(%s)",rsq);
		} else {
			mknb_assign (rinv,"__frsqrtes(%s)",rsq);
		}
	}

	/* Newton-Rhapson iteration step */    
	mknb_assign(rinv,"(0.5*%s*(3.0-((%s*%s)*%s)))",rinv,rsq,rinv,rinv);
	nflops += 5; /* 4 mult and one sub on the last line */
	
	if(mknb_options.ppc_invsqrt==2)
	{
		/* Older powerpc architectures need two iterations for single, 3 for double */
		mknb_assign(rinv,"(0.5*%s*(3.0-((%s*%s)*%s)))",rinv,rsq,rinv,rinv);
		nflops += 5;
	}
	
	if(mknb_double) {
		mknb_assign(rinv,"(0.5*%s*(3.0-((%s*%s)*%s)))",rinv,rsq,rinv,rinv);
		nflops += 5; /* 4 mult and one sub on the last line */
	}    
#ifdef IBM_FORTRAN_CPP
	fprintf(mknb_output,"#ifdef GMX_DOUBLE\n");
	mknb_assign(rinv,"(0.5*%s*(3.0-((%s*%s)*%s)))",rinv,rsq,rinv,rinv);
	fprintf(mknb_output,"#endif\n");
#endif
	return nflops;
}


int
software_invsqrt(char *rsq, char *rinv)
{
    int nflops = 0;

    mknb_comment("Gromacs software 1/sqrt");
    /* table lookup step */
    mknb_assign ( (mknb_fortran) ? "fval" : "bitpattern.fval", rsq);

    if(mknb_fortran) {
        mknb_assign("iexp","rshift(and(bval,expmask),expshift)");
        mknb_assign("addr","rshift(and(bval,or(fractmask,explsb)),fractshift)");
        mknb_assign("result","or(invsqrtexptab(iexp+1),invsqrtfracttab(addr+1))");
    } else {
        mknb_assign("iexp","(((bitpattern.bval)&expmask)>>expshift)");
        mknb_assign("addr","(((bitpattern.bval)&(fractmask|explsb))>>fractshift)");

        mknb_assign("result.bval",
                    "gmx_invsqrt_exptab[iexp] | gmx_invsqrt_fracttab[addr]");
        mknb_assign("lu", "result.fval");
    }

    /* Newton-Rhapson iteration step */
    mknb_assign(rinv,"(0.5*lu*(3.0-((%s*lu)*lu)))",rsq);
    nflops += 5; /* 4 mult and one sub on the last line */
    if(mknb_double) {
        mknb_assign(rinv,"(0.5*%s*(3.0-((%s*%s)*%s)))",rinv,rsq,rinv,rinv);
        nflops += 5; /* 4 mult and one sub on the last line */
    }
    return nflops;
}


int
mknb_load_inner_coordinates()
{
	int j,firstj;
	char tmp[255];

	/* TIP4P water doesnt have any coulomb interaction
	 * on atom 1, so we skip it if we dont do LJ
	 */
	firstj = ((mknb_func.vdw==MKNB_VDW_NO) && 
			  (mknb_func.water==MKNB_WATER_TIP4P_PAIR)) ? 2 : 1;
  
	mknb_comment("load j atom coordinates");
	for(j=firstj;j<=mknb_func.nj;j++) {
		sprintf(tmp,"j3+%d",3*(j-1));
		mknb_assign("jx%d","%s",j,mknb_array("pos",tmp));
		sprintf(tmp,"j3+%d",3*(j-1)+1);
		mknb_assign("jy%d","%s",j,mknb_array("pos",tmp));
		sprintf(tmp,"j3+%d",3*(j-1)+2);
		mknb_assign("jz%d","%s",j,mknb_array("pos",tmp));
	}
	/* Only assignment, no flops */
	return 0;
}

int
mknb_calc_distance()
{
	int i,j,firsti,firstj,nflops=0;
  
	/* TIP4P water doesnt have any coulomb interaction
	 * on atom 1, so we skip it if we dont do LJ
	 */
	if(mknb_func.vdw==MKNB_VDW_NO) {
		firsti = (mknb_func.water==MKNB_WATER_TIP4P_SINGLE || 
				  mknb_func.water==MKNB_WATER_TIP4P_PAIR) ? 2 : 1;
		firstj = (mknb_func.water==MKNB_WATER_TIP4P_PAIR) ? 2 : 1;
	} else {
		firsti = 1;
		firstj = 1;
	}

	mknb_comment("Calculate distance");
	for(i=firsti;i<=mknb_func.ni;i++)
		for(j=firstj;j<=mknb_func.nj;j++) {
			/* For TIP4p, site 1 never interacts with site 2,3,4 */
			if(mknb_func.water==MKNB_WATER_TIP4P_PAIR && 
			   ((i==1 && j>1) || (j==1 && i>1)))
				continue;
			mknb_assign("dx%d%d","ix%d - jx%d",i,j,i,j);
			mknb_assign("dy%d%d","iy%d - jy%d",i,j,i,j);
			mknb_assign("dz%d%d","iz%d - jz%d",i,j,i,j);
			mknb_assign("rsq%d%d","dx%d%d*dx%d%d+dy%d%d*dy%d%d+dz%d%d*dz%d%d",
						i,j,i,j,i,j,i,j,i,j,i,j,i,j);
			
			/* three subtractions, two adds and three multiplications */
			nflops += 8; 
		}
    return nflops;
}


int
mknb_prefetch_inner_forces()
{
	int j,firstj,nflops=0;
	char tmp[255];

	/* TIP4P water doesnt have any coulomb interaction
	 * on atom 1, so we skip it if we dont do LJ
	 */
	firstj = ((mknb_func.vdw==MKNB_VDW_NO) && 
			  (mknb_func.water==MKNB_WATER_TIP4P_PAIR)) ? 2 : 1;

	if(mknb_func.do_force) {
		mknb_comment("prefetch forces");

		for(j=firstj;j<=mknb_func.nj;j++) {
			sprintf(tmp,"j3+%d",3*(j-1));
			mknb_assign("fjx%d","%s",j,mknb_array("faction",tmp));
			sprintf(tmp,"j3+%d",3*(j-1)+1);
			mknb_assign("fjy%d","%s",j,mknb_array("faction",tmp));
			sprintf(tmp,"j3+%d",3*(j-1)+2);
			mknb_assign("fjz%d","%s",j,mknb_array("faction",tmp));
		}
	}
	return nflops;
}


int
mknb_load_inner_parameters(int i,int j)
{
	int nflops=0;
	char idx[255];
  
	/* For water-water functions all parameters are constant,
	 * so we only assign them to the right variable. In
	 * the other cases we also have to load the data.
	 */
	mknb_comment("Load parameters for j atom");
  
	/* Coulomb parameters */
	if(mknb_func.coul) {
		if(mknb_func.water==MKNB_WATER_NO) {
			if(mknb_func.coul==MKNB_COUL_GB) {
				/* Generalized born: load 1/sqrt(a) */
				mknb_assign("isaj",mknb_array("invsqrta","jnr"));
				mknb_assign("isaprod","isai*isaj");
				
				/* GB-PL */
				mknb_assign("qq","iq*%s",mknb_array("charge","jnr"));
				mknb_assign("vcoul","qq*rinv%d%d",i,j);
				
				if(mknb_func.do_force)
					mknb_assign("fscal","vcoul*rinv%d%d",i,j);
				
				/* Save a flop by multiplying qq with isa1*isa2 already here */
				mknb_assign("qq","isaprod*(-qq)");
				mknb_assign("gbscale","isaprod*gbtabscale");
				nflops+=4;
			} else {
				/* Load normal (non-GB) charges */
				mknb_assign("qq","iq*%s",mknb_array("charge","jnr"));
				nflops++;
			}
		} else if(mknb_func.water==MKNB_WATER_SPC_SINGLE) {
			sprintf(idx,"jnr+%d",j-1);
			if(i==1) 
				mknb_assign("jq", mknb_array("charge",idx));
			if(i==1 || i==2) {
				mknb_assign("qq","%s*jq", (i==1) ? "qO" : "qH");
				nflops++;
			}
		} else if(mknb_func.water==MKNB_WATER_TIP4P_SINGLE) {
			sprintf(idx,"jnr+%d",j-1);
			if(i==2)
				mknb_assign("jq", mknb_array("charge",idx));
			if(i==2 || i==4) {
				mknb_assign("qq","%s*jq", (i==4) ? "qM" : "qH");
				nflops++;
			}
		} else if(mknb_func.water==MKNB_WATER_SPC_PAIR) {
			if(i==1 && j==1)
				mknb_assign("qq","qqOO");
			else if(i>1 && j>1)
				mknb_assign("qq","qqHH");
			else
				mknb_assign("qq","qqOH");
		} else if(mknb_func.water==MKNB_WATER_TIP4P_PAIR) {
			if(i==1 || j==1) {
				/* Do nothing - site 1 is LJ-only in TIP4P */
			} else if(i==4 && j==4)
				mknb_assign("qq","qqMM");
			else if(i<4 && j<4)
				mknb_assign("qq","qqHH");
			else
				mknb_assign("qq","qqMH");
		}
	}

	/* VdW parameters */
	if(mknb_func.vdw) {
		if(mknb_func.water==MKNB_WATER_NO ||
		   ((mknb_func.water==MKNB_WATER_SPC_SINGLE || 
			 mknb_func.water==MKNB_WATER_TIP4P_SINGLE) && (i==1))) {

			mknb_assign("tj","nti+%d*%s%s",mknb_func.nvdw_parameters,
						mknb_array("type","jnr"), (mknb_fortran) ? "+1" : "");

			mknb_assign("c6",mknb_array("vdwparam","tj"));

			if(mknb_func.vdw==MKNB_VDW_BHAM) {
				mknb_assign("cexp1",mknb_array("vdwparam","tj+1"));
				mknb_assign("cexp2",mknb_array("vdwparam","tj+2"));
			} else {
				mknb_assign("c12",mknb_array("vdwparam","tj+1"));
			}
		}
		/* For water-water interactions the LJ parameters are constant
		 * and calculate outside the loops.
		 */
	}
	return nflops;
}


int
mknb_calc_square_root()
{
	int i,j,firsti,firstj,nflops=0;
	char tmp[255],tmp2[255];
  
	/* TIP4P water doesnt have any coulomb interaction
	 * on atom 1, so skip it completely if we dont do LJ
	 */
	if(mknb_func.vdw==MKNB_VDW_NO) {
		firsti = (mknb_func.water==MKNB_WATER_TIP4P_SINGLE || 
				  mknb_func.water==MKNB_WATER_TIP4P_PAIR) ? 2 : 1;
		firstj = (mknb_func.water==MKNB_WATER_TIP4P_PAIR) ? 2 : 1;
	} else {
		firsti = 1;
		firstj = 1;
	}
  
	mknb_comment("Calculate 1/r and 1/r2");
	for(i=firsti;i<=mknb_func.ni;i++) {
		for(j=firstj;j<=mknb_func.nj;j++) {
			/* For TIP4p, site 1 never interacts with site 2,3,4 */
			if(mknb_func.water==MKNB_WATER_TIP4P_PAIR && 
			   ((i==1 && j>1) || (j==1 && i>1)))
				continue;      
			/* For LJ-only interactions we only need 1/rsq.
			 * This is always true for the LJ-only site in TIP4P.
			 */
			if( mknb_func.vdw==MKNB_VDW_LJ &&
				(mknb_func.coul==MKNB_COUL_NO || 
				 ((i==1) && (mknb_func.water==MKNB_WATER_TIP4P_SINGLE || 
							 mknb_func.water==MKNB_WATER_TIP4P_PAIR)))) {
				mknb_assign("rinvsq","1.0/rsq%d%d",i,j,i,j);
				nflops+=4; /* Estimate 1/x to 4 flops (SSE value) */
			} else {
				if(mknb_options.software_invsqrt) {
					sprintf(tmp,"rsq%d%d",i,j);
					sprintf(tmp2,"rinv%d%d",i,j);
					nflops += software_invsqrt(tmp,tmp2);
				} else if (mknb_options.ppc_invsqrt) {
					sprintf(tmp,"rsq%d%d",i,j);
					sprintf(tmp2,"rinv%d%d",i,j);
					nflops += ppc_invsqrt(tmp,tmp2);
				} else {
					mknb_assign("rinv%d%d","1.0/sqrt(rsq%d%d)",i,j,i,j);
					/* Estimate 1/sqrt(x) to 5 flops in single, 10 in double */
					nflops += mknb_double ? 10 : 5;
				}
			}
		}
	}
	return nflops;
}

void mknb_innerloop()
{
	int i,j,firsti,firstj;
	int nflops = 0;
	int vdwsave,coulsave,read_from_mem,write_to_mem;
	char tmp[255],mem[255],var[255],rsq[255],rinv[255];

	
	mknb_start_loop("k", "nj0" ,"nj1");

	mknb_comment("Get j neighbor index, and coordinate index");
	mknb_assign("jnr", "%s%s", mknb_array("jjnr","k"), 
				(mknb_fortran) ? "+1" : "");
	mknb_assign("j3","3*jnr%s", (mknb_fortran) ? "-2" : "");
  
	/* Load j particle coordinates */
	nflops += mknb_load_inner_coordinates();

	nflops += mknb_calc_distance();
  
	/* calculate inverse square root, except when we only do LJ coulomb -
	 * in that case we only need 1/rsq and can use a faster division.
	 * This also applies to the first atom in TIP4P water, which
	 * only has LJ.   
	 */
	nflops += mknb_calc_square_root();

	if(mknb_func.do_force && mknb_options.prefetch_forces)
		nflops += mknb_prefetch_inner_forces();

	/* TIP4P water doesnt have any coulomb interaction
	 * on atom 1, so skip it completely if we dont do LJ
	 */
	if(mknb_func.vdw==MKNB_VDW_NO) {
		firsti = (mknb_func.water==MKNB_WATER_TIP4P_SINGLE || 
				  mknb_func.water==MKNB_WATER_TIP4P_PAIR) ? 2 : 1;
		firstj = (mknb_func.water==MKNB_WATER_TIP4P_PAIR) ? 2 : 1;
	} else {
		firsti = 1;
		firstj = 1;
	}
  
	/* Do each interaction and update forces directly. */

	for(i=firsti;i<=mknb_func.ni;i++) {
		for(j=firstj;j<=mknb_func.nj;j++) {
			/* For TIP4P, site 1 never interacts with site 2,3,4 */
			if(mknb_func.water==MKNB_WATER_TIP4P_PAIR && 
			   ((i==1 && j>1) || (j==1 && i>1)))
				continue;

			coulsave=mknb_func.coul;
			vdwsave=mknb_func.vdw;

			/* Do not do VdW for atoms 2,3,(4) in waters */
			if(i>1 || j>1) 
				mknb_func.vdw=MKNB_VDW_NO;
      
			/* Do not do coulomb for atom 1 on TIP4P waters */
			if(((mknb_func.water==MKNB_WATER_TIP4P_SINGLE || 
				 mknb_func.water==MKNB_WATER_TIP4P_PAIR) && i==1) ||
			   (mknb_func.water==MKNB_WATER_TIP4P_PAIR && j==1)) 
				mknb_func.coul=MKNB_COUL_NO;
         
			/* load j atom parameters, and calculate charge products, etc. */
			nflops += mknb_load_inner_parameters(i,j);

			sprintf(rsq,"rsq%d%d",i,j);
			sprintf(rinv,"rinv%d%d",i,j);

			nflops += mknb_calculate_interaction(rsq,rinv);

			mknb_func.coul=coulsave;
			mknb_func.vdw=vdwsave;
                
			/* now we have the scalar force - update the forces if necessary */
			if(mknb_func.do_force) {
				mknb_comment("Calculate temporary vectorial force");
				mknb_assign("tx","fscal*dx%d%d",i,j);
				mknb_assign("ty","fscal*dy%d%d",i,j);
				mknb_assign("tz","fscal*dz%d%d",i,j);
				nflops += 3;
				/* increment i force */
				mknb_comment("Increment i atom force");
				mknb_assign("fix%d","fix%d + tx",i,i);
				mknb_assign("fiy%d","fiy%d + ty",i,i);
				mknb_assign("fiz%d","fiz%d + tz",i,i);
				nflops += 3;
				/* decrement j force.
				 * The first time we access it we load from memory,
				 * unless we already prefetched the forces above. 
				 * The last time the result is stored to memory.
				 */
				mknb_comment("Decrement j atom force");
				/* Read from memory if we did not prefetch, if one of this is true:
				 * Non-water
				 * SPC_SINGLE, and i==1
				 * TIP4P_SINGLE, and i==1, or i==2 if no LJ
				 * SPC_PAIR, and i==1
				 * TIP4P_PAIR, and i==1, or i==2.
				 */
				read_from_mem = !mknb_options.prefetch_forces &&
					((mknb_func.water==MKNB_WATER_NO) ||
					 (mknb_func.water==MKNB_WATER_SPC_SINGLE && i==1) ||
					 (mknb_func.water==MKNB_WATER_TIP4P_SINGLE && 
					  (i==1 || (i==2 && mknb_func.vdw==MKNB_VDW_NO))) ||
					 (mknb_func.water==MKNB_WATER_SPC_PAIR && i==1) ||
					 (mknb_func.water==MKNB_WATER_TIP4P_PAIR && (i==1 || i==2)));

				/* Write the force directly to memory in these cases:
				 * Non-water
				 * SPC_SINGLE, and i==3
				 * TIP4P_SINGLE, and i==1 or i==4.
				 * SPC_PAIR, and i==3
				 * TIP4P_PAIR, and i==1, or i==4.
				 */
				write_to_mem =
					(mknb_func.water==MKNB_WATER_NO) ||
					(mknb_func.water==MKNB_WATER_SPC_SINGLE && i==3) ||
					(mknb_func.water==MKNB_WATER_TIP4P_SINGLE && i==4) ||
					(mknb_func.water==MKNB_WATER_SPC_PAIR && i==3) ||
					(mknb_func.water==MKNB_WATER_TIP4P_PAIR && (i==1 || i==4));
        
				sprintf(tmp,"j3+%d",3*(j-1));
				sprintf(mem,"%s",mknb_array("faction",tmp));
				sprintf(var,"fjx%d",j);

				mknb_assign( write_to_mem ? mem : var, "%s - tx",
							 read_from_mem ? mem : var);

				sprintf(tmp,"j3+%d",3*(j-1)+1);
				sprintf(mem,"%s",mknb_array("faction",tmp));
				sprintf(var,"fjy%d",j);

				mknb_assign( write_to_mem ? mem : var, "%s - ty",
							 read_from_mem ? mem : var);

				sprintf(tmp,"j3+%d",3*(j-1)+2);
				sprintf(mem,"%s",mknb_array("faction",tmp));
				sprintf(var,"fjz%d",j);

				mknb_assign( write_to_mem ? mem : var, "%s - tz",
							 read_from_mem ? mem : var);

				nflops += 3; /* Decrementing j forces */
			}
		}
	}
	sprintf(tmp,"Inner loop uses %d flops/iteration",nflops);
	mknb_comment(tmp);
	
	mknb_end_loop();
}

