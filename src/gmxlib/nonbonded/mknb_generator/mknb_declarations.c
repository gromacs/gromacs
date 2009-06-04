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
#include <stdlib.h>


#include <mknb_common.h>
#include <mknb_declarations.h>
#include <mknb_metacode.h>


void
mknb_function_header(char *funcname)
{

	/* Important - make sure the C and Fortran declarations match...
	 * The rest here is just pretty formatting that we could live without.
	 */
	if(!mknb_fortran) {

#define C_REAL  (mknb_double ? "double *" : "float *")

		fprintf(mknb_output,"void %s(\n",funcname); 
		fprintf(mknb_output,"%19s %-8s %6s p_nri,\n",       "", "int *",  "");
		fprintf(mknb_output,"%19s %-8s %6s iinr,\n",       "", "int *",   "");
		fprintf(mknb_output,"%19s %-8s %6s jindex,\n",     "", "int *",   "");
		fprintf(mknb_output,"%19s %-8s %6s jjnr,\n",       "", "int *",   "");
		fprintf(mknb_output,"%19s %-8s %6s shift,\n",      "", "int *",   "");
		fprintf(mknb_output,"%19s %-8s %6s shiftvec,\n",   "", C_REAL,    "");
		fprintf(mknb_output,"%19s %-8s %6s fshift,\n",     "", C_REAL,    "");
		fprintf(mknb_output,"%19s %-8s %6s gid,\n",        "", "int *",   "");
		fprintf(mknb_output,"%19s %-8s %6s pos,\n",        "", C_REAL,    "");
		fprintf(mknb_output,"%19s %-8s %6s faction,\n",    "", C_REAL,    "");
		fprintf(mknb_output,"%19s %-8s %6s charge,\n",     "", C_REAL,    "");
		fprintf(mknb_output,"%19s %-8s %6s p_facel,\n",     "", C_REAL,   "");
		fprintf(mknb_output,"%19s %-8s %6s p_krf,\n",       "", C_REAL,   "");
		fprintf(mknb_output,"%19s %-8s %6s p_crf,\n",       "", C_REAL,   "");
		fprintf(mknb_output,"%19s %-8s %6s Vc,\n",         "", C_REAL,    "");
		fprintf(mknb_output,"%19s %-8s %6s type,\n",       "", "int *",   "");
		fprintf(mknb_output,"%19s %-8s %6s p_ntype,\n",    "", "int *",   "");
		fprintf(mknb_output,"%19s %-8s %6s vdwparam,\n",   "", C_REAL,    "");
		fprintf(mknb_output,"%19s %-8s %6s Vvdw,\n",       "", C_REAL,    "");
		fprintf(mknb_output,"%19s %-8s %6s p_tabscale,\n",  "", C_REAL,   "");
		fprintf(mknb_output,"%19s %-8s %6s VFtab,\n",      "", C_REAL,    "");
		fprintf(mknb_output,"%19s %-8s %6s invsqrta,\n",   "", C_REAL,    "");
		fprintf(mknb_output,"%19s %-8s %6s dvda,\n",       "", C_REAL,    "");
		fprintf(mknb_output,"%19s %-8s %6s p_gbtabscale,\n","", C_REAL,   "");
		fprintf(mknb_output,"%19s %-8s %6s GBtab,\n",      "", C_REAL,    "");
		fprintf(mknb_output,"%19s %-8s %6s p_nthreads,\n",  "", "int *",  "");
		fprintf(mknb_output,"%19s %-8s %6s count,\n",      "", "int *",   "");
		fprintf(mknb_output,"%19s %-8s %6s mtx,\n",        "", "void *",  "");
		fprintf(mknb_output,"%19s %-8s %6s outeriter,\n",  "", "int *",   "");
		fprintf(mknb_output,"%19s %-8s %6s inneriter,\n",  "", "int *",   "");
		fprintf(mknb_output,"%19s %-8s %6s work)\n{",      "", C_REAL,    "");

#undef C_REAL

	} else {

		/* Fortran */

		fprintf(mknb_output,"      subroutine %s(\n",funcname);
		fprintf(mknb_output,"     &                          nri,\n");
		fprintf(mknb_output,"     &                          iinr,\n");
		fprintf(mknb_output,"     &                          jindex,\n");
		fprintf(mknb_output,"     &                          jjnr,\n");
		fprintf(mknb_output,"     &                          shift,\n");
		fprintf(mknb_output,"     &                          shiftvec,\n");
		fprintf(mknb_output,"     &                          fshift,\n");
		fprintf(mknb_output,"     &                          gid,\n");
		fprintf(mknb_output,"     &                          pos,\n");
		fprintf(mknb_output,"     &                          faction,\n");
		fprintf(mknb_output,"     &                          charge,\n");
		fprintf(mknb_output,"     &                          facel,\n");
		fprintf(mknb_output,"     &                          krf,\n");
		fprintf(mknb_output,"     &                          crf,\n");
		fprintf(mknb_output,"     &                          Vc,\n");
		fprintf(mknb_output,"     &                          type,\n");
		fprintf(mknb_output,"     &                          ntype,\n");
		fprintf(mknb_output,"     &                          vdwparam,\n");
		fprintf(mknb_output,"     &                          Vvdw,\n");
		fprintf(mknb_output,"     &                          tabscale,\n");
		fprintf(mknb_output,"     &                          VFtab,\n");
		fprintf(mknb_output,"     &                          invsqrta,\n");
		fprintf(mknb_output,"     &                          dvda,\n");
		fprintf(mknb_output,"     &                          gbtabscale,\n");
		fprintf(mknb_output,"     &                          GBtab,\n");
		fprintf(mknb_output,"     &                          nthreads,\n");
		fprintf(mknb_output,"     &                          count,\n");
		fprintf(mknb_output,"     &                          mtx,\n");
		fprintf(mknb_output,"     &                          outeriter,\n");
		fprintf(mknb_output,"     &                          inneriter,\n");
		fprintf(mknb_output,"     &                          work)\n");

		/* Declare fortran call arguments after header */

		mknb_declare_other("implicit","none");
		mknb_declare_int("nri,iinr(*),jindex(*),jjnr(*),shift(*)");
		mknb_declare_real("shiftvec(*),fshift(*),pos(*),faction(*)");
		mknb_declare_int("gid(*),type(*),ntype");
		mknb_declare_real("charge(*),facel,krf,crf,Vc(*),vdwparam(*)");
		mknb_declare_real("Vvdw(*),tabscale,VFtab(*)");
		mknb_declare_real("invsqrta(*),dvda(*),gbtabscale,GBtab(*)");

		/*
		 * mask the mutex pointer as an integer passed by
		 * reference when using fortran, or as placeholder
		 */
		mknb_declare_int("nthreads,count,mtx,outeriter,inneriter");

        /* Workspace */
        mknb_declare_real("work(*)");

	}
	fprintf(mknb_output,"\n");
}


void
mknb_finish_function()
{

	mknb_comment("Write outer/inner iteration count to pointers");
	mknb_assign( mknb_fortran ? "outeriter" : "*outeriter" , 
				 mknb_options.threads ? "nouter" : "nri");
	mknb_assign( mknb_fortran ? "inneriter" : "*inneriter" , 
				 mknb_options.threads ? "ninner" : "nj1");

	if(mknb_fortran) {
		mknb_code("return");
		mknb_code("end");
	} else {
		fprintf(mknb_output,"}");
	}
	fprintf(mknb_output,"\n\n\n");
}
  

void
mknb_declare_variables()
{
	int i,j,firsti,firstj;
	char buf[255],buf2[255];

	/* Never mind the if-statements when you are developing -
	 * they are just there to avoid tons of warnings about
	 * unused variables in the shipping code.
	 * The for-statements have the same function - we only
	 * declare what we need to make the generated code prettier.
	 */

	/* Scalar versions of arguments passed by reference */
	if(!mknb_fortran) {
		mknb_declare_int("nri,ntype,nthreads");
		mknb_declare_real("facel,krf,crf,tabscale,gbtabscale");
	}

	/* loop indices and group id*/
	mknb_declare_int("n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid");
	if(mknb_options.threads)
		mknb_declare_int("nn0,nn1,nouter,ninner");
	mknb_declare_real("shX,shY,shZ"); /* shift vectors */
	/* scalar force and vectorial force */
	if(mknb_func.do_force)             
		mknb_declare_real("fscal,tx,ty,tz");
               
	if((mknb_func.do_force && 
		(mknb_func.coul==MKNB_COUL_NORMAL || mknb_func.coul==MKNB_COUL_RF)) ||
	   mknb_func.vdw==MKNB_VDW_LJ || mknb_func.vdw==MKNB_VDW_BHAM) 
		mknb_declare_real("rinvsq");  
               
	if(mknb_func.coul) {
		if(mknb_func.water==MKNB_WATER_NO)
			mknb_declare_real("iq");
		else if(mknb_func.water==MKNB_WATER_SPC_SINGLE || 
				mknb_func.water==MKNB_WATER_TIP4P_SINGLE)
			mknb_declare_real("jq");
		mknb_declare_real("qq,vcoul,vctot");
	}
               
	if(mknb_func.vdw) {
		if(mknb_func.water!=MKNB_WATER_SPC_PAIR &&
		   mknb_func.water!=MKNB_WATER_TIP4P_PAIR)
			mknb_declare_int("nti");
		mknb_declare_int("tj");
		if(mknb_func.vdw!=MKNB_VDW_TAB)
			mknb_declare_real("rinvsix");
		mknb_declare_real("Vvdw6,Vvdwtot");
		if(mknb_func.vdw!=MKNB_VDW_BHAM)
			mknb_declare_real("Vvdw12");
	}
               
	if(mknb_func.coul==MKNB_COUL_TAB || mknb_func.coul==MKNB_COUL_GB || 
	   mknb_func.vdw==MKNB_VDW_TAB) {
		mknb_declare_real("r,rt,eps,eps2");
		mknb_declare_int("n0,nnn");
		mknb_declare_real("Y,F,Geps,Heps2,Fp,VV");
		if(mknb_func.do_force) {
			mknb_declare_real("FF");
			if(mknb_func.coul==MKNB_COUL_TAB || mknb_func.coul==MKNB_COUL_GB)
				mknb_declare_real("fijC");
			if(mknb_func.vdw==MKNB_VDW_TAB)
				mknb_declare_real("fijD,fijR");
		}
	}
	if(mknb_func.coul==MKNB_COUL_RF)
		mknb_declare_real("krsq");
	if(mknb_func.coul==MKNB_COUL_GB) {
		mknb_declare_real("isai,isaj,isaprod,gbscale,vgb");
		if(mknb_func.do_force) 
			mknb_declare_real("dvdasum,dvdatmp,dvdaj,fgb");
	}
    
	if(mknb_func.vdw==MKNB_VDW_BHAM) 
		mknb_declare_real("Vvdwexp,br");

	/* TIP4P water doesnt have any coulomb interaction
	 * on atom 1, so we skip it if we dont do LJ
	 */
	/* check for outer/i atom */
	firsti = ((mknb_func.vdw==MKNB_VDW_NO) &&
			  (mknb_func.water==MKNB_WATER_TIP4P_SINGLE || 
			   mknb_func.water==MKNB_WATER_TIP4P_PAIR)) ? 2 : 1;

	/* inner/j atom */
	firstj = ((mknb_func.vdw==MKNB_VDW_NO) && 
			  (mknb_func.water==MKNB_WATER_TIP4P_PAIR)) ? 2 : 1;
               
	/* i coordinates and forces */
	for(i=firsti;i<=mknb_func.ni;i++) {
		sprintf(buf,"ix%d,iy%d,iz%d",i,i,i);
		if(mknb_func.do_force) {
			sprintf(buf2,",fix%d,fiy%d,fiz%d",i,i,i);
			strcat(buf,buf2);
		}
		mknb_declare_real(buf);
	}
	/* j coordinates and forces */
	for(j=firstj;j<=mknb_func.nj;j++) {
		sprintf(buf,"jx%d,jy%d,jz%d",j,j,j);
		if(mknb_func.do_force && 
		   (mknb_func.water || mknb_options.prefetch_forces)) {
			if(!(mknb_func.water==MKNB_WATER_TIP4P_PAIR && j==1) || 
			   mknb_options.prefetch_forces) {
				sprintf(buf2,",fjx%d,fjy%d,fjz%d",j,j,j);
				strcat(buf,buf2);
			}
		}
		mknb_declare_real(buf);
	}
	/* i-j vectorial distance, rsq and rinv. */
	for(i=firsti;i<=mknb_func.ni;i++) {
		for(j=firstj;j<=mknb_func.nj;j++) {
			/* For TIP4p, site 1 never interacts with site 2,3,4 */
			if(mknb_func.water==MKNB_WATER_TIP4P_PAIR && 
			   ((i==1 && j>1) || (j==1 && i>1)))
				continue;
			sprintf(buf,"dx%d%d,dy%d%d,dz%d%d,rsq%d%d",
					i,j,i,j,i,j,i,j);
			if(mknb_func.coul || mknb_func.vdw!=MKNB_VDW_LJ) {
				if(!((mknb_func.water==MKNB_WATER_TIP4P_SINGLE || 
					  mknb_func.water==MKNB_WATER_TIP4P_PAIR) && 
					 i==1 && j==1 && mknb_func.vdw==MKNB_VDW_LJ)) {
					sprintf(buf2,",rinv%d%d",i,j);
					strcat(buf,buf2);
				}
			}
			mknb_declare_real(buf);
		}
	}
  
	/* The water charges and VdW parameters dont change,
	 * so we can determine them outside the mknb_func.
	 */
	if(mknb_func.water==MKNB_WATER_SPC_SINGLE)
		mknb_declare_real("qO,qH");
	else if(mknb_func.water==MKNB_WATER_TIP4P_SINGLE) 
		mknb_declare_real("qH,qM");
  
	if(mknb_func.water==MKNB_WATER_SPC_PAIR)
		mknb_declare_real("qO,qH,qqOO,qqOH,qqHH");
	if(mknb_func.water==MKNB_WATER_TIP4P_PAIR)
		mknb_declare_real("qH,qM,qqMM,qqMH,qqHH");

	if(mknb_func.vdw==MKNB_VDW_BHAM)
		mknb_declare_real("c6,cexp1,cexp2");
	else if(mknb_func.vdw!=MKNB_VDW_NO)
		mknb_declare_real("c6,c12");
  
	/* Variables needed for the inlined software inverse square root */
	if(mknb_options.software_invsqrt && 
	   (mknb_func.coul || mknb_func.vdw==MKNB_VDW_BHAM || 
		mknb_func.vdw==MKNB_VDW_TAB)) {
		mknb_declare_const_int("fractshift",12);
		mknb_declare_const_int("fractmask",8388607);
		mknb_declare_const_int("expshift",23);
		mknb_declare_const_int("expmask",2139095040);
		mknb_declare_const_int("explsb",8388608);
		mknb_declare_real4("lu");
		mknb_declare_int4("iexp,addr");
		/* To do bitwise manipulation of a FP number we need to move
		 * it back and forth between floating point and integer registers,
		 * without converting the actual data.
		 */
		if(mknb_fortran) {
			mknb_declare_int4("bval,result");
			mknb_declare_real4("fval");
			mknb_code("equivalence(bval,fval)");
			mknb_code("equivalence(result,lu)");
			mknb_declare_int4("invsqrtexptab,invsqrtfracttab");
			mknb_code("common /gmxinvsqrtdata/ invsqrtexptab(256),invsqrtfracttab(4096)");
		} else {
			mknb_declare_other("union { unsigned int bval; float fval; }",
						  "bitpattern,result");
		} 
	}
	fprintf(mknb_output,"\n");
}

void
mknb_initialize_data(void)
{   
	char buf[255];

	/* move arguments passed by reference to local scalars */
	if(!mknb_fortran) {
		mknb_assign("nri",        "*p_nri");
		mknb_assign("ntype",      "*p_ntype");
		mknb_assign("nthreads",   "*p_nthreads");
		mknb_assign("facel",      "*p_facel");
		mknb_assign("krf",        "*p_krf");
		mknb_assign("crf",        "*p_crf");
		mknb_assign("tabscale",   "*p_tabscale");
		if(mknb_func.coul==MKNB_COUL_GB)
    	{
			mknb_assign("gbtabscale", "*p_gbtabscale");
		}
	}	   

	/* assign the charge combinations for OO,OH and HH, 
     * or HH/HL/LL for TIP4P/TIP5P 
	 */
	/* we're always doing coulomb */
	if(mknb_func.water) {
		mknb_comment("Initialize water data");

		mknb_assign("ii", "%s%s", 
					mknb_array("iinr",(mknb_fortran) ? "1" : "0"), 
					(mknb_fortran) ? "+1" : "");
		switch(mknb_func.water) {
		case MKNB_WATER_SPC_SINGLE:
			mknb_assign("qO", "facel*%s", mknb_array("charge","ii"));
			mknb_assign("qH", "facel*%s", mknb_array("charge","ii+1"));
			break;
		case MKNB_WATER_TIP4P_SINGLE:
			mknb_assign("qH", "facel*%s", mknb_array("charge","ii+1"));
			mknb_assign("qM", "facel*%s", mknb_array("charge","ii+3"));
			break;
		case MKNB_WATER_SPC_PAIR:
			mknb_assign("qO", mknb_array("charge","ii"));
	 		mknb_assign("qH", mknb_array("charge","ii+1"));
			mknb_assign("qqOO","facel*qO*qO");
			mknb_assign("qqOH","facel*qO*qH");
			mknb_assign("qqHH","facel*qH*qH");
			break;
		case MKNB_WATER_TIP4P_PAIR:
			mknb_assign("qH", mknb_array("charge","ii+1"));
			mknb_assign("qM", mknb_array("charge","ii+3"));
			mknb_assign("qqMM","facel*qM*qM");
			mknb_assign("qqMH","facel*qM*qH");
			mknb_assign("qqHH","facel*qH*qH");
			break;
		default:
			printf("Error, unidentified water model (mknb_declarations.c)\n");
			exit(0);
		}

		if((mknb_func.water==MKNB_WATER_SPC_SINGLE || 
			mknb_func.water==MKNB_WATER_TIP4P_SINGLE) && mknb_func.vdw)
        {
			mknb_assign("nti","%d*ntype*%s",
                        mknb_func.nvdw_parameters,
                        mknb_array("type","ii"));
		}
        
		/* assign the nonbonded combination for the
		 * Oxygen-oxygen interactions 
		 */
		if((mknb_func.water==MKNB_WATER_SPC_PAIR || 
			mknb_func.water==MKNB_WATER_TIP4P_PAIR) && mknb_func.vdw) {
			sprintf(buf,"%d*(ntype+1)*%s%s",mknb_func.nvdw_parameters,
					mknb_array("type","ii"), (mknb_fortran) ? "+1" : "");
			mknb_assign("tj",buf);
			mknb_assign("c6",mknb_array("vdwparam","tj"));
			if(mknb_func.vdw==MKNB_VDW_BHAM) {
				mknb_assign("cexp1",mknb_array("vdwparam","tj+1"));
				mknb_assign("cexp2",mknb_array("vdwparam","tj+2"));
			} else
				mknb_assign("c12",mknb_array("vdwparam","tj+1"));
		}
		fprintf(mknb_output,"\n");
	}

	mknb_comment("Reset outer and inner iteration counters");

	if(mknb_options.threads)
	{
	    mknb_assign("nouter","0");
		mknb_assign("ninner","0");
	} 
	else
	{
		mknb_comment("Avoid compiler warning about unassigned variable");
		mknb_assign("nj1","0");
	}
}
 

