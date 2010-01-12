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

/* This program is NOT threadsafe, but it is only used to create
 * the nonbonded functions during the build process, so it will never be
 * executed by multiple threads.
 */

#include <mknb_common.h>
#include <mknb_declarations.h>
#include <mknb_outerloop.h>
#include <mknb_metacode.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* This program generates nonbonded interaction source code for the GROMACS
 * Molecular Dynamics package. The functions are nearly optimal, i.e.
 * there are no superfluous declarations or statements in the code.
 * The functions are generated in either C or Fortran 77. Note that some
 * special functions (free energy, generalized-born interactions) are
 * written outside this generator.
 *
 * C is somewhat more portable, but Fortran is faster on some machines.
 * There are also a lot of special options like prefetching, the software
 * version of 1/sqrt(x), thread synchronization, etc. In C we could handle
 * some of this with a lot of defines, but not all Fortran compilers
 * support preprocessing. It is also very error-prone to have 200-300
 * different functions to check - by generating them automatically any error
 * will probably show up in all functions and be fixed immediately.
 *
 */



void
mknb_write_function(void)
{
	char funcname[255];

	/* Give the routines names like nb_kernel123nf */

#ifdef IBM_FORTRAN_CPP
	sprintf(funcname,"pwr6kernel%d%d%d%s",
			mknb_func.coul,mknb_func.vdw,mknb_func.water,
			mknb_func.do_force ? "" : "nf");
#else
	sprintf(funcname,"nb_kernel%d%d%d%s",
			mknb_func.coul,mknb_func.vdw,mknb_func.water,
			mknb_func.do_force ? "" : "nf");
#endif
	
	/* Set variables we need when writing the code */
  
	/* ni is the number of atoms we work with in each 
	 * iteration of the outer loop.
	 */
	if(mknb_func.water==MKNB_WATER_SPC_SINGLE || 
	   mknb_func.water==MKNB_WATER_SPC_PAIR) {
		/* For all SPC/TIP3P loops we have 3 atoms in outer loop */
		mknb_func.ni=3; 
	} 
	else if(mknb_func.water==MKNB_WATER_TIP4P_SINGLE || 
			mknb_func.water==MKNB_WATER_TIP4P_PAIR) {
		/* For all TIP4P loops there are 4 atoms in outer loop */
		mknb_func.ni=4;
	}
	else {
		/* No water optimization */
		mknb_func.ni=1;
	}

	/* nj is the number of atoms we work with in each 
	 * iteration of the inner loop.
	 */
	if(mknb_func.water==MKNB_WATER_SPC_PAIR) {
		/* Water-water optimization for SPC/TIP3P */
		mknb_func.nj=3;
	} 
	else if(mknb_func.water==MKNB_WATER_TIP4P_PAIR) {
		/* Water-water optimization for TIP4P */
		mknb_func.nj=4;
	} 
	else {
		/* No water-water optimization.
		 * (i.e. non-water or water-other atoms kernel).
		 */
		mknb_func.nj=1;
	}

	/* LJ uses 2 parameters, Buckingham 3 */
	if(mknb_func.vdw==MKNB_VDW_BHAM)
		mknb_func.nvdw_parameters=3;
	else
		mknb_func.nvdw_parameters=2;
  
	/* Each table point needs 4 floats per interaction.
	 * Since the repulsive and dispersive parts of LJ are separate,
	 * this means we need 4 floats for coulomb-only, 8 for lj-only
	 * and 12 for LJ+coulomb tabulated interactions.
	 */
	mknb_func.table_element_size=0;

	if(mknb_func.coul==MKNB_COUL_TAB)
		mknb_func.table_element_size += 4;

	if(mknb_func.vdw==MKNB_VDW_TAB)
		mknb_func.table_element_size += 8;
  

	/* Cosmetics for C, important for Fortran */
	mknb_indent_level = 0;
  
	/* Write info about this particular kernel */
	if(!mknb_fortran)
		fprintf(mknb_output,"\n\n\n/*\n"
				" * Gromacs nonbonded kernel %s\n"
				" * Coulomb interaction:     %s\n"
				" * VdW interaction:         %s\n"
				" * water optimization:      %s\n"
				" * Calculate forces:        %s\n"
				" */\n",funcname, 
				mknb_coul_names[mknb_func.coul],
				mknb_vdw_names[mknb_func.vdw],
				mknb_water_names[mknb_func.water], 
				mknb_func.do_force ? "yes" : "no");
	else
		fprintf(mknb_output,"\n\n\nC\n"
				"C Gromacs nonbonded kernel %s\n"
				"C Coulomb interaction:     %s\n"
				"C VdW interaction:         %s\n"
				"C water optimization:      %s\n"
				"C Calculate forces:        %s\n"
				"C\n",funcname, 
				mknb_coul_names[mknb_func.coul],
				mknb_vdw_names[mknb_func.vdw],
				mknb_water_names[mknb_func.water], 
				mknb_func.do_force ? "yes" : "no");
  
	/* Write the function header and call parameters */
	mknb_function_header(funcname);

	/* Local variables */
	mknb_declare_variables();
  
	/* Initialize local data; e.g. charges for water-water */
	mknb_initialize_data();
  
	/* Start the loops over outer and inner atoms/groups */
	mknb_outerloop();
  
	mknb_finish_function();
}


void
mknb_write_file_header(void)
{

	/* Write info about an entire file 
	 * (usually two kernels, with and without forces).
	 */
	if(!mknb_fortran) {
		fprintf(mknb_output,"/*\n"
				" * Copyright (c) Erik Lindahl, David van der Spoel 2003\n"
				" * \n"
				" * This file is generated automatically at compile time\n"
				" * by the program mknb in the Gromacs distribution.\n"
				" *\n"
				" * Options used when generation this file:\n"
				" * Language:         c\n"
				" * Precision:        %s\n"
				" * Threads:          %s\n"
				" * Software invsqrt: %s\n"
				" * PowerPC invsqrt:  %s\n"
				" * Prefetch forces:  %s\n"
				" * Comments:         %s\n */\n",
				mknb_double ? "double" : "single",
				mknb_options.threads ? "yes" : "no",
				mknb_options.software_invsqrt ? "yes" : "no",
				mknb_options.ppc_invsqrt ? "yes" : "no",
				mknb_options.prefetch_forces ? "yes" : "no",
				mknb_keep_comments ? "yes" : "no");
	  
		fprintf(mknb_output,"#ifdef HAVE_CONFIG_H\n#include<config.h>\n#endif\n");

		if(mknb_options.threads)
                {
			/* thread_mpi.h must come before all other includes 
                           (except config.h) */
			fprintf(mknb_output,"#ifdef GMX_THREAD_SHM_FDECOMP\n"); 
			fprintf(mknb_output,"#include<thread_mpi.h>\n"); 
			fprintf(mknb_output,"#endif\n"); 
                }

		fprintf(mknb_output,"#include<math.h>\n");

		if(mknb_options.software_invsqrt)
			fprintf(mknb_output,"#include<vec.h>\n");

	} else {
		fprintf(mknb_output,
				"C\n"
				"C                This source code is part of\n"
				"C\n"
				"C                 G   R   O   M   A   C   S\n"
				"C\n"
				"C Copyright (c) 1991-2000, University of Groningen, The Netherlands.\n"
				"C Copyright (c) 2001-2009, The GROMACS Development Team\n"
				"C\n"
				"C Gromacs is a library for molecular simulation and trajectory analysis,\n"
				"C written by Erik Lindahl, David van der Spoel, Berk Hess, and others - for\n"
				"C a full list of developers and information, check out http://www.gromacs.org\n"
				"C\n"
				"C This program is free software; you can redistribute it and/or modify it under\n"
				"C the terms of the GNU Lesser General Public License as published by the Free\n"
				"C Software Foundation; either version 2 of the License, or (at your option) any\n" 
				"C later version.\n"
				"C As a special exception, you may use this file as part of a free software\n"
				"C library without restriction.  Specifically, if other files instantiate\n"
				"C templates or use macros or inline functions from this file, or you compile\n"
				"C this file and link it with other files to produce an executable, this\n"
				"C file does not by itself cause the resulting executable to be covered by\n"
				"C the GNU Lesser General Public License.\n"
				"C\n"
				"C In plain-speak: do not worry about classes/macros/templates either - only\n"
				"C changes to the library have to be LGPL, not an application linking with it.\n"
				"C\n"
				"C To help fund GROMACS development, we humbly ask that you cite\n"
				"C the papers people have written on it - you can find them on the website!\n"
				"C\n\n");
#ifdef IBM_FORTRAN_CPP
		fprintf(mknb_output,"#ifdef HAVE_CONFIG_H\n#  include<config.h>\n#endif\n\n");
		fprintf(mknb_output,"#ifdef GMX_DOUBLE\n");
		fprintf(mknb_output,"#  define gmxreal real*8\n");
		fprintf(mknb_output,"#else\n");
		fprintf(mknb_output,"#  define gmxreal real*4\n");
		fprintf(mknb_output,"#endif\n");
		
#endif
	}
}


int
main(int argc,char *argv[])
{
    char filename[255];
	int  i,nfiles;
	FILE *fp;
	char *ch;

	/* First set options to default values */
	mknb_fortran                  = 0; /* global variable in mknb_metacode.c */
	mknb_double                   = 0; /* global variable in mknb_metacode.c */
	mknb_keep_comments            = 0; /* global variable in mknb_metacode.c */
	mknb_options.threads          = 0; /* global variable in mknb.c */
	mknb_options.software_invsqrt = 0; /* global variable in mknb.c */
	mknb_options.ppc_invsqrt      = 0; /* global variable in mknb.c */
	mknb_options.prefetch_forces  = 0; /* global variable in mknb.c */

	fprintf(stderr,">>> Gromacs nonbonded kernel generator (-h for help)\n");

	/* Change options from command line arguments */
	for(i=1;i<argc;i++) {
		if(argv[i][0]!='-')
			fprintf(stderr,"Unknown option: %s\n",argv[i]);
		else if(argv[i][1]=='f') /* f as in fortran */
			mknb_fortran                  = 1;
		else if(argv[i][1]=='d') /* d as in double */
			mknb_double                   = 1;
		else if(argv[i][1]=='t') /* t as in threads */
			mknb_options.threads          = 1;
		else if(argv[i][1]=='s') /* s as in software_invsqrt */
			mknb_options.software_invsqrt = 1;
		else if(argv[i][1]=='p' && argv[i][2]=='p') /* pp as in ppc_invsqrt */
	        {
			/* Check if it says -ppc_invsqrt=1 , otherwise assume 2 */
			if( (ch=strchr(argv[i],'=')) != NULL && *(++ch)=='1')
			{
			    mknb_options.ppc_invsqrt = 1;
			}
			else
			{
			    mknb_options.ppc_invsqrt = 2;
			}
		}
		else if(argv[i][1]=='p') /* p as in prefetch_forces (but not pp) */
			mknb_options.prefetch_forces  = 1;
		else if(argv[i][1]=='c') /* c as in comments */
			mknb_keep_comments            = 1;
		else if(argv[i][1]=='h') { /* h as in help */
			fprintf(stderr,
					"Available options:\n"
					" -fortran           Write Fortran77 code instead of C\n"
					" -double            Use double precision iso. single\n"
					" -threads           Write kernels with thread support\n"
					" -software_invsqrt  Use Gromacs software for 1/sqrt(x)\n"
					" -ppc_invsqrt       Use PowerPC intrinsics for 1/sqrt(x)\n"
					"                    (even better: -ppc_invsqrt=1 for pwr4/ppc440/450)\n"
					" -prefetch_forces   Prefetch force memory in loops\n"
					" -comments          Write comments in output files\n\n");
			exit(0);
		}
	}


	fprintf(stderr,">>> Generating %s%s precision functions in %s.\n",
			(mknb_options.threads==1) ? "multithreaded " : "",
			(mknb_double) ? "double" : "single",
			(mknb_fortran) ? "Fortran77" : "C");
	if(mknb_options.software_invsqrt)
		fprintf(stderr,">>> Using Gromacs software version of 1/sqrt(x).\n");

	if(mknb_options.ppc_invsqrt>0)
	  fprintf(stderr,">>> Using PowerPC intrinsics for 1/sqrt(x) table (%d iterations for %s precision).\n",
		  mknb_double ? mknb_options.ppc_invsqrt+1 : mknb_options.ppc_invsqrt,
		  mknb_double ? "double" : "single");

	if(mknb_options.prefetch_forces)
		fprintf(stderr,">>> Prefetching forces in loops.\n");

 
	/* Start to write the nonbonded functions.
	 *
	 * To speed up the compile on SMP machines we 
	 * use one file for each combination of Coulomb & VdW
	 * interactions. Every file then contains the standard
	 * version of the function, and functions otimized for water.
	 * Finally, we also write versions of the functions without forces,
	 * which can be significantly faster for MC simulations.
	 */  

	nfiles=0;
	/* Loop over all combinations to construct nonbonded functions.
	 * mknb_func is a global variable structure in mknb.c.
	 */

	/* Coulomb interaction alteratives */
	for(mknb_func.coul=MKNB_COUL_NO;
		mknb_func.coul<MKNB_COUL_NR;
		mknb_func.coul = (enum mknb_coul)(mknb_func.coul + 1)) {

		/* VdW interaction alternatives */
		for(mknb_func.vdw=MKNB_VDW_NO;
			mknb_func.vdw<MKNB_VDW_NR;
			mknb_func.vdw = (enum mknb_vdw)(mknb_func.vdw + 1)) {

			/* Skip the case when we dont have any interaction at all */
			if(mknb_func.coul==MKNB_COUL_NO && mknb_func.vdw==MKNB_VDW_NO)
				continue;
			
			/* Water optimization alternatives */
			for(mknb_func.water=MKNB_WATER_NO;
				mknb_func.water<MKNB_WATER_NR;
				mknb_func.water = (enum mknb_water)(mknb_func.water + 1)) {

				/* Water optimization is useless without coulomb */
				if(mknb_func.coul==MKNB_COUL_NO && 
				   mknb_func.water!=MKNB_WATER_NO)
					continue;

				/* Only do Generalized-Born for non-water loops */
				if(mknb_func.coul==MKNB_COUL_GB && 
				   mknb_func.water!=MKNB_WATER_NO)
					continue;
	
				/* Open a new file for this function type */
#ifdef IBM_FORTRAN_CPP
				sprintf(filename,"pwr6kernel%d%d%d.F",
						mknb_func.coul,mknb_func.vdw,mknb_func.water);				
#else
				sprintf(filename,"nb_kernel%d%d%d_%s.%s",
						mknb_func.coul,mknb_func.vdw,mknb_func.water,
						(mknb_fortran) ? "f" : "c",
						(mknb_fortran) ? "f" : "c");
#endif
			
				/* Make sure that we can open it */
				if( (mknb_output = fopen(filename,"w")) == NULL ) {
					fprintf(stderr,"Error: Cannot open %s for writing.\n",
							filename);
					exit(1);
				}
        
				/* Write a header for the entire file */
				mknb_write_file_header();

				/* Write two kernels. First one with force, then without. */
				for(i=0;i<=1;i++) {
					mknb_func.do_force=(i==0);
	  
					mknb_write_function();
				}

				/* Close the file */
				fclose(mknb_output);

				/* Wrote one more without crashing - be happy! */
				nfiles++;

				/* Apparently we have 67 files in total now... */
				fprintf(stderr,"\rProgress: %2d%%",100*nfiles/67);
			}
		}
	}
	fprintf(stderr,"\n");

	/* Touch a stamp-file to show that kernels have been updated */
	fp = fopen("kernel-stamp","w");
	/* Just write something so the file is not empty, otherwise the
     * time-stamping will not work e.g. on AFS.
     */
	fprintf(fp,"kernel-stamp\n");
	fclose(fp);
	
	return 0;
} 


