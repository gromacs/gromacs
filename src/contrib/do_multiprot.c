/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sysstuff.h"
#include "typedefs.h"
#include "string2.h"
#include "strdb.h"
#include "macros.h"
#include "smalloc.h"
#include "mshift.h"
#include "statutil.h"
#include "copyrite.h"
#include "pdbio.h"
#include "gmx_fatal.h"
#include "xvgr.h"
#include "matio.h"
#include "index.h"
#include "gstat.h"
#include "tpxio.h"
#include "viewit.h"
#include "gbutil.h"
#include "vec.h"
#include "confio.h"
#include "gmxfio.h"

typedef struct {
    int resnr; 
    int count; 
} t_countres;

static void process_multiprot_output(const char *fn, real *rmsd, int *nres, rvec rotangles, 
				     rvec translation, bool bCountres, t_countres *countres)
{
    FILE       *mpoutput;
    char       line[256];
    char       *string;
    int        i=0,j=0,res;
    
    (*rmsd)=-1;
    (*nres)=0;
    mpoutput=ffopen (fn,"r");
    
    if (bCountres) {
	do {
	    fgets(line, 256, mpoutput);
	} while (strstr(line,"Match List") == NULL);
	fgets(line, 256, mpoutput);
	do {
	    string = strtok (line,".");
	    while (i<2 && string != NULL) {
		string = strtok (NULL,". ");
		i++;
	    }
	    i=0;
	    res=atoi(string);
	    if (res > 0) {
		while (countres[j].resnr!=res)
		    j++;
		countres[j].count++;
	    }
	    fgets(line, 256, mpoutput);
	} while (strstr(line,"End of Match List") == NULL);
	rewind(mpoutput);
    }
    do {
	fgets(line, 256, mpoutput);
    } while (strstr(line,"Trans : ") == NULL);
    
    string = strtok (line," :");
    string = strtok (NULL," ");
    while (i<3 && string != NULL) {
	string = strtok (NULL," ");
	rotangles[i]=atof(string);
	i++;
    }
    i=0;
    while (i<3 && string != NULL) {
	string = strtok (NULL," ");
	translation[i]=atof(string)/10;
	i++;
    }
    if (i!=3) {
	gmx_warning("Not enough values for rotation and translation vectors in the output of multiprot");
    }
    
    rotangles[YY]=rotangles[YY]*(-1); 
    
    while ((*rmsd) <0) {
	fgets(line, 256, mpoutput);
	if (strstr(line,"RMSD : ") != NULL) {
	    string = strtok (line,":");
	    string = strtok (NULL,":");
	    (*rmsd)=atof(string)/10;
	}
    }
    while (!(*nres)) {
	fgets(line,256, mpoutput);
	if (strstr(line,"Match List") != NULL) {
	    string = strtok (line,":");
	    string = strtok (NULL,":");
	    (*nres) = atoi(string);
	}
    }
    ffclose(mpoutput);
}

int main(int argc,char *argv[])
{
    const char *desc[] = {
	"[TT]do_multiprot[tt] ", 
	"reads a trajectory file and aligns it to a reference structure  ",
	"each time frame",
	"calling the multiprot program. This allows you to use a reference",
	"structure whose sequence is different than that of the protein in the ",
	"trajectory, since the alignment is based on the geometry, not sequence.",
	"The output of [TT]do_multiprot[tt] includes the rmsd and the number of residues",
	"on which it was calculated.",
	"[PAR]",
	"An aligned trajectory file is generated with the [TT]-ox[tt] option.[PAR]",
	"With the [TT]-cr[tt] option, the number of hits in the alignment is given",
	"per residue. This number can be between 0 and the number of frames, and",
	"indicates the structural conservation of this residue.[PAR]",
	"If you do not have the [TT]multiprot[tt] program, get it. [TT]do_multiprot[tt] assumes", 
	"that the [TT]multiprot[tt] executable is [TT]/usr/local/bin/multiprot[tt]. If this is ",
	"not the case, then you should set an environment variable [BB]MULTIPROT[bb]", 
	"pointing to the [TT]multiprot[tt] executable, e.g.: [PAR]",
	"[TT]setenv MULTIPROT /usr/MultiProtInstall/multiprot.Linux[tt][PAR]",
	"Note that at the current implementation only binary alignment (your",
	"molecule to a reference) is supported. In addition, note that the ",
	"by default [TT]multiprot[tt] aligns the two proteins on their C-alpha carbons.",
	"and that this depends on the [TT]multiprot[tt] parameters which are not dealt ",
	"with here. Thus, the C-alpha carbons is expected to give the same "
	"results as choosing the whole protein and will be slightly faster.[PAR]",
	"For information about [TT]multiprot[tt], see:",
	"http://bioinfo3d.cs.tau.ac.il/MultiProt/.[PAR]"
    };
    static bool bVerbose;
    t_pargs pa[] = {
	{ "-v",  FALSE, etBOOL, {&bVerbose},
	  "HIDDENGenerate miles of useless information" }
    };
  
    const char *bugs[] = { 
	"The program is very slow, since multiprot is run externally"
    };
  
    t_trxstatus *status;
    t_trxstatus *trxout=NULL;
    FILE        *tapein,*fo,*frc,*tmpf,*out=NULL,*fres=NULL;
    const char  *fnRef;
    const char  *fn="2_sol.res";
    t_topology  top;
    int         ePBC;
    t_atoms     *atoms,ratoms,useatoms;
    t_trxframe  fr;
    t_pdbinfo   p;
    int         nres,nres2,nr0;
    real        t;
    int         i,j,natoms,nratoms,nframe=0,model_nr=-1;
    int         cur_res,prev_res;
    int         nout;
    t_countres  *countres=NULL;
    matrix      box,rbox;
    int         gnx;
    char        *grpnm,*ss_str; 
    atom_id     *index;
    rvec        *xp,*x,*xr;
    char        pdbfile[32],refpdb[256],title[256],rtitle[256],filemode[5];
    char        out_title[256];
    char        multiprot[256],*mptr;
    int         ftp;
    int         outftp=-1;
    real        rmsd;
    bool        bTrjout,bCountres;
    const char  *TrjoutFile=NULL;
    output_env_t oenv;
    static rvec translation={0,0,0},rotangles={0,0,0};
    gmx_rmpbc_t gpbc=NULL;
    
    t_filenm   fnm[] = {
	{ efTRX, "-f",   NULL,      ffREAD },
	{ efTPS, NULL,   NULL,      ffREAD },
	{ efNDX, NULL,   NULL,      ffOPTRD },
	{ efSTX, "-r",   NULL     , ffREAD },
	{ efXVG, "-o",  "rmss",     ffWRITE },
	{ efXVG, "-rc", "rescount", ffWRITE},
	{ efXVG, "-cr", "countres", ffOPTWR},
	{ efTRX, "-ox", "aligned",  ffOPTWR }
    };
#define NFILE asize(fnm)
    
    CopyRight(stderr,argv[0]);
    parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW | PCA_TIME_UNIT | PCA_BE_NICE ,
		      NFILE,fnm, asize(pa),pa, asize(desc),desc,
		      asize(bugs),bugs,&oenv
	);
    fnRef=opt2fn("-r",NFILE,fnm);
    bTrjout = opt2bSet("-ox",NFILE,fnm);
    bCountres=  opt2bSet("-cr",NFILE,fnm);
    
    if (bTrjout) {
	TrjoutFile = opt2fn_null("-ox",NFILE,fnm);
    }
    
    read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xp,NULL,box,FALSE);
    gpbc = gmx_rmpbc_init(&top.idef,ePBC,top.atoms.nr,box);
    atoms=&(top.atoms);

    ftp=fn2ftp(fnRef);
 
    get_stx_coordnum(fnRef,&nratoms);
    init_t_atoms(&ratoms,nratoms,TRUE);  
    snew(xr,nratoms);
    read_stx_conf(fnRef,rtitle,&ratoms,xr,NULL,&ePBC,rbox);
    
    if (bVerbose) {
	fprintf(stderr,"Read %d atoms\n",atoms->nr); 
	fprintf(stderr,"Read %d reference atoms\n",ratoms.nr); 
    }
    if (bCountres) {
	snew(countres,ratoms.nres);
	j=0;
	cur_res=0;
	for (i=0;i<ratoms.nr;i++) {
	    prev_res=cur_res;
	    cur_res=ratoms.atom[i].resind;
	    if (cur_res != prev_res) {
		countres[j].resnr=cur_res;
		countres[j].count=0;
		j++;
	    }
	}
    }
    get_index(atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&gnx,&index,&grpnm);
    nres=0;
    nr0=-1;
    for(i=0; (i<gnx); i++) {
	if (atoms->atom[index[i]].resind != nr0) {
	    nr0=atoms->atom[index[i]].resind;
	    nres++;
	}
    }
    fprintf(stderr,"There are %d residues in your selected group\n",nres);
    
    strcpy(pdbfile,"ddXXXXXX");
    gmx_tmpnam(pdbfile);
    if ((tmpf = fopen(pdbfile,"w")) == NULL) {
	sprintf(pdbfile,"%ctmp%cfilterXXXXXX",DIR_SEPARATOR,DIR_SEPARATOR);
	gmx_tmpnam(pdbfile);
	if ((tmpf = fopen(pdbfile,"w")) == NULL) {
	    gmx_fatal(FARGS,"Can not open tmp file %s",pdbfile);
	}
    }
    else {
	ffclose(tmpf);
    }

    if (ftp != efPDB) {
	strcpy(refpdb,"ddXXXXXX");
	gmx_tmpnam(refpdb);
	strcat(refpdb,".pdb");
	write_sto_conf(refpdb,rtitle,&ratoms,xr,NULL,ePBC,rbox);
    }
    else {
	strcpy(refpdb,fnRef);
    }

    if ((mptr=getenv("MULTIPROT")) == NULL) {
	mptr="/usr/local/bin/multiprot";
    }
    if (!gmx_fexist(mptr)) {
	gmx_fatal(FARGS,"MULTIPROT executable (%s) does not exist (use setenv MULTIPROT)",
		  mptr);
    }
    sprintf (multiprot,"%s %s %s > /dev/null %s",
	     mptr, refpdb, pdbfile, "2> /dev/null");
    
    if (bVerbose)
	fprintf(stderr,"multiprot cmd='%s'\n",multiprot);
    
    if (!read_first_frame(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&fr,TRX_READ_X)) 
      	gmx_fatal(FARGS,"Could not read a frame from %s",ftp2fn(efTRX,NFILE,fnm));
    natoms = fr.natoms;

    if (bTrjout) {
	nout=natoms;
	/* open file now */
	outftp=fn2ftp(TrjoutFile);
	if (bVerbose)
	    fprintf(stderr,"Will write %s: %s\n",ftp2ext(ftp),ftp2desc(outftp));
	strcpy(filemode,"w");
	switch (outftp) {
	    case efXTC:
	    case efG87:
	    case efTRR:
	    case efTRJ:
		out=NULL;
		trxout = open_trx(TrjoutFile,filemode);
		break;
	    case efGRO:
	    case efG96:
	    case efPDB:
		/* Make atoms struct for output in GRO or PDB files */
		/* get memory for stuff to go in pdb file */
		init_t_atoms(&useatoms,nout,FALSE);
		sfree(useatoms.resinfo);
		useatoms.resinfo=atoms->resinfo;
		for(i=0;(i<nout);i++) {
		    useatoms.atomname[i]=atoms->atomname[i];
		    useatoms.atom[i]=atoms->atom[i];
		    useatoms.nres=max(useatoms.nres,useatoms.atom[i].resind+1);
		}
		useatoms.nr=nout;
		out=ffopen(TrjoutFile,filemode);
		break;
	}
	if (outftp == efG87)
	    fprintf(gmx_fio_getfp(trx_get_fileio(trxout)),"Generated by %s. #atoms=%d, a BOX is"
		    " stored in this file.\n",Program(),nout);   
    }
    
    if (natoms > atoms->nr) {
	gmx_fatal(FARGS,"\nTrajectory does not match topology!");
    }
    if (gnx > natoms) {
	gmx_fatal(FARGS,"\nTrajectory does not match selected group!");
    }

    fo = xvgropen(opt2fn("-o",NFILE,fnm),"RMSD","Time (ps)","RMSD (nm)",oenv);
    frc = xvgropen(opt2fn("-rc",NFILE,fnm),"Number of Residues in the alignment","Time (ps)","Residues",oenv);
    
    do {
	t = output_env_conv_time(oenv,fr.time);
	gmx_rmpbc(gpbc,natoms,fr.box,fr.x);
	tapein=ffopen(pdbfile,"w");
	write_pdbfile_indexed(tapein,NULL,atoms,fr.x,ePBC,fr.box,' ',-1,gnx,index,NULL,TRUE); 
	ffclose(tapein);
	system(multiprot);
	remove(pdbfile);
	process_multiprot_output(fn, &rmsd, &nres2,rotangles,translation,bCountres,countres);
	fprintf(fo,"%12.7f",t);
	fprintf(fo," %12.7f\n",rmsd);
	fprintf(frc,"%12.7f",t);
	fprintf(frc,"%12d\n",nres2);
	if (bTrjout) {
	    rotate_conf(natoms,fr.x,NULL,rotangles[XX],rotangles[YY],rotangles[ZZ]);
	    for(i=0; i<natoms; i++) {
		rvec_inc(fr.x[i],translation);
	    }
	    switch(outftp) {
		case efTRJ:
		case efTRR:
		case efG87:
		case efXTC:
		    write_trxframe(trxout,&fr,NULL);
		    break;
		case efGRO:
		case efG96:
		case efPDB:
		    sprintf(out_title,"Generated by do_multiprot : %s t= %g %s",
			    title,output_env_conv_time(oenv,fr.time),output_env_get_time_unit(oenv));
		    switch(outftp) {
			case efGRO: 
			    write_hconf_p(out,out_title,&useatoms,prec2ndec(fr.prec),
					  fr.x,NULL,fr.box);
			    break;
			case efPDB:
			    fprintf(out,"REMARK    GENERATED BY DO_MULTIPROT\n");
			    sprintf(out_title,"%s t= %g %s",title,output_env_conv_time(oenv,fr.time),output_env_get_time_unit(oenv));
			    /* if reading from pdb, we want to keep the original 
			       model numbering else we write the output frame
			       number plus one, because model 0 is not allowed in pdb */
			    if (ftp==efPDB && fr.step > model_nr) {
				model_nr = fr.step;
			    }
			    else {
				model_nr++;
			    }
			    write_pdbfile(out,out_title,&useatoms,fr.x,ePBC,fr.box,' ',model_nr,NULL,TRUE);
			    break;
			case efG96:
			    fr.title = out_title;
			    fr.bTitle = (nframe == 0);
			    fr.bAtoms = FALSE;
			    fr.bStep = TRUE;
			    fr.bTime = TRUE;
			    write_g96_conf(out,&fr,-1,NULL);
		    }
		    break;
	    }
	}
	nframe++;
    } while(read_next_frame(oenv,status,&fr));
    if (bCountres) {
	fres=  xvgropen(opt2fn("-cr",NFILE,fnm),"Number of frames in which the residues are aligned to","Residue","Number",oenv);
	for (i=0;i<ratoms.nres;i++) {
	    fprintf(fres,"%10d  %12d\n",countres[i].resnr,countres[i].count);
	}
	ffclose(fres);
    }
    ffclose(fo);
    ffclose(frc);
    fprintf(stderr,"\n");
    close_trj(status);
    if (trxout != NULL) {
	close_trx(trxout);
    }
    else if (out != NULL) {
	ffclose(out);
    }
    view_all(oenv,NFILE, fnm);
    sfree(xr);
    if (bCountres) {
	sfree(countres);
    }
    free_t_atoms(&ratoms,TRUE);
    if (bTrjout) {
	if (outftp==efPDB || outftp==efGRO || outftp==efG96) {
	    free_t_atoms(&useatoms,TRUE);
	}
    }
    thanx(stderr);
    return 0;
}
