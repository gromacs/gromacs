/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: gentop.c,v 1.26 2009/05/20 10:48:03 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include "maths.h"
#include "macros.h"
#include "copyrite.h"
#include "bondf.h"
#include "string2.h"
#include "smalloc.h"
#include "strdb.h"
#include "sysstuff.h"
#include "confio.h"
#include "physics.h"
#include "statutil.h"
#include "vec.h"
#include "random.h"
#include "3dview.h"
#include "txtdump.h"
#include "readinp.h"
#include "names.h"
#include "filenm.h"
#include "pbc.h"
#include "pdbio.h"
#include "vec.h"
#include "gmx_random.h"
#include "pdb2top.h"
#include "toputil.h"
#include "slater_integrals.h"
#include "gpp_atomtype.h"
#include "gmx_resp.h"
#include "gentop_vsite.h"
#include "gentop_nm2type.h"
#include "gentop_core.h"
#include "gentop_qgen.h"
#include "atomprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "poldata.h"
#include "poldata_xml.h"

enum { edihNo, edihOne, edihAll, edihNR };

static int get_option(const char **opts)
{
    int val = 0;
  
    if (!opts)
        return NOTSET;
    if (opts[val] != NULL)
        for(val=1; (opts[val] != NULL); val++)
            if (strcasecmp(opts[0],opts[val]) == 0)
                break;
    if (opts[val] == NULL)
        val = 0;
    else
        val--;
    
    return val;
}

static void clean_pdb_names(t_atoms *atoms,t_symtab *tab)
{
    int  i,changed;
    char *ptr,buf[128];
  
    for(i=0; (i<atoms->nr); i++) {
        changed = 0;
        strncpy(buf,*(atoms->atomname[i]),sizeof(buf));
        while ((ptr = strchr(buf,' ')) != NULL) {
            *ptr = '_';
            changed = 1;
        }
        if (changed)
            atoms->atomname[i] = put_symtab(tab,buf);
    }
}

static void print_qpol(t_atoms *atoms,char **smnames,gmx_poldata_t pd)
{
    int    i,np;
    double poltot,pol,sigpol,qtot,sptot;
    char   *gt_type;
    
    poltot = 0;
    sptot  = 0;
    np     = 0;
    for(i=0; (i<atoms->nr); i++) 
    {
        gt_type = gmx_poldata_get_gt_type(pd,smnames[i]);  
        if ((NULL != gt_type) && 
            gmx_poldata_gt_type_polarizability(pd,gt_type,&pol,&sigpol))
        {
            np++;
            poltot += pol;
            sptot  += sqr(sigpol);
        }
    }
    printf("Polarizability is %g +/- %g A^3.\n",poltot,sqrt(sptot/atoms->nr));
}

void put_in_box(int natom,matrix box,rvec x[],real dbox)
{
    int  i,m;
    real xx,yy,zz;
    rvec invxyz;
    rvec xmin,xmax,xcom,xcob;
  
    clear_rvec(xcom);
    copy_rvec(x[0],xmin);
    copy_rvec(x[0],xmax);
    for(i=0; (i<natom); i++)
    {
        rvec_inc(xcom,x[i]);
        for(m=0; (m<DIM); m++)
        {
            if (xmin[m] > x[i][m])
                xmin[m] = x[i][m];
            else if (xmax[m] < x[i][m])
                xmax[m] = x[i][m];
        }
    }
    for(m=0; (m<DIM); m++)
    {
        xcom[m] /= natom;
        box[m][m] = (dbox+xmax[m]-xmin[m]);
        xcob[m] = box[m][m]/2;
    }
}

void write_zeta_q(FILE *fp,gentop_qgen_t qgen,
                  t_atoms *atoms,gmx_poldata_t pd,int iModel)
{
    int    i,j,k,nz,row;
    double zeta,q,qtot;
    gmx_bool   bAtom;
    
    if (NULL == qgen)
        return;
        
    fprintf(fp,"[ zeta_q ]\n");
    fprintf(fp,"; i     type    nq  row       zeta          q\n");
    k = -1;
    for(i=0; (i<atoms->nr); i++)
    {
        bAtom = (atoms->atom[i].ptype == eptAtom);
        if (bAtom)
            k++;
        if (k == -1)
            gmx_fatal(FARGS,"The first atom must be a real atom, not a shell");
        nz = gentop_qgen_get_nzeta(qgen,k);
        if (nz != NOTSET)
        {
            fprintf(fp,"%5d %6s %5d",i+1,get_eemtype_name(iModel),(bAtom) ? nz-1 : 1);
            qtot = 0;
            for(j=(bAtom ? 0 : nz-1); (j<(bAtom ? nz-1 : nz)); j++)
            {
                row   = gentop_qgen_get_row(qgen,k,j);
                q     = gentop_qgen_get_q(qgen,k,j);
                zeta  = gentop_qgen_get_zeta(qgen,k,j);
                if ((row != NOTSET) && (q != NOTSET) && (zeta != NOTSET)) 
                {
                    qtot += q;
                    fprintf(fp,"%5d %10g %10g",row,zeta,q);
                }
            }
            atoms->atom[i].q = qtot;
            fprintf(fp,"\n");
        }
    }
    fprintf(fp,"\n");
}

void write_zeta_q2(gentop_qgen_t qgen,gpp_atomtype_t atype,
                   t_atoms *atoms,gmx_poldata_t pd,int iModel)
{
    FILE   *fp;
    int    i,j,k,nz,row;
    double zeta,q,qtot;
    gmx_bool   bAtom;
    
    if (NULL == qgen)
        return;

    fp = fopen("zeta_q.txt","w");        
    k = -1;
    for(i=0; (i<atoms->nr); i++)
    {
        bAtom = (atoms->atom[i].ptype == eptAtom);
        if (bAtom)
            k++;
        if (k == -1)
            gmx_fatal(FARGS,"The first atom must be a real atom, not a shell");
        nz = gentop_qgen_get_nzeta(qgen,k);
        if (nz != NOTSET)
        {
            fprintf(fp,"%6s  %5s  %5d",get_eemtype_name(iModel),
                    get_atomtype_name(atoms->atom[i].type,atype),
                    (bAtom) ? nz-1 : 1);
            qtot = 0;
            for(j=(bAtom ? 0 : nz-1); (j<(bAtom ? nz-1 : nz)); j++)
            {
                row   = gentop_qgen_get_row(qgen,k,j);
                q     = gentop_qgen_get_q(qgen,k,j);
                zeta  = gentop_qgen_get_zeta(qgen,k,j);
                if ((row != NOTSET) && (q != NOTSET) && (zeta != NOTSET)) 
                {
                    qtot += q;
                    fprintf(fp,"%5d %10g %10g",row,zeta,q);
                }
            }
            atoms->atom[i].q = qtot;
            fprintf(fp,"\n");
        }
    }
    fprintf(fp,"\n");
    fclose(fp);
}

static void unique_atomnames(t_atoms *atoms)
{
    fprintf(stderr,"WARNING: Generating unique atomnames is not implemented yet.\n");
}

int main(int argc, char *argv[])
{
    static const char *desc[] = {
        "gentop generates a topology from molecular coordinates",
        "either from a file, from a database, or from a gaussian log file.",
        "The program assumes all hydrogens are present when defining",
        "the hybridization from the atom name and the number of bonds.",
        "The program can also make an rtp entry, which you can then add",
        "to the rtp database.[PAR]",
        "If the [TT]-oi[tt] option is set an [TT]itp[tt] file will be generated",
        "instead of a [TT]top[tt] file.",
        "When [TT]-param[tt] is set, equilibrium distances and angles",
        "and force constants will be printed in the topology for all",
        "interactions. The equilibrium distances and angles are taken",
        "from the input coordinates, the force constant are set with",
        "command line options.",
        "With the [TT]-db molecule[tt] option a file is extracted from the",
        "database from one of the specified QM calculations (given with [TT]-lot[tt]).",
        "An alternative to the system-wide database [TT]molprops.dat[tt]",
        "can be passed along using the [TT]-mpdb[tt] flag.[PAR]",
        "If the flag [TT]-qgen[tt] is given, charges will be generated using the",
        "specified algorithm. Without the flag the charges from the QM calculation",
        "will be used.",
        "The only supported force field for this tool is Alexandria.[PAR]",
        "oplsaa OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)[PAR]",
        "The corresponding data files can be found in the library directory",
        "with names like ffXXXX.YYY. Check chapter 5 of the manual for more",
        "information about file formats. By default the forcefield selection",
        "is interactive, but you can use the [TT]-ff[tt] option to specify",
        "one of the short names above on the command line instead. In that",
        "case gentop just looks for the corresponding file.[PAR]",
    };
    const char *bugs[] = {
        "No force constants for impropers are generated"
    };
    FILE       *fp;
    t_params   *plist;
    t_excls    *excls;
    t_atoms    *atoms;       /* list with all atoms */
    gpp_atomtype_t atype;
    output_env_t oenv;
    t_mols     mymol;
    gmx_atomprop_t aps;
    gmx_poldata_t  pd;
    char       title[STRLEN];
    rvec       *x=NULL;        /* coordinates? */
    int        *nbonds;
    char       **smnames;
    int        *cgnr;
    int        bts[] = { 1,1,3,1 };
    int        nalloc,nelec,ePBC;
    matrix     box;          /* box length matrix */
    t_pbc      pbc;
    int        natoms;       /* number of atoms in one molecule  */
    int        nres;         /* number of molecules? */
    int        i,j,k,l,m,ndih,dih,cgtp,eQGEN,*symmetric_charges=NULL;
    real       mu;
    gmx_bool   bPDB,bRTP,bTOP;
    t_symtab   symtab;
    int        nqa=0,iModel;
    real       cutoff,qtot,mtot;
    char       *gentop_version = (char *)"v0.99b";
    const char *fn,*xmlf;
    char       rtp[STRLEN],forcefield[STRLEN],ffdir[STRLEN];
    char       ffname[STRLEN],suffix[STRLEN],buf[STRLEN],gentopdat[STRLEN];
    gmx_conect gc = NULL;
    const char *potfn,*rhofn,*hisfn,*difffn,*diffhistfn,*reffn;
    gmx_resp_t gr = NULL;
    gmx_resp_t grref = NULL;
    gentop_vsite_t gvt = NULL;
    gentop_qgen_t  qqgen = NULL;
    char qgen_msg[STRLEN];
    
    t_filenm fnm[] = {
        { efSTX, "-f",    "conf", ffOPTRD },
        { efTOP, "-o",    "out",  ffOPTWR },
        { efITP, "-oi",   "out",  ffOPTWR },
        { efSTO, "-c",    "out",  ffWRITE },
        { efRTP, "-r",    "out",  ffOPTWR },
        { efLOG, "-g03",  "gauss",  ffOPTRD },
        { efNDX, "-n",    "renum", ffOPTWR },
        { efDAT, "-q",    "qout", ffOPTWR },
        { efDAT, "-x",    "mol",  ffOPTWR },
        { efDAT, "-mpdb", "molprops", ffOPTRD },
        { efCUB, "-pot",  "potential", ffOPTWR },
        { efCUB, "-ref",  "refpot", ffOPTRD },
        { efCUB, "-diff", "diffpot", ffOPTWR },
        { efCUB, "-rho",  "density", ffOPTWR },
        { efXVG, "-diffhist", "diffpot", ffOPTWR },
        { efXVG, "-his",  "pot-histo", ffOPTWR },
        { efXVG, "-pc",   "pot-comp", ffOPTWR },
        { efPDB, "-pdbdiff", "pdbdiff", ffOPTWR }
    };
#define NFILE asize(fnm)
    static real scale = 1.1, kb = 4e5,kt = 400,kp = 5;
    static real btol=0.2,qtol=1e-6,zmin=5,zmax=100,delta_z=-1;
    static real hfac=0,qweight=1e-3,bhyper=0.1;
    static real th_toler=170,ph_toler=5,watoms=0,spacing=0.1;
    static real dbox=0.370424,pfac=1,epsr=1;
    static int  qtotref=0,nexcl = 2;
    static int  maxiter=25000,maxcycle=1;
    static int  idihtp=1,pdihtp=3,nmol=1;
    static real rDecrZeta = -1;
    static gmx_bool bRemoveDih=FALSE,bQsym=TRUE,bZatype=TRUE,bFitCube=FALSE;
    static gmx_bool bParam=FALSE,bH14=TRUE,bRound=TRUE,bITP,bAddShells=FALSE;
    static gmx_bool bPairs = TRUE, bPBC = TRUE, bResp = FALSE;
    static gmx_bool bUsePDBcharge = FALSE,bVerbose=FALSE,bAXpRESP=FALSE;
    static gmx_bool bCONECT=FALSE,bRandZeta=FALSE,bFitZeta=TRUE,bEntropy=FALSE;
    static gmx_bool bGenVSites=FALSE,bSkipVSites=TRUE,bUnique=FALSE;
    static char *molnm = "",*dbname = "", *symm_string = "";
    static const char *cqgen[] = { NULL, "None", "Yang", "Bultinck", "Rappe", 
                                   "AXp", "AXs", "AXg", "ESP", "RESP", NULL };
    static const char *dihopt[] = { NULL, "No", "Single", "All", NULL };
    static const char *cgopt[] = { NULL, "Atom", "Group", "Neutral", NULL };
    static const char *lot = "B3LYP/aug-cc-pVTZ",*cat="Other";
    static const char *dzatoms = "";
    static const char *ff = "select";
    t_pargs pa[] = {
        { "-v",      FALSE, etBOOL, {&bVerbose},
          "Generate verbose output in the top file and on terminal." },
        { "-ff",     FALSE, etSTR,  {&ff},
          "Force field, interactive by default. Use -h for information." },
        { "-db",     FALSE, etSTR,  {&dbname},
          "Read a molecule from the database rather than from a file" },
        { "-lot",    FALSE, etSTR,  {&lot},
          "Use this method and level of theory when selecting coordinates and charges" },
        { "-nexcl", FALSE, etINT,  {&nexcl},
          "Number of exclusions" },
        { "-H14",    FALSE, etBOOL, {&bH14}, 
          "Use 3rd neighbour interactions for hydrogen atoms" },
        { "-dih",    FALSE, etSTR,  {dihopt}, 
          "Which proper dihedrals to generate: none, one per rotatable bond, or all possible." },
        { "-remdih", FALSE, etBOOL, {&bRemoveDih}, 
          "Remove dihedrals on the same bond as an improper" },
        { "-pairs",  FALSE, etBOOL, {&bPairs},
          "Output 1-4 interactions (pairs) in topology file" },
        { "-name",   FALSE, etSTR,  {&molnm},
          "Name of your molecule" },
        { "-pbc",    FALSE, etBOOL, {&bPBC},
          "Use periodic boundary conditions." },
        { "-conect", FALSE, etBOOL, {&bCONECT},
          "Use CONECT records in an input pdb file to signify bonds" },
        { "-unique", FALSE, etBOOL, {&bUnique},
          "Make atom names unique" },
        { "-genvsites", FALSE, etBOOL, {&bGenVSites},
          "[HIDDEN]Generate virtual sites for linear groups. Check and double check." },
        { "-skipvsites", FALSE, etBOOL, {&bSkipVSites},
          "[HIDDEN]Skip virtual sites in the input file" },
        { "-pdbq",  FALSE, etBOOL, {&bUsePDBcharge},
          "[HIDDEN]Use the B-factor supplied in a pdb file for the atomic charges" },
        { "-btol",  FALSE, etREAL, {&btol},
          "Relative tolerance for determining whether two atoms are bonded." },
        { "-epsr", FALSE, etREAL, {&epsr},
          "Relative dielectric constant to account for intramolecular polarization. Should be >= 1." },
        { "-spacing", FALSE, etREAL, {&spacing},
          "Spacing of grid points for computing the potential (not used when a reference file is read)." },
        { "-dbox", FALSE, etREAL, {&dbox},
          "[HIDDEN]Extra space around the molecule when generating an ESP output file with the [TT]-pot[tt] option. The strange default value corresponds to 0.7 a.u. that is sometimes used in other programs." },
        { "-axpresp", FALSE, etBOOL, {&bAXpRESP}, 
          "Turn on RESP features for AXp fitting" },
        { "-qweight", FALSE, etREAL, {&qweight},
          "Restraining force constant for the RESP algorithm (AXp only, and with [TT]-axpresp[tt])." },
        { "-bhyper", FALSE, etREAL, {&bhyper},
          "Hyperbolic term for the RESP algorithm (AXp only), and with [TT]-axpresp[tt])." },
        { "-entropy", FALSE, etBOOL, {&bEntropy},
          "Use maximum entropy criterion for optimizing to ESP data rather than good ol' RMS" },
        { "-fitcube", FALSE, etBOOL, {&bFitCube},
          "Fit to the potential in the cube file rather than the log file. This typically gives incorrect results if it converges at all, because points close to the atoms are taken into account on equal footing with points further away." },
        { "-zmin",  FALSE, etREAL, {&zmin},
          "Minimum allowed zeta (1/nm) when fitting models containing gaussian or Slater charges to the ESP" },
        { "-zmax",  FALSE, etREAL, {&zmax},
          "Maximum allowed zeta (1/nm) when fitting models containing gaussian or Slater charges to the ESP" },
        { "-deltaz", FALSE, etREAL, {&delta_z},
          "Maximum allowed deviation from the starting value of zeta. If this option is set then both zmin and zmax will be ignored. A reasonable value would be 10/nm." },
        { "-dzatoms", FALSE, etSTR, {&dzatoms},
          "List of atomtypes for which the fitting is restrained by the -deltaz option." },
        { "-zatype", FALSE, etBOOL, {&bZatype},
          "Use the same zeta for each atom with the same atomtype in a molecule when fitting gaussian or Slater charges to the ESP" },
        { "-decrzeta", FALSE, etREAL, {&rDecrZeta},
          "Generate decreasing zeta with increasing row numbers for atoms that have multiple distributed charges. In this manner the 1S electrons are closer to the nucleus than 2S electrons and so on. If this number is < 0, nothing is done, otherwise a penalty is imposed in fitting if the Z2-Z1 < this number." },
        { "-randzeta", FALSE, etBOOL, {&bRandZeta},
          "Use random zeta values within the zmin zmax interval when optimizing against Gaussian ESP data. If FALSE the initial values from the gentop.dat file will be used." },
        { "-fitzeta", FALSE, etBOOL, {&bFitZeta},
          "Controls whether or not the Gaussian/Slater widths are optimized when fitting to a QM computed ESP" },
        { "-pfac",   FALSE, etREAL, {&pfac},
          "Factor for weighing penalty function for e.g. [TT]-decrzeta[tt] option." },
        { "-watoms", FALSE, etREAL, {&watoms},
          "Weight for the atoms when fitting the charges to the electrostatic potential. The potential on atoms is usually two orders of magnitude larger than on other points (and negative). For point charges or single smeared charges use 0. For point+smeared charges 1 is recommended." },
        { "-param", FALSE, etBOOL, {&bParam},
          "Print parameters in the output" },
        { "-round",  FALSE, etBOOL, {&bRound},
          "Round off measured values" },
        { "-qgen",   FALSE, etENUM, {cqgen},
          "Algorithm used for charge generation" },
        { "-shells", FALSE, etBOOL, {&bAddShells},
          "Add shell particles to the topology" },
        { "-qtol",   FALSE, etREAL, {&qtol},
          "Tolerance for assigning charge generation algorithm" },
        { "-qtot",   FALSE, etINT,  {&qtotref},
          "Net charge on molecule when generating a charge" },
        { "-maxiter",FALSE, etINT, {&maxiter},
          "Max number of iterations for charge generation algorithm" },
        { "-maxcycle", FALSE, etINT, {&maxcycle},
          "Max number of tries for optimizing the charges. The trial with lowest chi2 will be used for generating a topology. Will be turned off if randzeta is No." },
        { "-hfac",    FALSE, etREAL, {&hfac},
          "Fudge factor for AXx algorithms that modulates J00 for hydrogen atoms by multiplying it by (1 + hfac*qH). This hack is originally due to Rappe & Goddard." },
        { "-qsymm",  FALSE, etBOOL, {&bQsym},
          "Symmetrize the charges on symmetric groups, e.g. CH3, NH2." },
        { "-symm",   FALSE, etSTR, {&symm_string},
          "Use the order given here for symmetrizing, e.g. when specifying [TT]-symm '0 1 0'[tt] for a water molecule (H-O-H) the hydrogens will have obtain the same charge. For simple groups, like methyl this is done automatically, but higher symmetry is not detected by the program. The numbers should correspond to atom numbers ,minus 1, and point to either the atom itself or to a previous atom." },
        { "-cgsort", FALSE, etSTR, {cgopt},
          "Option for assembling charge groups: based on Group (default, e.g. CH3 groups are kept together), Atom, or Neutral sections. If the order of atoms is changed an index file is written in order to facilitate changing the order in old files." },
        { "-nmolsort", FALSE, etINT, {&nmol},
          "[HIDDEN]Number of molecules to output to the index file in case of sorting. This is a convenience option to reorder trajectories for use with a new force field." },
        { "-th_toler",FALSE, etREAL, {&th_toler},
          "If bond angles are larger than this value the group will be treated as a linear one and a virtual site will be created to keep the group linear" },
        { "-ph_toler", FALSE, etREAL, {&ph_toler},
          "If dihedral angles are less than this (in absolute value) the atoms will be treated as a planar group with an improper dihedral being added to keep the group planar" },
        { "-kb",    FALSE, etREAL, {&kb},
          "Bonded force constant (kJ/mol/nm^2)" },
        { "-kt",    FALSE, etREAL, {&kt},
          "Angle force constant (kJ/mol/rad^2)" },
        { "-kp",    FALSE, etREAL, {&kp},
          "Dihedral angle force constant (kJ/mol/rad^2)" },
        { "-pdihtp", FALSE, etINT, {&pdihtp},
          "Type to use for proper dihedrals (see GROMACS manual)" },
        { "-idihtp", FALSE, etINT, {&idihtp},
          "Type to use for improper dihedrals (see GROMACS manual)" }
    };
  
    CopyRight(stdout,argv[0]);

    parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
                      asize(desc),desc,asize(bugs),bugs,&oenv);
    
    /* Force field selection, interactive or direct */
    choose_ff(strcmp(ff,"select") == 0 ? NULL : ff,
              forcefield,sizeof(forcefield),
              ffdir,sizeof(ffdir));

    if (strlen(forcefield) > 0) 
    {
        strcpy(ffname,forcefield);
        ffname[0] = toupper(ffname[0]);
    } 
    else 
    {
        gmx_fatal(FARGS,"Empty forcefield string");
    }
    sprintf(gentopdat,"%s/gentop.dat",ffdir);
    printf("\nUsing the %s force field file %s\n\n",
           ffname,gentopdat);
    
    /* Check the options */
    bRTP = opt2bSet("-r",NFILE,fnm);
    bITP = opt2bSet("-oi",NFILE,fnm);
    bTOP = TRUE;
  
    if (!bRandZeta)
        maxcycle = 1;
        
    if (!bRTP && !bTOP)
        gmx_fatal(FARGS,"Specify at least one output file");
  
    if ((btol < 0) || (btol > 1)) 
        gmx_fatal(FARGS,"Bond tolerance should be between 0 and 1 (not %g)",btol);
    if ((qtol < 0) || (qtol > 1)) 
        gmx_fatal(FARGS,"Charge tolerance should be between 0 and 1 (not %g)",qtol);
    /* Check command line options of type enum */
    dih  = get_option(dihopt);
    cgtp = get_option(cgopt);
    if ((iModel = name2eemtype(cqgen[0])) == -1)
        gmx_fatal(FARGS,"Invalid model %s. How could you!\n",cqgen[0]);
    
    /* Set bts for topology output */
    bts[2] = pdihtp;
    bts[3] = idihtp;
    
    /* Read standard atom properties */
    aps = gmx_atomprop_init();
  
    /* Read polarization stuff */
    if ((pd = gmx_poldata_read(gentopdat,aps)) == NULL)
        gmx_fatal(FARGS,"Can not read the force field information. File missing or incorrect.");
    if (bVerbose)
        printf("Reading force field information. There are %d atomtypes.\n",
               gmx_poldata_get_natypes(pd));
  
    /* Init parameter lists */
    snew(plist,F_NRE);
    init_plist(plist);
    snew(atoms,1);
    open_symtab(&symtab);
    reffn  = opt2fn_null("-ref",NFILE,fnm);

    if (strlen(dbname) > 0) 
    {
        gmx_molprop_t *mp;
        int np;
    
        mymol.name = strdup(dbname);
        mymol.nr   = 1;
        if (bVerbose)
            printf("Reading molecule database.\n");
        mp = gmx_molprops_read(opt2fn_null("-mpdb",NFILE,fnm),&np);
        for(i=0; (i<np); i++) 
        {
            if (strcasecmp(dbname,gmx_molprop_get_molname(mp[i])) == 0)
                break;
        }
        if (i == np)
            gmx_fatal(FARGS,"Molecule %s not found in database",dbname);
        if (molprop_2_atoms(mp[i],aps,&symtab,lot,atoms,get_eemtype_name(iModel),
                            &x) == 0)
            gmx_fatal(FARGS,"Could not convert molprop to atoms structure");
        ePBC = epbcNONE;
        clear_mat(box);
        put_in_box(atoms->nr,box,x,dbox);
    }
    else if (opt2bSet("-g03",NFILE,fnm))
    {
        gr = gmx_resp_init(pd,iModel,bAXpRESP,qweight,bhyper,qtotref,
                           zmin,zmax,delta_z,
                           bZatype,watoms,rDecrZeta,bRandZeta,pfac,bFitZeta,
                           bEntropy,dzatoms);
        gmx_resp_read_log(gr,aps,pd,opt2fn("-g03",NFILE,fnm));
        
        gmx_resp_get_atom_info(gr,atoms,&symtab,&x);
        if ((NULL != reffn) && bFitCube)
            gmx_resp_read_cube(gr,reffn,TRUE);
        ePBC = epbcNONE;
        clear_mat(box);
        put_in_box(atoms->nr,box,x,dbox);
        if ((NULL == molnm) || (strlen(molnm) == 0))
        {
            char *ptr = gmx_resp_get_stoichiometry(gr);
            if (NULL == ptr)
                mymol.name = strdup("BOE");
            else
                mymol.name = strdup(ptr);
        }
        else
            mymol.name = strdup(molnm);
        mymol.nr   = 1;
    }
    else 
    {
        if ((NULL == molnm) || (strlen(molnm) == 0))
        {
            mymol.name = strdup("BOE");
            molnm = strdup(mymol.name);
        }
        else
            mymol.name = strdup(molnm);
        mymol.nr   = 1;
        /* Read coordinates */
        get_stx_coordnum(opt2fn("-f",NFILE,fnm),&natoms); 
    
        /* make space for all the atoms */
        init_t_atoms(atoms,natoms,TRUE);
        snew(x,natoms);              
    
        bPDB = (fn2ftp(opt2fn("-f",NFILE,fnm)) == efPDB);
        if (bPDB) 
        {
            if (bCONECT)
                gc = gmx_conect_init();
            read_pdb_conf(opt2fn("-f",NFILE,fnm),title,atoms,x,
                          &ePBC,box,FALSE,gc);
        }
        else
        {
            bCONECT = FALSE;
            read_stx_conf(opt2fn("-f",NFILE,fnm),title,atoms,x,NULL,
                          &ePBC,box);
        }
        if (bSkipVSites) 
        {
            for(i=0; (i<atoms->nr); i++) {
                if (strcmp(*atoms->atomname[i],"ML") == 0) {
                    for(j=i; (j<atoms->nr-1); j++) {
                        atoms->atom[j] = atoms->atom[j+1];
                        atoms->atomname[j] = atoms->atomname[j+1];
                        if (NULL != atoms->atomtype)
                            atoms->atomtype[j] = atoms->atomtype[j+1];
                        if (NULL != atoms->atomtypeB)
                            atoms->atomtypeB[j] = atoms->atomtypeB[j+1];
                    }
                    atoms->nr--;
                }
            }
        }
        get_pdb_atomnumber(atoms,aps);
        
        /* Check input consistency */  
        nelec = 0;
        for(i=0; (i<atoms->nr); i++)
            nelec += atoms->atom[i].atomnumber;
        if (((nelec % 2) == 1) && (((qtotref+32) % 2) != 1))
            gmx_fatal(FARGS,"Number of electrons in %s is %d, but your charge is %d",
                      molnm,nelec,qtotref);
        clean_pdb_names(atoms,&symtab);
        if (bCONECT && (NULL != debug))
            gmx_conect_dump(debug,gc);
        /*for(i=0; (i<atoms->nr); i++)
          t_atoms_set_resinfo(atoms,i,&symtab,mymol.name,1,' ',' ');
        */
    }
    set_pbc(&pbc,ePBC,box);

    if (bVerbose)  
        printf("Generating bonds from distances...\n");
    snew(nbonds,atoms->nr);
    snew(smnames,atoms->nr);
    mk_bonds(pd,atoms,x,gc,plist,nbonds,bH14,(dih == edihAll),bRemoveDih,
             nexcl,&excls,bPBC,box,aps,btol);

    /* Setting the atom types: this depends on the bonding */
    gvt = gentop_vsite_init(egvtALL);
    if ((atype = set_atom_type(stderr,molnm,&symtab,atoms,&(plist[F_BONDS]),
                               nbonds,smnames,pd,aps,x,&pbc,th_toler,ph_toler,
                               gvt)) == NULL) 
        gmx_fatal(FARGS,"Can not find all atomtypes. Better luck next time!");
    
    if (NULL != gr)
        gmx_resp_update_atomtypes(gr,atoms);
    
    if ((NULL != debug)) 
        dump_hybridization(debug,atoms,nbonds);
    sfree(nbonds);
  
    /* Check which algorithm to use for charge generation */
    bQsym = bQsym || (opt2parg_bSet("-symm",asize(pa),pa));
    symmetric_charges = symmetrize_charges(bQsym,atoms,&(plist[F_BONDS]),
                                           pd,aps,symm_string);
    
    if (NULL != gr)
    {    
        gmx_resp_add_atom_info(gr,atoms,pd);
        gmx_resp_add_atom_symmetry(gr,pd,symmetric_charges);
        gmx_resp_summary(stdout,gr,symmetric_charges);
    }
        
    strcpy(qgen_msg,"");
    if (iModel == eqgNone)
    {
        printf("Using charges from gentop.dat\n");
        for(i=0; (i<atoms->nr); i++)
        {
            char *qq = gmx_poldata_get_charge(pd,smnames[i]);
            double q;
            if (NULL != qq)
                sscanf(qq,"%lf",&q);
            else
                q = 0;
            atoms->atom[i].q  = atoms->atom[i].qB = q;
        }
        eQGEN = eQGEN_OK;
    }
    else if (((iModel == eqgESP) || (iModel == eqgRESP)) &&
             opt2parg_bSet("-db",asize(pa),pa))
    {
        printf("Using %s charges\n",get_eemtype_name(iModel));
        eQGEN = eQGEN_OK;
    }
    else 
    {
        qqgen = gentop_qgen_init(pd,atoms,aps,x,iModel,hfac,qtotref,epsr);

        if (NULL == qqgen)
            gmx_fatal(FARGS,"Can not generate charges for %s. Probably due to issues with atomtype detection or support.\n",molnm);
        eQGEN = generate_charges(bVerbose ? stdout : NULL,
                                 qqgen,gr,molnm,pd,atoms,x,qtol,
                                 maxiter,maxcycle,aps,hfac);
        qgen_message(qqgen,sizeof(qgen_msg),qgen_msg,gr);
    }
    if (NULL != gr) 
    {
        /* This has to be done before the grid is f*cked up by 
           writing a cube file */
        gmx_resp_potcomp(gr,opt2fn_null("-pc",NFILE,fnm),
                         opt2fn_null("-pdbdiff",NFILE,fnm),oenv);
        /* Dump potential if needed */
        potfn  = opt2fn_null("-pot",NFILE,fnm);
        rhofn  = opt2fn_null("-rho",NFILE,fnm);
        hisfn  = opt2fn_null("-his",NFILE,fnm);
        difffn = opt2fn_null("-diff",NFILE,fnm);
        diffhistfn = opt2fn_null("-diffhist",NFILE,fnm);
        if ((NULL != potfn) || (NULL != hisfn) || (NULL != rhofn) || 
            ((NULL != difffn) && (NULL != reffn)))
        {
            char buf[256];
            
            sprintf(buf,"Potential generated by %s based on %s charges",
                    gentop_version,
                    get_eemtype_name(iModel));
                    
            if (NULL != difffn)
            {
                grref = gmx_resp_init(pd,iModel,bAXpRESP,qweight,bhyper,
                                      qtotref,zmin,zmax,delta_z,bZatype,watoms,
                                      rDecrZeta,bRandZeta,pfac,bFitZeta,
                                      bEntropy,dzatoms);
                gmx_resp_add_atom_info(grref,atoms,pd);
                gmx_resp_add_atom_symmetry(grref,pd,symmetric_charges);
                gmx_resp_read_cube(grref,opt2fn("-ref",NFILE,fnm),FALSE);
                gmx_resp_copy_grid(gr,grref);
            }
            else 
            {
                gmx_resp_make_grid(gr,spacing,box,x);
            }
            if (NULL != rhofn)
            {
                sprintf(buf,"Electron density generated by %s based on %s charges",
                        gentop_version,get_eemtype_name(iModel));
                gmx_resp_calc_rho(gr);
                gmx_resp_write_rho(gr,rhofn,buf);
            }
            sprintf(buf,"Potential generated by %s based on %s charges",
                    gentop_version,get_eemtype_name(iModel));
            if (NULL != potfn)
            {
                gmx_resp_calc_pot(gr);
                gmx_resp_write_cube(gr,potfn,buf);
            }
            if (NULL != hisfn)
            {
                gmx_resp_write_histo(gr,hisfn,buf,oenv);
            }
            if ((NULL != difffn) || (NULL != diffhistfn))
            {
                sprintf(buf,"Potential difference generated by %s based on %s charges",
                        gentop_version,
                        get_eemtype_name(iModel));
                gmx_resp_write_diff_cube(grref,gr,difffn,diffhistfn,buf,oenv,0);
                gmx_resp_destroy(grref);
            }
        }
        gmx_resp_destroy(gr);
    }
    /* Check whether our charges are OK, quit otherwise */
    if (eQGEN != eQGEN_OK) 
    {
        fprintf(stderr,"Not generating topology and coordinate files\n");
    }
    else 
    {
        int  anr;
        
        anr = atoms->nr;    
        gentop_vsite_generate_special(gvt,bGenVSites,atoms,&x,plist,
                                      &symtab,atype,&excls);
        if (atoms->nr > anr) 
        {
            srenew(smnames,atoms->nr);
            for(i=anr; (i<atoms->nr); i++) {
                smnames[i] = strdup("ML");
            }
        }
    
        if (!bPairs)
            plist[F_LJ14].nr = 0;
        if (dih == edihNo)
            plist[F_PDIHS].nr = 0;
    
        if (bAddShells)
        {
            nalloc = atoms->nr*2+2;
            srenew(x,nalloc);
            srenew(smnames,nalloc);
            srenew(excls,nalloc);
            add_shells(pd,nalloc,atoms,atype,plist,x,&symtab,&excls,smnames);
        }
        mu = calc_dip(atoms,x);
        
        calc_angles_dihs(&plist[F_ANGLES],&plist[F_PDIHS],x,bPBC,box);
    
        if ((cgnr = generate_charge_groups(cgtp,atoms,
                                           &plist[F_BONDS],&plist[F_POLARIZATION],
                                           bUsePDBcharge,
                                           &qtot,&mtot)) == NULL)
            gmx_fatal(FARGS,"Error generating charge groups");
        if (cgtp != ecgAtom)
            sort_on_charge_groups(cgnr,atoms,plist,x,excls,smnames,
                                  opt2fn("-n",NFILE,fnm),nmol);
    
        if (bVerbose) 
        {
            printf("There are %4d proper dihedrals, %4d impropers\n"
                   "          %4d angles, %4d linear angles\n"
                   "          %4d pairs, %4d bonds, %4d atoms\n"
                   "          %4d polarizations\n",
                   plist[F_PDIHS].nr,  plist[F_IDIHS].nr, 
                   plist[F_ANGLES].nr, plist[F_LINEAR_ANGLES].nr,
                   plist[F_LJ14].nr,   plist[F_BONDS].nr,atoms->nr,
                   plist[F_POLARIZATION].nr);
        }
        printf("Total charge is %g, total mass is %g, dipole is %g D\n",
               qtot,mtot,mu);
        reset_q(atoms);
        print_qpol(atoms,smnames,pd);

        snew(atoms->atomtype,atoms->nr);
        for(i=0; (i<atoms->nr); i++)
            atoms->atomtype[i] = put_symtab(&symtab,
                                            get_atomtype_name(atoms->atom[i].type,atype));
    
        set_force_const(plist,kb,kt,kp,bRound,bParam);
        
        if (bUnique)
            unique_atomnames(atoms);
            
        if (bTOP) 
        {    
            /* Write topology file */
            if (bITP)
                fn = ftp2fn(efITP,NFILE,fnm);
            else
                fn = ftp2fn(efTOP,NFILE,fnm);
            fp = ffopen(fn,"w");
            if (!bITP)
                print_top_header(fp,fn,gentop_version,bITP,forcefield,
                                 1.0,qgen_msg);

            if (strlen(dbname) > 0)
            {
                fprintf(fp,"; This topology was generated based on calculations\n");
                fprintf(fp,"; due to Paul J. van Maaren and David van der Spoel\n");
                fprintf(fp,"; (in preparation)\n");
            }
            write_top(fp,NULL,mymol.name,atoms,FALSE,bts,plist,excls,atype,cgnr,nexcl);
            if (bAddShells)
            {
                /* write_zeta_q(fp,qqgen,atoms,pd,iModel);*/
                write_zeta_q2(qqgen,atype,atoms,pd,iModel);
            }
            if (!bITP)
                print_top_mols(fp,mymol.name,forcefield,NULL,0,NULL,1,&mymol);
            fclose(fp);
        }
        if (iModel != eqgNone)
            gentop_qgen_done(atoms,qqgen);
        if (bRTP) 
        {
            /* Write force field component */
            print_rtp((char *)ftp2fn(efRTP,NFILE,fnm),gentop_version,
                      atoms,plist,cgnr,asize(bts),bts);
        }
        /* Write coordinates */ 
        if ((xmlf = opt2fn_null("-x",NFILE,fnm)) != NULL) 
        {
            gmx_molprop_t mp;
            mp = atoms_2_molprop(molnm,atoms,smnames,aps,pd,TRUE,th_toler,ph_toler);
            gmx_molprops_write(xmlf,1,&mp,0);
            gmx_molprop_delete(mp);
        }
    
        close_symtab(&symtab);
    }
    sprintf(title,"%s processed by %s",molnm,ShortProgram());
    write_sto_conf(opt2fn("-c",NFILE,fnm),title,atoms,x,NULL,ePBC,box);
    
  
    thanx(stderr);
  
    return 0;
}
