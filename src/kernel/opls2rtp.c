/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
static char *SRCID_opls2rtp_c = "$Id$";
/**********************************************************************
 * Program: opls2rtp.c                                                *
 * Usage: opls2rtp                                                    *
 *      input files: opls-parameters -> a file with OPLS forcefield   *
 *            parameters in the form:                                 *
 * TYPE ATOMNUMBER ATOMTYPE CHARGE SIGMA EPSILON COMMENT              *
 * No comments are allowed in the file. sigma and epsilon are in A and*
 * kcal/mol. They will be converted to nm and kJ/mol                  *
 *                : opls.rtp -> a file with residue dbase             *
 *            entries in gromacs format, except that the types have   *
 * been substituted by OPLS numbers, bonds, dihedrals, charges, and   *
 * chargegroups have been left out. (`hacked' rtp file)               *
 *                                                                    *
 *                                                                    *
 *      output files: opls.atp with entries OPLSTYPE MASS ; COMMENT   *
 *                    opls-rtp.rtp the residue dbase file             *
 *                   opls-nb.itp the nonbonded interaction file       *
 * the file opls.hdb has to be modified outside this program, and     *
 * opls-bon.itp has to be created using the table opls2gromacs.tab    *
 * and a seperate Perl script.                                        *
 **********************************************************************
 * author: Peter Tieleman. tieleman@chem.rug.nl  September 1994       *
 **********************************************************************/
#define ATP_FILE "opls.atp"
#define ITP_NB_FILE "opls-nb.itp"
#define OPLS_FILE "opls-parameters"

#include <string.h>
#include <stdio.h>
#include "symtab.h"
#include "string2.h"
#include "topio.h"
#include "smalloc.h" /* mem-allocation routines like snew */
#include "resall.h"  /* routines for working with residue dbase files */
#include "copyrite.h" /* copyright message and such */
#include "futil.h"   /* Gromacs file utils, ffopen and such */
#include <math.h>    /* for fabs in the function modify_atoms */
#include <stdlib.h>

typedef struct {    /* struct for atoms in residue dbase. Name, OPLS-code */
  char *aname;      /* the 'official name' of the atom */
  int  op_code;     /* the OPLS-code (int between 1 and 300something) */
} t_opatom;

typedef struct {    /* struct for one residue in residue dbase file */
  char     *resname; /* name of the residue (ex. ARG, ALA, LYS) */
  int      nat;      /* number of atoms in residue */
  t_opatom *op;      /* array of atoms. op[i].op_code and op[i].aname */
} t_opres;           

struct opls_entry {
  int opls_nr; /* atom type. This is an atom number in the OPLS system */
  int atom_number; /* atom number in periodic system of elements */
  char opls_type[10]; /* type in OPLS system */
  float opls_charge; /* charge in opls_system */
  float opls_sigma; /* LJ-parameter sigma (in Angstrom) */
  float opls_epsilon; /* LJ-parameter epsilon */
  char *opls_comment; /* comment in opls-parameters file */
};

struct opls_entry_list{
  struct opls_entry opls_data_entry;
  struct opls_entry_list *next_entry_ptr;
};

struct opls_entry_list *first_entry_ptr = NULL; 
/* pointer to first element of opls-data struct from input file */

/***********************************************************************
 * read_opatoms: reads a hacked residue dbase file.                    *
 * arguments: *fn: name of hacked residue dbase file to be read        *
 *            **opa: array for residues (contains residues after the   *
 * function is done                                                    *
 **********************************************************************/
int read_opatoms(char *fn,t_opres **opa)
{
  FILE    *in;  /* input file */
  char    buf[STRLEN+1]; /* buffer for one line from the input file */
  char    resnm[20]; /* name of residue */
  char    anm[20];   /* 'official' atomname */
  t_opres *opr;      /* pointer to residue */
  int     nopa,nat,type; /* number of residues, number of atoms, OPLS-type */
  int     i,j; /* useless indices */
  
  in=ffopen(fn,"r");     /* Gromacs version of ffopen. Performs some checks */
  fgets2(buf,STRLEN,in); /* Gromacs version of fgets. Performs some checks */
  sscanf(buf,"%d",&nopa); /* number of residues in hacked dbase file */
  snew((*opa),nopa);      /* allocate nopa * sizeof(*opa) using calloc */
  for(i=0; (i<nopa); i++) { /* iterate over all residues */
    opr=&((*opa)[i]);      /* opr points to the i-th element of the res.array*/
    fgets2(buf,STRLEN,in); /* get next line */
    sscanf(buf,"%s",resnm);  /* get name of residue */
    opr->resname=strdup(resnm); /* put name of residue in residue structure */
    fgets2(buf,STRLEN,in);   /* get number of atoms  */
    sscanf(buf,"%d",&nat);
    opr->nat=nat;            /* put number of atoms in residue structure */
    snew(opr->op,nat);       /* allocate memory for all atoms */
    for(j=0; (j<nat); j++) {  
      fgets2(buf,STRLEN,in);  /* for each atom, read name and number */
      sscanf(buf,"%s%d",anm,&type);
      opr->op[j].aname=strdup(anm); /* put this data in residue structure */
      opr->op[j].op_code=type;
    }
  }
  
  fclose(in);       /* close file */
  
  return nopa;      /* return number of residues read */
}


/*****************************************************************
  * pr_atoms: Prints the array read in with read_opatoms         *
 ****************************************************************/
void pr_opatoms(FILE *out,int nopa,t_opres opa[])
{
  int i,j;
  
  fprintf(out,"%d\n",nopa);
  for(i=0; (i<nopa); i++) {
    fprintf(out,"%s\n%5d\n",opa[i].resname,opa[i].nat);
    for(j=0; (j<opa[i].nat); j++)
      fprintf(out,"%10s  %5d\n",opa[i].op[j].aname,opa[i].op[j].op_code);
  }
}


/*****************************************************************
 * modify_atoms: modifies the atoms from the origal rtp file,    *
 * putting the modified residues in a new array with residues,   *
 *  *opls_rtp. It returns this array                             *
******************************************************************/
t_restp * modify_atoms(int nres,t_restp rtp[], int *nwnres)
{
  int find_res_in_opls();

  int nopres; /* number of residues in opls.rtp */
  int i,j;
  int opls_set = 0; /* if 1, don't put residue in the new array with mod. res*/
  t_opres *opa;  /* array with opls residues */
  float nwcharge; /* OPLS-charge */
  float total_charge; /* total charge for chargegroup. Must add up to -1,0,1 */
  int charge_group; /* the new chargegroup for an atom */
  int nwtype; /* OPLS type */
  int status; /* 1: res,atom found, 0: res,atom not found */  
  t_restp *opls_rtp; /* array with the new opls-residues */

  snew(opls_rtp,nres); /* allocate memory for opls-rtp.rtp */
  nopres=read_opatoms("opls.rtp",&opa);
  for(i=0;i<nres;i++){           /* loop over all residues */
    total_charge = 0; 
    charge_group = 0;
    opls_set = 0;
    
    for(j=0;j<rtp[i].natom;j++){ /* loop over all atoms in residue */
      status = find_res_in_opls(rtp[i].resname,
				*(rtp[i].atomname[j]),
				&nwcharge,&nwtype,opa,nopres);
      if(status){ 
	rtp[i].atom[j].type = nwtype-1; /* assign atom OPLS type */
	rtp[i].atom[j].q = nwcharge;/* assign atom OPLS charge */
	rtp[i].cgnr[j] = charge_group; /*charge_group;*/
	total_charge = total_charge + nwcharge;
	if((fabs(total_charge - 0) < 1E-5)
	   || (fabs(total_charge - 1) < 1E-5)
	   || (fabs(total_charge + 1) < 1E-5)){
	  charge_group++; total_charge = 0;
	}
	if(!opls_set) {  /* add residue to the new array with mod. residues */
	  opls_rtp[i] = rtp[i];
	  opls_set = 1; /* only add a residue once, not once for each atom */
	}
      } /* endif status */
    } 
  }
  *nwnres = nres;
  pr_opatoms(stdout,nopres,opa); 
  return opls_rtp;
}

/*****************************************************************
 * find_res_in_opls: Searches for 'atom_name' in residue         *
 * 'residue_name' in the hacked rtp file (*opa). If it finds it, *
 * it changes the type and charge to the opls values and return1 *
******************************************************************/
int find_res_in_opls(char *residue_name, 
		     char *atom_name,
		     float *charge,
		     int *nwtype,
		      t_opres opa[],
		      int numres){

  int i = 0; int j = 0; 
  float find_opls_charge();
/*  fprintf(stdout,"looking for %s and %s\n",residue_name,atom_name); */
  while(i<numres){
    if (!strcmp(opa[i].resname,residue_name)) {
      fprintf(stdout,"found residue %s and atom %s in opls file\n",opa[i].resname,atom_name); 
      for(j=0;j<opa[i].nat;j++){
	if(strcmp(opa[i].op[j].aname,atom_name)==0){
	  *nwtype = opa[i].op[j].op_code;
	  *charge = find_opls_charge(opa[i].op[j].op_code);
	  return 1;
	}
      }
    }
    i++;
  }
  printf("Couldn't find residue %s\n",residue_name);
  return 0;
}

/*****************************************************************
 * find_opls_charge: given a opls_code, get the corresponding    *
 * charge. This information is available in the linked list of   *
 * entries from the opls-parameters file                         *
******************************************************************/
float find_opls_charge(int opls_code){
  struct opls_entry_list *search_entry_ptr;
  struct opls_entry entry_to_search;
 
  search_entry_ptr = first_entry_ptr;
  while(search_entry_ptr != NULL){
    entry_to_search = search_entry_ptr->opls_data_entry;
    if(opls_code == entry_to_search.opls_nr){ 
      return entry_to_search.opls_charge;
    }
    search_entry_ptr = search_entry_ptr->next_entry_ptr;
  }
  fprintf(stderr,"Error: could not find opls-charge. Send hatemail to spoel@chem.rug.nl\n");
return 10; /* nice high charge, will draw attention when looked at */
}

/*****************************************************************
 * read_opls_file: Read the file opls-parameters and build a list*
 * of the entries in that file.                                  *
******************************************************************/
void read_opls_file(){
  FILE *input_file;
  char line[100];
  char *eof_ptr;
  char *comment_ptr;
  int offset = 43; /* comments start at 43 */
  struct opls_entry opls_line;
  struct opls_entry_list *new_entry_ptr = NULL;
  struct opls_entry_list *previous_entry_ptr = NULL;
  int scan_count,i;
  
  input_file = fopen(OPLS_FILE, "r");
  while(1){
    eof_ptr = fgets2(line, sizeof(line), input_file);
    if(eof_ptr == NULL) break;
    scan_count = sscanf(line,"%d %d %s %f %f %f", &opls_line.opls_nr, \
			&opls_line.atom_number, opls_line.opls_type, \
			&opls_line.opls_charge, &opls_line.opls_sigma, \
			&opls_line.opls_epsilon);
    comment_ptr = line + offset;
    opls_line.opls_comment = strdup(comment_ptr);
    for(i=0;i<80;i++)line[i]=' ';
    new_entry_ptr = (struct opls_entry_list *)malloc(sizeof(struct opls_entry_list));
    new_entry_ptr->opls_data_entry = opls_line;
    if (first_entry_ptr == NULL) {
      first_entry_ptr = new_entry_ptr;
      previous_entry_ptr = first_entry_ptr;
    }
    previous_entry_ptr->next_entry_ptr = new_entry_ptr;
    previous_entry_ptr = new_entry_ptr;
  }
}

/*****************************************************************
 * print_opls_file: Print out the list of entries from the opls- *
 * parameter file as read in by read_opls_file                   *
 *****************************************************************/
void print_opls_file(){
  struct opls_entry_list *print_entry_ptr;
  struct opls_entry entry_to_print;
  print_entry_ptr = first_entry_ptr;
  while(print_entry_ptr != NULL){
    entry_to_print = print_entry_ptr->opls_data_entry;
    printf("type: %d number: %d\n", entry_to_print.opls_nr
, entry_to_print.atom_number);
    print_entry_ptr = print_entry_ptr->next_entry_ptr;
  }
}

/***************************************************************** 
 * write_atom_type_file: make the files opls.atp and opls-nb.itp *
 * and write them to disk.                                       *
 *****************************************************************/
void write_atomtype_file(int nr_atom){
  struct opls_entry_list *print_entry_ptr;
  struct opls_entry entry_to_print;
  FILE *output_file; /* the atomic type file */
  FILE *nb_itp_file; /* the nonbonded interaction file */
  int nr=0;
  float sigma; /* sigma in nm, calculated from opls_sigma in angstrom */
  float epsilon; /* epsilon in kJ/mol, from opls_epsilon in kCal/mol */
  
  float atomic_weights[] = {0, 1.00800, 4.0260, 6.941, 9.01218, 10.81, 
    12.011, 14.0067, 15.9994, 18.998403, 20.179, 22.98977, 24.305, 
  26.98154, 28.0855, 30.97376, 32.06, 35.453, 39.948, 39.0983, 
  40.08, 44.9559, 47.88, 50.9415, 51.996, 54.9380, 55.847, 58.9332, 
  58.69, 63.546, 65.38, 69.72, 72.59, 74.9216, 78.96, 79.904, 83.80, 
  85.4678, 87.62, 88.9059, 91.22, 92.9064, 95.94, 98, 101.07, 102.9055, 
  106.42, 107.8682, 112.41, 114.82, 118.69, 121.75, 127.60, 126.9045, 131.29,
  132.9054,137.33};

  output_file = fopen(ATP_FILE,"w");
  nb_itp_file = fopen(ITP_NB_FILE,"w");

  print_entry_ptr = first_entry_ptr;
  
  fprintf(nb_itp_file,"[ atomtypes ]\n");
  fprintf(nb_itp_file,";type\tmass\t\tcharge\t\tsigma(nm)\tepsilon(kj/mol)\tcomment\n");
  fprintf(output_file,"%d ;Type\tMass\tComment\n",nr_atom);  

  while(print_entry_ptr != NULL){
    entry_to_print = print_entry_ptr->opls_data_entry;
    nr = entry_to_print.atom_number;
    if((nr<1) || (nr>56)){
      fprintf(stdout,"Error: atomic number too high: %d. Skipping entry\n",nr);
    } else {
      fprintf(output_file,"OP%d\t %f\t; %s\n", entry_to_print.opls_nr
	      , atomic_weights[entry_to_print.atom_number], entry_to_print.opls_comment);
      sigma = entry_to_print.opls_sigma / 10;
      epsilon = entry_to_print.opls_epsilon * 4.1868;
      fprintf(nb_itp_file,"OP%d\t%f\t%f\t%f\t%f ;%s\n",entry_to_print.opls_nr,
	      atomic_weights[entry_to_print.atom_number],
	 entry_to_print.opls_charge,sigma,epsilon,entry_to_print.opls_comment);
    }
    
    print_entry_ptr = print_entry_ptr->next_entry_ptr;
  }
  fclose(nb_itp_file);
  fclose(output_file);
}

/*****************************************************************
 * main                                                          *
 *****************************************************************/
int main(int argc,char *argv[])
{
  int        nres;
  t_restp    *rtp;
  t_resbond  *rb;
  t_resang   *ra;
  t_idihres  *ires;
  t_atomtype *atype,*atype_opls;
  t_symtab   tab;
  char       *ff="rt37AM";
  int nr_atoms;   /* number of atoms in the atp-file */
  t_restp    *opls_rtp; /* the new array with opls-residues */
  int nwnres; /* number of residues in array with opls-residues */
  FILE       *aa; 

  aa = ffopen("opls-rtp.rtp","w");
  CopyRight(stdout,argv[0]);
  open_symtab(&tab);
  atype=read_atype(ff,&tab);
  atype_opls=read_atype("opls",&tab); /* get opls.atp for print_resall */
  nr_atoms = atype_opls->nr;
  nres=read_resall(ff,&rtp,&rb,&ra,&ires,atype,&tab);
  read_opls_file();
#ifdef DEBUG
  print_opls_file();
#endif 
  write_atomtype_file(nr_atoms); 
  opls_rtp = modify_atoms(nres,rtp, &nwnres);
  printf("Number of residues in new rtp file: %d",nwnres);
  print_resall(aa,nwnres,opls_rtp,rb,ra,ires,atype_opls); 
  fclose(aa);
  thanx(stderr);
  
  printf("hallo! Don't Panic!\n");
  
  return 0;
}
 
