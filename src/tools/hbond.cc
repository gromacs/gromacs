/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * GRoups of Organic Molecules in ACtion for Science
 */
static char *SRCID_hbond_cc = "$Id$";

#include "hbond.h"
#include <pbc.h>
#include <physics.h>
#include <vec.h>

static atom_id *w;
static int     nrw;

real pbc_distance2(matrix box,rvec v_1,rvec v_2)
{
  rvec dx;
  pbc_dx(box,v_1,v_2,dx);
  return (iprod(dx,dx));
}


atom_id find_nearest_atom(atom_id atom,atom_id *water,int nrwater)
{
  int i;
  real min_distance2=1000;
  atom_id dummy=0;
  for(i=0;(i<nrwater);i++)
    if ((water[i]!=atom)&&(pbc_distance2(box,x[water[i]],x[atom])<min_distance2)) {
      min_distance2=pbc_distance2(box,x[water[i]],x[atom]);
      dummy=water[i];
    }
  return(dummy);
}

void find_atoms(atom_id atom,atom_id **fw, int& nrfw)
{
  int i;
  nrfw=0;
  for(i=0;(i<nrw);i++)
    if (pbc_distance2(box,x[w[i]],x[atom])<rcut) {
      (*fw) = (atom_id *)realloc(*fw,++nrfw*sizeof(atom_id));
      (*fw)[nrfw-1]=w[i];
    }
  if (nrfw==0) {
    (*fw)=(atom_id *)malloc(sizeof(atom_id));
    (*fw)[0]=find_nearest_atom(atom,w,nrw);
    nrfw=1;
  }
}

void water_list_init(atom_id *wp,int nrwp)
{
  w = wp;
  nrw = nrwp;
}

Hbond::Hbond(atom_id dd,atom_id aa,atom_id hh)
{
  d = dd;
  a = aa;
  h = hh;
  freq=0;
  type=INTER;
}

Hbond::Hbond(Hbond *that)
{
  d = that->d;
  a = that->a;
  h = that->h;
  freq=0;
  type=INTER;
}

Hbond::~Hbond()
{
}

real Hbond::distance2()
{
  rvec dx;
  pbc_dx(box,x[a],x[d],dx);
  return (iprod(dx,dx));
}

real Hbond::angle()
{
  rvec v_1,v_2;
  pbc_dx(box,x[d],x[h],v_1);
  pbc_dx(box,x[d],x[a],v_2);
  return acos(cos_angle(v_1,v_2));
}

#ifdef ENERGY 
real Hbond::energy()
{
  /* calculate coulomb energy */
  real qener;
  rvec dx;
    
  /* donor and hydrogen */
  pbc_dx(box,x[d],x[h],dx);
  qener = top->atoms.atom[d].q * top->atoms.atom[h].q / norm(dx);

  /* donor and acceptor */
  pbc_dx(box,x[d],x[a],dx);
  qener += top->atoms.atom[d].q * top->atoms.atom[a].q / norm(dx);

  /* hydrogen and acceptor */
  pbc_dx(box,x[h],x[a],dx);
  qener += top->atoms.atom[h].q * top->atoms.atom[a].q / norm(dx);

  return qener;
}
#endif

bool Hbond::exist()
{
  if ((distance2()<rcut2)&&(angle()<alfcut))
    return TRUE;
  else 
    return FALSE;
} 

void Hbond::print(FILE *output)
{
  fprintf(output,"%5d%-5.5s%5.5s%5d--",
	  top->atoms.atom[d].resnr+1,
	  *top->atoms.resname[top->atoms.atom[d].resnr],
	  *top->atoms.atomname[d],
	  d+1);
  fprintf(output,"%5d%-5.5s%5.5s%5d--",
	  top->atoms.atom[h].resnr+1,
	  *top->atoms.resname[top->atoms.atom[h].resnr],
	  *top->atoms.atomname[h],
	  h+1);
  fprintf(output,"%5d%-5.5s%5.5s%5d",
	  top->atoms.atom[a].resnr+1,
	  *top->atoms.resname[top->atoms.atom[a].resnr],
	  *top->atoms.atomname[a],
	  a+1);
  fprintf(output,"%8.3f %8.3f %8d",first < 0.0 ? 0.0 : first,last,freq);
}

void Hbond::analyse_init()
{
  freq=0;
  first=-1;
  last=0;

  /* malloc distance array for selected hbonds */
  if ((mode==SELECTED)||(mode==INSERT)) 
    dist=(real *)malloc(nr_frames * sizeof(real));

  if (mode==INSERT) {
    d_dist=(real *)malloc(nr_frames * sizeof(real));
    a_dist=(real *)malloc(nr_frames * sizeof(real));
    h_dist=(real *)malloc(nr_frames * sizeof(real));
    waid  =(atom_id *)malloc(nr_frames * sizeof(atom_id));
  }
  
}

void Hbond::analyse()
{
  if ( exist() ) {
    freq++;
    if ( first < 0 )
      first=this_time;
    last=this_time;
  }

  /* update selected distances */
  if ((mode==SELECTED)||(mode==INSERT)) {
    dist[this_frame]=sqrt(distance2());
  }

  if (mode==INSERT) {
    atom_id *water_index=NULL;
    int water_size=0;
    find_atoms(a,&water_index,water_size);
    waid[this_frame]=find_nearest_atom(d,water_index,water_size);
    d_dist[this_frame]=sqrt(pbc_distance2(box,x[waid[this_frame]],x[d]));
    a_dist[this_frame]=sqrt(pbc_distance2(box,x[waid[this_frame]],x[a]));
    h_dist[this_frame]=sqrt(pbc_distance2(box,x[waid[this_frame]],x[h]));
    dist[this_frame]=sqrt(pbc_distance2(box,x[d],x[a]));
  }
}
 
bool Hbond::insert()
{
  if (mode!=INSERT) {
    fprintf(stderr,"error in program\n Please contact gromacs@chem.rug.nl\n");
    exit(0);
  }
  
  if ((a_dist[this_frame]<rcut) && (d_dist[this_frame]<rcut))
    return TRUE;
  else
    return FALSE;
  
}

void Hbond::settype(atom_id left,atom_id right)
{
  int ra,rd;
  if ((d>=left) && (d<right) && (a>=left) && (a<right)) {
    type=INTRA;
    /* Check for helix stuff */
    rd=top->atoms.atom[d].resnr;
    ra=top->atoms.atom[a].resnr;
    switch (abs(rd-ra)) {
    case 3:
      hhb=NN3;
      break;
    case 4:
      hhb=NN4;
      break;
    case 5:
      hhb=NN5;
      break;
    default:
      hhb=NN0;
      break;
    }
  } 
}

void Hbond::ndx(t_type t,FILE *fp)
{
  if ((type==t)||(t==ALL))
    fprintf(fp," %5d %5d %5d\n",d,h,a);  
}

int Hbond::compare(Hbond *that)
{
  if (this->d<that->d)
    return(-1);
  if (this->d>that->d)
    return(1);
  if (this->a<that->a)
    return(-1);
  if (this->a>that->a)
    return(1);
  if (this->h<that->h)
    return(-1);
  if (this->h>that->h)
    return(1);
  return(0);

}

char *Hbond::hb_name(t_atoms *atoms)
{
  static char buf[256];
  int    nd,na;
  
  nd=atoms->atom[d].resnr;
  na=atoms->atom[a].resnr;
  sprintf(buf,"%s%d:%s%s-%s%d:%s",*(atoms->resname[nd]),nd+1,
	  *(atoms->atomname[d]),*(atoms->atomname[h]),
	  *(atoms->resname[na]),na+1,*(atoms->atomname[a]));
  
  return buf;
}


/* water insertion functions */
void Hbond::insert_print(int t,FILE *fp)
{
  fprintf(fp,"%8.3f %8.3f %8.3f %5d",dist[t],d_dist[t],a_dist[t],waid[t]+1);
}

void Hbond::selected_print(int t,FILE *fp)
{
  fprintf(fp," %8.3f ",dist[t]);
}








