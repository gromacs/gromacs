#ifndef _pp2shift_h
#define _pp2shift_h
	
enum { edPhi=0, edPsi, edChi1, edChi2, edChi3, edChi4, edChi5, edChi6, edMax };

#define NHISTO 360
#define NONCHI 2
#define MAXCHI edMax-NONCHI

typedef struct {
  int H,N,C,O,Cn[MAXCHI+3];
} t_dihatms; /* Cn[0]=N, Cn[1]=Ca, Cn[2]=Cb etc. */

typedef struct {
  char name[12];
  int  resnr;
  int  index;       /* Index for amino acids (histograms) */
  int  j0[edMax];   /* Index in dih array (phi angle is first...) */
  t_dihatms  atm;
  int  b[edMax];
  int  ntr[edMax];
  real S2[edMax];
} t_dlist;

extern void do_pp2shifts(char *fn,int nframes,
			 int nlist,t_dlist dlist[],real **dih);

extern bool has_dihedral(int Dih,t_dlist *dl);

#endif
