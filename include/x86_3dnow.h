#ifndef _x86_3dnow_h
#define _x86_3dnow_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if (defined USE_3DNOW && !defined DOUBLE)

void check3dnow();
void vecinvsqrt_3dnow(float in[],float out[],int n);
void vecrecip_3dnow(float in[],float out[],int n);


void inl0100_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],int type[],int ntype,float nbfp[],
		   float Vnb[]);
void inl0110_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],int type[],int ntype,float nbfp[],
		   float Vnb[], int nsatoms[]);
void inl0300_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],int type[],int ntype,float nbfp[],
		   float Vnb[],float tabscale,float VFtab[]);
void inl0310_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],int type[],int ntype,float nbfp[],
		   float Vnb[],float tabscale,float VFtab[], int nsatoms[]);
void inl1000_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
                   float shiftvec[],float fshift[],int gid[],float pos[],
                   float faction[],float charge[],float facel,float Vc[]);
void inl1010_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
                   float shiftvec[],float fshift[],int gid[],float pos[],
                   float faction[],float charge[],float facel, float Vc[],
		   int nsatoms[]);
void inl1020_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
                   float shiftvec[],float fshift[],int gid[],float pos[],
                   float faction[],float charge[],float facel,float Vc[]);
void inl1030_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
                   float shiftvec[],float fshift[],int gid[],float pos[],
                   float faction[],float charge[],float facel,float Vc[]);
void inl1100_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[]);
void inl1110_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[],
		   int nsatoms[]);
void inl1120_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[]);
void inl1130_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[]);
void inl3000_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   float tabscale,float VFtab[]);
void inl3010_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   float tabscale,float VFtab[], int nsatoms[]);
void inl3020_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   float tabscale,float VFtab[]);
void inl3030_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   float tabscale,float VFtab[]);
void inl3100_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[],
		   float tabscale, float VFtab[]);
void inl3110_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[],
		   float tabscale, float VFtab[], int nsatoms[]);
void inl3120_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[],
		   float tabscale, float VFtab[]);
void inl3130_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[],
		   float tabscale, float VFtab[]);
void inl3300_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[],
		   float tabscale,float VFtab[]);
void inl3310_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[],
		   float tabscale,float VFtab[], int nsatoms[]);
void inl3320_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[],
		   float tabscale,float VFtab[]);
void inl3330_3dnow(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
	           float shiftvec[],float fshift[],int gid[],float pos[],
	           float faction[],float charge[],float facel,float Vc[],
		   int type[],int ntype,float nbfp[],float Vnb[],
		   float tabscale,float VFtab[]);

 
#endif
#endif

