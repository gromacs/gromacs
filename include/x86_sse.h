#ifndef _x86_sse_h
#define _x86_sse_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if (defined USE_SSE && !defined DOUBLE)

void checksse();
void vecinvsqrt_sse(float in[],float out[],int n);
void vecrecip_sse(float in[],float out[],int n);

void inl0100_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],int type[],int ntype,float nbfp[],
		 float Vnb[]);
void inl0110_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],int type[],int ntype,float nbfp[],
		 float Vnb[], int nsatoms[]);
void inl0300_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],int type[],int ntype,float nbfp[],
		 float Vnb[],float tabscale,float VFtab[]);
void inl0310_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],int type[],int ntype,float nbfp[],
		 float Vnb[],float tabscale,float VFtab[], int nsatoms[]);
void inl1000_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[]);
void inl1010_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel, float Vc[],
		 int nsatoms[]);
void inl1020_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[]);
void inl1030_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[]);
void inl1100_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 int type[],int ntype,float nbfp[],float Vnb[]);
void inl2000_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 float krf, float crf);
void inl2100_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 float krf, float crf, int type[],int ntype,
		 float nbfp[],float Vnb[]);
void inl1110_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 int type[],int ntype,float nbfp[],float Vnb[],
		 int nsatoms[]);
void inl1120_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 int type[],int ntype,float nbfp[],float Vnb[]);
void inl2020_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 float krf, float crf);
void inl2120_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 float krf, float crf, int type[],int ntype,
		 float nbfp[],float Vnb[]);
void inl1130_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 int type[],int ntype,float nbfp[],float Vnb[]);
void inl2030_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 float krf, float crf);
void inl2130_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 float krf, float crf, int type[],int ntype,
		 float nbfp[],float Vnb[]);
void inl3000_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 float tabscale,float VFtab[]); 
void inl3010_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 float tabscale,float VFtab[], int nsatoms[]);
void inl3020_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 float tabscale,float VFtab[]);
void inl3030_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 float tabscale,float VFtab[]);
void inl3100_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 int type[],int ntype,float nbfp[],float Vnb[],
		 float tabscale, float VFtab[]);
void inl3110_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 int type[],int ntype,float nbfp[],float Vnb[],
		 float tabscale, float VFtab[], int nsatoms[]);
void inl3120_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 int type[],int ntype,float nbfp[],float Vnb[],
		 float tabscale, float VFtab[]);
void inl3130_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 int type[],int ntype,float nbfp[],float Vnb[],
		 float tabscale, float VFtab[]);
void inl3300_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 int type[],int ntype,float nbfp[],float Vnb[],
		 float tabscale,float VFtab[]);
void inl3310_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 int type[],int ntype,float nbfp[],float Vnb[],
		 float tabscale,float VFtab[], int nsatoms[]);
void inl3320_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 int type[],int ntype,float nbfp[],float Vnb[],
		 float tabscale,float VFtab[]);
void inl3330_sse(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		 float shiftvec[],float fshift[],int gid[],float pos[],
		 float faction[],float charge[],float facel,float Vc[],
		 int type[],int ntype,float nbfp[],float Vnb[],
		 float tabscale,float VFtab[]);

#endif
#endif

 
