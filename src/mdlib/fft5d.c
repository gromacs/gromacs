/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fft5d.h"
#include <float.h>
#include <math.h>
#include <assert.h>

static double fft5d_fmax(double a, double b){
	return (a>b)?a:b;
}
static double fft5d_fmin(double a, double b){
	return (a<b)?a:b;
}

/* largest factor smaller than sqrt */
static int lfactor(int z) {  
	int i;
	for (i=sqrt(z);;i--)
		if (z%i==0) return i;
}

/* largest factor */
static int l2factor(int z) {  
	int i;
	if (z==1) return 1;
	for (i=z/2;;i--)
		if (z%i==0) return i;
}

/* largest prime factor: WARNING: slow recursion, only use for small numbers */
static int lpfactor(int z) {
	int f = l2factor(z);
	if (f==1) return z;
	return fft5d_fmax(lpfactor(f),lpfactor(z/f));
}


/* NxMxK the size of the data
 * comm communicator to use for fft5d
 * P0 number of processor in 1st axes (can be null for automatic)
 * lin is allocated by fft5d because size of array is only known after planning phase */

fft5d_plan fft5d_plan_3d(int NG, int MG, int KG, MPI_Comm comm[2], fft5d_flags flags, fft5d_type** rlin, fft5d_type** rlout, FILE* debug) {
    int P[2],bMaster,prank[2];
	/* comm, prank and P are in the order of the decomposition (plan->cart is in the order of transposes) */
    if (gmx_parallel_env_initialized() && comm[0] != NULL)
    {
	    MPI_Comm_size(comm[0],&P[0]);
	    MPI_Comm_rank(comm[0],&prank[0]);
    }
    else
    {
        P[0] = 0;
        prank[0] = 0;
    }
    if (gmx_parallel_env_initialized() && comm[1] != NULL)
    {
	    MPI_Comm_size(comm[1],&P[1]);
	    MPI_Comm_rank(comm[1],&prank[1]);
    }
    else
    {
        P[1] = 0;
        prank[1] = 0;
    }
   
	bMaster=(prank[0]==0&&prank[1]==0);
	
	if (bMaster && debug) fprintf(debug,"FFT5D: Using %dx%d processor grid\n",P[0],P[1]);
	
	
	if (bMaster) {
   	        if (debug) 
		    fprintf(debug,"FFT5D: N: %d, M: %d, K: %d, P: %dx%d, real2complex: %d, backward: %d, order yz: %d, debug %d\n",
			    NG,MG,KG,P[0],P[1],(flags&FFT5D_REALCOMPLEX)>0,(flags&FFT5D_BACKWARD)>0,(flags&FFT5D_ORDER_YZ)>0,(flags&FFT5D_DEBUG)>0);
		if (fft5d_fmax(fft5d_fmax(lpfactor(NG),lpfactor(MG)),lpfactor(KG))>7) {
			printf("WARNING: FFT very slow with prime factors larger 7\n");
			printf("Change FFT size or in case you cannot change it look at\n");
			printf("http://www.fftw.org/fftw3_doc/Generating-your-own-code.html\n");
		}
	}
	
	if (NG==0 || MG==0 || KG==0) {
		if (bMaster) printf("FFT5D: FATAL: Datasize cannot be zero in any dimension\n");
		return 0;
	}

	int rNG=NG,rMG=MG,rKG=KG;
	
	if (flags&FFT5D_REALCOMPLEX) {
		if (!(flags&FFT5D_BACKWARD)) NG = NG/2+1;
		else {
			if (!(flags&FFT5D_ORDER_YZ)) MG=MG/2+1;
			else KG=KG/2+1;
		}
	}
	
	/* setting it to the length of the current processor would allow also other 
       distribution than last processor less. And would not require max() in 
       loops */
	int N0=ceil((double)NG/P[0]),N1=ceil((double)NG/P[1]); 
	int M0=ceil((double)MG/P[0]),M1=ceil((double)MG/P[1]); 
	int K0=ceil((double)KG/P[0]),K1=ceil((double)KG/P[1]); 

	/* for step 1-3 the local N,M,K sizes of the transposed system
	   C: contiguous dimension, and nP: number of processor in subcommunicator 
       for that step */
	int N[3],M[3],K[3],C[3],rC[3],nP[2];
	
	
	M[0] = M0;
	K[0] = K1;
	C[0] = NG;
	rC[0] = rNG;
	if (!(flags&FFT5D_ORDER_YZ)) {
		N[0] = N1;
		nP[0] = P[1];
		C[1] = KG;
		rC[1] =rKG;
		N[1] = K0;
		M[1] = M0;
		K[1] = N1;
		nP[1] = P[0];
		C[2] = MG;
		rC[2] = rMG;
		M[2] = K0;
		K[2] = N1;
	} else {
		N[0] = N0;
		nP[0] = P[0];
		C[1] = MG;
		rC[1] =rMG;
		N[1] = M1;
		M[1] = N0;
		K[1] = K1;
		nP[1] = P[1];
		C[2] = KG;
		rC[2] = rKG;
		M[2] = N0;
		K[2] = M1;
	}
	
	/*
      Difference between x-y-z regarding 2d decomposition is whether they are 
      distributed along axis 1, 2 or both 
    */
	
	/* int lsize = fmax(N[0]*M[0]*K[0]*nP[0],N[1]*M[1]*K[1]*nP[1]); */
	int lsize = fft5d_fmax(N[0]*M[0]*K[0]*nP[0],fft5d_fmax(N[1]*M[1]*K[1]*nP[1],C[2]*M[2]*K[2])); 
	/* int lsize = fmax(C[0]*M[0]*K[0],fmax(C[1]*M[1]*K[1],C[2]*M[2]*K[2])); */
	fft5d_type* lin,*lout;
	if (!(flags&FFT5D_NOMALLOC)) { 
        /* local in	*/
		lin = (fft5d_type*)FFTW(malloc)(sizeof(fft5d_type) * lsize); 
        /* local output */
		lout = (fft5d_type*)FFTW(malloc)(sizeof(fft5d_type) * lsize); 
	} else {
		lin = *rlin;
		lout = *rlout;
	}
	
	int fftwflags=FFTW_DESTROY_INPUT;
	if (!(flags&FFT5D_NOMEASURE)) fftwflags|=FFTW_MEASURE;
	fft5d_type* output=lout;
	fft5d_plan plan = (fft5d_plan)malloc(sizeof(struct fft5d_plan_t));
	int s;
	for (s=0;s<3;s++) {
		if ((flags&FFT5D_INPLACE) && s==2) {
			output=lin;
			fftwflags&=~FFTW_DESTROY_INPUT;
		}
		if ((flags&FFT5D_REALCOMPLEX) && !(flags&FFT5D_BACKWARD) && s==0) {
			plan->p1d[s] = FFTW(plan_many_dft_r2c)(1, &rC[s], M[s]*K[s],   
                    (fft5d_rtype*)lin, &rC[s], 1,   C[s]*2, /* why *2 */
					(FFTW(complex)*)output, &C[s], 1,   C[s], fftwflags);
		} else if ((flags&FFT5D_REALCOMPLEX) && (flags&FFT5D_BACKWARD) && s==2) {
			plan->p1d[s] = FFTW(plan_many_dft_c2r)(1, &rC[s], M[s]*K[s],   
					(FFTW(complex)*)lin, &C[s], 1,   C[s], 
					(fft5d_rtype*)output, &rC[s], 1,   C[s]*2, fftwflags);
		} else {
			plan->p1d[s] = FFTW(plan_many_dft)(1, &C[s], M[s]*K[s],   
					(FFTW(complex)*)lin, &C[s], 1,   C[s], 
					(FFTW(complex)*)output, &C[s], 1,   C[s], (flags&FFT5D_BACKWARD)?1:-1, fftwflags);
		}
	}

	if ((flags&FFT5D_ORDER_YZ)) { /*plan->cart is in the order of transposes */
	    plan->cart[0]=comm[0]; plan->cart[1]=comm[1];
	} else {
	    plan->cart[1]=comm[0]; plan->cart[0]=comm[1];
	}
#ifdef FFT5D_MPI_TRANSPOSE
	for (s=0;s<2;s++) {
		plan->mpip[s] = FFTW(mpi_plan_many_transpose)(nP[s], nP[s], N[s]*K[s]*M[s]*2, 1, 1, (fft5d_rtype*)lin, (fft5d_rtype*)lout, plan->comm[s], FFTW_PATIENT);
	}
#endif 

	
	plan->lin=lin;
	plan->lout=lout;
	
	plan->NG=NG;plan->MG=MG;plan->KG=KG;
	
	for (s=0;s<3;s++) {
		plan->N[s]=N[s];plan->M[s]=M[s];plan->K[s]=K[s];plan->C[s]=C[s];plan->rC[s]=rC[s];
	}
	for (s=0;s<2;s++) {
		plan->P[s]=nP[s];plan->coor[s]=prank[s];
	}
	
/*	plan->fftorder=fftorder;
	plan->direction=direction;	
	plan->realcomplex=realcomplex;
*/
	plan->flags=flags;
	*rlin=lin;
	*rlout=lout;
	return plan;
}


enum order {
	XYZ,
	XZY,
	YXZ,
	YZX,
	ZXY,
	ZYX
};



/*here x,y,z and N,M,K is in rotated coordinate system!!
  x (and N) is mayor (consecutive) dimension, y (M) middle and z (K) major
  N,M,K is size of local data!
  NG, MG, KG is size of global data*/
static void splitaxes(fft5d_type* lin,const fft5d_type* lout,int N,int M,int K,int P,int NG) {
	int x,y,z,i;
	for (i=0;i<P;i++) { /*index cube along long axis*/
		for (z=0;z<K;z++) { /*3. z l*/
			for (y=0;y<M;y++) { /*2. y k*/
				for (x=0;x<fft5d_fmin(N,NG-N*i);x++) { /*1. x j*/
					lin[x+y*N+z*M*N+i*M*N*K]=lout[(i*N+x)+y*NG+z*NG*M];
				}
			}
		}
	}
}

/*make axis contiguous again (after AllToAll) and also do local transpose*/
/*transpose mayor and major dimension
  variables see above
  the major, middle, minor order is only correct for x,y,z (N,M,K) for the input
  N,M,K local dimensions
  KG global size*/
static void joinAxesTrans13(fft5d_type* lin,const fft5d_type* lout,int N,int M,int K,int P,int KG) {
	int i,x,y,z;
	for (i=0;i<P;i++) { /*index cube along long axis*/
		for (x=0;x<N;x++) { /*1.j*/
			for (z=0;z<fft5d_fmin(K,KG-K*i);z++) { /*3.l*/
				for (y=0;y<M;y++) { /*2.k*/
					lin[(i*K+z)+y*KG+x*KG*M]=lout[x+y*N+z*M*N+i*M*N*K];
				}
			}
		}
	}
}

/*make axis contiguous again (after AllToAll) and also do local transpose
  tranpose mayor and middle dimension
  variables see above
  the minor, middle, major order is only correct for x,y,z (N,M,K) for the input
  N,M,K local size
  MG, global size*/
static void joinAxesTrans12(fft5d_type* lin,const fft5d_type* lout,int N,int M,int K,int P,int MG) {
	int i,z,y,x;
	for (i=0;i<P;i++) { /*index cube along long axis*/
		for (z=0;z<K;z++) { 
			for (x=0;x<N;x++) { 
				for (y=0;y<fft5d_fmin(M,MG-M*i);y++) { 
					lin[(i*M+y)+x*MG+z*MG*N]=lout[x+y*N+z*M*N+i*M*N*K];
				}
			}
		}
	}
}


static void rotate(int x[]) {
	int t=x[0];
/*	x[0]=x[2];
	x[2]=x[1];
    x[1]=t;*/
	x[0]=x[1];
	x[1]=x[2];
	x[2]=t;
}

/*compute the offset to compare or print transposed local data in original input coordinates
  xo matrix offset, xl dimension length, xc decomposition offset 
  s: step in computation = number of transposes*/
static void compute_offsets(fft5d_plan plan, int xo[], int xl[], int xc[], int NG[], int s) {
/*	int direction = plan->direction;
	int fftorder = plan->fftorder;*/
	
	int o;
	
	NG[0]=plan->NG;NG[1]=plan->MG;NG[2]=plan->KG;

	if (!(plan->flags&FFT5D_ORDER_YZ)) {
		switch (s) {
		case 0: o=XYZ; break;
		case 1: o=ZYX; break;
		case 2: o=YZX; break;
		default: assert(0);
		}
	} else {
		switch (s) {
		case 0: o=XYZ; break;
		case 1: o=YXZ; break;
		case 2: o=ZXY; break;
		default: assert(0);
		}
	}
 
	int pos[3],i;
	
	switch (o) {
		case XYZ:pos[0]=1;pos[1]=2;pos[2]=3;break;
		case XZY:pos[0]=1;pos[1]=3;pos[2]=2;break;
		case YXZ:pos[0]=2;pos[1]=1;pos[2]=3;break;
		case YZX:pos[0]=3;pos[1]=1;pos[2]=2;break;
		case ZXY:pos[0]=2;pos[1]=3;pos[2]=1;break;
		case ZYX:pos[0]=3;pos[1]=2;pos[2]=1;break;
	}
	/*if (debug) printf("pos: %d %d %d\n",pos[0],pos[1],pos[2]);*/
	int *M=plan->M, *K=plan->K, *C=plan->C, *rC=plan->rC;
	int *coor=plan->coor;
	
	/*xo, xl give offset and length in local transposed coordinate system
      for 0(/1/2): x(/y/z) in original coordinate system*/
	for (i=0;i<3;i++) {
		switch (pos[i]) {
		case 1: xo[i]=1;         xc[i]=0;            xl[i]=C[s];break;
		case 2: xo[i]=C[s];      xc[i]=coor[0]*M[s]; xl[i]=fft5d_fmin(M[s],NG[i]-coor[0]*M[s]);break;
		case 3: xo[i]=C[s]*M[s]; xc[i]=coor[1]*K[s]; xl[i]=fft5d_fmin(K[s],NG[i]-coor[1]*K[s]);break;
		}
	}
    /*input order is different for test programm to match FFTW order 
      (important for complex to real)*/
	if (plan->flags&FFT5D_BACKWARD) {
		rotate(xo);
		rotate(xl);
		rotate(xc);
		rotate(NG);
		if (plan->flags&FFT5D_ORDER_YZ) {
			rotate(xo);
			rotate(xl);
			rotate(xc);
			rotate(NG);			
		}
	}
	if (plan->flags&FFT5D_REALCOMPLEX && ((!(plan->flags&FFT5D_BACKWARD) && s==0) || (plan->flags&FFT5D_BACKWARD && s==2))) {
		xl[0] = rC[s];
	}
}

/*N, M, K not used anymore*/
static void print_localdata(const fft5d_type* lin, const char* txt, int N,int M,int K, int s, fft5d_plan plan) {
	int x,y,z,l;
	int *coor = plan->coor;
	int xo[3],xl[3],xc[3],NG[3];		
	compute_offsets(plan,xo,xl,xc,NG,s);
	int ll=(plan->flags&FFT5D_REALCOMPLEX)?1:2;
	printf(txt,coor[0],coor[1],s);
	/*printf("xo: %d %d %d, xl: %d %d %d\n",xo[0],xo[1],xo[2],xl[0],xl[1],xl[2]);*/
	for (z=0;z<xl[2];z++) {
		for(y=0;y<xl[1];y++) {
			printf("%d %d: ",coor[0],coor[1]);
			for (x=0;x<xl[0];x++) {
				for (l=0;l<ll;l++) {
					printf("%f ",((fft5d_rtype*)lin)[(z*xo[2]+y*xo[1])*2+(x*xo[0])*ll+l]);
				}
				printf(",");
			}
			printf("\n");
		}
	}
}

void fft5d_execute(fft5d_plan plan,fft5d_time times) {
	fft5d_type *lin = plan->lin;
	fft5d_type *lout = plan->lout;

	FFTW(plan) *p1d=plan->p1d;
#ifdef FFT5D_MPI_TRANSPOSE
	FFTW(plan) *mpip=plan->mpip;
#else
	MPI_Comm *cart=plan->cart;
#endif
	double time_fft=0,time_local=0,time_mpi[2]={0},time;	
	int *N=plan->N,*M=plan->M,*K=plan->K,*C=plan->C,*P=plan->P;
	
	

	/*lin: x,y,z*/
	int s=0;
	if (plan->flags&FFT5D_DEBUG) print_localdata(lin, "%d %d: copy in lin\n", C[0], M[0], K[0], s, plan);
	for (s=0;s<2;s++) {
		time=MPI_Wtime();
		FFTW(execute)(p1d[s]); /*in:lin out:lout*/
		time_fft+=MPI_Wtime()-time;
	
		if (plan->flags&FFT5D_DEBUG) print_localdata(lout, "%d %d: FFT %d\n", C[s], M[s], K[s], s, plan);
		
		time=MPI_Wtime(); 
		/*prepare for AllToAll
          1. (most outer) axes (x) is split into P[s] parts of size N[s] 
             for sending*/
		splitaxes(lin,lout,N[s],M[s],K[s],P[s],C[s]);
		time_local+=MPI_Wtime()-time;
		
		/*send, recv*/
		time=MPI_Wtime();
#ifdef GMX_MPI
        if (gmx_parallel_env_initialized() && cart[s] != NULL)
        {
#ifdef FFT5D_MPI_TRANSPOSE
			FFTW(execute)(mpip[s]);
#else
            MPI_Alltoall(lin,N[s]*M[s]*K[s]*sizeof(fft5d_type)/sizeof(fft5d_rtype),FFT5D_MPI_RTYPE,lout,N[s]*M[s]*K[s]*sizeof(fft5d_type)/sizeof(fft5d_rtype),FFT5D_MPI_RTYPE,cart[s]);
#endif
        }
        else
#endif
        {
            memcpy(lin,lout,N[s]*M[s]*K[s]*sizeof(fft5d_type));
        }
		time_mpi[s]=MPI_Wtime()-time;
	
		time=MPI_Wtime();
		/*bring back in matrix form (could be avoided by storing blocks as eleftheriou, really?)
          thus make  new 1. axes contiguos
          also local transpose 1and 2/3 */
		if ((s==0 && !(plan->flags&FFT5D_ORDER_YZ)) || (s==1 && (plan->flags&FFT5D_ORDER_YZ))) 
			joinAxesTrans13(lin,lout,N[s],M[s],K[s],P[s],C[s+1]);
		else 
			joinAxesTrans12(lin,lout,N[s],M[s],K[s],P[s],C[s+1]);	
		time_local+=MPI_Wtime()-time;
	
		if (plan->flags&FFT5D_DEBUG) print_localdata(lin, "%d %d: tranposed %d\n", C[s+1], M[s+1], K[s+1], s+1, plan);
				
		/*if (debug) print_localdata(lin, "%d %d: transposed x-z\n", N1, M0, K, ZYX, coor);*/
	}	
	
	time=MPI_Wtime();
	FFTW(execute)(p1d[2]);
	time_fft+=MPI_Wtime()-time;
	if (plan->flags&FFT5D_DEBUG) print_localdata(lout, "%d %d: FFT %d\n", C[s], M[s], K[s], s, plan);
	/*if (debug) print_localdata(lout, "%d %d: FFT in y\n", N1, M, K0, YZX, coor);*/
	
	if (times!=0) {
		times->fft+=time_fft;
		times->local+=time_local;
		times->mpi2+=time_mpi[1];
		times->mpi1+=time_mpi[0];
	}
}

void fft5d_destroy(fft5d_plan plan) {
	int s;
	for (s=0;s<3;s++)
		FFTW(destroy_plan)(plan->p1d[s]);
	
#ifdef FFT5D_MPI_TRANSPOSE
	for (s=0;s<2;s++)	
		FFTW(destroy_plan)(plan->mpip[s]);
#endif
	if (!(plan->flags&FFT5D_NOMALLOC)) { 
		FFTW(free)(plan->lin);
		FFTW(free)(plan->lout);
	}
	free(plan);
	
}

/*TODO better than direct access of plan? enough data?
  here 0,1 reference divided by which processor grid dimension (not FFT step!)*/
void fft5d_local_size(fft5d_plan plan,int* N1,int* M0,int* K0,int* K1,int** coor) {
	*N1=plan->N[0];
	*M0=plan->M[0];
	*K1=plan->K[0];
	*K0=plan->N[1];
	
	*coor=plan->coor;
}


/*same as fft5d_plan_3d but with cartesian coordinator and automatic splitting 
  of processor dimensions*/
fft5d_plan fft5d_plan_3d_cart(int NG, int MG, int KG, MPI_Comm comm, int P0, fft5d_flags flags, fft5d_type** rlin, fft5d_type** rlout, FILE* debug) {
	int size,prank;
	MPI_Comm_size(comm,&size);
	MPI_Comm_rank(comm,&prank);
	if (P0==0) P0 = lfactor(size);
	if (size%P0!=0) {
		if (prank==0) printf("FFT5D: WARNING: Number of processors %d not evenly dividable by %d\n",size,P0);
		P0 = lfactor(size);
	}
		
	int P[] = {P0,size/P0}; /*number of processors in the two dimensions*/
	
	/*Difference between x-y-z regarding 2d decomposition is whether they are 
      distributed along axis 1, 2 or both*/
	
	int coor[2];
	
	int wrap[]={0,0};
	MPI_Comm gcart;
	MPI_Cart_create(comm,2,P,wrap,1,&gcart); /*parameter 4: value 1: reorder*/
	MPI_Cart_get(gcart,2,P,wrap,coor); 
	int rdim1[] = {0,1}, rdim2[] = {1,0};
	MPI_Comm cart[2];
	MPI_Cart_sub(gcart, rdim1 , &cart[0]);
	MPI_Cart_sub(gcart, rdim2 , &cart[1]);

	return fft5d_plan_3d(NG, MG, KG, cart, flags, rlin, rlout, debug); 
}



/*prints in original coordinate system of data (as the input to FFT)*/
void fft5d_compare_data(const fft5d_type* lin, const fft5d_type* in, fft5d_plan plan, int bothLocal, int normalize) {
	int xo[3],xl[3],xc[3],NG[3];
	int x,y,z,l;
	int *coor = plan->coor;
	int ll=2; /*compare ll values per element (has to be 2 for complex)*/
	if (plan->flags&FFT5D_REALCOMPLEX && plan->flags&FFT5D_BACKWARD) 
		ll=1;
		
	compute_offsets(plan,xo,xl,xc,NG,2);
	if (plan->flags&FFT5D_DEBUG) printf("Compare2\n");
	for (z=0;z<xl[2];z++) {
		for(y=0;y<xl[1];y++) {
			if (plan->flags&FFT5D_DEBUG) printf("%d %d: ",coor[0],coor[1]);
			for (x=0;x<xl[0];x++) {
			    for (l=0;l<ll;l++) { /*loop over real/complex parts*/
					fft5d_rtype a,b;
					a=((fft5d_rtype*)lin)[(z*xo[2]+y*xo[1])*2+x*xo[0]*ll+l];
					if (normalize) a/=plan->rC[0]*plan->rC[1]*plan->rC[2];
					if (!bothLocal) 
						b=((fft5d_rtype*)in)[((z+xc[2])*NG[0]*NG[1]+(y+xc[1])*NG[0])*2+(x+xc[0])*ll+l];
					else 
						b=((fft5d_rtype*)in)[(z*xo[2]+y*xo[1])*2+x*xo[0]*ll+l];
					if (plan->flags&FFT5D_DEBUG) {
						printf("%f %f, ",a,b);
					} else {
						if (fabs(a-b)>2*NG[0]*NG[1]*NG[2]*FFT5D_EPS) {
							printf("result incorrect on %d,%d at %d,%d,%d: FFT5D:%f reference:%f\n",coor[0],coor[1],x,y,z,a,b);
						}
/*						assert(fabs(a-b)<2*NG[0]*NG[1]*NG[2]*FFT5D_EPS);*/
					}
				}
				if (plan->flags&FFT5D_DEBUG) printf(",");
			}
			if (plan->flags&FFT5D_DEBUG) printf("\n");
		}
	}
	
}

