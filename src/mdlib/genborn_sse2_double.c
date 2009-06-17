#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>

#include "typedefs.h"
#include "smalloc.h"
#include "genborn.h"
#include "vec.h"
#include "grompp.h"
#include "pdbio.h"
#include "names.h"
#include "physics.h"
#include "domdec.h"
#include "partdec.h"
#include "network.h"
#include "gmx_fatal.h"
#include "mtop_util.h"
#include "genborn.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREAD_MPI
#include "thread_mpi.h"
#endif

/* Only compile this file if SSE2 intrinsics are available */
#if ( (defined(GMX_IA32_SSE2) || defined(GMX_X86_64_SSE2) || defined(GMX_SSE2)) && defined(GMX_DOUBLE) )
#include <xmmintrin.h>
#include <emmintrin.h>

#if (defined (_MSC_VER) || defined(__INTEL_COMPILER))
#define gmx_castsi128_pd(a) _mm_castsi128_pd(a)
#define gmx_castpd_si128(a) _mm_castpd_si128(a)
#elif defined(__GNUC__)
#define gmx_castsi128_pd(a) ((__m128d)(a))
#define gmx_castpd_si128(a) ((__m128i)(a))
#else
static __m128d gmx_castsi128_pd(__m128i a) { return *(__m128d *) &a; } 
static __m128i gmx_castpd_si128(__m128d a) { return *(__m128i *) &a; } 
#endif



/* Still parameters - make sure to edit in genborn.c too if you change these! */
#define STILL_P1  0.073*0.1              /* length        */
#define STILL_P2  0.921*0.1*CAL2JOULE    /* energy*length */
#define STILL_P3  6.211*0.1*CAL2JOULE    /* energy*length */
#define STILL_P4  15.236*0.1*CAL2JOULE
#define STILL_P5  1.254 

#define STILL_P5INV (1.0/STILL_P5)
#define STILL_PIP5  (M_PI*STILL_P5)

#ifdef _MSC_VER /* visual c++ */
# define ALIGN16_BEG __declspec(align(16))
# define ALIGN16_END 
#else /* gcc or icc */
# define ALIGN16_BEG
# define ALIGN16_END __attribute__((aligned(16)))
#endif

#define _PS_CONST(Name, Val)                                            \
static const ALIGN16_BEG float _ps_##Name[4] ALIGN16_END = { Val, Val, Val, Val }
#define _PI32_CONST(Name, Val)                                            \
static const ALIGN16_BEG int _pi32_##Name[4] ALIGN16_END = { Val, Val, Val, Val }
#define _PS_CONST_TYPE(Name, Type, Val)                                 \
static const ALIGN16_BEG Type _ps_##Name[4] ALIGN16_END = { Val, Val, Val, Val }

/* Still parameters - make sure to edit in genborn.c too if you change these! */
#define STILL_P1  0.073*0.1              /* length        */
#define STILL_P2  0.921*0.1*CAL2JOULE    /* energy*length */
#define STILL_P3  6.211*0.1*CAL2JOULE    /* energy*length */
#define STILL_P4  15.236*0.1*CAL2JOULE
#define STILL_P5  1.254 

#define STILL_P5INV (1.0/STILL_P5)
#define STILL_PIP5  (M_PI*STILL_P5)



typedef
union 
{
	__m128 sse;
	float  f[4];
} my_sse_t;

typedef union xmm_mm_union {
	__m128 xmm;
	__m64 mm[2];
} xmm_mm_union;


__m128d log_pd(__m128d x)
{
	const __m128i exp_mask   = _mm_set_epi32(0x7FF00000,0,0x7FF00000,0);
	const __m128i exp_bias   = _mm_set_epi32(0,1023,0,1023);
	
	const __m128d const_loge = _mm_set1_pd(0.69314718055994529);
	const __m128d const_one  = _mm_set1_pd(1.0);
	const __m128d const_two  = _mm_set1_pd(2.0);
	/* Almost full single precision accuracy (~20 bits worst case) */
	const __m128d P0      = _mm_set1_pd(6.108179944792157749153050);
	const __m128d P1      = _mm_set1_pd(52.43691313715523327631139);
	const __m128d P2      = _mm_set1_pd(71.53664010795613671168440);
	const __m128d P3      = _mm_set1_pd(18.92097516931559547548485);
	const __m128d P4      = _mm_set1_pd(0.3504714784635941984522153);
	const __m128d P5      = _mm_set1_pd(-0.007105890734229368515879);
	const __m128d Q1      = _mm_set1_pd(17.73314231909420567454406);
	const __m128d Q2      = _mm_set1_pd(48.82373085428713023213363);
	const __m128d Q3      = _mm_set1_pd(31.65945943354513166309101);
	const __m128d Q4      = _mm_set1_pd(4.302477477108162270199051);
	
	__m128d xmm0,xmm1,xmm2,xmm3, xmm4;
	__m128i xmmi,xmmj;
	
	xmmi = gmx_castpd_si128(x);	
	xmm1 = _mm_cvtepi32_pd(_mm_shuffle_epi32(_mm_sub_epi64(_mm_srli_epi64(_mm_and_si128(xmmi, exp_mask), 52), exp_bias),_MM_SHUFFLE(3,1,2,0)));
	xmm0 = _mm_or_pd(gmx_castsi128_pd(_mm_andnot_si128(exp_mask, xmmi)), const_one);
	
	xmm2  = _mm_mul_pd(P5,xmm0);
	xmm2  = _mm_add_pd(xmm2,P4);
	xmm2  = _mm_mul_pd(xmm2,xmm0);
	xmm2  = _mm_add_pd(xmm2,P3);
	xmm2  = _mm_mul_pd(xmm2,xmm0);
	xmm2  = _mm_add_pd(xmm2,P2);
	xmm2  = _mm_mul_pd(xmm2,xmm0);
	xmm2  = _mm_add_pd(xmm2,P1);
	xmm2  = _mm_mul_pd(xmm2,xmm0);
	xmm2  = _mm_add_pd(xmm2,P0);
	
	xmm3  = _mm_mul_pd(Q4,xmm0);
	xmm3  = _mm_add_pd(xmm3,Q3);
	xmm3  = _mm_mul_pd(xmm3,xmm0);
	xmm3  = _mm_add_pd(xmm3,Q2);
	xmm3  = _mm_mul_pd(xmm3,xmm0);
	xmm3  = _mm_add_pd(xmm3,Q1);
	xmm3  = _mm_mul_pd(xmm3,xmm0);
	xmm3  = _mm_add_pd(xmm3,const_one);
	
	/* xmm4=1.0/xmm3 */
	xmm4 = _mm_cvtps_pd(_mm_rcp_ps(_mm_cvtpd_ps(xmm3)));
	xmm4 = _mm_mul_pd(xmm4,_mm_sub_pd(const_two,_mm_mul_pd(xmm3,xmm4)));
	xmm4 = _mm_mul_pd(xmm4,_mm_sub_pd(const_two,_mm_mul_pd(xmm3,xmm4)));
	xmm2 = _mm_mul_pd(xmm2,xmm4);
	
	xmm0  = _mm_sub_pd(xmm0, const_one);
	xmm0  = _mm_mul_pd(xmm0,xmm2);
	
	xmm0  = _mm_add_pd(xmm0,xmm1);
	
    return _mm_mul_pd(xmm0, const_loge);
}

/* This exp() routine provides accuracy of 10E-9 to 10E-11.
 * The polynomial minimax coefficients are actually accurate to 10E-14,
 * but we lose some accuracy in the polynomial evaluation.
 */
__m128d exp_pd(__m128d x)
{
    const __m128d lim1   = _mm_set1_pd(1025.0);   /* 1025.00000e+0d */
    const __m128d lim2   = _mm_set1_pd(-1022.99999999);   /* -1022.99999e+0f */
	
    const __m128i base   = _mm_set_epi32(0,0,1023,1023);
	const __m128d half   = _mm_set1_pd(0.5); 
	const __m128d log2e  = _mm_set1_pd(1.4426950408889634);
	
    const __m128d exp_P0 = _mm_set1_pd(1.00000000000001276211229749);
    const __m128d exp_P1 = _mm_set1_pd(6.931471805598709708635169E-1);
    const __m128d exp_P2 = _mm_set1_pd(2.402265069564965287455972E-1);
    const __m128d exp_P3 = _mm_set1_pd(5.550410866868561155599683E-2);
    const __m128d exp_P4 = _mm_set1_pd(9.618129192067246128919915E-3);
    const __m128d exp_P5 = _mm_set1_pd(1.333355761760444302342084E-3);
    const __m128d exp_P6 = _mm_set1_pd(1.540343494807179111289781E-4);
    const __m128d exp_P7 = _mm_set1_pd(1.525298483865349629325421E-5);
    const __m128d exp_P8 = _mm_set1_pd(1.325940560934510417818578E-6);
	const __m128d exp_P9 = _mm_set1_pd(1.015033670529892589443421E-7);
	
	__m128d xmm0,xmm1;
	__m128i xmmi;
	
	xmm0 = _mm_mul_pd(x,log2e);
	xmm0 = _mm_min_pd(xmm0,lim1);
	xmm0 = _mm_max_pd(xmm0,lim2);
	xmm1 = _mm_sub_pd(xmm0,half);
	
	xmmi = _mm_cvtpd_epi32(xmm1);
	xmm1 = _mm_cvtepi32_pd(xmmi);
	
	xmmi = _mm_add_epi32(xmmi,base);
	xmmi = _mm_shuffle_epi32(xmmi,_MM_SHUFFLE(3,1,2,0));
	xmmi = _mm_slli_epi64(xmmi,52);
	
	xmm0 = _mm_sub_pd(xmm0,xmm1);
	
	xmm1 = _mm_mul_pd(exp_P9,xmm0);
	xmm1 = _mm_add_pd(xmm1,exp_P8);
	xmm1 = _mm_mul_pd(xmm1,xmm0);
	xmm1 = _mm_add_pd(xmm1,exp_P7);
	xmm1 = _mm_mul_pd(xmm1,xmm0);
	xmm1 = _mm_add_pd(xmm1,exp_P6);
	xmm1 = _mm_mul_pd(xmm1,xmm0);
	xmm1 = _mm_add_pd(xmm1,exp_P5);
	xmm1 = _mm_mul_pd(xmm1,xmm0);
	xmm1 = _mm_add_pd(xmm1,exp_P4);
	xmm1 = _mm_mul_pd(xmm1,xmm0);
	xmm1 = _mm_add_pd(xmm1,exp_P3);
	xmm1 = _mm_mul_pd(xmm1,xmm0);
	xmm1 = _mm_add_pd(xmm1,exp_P2);
	xmm1 = _mm_mul_pd(xmm1,xmm0);
	xmm1 = _mm_add_pd(xmm1,exp_P1);
	xmm1 = _mm_mul_pd(xmm1,xmm0);
	xmm1 = _mm_add_pd(xmm1,exp_P0);
	xmm1 = _mm_mul_pd(xmm1,gmx_castsi128_pd(xmmi));
	
    return xmm1;
}


static inline __m128d
my_invrsq_pd(__m128d x)
{
	const __m128d three = {3.0, 3.0};
	const __m128d half  = {0.5, 0.5};
	
	__m128  t  = _mm_rsqrt_ps(_mm_cvtpd_ps(x)); /* Convert to single precision and do _mm_rsqrt_ps() */
	__m128d t1 = _mm_cvtps_pd(t); /* Convert back to double precision */
	
	/* First Newton-Rapson step, accuracy is now 24 bits */
	__m128d t2 = _mm_mul_pd(half,_mm_mul_pd(t1,_mm_sub_pd(three,_mm_mul_pd(x,_mm_mul_pd(t1,t1)))));
	
	/* Return second Newton-Rapson step, accuracy 48 bits */
	return _mm_mul_pd(half,_mm_mul_pd(t2,_mm_sub_pd(three,_mm_mul_pd(x,_mm_mul_pd(t2,t2)))));
}

static inline __m128d
my_inv_pd(__m128d x)
{
	const __m128d two = {2.0, 2.0};
	
	__m128  t  = _mm_rcp_ps(_mm_cvtpd_ps(x));
	__m128d t1 = _mm_cvtps_pd(t);
	__m128d t2 = _mm_mul_pd(t1,_mm_sub_pd(two,_mm_mul_pd(t1,x)));
	
	return _mm_mul_pd(t2,_mm_sub_pd(two,_mm_mul_pd(t2,x)));
}

#if 0
int
calc_gb_rad_still_sse2_double(t_commrec *cr, t_forcerec *fr,int natoms, gmx_mtop_t *mtop,
							  const t_atomtypes *atype, double *x, t_nblist *nl, gmx_genborn_t *born, t_mdatoms *md)
{
	int i,k,n,ai,ai3,aj1,aj2,aj13,aj23;
	int at0,at1,nj0,nj1,offset,taj1,taj2;

	double factor,gpi_ai,gpi_tmp,gpi2;
	double *sum_gpi;
	
	__m128d ix,iy,iz,jx,jy,jz,dx,dy,dz;
	__m128d t1,t2,t3,rsq11,rinv,rinv2,rinv4,rinv6;
	__m128d ratio,gpi,rai,raj,vaj,rvdw,mask_cmp;
	__m128d ccf,dccf,theta,cosq,term,sinq,res,prod;
	__m128d xmm1,xmm2,xmm3;
	__m128  tmp,sinq_single,cosq_single;
	
	const __m128d half  = {0.5, 0.5};
	const __m128d three = {3.0, 3.0};
	const __m128d one   = {1.0, 1.0};
	const __m128d two   = {2.0, 2.0};
	const __m128d zero  = {0.0, 0.0};
	const __m128d four  = {4.0, 4.0};
	
	const __m128d p5inv  = {STILL_P5INV, STILL_P5INV};
	const __m128d pip5   = {STILL_PIP5,  STILL_PIP5};
	const __m128d p4     = {STILL_P4,    STILL_P4};
	
	factor              = 0.5 * ONE_4PI_EPS0;
	sum_gpi             = born->work;
	n                   = 0;
	
	for(i=0;i<nl->nri;i++)
	{
		ai     = nl->iinr[i];
		ai3    = ai * 3;
		nj0    = nl->jindex[ai];
		nj1    = nl->jindex[ai+1];
		
		offset = (nj1-nj0)%2;
		
		/* Polarization energy for atom ai */
		gpi_ai = born->gpol[ai];
		gpi    = _mm_setzero_pd();
	
		/* Load particle ai coordinates */
		ix     = _mm_load1_pd(x+ai3);
		iy     = _mm_load1_pd(x+ai3+1);
		iz     = _mm_load1_pd(x+ai3+2);
		
		/* Load particle ai gb_radius */
		rai    = _mm_set1_pd(mtop->atomtypes.gb_radius[md->typeA[ai]]);
							 
		if(PAR(cr))
		{
			sum_gpi[ai] = 0;
		}
		
		for(k=nj0;k<nj1-offset;k+=2)
		{
			aj1      = nl->jjnr[k];
			aj2      = nl->jjnr[k+1];
			
			aj13     = aj1 * 3;
			aj23     = aj2 * 3;
			
			taj1     = md->typeA[aj1];
			taj2     = md->typeA[aj2];
			
			/* Load particle aj1-2 coordinates */
			xmm1     = _mm_loadu_pd(x+aj13);
			xmm2     = _mm_loadu_pd(x+aj23);
			jx       = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(0,0));
			jy       = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(1,1));
			
			jz       = _mm_loadl_pd(jz,x+aj13+2);
			jz       = _mm_loadh_pd(jz,x+aj23+2);
			
			dx       = _mm_sub_pd(ix,jx);
			dy       = _mm_sub_pd(iy,jy);
			dz       = _mm_sub_pd(iz,jz);
			
			rsq11    = _mm_add_pd( _mm_add_pd( _mm_mul_pd(dx,dx) , _mm_mul_pd(dy,dy) ) , _mm_mul_pd(dz,dz) );
			rinv     = my_invrsq_pd(rsq11);
			
			rinv2    = _mm_mul_pd(rinv,rinv);
			rinv4    = _mm_mul_pd(rinv2,rinv2);
			rinv6    = _mm_mul_pd(rinv4,rinv2);
			
			vaj      = _mm_loadl_pd(vaj,born->vsolv+aj1);
			vaj      = _mm_loadh_pd(vaj,born->vsolv+aj2);
			
			raj      = _mm_loadl_pd(raj, mtop->atomtypes.gb_radius+taj1);
			raj      = _mm_loadh_pd(raj, mtop->atomtypes.gb_radius+taj2);
			
			rvdw     = _mm_add_pd(rai,raj);
			rvdw     = _mm_mul_pd(rvdw,rvdw);
			ratio    = _mm_div_pd(rsq11,rvdw);
			
			mask_cmp = _mm_cmpgt_pd(ratio,p5inv);
								  
			switch(_mm_movemask_pd(mask_cmp))
			{
				case 0xF:
					ccf  = one;
					dccf = zero;
					break;
				default:
					theta = _mm_mul_pd(ratio,pip5);
				
					sincos_ps(_mm_cvtpd_ps(theta),&sinq_single,&cosq_single);
					sinq  = _mm_cvtps_pd(sinq_single);
					cosq  = _mm_cvtps_pd(cosq_single);
					
					term  = _mm_sub_pd(one,cosq);
					term  = _mm_mul_pd(half,term);
					ccf   = _mm_mul_pd(term,term);
					dccf  = _mm_mul_pd(two,term);
					dccf  = _mm_mul_pd(dccf,sinq);
					dccf  = _mm_mul_pd(dccf,pip5);
					dccf  = _mm_mul_pd(dccf,ratio);
					
					ccf   = _mm_or_pd(_mm_and_pd(mask_cmp,one)  ,_mm_andnot_pd(mask_cmp,ccf));
					dccf  = _mm_or_pd(_mm_and_pd(mask_cmp,zero) ,_mm_andnot_pd(mask_cmp,dccf));
			}
			
			prod      = _mm_mul_pd(p4,vaj);
			xmm2      = _mm_mul_pd(ccf,rinv4);
			xmm2      = _mm_mul_pd(xmm2,prod);
			gpi       = _mm_add_pd(gpi,xmm2);
			
			/* Chain rule terms */
			ccf       = _mm_mul_pd(four,ccf);
			xmm3      = _mm_sub_pd(ccf,dccf);
			xmm3      = _mm_mul_pd(xmm3,rinv6);
			xmm1      = _mm_mul_pd(xmm3,prod);
			
			_mm_storeu_pd(fr->dadx+n,xmm1);
			
			n         = n + 2;
				   
		}
			
		/* Deal with odd elements */
		if(offset!=0)
		{
			aj1  = nl->jjnr[k];
			aj13 = aj1 * 3;
			taj1 = md->typeA[aj1];
			
			jx    = _mm_load_sd(x+aj13);
			jy    = _mm_load_sd(x+aj13+1);
			jz    = _mm_load_sd(x+aj13+2);
			
			raj   = _mm_load_sd(mtop->atomtypes.gb_radius+taj1);
			vaj   = _mm_load_sd(born->vsolv+aj1);
						
			dx    = _mm_sub_sd(ix,jx);
			dy    = _mm_sub_pd(iy,jy);
			dz    = _mm_sub_pd(iz,jz);
			
			rsq11 = _mm_add_sd( _mm_add_sd( _mm_mul_sd(dx,dx) , _mm_mul_sd(dy,dy) ) , _mm_mul_sd(dz,dz) );
			rinv  = my_invrsq_pd(rsq11);
			
			rinv2 = _mm_mul_sd(rinv,rinv);
			rinv4 = _mm_mul_sd(rinv2,rinv2);
			rinv6 = _mm_mul_sd(rinv4,rinv2);
			
			rvdw  = _mm_add_sd(rai,raj);
			rvdw  = _mm_mul_sd(rvdw,rvdw);
			ratio = _mm_div_sd(rsq11,rvdw);
			
			mask_cmp = _mm_cmpgt_sd(ratio,p5inv);
			
			switch(_mm_movemask_pd(mask_cmp))
			{
				case 0xF:
					ccf  = one;
					dccf = zero;
					break;
				default:
					theta = _mm_mul_sd(ratio,pip5);
					
					sincos_ps(_mm_cvtpd_ps(theta),&sinq_single,&cosq_single);
					sinq  = _mm_cvtps_pd(sinq_single);
					cosq  = _mm_cvtps_pd(cosq_single);
					
					term  = _mm_sub_sd(one,cosq);
					term  = _mm_mul_sd(half,term);
					ccf   = _mm_mul_sd(term,term);
					dccf  = _mm_mul_sd(two,term);
					dccf  = _mm_mul_sd(dccf,sinq);
					dccf  = _mm_mul_sd(dccf,pip5);
					dccf  = _mm_mul_sd(dccf,ratio);
					
					ccf   = _mm_or_pd(_mm_and_pd(mask_cmp,one)  ,_mm_andnot_pd(mask_cmp,ccf));
					dccf  = _mm_or_pd(_mm_and_pd(mask_cmp,zero) ,_mm_andnot_pd(mask_cmp,dccf));
			}
			
			prod      = _mm_mul_sd(p4,vaj);
			xmm2      = _mm_mul_sd(ccf,rinv4);
			xmm2      = _mm_mul_sd(xmm2,prod);
			gpi       = _mm_add_sd(gpi,xmm2);
			
			/* Chain rule terms */
			ccf       = _mm_mul_sd(four,ccf);
			xmm3      = _mm_sub_sd(ccf,dccf);
			xmm3      = _mm_mul_sd(xmm3,rinv6);
			xmm1      = _mm_mul_sd(xmm3,prod);
			
			_mm_store_sd(fr->dadx+n,xmm1);
			
			n         = n + 1;
		} /* End offset */

		/* Do end processing ... */
		xmm2  = _mm_unpacklo_pd(xmm2,gpi);
		gpi   = _mm_add_pd(gpi,xmm2);
		gpi   = _mm_shuffle_pd(gpi,gpi,_MM_SHUFFLE2(1,1));
			
		_mm_store_sd(&gpi_tmp,gpi);
		
		if(PAR(cr))
		{
			sum_gpi[ai] = gpi_tmp;
		}
		else
		{
			gpi_ai = gpi_ai + gpi_tmp;
			gpi2   = gpi_ai * gpi_ai;
			
			born->bRad[ai]   = factor * invsqrt(gpi2);
			fr->invsqrta[ai] = invsqrt(born->bRad[ai]);
		}
	}

	if(PARTDECOMP(cr))
	{
		pd_at_range(cr,&at0,&at1);
		gmx_sum(natoms,sum_gpi,cr);
		
		for(i=at0;i<at1;i++)
		{
			ai = i;
			ai = i;
			gpi_ai = born->gpol[ai];
			gpi_ai = gpi_ai + sum_gpi[ai];
			gpi2   = gpi_ai*gpi_ai;
			
			born->bRad[ai]=factor*invsqrt(gpi2);
			fr->invsqrta[ai]=invsqrt(born->bRad[ai]);
		}
		
		/* Communicate Born radii*/
		gb_pd_send(cr,born->bRad,md->nr);
		gb_pd_send(cr,fr->invsqrta,md->nr);
	}
	else if(DOMAINDECOMP(cr))
	{
		dd_atom_sum_real(cr->dd,sum_gpi);
		
		for(i=0;i<nl->nri;i++)
		{
			ai     = nl->iinr[i];
			gpi_ai = born->gpol[cr->dd->gatindex[ai]];
			gpi_ai = gpi_ai + sum_gpi[i];
			gpi2   = gpi_ai*gpi_ai;
			
			born->bRad[ai]   = factor*invsqrt(gpi2);
			fr->invsqrta[ai] = invsqrt(born->bRad[ai]);
		}
		
		/* Communicate Born radii */
		dd_atom_spread_real(cr->dd,born->bRad);
		dd_atom_spread_real(cr->dd,fr->invsqrta);
	}
	
	
	return 0;
}

int 
calc_gb_rad_hct_sse2_double(t_commrec *cr, t_forcerec *fr, int natoms, gmx_mtop_t *mtop, 
							const t_atomtypes *atype, double *x, t_nblist *nl, gmx_genborn_t *born, t_mdatoms *md)
{
	int i,k,n,ai,ai3,aj1,aj2,aj13,aj23,nj0,nj1,at0,at1,offset;
	double ri,rr,sum,sum_tmp,sum_ai,min_rad,rad,doffset;
	double *sum_mpi;
	
	__m128d ix,iy,iz,jx,jy,jz,dx,dy,dz;
	__m128d t1,t2,t3,rsq11,rinv,r,rai;
	__m128d rai_inv,sk,sk2,lij,dlij,duij;
	__m128d uij,lij2,uij2,lij3,uij3,diff2;
	__m128d lij_inv,sk2_inv,prod,log_term,tmp,tmp_sum;
	__m128d mask_cmp,mask_cmp2,mask_cmp3;
	__m128d xmm1,xmm2,xmm3,xmm4,xmm8;
	__m128  log_term_single;
	
	const __m128d neg   = {-1.0, -1.0};
	const __m128d zero  = {0.0, 0.0};
	const __m128d eigth = {0.125, 0.125};
	const __m128d qrtr  = {0.25, 0.25};
	const __m128d half  = {0.5, 0.5};
	const __m128d one   = {1.0, 1.0};
	const __m128d two   = {2.0, 2.0};
	const __m128d three = {3.0, 3.0};
	
	if(PAR(cr))
	{
		pd_at_range(cr,&at0,&at1);
	}
	else
	{
		at0=0;
		at1=natoms;
	}
	
	doffset = born->gb_doffset;
	sum_mpi = born->work;
	n       = 0;
		
	for(i=0;i<nl->nri;i++)
	{
		ai      = nl->iinr[i];
		ai3     = ai * 3;
		
		nj0     = nl->jindex[ai];
		nj1     = nl->jindex[ai+1];
		
		offset  = (nj1-nj0)%2;
		
		/* Load rai */
		rr      = mtop->atomtypes.gb_radius[md->typeA[ai]]-doffset;
		rai     = _mm_load1_pd(&rr);
		
		/* Load cumulative sums for born radii calculation */
		sum     = 1.0/rr;
		rai_inv = _mm_load1_pd(&sum);
		tmp_sum = zero;
		
		/* Load ai coordinates */
		ix       = _mm_load1_pd(x+ai3);
		iy       = _mm_load1_pd(x+ai3+1);
		iz       = _mm_load1_pd(x+ai3+2);
		
		if(PAR(cr))
		{
			/* Only have the master node do this, since we only want one value at one time */
			if(MASTER(cr))
				sum_mpi[ai]=sum;
			else
				sum_mpi[ai]=0;
		}
		
		for(k=nj0;k<nj1-offset;k+=2)
		{
			aj1  = nl->jjnr[k];
			aj2  = nl->jjnr[k+1];
			
			aj13 = aj1 * 3;
			aj23 = aj2 * 3;
			
			/* Load particle aj1-2 coordinates */
			xmm1     = _mm_loadu_pd(x+aj13);
			xmm2     = _mm_loadu_pd(x+aj23);
			jx       = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(0,0));
			jy       = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(1,1));
			
			jz         = _mm_loadl_pd(jz, x+aj13+2);
			jz         = _mm_loadh_pd(jz, x+aj23+2);
			
			dx         = _mm_sub_pd(ix,jx);
			dy         = _mm_sub_pd(iy,jy);
			dz         = _mm_sub_pd(iz,jz);
			
			rsq11      = _mm_add_pd( _mm_add_pd( _mm_mul_pd(dx,dx) , _mm_mul_pd(dy,dy) ) , _mm_mul_pd(dz,dz) );
			rinv       = my_invrsq_pd(rsq11);
			r          = _mm_mul_pd(rinv,rsq11);
			
			sk         = _mm_loadl_pd(sk,born->param+aj1);
			sk         = _mm_loadh_pd(sk,born->param+aj2);
			
			/* conditional mask for rai<dr+sk */
			xmm1       = _mm_add_pd(r,sk);
			mask_cmp   = _mm_cmplt_pd(rai,xmm1);
			
			/* conditional for rai>dr-sk, ends with mask_cmp2 */
			xmm2      = _mm_sub_pd(r,sk);
			xmm3      = my_inv_pd(xmm2);
			mask_cmp2 = _mm_cmpgt_pd(rai,xmm2);
			
			lij	      = _mm_or_pd(_mm_and_pd(mask_cmp2,rai_inv)  ,_mm_andnot_pd(mask_cmp2,xmm3)); /*conditional as a mask*/
			dlij      = _mm_or_pd(_mm_and_pd(mask_cmp2,zero) ,_mm_andnot_pd(mask_cmp2,one));
			
			uij       = my_inv_pd(xmm1);
			lij2      = _mm_mul_pd(lij,lij);
			lij3      = _mm_mul_pd(lij2,lij);
			uij2      = _mm_mul_pd(uij,uij);
			uij3      = _mm_mul_pd(uij2,uij);
			
			diff2     = _mm_sub_pd(uij2,lij2);
			
			lij_inv   = my_invrsq_pd(lij2);
			sk2       = _mm_mul_pd(sk,sk);
			sk2_inv   = _mm_mul_pd(sk2,rinv);
			prod      = _mm_mul_pd(qrtr,sk2_inv);
			
			log_term  = _mm_mul_pd(uij,lij_inv);
			log_term_single  = log2_ps(_mm_cvtpd_ps(log_term)); /* Use single precision log */
			log_term  = _mm_cvtps_pd(log_term_single);
			
			xmm1      = _mm_sub_pd(lij,uij);
			xmm2      = _mm_mul_pd(qrtr,r);
			xmm2      = _mm_mul_pd(xmm2,diff2);
			xmm1      = _mm_add_pd(xmm1,xmm2); 
			xmm2      = _mm_mul_pd(half,rinv); 
			xmm2      = _mm_mul_pd(xmm2,log_term); 
			xmm1      = _mm_add_pd(xmm1,xmm2); 
			xmm8      = _mm_mul_pd(neg,diff2); 
			xmm2      = _mm_mul_pd(xmm8,prod); 
			tmp       = _mm_add_pd(xmm1,xmm2); 
			
			/* contitional for rai<sk-dr */
			xmm3      = _mm_sub_pd(sk,r);
			mask_cmp3 = _mm_cmplt_pd(rai,xmm3); /* rai<sk-dr */
			
			xmm4    = _mm_sub_pd(rai_inv,lij);
			xmm4    = _mm_mul_pd(two,xmm4);
			xmm4    = _mm_add_pd(tmp,xmm4);
			
			tmp	    = _mm_or_pd(_mm_and_pd(mask_cmp3,xmm4)  ,_mm_andnot_pd(mask_cmp3,tmp)); /*conditional as a mask*/
			
			/* the tmp will now contain four partial values, that not all are to be used. Which */
			/* ones are governed by the mask_cmp mask. */
			tmp     = _mm_mul_pd(half,tmp); 
			tmp     = _mm_or_pd(_mm_and_pd(mask_cmp,tmp)  ,_mm_andnot_pd(mask_cmp,zero)); /*conditional as a mask*/
			tmp_sum = _mm_sub_pd(tmp_sum,tmp);
						
			xmm2   = _mm_mul_pd(half,lij2); 
			xmm3   = _mm_mul_pd(prod,lij3); 
			xmm2   = _mm_add_pd(xmm2,xmm3); 
			xmm3   = _mm_mul_pd(lij,rinv); 
			xmm4   = _mm_mul_pd(lij3,r); 
			xmm3   = _mm_add_pd(xmm3,xmm4); 
			xmm3   = _mm_mul_pd(qrtr,xmm3); 
			t1     = _mm_sub_pd(xmm2,xmm3); 
			
			xmm2   = _mm_mul_pd(half,uij2);
			xmm2   = _mm_mul_pd(neg,xmm2); 
			xmm3   = _mm_mul_pd(qrtr,sk2_inv);
			xmm3   = _mm_mul_pd(xmm3,uij3); 
			xmm2   = _mm_sub_pd(xmm2,xmm3); 
			xmm3   = _mm_mul_pd(uij,rinv); 
			xmm4   = _mm_mul_pd(uij3,r); 
			xmm3   = _mm_add_pd(xmm3,xmm4); 
			xmm3   = _mm_mul_pd(qrtr,xmm3); 
			t2     = _mm_add_pd(xmm2,xmm3); 
			
			xmm2   = _mm_mul_pd(sk2_inv,rinv);
			xmm2   = _mm_add_pd(one,xmm2); 
			xmm2   = _mm_mul_pd(eigth,xmm2);
			xmm2   = _mm_mul_pd(xmm2,xmm8); 
			xmm3   = _mm_mul_pd(log_term, rinv);
			xmm3   = _mm_mul_pd(xmm3,rinv); 
			xmm3   = _mm_mul_pd(qrtr,xmm3); 
			t3     = _mm_add_pd(xmm2,xmm3); 
			
			/* chain rule terms */
			xmm2   = _mm_mul_pd(dlij,t1); 
			xmm2   = _mm_add_pd(xmm2,t2);
			xmm2   = _mm_add_pd(xmm2,t3); 
			xmm2   = _mm_mul_pd(xmm2,rinv);
			
			_mm_storeu_pd(fr->dadx+n,xmm2);
			
			n      =  n + 2;
		}
		
		if(offset!=0)
		{
			aj1       = nl->jjnr[k];
			aj13      = aj1 * 3;
			
			jx        = _mm_load_sd(x+aj13);
			jy        = _mm_load_sd(x+aj13+1);
			jz        = _mm_load_sd(x+aj13+2);
			
			sk        = _mm_load_sd(born->param+aj1);
			
			dx        = _mm_sub_sd(ix,jx);
			dy        = _mm_sub_pd(iy,jy);
			dz        = _mm_sub_pd(iz,jz);
			
			rsq11     = _mm_add_sd( _mm_add_sd( _mm_mul_sd(dx,dx) , _mm_mul_sd(dy,dy) ) , _mm_mul_sd(dz,dz) );
			rinv      = my_invrsq_pd(rsq11);
			r         = _mm_mul_sd(rinv,rsq11);
			
			/* conditional mask for rai<dr+sk */
			xmm1      = _mm_add_sd(r,sk);
			mask_cmp  = _mm_cmplt_sd(rai,xmm1);
			
			/* conditional for rai>dr-sk, ends with mask_cmp2 */
			xmm2      = _mm_sub_sd(r,sk);
			xmm3      = my_inv_pd(xmm2);
			mask_cmp2 = _mm_cmpgt_sd(rai,xmm2);
			
			lij	      = _mm_or_pd(_mm_and_pd(mask_cmp2,rai_inv)  ,_mm_andnot_pd(mask_cmp2,xmm3)); /*conditional as a mask*/
			dlij      = _mm_or_pd(_mm_and_pd(mask_cmp2,zero) ,_mm_andnot_pd(mask_cmp2,one));
			
			uij       = my_inv_pd(xmm1);
			lij2      = _mm_mul_sd(lij,lij);
			lij3      = _mm_mul_sd(lij2,lij);
			uij2      = _mm_mul_sd(uij,uij);
			uij3      = _mm_mul_sd(uij2,uij);
			
			diff2     = _mm_sub_sd(uij2,lij2);
			
			lij_inv   = my_invrsq_pd(lij2);
			sk2       = _mm_mul_sd(sk,sk);
			sk2_inv   = _mm_mul_sd(sk2,rinv);
			prod      = _mm_mul_sd(qrtr,sk2_inv);
		
			log_term  = _mm_mul_sd(uij,lij_inv);
			log_term_single  = log2_ps(_mm_cvtpd_ps(log_term)); /* Use single precision log */
			log_term  = _mm_cvtps_pd(log_term_single);
			
			xmm1      = _mm_sub_sd(lij,uij);
			xmm2      = _mm_mul_sd(qrtr,r);
			xmm2      = _mm_mul_sd(xmm2,diff2);
			xmm1      = _mm_add_sd(xmm1,xmm2); 
			xmm2      = _mm_mul_sd(half,rinv); 
			xmm2      = _mm_mul_sd(xmm2,log_term); 
			xmm1      = _mm_add_sd(xmm1,xmm2); 
			xmm8      = _mm_mul_sd(neg,diff2); 
			xmm2      = _mm_mul_sd(xmm8,prod); 
			tmp       = _mm_add_sd(xmm1,xmm2); 
			
			/* contitional for rai<sk-dr */
			xmm3      = _mm_sub_sd(sk,r);
			mask_cmp3 = _mm_cmplt_sd(rai,xmm3); /* rai<sk-dr */
			
			xmm4    = _mm_sub_sd(rai_inv,lij);
			xmm4    = _mm_mul_sd(two,xmm4);
			xmm4    = _mm_add_sd(tmp,xmm4);
			
			tmp	    = _mm_or_pd(_mm_and_pd(mask_cmp3,xmm4)  ,_mm_andnot_pd(mask_cmp3,tmp)); /*conditional as a mask*/
			
			/* the tmp will now contain four partial values, that not all are to be used. Which */
			/* ones are governed by the mask_cmp mask. */
			tmp     = _mm_mul_pd(half,tmp); 
			tmp     = _mm_or_pd(_mm_and_pd(mask_cmp,tmp)  ,_mm_andnot_pd(mask_cmp,zero)); /*conditional as a mask*/
			tmp_sum = _mm_sub_sd(tmp_sum,tmp);
					
			xmm2   = _mm_mul_sd(half,lij2); 
			xmm3   = _mm_mul_sd(prod,lij3); 
			xmm2   = _mm_add_sd(xmm2,xmm3); 
			xmm3   = _mm_mul_sd(lij,rinv); 
			xmm4   = _mm_mul_sd(lij3,r); 
			xmm3   = _mm_add_sd(xmm3,xmm4); 
			xmm3   = _mm_mul_sd(qrtr,xmm3); 
			t1     = _mm_sub_sd(xmm2,xmm3); 
			
			xmm2   = _mm_mul_sd(half,uij2);
			xmm2   = _mm_mul_sd(neg,xmm2); 
			xmm3   = _mm_mul_sd(qrtr,sk2_inv);
			xmm3   = _mm_mul_sd(xmm3,uij3); 
			xmm2   = _mm_sub_sd(xmm2,xmm3); 
			xmm3   = _mm_mul_sd(uij,rinv); 
			xmm4   = _mm_mul_sd(uij3,r); 
			xmm3   = _mm_add_sd(xmm3,xmm4); 
			xmm3   = _mm_mul_sd(qrtr,xmm3); 
			t2     = _mm_add_sd(xmm2,xmm3); 
			
			xmm2   = _mm_mul_sd(sk2_inv,rinv);
			xmm2   = _mm_add_sd(one,xmm2); 
			xmm2   = _mm_mul_sd(eigth,xmm2);
			xmm2   = _mm_mul_sd(xmm2,xmm8); 
			xmm3   = _mm_mul_sd(log_term, rinv);
			xmm3   = _mm_mul_sd(xmm3,rinv); 
			xmm3   = _mm_mul_sd(qrtr,xmm3); 
			t3     = _mm_add_sd(xmm2,xmm3); 
			
			/* chain rule terms */
			xmm2   = _mm_mul_sd(dlij,t1); 
			xmm2   = _mm_add_sd(xmm2,t2);
			xmm2   = _mm_add_sd(xmm2,t3); 
			xmm2   = _mm_mul_sd(xmm2,rinv);
			
			_mm_store_sd(fr->dadx+n,xmm2);
			
			n      =  n + 1;
		}
		
		/* Do end processing ...  */
		tmp     = _mm_unpacklo_pd(tmp,tmp_sum);
		tmp_sum = _mm_add_pd(tmp_sum,tmp);
		tmp_sum = _mm_shuffle_pd(tmp_sum,tmp_sum,_MM_SHUFFLE2(1,1));
		
		_mm_store_sd(&sum_tmp,tmp_sum);
		
		if(PAR(cr))
		{
			sum_mpi[ai]=sum+sum_tmp;
		}
		else
		{
			sum_ai=sum+sum_tmp; 
			
			min_rad = rr + doffset;
			rad=1.0/sum_ai; 
			
			born->bRad[ai]=rad > min_rad ? rad : min_rad;
			fr->invsqrta[ai]=invsqrt(born->bRad[ai]);
		}
	}
	
	if(PAR(cr))
	{
		gmx_sum(natoms,sum_mpi,cr);
		
		for(i=at0;i<at1;i++)
		{
			ai      = i;
			min_rad = mtop->atomtypes.gb_radius[md->typeA[ai]]; 
			rad     = 1.0/sum_mpi[ai];
			
			born->bRad[ai]=rad > min_rad ? rad : min_rad;
			fr->invsqrta[ai]=invsqrt(born->bRad[ai]);
		}
		
		gb_pd_send(cr,born->bRad,md->nr);
		gb_pd_send(cr,fr->invsqrta,md->nr);
		
	}
	
	return 0;
}

int 
calc_gb_rad_obc_sse2_double(t_commrec *cr, t_forcerec * fr, int natoms, gmx_mtop_t *mtop,
							const t_atomtypes *atype, double *x, t_nblist *nl, gmx_genborn_t *born,t_mdatoms *md)
{
	int i,k,n,ai,ai3,aj1,aj2,aj13,aj23,nj0,nj1,at0,at1,offset;
	double ri,rr,sum,sum_tmp,sum_ai,sum_ai2,sum_ai3,min_rad,rad,doffset;
	double tsum,tchain,rr_inv,gbr;
	double *sum_mpi;
	
	__m128d ix,iy,iz,jx,jy,jz,dx,dy,dz;
	__m128d t1,t2,t3,rsq11,rinv,r,rai;
	__m128d rai_inv,sk,sk2,lij,dlij,duij;
	__m128d uij,lij2,uij2,lij3,uij3,diff2;
	__m128d lij_inv,sk2_inv,prod,log_term,tmp,tmp_sum;
	__m128d mask_cmp,mask_cmp2,mask_cmp3;
	__m128d xmm1,xmm2,xmm3,xmm4,xmm8;
	__m128 log_term_single;
	
	const __m128d neg   = {-1.0, -1.0};
	const __m128d zero  = {0.0, 0.0};
	const __m128d eigth = {0.125, 0.125};
	const __m128d qrtr  = {0.25, 0.25};
	const __m128d half  = {0.5, 0.5};
	const __m128d one   = {1.0, 1.0};
	const __m128d two   = {2.0, 2.0};
	const __m128d three = {3.0, 3.0};
	
	if(PAR(cr))
	{
		pd_at_range(cr,&at0,&at1);
	}
	else
	{
		at0=0;
		at1=natoms;
	}
	
	doffset = born->gb_doffset;
	sum_mpi = born->work;
	n       = 0;
	
	for(i=0;i<nl->nri;i++)
	{
		ai      = nl->iinr[i];
		ai3     = ai * 3;
		
		nj0     = nl->jindex[ai];
		nj1     = nl->jindex[ai+1];
		
		offset  = (nj1-nj0)%2;
		
		/* Load rai */
		rr      = mtop->atomtypes.gb_radius[md->typeA[ai]];
		ri      = rr - doffset;
		rai     = _mm_load1_pd(&ri);
		
		ri      = 1.0/ri;
		rai_inv = _mm_load1_pd(&ri);
		
		/* Load ai coordinates */
		ix       = _mm_load1_pd(x+ai3);
		iy       = _mm_load1_pd(x+ai3+1);
		iz       = _mm_load1_pd(x+ai3+2);
		
		tmp_sum  = _mm_setzero_pd();
		
		if(PAR(cr))
		{
			sum_mpi[ai] = 0;
		}
		
		for(k=nj0;k<nj1-offset;k+=2)
		{
			aj1  = nl->jjnr[k];
			aj2  = nl->jjnr[k+1];
			
			aj13 = aj1 * 3;
			aj23 = aj2 * 3;
			
			/* Load particle aj1-2 coordinates */
			xmm1     = _mm_loadu_pd(x+aj13);
			xmm2     = _mm_loadu_pd(x+aj23);
			jx       = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(0,0));
			jy       = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(1,1));
			
			jz         = _mm_loadl_pd(jz, x+aj13+2);
			jz         = _mm_loadh_pd(jz, x+aj23+2);
			
			dx         = _mm_sub_pd(ix,jx);
			dy         = _mm_sub_pd(iy,jy);
			dz         = _mm_sub_pd(iz,jz);
			
			rsq11      = _mm_add_pd( _mm_add_pd( _mm_mul_pd(dx,dx) , _mm_mul_pd(dy,dy) ) , _mm_mul_pd(dz,dz) );
			rinv       = my_invrsq_pd(rsq11);
			r          = _mm_mul_pd(rinv,rsq11);
			
			sk         = _mm_loadl_pd(sk,born->param+aj1);
			sk         = _mm_loadh_pd(sk,born->param+aj2);
			
			/* conditional mask for rai<dr+sk */
			xmm1       = _mm_add_pd(r,sk);
			mask_cmp   = _mm_cmplt_pd(rai,xmm1);
			
			/* conditional for rai>dr-sk, ends with mask_cmp2 */
			xmm2      = _mm_sub_pd(r,sk);
			xmm3      = my_inv_pd(xmm2);
			mask_cmp2 = _mm_cmpgt_pd(rai,xmm2);
			
			lij	      = _mm_or_pd(_mm_and_pd(mask_cmp2,rai_inv)  ,_mm_andnot_pd(mask_cmp2,xmm3)); /*conditional as a mask*/
			dlij      = _mm_or_pd(_mm_and_pd(mask_cmp2,zero) ,_mm_andnot_pd(mask_cmp2,one));
			
			uij       = my_inv_pd(xmm1);
			lij2      = _mm_mul_pd(lij,lij);
			lij3      = _mm_mul_pd(lij2,lij);
			uij2      = _mm_mul_pd(uij,uij);
			uij3      = _mm_mul_pd(uij2,uij);
			
			diff2     = _mm_sub_pd(uij2,lij2);
			
			lij_inv   = my_invrsq_pd(lij2);
			sk2       = _mm_mul_pd(sk,sk);
			sk2_inv   = _mm_mul_pd(sk2,rinv);
			prod      = _mm_mul_pd(qrtr,sk2_inv);
			
			log_term  = _mm_mul_pd(uij,lij_inv);
			log_term_single  = log2_ps(_mm_cvtpd_ps(log_term)); /* Use single precision log */
			log_term  = _mm_cvtps_pd(log_term_single);
			
			xmm1      = _mm_sub_pd(lij,uij);
			xmm2      = _mm_mul_pd(qrtr,r);
			xmm2      = _mm_mul_pd(xmm2,diff2);
			xmm1      = _mm_add_pd(xmm1,xmm2); 
			xmm2      = _mm_mul_pd(half,rinv); 
			xmm2      = _mm_mul_pd(xmm2,log_term); 
			xmm1      = _mm_add_pd(xmm1,xmm2); 
			xmm8      = _mm_mul_pd(neg,diff2); 
			xmm2      = _mm_mul_pd(xmm8,prod); 
			tmp       = _mm_add_pd(xmm1,xmm2); 
			
			/* contitional for rai<sk-dr */
			xmm3      = _mm_sub_pd(sk,r);
			mask_cmp3 = _mm_cmplt_pd(rai,xmm3); /* rai<sk-dr */
			
			xmm4    = _mm_sub_pd(rai_inv,lij);
			xmm4    = _mm_mul_pd(two,xmm4);
			xmm4    = _mm_add_pd(tmp,xmm4);
			
			tmp	    = _mm_or_pd(_mm_and_pd(mask_cmp3,xmm4)  ,_mm_andnot_pd(mask_cmp3,tmp)); /*conditional as a mask*/
			
			/* the tmp will now contain four partial values, that not all are to be used. Which */
			/* ones are governed by the mask_cmp mask. */
			tmp     = _mm_mul_pd(half,tmp); 
			tmp     = _mm_or_pd(_mm_and_pd(mask_cmp,tmp)  ,_mm_andnot_pd(mask_cmp,zero)); /*conditional as a mask*/
			tmp_sum = _mm_add_pd(tmp_sum,tmp);
			
			xmm2   = _mm_mul_pd(half,lij2); 
			xmm3   = _mm_mul_pd(prod,lij3); 
			xmm2   = _mm_add_pd(xmm2,xmm3); 
			xmm3   = _mm_mul_pd(lij,rinv); 
			xmm4   = _mm_mul_pd(lij3,r); 
			xmm3   = _mm_add_pd(xmm3,xmm4); 
			xmm3   = _mm_mul_pd(qrtr,xmm3); 
			t1     = _mm_sub_pd(xmm2,xmm3); 
			
			xmm2   = _mm_mul_pd(half,uij2);
			xmm2   = _mm_mul_pd(neg,xmm2); 
			xmm3   = _mm_mul_pd(qrtr,sk2_inv);
			xmm3   = _mm_mul_pd(xmm3,uij3); 
			xmm2   = _mm_sub_pd(xmm2,xmm3); 
			xmm3   = _mm_mul_pd(uij,rinv); 
			xmm4   = _mm_mul_pd(uij3,r); 
			xmm3   = _mm_add_pd(xmm3,xmm4); 
			xmm3   = _mm_mul_pd(qrtr,xmm3); 
			t2     = _mm_add_pd(xmm2,xmm3); 
			
			xmm2   = _mm_mul_pd(sk2_inv,rinv);
			xmm2   = _mm_add_pd(one,xmm2); 
			xmm2   = _mm_mul_pd(eigth,xmm2);
			xmm2   = _mm_mul_pd(xmm2,xmm8); 
			xmm3   = _mm_mul_pd(log_term, rinv);
			xmm3   = _mm_mul_pd(xmm3,rinv); 
			xmm3   = _mm_mul_pd(qrtr,xmm3); 
			t3     = _mm_add_pd(xmm2,xmm3); 
			
			/* chain rule terms */
			xmm2   = _mm_mul_pd(dlij,t1); 
			xmm2   = _mm_add_pd(xmm2,t2);
			xmm2   = _mm_add_pd(xmm2,t3); 
			xmm2   = _mm_mul_pd(xmm2,rinv);
			
			_mm_storeu_pd(fr->dadx+n,xmm2);
			
			n      =  n + 2;
		}
		if(offset!=0)
		{
			aj1       = nl->jjnr[k];
			aj13      = aj1 * 3;
			
			jx        = _mm_load_sd(x+aj13);
			jy        = _mm_load_sd(x+aj13+1);
			jz        = _mm_load_sd(x+aj13+2);
			
			sk        = _mm_load_sd(born->param+aj1);
			
			dx        = _mm_sub_sd(ix,jx);
			dy        = _mm_sub_pd(iy,jy);
			dz        = _mm_sub_pd(iz,jz);
			
			rsq11     = _mm_add_sd( _mm_add_sd( _mm_mul_sd(dx,dx) , _mm_mul_sd(dy,dy) ) , _mm_mul_sd(dz,dz) );
			rinv      = my_invrsq_pd(rsq11);
			r         = _mm_mul_sd(rinv,rsq11);
			
			/* conditional mask for rai<dr+sk */
			xmm1      = _mm_add_sd(r,sk);
			mask_cmp  = _mm_cmplt_sd(rai,xmm1);
			
			/* conditional for rai>dr-sk, ends with mask_cmp2 */
			xmm2      = _mm_sub_sd(r,sk);
			xmm3      = my_inv_pd(xmm2);
			mask_cmp2 = _mm_cmpgt_sd(rai,xmm2);
			
			lij	      = _mm_or_pd(_mm_and_pd(mask_cmp2,rai_inv)  ,_mm_andnot_pd(mask_cmp2,xmm3)); /*conditional as a mask*/
			dlij      = _mm_or_pd(_mm_and_pd(mask_cmp2,zero) ,_mm_andnot_pd(mask_cmp2,one));
			
			uij       = my_inv_pd(xmm1);
			lij2      = _mm_mul_sd(lij,lij);
			lij3      = _mm_mul_sd(lij2,lij);
			uij2      = _mm_mul_sd(uij,uij);
			uij3      = _mm_mul_sd(uij2,uij);
			
			diff2     = _mm_sub_sd(uij2,lij2);
			
			lij_inv   = my_invrsq_pd(lij2);
			sk2       = _mm_mul_sd(sk,sk);
			sk2_inv   = _mm_mul_sd(sk2,rinv);
			prod      = _mm_mul_sd(qrtr,sk2_inv);
			
			log_term  = _mm_mul_pd(uij,lij_inv);
			log_term_single  = log2_ps(_mm_cvtpd_ps(log_term)); /* Use single precision log */
			log_term  = _mm_cvtps_pd(log_term_single);
			
			xmm1      = _mm_sub_sd(lij,uij);
			xmm2      = _mm_mul_sd(qrtr,r);
			xmm2      = _mm_mul_sd(xmm2,diff2);
			xmm1      = _mm_add_sd(xmm1,xmm2); 
			xmm2      = _mm_mul_sd(half,rinv); 
			xmm2      = _mm_mul_sd(xmm2,log_term); 
			xmm1      = _mm_add_sd(xmm1,xmm2); 
			xmm8      = _mm_mul_sd(neg,diff2); 
			xmm2      = _mm_mul_sd(xmm8,prod); 
			tmp       = _mm_add_sd(xmm1,xmm2); 
			
			/* contitional for rai<sk-dr */
			xmm3      = _mm_sub_sd(sk,r);
			mask_cmp3 = _mm_cmplt_sd(rai,xmm3); /* rai<sk-dr */
			
			xmm4    = _mm_sub_sd(rai_inv,lij);
			xmm4    = _mm_mul_sd(two,xmm4);
			xmm4    = _mm_add_sd(tmp,xmm4);
			
			tmp	    = _mm_or_pd(_mm_and_pd(mask_cmp3,xmm4)  ,_mm_andnot_pd(mask_cmp3,tmp)); /*conditional as a mask*/
			
			/* the tmp will now contain four partial values, that not all are to be used. Which */
			/* ones are governed by the mask_cmp mask. */
			tmp     = _mm_mul_pd(half,tmp); 
			tmp     = _mm_or_pd(_mm_and_pd(mask_cmp,tmp)  ,_mm_andnot_pd(mask_cmp,zero)); /*conditional as a mask*/
			tmp_sum = _mm_add_sd(tmp_sum,tmp);
			
			xmm2   = _mm_mul_sd(half,lij2); 
			xmm3   = _mm_mul_sd(prod,lij3); 
			xmm2   = _mm_add_sd(xmm2,xmm3); 
			xmm3   = _mm_mul_sd(lij,rinv); 
			xmm4   = _mm_mul_sd(lij3,r); 
			xmm3   = _mm_add_sd(xmm3,xmm4); 
			xmm3   = _mm_mul_sd(qrtr,xmm3); 
			t1     = _mm_sub_sd(xmm2,xmm3); 
			
			xmm2   = _mm_mul_sd(half,uij2);
			xmm2   = _mm_mul_sd(neg,xmm2); 
			xmm3   = _mm_mul_sd(qrtr,sk2_inv);
			xmm3   = _mm_mul_sd(xmm3,uij3); 
			xmm2   = _mm_sub_sd(xmm2,xmm3); 
			xmm3   = _mm_mul_sd(uij,rinv); 
			xmm4   = _mm_mul_sd(uij3,r); 
			xmm3   = _mm_add_sd(xmm3,xmm4); 
			xmm3   = _mm_mul_sd(qrtr,xmm3); 
			t2     = _mm_add_sd(xmm2,xmm3); 
			
			xmm2   = _mm_mul_sd(sk2_inv,rinv);
			xmm2   = _mm_add_sd(one,xmm2); 
			xmm2   = _mm_mul_sd(eigth,xmm2);
			xmm2   = _mm_mul_sd(xmm2,xmm8); 
			xmm3   = _mm_mul_sd(log_term, rinv);
			xmm3   = _mm_mul_sd(xmm3,rinv); 
			xmm3   = _mm_mul_sd(qrtr,xmm3); 
			t3     = _mm_add_sd(xmm2,xmm3); 
			
			/* chain rule terms */
			xmm2   = _mm_mul_sd(dlij,t1); 
			xmm2   = _mm_add_sd(xmm2,t2);
			xmm2   = _mm_add_sd(xmm2,t3); 
			xmm2   = _mm_mul_sd(xmm2,rinv);
			
			_mm_store_sd(fr->dadx+n,xmm2);
			
			n      =  n + 1;
		}
		
		/* Do end processing ...  */
		tmp     = _mm_unpacklo_pd(tmp,tmp_sum);
		tmp_sum = _mm_add_pd(tmp_sum,tmp);
		tmp_sum = _mm_shuffle_pd(tmp_sum,tmp_sum,_MM_SHUFFLE2(1,1));
			
		_mm_store_sd(&sum_tmp,tmp_sum);
		
		if(PAR(cr))
		{
			sum_mpi[ai]=sum_tmp;
		}
		else
		{
			sum_ai=sum_tmp;
			
			/* calculate the born radii */
			sum_ai  = (rr-doffset) * sum_ai;
			sum_ai2 = sum_ai       * sum_ai;
			sum_ai3 = sum_ai2      * sum_ai;
			
			tsum    = tanh(born->obc_alpha*sum_ai-born->obc_beta*sum_ai2+born->obc_gamma*sum_ai3);
			born->bRad[ai] = ri - tsum/rr;
			born->bRad[ai] = 1.0/(born->bRad[ai]);
			
			fr->invsqrta[ai] = invsqrt(born->bRad[ai]);
			tchain = (rr-doffset)*(born->obc_alpha-2*born->obc_beta*sum_ai+3*born->obc_gamma*sum_ai2);
			born->drobc[ai] = (1.0 - tsum*tsum)*tchain*(1.0/rr);
		}
		
	}
	
	if(PAR(cr))
	{
		gmx_sum(natoms,sum_mpi,cr);
		
		for(i=at0;i<at1;i++)
		{
			ai      = i;
			rr      = mtop->atomtypes.gb_radius[md->typeA[ai]];
			rr_inv  = 1.0/rr;
			
			sum_ai  = sum_mpi[ai];
			sum_ai  = (rr-doffset) * sum_ai;
			sum_ai2 = sum_ai       * sum_ai;
			sum_ai3 = sum_ai2      * sum_ai;
			
			tsum    = tanh(born->obc_alpha*sum_ai-born->obc_beta*sum_ai2+born->obc_gamma*sum_ai3);
			gbr     = 1.0/(rr-doffset) - tsum*rr_inv;
			
			born->bRad[ai] = 1.0 / gbr;
			fr->invsqrta[ai]=invsqrt(born->bRad[ai]);
			
			tchain  = (rr-doffset) * (born->obc_alpha-2*born->obc_beta*sum_ai+3*born->obc_gamma*sum_ai2);
			born->drobc[ai] = (1.0-tsum*tsum)*tchain*rr_inv;
		}
		
		gb_pd_send(cr,born->bRad,md->nr);
		gb_pd_send(cr,fr->invsqrta,md->nr);
		gb_pd_send(cr,born->drobc,md->nr);
	}
	
	return 0;

}
#endif

int
calc_gb_chainrule_sse2_double(int natoms, t_nblist *nl, double *dadx, double *dvda, double *xd, double *f, int gb_algorithm, gmx_genborn_t *born)
{
	int i,k,n,ai,aj,ai3,aj1,aj2,aj13,aj23,aj4,nj0,nj1,offset;
	double rbi;
	double *rb;
	
	__m128d ix,iy,iz,jx,jy,jz,fix,fiy,fiz;
	__m128d dx,dy,dz,t1,t2,t3,dva,dax,fgb;
	__m128d xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;
	
	t1 = t2 = t3 = _mm_setzero_pd();
	rb = born->work;
	
	/* Loop to get the proper form for the Born radius term */
	if(gb_algorithm==egbSTILL) {
		for(i=0;i<natoms;i++)
		{
			rbi   = born->bRad[i];
			rb[i] = (2 * rbi * rbi * dvda[i])/ONE_4PI_EPS0;
		}
	}
	
	else if(gb_algorithm==egbHCT) {
		for(i=0;i<natoms;i++)
		{
			rbi   = born->bRad[i];
			rb[i] = rbi * rbi * dvda[i];
		}
	}
	
	else if(gb_algorithm==egbOBC) {
		for(i=0;i<natoms;i++)
		{
			rbi   = born->bRad[i];
			rb[i] = rbi * rbi * born->drobc[i] * dvda[i];
		}
	}
	
	n = 0;
	
	for(i=0;i<nl->nri;i++)
	{
		ai     = nl->iinr[i];
		ai3    = ai * 3;
		
		nj0    = nl->jindex[ai];
		nj1    = nl->jindex[ai+1];
		
		offset = (nj1-nj0)%2;
		
		/* Load particle ai coordinates */
		ix  = _mm_load1_pd(xd+ai3);
		iy  = _mm_load1_pd(xd+ai3+1);
		iz  = _mm_load1_pd(xd+ai3+2);
		
		/* Load particle ai dvda */
		dva = _mm_load1_pd(rb+ai);
		
		fix = _mm_setzero_pd();
		fiy = _mm_setzero_pd();
		fiz = _mm_setzero_pd();	
		
		for(k=nj0;k<nj1-offset;k+=2)
		{
			aj1 = nl->jjnr[k];
			aj2 = nl->jjnr[k+1];
			
			aj13 = aj1 * 3;
			aj23 = aj2 * 3;
			
			/* Load particle aj1-2 coordinates */
			xmm1 = _mm_loadu_pd(xd+aj13);
			xmm2 = _mm_loadu_pd(xd+aj23);
			
			xmm5 = _mm_load_sd(xd+aj13+2);
			xmm6 = _mm_load_sd(xd+aj23+2);
			
			jx   = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(0,0));
			jy   = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(1,1));
			jz   = _mm_shuffle_pd(xmm5,xmm6,_MM_SHUFFLE2(0,0));
					
			/* load derivative of born radii w.r.t. coordinates */
			dax  = _mm_loadu_pd(dadx+n);
			n    = n + 2;
			
			/* distances */
			dx   = _mm_sub_pd(ix,jx);
			dy   = _mm_sub_pd(iy,jy);
			dz   = _mm_sub_pd(iz,jz);

			/* scalar force */
			fgb  = _mm_mul_pd(dva,dax); 
		
			/* partial forces */
			t1   = _mm_mul_pd(fgb,dx); 
			t2   = _mm_mul_pd(fgb,dy); 
			t3   = _mm_mul_pd(fgb,dz); 
			
			/* update the i force */
			fix  = _mm_add_pd(fix,t1);
			fiy  = _mm_add_pd(fiy,t2);
			fiz  = _mm_add_pd(fiz,t3);	
			
			/* accumulate forces from memory */
			xmm1 = _mm_loadu_pd(f+aj13); /* fx1 fx2 */
			xmm2 = _mm_loadu_pd(f+aj23); /* fy1 fy2 */
			
			xmm5 = _mm_load1_pd(f+aj13+2); /* fz1 fz1 */
			xmm6 = _mm_load1_pd(f+aj23+2); /* fz2 fz2 */
			
			/* transpose */
			xmm7 = _mm_shuffle_pd(xmm5,xmm6,_MM_SHUFFLE2(0,0)); /* fz1 fz2 */
			xmm5 = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(0,0)); /* fx1 fx2 */
			xmm6 = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(1,1)); /* fy1 fy2 */
						
			/* subtract partial forces */
			xmm5 = _mm_sub_pd(xmm5, t1); 
			xmm6 = _mm_sub_pd(xmm6, t2); 
			xmm7 = _mm_sub_pd(xmm7, t3); 
			
			/* transpose back fx and fy */
			xmm1 = _mm_shuffle_pd(xmm5,xmm6,_MM_SHUFFLE2(0,0));
			xmm2 = _mm_shuffle_pd(xmm5,xmm6,_MM_SHUFFLE2(1,1));
			
			/* store the force, first fx's and fy's */
			_mm_storeu_pd(f+aj13,xmm1);
			_mm_storeu_pd(f+aj23,xmm2);
			
			/* then fz */
			_mm_storel_pd(f+aj13+2,xmm7);
			_mm_storeh_pd(f+aj23+2,xmm7);
		}

		if(offset!=0)
		{
			aj1  = nl->jjnr[k];
			aj13 = aj1 * 3;
			
			jx    = _mm_load_sd(xd+aj13);
			jy    = _mm_load_sd(xd+aj13+1);
			jz    = _mm_load_sd(xd+aj13+2);
			
			dax  = _mm_load_sd(dadx+n);
			n    = n + 1;
			
			dx   = _mm_sub_sd(ix,jx);
			dy   = _mm_sub_sd(iy,jy);
			dz   = _mm_sub_sd(iz,jz);
			
			fgb  = _mm_mul_sd(dva,dax); 
			
			t1   = _mm_mul_sd(fgb,dx); 
			t2   = _mm_mul_sd(fgb,dy); 
			t3   = _mm_mul_sd(fgb,dz); 
			
			fix  = _mm_add_sd(fix,t1);
			fiy  = _mm_add_sd(fiy,t2);
			fiz  = _mm_add_sd(fiz,t3);	
			
			/* accumulate forces from memory */
			xmm5 = _mm_load_sd(f+aj13);
			xmm6 = _mm_load_sd(f+aj13+1);
			xmm7 = _mm_load_sd(f+aj13+2); 
					
			/* subtract partial forces */
			xmm5 = _mm_sub_sd(xmm5, t1); 
			xmm6 = _mm_sub_sd(xmm6, t2); 
			xmm7 = _mm_sub_sd(xmm7, t3); 
					
			/* store the force */
			_mm_store_sd(f+aj13,xmm5);
			_mm_store_sd(f+aj13+1,xmm6);
			_mm_store_sd(f+aj13+2,xmm7);
		}
		
		/* fix/fiy/fiz now contain two partial terms, that all should be
		 * added to the i particle forces
		 */
		t1		 = _mm_unpacklo_pd(t1,fix);
		t2		 = _mm_unpacklo_pd(t2,fiy);
		t3		 = _mm_unpacklo_pd(t3,fiz);
		
		fix		 = _mm_add_pd(fix,t1);
		fiy		 = _mm_add_pd(fiy,t2);
		fiz		 = _mm_add_pd(fiz,t3);
		
		fix      = _mm_shuffle_pd(fix,fix,_MM_SHUFFLE2(1,1));
		fiy      = _mm_shuffle_pd(fiy,fiy,_MM_SHUFFLE2(1,1));
		fiz      = _mm_shuffle_pd(fiz,fiz,_MM_SHUFFLE2(1,1));
		
		/* Load i forces from memory */
		xmm1     = _mm_load_sd(f+ai3);
		xmm2     = _mm_load_sd(f+ai3+1);
		xmm3     = _mm_load_sd(f+ai3+2);
		
		/* Add to i force */
		fix      = _mm_add_sd(fix,xmm1);
		fiy      = _mm_add_sd(fiy,xmm2);
		fiz      = _mm_add_sd(fiz,xmm3);
		
		/* store i forces to memory */
		_mm_store_sd(f+ai3,fix);
		_mm_store_sd(f+ai3+1,fiy);
		_mm_store_sd(f+ai3+2,fiz);

	}
	
	return 0;
}

#else
/* keep compiler happy */
int genborn_sse_dummy;

#endif /* SSE2 intrinsics available */

