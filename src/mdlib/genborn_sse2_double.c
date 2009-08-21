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
#ifdef GMX_THREADS
#include "tmpi.h"
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

static inline void
sincos_sse2double(__m128d x, __m128d *sinval, __m128d *cosval)
{
    const __m128d two_over_pi = {2.0/M_PI,2.0/M_PI};
    const __m128d half        = {0.5,0.5};
    const __m128d one         = {1.0,1.0};
    const __m128i izero       = _mm_set1_epi32(0);
    const __m128i ione        = _mm_set1_epi32(1);
    const __m128i itwo        = _mm_set1_epi32(2);
    const __m128i ithree      = _mm_set1_epi32(3);
    const __m128d sincosd_kc1 = {(13176794.0 / 8388608.0),(13176794.0 / 8388608.0)};
    const __m128d sincosd_kc2 = {7.5497899548918821691639751442098584e-8,7.5497899548918821691639751442098584e-8};
    const __m128d sincosd_cc0 = {0.00000000206374484196,0.00000000206374484196};
    const __m128d sincosd_cc1 = {-0.00000027555365134677,-0.00000027555365134677};
    const __m128d sincosd_cc2 = {0.00002480157946764225,0.00002480157946764225};
    const __m128d sincosd_cc3 = {-0.00138888888730525966,-0.00138888888730525966};
    const __m128d sincosd_cc4 = {0.04166666666651986722,0.04166666666651986722};
    const __m128d sincosd_cc5 = {-0.49999999999999547304,-0.49999999999999547304};
    const __m128d sincosd_sc0 = {0.00000000015893606014,0.00000000015893606014};
    const __m128d sincosd_sc1 = {-0.00000002505069049138,-0.00000002505069049138};
    const __m128d sincosd_sc2 = {0.00000275573131527032,0.00000275573131527032};
    const __m128d sincosd_sc3 = {-0.00019841269827816117,-0.00019841269827816117};
    const __m128d sincosd_sc4 = {0.00833333333331908278,0.00833333333331908278};
    const __m128d sincosd_sc5 = {-0.16666666666666612594,-0.16666666666666612594};
    
    __m128d signbit           = (__m128d) _mm_set1_epi64x(0x8000000000000000ULL);
    __m128d tiny              = (__m128d) _mm_set1_epi64x(0x3e40000000000000ULL);
    
    __m128d xl,xl2,xl3,qd,absxl,p1,cx,sx,ts,tc,tsn,tcn;
    __m128i q;
    __m128i offsetSin,offsetCos;
    __m128d sinMask,cosMask,isTiny;
    __m128d ct0,ct1,ct2,ct3,ct4,ct5,ct6,st1,st2,st3,st4,st6;
    
    /* Rescale the angle to the range 0..4, and find which quadrant it is in */
    xl        = _mm_mul_pd(x,two_over_pi);
    
    /* q=integer part of xl, rounded _away_ from 0.0 */
    /* Add 0.5 away from 0.0 */
    xl        = _mm_add_pd(xl,_mm_or_pd(_mm_and_pd(xl,signbit),half));
	
    q         = _mm_cvttpd_epi32(xl);
    qd        = _mm_cvtepi32_pd(q);
    q         = _mm_shuffle_epi32(q,_MM_SHUFFLE(1,1,0,0));
	
    /* Compute offset based on quadrant the arg falls in */
    offsetSin   = _mm_and_si128(q,ithree);
    offsetCos   = _mm_add_epi32(offsetSin,ione);
    
    /* Remainder in range [-pi/4..pi/4] */
    p1 = _mm_mul_pd(qd,sincosd_kc1);
    xl = _mm_mul_pd(qd,sincosd_kc2);
    p1 = _mm_sub_pd(x,p1);    
    xl = _mm_sub_pd(p1,xl);
    
    absxl  = _mm_andnot_pd(signbit,xl);
    isTiny = _mm_cmpgt_pd(tiny,absxl);
    
    xl2    = _mm_mul_pd(xl,xl);
    xl3    = _mm_mul_pd(xl2,xl);
    
    ct0    = _mm_mul_pd(xl2,xl2);
    ct1    = _mm_mul_pd(sincosd_cc0,xl2);
    ct2    = _mm_mul_pd(sincosd_cc2,xl2);
    ct3    = _mm_mul_pd(sincosd_cc4,xl2);
    st1    = _mm_mul_pd(sincosd_sc0,xl2);
    st2    = _mm_mul_pd(sincosd_sc2,xl2);
    st3    = _mm_mul_pd(sincosd_sc4,xl2);
    ct1    = _mm_add_pd(ct1,sincosd_cc1);
    ct2    = _mm_add_pd(ct2,sincosd_cc3);
    ct3    = _mm_add_pd(ct3,sincosd_cc5);
    st1    = _mm_add_pd(st1,sincosd_sc1);
    st2    = _mm_add_pd(st2,sincosd_sc3);
    st3    = _mm_add_pd(st3,sincosd_sc5);
	
    ct4    = _mm_mul_pd(ct2,ct0);
    ct4    = _mm_add_pd(ct4,ct3);
    
    st4    = _mm_mul_pd(st2,ct0);
    st4    = _mm_add_pd(st4,st3);
    ct5    = _mm_mul_pd(ct0,ct0);
    
    ct6    = _mm_mul_pd(ct5,ct1);
    ct6    = _mm_add_pd(ct6,ct4);
    
    st6    = _mm_mul_pd(ct5,st1);
    st6    = _mm_add_pd(st6,st4);
    
    cx     = _mm_mul_pd(ct6,xl2);
    cx     = _mm_add_pd(cx,one);
    
    sx     = _mm_mul_pd(st6,xl3);
    sx     = _mm_add_pd(sx,xl);
    
    /* Small angle approximation, sin(tiny)=tiny, cos(tiny)=1.0 */
    sx     = _mm_or_pd( _mm_and_pd(isTiny,xl) , _mm_andnot_pd(isTiny,sx) );
    cx     = _mm_or_pd( _mm_and_pd(isTiny,one) , _mm_andnot_pd(isTiny,cx) );
	
    sinMask = (__m128d) _mm_cmpeq_epi32( _mm_and_si128(offsetSin,ione), izero);
    cosMask = (__m128d) _mm_cmpeq_epi32( _mm_and_si128(offsetCos,ione), izero);
    
    ts     = _mm_or_pd( _mm_and_pd(sinMask,sx) , _mm_andnot_pd(sinMask,cx) );
    tc     = _mm_or_pd( _mm_and_pd(cosMask,sx) , _mm_andnot_pd(cosMask,cx) );
	
    /* Flip the sign of the result when (offset mod 4) = 1 or 2 */
    sinMask = (__m128d) _mm_cmpeq_epi32( _mm_and_si128(offsetSin,itwo), izero);
    tsn    = _mm_xor_pd(signbit,ts);
    ts     = _mm_or_pd( _mm_and_pd(sinMask,ts) , _mm_andnot_pd(sinMask,tsn) );
	
    cosMask = (__m128d) _mm_cmpeq_epi32( _mm_and_si128(offsetCos,itwo), izero);
    tcn    = _mm_xor_pd(signbit,tc);
    tc     = _mm_or_pd( _mm_and_pd(cosMask,tc) , _mm_andnot_pd(cosMask,tcn) );
	
    *sinval = ts;
    *cosval = tc;
    
    return;
}

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


int
calc_gb_rad_still_sse2_double(t_commrec *cr, t_forcerec *fr,int natoms, gmx_localtop_t *top,
							  const t_atomtypes *atype, double *x, t_nblist *nl, gmx_genborn_t *born, t_mdatoms *md)
{
	int i,k,n,ai,ai3,aj1,aj2,aj13,aj23;
	int at0,at1,nj0,nj1,offset,taj1,taj2;

	double factor,gpi_ai,gpi_tmp,gpi2;
	
	__m128d ix,iy,iz,jx,jy,jz,dx,dy,dz;
	__m128d t1,t2,t3,rsq11,rinv,rinv2,rinv4,rinv6;
	__m128d ratio,gpi,rai,raj,vaj,rvdw,mask_cmp;
	__m128d ccf,dccf,theta,cosq,term,sinq,res,prod;
	__m128d xmm1,xmm2,xmm3,xmm4,xmm7,vai,prod_ai,icf4,icf6;
	
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
	n                   = 0;
		
	/* Keep the compiler happy */
	raj         = _mm_setzero_pd();
	vaj         = _mm_setzero_pd();
	jx          = _mm_setzero_pd();
	jy          = _mm_setzero_pd();
	jz          = _mm_setzero_pd();
	xmm2        = _mm_setzero_pd();
	xmm7        = _mm_setzero_pd();
		
	for(i=0;i<born->nr;i++)
	{
		born->gpol_still_work[i]=0;
	}
	
	for(i=0;i<nl->nri;i++)
	{
		ai      = nl->iinr[i];
		ai3     = ai * 3;
		nj0     = nl->jindex[ai];
		nj1     = nl->jindex[ai+1];
		
		offset  = (nj1-nj0)%2;
		
		/* Polarization energy for atom ai */
		gpi     = _mm_setzero_pd();
	
		/* Load particle ai coordinates */
		ix      = _mm_load1_pd(x+ai3);
		iy      = _mm_load1_pd(x+ai3+1);
		iz      = _mm_load1_pd(x+ai3+2);
		
		/* Load particle ai gb_radius */
		rai     = _mm_set1_pd(top->atomtypes.gb_radius[md->typeA[ai]]);
		vai     = _mm_set1_pd(born->vsolv[ai]);
		prod_ai = _mm_mul_pd(p4,vai);
				
		for(k=nj0;k<nj1-offset;k+=2)
		{
			aj1      = nl->jjnr[k];
			aj2      = nl->jjnr[k+1];
			
			aj13     = aj1 * 3;
			aj23     = aj2 * 3;
			
			taj1     = md->typeA[aj1];
			taj2     = md->typeA[aj2];
			
			/* Load particle aj1-2 coordinates and compute ai->aj distances */
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
			
			raj      = _mm_loadl_pd(raj, top->atomtypes.gb_radius+taj1);
			raj      = _mm_loadh_pd(raj, top->atomtypes.gb_radius+taj2);
			
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
					sincos_sse2double(theta,&sinq,&cosq);
					
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
			icf4      = _mm_mul_pd(ccf,rinv4);
			xmm2      = _mm_mul_pd(icf4,prod);
			xmm3      = _mm_mul_pd(icf4,prod_ai);
			gpi       = _mm_add_pd(gpi,xmm2);
			
			/* Load, subtract and store atom aj gpol energy */
			xmm7      = _mm_loadl_pd(xmm7,born->gpol_still_work+aj1);
			xmm7      = _mm_loadh_pd(xmm7,born->gpol_still_work+aj2);
			
			xmm3      = _mm_add_pd(xmm7,xmm3);
			
			_mm_storel_pd(born->gpol_still_work+aj1,xmm3);
			_mm_storeh_pd(born->gpol_still_work+aj2,xmm3);
			
			/* Chain rule terms */
			ccf       = _mm_mul_pd(four,ccf);
			xmm3      = _mm_sub_pd(ccf,dccf);
			icf6      = _mm_mul_pd(xmm3,rinv6);
			xmm1      = _mm_mul_pd(icf6,prod);
			xmm2      = _mm_mul_pd(icf6,prod_ai);
			
			/* As with single precision, we need to shift stuff around, to change the order of
			 * the interactions from ai->aj1, ai->aj2 to ai->aj1, aj1->ai, ai->aj2, aj2->ai etc,
			 to do 2 instead of 4 store operations
			 */
			xmm3      = _mm_unpacklo_pd(xmm1,xmm2);
			xmm4      = _mm_unpackhi_pd(xmm1,xmm2);
			
			_mm_storeu_pd(fr->dadx+n,xmm3);
			n         = n + 2;
			_mm_storeu_pd(fr->dadx+n,xmm4);
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
			
			raj   = _mm_load_sd(top->atomtypes.gb_radius+taj1);
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
					sincos_sse2double(theta,&sinq,&cosq);
										
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
			icf4      = _mm_mul_sd(ccf,rinv4);
			xmm2      = _mm_mul_sd(icf4,prod);
			xmm3      = _mm_mul_sd(icf4,prod_ai);
			gpi       = _mm_add_sd(gpi,xmm2);
			
			/* Load, subtract and store atom aj gpol energy */
			xmm7      = _mm_load_sd(born->gpol_still_work+aj1);
			xmm3      = _mm_add_sd(xmm7,xmm3);
			_mm_store_sd(born->gpol_still_work+aj1,xmm3);
			
			/* Chain rule terms */
			ccf       = _mm_mul_sd(four,ccf);
			xmm3      = _mm_sub_sd(ccf,dccf);
			icf6      = _mm_mul_sd(xmm3,rinv6);
			xmm1      = _mm_mul_sd(icf6,prod);
			xmm2      = _mm_mul_sd(icf6,prod_ai);
			
			/* Here we only have ai->aj1 and aj1->ai, so we can store directly */
			_mm_storel_pd(fr->dadx+n,xmm1);
			n         = n + 1;
			_mm_storel_pd(fr->dadx+n,xmm2);
			n         = n + 1;
		} /* End offset */

		/* Do end processing ... */
		xmm2  = _mm_unpacklo_pd(xmm2,gpi);
		gpi   = _mm_add_pd(gpi,xmm2);
		gpi   = _mm_shuffle_pd(gpi,gpi,_MM_SHUFFLE2(1,1));
			
		/* Load, add and store atom ai polarisation energy */
		xmm2 = _mm_load_sd(born->gpol_still_work+ai);
		gpi = _mm_add_sd(gpi,xmm2);
		_mm_store_sd(born->gpol_still_work+ai,gpi);
	}
	
	/* Parallell summations */
	if(PARTDECOMP(cr))
	{
		gmx_sum(natoms,born->gpol_still_work, cr);
	}
	else if(DOMAINDECOMP(cr))
	{
		dd_atom_sum_real(cr->dd, born->gpol_still_work);
	}
	
	/* Compute the radii */
	for(i=0;i<nl->nri;i++)
	{
		ai = nl->iinr[i];
		gpi_ai = born->gpol[ai] + born->gpol_still_work[ai];
		gpi2   = gpi_ai*gpi_ai;
		
		born->bRad[ai]=factor*invsqrt(gpi2);
		fr->invsqrta[ai]=invsqrt(born->bRad[ai]);
	}
	
	/* Extra (local) communication reqiured for DD */
	if(DOMAINDECOMP(cr))
	{
		dd_atom_spread_real(cr->dd, born->bRad);
		dd_atom_spread_real(cr->dd, fr->invsqrta);
	}
		
	return 0;
}

int 
calc_gb_rad_hct_sse2_double(t_commrec *cr, t_forcerec *fr, int natoms, gmx_localtop_t *top, 
							const t_atomtypes *atype, double *x, t_nblist *nl, gmx_genborn_t *born, t_mdatoms *md)
{
	int i,k,n,ai,ai3,aj1,aj2,aj13,aj23,nj0,nj1,at0,at1,offset;
	int p1, p2;
	double rr,sum,sum_tmp,min_rad,rad,doff;
	
	__m128d ix,iy,iz,jx,jy,jz,dx,dy,dz;
	__m128d t1,t2,t3,rsq11,rinv,r,rai;
	__m128d rai_inv,sk,sk2,lij,dlij,duij;
	__m128d uij,lij2,uij2,lij3,uij3,diff2;
	__m128d lij_inv,sk2_inv,prod,log_term,tmp,tmp_sum;
	__m128d mask_cmp,mask_cmp2,mask_cmp3;
	__m128d xmm1,xmm2,xmm3,xmm4,xmm7,xmm8,xmm9,doffset;
	__m128d sum_ai,chrule,chrule_ai,tmp_ai; 
	__m128d sk_ai, sk2_ai,raj,raj_inv;
	
	const __m128d neg   = {-1.0, -1.0};
	const __m128d zero  = {0.0, 0.0};
	const __m128d eigth = {0.125, 0.125};
	const __m128d qrtr  = {0.25, 0.25};
	const __m128d half  = {0.5, 0.5};
	const __m128d one   = {1.0, 1.0};
	const __m128d two   = {2.0, 2.0};
	const __m128d three = {3.0, 3.0};
	
	/* Keep the compiler happy */
	tmp_ai      = _mm_setzero_pd();
	tmp         = _mm_setzero_pd();
	xmm7        = _mm_setzero_pd();
	xmm8        = _mm_setzero_pd();
	raj_inv     = _mm_setzero_pd();
	raj         = _mm_setzero_pd();
	sk          = _mm_setzero_pd();
	jx          = _mm_setzero_pd();
	jy          = _mm_setzero_pd();
	jz          = _mm_setzero_pd();
	
	/* Set the dielectric offset */
	doff = born->gb_doffset;
	doffset = _mm_load1_pd(&doff);
	n       = 0;
	
	for(i=0;i<born->nr;i++)
	{
		born->gpol_hct_work[i] = 0;
	}
	
	for(i=0;i<nl->nri;i++)
	{
		ai      = nl->iinr[i];
		ai3     = ai * 3;
		
		nj0     = nl->jindex[ai];
		nj1     = nl->jindex[ai+1];
		
		offset  = (nj1-nj0)%2;
		
		/* Load rai */
		rr      = top->atomtypes.gb_radius[md->typeA[ai]]-doff;
		rai     = _mm_load1_pd(&rr);
		rr      = 1.0/rr;
		rai_inv = _mm_load1_pd(&rr);
		
		/* Zero out sums for polarisation energies */
		sum_ai  = _mm_setzero_pd();
		
		/* Load ai coordinates */
		ix       = _mm_load1_pd(x+ai3);
		iy       = _mm_load1_pd(x+ai3+1);
		iz       = _mm_load1_pd(x+ai3+2);
		
		sk_ai    = _mm_load1_pd(born->param+ai);
		sk2_ai   = _mm_mul_pd(sk_ai,sk_ai);
	
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
			
			/* Load aj1,aj2 raj */
			p1        = md->typeA[aj1];
			p2        = md->typeA[aj2];
			
			raj       = _mm_loadl_pd(raj,top->atomtypes.gb_radius+p1);
			raj       = _mm_loadh_pd(raj,top->atomtypes.gb_radius+p2);
			raj       = _mm_sub_pd(raj,doffset);
			
			/* Compute 1.0/raj */
			raj_inv   = my_inv_pd(raj);
			
			/* INTERACTION aj->ai STARS HERE */
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
			log_term  = log_pd(log_term);
			
			xmm1      = _mm_sub_pd(lij,uij);
			xmm2      = _mm_mul_pd(qrtr,r);
			xmm2      = _mm_mul_pd(xmm2,diff2);
			xmm1      = _mm_add_pd(xmm1,xmm2); 
			xmm2      = _mm_mul_pd(half,rinv); 
			xmm2      = _mm_mul_pd(xmm2,log_term); 
			xmm1      = _mm_add_pd(xmm1,xmm2); 
			xmm9      = _mm_mul_pd(neg,diff2); 
			xmm2      = _mm_mul_pd(xmm9,prod); 
			tmp_ai    = _mm_add_pd(xmm1,xmm2); 
			
			/* contitional for rai<sk-dr */
			xmm3      = _mm_sub_pd(sk,r);
			mask_cmp3 = _mm_cmplt_pd(rai,xmm3); /* rai<sk-dr */
			
			xmm4    = _mm_sub_pd(rai_inv,lij);
			xmm4    = _mm_mul_pd(two,xmm4);
			xmm4    = _mm_add_pd(tmp_ai,xmm4);
			
			tmp_ai	    = _mm_or_pd(_mm_and_pd(mask_cmp3,xmm4)  ,_mm_andnot_pd(mask_cmp3,tmp_ai)); /*conditional as a mask*/
			
			/* the tmp will now contain two partial values, that not all are to be used. Which */
			/* ones are governed by the mask_cmp mask. */
			tmp_ai     = _mm_mul_pd(half,tmp_ai); 
			tmp_ai     = _mm_or_pd(_mm_and_pd(mask_cmp,tmp_ai)  ,_mm_andnot_pd(mask_cmp,zero)); /*conditional as a mask*/
			sum_ai     = _mm_add_pd(sum_ai,tmp_ai);
			
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
			
			/* temporary storage of chain rule terms, since we have to compute
			 the reciprocal terms also before storing them */
			chrule = _mm_mul_pd(xmm2,rinv); 	
			
			/* INTERACTION ai->aj starts here */
			/* Conditional mask for raj<dr+sk_ai */
			xmm1   = _mm_add_pd(r,sk_ai);
			mask_cmp = _mm_cmplt_pd(raj,xmm1);
			
			/* Conditional for rai>dr-sk, ends with mask_cmp2 */
			xmm2  = _mm_sub_pd(r,sk_ai);
			xmm3  = my_inv_pd(xmm2);
			mask_cmp2 = _mm_cmpgt_pd(raj,xmm2);
			
			lij	      = _mm_or_pd(_mm_and_pd(mask_cmp2,raj_inv)  ,_mm_andnot_pd(mask_cmp2,xmm3)); /*conditional as a mask*/
			dlij      = _mm_or_pd(_mm_and_pd(mask_cmp2,zero) ,_mm_andnot_pd(mask_cmp2,one));
			
			uij       = my_inv_pd(xmm1);
			lij2      = _mm_mul_pd(lij,lij);
			lij3      = _mm_mul_pd(lij2,lij);
			uij2      = _mm_mul_pd(uij,uij);
			uij3      = _mm_mul_pd(uij2,uij);
			
			diff2     = _mm_sub_pd(uij2,lij2);
			
			lij_inv   = my_invrsq_pd(lij2);
			
			sk2       = sk2_ai;
			sk2_inv   = _mm_mul_pd(sk2,rinv);
			prod      = _mm_mul_pd(qrtr,sk2_inv);
			
			log_term  = _mm_mul_pd(uij,lij_inv);
			log_term = log_pd(log_term);
			
			xmm1      = _mm_sub_pd(lij,uij);
			xmm2      = _mm_mul_pd(qrtr,r);
			xmm2      = _mm_mul_pd(xmm2,diff2);
			xmm1      = _mm_add_pd(xmm1,xmm2); 
			xmm2      = _mm_mul_pd(half,rinv); 
			xmm2      = _mm_mul_pd(xmm2,log_term); 
			xmm1      = _mm_add_pd(xmm1,xmm2); 
			xmm9      = _mm_mul_pd(neg,diff2); 
			xmm2      = _mm_mul_pd(xmm9,prod); 
			tmp       = _mm_add_pd(xmm1,xmm2); 
			
			/* contitional for rai<sk_ai-dr */
			xmm3      = _mm_sub_pd(sk_ai,r);
			mask_cmp3 = _mm_cmplt_pd(raj,xmm3); /* rai<sk-dr */
			
			xmm4    = _mm_sub_pd(raj_inv,lij);
			xmm4    = _mm_mul_pd(two,xmm4);
			xmm4    = _mm_add_pd(tmp,xmm4);
			
			tmp	    = _mm_or_pd(_mm_and_pd(mask_cmp3,xmm4)  ,_mm_andnot_pd(mask_cmp3,tmp)); /*conditional as a mask*/
			
			/* the tmp will now contain two partial values, that not all are to be used. Which */
			/* ones are governed by the mask_cmp mask. */
			tmp     = _mm_mul_pd(half,tmp); 
			tmp     = _mm_or_pd(_mm_and_pd(mask_cmp,tmp)  ,_mm_andnot_pd(mask_cmp,zero)); /*conditional as a mask*/
			
			/* Load, add and store ai->aj pol energy */
			xmm7    = _mm_loadl_pd(xmm7,born->gpol_hct_work+aj1);
			xmm7    = _mm_loadh_pd(xmm7,born->gpol_hct_work+aj2);
			
			xmm7    = _mm_add_pd(xmm7,tmp);
			
			_mm_storel_pd(born->gpol_hct_work+aj1,xmm7);
			_mm_storeh_pd(born->gpol_hct_work+aj2,xmm7);
			
			/* Start chain rule terms */
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
			chrule_ai = _mm_mul_pd(xmm2,rinv);
			
			/* Store chain rule terms 
			 * same unpacking rationale as with Still above 
			 */
			xmm3 = _mm_unpacklo_pd(chrule, chrule_ai);
			xmm4 = _mm_unpackhi_pd(chrule, chrule_ai);
			
			_mm_storeu_pd(fr->dadx+n, xmm3);
			n = n + 2;
			_mm_storeu_pd(fr->dadx+n, xmm4);
			n = n + 2;
		}
		
		if(offset!=0)
		{
			aj1       = nl->jjnr[k];
			aj13      = aj1 * 3;
			p1        = md->typeA[aj1];
			
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
			
			/* Load raj */
			raj       = _mm_load_sd(top->atomtypes.gb_radius+p1);
			raj       = _mm_sub_sd(raj,doffset);
			raj_inv   = my_inv_pd(raj);
			
			/* OFFSET INTERATIONS aj->ai STARTS HERE */
			/* conditional mask for rai<dr+sk */
			xmm1      = _mm_add_sd(r,sk);
			mask_cmp  = _mm_cmplt_sd(rai,xmm1);
			
			/* conditional for rai>dr-sk, ends with mask_cmp2 */
			xmm2      = _mm_sub_sd(r,sk);
			xmm3      = my_inv_pd(xmm2);
			mask_cmp2 = _mm_cmpgt_pd(rai,xmm2);
			
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
			log_term = log_pd(log_term);
			
			xmm1      = _mm_sub_sd(lij,uij);
			xmm2      = _mm_mul_sd(qrtr,r);
			xmm2      = _mm_mul_sd(xmm2,diff2);
			xmm1      = _mm_add_sd(xmm1,xmm2); 
			xmm2      = _mm_mul_sd(half,rinv); 
			xmm2      = _mm_mul_sd(xmm2,log_term); 
			xmm1      = _mm_add_sd(xmm1,xmm2); 
			xmm9      = _mm_mul_sd(neg,diff2); 
			xmm2      = _mm_mul_sd(xmm9,prod); 
			
			tmp_ai    = _mm_add_sd(xmm1,xmm2); 
			
			/* contitional for rai<sk-dr */
			xmm3      = _mm_sub_sd(sk,r);
			mask_cmp3 = _mm_cmplt_sd(rai,xmm3); /* rai<sk-dr */
			
			xmm4    = _mm_sub_sd(rai_inv,lij);
			xmm4    = _mm_mul_sd(two,xmm4);
			xmm4    = _mm_add_sd(tmp_ai,xmm4);
			
			tmp_ai	    = _mm_or_pd(_mm_and_pd(mask_cmp3,xmm4)  ,_mm_andnot_pd(mask_cmp3,tmp_ai)); /*conditional as a mask*/
			
			/* the tmp will now contain two partial values, that not all are to be used. Which */
			/* ones are governed by the mask_cmp mask. */
			tmp_ai     = _mm_mul_pd(half,tmp_ai); 
			tmp_ai     = _mm_or_pd(_mm_and_pd(mask_cmp,tmp_ai)  ,_mm_andnot_pd(mask_cmp,zero)); /*conditional as a mask*/
			sum_ai     = _mm_add_sd(sum_ai,tmp_ai);
					
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
			
			/* temporary storage of chain rule terms, since we have to compute
			 the reciprocal terms also before storing them */
			chrule = _mm_mul_sd(xmm2,rinv);
			
			/* OFFSET INTERACTION ai->aj starts here */
			/* conditional mask for raj<dr+sk */
			xmm1      = _mm_add_sd(r,sk_ai);
			mask_cmp  = _mm_cmplt_sd(raj,xmm1);
			
			/* conditional for rai>dr-sk, ends with mask_cmp2 */
			xmm2      = _mm_sub_sd(r,sk_ai);
			xmm3      = my_inv_pd(xmm2);
			mask_cmp2 = _mm_cmpgt_pd(raj,xmm2);
			
			lij	      = _mm_or_pd(_mm_and_pd(mask_cmp2,raj_inv)  ,_mm_andnot_pd(mask_cmp2,xmm3)); /*conditional as a mask*/
			dlij      = _mm_or_pd(_mm_and_pd(mask_cmp2,zero) ,_mm_andnot_pd(mask_cmp2,one));
			
			uij       = my_inv_pd(xmm1);
			lij2      = _mm_mul_sd(lij,lij);
			lij3      = _mm_mul_sd(lij2,lij);
			uij2      = _mm_mul_sd(uij,uij);
			uij3      = _mm_mul_sd(uij2,uij);
			
			diff2     = _mm_sub_sd(uij2,lij2);
			
			lij_inv   = my_invrsq_pd(lij2);
			
			sk2       = sk2_ai;
			sk2_inv   = _mm_mul_sd(sk2,rinv);
			prod      = _mm_mul_sd(qrtr,sk2_inv);
			
			log_term  = _mm_mul_pd(uij,lij_inv);
			log_term = log_pd(log_term);
			
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
			xmm3      = _mm_sub_sd(sk_ai,r);
			mask_cmp3 = _mm_cmplt_sd(raj,xmm3); /* rai<sk-dr */
			
			xmm4    = _mm_sub_sd(raj_inv,lij);
			xmm4    = _mm_mul_sd(two,xmm4);
			xmm4    = _mm_add_sd(tmp,xmm4);
			
			tmp	    = _mm_or_pd(_mm_and_pd(mask_cmp3,xmm4)  ,_mm_andnot_pd(mask_cmp3,tmp)); /*conditional as a mask*/
			
			/* the tmp will now contain two partial values, that not all are to be used. Which */
			/* ones are governed by the mask_cmp mask. */
			tmp     = _mm_mul_pd(half,tmp); 
			tmp     = _mm_or_pd(_mm_and_pd(mask_cmp,tmp)  ,_mm_andnot_pd(mask_cmp,zero)); /*conditional as a mask*/
			
			/* Load, add and store gpol energy */
			xmm7    = _mm_load_sd(born->gpol_hct_work+aj1);
			xmm7    = _mm_add_sd(xmm7,tmp);
			_mm_store_sd(born->gpol_hct_work+aj1,xmm7);
		
			/* Start chain rule terms, t1 */
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
			chrule_ai = _mm_mul_sd(xmm2,rinv);
			
			_mm_store_sd(fr->dadx+n, chrule); 
			n = n + 1;
			_mm_store_sd(fr->dadx+n, chrule_ai);
			n = n + 1;
		}
		
		/* Do end processing ...  */
		tmp_ai = _mm_unpacklo_pd(tmp_ai,sum_ai);
		sum_ai = _mm_add_pd(sum_ai,tmp_ai);
		sum_ai = _mm_shuffle_pd(sum_ai,sum_ai,_MM_SHUFFLE2(1,1));
				
		/* Load, add and store atom ai polarisation energy */
		xmm2 = _mm_load_sd(born->gpol_hct_work+ai);
		sum_ai = _mm_add_sd(sum_ai,xmm2);
		_mm_store_sd(born->gpol_hct_work+ai,sum_ai);
	}
	
	/* Parallel summations */
	if(PARTDECOMP(cr))
	{
		gmx_sum(natoms, born->gpol_hct_work, cr);
	}
	else if(DOMAINDECOMP(cr))
	{
		dd_atom_sum_real(cr->dd, born->gpol_hct_work);
	}
	
	/* Compute the radii */
	for(i=0;i<nl->nri;i++)
	{
		ai      = nl->iinr[i];
		rr      = top->atomtypes.gb_radius[md->typeA[ai]]-doff; 
		sum     = 1.0/rr - born->gpol_hct_work[ai];
		min_rad = rr + doff;
		rad     = 1.0/sum;  
		
		born->bRad[ai]   = rad > min_rad ? rad : min_rad;
		fr->invsqrta[ai] = invsqrt(born->bRad[ai]);
	}
	
	/* Extra (local) communication required for DD */
	if(DOMAINDECOMP(cr))
	{
		dd_atom_spread_real(cr->dd, born->bRad);
		dd_atom_spread_real(cr->dd, fr->invsqrta);
	}
	
	return 0;
}

int 
calc_gb_rad_obc_sse2_double(t_commrec *cr, t_forcerec * fr, int natoms, gmx_localtop_t *top,
							const t_atomtypes *atype, double *x, t_nblist *nl, gmx_genborn_t *born,t_mdatoms *md)
{
	int i,k,n,ai,ai3,aj1,aj2,aj13,aj23,nj0,nj1,at0,at1,offset;
	int p1,p2,p3,p4;
	double rr,sum,sum_tmp,sum2,sum3,min_rad,rad,doff;
	double tsum,tchain,rr_inv,rr_inv2,gbr;
	
	__m128d ix,iy,iz,jx,jy,jz,dx,dy,dz;
	__m128d t1,t2,t3,rsq11,rinv,r,rai;
	__m128d rai_inv,sk,sk2,lij,dlij,duij;
	__m128d uij,lij2,uij2,lij3,uij3,diff2;
	__m128d lij_inv,sk2_inv,prod,log_term,tmp,tmp_sum;
	__m128d mask_cmp,mask_cmp2,mask_cmp3,doffset,raj,raj_inv;
	__m128d xmm1,xmm2,xmm3,xmm4,xmm7,xmm8,xmm9;
	__m128d sum_ai,chrule,chrule_ai,tmp_ai,sk_ai,sk2_ai;
	
	const __m128d neg   = {-1.0, -1.0};
	const __m128d zero  = {0.0, 0.0};
	const __m128d eigth = {0.125, 0.125};
	const __m128d qrtr  = {0.25, 0.25};
	const __m128d half  = {0.5, 0.5};
	const __m128d one   = {1.0, 1.0};
	const __m128d two   = {2.0, 2.0};
	const __m128d three = {3.0, 3.0};
	
	/* Keep the compiler happy */
	tmp_ai      = _mm_setzero_pd();
	tmp         = _mm_setzero_pd();
	raj_inv     = _mm_setzero_pd();
	raj         = _mm_setzero_pd();
	sk          = _mm_setzero_pd();
	jx          = _mm_setzero_pd();
	jy          = _mm_setzero_pd();
	jz          = _mm_setzero_pd();
	xmm7        = _mm_setzero_pd();
	xmm8        = _mm_setzero_pd();
	
	/* Set the dielectric offset */
	doff = born->gb_doffset;
	doffset = _mm_load1_pd(&doff);
	
	n       = 0;
	
	for(i=0;i<born->nr;i++)
	{
		born->gpol_hct_work[i] = 0;
	}

	for(i=0;i<nl->nri;i++)
	{
		ai      = nl->iinr[i];
		ai3     = ai * 3;
		
		nj0     = nl->jindex[ai];
		nj1     = nl->jindex[ai+1];
		
		offset  = (nj1-nj0)%2;
		
		/* Load rai */
		rr      = top->atomtypes.gb_radius[md->typeA[ai]]-doff;
		rai     = _mm_load1_pd(&rr);
		rr      = 1.0/rr;
		rai_inv = _mm_load1_pd(&rr);
		
		/* Load ai coordinates */
		ix       = _mm_load1_pd(x+ai3);
		iy       = _mm_load1_pd(x+ai3+1);
		iz       = _mm_load1_pd(x+ai3+2);
		
		/* Zero out sums for polarisation energies */
		sum_ai = _mm_setzero_pd();
		
		sk_ai  =  _mm_load1_pd(born->param+ai);
		sk2_ai = _mm_mul_pd(sk_ai,sk_ai);
				
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
			
			/* Load atom aj1,aj2 raj */
			p1         = md->typeA[aj1];
			p2         = md->typeA[aj2];
			
			raj        = _mm_loadl_pd(raj,top->atomtypes.gb_radius+p1);
			raj        = _mm_loadh_pd(raj,top->atomtypes.gb_radius+p2);
			raj        = _mm_sub_pd(raj,doffset);
			
			/* Compute 1.0/raj */
			raj_inv    = my_inv_pd(raj);
			
			sk         = _mm_loadl_pd(sk,born->param+aj1);
			sk         = _mm_loadh_pd(sk,born->param+aj2);
			
			/* INTERACTION aj->ai STARTS HERE */
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
			log_term = log_pd(log_term);
			
			xmm1      = _mm_sub_pd(lij,uij);
			xmm2      = _mm_mul_pd(qrtr,r);
			xmm2      = _mm_mul_pd(xmm2,diff2);
			xmm1      = _mm_add_pd(xmm1,xmm2); 
			xmm2      = _mm_mul_pd(half,rinv); 
			xmm2      = _mm_mul_pd(xmm2,log_term); 
			xmm1      = _mm_add_pd(xmm1,xmm2); 
			xmm9      = _mm_mul_pd(neg,diff2); 
			xmm2      = _mm_mul_pd(xmm9,prod); 
			tmp_ai    = _mm_add_pd(xmm1,xmm2); 
			
			/* contitional for rai<sk-dr */
			xmm3      = _mm_sub_pd(sk,r);
			mask_cmp3 = _mm_cmplt_pd(rai,xmm3); /* rai<sk-dr */
			
			xmm4    = _mm_sub_pd(rai_inv,lij);
			xmm4    = _mm_mul_pd(two,xmm4);
			xmm4    = _mm_add_pd(tmp_ai,xmm4);
			
			tmp_ai	    = _mm_or_pd(_mm_and_pd(mask_cmp3,xmm4)  ,_mm_andnot_pd(mask_cmp3,tmp_ai)); /*conditional as a mask*/
			
			/* the tmp will now contain two partial values, that not all are to be used. Which */
			/* ones are governed by the mask_cmp mask. */
			tmp_ai     = _mm_mul_pd(half,tmp_ai); 
			tmp_ai     = _mm_or_pd(_mm_and_pd(mask_cmp,tmp_ai)  ,_mm_andnot_pd(mask_cmp,zero)); /*conditional as a mask*/
			sum_ai     = _mm_add_pd(sum_ai,tmp_ai);
			
			/* Start the dadx chain rule terms */
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
			
			/* temporary storage of chain rule terms, since we have to compute
			 the reciprocal terms also before storing them */
			chrule   = _mm_mul_pd(xmm2,rinv);
			
			/* INTERACTION ai->aj STARTS HERE */
			/* conditional mask for raj<dr+sk_ai */
			xmm1       = _mm_add_pd(r,sk_ai);
			mask_cmp   = _mm_cmplt_pd(raj,xmm1);
			
			/* conditional for rai>dr-sk, ends with mask_cmp2 */
			xmm2      = _mm_sub_pd(r,sk_ai);
			xmm3      = my_inv_pd(xmm2);
			mask_cmp2 = _mm_cmpgt_pd(raj,xmm2);
			
			lij	      = _mm_or_pd(_mm_and_pd(mask_cmp2,raj_inv)  ,_mm_andnot_pd(mask_cmp2,xmm3)); /*conditional as a mask*/
			dlij      = _mm_or_pd(_mm_and_pd(mask_cmp2,zero) ,_mm_andnot_pd(mask_cmp2,one));
			
			uij       = my_inv_pd(xmm1);
			lij2      = _mm_mul_pd(lij,lij);
			lij3      = _mm_mul_pd(lij2,lij);
			uij2      = _mm_mul_pd(uij,uij);
			uij3      = _mm_mul_pd(uij2,uij);
			
			diff2     = _mm_sub_pd(uij2,lij2);
			
			lij_inv   = my_invrsq_pd(lij2);
			
			sk2       = sk2_ai;
			sk2_inv   = _mm_mul_pd(sk2,rinv);
			prod      = _mm_mul_pd(qrtr,sk2_inv);
			
			log_term  = _mm_mul_pd(uij,lij_inv);
			log_term = log_pd(log_term);
			
			xmm1      = _mm_sub_pd(lij,uij);
			xmm2      = _mm_mul_pd(qrtr,r);
			xmm2      = _mm_mul_pd(xmm2,diff2);
			xmm1      = _mm_add_pd(xmm1,xmm2); 
			xmm2      = _mm_mul_pd(half,rinv); 
			xmm2      = _mm_mul_pd(xmm2,log_term); 
			xmm1      = _mm_add_pd(xmm1,xmm2); 
			xmm9      = _mm_mul_pd(neg,diff2); 
			xmm2      = _mm_mul_pd(xmm9,prod); 
			tmp       = _mm_add_pd(xmm1,xmm2); 
			
			/* contitional for raj<sk_ai-dr */
			xmm3      = _mm_sub_pd(sk_ai,r);
			mask_cmp3 = _mm_cmplt_pd(raj,xmm3); /* rai<sk-dr */
			
			xmm4    = _mm_sub_pd(raj_inv,lij);
			xmm4    = _mm_mul_pd(two,xmm4);
			xmm4    = _mm_add_pd(tmp,xmm4);
			
			tmp	    = _mm_or_pd(_mm_and_pd(mask_cmp3,xmm4)  ,_mm_andnot_pd(mask_cmp3,tmp)); /*conditional as a mask*/
			
			/* the tmp will now contain two partial values, that not all are to be used. Which */
			/* ones are governed by the mask_cmp mask. */
			tmp     = _mm_mul_pd(half,tmp); 
			tmp     = _mm_or_pd(_mm_and_pd(mask_cmp,tmp)  ,_mm_andnot_pd(mask_cmp,zero)); /*conditional as a mask*/
			
			/* Load, add and store ai->aj pol energy */
			xmm7    = _mm_loadl_pd(xmm7,born->gpol_hct_work+aj1);
			xmm7    = _mm_loadh_pd(xmm7,born->gpol_hct_work+aj2);
			
			xmm7    = _mm_add_pd(xmm7,tmp);
			
			_mm_storel_pd(born->gpol_hct_work+aj1,xmm7);
			_mm_storeh_pd(born->gpol_hct_work+aj2,xmm7);
			
			/* Start dadx chain rule terms */
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
			chrule_ai   = _mm_mul_pd(xmm2,rinv);
			
			/* Store chain rule terms 
			 * same unpacking rationale as with Still above 
			 */
			xmm3 = _mm_unpacklo_pd(chrule, chrule_ai);
			xmm4 = _mm_unpackhi_pd(chrule, chrule_ai);
			
			_mm_storeu_pd(fr->dadx+n, xmm3);
			n = n + 2;
			_mm_storeu_pd(fr->dadx+n, xmm4);
			n = n + 2;
		}
		
		if(offset!=0)
		{
			aj1       = nl->jjnr[k];
			aj13      = aj1 * 3;
			p1        = md->typeA[aj1];
			
			jx        = _mm_load_sd(x+aj13);
			jy        = _mm_load_sd(x+aj13+1);
			jz        = _mm_load_sd(x+aj13+2);
			
			sk        = _mm_load_sd(born->param+aj1);
			
			/* Load raj */
			raj        = _mm_load_sd(top->atomtypes.gb_radius+p1);
			raj        = _mm_sub_sd(raj,doffset);
			raj_inv    = my_inv_pd(raj);
			
			dx        = _mm_sub_sd(ix,jx);
			dy        = _mm_sub_pd(iy,jy);
			dz        = _mm_sub_pd(iz,jz);
			
			rsq11     = _mm_add_sd( _mm_add_sd( _mm_mul_sd(dx,dx) , _mm_mul_sd(dy,dy) ) , _mm_mul_sd(dz,dz) );
			rinv      = my_invrsq_pd(rsq11);
			r         = _mm_mul_sd(rinv,rsq11);
			
			/* OFFSET INTERACTION aj->ai STARTS HERE */
			/* conditional mask for rai<dr+sk */
			xmm1      = _mm_add_sd(r,sk);
			mask_cmp  = _mm_cmplt_sd(rai,xmm1);
			
			/* conditional for rai>dr-sk, ends with mask_cmp2 */
			xmm2      = _mm_sub_sd(r,sk);
			xmm3      = my_inv_pd(xmm2);
			mask_cmp2 = _mm_cmpgt_pd(rai,xmm2);
			
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
			log_term = log_pd(log_term);
			
			xmm1      = _mm_sub_sd(lij,uij);
			xmm2      = _mm_mul_sd(qrtr,r);
			xmm2      = _mm_mul_sd(xmm2,diff2);
			xmm1      = _mm_add_sd(xmm1,xmm2); 
			xmm2      = _mm_mul_sd(half,rinv); 
			xmm2      = _mm_mul_sd(xmm2,log_term); 
			xmm1      = _mm_add_sd(xmm1,xmm2); 
			xmm9      = _mm_mul_sd(neg,diff2); 
			xmm2      = _mm_mul_sd(xmm9,prod); 
			tmp_ai    = _mm_add_sd(xmm1,xmm2); 
			
			/* contitional for rai<sk-dr */
			xmm3      = _mm_sub_sd(sk,r);
			mask_cmp3 = _mm_cmplt_sd(rai,xmm3); /* rai<sk-dr */
			
			xmm4    = _mm_sub_sd(rai_inv,lij);
			xmm4    = _mm_mul_sd(two,xmm4);
			xmm4    = _mm_add_sd(tmp_ai,xmm4);
			
			tmp_ai  = _mm_or_pd(_mm_and_pd(mask_cmp3,xmm4)  ,_mm_andnot_pd(mask_cmp3,tmp_ai)); /*conditional as a mask*/
			
			/* the tmp will now contain two partial values, that not all are to be used. Which */
			/* ones are governed by the mask_cmp mask. */
			tmp_ai = _mm_mul_pd(half,tmp_ai); 
			tmp_ai = _mm_or_pd(_mm_and_pd(mask_cmp,tmp_ai)  ,_mm_andnot_pd(mask_cmp,zero)); /*conditional as a mask*/
			sum_ai = _mm_add_sd(sum_ai,tmp_ai);
			
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
			chrule   = _mm_mul_sd(xmm2,rinv);
			
			/* OFFSET INTERACTION ai->aj STARTS HERE */
			/* conditional mask for raj<dr+sk_ai */
			xmm1      = _mm_add_sd(r,sk_ai);
			mask_cmp  = _mm_cmplt_sd(raj,xmm1);
			
			/* conditional for rai>dr-sk, ends with mask_cmp2 */
			xmm2      = _mm_sub_sd(r,sk_ai);
			xmm3      = my_inv_pd(xmm2);
			mask_cmp2 = _mm_cmpgt_pd(raj,xmm2);
			
			lij	      = _mm_or_pd(_mm_and_pd(mask_cmp2,raj_inv)  ,_mm_andnot_pd(mask_cmp2,xmm3)); /*conditional as a mask*/
			dlij      = _mm_or_pd(_mm_and_pd(mask_cmp2,zero) ,_mm_andnot_pd(mask_cmp2,one));
			
			uij       = my_inv_pd(xmm1);
			lij2      = _mm_mul_sd(lij,lij);
			lij3      = _mm_mul_sd(lij2,lij);
			uij2      = _mm_mul_sd(uij,uij);
			uij3      = _mm_mul_sd(uij2,uij);
			
			diff2     = _mm_sub_sd(uij2,lij2);
			
			lij_inv   = my_invrsq_pd(lij2);
			
			sk2       = sk2_ai;
			sk2_inv   = _mm_mul_sd(sk2,rinv);
			prod      = _mm_mul_sd(qrtr,sk2_inv);
			
			log_term  = _mm_mul_pd(uij,lij_inv);
			log_term = log_pd(log_term);
			
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
			
			/* contitional for raj<sk_ai-dr */
			xmm3      = _mm_sub_sd(sk_ai,r);
			mask_cmp3 = _mm_cmplt_sd(raj,xmm3); /* rai<sk-dr */
			
			xmm4    = _mm_sub_sd(raj_inv,lij);
			xmm4    = _mm_mul_sd(two,xmm4);
			xmm4    = _mm_add_sd(tmp,xmm4);
			
			tmp	    = _mm_or_pd(_mm_and_pd(mask_cmp3,xmm4)  ,_mm_andnot_pd(mask_cmp3,tmp)); /*conditional as a mask*/
			
			/* the tmp will now contain two partial values, that not all are to be used. Which */
			/* ones are governed by the mask_cmp mask. */
			tmp     = _mm_mul_pd(half,tmp); 
			tmp     = _mm_or_pd(_mm_and_pd(mask_cmp,tmp)  ,_mm_andnot_pd(mask_cmp,zero)); /*conditional as a mask*/
			
			/* Load, add and store gpol energy */
			xmm7    = _mm_load_sd(born->gpol_hct_work+aj1);
			xmm7    = _mm_add_sd(xmm7,tmp);
			_mm_store_sd(born->gpol_hct_work+aj1,xmm7);
			
			/* Start chain rule terms, t1 */
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
			chrule_ai   = _mm_mul_sd(xmm2,rinv);
			
			_mm_store_sd(fr->dadx+n, chrule); 
			n = n + 1;
			_mm_store_sd(fr->dadx+n, chrule_ai);
			n = n + 1;
		}
		
		/* Do end processing ...  */
		tmp_ai = _mm_unpacklo_pd(tmp_ai,sum_ai);
		sum_ai = _mm_add_pd(sum_ai,tmp_ai);
		sum_ai = _mm_shuffle_pd(sum_ai,sum_ai,_MM_SHUFFLE2(1,1));
		
		/* Load, add and store atom ai polarisation energy */
		xmm2 = _mm_load_sd(born->gpol_hct_work+ai);
		sum_ai = _mm_add_sd(sum_ai,xmm2);
		_mm_store_sd(born->gpol_hct_work+ai,sum_ai);
	}
	
	/* Parallel summations */
	if(PARTDECOMP(cr))
	{
		gmx_sum(natoms, born->gpol_hct_work, cr);
	}
	else if(DOMAINDECOMP(cr))
	{
		dd_atom_sum_real(cr->dd, born->gpol_hct_work);
	}
	
	/* Compute the radii */
	for(i=0;i<nl->nri;i++)
	{
		ai      = nl->iinr[i];
		rr      = top->atomtypes.gb_radius[md->typeA[ai]];
		rr_inv2 = 1.0/rr;
		rr      = rr-doff; 
		rr_inv  = 1.0/rr;
		sum     = rr * born->gpol_hct_work[ai];
		sum2    = sum  * sum;
		sum3    = sum2 * sum;
		
		tsum    = tanh(born->obc_alpha*sum-born->obc_beta*sum2+born->obc_gamma*sum3);
		born->bRad[ai] = rr_inv - tsum*rr_inv2;
		born->bRad[ai] = 1.0 / born->bRad[ai];
		
		fr->invsqrta[ai]=invsqrt(born->bRad[ai]);
		
		tchain  = rr * (born->obc_alpha-2*born->obc_beta*sum+3*born->obc_gamma*sum2);
		born->drobc[ai] = (1.0-tsum*tsum)*tchain*rr_inv2;
	}
	
	/* Extra (local) communication required for DD */
	if(DOMAINDECOMP(cr))
	{
		dd_atom_spread_real(cr->dd, born->bRad);
		dd_atom_spread_real(cr->dd, fr->invsqrta);
		dd_atom_spread_real(cr->dd, born->drobc);
	}
	
	
	return 0;

}


int
calc_gb_chainrule_sse2_double(int natoms, t_nblist *nl, double *dadx, double *dvda, double *xd, double *f, int gb_algorithm, gmx_genborn_t *born)
{
	int i,k,n,ai,aj,ai3,aj1,aj2,aj13,aj23,aj4,nj0,nj1,offset;
	double rbi;
	double *rb;
	
	__m128d ix,iy,iz,jx,jy,jz,fix,fiy,fiz;
	__m128d dx,dy,dz,t1,t2,t3,dva,dvaj,dax,dax_ai,fgb,fgb_ai;
	__m128d xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xmm8;
	
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
			
			/* Load dvda_j */
			dvaj = _mm_set_pd(rb[aj2],rb[aj1]);
			
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
			xmm7  = _mm_loadu_pd(dadx+n);
			n     = n + 2;
			xmm8  = _mm_loadu_pd(dadx+n);
			n     = n + 2;
			
			/* Shuffle to get the ai and aj components right */ 
			dax    = _mm_shuffle_pd(xmm7,xmm8,_MM_SHUFFLE2(2,0));
			dax_ai = _mm_shuffle_pd(xmm7,xmm8,_MM_SHUFFLE2(3,1));
			
			/* distances */
			dx   = _mm_sub_pd(ix,jx);
			dy   = _mm_sub_pd(iy,jy);
			dz   = _mm_sub_pd(iz,jz);

			/* scalar force */
			fgb  = _mm_mul_pd(dva,dax); 
			fgb_ai = _mm_mul_pd(dvaj,dax_ai);
			fgb = _mm_add_pd(fgb,fgb_ai);
		
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
			
			dvaj = _mm_load_sd(rb+aj1);
			
			aj13 = aj1 * 3;
			
			jx    = _mm_load_sd(xd+aj13);
			jy    = _mm_load_sd(xd+aj13+1);
			jz    = _mm_load_sd(xd+aj13+2);
			
			dax  = _mm_load_sd(dadx+n);
			n    = n + 1;
			dax_ai = _mm_load_sd(dadx+n);
			n    = n + 1;
			
			dx   = _mm_sub_sd(ix,jx);
			dy   = _mm_sub_sd(iy,jy);
			dz   = _mm_sub_sd(iz,jz);
			
			fgb  = _mm_mul_sd(dva,dax); 
			fgb_ai = _mm_mul_sd(dvaj,dax_ai);
			fgb  = _mm_add_pd(fgb,fgb_ai);
			
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

