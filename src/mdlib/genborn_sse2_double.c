#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef GMX_DOUBLE

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
#include "partdec.h"
#include "network.h"
#include "gmx_fatal.h"
#include "mtop_util.h"
#include "genborn.h"

#ifdef GMX_MPI
#include "mpi.h"
#endif

/* Only compile this file if SSE2 intrinsics are available */
#if ( defined(GMX_IA32_SSE2) || defined(GMX_X86_64_SSE2) )
#include <xmmintrin.h>
#include <emmintrin.h>

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



_PS_CONST(1  , 1.0f);
_PS_CONST(0p5, 0.5f);
/* the smallest non denormalized float number */
_PS_CONST_TYPE(min_norm_pos, int, 0x00800000);
_PS_CONST_TYPE(mant_mask, int, 0x7f800000);
_PS_CONST_TYPE(inv_mant_mask, int, ~0x7f800000);

_PS_CONST_TYPE(sign_mask, int, 0x80000000);
_PS_CONST_TYPE(inv_sign_mask, int, ~0x80000000);

_PI32_CONST(1, 1);
_PI32_CONST(inv1, ~1);
_PI32_CONST(2, 2);
_PI32_CONST(4, 4);
_PI32_CONST(0x7f, 0x7f);

_PS_CONST(cephes_SQRTHF, 0.707106781186547524);
_PS_CONST(cephes_log_p0, 7.0376836292E-2);
_PS_CONST(cephes_log_p1, - 1.1514610310E-1);
_PS_CONST(cephes_log_p2, 1.1676998740E-1);
_PS_CONST(cephes_log_p3, - 1.2420140846E-1);
_PS_CONST(cephes_log_p4, + 1.4249322787E-1);
_PS_CONST(cephes_log_p5, - 1.6668057665E-1);
_PS_CONST(cephes_log_p6, + 2.0000714765E-1);
_PS_CONST(cephes_log_p7, - 2.4999993993E-1);
_PS_CONST(cephes_log_p8, + 3.3333331174E-1);
_PS_CONST(cephes_log_q1, -2.12194440e-4);
_PS_CONST(cephes_log_q2, 0.693359375);

_PS_CONST(minus_cephes_DP1, -0.78515625);
_PS_CONST(minus_cephes_DP2, -2.4187564849853515625e-4);
_PS_CONST(minus_cephes_DP3, -3.77489497744594108e-8);
_PS_CONST(sincof_p0, -1.9515295891E-4);
_PS_CONST(sincof_p1,  8.3321608736E-3);
_PS_CONST(sincof_p2, -1.6666654611E-1);
_PS_CONST(coscof_p0,  2.443315711809948E-005);
_PS_CONST(coscof_p1, -1.388731625493765E-003);
_PS_CONST(coscof_p2,  4.166664568298827E-002);
_PS_CONST(cephes_FOPI, 1.27323954473516); /* 4 / M_PI */

_PS_CONST(exp_hi,	88.3762626647949f);
_PS_CONST(exp_lo,	-88.3762626647949f);

_PS_CONST(cephes_LOG2EF, 1.44269504088896341);
_PS_CONST(cephes_exp_C1, 0.693359375);
_PS_CONST(cephes_exp_C2, -2.12194440e-4);

_PS_CONST(cephes_exp_p0, 1.9875691500E-4);
_PS_CONST(cephes_exp_p1, 1.3981999507E-3);
_PS_CONST(cephes_exp_p2, 8.3334519073E-3);
_PS_CONST(cephes_exp_p3, 4.1665795894E-2);
_PS_CONST(cephes_exp_p4, 1.6666665459E-1);
_PS_CONST(cephes_exp_p5, 5.0000001201E-1);


#define COPY_XMM_TO_MM(xmm_, mm0_, mm1_) {          \
xmm_mm_union u; u.xmm = xmm_;                   \
mm0_ = u.mm[0];                                 \
mm1_ = u.mm[1];                                 \
}

#define COPY_MM_TO_XMM(mm0_, mm1_, xmm_) {                         \
xmm_mm_union u; u.mm[0]=mm0_; u.mm[1]=mm1_; xmm_ = u.xmm;      \
}

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

void sincos_ps(__m128 x, __m128 *s, __m128 *c) {
	__m128 xmm1, xmm2, xmm3 = _mm_setzero_ps(), sign_bit_sin, y;
	__m64 mm0, mm1, mm2, mm3, mm4, mm5;
	sign_bit_sin = x;
	/* take the absolute value */
	x = _mm_and_ps(x, *(__m128*)_ps_inv_sign_mask);
	/* extract the sign bit (upper one) */
	sign_bit_sin = _mm_and_ps(sign_bit_sin, *(__m128*)_ps_sign_mask);
	
	/* scale by 4/Pi */
	y = _mm_mul_ps(x, *(__m128*)_ps_cephes_FOPI);
    
	/* store the integer part of y in mm0:mm1 */
	xmm3 = _mm_movehl_ps(xmm3, y);
	mm2 = _mm_cvttps_pi32(y);
	mm3 = _mm_cvttps_pi32(xmm3);
	
	/* j=(j+1) & (~1) (see the cephes sources) */
	mm2 = _mm_add_pi32(mm2, *(__m64*)_pi32_1);
	mm3 = _mm_add_pi32(mm3, *(__m64*)_pi32_1);
	mm2 = _mm_and_si64(mm2, *(__m64*)_pi32_inv1);
	mm3 = _mm_and_si64(mm3, *(__m64*)_pi32_inv1);
	
	y = _mm_cvtpi32x2_ps(mm2, mm3);
	
	mm4 = mm2;
	mm5 = mm3;
	
	/* get the swap sign flag for the sine */
	mm0 = _mm_and_si64(mm2, *(__m64*)_pi32_4);
	mm1 = _mm_and_si64(mm3, *(__m64*)_pi32_4);
	mm0 = _mm_slli_pi32(mm0, 29);
	mm1 = _mm_slli_pi32(mm1, 29);
	__m128 swap_sign_bit_sin;
	COPY_MM_TO_XMM(mm0, mm1, swap_sign_bit_sin);
	
	/* get the polynom selection mask for the sine */
	
	mm2 = _mm_and_si64(mm2, *(__m64*)_pi32_2);
	mm3 = _mm_and_si64(mm3, *(__m64*)_pi32_2);
	mm2 = _mm_cmpeq_pi32(mm2, _mm_setzero_si64());
	mm3 = _mm_cmpeq_pi32(mm3, _mm_setzero_si64());
	__m128 poly_mask;
	COPY_MM_TO_XMM(mm2, mm3, poly_mask);
	
	/* The magic pass: "Extended precision modular arithmetic" 
     x = ((x - y * DP1) - y * DP2) - y * DP3; */
	xmm1 = *(__m128*)_ps_minus_cephes_DP1;
	xmm2 = *(__m128*)_ps_minus_cephes_DP2;
	xmm3 = *(__m128*)_ps_minus_cephes_DP3;
	xmm1 = _mm_mul_ps(y, xmm1);
	xmm2 = _mm_mul_ps(y, xmm2);
	xmm3 = _mm_mul_ps(y, xmm3);
	x = _mm_add_ps(x, xmm1);
	x = _mm_add_ps(x, xmm2);
	x = _mm_add_ps(x, xmm3);
	
	
	/* get the sign flag for the cosine */
	mm4 = _mm_sub_pi32(mm4, *(__m64*)_pi32_2);
	mm5 = _mm_sub_pi32(mm5, *(__m64*)_pi32_2);
	mm4 = _mm_andnot_si64(mm4, *(__m64*)_pi32_4);
	mm5 = _mm_andnot_si64(mm5, *(__m64*)_pi32_4);
	mm4 = _mm_slli_pi32(mm4, 29);
	mm5 = _mm_slli_pi32(mm5, 29);
	__m128 sign_bit_cos;
	COPY_MM_TO_XMM(mm4, mm5, sign_bit_cos);
	
	sign_bit_sin = _mm_xor_ps(sign_bit_sin, swap_sign_bit_sin);
	
	
	/* Evaluate the first polynom  (0 <= x <= Pi/4) */
	__m128 z = _mm_mul_ps(x,x);
	y = *(__m128*)_ps_coscof_p0;
	
	y = _mm_mul_ps(y, z);
	y = _mm_add_ps(y, *(__m128*)_ps_coscof_p1);
	y = _mm_mul_ps(y, z);
	y = _mm_add_ps(y, *(__m128*)_ps_coscof_p2);
	y = _mm_mul_ps(y, z);
	y = _mm_mul_ps(y, z);
	__m128 tmp = _mm_mul_ps(z, *(__m128*)_ps_0p5);
	y = _mm_sub_ps(y, tmp);
	y = _mm_add_ps(y, *(__m128*)_ps_1);
	
	/* Evaluate the second polynom  (Pi/4 <= x <= 0) */
	__m128 y2 = *(__m128*)_ps_sincof_p0;
	y2 = _mm_mul_ps(y2, z);
	y2 = _mm_add_ps(y2, *(__m128*)_ps_sincof_p1);
	y2 = _mm_mul_ps(y2, z);
	y2 = _mm_add_ps(y2, *(__m128*)_ps_sincof_p2);
	y2 = _mm_mul_ps(y2, z);
	y2 = _mm_mul_ps(y2, x);
	y2 = _mm_add_ps(y2, x);
	
	/* select the correct result from the two polynoms */  
	xmm3 = poly_mask;
	__m128 ysin2 = _mm_and_ps(xmm3, y2);
	__m128 ysin1 = _mm_andnot_ps(xmm3, y);
	y2 = _mm_sub_ps(y2,ysin2);
	y = _mm_sub_ps(y, ysin1);
	
	xmm1 = _mm_add_ps(ysin1,ysin2);
	xmm2 = _mm_add_ps(y,y2);
	
	/* update the sign */
	*s = _mm_xor_ps(xmm1, sign_bit_sin);
	*c = _mm_xor_ps(xmm2, sign_bit_cos);
	_mm_empty(); /* good-bye mmx */
}


__m128 log2_ps(__m128 x)
{
	__m128i exp   = _mm_set_epi32(0x7F800000, 0x7F800000, 0x7F800000, 0x7F800000);
	__m128i one   = _mm_set_epi32(0x3F800000, 0x3F800000, 0x3F800000, 0x3F800000); 
	__m128i off   = _mm_set_epi32(0x3FBF8000, 0x3FBF8000, 0x3FBF8000, 0x3FBF8000); 
	__m128i mant  = _mm_set_epi32(0x007FFFFF, 0x007FFFFF, 0x007FFFFF, 0x007FFFFF);
	__m128i sign  = _mm_set_epi32(0x80000000, 0x80000000, 0x80000000, 0x80000000);
	__m128i base  = _mm_set_epi32(0x43800000, 0x43800000, 0x43800000, 0x43800000);
	__m128i loge  = _mm_set_epi32(0x3F317218, 0x3F317218, 0x3F317218, 0x3F317218);
	
	__m128i D5     = _mm_set_epi32(0xBD0D0CC5, 0xBD0D0CC5, 0xBD0D0CC5, 0xBD0D0CC5);
	__m128i D4     = _mm_set_epi32(0x3EA2ECDD, 0x3EA2ECDD, 0x3EA2ECDD, 0x3EA2ECDD); 
	__m128i D3     = _mm_set_epi32(0xBF9dA2C9, 0xBF9dA2C9, 0xBF9dA2C9, 0xBF9dA2C9);
	__m128i D2     = _mm_set_epi32(0x4026537B, 0x4026537B, 0x4026537B, 0x4026537B);
	__m128i D1     = _mm_set_epi32(0xC054bFAD, 0xC054bFAD, 0xC054bFAD, 0xC054bFAD); 
	__m128i D0     = _mm_set_epi32(0x4047691A, 0x4047691A, 0x4047691A, 0x4047691A);
	
	__m128  xmm0,xmm1,xmm2;
	__m128i xmm1i;
	
	xmm0  = x;
	xmm1  = xmm0;
	xmm1  = _mm_and_ps(xmm1, (__m128) exp);
	xmm1 = (__m128) _mm_srli_epi32( (__m128i) xmm1,8); 
	
	xmm1  = _mm_or_ps(xmm1,(__m128) one);
	xmm1  = _mm_sub_ps(xmm1,(__m128) off);
	
	xmm1  = _mm_mul_ps(xmm1,(__m128) base);
	xmm0  = _mm_and_ps(xmm0,(__m128) mant);
	xmm0  = _mm_or_ps(xmm0,(__m128) one);
	
	xmm2  = _mm_mul_ps(xmm0, (__m128) D5);
	xmm2  = _mm_add_ps(xmm2, (__m128) D4);
	xmm2  = _mm_mul_ps(xmm2,xmm0);
	xmm2  = _mm_add_ps(xmm2, (__m128) D3);
	xmm2  = _mm_mul_ps(xmm2,xmm0);
	xmm2  = _mm_add_ps(xmm2, (__m128) D2);
	xmm2  = _mm_mul_ps(xmm2,xmm0);
	xmm2  = _mm_add_ps(xmm2, (__m128) D1);
	xmm2  = _mm_mul_ps(xmm2,xmm0);
	xmm2  = _mm_add_ps(xmm2, (__m128) D0);
	xmm0  = _mm_sub_ps(xmm0, (__m128) one);
	xmm0  = _mm_mul_ps(xmm0,xmm2);
	xmm1  = _mm_add_ps(xmm1,xmm0);
	
	
	x     = xmm1;
	x  = _mm_mul_ps(x,(__m128)loge);
	
    return x;
}

static inline __m128d
my_invrsq_pd(__m128d x)
{
	const __m128d three = (const __m128d) {3.0f, 3.0f};
	const __m128d half  = (const __m128d) {0.5f, 0.5f};
	
	__m128  t  = _mm_rsqrt_ps(_mm_cvtpd_ps(x)); /* Convert to single precision and do _mm_rsqrt_ps() */
	__m128d t1 = _mm_cvtps_pd(t); /* Convert back to double precision */
	
	/* First Newton-Rapson step, accuracy is now 24 bits */
	__m128d t2 = _mm_mul_pd(half,_mm_mul_pd(t1,_mm_sub_pd(three,_mm_mul_pd(x,_mm_mul_pd(t1,t1)))));
	
	/* Return second Newton-Rapson step, accuracy 48 bits */
	return (__m128d) _mm_mul_pd(half,_mm_mul_pd(t2,_mm_sub_pd(three,_mm_mul_pd(x,_mm_mul_pd(t2,t2)))));
}

static inline __m128d
my_inv_pd(__m128d x)
{
	const __m128d two = (const __m128d) {2.0f, 2.0f};
	
	__m128  t  = _mm_rcp_ps(_mm_cvtpd_ps(x));
	__m128d t1 = _mm_cvtps_pd(t);
	__m128d t2 = _mm_mul_pd(t1,_mm_sub_pd(two,_mm_mul_pd(t1,x)));
	
	return (__m128d) _mm_mul_pd(t2,_mm_sub_pd(two,_mm_mul_pd(t2,x)));
}

int
calc_gb_rad_still_sse2_double(t_commrec *cr, t_forcerec *fr,int natoms, gmx_mtop_t *mtop,
							  const t_atomtypes *atype, real *x, t_nblist *nl, gmx_genborn_t *born, t_mdatoms *md)
{
	int i,k,n,ai,ai3,aj1,aj2,aj13,aj23;
	int at0,at1,nj0,nj1,offset,taj1,taj2;

	real factor,gpi_ai,gpi_tmp,gpi2;
	real *sum_gpi;
	
	__m128d ix,iy,iz,jx,jy,jz,dx,dy,dz;
	__m128d t1,t2,t3,rsq11,rinv,rinv2,rinv4,rinv6;
	__m128d ratio,gpi,rai,raj,vaj,rvdw,mask_cmp;
	__m128d ccf,dccf,theta,cosq,term,sinq,res,prod;
	__m128d xmm1,xmm2,xmm3;
	__m128  tmp,sinq_single,cosq_single;
	
	const __m128d half  = {0.5f, 0.5f};
	const __m128d three = {3.0f, 3.0f};
	const __m128d one   = {1.0f, 1.0f};
	const __m128d two   = {2.0f, 2.0f};
	const __m128d zero  = {0.0f, 0.0f};
	const __m128d four  = {4.0f, 4.0f};
	
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
			sum_gpi[ai] = 0;
				
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
			sum_gpi[ai]=gpi_tmp;
		}
		else
		{
			gpi_ai = gpi_ai + gpi_tmp;
			gpi2   = gpi_ai * gpi_ai;
			
			born->bRad[ai]   = factor * invsqrt(gpi2);
			fr->invsqrta[ai] = invsqrt(born->bRad[ai]);
		}
	}

	if(PAR(cr))
	{
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
	
	
	return 0;
}

int 
calc_gb_rad_hct_sse2_double(t_commrec *cr, t_forcerec *fr, int natoms, gmx_mtop_t *mtop, 
							const t_atomtypes *atype, real *x, t_nblist *nl, gmx_genborn_t *born, t_mdatoms *md)
{
	int i,k,n,ai,ai3,aj1,aj2,aj13,aj23,nj0,nj1,at0,at1,offset;
	real ri,rr,sum,sum_tmp,sum_ai,min_rad,rad,doffset;
	real *sum_mpi;
	
	__m128d ix,iy,iz,jx,jy,jz,dx,dy,dz;
	__m128d t1,t2,t3,rsq11,rinv,r,rai;
	__m128d rai_inv,sk,sk2,lij,dlij,duij;
	__m128d uij,lij2,uij2,lij3,uij3,diff2;
	__m128d lij_inv,sk2_inv,prod,log_term,tmp,tmp_sum;
	__m128d mask_cmp,mask_cmp2,mask_cmp3;
	__m128d xmm1,xmm2,xmm3,xmm4,xmm8;
	__m128  log_term_single;
	
	const __m128d neg   = {-1.0f, -1.0f};
	const __m128d zero  = {0.0f, 0.0f};
	const __m128d eigth = {0.125f, 0.125f};
	const __m128d qrtr  = {0.25f, 0.25f};
	const __m128d half  = {0.5f, 0.5f};
	const __m128d one   = {1.0f, 1.0f};
	const __m128d two   = {2.0f, 2.0f};
	const __m128d three = {3.0f, 3.0f};
	
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
							const t_atomtypes *atype, real *x, t_nblist *nl, gmx_genborn_t *born,t_mdatoms *md)
{
	int i,k,n,ai,ai3,aj1,aj2,aj13,aj23,nj0,nj1,at0,at1,offset;
	real ri,rr,sum,sum_tmp,sum_ai,sum_ai2,sum_ai3,min_rad,rad,doffset;
	real tsum,tchain,rr_inv,gbr;
	real *sum_mpi;
	
	__m128d ix,iy,iz,jx,jy,jz,dx,dy,dz;
	__m128d t1,t2,t3,rsq11,rinv,r,rai;
	__m128d rai_inv,sk,sk2,lij,dlij,duij;
	__m128d uij,lij2,uij2,lij3,uij3,diff2;
	__m128d lij_inv,sk2_inv,prod,log_term,tmp,tmp_sum;
	__m128d mask_cmp,mask_cmp2,mask_cmp3;
	__m128d xmm1,xmm2,xmm3,xmm4,xmm8;
	__m128 log_term_single;
	
	const __m128d neg   = {-1.0f, -1.0f};
	const __m128d zero  = {0.0f, 0.0f};
	const __m128d eigth = {0.125f, 0.125f};
	const __m128d qrtr  = {0.25f, 0.25f};
	const __m128d half  = {0.5f, 0.5f};
	const __m128d one   = {1.0f, 1.0f};
	const __m128d two   = {2.0f, 2.0f};
	const __m128d three = {3.0f, 3.0f};
	
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

int
calc_gb_chainrule_sse2_double(int natoms, t_nblist *nl, real *dadx, real *dvda, real *xd, real *f, int gb_algorithm, gmx_genborn_t *born)
{
	int i,k,n,ai,aj,ai3,aj1,aj2,aj13,aj23,aj4,nj0,nj1,offset;
	real rbi;
	real rb[natoms+4];
	
	__m128d ix,iy,iz,jx,jy,jz,fix,fiy,fiz;
	__m128d dx,dy,dz,t1,t2,t3,dva,dax,fgb;
	__m128d xmm1,xmm2,xmm3,xmm4,xmm5;
	
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
			xmm1 = _mm_loadl_pd(xmm1,xd+aj13);
			xmm2 = _mm_loadl_pd(xmm2,xd+aj23);
			jx   = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(0,0));
			jy   = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(1,1));
			
			jz   = _mm_loadl_pd(jz, xd+aj13+2);
			jz   = _mm_loadh_pd(jz, xd+aj23+2);
			
			dax  = _mm_loadu_pd(dadx+n);
			n    = n + 2;
			
			dx   = _mm_sub_pd(ix,jx);
			dy   = _mm_sub_pd(iy,jy);
			dz   = _mm_sub_pd(iz,jz);

			fgb  = _mm_mul_pd(dva,dax); 
		
			t1   = _mm_mul_pd(fgb,dx); 
			t2   = _mm_mul_pd(fgb,dy); 
			t3   = _mm_mul_pd(fgb,dz); 
			
			fix  = _mm_add_pd(fix,t1);
			fiy  = _mm_add_pd(fiy,t2);
			fiz  = _mm_add_pd(fiz,t3);	
			
			/* accumulate the aj1-2 fx and fy forces from memory */
			xmm1 = _mm_loadl_pd(xmm1,f+aj13);
			xmm2 = _mm_loadl_pd(xmm2,f+aj23);
			xmm1 = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(0,0)); /* fx1 fx2 */
			xmm2 = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(1,1)); /* fy1 fy2 */
			
			xmm3   = _mm_loadl_pd(xmm3, f+aj13+2); 
			xmm3   = _mm_loadh_pd(xmm3, f+aj23+2); /* fz1 fz2 */
			
			/* subtract partial forces */
			xmm1 = _mm_sub_pd(xmm1, t1); 
			xmm2 = _mm_sub_pd(xmm2, t2); 
			xmm3 = _mm_sub_pd(xmm3, t3); 
			
			/* transpose back fx and fy */
			xmm4 = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(0,0));
			xmm5 = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(1,1));
			
			/* store the force, first fx's and fy's */
			_mm_storeu_pd(f+aj13,xmm4);
			_mm_storeu_pd(f+aj23,xmm5);
			
			/* then fz */
			_mm_storel_pd(f+aj13+2,xmm3);
			_mm_storeh_pd(f+aj23+2,xmm3);
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
			
			/* accumulate the aj1-2 fx and fy forces from memory */
			xmm1 = _mm_load_sd(f+aj13);
			xmm2 = _mm_load_sd(f+aj13+1);
			xmm3 = _mm_load_sd(f+aj13+2); 
					
			/* subtract partial forces */
			xmm1 = _mm_sub_sd(xmm1, t1); 
			xmm2 = _mm_sub_sd(xmm2, t2); 
			xmm3 = _mm_sub_sd(xmm3, t3); 
					
			/* store the force */
			_mm_store_sd(f+aj13,xmm1);
			_mm_store_sd(f+aj13,xmm2);
			_mm_store_sd(f+aj13,xmm3);
		}
		
		/* Do end processing ...  */
		
		
		/* ... until here */
		
	}

	return 0;
}

#else
/* keep compiler happy */
int genborn_sse_dummy;

#endif /* SSE2 intrinsics available */

#else
/* keep compiler happy */
int genborn_sse2_double_dummy;

#endif /* double precision available */
