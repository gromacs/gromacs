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
#include <gmx_sse2_double.h>
#include <xmmintrin.h>
#include <emmintrin.h>

#include "genborn_sse2_double.h"

int
calc_gb_rad_still_sse2_double(t_commrec *cr, t_forcerec *fr,int natoms, gmx_localtop_t *top,
							  const t_atomtypes *atype, double *x, t_nblist *nl, gmx_genborn_t *born, t_mdatoms *md)
{
	int i,k,n,ai,ai3,aj1,aj2,aj13,aj23;
	int at0,at1,nj0,nj1,offset,taj1,taj2;
	int shift;

	double factor,gpi_ai,gpi_tmp,gpi2;
	double shX,shY,shZ;
	
	__m128d ix,iy,iz,jx,jy,jz,dx,dy,dz,sX,sY,sZ;
	__m128d t1,t2,t3,rsq11,rinv,rinv2,rinv4,rinv6;
	__m128d ratio,gpi,rai,raj,vaj,rvdw,mask_cmp;
	__m128d ccf,dccf,theta,cosq,term,sinq,res,prod;
	__m128d xmm1,xmm2,xmm3,xmm4,xmm7,vai,prod_ai,icf4,icf6;
	
	const __m128d half  = _mm_set1_pd(0.5);
	const __m128d three = _mm_set1_pd(3.0);
	const __m128d one   = _mm_set1_pd(1.0);
	const __m128d two   = _mm_set1_pd(2.0);
	const __m128d zero  = _mm_set1_pd(0.0);
	const __m128d four  = _mm_set1_pd(4.0);
	
	const __m128d p5inv  = _mm_set1_pd(STILL_P5INV);
	const __m128d pip5   = _mm_set1_pd(STILL_PIP5);
	const __m128d p4     = _mm_set1_pd(STILL_P4);
	
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
		
		nj0     = nl->jindex[i];
		nj1     = nl->jindex[i+1];
		
		/* Load shifts for this list */
		shift   = nl->shift[i];
		shX     = fr->shift_vec[shift][0];
		shY     = fr->shift_vec[shift][1];
		shZ     = fr->shift_vec[shift][2];
		
		/* Splat the shifts */
		sX = _mm_load1_pd(&shX);
		sY = _mm_load1_pd(&shY);
		sZ = _mm_load1_pd(&shZ);
		
		offset  = (nj1-nj0)%2;
		
		/* Polarization energy for atom ai */
		gpi     = _mm_setzero_pd();
	
		/* Load particle ai coordinates and add shifts */
		ix      = _mm_load1_pd(x+ai3);
		iy      = _mm_load1_pd(x+ai3+1);
		iz      = _mm_load1_pd(x+ai3+2);
		
		ix      = _mm_add_pd(sX,ix);
		iy      = _mm_add_pd(sY,iy);
		iz      = _mm_add_pd(sZ,iz);
		
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
			rinv     = gmx_mm_invsqrt_pd(rsq11);
			
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
					gmx_mm_sincos_pd(theta,&sinq,&cosq);
					
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
			rinv  = gmx_mm_invsqrt_pd(rsq11);
			
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
					gmx_mm_sincos_pd(theta,&sinq,&cosq);
										
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
	
	/* Parallel summations */
	if(PARTDECOMP(cr))
	{
		gmx_sum(natoms,born->gpol_still_work, cr);
	}
	else if(DOMAINDECOMP(cr))
	{
		dd_atom_sum_real(cr->dd, born->gpol_still_work);
	}
	
	/* Compute the radii */
	for(i=0;i<fr->natoms_force;i++) /* PELA born->nr */
	{
		if(born->use[i] != 0)
		{
			gpi_ai = born->gpol[i] + born->gpol_still_work[i];
			gpi2   = gpi_ai*gpi_ai;
			
			born->bRad[i]=factor*gmx_invsqrt(gpi2);
			fr->invsqrta[i]=gmx_invsqrt(born->bRad[i]);
		}
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
	int p1, p2, shift;
	double rr,sum,sum_tmp,min_rad,rad,doff;
	double shX,shY,shZ;
	
	__m128d ix,iy,iz,jx,jy,jz,dx,dy,dz,sX,sY,sZ;
	__m128d t1,t2,t3,rsq11,rinv,r,rai;
	__m128d rai_inv,sk,sk2,lij,dlij,duij;
	__m128d uij,lij2,uij2,lij3,uij3,diff2;
	__m128d lij_inv,sk2_inv,prod,log_term,tmp,tmp_sum;
	__m128d mask_cmp,mask_cmp2,mask_cmp3;
	__m128d xmm1,xmm2,xmm3,xmm4,xmm7,xmm8,xmm9,doffset;
	__m128d sum_ai,chrule,chrule_ai,tmp_ai; 
	__m128d sk_ai, sk2_ai,raj,raj_inv;
	
	const __m128d neg   = _mm_set1_pd(-1.0);
	const __m128d zero  = _mm_set1_pd(0.0);
	const __m128d eigth = _mm_set1_pd(0.125);
	const __m128d qrtr  = _mm_set1_pd(0.25);
	const __m128d half  = _mm_set1_pd(0.5);
	const __m128d one   = _mm_set1_pd(1.0);
	const __m128d two   = _mm_set1_pd(2.0);
	const __m128d three = _mm_set1_pd(3.0);
	
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
		
		nj0     = nl->jindex[i];
		nj1     = nl->jindex[i+1];
		
		/* Load shifts for this list */
		shift   = nl->shift[i];
		shX     = fr->shift_vec[shift][0];
		shY     = fr->shift_vec[shift][1];
		shZ     = fr->shift_vec[shift][2];
		
		/* Splat the shifts */
		sX = _mm_load1_pd(&shX);
		sY = _mm_load1_pd(&shY);
		sZ = _mm_load1_pd(&shZ);
		
		offset  = (nj1-nj0)%2;
		
		/* Load rai */
		rr      = top->atomtypes.gb_radius[md->typeA[ai]]-doff;
		rai     = _mm_load1_pd(&rr);
		rr      = 1.0/rr;
		rai_inv = _mm_load1_pd(&rr);
		
		/* Zero out sums for polarisation energies */
		sum_ai  = _mm_setzero_pd();
		
		/* Load ai coordinates and add shifts */
		ix       = _mm_load1_pd(x+ai3);
		iy       = _mm_load1_pd(x+ai3+1);
		iz       = _mm_load1_pd(x+ai3+2);
		
		ix      = _mm_add_pd(sX,ix);
		iy      = _mm_add_pd(sY,iy);
		iz      = _mm_add_pd(sZ,iz);
		
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
			rinv       = gmx_mm_invsqrt_pd(rsq11);
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
			raj_inv   = gmx_mm_inv_pd(raj);
			
			/* INTERACTION aj->ai STARS HERE */
			/* conditional mask for rai<dr+sk */
			xmm1       = _mm_add_pd(r,sk);
			mask_cmp   = _mm_cmplt_pd(rai,xmm1);
			
			/* conditional for rai>dr-sk, ends with mask_cmp2 */
			xmm2      = _mm_sub_pd(r,sk);
			xmm3      = gmx_mm_inv_pd(xmm2);
			mask_cmp2 = _mm_cmpgt_pd(rai,xmm2);
			
			lij	      = _mm_or_pd(_mm_and_pd(mask_cmp2,rai_inv)  ,_mm_andnot_pd(mask_cmp2,xmm3)); /*conditional as a mask*/
			dlij      = _mm_or_pd(_mm_and_pd(mask_cmp2,zero) ,_mm_andnot_pd(mask_cmp2,one));
			
			uij       = gmx_mm_inv_pd(xmm1);
			lij2      = _mm_mul_pd(lij,lij);
			lij3      = _mm_mul_pd(lij2,lij);
			uij2      = _mm_mul_pd(uij,uij);
			uij3      = _mm_mul_pd(uij2,uij);
			
			diff2     = _mm_sub_pd(uij2,lij2);
			
			lij_inv   = gmx_mm_invsqrt_pd(lij2);
			sk2       = _mm_mul_pd(sk,sk);
			sk2_inv   = _mm_mul_pd(sk2,rinv);
			prod      = _mm_mul_pd(qrtr,sk2_inv);
			
			log_term  = _mm_mul_pd(uij,lij_inv);
			log_term  = gmx_mm_log_pd(log_term);
			
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
			xmm2   = _mm_mul_pd(xmm2,xmm9); 
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
			xmm3  = gmx_mm_inv_pd(xmm2);
			mask_cmp2 = _mm_cmpgt_pd(raj,xmm2);
			
			lij	      = _mm_or_pd(_mm_and_pd(mask_cmp2,raj_inv)  ,_mm_andnot_pd(mask_cmp2,xmm3)); /*conditional as a mask*/
			dlij      = _mm_or_pd(_mm_and_pd(mask_cmp2,zero) ,_mm_andnot_pd(mask_cmp2,one));
			
			uij       = gmx_mm_inv_pd(xmm1);
			lij2      = _mm_mul_pd(lij,lij);
			lij3      = _mm_mul_pd(lij2,lij);
			uij2      = _mm_mul_pd(uij,uij);
			uij3      = _mm_mul_pd(uij2,uij);
			
			diff2     = _mm_sub_pd(uij2,lij2);
			
			lij_inv   = gmx_mm_invsqrt_pd(lij2);
			
			sk2       = sk2_ai;
			sk2_inv   = _mm_mul_pd(sk2,rinv);
			prod      = _mm_mul_pd(qrtr,sk2_inv);
			
			log_term  = _mm_mul_pd(uij,lij_inv);
			log_term  = gmx_mm_log_pd(log_term);
			
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
			xmm2   = _mm_mul_pd(xmm2,xmm9); 
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
			rinv      = gmx_mm_invsqrt_pd(rsq11);
			r         = _mm_mul_sd(rinv,rsq11);
			
			/* Load raj */
			raj       = _mm_load_sd(top->atomtypes.gb_radius+p1);
			raj       = _mm_sub_sd(raj,doffset);
			raj_inv   = gmx_mm_inv_pd(raj);
			
			/* OFFSET INTERATIONS aj->ai STARTS HERE */
			/* conditional mask for rai<dr+sk */
			xmm1      = _mm_add_sd(r,sk);
			mask_cmp  = _mm_cmplt_sd(rai,xmm1);
			
			/* conditional for rai>dr-sk, ends with mask_cmp2 */
			xmm2      = _mm_sub_sd(r,sk);
			xmm3      = gmx_mm_inv_pd(xmm2);
			mask_cmp2 = _mm_cmpgt_pd(rai,xmm2);
			
			lij	      = _mm_or_pd(_mm_and_pd(mask_cmp2,rai_inv)  ,_mm_andnot_pd(mask_cmp2,xmm3)); /*conditional as a mask*/
			dlij      = _mm_or_pd(_mm_and_pd(mask_cmp2,zero) ,_mm_andnot_pd(mask_cmp2,one));
			
			uij       = gmx_mm_inv_pd(xmm1);
			lij2      = _mm_mul_sd(lij,lij);
			lij3      = _mm_mul_sd(lij2,lij);
			uij2      = _mm_mul_sd(uij,uij);
			uij3      = _mm_mul_sd(uij2,uij);
			
			diff2     = _mm_sub_sd(uij2,lij2);
			
			lij_inv   = gmx_mm_invsqrt_pd(lij2);
			sk2       = _mm_mul_sd(sk,sk);
			sk2_inv   = _mm_mul_sd(sk2,rinv);
			prod      = _mm_mul_sd(qrtr,sk2_inv);
		
			log_term  = _mm_mul_pd(uij,lij_inv);
			log_term  = gmx_mm_log_pd(log_term);
			
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
			xmm2   = _mm_mul_sd(xmm2,xmm9); 
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
			xmm3      = gmx_mm_inv_pd(xmm2);
			mask_cmp2 = _mm_cmpgt_pd(raj,xmm2);
			
			lij	      = _mm_or_pd(_mm_and_pd(mask_cmp2,raj_inv)  ,_mm_andnot_pd(mask_cmp2,xmm3)); /*conditional as a mask*/
			dlij      = _mm_or_pd(_mm_and_pd(mask_cmp2,zero) ,_mm_andnot_pd(mask_cmp2,one));
			
			uij       = gmx_mm_inv_pd(xmm1);
			lij2      = _mm_mul_sd(lij,lij);
			lij3      = _mm_mul_sd(lij2,lij);
			uij2      = _mm_mul_sd(uij,uij);
			uij3      = _mm_mul_sd(uij2,uij);
			
			diff2     = _mm_sub_sd(uij2,lij2);
			
			lij_inv   = gmx_mm_invsqrt_pd(lij2);
			
			sk2       = sk2_ai;
			sk2_inv   = _mm_mul_sd(sk2,rinv);
			prod      = _mm_mul_sd(qrtr,sk2_inv);
			
			log_term  = _mm_mul_pd(uij,lij_inv);
			log_term  = gmx_mm_log_pd(log_term);
			
			xmm1      = _mm_sub_sd(lij,uij);
			xmm2      = _mm_mul_sd(qrtr,r);
			xmm2      = _mm_mul_sd(xmm2,diff2);
			xmm1      = _mm_add_sd(xmm1,xmm2); 
			xmm2      = _mm_mul_sd(half,rinv); 
			xmm2      = _mm_mul_sd(xmm2,log_term); 
			xmm1      = _mm_add_sd(xmm1,xmm2); 
			xmm9      = _mm_mul_sd(neg,diff2); 
			xmm2      = _mm_mul_sd(xmm9,prod); 
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
			xmm2   = _mm_mul_sd(xmm2,xmm9); 
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
	for(i=0;i<fr->natoms_force;i++) /* PELA born->nr */
	{
		if(born->use[i] != 0)
		{
			rr      = top->atomtypes.gb_radius[md->typeA[i]]-doff; 
			sum     = 1.0/rr - born->gpol_hct_work[i];
			min_rad = rr + doff;
			rad     = 1.0/sum;  
			
			born->bRad[i]   = rad > min_rad ? rad : min_rad;
			fr->invsqrta[i] = gmx_invsqrt(born->bRad[i]);
		}
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
calc_gb_rad_hct_obc_sse2_double(t_commrec *cr, t_forcerec * fr, int natoms, gmx_localtop_t *top,
                                const t_atomtypes *atype, double *x, t_nblist *nl, gmx_genborn_t *born,t_mdatoms *md, int gb_algorithm)
{
	int i,k,n,ai,ai3,aj1,aj2,aj13,aj23,nj0,nj1,at0,at1,offset;
	int p1,p2,p3,p4,shift;
	double rr,sum,sum_tmp,sum2,sum3,min_rad,rad,doff;
	double tsum,tchain,rr_inv,rr_inv2,gbr;
	double shX,shY,shZ;
	
	__m128d ix,iy,iz,jx,jy,jz,dx,dy,dz,sX,sY,sZ;
	__m128d t1,t2,t3,rsq11,rinv,r,rai;
	__m128d rai_inv,sk,sk2,lij,dlij,duij;
	__m128d uij,lij2,uij2,lij3,uij3,diff2;
	__m128d lij_inv,sk2_inv,prod,log_term,tmp,tmp_sum;
	__m128d mask_cmp,mask_cmp2,mask_cmp3,doffset,raj,raj_inv;
	__m128d xmm1,xmm2,xmm3,xmm4,xmm7,xmm8,xmm9;
	__m128d sum_ai,chrule,chrule_ai,tmp_ai,sk_ai,sk2_ai;
	
	const __m128d neg   = _mm_set1_pd(-1.0);
	const __m128d zero  = _mm_set1_pd(0.0);
	const __m128d eigth = _mm_set1_pd(0.125);
	const __m128d qrtr  = _mm_set1_pd(0.25);
	const __m128d half  = _mm_set1_pd(0.5);
	const __m128d one   = _mm_set1_pd(1.0);
	const __m128d two   = _mm_set1_pd(2.0);
    const __m128d three = _mm_set1_pd(3.0);
	
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
		
		nj0     = nl->jindex[i];
		nj1     = nl->jindex[i+1];
		
		/* Load shifts for this list */
		shift   = nl->shift[i];
		shX     = fr->shift_vec[shift][0];
		shY     = fr->shift_vec[shift][1];
		shZ     = fr->shift_vec[shift][2];
		
		/* Splat the shifts */
		sX = _mm_load1_pd(&shX);
		sY = _mm_load1_pd(&shY);
		sZ = _mm_load1_pd(&shZ);
		
		offset  = (nj1-nj0)%2;
		
		/* Load rai */
		rr      = top->atomtypes.gb_radius[md->typeA[ai]]-doff;
		rai     = _mm_load1_pd(&rr);
		rr      = 1.0/rr;
		rai_inv = _mm_load1_pd(&rr);
		
		/* Load ai coordinates and add shifts */
		ix       = _mm_load1_pd(x+ai3);
		iy       = _mm_load1_pd(x+ai3+1);
		iz       = _mm_load1_pd(x+ai3+2);
		
		ix      = _mm_add_pd(sX,ix);
		iy      = _mm_add_pd(sY,iy);
		iz      = _mm_add_pd(sZ,iz);
		
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
			rinv       = gmx_mm_invsqrt_pd(rsq11);
			r          = _mm_mul_pd(rinv,rsq11);
			
			/* Load atom aj1,aj2 raj */
			p1         = md->typeA[aj1];
			p2         = md->typeA[aj2];
			
			raj        = _mm_loadl_pd(raj,top->atomtypes.gb_radius+p1);
			raj        = _mm_loadh_pd(raj,top->atomtypes.gb_radius+p2);
			raj        = _mm_sub_pd(raj,doffset);
			
			/* Compute 1.0/raj */
			raj_inv    = gmx_mm_inv_pd(raj);
			
			sk         = _mm_loadl_pd(sk,born->param+aj1);
			sk         = _mm_loadh_pd(sk,born->param+aj2);
			
			/* INTERACTION aj->ai STARTS HERE */
			/* conditional mask for rai<dr+sk */
			xmm1       = _mm_add_pd(r,sk);
			mask_cmp   = _mm_cmplt_pd(rai,xmm1);
			
			/* conditional for rai>dr-sk, ends with mask_cmp2 */
			xmm2      = _mm_sub_pd(r,sk);
			xmm3      = gmx_mm_inv_pd(xmm2);
			mask_cmp2 = _mm_cmpgt_pd(rai,xmm2);
			
			lij	      = _mm_or_pd(_mm_and_pd(mask_cmp2,rai_inv)  ,_mm_andnot_pd(mask_cmp2,xmm3)); /*conditional as a mask*/
			dlij      = _mm_or_pd(_mm_and_pd(mask_cmp2,zero) ,_mm_andnot_pd(mask_cmp2,one));
			
			uij       = gmx_mm_inv_pd(xmm1);
			lij2      = _mm_mul_pd(lij,lij);
			lij3      = _mm_mul_pd(lij2,lij);
			uij2      = _mm_mul_pd(uij,uij);
			uij3      = _mm_mul_pd(uij2,uij);
			
			diff2     = _mm_sub_pd(uij2,lij2);
			
			lij_inv   = gmx_mm_invsqrt_pd(lij2);
			sk2       = _mm_mul_pd(sk,sk);
			sk2_inv   = _mm_mul_pd(sk2,rinv);
			prod      = _mm_mul_pd(qrtr,sk2_inv);
			
			log_term  = _mm_mul_pd(uij,lij_inv);
			log_term  = gmx_mm_log_pd(log_term);
			
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
			xmm2   = _mm_mul_pd(xmm2,xmm9); 
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
			xmm3      = gmx_mm_inv_pd(xmm2);
			mask_cmp2 = _mm_cmpgt_pd(raj,xmm2);
			
			lij	      = _mm_or_pd(_mm_and_pd(mask_cmp2,raj_inv)  ,_mm_andnot_pd(mask_cmp2,xmm3)); /*conditional as a mask*/
			dlij      = _mm_or_pd(_mm_and_pd(mask_cmp2,zero) ,_mm_andnot_pd(mask_cmp2,one));
			
			uij       = gmx_mm_inv_pd(xmm1);
			lij2      = _mm_mul_pd(lij,lij);
			lij3      = _mm_mul_pd(lij2,lij);
			uij2      = _mm_mul_pd(uij,uij);
			uij3      = _mm_mul_pd(uij2,uij);
			
			diff2     = _mm_sub_pd(uij2,lij2);
			
			lij_inv   = gmx_mm_invsqrt_pd(lij2);
			
			sk2       = sk2_ai;
			sk2_inv   = _mm_mul_pd(sk2,rinv);
			prod      = _mm_mul_pd(qrtr,sk2_inv);
			
			log_term  = _mm_mul_pd(uij,lij_inv);
			log_term  = gmx_mm_log_pd(log_term);
			
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
			xmm2   = _mm_mul_pd(xmm2,xmm9); 
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
			raj_inv    = gmx_mm_inv_pd(raj);
			
			dx        = _mm_sub_sd(ix,jx);
			dy        = _mm_sub_pd(iy,jy);
			dz        = _mm_sub_pd(iz,jz);
			
			rsq11     = _mm_add_sd( _mm_add_sd( _mm_mul_sd(dx,dx) , _mm_mul_sd(dy,dy) ) , _mm_mul_sd(dz,dz) );
			rinv      = gmx_mm_invsqrt_pd(rsq11);
			r         = _mm_mul_sd(rinv,rsq11);
			
			/* OFFSET INTERACTION aj->ai STARTS HERE */
			/* conditional mask for rai<dr+sk */
			xmm1      = _mm_add_sd(r,sk);
			mask_cmp  = _mm_cmplt_sd(rai,xmm1);
			
			/* conditional for rai>dr-sk, ends with mask_cmp2 */
			xmm2      = _mm_sub_sd(r,sk);
			xmm3      = gmx_mm_inv_pd(xmm2);
			mask_cmp2 = _mm_cmpgt_pd(rai,xmm2);
			
			lij	      = _mm_or_pd(_mm_and_pd(mask_cmp2,rai_inv)  ,_mm_andnot_pd(mask_cmp2,xmm3)); /*conditional as a mask*/
			dlij      = _mm_or_pd(_mm_and_pd(mask_cmp2,zero) ,_mm_andnot_pd(mask_cmp2,one));
			
			uij       = gmx_mm_inv_pd(xmm1);
			lij2      = _mm_mul_sd(lij,lij);
			lij3      = _mm_mul_sd(lij2,lij);
			uij2      = _mm_mul_sd(uij,uij);
			uij3      = _mm_mul_sd(uij2,uij);
			
			diff2     = _mm_sub_sd(uij2,lij2);
			
			lij_inv   = gmx_mm_invsqrt_pd(lij2);
			sk2       = _mm_mul_sd(sk,sk);
			sk2_inv   = _mm_mul_sd(sk2,rinv);
			prod      = _mm_mul_sd(qrtr,sk2_inv);
			
			log_term  = _mm_mul_pd(uij,lij_inv);
			log_term  = gmx_mm_log_pd(log_term);
			
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
			xmm2   = _mm_mul_sd(xmm2,xmm9); 
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
			xmm3      = gmx_mm_inv_pd(xmm2);
			mask_cmp2 = _mm_cmpgt_pd(raj,xmm2);
			
			lij	      = _mm_or_pd(_mm_and_pd(mask_cmp2,raj_inv)  ,_mm_andnot_pd(mask_cmp2,xmm3)); /*conditional as a mask*/
			dlij      = _mm_or_pd(_mm_and_pd(mask_cmp2,zero) ,_mm_andnot_pd(mask_cmp2,one));
			
			uij       = gmx_mm_inv_pd(xmm1);
			lij2      = _mm_mul_sd(lij,lij);
			lij3      = _mm_mul_sd(lij2,lij);
			uij2      = _mm_mul_sd(uij,uij);
			uij3      = _mm_mul_sd(uij2,uij);
			
			diff2     = _mm_sub_sd(uij2,lij2);
			
			lij_inv   = gmx_mm_invsqrt_pd(lij2);
			
			sk2       = sk2_ai;
			sk2_inv   = _mm_mul_sd(sk2,rinv);
			prod      = _mm_mul_sd(qrtr,sk2_inv);
			
			log_term  = _mm_mul_pd(uij,lij_inv);
			log_term  = gmx_mm_log_pd(log_term);
			
			xmm1      = _mm_sub_sd(lij,uij);
			xmm2      = _mm_mul_sd(qrtr,r);
			xmm2      = _mm_mul_sd(xmm2,diff2);
			xmm1      = _mm_add_sd(xmm1,xmm2); 
			xmm2      = _mm_mul_sd(half,rinv); 
			xmm2      = _mm_mul_sd(xmm2,log_term); 
			xmm1      = _mm_add_sd(xmm1,xmm2); 
			xmm9      = _mm_mul_sd(neg,diff2); 
			xmm2      = _mm_mul_sd(xmm9,prod); 
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
			xmm2   = _mm_mul_sd(xmm2,xmm9); 
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
	for(i=0;i<fr->natoms_force;i++) /* PELA born->nr */
	{
		if(born->use[i] != 0)
		{
			rr      = top->atomtypes.gb_radius[md->typeA[i]];
			rr_inv2 = 1.0/rr;
			rr      = rr-doff; 
			rr_inv  = 1.0/rr;
			sum     = rr * born->gpol_hct_work[i];
			sum2    = sum  * sum;
			sum3    = sum2 * sum;
			
			tsum    = tanh(born->obc_alpha*sum-born->obc_beta*sum2+born->obc_gamma*sum3);
			born->bRad[i] = rr_inv - tsum*rr_inv2;
			born->bRad[i] = 1.0 / born->bRad[i];
			
			fr->invsqrta[i] = gmx_invsqrt(born->bRad[i]);
			
			tchain  = rr * (born->obc_alpha-2*born->obc_beta*sum+3*born->obc_gamma*sum2);
			born->drobc[i] = (1.0-tsum*tsum)*tchain*rr_inv2;
		}
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
calc_gb_chainrule_sse2_double(int natoms, t_nblist *nl, double *dadx, double *dvda, double *xd, double *f, 
							  double *fshift, double *shift_vec, int gb_algorithm, gmx_genborn_t *born)
{
	int i,k,n,ai,aj,ai3,aj1,aj2,aj13,aj23,aj4,nj0,nj1,offset;
	int shift;
	double rbi,shX,shY,shZ;
	double *rb;
	
	__m128d ix,iy,iz,jx,jy,jz,fix,fiy,fiz,sX,sY,sZ;
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
		
		nj0    = nl->jindex[i];
		nj1    = nl->jindex[i+1];
		
		/* Load shifts for this list */
		shift   = 3*nl->shift[i];
		shX     = shift_vec[shift+0];
		shY     = shift_vec[shift+1];
		shZ     = shift_vec[shift+2];
		
		/* Splat the shifts */
		sX = _mm_load1_pd(&shX);
		sY = _mm_load1_pd(&shY);
		sZ = _mm_load1_pd(&shZ);
		
		offset = (nj1-nj0)%2;
		
		/* Load particle ai coordinates and add shifts */
		ix  = _mm_load1_pd(xd+ai3);
		iy  = _mm_load1_pd(xd+ai3+1);
		iz  = _mm_load1_pd(xd+ai3+2);
		
		ix      = _mm_add_pd(sX,ix);
		iy      = _mm_add_pd(sY,iy);
		iz      = _mm_add_pd(sZ,iz);
		
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
		
		/* load, add and store i forces */
		xmm1     = _mm_load_sd(f+ai3);
		xmm2     = _mm_load_sd(f+ai3+1);
		xmm3     = _mm_load_sd(f+ai3+2);
		
		fix      = _mm_add_sd(fix,xmm1);
		fiy      = _mm_add_sd(fiy,xmm2);
		fiz      = _mm_add_sd(fiz,xmm3);
		
		_mm_store_sd(f+ai3,fix);
		_mm_store_sd(f+ai3+1,fiy);
		_mm_store_sd(f+ai3+2,fiz);
		
		/* load, add and store i shift forces */
		xmm1     = _mm_load_sd(fshift+shift);
		xmm2     = _mm_load_sd(fshift+shift+1);
		xmm3     = _mm_load_sd(fshift+shift+2);
		
		fix      = _mm_add_sd(fix,xmm1);
		fiy      = _mm_add_sd(fiy,xmm2);
		fiz      = _mm_add_sd(fiz,xmm3);
		
		_mm_store_sd(fshift+shift,fix);
		_mm_store_sd(fshift+shift+1,fiy);
		_mm_store_sd(fshift+shift+2,fiz);

	}
	
	return 0;
}

#else
/* keep compiler happy */
int genborn_sse_dummy;

#endif /* SSE2 intrinsics available */

