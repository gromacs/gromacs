/*
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations 
 *
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org 
 * 
 * And Hey:
 * GROup of MAchos and Cynical Suckers
 */	

/* This file contains a subset of the gromacs innerloops
 * manually written in assembly to optimize performance
 * on AMD extended 3DNow-enabled processors like Athlon 
 * and later generations. 
 * Erik Lindahl, 2000-2001, erik@theophys.kth.se
 *
 * We use intel syntax for portability. There are probably some GNU-specific
 * things here, but they are easy to fix.
 */
	
.intel_syntax noprefix

.text
	
mm_two:	
	.long 0x40000000
	.long 0x40000000
mm_six:	
	.long 0x40c00000
	.long 0x40c00000
mm_twelve:	
	.long 0x41400000
	.long 0x41400000

	.align 4

.globl check3dnow  /* try to issue an Extended 3DNow instruction */
	.type check3dnow,@function
check3dnow:	
	femms
	pswapd mm0,mm0
	femms
	ret

			
.globl vecrecip_3dnow
	.type vecrecip_3dnow,@function
vecrecip_3dnow:	
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx

	mov eax, [ebp + 8]
	mov ebx, [ebp + 12]	
	mov ecx, [ebp + 16]
        mov edx, ecx
        shr ecx, 2 
        jecxz .vecrecip_tail
        emms	
.vecrecip_mainloop:	
        movq mm0,[eax]
	add eax,  8
        pfrcp mm1,mm0
	movq mm4,[eax]
	pswapd mm0,mm0
	add eax,  8 
        pfrcp mm2,mm0
	pswapd mm0,mm0
        pfrcp mm5,mm4
	pswapd mm4,mm4	
	punpckldq mm1,mm2
	pfrcp mm6,mm4
	pswapd mm4,mm4
	pfrcpit1 mm0,mm1
	punpckldq mm5,mm6	
	pfrcpit2 mm0,mm1
        movq [ebx],mm0
	pfrcpit1 mm4,mm5
	add ebx,  8
	pfrcpit2 mm4,mm5	
        movq [ebx],mm4
	add ebx,  8	
        dec ecx
        jecxz .vecrecip_tail
        jmp short .vecrecip_mainloop
.vecrecip_tail:
        mov ecx,edx
        and ecx,3
        jecxz .vecrecip_end
.vecrecip_tailloop:	
        movd mm0,[eax]
	add eax,  4
        pfrcp mm1,mm0
        pfrcpit1 mm0,mm1
        pfrcpit2 mm0,mm1
        movd [ebx],mm0	
	add ebx,  4
	dec ecx
	jecxz .vecrecip_end
	jmp short .vecrecip_tailloop
.vecrecip_end:	
	emms
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret
	

.globl vecinvsqrt_3dnow
	.type vecinvsqrt_3dnow,@function
vecinvsqrt_3dnow:	
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx

	mov eax, [ebp + 8]
	mov ebx, [ebp + 12]	
	mov ecx, [ebp + 16]
        mov edx, ecx
        shr ecx, 2 
        jecxz .vecinvsqrt_tail
        emms	
.vecinvsqrt_mainloop:	
        movq mm0,[eax]
	add eax,  8
        pfrsqrt mm1,mm0
	movq mm4,[eax]
	pswapd mm0,mm0
	add eax,  8
        pfrsqrt mm2,mm0
	pswapd mm0,mm0
        pfrsqrt mm5,mm4
	pswapd mm4,mm4	
	punpckldq mm1,mm2
	pfrsqrt mm6,mm4
	movq mm3,mm1
	pswapd mm4,mm4
	pfmul mm1,mm1
	punpckldq mm5,mm6	
	pfrsqit1 mm1,mm0
	movq mm7,mm5	
	pfrcpit2 mm1,mm3
	pfmul mm5,mm5
        movq [ebx],mm1
	pfrsqit1 mm5,mm4
	add ebx,  8
	pfrcpit2 mm5,mm7	
        movq [ebx],mm5
	add ebx,  8	
        dec ecx
        jecxz .vecinvsqrt_tail
        jmp short .vecinvsqrt_mainloop
.vecinvsqrt_tail:
        mov ecx,edx
        and ecx,3
        jecxz .vecinvsqrt_end
.vecinvsqrt_tailloop:	
        movd mm0,[eax]
	add eax,  4
        pfrsqrt mm1,mm0
        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0
        pfrcpit2 mm1,mm2
        movd [ebx],mm1		
	add ebx,  4
	dec ecx
	jecxz .vecinvsqrt_end
	jmp short .vecinvsqrt_tailloop
.vecinvsqrt_end:	
	emms
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret
	

.globl inl0100_3dnow
	.type inl0100_3dnow,@function
inl0100_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40
.equ		_faction,	44
.equ		_type,		48
.equ		_ntype,		52
.equ		_nbfp,		56
.equ		_Vnb,		60
	/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_ix,		8
.equ		_iy,		12
.equ		_iz,		16
.equ		_vnbtot,	20  
.equ		_c6,		28  
.equ		_c12,		36  
.equ		_six,		44  
.equ		_twelve,	52  
.equ		_ntia,		60
.equ		_innerjjnr,	64
.equ		_innerk,	68		
.equ		_fix,		72
.equ		_fiy,		76
.equ		_fiz,		80
.equ		_dx1,		84
.equ		_dy1,		88
.equ		_dz1,		92
.equ		_dx2,		96
.equ		_dy2,		100
.equ		_dz2,		104						
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 108		/* local stack space */
	femms
	/* move data to local stack */ 
	movq  mm0, [mm_six]
	movq  mm1, [mm_twelve]
	movq  [esp + _six ], mm0
	movq  [esp + _twelve ], mm1
	/* assume we have at least one i particle - start directly */	
.i0100_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift], 4		/* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm0, [eax + ebx*4]	/* move shX/shY to mm0 and shZ to mm1. */
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx =ii */

	mov   edx, [ebp + _type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + _ntype]
	shl   edx, 1
	mov   [esp + _ntia], edx

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	
	pfadd mm0, [eax + ebx*4]        /* ix = shX + posX (and iy too) */
	movd  mm3, [eax + ebx*4 + 8]    /* cant use direct memory add for 4 bytes (iz) */
	mov   [esp + _ii3], ebx
	pfadd mm1, mm3
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	
				
	/* clear total potential and i forces */
	pxor  mm7,mm7
	movq  [esp + _vnbtot], mm7
	movq  [esp + _fix],    mm7
	movd  [esp + _fiz],    mm7

	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */

	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]	
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr], eax     /*  pointer to jjnr[nj0] */
	sub   edx,  2
	mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i0100_unroll_loop
	jmp   .i0100_finish_inner
.i0100_unroll_loop:
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */	

	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	mov ecx, [ecx + ebx*4]           /* type [jnr2] */

	mov esi, [ebp + _nbfp]		/* base of nbfp */ 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	add ecx, [esp + _ntia]

	movq mm5, [esi + edx*4]		/* mm5 = 1st c6 / c12 */		
	movq mm7, [esi + ecx*4]		/* mm7 = 2nd c6 / c12 */	
	movq mm6,mm5			
	punpckldq mm5,mm7		/* mm5 = 1st c6 / 2nd c6 */
	punpckhdq mm6,mm7		/* mm6 = 1st c12 / 2nd c12 */
	movq [esp + _c6], mm5
	movq [esp + _c12], mm6

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrcp mm0, mm4	                 /* lookup reciprocal seed */ 
        pfrcp mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs. */
                  	        	/* amd 3dnow N-R iteration to get full precision. */
        pfrcpit1 mm4,mm0				
        pfrcpit2 mm4,mm0	
	/* mm4 now contains invsq,
	 * do potential and fscal
	 */
	movq  mm0, mm4
	pfmul mm4, mm0
	pfmul mm4, mm0             	/* mm4=rinvsix */
	movq  mm5, mm4	
	pfmul mm5, mm5	                /* mm5=rinvtwelve */

	pfmul mm5, [esp + _c12]
	pfmul mm4, [esp + _c6]	
	movq mm6, mm5	/* mm6 is vnb12-vnb6 */ 
	pfsub mm6, mm4

	pfmul mm4, [esp + _six]

	pfmul mm5, [esp + _twelve]
	pfsub mm5,mm4
 	pfmul mm0, mm5        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	/* update vnbtot */ 
	pfadd mm6, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm6       /* store the sum */

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i0100_finish_inner
	jmp   .i0100_unroll_loop
.i0100_finish_inner:	
	and [esp + _innerk],  1
	jnz  .i0100_single_inner
	jmp  .i0100_updateouterdata		
.i0100_single_inner:
	/* a single j particle iteration here - compare with the unrolled code for comments */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov esi, [ebp + _nbfp]
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	/* type [jnr1] */
	shl edx, 1
	add edx, [esp + _ntia]	        /* tja = ntia + 2*type */
	movd mm5, [esi + edx*4]		/* mm5 = 1st c6 */		
	movq [esp + _c6], mm5
	movd mm5, [esi + edx*4 + 4]	/* mm5 = 1st c12 */		
	movq [esp + _c12], mm5

	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm4=rsq */
	
        pfrcp mm0,mm4
        pfrcpit1 mm4,mm0				
        pfrcpit2 mm4,mm0	/* mm4=invsq */
	/* calculate potentials and scalar force */
	movq  mm0, mm4

	pfmul mm4, mm0
	pfmul mm4, mm0             	/* mm4=rinvsix */
	movq  mm5, mm4	
	pfmul mm5, mm5	                /* mm5=rinvtwelve */

	pfmul mm5, [esp + _c12]
	pfmul mm4, [esp + _c6]	
	movq mm6, mm5	/* mm6 is vnb12-vnb6 */
	pfsub mm6, mm4

	pfmul mm4, [esp + _six]

	pfmul mm5, [esp + _twelve]
	pfsub mm5, mm4
 	pfmul mm0, mm5        /* mm0 is total fscal now */

	/* update vnbtot */
	pfadd mm6, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm6       /* store the sum */  

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i0100_updateouterdata:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vnbtot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vnb[gid] */

	/* finish if last */
	mov   ecx, [ebp + _nri]
	dec ecx
	jecxz .i0100_end
	/* not last, iterate once more! */
	mov [ebp + _nri], ecx
	jmp .i0100_outer
.i0100_end:
	femms
	add esp, 108
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret
	



		
		
.globl inl0110_3dnow
	.type inl0110_3dnow,@function
inl0110_3dnow:	 
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_type,		48
.equ		_ntype,		52
.equ		_nbfp,		56	
.equ		_Vnb,		60				
.equ		_nsatoms,	64		
	/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_shX,		8
.equ		_shY,		12 
.equ		_shZ,		16	
.equ		_ix,		20
.equ		_iy,		24
.equ		_iz,		28	
.equ		_vnbtot,	32 
.equ		_c6,		40 
.equ		_c12,		48 
.equ		_six,		56 
.equ		_twelve,	64 
.equ		_ntia,		72	
.equ		_innerjjnr0,	76
.equ		_innerk0,	80		
.equ		_innerjjnr,	84
.equ		_innerk,	88	
.equ		_fix,		92
.equ		_fiy,		96
.equ		_fiz,		100
.equ		_dx1,		104
.equ		_dy1,		108
.equ		_dz1,		112
.equ		_dx2,		116
.equ		_dy2,		120
.equ		_dz2,		124					
.equ		_nsvdwc,	128
.equ		_nscoul,	132
.equ		_nsvdw,		136
.equ		_solnr,		140		
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 144		/* local stack space */
	femms
	movq  mm0, [mm_six]
	movq  mm1, [mm_twelve]
	movq  [esp + _six],    mm0
	movq  [esp + _twelve], mm1
	/* assume we have at least one i particle - start directly */		
.i0110_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm0, [eax + ebx*4]	/* move shX/shY to mm0 and shZ to mm1 */
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + _shX], mm0
	movd  [esp + _shZ], mm1

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   eax, [ebp + _nsatoms]
	add   [ebp + _nsatoms],  12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + _nsvdwc], edx
	mov   [esp + _nscoul], eax
	mov   [esp + _nsvdw], ecx
		
	/* clear potential */
	pxor  mm7,mm7
	movq  [esp + _vnbtot], mm7
	mov   [esp + _solnr],  ebx
	
	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr0], eax     /* pointer to jjnr[nj0] */

	mov   [esp + _innerk0], edx        /* number of innerloop atoms */
	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]
	
	mov   ecx, [esp + _nsvdwc]
	cmp   ecx,  0
	jnz   .i0110_mno_vdwc
	jmp   .i0110_testvdw
.i0110_mno_vdwc:
	mov   ebx, [esp + _solnr]
	inc   dword ptr [esp + _solnr]

	mov   edx, [ebp + _type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + _ntype]
	shl   edx, 1
	mov   [esp + _ntia], edx	

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	mov   [esp + _ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + _shX]
	pfadd mm1, [esp + _shZ]
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	

	/* clear forces */
	pxor  mm7,mm7
	movq  [esp + _fix],   mm7
	movd  [esp + _fiz],   mm7

	mov   ecx, [esp + _innerjjnr0]
	mov   [esp + _innerjjnr], ecx
	mov   edx, [esp + _innerk0]
        sub   edx,  2
        mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i0110_unroll_vdwc_loop
	jmp   .i0110_finish_vdwc_inner
.i0110_unroll_vdwc_loop:	
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */	

	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	mov ecx, [ecx + ebx*4]           /* type [jnr2] */

	mov esi, [ebp + _nbfp]		/* base of nbfp */ 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	add ecx, [esp + _ntia]

	movq mm5, [esi + edx*4]		/* mm5 = 1st c6 / c12 */		
	movq mm7, [esi + ecx*4]		/* mm7 = 2nd c6 / c12 */	
	movq mm6,mm5			
	punpckldq mm5,mm7		/* mm5 = 1st c6 / 2nd c6 */
	punpckhdq mm6,mm7		/* mm6 = 1st c12 / 2nd c12 */
	movq [esp + _c6], mm5
	movq [esp + _c12], mm6

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrcp mm0, mm4	                 /* lookup reciprocal seed */ 
        pfrcp mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs */
                  	        	/* amd 3dnow N-R iteration to get full precision */
        pfrcpit1 mm4,mm0				
        pfrcpit2 mm4,mm0	
	/* mm4 now contains invsq,
	 * do potential and fscal
	 */
	movq  mm0, mm4
	pfmul mm4, mm0
	pfmul mm4, mm0             	/* mm4=rinvsix */
	movq  mm5, mm4	
	pfmul mm5, mm5	                /* mm5=rinvtwelve */

	pfmul mm5, [esp + _c12]
	pfmul mm4, [esp + _c6]	
	movq mm6, mm5	/* mm6 is vnb12-vnb6 */ 
	pfsub mm6, mm4

	pfmul mm4, [esp + _six]

	pfmul mm5, [esp + _twelve]
	pfsub mm5,mm4
 	pfmul mm0, mm5        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	/* update vnbtot */
	pfadd mm6, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm6       /* store the sum */      

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i0110_finish_vdwc_inner
	jmp   .i0110_unroll_vdwc_loop
.i0110_finish_vdwc_inner:	
	and [esp + _innerk],  1
	jnz  .i0110_single_vdwc_inner
	jmp  .i0110_updateouterdata_vdwc		
.i0110_single_vdwc_inner:	
	/* a single j particle iteration here - compare with the unrolled code for comments */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov esi, [ebp + _nbfp]
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	shl edx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	movd mm5, [esi + edx*4]		/* mm5 = 1st c6 */ 		
	movq [esp + _c6], mm5
	movd mm5, [esi + edx*4 + 4]	/* mm5 = 1st c12 */ 		
	movq [esp + _c12], mm5

	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm4=rsq */
	
        pfrcp mm0,mm4
        pfrcpit1 mm4,mm0				
        pfrcpit2 mm4,mm0	/* mm4=invsq */ 
	/* calculate potentials and scalar force */
	movq  mm0, mm4

	pfmul mm4, mm0
	pfmul mm4, mm0             	/* mm4=rinvsix */
	movq  mm5, mm4	
	pfmul mm5, mm5	                /* mm5=rinvtwelve */

	pfmul mm5, [esp + _c12]
	pfmul mm4, [esp + _c6]	
	movq mm6, mm5	/* mm6 is vnb12-vnb6 */ 
	pfsub mm6, mm4

	pfmul mm4, [esp + _six]

	pfmul mm5, [esp + _twelve]
	pfsub mm5, mm4
 	pfmul mm0, mm5        /* mm0 is total fscal now */

	/* update vnbtot */
	pfadd mm6, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm6       /* store the sum */      

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i0110_updateouterdata_vdwc:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	/* loop back to mno */
	dec dword ptr [esp + _nsvdwc]
	jz  .i0110_testvdw
	jmp .i0110_mno_vdwc
.i0110_testvdw:	
	mov  ebx,  [esp + _nscoul]
	add  [esp + _solnr],  ebx

	mov  ecx, [esp + _nsvdw]
	cmp  ecx,  0
	jnz  .i0110_mno_vdw
	jmp  .i0110_last_mno
.i0110_mno_vdw:
	mov   ebx,  [esp + _solnr]
	inc   dword ptr [esp + _solnr]

	mov   edx, [ebp + _type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + _ntype]
	shl   edx, 1
	mov   [esp + _ntia], edx	

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	mov   [esp + _ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + _shX]
	pfadd mm1, [esp + _shZ]
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	

	/* clear forces */
	pxor  mm7,mm7
	movq  [esp + _fix],   mm7
	movd  [esp + _fiz],   mm7

	mov   ecx, [esp + _innerjjnr0]
	mov   [esp + _innerjjnr], ecx
	mov   edx, [esp + _innerk0]
        sub   edx,  2
        mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i0110_unroll_vdw_loop
	jmp   .i0110_finish_vdw_inner
.i0110_unroll_vdw_loop:	
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */	

	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	mov ecx, [ecx + ebx*4]           /* type [jnr2] */

	mov esi, [ebp + _nbfp]		/* base of nbfp */ 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	add ecx, [esp + _ntia]

	movq mm5, [esi + edx*4]		/* mm5 = 1st c6 / c12 */		
	movq mm7, [esi + ecx*4]		/* mm7 = 2nd c6 / c12 */	
	movq mm6,mm5			
	punpckldq mm5,mm7		/* mm5 = 1st c6 / 2nd c6 */
	punpckhdq mm6,mm7		/* mm6 = 1st c12 / 2nd c12 */
	movq [esp + _c6], mm5
	movq [esp + _c12], mm6

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrcp mm0, mm4	                 /* lookup reciprocal seed */ 
        pfrcp mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs */
                  	        	/* amd 3dnow N-R iteration to get full precision */
        pfrcpit1 mm4,mm0				
        pfrcpit2 mm4,mm0	
	/* mm4 now contains invsq,
	 * do potential and fscal
	 */
	movq  mm0, mm4
	pfmul mm4, mm0
	pfmul mm4, mm0             	/* mm4=rinvsix */
	movq  mm5, mm4	
	pfmul mm5, mm5	                /* mm5=rinvtwelve */

	pfmul mm5, [esp + _c12]
	pfmul mm4, [esp + _c6]	
	movq mm6, mm5	/* mm6 is vnb12-vnb6 */ 
	pfsub mm6, mm4

	pfmul mm4, [esp + _six]

	pfmul mm5, [esp + _twelve]
	pfsub mm5,mm4
 	pfmul mm0, mm5        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	/* update vnbtot */
	pfadd mm6, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm6       /* store the sum */      

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i0110_finish_vdw_inner
	jmp   .i0110_unroll_vdw_loop
.i0110_finish_vdw_inner:	
	and [esp + _innerk],  1
	jnz  .i0110_single_vdw_inner
	jmp  .i0110_updateouterdata_vdw		
.i0110_single_vdw_inner:	
	/* a single j particle iteration here - compare with the unrolled code for comments */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov esi, [ebp + _nbfp]
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	shl edx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	movd mm5, [esi + edx*4]		/* mm5 = 1st c6 */ 		
	movq [esp + _c6], mm5
	movd mm5, [esi + edx*4 + 4]	/* mm5 = 1st c12 */ 		
	movq [esp + _c12], mm5

	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm4=rsq */
	
        pfrcp mm0,mm4
        pfrcpit1 mm4,mm0				
        pfrcpit2 mm4,mm0	/* mm4=invsq */ 
	/* calculate potentials and scalar force */
	movq  mm0, mm4

	pfmul mm4, mm0
	pfmul mm4, mm0             	/* mm4=rinvsix */
	movq  mm5, mm4	
	pfmul mm5, mm5	                /* mm5=rinvtwelve */

	pfmul mm5, [esp + _c12]
	pfmul mm4, [esp + _c6]	
	movq mm6, mm5	/* mm6 is vnb12-vnb6 */ 
	pfsub mm6, mm4

	pfmul mm4, [esp + _six]

	pfmul mm5, [esp + _twelve]
	pfsub mm5, mm4
 	pfmul mm0, mm5        /* mm0 is total fscal now */

	/* update vnbtot */
	pfadd mm6, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm6       /* store the sum */      

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i0110_updateouterdata_vdw:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	/* loop back to mno */
	dec dword ptr [esp + _nsvdw]
	jz  .i0110_last_mno
	jmp .i0110_mno_vdw
	
.i0110_last_mno:	
	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vnbtot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */
	/* finish if last */
	mov   ecx, [ebp + _nri]
	dec ecx
	jecxz .i0110_end
	/* not last, iterate once more! */
	mov [ebp + _nri], ecx
	jmp .i0110_outer
.i0110_end:
	femms
	add esp, 144
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret



.globl inl0300_3dnow
	.type inl0300_3dnow,@function
inl0300_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_type,		48
.equ		_ntype,		52
.equ		_nbfp,		56	
.equ		_Vnb,		60
.equ		_tabscale,	64
.equ		_VFtab,		68
	/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_ix,		8
.equ		_iy,		12
.equ		_iz,		16
.equ		_vnbtot,	20 
.equ		_c6,		28 
.equ		_c12,		36 
.equ		_two,		44 
.equ		_n1,		52 
.equ		_tsc,		60 
.equ		_ntia,		68
.equ		_innerjjnr,	72
.equ		_innerk,	76		
.equ		_fix,		80
.equ		_fiy,		84
.equ		_fiz,		88
.equ		_dx1,		92
.equ		_dy1,		96
.equ		_dz1,		100
.equ		_dx2,		104
.equ		_dy2,		108
.equ		_dz2,		112						
        push ebp
        mov ebp,esp
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 116		/* local stack space */
	femms
	/* move data to local stack */ 
	movq  mm0, [mm_two]
	movd  mm3, [ebp + _tabscale]
	movq  [esp + _two],    mm0
	punpckldq mm3,mm3
	movq  [esp + _tsc], mm3	
	/* assume we have at least one i particle - start directly */	
.i0300_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm0, [eax + ebx*4]	/* move shX/shY to mm0 and shZ to mm1 */
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   edx, [ebp + _type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + _ntype]
	shl   edx, 1
	mov   [esp + _ntia], edx

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	
	pfadd mm0, [eax + ebx*4]        /* ix = shX + posX (and iy too) */
	movd  mm3, [eax + ebx*4 + 8]    /* cant use direct memory add for 4 bytes (iz) */
	mov   [esp + _ii3], ebx
	pfadd mm1, mm3
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	
				
	/* clear total potential and i forces */
	pxor  mm7,mm7
	movq  [esp + _vnbtot], mm7
	movq  [esp + _fix],    mm7
	movd  [esp + _fiz],    mm7

	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */

	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]	
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr], eax     /* pointer to jjnr[nj0] */
	sub   edx,  2
	mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i0300_unroll_loop
	jmp   .i0300_finish_inner
.i0300_unroll_loop:
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */ 
	
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	mov ecx, [ecx + ebx*4]           /* type [jnr2] */

	mov esi, [ebp + _nbfp]		/* base of nbfp */ 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */ 
	add ecx, [esp + _ntia]

	movq mm5, [esi + edx*4]		/* mm5 = 1st c6 / c12 */		
	movq mm7, [esi + ecx*4]		/* mm7 = 2nd c6 / c12 */	
	movq mm6,mm5			
	punpckldq mm5,mm7		/* mm5 = 1st c6 / 2nd c6 */
	punpckhdq mm6,mm7		/* mm6 = 1st c12 / 2nd c12 */
	movq [esp + _c6], mm5
	movq [esp + _c12], mm6

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrsqrt mm0, mm4	         /* lookup inverse square root seed */
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs */
        movq mm2,mm0	        	/* amd 3dnow N-R iteration to get full precision */
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r */
	/* do potential and fscal */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */
	
	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	mov edx, [ebp + _VFtab]
	/* dispersion table */
	mov ecx, [esp + _n1]
	shl ecx, 3
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */	

	movq mm4, [esp + _c6]
	pfmul mm7, mm4	/* fijD */
	pfmul mm5, mm4	/* vnb6 */           
	movq mm3, mm7	/* add to fscal */

	/* update vnbtot to release mm5! */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* repulsion table */
	mov ecx, [esp + _n1]
	shl ecx, 3
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	movq mm6, [esp + _c12]
	pfmul mm7, mm6	/* fijR */
	pfmul mm5, mm6	/* vnb12 */
	pfadd mm3, mm7	/* total fscal fijD+ fijR */

	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm1, [esp + _tsc]
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	/* update vnbtot */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i0300_finish_inner
	jmp   .i0300_unroll_loop
.i0300_finish_inner:	
	and [esp + _innerk],  1
	jnz  .i0300_single_inner
	jmp  .i0300_updateouterdata		
.i0300_single_inner:
	/* a single j particle iteration here - compare with the unrolled code for comments */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov esi, [ebp + _nbfp]
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	shl edx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	movd mm5, [esi + edx*4]		/* mm5 = 1st c6 */ 		
	movq [esp + _c6], mm5
	movd mm5, [esi + edx*4 + 4]	/* mm5 = 1st c12 */ 		
	movq [esp + _c12], mm5

	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm0=rsq */
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	/* mm1=invsqrt */
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r */

	/* calculate potentials and scalar force */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */

	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 3
	/* dispersion table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */	

	movq mm4, [esp + _c6]
	pfmul mm7, mm4	/* fijD */
	pfmul mm5, mm4	/* vnb6 */           
	movq mm3, mm7	/* add to fscal */

	/* update vnbtot to release mm5! */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* repulsion table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	movq mm6, [esp + _c12]
	pfmul mm7, mm6	/* fijR */
	pfmul mm5, mm6	/* vnb12 */
	pfadd mm3, mm7	/* total fscal fijC+ fijD+ fijR */

	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + _tsc]
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	/* update vnbtot */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i0300_updateouterdata:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vnbtot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vnb[gid] */

	/* finish if last */
	mov   ecx, [ebp + _nri]
	dec ecx
	jecxz .i0300_end
	/* not last, iterate once more! */
	mov [ebp + _nri], ecx
	jmp .i0300_outer
.i0300_end:
	femms
	add esp, 116
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret


			
	
.globl inl0310_3dnow
	.type inl0310_3dnow,@function
inl0310_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_type,		48
.equ		_ntype,		52
.equ		_nbfp,		56	
.equ		_Vnb,		60
.equ		_tabscale,	64
.equ		_VFtab,		68
.equ		_nsatoms,	72		
	/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_shX,		8
.equ		_shY,		12 
.equ		_shZ,		16	
.equ		_ix,		20
.equ		_iy,		24
.equ		_iz,		28	
.equ		_vnbtot,	32 
.equ		_c6,		40 
.equ		_c12,		48 
.equ		_two,		56 
.equ		_n1,		64 
.equ		_tsc,		72 
.equ		_ntia,		80	
.equ		_innerjjnr0,	84
.equ		_innerk0,	88		
.equ		_innerjjnr,	92
.equ		_innerk,	96	
.equ		_fix,		100
.equ		_fiy,		104
.equ		_fiz,		108
.equ		_dx1,		112
.equ		_dy1,		116
.equ		_dz1,		120
.equ		_dx2,		124
.equ		_dy2,		128
.equ		_dz2,		132								
.equ		_nsvdwc,	136
.equ		_nscoul,	140
.equ		_nsvdw,		144
.equ		_solnr,		148		
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 152		/* local stack space */
	femms
	movq  mm0, [mm_two]
	movd  mm3, [ebp + _tabscale]
	movq  [esp + _two],    mm0
	punpckldq mm3,mm3
	movq  [esp + _tsc], mm3	
	
	/* assume we have at least one i particle - start directly */		
.i0310_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm0, [eax + ebx*4]	/* move shX/shY to mm0 and shZ to mm1 */
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + _shX], mm0
	movd  [esp + _shZ], mm1

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   eax, [ebp + _nsatoms]
	add   [ebp + _nsatoms],  12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + _nsvdwc], edx
	mov   [esp + _nscoul], eax
	mov   [esp + _nsvdw], ecx
		
	/* clear potential */
	pxor  mm7,mm7
	movq  [esp + _vnbtot], mm7
	mov   [esp + _solnr],  ebx
	
	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr0], eax     /* pointer to jjnr[nj0] */

	mov   [esp + _innerk0], edx        /* number of innerloop atoms */
	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]
	
	mov   ecx, [esp + _nsvdwc]
	cmp   ecx,  0
	jnz   .i0310_mno_vdwc
	jmp   .i0310_testvdw
.i0310_mno_vdwc:
	mov   ebx,  [esp + _solnr]
	inc   dword ptr [esp + _solnr]

	mov   edx, [ebp + _type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + _ntype]
	shl   edx, 1
	mov   [esp + _ntia], edx	

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	mov   [esp + _ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + _shX]
	pfadd mm1, [esp + _shZ]
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	

	/* clear forces */
	pxor  mm7,mm7
	movq  [esp + _fix],   mm7
	movd  [esp + _fiz],   mm7

	mov   ecx, [esp + _innerjjnr0]
	mov   [esp + _innerjjnr], ecx
	mov   edx, [esp + _innerk0]
        sub   edx,  2
        mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i0310_unroll_vdwc_loop
	jmp   .i0310_finish_vdwc_inner
.i0310_unroll_vdwc_loop:	
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */
	
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	mov ecx, [ecx + ebx*4]           /* type [jnr2] */

	mov esi, [ebp + _nbfp]		/* base of nbfp */ 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	add ecx, [esp + _ntia]

	movq mm5, [esi + edx*4]		/* mm5 = 1st c6 / c12 */		
	movq mm7, [esi + ecx*4]		/* mm7 = 2nd c6 / c12 */	
	movq mm6,mm5			
	punpckldq mm5,mm7		/* mm5 = 1st c6 / 2nd c6 */
	punpckhdq mm6,mm7		/* mm6 = 1st c12 / 2nd c12 */
	movq [esp + _c6], mm5
	movq [esp + _c12], mm6

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrsqrt mm0, mm4	         /* lookup inverse square root seed */
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs. */
        movq mm2,mm0	        	/* amd 3dnow N-R iteration to get full precision */
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r */
	/* do potential and fscal */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */
	
	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	mov edx, [ebp + _VFtab]
	/* dispersion table */
	mov ecx, [esp + _n1]
	shl ecx, 3
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */	

	movq mm4, [esp + _c6]
	pfmul mm7, mm4	/* fijD */
	pfmul mm5, mm4	/* vnb6 */           
	movq mm3, mm7	/* add to fscal */

	/* update vnbtot to release mm5! */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* repulsion table */
	mov ecx, [esp + _n1]
	shl ecx, 3
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	movq mm6, [esp + _c12]
	pfmul mm7, mm6	/* fijR */
	pfmul mm5, mm6	/* vnb12 */
	pfadd mm3, mm7	/* total fscal fijD+ fijR */

	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm1, [esp + _tsc]
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	/* update vnbtot */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
		
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i0310_finish_vdwc_inner
	jmp   .i0310_unroll_vdwc_loop
.i0310_finish_vdwc_inner:	
	and [esp + _innerk],  1
	jnz  .i0310_single_vdwc_inner
	jmp  .i0310_updateouterdata_vdwc		
.i0310_single_vdwc_inner:	
	/* a single j particle iteration here - compare with the unrolled code for comments */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov esi, [ebp + _nbfp]
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	shl edx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	movd mm5, [esi + edx*4]		/* mm5 = 1st c6 */ 		
	movq [esp + _c6], mm5
	movd mm5, [esi + edx*4 + 4]	/* mm5 = 1st c12 */ 		
	movq [esp + _c12], mm5

	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm0=rsq */
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	/* mm1=invsqrt */
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r */

	/* calculate potentials and scalar force */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */

	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 3
	/* dispersion table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */	

	movq mm4, [esp + _c6]
	pfmul mm7, mm4	/* fijD */
	pfmul mm5, mm4	/* vnb6 */           
	movq mm3, mm7	/* add to fscal */

	/* update vnbtot to release mm5! */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* repulsion table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	movq mm6, [esp + _c12]
	pfmul mm7, mm6	/* fijR */
	pfmul mm5, mm6	/* vnb12 */
	pfadd mm3, mm7	/* total fscal fijC+ fijD+ fijR */

	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + _tsc]
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	/* update vnbtot */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i0310_updateouterdata_vdwc:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	/* loop back to mno */
	dec dword ptr [esp + _nsvdwc]
	jz  .i0310_testvdw
	jmp .i0310_mno_vdwc
.i0310_testvdw:	
	mov  ebx,  [esp + _nscoul]
	add  [esp + _solnr],  ebx

	mov  ecx, [esp + _nsvdw]
	cmp  ecx,  0
	jnz  .i0310_mno_vdw
	jmp  .i0310_last_mno
.i0310_mno_vdw:
	mov   ebx,  [esp + _solnr]
	inc   dword ptr [esp + _solnr]

	mov   edx, [ebp + _type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + _ntype]
	shl   edx, 1
	mov   [esp + _ntia], edx	

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	mov   [esp + _ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + _shX]
	pfadd mm1, [esp + _shZ]
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	

	/* clear forces */
	pxor  mm7,mm7
	movq  [esp + _fix],   mm7
	movd  [esp + _fiz],   mm7

	mov   ecx, [esp + _innerjjnr0]
	mov   [esp + _innerjjnr], ecx
	mov   edx, [esp + _innerk0]
        sub   edx,  2
        mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i0310_unroll_vdw_loop
	jmp   .i0310_finish_vdw_inner
.i0310_unroll_vdw_loop:	
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */
	
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	mov ecx, [ecx + ebx*4]           /* type [jnr2] */

	mov esi, [ebp + _nbfp]		/* base of nbfp */ 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	add ecx, [esp + _ntia]

	movq mm5, [esi + edx*4]		/* mm5 = 1st c6 / c12 */		
	movq mm7, [esi + ecx*4]		/* mm7 = 2nd c6 / c12 */	
	movq mm6,mm5			
	punpckldq mm5,mm7		/* mm5 = 1st c6 / 2nd c6 */
	punpckhdq mm6,mm7		/* mm6 = 1st c12 / 2nd c12 */
	movq [esp + _c6], mm5
	movq [esp + _c12], mm6

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrsqrt mm0, mm4	         /* lookup inverse square root seed */
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs */
        movq mm2,mm0	        	/* amd 3dnow N-R iteration to get full precision */
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r */
	/* do potential and fscal */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */
	
	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	mov edx, [ebp + _VFtab]
	/* dispersion table */
	mov ecx, [esp + _n1]
	shl ecx, 3
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */	

	movq mm4, [esp + _c6]
	pfmul mm7, mm4	/* fijD */
	pfmul mm5, mm4	/* vnb6 */           
	movq mm3, mm7	/* add to fscal */

	/* update vnbtot to release mm5! */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* repulsion table */
	mov ecx, [esp + _n1]
	shl ecx, 3
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	movq mm6, [esp + _c12]
	pfmul mm7, mm6	/* fijR */
	pfmul mm5, mm6	/* vnb12 */
	pfadd mm3, mm7	/* total fscal fijD+ fijR */

	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm1, [esp + _tsc]
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	/* update vnbtot */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i0310_finish_vdw_inner
	jmp   .i0310_unroll_vdw_loop
.i0310_finish_vdw_inner:	
	and [esp + _innerk],  1
	jnz  .i0310_single_vdw_inner
	jmp  .i0310_updateouterdata_vdw		
.i0310_single_vdw_inner:	
	/* a single j particle iteration here - compare with the unrolled code for comments */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov esi, [ebp + _nbfp]
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	shl edx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	movd mm5, [esi + edx*4]		/* mm5 = 1st c6 */ 		
	movq [esp + _c6], mm5
	movd mm5, [esi + edx*4 + 4]	/* mm5 = 1st c12 */ 		
	movq [esp + _c12], mm5

	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm0=rsq */
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	/* mm1=invsqrt */
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r */

	/* calculate potentials and scalar force */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */

	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 3
	/* dispersion table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */	

	movq mm4, [esp + _c6]
	pfmul mm7, mm4	/* fijD */
	pfmul mm5, mm4	/* vnb6 */           
	movq mm3, mm7	/* add to fscal */

	/* update vnbtot to release mm5! */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* repulsion table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	movq mm6, [esp + _c12]
	pfmul mm7, mm6	/* fijR */
	pfmul mm5, mm6	/* vnb12 */
	pfadd mm3, mm7	/* total fscal fijC+ fijD+ fijR */

	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + _tsc]
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	/* update vnbtot */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i0310_updateouterdata_vdw:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	/* loop back to mno */
	dec dword ptr [esp + _nsvdw]
	jz  .i0310_last_mno
	jmp .i0310_mno_vdw
	
.i0310_last_mno:	
	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vnbtot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */
	/* finish if last */
	mov   ecx, [ebp + _nri]
	dec ecx
	jecxz .i0310_end
	/* not last, iterate once more! */
	mov [ebp + _nri], ecx
	jmp .i0310_outer
.i0310_end:
	femms
	add esp, 152
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret


.globl inl1000_3dnow
	.type inl1000_3dnow,@function
inl1000_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56			
	/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_ix,		8
.equ		_iy,		12
.equ		_iz,		16
.equ		_iq,		20		
.equ		_vctot,		28 
.equ		_innerjjnr,	36
.equ		_innerk,	40		
.equ		_fix,		44
.equ		_fiy,		48
.equ		_fiz,		52
.equ		_dx1,		56
.equ		_dy1,		60
.equ		_dz1,		64
.equ		_dx2,		68
.equ		_dy2,		72
.equ		_dz2,		76									
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 80		/* 80 bytes local stack space */
	femms
	/* assume we have at least one i particle - start directly */	
.i1000_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm0, [eax + ebx*4]	/* move shX/shY to mm0 and shZ to mm1 */
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   edx, [ebp + _charge]
	movd  mm2, [edx + ebx*4]        /* mm2=charge[ii] */
	pfmul mm2, [ebp + _facel]
	punpckldq mm2,mm2	        /* spread to both halves */
	movq  [esp + _iq], mm2	        /* iq =facel*charge[ii] */

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	
	pfadd mm0, [eax + ebx*4]        /* ix = shX + posX (and iy too) */
	movd  mm3, [eax + ebx*4 + 8]    /* cant use direct memory add for 4 bytes (iz) */
	mov   [esp + _ii3], ebx
	pfadd mm1, mm3
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	
				
	/* clear vctot and i forces */
	pxor  mm7,mm7
	movq  [esp + _vctot], mm7
	movq  [esp + _fix],   mm7
	movd  [esp + _fiz],   mm7

	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */

	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]	
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr], eax     /* pointer to jjnr[nj0] */
	sub   edx,  2
	mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i1000_unroll_loop
	jmp   .i1000_finish_inner
.i1000_unroll_loop:	
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */
	
	mov ecx, [ebp + _charge]        /* base of charge[] */
	movq mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]	         /* charge[jnr1] */
	movd mm7, [ecx + ebx*4]  	 /* charge[jnr2] */
	punpckldq mm3,mm7	         /* move charge 2 to high part of mm3 */
	pfmul mm3,mm5		         /* mm3 now has qq for both particles */

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]	

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrsqrt mm0, mm4	         /* lookup inverse square root seed */
        pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs */
        movq mm2,mm0	        	/* amd 3dnow N-R iteration to get full precision */
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	/* mm0 now contains invsq, and mm1 invsqrt
	 * do potential and fscal
	 */
	prefetchw [esp + _dx1]	/* prefetch i forces to cache */
	
	pfmul mm3,mm1		/* 6 has both vcoul */
	pfmul mm0,mm3		/* 0 has both fscal */

	/* update vctot */

	pfadd mm3, [esp + _vctot]      /* add the earlier value */ 
	movq [esp + _vctot], mm3       /* store the sum */
	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1
	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]
	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i1000_finish_inner
	jmp   .i1000_unroll_loop
.i1000_finish_inner:	
	and [esp + _innerk],  1
	jnz  .i1000_single_inner
	jmp  .i1000_updateouterdata		
.i1000_single_inner:	
	/* a single j particle iteration here - compare with the unrolled code for comments */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov ecx, [ebp + _charge]
	movd mm6, [esp + _iq]
	movd mm7, [ecx + eax*4]
	pfmul mm6, mm7	  	/* mm6=qq */
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm2, [esi + eax*4]
	movd  mm3, [esi + eax*4 + 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq  [esp + _dx1], mm0
	pfmul mm0,mm0
	movd  [esp + _dz1], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfacc mm0, mm1		/* mm0=rsq */
	
        pfrsqrt mm1,mm0
        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	movq  mm4, mm1
	pfmul mm4, mm4		/* mm4=invsq */
	/* calculate potential and scalar force */
	pfmul mm6, mm1		/* mm6=vcoul */
	pfmul mm4, mm6		/* mm4=fscalar */ 
	/* update vctot */
	movq mm5, [esp + _vctot]
	pfadd mm5, mm6
	movq [esp + _vctot], mm5
	/* spread fscalar to both positions */
	punpckldq mm4,mm4
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm0,  [esp + _dx1]
	movd mm1,  [esp + _dz1]
	pfmul mm0, mm4
	pfmul mm1, mm4
	/* update i particle force */
	movq mm2,  [esp + _fix]
	movd mm3,  [esp + _fiz]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fix], mm2
	movd [esp + _fiz], mm3
	/* update j particle force */
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax *4+ 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	/* done! */
.i1000_updateouterdata:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */
	/* finish if last */
	mov   ecx, [ebp + _nri]
	dec ecx
	jecxz .i1000_end
	/* not last, iterate once more! */
	mov [ebp + _nri], ecx
	jmp .i1000_outer
.i1000_end:
	femms
	add esp, 80
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret


.globl inl1010_3dnow
	.type inl1010_3dnow,@function
inl1010_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56
.equ		_nsatoms,	60		
	/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_shX,		8
.equ		_shY,		12 
.equ		_shZ,		16	
.equ		_ix,		20
.equ		_iy,		24
.equ		_iz,		28	
.equ		_iq,		32 
.equ		_vctot,		40 
.equ		_innerjjnr0,	48
.equ		_innerk0,	52		
.equ		_innerjjnr,	56
.equ		_innerk,	60		
.equ		_fix,		64
.equ		_fiy,		68
.equ		_fiz,		72
.equ		_dx1,		76
.equ		_dy1,		80
.equ		_dz1,		84
.equ		_dx2,		88
.equ		_dy2,		92
.equ		_dz2,		96
.equ		_nscoul,	100
.equ		_solnr,		104		
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 108		/* local stack space */
	femms
	/* assume we have at least one i particle - start directly */	
	add   [ebp + _nsatoms],  8

.i1010_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm0, [eax + ebx*4]	/* move shX/shY to mm0 and shZ to mm1 */
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + _shX], mm0
	movd  [esp + _shZ], mm1

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   eax, [ebp + _nsatoms]
	mov   ecx, [eax]
	add   [ebp + _nsatoms],  12
	mov   [esp + _nscoul], ecx
		
	/* clear potential */
	pxor  mm7,mm7
	movq  [esp + _vctot], mm7
	mov   [esp + _solnr], ebx
	
	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr0], eax     /* pointer to jjnr[nj0] */

	mov   [esp + _innerk0], edx        /* number of innerloop atoms */
	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]

	mov   ecx, [esp + _nscoul]
	cmp   ecx,  0
	jnz   .i1010_mno_coul
	jmp   .i1010_last_mno
.i1010_mno_coul:				
	mov   ebx,  [esp + _solnr]
	inc   dword ptr [esp + _solnr]
	mov   edx, [ebp + _charge]
	movd  mm2, [edx + ebx*4]        /* mm2=charge[ii] */
	pfmul mm2, [ebp + _facel]
	punpckldq mm2,mm2	        /* spread to both halves */
	movq  [esp + _iq], mm2	        /* iq =facel*charge[ii] */

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	mov   [esp + _ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + _shX]
	pfadd mm1, [esp + _shZ]
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	

	/* clear forces */
	pxor  mm7,mm7
	movq  [esp + _fix],   mm7
	movd  [esp + _fiz],   mm7

	mov   ecx, [esp + _innerjjnr0]
	mov   [esp + _innerjjnr], ecx
	mov   edx, [esp + _innerk0]
        sub   edx,  2
        mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i1010_unroll_coul_loop
	jmp   .i1010_finish_coul_inner
.i1010_unroll_coul_loop:	
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */
	
	mov ecx, [ebp + _charge]        /* base of charge[] */
	movq mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]	         /* charge[jnr1] */
	movd mm7, [ecx + ebx*4]  	 /* charge[jnr2] */
	punpckldq mm3,mm7	         /* move charge 2 to high part of mm3 */
	pfmul mm3,mm5		         /* mm3 now has qq for both particles */

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]	

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrsqrt mm0, mm4	         /* lookup inverse square root seed */
        pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs */
        movq mm2,mm0	        	/* amd 3dnow N-R iteration to get full precision */
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	/* mm0 now contains invsq, and mm1 invsqrt */
	/* do potential and fscal */
	prefetchw [esp + _dx1]	/* prefetch i forces to cache */
	
	pfmul mm3,mm1		/* 6 has both vcoul */
	pfmul mm0,mm3		/* 0 has both fscal */

	/* update vctot */

	pfadd mm3, [esp + _vctot]      /* add the earlier value */ 
	movq [esp + _vctot], mm3       /* store the sum */
	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1
	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]
	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i1010_finish_coul_inner
	jmp   .i1010_unroll_coul_loop
.i1010_finish_coul_inner:	
	and [esp + _innerk],  1
	jnz  .i1010_single_coul_inner
	jmp  .i1010_updateouterdata_coul		
.i1010_single_coul_inner:	
	/* a single j particle iteration here - compare with the unrolled code for comments */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov ecx, [ebp + _charge]
	movd mm6, [esp + _iq]
	movd mm7, [ecx + eax*4]
	pfmul mm6, mm7	  	/* mm6=qq */
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm2, [esi + eax*4]
	movd  mm3, [esi + eax*4 + 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq  [esp + _dx1], mm0
	pfmul mm0,mm0
	movd  [esp + _dz1], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfacc mm0, mm1		/* mm0=rsq */
	
        pfrsqrt mm1,mm0
        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	movq  mm4, mm1
	pfmul mm4, mm4		/* mm4=invsq */
	/* calculate potential and scalar force */
	pfmul mm6, mm1		/* mm6=vcoul */
	pfmul mm4, mm6		/* mm4=fscalar */ 
	/* update vctot */
	movq mm5, [esp + _vctot]
	pfadd mm5, mm6
	movq [esp + _vctot], mm5
	/* spread fscalar to both positions */
	punpckldq mm4,mm4
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm0,  [esp + _dx1]
	movd mm1,  [esp + _dz1]
	pfmul mm0, mm4
	pfmul mm1, mm4
	/* update i particle force */
	movq mm2,  [esp + _fix]
	movd mm3,  [esp + _fiz]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fix], mm2
	movd [esp + _fiz], mm3
	/* update j particle force */
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax *4+ 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	/* done! */
.i1010_updateouterdata_coul:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	/* loop back to mno */
	dec dword ptr [esp + _nscoul]
	jz  .i1010_last_mno
	jmp .i1010_mno_coul
.i1010_last_mno:	
	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */
	/* finish if last */
	mov   ecx, [ebp + _nri]
	dec ecx
	jecxz .i1010_end
	/* not last, iterate once more! */
	mov [ebp + _nri], ecx
	jmp .i1010_outer
.i1010_end:
	femms
	add esp, 108
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret

			
.globl inl1020_3dnow
	.type inl1020_3dnow,@function
inl1020_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56			
			/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_ixO,		8
.equ		_iyO,		12
.equ		_izO,		16	
.equ		_ixH,		20 
.equ		_iyH,		28 
.equ		_izH,		36 
.equ		_iqO,		44	
.equ		_iqH,		52		
.equ		_vctot,		60 
.equ		_innerjjnr,	68
.equ		_innerk,	72		
.equ		_fixO,		76 
.equ		_fiyO,		80
.equ		_fizO,		84
.equ		_fixH,		88
.equ		_fiyH,		96
.equ		_fizH,		104         
.equ		_dxO,		112
.equ		_dyO,		116
.equ		_dzO,		120
.equ		_dxH,		124         
.equ		_dyH,		132         
.equ		_dzH,		140         
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 148		/* local stack space */
	femms
	/* assume we have at least one i particle - start directly */	

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   edx, [ebp + _charge]
	movd  mm1, [ebp + _facel]
	movd  mm2, [edx + ebx*4]        /* mm2=charge[i.i0] */
	pfmul mm2, mm1		
	movq  [esp + _iqO], mm2	        /* iqO = facel*charge[ii] */
	
	movd  mm2, [edx + ebx*4 + 4]    /* mm2=charge[i.i0+1] */
	pfmul mm2, mm1
	punpckldq mm2,mm2	        /* spread to both halves */
	movq  [esp + _iqH], mm2	        /* iqH = facel*charge[i.i0+1] */
.i1020_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm5, [eax + ebx*4]	/* move shX/shY to mm5 and shZ to mm6 */
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        /* also expand shX,Y,Z in mm0--mm2 */
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	
	pfadd mm5, [eax + ebx*4]        /* ix = shX + posX (and iy too) */
	movd  mm7, [eax + ebx*4 + 8]    /* cant use direct memory add for 4 bytes (iz) */
	mov   [esp + _ii3], ebx	        /* (use mm7 as temp storage for iz) */
	pfadd mm6, mm7
	movq  [esp + _ixO], mm5	
	movq  [esp + _izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] /* coords of H1 in low mm3-mm5, H2 in high */
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + _ixH], mm0	
	movq [esp + _iyH], mm1	
	movq [esp + _izH], mm2	
					
	/* clear vctot and i forces */
	pxor  mm7,mm7
	movq  [esp + _vctot], mm7
	movq  [esp + _fixO],   mm7
	movd  [esp + _fizO],   mm7
	movq  [esp + _fixH],   mm7
	movq  [esp + _fiyH],   mm7
	movq  [esp + _fizH],   mm7

	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */
	mov   [esp + _innerk], edx        /* number of innerloop atoms */

	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]	
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr], eax     /* pointer to jjnr[nj0] */
.i1020_inner_loop:	
	/* a single j particle iteration here - compare with the unrolled code for comments */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */
        add   [esp + _innerjjnr],  4 /* advance pointer */
	prefetch [ecx + 16]	   /* prefetch data - trial and error says 16 is best */

	mov ecx, [ebp + _charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + _iqO]
	pfmul mm7, [esp + _iqH]	/* mm6=qqO, mm7=qqH */
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	movq  mm4, mm1
	pfmul mm4, mm4		/* mm4=invsq */
	/* calculate potential and scalar force */
	pfmul mm6, mm1		/* mm6=vcoul */
	pfmul mm4, mm6		/* mm4=fscalar */ 

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	/* seeds are in mm5 now, and rsq in mm3 */

	movq mm2, mm5
	pfmul mm5,mm5
        pfrsqit1 mm5,mm3				
        pfrcpit2 mm5,mm2	/* mm5=invsqrt */
	movq mm3,mm5
	pfmul mm3,mm3		/* mm3=invsq */
	pfmul mm7, mm5		/* mm7=vcoul */
	pfmul mm3, mm7		/* mm3=fscal for the two H's */

	/* update vctot */
	pfadd mm7, mm6
	pfadd mm7, [esp + _vctot]
	movq [esp + _vctot], mm7
	
	/* spread oxygen fscalar to both positions */
	punpckldq mm4,mm4
	/* calc vectorial force for O */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4
	
	/* pack j forces from H in the same form as the oxygen force */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	
	/*  done  - one more? */
	dec dword ptr [esp + _innerk]
	jz  .i1020_updateouterdata
	jmp .i1020_inner_loop
.i1020_updateouterdata:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment iO force */ 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + _fixH]
	movq  mm3, [esp + _fiyH]
	movq  mm1, [esp + _fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	/* mm0(l)=fxH1, mm0(h)=fyH1 */
	punpckhdq mm2, mm3	/* mm2(l)=fxH2, mm2(h)=fyH2 */
	movq mm3, mm1
	pswapd mm3,mm3		
	/* mm1 is fzH1 */
	/* mm3 is fzH2 */
	
	movq  mm6, [edi + ecx*4 + 12]       /* increment iH1 force */ 
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       /* increment iH2 force */
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */

	/* finish if last */
	dec dword ptr [ebp + _nri]
	jz  .i1020_end
	/* not last, iterate once more! */
	jmp .i1020_outer
.i1020_end:
	femms
	add esp, 148
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret


.globl inl1030_3dnow
	.type inl1030_3dnow,@function
inl1030_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56						
			/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_ixO,		8
.equ		_iyO,		12
.equ		_izO,		16	
.equ		_ixH,		20
.equ		_iyH,		28
.equ		_izH,		36
.equ		_qqOO,		44		
.equ		_qqOH,		52		
.equ		_qqHH,		60     	
.equ		_vctot,		68
.equ		_innerjjnr,	76
.equ		_innerk,	80		
.equ		_fixO,		84 
.equ		_fiyO,		88
.equ		_fizO,		92
.equ		_fixH,		96
.equ		_fiyH,		104         
.equ		_fizH,		112         
.equ		_dxO,		120
.equ		_dyO,		124
.equ		_dzO,		128
.equ		_dxH,		132         
.equ		_dyH,		140         
.equ		_dzH,		148         
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 156		/* local stack space */
	femms
	/* assume we have at least one i particle - start directly */	

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   edx, [ebp + _charge]
	movd  mm1, [ebp + _facel]	/* mm1=facel */
	movd  mm2, [edx + ebx*4]        /* mm2=charge[i.i0] (O) */
	movd  mm3, [edx + ebx*4 + 4]    /* mm2=charge[i.i0+1] (H) */ 
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			/* mm4=qqOO*facel */
	pfmul mm5, mm3			/* mm5=qqOH*facel */
	pfmul mm6, mm3			/* mm6=qqHH*facel */
	punpckldq mm5,mm5	        /* spread to both halves */
	punpckldq mm6,mm6	        /* spread to both halves */
	movq  [esp + _qqOO], mm4
	movq  [esp + _qqOH], mm5
	movq  [esp + _qqHH], mm6
.i1030_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm5, [eax + ebx*4]	/* move shX/shY to mm5 and shZ to mm6 */
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        /* also expand shX,Y,Z in mm0--mm2 */
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */

	pfadd mm5, [eax + ebx*4]        /* ix = shX + posX (and iy too) */
	movd  mm7, [eax + ebx*4 + 8]    /* cant use direct memory add for 4 bytes (iz) */
	mov   [esp + _ii3], ebx	        /* (use mm7 as temp storage for iz) */
	pfadd mm6, mm7
	movq  [esp + _ixO], mm5	
	movq  [esp + _izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] /* coords of H1 in low mm3-mm5, H2 in high */
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + _ixH], mm0	
	movq [esp + _iyH], mm1	
	movq [esp + _izH], mm2	

	/* clear vctot and i forces */
	pxor  mm7,mm7
	movq  [esp + _vctot], mm7
	movq  [esp + _fixO],  mm7
	movq  [esp + _fizO],  mm7
	movq  [esp + _fixH],  mm7
	movq  [esp + _fiyH],  mm7
	movq  [esp + _fizH],  mm7

	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */
	mov   [esp + _innerk], edx        /* number of innerloop atoms */

	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]	
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr], eax     /* pointer to jjnr[nj0] */
.i1030_inner_loop:
	/* a single j particle iteration here - compare with the unrolled code for comments */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */
        add   [esp + _innerjjnr],  4 /* advance pointer */

	movd  mm6, [esp + _qqOO]
	movq  mm7, [esp + _qqOH]

	lea   eax, [eax + eax*2]
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	movq  mm4, mm1
	pfmul mm4, mm4		/* mm4=invsq */
	/* calculate potential and scalar force */
	pfmul mm6, mm1		/* mm6=vcoul */
	pfmul mm4, mm6		/* mm4=fscalar */ 

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	/* seeds are in mm5 now, and rsq in mm3 */

	movq mm2, mm5
	pfmul mm5,mm5
        pfrsqit1 mm5,mm3				
        pfrcpit2 mm5,mm2	/* mm5=invsqrt */
	movq mm3,mm5
	pfmul mm3,mm3		/* mm3=invsq */
	pfmul mm7, mm5		/* mm7=vcoul */
	pfmul mm3, mm7		/* mm3=fscal for the two H's */

	/* update vctot */
	pfadd mm7, mm6
	pfadd mm7, [esp + _vctot]
	movq [esp + _vctot], mm7
	
	/* spread oxygen fscalar to both positions */
	punpckldq mm4,mm4
	/* calc vectorial force for O */
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4
	
	/* pack j forces from H in the same form as the oxygen force */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3

	/* interactions with j H1 */
	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	movd mm6, [esp + _qqOH]
	movq mm7, [esp + _qqHH]
	
	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	movq  mm4, mm1
	pfmul mm4, mm4		/* mm4=invsq */
	/* calculate potential and scalar force */
	pfmul mm6, mm1		/* mm6=vcoul */
	pfmul mm4, mm6		/* mm4=fscalar */ 

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	/* seeds are in mm5 now, and rsq in mm3 */

	movq mm2, mm5
	pfmul mm5,mm5
        pfrsqit1 mm5,mm3				
        pfrcpit2 mm5,mm2	/* mm5=invsqrt */
	movq mm3,mm5
	pfmul mm3,mm3		/* mm3=invsq */
	pfmul mm7, mm5		/* mm7=vcoul */
	pfmul mm3, mm7		/* mm3=fscal for the two H's */

	/* update vctot */
	pfadd mm7, mm6
	pfadd mm7, [esp + _vctot]
	movq [esp + _vctot], mm7
	
	/* spread oxygen fscalar to both positions */
	punpckldq mm4,mm4
	/* calc vectorial force for O */
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4
	
	/* pack j forces from H in the same form as the oxygen force */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4 + 12]
	movd mm3,  [edi + eax*4 + 20]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 12], mm2
	movd [edi + eax*4 + 20], mm3

	/* interactions with j H2 */
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	movd mm6, [esp + _qqOH]
	movq mm7, [esp + _qqHH]

	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	movq  mm4, mm1
	pfmul mm4, mm4		/* mm4=invsq */
	/* calculate potential and scalar force */
	pfmul mm6, mm1		/* mm6=vcoul */
	pfmul mm4, mm6		/* mm4=fscalar */ 

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	/* seeds are in mm5 now, and rsq in mm3 */

	movq mm2, mm5
	pfmul mm5,mm5
        pfrsqit1 mm5,mm3				
        pfrcpit2 mm5,mm2	/* mm5=invsqrt */
	movq mm3,mm5
	pfmul mm3,mm3		/* mm3=invsq */
	pfmul mm7, mm5		/* mm7=vcoul */
	pfmul mm3, mm7		/* mm3=fscal for the two H's */

	/* update vctot */
	pfadd mm7, mm6
	pfadd mm7, [esp + _vctot]
	movq [esp + _vctot], mm7
	
	/* spread oxygen fscalar to both positions */
	punpckldq mm4,mm4
	/* calc vectorial force for O */
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4	

	/* pack j forces from H in the same form as the oxygen force */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4 + 24]
	movd mm3,  [edi + eax*4 + 32]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 24], mm2
	movd [edi + eax*4 + 32], mm3
	
	/*  done  - one more? */
	dec dword ptr [esp + _innerk]
	jz  .i1030_updateouterdata
	jmp .i1030_inner_loop	
.i1030_updateouterdata:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment iO force */ 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + _fixH]
	movq  mm3, [esp + _fiyH]
	movq  mm1, [esp + _fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	/* mm0(l)=fxH1, mm0(h)=fyH1 */
	punpckhdq mm2, mm3	/* mm2(l)=fxH2, mm2(h)=fyH2 */
	movq mm3, mm1
	pswapd mm3,mm3		
	/* mm1 is fzH1 */
	/* mm3 is fzH2 */

	movq  mm6, [edi + ecx*4 + 12]       /* increment iH1 force */ 
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       /* increment iH2 force */
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */

	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */
	/* finish if last */
	dec dword ptr [ebp + _nri]
	jz  .i1030_end
	/* not last, iterate once more! */
	jmp .i1030_outer
.i1030_end:
	femms
	add esp, 156
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret


.globl inl1100_3dnow
	.type inl1100_3dnow,@function
inl1100_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56			
.equ		_type,		60
.equ		_ntype,		64
.equ		_nbfp,		68	
.equ		_Vnb,		72	
	/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_ix,		8
.equ		_iy,		12
.equ		_iz,		16
.equ		_iq,		20 
.equ		_vctot,		28 
.equ		_vnbtot,	36 
.equ		_c6,		44 
.equ		_c12,		52 
.equ		_six,		60 
.equ		_twelve,	68 
.equ		_ntia,		76
.equ		_innerjjnr,	80
.equ		_innerk,	84		
.equ		_fix,		88
.equ		_fiy,		92
.equ		_fiz,		96
.equ		_dx1,		100
.equ		_dy1,		104
.equ		_dz1,		108
.equ		_dx2,		112
.equ		_dy2,		116
.equ		_dz2,		120						
	push ebp
	mov ebp,esp
	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 124		/* local stack space */
	femms
	/* move data to local stack */ 
	movq  mm0, [mm_six]
	movq  mm1, [mm_twelve]
	movq  [esp + _six],    mm0
	movq  [esp + _twelve], mm1
	/* assume we have at least one i particle - start directly */	
.i1100_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm0, [eax + ebx*4]	/* move shX/shY to mm0 and shZ to mm1 */
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   edx, [ebp + _charge]
	movd  mm2, [edx + ebx*4]        /* mm2=charge[ii] */
	pfmul mm2, [ebp + _facel]
	punpckldq mm2,mm2	        /* spread to both halves */
	movq  [esp + _iq], mm2	        /* iq =facel*charge[ii] */

	mov   edx, [ebp + _type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + _ntype]
	shl   edx, 1
	mov   [esp + _ntia], edx

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	
	pfadd mm0, [eax + ebx*4]        /* ix = shX + posX (and iy too) */
	movd  mm3, [eax + ebx*4 + 8]    /* cant use direct memory add for 4 bytes (iz) */
	mov   [esp + _ii3], ebx
	pfadd mm1, mm3
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	
				
	/* clear total potential and i forces */
	pxor  mm7,mm7
	movq  [esp + _vctot],  mm7
	movq  [esp + _vnbtot], mm7
	movq  [esp + _fix],    mm7
	movd  [esp + _fiz],    mm7

	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */

	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]	
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr], eax     /* pointer to jjnr[nj0] */
	sub   edx,  2
	mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i1100_unroll_loop
	jmp   .i1100_finish_inner
.i1100_unroll_loop:
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */
	
	mov ecx, [ebp + _charge]        /* base of charge[] */
	movq mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]	         /* charge[jnr1] */
        punpckldq mm3, [ecx + ebx*4]     /* move charge 2 to high part of mm3 */
	pfmul mm3,mm5		         /* mm3 now has qq for both particles */

	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	mov ecx, [ecx + ebx*4]           /* type [jnr2] */

	mov esi, [ebp + _nbfp]		/* base of nbfp */ 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	add ecx, [esp + _ntia]

	movq mm5, [esi + edx*4]		/* mm5 = 1st c6 / c12 */		
	movq mm7, [esi + ecx*4]		/* mm7 = 2nd c6 / c12 */	
	movq mm6,mm5			
	punpckldq mm5,mm7		/* mm5 = 1st c6 / 2nd c6 */
	punpckhdq mm6,mm7		/* mm6 = 1st c12 / 2nd c12 */
	movq [esp + _c6], mm5
	movq [esp + _c12], mm6

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrsqrt mm0, mm4	         /* lookup inverse square root seed */
        pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs */
        movq mm2,mm0	        	/* amd 3dnow N-R iteration to get full precision */
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	/* mm0 now contains invsq, and mm1 invsqrt */
	/* do potential and fscal */
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	/* mm4=rinvsix */
	movq  mm5, mm4	
	pfmul mm5, mm5	                /* mm5=rinvtwelve */

	pfmul mm3, mm1		/* mm3 has vcoul for both interactions */
	movq  mm7, mm3	        /* use mm7 for sum to make fscal */ 

	pfmul mm5, [esp + _c12]
	pfmul mm4, [esp + _c6]	
	movq mm6, mm5	/* mm6 is vnb12-vnb6 */ 
	pfsub mm6, mm4

	pfmul mm4, [esp + _six]

	pfmul mm5, [esp + _twelve]
	pfsub mm7,mm4
	pfadd mm7, mm5
 	pfmul mm0, mm7        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* update vctot */
	pfadd mm3, [esp + _vctot]      /* add the earlier value */
	movq [esp + _vctot], mm3       /* store the sum */      

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	/* update vnbtot */
	pfadd mm6, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm6       /* store the sum */      

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i1100_finish_inner
	jmp   .i1100_unroll_loop
.i1100_finish_inner:	
	and [esp + _innerk],  1
	jnz  .i1100_single_inner
	jmp  .i1100_updateouterdata		
.i1100_single_inner:
	/* a single j particle iteration here - compare with the unrolled code for comments */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov ecx, [ebp + _charge]
	movd mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	/* mm3=qq */

	mov esi, [ebp + _nbfp]
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	shl edx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	movd mm5, [esi + edx*4]		/* mm5 = 1st c6 */ 		
	movq [esp + _c6], mm5
	movd mm5, [esi + edx*4 + 4]	/* mm5 = 1st c12 */ 		
	movq [esp + _c12], mm5


	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm0=rsq */
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	/* mm1=invsqrt */
	movq  mm1, mm0
	pfmul mm0, mm0		/* mm0=invsq */
	/* calculate potentials and scalar force */
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	/* mm4=rinvsix */
	movq  mm5, mm4	
	pfmul mm5, mm5	                /* mm5=rinvtwelve */

	pfmul mm3, mm1		/* mm3 has vcoul for both interactions */
	movq  mm7, mm3	        /* use mm7 for sum to make fscal */ 

	pfmul mm5, [esp + _c12]
	pfmul mm4, [esp + _c6]	
	movq mm6, mm5	/* mm6 is vnb12-vnb6 */ 
	pfsub mm6, mm4

	pfmul mm4, [esp + _six]

	pfmul mm5, [esp + _twelve]
	pfsub mm7,mm4
	pfadd mm7, mm5
 	pfmul mm0, mm7        /* mm0 is total fscal now */

	/* update vctot */
	pfadd mm3, [esp + _vctot]
	movq [esp + _vctot], mm3

	/* update vnbtot */
	pfadd mm6, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm6       /* store the sum */      

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i1100_updateouterdata:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */

	movq  mm7, [esp + _vnbtot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vnb[gid] */

	/* finish if last */
	mov   ecx, [ebp + _nri]
	dec ecx
	jecxz .i1100_end
	/* not last, iterate once more! */
	mov [ebp + _nri], ecx
	jmp .i1100_outer
.i1100_end:
	femms
	add esp, 124
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret

	



.globl inl1110_3dnow
	.type inl1110_3dnow,@function
inl1110_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56	
.equ		_type,		60
.equ		_ntype,		64
.equ		_nbfp,		68	
.equ		_Vnb,		72				
.equ		_nsatoms,	76		
	/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_shX,		8
.equ		_shY,		12 
.equ		_shZ,		16	
.equ		_ix,		20
.equ		_iy,		24
.equ		_iz,		28	
.equ		_iq,		32		 
.equ		_vctot,		40 
.equ		_vnbtot,	48 
.equ		_c6,		56 
.equ		_c12,		64 
.equ		_six,		72 
.equ		_twelve,	80 
.equ		_ntia,		88	
.equ		_innerjjnr0,	92
.equ		_innerk0,	96		
.equ		_innerjjnr,	100
.equ		_innerk,	104	
.equ		_fix,		108
.equ		_fiy,		112
.equ		_fiz,		116
.equ		_dx1,		120
.equ		_dy1,		124
.equ		_dz1,		128
.equ		_dx2,		132
.equ		_dy2,		136
.equ		_dz2,		140								
.equ		_nsvdwc,	144
.equ		_nscoul,	148
.equ		_nsvdw,		152
.equ		_solnr,		156		
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 160		/* local stack space */
	femms
	movq  mm0, [mm_six]
	movq  mm1, [mm_twelve]
	movq  [esp + _six],    mm0
	movq  [esp + _twelve], mm1
	/* assume we have at least one i particle - start directly */		
.i1110_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm0, [eax + ebx*4]	/* move shX/shY to mm0 and shZ to mm1 */
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + _shX], mm0
	movd  [esp + _shZ], mm1

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   eax, [ebp + _nsatoms]
	add   [ebp + _nsatoms],  12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + _nsvdwc], edx
	mov   [esp + _nscoul], eax
	mov   [esp + _nsvdw], ecx
		
	/* clear potential */
	pxor  mm7,mm7
	movq  [esp + _vctot],  mm7
	movq  [esp + _vnbtot], mm7
	mov   [esp + _solnr],  ebx
	
	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr0], eax     /* pointer to jjnr[nj0] */

	mov   [esp + _innerk0], edx        /* number of innerloop atoms */
	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]
	
	mov   ecx, [esp + _nsvdwc]
	cmp   ecx,  0
	jnz   .i1110_mno_vdwc
	jmp   .i1110_testcoul
.i1110_mno_vdwc:
	mov   ebx,  [esp + _solnr]
	inc   dword ptr [esp + _solnr]
	mov   edx, [ebp + _charge]
	movd  mm2, [edx + ebx*4]        /* mm2=charge[ii] */
	pfmul mm2, [ebp + _facel]
	punpckldq mm2,mm2	        /* spread to both halves */
	movq  [esp + _iq], mm2	        /* iq =facel*charge[ii] */

	mov   edx, [ebp + _type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + _ntype]
	shl   edx, 1
	mov   [esp + _ntia], edx	

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	mov   [esp + _ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + _shX]
	pfadd mm1, [esp + _shZ]
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	

	/* clear forces */
	pxor  mm7,mm7
	movq  [esp + _fix],   mm7
	movd  [esp + _fiz],   mm7

	mov   ecx, [esp + _innerjjnr0]
	mov   [esp + _innerjjnr], ecx
	mov   edx, [esp + _innerk0]
        sub   edx,  2
        mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i1110_unroll_vdwc_loop
	jmp   .i1110_finish_vdwc_inner
.i1110_unroll_vdwc_loop:	
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */
	
	mov ecx, [ebp + _charge]        /* base of charge[] */
	movq mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]	         /* charge[jnr1] */
        punpckldq mm3, [ecx + ebx*4]     /* move charge 2 to high part of mm3 */
	pfmul mm3,mm5		         /* mm3 now has qq for both particles */

	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	mov ecx, [ecx + ebx*4]           /* type [jnr2] */

	mov esi, [ebp + _nbfp]		/* base of nbfp */ 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	add ecx, [esp + _ntia]

	movq mm5, [esi + edx*4]		/* mm5 = 1st c6 / c12 */		
	movq mm7, [esi + ecx*4]		/* mm7 = 2nd c6 / c12 */	
	movq mm6,mm5			
	punpckldq mm5,mm7		/* mm5 = 1st c6 / 2nd c6 */
	punpckhdq mm6,mm7		/* mm6 = 1st c12 / 2nd c12 */
	movq [esp + _c6], mm5
	movq [esp + _c12], mm6

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrsqrt mm0, mm4	         /* lookup inverse square root seed */
        pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs */
        movq mm2,mm0	        	/* amd 3dnow N-R iteration to get full precision */
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	/* mm0 now contains invsq, and mm1 invsqrt */
	/* do potential and fscal */
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	/* mm4=rinvsix */
	movq  mm5, mm4	
	pfmul mm5, mm5	                /* mm5=rinvtwelve */

	pfmul mm3, mm1		/* mm3 has vcoul for both interactions */
	movq  mm7, mm3	        /* use mm7 for sum to make fscal */ 

	pfmul mm5, [esp + _c12]
	pfmul mm4, [esp + _c6]	
	movq mm6, mm5	/* mm6 is vnb12-vnb6 */ 
	pfsub mm6, mm4

	pfmul mm4, [esp + _six]

	pfmul mm5, [esp + _twelve]
	pfsub mm7,mm4
	pfadd mm7, mm5
 	pfmul mm0, mm7        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* update vctot */
	pfadd mm3, [esp + _vctot]      /* add the earlier value */
	movq [esp + _vctot], mm3       /* store the sum */      

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	/* update vnbtot */
	pfadd mm6, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm6       /* store the sum */      

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i1110_finish_vdwc_inner
	jmp   .i1110_unroll_vdwc_loop
.i1110_finish_vdwc_inner:	
	and [esp + _innerk],  1
	jnz  .i1110_single_vdwc_inner
	jmp  .i1110_updateouterdata_vdwc		
.i1110_single_vdwc_inner:	
	/* a single j particle iteration here - compare with the unrolled code for comments */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov ecx, [ebp + _charge]
	movd mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	/* mm3=qq */

	mov esi, [ebp + _nbfp]
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	shl edx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	movd mm5, [esi + edx*4]		/* mm5 = 1st c6 */ 		
	movq [esp + _c6], mm5
	movd mm5, [esi + edx*4 + 4]	/* mm5 = 1st c12 */ 		
	movq [esp + _c12], mm5


	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm0=rsq */
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	/* mm1=invsqrt */
	movq  mm1, mm0
	pfmul mm0, mm0		/* mm0=invsq */
	/* calculate potentials and scalar force */
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	/* mm4=rinvsix */
	movq  mm5, mm4	
	pfmul mm5, mm5	                /* mm5=rinvtwelve */

	pfmul mm3, mm1		/* mm3 has vcoul for both interactions */
	movq  mm7, mm3	        /* use mm7 for sum to make fscal */ 

	pfmul mm5, [esp + _c12]
	pfmul mm4, [esp + _c6]	
	movq mm6, mm5	/* mm6 is vnb12-vnb6 */ 
	pfsub mm6, mm4

	pfmul mm4, [esp + _six]

	pfmul mm5, [esp + _twelve]
	pfsub mm7,mm4
	pfadd mm7, mm5
 	pfmul mm0, mm7        /* mm0 is total fscal now */

	/* update vctot */
	pfadd mm3, [esp + _vctot]
	movq [esp + _vctot], mm3

	/* update vnbtot */
	pfadd mm6, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm6       /* store the sum */      

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i1110_updateouterdata_vdwc:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	/* loop back to mno */
	dec  dword ptr [esp + _nsvdwc]
	jz  .i1110_testcoul
	jmp .i1110_mno_vdwc
.i1110_testcoul:	
	mov  ecx, [esp + _nscoul]
	cmp  ecx,  0
	jnz  .i1110_mno_coul
	jmp  .i1110_testvdw
.i1110_mno_coul:
	mov   ebx,  [esp + _solnr]
	inc   dword ptr [esp + _solnr]
	mov   edx, [ebp + _charge]
	movd  mm2, [edx + ebx*4]        /* mm2=charge[ii] */
	pfmul mm2, [ebp + _facel]
	punpckldq mm2,mm2	        /* spread to both halves */
	movq  [esp + _iq], mm2	        /* iq =facel*charge[ii] */

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	mov   [esp + _ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + _shX]
	pfadd mm1, [esp + _shZ]
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	

	/* clear forces */
	pxor  mm7,mm7
	movq  [esp + _fix],   mm7
	movd  [esp + _fiz],   mm7

	mov   ecx, [esp + _innerjjnr0]
	mov   [esp + _innerjjnr], ecx
	mov   edx, [esp + _innerk0]
        sub   edx,  2
        mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i1110_unroll_coul_loop
	jmp   .i1110_finish_coul_inner
.i1110_unroll_coul_loop:	
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */
	
	mov ecx, [ebp + _charge]        /* base of charge[] */
	movq mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]	         /* charge[jnr1] */
	movd mm7, [ecx + ebx*4]  	 /* charge[jnr2] */
	punpckldq mm3,mm7	         /* move charge 2 to high part of mm3 */
	pfmul mm3,mm5		         /* mm3 now has qq for both particles */

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]	

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrsqrt mm0, mm4	         /* lookup inverse square root seed */
        pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs */
        movq mm2,mm0	        	/* amd 3dnow N-R iteration to get full precision */
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	/* mm0 now contains invsq, and mm1 invsqrt */
	/* do potential and fscal */
	prefetchw [esp + _dx1]	/* prefetch i forces to cache */
	
	pfmul mm3,mm1		/* 6 has both vcoul */
	pfmul mm0,mm3		/* 0 has both fscal */

	/* update vctot */

	pfadd mm3, [esp + _vctot]      /* add the earlier value */ 
	movq [esp + _vctot], mm3       /* store the sum */
	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1
	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]
	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i1110_finish_coul_inner
	jmp   .i1110_unroll_coul_loop
.i1110_finish_coul_inner:	
	and [esp + _innerk],  1
	jnz  .i1110_single_coul_inner
	jmp  .i1110_updateouterdata_coul		
.i1110_single_coul_inner:	
	/* a single j particle iteration here - compare with the unrolled code for comments */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov ecx, [ebp + _charge]
	movd mm6, [esp + _iq]
	movd mm7, [ecx + eax*4]
	pfmul mm6, mm7	  	/* mm6=qq */
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm2, [esi + eax*4]
	movd  mm3, [esi + eax*4 + 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq  [esp + _dx1], mm0
	pfmul mm0,mm0
	movd  [esp + _dz1], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfacc mm0, mm1		/* mm0=rsq */
	
        pfrsqrt mm1,mm0
        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	movq  mm4, mm1
	pfmul mm4, mm4		/* mm4=invsq */
	/* calculate potential and scalar force */
	pfmul mm6, mm1		/* mm6=vcoul */
	pfmul mm4, mm6		/* mm4=fscalar */ 
	/* update vctot */
	movq mm5, [esp + _vctot]
	pfadd mm5, mm6
	movq [esp + _vctot], mm5
	/* spread fscalar to both positions */
	punpckldq mm4,mm4
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm0,  [esp + _dx1]
	movd mm1,  [esp + _dz1]
	pfmul mm0, mm4
	pfmul mm1, mm4
	/* update i particle force */
	movq mm2,  [esp + _fix]
	movd mm3,  [esp + _fiz]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fix], mm2
	movd [esp + _fiz], mm3
	/* update j particle force */
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax *4+ 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	/* done! */
.i1110_updateouterdata_coul:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	/* loop back to mno */
	dec dword ptr [esp + _nscoul]
	jz  .i1110_testvdw
	jmp .i1110_mno_coul
.i1110_testvdw:	
	mov  ecx, [esp + _nsvdw]
	cmp  ecx,  0
	jnz  .i1110_mno_vdw
	jmp  .i1110_last_mno
.i1110_mno_vdw:
	mov   ebx,  [esp + _solnr]
	inc   dword ptr [esp + _solnr]

	mov   edx, [ebp + _type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + _ntype]
	shl   edx, 1
	mov   [esp + _ntia], edx	

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	mov   [esp + _ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + _shX]
	pfadd mm1, [esp + _shZ]
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	

	/* clear forces */
	pxor  mm7,mm7
	movq  [esp + _fix],   mm7
	movd  [esp + _fiz],   mm7

	mov   ecx, [esp + _innerjjnr0]
	mov   [esp + _innerjjnr], ecx
	mov   edx, [esp + _innerk0]
        sub   edx,  2
        mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i1110_unroll_vdw_loop
	jmp   .i1110_finish_vdw_inner
.i1110_unroll_vdw_loop:	
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */	

	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	mov ecx, [ecx + ebx*4]           /* type [jnr2] */

	mov esi, [ebp + _nbfp]		/* base of nbfp */ 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	add ecx, [esp + _ntia]

	movq mm5, [esi + edx*4]		/* mm5 = 1st c6 / c12 */		
	movq mm7, [esi + ecx*4]		/* mm7 = 2nd c6 / c12 */	
	movq mm6,mm5			
	punpckldq mm5,mm7		/* mm5 = 1st c6 / 2nd c6 */
	punpckhdq mm6,mm7		/* mm6 = 1st c12 / 2nd c12 */
	movq [esp + _c6], mm5
	movq [esp + _c12], mm6

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrsqrt mm0, mm4	         /* lookup inverse square root seed */
        pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs */
        movq mm2,mm0	        	/* amd 3dnow N-R iteration to get full precision */
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	/* mm0 now contains invsq, and mm1 invsqrt */
	/* do potential and fscal */
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	/* mm4=rinvsix */
	movq  mm5, mm4	
	pfmul mm5, mm5	                /* mm5=rinvtwelve */

	pfmul mm5, [esp + _c12]
	pfmul mm4, [esp + _c6]	
	movq mm6, mm5	/* mm6 is vnb12-vnb6 */ 
	pfsub mm6, mm4

	pfmul mm4, [esp + _six]

	pfmul mm5, [esp + _twelve]
	movq  mm7, mm5
	pfsub mm7,mm4
 	pfmul mm0, mm7        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	/* update vnbtot */
	pfadd mm6, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm6       /* store the sum */      

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i1110_finish_vdw_inner
	jmp   .i1110_unroll_vdw_loop
.i1110_finish_vdw_inner:	
	and [esp + _innerk],  1
	jnz  .i1110_single_vdw_inner
	jmp  .i1110_updateouterdata_vdw		
.i1110_single_vdw_inner:	
	/* a single j particle iteration here - compare with the unrolled code for comments */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov esi, [ebp + _nbfp]
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	shl edx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	movd mm5, [esi + edx*4]		/* mm5 = 1st c6 */ 		
	movq [esp + _c6], mm5
	movd mm5, [esi + edx*4 + 4]	/* mm5 = 1st c12 */ 		
	movq [esp + _c12], mm5


	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm0=rsq */
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	/* mm1=invsqrt */
	movq  mm1, mm0
	pfmul mm0, mm0		/* mm0=invsq */
	/* calculate potentials and scalar force */
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	/* mm4=rinvsix */
	movq  mm5, mm4	
	pfmul mm5, mm5	                /* mm5=rinvtwelve */

	pfmul mm5, [esp + _c12]
	pfmul mm4, [esp + _c6]	
	movq mm6, mm5	/* mm6 is vnb12-vnb6 */ 
	pfsub mm6, mm4

	pfmul mm4, [esp + _six]

	pfmul mm5, [esp + _twelve]
	movq  mm7, mm5
	pfsub mm7,mm4
 	pfmul mm0, mm7        /* mm0 is total fscal now */

	/* update vnbtot */
	pfadd mm6, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm6       /* store the sum */      

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i1110_updateouterdata_vdw:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	/* loop back to mno */
	dec dword ptr [esp + _nsvdw]
	jz  .i1110_last_mno
	jmp .i1110_mno_vdw
	
.i1110_last_mno:	
	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */

	movq  mm7, [esp + _vnbtot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */
	/* finish if last */
	mov   ecx, [ebp + _nri]
	dec ecx
	jecxz .i1110_end
	/* not last, iterate once more! */
	mov [ebp + _nri], ecx
	jmp .i1110_outer
.i1110_end:
	femms
	add esp, 160
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret



.globl inl1120_3dnow
	.type inl1120_3dnow,@function
inl1120_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56	
.equ		_type,		60
.equ		_ntype,		64
.equ		_nbfp,		68	
.equ		_Vnb,		72				
			/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_ixO,		8
.equ		_iyO,		12
.equ		_izO,		16	
.equ		_ixH,		20  
.equ		_iyH,		28  
.equ		_izH,		36  
.equ		_iqO,		44  
.equ		_iqH,		52  
.equ		_vctot,		60  
.equ		_vnbtot,		68  
.equ		_c6,		76  
.equ		_c12,		84  
.equ		_six,		92  
.equ		_twelve,		100 
.equ		_ntia,		108 
.equ		_innerjjnr,		116
.equ		_innerk,		120	
.equ		_fixO,		124
.equ		_fiyO,		128
.equ		_fizO,		132
.equ		_fixH,		136  
.equ		_fiyH,		144  
.equ		_fizH,		152  
.equ		_dxO,		160
.equ		_dyO,		164
.equ		_dzO,		168
.equ		_dxH,		172  
.equ		_dyH,		180  
.equ		_dzH,		188  
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 196		/* local stack space */
	femms
	/* assume we have at least one i particle - start directly */	

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   edx, [ebp + _charge]
	movd  mm1, [ebp + _facel]
	movd  mm2, [edx + ebx*4]        /* mm2=charge[i.i0] */
	pfmul mm2, mm1		
	movq  [esp + _iqO], mm2	        /* iqO = facel*charge[ii] */
	
	movd  mm2, [edx + ebx*4 + 4]    /* mm2=charge[i.i0+1] */
	pfmul mm2, mm1
	punpckldq mm2,mm2	        /* spread to both halves */
	movq  [esp + _iqH], mm2	        /* iqH = facel*charge[i.i0+1] */

	mov   edx, [ebp + _type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	imul  ecx, [ebp + _ntype]      /* ecx = ntia = 2*ntype*type[i.i0] */ 
	mov   [esp + _ntia], ecx
	
	movq  mm3, [mm_six]
	movq  mm4, [mm_twelve]
	movq  [esp + _six],    mm3
	movq  [esp + _twelve], mm4  
.i1120_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm5, [eax + ebx*4]	/* move shX/shY to mm5 and shZ to mm6. */
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        /* also expand shX,Y,Z in mm0--mm2. */
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	
	pfadd mm5, [eax + ebx*4]        /* ix = shX + posX (and iy too) */
	movd  mm7, [eax + ebx*4 + 8]    /* cant use direct memory add for 4 bytes (iz) */
	mov   [esp + _ii3], ebx	        /* (use mm7 as temp. storage for iz.) */
	pfadd mm6, mm7
	movq  [esp + _ixO], mm5	
	movq  [esp + _izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] /* coords of H1 in low mm3-mm5, H2 in high */
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + _ixH], mm0	
	movq [esp + _iyH], mm1	
	movq [esp + _izH], mm2	
					
	/* clear vctot and i forces */
	pxor  mm7,mm7
	movq  [esp + _vctot], mm7
	movq  [esp + _vnbtot], mm7
	movq  [esp + _fixO],   mm7
	movd  [esp + _fizO],   mm7
	movq  [esp + _fixH],   mm7
	movq  [esp + _fiyH],   mm7
	movq  [esp + _fizH],   mm7

	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */
	mov   [esp + _innerk], edx        /* number of innerloop atoms */

	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]	
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr], eax     /* pointer to jjnr[nj0] */
.i1120_inner_loop:	
	/* a single j particle iteration here - compare with the unrolled code for comments. */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */
        add   [esp + _innerjjnr],  4 /* advance pointer */
	prefetch [ecx + 16]	   /* prefetch data - trial and error says 16 is best */

	mov ecx, [ebp + _charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + _iqO]
	pfmul mm7, [esp + _iqH]	/* mm6=qqO, mm7=qqH */

	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr] */
	mov ecx, [ebp + _nbfp]
	shl edx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	movd mm5, [ecx + edx*4]		/* mm5 = 1st c6 */ 		
	movq [esp + _c6], mm5
	movd mm5, [ecx + edx*4 + 4]	/* mm5 = 1st c12 */ 		
	movq [esp + _c12], mm5	
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	movq  mm4, mm1
	pfmul mm4, mm4		/* mm4=invsq */

	movq  mm0, mm4
	pfmul mm0, mm4
	pfmul mm0, mm4		/* mm0=rinvsix */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm2=rintwelve */
	
	/* calculate potential and scalar force */
	pfmul mm6, mm1		/* mm6=vcoul */
	movq  mm1, mm6		/* use mm1 for fscal sum */

	/* LJ for the oxygen */
	pfmul mm0, [esp + _c6]	 
	pfmul mm2, [esp + _c12]	 

	/* calc nb potential */
	movq mm5, mm2
	pfsub mm5, mm0

	/* calc nb force */
	pfmul mm0, [esp + _six]
	pfmul mm2, [esp + _twelve]
	
	/* increment scalar force */
	pfsub mm1, mm0
	pfadd mm1, mm2
	pfmul mm4, mm1		/* total scalar force on oxygen. */
	
	/* update nb potential */
	pfadd mm5, [esp + _vnbtot]
	movq [esp + _vnbtot], mm5
	
	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	/* seeds are in mm5 now, and rsq in mm3. */

	movq mm2, mm5
	pfmul mm5,mm5
        pfrsqit1 mm5,mm3				
        pfrcpit2 mm5,mm2	/* mm5=invsqrt */
	movq mm3,mm5
	pfmul mm3,mm3		/* mm3=invsq */
	pfmul mm7, mm5		/* mm7=vcoul */
	pfmul mm3, mm7		/* mm3=fscal for the two H's. */

	/* update vctot */
	pfadd mm7, mm6
	pfadd mm7, [esp + _vctot]
	movq [esp + _vctot], mm7
	
	/* spread oxygen fscalar to both positions */
	punpckldq mm4,mm4
	/* calc vectorial force for O */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4
	
	/* pack j forces from H in the same form as the oxygen force. */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle. */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	
	/*  done  - one more? */
	dec dword ptr [esp + _innerk]
	jz  .i1120_updateouterdata
	jmp .i1120_inner_loop
.i1120_updateouterdata:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment iO force */ 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + _fixH]
	movq  mm3, [esp + _fiyH]
	movq  mm1, [esp + _fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	/* mm0(l)=fxH1, mm0(h)=fyH1 */
	punpckhdq mm2, mm3	/* mm2(l)=fxH2, mm2(h)=fyH2 */
	movq mm3, mm1
	pswapd mm3,mm3		
	/* mm1 is fzH1 */
	/* mm3 is fzH2 */
	
	movq  mm6, [edi + ecx*4 + 12]       /* increment iH1 force */ 
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       /* increment iH2 force */
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */

	movq  mm7, [esp + _vnbtot]     
	pfacc mm7,mm7	              /* same for Vnb */
	
	mov   eax, [ebp + _Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vnb[gid] */
	/* finish if last */
	dec dword ptr [ebp + _nri]
	jz  .i1120_end
	/* not last, iterate once more! */
	jmp .i1120_outer
.i1120_end:
	femms
	add esp, 196
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret

	

.globl inl1130_3dnow
	.type inl1130_3dnow,@function
inl1130_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56						
.equ		_type,		60
.equ		_ntype,		64
.equ		_nbfp,		68	
.equ		_Vnb,		72			
			/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_ixO,		8
.equ		_iyO,		12
.equ		_izO,		16	
.equ		_ixH,		20  
.equ		_iyH,		28  
.equ		_izH,		36  
.equ		_qqOO,		44  
.equ		_qqOH,		52  
.equ		_qqHH,		60  
.equ		_c6,		68  
.equ		_c12,		76  
.equ		_six,		84  
.equ		_twelve,	92  
.equ		_vctot,		100 
.equ		_vnbtot,	108 
.equ		_innerjjnr,	116
.equ		_innerk,	120	
.equ		_fixO,		124
.equ		_fiyO,		128
.equ		_fizO,		132
.equ		_fixH,		136 
.equ		_fiyH,		144 
.equ		_fizH,		152 
.equ		_dxO,		160
.equ		_dyO,		164
.equ		_dzO,		168
.equ		_dxH,		172 
.equ		_dyH,		180 
.equ		_dzH,		188 
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 196		/* local stack space */
	femms
	/* assume we have at least one i particle - start directly */	

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   edx, [ebp + _charge]
	movd  mm1, [ebp + _facel]	/* mm1=facel */
	movd  mm2, [edx + ebx*4]        /* mm2=charge[i.i0] (O) */
	movd  mm3, [edx + ebx*4 + 4]    /* mm2=charge[i.i0+1] (H) */ 
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			/* mm4=qqOO*facel */
	pfmul mm5, mm3			/* mm5=qqOH*facel */
	pfmul mm6, mm3			/* mm6=qqHH*facel */
	punpckldq mm5,mm5	        /* spread to both halves */
	punpckldq mm6,mm6	        /* spread to both halves */
	movq  [esp + _qqOO], mm4
	movq  [esp + _qqOH], mm5
	movq  [esp + _qqHH], mm6
	mov   edx, [ebp + _type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [ebp + _ntype]
	add   edx, ecx
	mov   eax, [ebp + _nbfp]
	movd  mm0, [eax + edx*4]          
	movd  mm1, [eax + edx*4 + 4]
	movq  [esp + _c6], mm0
	movq  [esp + _c12], mm1
	movq  mm2, [mm_six]
	movq  mm3, [mm_twelve]
	movq  [esp + _six], mm2
	movq  [esp + _twelve], mm3
.i1130_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm5, [eax + ebx*4]	/* move shX/shY to mm5 and shZ to mm6. */
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        /* also expand shX,Y,Z in mm0--mm2. */
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */

	pfadd mm5, [eax + ebx*4]        /* ix = shX + posX (and iy too) */
	movd  mm7, [eax + ebx*4 + 8]    /* cant use direct memory add for 4 bytes (iz) */
	mov   [esp + _ii3], ebx	        /* (use mm7 as temp. storage for iz.) */
	pfadd mm6, mm7
	movq  [esp + _ixO], mm5	
	movq  [esp + _izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] /* coords of H1 in low mm3-mm5, H2 in high */
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + _ixH], mm0	
	movq [esp + _iyH], mm1	
	movq [esp + _izH], mm2	

	/* clear vctot and i forces */
	pxor  mm7,mm7
	movq  [esp + _vctot], mm7
	movq  [esp + _vnbtot], mm7
	movq  [esp + _fixO],  mm7
	movq  [esp + _fizO],  mm7
	movq  [esp + _fixH],  mm7
	movq  [esp + _fiyH],  mm7
	movq  [esp + _fizH],  mm7

	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */
	mov   [esp + _innerk], edx        /* number of innerloop atoms */

	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]	
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr], eax     /* pointer to jjnr[nj0] */
.i1130_inner_loop:
	/* a single j particle iteration here - compare with the unrolled code for comments. */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */
        add   [esp + _innerjjnr],  4 /* advance pointer */

	movd  mm6, [esp + _qqOO]
	movq  mm7, [esp + _qqOH]

	lea   eax, [eax + eax*2]
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */ OO
	movq  mm4, mm1
	pfmul mm4, mm4		/* mm4=invsq */ OO

	movq mm2, mm4
	pfmul mm2, mm4
	pfmul mm2, mm4
	movq mm0, mm2
	pfmul mm0,mm0
	pfmul mm2, [esp + _c6]
	pfmul mm0, [esp + _c12]
	movq mm5, mm0
	pfsub mm5, mm2		/* vnb */

	pfmul mm2, [esp + _six]
	pfmul mm0, [esp + _twelve]

	pfsub mm0, mm2
	
	/* calculate potential and scalar force */
	pfmul mm6, mm1		/* mm6=vcoul */
	pfadd mm0, mm6
	pfmul mm4, mm0		/* mm4=fscalar */ 

	/* update nb potential */
	pfadd mm5, [esp + _vnbtot]
	movq [esp + _vnbtot], mm5

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	/* seeds are in mm5 now, and rsq in mm3 */

	movq mm2, mm5
	pfmul mm5,mm5
        pfrsqit1 mm5,mm3				
        pfrcpit2 mm5,mm2	/* mm5=invsqrt */
	movq mm3,mm5
	pfmul mm3,mm3		/* mm3=invsq */
	pfmul mm7, mm5		/* mm7=vcoul */
	pfmul mm3, mm7		/* mm3=fscal for the two H's. */

	/* update vctot */
	pfadd mm7, mm6
	pfadd mm7, [esp + _vctot]
	movq [esp + _vctot], mm7
	
	/* spread oxygen fscalar to both positions */
	punpckldq mm4,mm4
	/* calc vectorial force for O */
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4
	
	/* pack j forces from H in the same form as the oxygen force. */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle. */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3

	/* interactions with j H1 */
	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	movd mm6, [esp + _qqOH]
	movq mm7, [esp + _qqHH]
	
	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	movq  mm4, mm1
	pfmul mm4, mm4		/* mm4=invsq */
	/* calculate potential and scalar force */
	pfmul mm6, mm1		/* mm6=vcoul */
	pfmul mm4, mm6		/* mm4=fscalar */ 

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	/* seeds are in mm5 now, and rsq in mm3 */

	movq mm2, mm5
	pfmul mm5,mm5
        pfrsqit1 mm5,mm3				
        pfrcpit2 mm5,mm2	/* mm5=invsqrt */
	movq mm3,mm5
	pfmul mm3,mm3		/* mm3=invsq */
	pfmul mm7, mm5		/* mm7=vcoul */
	pfmul mm3, mm7		/* mm3=fscal for the two H's. */

	/* update vctot */
	pfadd mm7, mm6
	pfadd mm7, [esp + _vctot]
	movq [esp + _vctot], mm7
	
	/* spread oxygen fscalar to both positions */
	punpckldq mm4,mm4
	/* calc vectorial force for O */
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4
	
	/* pack j forces from H in the same form as the oxygen force. */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle. */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4 + 12]
	movd mm3,  [edi + eax*4 + 20]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 12], mm2
	movd [edi + eax*4 + 20], mm3

	/* interactions with j H2 */
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	movd mm6, [esp + _qqOH]
	movq mm7, [esp + _qqHH]

	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	movq  mm4, mm1
	pfmul mm4, mm4		/* mm4=invsq */
	/* calculate potential and scalar force */
	pfmul mm6, mm1		/* mm6=vcoul */
	pfmul mm4, mm6		/* mm4=fscalar */ 

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	/* seeds are in mm5 now, and rsq in mm3. */

	movq mm2, mm5
	pfmul mm5,mm5
        pfrsqit1 mm5,mm3				
        pfrcpit2 mm5,mm2	/* mm5=invsqrt */
	movq mm3,mm5
	pfmul mm3,mm3		/* mm3=invsq */
	pfmul mm7, mm5		/* mm7=vcoul */
	pfmul mm3, mm7		/* mm3=fscal for the two H's. */

	/* update vctot */
	pfadd mm7, mm6
	pfadd mm7, [esp + _vctot]
	movq [esp + _vctot], mm7
	
	/* spread oxygen fscalar to both positions */
	punpckldq mm4,mm4
	/* calc vectorial force for O */
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4	

	/* pack j forces from H in the same form as the oxygen force. */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle. */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4 + 24]
	movd mm3,  [edi + eax*4 + 32]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 24], mm2
	movd [edi + eax*4 + 32], mm3
	
	/*  done  - one more? */
	dec dword ptr [esp + _innerk]
	jz  .i1130_updateouterdata
	jmp .i1130_inner_loop	
.i1130_updateouterdata:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment iO force */ 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + _fixH]
	movq  mm3, [esp + _fiyH]
	movq  mm1, [esp + _fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	/* mm0(l)=fxH1, mm0(h)=fyH1 */
	punpckhdq mm2, mm3	/* mm2(l)=fxH2, mm2(h)=fyH2 */
	movq mm3, mm1
	pswapd mm3,mm3		
	/* mm1 is fzH1 */
	/* mm3 is fzH2 */

	movq  mm6, [edi + ecx*4 + 12]       /* increment iH1 force */ 
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       /* increment iH2 force */
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */

	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */

	movq  mm7, [esp + _vnbtot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */

	mov   eax, [ebp + _Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vnbtot[gid] */
	/* finish if last */
	dec dword ptr [ebp + _nri]
	jz  .i1130_end
	/* not last, iterate once more! */
	jmp .i1130_outer
.i1130_end:
	femms
	add esp, 196
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret




.globl inl3000_3dnow
	.type inl3000_3dnow,@function
inl3000_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56			
.equ		_tabscale,	60
.equ		_VFtab,		64
	/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_ix,		8
.equ		_iy,		12
.equ		_iz,		16
.equ		_iq,		20 
.equ		_vctot,		28 
.equ		_two,		36 
.equ		_n1,		44 
.equ		_tsc,		52 
.equ		_ntia,		60
.equ		_innerjjnr,	64
.equ		_innerk,	68		
.equ		_fix,		72
.equ		_fiy,		76
.equ		_fiz,		80
.equ		_dx1,		84
.equ		_dy1,		88
.equ		_dz1,		92
.equ		_dx2,		96
.equ		_dy2,		100
.equ		_dz2,		104						
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 108		/* local stack space */
	femms
	/* move data to local stack */ 
	movq  mm0, [mm_two]
	movd  mm3, [ebp + _tabscale]
	movq  [esp + _two],    mm0
	punpckldq mm3,mm3
	movq  [esp + _tsc], mm3	
	/* assume we have at least one i particle - start directly */	
.i3000_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm0, [eax + ebx*4]	/* move shX/shY to mm0 and shZ to mm1 */
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   edx, [ebp + _charge]
	movd  mm2, [edx + ebx*4]        /* mm2=charge[ii] */
	pfmul mm2, [ebp + _facel]
	punpckldq mm2,mm2	        /* spread to both halves */
	movq  [esp + _iq], mm2	        /* iq =facel*charge[ii] */

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	
	pfadd mm0, [eax + ebx*4]        /* ix = shX + posX (and iy too) */
	movd  mm3, [eax + ebx*4 + 8]    /* cant use direct memory add for 4 bytes (iz) */
	mov   [esp + _ii3], ebx
	pfadd mm1, mm3
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	
				
	/* clear total potential and i forces */
	pxor  mm7,mm7
	movq  [esp + _vctot],  mm7
	movq  [esp + _fix],    mm7
	movd  [esp + _fiz],    mm7

	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */

	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]	
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr], eax     /* pointer to jjnr[nj0] */
	sub   edx,  2
	mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i3000_unroll_loop
	jmp   .i3000_finish_inner
.i3000_unroll_loop:
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */
	
	mov ecx, [ebp + _charge]        /* base of charge[] */
	movq mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]	         /* charge[jnr1] */
        punpckldq mm3, [ecx + ebx*4]     /* move charge 2 to high part of mm3 */
	pfmul mm3,mm5		         /* mm3 now has qq for both particles */

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrsqrt mm0, mm4	         /* lookup inverse square root seed */
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs. */
        movq mm2,mm0	        	/* amd 3dnow N-R iteration to get full precision. */
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r. */
	/* do potential and fscal */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */
	
	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* coulomb table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */

	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, mm3	/* vcoul=qq*VV */
	pfmul mm3, mm7	/* fijC=FF*qq */ 

	/* at this point mm5 contains vcoul and mm3 fijC. */
	/* increment vcoul - then we can get rid of mm5. */
	/* update vctot */
	pfadd mm5, [esp + _vctot]      /* add the earlier value */
	movq [esp + _vctot], mm5       /* store the sum */      

	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3
	pfmul mm1, [esp + _tsc]	
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i3000_finish_inner
	jmp   .i3000_unroll_loop
.i3000_finish_inner:	
	and [esp + _innerk],  1
	jnz  .i3000_single_inner
	jmp  .i3000_updateouterdata		
.i3000_single_inner:
	/* a single j particle iteration here - compare with the unrolled code for comments. */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov ecx, [ebp + _charge]
	movd mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	/* mm3=qq */

	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm0=rsq */
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	/* mm1=invsqrt */
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r. */

	/* calculate potentials and scalar force */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */

	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */

	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, mm3	/* vcoul=qq*VV */
	pfmul mm3, mm7	/* fijC=FF*qq */ 

	/* at this point mm5 contains vcoul and mm3 fijC */
	/* increment vcoul - then we can get rid of mm5 */
	/* update vctot */
	pfadd mm5, [esp + _vctot]      /* add the earlier value */
	movq [esp + _vctot], mm5       /* store the sum */      
	
	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + _tsc]
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i3000_updateouterdata:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */

	/* finish if last */
	mov   ecx, [ebp + _nri]
	dec ecx
	jecxz .i3000_end
	/* not last, iterate once more! */
	mov [ebp + _nri], ecx
	jmp .i3000_outer
.i3000_end:
	femms
	add esp, 108
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret



	
.globl inl3010_3dnow
	.type inl3010_3dnow,@function
inl3010_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56
.equ		_tabscale,	60		
.equ		_VFtab,		64
.equ		_nsatoms,	68		
	/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_shX,		8
.equ		_shY,		12 
.equ		_shZ,		16	
.equ		_ix,		20
.equ		_iy,		24
.equ		_iz,		28	
.equ		_iq,		32 
.equ		_vctot,		40 
.equ		_two,		48 
.equ		_n1,		56 
.equ		_tsc,		64 			
.equ		_innerjjnr0,	72
.equ		_innerk0,	76		
.equ		_innerjjnr,	80
.equ		_innerk,	84		
.equ		_fix,		88
.equ		_fiy,		92
.equ		_fiz,		96
.equ		_dx1,		100
.equ		_dy1,		104
.equ		_dz1,		108
.equ		_dx2,		112
.equ		_dy2,		116
.equ		_dz2,		120								
.equ		_nscoul,	124
.equ		_solnr,		128		
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 132		/* local stack space */
	femms
	
	add   [ebp + _nsatoms],  8
	movq  mm2, [mm_two]
	movq  [esp + _two], mm2
	movd  mm3, [ebp + _tabscale]
	punpckldq mm3,mm3
	movq  [esp + _tsc], mm3
	
	/* assume we have at least one i particle - start directly */		
.i3010_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm0, [eax + ebx*4]	/* move shX/shY to mm0 and shZ to mm1 */
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + _shX], mm0
	movd  [esp + _shZ], mm1

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   eax, [ebp + _nsatoms]
	mov   ecx, [eax]
	add   [ebp + _nsatoms],  12
	mov   [esp + _nscoul], ecx
		
	/* clear potential */
	pxor  mm7,mm7
	movq  [esp + _vctot], mm7
	mov   [esp + _solnr], ebx
	
	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr0], eax     /* pointer to jjnr[nj0] */

	mov   [esp + _innerk0], edx        /* number of innerloop atoms */
	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]
	mov   ecx, [esp + _nscoul]
	cmp   ecx,  0
	jnz  .i3010_mno_coul
	jmp  .i3010_last_mno
.i3010_mno_coul:				
	mov   ebx,  [esp + _solnr]
	inc   dword ptr [esp + _solnr]
	mov   edx, [ebp + _charge]
	movd  mm2, [edx + ebx*4]        /* mm2=charge[ii] */
	pfmul mm2, [ebp + _facel]
	punpckldq mm2,mm2	        /* spread to both halves */
	movq  [esp + _iq], mm2	        /* iq =facel*charge[ii] */

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	mov   [esp + _ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + _shX]
	pfadd mm1, [esp + _shZ]
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	

	/* clear forces */
	pxor  mm7,mm7
	movq  [esp + _fix],   mm7
	movd  [esp + _fiz],   mm7

	mov   ecx, [esp + _innerjjnr0]
	mov   [esp + _innerjjnr], ecx
	mov   edx, [esp + _innerk0]
        sub   edx,  2
        mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i3010_unroll_coul_loop
	jmp   .i3010_finish_coul_inner
.i3010_unroll_coul_loop:	
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */
	
	mov ecx, [ebp + _charge]        /* base of charge[] */
	movq mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]	         /* charge[jnr1] */
        punpckldq mm3, [ecx + ebx*4]     /* move charge 2 to high part of mm3 */
	pfmul mm3,mm5		         /* mm3 now has qq for both particles */

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrsqrt mm0, mm4	         /* lookup inverse square root seed */
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs. */
        movq mm2,mm0	        	/* amd 3dnow N-R iteration to get full precision. */
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r. */
	/* do potential and fscal */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */
	
	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* coulomb table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */

	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, mm3	/* vcoul=qq*VV */
	pfmul mm3, mm7	/* fijC=FF*qq */ 

	/* at this point mm5 contains vcoul and mm3 fijC */
	/* increment vcoul - then we can get rid of mm5 */
	/* update vctot */
	pfadd mm5, [esp + _vctot]      /* add the earlier value */
	movq [esp + _vctot], mm5       /* store the sum */      

	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3
	pfmul mm1, [esp + _tsc]	
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i3010_finish_coul_inner
	jmp   .i3010_unroll_coul_loop
.i3010_finish_coul_inner:	
	and [esp + _innerk],  1
	jnz  .i3010_single_coul_inner
	jmp  .i3010_updateouterdata_coul		
.i3010_single_coul_inner:	
	/* a single j particle iteration here - compare with the unrolled code for comments. */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov ecx, [ebp + _charge]
	movd mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	/* mm3=qq */

	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm0=rsq */
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	/* mm1=invsqrt */
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r. */

	/* calculate potentials and scalar force */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */

	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */

	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, mm3	/* vcoul=qq*VV */
	pfmul mm3, mm7	/* fijC=FF*qq */ 

	/* at this point mm5 contains vcoul and mm3 fijC */
	/* increment vcoul - then we can get rid of mm5 */
	/* update vctot */
	pfadd mm5, [esp + _vctot]      /* add the earlier value */
	movq [esp + _vctot], mm5       /* store the sum */      
	
	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + _tsc]
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i3010_updateouterdata_coul:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	/* loop back to mno */
	dec dword ptr [esp + _nscoul]
	jz  .i3010_last_mno
	jmp .i3010_mno_coul
.i3010_last_mno:	
	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */
	/* finish if last */
	mov   ecx, [ebp + _nri]
	dec ecx
	jecxz .i3010_end
	/* not last, iterate once more! */
	mov [ebp + _nri], ecx
	jmp .i3010_outer
.i3010_end:
	femms
	add esp, 132
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret

	


.globl inl3020_3dnow
	.type inl3020_3dnow,@function
inl3020_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56			
.equ		_tabscale,	60
.equ		_VFtab,		64
			/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_ixO,		8
.equ		_iyO,		12
.equ		_izO,		16	
.equ		_ixH,		20  
.equ		_iyH,		28  
.equ		_izH,		36  
.equ		_iqO,		44  
.equ		_iqH,		52  
.equ		_qqO,		60  
.equ		_qqH,		68  
.equ		_vctot,		76  
.equ		_two,		84  
.equ		_n1,		92  
.equ		_tsc,		100 
.equ		_innerjjnr,	108
.equ		_innerk,	112	
.equ		_fixO,		116
.equ		_fiyO,		120
.equ		_fizO,		124
.equ		_fixH,		128 
.equ		_fiyH,		136 
.equ		_fizH,		144 
.equ		_dxO,		152
.equ		_dyO,		156
.equ		_dzO,		160
.equ		_dxH,		164 
.equ		_dyH,		172 
.equ		_dzH,		180 
.equ		_tmprsqH,	188 
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 196		/* local stack space */
	femms

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   edx, [ebp + _charge]
	movd  mm1, [ebp + _facel]
	movd  mm2, [edx + ebx*4]        /* mm2=charge[i.i0] */
	pfmul mm2, mm1		
	movq  [esp + _iqO], mm2	        /* iqO = facel*charge[ii] */
	
	movd  mm2, [edx + ebx*4 + 4]    /* mm2=charge[i.i0+1] */
	pfmul mm2, mm1
	punpckldq mm2,mm2	        /* spread to both halves */
	movq  [esp + _iqH], mm2	        /* iqH = facel*charge[i.i0+1] */

	movq  mm3, [mm_two]
	movd  mm4, [ebp + _tabscale]
	punpckldq mm4,mm4	        /* spread to both halves */
	movq  [esp + _two],    mm3
	movq  [esp + _tsc], mm4	      
	/* assume we have at least one i particle - start directly */	 
.i3020_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm5, [eax + ebx*4]	/* move shX/shY to mm5 and shZ to mm6. */
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        /* also expand shX,Y,Z in mm0--mm2. */
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	
	pfadd mm5, [eax + ebx*4]        /* ix = shX + posX (and iy too) */
	movd  mm7, [eax + ebx*4 + 8]    /* cant use direct memory add for 4 bytes (iz) */
	mov   [esp + _ii3], ebx	        /* (use mm7 as temp. storage for iz.) */
	pfadd mm6, mm7
	movq  [esp + _ixO], mm5	
	movq  [esp + _izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] /* coords of H1 in low mm3-mm5, H2 in high */
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + _ixH], mm0	
	movq [esp + _iyH], mm1	
	movq [esp + _izH], mm2	
					
	/* clear vctot and i forces */
	pxor  mm7,mm7
	movq  [esp + _vctot], mm7
	movq  [esp + _fixO],   mm7
	movd  [esp + _fizO],   mm7
	movq  [esp + _fixH],   mm7
	movq  [esp + _fiyH],   mm7
	movq  [esp + _fizH],   mm7

	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */
	mov   [esp + _innerk], edx        

	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]	
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr], eax     /* pointer to jjnr[nj0] */
.i3020_inner_loop:	
	/* a single j particle iteration */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	 /* eax=jnr offset */
        add   [esp + _innerjjnr],  4 /* advance pointer */
	prefetch [ecx + 16]	   /* prefetch data - trial and error says 16 is best */

	mov ecx, [ebp + _charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + _iqO]
	pfmul mm7, [esp + _iqH]	 /* mm6=qqO, mm7=qqH */
	movd [esp + _qqO], mm6
	movq [esp + _qqH], mm7
		
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */
	movq [esp + _tmprsqH], mm3
	
        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */

	pfmul mm0, mm1		/* mm0=r */

	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */

	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqO]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqO]	/* fijC=qq*FF */
	/* update vctot directly, use mm3 for fscal sum. */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	movq mm3, mm7	

	/* change sign of fscal and multiply with rinv */ 
        pxor mm0,mm0
	pfsubr mm3, mm0	
	pfmul mm3, [esp + _tsc]
 	pfmul mm3, mm1        /* mm3 is total fscal (for the oxygen) now */	
	
	/* Ready with the oxygen - potential is updated, fscal is in mm3. */
	/* now do the two hydrogens. */
	 
	movq mm0, [esp + _tmprsqH] /* mm0=r */sqH

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	/* seeds are in mm1 now, and rsq in mm0. */

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	
	pfmul mm0,mm1		/* mm0=r */
	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */

	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqH]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqH]	/* fijC=qq*FF */

	/* update vctot */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	
	/* change sign of fijC and multiply by rinv */
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + _tsc]
 	pfmul mm4, mm1        /* mm4 is total fscal (for the hydrogens) now */	

	/* spread oxygen fscalar to both positions */
	punpckldq mm3,mm3
	/* calc vectorial force for O */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4
	
	/* pack j forces from H in the same form as the oxygen force. */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle. */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 + 8], mm3
	
	/*  done  - one more? */
	dec dword ptr [esp + _innerk]
	jz  .i3020_updateouterdata
	jmp .i3020_inner_loop
.i3020_updateouterdata:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment iO force */ 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + _fixH]
	movq  mm3, [esp + _fiyH]
	movq  mm1, [esp + _fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	/* mm0(l)=fxH1, mm0(h)=fyH1 */
	punpckhdq mm2, mm3	/* mm2(l)=fxH2, mm2(h)=fyH2 */
	movq mm3, mm1
	pswapd mm3, mm3	
	/* mm1 is fzH1 */
	/* mm3 is fzH2 */
	
	movq  mm6, [edi + ecx*4 + 12]       /* increment iH1 force */ 
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       /* increment iH2 force */
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */

	/* finish if last */
	dec dword ptr [ebp + _nri]
	jz  .i3020_end
	/* not last, iterate once more! */
	jmp .i3020_outer
.i3020_end:
	femms
	add esp, 196
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret

	

.globl inl3030_3dnow
	.type inl3030_3dnow,@function
inl3030_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56			
.equ		_tabscale,	60
.equ		_VFtab,		64
			/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_ixO,		8
.equ		_iyO,		12
.equ		_izO,		16	
.equ		_ixH,		20  
.equ		_iyH,		28  
.equ		_izH,		36  
.equ		_qqOO,		44  
.equ		_qqOH,		52  
.equ		_qqHH,		60  
.equ		_two,		68  
.equ		_n1,		76  
.equ		_tsc,		84  
.equ		_vctot,		92  
.equ		_innerjjnr,	100
.equ		_innerk,	104	
.equ		_fixO,		108
.equ		_fiyO,		112
.equ		_fizO,		116
.equ		_fixH,		120 
.equ		_fiyH,		128 
.equ		_fizH,		136 
.equ		_dxO,		144
.equ		_dyO,		148
.equ		_dzO,		152
.equ		_dxH,		156 
.equ		_dyH,		164 
.equ		_dzH,		172 
.equ		_tmprsqH,	180 
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 188		/* local stack space */
	femms
	/* assume we have at least one i particle - start directly */	

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   edx, [ebp + _charge]
	movd  mm1, [ebp + _facel]	/* mm1=facel */
	movd  mm2, [edx + ebx*4]        /* mm2=charge[i.i0] (O) */
	movd  mm3, [edx + ebx*4 + 4]    /* mm2=charge[i.i0+1] (H) */ 
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			/* mm4=qqOO*facel */
	pfmul mm5, mm3			/* mm5=qqOH*facel */
	pfmul mm6, mm3			/* mm6=qqHH*facel */
	punpckldq mm5,mm5	        /* spread to both halves */
	punpckldq mm6,mm6	        /* spread to both halves */
	movq  [esp + _qqOO], mm4
	movq  [esp + _qqOH], mm5
	movq  [esp + _qqHH], mm6
	movq  mm2, [mm_two]
	movq  [esp + _two], mm2
	movd  mm3, [ebp + _tabscale]
	punpckldq mm3,mm3
	movq  [esp + _tsc], mm3
.i3030_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm5, [eax + ebx*4]	/* move shX/shY to mm5 and shZ to mm6. */
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        /* also expand shX,Y,Z in mm0--mm2. */
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */

	pfadd mm5, [eax + ebx*4]        /* ix = shX + posX (and iy too) */
	movd  mm7, [eax + ebx*4 + 8]    /* cant use direct memory add for 4 bytes (iz) */
	mov   [esp + _ii3], ebx	        /* (use mm7 as temp. storage for iz.) */
	pfadd mm6, mm7
	movq  [esp + _ixO], mm5	
	movq  [esp + _izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] /* coords of H1 in low mm3-mm5, H2 in high */
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + _ixH], mm0	
	movq [esp + _iyH], mm1	
	movq [esp + _izH], mm2	

	/* clear vctot and i forces */
	pxor  mm7,mm7
	movq  [esp + _vctot], mm7
	movq  [esp + _fixO],  mm7
	movq  [esp + _fizO],  mm7
	movq  [esp + _fixH],  mm7
	movq  [esp + _fiyH],  mm7
	movq  [esp + _fizH],  mm7

	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */
	mov   [esp + _innerk], edx        /* number of innerloop atoms */

	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]	
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr], eax     /* pointer to jjnr[nj0] */
.i3030_inner_loop:
	/* a single j particle iteration here - compare with the unrolled code for comments. */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */
        add   [esp + _innerjjnr],  4 /* advance pointer */

	lea   eax, [eax + eax*2]

	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */
	movq [esp + _tmprsqH], mm3

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */ OO
	pfmul mm0, mm1		/* mm0=rsq */ OO

	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */

	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2

	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqOO]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqOO]	/* fijC=qq*FF */

	/* update vctot directly, use mm3 for fscal sum. */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	movq mm3, mm7

	/* change sign of fscal and multiply with rinv */ 
        pxor mm0,mm0
	pfsubr mm3, mm0	
	pfmul mm3, [esp + _tsc]
 	pfmul mm3, mm1        /* mm3 is total fscal (for the oxygen) now */
	
	/* Ready with the oxygen - potential is updated, fscal is in mm3. */
	/* time for hydrogens! */

	movq mm0, [esp + _tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	/* seeds are in mm1 now, and rsq in mm0. */

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	
	pfmul mm0,mm1		/* mm0=r */
	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqOH]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqOH]	/* fijC=qq*FF */
	/* update vctot */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	
	/* change sign of fijC and multiply by rinv */
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + _tsc]
 	pfmul mm4, mm1        /* mm4 is total fscal (for the hydrogens) now */	
	
	/* spread oxygen fscalar to both positions */
	punpckldq mm3,mm3
	/* calc vectorial force for O */
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4
	
	/* pack j forces from H in the same form as the oxygen force. */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle. */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3

	/* interactions with j H1 */

	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */
	movq [esp + _tmprsqH], mm3

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	pfmul mm0, mm1		/* mm0=rsq */ 
	
	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */

	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2

	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqOH]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqOH]	/* fijC=qq*FF */

	/* update vctot  directly, force is moved to mm3 */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	pxor mm3, mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + _tsc]
	pfmul mm3, mm1        /* mm3 is total fscal (for the oxygen) now */

	movq mm0, [esp + _tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	/* seeds are in mm1 now, and rsq in mm0. */

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	
	pfmul mm0,mm1		/* mm0=r */
	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqHH]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqHH]	/* fijC=qq*FF */
	/* update vctot */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	
	/* change sign of fijC and multiply by rinv */
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + _tsc]
 	pfmul mm4, mm1        /* mm4 is total fscal (for the hydrogens) now */		

	/* spread oxygen fscalar to both positions */
	punpckldq mm3,mm3
	/* calc vectorial force for O */
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4
	
	/* pack j forces from H in the same form as the oxygen force. */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle. */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4 + 12]
	movd mm3,  [edi + eax*4 + 20]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 12], mm2
	movd [edi + eax*4 + 20], mm3

	/* interactions with j H2 */
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */
	movq [esp + _tmprsqH], mm3

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	pfmul mm0, mm1

	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */

	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2

	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqOH]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqOH]	/* fijC=qq*FF */

	/* update vctot directly, use mm3 for fscal sum. */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	pxor mm3,mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + _tsc]
 	pfmul mm3, mm1        /* mm3 is total fscal (for the oxygen) now */

	movq mm0, [esp + _tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	/* seeds are in mm1 now, and rsq in mm0. */

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	
	pfmul mm0,mm1		/* mm0=r */
	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqHH]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqHH]	/* fijC=qq*FF */
	/* update vctot */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	
	/* change sign of fijC and multiply by rinv */
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + _tsc]
 	pfmul mm4, mm1        /* mm4 is total fscal (for the hydrogens) now */	

	/* spread oxygen fscalar to both positions */
	punpckldq mm3,mm3
	/* calc vectorial force for O */
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4	

	/* pack j forces from H in the same form as the oxygen force. */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle. */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4 + 24]
	movd mm3,  [edi + eax*4 + 32]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 24], mm2
	movd [edi + eax*4 + 32], mm3
	
	/*  done  - one more? */
	dec dword ptr [esp + _innerk]
	jz  .i3030_updateouterdata
	jmp .i3030_inner_loop	
.i3030_updateouterdata:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment iO force */ 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + _fixH]
	movq  mm3, [esp + _fiyH]
	movq  mm1, [esp + _fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	/* mm0(l)=fxH1, mm0(h)=fyH1 */
	punpckhdq mm2, mm3	/* mm2(l)=fxH2, mm2(h)=fyH2 */
	movq mm3, mm1
	pswapd mm3,mm3		
	/* mm1 is fzH1 */
	/* mm3 is fzH2 */

	movq  mm6, [edi + ecx*4 + 12]       /* increment iH1 force */ 
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       /* increment iH2 force */
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */

	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */

	/* finish if last */
	dec dword ptr [ebp + _nri]
	jz  .i3030_end
	/* not last, iterate once more! */
	jmp .i3030_outer
.i3030_end:
	femms
	add esp, 188
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret




.globl inl3100_3dnow
	.type inl3100_3dnow,@function
inl3100_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56			
.equ		_type,		60
.equ		_ntype,		64
.equ		_nbfp,		68	
.equ		_Vnb,		72
.equ		_tabscale,	76
.equ		_VFtab,		80
	/* stack offsets for local variables */
.equ		_is3,		0 
.equ		_ii3,		4
.equ		_ix,		8
.equ		_iy,		12
.equ		_iz,		16
.equ		_iq,		20 
.equ		_vctot,		28 
.equ		_vnbtot,	36 
.equ		_c6,		44 
.equ		_c12,		52 
.equ		_six,		60 
.equ		_twelve,	68 
.equ		_two,		76 
.equ		_n1,		84 
.equ		_tsc,		92 
.equ		_ntia,		100
.equ		_innerjjnr,	104
.equ		_innerk,	108	
.equ		_fix,		112
.equ		_fiy,		116
.equ		_fiz,		120
.equ		_dx1,		124
.equ		_dy1,		128
.equ		_dz1,		132
.equ		_dx2,		136
.equ		_dy2,		140
.equ		_dz2,		144						
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 148		/* local stack space */
	femms
	/* move data to local stack */ 
	movq  mm0, [mm_two]
	movq  mm1, [mm_six]
	movq  mm2, [mm_twelve]
	movd  mm3, [ebp + _tabscale]
	movq  [esp + _two],    mm0
	movq  [esp + _six],    mm1
	movq  [esp + _twelve],    mm2
	punpckldq mm3,mm3
	movq  [esp + _tsc], mm3	
	/* assume we have at least one i particle - start directly */	
.i3100_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm0, [eax + ebx*4]	/* move shX/shY to mm0 and shZ to mm1 */
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   edx, [ebp + _charge]
	movd  mm2, [edx + ebx*4]        /* mm2=charge[ii] */
	pfmul mm2, [ebp + _facel]
	punpckldq mm2,mm2	        /* spread to both halves */
	movq  [esp + _iq], mm2	        /* iq =facel*charge[ii] */

	mov   edx, [ebp + _type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + _ntype]
	shl   edx, 1
	mov   [esp + _ntia], edx

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	
	pfadd mm0, [eax + ebx*4]        /* ix = shX + posX (and iy too) */
	movd  mm3, [eax + ebx*4 + 8]    /* cant use direct memory add for 4 bytes (iz) */
	mov   [esp + _ii3], ebx
	pfadd mm1, mm3
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	
				
	/* clear total potential and i forces */
	pxor  mm7,mm7
	movq  [esp + _vctot],  mm7
	movq  [esp + _vnbtot], mm7
	movq  [esp + _fix],    mm7
	movd  [esp + _fiz],    mm7

	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */

	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]	
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr], eax     /* pointer to jjnr[nj0] */
	sub   edx,  2
	mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i3100_unroll_loop
	jmp   .i3100_finish_inner
.i3100_unroll_loop:
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */
	
	mov ecx, [ebp + _charge]        /* base of charge[] */
	movq mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]	         /* charge[jnr1] */
        punpckldq mm3, [ecx + ebx*4]     /* move charge 2 to high part of mm3 */
	pfmul mm3,mm5		         /* mm3 now has qq for both particles */

	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	mov ecx, [ecx + ebx*4]           /* type [jnr2] */

	mov esi, [ebp + _nbfp]		/* base of nbfp */ 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	add ecx, [esp + _ntia]

	movq mm5, [esi + edx*4]		/* mm5 = 1st c6 / c12 */		
	movq mm7, [esi + ecx*4]		/* mm7 = 2nd c6 / c12 */	
	movq mm6,mm5			
	punpckldq mm5,mm7		/* mm5 = 1st c6 / 2nd c6 */
	punpckhdq mm6,mm7		/* mm6 = 1st c12 / 2nd c12 */
	movq [esp + _c6], mm5
	movq [esp + _c12], mm6

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrsqrt mm0, mm4	         /* lookup inverse square root seed */
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs. */
        movq mm2,mm0	        	/* amd 3dnow N-R iteration to get full precision. */
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r. */
	/* do potential and fscal */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */
	
	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* coulomb table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */

	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, mm3	/* vcoul=qq*VV */
	pfmul mm3, mm7	/* fijC=FF*qq */ 

	movq mm1, mm0
	pfmul mm1,mm1 	/* mm1=invsq */
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	/* mm2=rinvsix */
	movq  mm1,mm2
	pfmul mm1,mm1	/* mm1=rinvtwelve */
	
	pfmul mm3, [esp + _tsc]
	
	pfmul mm1, [esp + _c12]

	pfmul mm2, [esp + _c6]

	movq mm4, mm1
	pfsub mm4, mm2	/* mm4 = vnb12-vnb6 */

	pfmul mm2, [esp + _six]
	pfmul mm1, [esp + _twelve]

	pfsub mm1, mm2
	pfmul mm1, mm0	/* mm1=	(12*vnb12-6*vnb6)*rinv11 */

	pfsub mm1, mm3

	/* update vctot */
	pfadd mm5, [esp + _vctot]      /* add the earlier value */
	movq [esp + _vctot], mm5       /* store the sum */      

 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	/* update vnbtot */
	pfadd mm4, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm4       /* store the sum */      

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i3100_finish_inner
	jmp   .i3100_unroll_loop
.i3100_finish_inner:	
	and [esp + _innerk],  1
	jnz  .i3100_single_inner
	jmp  .i3100_updateouterdata		
.i3100_single_inner:
	/* a single j particle iteration here - compare with the unrolled code for comments. */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov ecx, [ebp + _charge]
	movd mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	/* mm3=qq */

	mov esi, [ebp + _nbfp]
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	shl edx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	movd mm5, [esi + edx*4]		/* mm5 = 1st c6 */ 		
	movq [esp + _c6], mm5
	movd mm5, [esi + edx*4 + 4]	/* mm5 = 1st c12 */ 		
	movq [esp + _c12], mm5


	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm4=rsq */
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	/* mm1=invsqrt */
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r. */
	/* calculate potentials and scalar force */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */

	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */

	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, mm3	/* vcoul=qq*VV */
	pfmul mm3, mm7	/* fijC=FF*qq */ 
	
	/* at this point mm5 contains vcoul and mm3 fijC */

	movq mm1, mm0
	pfmul mm1,mm1 	/* mm1=invsq */
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	/* mm2=rinvsix */
	movq  mm1,mm2
	pfmul mm1,mm1	/* mm1=rinvtwelve */
	
	pfmul mm3, [esp + _tsc]
	
	pfmul mm1, [esp + _c12]

	pfmul mm2, [esp + _c6]

	movq mm4, mm1
	pfsub mm4, mm2	/* mm4 = vnb12-vnb6 */

	pfmul mm2, [esp + _six]
	pfmul mm1, [esp + _twelve]

	pfsub mm1, mm2
	pfmul mm1, mm0	/* mm1=	(12*vnb12-6*vnb6)*rinv11 */

	pfsub mm1, mm3

	/* update vctot */
	pfadd mm5, [esp + _vctot]      /* add the earlier value */
	movq [esp + _vctot], mm5       /* store the sum */      

 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]

	/* update vnbtot */
	pfadd mm4, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm4       /* store the sum */      

	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i3100_updateouterdata:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix] 
	pfadd mm7, [esp + _fiz] 
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */
 
	movq  mm7, [esp + _vnbtot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vnb] 
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vnb[gid] */

	/* finish if last */
	mov   ecx, [ebp + _nri]
	dec ecx
	jecxz .i3100_end
	/* not last, iterate once more! */
	mov [ebp + _nri], ecx
	jmp .i3100_outer
.i3100_end:
	femms
	add esp, 148
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret







.globl inl3110_3dnow
	.type inl3110_3dnow,@function
inl3110_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56			
.equ		_type,		60
.equ		_ntype,		64
.equ		_nbfp,		68	
.equ		_Vnb,		72
.equ		_tabscale,	76
.equ		_VFtab,		80
.equ		_nsatoms,	84	
	/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_shX,		8
.equ		_shY,		12 
.equ		_shZ,		16	
.equ		_ix,		20
.equ		_iy,		24
.equ		_iz,		28	
.equ		_iq,		32  
.equ		_vctot,		40  
.equ		_vnbtot,	48  
.equ		_c6,		56  
.equ		_c12,		64  
.equ		_six,		72  
.equ		_twelve,	80  
.equ		_two,		88  
.equ		_n1,		96  
.equ		_tsc,		104 
.equ		_ntia,		112
.equ		_innerjjnr0,	116
.equ		_innerk0,	120	
.equ		_innerjjnr,	124
.equ		_innerk,	128	
.equ		_fix,		132
.equ		_fiy,		136
.equ		_fiz,		140
.equ		_dx1,		144
.equ		_dy1,		148
.equ		_dz1,		152
.equ		_dx2,		156
.equ		_dy2,		160
.equ		_dz2,		164								
.equ		_nsvdwc,	168
.equ		_nscoul,	172
.equ		_nsvdw,		176
.equ		_solnr,		180		
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 184		/* local stack space */
	femms
	movq  mm0, [mm_six]
	movq  mm1, [mm_twelve]
	movq  [esp + _six],    mm0
	movq  [esp + _twelve], mm1
	movq  mm2, [mm_two]
	movd  mm3, [ebp + _tabscale]
	movq  [esp + _two],    mm2
	punpckldq mm3,mm3
	movq  [esp + _tsc], mm3	
	/* assume we have at least one i particle - start directly */		
.i3110_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm0, [eax + ebx*4]	/* move shX/shY to mm0 and shZ to mm1 */
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + _shX], mm0
	movd  [esp + _shZ], mm1

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   eax, [ebp + _nsatoms]
	add   [ebp + _nsatoms],  12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + _nsvdwc], edx
	mov   [esp + _nscoul], eax
	mov   [esp + _nsvdw], ecx
		
	/* clear potential */
	pxor  mm7,mm7
	movq  [esp + _vctot],  mm7
	movq  [esp + _vnbtot], mm7
	mov   [esp + _solnr],  ebx
	
	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr0], eax     /* pointer to jjnr[nj0] */

	mov   [esp + _innerk0], edx        /* number of innerloop atoms */
	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]
	
	mov   ecx, [esp + _nsvdwc]
	cmp   ecx,  0
	jnz   .i3110_mno_vdwc
	jmp   .i3110_testcoul
.i3110_mno_vdwc:
	mov   ebx,  [esp + _solnr]
	inc   dword ptr [esp + _solnr]
	mov   edx, [ebp + _charge]
	movd  mm2, [edx + ebx*4]        /* mm2=charge[ii] */
	pfmul mm2, [ebp + _facel]
	punpckldq mm2,mm2	        /* spread to both halves */
	movq  [esp + _iq], mm2	        /* iq =facel*charge[ii] */

	mov   edx, [ebp + _type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + _ntype]
	shl   edx, 1
	mov   [esp + _ntia], edx	

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	mov   [esp + _ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + _shX]
	pfadd mm1, [esp + _shZ]
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	

	/* clear forces */
	pxor  mm7,mm7
	movq  [esp + _fix],   mm7
	movd  [esp + _fiz],   mm7

	mov   ecx, [esp + _innerjjnr0]
	mov   [esp + _innerjjnr], ecx
	mov   edx, [esp + _innerk0]
        sub   edx,  2
        mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i3110_unroll_vdwc_loop
	jmp   .i3110_finish_vdwc_inner
.i3110_unroll_vdwc_loop:	
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */
	
	mov ecx, [ebp + _charge]        /* base of charge[] */
	movq mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]	         /* charge[jnr1] */
        punpckldq mm3, [ecx + ebx*4]     /* move charge 2 to high part of mm3 */
	pfmul mm3,mm5		         /* mm3 now has qq for both particles */

	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	mov ecx, [ecx + ebx*4]           /* type [jnr2] */

	mov esi, [ebp + _nbfp]		/* base of nbfp */ 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	add ecx, [esp + _ntia]

	movq mm5, [esi + edx*4]		/* mm5 = 1st c6 / c12 */		
	movq mm7, [esi + ecx*4]		/* mm7 = 2nd c6 / c12 */	
	movq mm6,mm5			
	punpckldq mm5,mm7		/* mm5 = 1st c6 / 2nd c6 */
	punpckhdq mm6,mm7		/* mm6 = 1st c12 / 2nd c12 */
	movq [esp + _c6], mm5
	movq [esp + _c12], mm6

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrsqrt mm0, mm4	         /* lookup inverse square root seed */
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs. */
        movq mm2,mm0	        	/* amd 3dnow N-R iteration to get full precision. */
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r. */
	/* do potential and fscal */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */
	
	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* coulomb table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */

	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, mm3	/* vcoul=qq*VV */
	pfmul mm3, mm7	/* fijC=FF*qq */ 

	movq mm1, mm0
	pfmul mm1,mm1 	/* mm1=invsq */
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	/* mm2=rinvsix */
	movq  mm1,mm2
	pfmul mm1,mm1	/* mm1=rinvtwelve */
	
	pfmul mm3, [esp + _tsc]
	
	pfmul mm1, [esp + _c12]

	pfmul mm2, [esp + _c6]

	movq mm4, mm1
	pfsub mm4, mm2	/* mm4 = vnb12-vnb6 */

	pfmul mm2, [esp + _six]
	pfmul mm1, [esp + _twelve]

	pfsub mm1, mm2
	pfmul mm1, mm0	/* mm1=	(12*vnb12-6*vnb6)*rinv11 */

	pfsub mm1, mm3

	/* update vctot */
	pfadd mm5, [esp + _vctot]      /* add the earlier value */
	movq [esp + _vctot], mm5       /* store the sum */      

 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	/* update vnbtot */
	pfadd mm4, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm4       /* store the sum */      

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7	
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i3110_finish_vdwc_inner
	jmp   .i3110_unroll_vdwc_loop
.i3110_finish_vdwc_inner:	
	and [esp + _innerk],  1
	jnz  .i3110_single_vdwc_inner
	jmp  .i3110_updateouterdata_vdwc		
.i3110_single_vdwc_inner:	
	/* a single j particle iteration here - compare with the unrolled code for comments. */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov ecx, [ebp + _charge]
	movd mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	/* mm3=qq */

	mov esi, [ebp + _nbfp]
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	shl edx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	movd mm5, [esi + edx*4]		/* mm5 = 1st c6 */ 		
	movq [esp + _c6], mm5
	movd mm5, [esi + edx*4 + 4]	/* mm5 = 1st c12 */ 		
	movq [esp + _c12], mm5


	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm4=rsq */
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	/* mm1=invsqrt */
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r. */
	/* calculate potentials and scalar force */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */

	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */

	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, mm3	/* vcoul=qq*VV */
	pfmul mm3, mm7	/* fijC=FF*qq */ 

	movq mm1, mm0
	pfmul mm1,mm1 	/* mm1=invsq */
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	/* mm2=rinvsix */
	movq  mm1,mm2
	pfmul mm1,mm1	/* mm1=rinvtwelve */
	
	pfmul mm3, [esp + _tsc]
	
	pfmul mm1, [esp + _c12]

	pfmul mm2, [esp + _c6]

	movq mm4, mm1
	pfsub mm4, mm2	/* mm4 = vnb12-vnb6 */

	pfmul mm2, [esp + _six]
	pfmul mm1, [esp + _twelve]

	pfsub mm1, mm2
	pfmul mm1, mm0	/* mm1=	(12*vnb12-6*vnb6)*rinv11 */

	pfsub mm1, mm3

	/* update vctot */
	pfadd mm5, [esp + _vctot]      /* add the earlier value */
	movq [esp + _vctot], mm5       /* store the sum */      

 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]

	/* update vnbtot */
	pfadd mm4, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm4       /* store the sum */      

	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i3110_updateouterdata_vdwc:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	/* loop back to mno */
	dec dword ptr [esp + _nsvdwc]
	jz  .i3110_testcoul
	jmp .i3110_mno_vdwc
.i3110_testcoul:	
	mov  ecx, [esp + _nscoul]
	cmp  ecx,  0
	jnz  .i3110_mno_coul
	jmp  .i3110_testvdw
.i3110_mno_coul:
	mov   ebx,  [esp + _solnr]
	inc   dword ptr [esp + _solnr]
	mov   edx, [ebp + _charge]
	movd  mm2, [edx + ebx*4]        /* mm2=charge[ii] */
	pfmul mm2, [ebp + _facel]
	punpckldq mm2,mm2	        /* spread to both halves */
	movq  [esp + _iq], mm2	        /* iq =facel*charge[ii] */

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	mov   [esp + _ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + _shX]
	pfadd mm1, [esp + _shZ]
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	

	/* clear forces */
	pxor  mm7,mm7
	movq  [esp + _fix],   mm7
	movd  [esp + _fiz],   mm7

	mov   ecx, [esp + _innerjjnr0]
	mov   [esp + _innerjjnr], ecx
	mov   edx, [esp + _innerk0]
        sub   edx,  2
        mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i3110_unroll_coul_loop
	jmp   .i3110_finish_coul_inner
.i3110_unroll_coul_loop:	
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */
	
	mov ecx, [ebp + _charge]        /* base of charge[] */
	movq mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]	         /* charge[jnr1] */
        punpckldq mm3, [ecx + ebx*4]     /* move charge 2 to high part of mm3 */
	pfmul mm3,mm5		         /* mm3 now has qq for both particles */

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrsqrt mm0, mm4	         /* lookup inverse square root seed */
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs. */
        movq mm2,mm0	        	/* amd 3dnow N-R iteration to get full precision. */
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r. */
	/* do potential and fscal */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */
	
	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* coulomb table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */

	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, mm3	/* vcoul=qq*VV */
	pfmul mm3, mm7	/* fijC=FF*qq */ 

	/* at this point mm5 contains vcoul and mm3 fijC */
	/* increment vcoul - then we can get rid of mm5 */
	/* update vctot */
	pfadd mm5, [esp + _vctot]      /* add the earlier value */
	movq [esp + _vctot], mm5       /* store the sum */      

	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3
	pfmul mm1, [esp + _tsc]	
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i3110_finish_coul_inner
	jmp   .i3110_unroll_coul_loop
.i3110_finish_coul_inner:	
	and [esp + _innerk],  1
	jnz  .i3110_single_coul_inner
	jmp  .i3110_updateouterdata_coul		
.i3110_single_coul_inner:	
	/* a single j particle iteration here - compare with the unrolled code for comments. */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov ecx, [ebp + _charge]
	movd mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	/* mm3=qq */

	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm0=rsq */
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	/* mm1=invsqrt */
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r. */

	/* calculate potentials and scalar force */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */

	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */

	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, mm3	/* vcoul=qq*VV */
	pfmul mm3, mm7	/* fijC=FF*qq */ 

	/* at this point mm5 contains vcoul and mm3 fijC */
	/* increment vcoul - then we can get rid of mm5 */
	/* update vctot */
	pfadd mm5, [esp + _vctot]      /* add the earlier value */
	movq [esp + _vctot], mm5       /* store the sum */      
	
	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + _tsc]
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i3110_updateouterdata_coul:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	/* loop back to mno */
	dec dword ptr [esp + _nscoul]
	jz  .i3110_testvdw
	jmp .i3110_mno_coul
.i3110_testvdw:	
	mov  ecx, [esp + _nsvdw]
	cmp  ecx,  0
	jnz  .i3110_mno_vdw
	jmp  .i3110_last_mno
.i3110_mno_vdw:
	mov   ebx,  [esp + _solnr]
	inc   dword ptr [esp + _solnr]

	mov   edx, [ebp + _type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + _ntype]
	shl   edx, 1
	mov   [esp + _ntia], edx	

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	mov   [esp + _ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + _shX]
	pfadd mm1, [esp + _shZ]
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	

	/* clear forces */
	pxor  mm7,mm7
	movq  [esp + _fix],   mm7
	movd  [esp + _fiz],   mm7

	mov   ecx, [esp + _innerjjnr0]
	mov   [esp + _innerjjnr], ecx
	mov   edx, [esp + _innerk0]
        sub   edx,  2
        mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i3110_unroll_vdw_loop
	jmp   .i3110_finish_vdw_inner
.i3110_unroll_vdw_loop:	
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */	

	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	mov ecx, [ecx + ebx*4]           /* type [jnr2] */

	mov esi, [ebp + _nbfp]		/* base of nbfp */ 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	add ecx, [esp + _ntia]

	movq mm5, [esi + edx*4]		/* mm5 = 1st c6 / c12 */		
	movq mm7, [esi + ecx*4]		/* mm7 = 2nd c6 / c12 */	
	movq mm6,mm5			
	punpckldq mm5,mm7		/* mm5 = 1st c6 / 2nd c6 */
	punpckhdq mm6,mm7		/* mm6 = 1st c12 / 2nd c12 */
	movq [esp + _c6], mm5
	movq [esp + _c12], mm6

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrcp mm0, mm4	                 /* lookup reciprocal seed */ 
        pfrcp mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs. */
                  	        	/* amd 3dnow N-R iteration to get full precision. */
        pfrcpit1 mm4,mm0				
        pfrcpit2 mm4,mm0	
	/* mm4 now contains invsq,
	 * do potential and fscal 
	 */
	movq  mm0, mm4
	pfmul mm4, mm0
	pfmul mm4, mm0             	/* mm4=rinvsix */
	movq  mm5, mm4	
	pfmul mm5, mm5	                /* mm5=rinvtwelve */

	pfmul mm5, [esp + _c12]
	pfmul mm4, [esp + _c6]	
	movq mm6, mm5	/* mm6 is vnb12-vnb6 */ 
	pfsub mm6, mm4

	pfmul mm4, [esp + _six]

	pfmul mm5, [esp + _twelve]
	pfsub mm5,mm4
 	pfmul mm0, mm5        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	/* update vnbtot */
	pfadd mm6, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm6       /* store the sum */      

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i3110_finish_vdw_inner
	jmp   .i3110_unroll_vdw_loop
.i3110_finish_vdw_inner:	
	and [esp + _innerk],  1
	jnz  .i3110_single_vdw_inner
	jmp  .i3110_updateouterdata_vdw		
.i3110_single_vdw_inner:	
	/* a single j particle iteration here - compare with the unrolled code for comments. */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov esi, [ebp + _nbfp]
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	shl edx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	movd mm5, [esi + edx*4]		/* mm5 = 1st c6 */ 		
	movq [esp + _c6], mm5
	movd mm5, [esi + edx*4 + 4]	/* mm5 = 1st c12 */ 		
	movq [esp + _c12], mm5

	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm4=rsq */
	
        pfrcp mm0,mm4
        pfrcpit1 mm4,mm0				
        pfrcpit2 mm4,mm0	/* mm4=invsq */ 
	/* calculate potentials and scalar force */
	movq  mm0, mm4

	pfmul mm4, mm0
	pfmul mm4, mm0             	/* mm4=rinvsix */
	movq  mm5, mm4	
	pfmul mm5, mm5	                /* mm5=rinvtwelve */

	pfmul mm5, [esp + _c12]
	pfmul mm4, [esp + _c6]	
	movq mm6, mm5	/* mm6 is vnb12-vnb6 */ 
	pfsub mm6, mm4

	pfmul mm4, [esp + _six]

	pfmul mm5, [esp + _twelve]
	pfsub mm5, mm4
 	pfmul mm0, mm5        /* mm0 is total fscal now */

	/* update vnbtot */
	pfadd mm6, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm6       /* store the sum */      

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i3110_updateouterdata_vdw:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	/* loop back to mno */
	dec dword ptr [esp + _nsvdw]
	jz  .i3110_last_mno
	jmp .i3110_mno_vdw
	
.i3110_last_mno:	
	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */

	movq  mm7, [esp + _vnbtot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */
	/* finish if last */
	mov   ecx, [ebp + _nri]
	dec ecx
	jecxz .i3110_end
	/* not last, iterate once more! */
	mov [ebp + _nri], ecx
	jmp .i3110_outer
.i3110_end:
	femms
	add esp, 184
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret


	

.globl inl3120_3dnow
	.type inl3120_3dnow,@function
inl3120_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56			
.equ		_type,		60
.equ		_ntype,		64
.equ		_nbfp,		68	
.equ		_Vnb,		72
.equ		_tabscale,	76
.equ		_VFtab,		80
			/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_ixO,		8
.equ		_iyO,		12
.equ		_izO,		16	
.equ		_ixH,		20  
.equ		_iyH,		28  
.equ		_izH,		36  
.equ		_iqO,		44  
.equ		_iqH,		52  
.equ		_qqO,		60  
.equ		_qqH,		68  
.equ		_vctot,		76  
.equ		_vnbtot,	84  
.equ		_c6,		92  
.equ		_c12,		100 
.equ		_six,		108 
.equ		_twelve,	116 
.equ		_two,		124 
.equ		_n1,		132 
.equ		_tsc,		140 
.equ		_ntia,		148 
.equ		_innerjjnr,	156
.equ		_innerk,	160	
.equ		_fixO,		164
.equ		_fiyO,		168
.equ		_fizO,		172
.equ		_fixH,		176 
.equ		_fiyH,		184 
.equ		_fizH,		192 
.equ		_dxO,		200
.equ		_dyO,		204
.equ		_dzO,		208
.equ		_dxH,		212 
.equ		_dyH,		220 
.equ		_dzH,		228 
.equ		_tmprsqH,	236 
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 244		/* local stack space */
	femms

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   edx, [ebp + _charge]
	movd  mm1, [ebp + _facel]
	movd  mm2, [edx + ebx*4]        /* mm2=charge[i.i0] */
	pfmul mm2, mm1
	movq  [esp + _iqO], mm2	        /* iqO = facel*charge[ii] */
	
	movd  mm2, [edx + ebx*4 + 4]    /* mm2=charge[i.i0+1] */
	pfmul mm2, mm1
	punpckldq mm2,mm2	        /* spread to both halves */
	movq  [esp + _iqH], mm2	        /* iqH = facel*charge[i.i0+1] */

	mov   edx, [ebp + _type] 	
	mov   edx, [edx + ebx*4]
	shl   edx, 1
	mov   ecx, edx		        
	imul  ecx, [ebp + _ntype]      /* ecx = ntia = 2*ntype*type[i.i0] */ 
	mov   [esp + _ntia], ecx
	 	
	movq  mm3, [mm_two]
	movq  mm4, [mm_six]
	movq  mm5, [mm_twelve]
	movq  mm6, [ebp + _tabscale]
	punpckldq mm6,mm6	        /* spread to both halves */
	movq  [esp + _two], mm3
	movq  [esp + _six], mm4
	movq  [esp + _twelve], mm5
	movq  [esp + _tsc], mm6	      
 	/* assume we have at least one i particle - start directly */	
.i3120_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm5, [eax + ebx*4]	/* move shX/shY to mm5 and shZ to mm6. */
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        /* also expand shX,Y,Z in mm0--mm2. */
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	
	pfadd mm5, [eax + ebx*4]        /* ix = shX + posX (and iy too) */
	movd  mm7, [eax + ebx*4 + 8]    /* cant use direct memory add for 4 bytes (iz) */
	mov   [esp + _ii3], ebx	        /* (use mm7 as temp. storage for iz.) */
	pfadd mm6, mm7
	movq  [esp + _ixO], mm5	
	movq  [esp + _izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] /* coords of H1 in low mm3-mm5, H2 in high */
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + _ixH], mm0	
	movq [esp + _iyH], mm1	
	movq [esp + _izH], mm2	
					
	/* clear vctot and i forces */
	pxor  mm7,mm7
	movq  [esp + _vctot], mm7
	movq  [esp + _vnbtot], mm7
	movq  [esp + _fixO],   mm7
	movd  [esp + _fizO],   mm7
	movq  [esp + _fixH],   mm7
	movq  [esp + _fiyH],   mm7
	movq  [esp + _fizH],   mm7

	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */
	mov   [esp + _innerk], edx        

	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]	
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr], eax     /* pointer to jjnr[nj0] */
.i3120_inner_loop:	
	/* a single j particle iteration here - compare with the unrolled code for comments. */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */
        add   [esp + _innerjjnr],  4 /* advance pointer */
	prefetch [ecx + 16]	   /* prefetch data - trial and error says 16 is best */

	mov ecx, [ebp + _charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + _iqO]
	pfmul mm7, [esp + _iqH]	/* mm6=qqO, mm7=qqH */
	movd [esp + _qqO], mm6
	movq [esp + _qqH], mm7

	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr] */
	mov ecx, [ebp + _nbfp]
	shl edx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	movd mm5, [ecx + edx*4]		/* mm5 = 1st c6 */ 		
	movq [esp + _c6], mm5
	movd mm5, [ecx + edx*4 + 4]	/* mm5 = 1st c12 */ 		
	movq [esp + _c12], mm5	
			
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */
	movq [esp + _tmprsqH], mm3
	
        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */

	pfmul mm0, mm1		/* mm0=r */

	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */

	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqO]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqO]	/* fijC=qq*FF */
	/* update vctot directly, use mm3 for fscal sum. */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5

	movq mm3, mm7
	pfmul mm3, [esp + _tsc]
	
	/* nontabulated LJ - mm1 is invsqrt. - keep mm1! */
	movq mm0, mm1
	pfmul mm0, mm0		/* mm0 is invsq */
	movq mm2, mm0
	pfmul mm2, mm0
	pfmul mm2, mm0		/* mm2 = rinvsix */
	movq mm4, mm2
	pfmul mm4, mm4		/* mm4=rinvtwelve */

	pfmul mm4, [esp + _c12]
	pfmul mm2, [esp + _c6]
	movq mm5, mm4
	pfsub mm5, mm2		/* mm5=vnb12-vnb6 */

	pfmul mm2, [esp + _six]
	pfmul mm4, [esp + _twelve]
	pfsub mm4, mm2
	pfmul mm4, mm1        /* mm4=(12*vnb12-6*vnb6)*rinv11 */

	pfsubr mm3, mm4 
 	pfmul mm3, mm1        /* mm3 is total fscal (for the oxygen) now */
	
	/* update vnbtot */ 
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      
	
	/* Ready with the oxygen - potential is updated, fscal is in mm3. */
	/* now do the two hydrogens. */
	movq mm0, [esp + _tmprsqH] /* mm0=r */sqH

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	/* seeds are in mm1 now, and rsq in mm0. */

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	
	pfmul mm0,mm1		/* mm0=r */
	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqH]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqH]	/* fijC=qq*FF */
	/* update vctot */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	
	/* change sign of fijC and multiply by rinv */
        pxor mm4,mm4
	pfsub mm4, mm7
	pfmul mm4, [esp + _tsc]	
 	pfmul mm4, mm1        /* mm4 is total fscal (for the hydrogens) now */	
	
	/* spread oxygen fscalar to both positions */
	punpckldq mm3,mm3
	/* calc vectorial force for O */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4
	
	/* pack j forces from H in the same form as the oxygen force. */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle. */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	
	/*  done  - one more? */
	dec dword ptr [esp + _innerk]
	jz  .i3120_updateouterdata
	jmp .i3120_inner_loop
.i3120_updateouterdata:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment iO force */ 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + _fixH]
	movq  mm3, [esp + _fiyH]
	movq  mm1, [esp + _fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	/* mm0(l)=fxH1, mm0(h)=fyH1 */
	punpckhdq mm2, mm3	/* mm2(l)=fxH2, mm2(h)=fyH2 */
	movq mm3, mm1
	pswapd mm3,mm3		
	/* mm1 is fzH1 */
	/* mm3 is fzH2 */
	
	movq  mm6, [edi + ecx*4 + 12]       /* increment iH1 force */ 
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       /* increment iH2 force */
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */

	movq  mm7, [esp + _vnbtot]     
	pfacc mm7,mm7	              /* same for Vnb */
	
	mov   eax, [ebp + _Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vnb[gid] */
	/* finish if last */
	dec dword ptr [ebp + _nri]
	jz  .i3120_end
	/* not last, iterate once more! */
	jmp .i3120_outer
.i3120_end:
	femms
	add esp, 244
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret





.globl inl3130_3dnow
	.type inl3130_3dnow,@function
inl3130_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56			
.equ		_type,		60
.equ		_ntype,		64
.equ		_nbfp,		68	
.equ		_Vnb,		72
.equ		_tabscale,	76
.equ		_VFtab,		80
			/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_ixO,		8
.equ		_iyO,		12
.equ		_izO,		16	
.equ		_ixH,		20  
.equ		_iyH,		28  
.equ		_izH,		36  
.equ		_qqOO,		44  
.equ		_qqOH,		52  
.equ		_qqHH,		60  
.equ		_c6,		68  
.equ		_c12,		76  
.equ		_six,		84  
.equ		_twelve,	92  
.equ		_two,		100 
.equ		_n1,		108 
.equ		_tsc,		116 
.equ		_vctot,		124 
.equ		_vnbtot,	132 
.equ		_innerjjnr,	140
.equ		_innerk,	144	
.equ		_fixO,		148
.equ		_fiyO,		152
.equ		_fizO,		156
.equ		_fixH,		160 
.equ		_fiyH,		168 
.equ		_fizH,		176 
.equ		_dxO,		184
.equ		_dyO,		188
.equ		_dzO,		192
.equ		_dxH,		200 
.equ		_dyH,		208 
.equ		_dzH,		216 
.equ		_tmprsqH,	224 
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 232		/* local stack space */
	femms
	/* assume we have at least one i particle - start directly */	

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   edx, [ebp + _charge]
	movd  mm1, [ebp + _facel]	/* mm1=facel */
	movd  mm2, [edx + ebx*4]        /* mm2=charge[i.i0] (O) */
	movd  mm3, [edx + ebx*4 + 4]    /* mm2=charge[i.i0+1] (H) */ 
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			/* mm4=qqOO*facel */
	pfmul mm5, mm3			/* mm5=qqOH*facel */
	pfmul mm6, mm3			/* mm6=qqHH*facel */
	punpckldq mm5,mm5	        /* spread to both halves */
	punpckldq mm6,mm6	        /* spread to both halves */
	movq  [esp + _qqOO], mm4
	movq  [esp + _qqOH], mm5
	movq  [esp + _qqHH], mm6
	mov   edx, [ebp + _type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [ebp + _ntype]
	add   edx, ecx
	mov   eax, [ebp + _nbfp]
	movd  mm0, [eax + edx*4]
	movd  mm1, [eax + edx*4 + 4]
	movq  [esp + _c6], mm0
	movq  [esp + _c12], mm1
	movq  mm2, [mm_two]
	movq  mm3, [mm_six]
	movq  mm4, [mm_twelve]
	movq  [esp + _two], mm2
	movq  [esp + _six], mm3
	movq  [esp + _twelve], mm4
	movd  mm5, [ebp + _tabscale]
	punpckldq mm5,mm5
	movq  [esp + _tsc], mm5
.i3130_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm5, [eax + ebx*4]	/* move shX/shY to mm5 and shZ to mm6. */
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        /* also expand shX,Y,Z in mm0--mm2. */
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */

	pfadd mm5, [eax + ebx*4]        /* ix = shX + posX (and iy too) */
	movd  mm7, [eax + ebx*4 + 8]    /* cant use direct memory add for 4 bytes (iz) */
	mov   [esp + _ii3], ebx	        /* (use mm7 as temp. storage for iz.) */
	pfadd mm6, mm7
	movq  [esp + _ixO], mm5	
	movq  [esp + _izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] /* coords of H1 in low mm3-mm5, H2 in high */
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + _ixH], mm0	
	movq [esp + _iyH], mm1	
	movq [esp + _izH], mm2	

	/* clear vctot and i forces */
	pxor  mm7,mm7
	movq  [esp + _vctot], mm7
	movq  [esp + _vnbtot], mm7
	movq  [esp + _fixO],  mm7
	movq  [esp + _fizO],  mm7
	movq  [esp + _fixH],  mm7
	movq  [esp + _fiyH],  mm7
	movq  [esp + _fizH],  mm7

	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */
	mov   [esp + _innerk], edx        /* number of innerloop atoms */

	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]	
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr], eax     /* pointer to jjnr[nj0] */
.i3130_inner_loop:
	/* a single j particle iteration here - compare with the unrolled code for comments. */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */
        add   [esp + _innerjjnr],  4 /* advance pointer */

	lea   eax, [eax + eax*2]

	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */
	movq [esp + _tmprsqH], mm3

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */ OO
	pfmul mm0, mm1		/* mm0=rsq */ OO

	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */

	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2

	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqOO]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqOO]	/* fijC=qq*FF */

	/* update vctot directly, use mm3 for fscal sum. */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	movq mm3, mm7
	pfmul mm3, [esp + _tsc]
	
	movq mm5, mm1
	pfmul mm5,mm5
	movq mm4, mm5
	pfmul mm4,mm5
	pfmul mm4,mm5
	movq mm5, mm4
	pfmul mm5,mm5	/* mm4=rinvsix, mm5=rinvtwelve */

	pfmul mm4, [esp + _c6]
	pfmul mm5, [esp + _c12]
	movq mm6,mm5
	pfsub mm6,mm4

	pfmul mm4, [esp + _six]
	pfmul mm5, [esp + _twelve]
	pfsub mm5,mm4
	pfmul mm5, mm1
	pfsubr mm3, mm5

 	pfmul mm3, mm1        /* mm3 is total fscal (for the oxygen) now */
	
	/* update vnbtot */ 
	pfadd mm6, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm6       /* store the sum */      
	
	/* Ready with the oxygen - potential is updated, fscal is in mm3. */
	/* time for hydrogens! */

	movq mm0, [esp + _tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	/* seeds are in mm1 now, and rsq in mm0. */

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	
	pfmul mm0,mm1		/* mm0=r */
	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqOH]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqOH]	/* fijC=qq*FF */
	/* update vctot */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	
	/* change sign of fijC and multiply by rinv */
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + _tsc]
 	pfmul mm4, mm1        /* mm4 is total fscal (for the hydrogens) now */	
	
	/* spread oxygen fscalar to both positions */
	punpckldq mm3,mm3
	/* calc vectorial force for O */
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4
	
	/* pack j forces from H in the same form as the oxygen force. */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle. */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3

	/* interactions with j H1 */

	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */
	movq [esp + _tmprsqH], mm3

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	pfmul mm0, mm1		/* mm0=rsq */ 
	
	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */

	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2

	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqOH]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqOH]	/* fijC=qq*FF */

	/* update vctot  directly, force is moved to mm3 */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	pxor mm3, mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + _tsc]
	pfmul mm3, mm1        /* mm3 is total fscal (for the oxygen) now */

	movq mm0, [esp + _tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	/* seeds are in mm1 now, and rsq in mm0. */

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	
	pfmul mm0,mm1		/* mm0=r */
	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqHH]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqHH]	/* fijC=qq*FF */
	/* update vctot */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	
	/* change sign of fijC and multiply by rinv */
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + _tsc]
 	pfmul mm4, mm1        /* mm4 is total fscal (for the hydrogens) now */		

	/* spread oxygen fscalar to both positions */
	punpckldq mm3,mm3
	/* calc vectorial force for O */
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4
	
	/* pack j forces from H in the same form as the oxygen force. */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle. */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4 + 12]
	movd mm3,  [edi + eax*4 + 20]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 12], mm2
	movd [edi + eax*4 + 20], mm3

	/* interactions with j H2 */
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */
	movq [esp + _tmprsqH], mm3

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	pfmul mm0, mm1

	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */

	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2

	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqOH]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqOH]	/* fijC=qq*FF */

	/* update vctot directly, use mm3 for fscal sum. */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	pxor mm3,mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + _tsc]
 	pfmul mm3, mm1        /* mm3 is total fscal (for the oxygen) now */

	movq mm0, [esp + _tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	/* seeds are in mm1 now, and rsq in mm0. */

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	
	pfmul mm0,mm1		/* mm0=r */
	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	shl ecx, 2
	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqHH]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqHH]	/* fijC=qq*FF */
	/* update vctot */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	
	/* change sign of fijC and multiply by rinv */
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + _tsc]
 	pfmul mm4, mm1        /* mm4 is total fscal (for the hydrogens) now */	

	/* spread oxygen fscalar to both positions */
	punpckldq mm3,mm3
	/* calc vectorial force for O */
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4	

	/* pack j forces from H in the same form as the oxygen force. */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle. */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4 + 24]
	movd mm3,  [edi + eax*4 + 32]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 24], mm2
	movd [edi + eax*4 + 32], mm3
	
	/*  done  - one more? */
	dec dword ptr [esp + _innerk]
	jz  .i3130_updateouterdata
	jmp .i3130_inner_loop	
.i3130_updateouterdata:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment iO force */ 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + _fixH]
	movq  mm3, [esp + _fiyH]
	movq  mm1, [esp + _fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	/* mm0(l)=fxH1, mm0(h)=fyH1 */
	punpckhdq mm2, mm3	/* mm2(l)=fxH2, mm2(h)=fyH2 */
	movq mm3, mm1
	pswapd mm3,mm3		
	/* mm1 is fzH1 */
	/* mm3 is fzH2 */

	movq  mm6, [edi + ecx*4 + 12]       /* increment iH1 force */ 
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       /* increment iH2 force */
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */

	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */

	movq  mm7, [esp + _vnbtot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */

	mov   eax, [ebp + _Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vnbtot[gid] */
	/* finish if last */
	dec dword ptr [ebp + _nri]
	jz  .i3130_end
	/* not last, iterate once more! */
	jmp .i3130_outer
.i3130_end:
	femms
	add esp, 232
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret


.globl inl3300_3dnow
	.type inl3300_3dnow,@function
inl3300_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56			
.equ		_type,		60
.equ		_ntype,		64
.equ		_nbfp,		68	
.equ		_Vnb,		72
.equ		_tabscale,	76
.equ		_VFtab,		80
	/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_ix,		8
.equ		_iy,		12
.equ		_iz,		16
.equ		_iq,		20  
.equ		_vctot,		28  
.equ		_vnbtot,	36  
.equ		_c6,		44  
.equ		_c12,		52  
.equ		_two,		60  
.equ		_n1,		68  
.equ		_tsc,		76  
.equ		_ntia,		84
.equ		_innerjjnr,	88
.equ		_innerk,	92		
.equ		_fix,		96
.equ		_fiy,		100
.equ		_fiz,		104
.equ		_dx1,		108
.equ		_dy1,		112
.equ		_dz1,		116
.equ		_dx2,		120
.equ		_dy2,		124
.equ		_dz2,		128						
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 132		/* local stack space */
	femms
	/* move data to local stack */ 
	movq  mm0, [mm_two]
	movd  mm3, [ebp + _tabscale]
	movq  [esp + _two],    mm0
	punpckldq mm3,mm3
	movq  [esp + _tsc], mm3	
	/* assume we have at least one i particle - start directly */	
.i3300_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm0, [eax + ebx*4]	/* move shX/shY to mm0 and shZ to mm1 */
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   edx, [ebp + _charge]
	movd  mm2, [edx + ebx*4]        /* mm2=charge[ii] */
	pfmul mm2, [ebp + _facel]
	punpckldq mm2,mm2	        /* spread to both halves */
	movq  [esp + _iq], mm2	        /* iq =facel*charge[ii] */

	mov   edx, [ebp + _type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + _ntype]
	shl   edx, 1
	mov   [esp + _ntia], edx

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	
	pfadd mm0, [eax + ebx*4]        /* ix = shX + posX (and iy too) */
	movd  mm3, [eax + ebx*4 + 8]    /* cant use direct memory add for 4 bytes (iz) */
	mov   [esp + _ii3], ebx
	pfadd mm1, mm3
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	
				
	/* clear total potential and i forces */
	pxor  mm7,mm7
	movq  [esp + _vctot],  mm7
	movq  [esp + _vnbtot], mm7
	movq  [esp + _fix],    mm7
	movd  [esp + _fiz],    mm7

	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */

	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]	
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr], eax     /* pointer to jjnr[nj0] */
	sub   edx,  2
	mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i3300_unroll_loop
	jmp   .i3300_finish_inner
.i3300_unroll_loop:
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */
	
	mov ecx, [ebp + _charge]        /* base of charge[] */
	movq mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]	         /* charge[jnr1] */
        punpckldq mm3, [ecx + ebx*4]     /* move charge 2 to high part of mm3 */
	pfmul mm3,mm5		         /* mm3 now has qq for both particles */

	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	mov ecx, [ecx + ebx*4]           /* type [jnr2] */

	mov esi, [ebp + _nbfp]		/* base of nbfp */ 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	add ecx, [esp + _ntia]

	movq mm5, [esi + edx*4]		/* mm5 = 1st c6 / c12 */		
	movq mm7, [esi + ecx*4]		/* mm7 = 2nd c6 / c12 */	
	movq mm6,mm5			
	punpckldq mm5,mm7		/* mm5 = 1st c6 / 2nd c6 */
	punpckhdq mm6,mm7		/* mm6 = 1st c12 / 2nd c12 */
	movq [esp + _c6], mm5
	movq [esp + _c12], mm6

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrsqrt mm0, mm4	         /* lookup inverse square root seed */
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs. */
        movq mm2,mm0	        	/* amd 3dnow N-R iteration to get full precision. */
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r. */
	/* do potential and fscal */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */
	
	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */

	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, mm3	/* vcoul=qq*VV */
	pfmul mm3, mm7	/* fijC=FF*qq */ 

	/* at this point mm5 contains vcoul and mm3 fijC */
	/* increment vcoul - then we can get rid of mm5 */
	/* update vctot */
	pfadd mm5, [esp + _vctot]      /* add the earlier value */
	movq [esp + _vctot], mm5       /* store the sum */      

	/* dispersion table */
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + _n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */	

	movq mm4, [esp + _c6]
	pfmul mm7, mm4	/* fijD */
	pfmul mm5, mm4	/* vnb6 */           
	pfadd mm3, mm7	/* add to fscal */

	/* update vnbtot to release mm5! */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* repulsion table */
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]
	mov ecx, [esp + _n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 32]
	punpckldq mm5, [edx + ecx*4 + 36]
	punpckldq mm6, [edx + ecx*4 + 40]
	punpckldq mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	movq mm6, [esp + _c12]
	pfmul mm7, mm6	/* fijR */
	pfmul mm5, mm6	/* vnb12 */
	pfadd mm3, mm7	/* total fscal fijC+ fijD+ fijR */

	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + _tsc]
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	/* update vnbtot */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i3300_finish_inner
	jmp   .i3300_unroll_loop
.i3300_finish_inner:	
	and [esp + _innerk],  1
	jnz  .i3300_single_inner
	jmp  .i3300_updateouterdata		
.i3300_single_inner:
	/* a single j particle iteration here - compare with the unrolled code for comments. */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov ecx, [ebp + _charge]
	movd mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	/* mm3=qq */

	mov esi, [ebp + _nbfp]
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	shl edx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	movd mm5, [esi + edx*4]		/* mm5 = 1st c6 */ 		
	movq [esp + _c6], mm5
	movd mm5, [esi + edx*4 + 4]	/* mm5 = 1st c12 */ 		
	movq [esp + _c12], mm5

	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm0=rsq */
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	/* mm1=invsqrt */
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r. */

	/* calculate potentials and scalar force */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */

	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */

	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, mm3	/* vcoul=qq*VV */
	pfmul mm3, mm7	/* fijC=FF*qq */ 

	/* at this point mm5 contains vcoul and mm3 fijC */
	/* increment vcoul - then we can get rid of mm5 */
	/* update vctot */
	pfadd mm5, [esp + _vctot]      /* add the earlier value */
	movq [esp + _vctot], mm5       /* store the sum */      
	
	/* dispersion table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */	

	movq mm4, [esp + _c6]
	pfmul mm7, mm4	/* fijD */
	pfmul mm5, mm4	/* vnb6 */           
	pfadd mm3, mm7	/* add to fscal */

	/* update vnbtot to release mm5! */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* repulsion table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	movq mm6, [esp + _c12]
	pfmul mm7, mm6	/* fijR */
	pfmul mm5, mm6	/* vnb12 */
	pfadd mm3, mm7	/* total fscal fijC+ fijD+ fijR */

	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + _tsc]
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	/* update vnbtot */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i3300_updateouterdata:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */

	movq  mm7, [esp + _vnbtot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vnb[gid] */

	/* finish if last */
	mov   ecx, [ebp + _nri]
	dec ecx
	jecxz .i3300_end
	/* not last, iterate once more! */
	mov [ebp + _nri], ecx
	jmp .i3300_outer
.i3300_end:
	femms
	add esp, 132
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret





.globl inl3310_3dnow
	.type inl3310_3dnow,@function
inl3310_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56			
.equ		_type,		60
.equ		_ntype,		64
.equ		_nbfp,		68	
.equ		_Vnb,		72
.equ		_tabscale,	76
.equ		_VFtab,		80
.equ		_nsatoms,	84		
	/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_shX,		8
.equ		_shY,		12 
.equ		_shZ,		16	
.equ		_ix,		20
.equ		_iy,		24
.equ		_iz,		28	
.equ		_iq,		32  
.equ		_vctot,		40  
.equ		_vnbtot,	48  
.equ		_c6,		56  
.equ		_c12,		64  
.equ		_two,		72  
.equ		_n1,		80  
.equ		_tsc,		88  
.equ		_ntia,		96	
.equ		_innerjjnr0,	100
.equ		_innerk0,	104	
.equ		_innerjjnr,	108
.equ		_innerk,	112	
.equ		_fix,		116
.equ		_fiy,		120
.equ		_fiz,		124
.equ		_dx1,		128
.equ		_dy1,		132
.equ		_dz1,		136
.equ		_dx2,		140
.equ		_dy2,		144
.equ		_dz2,		148								
.equ		_nsvdwc,	152
.equ		_nscoul,	156
.equ		_nsvdw,		160
.equ		_solnr,		164		
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 168		/* local stack space */
	femms
	movq  mm0, [mm_two]
	movd  mm3, [ebp + _tabscale]
	movq  [esp + _two],    mm0
	punpckldq mm3,mm3
	movq  [esp + _tsc], mm3	
	/* assume we have at least one i particle - start directly */		
.i3310_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm0, [eax + ebx*4]	/* move shX/shY to mm0 and shZ to mm1 */
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + _shX], mm0
	movd  [esp + _shZ], mm1

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   eax, [ebp + _nsatoms]
	add   [ebp + _nsatoms],  12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + _nsvdwc], edx
	mov   [esp + _nscoul], eax
	mov   [esp + _nsvdw], ecx
		
	/* clear potential */
	pxor  mm7,mm7
	movq  [esp + _vctot],  mm7
	movq  [esp + _vnbtot], mm7
	mov   [esp + _solnr],  ebx
	
	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr0], eax     /* pointer to jjnr[nj0] */

	mov   [esp + _innerk0], edx        /* number of innerloop atoms */
	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]
	
	mov   ecx, [esp + _nsvdwc]
	cmp   ecx,  0
	jnz   .i3310_mno_vdwc
	jmp   .i3310_testcoul
.i3310_mno_vdwc:
	mov   ebx,  [esp + _solnr]
	inc   dword ptr [esp + _solnr]
	mov   edx, [ebp + _charge]
	movd  mm2, [edx + ebx*4]        /* mm2=charge[ii] */
	pfmul mm2, [ebp + _facel]
	punpckldq mm2,mm2	        /* spread to both halves */
	movq  [esp + _iq], mm2	        /* iq =facel*charge[ii] */

	mov   edx, [ebp + _type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + _ntype]
	shl   edx, 1
	mov   [esp + _ntia], edx	

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	mov   [esp + _ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + _shX]
	pfadd mm1, [esp + _shZ]
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	

	/* clear forces */
	pxor  mm7,mm7
	movq  [esp + _fix],   mm7
	movd  [esp + _fiz],   mm7

	mov   ecx, [esp + _innerjjnr0]
	mov   [esp + _innerjjnr], ecx
	mov   edx, [esp + _innerk0]
        sub   edx,  2
        mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i3310_unroll_vdwc_loop
	jmp   .i3310_finish_vdwc_inner
.i3310_unroll_vdwc_loop:	
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */
	
	mov ecx, [ebp + _charge]        /* base of charge[] */
	movq mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]	         /* charge[jnr1] */
        punpckldq mm3, [ecx + ebx*4]     /* move charge 2 to high part of mm3 */
	pfmul mm3,mm5		         /* mm3 now has qq for both particles */

	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	mov ecx, [ecx + ebx*4]           /* type [jnr2] */

	mov esi, [ebp + _nbfp]		/* base of nbfp */ 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	add ecx, [esp + _ntia]

	movq mm5, [esi + edx*4]		/* mm5 = 1st c6 / c12 */		
	movq mm7, [esi + ecx*4]		/* mm7 = 2nd c6 / c12 */	
	movq mm6,mm5			
	punpckldq mm5,mm7		/* mm5 = 1st c6 / 2nd c6 */
	punpckhdq mm6,mm7		/* mm6 = 1st c12 / 2nd c12 */
	movq [esp + _c6], mm5
	movq [esp + _c12], mm6

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrsqrt mm0, mm4	         /* lookup inverse square root seed */
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs. */
        movq mm2,mm0	        	/* amd 3dnow N-R iteration to get full precision. */
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r. */
	/* do potential and fscal */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */
	
	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */

	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, mm3	/* vcoul=qq*VV */
	pfmul mm3, mm7	/* fijC=FF*qq */ 

	/* at this point mm5 contains vcoul and mm3 fijC */
	/* increment vcoul - then we can get rid of mm5 */
	/* update vctot */
	pfadd mm5, [esp + _vctot]      /* add the earlier value */
	movq [esp + _vctot], mm5       /* store the sum */      

	/* dispersion table */
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + _n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */	

	movq mm4, [esp + _c6]
	pfmul mm7, mm4	/* fijD */
	pfmul mm5, mm4	/* vnb6 */           
	pfadd mm3, mm7	/* add to fscal */

	/* update vnbtot to release mm5! */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* repulsion table */
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]
	mov ecx, [esp + _n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 32]
	punpckldq mm5, [edx + ecx*4 + 36]
	punpckldq mm6, [edx + ecx*4 + 40]
	punpckldq mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	movq mm6, [esp + _c12]
	pfmul mm7, mm6	/* fijR */
	pfmul mm5, mm6	/* vnb12 */
	pfadd mm3, mm7	/* total fscal fijC+ fijD+ fijR */

	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + _tsc]
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	/* update vnbtot */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i3310_finish_vdwc_inner
	jmp   .i3310_unroll_vdwc_loop
.i3310_finish_vdwc_inner:	
	and [esp + _innerk],  1
	jnz  .i3310_single_vdwc_inner
	jmp  .i3310_updateouterdata_vdwc		
.i3310_single_vdwc_inner:	
	/* a single j particle iteration here - compare with the unrolled code for comments. */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov ecx, [ebp + _charge]
	movd mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	/* mm3=qq */

	mov esi, [ebp + _nbfp]
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	shl edx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	movd mm5, [esi + edx*4]		/* mm5 = 1st c6 */ 		
	movq [esp + _c6], mm5
	movd mm5, [esi + edx*4 + 4]	/* mm5 = 1st c12 */ 		
	movq [esp + _c12], mm5

	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm0=rsq */
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	/* mm1=invsqrt */
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r. */

	/* calculate potentials and scalar force */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */

	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */

	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, mm3	/* vcoul=qq*VV */
	pfmul mm3, mm7	/* fijC=FF*qq */ 

	/* at this point mm5 contains vcoul and mm3 fijC */
	/* increment vcoul - then we can get rid of mm5 */
	/* update vctot */
	pfadd mm5, [esp + _vctot]      /* add the earlier value */
	movq [esp + _vctot], mm5       /* store the sum */      
	
	/* dispersion table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */	

	movq mm4, [esp + _c6]
	pfmul mm7, mm4	/* fijD */
	pfmul mm5, mm4	/* vnb6 */           
	pfadd mm3, mm7	/* add to fscal */

	/* update vnbtot to release mm5! */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* repulsion table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	movq mm6, [esp + _c12]
	pfmul mm7, mm6	/* fijR */
	pfmul mm5, mm6	/* vnb12 */
	pfadd mm3, mm7	/* total fscal fijC+ fijD+ fijR */

	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + _tsc]
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	/* update vnbtot */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i3310_updateouterdata_vdwc:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	/* loop back to mno */
	dec dword ptr [esp + _nsvdwc]
	jz  .i3310_testcoul
	jmp .i3310_mno_vdwc
.i3310_testcoul:	
	mov  ecx, [esp + _nscoul]
	cmp  ecx,  0
	jnz  .i3310_mno_coul
	jmp  .i3310_testvdw
.i3310_mno_coul:
	mov   ebx,  [esp + _solnr]
	inc   dword ptr [esp + _solnr]
	mov   edx, [ebp + _charge]
	movd  mm2, [edx + ebx*4]        /* mm2=charge[ii] */
	pfmul mm2, [ebp + _facel]
	punpckldq mm2,mm2	        /* spread to both halves */
	movq  [esp + _iq], mm2	        /* iq =facel*charge[ii] */

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	mov   [esp + _ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + _shX]
	pfadd mm1, [esp + _shZ]
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	

	/* clear forces */
	pxor  mm7,mm7
	movq  [esp + _fix],   mm7
	movd  [esp + _fiz],   mm7

	mov   ecx, [esp + _innerjjnr0]
	mov   [esp + _innerjjnr], ecx
	mov   edx, [esp + _innerk0]
        sub   edx,  2
        mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i3310_unroll_coul_loop
	jmp   .i3310_finish_coul_inner
.i3310_unroll_coul_loop:	
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */
	
	mov ecx, [ebp + _charge]        /* base of charge[] */
	movq mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]	         /* charge[jnr1] */
        punpckldq mm3, [ecx + ebx*4]     /* move charge 2 to high part of mm3 */
	pfmul mm3,mm5		         /* mm3 now has qq for both particles */

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 /* dr = ir - jr */ 
	pfsubr mm7,mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6,mm6	                 /* square dx,dy,dz */
	pfmul mm7,mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrsqrt mm0, mm4	         /* lookup inverse square root seed */
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	/* now 4 has rsq and 0 the seed for both pairs. */
        movq mm2,mm0	        	/* amd 3dnow N-R iteration to get full precision. */
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r. */
	/* do potential and fscal */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */
	
	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	/* coulomb table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */

	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, mm3	/* vcoul=qq*VV */
	pfmul mm3, mm7	/* fijC=FF*qq */ 

	/* at this point mm5 contains vcoul and mm3 fijC */
	/* increment vcoul - then we can get rid of mm5 */
	/* update vctot */
	pfadd mm5, [esp + _vctot]      /* add the earlier value */
	movq [esp + _vctot], mm5       /* store the sum */      

	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3
	pfmul mm1, [esp + _tsc]	
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i3310_finish_coul_inner
	jmp   .i3310_unroll_coul_loop
.i3310_finish_coul_inner:	
	and [esp + _innerk],  1
	jnz  .i3310_single_coul_inner
	jmp  .i3310_updateouterdata_coul		
.i3310_single_coul_inner:	
	/* a single j particle iteration here - compare with the unrolled code for comments. */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov ecx, [ebp + _charge]
	movd mm5, [esp + _iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	/* mm3=qq */

	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm0=rsq */
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	/* mm1=invsqrt */
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r. */

	/* calculate potentials and scalar force */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */

	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */

	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, mm3	/* vcoul=qq*VV */
	pfmul mm3, mm7	/* fijC=FF*qq */ 

	/* at this point mm5 contains vcoul and mm3 fijC */
	/* increment vcoul - then we can get rid of mm5 */
	/* update vctot */
	pfadd mm5, [esp + _vctot]      /* add the earlier value */
	movq [esp + _vctot], mm5       /* store the sum */      
	
	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + _tsc]
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i3310_updateouterdata_coul:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	/* loop back to mno */
	dec dword ptr [esp + _nscoul]
	jz  .i3310_testvdw
	jmp .i3310_mno_coul
.i3310_testvdw:	
	mov  ecx, [esp + _nsvdw]
	cmp  ecx,  0
	jnz  .i3310_mno_vdw
	jmp  .i3310_last_mno
.i3310_mno_vdw:
	mov   ebx,  [esp + _solnr]
	inc   dword ptr [esp + _solnr]

	mov   edx, [ebp + _type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + _ntype]
	shl   edx, 1
	mov   [esp + _ntia], edx	

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	mov   [esp + _ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + _shX]
	pfadd mm1, [esp + _shZ]
	movq  [esp + _ix], mm0	
	movd  [esp + _iz], mm1	

	/* clear forces */
	pxor  mm7,mm7
	movq  [esp + _fix],   mm7
	movd  [esp + _fiz],   mm7

	mov   ecx, [esp + _innerjjnr0]
	mov   [esp + _innerjjnr], ecx
	mov   edx, [esp + _innerk0]
        sub   edx,  2
        mov   [esp + _innerk], edx        /* number of innerloop atoms */
	jge   .i3310_unroll_vdw_loop
	jmp   .i3310_finish_vdw_inner
.i3310_unroll_vdw_loop:	
	/* paired innerloop starts here */
	mov   ecx, [esp + _innerjjnr]     /* pointer to jjnr[k] */
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             /* eax/ebx=jnr */
	add   [esp + _innerjjnr],  8 /* advance pointer (unrolled 2) */
	prefetch [ecx + 16]	         /* prefetch data - trial and error says 16 is best */
	
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	mov ecx, [ecx + ebx*4]           /* type [jnr2] */

	mov esi, [ebp + _nbfp]		/* base of nbfp */ 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	add ecx, [esp + _ntia]

	movq mm5, [esi + edx*4]		/* mm5 = 1st c6 / c12 */		
	movq mm7, [esi + ecx*4]		/* mm7 = 2nd c6 / c12 */	
	movq mm6, mm5			
	punpckldq mm5, mm7		/* mm5 = 1st c6 / 2nd c6 */
	punpckhdq mm6, mm7		/* mm6 = 1st c12 / 2nd c12 */
	movq [esp + _c6], mm5
	movq [esp + _c12], mm6

	lea   eax, [eax + eax*2]         /* replace jnr with j3 */
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + _pos]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]	 	
	movq  mm4, [esi + eax*4]         /* fetch first j coordinates */
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         /* dr = ir - jr */ 
	pfsubr mm5,mm1
	movq  [esp + _dx1], mm4	         /* store dr */
	movd  [esp + _dz1], mm5
	pfmul mm4,mm4	                 /* square dx,dy,dz */		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm4, mm5		         /* first rsq in lower mm4 */

	movq  mm6, [esi + ebx*4]         /* fetch second j coordinates */ 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6, mm0	                 /* dr = ir - jr */ 
	pfsubr mm7, mm1
	movq  [esp + _dx2], mm6	         /* store dr */
	movd  [esp + _dz2], mm7
	pfmul mm6, mm6	                 /* square dx,dy,dz */
	pfmul mm7, mm7
	pfacc mm6, mm7		         /* accumulate to get dx*dx+ dy*dy+ dz*dz */
	pfacc mm6, mm7	                 /* second rsq in lower mm6 */

        pfrsqrt mm0, mm4	         /* lookup inverse square root seed */
        pfrsqrt mm1, mm6
 

	punpckldq mm0, mm1
	punpckldq mm4, mm6        	/* now 4 has rsq and 0 the seed for both pairs. */
        movq mm2, mm0	        	/* amd 3dnow N-R iteration to get full precision. */
	pfmul mm0, mm0
        pfrsqit1 mm0, mm4				
        pfrcpit2 mm0, mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r. */
	/* do potential and fscal */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4, mm1
	movq [esp + _n1], mm4
	pi2fd mm4, mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 is n0 */
	
	movq mm2, mm1
	pfmul mm2, mm2	/* mm1 is eps, mm2 is eps2 */
	
	mov edx, [ebp + _VFtab]
	/* dispersion table */
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]	
	shl ecx, 2
	/* load all the table values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */	

	movq mm4, [esp + _c6]
	pfmul mm7, mm4	/* fijD */
	pfmul mm5, mm4	/* vnb6 */           
	movq mm3, mm7	/* add to fscal */

	/* update vnbtot to release mm5! */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* repulsion table */
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + _n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	movq mm6, [esp + _c12]
	pfmul mm7, mm6	/* fijR */Toggle
	pfmul mm5, mm6	/* vnb12 */
	pfadd mm3, mm7	/* total fscal fijD+ fijR */

	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm1, [esp + _tsc]
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	prefetchw [esp + _dx1]	/* prefetch i forces to cache */

	/* spread fscalar to both positions */
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	/* calc vector force */
	prefetchw [edi + eax*4]	/* prefetch the 1st faction to cache */
	movq mm2,  [esp + _dx1]	/* fetch dr */
	movd mm3,  [esp + _dz1]

	/* update vnbtot */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	prefetchw [edi + ebx*4]	/* prefetch the 2nd faction to cache */
	pfmul mm2, mm0		/* mult by fs */ 
	pfmul mm3, mm0

	movq mm4,  [esp + _dx2] 	/* fetch dr */
	movd mm5,  [esp + _dz2]
	pfmul mm4, mm1   	/* mult by fs */ 
	pfmul mm5, mm1
	/* update i forces */

	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j forces */

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	/* should we do one more iteration? */
	sub   [esp + _innerk],  2
	jl    .i3310_finish_vdw_inner
	jmp   .i3310_unroll_vdw_loop
.i3310_finish_vdw_inner:	
	and [esp + _innerk],  1
	jnz  .i3310_single_vdw_inner
	jmp  .i3310_updateouterdata_vdw		
.i3310_single_vdw_inner:	
	/* a single j particle iteration here - compare with the unrolled code for comments. */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */

	mov esi, [ebp + _nbfp]
	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr1] */
	shl edx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	movd mm5, [esi + edx*4]		/* mm5 = 1st c6 */ 		
	movq [esp + _c6], mm5
	movd mm5, [esi + edx*4 + 4]	/* mm5 = 1st c12 */ 		
	movq [esp + _c12], mm5

	mov   esi, [ebp + _pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + _ix]
	movd  mm1, [esp + _iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + _dx1], mm4
	pfmul mm4,mm4
	movd  [esp + _dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		/* mm0=rsq */
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	/* mm1=invsqrt */
	pfmul mm4, mm0
	movq mm1, mm4
	/* mm0 is invsqrt, and mm1 r. */

	/* calculate potentials and scalar force */
	pfmul mm1, [esp + _tsc]	/* mm1=rt */
	pf2iw mm4,mm1
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   /* now mm1 is eps and mm4 n0. */

	movq mm2,mm1
	pfmul mm2,mm2	/* mm1 is eps, mm2 is eps2 */
	
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]	
	shl ecx, 2
	/* dispersion table
	 * load all the table values we need
	 */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */	

	movq mm4, [esp + _c6]
	pfmul mm7, mm4	/* fijD */
	pfmul mm5, mm4	/* vnb6 */           
	movq mm3, mm7	/* add to fscal */

	/* update vnbtot to release mm5! */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* repulsion table
	 * load all the table values we need
	 */
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm1  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	movq mm6, [esp + _c12]
	pfmul mm7, mm6	/* fijR */
	pfmul mm5, mm6	/* vnb12 */
	pfadd mm3, mm7	/* total fscal fijC+ fijD+ fijR */

	/* change sign of mm3 */
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + _tsc]
 	pfmul mm0, mm1        /* mm0 is total fscal now */	

	/* update vnbtot */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* spread fscalar to both positions */
	punpckldq mm0,mm0
	/* calc vectorial force */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm2,  [esp + _dx1]
	movd mm3,  [esp + _dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	/* update i particle force */
	movq mm0,  [esp + _fix]
	movd mm1,  [esp + _fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + _fix], mm0
	movd [esp + _fiz], mm1
	/* update j particle force */
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	/* done! */
.i3310_updateouterdata_vdw:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment i force */
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fix]
	pfadd mm7, [esp + _fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	/* loop back to mno */
	dec dword ptr [esp + _nsvdw]
	jz  .i3310_last_mno
	jmp .i3310_mno_vdw
	
.i3310_last_mno:	
	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */

	movq  mm7, [esp + _vnbtot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */
	/* finish if last */
	mov   ecx, [ebp + _nri]
	dec ecx
	jecxz .i3310_end
	/* not last, iterate once more! */
	mov [ebp + _nri], ecx
	jmp .i3310_outer
.i3310_end:
	femms
	add esp, 168
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret


.globl inl3320_3dnow
	.type inl3320_3dnow,@function
inl3320_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56			
.equ		_type,		60
.equ		_ntype,		64
.equ		_nbfp,		68	
.equ		_Vnb,		72
.equ		_tabscale,	76
.equ		_VFtab,		80
			/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_ixO,		8
.equ		_iyO,		12
.equ		_izO,		16	
.equ		_ixH,		20  
.equ		_iyH,		28  
.equ		_izH,		36  
.equ		_iqO,		44  
.equ		_iqH,		52  
.equ		_qqO,		60  
.equ		_qqH,		68  
.equ		_vctot,		76  
.equ		_vnbtot,	84  
.equ		_c6,		92  
.equ		_c12,		100 
.equ		_two,		108 
.equ		_n1,		116 
.equ		_tsc,		124 
.equ		_ntia,		132 
.equ		_innerjjnr,	140
.equ		_innerk,	144	
.equ		_fixO,		148
.equ		_fiyO,		152
.equ		_fizO,		156
.equ		_fixH,		160 
.equ		_fiyH,		168 
.equ		_fizH,		176 
.equ		_dxO,		184
.equ		_dyO,		188
.equ		_dzO,		192
.equ		_dxH,		196 
.equ		_dyH,		204 
.equ		_dzH,		212 
.equ		_tmprsqH,	220 
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 228		/* local stack space */
	femms

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   edx, [ebp + _charge]
	movd  mm1, [ebp + _facel]
	movd  mm2, [edx + ebx*4]        /* mm2=charge[i.i0] */
	pfmul mm2, mm1		
	movq  [esp + _iqO], mm2	        /* iqO = facel*charge[ii] */
	
	movd  mm2, [edx + ebx*4 + 4]    /* mm2=charge[i.i0+1] */
	pfmul mm2, mm1
	punpckldq mm2,mm2	        /* spread to both halves */
	movq  [esp + _iqH], mm2	        /* iqH = facel*charge[i.i0+1] */

	mov   edx, [ebp + _type] 	
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1		        
	imul  ecx, [ebp + _ntype]      /* ecx = ntia = 2*ntype*type[i.i0] */ 
	mov   [esp + _ntia], ecx
	 	
	movq  mm3, [mm_two]
	movq  mm4, [ebp + _tabscale]
	punpckldq mm4,mm4	        /* spread to both halves */
	movq  [esp + _two],    mm3
	movq  [esp + _tsc], mm4	      
	/* assume we have at least one i particle - start directly */	 
.i3320_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm5, [eax + ebx*4]	/* move shX/shY to mm5 and shZ to mm6. */
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        /* also expand shX,Y,Z in mm0--mm2. */
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */
	
	pfadd mm5, [eax + ebx*4]        /* ix = shX + posX (and iy too) */
	movd  mm7, [eax + ebx*4 + 8]    /* cant use direct memory add for 4 bytes (iz) */
	mov   [esp + _ii3], ebx	        /* (use mm7 as temp. storage for iz.) */
	pfadd mm6, mm7
	movq  [esp + _ixO], mm5	
	movq  [esp + _izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] /* coords of H1 in low mm3-mm5, H2 in high */
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + _ixH], mm0	
	movq [esp + _iyH], mm1	
	movq [esp + _izH], mm2	
					
	/* clear vctot and i forces */
	pxor  mm7,mm7
	movq  [esp + _vctot], mm7
	movq  [esp + _vnbtot], mm7
	movq  [esp + _fixO],   mm7
	movd  [esp + _fizO],   mm7
	movq  [esp + _fixH],   mm7
	movq  [esp + _fiyH],   mm7
	movq  [esp + _fizH],   mm7

	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */
	mov   [esp + _innerk], edx        

	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]	
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr], eax     /* pointer to jjnr[nj0] */
.i3320_inner_loop:	
	/* a single j particle iteration here - compare with the unrolled code for comments. */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */
        add   [esp + _innerjjnr],  4 /* advance pointer */
	prefetch [ecx + 16]	   /* prefetch data - trial and error says 16 is best */

	mov ecx, [ebp + _charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + _iqO]
	pfmul mm7, [esp + _iqH]	/* mm6=qqO, mm7=qqH */
	movd [esp + _qqO], mm6
	movq [esp + _qqH], mm7

	mov ecx, [ebp + _type]
	mov edx, [ecx + eax*4]        	 /* type [jnr] */
	mov ecx, [ebp + _nbfp]
	shl edx, 1
	add edx, [esp + _ntia]	         /* tja = ntia + 2*type */
	movd mm5, [ecx + edx*4]		/* mm5 = 1st c6 */ 		
	movq [esp + _c6], mm5
	movd mm5, [ecx + edx*4 + 4]	/* mm5 = 1st c12 */ 		
	movq [esp + _c12], mm5	
			
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */
	movq [esp + _tmprsqH], mm3
	
        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */

	pfmul mm0, mm1		/* mm0=r */

	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */

	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqO]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqO]	/* fijC=qq*FF */

	/* update vctot directly, use mm3 for fscal sum. */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	movq mm3, mm7
	
	/* dispersion table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */	

	movq mm4, [esp + _c6]
	pfmul mm7, mm4	/* fijD */
	pfmul mm5, mm4	/* vnb6 */           
	pfadd mm3, mm7	/* add to fscal */ 

	/* update vnbtot to release mm5! */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* repulsion table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	movq mm6, [esp + _c12]
	pfmul mm7, mm6	/* fijR */
	pfmul mm5, mm6	/* vnb12 */
	pfadd mm3, mm7	/* total fscal fijC+ fijD+ fijR */

	/* change sign of fscal and multiply with rinv */ 
        pxor mm0,mm0
	pfsubr mm3, mm0	
	pfmul mm3, [esp + _tsc]
 	pfmul mm3, mm1        /* mm3 is total fscal (for the oxygen) now */
	
	/* update vnbtot */ 
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      
	
	/* Ready with the oxygen - potential is updated, fscal is in mm3. */
	/* now do the two hydrogens. */
	movq mm0, [esp + _tmprsqH] /* mm0=r */sqH

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	/* seeds are in mm1 now, and rsq in mm0. */

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	
	pfmul mm0,mm1		/* mm0=r */
	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqH]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqH]	/* fijC=qq*FF */
	/* update vctot */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	
	/* change sign of fijC and multiply by rinv */
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + _tsc]
 	pfmul mm4, mm1        /* mm4 is total fscal (for the hydrogens) now */	
	
	/* spread oxygen fscalar to both positions */
	punpckldq mm3,mm3
	/* calc vectorial force for O */
	prefetchw [edi + eax*4]	/* prefetch faction to cache */ 
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4
	
	/* pack j forces from H in the same form as the oxygen force. */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle. */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	
	/*  done  - one more? */
	dec dword ptr [esp + _innerk]
	jz  .i3320_updateouterdata
	jmp .i3320_inner_loop
.i3320_updateouterdata:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment iO force */ 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + _fixH]
	movq  mm3, [esp + _fiyH]
	movq  mm1, [esp + _fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	/* mm0(l)=fxH1, mm0(h)=fyH1 */
	punpckhdq mm2, mm3	/* mm2(l)=fxH2, mm2(h)=fyH2 */
	movq mm3, mm1
	pswapd mm3,mm3		
	/* mm1 is fzH1 */
	/* mm3 is fzH2 */
	
	movq  mm6, [edi + ecx*4 + 12]       /* increment iH1 force */ 
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       /* increment iH2 force */
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */
	
	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */

	movq  mm7, [esp + _vnbtot]     
	pfacc mm7,mm7	              /* same for Vnb */
	
	mov   eax, [ebp + _Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vnb[gid] */
	/* finish if last */
	dec dword ptr [ebp + _nri]
	jz  .i3320_end
	/* not last, iterate once more! */
	jmp .i3320_outer
.i3320_end:
	femms
	add esp, 228
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret

	

.globl inl3330_3dnow
	.type inl3330_3dnow,@function
inl3330_3dnow:	
.equ		_nri,		8
.equ		_iinr,		12
.equ		_jindex,	16
.equ		_jjnr,		20
.equ		_shift,		24
.equ		_shiftvec,	28
.equ		_fshift,	32
.equ		_gid,		36
.equ		_pos,		40		
.equ		_faction,	44
.equ		_charge,	48
.equ		_facel,		52
.equ		_Vc,		56			
.equ		_type,		60
.equ		_ntype,		64
.equ		_nbfp,		68	
.equ		_Vnb,		72
.equ		_tabscale,	76
.equ		_VFtab,		80
			/* stack offsets for local variables */
.equ		_is3,		0
.equ		_ii3,		4
.equ		_ixO,		8
.equ		_iyO,		12
.equ		_izO,		16	
.equ		_ixH,		20  
.equ		_iyH,		28  
.equ		_izH,		36  
.equ		_qqOO,		44  
.equ		_qqOH,		52  
.equ		_qqHH,		60  
.equ		_c6,		68  
.equ		_c12,		76  
.equ		_two,		84  
.equ		_n1,		92  
.equ		_tsc,		100 
.equ		_vctot,		108 
.equ		_vnbtot,	116 
.equ		_innerjjnr,	124
.equ		_innerk,	128	
.equ		_fixO,		132
.equ		_fiyO,		136
.equ		_fizO,		140
.equ		_fixH,		144 
.equ		_fiyH,		152 
.equ		_fizH,		160 
.equ		_dxO,		168
.equ		_dyO,		172
.equ		_dzO,		176
.equ		_dxH,		180  
.equ		_dyH,		188  
.equ		_dzH,		196  
.equ		_tmprsqH,	204  
	push ebp
	mov ebp,esp	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 212		/* local stack space */
	femms
	/* assume we have at least one i particle - start directly */	

	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	mov   ebx, [ecx]	        /* ebx=ii */

	mov   edx, [ebp + _charge]
	movd  mm1, [ebp + _facel]	/* mm1=facel */
	movd  mm2, [edx + ebx*4]        /* mm2=charge[i.i0] (O) */
	movd  mm3, [edx + ebx*4 + 4]    /* mm2=charge[i.i0+1] (H) */ 
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			/* mm4=qqOO*facel */
	pfmul mm5, mm3			/* mm5=qqOH*facel */
	pfmul mm6, mm3			/* mm6=qqHH*facel */
	punpckldq mm5,mm5	        /* spread to both halves */
	punpckldq mm6,mm6	        /* spread to both halves */
	movq  [esp + _qqOO], mm4
	movq  [esp + _qqOH], mm5
	movq  [esp + _qqHH], mm6
	mov   edx, [ebp + _type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [ebp + _ntype]
	add   edx, ecx
	mov   eax, [ebp + _nbfp]
	movd  mm0, [eax + edx*4]
	movd  mm1, [eax + edx*4 + 4]
	movq  [esp + _c6], mm0
	movq  [esp + _c12], mm1
	movq  mm2, [mm_two]
	movq  [esp + _two], mm2
	movd  mm3, [ebp + _tabscale]
	punpckldq mm3,mm3
	movq  [esp + _tsc], mm3
.i3330_outer:
	mov   eax, [ebp + _shift]      /* eax = pointer into shift[] */
	mov   ebx, [eax]		/* ebx=shift[n] */
	add   [ebp + _shift],  4  /* advance pointer one step */
	
	lea   ebx, [ebx + ebx*2]        /* ebx=3*is */
	mov   [esp + _is3],ebx    	/* store is3 */

	mov   eax, [ebp + _shiftvec]   /* eax = base of shiftvec[] */
	
	movq  mm5, [eax + ebx*4]	/* move shX/shY to mm5 and shZ to mm6. */
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        /* also expand shX,Y,Z in mm0--mm2. */
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + _iinr]       /* ecx = pointer into iinr[] */	
	add   [ebp + _iinr],  4   /* advance pointer */
	mov   ebx, [ecx]	        /* ebx=ii */

	lea   ebx, [ebx + ebx*2]	/* ebx = 3*ii=ii3 */
	mov   eax, [ebp + _pos]        /* eax = base of pos[] */

	pfadd mm5, [eax + ebx*4]        /* ix = shX + posX (and iy too) */
	movd  mm7, [eax + ebx*4 + 8]    /* cant use direct memory add for 4 bytes (iz) */
	mov   [esp + _ii3], ebx	        /* (use mm7 as temp. storage for iz.) */
	pfadd mm6, mm7
	movq  [esp + _ixO], mm5	
	movq  [esp + _izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] /* coords of H1 in low mm3-mm5, H2 in high */
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + _ixH], mm0	
	movq [esp + _iyH], mm1	
	movq [esp + _izH], mm2	

	/* clear vctot and i forces */
	pxor  mm7,mm7
	movq  [esp + _vctot], mm7
	movq  [esp + _vnbtot], mm7
	movq  [esp + _fixO],  mm7
	movq  [esp + _fizO],  mm7
	movq  [esp + _fixH],  mm7
	movq  [esp + _fiyH],  mm7
	movq  [esp + _fizH],  mm7

	mov   eax, [ebp + _jindex]
	mov   ecx, [eax]	         /* jindex[n] */
	mov   edx, [eax + 4]	         /* jindex[n+1] */
	add   [ebp + _jindex],  4
	sub   edx, ecx                   /* number of innerloop atoms */
	mov   [esp + _innerk], edx        

	mov   esi, [ebp + _pos]
	mov   edi, [ebp + _faction]	
	mov   eax, [ebp + _jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + _innerjjnr], eax     /* pointer to jjnr[nj0] */
.i3330_inner_loop:
	/* a single j particle iteration here - compare with the unrolled code for comments. */
	mov   eax, [esp + _innerjjnr]
	mov   eax, [eax]	/* eax=jnr offset */
        add   [esp + _innerjjnr],  4 /* advance pointer */

	lea   eax, [eax + eax*2]

	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */
	movq [esp + _tmprsqH], mm3

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */ OO
	pfmul mm0, mm1		/* mm0=rsq */ OO

	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */

	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2

	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqOO]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqOO]	/* fijC=qq*FF */

	/* update vctot directly, use mm3 for fscal sum. */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	movq mm3, mm7

	/* dispersion table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */	

	movq mm4, [esp + _c6]
	pfmul mm7, mm4	/* fijD */
	pfmul mm5, mm4	/* vnb6 */           
	pfadd mm3, mm7	/* add to fscal */ 

	/* update vnbtot to release mm5! */
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      

	/* repulsion table */
	/* load all the table values we need */
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */
	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */
	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	movq mm6, [esp + _c12]
	pfmul mm7, mm6	/* fijR */
	pfmul mm5, mm6	/* vnb12 */
	pfadd mm3, mm7	/* total fscal fijC+ fijD+ fijR */

	/* change sign of fscal and multiply with rinv */ 
        pxor mm0,mm0
	pfsubr mm3, mm0	
	pfmul mm3, [esp + _tsc]
 	pfmul mm3, mm1        /* mm3 is total fscal (for the oxygen) now */
	
	/* update vnbtot */ 
	pfadd mm5, [esp + _vnbtot]      /* add the earlier value */
	movq [esp + _vnbtot], mm5       /* store the sum */      
	
	/* Ready with the oxygen - potential is updated, fscal is in mm3.
	 * time for hydrogens!
         */
	
	movq mm0, [esp + _tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	/* seeds are in mm1 now, and rsq in mm0. */

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	
	pfmul mm0,mm1		/* mm0=r */
	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqOH]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqOH]	/* fijC=qq*FF */
	/* update vctot */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	
	/* change sign of fijC and multiply by rinv */
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + _tsc]
 	pfmul mm4, mm1        /* mm4 is total fscal (for the hydrogens) now */	
	
	/* spread oxygen fscalar to both positions */
	punpckldq mm3,mm3
	/* calc vectorial force for O */
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4
	
	/* pack j forces from H in the same form as the oxygen force. */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle. */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3

	/* interactions with j H1 */

	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */
	movq [esp + _tmprsqH], mm3

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	pfmul mm0, mm1		/* mm0=rsq */ 
	
	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */

	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2

	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqOH]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqOH]	/* fijC=qq*FF */

	/* update vctot directly, force is moved to mm3. */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	pxor mm3, mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + _tsc]
	pfmul mm3, mm1        /* mm3 is total fscal (for the oxygen) now */

	movq mm0, [esp + _tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	/* seeds are in mm1 now, and rsq in mm0. */

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	
	pfmul mm0,mm1		/* mm0=r */
	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqHH]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqHH]	/* fijC=qq*FF */
	/* update vctot */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	
	/* change sign of fijC and multiply by rinv */
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + _tsc]
 	pfmul mm4, mm1        /* mm4 is total fscal (for the hydrogens) now */		

	/* spread oxygen fscalar to both positions */
	punpckldq mm3,mm3
	/* calc vectorial force for O */
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4
	
	/* pack j forces from H in the same form as the oxygen force. */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle. */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4 + 12]
	movd mm3,  [edi + eax*4 + 20]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 12], mm2
	movd [edi + eax*4 + 20], mm3

	/* interactions with j H2 */
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	/* copy & expand to mm2-mm4 for the H interactions */
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	pfsubr mm0, [esp + _ixO]
	pfsubr mm1, [esp + _izO]
		
	movq  [esp + _dxO], mm0
	pfmul mm0,mm0
	movd  [esp + _dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		/* mm0=rsqO */
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  /* mm2-mm4 is jx-jz */
	pfsubr mm2, [esp + _ixH]
	pfsubr mm3, [esp + _iyH]
	pfsubr mm4, [esp + _izH] /* mm2-mm4 is dxH-dzH */
	
	movq [esp + _dxH], mm2
	movq [esp + _dyH], mm3
	movq [esp + _dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		/* mm3=rsqH */
	movq [esp + _tmprsqH], mm3

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	pfmul mm0, mm1

	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movd [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */

	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2

	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqOH]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqOH]	/* fijC=qq*FF */

	/* update vctot directly, use mm3 for fscal sum */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	pxor mm3,mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + _tsc]
 	pfmul mm3, mm1         /* mm3 is total fscal (for the oxygen) now */	

	movq mm0, [esp + _tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	/* seeds are in mm1 now, and rsq in mm0. */

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	/* mm1=invsqrt */
	
	pfmul mm0,mm1		/* mm0=r */
	pfmul mm0, [esp + _tsc]
	pf2iw mm4, mm0
	movq [esp + _n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   /* now mm0 is eps and mm4 n0 */
	movq  mm2, mm0
	pfmul mm2, mm2		/* mm0 is eps, mm2 eps2 */
	
	/* coulomb table */
	mov edx, [ebp + _VFtab]
	mov ecx, [esp + _n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	/* load all values we need */
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + _n1 + 4]/* mm5 = Fp */
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  /* mm6 = Geps */		
	pfmul mm7, mm2	/* mm7 = Heps2 */
	
	pfadd mm5, mm6
	pfadd mm5, mm7	/* mm5 = Fp */

	pfmul mm7, [esp + _two]	/* two*Heps2 */
	pfadd mm7, mm6
	pfadd mm7, mm5	/* mm7=FF */

	pfmul mm5, mm0  /* mm5=eps*Fp */
	pfadd mm5, mm4	/*  mm5= VV */

	pfmul mm5, [esp + _qqHH]	/* vcoul=qq*VV */
	pfmul mm7, [esp + _qqHH]	/* fijC=qq*FF */
	/* update vctot */
	pfadd mm5, [esp + _vctot]
	movq [esp + _vctot], mm5
	
	/* change sign of fijC and multiply by rinv */
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + _tsc]
 	pfmul mm4, mm1        /* mm4 is total fscal (for the hydrogens) now */	

	/* spread oxygen fscalar to both positions */
	punpckldq mm3,mm3
	/* calc vectorial force for O */
	movq mm0,  [esp + _dxO]
	movd mm1,  [esp + _dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	/* calc vectorial force for H's */
	movq mm5, [esp + _dxH]
	movq mm6, [esp + _dyH]
	movq mm7, [esp + _dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	/* update iO particle force */
	movq mm2,  [esp + _fixO]
	movd mm3,  [esp + _fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + _fixO], mm2
	movd [esp + _fizO], mm3

	/* update iH forces */
	movq mm2, [esp + _fixH]
	movq mm3, [esp + _fiyH]
	movq mm4, [esp + _fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + _fixH], mm2
	movq [esp + _fiyH], mm3
	movq [esp + _fizH], mm4	

	/* pack j forces from H in the same form as the oxygen force. */
	pfacc mm5, mm6		/* mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) */
	pfacc mm7, mm7		/* mm7(l)=fjz(H1+ h2) */
	
	pfadd mm0, mm5		/* add up total force on j particle. */ 
	pfadd mm1, mm7

	/* update j particle force */
	movq mm2,  [edi + eax*4 + 24]
	movd mm3,  [edi + eax*4 + 32]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 24], mm2
	movd [edi + eax*4 + 32], mm3
	
	/*  done  - one more? */
	dec dword ptr [esp + _innerk]
	jz  .i3330_updateouterdata
	jmp .i3330_inner_loop	
.i3330_updateouterdata:	
	mov   ecx, [esp + _ii3]

	movq  mm6, [edi + ecx*4]       /* increment iO force */ 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + _fixH]
	movq  mm3, [esp + _fiyH]
	movq  mm1, [esp + _fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	/* mm0(l)=fxH1, mm0(h)=fyH1 */
	punpckhdq mm2, mm3	/* mm2(l)=fxH2, mm2(h)=fyH2 */
	movq mm3, mm1
	pswapd mm3,mm3		
	/* mm1 is fzH1 */
	/* mm3 is fzH2 */

	movq  mm6, [edi + ecx*4 + 12]       /* increment iH1 force */ 
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       /* increment iH2 force */
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + _fshift]    /* increment fshift force */
	mov   edx, [esp + _is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + _fixO]
	pfadd mm7, [esp + _fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + _gid]      /* get group index for this i particle */
	mov   edx, [edx]
	add   [ebp + _gid],  4  /* advance pointer */

	movq  mm7, [esp + _vctot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */

	mov   eax, [ebp + _Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vc[gid] */

	movq  mm7, [esp + _vnbtot]     
	pfacc mm7,mm7	              /* get and sum the two parts of total potential */

	mov   eax, [ebp + _Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              /* increment vnbtot[gid] */
	/* finish if last */
	dec dword ptr [ebp + _nri]
	jz  .i3330_end
	/* not last, iterate once more! */
	jmp .i3330_outer
.i3330_end:
	femms
	add esp, 212
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	leave
	ret
 
