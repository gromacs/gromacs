;#
;# $Id$
;#
;# Gromacs 4.0                         Copyright (c) 1991-2003 
;# David van der Spoel, Erik Lindahl
;#
;# This program is free software; you can redistribute it and/or
;# modify it under the terms of the GNU General Public License
;# as published by the Free Software Foundation; either version 2
;# of the License, or (at your option) any later version.
;#
;# To help us fund GROMACS development, we humbly ask that you cite
;# the research papers on the package. Check out http://www.gromacs.org
;# 
;# And Hey:
;# Gnomes, ROck Monsters And Chili Sauce
;#

;# These files require GNU binutils 2.10 or later, since we
;# use intel syntax for portability, or a recent version 
;# of NASM that understands Extended 3DNow and SSE2 instructions.
;# (NASM is normally only used with MS Visual C++).
;# Since NASM and gnu as disagree on some definitions and use 
;# completely different preprocessing options I have to introduce a
;# trick: NASM uses ';' for comments, while gnu as uses '#' on x86.
;# Gnu as treats ';' as a line break, i.e. ignores it. This is the
;# reason why all comments need both symbols...
;# The source is written for GNU as, with intel syntax. When you use
;# NASM we redefine a couple of things. The false if-statement around 
;# the following code is seen by GNU as, but NASM doesn't see it, so 
;# the code inside is read by NASM but not gcc.

; .if 0    # block below only read by NASM
%define .section	section
%define .long		dd
%define .align		align
%define .globl		global
;# NASM only wants 'dword', not 'dword ptr'.
%define ptr
%macro .equiv 2
   %1 equ %2
%endmacro
; .endif                   # End of NASM-specific block
; .intel_syntax noprefix   # Line only read by gnu as




.globl nb_kernel030_ia32_3dnow
.globl _nb_kernel030_ia32_3dnow
nb_kernel030_ia32_3dnow:	
_nb_kernel030_ia32_3dnow:	
.equiv		nb030_p_nri,		8
.equiv		nb030_iinr,		12
.equiv		nb030_jindex,		16
.equiv		nb030_jjnr,		20
.equiv		nb030_shift,		24
.equiv		nb030_shiftvec,		28
.equiv		nb030_fshift,		32
.equiv		nb030_gid,		36
.equiv		nb030_pos,		40		
.equiv		nb030_faction,		44
.equiv		nb030_charge,		48
.equiv		nb030_p_facel,		52
.equiv		nb030_p_krf,		56	
.equiv		nb030_p_crf,		60	
.equiv		nb030_Vc,		64	
.equiv		nb030_type,		68
.equiv		nb030_p_ntype,		72
.equiv		nb030_vdwparam,		76	
.equiv		nb030_Vvdw,		80	
.equiv		nb030_p_tabscale,	84	
.equiv		nb030_VFtab,		88
.equiv		nb030_invsqrta,		92	
.equiv		nb030_dvda,		96
.equiv          nb030_p_gbtabscale,     100
.equiv          nb030_GBtab,            104
.equiv          nb030_p_nthreads,       108
.equiv          nb030_count,            112
.equiv 	        nb030_mtx,		116
.equiv 	        nb030_outeriter,	120
.equiv 	        nb030_inneriter,	124
.equiv 	        nb030_work,     	128
	;# stack offsets for local variables 
.equiv		nb030_is3,		0
.equiv		nb030_ii3,		4
.equiv		nb030_ix,		8
.equiv		nb030_iy,		12
.equiv		nb030_iz,		16
.equiv		nb030_Vvdwtot,		20 
.equiv		nb030_c6,		28 
.equiv		nb030_c12,		36 
.equiv		nb030_two,		44 
.equiv		nb030_n1,		52 
.equiv		nb030_tsc,		60 
.equiv		nb030_ntia,		68
.equiv		nb030_innerjjnr,	72
.equiv		nb030_innerk,		76		
.equiv		nb030_fix,		80
.equiv		nb030_fiy,		84
.equiv		nb030_fiz,		88
.equiv		nb030_dx1,		92
.equiv		nb030_dy1,		96
.equiv		nb030_dz1,		100
.equiv		nb030_dx2,		104
.equiv		nb030_dy2,		108
.equiv		nb030_dz2,		112						
.equiv          nb030_n,                116 ;# idx for outer loop
.equiv          nb030_nn1,              120 ;# number of outer iterations
.equiv		nb030_nri,		124
.equiv		nb030_ntype,		128
.equiv		nb030_nouter,		132
.equiv		nb030_ninner,		136

    	push ebp
    	mov ebp,esp
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 140		;# local stack space 
	femms
	;# move data to local stack  
	mov eax, 0x40000000
        mov [esp + nb030_two], eax
        mov [esp + nb030_two+4], eax

	mov ecx, [ebp + nb030_p_nri]
	mov edx, [ebp + nb030_p_ntype]
	mov esi, [ebp + nb030_p_tabscale]
	mov ecx, [ecx]
	mov edx, [edx]
	mov [esp + nb030_nri], ecx
	mov [esp + nb030_ntype], edx

	movd  mm3, [esi]
	punpckldq mm3,mm3
	movq  [esp + nb030_tsc], mm3	

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb030_nouter], eax
	mov [esp + nb030_ninner], eax

.nb030_threadloop:
        mov   esi, [ebp + nb030_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb030_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb030_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb030_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb030_n], eax
        mov [esp + nb030_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb030_outerstart
        jmp .nb030_end

.nb030_outerstart:	
	;# ebx contains number of outer iterations
	add ebx, [esp + nb030_nouter]
        mov [esp + nb030_nouter], ebx
	

.nb030_outer:
	mov   eax, [ebp + nb030_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb030_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb030_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + nb030_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx=ii 

	mov   edx, [ebp + nb030_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [esp + nb030_ntype]
	shl   edx, 1
	mov   [esp + nb030_ntia], edx

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb030_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb030_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + nb030_ix], mm0	
	movd  [esp + nb030_iz], mm1	
				
	;# clear total potential and i forces 
	pxor  mm7,mm7
	movq  [esp + nb030_Vvdwtot], mm7
	movq  [esp + nb030_fix],    mm7
	movd  [esp + nb030_fiz],    mm7

	mov   eax, [ebp + nb030_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb030_pos]
	mov   edi, [ebp + nb030_faction]	
	mov   eax, [ebp + nb030_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb030_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb030_ninner]
	mov   [esp + nb030_ninner], ecx
	mov   [esp + nb030_innerk], edx    ;# number of innerloop atoms 
	add   edx, 0
	jge   .nb030_unroll_loop
	jmp   .nb030_finish_inner
.nb030_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + nb030_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + nb030_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best  
	
	mov ecx, [ebp + nb030_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + nb030_vdwparam]		;# base of vdwparam  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + nb030_ntia]	     ;# tja = ntia + 2*type  
	add ecx, [esp + nb030_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + nb030_c6], mm5
	movq [esp + nb030_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + nb030_pos]

	movq  mm0, [esp + nb030_ix]
	movd  mm1, [esp + nb030_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + nb030_dx1], mm4	     ;# store dr 
	movd  [esp + nb030_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + nb030_dx2], mm6	     ;# store dr 
	movd  [esp + nb030_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r 
	;# do potential and fscal 
	pfmul mm1, [esp + nb030_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + nb030_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + nb030_VFtab]
	;# dispersion table 
	mov ecx, [esp + nb030_n1]
	shl ecx, 3
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + nb030_n1+4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + nb030_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + nb030_c6]
	pfmul mm7, mm4	;# fijD 
	pfmul mm5, mm4	;# Vvdw6            
	movq mm3, mm7	;# add to fscal 

	;# update Vvdwtot to release mm5! 
	pfadd mm5, [esp + nb030_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb030_Vvdwtot], mm5       ;# store the sum       

	;# repulsion table 
	mov ecx, [esp + nb030_n1]
	shl ecx, 3
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + nb030_n1+4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + nb030_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + nb030_c12]
	pfmul mm7, mm6	;# fijR 
	pfmul mm5, mm6	;# Vvdw12 
	pfadd mm3, mm7	;# total fscal fijD+ fijR 

	;# change sign of mm3 
    pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm1, [esp + nb030_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	prefetchw [esp + nb030_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + nb030_dx1]	;# fetch dr 
	movd mm3,  [esp + nb030_dz1]

	;# update Vvdwtot 
	pfadd mm5, [esp + nb030_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb030_Vvdwtot], mm5       ;# store the sum       

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + nb030_dx2] 	;# fetch dr 
	movd mm5,  [esp + nb030_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + nb030_fix]
	movd mm1,  [esp + nb030_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + nb030_fix], mm0
	movd [esp + nb030_fiz], mm1
	;# update j forces 

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
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb030_innerk],  2
	jl    .nb030_finish_inner
	jmp   .nb030_unroll_loop
.nb030_finish_inner:	
	and dword ptr [esp + nb030_innerk],  1
	jnz  .nb030_single_inner
	jmp  .nb030_updateouterdata		
.nb030_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + nb030_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + nb030_vdwparam]
	mov ecx, [ebp + nb030_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + nb030_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + nb030_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + nb030_c12], mm5

	mov   esi, [ebp + nb030_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + nb030_ix]
	movd  mm1, [esp + nb030_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + nb030_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + nb030_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + nb030_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + nb030_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + nb030_VFtab]
	mov ecx, [esp + nb030_n1]
	shl ecx, 3
	;# dispersion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + nb030_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + nb030_c6]
	pfmul mm7, mm4	;# fijD 
	pfmul mm5, mm4	;# Vvdw6            
	movq mm3, mm7	;# add to fscal 

	;# update Vvdwtot to release mm5! 
	pfadd mm5, [esp + nb030_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb030_Vvdwtot], mm5       ;# store the sum       

	;# repulsion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + nb030_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + nb030_c12]
	pfmul mm7, mm6	;# fijR 
	pfmul mm5, mm6	;# Vvdw12 
	pfadd mm3, mm7	;# total fscal fijC+ fijD+ fijR 

	;# change sign of mm3 
    pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + nb030_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	;# update Vvdwtot 
	pfadd mm5, [esp + nb030_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb030_Vvdwtot], mm5       ;# store the sum       

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + nb030_dx1]
	movd mm3,  [esp + nb030_dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + nb030_fix]
	movd mm1,  [esp + nb030_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + nb030_fix], mm0
	movd [esp + nb030_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.nb030_updateouterdata:	
	mov   ecx, [esp + nb030_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + nb030_fix]
	pfadd mm7, [esp + nb030_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + nb030_fshift]    ;# increment fshift force 
	mov   edx, [esp + nb030_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + nb030_fix]
	pfadd mm7, [esp + nb030_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# get n from stack
	mov esi, [esp + nb030_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb030_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb030_Vvdwtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb030_Vvdw]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment Vvdw[gid] 

       ;# finish if last 
        mov ecx, [esp + nb030_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb030_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb030_n], esi
        jmp .nb030_outer
.nb030_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb030_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb030_end
        ;# non-zero, do one more workunit
        jmp   .nb030_threadloop
.nb030_end:
	femms

	mov eax, [esp + nb030_nouter] 	
	mov ebx, [esp + nb030_ninner]
	mov ecx, [ebp + nb030_outeriter]
	mov edx, [ebp + nb030_inneriter]
	mov [ecx], eax
	mov [edx], ebx
	
	add esp, 140
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


			

	

.globl nb_kernel030nf_ia32_3dnow
.globl _nb_kernel030nf_ia32_3dnow
nb_kernel030nf_ia32_3dnow:	
_nb_kernel030nf_ia32_3dnow:	
.equiv		nb030nf_p_nri,		8
.equiv		nb030nf_iinr,		12
.equiv		nb030nf_jindex,		16
.equiv		nb030nf_jjnr,		20
.equiv		nb030nf_shift,		24
.equiv		nb030nf_shiftvec,	28
.equiv		nb030nf_fshift,		32
.equiv		nb030nf_gid,		36
.equiv		nb030nf_pos,		40		
.equiv		nb030nf_faction,	44
.equiv		nb030nf_charge,		48
.equiv		nb030nf_p_facel,	52
.equiv		nb030nf_p_krf,		56	
.equiv		nb030nf_p_crf,		60	
.equiv		nb030nf_Vc,		64	
.equiv		nb030nf_type,		68
.equiv		nb030nf_p_ntype,	72
.equiv		nb030nf_vdwparam,	76	
.equiv		nb030nf_Vvdw,		80	
.equiv		nb030nf_p_tabscale,	84	
.equiv		nb030nf_VFtab,		88
.equiv		nb030nf_invsqrta,	92	
.equiv		nb030nf_dvda,		96
.equiv          nb030nf_p_gbtabscale,   100
.equiv          nb030nf_GBtab,          104
.equiv          nb030nf_p_nthreads,     108
.equiv          nb030nf_count,          112
.equiv          nb030nf_mtx,            116
.equiv          nb030nf_outeriter,      120
.equiv          nb030nf_inneriter,      124
.equiv          nb030nf_work,           128
	;# stack offsets for local variables 
.equiv		nb030nf_is3,		0
.equiv		nb030nf_ii3,		4
.equiv		nb030nf_ix,		8
.equiv		nb030nf_iy,		12
.equiv		nb030nf_iz,		16
.equiv		nb030nf_Vvdwtot,	20 
.equiv		nb030nf_c6,		28 
.equiv		nb030nf_c12,		36
.equiv		nb030nf_n1,		44 
.equiv		nb030nf_tsc,		52 
.equiv		nb030nf_ntia,		60
.equiv		nb030nf_innerjjnr,	64
.equiv		nb030nf_innerk,		68
.equiv          nb030nf_n,              72 ;# idx for outer loop
.equiv          nb030nf_nn1,            76 ;# number of outer iterations
.equiv		nb030nf_nri,		80
.equiv		nb030nf_ntype,		84
.equiv		nb030nf_nouter,		88
.equiv		nb030nf_ninner,		92

    	push ebp
    	mov ebp,esp
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 96		;# local stack space 
	femms
	;# move data to local stack  
	mov ecx, [ebp + nb030nf_p_nri]
	mov edx, [ebp + nb030nf_p_ntype]
	mov esi, [ebp + nb030nf_p_tabscale]
	mov ecx, [ecx]
	mov edx, [edx]
	mov [esp + nb030nf_nri], ecx
	mov [esp + nb030nf_ntype], edx

	movd  mm3, [esi]
	punpckldq mm3,mm3
	movq  [esp + nb030nf_tsc], mm3	

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb030nf_nouter], eax
	mov [esp + nb030nf_ninner], eax

.nb030nf_threadloop:
        mov   esi, [ebp + nb030nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb030nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb030nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb030nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb030nf_n], eax
        mov [esp + nb030nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb030nf_outerstart
        jmp .nb030nf_end

.nb030nf_outerstart:	
	;# ebx contains number of outer iterations
	add ebx, [esp + nb030nf_nouter]
        mov [esp + nb030nf_nouter], ebx
	
.nb030nf_outer:
	mov   eax, [ebp + nb030nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb030nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb030nf_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + nb030nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx=ii 

	mov   edx, [ebp + nb030nf_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [esp + nb030nf_ntype]
	shl   edx, 1
	mov   [esp + nb030nf_ntia], edx

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb030nf_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb030nf_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + nb030nf_ix], mm0	
	movd  [esp + nb030nf_iz], mm1	
				
	;# clear total potential 
	pxor  mm7,mm7
	movq  [esp + nb030nf_Vvdwtot], mm7

	mov   eax, [ebp + nb030nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb030nf_pos]
	mov   eax, [ebp + nb030nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb030nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb030nf_ninner]
	mov   [esp + nb030nf_ninner], ecx
	mov   [esp + nb030nf_innerk], edx    ;# number of innerloop atoms 
	add   edx, 0
	jge   .nb030nf_unroll_loop
	jmp   .nb030nf_finish_inner
.nb030nf_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + nb030nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + nb030nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best  
	
	mov ecx, [ebp + nb030nf_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + nb030nf_vdwparam]		;# base of vdwparam  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + nb030nf_ntia]	     ;# tja = ntia + 2*type  
	add ecx, [esp + nb030nf_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + nb030nf_c6], mm5
	movq [esp + nb030nf_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + nb030nf_pos]

	movq  mm0, [esp + nb030nf_ix]
	movd  mm1, [esp + nb030nf_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r 
	;# do potential and fscal 
	pfmul mm1, [esp + nb030nf_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + nb030nf_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + nb030nf_VFtab]
	;# dispersion table 
	mov ecx, [esp + nb030nf_n1]
	shl ecx, 3
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + nb030nf_n1+4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + nb030nf_c6]
	pfmul mm5, mm4	;# Vvdw6        
	;# update Vvdwtot to release mm5! 
	pfadd mm5, [esp + nb030nf_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb030nf_Vvdwtot], mm5       ;# store the sum       

	;# repulsion table 
	mov ecx, [esp + nb030nf_n1]
	shl ecx, 3
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + nb030nf_n1+4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + nb030nf_c12]
	pfmul mm5, mm6	;# Vvdw12 
	;# update Vvdwtot 
	pfadd mm5, [esp + nb030nf_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb030nf_Vvdwtot], mm5       ;# store the sum       
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb030nf_innerk],  2
	jl    .nb030nf_finish_inner
	jmp   .nb030nf_unroll_loop
.nb030nf_finish_inner:	
	and dword ptr [esp + nb030nf_innerk],  1
	jnz  .nb030nf_single_inner
	jmp  .nb030nf_updateouterdata		
.nb030nf_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + nb030nf_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + nb030nf_vdwparam]
	mov ecx, [ebp + nb030nf_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + nb030nf_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + nb030nf_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + nb030nf_c12], mm5

	mov   esi, [ebp + nb030nf_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + nb030nf_ix]
	movd  mm1, [esp + nb030nf_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + nb030nf_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + nb030nf_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + nb030nf_VFtab]
	mov ecx, [esp + nb030nf_n1]
	shl ecx, 3
	;# dispersion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + nb030nf_c6]
	pfmul mm5, mm4	;# Vvdw6            
	;# update Vvdwtot to release mm5! 
	pfadd mm5, [esp + nb030nf_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb030nf_Vvdwtot], mm5       ;# store the sum       

	;# repulsion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + nb030nf_c12]
	pfmul mm5, mm6	;# Vvdw12 
	;# update Vvdwtot 
	pfadd mm5, [esp + nb030nf_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb030nf_Vvdwtot], mm5       ;# store the sum       

.nb030nf_updateouterdata:	
	;# get n from stack
	mov esi, [esp + nb030nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb030nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb030nf_Vvdwtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb030nf_Vvdw]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment Vvdw[gid] 

       ;# finish if last 
        mov ecx, [esp + nb030nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb030nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb030nf_n], esi
        jmp .nb030nf_outer
.nb030nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb030nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb030nf_end
        ;# non-zero, do one more workunit
        jmp   .nb030nf_threadloop
.nb030nf_end:
	femms

	mov eax, [esp + nb030nf_nouter] 	
	mov ebx, [esp + nb030nf_ninner]
	mov ecx, [ebp + nb030nf_outeriter]
	mov edx, [ebp + nb030nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx
	
	add esp, 96
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret
