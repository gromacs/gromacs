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





	
	

.globl nb_kernel310_ia32_3dnow
.globl _nb_kernel310_ia32_3dnow
nb_kernel310_ia32_3dnow:	
_nb_kernel310_ia32_3dnow:	
.equiv		nb310_p_nri,		8
.equiv		nb310_iinr,		12
.equiv		nb310_jindex,		16
.equiv		nb310_jjnr,		20
.equiv		nb310_shift,		24
.equiv		nb310_shiftvec,		28
.equiv		nb310_fshift,		32
.equiv		nb310_gid,		36
.equiv		nb310_pos,		40		
.equiv		nb310_faction,		44
.equiv		nb310_charge,		48
.equiv		nb310_p_facel,		52
.equiv		nb310_p_krf,		56	
.equiv		nb310_p_crf,		60	
.equiv		nb310_Vc,		64	
.equiv		nb310_type,		68
.equiv		nb310_p_ntype,		72
.equiv		nb310_vdwparam,		76	
.equiv		nb310_Vvdw,		80	
.equiv		nb310_p_tabscale,	84	
.equiv		nb310_VFtab,		88
.equiv		nb310_invsqrta,		92	
.equiv		nb310_dvda,		96
.equiv          nb310_p_gbtabscale,     100
.equiv          nb310_GBtab,            104
.equiv          nb310_p_nthreads,       108
.equiv          nb310_count,            112
.equiv          nb310_mtx,              116
.equiv          nb310_outeriter,        120
.equiv          nb310_inneriter,        124
.equiv          nb310_work,             128
	;# stack offsets for local variables 
.equiv		nb310_is3,		0 
.equiv		nb310_ii3,		4
.equiv		nb310_ix,		8
.equiv		nb310_iy,		12
.equiv		nb310_iz,		16
.equiv		nb310_iq,		20 
.equiv		nb310_vctot,		28 
.equiv		nb310_Vvdwtot,		36 
.equiv		nb310_c6,		44 
.equiv		nb310_c12,		52 
.equiv		nb310_six,		60 
.equiv		nb310_twelve,		68 
.equiv		nb310_two,		76 
.equiv		nb310_n1,		84 
.equiv		nb310_tsc,		92 
.equiv		nb310_ntia,		100
.equiv		nb310_innerjjnr,	104
.equiv		nb310_innerk,		108	
.equiv		nb310_fix,		112
.equiv		nb310_fiy,		116
.equiv		nb310_fiz,		120
.equiv		nb310_dx1,		124
.equiv		nb310_dy1,		128
.equiv		nb310_dz1,		132
.equiv		nb310_dx2,		136
.equiv		nb310_dy2,		140
.equiv		nb310_dz2,		144
.equiv          nb310_n,                148 ;# idx for outer loop
.equiv          nb310_nn1,              152 ;# number of outer iterations
.equiv          nb310_nri,              156
.equiv          nb310_facel,            160
.equiv          nb310_ntype,            164
.equiv          nb310_nouter,           168
.equiv          nb310_ninner,           172
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 176		;# local stack space 
	femms
	;# move data to local stack
	mov ecx, [ebp + nb310_p_nri]
	mov edx, [ebp + nb310_p_ntype]
	mov esi, [ebp + nb310_p_facel]
	mov edi, [ebp + nb310_p_tabscale]
	mov ecx, [ecx]
	mov edx, [edx]
	mov esi, [esi]
	mov [esp + nb310_nri], ecx
	mov [esp + nb310_ntype], edx
	mov [esp + nb310_facel], esi

	movd  mm3, [edi]
	punpckldq mm3,mm3
	movq  [esp + nb310_tsc], mm3	
	mov eax, 0x40000000    ;# 2.0
	mov [esp + nb310_two], eax
	mov [esp + nb310_two + 4], eax
	mov ebx, 0x40c00000   ;# 6.0
	mov [esp + nb310_six], ebx
	mov [esp + nb310_six + 4], ebx
	mov ecx, 0x41400000    ;# 12.0
	mov [esp + nb310_twelve], ecx
	mov [esp + nb310_twelve + 4], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb310_nouter], eax
	mov [esp + nb310_ninner], eax
	
.nb310_threadloop:
        mov   esi, [ebp + nb310_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb310_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [esi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb310_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb310_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb310_n], eax
        mov [esp + nb310_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb310_outerstart
        jmp .nb310_end

.nb310_outerstart:	
	;# ebx contains number of outer iterations
	add ebx, [esp + nb310_nouter]
        mov [esp + nb310_nouter], ebx
	
.nb310_outer:
	mov   eax, [ebp + nb310_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb310_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb310_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + nb310_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx=ii 

	mov   edx, [ebp + nb310_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [esp + nb310_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + nb310_iq], mm2	    ;# iq =facel*charge[ii] 

	mov   edx, [ebp + nb310_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [esp + nb310_ntype]
	shl   edx, 1
	mov   [esp + nb310_ntia], edx

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb310_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb310_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + nb310_ix], mm0	
	movd  [esp + nb310_iz], mm1	
				
	;# clear total potential and i forces 
	pxor  mm7,mm7
	movq  [esp + nb310_vctot],  mm7
	movq  [esp + nb310_Vvdwtot], mm7
	movq  [esp + nb310_fix],    mm7
	movd  [esp + nb310_fiz],    mm7

	mov   eax, [ebp + nb310_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb310_pos]
	mov   edi, [ebp + nb310_faction]	
	mov   eax, [ebp + nb310_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb310_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb310_ninner]
	mov   [esp + nb310_ninner], ecx
	mov   [esp + nb310_innerk], edx    ;# number of innerloop atoms 
	add   edx, 0
	jge   .nb310_unroll_loop
	jmp   .nb310_finish_inner
.nb310_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + nb310_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + nb310_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + nb310_charge]    ;# base of charge[] 
	movq mm5, [esp + nb310_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    	punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	mov ecx, [ebp + nb310_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + nb310_vdwparam]		;# base of vdwparam  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + nb310_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + nb310_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + nb310_c6], mm5
	movq [esp + nb310_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + nb310_pos]

	movq  mm0, [esp + nb310_ix]
	movd  mm1, [esp + nb310_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + nb310_dx1], mm4	     ;# store dr 
	movd  [esp + nb310_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + nb310_dx2], mm6	     ;# store dr 
	movd  [esp + nb310_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    	pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    	pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    	movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0,mm0
    	pfrsqit1 mm0,mm4				
    	pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + nb310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + nb310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + nb310_VFtab]
	mov ecx, [esp + nb310_n1]
	shl ecx, 2
	;# coulomb table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + nb310_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + nb310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  

	movq mm1, mm0
	pfmul mm1,mm1 	;# mm1=invsq 
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	;# mm2=rinvsix 
	movq  mm1,mm2
	pfmul mm1,mm1	;# mm1=rinvtwelve 
	
	pfmul mm3, [esp + nb310_tsc]
	
	pfmul mm1, [esp + nb310_c12]

	pfmul mm2, [esp + nb310_c6]

	movq mm4, mm1
	pfsub mm4, mm2	;# mm4 = Vvdw12-Vvdw6 

	pfmul mm2, [esp + nb310_six]
	pfmul mm1, [esp + nb310_twelve]

	pfsub mm1, mm2
	pfmul mm1, mm0	;# mm1=	(12*Vvdw12-6*Vvdw6)*rinv11 

	pfsub mm1, mm3

	;# update vctot 
	pfadd mm5, [esp + nb310_vctot]      ;# add the earlier value 
	movq [esp + nb310_vctot], mm5       ;# store the sum       

 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	prefetchw [esp + nb310_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + nb310_dx1]	;# fetch dr 
	movd mm3,  [esp + nb310_dz1]

	;# update Vvdwtot 
	pfadd mm4, [esp + nb310_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb310_Vvdwtot], mm4       ;# store the sum       

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + nb310_dx2] 	;# fetch dr 
	movd mm5,  [esp + nb310_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + nb310_fix]
	movd mm1,  [esp + nb310_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + nb310_fix], mm0
	movd [esp + nb310_fiz], mm1
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
	sub dword ptr [esp + nb310_innerk],  2
	jl    .nb310_finish_inner
	jmp   .nb310_unroll_loop
.nb310_finish_inner:	
	and dword ptr [esp + nb310_innerk],  1
	jnz  .nb310_single_inner
	jmp  .nb310_updateouterdata		
.nb310_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + nb310_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + nb310_charge]
	movd mm5, [esp + nb310_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov esi, [ebp + nb310_vdwparam]
	mov ecx, [ebp + nb310_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + nb310_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + nb310_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + nb310_c12], mm5


	mov   esi, [ebp + nb310_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + nb310_ix]
	movd  mm1, [esp + nb310_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + nb310_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + nb310_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm4=rsq 
	
    	pfrsqrt mm0,mm4
    	movq mm2,mm0
    	pfmul mm0,mm0
    	pfrsqit1 mm0,mm4				
    	pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# calculate potentials and scalar force 
	pfmul mm1, [esp + nb310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + nb310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + nb310_VFtab]
	mov ecx, [esp + nb310_n1]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + nb310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  
	
	;# at this point mm5 contains vcoul and mm3 fijC 

	movq mm1, mm0
	pfmul mm1,mm1 	;# mm1=invsq 
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	;# mm2=rinvsix 
	movq  mm1,mm2
	pfmul mm1,mm1	;# mm1=rinvtwelve 
	
	pfmul mm3, [esp + nb310_tsc]
	
	pfmul mm1, [esp + nb310_c12]

	pfmul mm2, [esp + nb310_c6]

	movq mm4, mm1
	pfsub mm4, mm2	;# mm4 = Vvdw12-Vvdw6 

	pfmul mm2, [esp + nb310_six]
	pfmul mm1, [esp + nb310_twelve]

	pfsub mm1, mm2
	pfmul mm1, mm0	;# mm1=	(12*Vvdw12-6*Vvdw6)*rinv11 

	pfsub mm1, mm3

	;# update vctot 
	pfadd mm5, [esp + nb310_vctot]      ;# add the earlier value 
	movq [esp + nb310_vctot], mm5       ;# store the sum       

 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + nb310_dx1]
	movd mm3,  [esp + nb310_dz1]

	;# update Vvdwtot 
	pfadd mm4, [esp + nb310_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb310_Vvdwtot], mm4       ;# store the sum       

	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + nb310_fix]
	movd mm1,  [esp + nb310_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + nb310_fix], mm0
	movd [esp + nb310_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.nb310_updateouterdata:	
	mov   ecx, [esp + nb310_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + nb310_fix]
	pfadd mm7, [esp + nb310_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + nb310_fshift]    ;# increment fshift force 
	mov   edx, [esp + nb310_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + nb310_fix] 
	pfadd mm7, [esp + nb310_fiz] 
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# get n from stack
	mov esi, [esp + nb310_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb310_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb310_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb310_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
 
	movq  mm7, [esp + nb310_Vvdwtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb310_Vvdw] 
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment Vvdw[gid] 

       	;# finish if last 
        mov ecx, [esp + nb310_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb310_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb310_n], esi
        jmp .nb310_outer
.nb310_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb310_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb310_end
        ;# non-zero, do one more workunit
        jmp   .nb310_threadloop
.nb310_end:
	femms

	mov eax, [esp + nb310_nouter] 	
	mov ebx, [esp + nb310_ninner]
	mov ecx, [ebp + nb310_outeriter]
	mov edx, [ebp + nb310_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	add esp, 176
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret







.globl nb_kernel310nf_ia32_3dnow
.globl _nb_kernel310nf_ia32_3dnow
nb_kernel310nf_ia32_3dnow:	
_nb_kernel310nf_ia32_3dnow:	
.equiv		nb310nf_p_nri,		8
.equiv		nb310nf_iinr,		12
.equiv		nb310nf_jindex,		16
.equiv		nb310nf_jjnr,		20
.equiv		nb310nf_shift,		24
.equiv		nb310nf_shiftvec,	28
.equiv		nb310nf_fshift,		32
.equiv		nb310nf_gid,		36
.equiv		nb310nf_pos,		40		
.equiv		nb310nf_faction,	44
.equiv		nb310nf_charge,		48
.equiv		nb310nf_p_facel,	52
.equiv		nb310nf_p_krf,		56	
.equiv		nb310nf_p_crf,		60	
.equiv		nb310nf_Vc,		64	
.equiv		nb310nf_type,		68
.equiv		nb310nf_p_ntype,	72
.equiv		nb310nf_vdwparam,	76	
.equiv		nb310nf_Vvdw,		80	
.equiv		nb310nf_p_tabscale,	84	
.equiv		nb310nf_VFtab,		88
.equiv		nb310nf_invsqrta,	92	
.equiv		nb310nf_dvda,		96
.equiv          nb310nf_p_gbtabscale,   100
.equiv          nb310nf_GBtab,          104
.equiv          nb310nf_p_nthreads,     108
.equiv          nb310nf_count,          112
.equiv          nb310nf_mtx,            116
.equiv          nb310nf_outeriter,      120
.equiv          nb310nf_inneriter,      124
.equiv          nb310nf_work,           128
	;# stack offsets for local variables 
.equiv		nb310nf_is3,		0 
.equiv		nb310nf_ii3,		4
.equiv		nb310nf_ix,		8
.equiv		nb310nf_iy,		12
.equiv		nb310nf_iz,		16
.equiv		nb310nf_iq,		20 
.equiv		nb310nf_vctot,		28 
.equiv		nb310nf_Vvdwtot,	36 
.equiv		nb310nf_c6,		44 
.equiv		nb310nf_c12,		52
.equiv		nb310nf_n1,		60 
.equiv		nb310nf_tsc,		68 
.equiv		nb310nf_ntia,		76
.equiv		nb310nf_innerjjnr,	80
.equiv		nb310nf_innerk,		84				
.equiv          nb310nf_n,              88 ;# idx for outer loop
.equiv          nb310nf_nn1,            92 ;# number of outer iterations
.equiv          nb310nf_nri,            96
.equiv          nb310nf_facel,          100
.equiv          nb310nf_ntype,          104
.equiv          nb310nf_nouter,         108
.equiv          nb310nf_ninner,         112
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 116		;# local stack space 
	femms
	;# move data to local stack  
	mov ecx, [ebp + nb310nf_p_nri]
	mov edx, [ebp + nb310nf_p_ntype]
	mov esi, [ebp + nb310nf_p_facel]
	mov edi, [ebp + nb310nf_p_tabscale]
	mov ecx, [ecx]
	mov edx, [edx]
	mov esi, [esi]
	mov [esp + nb310nf_nri], ecx
	mov [esp + nb310nf_ntype], edx
	mov [esp + nb310nf_facel], esi

	movd  mm3, [edi]
	punpckldq mm3,mm3
	movq  [esp + nb310nf_tsc], mm3	

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb310nf_nouter], eax
	mov [esp + nb310nf_ninner], eax
	
.nb310nf_threadloop:
        mov   esi, [ebp + nb310nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb310nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [esi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb310nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb310nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb310nf_n], eax
        mov [esp + nb310nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb310nf_outerstart
        jmp .nb310nf_end

.nb310nf_outerstart:	
	;# ebx contains number of outer iterations
	add ebx, [esp + nb310nf_nouter]
        mov [esp + nb310nf_nouter], ebx
	
.nb310nf_outer:
	mov   eax, [ebp + nb310nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 

	mov   eax, [ebp + nb310nf_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [esp + nb310nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx=ii 

	mov   edx, [ebp + nb310nf_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [esp + nb310nf_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + nb310nf_iq], mm2	    ;# iq =facel*charge[ii] 

	mov   edx, [ebp + nb310nf_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [esp + nb310nf_ntype]
	shl   edx, 1
	mov   [esp + nb310nf_ntia], edx

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb310nf_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	pfadd mm1, mm3
	movq  [esp + nb310nf_ix], mm0	
	movd  [esp + nb310nf_iz], mm1	
				
	;# clear total potential and i forces 
	pxor  mm7,mm7
	movq  [esp + nb310nf_vctot],  mm7
	movq  [esp + nb310nf_Vvdwtot], mm7

	mov   eax, [ebp + nb310nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb310nf_pos]
	mov   eax, [ebp + nb310nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb310nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb310nf_ninner]
	mov   [esp + nb310nf_ninner], ecx
	mov   [esp + nb310nf_innerk], edx    ;# number of innerloop atoms 
	add   edx, 0
	jmp .nb310nf_finish_inner
	jge   .nb310nf_unroll_loop
	jmp   .nb310nf_finish_inner
.nb310nf_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + nb310nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + nb310nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + nb310nf_charge]    ;# base of charge[] 
	movq mm5, [esp + nb310nf_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    	punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	mov ecx, [ebp + nb310nf_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + nb310nf_vdwparam]		;# base of vdwparam  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + nb310nf_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + nb310nf_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + nb310nf_c6], mm5
	movq [esp + nb310nf_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + nb310nf_pos]

	movq  mm0, [esp + nb310nf_ix]
	movd  mm1, [esp + nb310nf_iz]	 	
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
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    	movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0,mm0
    	pfrsqit1 mm0,mm4				
    	pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + nb310nf_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + nb310nf_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + nb310nf_VFtab]
	mov ecx, [esp + nb310nf_n1]
	shl ecx, 2
	;# coulomb table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + nb310nf_n1 + 4]
	shl ecx, 2
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

	pfmul mm5, mm3	;# vcoul=qq*VV 

	movq mm1, mm0
	pfmul mm1,mm1 	;# mm1=invsq 
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	;# mm2=rinvsix 
	movq  mm1,mm2
	pfmul mm1,mm1	;# mm1=rinvtwelve 
	
	pfmul mm3, [esp + nb310nf_tsc]
	
	pfmul mm1, [esp + nb310nf_c12]

	pfmul mm2, [esp + nb310nf_c6]

	movq mm4, mm1
	pfsub mm4, mm2	;# mm4 = Vvdw12-Vvdw6 
	;# update vctot 
	pfadd mm5, [esp + nb310nf_vctot]      ;# add the earlier value 
	movq [esp + nb310nf_vctot], mm5       ;# store the sum       
	;# update Vvdwtot 
	pfadd mm4, [esp + nb310nf_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb310nf_Vvdwtot], mm4       ;# store the sum       
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb310nf_innerk],  2
	jl    .nb310nf_finish_inner
	jmp   .nb310nf_unroll_loop
.nb310nf_finish_inner:	
	and dword ptr [esp + nb310nf_innerk],  1
	jnz  .nb310nf_single_inner
	jmp  .nb310nf_updateouterdata		
.nb310nf_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + nb310nf_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + nb310nf_charge]
	movd mm5, [esp + nb310nf_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov esi, [ebp + nb310nf_vdwparam]
	mov ecx, [ebp + nb310nf_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + nb310nf_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + nb310nf_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + nb310nf_c12], mm5


	mov   esi, [ebp + nb310nf_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + nb310nf_ix]
	movd  mm1, [esp + nb310nf_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm4=rsq 
	
    	pfrsqrt mm0,mm4
    	movq mm2,mm0
    	pfmul mm0,mm0
    	pfrsqit1 mm0,mm4				
    	pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# calculate potentials and scalar force 
	pfmul mm1, [esp + nb310nf_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + nb310nf_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + nb310nf_VFtab]
	mov ecx, [esp + nb310nf_n1]
	shl ecx, 2
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

	pfmul mm5, mm3	;# vcoul=qq*VV 
	;# at this point mm5 contains vcoul 

	movq mm1, mm0
	pfmul mm1,mm1 	;# mm1=invsq 
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	;# mm2=rinvsix 
	movq  mm1,mm2
	pfmul mm1,mm1	;# mm1=rinvtwelve 
	
	pfmul mm3, [esp + nb310nf_tsc]
	
	pfmul mm1, [esp + nb310nf_c12]

	pfmul mm2, [esp + nb310nf_c6]

	movq mm4, mm1
	pfsub mm4, mm2	;# mm4 = Vvdw12-Vvdw6 
	;# update vctot 
	pfadd mm5, [esp + nb310nf_vctot]      ;# add the earlier value 
	movq [esp + nb310nf_vctot], mm5       ;# store the sum       
	;# update Vvdwtot 
	pfadd mm4, [esp + nb310nf_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb310nf_Vvdwtot], mm4       ;# store the sum       

.nb310nf_updateouterdata:	
	;# get n from stack
	mov esi, [esp + nb310nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb310nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb310nf_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb310nf_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
 
	movq  mm7, [esp + nb310nf_Vvdwtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb310nf_Vvdw] 
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment Vvdw[gid] 

       	;# finish if last 
        mov ecx, [esp + nb310nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb310nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb310nf_n], esi
        jmp .nb310nf_outer
.nb310nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb310nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb310nf_end
        ;# non-zero, do one more workunit
        jmp   .nb310nf_threadloop
.nb310nf_end:
	femms

	mov eax, [esp + nb310nf_nouter] 	
	mov ebx, [esp + nb310nf_ninner]
	mov ecx, [ebp + nb310nf_outeriter]
	mov edx, [ebp + nb310nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	add esp, 116
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret

