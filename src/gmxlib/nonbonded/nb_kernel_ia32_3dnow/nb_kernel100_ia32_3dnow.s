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




	

.globl nb_kernel100_ia32_3dnow
.globl _nb_kernel100_ia32_3dnow
nb_kernel100_ia32_3dnow:	
_nb_kernel100_ia32_3dnow:	
.equiv		nb100_p_nri,		8
.equiv		nb100_iinr,		12
.equiv		nb100_jindex,		16
.equiv		nb100_jjnr,		20
.equiv		nb100_shift,		24
.equiv		nb100_shiftvec,		28
.equiv		nb100_fshift,		32
.equiv		nb100_gid,		36
.equiv		nb100_pos,		40		
.equiv		nb100_faction,		44
.equiv		nb100_charge,		48
.equiv		nb100_p_facel,		52
.equiv		nb100_p_krf,		56	
.equiv		nb100_p_crf,		60	
.equiv		nb100_Vc,		64	
.equiv		nb100_type,		68
.equiv		nb100_p_ntype,		72
.equiv		nb100_vdwparam,		76	
.equiv		nb100_Vvdw,		80	
.equiv		nb100_p_tabscale,	84	
.equiv		nb100_VFtab,		88
.equiv		nb100_invsqrta,		92	
.equiv		nb100_dvda,		96
.equiv          nb100_p_gbtabscale,     100
.equiv          nb100_GBtab,            104
.equiv          nb100_p_nthreads,       108
.equiv          nb100_count,            112
.equiv          nb100_mtx,              116
.equiv          nb100_outeriter,        120
.equiv          nb100_inneriter,        124
.equiv          nb100_work,             128
	;# stack offsets for local variables 
.equiv		nb100_is3,		0
.equiv		nb100_ii3,		4
.equiv		nb100_ix,		8
.equiv		nb100_iy,		12
.equiv		nb100_iz,		16
.equiv		nb100_iq,		20		
.equiv		nb100_vctot,		28 
.equiv		nb100_innerjjnr,	36
.equiv		nb100_innerk,		40		
.equiv		nb100_fix,		44
.equiv		nb100_fiy,		48
.equiv		nb100_fiz,		52
.equiv		nb100_dx1,		56
.equiv		nb100_dy1,		60
.equiv		nb100_dz1,		64
.equiv		nb100_dx2,		68
.equiv		nb100_dy2,		72
.equiv		nb100_dz2,		76	
.equiv          nb100_n,                80 ;# idx for outer loop
.equiv          nb100_nn1,              84 ;# number of outer iterations
.equiv          nb100_nri,              88
.equiv          nb100_facel,            92
.equiv          nb100_nouter,           96
.equiv          nb100_ninner,          100
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 104		;# local stack space 
	femms

	mov ecx, [ebp + nb100_p_nri]
	mov esi, [ebp + nb100_p_facel]
	mov ecx, [ecx]
	mov esi, [esi]
	mov [esp + nb100_nri], ecx
	mov [esp + nb100_facel], esi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb100_nouter], eax
	mov [esp + nb100_ninner], eax

.nb100_threadloop:
        mov   esi, [ebp + nb100_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb100_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [esi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb100_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb100_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb100_n], eax
        mov [esp + nb100_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi

        jg  .nb100_outerstart
        jmp .nb100_end

.nb100_outerstart:	
	;# ebx contains number of outer iterations
	add ebx, [esp + nb100_nouter]
        mov [esp + nb100_nouter], ebx
	
.nb100_outer:
	mov   eax, [ebp + nb100_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb100_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb100_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + nb100_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx=ii 

	mov   edx, [ebp + nb100_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [esp + nb100_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + nb100_iq], mm2	    ;# iq =facel*charge[ii] 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb100_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb100_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + nb100_ix], mm0	
	movd  [esp + nb100_iz], mm1	
				
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + nb100_vctot], mm7
	movq  [esp + nb100_fix],   mm7
	movd  [esp + nb100_fiz],   mm7

	mov   eax, [ebp + nb100_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb100_pos]
	mov   edi, [ebp + nb100_faction]	
	mov   eax, [ebp + nb100_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb100_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb100_ninner]
	mov   [esp + nb100_ninner], ecx
	mov   [esp + nb100_innerk], edx    ;# number of innerloop atoms 
	add   edx, 0
	jge   .nb100_unroll_loop
	jmp   .nb100_finish_inner
.nb100_unroll_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + nb100_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + nb100_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + nb100_charge]    ;# base of charge[] 
	movq mm5, [esp + nb100_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
	movd mm7, [ecx + ebx*4]  	 ;# charge[jnr2] 
	punpckldq mm3,mm7	     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	movq  mm0, [esp + nb100_ix]
	movd  mm1, [esp + nb100_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + nb100_dx1], mm4	     ;# store dr 
	movd  [esp + nb100_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + nb100_dx2], mm6	     ;# store dr 
	movd  [esp + nb100_dz2], mm7
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
	movq mm1,mm0
	pfmul mm0,mm0
	;# mm0 now contains invsq, and mm1 invsqrt
	 ;# do potential and fscal
	 
	prefetchw [esp + nb100_dx1]	;# prefetch i forces to cache 
	
	pfmul mm3,mm1		;# mm3 has both vcoul 
	pfmul mm0,mm3		;# mm0 has both fscal 

	;# update vctot 

	pfadd mm3, [esp + nb100_vctot]      ;# add the earlier value  
	movq [esp + nb100_vctot], mm3       ;# store the sum 
	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1
	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + nb100_dx1]	;# fetch dr 
	movd mm3,  [esp + nb100_dz1]
	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + nb100_dx2] 	;# fetch dr 
	movd mm5,  [esp + nb100_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + nb100_fix]
	movd mm1,  [esp + nb100_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + nb100_fix], mm0
	movd [esp + nb100_fiz], mm1
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
	sub dword ptr [esp + nb100_innerk],  2
	jl    .nb100_finish_inner
	jmp   .nb100_unroll_loop
.nb100_finish_inner:	
	and dword ptr [esp + nb100_innerk],  1
	jnz  .nb100_single_inner
	jmp  .nb100_updateouterdata		
.nb100_single_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + nb100_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + nb100_charge]
	pxor mm6, mm6
	movd mm6, [esp + nb100_iq]
	movd mm7, [ecx + eax*4]
	pfmul mm6, mm7	  	;# mm6=qq 
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esp + nb100_ix]
	movd  mm1, [esp + nb100_iz]
	movq  mm2, [esi + eax*4]
	movd  mm3, [esi + eax*4 + 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq  [esp + nb100_dx1], mm0
	pfmul mm0,mm0
	movd  [esp + nb100_dz1], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfacc mm0, mm1		;# mm0=rsq 
	
	pfrsqrt mm1,mm0
	movq mm2,mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	movq  mm4, mm1
	pfmul mm4, mm4		;# mm4=invsq 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	pfmul mm4, mm6		;# mm4=fscalar  
	
	;# update vctot 
	pfadd mm6, [esp + nb100_vctot]
	movq [esp + nb100_vctot], mm6
	;# spread fscalar to both positions 
	punpckldq mm4,mm4
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm0,  [esp + nb100_dx1]
	movd mm1,  [esp + nb100_dz1]
	pfmul mm0, mm4
	pfmul mm1, mm4
	;# update i particle force 
	movq mm2,  [esp + nb100_fix]
	movd mm3,  [esp + nb100_fiz]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + nb100_fix], mm2
	movd [esp + nb100_fiz], mm3
	;# update j particle force 
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax *4+ 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	;# done! 
.nb100_updateouterdata:	
	mov   ecx, [esp + nb100_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + nb100_fix]
	pfadd mm7, [esp + nb100_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + nb100_fshift]    ;# increment fshift force 
	mov   edx, [esp + nb100_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + nb100_fix]
	pfadd mm7, [esp + nb100_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# get n from stack
	mov esi, [esp + nb100_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb100_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb100_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb100_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
       	;# finish if last 
        mov ecx, [esp + nb100_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb100_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb100_n], esi
        jmp .nb100_outer
.nb100_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb100_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb100_end
        ;# non-zero, do one more workunit
        jmp   .nb100_threadloop
.nb100_end:
	femms

	mov eax, [esp + nb100_nouter] 	
	mov ebx, [esp + nb100_ninner]
	mov ecx, [ebp + nb100_outeriter]
	mov edx, [ebp + nb100_inneriter]
	mov [ecx], eax
	mov [edx], ebx
	
	add esp, 104
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret








	

.globl nb_kernel100nf_ia32_3dnow
.globl _nb_kernel100nf_ia32_3dnow
nb_kernel100nf_ia32_3dnow:	
_nb_kernel100nf_ia32_3dnow:	
.equiv		nb100nf_p_nri,		8
.equiv		nb100nf_iinr,		12
.equiv		nb100nf_jindex,		16
.equiv		nb100nf_jjnr,		20
.equiv		nb100nf_shift,		24
.equiv		nb100nf_shiftvec,	28
.equiv		nb100nf_fshift,		32
.equiv		nb100nf_gid,		36
.equiv		nb100nf_pos,		40		
.equiv		nb100nf_faction,	44
.equiv		nb100nf_charge,		48
.equiv		nb100nf_p_facel,	52
.equiv		nb100nf_p_krf,		56	
.equiv		nb100nf_p_crf,		60	
.equiv		nb100nf_Vc,		64	
.equiv		nb100nf_type,		68
.equiv		nb100nf_p_ntype,	72
.equiv		nb100nf_vdwparam,	76	
.equiv		nb100nf_Vvdw,		80	
.equiv		nb100nf_p_tabscale,	84	
.equiv		nb100nf_VFtab,		88
.equiv		nb100nf_invsqrta,	92	
.equiv		nb100nf_dvda,		96
.equiv          nb100nf_p_gbtabscale,   100
.equiv          nb100nf_GBtab,          104
.equiv          nb100nf_p_nthreads,     108
.equiv          nb100nf_count,          112
.equiv          nb100nf_mtx,            116
.equiv          nb100nf_outeriter,      120
.equiv          nb100nf_inneriter,      124
.equiv          nb100nf_work,           128
	;# stack offsets for local variables 
.equiv		nb100nf_is3,		0
.equiv		nb100nf_ii3,		4
.equiv		nb100nf_ix,		8
.equiv		nb100nf_iy,		12
.equiv		nb100nf_iz,		16
.equiv		nb100nf_iq,		20		
.equiv		nb100nf_vctot,		28 
.equiv		nb100nf_innerjjnr,	36
.equiv		nb100nf_innerk,		40	
.equiv          nb100nf_n,              44 ;# idx for outer loop
.equiv          nb100nf_nn1,            48 ;# number of outer iterations
.equiv          nb100nf_nri,            52
.equiv          nb100nf_facel,          56
.equiv          nb100nf_nouter,         60
.equiv          nb100nf_ninner,         64
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 68		;# local stack space 
	femms

	mov ecx, [ebp + nb100nf_p_nri]
	mov esi, [ebp + nb100nf_p_facel]
	mov ecx, [ecx]
	mov esi, [esi]
	mov [esp + nb100nf_nri], ecx
	mov [esp + nb100nf_facel], esi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb100nf_nouter], eax
	mov [esp + nb100nf_ninner], eax

.nb100nf_threadloop:
        mov   esi, [ebp + nb100nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb100nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [esi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb100nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb100nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb100nf_n], eax
        mov [esp + nb100nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi

        jg  .nb100nf_outerstart
        jmp .nb100nf_end

.nb100nf_outerstart:	
	;# ebx contains number of outer iterations
	add ebx, [esp + nb100nf_nouter]
        mov [esp + nb100nf_nouter], ebx
	
.nb100nf_outer:
	mov   eax, [ebp + nb100nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb100nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb100nf_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + nb100nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx=ii 

	mov   edx, [ebp + nb100nf_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [esp + nb100nf_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + nb100nf_iq], mm2	    ;# iq =facel*charge[ii] 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb100nf_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb100nf_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + nb100nf_ix], mm0	
	movd  [esp + nb100nf_iz], mm1	
				
	;# clear vctot
	pxor  mm7,mm7
	movq  [esp + nb100nf_vctot], mm7

	mov   eax, [ebp + nb100nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb100nf_pos]
	mov   edi, [ebp + nb100nf_faction]	
	mov   eax, [ebp + nb100nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb100nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb100nf_ninner]
	mov   [esp + nb100nf_ninner], ecx
	mov   [esp + nb100nf_innerk], edx    ;# number of innerloop atoms 
	add   edx, 0
	jge   .nb100nf_unroll_loop
	jmp   .nb100nf_finish_inner
.nb100nf_unroll_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + nb100nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + nb100nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + nb100nf_charge]    ;# base of charge[] 
	movq mm5, [esp + nb100nf_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
	movd mm7, [ecx + ebx*4]  	 ;# charge[jnr2] 
	punpckldq mm3,mm7	     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	movq  mm0, [esp + nb100nf_ix]
	movd  mm1, [esp + nb100nf_iz]	 	
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
	;# mm0=invsqrt
	 ;# do potential and fscal
	 
	pfmul mm3,mm0		;# mm3 has both vcoul 
	;# update vctot 
	pfadd mm3, [esp + nb100nf_vctot]      ;# add the earlier value  
	movq [esp + nb100nf_vctot], mm3       ;# store the sum 
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb100nf_innerk],  2
	jl    .nb100nf_finish_inner
	jmp   .nb100nf_unroll_loop
.nb100nf_finish_inner:	
	and dword ptr [esp + nb100nf_innerk],  1
	jnz  .nb100nf_single_inner
	jmp  .nb100nf_updateouterdata		
.nb100nf_single_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + nb100nf_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + nb100nf_charge]
	pxor mm6, mm6
	movd mm6, [esp + nb100nf_iq]
	movd mm7, [ecx + eax*4]
	pfmul mm6, mm7	  	;# mm6=qq 
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esp + nb100nf_ix]
	movd  mm1, [esp + nb100nf_iz]
	movq  mm2, [esi + eax*4]
	movd  mm3, [esi + eax*4 + 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfacc mm0, mm1		;# mm0=rsq 
	
	pfrsqrt mm1,mm0
	movq mm2,mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	;# calculate potential 
	pfmul mm6, mm1		;# mm6=vcoul 
	;# update vctot 
	pfadd mm6, [esp + nb100nf_vctot]
	movq [esp + nb100nf_vctot], mm6
	;# done! 
.nb100nf_updateouterdata:	
	;# get n from stack
	mov esi, [esp + nb100nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb100nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb100nf_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb100nf_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
       	;# finish if last 
        mov ecx, [esp + nb100nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb100nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb100nf_n], esi
        jmp .nb100nf_outer
.nb100nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb100nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb100nf_end
        ;# non-zero, do one more workunit
        jmp   .nb100nf_threadloop
.nb100nf_end:
	femms

	mov eax, [esp + nb100nf_nouter] 	
	mov ebx, [esp + nb100nf_ninner]
	mov ecx, [ebp + nb100nf_outeriter]
	mov edx, [ebp + nb100nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx
	
	add esp, 68
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret





