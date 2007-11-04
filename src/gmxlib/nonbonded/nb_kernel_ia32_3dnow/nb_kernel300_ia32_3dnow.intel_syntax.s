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



.globl nb_kernel300_ia32_3dnow
.globl _nb_kernel300_ia32_3dnow
nb_kernel300_ia32_3dnow:	
_nb_kernel300_ia32_3dnow:	
.equiv		nb300_p_nri,		8
.equiv		nb300_iinr,		12
.equiv		nb300_jindex,		16
.equiv		nb300_jjnr,		20
.equiv		nb300_shift,		24
.equiv		nb300_shiftvec,		28
.equiv		nb300_fshift,		32
.equiv		nb300_gid,		36
.equiv		nb300_pos,		40		
.equiv		nb300_faction,		44
.equiv		nb300_charge,		48
.equiv		nb300_p_facel,		52
.equiv		nb300_p_krf,		56	
.equiv		nb300_p_crf,		60	
.equiv		nb300_Vc,		64	
.equiv		nb300_type,		68
.equiv		nb300_p_ntype,		72
.equiv		nb300_vdwparam,		76	
.equiv		nb300_Vvdw,		80	
.equiv		nb300_p_tabscale,	84	
.equiv		nb300_VFtab,		88
.equiv		nb300_invsqrta,		92	
.equiv		nb300_dvda,		96
.equiv          nb300_p_gbtabscale,     100
.equiv          nb300_GBtab,            104
.equiv          nb300_p_nthreads,       108
.equiv          nb300_count,            112
.equiv          nb300_mtx,              116
.equiv          nb300_outeriter,        120
.equiv          nb300_inneriter,        124
.equiv          nb300_work,             128
	;# stack offsets for local variables 
.equiv		nb300_is3,		0
.equiv		nb300_ii3,		4
.equiv		nb300_ix,		8
.equiv		nb300_iy,		12
.equiv		nb300_iz,		16
.equiv		nb300_iq,		20 
.equiv		nb300_vctot,		28 
.equiv		nb300_two,		36 
.equiv		nb300_n1,		44 
.equiv		nb300_tsc,		52 
.equiv		nb300_ntia,		60
.equiv		nb300_innerjjnr,	64
.equiv		nb300_innerk,		68		
.equiv		nb300_fix,		72
.equiv		nb300_fiy,		76
.equiv		nb300_fiz,		80
.equiv		nb300_dx1,		84
.equiv		nb300_dy1,		88
.equiv		nb300_dz1,		92
.equiv		nb300_dx2,		96
.equiv		nb300_dy2,		100
.equiv		nb300_dz2,		104						
.equiv          nb300_n,                108 ;# idx for outer loop
.equiv          nb300_nn1,              112 ;# number of outer iterations
.equiv          nb300_nri,              116
.equiv          nb300_facel,            120
.equiv          nb300_nouter,           124
.equiv          nb300_ninner,           128
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 132		;# local stack space 
	femms
	;# move data to local stack  
	mov ecx, [ebp + nb300_p_nri]
	mov esi, [ebp + nb300_p_facel]
	mov edi, [ebp + nb300_p_tabscale]
	mov ecx, [ecx]
	mov esi, [esi]
	mov [esp + nb300_nri], ecx
	mov [esp + nb300_facel], esi

	movd  mm3, [edi]
	punpckldq mm3,mm3
	movq  [esp + nb300_tsc], mm3	
	mov eax, 0x40000000
	mov [esp + nb300_two], eax
	mov [esp + nb300_two+4], eax

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb300_nouter], eax
	mov [esp + nb300_ninner], eax

.nb300_threadloop:
        mov   esi, [ebp + nb300_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb300_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb300_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb300_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb300_n], eax
        mov [esp + nb300_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb300_outerstart
        jmp .nb300_end

.nb300_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb300_nouter]
        mov [esp + nb300_nouter], ebx
	
.nb300_outer:
	mov   eax, [ebp + nb300_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb300_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb300_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + nb300_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx=ii 

	mov   edx, [ebp + nb300_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [esp + nb300_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + nb300_iq], mm2	    ;# iq =facel*charge[ii] 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb300_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb300_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + nb300_ix], mm0	
	movd  [esp + nb300_iz], mm1	
				
	;# clear total potential and i forces 
	pxor  mm7,mm7
	movq  [esp + nb300_vctot],  mm7
	movq  [esp + nb300_fix],    mm7
	movd  [esp + nb300_fiz],    mm7

	mov   eax, [ebp + nb300_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb300_pos]
	mov   edi, [ebp + nb300_faction]	
	mov   eax, [ebp + nb300_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb300_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb300_ninner]
	mov   [esp + nb300_ninner], ecx
	mov   [esp + nb300_innerk], edx    ;# number of innerloop atoms 
	add   edx, 0
	jge   .nb300_unroll_loop
	jmp   .nb300_finish_inner
.nb300_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + nb300_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + nb300_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + nb300_charge]    ;# base of charge[] 
	movq mm5, [esp + nb300_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + nb300_pos]

	movq  mm0, [esp + nb300_ix]
	movd  mm1, [esp + nb300_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + nb300_dx1], mm4	     ;# store dr 
	movd  [esp + nb300_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + nb300_dx2], mm6	     ;# store dr 
	movd  [esp + nb300_dz2], mm7
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
	pfmul mm1, [esp + nb300_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + nb300_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + nb300_VFtab]
	mov ecx, [esp + nb300_n1]
	shl ecx, 2
	;# coulomb table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + nb300_n1+4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + nb300_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  

	;# at this point mm5 contains vcoul and mm3 fijC. 
	;# increment vcoul - then we can get rid of mm5. 
	;# update vctot 
	pfadd mm5, [esp + nb300_vctot]      ;# add the earlier value 
	movq [esp + nb300_vctot], mm5       ;# store the sum       

	;# change sign of mm3 
    pxor mm1,mm1
	pfsub mm1, mm3
	pfmul mm1, [esp + nb300_tsc]	
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	prefetchw [esp + nb300_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + nb300_dx1]	;# fetch dr 
	movd mm3,  [esp + nb300_dz1]

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + nb300_dx2] 	;# fetch dr 
	movd mm5,  [esp + nb300_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + nb300_fix]
	movd mm1,  [esp + nb300_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + nb300_fix], mm0
	movd [esp + nb300_fiz], mm1
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
	sub dword ptr [esp + nb300_innerk],  2
	jl    .nb300_finish_inner
	jmp   .nb300_unroll_loop
.nb300_finish_inner:	
	and dword ptr [esp + nb300_innerk],  1
	jnz  .nb300_single_inner
	jmp  .nb300_updateouterdata		
.nb300_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + nb300_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + nb300_charge]
	movd mm5, [esp + nb300_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov   esi, [ebp + nb300_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + nb300_ix]
	movd  mm1, [esp + nb300_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + nb300_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + nb300_dz1], mm5	
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
	;# mm0 is invsqrt, and mm1 r. 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + nb300_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + nb300_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + nb300_VFtab]
	mov ecx, [esp + nb300_n1]
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

	pfmul mm7, [esp + nb300_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  

	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + nb300_vctot]      ;# add the earlier value 
	movq [esp + nb300_vctot], mm5       ;# store the sum       
	
	;# change sign of mm3 
    pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + nb300_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + nb300_dx1]
	movd mm3,  [esp + nb300_dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + nb300_fix]
	movd mm1,  [esp + nb300_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + nb300_fix], mm0
	movd [esp + nb300_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.nb300_updateouterdata:	
	mov   ecx, [esp + nb300_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + nb300_fix]
	pfadd mm7, [esp + nb300_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + nb300_fshift]    ;# increment fshift force 
	mov   edx, [esp + nb300_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + nb300_fix]
	pfadd mm7, [esp + nb300_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# get n from stack
	mov esi, [esp + nb300_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb300_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb300_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb300_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

       	;# finish if last 
        mov ecx, [esp + nb300_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb300_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb300_n], esi
        jmp .nb300_outer
.nb300_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb300_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb300_end
        ;# non-zero, do one more workunit
        jmp   .nb300_threadloop
.nb300_end:
	femms

	mov eax, [esp + nb300_nouter] 	
	mov ebx, [esp + nb300_ninner]
	mov ecx, [ebp + nb300_outeriter]
	mov edx, [ebp + nb300_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	add esp, 132
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret



.globl nb_kernel300nf_ia32_3dnow
.globl _nb_kernel300nf_ia32_3dnow
nb_kernel300nf_ia32_3dnow:	
_nb_kernel300nf_ia32_3dnow:	
.equiv		nb300nf_p_nri,		8
.equiv		nb300nf_iinr,		12
.equiv		nb300nf_jindex,		16
.equiv		nb300nf_jjnr,		20
.equiv		nb300nf_shift,		24
.equiv		nb300nf_shiftvec,	28
.equiv		nb300nf_fshift,		32
.equiv		nb300nf_gid,		36
.equiv		nb300nf_pos,		40		
.equiv		nb300nf_faction,	44
.equiv		nb300nf_charge,		48
.equiv		nb300nf_p_facel,	52
.equiv		nb300nf_p_krf,		56	
.equiv		nb300nf_p_crf,		60	
.equiv		nb300nf_Vc,		64	
.equiv		nb300nf_type,		68
.equiv		nb300nf_p_ntype,	72
.equiv		nb300nf_vdwparam,	76	
.equiv		nb300nf_Vvdw,		80	
.equiv		nb300nf_p_tabscale,	84	
.equiv		nb300nf_VFtab,		88
.equiv		nb300nf_invsqrta,	92	
.equiv		nb300nf_dvda,		96
.equiv          nb300nf_p_gbtabscale,   100
.equiv          nb300nf_GBtab,          104
.equiv          nb300nf_p_nthreads,     108
.equiv          nb300nf_count,          112
.equiv          nb300nf_mtx,            116
.equiv          nb300nf_outeriter,      120
.equiv          nb300nf_inneriter,      124
.equiv          nb300nf_work,           128
	;# stack offsets for local variables 
.equiv		nb300nf_is3,		0
.equiv		nb300nf_ii3,		4
.equiv		nb300nf_ix,		8
.equiv		nb300nf_iy,		12
.equiv		nb300nf_iz,		16
.equiv		nb300nf_iq,		20 
.equiv		nb300nf_vctot,		28 
.equiv		nb300nf_n1,		36
.equiv		nb300nf_tsc,		44 
.equiv		nb300nf_ntia,		52
.equiv		nb300nf_innerjjnr,	56
.equiv		nb300nf_innerk,		60						
.equiv          nb300nf_n,              64 ;# idx for outer loop
.equiv          nb300nf_nn1,            68 ;# number of outer iterations
.equiv          nb300nf_nri,            72
.equiv          nb300nf_facel,          76
.equiv          nb300nf_nouter,         80
.equiv          nb300nf_ninner,         84
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 88		;# local stack space 
	femms
	;# move data to local stack  
	mov ecx, [ebp + nb300nf_p_nri]
	mov esi, [ebp + nb300nf_p_facel]
	mov edi, [ebp + nb300nf_p_tabscale]
	mov ecx, [ecx]
	mov esi, [esi]
	mov [esp + nb300nf_nri], ecx
	mov [esp + nb300nf_facel], esi

	movd  mm3, [edi]
	punpckldq mm3,mm3
	movq  [esp + nb300nf_tsc], mm3	

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb300nf_nouter], eax
	mov [esp + nb300nf_ninner], eax

.nb300nf_threadloop:
        mov   esi, [ebp + nb300nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb300nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb300nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb300nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb300nf_n], eax
        mov [esp + nb300nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb300nf_outerstart
        jmp .nb300nf_end

.nb300nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb300nf_nouter]
        mov [esp + nb300nf_nouter], ebx
	
.nb300nf_outer:
	mov   eax, [ebp + nb300nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb300nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb300nf_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + nb300nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx=ii 

	mov   edx, [ebp + nb300nf_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [esp + nb300nf_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + nb300nf_iq], mm2	    ;# iq =facel*charge[ii] 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb300nf_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb300nf_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + nb300nf_ix], mm0	
	movd  [esp + nb300nf_iz], mm1	
				
	;# clear total potential and i forces 
	pxor  mm7,mm7
	movq  [esp + nb300nf_vctot],  mm7

	mov   eax, [ebp + nb300nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb300nf_pos]	
	mov   eax, [ebp + nb300nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb300nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb300nf_ninner]
	mov   [esp + nb300nf_ninner], ecx
	mov   [esp + nb300nf_innerk], edx    ;# number of innerloop atoms 
	add   edx, 0
	jge   .nb300nf_unroll_loop
	jmp   .nb300nf_finish_inner
.nb300nf_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + nb300nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + nb300nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + nb300nf_charge]    ;# base of charge[] 
	movq mm5, [esp + nb300nf_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + nb300nf_pos]

	movq  mm0, [esp + nb300nf_ix]
	movd  mm1, [esp + nb300nf_iz]	 	
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
	pfmul mm1, [esp + nb300nf_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + nb300nf_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + nb300nf_VFtab]
	mov ecx, [esp + nb300nf_n1]
	shl ecx, 2
	;# coulomb table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + nb300nf_n1+4]
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
	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5. 
	;# update vctot 
	pfadd mm5, [esp + nb300nf_vctot]      ;# add the earlier value 
	movq [esp + nb300nf_vctot], mm5       ;# store the sum       
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb300nf_innerk],  2
	jl    .nb300nf_finish_inner
	jmp   .nb300nf_unroll_loop
.nb300nf_finish_inner:	
	and dword ptr [esp + nb300nf_innerk],  1
	jnz  .nb300nf_single_inner
	jmp  .nb300nf_updateouterdata		
.nb300nf_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + nb300nf_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + nb300nf_charge]
	movd mm5, [esp + nb300nf_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov   esi, [ebp + nb300nf_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + nb300nf_ix]
	movd  mm1, [esp + nb300nf_iz]
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
	;# mm0 is invsqrt, and mm1 r. 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + nb300nf_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + nb300nf_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + nb300nf_VFtab]
	mov ecx, [esp + nb300nf_n1]
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
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + nb300nf_vctot]      ;# add the earlier value 
	movq [esp + nb300nf_vctot], mm5       ;# store the sum       
	
.nb300nf_updateouterdata:	
	;# get n from stack
	mov esi, [esp + nb300nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb300nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb300nf_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb300nf_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

       	;# finish if last 
        mov ecx, [esp + nb300nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb300nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb300nf_n], esi
        jmp .nb300nf_outer
.nb300nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb300nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb300nf_end
        ;# non-zero, do one more workunit
        jmp   .nb300nf_threadloop
.nb300nf_end:
	femms

	mov eax, [esp + nb300nf_nouter] 	
	mov ebx, [esp + nb300nf_ninner]
	mov ecx, [ebp + nb300nf_outeriter]
	mov edx, [ebp + nb300nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	add esp, 88
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret



