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

mm_six:
        .long 0x40c00000
        .long 0x40c00000
mm_twelve:
        .long 0x41400000
        .long 0x41400000



.globl nb_kernel010_ia32_3dnow
.globl _nb_kernel010_ia32_3dnow
nb_kernel010_ia32_3dnow:	
_nb_kernel010_ia32_3dnow:	
.equiv		nb010_p_nri,		8
.equiv		nb010_iinr,		12
.equiv		nb010_jindex,		16
.equiv		nb010_jjnr,		20
.equiv		nb010_shift,		24
.equiv		nb010_shiftvec,		28
.equiv		nb010_fshift,		32
.equiv		nb010_gid,		36
.equiv		nb010_pos,		40		
.equiv		nb010_faction,		44
.equiv		nb010_charge,		48
.equiv		nb010_p_facel,		52
.equiv		nb010_p_krf,		56	
.equiv		nb010_p_crf,		60	
.equiv		nb010_Vc,		64	
.equiv		nb010_type,		68
.equiv		nb010_p_ntype,		72
.equiv		nb010_vdwparam,		76	
.equiv		nb010_Vvdw,		80	
.equiv		nb010_p_tabscale,	84	
.equiv		nb010_VFtab,		88
.equiv		nb010_invsqrta,		92	
.equiv		nb010_dvda,		96
.equiv          nb010_p_gbtabscale,     100
.equiv          nb010_GBtab,            104
.equiv          nb010_p_nthreads,       108
.equiv          nb010_count,            112
.equiv          nb010_mtx,              116
.equiv          nb010_outeriter,        120
.equiv          nb010_inneriter,        124
.equiv          nb010_work,             128
	;# stack offsets for local variables 
.equiv		nb010_is3,		0
.equiv		nb010_ii3,		4
.equiv		nb010_ix,		8
.equiv		nb010_iy,		12
.equiv		nb010_iz,		16
.equiv		nb010_Vvdwtot,		20  
.equiv		nb010_c6,		28  
.equiv		nb010_c12,		36  
.equiv		nb010_six,		44  
.equiv		nb010_twelve,		52  
.equiv		nb010_ntia,		60
.equiv		nb010_innerjjnr,	64
.equiv		nb010_innerk,		68		
.equiv		nb010_fix,		72
.equiv		nb010_fiy,		76
.equiv		nb010_fiz,		80
.equiv		nb010_dx1,		84
.equiv		nb010_dy1,		88
.equiv		nb010_dz1,		92
.equiv		nb010_dx2,		96
.equiv		nb010_dy2,		100
.equiv		nb010_dz2,		104			
.equiv          nb010_n,                108 ;# idx for outer loop
.equiv          nb010_nn1,              112 ;# number of outer iterations
.equiv          nb010_nri,              116
.equiv          nb010_ntype,            120
.equiv          nb010_nouter,           124
.equiv          nb010_ninner,           128

	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 132						;# local stack space 
	femms

	mov eax, 0x40c00000 ;# fp 6.0
	mov ebx, 0x41400000 ;# fp 12.0
	
        mov [esp + nb010_six], eax
        mov [esp + nb010_six + 4], eax
        mov [esp + nb010_twelve], ebx
        mov [esp + nb010_twelve +4], ebx

	mov ecx, [ebp + nb010_p_nri]
	mov edx, [ebp + nb010_p_ntype]
	mov ecx, [ecx]
	mov edx, [edx]
	mov [esp + nb010_nri], ecx
	mov [esp + nb010_ntype], edx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb010_nouter], eax
	mov [esp + nb010_ninner], eax

.nb010_threadloop:
        mov   esi, [ebp + nb010_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb010_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 10                           ;# ebx=nn1=nn0+10
        lock cmpxchg [esi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb010_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb010_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb010_n], eax
        mov [esp + nb010_nn1], ebx
	sub ebx, eax           			;# calc number of outer lists in ecx
	mov esi, eax				;# copy n to esi
	
        jg  .nb010_outerstart
        jmp .nb010_end

.nb010_outerstart:	
	;# ebx contains number of outer iterations
	add ebx, [esp + nb010_nouter]
        mov [esp + nb010_nouter], ebx
	
.nb010_outer:
		
	mov   eax, [ebp + nb010_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]			;# ebx=3*is 
	mov   [esp + nb010_is3],ebx    		;# store is3 

	mov   eax, [ebp + nb010_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]			;# move shX/shY to mm0 and shZ to mm1. 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + nb010_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]					;# ebx =ii 

	mov   edx, [ebp + nb010_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [esp + nb010_ntype]
	shl   edx, 1
	mov   [esp + nb010_ntia], edx

	lea   ebx, [ebx + ebx*2]			;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb010_pos]		;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]			;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]		;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb010_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + nb010_ix], mm0	
	movd  [esp + nb010_iz], mm1	
				
	;# clear total potential and i forces 
	pxor  mm7,mm7
	movq  [esp + nb010_Vvdwtot], mm7
	movq  [esp + nb010_fix],    mm7
	movd  [esp + nb010_fiz],    mm7

	mov   eax, [ebp + nb010_jindex]
	mov   ecx, [eax + esi*4]					;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]				;# jindex[n+1] 
	sub   edx, ecx						;# number of innerloop atoms 

	mov   esi, [ebp + nb010_pos]
	mov   edi, [ebp + nb010_faction]	
	mov   eax, [ebp + nb010_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb010_innerjjnr], eax     ;#  pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb010_ninner]
	mov   [esp + nb010_ninner], ecx
	mov   [esp + nb010_innerk], edx    ;# number of innerloop atoms 
	add   edx, 0
	jge   .nb010_unroll_loop
	jmp   .nb010_finish_inner
.nb010_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + nb010_innerjjnr]	;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]				;# eax/ebx=jnr 
	add   dword ptr [esp + nb010_innerjjnr],  8	;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]				;# prefetch data - trial and error says 16 is best 	

	mov ecx, [ebp + nb010_type]
	mov edx, [ecx + eax*4]        		;# type [jnr1] 
	mov ecx, [ecx + ebx*4]				;# type [jnr2] 

	mov esi, [ebp + nb010_vdwparam]			;# base of vdwparam  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + nb010_ntia]			;# tja = ntia + 2*type 
	add ecx, [esp + nb010_ntia]

	movq mm5, [esi + edx*4]				;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]				;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7					;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7					;# mm6 = 1st c12 / 2nd c12 
	movq [esp + nb010_c6], mm5
	movq [esp + nb010_c12], mm6

	lea   eax, [eax + eax*2]			;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + nb010_pos]

	movq  mm0, [esp + nb010_ix]
	movd  mm1, [esp + nb010_iz]	 	
	movq  mm4, [esi + eax*4]			;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0						;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + nb010_dx1], mm4	    ;# store dr 
	movd  [esp + nb010_dz1], mm5
	pfmul mm4,mm4						;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5						;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5						;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]			;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0						;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + nb010_dx2], mm6	    ;# store dr 
	movd  [esp + nb010_dz2], mm7
	pfmul mm6,mm6						;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7						;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7						;# second rsq in lower mm6 

    	pfrcp mm0, mm4						;# lookup reciprocal seed  
    	pfrcp mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        			;# now 4 has rsq and 0 the seed for both pairs. 
                  	        			;# amd 3dnow N-R iteration to get full precision. 
    	pfrcpit1 mm4,mm0				
    	pfrcpit2 mm4,mm0	
	;# mm4 now contains invsq,
	;# do potential and fscal
	 
	movq  mm0, mm4
	pfmul mm4, mm0
	pfmul mm4, mm0             			;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5						;# mm5=rinvtwelve 

	pfmul mm5, [esp + nb010_c12]
	pfmul mm4, [esp + nb010_c6]	
	movq mm6, mm5						;# mm6 is Vvdw12-Vvdw6  
	pfsub mm6, mm4

	pfmul mm4, [esp + nb010_six]

	pfmul mm5, [esp + nb010_twelve]
	pfsub mm5,mm4
 	pfmul mm0, mm5						;# mm0 is total fscal now 	

	prefetchw [esp + nb010_dx1]			;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]				;# prefetch the 1st faction to cache 
	movq mm2,  [esp + nb010_dx1]		;# fetch dr 
	movd mm3,  [esp + nb010_dz1]

	;# update Vvdwtot  
	pfadd mm6, [esp + nb010_Vvdwtot]     ;# add the earlier value 
	movq [esp + nb010_Vvdwtot], mm6      ;# store the sum 

	prefetchw [edi + ebx*4]				;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0						;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + nb010_dx2] 		;# fetch dr 
	movd mm5,  [esp + nb010_dz2]
	pfmul mm4, mm1   					;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + nb010_fix]
	movd mm1,  [esp + nb010_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + nb010_fix], mm0
	movd [esp + nb010_fiz], mm1
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
	sub  dword ptr [esp + nb010_innerk],  2
	jl    .nb010_finish_inner
	jmp   .nb010_unroll_loop
.nb010_finish_inner:	
	and dword ptr [esp + nb010_innerk],  1
	jnz  .nb010_single_inner
	jmp  .nb010_updateouterdata		
.nb010_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + nb010_innerjjnr]
	mov   eax, [eax]					;# eax=jnr offset 

	mov esi, [ebp + nb010_vdwparam]
	mov ecx, [ebp + nb010_type]
	mov edx, [ecx + eax*4]        		;# type [jnr1] 
	shl edx, 1
	add edx, [esp + nb010_ntia]			;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]				;# mm5 = 1st c6 		
	movq [esp + nb010_c6], mm5
	movd mm5, [esi + edx*4 + 4]			;# mm5 = 1st c12 		
	movq [esp + nb010_c12], mm5

	mov   esi, [ebp + nb010_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + nb010_ix]
	movd  mm1, [esp + nb010_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + nb010_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + nb010_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5						;# mm4=rsq 
	
    	pfrcp mm0,mm4
    	pfrcpit1 mm4,mm0				
    	pfrcpit2 mm4,mm0					;# mm4=invsq 
	;# calculate potentials and scalar force 
	movq  mm0, mm4

	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + nb010_c12]
	pfmul mm4, [esp + nb010_c6]	
	movq mm6, mm5	;# mm6 is Vvdw12-Vvdw6 
	pfsub mm6, mm4

	pfmul mm4, [esp + nb010_six]

	pfmul mm5, [esp + nb010_twelve]
	pfsub mm5, mm4
 	pfmul mm0, mm5    ;# mm0 is total fscal now 

	;# update Vvdwtot 
	pfadd mm6, [esp + nb010_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb010_Vvdwtot], mm6       ;# store the sum   

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache 
	movq mm2,  [esp + nb010_dx1]
	movd mm3,  [esp + nb010_dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + nb010_fix]
	movd mm1,  [esp + nb010_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + nb010_fix], mm0
	movd [esp + nb010_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.nb010_updateouterdata:	
	mov   ecx, [esp + nb010_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + nb010_fix]
	pfadd mm7, [esp + nb010_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + nb010_fshift]    ;# increment fshift force 
	mov   edx, [esp + nb010_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + nb010_fix]
	pfadd mm7, [esp + nb010_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# get n from stack
	mov esi, [esp + nb010_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb010_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb010_Vvdwtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb010_Vvdw]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment Vvdw[gid] 

       ;# finish if last 
        mov ecx, [esp + nb010_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb010_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb010_n], esi
        jmp .nb010_outer
.nb010_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb010_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb010_end
        ;# non-zero, do one more workunit
        jmp   .nb010_threadloop
.nb010_end:
	femms

	mov eax, [esp + nb010_nouter] 	
	mov ebx, [esp + nb010_ninner]
	mov ecx, [ebp + nb010_outeriter]
	mov edx, [ebp + nb010_inneriter]
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
	




.globl nb_kernel010nf_ia32_3dnow
.globl _nb_kernel010nf_ia32_3dnow
nb_kernel010nf_ia32_3dnow:	
_nb_kernel010nf_ia32_3dnow:	
.equiv		nb010nf_p_nri,		8
.equiv		nb010nf_iinr,		12
.equiv		nb010nf_jindex,		16
.equiv		nb010nf_jjnr,		20
.equiv		nb010nf_shift,		24
.equiv		nb010nf_shiftvec,	28
.equiv		nb010nf_fshift,		32
.equiv		nb010nf_gid,		36
.equiv		nb010nf_pos,		40		
.equiv		nb010nf_faction,	44
.equiv		nb010nf_charge,		48
.equiv		nb010nf_p_facel,	52
.equiv		nb010nf_p_krf,		56	
.equiv		nb010nf_p_crf,		60	
.equiv		nb010nf_Vc,		64	
.equiv		nb010nf_type,		68
.equiv		nb010nf_p_ntype,	72
.equiv		nb010nf_vdwparam,	76	
.equiv		nb010nf_Vvdw,		80	
.equiv		nb010nf_p_tabscale,	84	
.equiv		nb010nf_VFtab,		88
.equiv		nb010nf_invsqrta,	92	
.equiv		nb010nf_dvda,		96
.equiv          nb010nf_p_gbtabscale,   100
.equiv          nb010nf_GBtab,          104
.equiv          nb010nf_p_nthreads,     108
.equiv          nb010nf_count,          112
.equiv          nb010nf_mtx,            116
.equiv          nb010nf_outeriter,      120
.equiv          nb010nf_inneriter,      124
.equiv          nb010nf_work,           128
	;# stack offsets for local variables 
.equiv		nb010nf_is3,		0
.equiv		nb010nf_ii3,		4
.equiv		nb010nf_ix,		8
.equiv		nb010nf_iy,		12
.equiv		nb010nf_iz,		16
.equiv		nb010nf_Vvdwtot,	20  
.equiv		nb010nf_c6,		28  
.equiv		nb010nf_c12,		36  
.equiv		nb010nf_ntia,		44
.equiv		nb010nf_innerjjnr,	48
.equiv		nb010nf_innerk,		52	
.equiv          nb010nf_n,              56 ;# idx for outer loop
.equiv          nb010nf_nn1,            60 ;# number of outer iterations
.equiv          nb010nf_nri,            64
.equiv          nb010nf_ntype,          68
.equiv          nb010nf_nouter,         72
.equiv          nb010nf_ninner,         76
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx 
    	push edx
	push esi
	push edi
	sub esp, 80		;# local stack space 
	femms

	mov ecx, [ebp + nb010nf_p_nri]
	mov edx, [ebp + nb010nf_p_ntype]
	mov ecx, [ecx]
	mov edx, [edx]
	mov [esp + nb010nf_nri], ecx
	mov [esp + nb010nf_ntype], edx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb010nf_nouter], eax
	mov [esp + nb010nf_ninner], eax

.nb010nf_threadloop:
        mov   esi, [ebp + nb010nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb010nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 10                           ;# ebx=nn1=nn0+10
        lock cmpxchg [esi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb010nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb010nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb010nf_n], eax
        mov [esp + nb010nf_nn1], ebx

        sub ebx, eax		; # calc number of outer lists in ecx
        mov esi, eax	; # copy n to esi

        jg  .nb010nf_outerstart
        jmp .nb010nf_end

.nb010nf_outerstart:
	;; # ebx contains number of outer iterations
        add ebx, [esp + nb010nf_nouter]
        mov [esp + nb010nf_nouter], ebx
	
.nb010nf_outer:
	mov   eax, [ebp + nb010nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb010nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb010nf_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1. 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + nb010nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb010nf_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [esp + nb010nf_ntype]
	shl   edx, 1
	mov   [esp + nb010nf_ntia], edx

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb010nf_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb010nf_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + nb010nf_ix], mm0	
	movd  [esp + nb010nf_iz], mm1	
				
	;# clear total potential 
	pxor  mm7,mm7
	movq  [esp + nb010nf_Vvdwtot], mm7

	mov   eax, [ebp + nb010nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb010nf_pos]	
	mov   eax, [ebp + nb010nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb010nf_innerjjnr], eax     ;#  pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb010nf_ninner]
	mov   [esp + nb010nf_ninner], ecx
	mov   [esp + nb010nf_innerk], edx    ;# number of innerloop atoms 
	add   edx, 0
	jge   .nb010nf_unroll_loop
	jmp   .nb010nf_finish_inner
.nb010nf_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + nb010nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + nb010nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 	

	mov ecx, [ebp + nb010nf_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + nb010nf_vdwparam]		;# base of vdwparam  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + nb010nf_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + nb010nf_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + nb010nf_c6], mm5
	movq [esp + nb010nf_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + nb010nf_pos]

	movq  mm0, [esp + nb010nf_ix]
	movd  mm1, [esp + nb010nf_iz]	 	
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

    	pfrcp mm0, mm4	             ;# lookup reciprocal seed  
    	pfrcp mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
                  	        	;# amd 3dnow N-R iteration to get full precision. 
    	pfrcpit1 mm4,mm0				
    	pfrcpit2 mm4,mm0	
	;# mm4 now contains invsq,
	 ;# do potential and fscal
	 
	movq  mm0, mm4
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + nb010nf_c12]
	pfmul mm4, [esp + nb010nf_c6]	
	movq mm6, mm5	;# mm6 is Vvdw12-Vvdw6  
	pfsub mm6, mm4
	;# update Vvdwtot  
	pfadd mm6, [esp + nb010nf_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb010nf_Vvdwtot], mm6       ;# store the sum 
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb010nf_innerk],  2
	jl    .nb010nf_finish_inner
	jmp   .nb010nf_unroll_loop
.nb010nf_finish_inner:	
	and dword ptr [esp + nb010nf_innerk],  1
	jnz  .nb010nf_single_inner
	jmp  .nb010nf_updateouterdata
.nb010nf_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + nb010nf_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + nb010nf_vdwparam]
	mov ecx, [ebp + nb010nf_type]
	mov edx, [ecx + eax*4]        	;# type [jnr1] 
	shl edx, 1
	add edx, [esp + nb010nf_ntia]	    ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6 		
	movq [esp + nb010nf_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12 		
	movq [esp + nb010nf_c12], mm5

	mov   esi, [ebp + nb010nf_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + nb010nf_ix]
	movd  mm1, [esp + nb010nf_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm4=rsq 
	
    	pfrcp mm0,mm4
    	pfrcpit1 mm4,mm0				
    	pfrcpit2 mm4,mm0	;# mm4=invsq 
	;# calculate potentials and scalar force 
	movq  mm0, mm4

	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + nb010nf_c12]
	pfmul mm4, [esp + nb010nf_c6]	
	movq mm6, mm5	;# mm6 is Vvdw12-Vvdw6 
	pfsub mm6, mm4
	;# update Vvdwtot 
	pfadd mm6, [esp + nb010nf_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb010nf_Vvdwtot], mm6       ;# store the sum   

.nb010nf_updateouterdata:	
	;# get n from stack
	mov esi, [esp + nb010nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb010nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb010nf_Vvdwtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb010nf_Vvdw]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment Vvdw[gid] 

       ;# finish if last 
        mov ecx, [esp + nb010nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb010nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb010nf_n], esi
        jmp .nb010nf_outer
.nb010nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb010nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb010nf_end
        ;# non-zero, do one more workunit
        jmp   .nb010nf_threadloop
.nb010nf_end:
	femms

        mov eax, [esp + nb010nf_nouter]
	mov ebx, [esp + nb010nf_ninner]
	mov ecx, [ebp + nb010nf_outeriter]
	mov edx, [ebp + nb010nf_inneriter]
	mov [ecx], eax
        mov [edx], ebx
		
	add esp, 80
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret
	



