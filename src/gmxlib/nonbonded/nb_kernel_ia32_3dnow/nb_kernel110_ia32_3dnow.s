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



.globl nb_kernel110_ia32_3dnow
.globl _nb_kernel110_ia32_3dnow
nb_kernel110_ia32_3dnow:	
_nb_kernel110_ia32_3dnow:	
.equiv		nb110_p_nri,		8
.equiv		nb110_iinr,		12
.equiv		nb110_jindex,		16
.equiv		nb110_jjnr,		20
.equiv		nb110_shift,		24
.equiv		nb110_shiftvec,		28
.equiv		nb110_fshift,		32
.equiv		nb110_gid,		36
.equiv		nb110_pos,		40		
.equiv		nb110_faction,		44
.equiv		nb110_charge,		48
.equiv		nb110_p_facel,		52
.equiv		nb110_p_krf,		56	
.equiv		nb110_p_crf,		60	
.equiv		nb110_Vc,		64	
.equiv		nb110_type,		68
.equiv		nb110_p_ntype,		72
.equiv		nb110_vdwparam,		76	
.equiv		nb110_Vvdw,		80	
.equiv		nb110_p_tabscale,	84	
.equiv		nb110_VFtab,		88
.equiv		nb110_invsqrta,		92	
.equiv		nb110_dvda,		96
.equiv          nb110_p_gbtabscale,     100
.equiv          nb110_GBtab,            104
.equiv          nb110_p_nthreads,       108
.equiv          nb110_count,            112
.equiv          nb110_mtx,              116
.equiv          nb110_outeriter,        120
.equiv          nb110_inneriter,        124
.equiv          nb110_work,             128
	;# stack offsets for local variables 
.equiv		nb110_is3,		0
.equiv		nb110_ii3,		4
.equiv		nb110_ix,		8
.equiv		nb110_iy,		12
.equiv		nb110_iz,		16
.equiv		nb110_iq,		20 
.equiv		nb110_vctot,		28 
.equiv		nb110_Vvdwtot,		36 
.equiv		nb110_c6,		44 
.equiv		nb110_c12,		52 
.equiv		nb110_six,		60 
.equiv		nb110_twelve,		68 
.equiv		nb110_ntia,		76
.equiv		nb110_innerjjnr,	80
.equiv		nb110_innerk,		84		
.equiv		nb110_fix,		88
.equiv		nb110_fiy,		92
.equiv		nb110_fiz,		96
.equiv		nb110_dx1,		100
.equiv		nb110_dy1,		104
.equiv		nb110_dz1,		108
.equiv		nb110_dx2,		112
.equiv		nb110_dy2,		116
.equiv		nb110_dz2,		120						
.equiv          nb110_n,                124 ;# idx for outer loop
.equiv          nb110_nn1,              128 ;# number of outer iterations
.equiv          nb110_nri,              132
.equiv          nb110_facel,            136
.equiv          nb110_ntype,            140
.equiv          nb110_nouter,           144
.equiv          nb110_ninner,           148

	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 152		;# local stack space 
	femms
	;# move data to local stack  
	mov ecx, [ebp + nb110_p_nri]
	mov edx, [ebp + nb110_p_ntype]
	mov esi, [ebp + nb110_p_facel]
	mov ecx, [ecx]
	mov edx, [edx]
	mov esi, [esi]
	mov [esp + nb110_nri], ecx
	mov [esp + nb110_ntype], edx
	mov [esp + nb110_facel], esi

	mov eax, 0x40c00000 ;# fp 6.0
	mov ebx, 0x41400000 ;# fp 12.0
	
        mov [esp + nb110_six], eax
        mov [esp + nb110_six + 4], eax
        mov [esp + nb110_twelve], ebx
        mov [esp + nb110_twelve +4], ebx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb110_nouter], eax
	mov [esp + nb110_ninner], eax

.nb110_threadloop:
        mov   esi, [ebp + nb110_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb110_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [esi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb110_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb110_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb110_n], eax
        mov [esp + nb110_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi

        jg  .nb110_outerstart
        jmp .nb110_end

.nb110_outerstart:	
	;# ebx contains number of outer iterations
	add ebx, [esp + nb110_nouter]
        mov [esp + nb110_nouter], ebx
	
.nb110_outer:
	mov   eax, [ebp + nb110_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb110_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb110_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + nb110_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx=ii 

	mov   edx, [ebp + nb110_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [esp + nb110_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + nb110_iq], mm2	    ;# iq =facel*charge[ii] 

	mov   edx, [ebp + nb110_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [esp + nb110_ntype]
	shl   edx, 1
	mov   [esp + nb110_ntia], edx

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb110_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb110_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + nb110_ix], mm0	
	movd  [esp + nb110_iz], mm1	
				
	;# clear total potential and i forces 
	pxor  mm7,mm7
	movq  [esp + nb110_vctot],  mm7
	movq  [esp + nb110_Vvdwtot], mm7
	movq  [esp + nb110_fix],    mm7
	movd  [esp + nb110_fiz],    mm7

	mov   eax, [ebp + nb110_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb110_pos]
	mov   edi, [ebp + nb110_faction]	
	mov   eax, [ebp + nb110_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb110_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb110_ninner]
	mov   [esp + nb110_ninner], ecx
	mov   [esp + nb110_innerk], edx    ;# number of innerloop atoms 
	add   edx, 0
	jge   .nb110_unroll_loop
	jmp   .nb110_finish_inner
.nb110_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + nb110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + nb110_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + nb110_charge]    ;# base of charge[] 
	movq mm5, [esp + nb110_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	mov ecx, [ebp + nb110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + nb110_vdwparam]		;# base of vdwparam  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + nb110_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + nb110_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + nb110_c6], mm5
	movq [esp + nb110_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + nb110_pos]

	movq  mm0, [esp + nb110_ix]
	movd  mm1, [esp + nb110_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + nb110_dx1], mm4	     ;# store dr 
	movd  [esp + nb110_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + nb110_dx2], mm6	     ;# store dr 
	movd  [esp + nb110_dz2], mm7
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
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm3, mm1		;# mm3 has vcoul for both interactions 
	movq  mm7, mm3	    ;# use mm7 for sum to make fscal  

	pfmul mm5, [esp + nb110_c12]
	pfmul mm4, [esp + nb110_c6]	
	movq mm6, mm5	;# mm6 is Vvdw12-Vvdw6  
	pfsub mm6, mm4

	pfmul mm4, [esp + nb110_six]

	pfmul mm5, [esp + nb110_twelve]
	pfsub mm7,mm4
	pfadd mm7, mm5
 	pfmul mm0, mm7    ;# mm0 is total fscal now 	

	prefetchw [esp + nb110_dx1]	;# prefetch i forces to cache 

	;# update vctot 
	pfadd mm3, [esp + nb110_vctot]      ;# add the earlier value 
	movq [esp + nb110_vctot], mm3       ;# store the sum       

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + nb110_dx1]	;# fetch dr 
	movd mm3,  [esp + nb110_dz1]

	;# update Vvdwtot 
	pfadd mm6, [esp + nb110_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb110_Vvdwtot], mm6       ;# store the sum       

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + nb110_dx2] 	;# fetch dr 
	movd mm5,  [esp + nb110_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + nb110_fix]
	movd mm1,  [esp + nb110_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + nb110_fix], mm0
	movd [esp + nb110_fiz], mm1
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
	sub dword ptr [esp + nb110_innerk],  2
	jl    .nb110_finish_inner
	jmp   .nb110_unroll_loop
.nb110_finish_inner:	
	and dword ptr [esp + nb110_innerk],  1
	jnz  .nb110_single_inner
	jmp  .nb110_updateouterdata		
.nb110_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + nb110_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + nb110_charge]
	movd mm5, [esp + nb110_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov esi, [ebp + nb110_vdwparam]
	mov ecx, [ebp + nb110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + nb110_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + nb110_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + nb110_c12], mm5


	mov   esi, [ebp + nb110_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + nb110_ix]
	movd  mm1, [esp + nb110_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + nb110_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + nb110_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	movq  mm1, mm0
	pfmul mm0, mm0		;# mm0=invsq 
	;# calculate potentials and scalar force 
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm3, mm1		;# mm3 has vcoul for both interactions 
	movq  mm7, mm3	    ;# use mm7 for sum to make fscal  

	pfmul mm5, [esp + nb110_c12]
	pfmul mm4, [esp + nb110_c6]	
	movq mm6, mm5	;# mm6 is Vvdw12-Vvdw6  
	pfsub mm6, mm4

	pfmul mm4, [esp + nb110_six]

	pfmul mm5, [esp + nb110_twelve]
	pfsub mm7,mm4
	pfadd mm7, mm5
 	pfmul mm0, mm7    ;# mm0 is total fscal now 

	;# update vctot 
	pfadd mm3, [esp + nb110_vctot]
	movq [esp + nb110_vctot], mm3

	;# update Vvdwtot 
	pfadd mm6, [esp + nb110_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb110_Vvdwtot], mm6       ;# store the sum       

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + nb110_dx1]
	movd mm3,  [esp + nb110_dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + nb110_fix]
	movd mm1,  [esp + nb110_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + nb110_fix], mm0
	movd [esp + nb110_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.nb110_updateouterdata:	
	mov   ecx, [esp + nb110_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + nb110_fix]
	pfadd mm7, [esp + nb110_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + nb110_fshift]    ;# increment fshift force 
	mov   edx, [esp + nb110_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + nb110_fix]
	pfadd mm7, [esp + nb110_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# get n from stack
	mov esi, [esp + nb110_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb110_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb110_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb110_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + nb110_Vvdwtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb110_Vvdw]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment Vvdw[gid] 

       ;# finish if last 
        mov ecx, [esp + nb110_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb110_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb110_n], esi
        jmp .nb110_outer
.nb110_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb110_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb110_end
        ;# non-zero, do one more workunit
        jmp   .nb110_threadloop
.nb110_end:
	femms

	mov eax, [esp + nb110_nouter] 	
	mov ebx, [esp + nb110_ninner]
	mov ecx, [ebp + nb110_outeriter]
	mov edx, [ebp + nb110_inneriter]
	mov [ecx], eax
	mov [edx], ebx
	
	add esp, 152
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret




.globl nb_kernel110nf_ia32_3dnow
.globl _nb_kernel110nf_ia32_3dnow
nb_kernel110nf_ia32_3dnow:	
_nb_kernel110nf_ia32_3dnow:	
.equiv		nb110nf_p_nri,		8
.equiv		nb110nf_iinr,		12
.equiv		nb110nf_jindex,		16
.equiv		nb110nf_jjnr,		20
.equiv		nb110nf_shift,		24
.equiv		nb110nf_shiftvec,	28
.equiv		nb110nf_fshift,		32
.equiv		nb110nf_gid,		36
.equiv		nb110nf_pos,		40		
.equiv		nb110nf_faction,	44
.equiv		nb110nf_charge,		48
.equiv		nb110nf_p_facel,	52
.equiv		nb110nf_p_krf,		56	
.equiv		nb110nf_p_crf,		60	
.equiv		nb110nf_Vc,		64	
.equiv		nb110nf_type,		68
.equiv		nb110nf_p_ntype,	72
.equiv		nb110nf_vdwparam,	76	
.equiv		nb110nf_Vvdw,		80	
.equiv		nb110nf_p_tabscale,	84	
.equiv		nb110nf_VFtab,		88
.equiv		nb110nf_invsqrta,	92	
.equiv		nb110nf_dvda,		96
.equiv          nb110nf_p_gbtabscale,   100
.equiv          nb110nf_GBtab,          104
.equiv          nb110nf_p_nthreads,     108
.equiv          nb110nf_count,          112
.equiv          nb110nf_mtx,            116
.equiv          nb110nf_outeriter,      120
.equiv          nb110nf_inneriter,      124
.equiv          nb110nf_work,           128
	;# stack offsets for local variables 
.equiv		nb110nf_is3,		0
.equiv		nb110nf_ii3,		4
.equiv		nb110nf_ix,		8
.equiv		nb110nf_iy,		12
.equiv		nb110nf_iz,		16
.equiv		nb110nf_iq,		20 
.equiv		nb110nf_vctot,		28 
.equiv		nb110nf_Vvdwtot,	36 
.equiv		nb110nf_c6,		44 
.equiv		nb110nf_c12,		52
.equiv		nb110nf_ntia,		60
.equiv		nb110nf_innerjjnr,	64
.equiv		nb110nf_innerk,		68					
.equiv          nb110nf_n,              72 ;# idx for outer loop
.equiv          nb110nf_nn1,            76 ;# number of outer iterations
.equiv          nb110nf_nri,              80
.equiv          nb110nf_facel,            84
.equiv          nb110nf_ntype,            88
.equiv          nb110nf_nouter,           92
.equiv          nb110nf_ninner,           96
	push ebp
	mov ebp,esp
	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 100		;# local stack space 
	femms
	mov ecx, [ebp + nb110nf_p_nri]
	mov edx, [ebp + nb110nf_p_ntype]
	mov esi, [ebp + nb110nf_p_facel]
	mov ecx, [ecx]
	mov edx, [edx]
	mov esi, [esi]
	mov [esp + nb110nf_nri], ecx
	mov [esp + nb110nf_ntype], edx
	mov [esp + nb110nf_facel], esi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb110nf_nouter], eax
	mov [esp + nb110nf_ninner], eax

.nb110nf_threadloop:
        mov   esi, [ebp + nb110nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb110nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [esi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb110nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb110nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb110nf_n], eax
        mov [esp + nb110nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb110nf_outerstart
        jmp .nb110nf_end

.nb110nf_outerstart:	
	;# ebx contains number of outer iterations
	add ebx, [esp + nb110nf_nouter]
        mov [esp + nb110nf_nouter], ebx
	
.nb110nf_outer:
	mov   eax, [ebp + nb110nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb110nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb110nf_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + nb110nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx=ii 

	mov   edx, [ebp + nb110nf_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [esp + nb110nf_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + nb110nf_iq], mm2	    ;# iq =facel*charge[ii] 

	mov   edx, [ebp + nb110nf_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [esp + nb110nf_ntype]
	shl   edx, 1
	mov   [esp + nb110nf_ntia], edx

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb110nf_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb110nf_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + nb110nf_ix], mm0	
	movd  [esp + nb110nf_iz], mm1	
				
	;# clear total potential and i forces 
	pxor  mm7,mm7
	movq  [esp + nb110nf_vctot],  mm7
	movq  [esp + nb110nf_Vvdwtot], mm7

	mov   eax, [ebp + nb110nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb110nf_pos]
	mov   eax, [ebp + nb110nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb110nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb110nf_ninner]
	mov   [esp + nb110nf_ninner], ecx
	mov   [esp + nb110nf_innerk], edx    ;# number of innerloop atoms 
	add   edx, 0
	jge   .nb110nf_unroll_loop
	jmp   .nb110nf_finish_inner
.nb110nf_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + nb110nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + nb110nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + nb110nf_charge]    ;# base of charge[] 
	movq mm5, [esp + nb110nf_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
        punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	mov ecx, [ebp + nb110nf_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + nb110nf_vdwparam]		;# base of vdwparam  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + nb110nf_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + nb110nf_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + nb110nf_c6], mm5
	movq [esp + nb110nf_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + nb110nf_pos]

	movq  mm0, [esp + nb110nf_ix]
	movd  mm1, [esp + nb110nf_iz]	 	
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
	movq mm1,mm0
	pfmul mm0,mm0
	;# mm0 now contains invsq, and mm1 invsqrt 
	;# do potential and fscal 
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm3, mm1		;# mm3 has vcoul for both interactions 
	pfmul mm5, [esp + nb110nf_c12]
	pfmul mm4, [esp + nb110nf_c6]	
	movq mm6, mm5	;# mm6 is Vvdw12-Vvdw6  
	pfsub mm6, mm4
	;# update vctot 
	pfadd mm3, [esp + nb110nf_vctot]      ;# add the earlier value 
	movq [esp + nb110nf_vctot], mm3       ;# store the sum       
	;# update Vvdwtot 
	pfadd mm6, [esp + nb110nf_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb110nf_Vvdwtot], mm6       ;# store the sum       
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb110nf_innerk],  2
	jl    .nb110nf_finish_inner
	jmp   .nb110nf_unroll_loop
.nb110nf_finish_inner:	
	and dword ptr [esp + nb110nf_innerk],  1
	jnz  .nb110nf_single_inner
	jmp  .nb110nf_updateouterdata		
.nb110nf_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + nb110nf_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + nb110nf_charge]
	movd mm5, [esp + nb110nf_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov esi, [ebp + nb110nf_vdwparam]
	mov ecx, [ebp + nb110nf_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + nb110nf_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + nb110nf_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + nb110nf_c12], mm5


	mov   esi, [ebp + nb110nf_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + nb110nf_ix]
	movd  mm1, [esp + nb110nf_iz]
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
	movq  mm1, mm0
	pfmul mm0, mm0		;# mm0=invsq 
	;# calculate potentials and scalar force 
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm3, mm1		;# mm3 has vcoul for both interactions 
	pfmul mm5, [esp + nb110nf_c12]
	pfmul mm4, [esp + nb110nf_c6]	
	movq mm6, mm5	;# mm6 is Vvdw12-Vvdw6  
	pfsub mm6, mm4
	;# update vctot 
	pfadd mm3, [esp + nb110nf_vctot]
	movq [esp + nb110nf_vctot], mm3
	;# update Vvdwtot 
	pfadd mm6, [esp + nb110nf_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb110nf_Vvdwtot], mm6       ;# store the sum       

.nb110nf_updateouterdata:	
	;# get n from stack
	mov esi, [esp + nb110nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb110nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb110nf_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb110nf_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + nb110nf_Vvdwtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb110nf_Vvdw]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment Vvdw[gid] 

       	;# finish if last 
        mov ecx, [esp + nb110nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb110nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb110nf_n], esi
        jmp .nb110nf_outer
.nb110nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb110nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb110nf_end
        ;# non-zero, do one more workunit
        jmp   .nb110nf_threadloop
.nb110nf_end:
	femms

	mov eax, [esp + nb110nf_nouter] 	
	mov ebx, [esp + nb110nf_ninner]
	mov ecx, [ebp + nb110nf_outeriter]
	mov edx, [ebp + nb110nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	add esp, 100
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret
