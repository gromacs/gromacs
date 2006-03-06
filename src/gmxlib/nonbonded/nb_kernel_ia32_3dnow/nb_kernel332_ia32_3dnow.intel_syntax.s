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



.globl nb_kernel332_ia32_3dnow
.globl _nb_kernel332_ia32_3dnow
nb_kernel332_ia32_3dnow:	
_nb_kernel332_ia32_3dnow:	
.equiv		nb332_p_nri,		8
.equiv		nb332_iinr,		12
.equiv		nb332_jindex,		16
.equiv		nb332_jjnr,		20
.equiv		nb332_shift,		24
.equiv		nb332_shiftvec,		28
.equiv		nb332_fshift,		32
.equiv		nb332_gid,		36
.equiv		nb332_pos,		40		
.equiv		nb332_faction,		44
.equiv		nb332_charge,		48
.equiv		nb332_p_facel,		52
.equiv		nb332_p_krf,		56	
.equiv		nb332_p_crf,		60	
.equiv		nb332_Vc,		64	
.equiv		nb332_type,		68
.equiv		nb332_p_ntype,		72
.equiv		nb332_vdwparam,		76	
.equiv		nb332_Vvdw,		80	
.equiv		nb332_p_tabscale,	84	
.equiv		nb332_VFtab,		88
.equiv		nb332_invsqrta,		92	
.equiv		nb332_dvda,		96
.equiv          nb332_p_gbtabscale,     100
.equiv          nb332_GBtab,            104
.equiv          nb332_p_nthreads,       108
.equiv          nb332_count,            112
.equiv          nb332_mtx,              116
.equiv          nb332_outeriter,        120
.equiv          nb332_inneriter,        124
.equiv          nb332_work,             128
			;# stack offsets for local variables 
.equiv		nb332_is3,		0
.equiv		nb332_ii3,		4
.equiv		nb332_ixO,		8
.equiv		nb332_iyO,		12
.equiv		nb332_izO,		16	
.equiv		nb332_ixH,		20  
.equiv		nb332_iyH,		28  
.equiv		nb332_izH,		36  
.equiv		nb332_qqOO,		44  
.equiv		nb332_qqOH,		52  
.equiv		nb332_qqHH,		60  
.equiv		nb332_c6,		68  
.equiv		nb332_c12,		76  
.equiv		nb332_two,		84  
.equiv		nb332_n1,		92  
.equiv		nb332_tsc,		100 
.equiv		nb332_vctot,		108 
.equiv		nb332_Vvdwtot,		116 
.equiv		nb332_innerjjnr,	124
.equiv		nb332_innerk,		128	
.equiv		nb332_fixO,		132
.equiv		nb332_fiyO,		136
.equiv		nb332_fizO,		140
.equiv		nb332_fixH,		144 
.equiv		nb332_fiyH,		152 
.equiv		nb332_fizH,		160 
.equiv		nb332_dxO,		168
.equiv		nb332_dyO,		172
.equiv		nb332_dzO,		176
.equiv		nb332_dxH,		180  
.equiv		nb332_dyH,		188  
.equiv		nb332_dzH,		196  
.equiv		nb332_tmprsqH,		204  
.equiv          nb332_n,                212 ;# idx for outer loop
.equiv          nb332_nn1,              216 ;# number of outer iterations
.equiv          nb332_nri,              220
.equiv          nb332_nouter,           224
.equiv          nb332_ninner,           228
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 232		;# local stack space 
	femms
	mov ecx, [ebp + nb332_p_nri]
	mov esi, [ebp + nb332_p_facel]
	mov edi, [ebp + nb332_p_tabscale]
	mov ecx, [ecx]
	mov [esp + nb332_nri], ecx
	movd  mm1, [esi] 	;# facel
	mov esi, [ebp + nb332_p_ntype]

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb332_nouter], eax
	mov [esp + nb332_ninner], eax

	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + nb332_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + nb332_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] (O) 
	movd  mm3, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] (H)  
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			;# mm4=qqOO*facel 
	pfmul mm5, mm3			;# mm5=qqOH*facel 
	pfmul mm6, mm3			;# mm6=qqHH*facel 
	punpckldq mm5,mm5	    ;# spread to both halves 
	punpckldq mm6,mm6	    ;# spread to both halves 
	movq  [esp + nb332_qqOO], mm4
	movq  [esp + nb332_qqOH], mm5
	movq  [esp + nb332_qqHH], mm6
	mov   edx, [ebp + nb332_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [esi]
	add   edx, ecx
	mov   eax, [ebp + nb332_vdwparam]
	movd  mm0, [eax + edx*4]
	movd  mm1, [eax + edx*4 + 4]
	movq  [esp + nb332_c6], mm0
	movq  [esp + nb332_c12], mm1
	movd  mm3, [edi]
	punpckldq mm3,mm3
	movq  [esp + nb332_tsc], mm3
	mov eax, 0x40000000
	mov [esp + nb332_two], eax
	mov [esp + nb332_two+4], eax
.nb332_threadloop:
        mov   esi, [ebp + nb332_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb332_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb332_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb332_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb332_n], eax
        mov [esp + nb332_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb332_outerstart
        jmp .nb332_end

.nb332_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb332_nouter]
        mov [esp + nb332_nouter], ebx
	
.nb332_outer:
	mov   eax, [ebp + nb332_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb332_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb332_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + nb332_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb332_pos]    ;# eax = base of pos[] 

	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb332_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + nb332_ixO], mm5	
	movq  [esp + nb332_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + nb332_ixH], mm0	
	movq [esp + nb332_iyH], mm1	
	movq [esp + nb332_izH], mm2	

	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + nb332_vctot], mm7
	movq  [esp + nb332_Vvdwtot], mm7
	movq  [esp + nb332_fixO],  mm7
	movq  [esp + nb332_fizO],  mm7
	movq  [esp + nb332_fixH],  mm7
	movq  [esp + nb332_fiyH],  mm7
	movq  [esp + nb332_fizH],  mm7

	mov   eax, [ebp + nb332_jindex]
	mov   ecx, [eax+esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + nb332_innerk], edx        
	add   edx, [esp + nb332_ninner]
	mov   [esp + nb332_ninner], edx

	mov   esi, [ebp + nb332_pos]
	mov   edi, [ebp + nb332_faction]	
	mov   eax, [ebp + nb332_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb332_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.nb332_inner_loop:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + nb332_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + nb332_innerjjnr],  4 ;# advance pointer 

	lea   eax, [eax + eax*2]

	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + nb332_ixO]
	pfsubr mm1, [esp + nb332_izO]
		
	movq  [esp + nb332_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + nb332_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb332_ixH]
	pfsubr mm3, [esp + nb332_iyH]
	pfsubr mm4, [esp + nb332_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + nb332_dxH], mm2
	movq [esp + nb332_dyH], mm3
	movq [esp + nb332_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + nb332_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt  
	pfmul mm0, mm1		;# mm0=rsq  

	pfmul mm0, [esp + nb332_tsc]
	pf2iw mm4, mm0
	movd [esp + nb332_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + nb332_VFtab]
	mov ecx, [esp + nb332_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + nb332_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + nb332_qqOO]	;# vcoul=qq*VV 
	pfmul mm7, [esp + nb332_qqOO]	;# fijC=qq*FF 

	;# update vctot directly, use mm3 for fscal sum. 
	pfadd mm5, [esp + nb332_vctot]
	movq [esp + nb332_vctot], mm5
	movq mm3, mm7

	;# dispersion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + nb332_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + nb332_c6]
	pfmul mm7, mm4	;# fijD 
	pfmul mm5, mm4	;# Vvdw6            
	pfadd mm3, mm7	;# add to fscal  

	;# update Vvdwtot to release mm5! 
	pfadd mm5, [esp + nb332_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb332_Vvdwtot], mm5       ;# store the sum       

	;# repulsion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + nb332_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + nb332_c12]
	pfmul mm7, mm6	;# fijR 
	pfmul mm5, mm6	;# Vvdw12 
	pfadd mm3, mm7	;# total fscal fijC+ fijD+ fijR 

	;# change sign of fscal and multiply with rinv  
    pxor mm0,mm0
	pfsubr mm3, mm0	
	pfmul mm3, [esp + nb332_tsc]
 	pfmul mm3, mm1    ;# mm3 is total fscal (for the oxygen) now 
	
	;# update Vvdwtot  
	pfadd mm5, [esp + nb332_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb332_Vvdwtot], mm5       ;# store the sum       
	
	;# Ready with the oxygen - potential is updated, fscal is in mm3.
	 ;# time for hydrogens!
         
	
	movq mm0, [esp + nb332_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + nb332_tsc]
	pf2iw mm4, mm0
	movq [esp + nb332_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + nb332_VFtab]
	mov ecx, [esp + nb332_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + nb332_n1+4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + nb332_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + nb332_qqOH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + nb332_qqOH]	;# fijC=qq*FF 
	;# update vctot 
	pfadd mm5, [esp + nb332_vctot]
	movq [esp + nb332_vctot], mm5
	
	;# change sign of fijC and multiply by rinv 
    pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + nb332_tsc]
 	pfmul mm4, mm1    ;# mm4 is total fscal (for the hydrogens) now 	
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm3,mm3
	;# calc vectorial force for O 
	movq mm0,  [esp + nb332_dxO]
	movd mm1,  [esp + nb332_dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;# calc vectorial force for H's 
	movq mm5, [esp + nb332_dxH]
	movq mm6, [esp + nb332_dyH]
	movq mm7, [esp + nb332_dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;# update iO particle force 
	movq mm2,  [esp + nb332_fixO]
	movd mm3,  [esp + nb332_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + nb332_fixO], mm2
	movd [esp + nb332_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + nb332_fixH]
	movq mm3, [esp + nb332_fiyH]
	movq mm4, [esp + nb332_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + nb332_fixH], mm2
	movq [esp + nb332_fiyH], mm3
	movq [esp + nb332_fizH], mm4
	
	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3

	;# interactions with j H1 

	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + nb332_ixO]
	pfsubr mm1, [esp + nb332_izO]
		
	movq  [esp + nb332_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + nb332_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb332_ixH]
	pfsubr mm3, [esp + nb332_iyH]
	pfsubr mm4, [esp + nb332_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + nb332_dxH], mm2
	movq [esp + nb332_dyH], mm3
	movq [esp + nb332_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + nb332_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	pfmul mm0, mm1		;# mm0=rsq  
	
	pfmul mm0, [esp + nb332_tsc]
	pf2iw mm4, mm0
	movd [esp + nb332_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + nb332_VFtab]
	mov ecx, [esp + nb332_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + nb332_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + nb332_qqOH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + nb332_qqOH]	;# fijC=qq*FF 

	;# update vctot directly, force is moved to mm3. 
	pfadd mm5, [esp + nb332_vctot]
	movq [esp + nb332_vctot], mm5
	pxor mm3, mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + nb332_tsc]
	pfmul mm3, mm1    ;# mm3 is total fscal (for the oxygen) now 

	movq mm0, [esp + nb332_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + nb332_tsc]
	pf2iw mm4, mm0
	movq [esp + nb332_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + nb332_VFtab]
	mov ecx, [esp + nb332_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + nb332_n1+4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + nb332_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + nb332_qqHH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + nb332_qqHH]	;# fijC=qq*FF 
	;# update vctot 
	pfadd mm5, [esp + nb332_vctot]
	movq [esp + nb332_vctot], mm5
	
	;# change sign of fijC and multiply by rinv 
    pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + nb332_tsc]
 	pfmul mm4, mm1    ;# mm4 is total fscal (for the hydrogens) now 		

	;# spread oxygen fscalar to both positions 
	punpckldq mm3,mm3
	;# calc vectorial force for O 
	movq mm0,  [esp + nb332_dxO]
	movd mm1,  [esp + nb332_dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;# calc vectorial force for H's 
	movq mm5, [esp + nb332_dxH]
	movq mm6, [esp + nb332_dyH]
	movq mm7, [esp + nb332_dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;# update iO particle force 
	movq mm2,  [esp + nb332_fixO]
	movd mm3,  [esp + nb332_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + nb332_fixO], mm2
	movd [esp + nb332_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + nb332_fixH]
	movq mm3, [esp + nb332_fiyH]
	movq mm4, [esp + nb332_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + nb332_fixH], mm2
	movq [esp + nb332_fiyH], mm3
	movq [esp + nb332_fizH], mm4
	
	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4 + 12]
	movd mm3,  [edi + eax*4 + 20]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 12], mm2
	movd [edi + eax*4 + 20], mm3

	;# interactions with j H2 
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	pfsubr mm0, [esp + nb332_ixO]
	pfsubr mm1, [esp + nb332_izO]
		
	movq  [esp + nb332_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + nb332_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb332_ixH]
	pfsubr mm3, [esp + nb332_iyH]
	pfsubr mm4, [esp + nb332_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + nb332_dxH], mm2
	movq [esp + nb332_dyH], mm3
	movq [esp + nb332_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + nb332_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	pfmul mm0, mm1

	pfmul mm0, [esp + nb332_tsc]
	pf2iw mm4, mm0
	movd [esp + nb332_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + nb332_VFtab]
	mov ecx, [esp + nb332_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + nb332_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + nb332_qqOH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + nb332_qqOH]	;# fijC=qq*FF 

	;# update vctot directly, use mm3 for fscal sum 
	pfadd mm5, [esp + nb332_vctot]
	movq [esp + nb332_vctot], mm5
	pxor mm3,mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + nb332_tsc]
 	pfmul mm3, mm1     ;# mm3 is total fscal (for the oxygen) now 	

	movq mm0, [esp + nb332_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + nb332_tsc]
	pf2iw mm4, mm0
	movq [esp + nb332_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + nb332_VFtab]
	mov ecx, [esp + nb332_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + nb332_n1+4];# mm5 = Fp 
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + nb332_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + nb332_qqHH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + nb332_qqHH]	;# fijC=qq*FF 
	;# update vctot 
	pfadd mm5, [esp + nb332_vctot]
	movq [esp + nb332_vctot], mm5
	
	;# change sign of fijC and multiply by rinv 
    pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + nb332_tsc]
 	pfmul mm4, mm1    ;# mm4 is total fscal (for the hydrogens) now 	

	;# spread oxygen fscalar to both positions 
	punpckldq mm3,mm3
	;# calc vectorial force for O 
	movq mm0,  [esp + nb332_dxO]
	movd mm1,  [esp + nb332_dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;# calc vectorial force for H's 
	movq mm5, [esp + nb332_dxH]
	movq mm6, [esp + nb332_dyH]
	movq mm7, [esp + nb332_dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;# update iO particle force 
	movq mm2,  [esp + nb332_fixO]
	movd mm3,  [esp + nb332_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + nb332_fixO], mm2
	movd [esp + nb332_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + nb332_fixH]
	movq mm3, [esp + nb332_fiyH]
	movq mm4, [esp + nb332_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + nb332_fixH], mm2
	movq [esp + nb332_fiyH], mm3
	movq [esp + nb332_fizH], mm4	

	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4 + 24]
	movd mm3,  [edi + eax*4 + 32]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 24], mm2
	movd [edi + eax*4 + 32], mm3
	
	;#  done  - one more? 
	dec dword ptr [esp + nb332_innerk]
	jz  .nb332_updateouterdata
	jmp .nb332_inner_loop	
.nb332_updateouterdata:	
	mov   ecx, [esp + nb332_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment iO force  
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + nb332_fixO]
	pfadd mm7, [esp + nb332_fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + nb332_fixH]
	movq  mm3, [esp + nb332_fiyH]
	movq  mm1, [esp + nb332_fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;# mm0(l)=fxH1, mm0(h)=fyH1 
	punpckhdq mm2, mm3	;# mm2(l)=fxH2, mm2(h)=fyH2 
	movq mm3, mm1
	pswapd mm3,mm3		
	;# mm1 is fzH1 
	;# mm3 is fzH2 

	movq  mm6, [edi + ecx*4 + 12]       ;# increment iH1 force  
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;# increment iH2 force 
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + nb332_fshift]    ;# increment fshift force 
	mov   edx, [esp + nb332_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + nb332_fixO]
	pfadd mm7, [esp + nb332_fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	;# get n from stack
	mov esi, [esp + nb332_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb332_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb332_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + nb332_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + nb332_Vvdwtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + nb332_Vvdw]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment Vvdwtot[gid] 
       	;# finish if last 
        mov ecx, [esp + nb332_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb332_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb332_n], esi
        jmp .nb332_outer
.nb332_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb332_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb332_end
        ;# non-zero, do one more workunit
        jmp   .nb332_threadloop
.nb332_end:
	femms
	mov eax, [esp + nb332_nouter] 	
	mov ebx, [esp + nb332_ninner]
	mov ecx, [ebp + nb332_outeriter]
	mov edx, [ebp + nb332_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	add esp, 232
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret



.globl nb_kernel332nf_ia32_3dnow
.globl _nb_kernel332nf_ia32_3dnow
nb_kernel332nf_ia32_3dnow:	
_nb_kernel332nf_ia32_3dnow:	
.equiv		nb332nf_p_nri,		8
.equiv		nb332nf_iinr,		12
.equiv		nb332nf_jindex,		16
.equiv		nb332nf_jjnr,		20
.equiv		nb332nf_shift,		24
.equiv		nb332nf_shiftvec,	28
.equiv		nb332nf_fshift,		32
.equiv		nb332nf_gid,		36
.equiv		nb332nf_pos,		40		
.equiv		nb332nf_faction,	44
.equiv		nb332nf_charge,		48
.equiv		nb332nf_p_facel,	52
.equiv		nb332nf_p_krf,		56	
.equiv		nb332nf_p_crf,		60	
.equiv		nb332nf_Vc,		64	
.equiv		nb332nf_type,		68
.equiv		nb332nf_p_ntype,	72
.equiv		nb332nf_vdwparam,	76	
.equiv		nb332nf_Vvdw,		80	
.equiv		nb332nf_p_tabscale,	84	
.equiv		nb332nf_VFtab,		88
.equiv		nb332nf_invsqrta,	92	
.equiv		nb332nf_dvda,		96
.equiv          nb332nf_p_gbtabscale,   100
.equiv          nb332nf_GBtab,          104
.equiv          nb332nf_p_nthreads,     108
.equiv          nb332nf_count,          112
.equiv          nb332nf_mtx,            116
.equiv          nb332nf_outeriter,      120
.equiv          nb332nf_inneriter,      124
.equiv          nb332nf_work,           128
			;# stack offsets for local variables 
.equiv		nb332nf_is3,		0
.equiv		nb332nf_ii3,		4
.equiv		nb332nf_ixO,		8
.equiv		nb332nf_iyO,		12
.equiv		nb332nf_izO,		16	
.equiv		nb332nf_ixH,		20  
.equiv		nb332nf_iyH,		28  
.equiv		nb332nf_izH,		36  
.equiv		nb332nf_qqOO,		44  
.equiv		nb332nf_qqOH,		52  
.equiv		nb332nf_qqHH,		60  
.equiv		nb332nf_c6,		68  
.equiv		nb332nf_c12,		76 
.equiv		nb332nf_n1,		84  
.equiv		nb332nf_tsc,		92 
.equiv		nb332nf_vctot,		100 
.equiv		nb332nf_Vvdwtot,	108 
.equiv		nb332nf_innerjjnr,	116
.equiv		nb332nf_innerk,		120
.equiv		nb332nf_tmprsqH,	124  
.equiv          nb332nf_n,              132 ;# idx for outer loop
.equiv          nb332nf_nn1,            136 ;# number of outer iterations
.equiv          nb332nf_nri,            140
.equiv          nb332nf_nouter,         144
.equiv          nb332nf_ninner,         148
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
	mov ecx, [ebp + nb332nf_p_nri]
	mov esi, [ebp + nb332nf_p_facel]
	mov edi, [ebp + nb332nf_p_tabscale]
	mov ecx, [ecx]
	mov [esp + nb332nf_nri], ecx
	movd  mm1, [esi] 	;# facel
	mov esi, [ebp + nb332nf_p_ntype]

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb332nf_nouter], eax
	mov [esp + nb332nf_ninner], eax

	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + nb332nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + nb332nf_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] (O) 
	movd  mm3, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] (H)  
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			;# mm4=qqOO*facel 
	pfmul mm5, mm3			;# mm5=qqOH*facel 
	pfmul mm6, mm3			;# mm6=qqHH*facel 
	punpckldq mm5,mm5	    ;# spread to both halves 
	punpckldq mm6,mm6	    ;# spread to both halves 
	movq  [esp + nb332nf_qqOO], mm4
	movq  [esp + nb332nf_qqOH], mm5
	movq  [esp + nb332nf_qqHH], mm6
	mov   edx, [ebp + nb332nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [esi]
	add   edx, ecx
	mov   eax, [ebp + nb332nf_vdwparam]
	movd  mm0, [eax + edx*4]
	movd  mm1, [eax + edx*4 + 4]
	movq  [esp + nb332nf_c6], mm0
	movq  [esp + nb332nf_c12], mm1
	movd  mm3, [edi]
	punpckldq mm3,mm3
	movq  [esp + nb332nf_tsc], mm3
.nb332nf_threadloop:
        mov   esi, [ebp + nb332nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb332nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb332nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb332nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb332nf_n], eax
        mov [esp + nb332nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb332nf_outerstart
        jmp .nb332nf_end

.nb332nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb332nf_nouter]
        mov [esp + nb332nf_nouter], ebx
	
.nb332nf_outer:
	mov   eax, [ebp + nb332nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb332nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb332nf_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + nb332nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb332nf_pos]    ;# eax = base of pos[] 

	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb332nf_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + nb332nf_ixO], mm5	
	movq  [esp + nb332nf_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + nb332nf_ixH], mm0	
	movq [esp + nb332nf_iyH], mm1	
	movq [esp + nb332nf_izH], mm2	

	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + nb332nf_vctot], mm7
	movq  [esp + nb332nf_Vvdwtot], mm7

	mov   eax, [ebp + nb332nf_jindex]
	mov   ecx, [eax+esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + nb332nf_innerk], edx        
	add   edx, [esp + nb332nf_ninner]
	mov   [esp + nb332nf_ninner], edx

	mov   esi, [ebp + nb332nf_pos]	
	mov   eax, [ebp + nb332nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb332nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.nb332nf_inner_loop:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + nb332nf_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + nb332nf_innerjjnr],  4 ;# advance pointer 

	lea   eax, [eax + eax*2]

	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + nb332nf_ixO]
	pfsubr mm1, [esp + nb332nf_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb332nf_ixH]
	pfsubr mm3, [esp + nb332nf_iyH]
	pfsubr mm4, [esp + nb332nf_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + nb332nf_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt  
	pfmul mm0, mm1		;# mm0=rsq  

	pfmul mm0, [esp + nb332nf_tsc]
	pf2iw mm4, mm0
	movd [esp + nb332nf_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + nb332nf_VFtab]
	mov ecx, [esp + nb332nf_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + nb332nf_qqOO]	;# vcoul=qq*VV 
	;# update vctot directly, use mm3 for fscal sum. 
	pfadd mm5, [esp + nb332nf_vctot]
	movq [esp + nb332nf_vctot], mm5
	
	;# dispersion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + nb332nf_c6]
	pfmul mm5, mm4	;# Vvdw6            
	;# update Vvdwtot to release mm5! 
	pfadd mm5, [esp + nb332nf_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb332nf_Vvdwtot], mm5       ;# store the sum       

	;# repulsion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + nb332nf_c12]
	pfmul mm5, mm6	;# Vvdw12 
	;# change sign of fscal and multiply with rinv  
 	;# update Vvdwtot  
	pfadd mm5, [esp + nb332nf_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb332nf_Vvdwtot], mm5       ;# store the sum       
	
	;# Ready with the oxygen - time for hydrogens 
	
	movq mm0, [esp + nb332nf_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + nb332nf_tsc]
	pf2iw mm4, mm0
	movq [esp + nb332nf_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + nb332nf_VFtab]
	mov ecx, [esp + nb332nf_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + nb332nf_n1+4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + nb332nf_qqOH]	;# vcoul=qq*VV 
	;# update vctot 
	pfadd mm5, [esp + nb332nf_vctot]
	movq [esp + nb332nf_vctot], mm5

	;# interactions with j H1 
	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + nb332nf_ixO]
	pfsubr mm1, [esp + nb332nf_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb332nf_ixH]
	pfsubr mm3, [esp + nb332nf_iyH]
	pfsubr mm4, [esp + nb332nf_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + nb332nf_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	pfmul mm0, mm1		;# mm0=rsq  
	
	pfmul mm0, [esp + nb332nf_tsc]
	pf2iw mm4, mm0
	movd [esp + nb332nf_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + nb332nf_VFtab]
	mov ecx, [esp + nb332nf_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + nb332nf_qqOH]	;# vcoul=qq*VV 
	;# update vctot directly 
	pfadd mm5, [esp + nb332nf_vctot]
	movq [esp + nb332nf_vctot], mm5
	
	movq mm0, [esp + nb332nf_tmprsqH]
	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + nb332nf_tsc]
	pf2iw mm4, mm0
	movq [esp + nb332nf_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + nb332nf_VFtab]
	mov ecx, [esp + nb332nf_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + nb332nf_n1+4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + nb332nf_qqHH]	;# vcoul=qq*VV 
	;# update vctot 
	pfadd mm5, [esp + nb332nf_vctot]
	movq [esp + nb332nf_vctot], mm5
	
	;# interactions with j H2 
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	pfsubr mm0, [esp + nb332nf_ixO]
	pfsubr mm1, [esp + nb332nf_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb332nf_ixH]
	pfsubr mm3, [esp + nb332nf_iyH]
	pfsubr mm4, [esp + nb332nf_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + nb332nf_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	pfmul mm0, mm1

	pfmul mm0, [esp + nb332nf_tsc]
	pf2iw mm4, mm0
	movd [esp + nb332nf_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + nb332nf_VFtab]
	mov ecx, [esp + nb332nf_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + nb332nf_qqOH]	;# vcoul=qq*VV 
	;# update vctot directly 
	pfadd mm5, [esp + nb332nf_vctot]
	movq [esp + nb332nf_vctot], mm5
	
	movq mm0, [esp + nb332nf_tmprsqH]
	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + nb332nf_tsc]
	pf2iw mm4, mm0
	movq [esp + nb332nf_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + nb332nf_VFtab]
	mov ecx, [esp + nb332nf_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + nb332nf_n1+4] ;# mm5 = Fp 
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + nb332nf_qqHH]	;# vcoul=qq*VV 
	;# update vctot 
	pfadd mm5, [esp + nb332nf_vctot]
	movq [esp + nb332nf_vctot], mm5
		
	;#  done  - one more? 
	dec dword ptr [esp + nb332nf_innerk]
	jz  .nb332nf_updateouterdata
	jmp .nb332nf_inner_loop	
.nb332nf_updateouterdata:	
	;# get n from stack
	mov esi, [esp + nb332nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb332nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb332nf_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + nb332nf_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + nb332nf_Vvdwtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + nb332nf_Vvdw]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment Vvdwtot[gid] 
       	;# finish if last 
        mov ecx, [esp + nb332nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb332nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb332nf_n], esi
        jmp .nb332nf_outer
.nb332nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb332nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb332nf_end
        ;# non-zero, do one more workunit
        jmp   .nb332nf_threadloop
.nb332nf_end:
	femms
	mov eax, [esp + nb332nf_nouter] 	
	mov ebx, [esp + nb332nf_ninner]
	mov ecx, [ebp + nb332nf_outeriter]
	mov edx, [ebp + nb332nf_inneriter]
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
 

