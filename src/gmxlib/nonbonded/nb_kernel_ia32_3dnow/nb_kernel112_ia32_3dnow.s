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



	

.globl nb_kernel112_ia32_3dnow
.globl _nb_kernel112_ia32_3dnow
nb_kernel112_ia32_3dnow:	
_nb_kernel112_ia32_3dnow:	
.equiv		nb112_p_nri,		8
.equiv		nb112_iinr,		12
.equiv		nb112_jindex,		16
.equiv		nb112_jjnr,		20
.equiv		nb112_shift,		24
.equiv		nb112_shiftvec,		28
.equiv		nb112_fshift,		32
.equiv		nb112_gid,		36
.equiv		nb112_pos,		40		
.equiv		nb112_faction,		44
.equiv		nb112_charge,		48
.equiv		nb112_p_facel,		52
.equiv		nb112_p_krf,		56	
.equiv		nb112_p_crf,		60	
.equiv		nb112_Vc,		64	
.equiv		nb112_type,		68
.equiv		nb112_p_ntype,		72
.equiv		nb112_vdwparam,		76	
.equiv		nb112_Vvdw,		80	
.equiv		nb112_p_tabscale,	84	
.equiv		nb112_VFtab,		88
.equiv		nb112_invsqrta,		92	
.equiv		nb112_dvda,		96
.equiv          nb112_p_gbtabscale,     100
.equiv          nb112_GBtab,            104
.equiv          nb112_p_nthreads,       108
.equiv          nb112_count,            112
.equiv          nb112_mtx,              116
.equiv          nb112_outeriter,        120
.equiv          nb112_inneriter,        124
.equiv          nb112_work,             128
			;# stack offsets for local variables 
.equiv		nb112_is3,		0
.equiv		nb112_ii3,		4
.equiv		nb112_ixO,		8
.equiv		nb112_iyO,		12
.equiv		nb112_izO,		16	
.equiv		nb112_ixH,		20  
.equiv		nb112_iyH,		28  
.equiv		nb112_izH,		36  
.equiv		nb112_qqOO,		44  
.equiv		nb112_qqOH,		52  
.equiv		nb112_qqHH,		60  
.equiv		nb112_c6,		68  
.equiv		nb112_c12,		76  
.equiv		nb112_six,		84  
.equiv		nb112_twelve,		92  
.equiv		nb112_vctot,		100 
.equiv		nb112_Vvdwtot,		108 
.equiv		nb112_innerjjnr,	116
.equiv		nb112_innerk,		120	
.equiv		nb112_fixO,		124
.equiv		nb112_fiyO,		128
.equiv		nb112_fizO,		132
.equiv		nb112_fixH,		136 
.equiv		nb112_fiyH,		144 
.equiv		nb112_fizH,		152 
.equiv		nb112_dxO,		160
.equiv		nb112_dyO,		164
.equiv		nb112_dzO,		168
.equiv		nb112_dxH,		172 
.equiv		nb112_dyH,		180 
.equiv		nb112_dzH,		188 
.equiv          nb112_n,                196 ;# idx for outer loop
.equiv          nb112_nn1,              200 ;# number of outer iterations
.equiv          nb112_nri,              204
.equiv          nb112_ntype,            208
.equiv          nb112_nouter,           212
.equiv          nb112_ninner,           216
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 220		;# local stack space 
	femms

	mov ecx, [ebp + nb112_p_nri]
	mov edx, [ebp + nb112_p_ntype]
	mov esi, [ebp + nb112_p_facel]
	mov ecx, [ecx]
	mov edx, [edx]
	mov [esp + nb112_nri], ecx
	mov [esp + nb112_ntype], edx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb112_nouter], eax
	mov [esp + nb112_ninner], eax

	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + nb112_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	            ;# ebx=ii 

	mov   edx, [ebp + nb112_charge]
	movd  mm1, [esi]		    ;# mm1=facel 
	movd  mm2, [edx + ebx*4]            ;# mm2=charge[ii0] (O) 
	movd  mm3, [edx + ebx*4 + 4]        ;# mm2=charge[ii0+1] (H)  
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
	movq  [esp + nb112_qqOO], mm4
	movq  [esp + nb112_qqOH], mm5
	movq  [esp + nb112_qqHH], mm6
	mov   edx, [ebp + nb112_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [esp + nb112_ntype]
	add   edx, ecx
	mov   eax, [ebp + nb112_vdwparam]
	movd  mm0, [eax + edx*4]          
	movd  mm1, [eax + edx*4 + 4]
	movq  [esp + nb112_c6], mm0
	movq  [esp + nb112_c12], mm1
	;# move data to local stack  
	mov eax, 0x40c00000 ;# fp 6.0
	mov ebx, 0x41400000 ;# fp 12.0
		
        mov [esp + nb112_six], eax
        mov [esp + nb112_six + 4], eax
        mov [esp + nb112_twelve], ebx
        mov [esp + nb112_twelve +4], ebx

.nb112_threadloop:
        mov   esi, [ebp + nb112_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb112_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [esi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb112_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb112_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb112_n], eax
        mov [esp + nb112_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb112_outerstart
        jmp .nb112_end

.nb112_outerstart:	
	;# ebx contains number of outer iterations
	add ebx, [esp + nb112_nouter]
        mov [esp + nb112_nouter], ebx
	
.nb112_outer:
	mov   eax, [ebp + nb112_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb112_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb112_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + nb112_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb112_pos]    ;# eax = base of pos[] 

	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb112_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + nb112_ixO], mm5	
	movq  [esp + nb112_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + nb112_ixH], mm0	
	movq [esp + nb112_iyH], mm1	
	movq [esp + nb112_izH], mm2	

	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + nb112_vctot], mm7
	movq  [esp + nb112_Vvdwtot], mm7
	movq  [esp + nb112_fixO],  mm7
	movq  [esp + nb112_fizO],  mm7
	movq  [esp + nb112_fixH],  mm7
	movq  [esp + nb112_fiyH],  mm7
	movq  [esp + nb112_fizH],  mm7

	mov   eax, [ebp + nb112_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + nb112_innerk], edx    ;# number of innerloop atoms 
	add   edx, [esp + nb112_ninner]
	mov   [esp + nb112_ninner], edx

	mov   esi, [ebp + nb112_pos]
	mov   edi, [ebp + nb112_faction]	
	mov   eax, [ebp + nb112_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb112_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.nb112_inner_loop:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + nb112_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    	add dword ptr [esp + nb112_innerjjnr],  4 ;# advance pointer 

	movd  mm6, [esp + nb112_qqOO]
	movq  mm7, [esp + nb112_qqOH]

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
	
	pfsubr mm0, [esp + nb112_ixO]
	pfsubr mm1, [esp + nb112_izO]
		
	movq  [esp + nb112_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + nb112_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb112_ixH]
	pfsubr mm3, [esp + nb112_iyH]
	pfsubr mm4, [esp + nb112_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + nb112_dxH], mm2
	movq [esp + nb112_dyH], mm3
	movq [esp + nb112_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

    	pfrsqrt mm1,mm0

    	movq mm2,mm1
    	pfmul mm1,mm1
    	pfrsqit1 mm1,mm0				
    	pfrcpit2 mm1,mm2	;# mm1=invsqrt  
	movq  mm4, mm1
	pfmul mm4, mm4		;# mm4=invsq  

	movq mm2, mm4
	pfmul mm2, mm4
	pfmul mm2, mm4
	movq mm0, mm2
	pfmul mm0,mm0
	pfmul mm2, [esp + nb112_c6]
	pfmul mm0, [esp + nb112_c12]
	movq mm5, mm0
	pfsub mm5, mm2		;# Vvdw 

	pfmul mm2, [esp + nb112_six]
	pfmul mm0, [esp + nb112_twelve]

	pfsub mm0, mm2
	
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	pfadd mm0, mm6
	pfmul mm4, mm0		;# mm4=fscalar  

	;# update nb potential 
	pfadd mm5, [esp + nb112_Vvdwtot]
	movq [esp + nb112_Vvdwtot], mm5

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3 

	movq mm2, mm5
	pfmul mm5,mm5
    	pfrsqit1 mm5,mm3				
    	pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	movq mm3,mm5
	pfmul mm3,mm3		;# mm3=invsq 
	pfmul mm7, mm5		;# mm7=vcoul 
	pfmul mm3, mm7		;# mm3=fscal for the two H's. 

	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + nb112_vctot]
	movq [esp + nb112_vctot], mm7
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm4,mm4
	;# calc vectorial force for O 
	movq mm0,  [esp + nb112_dxO]
	movd mm1,  [esp + nb112_dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;# calc vectorial force for H's 
	movq mm5, [esp + nb112_dxH]
	movq mm6, [esp + nb112_dyH]
	movq mm7, [esp + nb112_dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;# update iO particle force 
	movq mm2,  [esp + nb112_fixO]
	movd mm3,  [esp + nb112_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + nb112_fixO], mm2
	movd [esp + nb112_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + nb112_fixH]
	movq mm3, [esp + nb112_fiyH]
	movq mm4, [esp + nb112_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + nb112_fixH], mm2
	movq [esp + nb112_fiyH], mm3
	movq [esp + nb112_fizH], mm4
	
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
	
	movd mm6, [esp + nb112_qqOH]
	movq mm7, [esp + nb112_qqHH]
	
	pfsubr mm0, [esp + nb112_ixO]
	pfsubr mm1, [esp + nb112_izO]
		
	movq  [esp + nb112_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + nb112_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb112_ixH]
	pfsubr mm3, [esp + nb112_iyH]
	pfsubr mm4, [esp + nb112_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + nb112_dxH], mm2
	movq [esp + nb112_dyH], mm3
	movq [esp + nb112_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

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

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3 

	movq mm2, mm5
	pfmul mm5,mm5
    	pfrsqit1 mm5,mm3				
    	pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	movq mm3,mm5
	pfmul mm3,mm3		;# mm3=invsq 
	pfmul mm7, mm5		;# mm7=vcoul 
	pfmul mm3, mm7		;# mm3=fscal for the two H's. 

	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + nb112_vctot]
	movq [esp + nb112_vctot], mm7
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm4,mm4
	;# calc vectorial force for O 
	movq mm0,  [esp + nb112_dxO]
	movd mm1,  [esp + nb112_dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;# calc vectorial force for H's 
	movq mm5, [esp + nb112_dxH]
	movq mm6, [esp + nb112_dyH]
	movq mm7, [esp + nb112_dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;# update iO particle force 
	movq mm2,  [esp + nb112_fixO]
	movd mm3,  [esp + nb112_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + nb112_fixO], mm2
	movd [esp + nb112_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + nb112_fixH]
	movq mm3, [esp + nb112_fiyH]
	movq mm4, [esp + nb112_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + nb112_fixH], mm2
	movq [esp + nb112_fiyH], mm3
	movq [esp + nb112_fizH], mm4
	
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

	movd mm6, [esp + nb112_qqOH]
	movq mm7, [esp + nb112_qqHH]

	pfsubr mm0, [esp + nb112_ixO]
	pfsubr mm1, [esp + nb112_izO]
		
	movq  [esp + nb112_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + nb112_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb112_ixH]
	pfsubr mm3, [esp + nb112_iyH]
	pfsubr mm4, [esp + nb112_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + nb112_dxH], mm2
	movq [esp + nb112_dyH], mm3
	movq [esp + nb112_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

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

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3. 

	movq mm2, mm5
	pfmul mm5,mm5
    	pfrsqit1 mm5,mm3				
    	pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	movq mm3,mm5
	pfmul mm3,mm3		;# mm3=invsq 
	pfmul mm7, mm5		;# mm7=vcoul 
	pfmul mm3, mm7		;# mm3=fscal for the two H's. 

	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + nb112_vctot]
	movq [esp + nb112_vctot], mm7
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm4,mm4
	;# calc vectorial force for O 
	movq mm0,  [esp + nb112_dxO]
	movd mm1,  [esp + nb112_dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;# calc vectorial force for H's 
	movq mm5, [esp + nb112_dxH]
	movq mm6, [esp + nb112_dyH]
	movq mm7, [esp + nb112_dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;# update iO particle force 
	movq mm2,  [esp + nb112_fixO]
	movd mm3,  [esp + nb112_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + nb112_fixO], mm2
	movd [esp + nb112_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + nb112_fixH]
	movq mm3, [esp + nb112_fiyH]
	movq mm4, [esp + nb112_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + nb112_fixH], mm2
	movq [esp + nb112_fiyH], mm3
	movq [esp + nb112_fizH], mm4	

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
	dec dword ptr [esp + nb112_innerk]
	jz  .nb112_updateouterdata
	jmp .nb112_inner_loop	
.nb112_updateouterdata:	
	mov   ecx, [esp + nb112_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment iO force  
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + nb112_fixO]
	pfadd mm7, [esp + nb112_fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + nb112_fixH]
	movq  mm3, [esp + nb112_fiyH]
	movq  mm1, [esp + nb112_fizH]
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

	
	mov   ebx, [ebp + nb112_fshift]    ;# increment fshift force 
	mov   edx, [esp + nb112_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + nb112_fixO]
	pfadd mm7, [esp + nb112_fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	;# get n from stack
	mov esi, [esp + nb112_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb112_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb112_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + nb112_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + nb112_Vvdwtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + nb112_Vvdw]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment Vvdwtot[gid] 
       ;# finish if last 
        mov ecx, [esp + nb112_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb112_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb112_n], esi
        jmp .nb112_outer
.nb112_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb112_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb112_end
        ;# non-zero, do one more workunit
        jmp   .nb112_threadloop
.nb112_end:
	femms

	mov eax, [esp + nb112_nouter] 	
	mov ebx, [esp + nb112_ninner]
	mov ecx, [ebp + nb112_outeriter]
	mov edx, [ebp + nb112_inneriter]
	mov [ecx], eax
	mov [edx], ebx
	
	add esp, 220
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


.globl nb_kernel112nf_ia32_3dnow
.globl _nb_kernel112nf_ia32_3dnow
nb_kernel112nf_ia32_3dnow:	
_nb_kernel112nf_ia32_3dnow:	
.equiv		nb112nf_p_nri,		8
.equiv		nb112nf_iinr,		12
.equiv		nb112nf_jindex,		16
.equiv		nb112nf_jjnr,		20
.equiv		nb112nf_shift,		24
.equiv		nb112nf_shiftvec,	28
.equiv		nb112nf_fshift,		32
.equiv		nb112nf_gid,		36
.equiv		nb112nf_pos,		40		
.equiv		nb112nf_faction,	44
.equiv		nb112nf_charge,		48
.equiv		nb112nf_p_facel,		52
.equiv		nb112nf_p_krf,		56	
.equiv		nb112nf_p_crf,		60	
.equiv		nb112nf_Vc,		64	
.equiv		nb112nf_type,		68
.equiv		nb112nf_p_ntype,	72
.equiv		nb112nf_vdwparam,	76	
.equiv		nb112nf_Vvdw,		80	
.equiv		nb112nf_p_tabscale,	84	
.equiv		nb112nf_VFtab,		88
.equiv		nb112nf_invsqrta,	92	
.equiv		nb112nf_dvda,		96
.equiv          nb112nf_p_gbtabscale,   100
.equiv          nb112nf_GBtab,          104
.equiv          nb112nf_p_nthreads,     108
.equiv          nb112nf_count,          112
.equiv          nb112nf_mtx,            116
.equiv          nb112nf_outeriter,      120
.equiv          nb112nf_inneriter,      124
.equiv          nb112nf_work,           128
			;# stack offsets for local variables 
.equiv		nb112nf_is3,		0
.equiv		nb112nf_ii3,		4
.equiv		nb112nf_ixO,		8
.equiv		nb112nf_iyO,		12
.equiv		nb112nf_izO,		16	
.equiv		nb112nf_ixH,		20  
.equiv		nb112nf_iyH,		28  
.equiv		nb112nf_izH,		36  
.equiv		nb112nf_qqOO,		44  
.equiv		nb112nf_qqOH,		52  
.equiv		nb112nf_qqHH,		60  
.equiv		nb112nf_c6,		68  
.equiv		nb112nf_c12,		76 
.equiv		nb112nf_vctot,		84 
.equiv		nb112nf_Vvdwtot,	92 
.equiv		nb112nf_innerjjnr,	100
.equiv		nb112nf_innerk,		104
.equiv          nb112nf_n,              108 ;# idx for outer loop
.equiv          nb112nf_nn1,            112 ;# number of outer iterations
.equiv          nb112nf_nri,            116
.equiv          nb112nf_ntype,          120
.equiv          nb112nf_nouter,         124
.equiv          nb112nf_ninner,         128
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
	mov ecx, [ebp + nb112nf_p_nri]
	mov edx, [ebp + nb112nf_p_ntype]
	mov esi, [ebp + nb112nf_p_facel]
	mov ecx, [ecx]
	mov edx, [edx]
	mov [esp + nb112nf_nri], ecx
	mov [esp + nb112nf_ntype], edx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb112nf_nouter], eax
	mov [esp + nb112nf_ninner], eax

	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + nb112nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + nb112nf_charge]
	movd  mm1, [esi]		;# mm1=facel 
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
	movq  [esp + nb112nf_qqOO], mm4
	movq  [esp + nb112nf_qqOH], mm5
	movq  [esp + nb112nf_qqHH], mm6
	mov   edx, [ebp + nb112nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [esp + nb112nf_ntype]
	add   edx, ecx
	mov   eax, [ebp + nb112nf_vdwparam]
	movd  mm0, [eax + edx*4]          
	movd  mm1, [eax + edx*4 + 4]
	movq  [esp + nb112nf_c6], mm0
	movq  [esp + nb112nf_c12], mm1
	
.nb112nf_threadloop:
        mov   esi, [ebp + nb112nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb112nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [esi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb112nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb112nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb112nf_n], eax
        mov [esp + nb112nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb112nf_outerstart
        jmp .nb112nf_end

.nb112nf_outerstart:	
	;# ebx contains number of outer iterations
	add ebx, [esp + nb112nf_nouter]
        mov [esp + nb112nf_nouter], ebx
	
.nb112nf_outer:
	mov   eax, [ebp + nb112nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb112nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb112nf_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + nb112nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb112nf_pos]    ;# eax = base of pos[] 

	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb112nf_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + nb112nf_ixO], mm5	
	movq  [esp + nb112nf_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + nb112nf_ixH], mm0	
	movq [esp + nb112nf_iyH], mm1	
	movq [esp + nb112nf_izH], mm2	

	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + nb112nf_vctot], mm7
	movq  [esp + nb112nf_Vvdwtot], mm7

	mov   eax, [ebp + nb112nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + nb112nf_innerk], edx    ;# number of innerloop atoms 
	add   edx, [esp + nb112nf_ninner]
	mov   [esp + nb112nf_ninner], edx

	mov   esi, [ebp + nb112nf_pos]
	mov   eax, [ebp + nb112nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb112nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.nb112nf_inner_loop:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + nb112nf_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    	add dword ptr [esp + nb112nf_innerjjnr],  4 ;# advance pointer 

	movd  mm6, [esp + nb112nf_qqOO]
	movq  mm7, [esp + nb112nf_qqOH]

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
	
	pfsubr mm0, [esp + nb112nf_ixO]
	pfsubr mm1, [esp + nb112nf_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb112nf_ixH]
	pfsubr mm3, [esp + nb112nf_iyH]
	pfsubr mm4, [esp + nb112nf_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

    	pfrsqrt mm1,mm0

    	movq mm2,mm1
    	pfmul mm1,mm1
    	pfrsqit1 mm1,mm0				
    	pfrcpit2 mm1,mm2	;# mm1=invsqrt  
	movq  mm4, mm1
	pfmul mm4, mm4		;# mm4=invsq  

	movq mm2, mm4
	pfmul mm2, mm4
	pfmul mm2, mm4
	movq mm0, mm2
	pfmul mm0,mm0
	pfmul mm2, [esp + nb112nf_c6]
	pfmul mm0, [esp + nb112nf_c12]
	movq mm5, mm0
	pfsub mm5, mm2		;# Vvdw 

	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	;# update nb potential 
	pfadd mm5, [esp + nb112nf_Vvdwtot]
	movq [esp + nb112nf_Vvdwtot], mm5

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3 

	movq mm2, mm5
	pfmul mm5,mm5
    	pfrsqit1 mm5,mm3				
    	pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	pfmul mm7, mm5		;# mm7=vcoul 
	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + nb112nf_vctot]
	movq [esp + nb112nf_vctot], mm7
	
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
	
	movd mm6, [esp + nb112nf_qqOH]
	movq mm7, [esp + nb112nf_qqHH]
	
	pfsubr mm0, [esp + nb112nf_ixO]
	pfsubr mm1, [esp + nb112nf_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb112nf_ixH]
	pfsubr mm3, [esp + nb112nf_iyH]
	pfsubr mm4, [esp + nb112nf_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

    	pfrsqrt mm1,mm0

    	movq mm2,mm1
    	pfmul mm1,mm1
    	pfrsqit1 mm1,mm0				
    	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3 

	movq mm2, mm5
	pfmul mm5,mm5
    	pfrsqit1 mm5,mm3				
    	pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	pfmul mm7, mm5		;# mm7=vcoul 
	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + nb112nf_vctot]
	movq [esp + nb112nf_vctot], mm7
	
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

	movd mm6, [esp + nb112nf_qqOH]
	movq mm7, [esp + nb112nf_qqHH]

	pfsubr mm0, [esp + nb112nf_ixO]
	pfsubr mm1, [esp + nb112nf_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb112nf_ixH]
	pfsubr mm3, [esp + nb112nf_iyH]
	pfsubr mm4, [esp + nb112nf_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

    	pfrsqrt mm1,mm0

    	movq mm2,mm1
    	pfmul mm1,mm1
    	pfrsqit1 mm1,mm0				
    	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3. 

	movq mm2, mm5
	pfmul mm5,mm5
    	pfrsqit1 mm5,mm3				
    	pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	pfmul mm7, mm5		;# mm7=vcoul 

	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + nb112nf_vctot]
	movq [esp + nb112nf_vctot], mm7
		
	;#  done  - one more? 
	dec dword ptr [esp + nb112nf_innerk]
	jz  .nb112nf_updateouterdata
	jmp .nb112nf_inner_loop	
.nb112nf_updateouterdata:	
	;# get n from stack
	mov esi, [esp + nb112nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb112nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb112nf_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + nb112nf_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + nb112nf_Vvdwtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + nb112nf_Vvdw]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment Vvdwtot[gid] 
       	;# finish if last 
        mov ecx, [esp + nb112nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb112nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb112nf_n], esi
        jmp .nb112nf_outer
.nb112nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb112nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb112nf_end
        ;# non-zero, do one more workunit
        jmp   .nb112nf_threadloop
.nb112nf_end:
	femms

	mov eax, [esp + nb112nf_nouter] 	
	mov ebx, [esp + nb112nf_ninner]
	mov ecx, [ebp + nb112nf_outeriter]
	mov edx, [ebp + nb112nf_inneriter]
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

