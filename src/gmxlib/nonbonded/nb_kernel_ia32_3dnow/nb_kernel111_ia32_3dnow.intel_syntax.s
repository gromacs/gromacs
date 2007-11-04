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





.globl nb_kernel111_ia32_3dnow
.globl _nb_kernel111_ia32_3dnow
nb_kernel111_ia32_3dnow:	
_nb_kernel111_ia32_3dnow:	
.equiv		nb111_p_nri,		8
.equiv		nb111_iinr,		12
.equiv		nb111_jindex,		16
.equiv		nb111_jjnr,		20
.equiv		nb111_shift,		24
.equiv		nb111_shiftvec,		28
.equiv		nb111_fshift,		32
.equiv		nb111_gid,		36
.equiv		nb111_pos,		40		
.equiv		nb111_faction,		44
.equiv		nb111_charge,		48
.equiv		nb111_p_facel,		52
.equiv		nb111_p_krf,		56	
.equiv		nb111_p_crf,		60	
.equiv		nb111_Vc,		64	
.equiv		nb111_type,		68
.equiv		nb111_p_ntype,		72
.equiv		nb111_vdwparam,		76	
.equiv		nb111_Vvdw,		80	
.equiv		nb111_p_tabscale,	84	
.equiv		nb111_VFtab,		88
.equiv		nb111_invsqrta,		92	
.equiv		nb111_dvda,		96
.equiv          nb111_p_gbtabscale,     100
.equiv          nb111_GBtab,            104
.equiv          nb111_p_nthreads,       108
.equiv          nb111_count,            112
.equiv          nb111_mtx,              116
.equiv          nb111_outeriter,        120
.equiv          nb111_inneriter,        124
.equiv          nb111_work,             128
			;# stack offsets for local variables 
.equiv		nb111_is3,		0
.equiv		nb111_ii3,		4
.equiv		nb111_ixO,		8
.equiv		nb111_iyO,		12
.equiv		nb111_izO,		16	
.equiv		nb111_ixH,		20  
.equiv		nb111_iyH,		28  
.equiv		nb111_izH,		36  
.equiv		nb111_iqO,		44  
.equiv		nb111_iqH,		52  
.equiv		nb111_vctot,		60  
.equiv		nb111_Vvdwtot,		68  
.equiv		nb111_c6,		76  
.equiv		nb111_c12,		84  
.equiv		nb111_six,		92  
.equiv		nb111_twelve,		100 
.equiv		nb111_ntia,		108 
.equiv		nb111_innerjjnr,	116
.equiv		nb111_innerk,		120	
.equiv		nb111_fixO,		124
.equiv		nb111_fiyO,		128
.equiv		nb111_fizO,		132
.equiv		nb111_fixH,		136  
.equiv		nb111_fiyH,		144  
.equiv		nb111_fizH,		152  
.equiv		nb111_dxO,		160
.equiv		nb111_dyO,		164
.equiv		nb111_dzO,		168
.equiv		nb111_dxH,		172  
.equiv		nb111_dyH,		180  
.equiv		nb111_dzH,		188  
.equiv          nb111_n,                196 ;# idx for outer loop
.equiv          nb111_nn1,              200 ;# number of outer iterations
.equiv          nb111_nri,              204
.equiv          nb111_ntype,            208
.equiv          nb111_nouter,           212
.equiv          nb111_ninner,           216
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
	mov ecx, [ebp + nb111_p_nri]
	mov edx, [ebp + nb111_p_ntype]
	mov esi, [ebp + nb111_p_facel]
	mov ecx, [ecx]
	mov edx, [edx]
	mov [esp + nb111_nri], ecx
	mov [esp + nb111_ntype], edx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb111_nouter], eax
	mov [esp + nb111_ninner], eax

	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + nb111_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + nb111_charge]
	movd  mm1, [esi]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] 
	pfmul mm2, mm1		
	movq  [esp + nb111_iqO], mm2	    ;# iqO = facel*charge[ii] 
	
	movd  mm2, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] 
	pfmul mm2, mm1
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + nb111_iqH], mm2	    ;# iqH = facel*charge[ii0+1] 

	mov   edx, [ebp + nb111_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	imul  ecx, [esp + nb111_ntype]      ;# ecx = ntia = 2*ntype*type[ii0]  
	mov   [esp + nb111_ntia], ecx
	
	;# move data to local stack  
	mov eax, 0x40c00000 ;# fp 6.0
	mov ebx, 0x41400000 ;# fp 12.0
	
        mov [esp + nb111_six], eax
        mov [esp + nb111_six+4], eax
        mov [esp + nb111_twelve], ebx
        mov [esp + nb111_twelve+4], ebx
.nb111_threadloop:
        mov   esi, [ebp + nb111_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb111_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb111_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb111_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb111_n], eax
        mov [esp + nb111_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb111_outerstart
        jmp .nb111_end

.nb111_outerstart:	
	;# ebx contains number of outer iterations
	add ebx, [esp + nb111_nouter]
        mov [esp + nb111_nouter], ebx
	
.nb111_outer:
	mov   eax, [ebp + nb111_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb111_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb111_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + nb111_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb111_pos]    ;# eax = base of pos[] 
	
	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb111_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + nb111_ixO], mm5	
	movq  [esp + nb111_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + nb111_ixH], mm0	
	movq [esp + nb111_iyH], mm1	
	movq [esp + nb111_izH], mm2	
					
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + nb111_vctot], mm7
	movq  [esp + nb111_Vvdwtot], mm7
	movq  [esp + nb111_fixO],   mm7
	movd  [esp + nb111_fizO],   mm7
	movq  [esp + nb111_fixH],   mm7
	movq  [esp + nb111_fiyH],   mm7
	movq  [esp + nb111_fizH],   mm7

	mov   eax, [ebp + nb111_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + nb111_innerk], edx    ;# number of innerloop atoms 
	add   edx, [esp + nb111_ninner]
	mov   [esp + nb111_ninner], edx

	mov   esi, [ebp + nb111_pos]
	mov   edi, [ebp + nb111_faction]	
	mov   eax, [ebp + nb111_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb111_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.nb111_inner_loop:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + nb111_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    	add dword ptr [esp + nb111_innerjjnr],  4 ;# advance pointer 
	prefetch [ecx + 16]	   ;# prefetch data - trial and error says 16 is best 

	mov ecx, [ebp + nb111_charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + nb111_iqO]
	pfmul mm7, [esp + nb111_iqH]	;# mm6=qqO, mm7=qqH 

	mov ecx, [ebp + nb111_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr] 
	mov ecx, [ebp + nb111_vdwparam]
	shl edx, 1
	add edx, [esp + nb111_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [ecx + edx*4]		;# mm5 = 1st c6  		
	movq [esp + nb111_c6], mm5
	movd mm5, [ecx + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + nb111_c12], mm5	
	
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
	
	pfsubr mm0, [esp + nb111_ixO]
	pfsubr mm1, [esp + nb111_izO]
		
	movq  [esp + nb111_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + nb111_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb111_ixH]
	pfsubr mm3, [esp + nb111_iyH]
	pfsubr mm4, [esp + nb111_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + nb111_dxH], mm2
	movq [esp + nb111_dyH], mm3
	movq [esp + nb111_dzH], mm4
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

	movq  mm0, mm4
	pfmul mm0, mm4
	pfmul mm0, mm4		;# mm0=rinvsix 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm2=rintwelve 
	
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	movq  mm1, mm6		;# use mm1 for fscal sum 

	;# LJ for the oxygen 
	pfmul mm0, [esp + nb111_c6]	 
	pfmul mm2, [esp + nb111_c12]	 

	;# calc nb potential 
	movq mm5, mm2
	pfsub mm5, mm0

	;# calc nb force 
	pfmul mm0, [esp + nb111_six]
	pfmul mm2, [esp + nb111_twelve]
	
	;# increment scalar force 
	pfsub mm1, mm0
	pfadd mm1, mm2
	pfmul mm4, mm1		;# total scalar force on oxygen. 
	
	;# update nb potential 
	pfadd mm5, [esp + nb111_Vvdwtot]
	movq [esp + nb111_Vvdwtot], mm5
	
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
	pfadd mm7, [esp + nb111_vctot]
	movq [esp + nb111_vctot], mm7
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm4,mm4
	;# calc vectorial force for O 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm0,  [esp + nb111_dxO]
	movd mm1,  [esp + nb111_dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;# calc vectorial force for H's 
	movq mm5, [esp + nb111_dxH]
	movq mm6, [esp + nb111_dyH]
	movq mm7, [esp + nb111_dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;# update iO particle force 
	movq mm2,  [esp + nb111_fixO]
	movd mm3,  [esp + nb111_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + nb111_fixO], mm2
	movd [esp + nb111_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + nb111_fixH]
	movq mm3, [esp + nb111_fiyH]
	movq mm4, [esp + nb111_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + nb111_fixH], mm2
	movq [esp + nb111_fiyH], mm3
	movq [esp + nb111_fizH], mm4
	
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
	
	;#  done  - one more? 
	dec dword ptr [esp + nb111_innerk]
	jz  .nb111_updateouterdata
	jmp .nb111_inner_loop
.nb111_updateouterdata:	
	mov   ecx, [esp + nb111_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment iO force  
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + nb111_fixO]
	pfadd mm7, [esp + nb111_fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + nb111_fixH]
	movq  mm3, [esp + nb111_fiyH]
	movq  mm1, [esp + nb111_fizH]
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

	
	mov   ebx, [ebp + nb111_fshift]    ;# increment fshift force 
	mov   edx, [esp + nb111_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + nb111_fixO]
	pfadd mm7, [esp + nb111_fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	;# get n from stack
	mov esi, [esp + nb111_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb111_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb111_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb111_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + nb111_Vvdwtot]     
	pfacc mm7,mm7	          ;# same for Vvdw 
	
	mov   eax, [ebp + nb111_Vvdw]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment Vvdw[gid] 
       ;# finish if last 
        mov ecx, [esp + nb111_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb111_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb111_n], esi
        jmp .nb111_outer
.nb111_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb111_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb111_end
        ;# non-zero, do one more workunit
        jmp   .nb111_threadloop
.nb111_end:
	femms

	mov eax, [esp + nb111_nouter] 	
	mov ebx, [esp + nb111_ninner]
	mov ecx, [ebp + nb111_outeriter]
	mov edx, [ebp + nb111_inneriter]
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




.globl nb_kernel111nf_ia32_3dnow
.globl _nb_kernel111nf_ia32_3dnow
nb_kernel111nf_ia32_3dnow:	
_nb_kernel111nf_ia32_3dnow:	
.equiv		nb111nf_p_nri,		8
.equiv		nb111nf_iinr,		12
.equiv		nb111nf_jindex,		16
.equiv		nb111nf_jjnr,		20
.equiv		nb111nf_shift,		24
.equiv		nb111nf_shiftvec,	28
.equiv		nb111nf_fshift,		32
.equiv		nb111nf_gid,		36
.equiv		nb111nf_pos,		40		
.equiv		nb111nf_faction,	44
.equiv		nb111nf_charge,		48
.equiv		nb111nf_p_facel,	52
.equiv		nb111nf_p_krf,		56	
.equiv		nb111nf_p_crf,		60	
.equiv		nb111nf_Vc,		64	
.equiv		nb111nf_type,		68
.equiv		nb111nf_p_ntype,	72
.equiv		nb111nf_vdwparam,	76	
.equiv		nb111nf_Vvdw,		80	
.equiv		nb111nf_p_tabscale,	84	
.equiv		nb111nf_VFtab,		88
.equiv		nb111nf_invsqrta,	92	
.equiv		nb111nf_dvda,		96
.equiv          nb111nf_p_gbtabscale,   100
.equiv          nb111nf_GBtab,          104
.equiv          nb111nf_p_nthreads,     108
.equiv          nb111nf_count,          112
.equiv          nb111nf_mtx,            116
.equiv          nb111nf_outeriter,      120
.equiv          nb111nf_inneriter,      124
.equiv          nb111nf_work,           128
			;# stack offsets for local variables 
.equiv		nb111nf_is3,		0
.equiv		nb111nf_ii3,		4
.equiv		nb111nf_ixO,		8
.equiv		nb111nf_iyO,		12
.equiv		nb111nf_izO,		16	
.equiv		nb111nf_ixH,		20  
.equiv		nb111nf_iyH,		28  
.equiv		nb111nf_izH,		36  
.equiv		nb111nf_iqO,		44  
.equiv		nb111nf_iqH,		52  
.equiv		nb111nf_vctot,		60  
.equiv		nb111nf_Vvdwtot,	68  
.equiv		nb111nf_c6,		76  
.equiv		nb111nf_c12,		84  
.equiv		nb111nf_ntia,		92
.equiv		nb111nf_innerjjnr,	96
.equiv		nb111nf_innerk,		100
.equiv          nb111nf_n,              104 ;# idx for outer loop
.equiv          nb111nf_nn1,            108 ;# number of outer iterations
.equiv          nb111nf_nri,            112
.equiv          nb111nf_ntype,          116
.equiv          nb111nf_nouter,         120
.equiv          nb111nf_ninner,         124
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 128		;# local stack space 
	femms
	;# zero iteration counters
	mov ecx, [ebp + nb111nf_p_nri]
	mov edx, [ebp + nb111nf_p_ntype]
	mov esi, [ebp + nb111nf_p_facel]
	mov ecx, [ecx]
	mov edx, [edx]
	mov [esp + nb111nf_nri], ecx
	mov [esp + nb111nf_ntype], edx

	mov eax, 0
	mov [esp + nb111nf_nouter], eax
	mov [esp + nb111nf_ninner], eax

	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + nb111nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + nb111nf_charge]
	movd  mm1, [esi]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] 
	pfmul mm2, mm1		
	movq  [esp + nb111nf_iqO], mm2	    ;# iqO = facel*charge[ii] 
	
	movd  mm2, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] 
	pfmul mm2, mm1
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + nb111nf_iqH], mm2	    ;# iqH = facel*charge[ii0+1] 

	mov   edx, [ebp + nb111nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	imul  ecx, [esp + nb111nf_ntype]      ;# ecx = ntia = 2*ntype*type[ii0]  
	mov   [esp + nb111nf_ntia], ecx
	
.nb111nf_threadloop:
        mov   esi, [ebp + nb111nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb111nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb111nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb111nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb111nf_n], eax
        mov [esp + nb111nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb111nf_outerstart
        jmp .nb111nf_end

.nb111nf_outerstart:	
	;# ebx contains number of outer iterations
	add ebx, [esp + nb111nf_nouter]
        mov [esp + nb111nf_nouter], ebx
	
.nb111nf_outer:
	mov   eax, [ebp + nb111nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb111nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb111nf_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + nb111nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb111nf_pos]    ;# eax = base of pos[] 
	
	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb111nf_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + nb111nf_ixO], mm5	
	movq  [esp + nb111nf_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + nb111nf_ixH], mm0	
	movq [esp + nb111nf_iyH], mm1	
	movq [esp + nb111nf_izH], mm2	
					
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + nb111nf_vctot], mm7
	movq  [esp + nb111nf_Vvdwtot], mm7

	mov   eax, [ebp + nb111nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + nb111nf_innerk], edx    ;# number of innerloop atoms 
	add   edx, [esp + nb111nf_ninner]
	mov   [esp + nb111nf_ninner], edx

	mov   esi, [ebp + nb111nf_pos]
	mov   eax, [ebp + nb111nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb111nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.nb111nf_inner_loop:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + nb111nf_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    	add dword ptr [esp + nb111nf_innerjjnr],  4 ;# advance pointer 
	prefetch [ecx + 16]	   ;# prefetch data - trial and error says 16 is best 

	mov ecx, [ebp + nb111nf_charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + nb111nf_iqO]
	pfmul mm7, [esp + nb111nf_iqH]	;# mm6=qqO, mm7=qqH 

	mov ecx, [ebp + nb111nf_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr] 
	mov ecx, [ebp + nb111nf_vdwparam]
	shl edx, 1
	add edx, [esp + nb111nf_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [ecx + edx*4]		;# mm5 = 1st c6  		
	movq [esp + nb111nf_c6], mm5
	movd mm5, [ecx + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + nb111nf_c12], mm5	
	
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
	
	pfsubr mm0, [esp + nb111nf_ixO]
	pfsubr mm1, [esp + nb111nf_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb111nf_ixH]
	pfsubr mm3, [esp + nb111nf_iyH]
	pfsubr mm4, [esp + nb111nf_izH] ;# mm2-mm4 is dxH-dzH 
	
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

	movq  mm0, mm4
	pfmul mm0, mm4
	pfmul mm0, mm4		;# mm0=rinvsix 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm2=rintwelve 
	
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	movq  mm1, mm6		;# use mm1 for fscal sum 

	;# LJ for the oxygen 
	pfmul mm0, [esp + nb111nf_c6]	 
	pfmul mm2, [esp + nb111nf_c12]	 

	;# calc nb potential 
	pfsub mm2, mm0
	;# update nb potential 
	pfadd mm2, [esp + nb111nf_Vvdwtot]
	movq [esp + nb111nf_Vvdwtot], mm2
	
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
	pfadd mm7, [esp + nb111nf_vctot]
	movq [esp + nb111nf_vctot], mm7
		
	;#  done  - one more? 
	dec dword ptr [esp + nb111nf_innerk]
	jz  .nb111nf_updateouterdata
	jmp .nb111nf_inner_loop
.nb111nf_updateouterdata:	
	;# get n from stack
	mov esi, [esp + nb111nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb111nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb111nf_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb111nf_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + nb111nf_Vvdwtot]     
	pfacc mm7,mm7	          ;# same for Vvdw 
	
	mov   eax, [ebp + nb111nf_Vvdw]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment Vvdw[gid] 
      	;# finish if last 
        mov ecx, [esp + nb111nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb111nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb111nf_n], esi
        jmp .nb111nf_outer
.nb111nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb111nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb111nf_end
        ;# non-zero, do one more workunit
        jmp   .nb111nf_threadloop
.nb111nf_end:
	femms

	mov eax, [esp + nb111nf_nouter] 	
	mov ebx, [esp + nb111nf_ninner]
	mov ecx, [ebp + nb111nf_outeriter]
	mov edx, [ebp + nb111nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx
	
	add esp, 128
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret

	
