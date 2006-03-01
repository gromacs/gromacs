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


			
.globl nb_kernel101_ia32_3dnow
.globl _nb_kernel101_ia32_3dnow
nb_kernel101_ia32_3dnow:	
_nb_kernel101_ia32_3dnow:	
.equiv		nb101_p_nri,		8
.equiv		nb101_iinr,		12
.equiv		nb101_jindex,		16
.equiv		nb101_jjnr,		20
.equiv		nb101_shift,		24
.equiv		nb101_shiftvec,		28
.equiv		nb101_fshift,		32
.equiv		nb101_gid,		36
.equiv		nb101_pos,		40		
.equiv		nb101_faction,		44
.equiv		nb101_charge,		48
.equiv		nb101_p_facel,		52
.equiv		nb101_p_krf,		56	
.equiv		nb101_p_crf,		60	
.equiv		nb101_Vc,		64	
.equiv		nb101_type,		68
.equiv		nb101_p_ntype,		72
.equiv		nb101_vdwparam,		76	
.equiv		nb101_Vvdw,		80	
.equiv		nb101_p_tabscale,	84	
.equiv		nb101_VFtab,		88
.equiv		nb101_invsqrta,		92	
.equiv		nb101_dvda,		96
.equiv          nb101_p_gbtabscale,     100
.equiv          nb101_GBtab,            104
.equiv          nb101_p_nthreads,       108
.equiv          nb101_count,            112
.equiv          nb101_mtx,              116
.equiv          nb101_outeriter,        120
.equiv          nb101_inneriter,        124
.equiv          nb101_work,             128
			;# stack offsets for local variables 
.equiv		nb101_is3,		0
.equiv		nb101_ii3,		4
.equiv		nb101_ixO,		8
.equiv		nb101_iyO,		12
.equiv		nb101_izO,		16	
.equiv		nb101_ixH,		20 
.equiv		nb101_iyH,		28 
.equiv		nb101_izH,		36 
.equiv		nb101_iqO,		44	
.equiv		nb101_iqH,		52		
.equiv		nb101_vctot,		60 
.equiv		nb101_innerjjnr,	68
.equiv		nb101_innerk,		72		
.equiv		nb101_fixO,		76 
.equiv		nb101_fiyO,		80
.equiv		nb101_fizO,		84
.equiv		nb101_fixH,		88
.equiv		nb101_fiyH,		96
.equiv		nb101_fizH,		104         
.equiv		nb101_dxO,		112
.equiv		nb101_dyO,		116
.equiv		nb101_dzO,		120
.equiv		nb101_dxH,		124         
.equiv		nb101_dyH,		132         
.equiv		nb101_dzH,		140         
.equiv          nb101_n,                148 ;# idx for outer loop
.equiv          nb101_nn1,              152 ;# number of outer iterations
.equiv          nb101_nri,              156
.equiv          nb101_nouter,           160
.equiv          nb101_ninner,           164
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 168		;# local stack space 
	femms

	mov ecx, [ebp + nb101_p_nri]
	mov esi, [ebp + nb101_p_facel]
	mov ecx, [ecx]
	mov [esp + nb101_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb101_nouter], eax
	mov [esp + nb101_ninner], eax

	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + nb101_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + nb101_charge]
	movd  mm1, [esi]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] 
	pfmul mm2, mm1		
	movq  [esp + nb101_iqO], mm2	    ;# iqO = facel*charge[ii] 
	
	movd  mm2, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] 
	pfmul mm2, mm1
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + nb101_iqH], mm2	    ;# iqH = facel*charge[ii0+1] 
.nb101_threadloop:
        mov   esi, [ebp + nb101_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb101_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb101_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb101_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb101_n], eax
        mov [esp + nb101_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb101_outerstart
        jmp .nb101_end

.nb101_outerstart:	
	;# ebx contains number of outer iterations
	add ebx, [esp + nb101_nouter]
        mov [esp + nb101_nouter], ebx
	
.nb101_outer:
	mov   eax, [ebp + nb101_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb101_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb101_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + nb101_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb101_pos]    ;# eax = base of pos[] 
	
	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb101_ii3], ebx	    ;# (use mm7 as temp storage for iz) 
	pfadd mm6, mm7
	movq  [esp + nb101_ixO], mm5	
	movq  [esp + nb101_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + nb101_ixH], mm0	
	movq [esp + nb101_iyH], mm1	
	movq [esp + nb101_izH], mm2	
					
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + nb101_vctot], mm7
	movq  [esp + nb101_fixO],   mm7
	movd  [esp + nb101_fizO],   mm7
	movq  [esp + nb101_fixH],   mm7
	movq  [esp + nb101_fiyH],   mm7
	movq  [esp + nb101_fizH],   mm7

	mov   eax, [ebp + nb101_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + nb101_innerk], edx    ;# number of innerloop atoms 
	add   edx, [esp + nb101_ninner]
	mov   [esp + nb101_ninner], edx

	mov   esi, [ebp + nb101_pos]
	mov   edi, [ebp + nb101_faction]	
	mov   eax, [ebp + nb101_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb101_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.nb101_inner_loop:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + nb101_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    	add dword ptr [esp + nb101_innerjjnr],  4 ;# advance pointer 
	prefetch [ecx + 16]	   ;# prefetch data - trial and error says 16 is best 

	mov ecx, [ebp + nb101_charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + nb101_iqO]
	pfmul mm7, [esp + nb101_iqH]	;# mm6=qqO, mm7=qqH 
	
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
	
	pfsubr mm0, [esp + nb101_ixO]
	pfsubr mm1, [esp + nb101_izO]
		
	movq  [esp + nb101_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + nb101_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb101_ixH]
	pfsubr mm3, [esp + nb101_iyH]
	pfsubr mm4, [esp + nb101_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + nb101_dxH], mm2
	movq [esp + nb101_dyH], mm3
	movq [esp + nb101_dzH], mm4
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
	pfmul mm3, mm7		;# mm3=fscal for the two H's 

	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + nb101_vctot]
	movq [esp + nb101_vctot], mm7
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm4,mm4
	;# calc vectorial force for O 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm0,  [esp + nb101_dxO]
	movd mm1,  [esp + nb101_dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;# calc vectorial force for H's 
	movq mm5, [esp + nb101_dxH]
	movq mm6, [esp + nb101_dyH]
	movq mm7, [esp + nb101_dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;# update iO particle force 
	movq mm2,  [esp + nb101_fixO]
	movd mm3,  [esp + nb101_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + nb101_fixO], mm2
	movd [esp + nb101_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + nb101_fixH]
	movq mm3, [esp + nb101_fiyH]
	movq mm4, [esp + nb101_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + nb101_fixH], mm2
	movq [esp + nb101_fiyH], mm3
	movq [esp + nb101_fizH], mm4
	
	;# pack j forces from H in the same form as the oxygen force 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	
	;#  done  - one more? 
	dec dword ptr [esp + nb101_innerk]
	jz  .nb101_updateouterdata
	jmp .nb101_inner_loop
.nb101_updateouterdata:	
	mov   ecx, [esp + nb101_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment iO force  
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + nb101_fixO]
	pfadd mm7, [esp + nb101_fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + nb101_fixH]
	movq  mm3, [esp + nb101_fiyH]
	movq  mm1, [esp + nb101_fizH]
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

	
	mov   ebx, [ebp + nb101_fshift]    ;# increment fshift force 
	mov   edx, [esp + nb101_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + nb101_fixO]
	pfadd mm7, [esp + nb101_fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	;# get n from stack
	mov esi, [esp + nb101_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb101_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb101_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb101_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

       ;# finish if last 
        mov ecx, [esp + nb101_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb101_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb101_n], esi
        jmp .nb101_outer
.nb101_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb101_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb101_end
        ;# non-zero, do one more workunit
        jmp   .nb101_threadloop
.nb101_end:
	femms
	mov eax, [esp + nb101_nouter] 	
	mov ebx, [esp + nb101_ninner]
	mov ecx, [ebp + nb101_outeriter]
	mov edx, [ebp + nb101_inneriter]
	mov [ecx], eax
	mov [edx], ebx
	
	add esp, 168
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


			

			
.globl nb_kernel101nf_ia32_3dnow
.globl _nb_kernel101nf_ia32_3dnow
nb_kernel101nf_ia32_3dnow:	
_nb_kernel101nf_ia32_3dnow:	
.equiv		nb101nf_p_nri,		8
.equiv		nb101nf_iinr,		12
.equiv		nb101nf_jindex,		16
.equiv		nb101nf_jjnr,		20
.equiv		nb101nf_shift,		24
.equiv		nb101nf_shiftvec,	28
.equiv		nb101nf_fshift,		32
.equiv		nb101nf_gid,		36
.equiv		nb101nf_pos,		40		
.equiv		nb101nf_faction,	44
.equiv		nb101nf_charge,		48
.equiv		nb101nf_p_facel,	52
.equiv		nb101nf_p_krf,		56	
.equiv		nb101nf_p_crf,		60	
.equiv		nb101nf_Vc,		64	
.equiv		nb101nf_type,		68
.equiv		nb101nf_p_ntype,	72
.equiv		nb101nf_vdwparam,	76	
.equiv		nb101nf_Vvdw,		80	
.equiv		nb101nf_p_tabscale,	84	
.equiv		nb101nf_VFtab,		88
.equiv		nb101nf_invsqrta,	92	
.equiv		nb101nf_dvda,		96
.equiv          nb101nf_p_gbtabscale,   100
.equiv          nb101nf_GBtab,          104
.equiv          nb101nf_p_nthreads,     108
.equiv          nb101nf_count,          112
.equiv          nb101nf_mtx,            116
.equiv          nb101nf_outeriter,      120
.equiv          nb101nf_inneriter,      124
.equiv          nb101nf_work,           128
			;# stack offsets for local variables 
.equiv		nb101nf_is3,		0
.equiv		nb101nf_ii3,		4
.equiv		nb101nf_ixO,		8
.equiv		nb101nf_iyO,		12
.equiv		nb101nf_izO,		16	
.equiv		nb101nf_ixH,		20 
.equiv		nb101nf_iyH,		28 
.equiv		nb101nf_izH,		36 
.equiv		nb101nf_iqO,		44	
.equiv		nb101nf_iqH,		52		
.equiv		nb101nf_vctot,		60 
.equiv		nb101nf_innerjjnr,	68
.equiv		nb101nf_innerk,		72   
.equiv          nb101nf_n,              76 ;# idx for outer loop
.equiv          nb101nf_nn1,            80 ;# number of outer iterations
.equiv          nb101nf_nri,            84
.equiv          nb101nf_nouter,         88
.equiv          nb101nf_ninner,         92
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

	mov ecx, [ebp + nb101nf_p_nri]
	mov esi, [ebp + nb101nf_p_facel]
	mov ecx, [ecx]
	mov [esp + nb101nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb101nf_nouter], eax
	mov [esp + nb101nf_ninner], eax

	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + nb101nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + nb101nf_charge]
	movd  mm1, [esi]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] 
	pfmul mm2, mm1		
	movq  [esp + nb101nf_iqO], mm2	    ;# iqO = facel*charge[ii] 
	
	movd  mm2, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] 
	pfmul mm2, mm1
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + nb101nf_iqH], mm2	    ;# iqH = facel*charge[ii0+1] 
.nb101nf_threadloop:
        mov   esi, [ebp + nb101nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb101nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb101nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb101nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb101nf_n], eax
        mov [esp + nb101nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb101nf_outerstart
        jmp .nb101nf_end

.nb101nf_outerstart:	
	;# ebx contains number of outer iterations
	add ebx, [esp + nb101nf_nouter]
        mov [esp + nb101nf_nouter], ebx
	
.nb101nf_outer:
	mov   eax, [ebp + nb101nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb101nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb101nf_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + nb101nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb101nf_pos]    ;# eax = base of pos[] 
	
	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb101nf_ii3], ebx	    ;# (use mm7 as temp storage for iz) 
	pfadd mm6, mm7
	movq  [esp + nb101nf_ixO], mm5	
	movq  [esp + nb101nf_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] 
	;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + nb101nf_ixH], mm0	
	movq [esp + nb101nf_iyH], mm1	
	movq [esp + nb101nf_izH], mm2	
					
	;# clear vctot 
	pxor  mm7,mm7
	movq  [esp + nb101nf_vctot], mm7

	mov   eax, [ebp + nb101nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + nb101nf_innerk], edx    ;# number of innerloop atoms 
	add   edx, [esp + nb101nf_ninner]
	mov   [esp + nb101nf_ninner], edx

	mov   esi, [ebp + nb101nf_pos]
	mov   eax, [ebp + nb101nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb101nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.nb101nf_inner_loop:	
	;# a single j particle iteration here 
	mov   eax, [esp + nb101nf_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    	add dword ptr [esp + nb101nf_innerjjnr],  4 ;# advance pointer 
	;# prefetch data - trial and error says 16 is best 
	prefetch [ecx + 16]	   

	mov ecx, [ebp + nb101nf_charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + nb101nf_iqO]
	pfmul mm7, [esp + nb101nf_iqH]	;# mm6=qqO, mm7=qqH 
	
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
	
	pfsubr mm0, [esp + nb101nf_ixO]
	pfsubr mm1, [esp + nb101nf_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb101nf_ixH]
	pfsubr mm3, [esp + nb101nf_iyH]
	pfsubr mm4, [esp + nb101nf_izH] ;# mm2-mm4 is dxH-dzH 
	
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
	pfadd mm7, [esp + nb101nf_vctot]
	movq [esp + nb101nf_vctot], mm7
	
	;#  done  - one more? 
	dec dword ptr [esp + nb101nf_innerk]
	jz  .nb101nf_updateouterdata
	jmp .nb101nf_inner_loop
.nb101nf_updateouterdata:	
	;# get n from stack
	mov esi, [esp + nb101nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb101nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb101nf_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb101nf_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

       ;# finish if last 
        mov ecx, [esp + nb101nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb101nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb101nf_n], esi
        jmp .nb101nf_outer
.nb101nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb101nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb101nf_end
        ;# non-zero, do one more workunit
        jmp   .nb101nf_threadloop
.nb101nf_end:
	femms
	mov eax, [esp + nb101nf_nouter] 	
	mov ebx, [esp + nb101nf_ninner]
	mov ecx, [ebp + nb101nf_outeriter]
	mov edx, [ebp + nb101nf_inneriter]
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


			
