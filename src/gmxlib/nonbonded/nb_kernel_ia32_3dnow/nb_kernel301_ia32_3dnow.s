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




.globl nb_kernel301_ia32_3dnow
.globl _nb_kernel301_ia32_3dnow
nb_kernel301_ia32_3dnow:	
_nb_kernel301_ia32_3dnow:	
.equiv		nb301_p_nri,		8
.equiv		nb301_iinr,		12
.equiv		nb301_jindex,		16
.equiv		nb301_jjnr,		20
.equiv		nb301_shift,		24
.equiv		nb301_shiftvec,		28
.equiv		nb301_fshift,		32
.equiv		nb301_gid,		36
.equiv		nb301_pos,		40		
.equiv		nb301_faction,		44
.equiv		nb301_charge,		48
.equiv		nb301_p_facel,		52
.equiv		nb301_p_krf,		56	
.equiv		nb301_p_crf,		60	
.equiv		nb301_Vc,		64	
.equiv		nb301_type,		68
.equiv		nb301_p_ntype,		72
.equiv		nb301_vdwparam,		76	
.equiv		nb301_Vvdw,		80	
.equiv		nb301_p_tabscale,	84	
.equiv		nb301_VFtab,		88
.equiv		nb301_invsqrta,		92	
.equiv		nb301_dvda,		96
.equiv          nb301_p_gbtabscale,     100
.equiv          nb301_GBtab,            104
.equiv          nb301_p_nthreads,       108
.equiv          nb301_count,            112
.equiv          nb301_mtx,              116
.equiv          nb301_outeriter,        120
.equiv          nb301_inneriter,        124
.equiv          nb301_work,             128
			;# stack offsets for local variables 
.equiv		nb301_is3,		0
.equiv		nb301_ii3,		4
.equiv		nb301_ixO,		8
.equiv		nb301_iyO,		12
.equiv		nb301_izO,		16	
.equiv		nb301_ixH,		20  
.equiv		nb301_iyH,		28  
.equiv		nb301_izH,		36  
.equiv		nb301_iqO,		44  
.equiv		nb301_iqH,		52  
.equiv		nb301_qqO,		60  
.equiv		nb301_qqH,		68  
.equiv		nb301_vctot,		76  
.equiv		nb301_two,		84  
.equiv		nb301_n1,		92  
.equiv		nb301_tsc,		100 
.equiv		nb301_innerjjnr,	108
.equiv		nb301_innerk,		112	
.equiv		nb301_fixO,		116
.equiv		nb301_fiyO,		120
.equiv		nb301_fizO,		124
.equiv		nb301_fixH,		128 
.equiv		nb301_fiyH,		136 
.equiv		nb301_fizH,		144 
.equiv		nb301_dxO,		152
.equiv		nb301_dyO,		156
.equiv		nb301_dzO,		160
.equiv		nb301_dxH,		164 
.equiv		nb301_dyH,		172 
.equiv		nb301_dzH,		180 
.equiv		nb301_tmprsqH,		188 
.equiv          nb301_n,                196 ;# idx for outer loop
.equiv          nb301_nn1,              200 ;# number of outer iterations
.equiv          nb301_nri,              204
.equiv          nb301_nouter,           208
.equiv          nb301_ninner,           212
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 216		;# local stack space 
	femms

	mov ecx, [ebp + nb301_p_nri]
	mov esi, [ebp + nb301_p_facel]
	mov edi, [ebp + nb301_p_tabscale]
	mov ecx, [ecx]
	mov [esp + nb301_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb301_nouter], eax
	mov [esp + nb301_ninner], eax

	mov   ecx, [ebp + nb301_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + nb301_charge]
	movd  mm1, [esi]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] 
	pfmul mm2, mm1		
	movq  [esp + nb301_iqO], mm2	    ;# iqO = facel*charge[ii] 
	
	movd  mm2, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] 
	pfmul mm2, mm1
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + nb301_iqH], mm2	    ;# iqH = facel*charge[ii0+1] 

	movd  mm4, [edi]
	punpckldq mm4,mm4	    ;# spread to both halves 
	movq  [esp + nb301_tsc], mm4	      
	mov eax, 0x40000000
	mov [esp + nb301_two], eax
	mov [esp + nb301_two + 4], eax
.nb301_threadloop:
        mov   esi, [ebp + nb301_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb301_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [esi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb301_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb301_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb301_n], eax
        mov [esp + nb301_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb301_outerstart
        jmp .nb301_end

.nb301_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb301_nouter]
        mov [esp + nb301_nouter], ebx
	
.nb301_outer:
	mov   eax, [ebp + nb301_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb301_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb301_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + nb301_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb301_pos]    ;# eax = base of pos[] 
	
	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb301_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + nb301_ixO], mm5	
	movq  [esp + nb301_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + nb301_ixH], mm0	
	movq [esp + nb301_iyH], mm1	
	movq [esp + nb301_izH], mm2	
					
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + nb301_vctot], mm7
	movq  [esp + nb301_fixO],   mm7
	movd  [esp + nb301_fizO],   mm7
	movq  [esp + nb301_fixH],   mm7
	movq  [esp + nb301_fiyH],   mm7
	movq  [esp + nb301_fizH],   mm7

	mov   eax, [ebp + nb301_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + nb301_innerk], edx             
	add   edx, [esp + nb301_ninner]
	mov   [esp + nb301_ninner], edx

	mov   esi, [ebp + nb301_pos]
	mov   edi, [ebp + nb301_faction]	
	mov   eax, [ebp + nb301_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb301_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.nb301_inner_loop:	
	;# a single j particle iteration 
	mov   eax, [esp + nb301_innerjjnr]
	mov   eax, [eax]	 ;# eax=jnr offset 
    	add dword ptr [esp + nb301_innerjjnr],  4 ;# advance pointer 
	prefetch [ecx + 16]	   ;# prefetch data - trial and error says 16 is best 

	mov ecx, [ebp + nb301_charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + nb301_iqO]
	pfmul mm7, [esp + nb301_iqH]	 ;# mm6=qqO, mm7=qqH 
	movd [esp + nb301_qqO], mm6
	movq [esp + nb301_qqH], mm7
		
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
	
	pfsubr mm0, [esp + nb301_ixO]
	pfsubr mm1, [esp + nb301_izO]
		
	movq  [esp + nb301_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + nb301_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb301_ixH]
	pfsubr mm3, [esp + nb301_iyH]
	pfsubr mm4, [esp + nb301_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + nb301_dxH], mm2
	movq [esp + nb301_dyH], mm3
	movq [esp + nb301_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + nb301_tmprsqH], mm3
	
    	pfrsqrt mm1,mm0

    	movq mm2,mm1
    	pfmul mm1,mm1
    	pfrsqit1 mm1,mm0				
    	pfrcpit2 mm1,mm2	;# mm1=invsqrt 

	pfmul mm0, mm1		;# mm0=r 

	pfmul mm0, [esp + nb301_tsc]
	pf2iw mm4, mm0
	movd [esp + nb301_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + nb301_VFtab]
	mov ecx, [esp + nb301_n1]
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

	pfmul mm7, [esp + nb301_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + nb301_qqO]	;# vcoul=qq*VV 
	pfmul mm7, [esp + nb301_qqO]	;# fijC=qq*FF 
	;# update vctot directly, use mm3 for fscal sum. 
	pfadd mm5, [esp + nb301_vctot]
	movq [esp + nb301_vctot], mm5
	movq mm3, mm7	

	;# change sign of fscal and multiply with rinv  
    	pxor mm0,mm0
	pfsubr mm3, mm0	
	pfmul mm3, [esp + nb301_tsc]
 	pfmul mm3, mm1    ;# mm3 is total fscal (for the oxygen) now 	
	
	;# Ready with the oxygen - potential is updated, fscal is in mm3. 
	;# now do the two hydrogens. 
	 
	movq mm0, [esp + nb301_tmprsqH] ;# mm0=rsqH 

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
	pfmul mm0, [esp + nb301_tsc]
	pf2iw mm4, mm0
	movq [esp + nb301_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + nb301_VFtab]
	mov ecx, [esp + nb301_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + nb301_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + nb301_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + nb301_qqH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + nb301_qqH]	;# fijC=qq*FF 

	;# update vctot 
	pfadd mm5, [esp + nb301_vctot]
	movq [esp + nb301_vctot], mm5
	
	;# change sign of fijC and multiply by rinv 
    	pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + nb301_tsc]
 	pfmul mm4, mm1    ;# mm4 is total fscal (for the hydrogens) now 	

	;# spread oxygen fscalar to both positions 
	punpckldq mm3,mm3
	;# calc vectorial force for O 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm0,  [esp + nb301_dxO]
	movd mm1,  [esp + nb301_dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;# calc vectorial force for H's 
	movq mm5, [esp + nb301_dxH]
	movq mm6, [esp + nb301_dyH]
	movq mm7, [esp + nb301_dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;# update iO particle force 
	movq mm2,  [esp + nb301_fixO]
	movd mm3,  [esp + nb301_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + nb301_fixO], mm2
	movd [esp + nb301_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + nb301_fixH]
	movq mm3, [esp + nb301_fiyH]
	movq mm4, [esp + nb301_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + nb301_fixH], mm2
	movq [esp + nb301_fiyH], mm3
	movq [esp + nb301_fizH], mm4
	
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
	movd [edi + eax*4 + 8], mm3
	
	;#  done  - one more? 
	dec dword ptr [esp + nb301_innerk]
	jz  .nb301_updateouterdata
	jmp .nb301_inner_loop
.nb301_updateouterdata:	
	mov   ecx, [esp + nb301_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment iO force  
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + nb301_fixO]
	pfadd mm7, [esp + nb301_fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + nb301_fixH]
	movq  mm3, [esp + nb301_fiyH]
	movq  mm1, [esp + nb301_fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;# mm0(l)=fxH1, mm0(h)=fyH1 
	punpckhdq mm2, mm3	;# mm2(l)=fxH2, mm2(h)=fyH2 
	movq mm3, mm1
	pswapd mm3, mm3	
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

	
	mov   ebx, [ebp + nb301_fshift]    ;# increment fshift force 
	mov   edx, [esp + nb301_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + nb301_fixO]
	pfadd mm7, [esp + nb301_fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	;# get n from stack
	mov esi, [esp + nb301_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb301_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb301_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb301_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

       	;# finish if last 
        mov ecx, [esp + nb301_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb301_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb301_n], esi
        jmp .nb301_outer
.nb301_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb301_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb301_end
        ;# non-zero, do one more workunit
        jmp   .nb301_threadloop
.nb301_end:
	femms

	mov eax, [esp + nb301_nouter] 	
	mov ebx, [esp + nb301_ninner]
	mov ecx, [ebp + nb301_outeriter]
	mov edx, [ebp + nb301_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	add esp, 216
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret




.globl nb_kernel301nf_ia32_3dnow
.globl _nb_kernel301nf_ia32_3dnow
nb_kernel301nf_ia32_3dnow:	
_nb_kernel301nf_ia32_3dnow:	
.equiv		nb301nf_p_nri,		8
.equiv		nb301nf_iinr,		12
.equiv		nb301nf_jindex,		16
.equiv		nb301nf_jjnr,		20
.equiv		nb301nf_shift,		24
.equiv		nb301nf_shiftvec,	28
.equiv		nb301nf_fshift,		32
.equiv		nb301nf_gid,		36
.equiv		nb301nf_pos,		40		
.equiv		nb301nf_faction,	44
.equiv		nb301nf_charge,		48
.equiv		nb301nf_p_facel,	52
.equiv		nb301nf_p_krf,		56	
.equiv		nb301nf_p_crf,		60	
.equiv		nb301nf_Vc,		64	
.equiv		nb301nf_type,		68
.equiv		nb301nf_p_ntype,	72
.equiv		nb301nf_vdwparam,	76	
.equiv		nb301nf_Vvdw,		80	
.equiv		nb301nf_p_tabscale,	84	
.equiv		nb301nf_VFtab,		88
.equiv		nb301nf_invsqrta,	92	
.equiv		nb301nf_dvda,		96
.equiv          nb301nf_p_gbtabscale,   100
.equiv          nb301nf_GBtab,          104
.equiv          nb301nf_p_nthreads,     108
.equiv          nb301nf_count,          112
.equiv          nb301nf_mtx,            116
.equiv          nb301nf_outeriter,      120
.equiv          nb301nf_inneriter,      124
.equiv          nb301nf_work,           128
			;# stack offsets for local variables 
.equiv		nb301nf_is3,		0
.equiv		nb301nf_ii3,		4
.equiv		nb301nf_ixO,		8
.equiv		nb301nf_iyO,		12
.equiv		nb301nf_izO,		16	
.equiv		nb301nf_ixH,		20  
.equiv		nb301nf_iyH,		28  
.equiv		nb301nf_izH,		36  
.equiv		nb301nf_iqO,		44  
.equiv		nb301nf_iqH,		52  
.equiv		nb301nf_qqO,		60  
.equiv		nb301nf_qqH,		68  
.equiv		nb301nf_vctot,		76  
.equiv		nb301nf_n1,		84  
.equiv		nb301nf_tsc,		92 
.equiv		nb301nf_innerjjnr,	100
.equiv		nb301nf_innerk,		104
.equiv		nb301nf_tmprsqH,	108 
.equiv          nb301nf_n,              116 ;# idx for outer loop
.equiv          nb301nf_nn1,            120 ;# number of outer iterations
.equiv          nb301nf_nri,              124
.equiv          nb301nf_nouter,           128
.equiv          nb301nf_ninner,           132
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 136		;# local stack space 
	femms

	mov ecx, [ebp + nb301nf_p_nri]
	mov esi, [ebp + nb301nf_p_facel]
	mov edi, [ebp + nb301nf_p_tabscale]
	mov ecx, [ecx]
	mov [esp + nb301nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb301nf_nouter], eax
	mov [esp + nb301nf_ninner], eax

	mov   ecx, [ebp + nb301nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + nb301nf_charge]
	movd  mm1, [esi]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] 
	pfmul mm2, mm1		
	movq  [esp + nb301nf_iqO], mm2	    ;# iqO = facel*charge[ii] 
	
	movd  mm2, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] 
	pfmul mm2, mm1
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + nb301nf_iqH], mm2	    ;# iqH = facel*charge[ii0+1] 

	movd  mm4, [edi]
	punpckldq mm4,mm4	    ;# spread to both halves 
	movq  [esp + nb301nf_tsc], mm4	      
	;# assume we have at least one i particle - start directly 	 
.nb301nf_threadloop:
        mov   esi, [ebp + nb301nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb301nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [esi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb301nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb301nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb301nf_n], eax
        mov [esp + nb301nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb301nf_outerstart
        jmp .nb301nf_end

.nb301nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb301nf_nouter]
        mov [esp + nb301nf_nouter], ebx
	
.nb301nf_outer:
	mov   eax, [ebp + nb301nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb301nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb301nf_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + nb301nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb301nf_pos]    ;# eax = base of pos[] 
	
	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb301nf_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + nb301nf_ixO], mm5	
	movq  [esp + nb301nf_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + nb301nf_ixH], mm0	
	movq [esp + nb301nf_iyH], mm1	
	movq [esp + nb301nf_izH], mm2	
					
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + nb301nf_vctot], mm7

	mov   eax, [ebp + nb301nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + nb301nf_innerk], edx        
	add   edx, [esp + nb301nf_ninner]
	mov   [esp + nb301nf_ninner], edx

	mov   esi, [ebp + nb301nf_pos]	
	mov   eax, [ebp + nb301nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb301nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.nb301nf_inner_loop:	
	;# a single j particle iteration 
	mov   eax, [esp + nb301nf_innerjjnr]
	mov   eax, [eax]	 ;# eax=jnr offset 
    	add dword ptr [esp + nb301nf_innerjjnr],  4 ;# advance pointer 
	prefetch [ecx + 16]	   ;# prefetch data - trial and error says 16 is best 

	mov ecx, [ebp + nb301nf_charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + nb301nf_iqO]
	pfmul mm7, [esp + nb301nf_iqH]	 ;# mm6=qqO, mm7=qqH 
	movd [esp + nb301nf_qqO], mm6
	movq [esp + nb301nf_qqH], mm7
		
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
	
	pfsubr mm0, [esp + nb301nf_ixO]
	pfsubr mm1, [esp + nb301nf_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb301nf_ixH]
	pfsubr mm3, [esp + nb301nf_iyH]
	pfsubr mm4, [esp + nb301nf_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + nb301nf_tmprsqH], mm3
	
    	pfrsqrt mm1,mm0

    	movq mm2,mm1
    	pfmul mm1,mm1
    	pfrsqit1 mm1,mm0				
    	pfrcpit2 mm1,mm2	;# mm1=invsqrt 

	pfmul mm0, mm1		;# mm0=r 

	pfmul mm0, [esp + nb301nf_tsc]
	pf2iw mm4, mm0
	movd [esp + nb301nf_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + nb301nf_VFtab]
	mov ecx, [esp + nb301nf_n1]
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

	pfmul mm5, [esp + nb301nf_qqO]	;# vcoul=qq*VV 
	;# update vctot directly 
	pfadd mm5, [esp + nb301nf_vctot]
	movq [esp + nb301nf_vctot], mm5
	
	;# now do the two hydrogens. 
	movq mm0, [esp + nb301nf_tmprsqH] ;# mm0=rsqH 

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
	pfmul mm0, [esp + nb301nf_tsc]
	pf2iw mm4, mm0
	movq [esp + nb301nf_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + nb301nf_VFtab]
	mov ecx, [esp + nb301nf_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + nb301nf_n1 + 4]
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

	pfmul mm5, [esp + nb301nf_qqH]	;# vcoul=qq*VV 

	;# update vctot 
	pfadd mm5, [esp + nb301nf_vctot]
	movq [esp + nb301nf_vctot], mm5
		
	;#  done  - one more? 
	dec dword ptr [esp + nb301nf_innerk]
	jz  .nb301nf_updateouterdata
	jmp .nb301nf_inner_loop
.nb301nf_updateouterdata:	
	;# get n from stack
	mov esi, [esp + nb301nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb301nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb301nf_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb301nf_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

       	;# finish if last 
        mov ecx, [esp + nb301nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb301nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb301nf_n], esi
        jmp .nb301nf_outer
.nb301nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb301nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb301nf_end
        ;# non-zero, do one more workunit
        jmp   .nb301nf_threadloop
.nb301nf_end:
	femms

	mov eax, [esp + nb301nf_nouter] 	
	mov ebx, [esp + nb301nf_ninner]
	mov ecx, [ebp + nb301nf_outeriter]
	mov edx, [ebp + nb301nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	add esp, 136
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret

	
