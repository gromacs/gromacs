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




.globl nb_kernel311_ia32_3dnow
.globl _nb_kernel311_ia32_3dnow
nb_kernel311_ia32_3dnow:	
_nb_kernel311_ia32_3dnow:	
.equiv		nb311_p_nri,		8
.equiv		nb311_iinr,		12
.equiv		nb311_jindex,		16
.equiv		nb311_jjnr,		20
.equiv		nb311_shift,		24
.equiv		nb311_shiftvec,		28
.equiv		nb311_fshift,		32
.equiv		nb311_gid,		36
.equiv		nb311_pos,		40		
.equiv		nb311_faction,		44
.equiv		nb311_charge,		48
.equiv		nb311_p_facel,		52
.equiv		nb311_p_krf,		56	
.equiv		nb311_p_crf,		60	
.equiv		nb311_Vc,		64	
.equiv		nb311_type,		68
.equiv		nb311_p_ntype,		72
.equiv		nb311_vdwparam,		76	
.equiv		nb311_Vvdw,		80	
.equiv		nb311_p_tabscale,	84	
.equiv		nb311_VFtab,		88
.equiv		nb311_invsqrta,		92	
.equiv		nb311_dvda,		96
.equiv          nb311_p_gbtabscale,     100
.equiv          nb311_GBtab,            104
.equiv          nb311_p_nthreads,       108
.equiv          nb311_count,            112
.equiv          nb311_mtx,              116
.equiv          nb311_outeriter,        120
.equiv          nb311_inneriter,        124
.equiv          nb311_work,             128
			;# stack offsets for local variables 
.equiv		nb311_is3,		0
.equiv		nb311_ii3,		4
.equiv		nb311_ixO,		8
.equiv		nb311_iyO,		12
.equiv		nb311_izO,		16	
.equiv		nb311_ixH,		20  
.equiv		nb311_iyH,		28  
.equiv		nb311_izH,		36  
.equiv		nb311_iqO,		44  
.equiv		nb311_iqH,		52  
.equiv		nb311_qqO,		60  
.equiv		nb311_qqH,		68  
.equiv		nb311_vctot,		76  
.equiv		nb311_Vvdwtot,		84  
.equiv		nb311_c6,		92  
.equiv		nb311_c12,		100 
.equiv		nb311_six,		108 
.equiv		nb311_twelve,		116 
.equiv		nb311_two,		124 
.equiv		nb311_n1,		132 
.equiv		nb311_tsc,		140 
.equiv		nb311_ntia,		148 
.equiv		nb311_innerjjnr,	156
.equiv		nb311_innerk,		160	
.equiv		nb311_fixO,		164
.equiv		nb311_fiyO,		168
.equiv		nb311_fizO,		172
.equiv		nb311_fixH,		176 
.equiv		nb311_fiyH,		184 
.equiv		nb311_fizH,		192 
.equiv		nb311_dxO,		200
.equiv		nb311_dyO,		204
.equiv		nb311_dzO,		208
.equiv		nb311_dxH,		212 
.equiv		nb311_dyH,		220 
.equiv		nb311_dzH,		228 
.equiv		nb311_tmprsqH,		236 
.equiv          nb311_n,                244 ;# idx for outer loop
.equiv          nb311_nn1,              248 ;# number of outer iterations
.equiv          nb311_nri,              252
.equiv          nb311_nouter,           256
.equiv          nb311_ninner,           260
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 264		;# local stack space 
	femms

	mov ecx, [ebp + nb311_p_nri]
	mov esi, [ebp + nb311_p_facel]
	mov edi, [ebp + nb311_p_tabscale]
	mov ecx, [ecx]
	mov [esp + nb311_nri], ecx
	movd  mm1, [esi] 	;# facel
	mov esi, [ebp + nb311_p_ntype]

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb311_nouter], eax
	mov [esp + nb311_ninner], eax

	mov   ecx, [ebp + nb311_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + nb311_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] 
	pfmul mm2, mm1
	movq  [esp + nb311_iqO], mm2	    ;# iqO = facel*charge[ii] 
	
	movd  mm2, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] 
	pfmul mm2, mm1
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + nb311_iqH], mm2	    ;# iqH = facel*charge[ii0+1] 

	mov   edx, [ebp + nb311_type] 	
	mov   edx, [edx + ebx*4]
	shl   edx, 1
	mov   ecx, edx		        
	imul  ecx, [esi]      ;# ecx = ntia = 2*ntype*type[ii0]  
	mov   [esp + nb311_ntia], ecx
	 	
	movq  mm6, [edi]
	punpckldq mm6,mm6	    ;# spread to both halves 
	movq  [esp + nb311_tsc], mm6	      
	mov eax, 0x40000000    ;# 2.0
	mov [esp + nb311_two], eax
	mov [esp + nb311_two+4], eax
	mov ebx, 0x40c00000   ;# 6.0
	mov [esp + nb311_six], ebx
	mov [esp + nb311_six+4], ebx
	mov ecx, 0x41400000    ;# 12.0
	mov [esp + nb311_twelve], ecx
	mov [esp + nb311_twelve+4], ecx
	
.nb311_threadloop:
        mov   esi, [ebp + nb311_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb311_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb311_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb311_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb311_n], eax
        mov [esp + nb311_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb311_outerstart
        jmp .nb311_end

.nb311_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb311_nouter]
        mov [esp + nb311_nouter], ebx
	
.nb311_outer:
	mov   eax, [ebp + nb311_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb311_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb311_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + nb311_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb311_pos]    ;# eax = base of pos[] 
	
	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb311_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + nb311_ixO], mm5	
	movq  [esp + nb311_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + nb311_ixH], mm0	
	movq [esp + nb311_iyH], mm1	
	movq [esp + nb311_izH], mm2	
					
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + nb311_vctot], mm7
	movq  [esp + nb311_Vvdwtot], mm7
	movq  [esp + nb311_fixO],   mm7
	movd  [esp + nb311_fizO],   mm7
	movq  [esp + nb311_fixH],   mm7
	movq  [esp + nb311_fiyH],   mm7
	movq  [esp + nb311_fizH],   mm7

	mov   eax, [ebp + nb311_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + nb311_innerk], edx        
	add   edx, [esp + nb311_ninner]
	mov   [esp + nb311_ninner], edx

	mov   esi, [ebp + nb311_pos]
	mov   edi, [ebp + nb311_faction]	
	mov   eax, [ebp + nb311_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb311_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.nb311_inner_loop:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + nb311_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + nb311_innerjjnr],  4 ;# advance pointer 
	prefetch [ecx + 16]	   ;# prefetch data - trial and error says 16 is best 

	mov ecx, [ebp + nb311_charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + nb311_iqO]
	pfmul mm7, [esp + nb311_iqH]	;# mm6=qqO, mm7=qqH 
	movd [esp + nb311_qqO], mm6
	movq [esp + nb311_qqH], mm7

	mov ecx, [ebp + nb311_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr] 
	mov ecx, [ebp + nb311_vdwparam]
	shl edx, 1
	add edx, [esp + nb311_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [ecx + edx*4]		;# mm5 = 1st c6  		
	movq [esp + nb311_c6], mm5
	movd mm5, [ecx + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + nb311_c12], mm5	
			
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
	
	pfsubr mm0, [esp + nb311_ixO]
	pfsubr mm1, [esp + nb311_izO]
		
	movq  [esp + nb311_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + nb311_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb311_ixH]
	pfsubr mm3, [esp + nb311_iyH]
	pfsubr mm4, [esp + nb311_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + nb311_dxH], mm2
	movq [esp + nb311_dyH], mm3
	movq [esp + nb311_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + nb311_tmprsqH], mm3
	
    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 

	pfmul mm0, mm1		;# mm0=r 

	pfmul mm0, [esp + nb311_tsc]
	pf2iw mm4, mm0
	movd [esp + nb311_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + nb311_VFtab]
	mov ecx, [esp + nb311_n1]
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

	pfmul mm7, [esp + nb311_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + nb311_qqO]	;# vcoul=qq*VV 
	pfmul mm7, [esp + nb311_qqO]	;# fijC=qq*FF 
	;# update vctot directly, use mm3 for fscal sum. 
	pfadd mm5, [esp + nb311_vctot]
	movq [esp + nb311_vctot], mm5

	movq mm3, mm7
	pfmul mm3, [esp + nb311_tsc]
	
	;# nontabulated LJ - mm1 is invsqrt. - keep mm1! 
	movq mm0, mm1
	pfmul mm0, mm0		;# mm0 is invsq 
	movq mm2, mm0
	pfmul mm2, mm0
	pfmul mm2, mm0		;# mm2 = rinvsix 
	movq mm4, mm2
	pfmul mm4, mm4		;# mm4=rinvtwelve 

	pfmul mm4, [esp + nb311_c12]
	pfmul mm2, [esp + nb311_c6]
	movq mm5, mm4
	pfsub mm5, mm2		;# mm5=Vvdw12-Vvdw6 

	pfmul mm2, [esp + nb311_six]
	pfmul mm4, [esp + nb311_twelve]
	pfsub mm4, mm2
	pfmul mm4, mm1    ;# mm4=(12*Vvdw12-6*Vvdw6)*rinv11 

	pfsubr mm3, mm4 
 	pfmul mm3, mm1    ;# mm3 is total fscal (for the oxygen) now 
	
	;# update Vvdwtot  
	pfadd mm5, [esp + nb311_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb311_Vvdwtot], mm5       ;# store the sum       
	
	;# Ready with the oxygen - potential is updated, fscal is in mm3. 
	;# now do the two hydrogens. 
	movq mm0, [esp + nb311_tmprsqH] ;# mm0=rsqH 

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
	pfmul mm0, [esp + nb311_tsc]
	pf2iw mm4, mm0
	movq [esp + nb311_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + nb311_VFtab]
	mov ecx, [esp + nb311_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + nb311_n1+4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + nb311_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + nb311_qqH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + nb311_qqH]	;# fijC=qq*FF 
	;# update vctot 
	pfadd mm5, [esp + nb311_vctot]
	movq [esp + nb311_vctot], mm5
	
	;# change sign of fijC and multiply by rinv 
    pxor mm4,mm4
	pfsub mm4, mm7
	pfmul mm4, [esp + nb311_tsc]	
 	pfmul mm4, mm1    ;# mm4 is total fscal (for the hydrogens) now 	
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm3,mm3
	;# calc vectorial force for O 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm0,  [esp + nb311_dxO]
	movd mm1,  [esp + nb311_dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;# calc vectorial force for H's 
	movq mm5, [esp + nb311_dxH]
	movq mm6, [esp + nb311_dyH]
	movq mm7, [esp + nb311_dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;# update iO particle force 
	movq mm2,  [esp + nb311_fixO]
	movd mm3,  [esp + nb311_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + nb311_fixO], mm2
	movd [esp + nb311_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + nb311_fixH]
	movq mm3, [esp + nb311_fiyH]
	movq mm4, [esp + nb311_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + nb311_fixH], mm2
	movq [esp + nb311_fiyH], mm3
	movq [esp + nb311_fizH], mm4
	
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
	dec dword ptr [esp + nb311_innerk]
	jz  .nb311_updateouterdata
	jmp .nb311_inner_loop
.nb311_updateouterdata:	
	mov   ecx, [esp + nb311_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment iO force  
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + nb311_fixO]
	pfadd mm7, [esp + nb311_fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + nb311_fixH]
	movq  mm3, [esp + nb311_fiyH]
	movq  mm1, [esp + nb311_fizH]
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

	
	mov   ebx, [ebp + nb311_fshift]    ;# increment fshift force 
	mov   edx, [esp + nb311_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + nb311_fixO]
	pfadd mm7, [esp + nb311_fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	;# get n from stack
	mov esi, [esp + nb311_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb311_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb311_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb311_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + nb311_Vvdwtot]     
	pfacc mm7,mm7	          ;# same for Vvdw 
	
	mov   eax, [ebp + nb311_Vvdw]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment Vvdw[gid] 
       	;# finish if last 
        mov ecx, [esp + nb311_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb311_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb311_n], esi
        jmp .nb311_outer
.nb311_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb311_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb311_end
        ;# non-zero, do one more workunit
        jmp   .nb311_threadloop
.nb311_end:
	femms
	mov eax, [esp + nb311_nouter] 	
	mov ebx, [esp + nb311_ninner]
	mov ecx, [ebp + nb311_outeriter]
	mov edx, [ebp + nb311_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	add esp, 264
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


.globl nb_kernel311nf_ia32_3dnow
.globl _nb_kernel311nf_ia32_3dnow
nb_kernel311nf_ia32_3dnow:	
_nb_kernel311nf_ia32_3dnow:	
.equiv		nb311nf_p_nri,		8
.equiv		nb311nf_iinr,		12
.equiv		nb311nf_jindex,		16
.equiv		nb311nf_jjnr,		20
.equiv		nb311nf_shift,		24
.equiv		nb311nf_shiftvec,	28
.equiv		nb311nf_fshift,		32
.equiv		nb311nf_gid,		36
.equiv		nb311nf_pos,		40		
.equiv		nb311nf_faction,	44
.equiv		nb311nf_charge,		48
.equiv		nb311nf_p_facel,	52
.equiv		nb311nf_p_krf,		56	
.equiv		nb311nf_p_crf,		60	
.equiv		nb311nf_Vc,		64	
.equiv		nb311nf_type,		68
.equiv		nb311nf_p_ntype,	72
.equiv		nb311nf_vdwparam,	76	
.equiv		nb311nf_Vvdw,		80	
.equiv		nb311nf_p_tabscale,	84	
.equiv		nb311nf_VFtab,		88
.equiv		nb311nf_invsqrta,	92	
.equiv		nb311nf_dvda,		96
.equiv          nb311nf_p_gbtabscale,   100
.equiv          nb311nf_GBtab,          104
.equiv          nb311nf_p_nthreads,     108
.equiv          nb311nf_count,          112
.equiv          nb311nf_mtx,            116
.equiv          nb311nf_outeriter,      120
.equiv          nb311nf_inneriter,      124
.equiv          nb311nf_work,           128
			;# stack offsets for local variables 
.equiv		nb311nf_is3,		0
.equiv		nb311nf_ii3,		4
.equiv		nb311nf_ixO,		8
.equiv		nb311nf_iyO,		12
.equiv		nb311nf_izO,		16	
.equiv		nb311nf_ixH,		20  
.equiv		nb311nf_iyH,		28  
.equiv		nb311nf_izH,		36  
.equiv		nb311nf_iqO,		44  
.equiv		nb311nf_iqH,		52  
.equiv		nb311nf_qqO,		60  
.equiv		nb311nf_qqH,		68  
.equiv		nb311nf_vctot,		76  
.equiv		nb311nf_Vvdwtot,	84  
.equiv		nb311nf_c6,		92  
.equiv		nb311nf_c12,		100 
.equiv		nb311nf_n1,		108
.equiv		nb311nf_tsc,		116 
.equiv		nb311nf_ntia,		124 
.equiv		nb311nf_innerjjnr,	128
.equiv		nb311nf_innerk,		132
.equiv		nb311nf_tmprsqH,	136 
.equiv          nb311nf_n,              144 ;# idx for outer loop
.equiv          nb311nf_nn1,            148 ;# number of outer iterations
.equiv          nb311nf_nri,            152
.equiv          nb311nf_nouter,         156
.equiv          nb311nf_ninner,         160
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 164		;# local stack space 
	femms

	mov ecx, [ebp + nb311nf_p_nri]
	mov esi, [ebp + nb311nf_p_facel]
	mov edi, [ebp + nb311nf_p_tabscale]
	mov ecx, [ecx]
	mov [esp + nb311nf_nri], ecx
	movd  mm1, [esi] 	;# facel
	mov esi, [ebp + nb311nf_p_ntype]

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb311nf_nouter], eax
	mov [esp + nb311nf_ninner], eax


	mov   ecx, [ebp + nb311nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + nb311nf_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] 
	pfmul mm2, mm1
	movq  [esp + nb311nf_iqO], mm2	    ;# iqO = facel*charge[ii] 
	
	movd  mm2, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] 
	pfmul mm2, mm1
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + nb311nf_iqH], mm2	    ;# iqH = facel*charge[ii0+1] 

	mov   edx, [ebp + nb311nf_type] 	
	mov   edx, [edx + ebx*4]
	shl   edx, 1
	mov   ecx, edx		        
	imul  ecx, [esi]      ;# ecx = ntia = 2*ntype*type[ii0]  
	mov   [esp + nb311nf_ntia], ecx
	 	
	movq  mm6, [edi]
	punpckldq mm6,mm6	    ;# spread to both halves 
	movq  [esp + nb311nf_tsc], mm6	      
.nb311nf_threadloop:
        mov   esi, [ebp + nb311nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb311nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb311nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb311nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb311nf_n], eax
        mov [esp + nb311nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb311nf_outerstart
        jmp .nb311nf_end

.nb311nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb311nf_nouter]
        mov [esp + nb311nf_nouter], ebx
	
.nb311nf_outer:
	mov   eax, [ebp + nb311nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb311nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb311nf_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + nb311nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb311nf_pos]    ;# eax = base of pos[] 
	
	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + nb311nf_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + nb311nf_ixO], mm5	
	movq  [esp + nb311nf_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + nb311nf_ixH], mm0	
	movq [esp + nb311nf_iyH], mm1	
	movq [esp + nb311nf_izH], mm2	
					
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + nb311nf_vctot], mm7
	movq  [esp + nb311nf_Vvdwtot], mm7

	mov   eax, [ebp + nb311nf_jindex]
	mov   ecx, [eax +esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + nb311nf_innerk], edx        
	add   edx, [esp + nb311nf_ninner]
	mov   [esp + nb311nf_ninner], edx

	mov   esi, [ebp + nb311nf_pos]	
	mov   eax, [ebp + nb311nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb311nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.nb311nf_inner_loop:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + nb311nf_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    	add dword ptr [esp + nb311nf_innerjjnr],  4 ;# advance pointer 
	prefetch [ecx + 16]	   ;# prefetch data - trial and error says 16 is best 

	mov ecx, [ebp + nb311nf_charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + nb311nf_iqO]
	pfmul mm7, [esp + nb311nf_iqH]	;# mm6=qqO, mm7=qqH 
	movd [esp + nb311nf_qqO], mm6
	movq [esp + nb311nf_qqH], mm7

	mov ecx, [ebp + nb311nf_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr] 
	mov ecx, [ebp + nb311nf_vdwparam]
	shl edx, 1
	add edx, [esp + nb311nf_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [ecx + edx*4]		;# mm5 = 1st c6  		
	movq [esp + nb311nf_c6], mm5
	movd mm5, [ecx + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + nb311nf_c12], mm5	
			
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
	
	pfsubr mm0, [esp + nb311nf_ixO]
	pfsubr mm1, [esp + nb311nf_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + nb311nf_ixH]
	pfsubr mm3, [esp + nb311nf_iyH]
	pfsubr mm4, [esp + nb311nf_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + nb311nf_tmprsqH], mm3
	
    	pfrsqrt mm1,mm0

    	movq mm2,mm1
    	pfmul mm1,mm1
    	pfrsqit1 mm1,mm0				
    	pfrcpit2 mm1,mm2	;# mm1=invsqrt 

	pfmul mm0, mm1		;# mm0=r 

	pfmul mm0, [esp + nb311nf_tsc]
	pf2iw mm4, mm0
	movd [esp + nb311nf_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + nb311nf_VFtab]
	mov ecx, [esp + nb311nf_n1]
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

	pfmul mm5, [esp + nb311nf_qqO]	;# vcoul=qq*VV 
	;# update vctot directly 
	pfadd mm5, [esp + nb311nf_vctot]
	movq [esp + nb311nf_vctot], mm5
	
	;# nontabulated LJ - mm1 is invsqrt. - keep mm1! 
	movq mm0, mm1
	pfmul mm0, mm0		;# mm0 is invsq 
	movq mm2, mm0
	pfmul mm2, mm0
	pfmul mm2, mm0		;# mm2 = rinvsix 
	movq mm4, mm2
	pfmul mm4, mm4		;# mm4=rinvtwelve 

	pfmul mm4, [esp + nb311nf_c12]
	pfmul mm2, [esp + nb311nf_c6]
	pfsub mm4, mm2		;# mm4=Vvdw12-Vvdw6 

	;# update Vvdwtot  
	pfadd mm4, [esp + nb311nf_Vvdwtot]      ;# add the earlier value 
	movq [esp + nb311nf_Vvdwtot], mm4       ;# store the sum       

		;# now do the two hydrogens. 
	movq mm0, [esp + nb311nf_tmprsqH] ;# mm0=rsqH 

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
	pfmul mm0, [esp + nb311nf_tsc]
	pf2iw mm4, mm0
	movq [esp + nb311nf_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + nb311nf_VFtab]
	mov ecx, [esp + nb311nf_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + nb311nf_n1+4]
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

	pfmul mm5, [esp + nb311nf_qqH]	;# vcoul=qq*VV 
	;# update vctot 
	pfadd mm5, [esp + nb311nf_vctot]
	movq [esp + nb311nf_vctot], mm5

	;#  done  - one more? 
	dec dword ptr [esp + nb311nf_innerk]
	jz  .nb311nf_updateouterdata
	jmp .nb311nf_inner_loop
.nb311nf_updateouterdata:	
	;# get n from stack
	mov esi, [esp + nb311nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb311nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movq  mm7, [esp + nb311nf_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + nb311nf_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + nb311nf_Vvdwtot]     
	pfacc mm7,mm7	          ;# same for Vvdw 
	
	mov   eax, [ebp + nb311nf_Vvdw]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment Vvdw[gid] 
       	;# finish if last 
        mov ecx, [esp + nb311nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb311nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb311nf_n], esi
        jmp .nb311nf_outer
.nb311nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb311nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb311nf_end
        ;# non-zero, do one more workunit
        jmp   .nb311nf_threadloop
.nb311nf_end:
	femms
	mov eax, [esp + nb311nf_nouter] 	
	mov ebx, [esp + nb311nf_ninner]
	mov ecx, [ebp + nb311nf_outeriter]
	mov edx, [ebp + nb311nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	add esp, 164
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret
