;#
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
%macro .equiv                  2
   %1 equ %2
%endmacro
; .endif                   # End of NASM-specific block
; .intel_syntax noprefix   # Line only read by gnu as

.section .text

	
.globl nb_kernel233_ia32_sse
.globl _nb_kernel233_ia32_sse
nb_kernel233_ia32_sse:	
_nb_kernel233_ia32_sse:	
.equiv          nb233_p_nri,            8
.equiv          nb233_iinr,             12
.equiv          nb233_jindex,           16
.equiv          nb233_jjnr,             20
.equiv          nb233_shift,            24
.equiv          nb233_shiftvec,         28
.equiv          nb233_fshift,           32
.equiv          nb233_gid,              36
.equiv          nb233_pos,              40
.equiv          nb233_faction,          44
.equiv          nb233_charge,           48
.equiv          nb233_p_facel,          52
.equiv          nb233_argkrf,           56
.equiv          nb233_argcrf,           60
.equiv          nb233_Vc,               64
.equiv          nb233_type,             68
.equiv          nb233_p_ntype,          72
.equiv          nb233_vdwparam,         76
.equiv          nb233_Vvdw,             80
.equiv          nb233_p_tabscale,       84
.equiv          nb233_VFtab,            88
.equiv          nb233_invsqrta,         92
.equiv          nb233_dvda,             96
.equiv          nb233_p_gbtabscale,     100
.equiv          nb233_GBtab,            104
.equiv          nb233_p_nthreads,       108
.equiv          nb233_count,            112
.equiv          nb233_mtx,              116
.equiv          nb233_outeriter,        120
.equiv          nb233_inneriter,        124
.equiv          nb233_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb233_ixO,              0
.equiv          nb233_iyO,              16
.equiv          nb233_izO,              32
.equiv          nb233_ixH1,             48
.equiv          nb233_iyH1,             64
.equiv          nb233_izH1,             80
.equiv          nb233_ixH2,             96
.equiv          nb233_iyH2,             112
.equiv          nb233_izH2,             128
.equiv          nb233_ixM,              144
.equiv          nb233_iyM,              160
.equiv          nb233_izM,              176
.equiv          nb233_iqM,              192
.equiv          nb233_iqH,              208
.equiv          nb233_dxO,              224
.equiv          nb233_dyO,              240
.equiv          nb233_dzO,              256
.equiv          nb233_dxH1,             272
.equiv          nb233_dyH1,             288
.equiv          nb233_dzH1,             304
.equiv          nb233_dxH2,             320
.equiv          nb233_dyH2,             336
.equiv          nb233_dzH2,             352
.equiv          nb233_dxM,              368
.equiv          nb233_dyM,              384
.equiv          nb233_dzM,              400
.equiv          nb233_qqM,              416
.equiv          nb233_qqH,              432
.equiv          nb233_rinvH1,           448
.equiv          nb233_rinvH2,           464
.equiv          nb233_rinvM,            480
.equiv          nb233_two,              496
.equiv          nb233_c6,               512
.equiv          nb233_c12,              528
.equiv          nb233_tsc,              544
.equiv          nb233_fstmp,            560
.equiv          nb233_krf,              576
.equiv          nb233_crf,              592
.equiv          nb233_krsqH1,           608
.equiv          nb233_krsqH2,           624
.equiv          nb233_krsqM,            640
.equiv          nb233_vctot,            656
.equiv          nb233_Vvdwtot,          672
.equiv          nb233_fixO,             688
.equiv          nb233_fiyO,             704
.equiv          nb233_fizO,             720
.equiv          nb233_fixH1,            736
.equiv          nb233_fiyH1,            752
.equiv          nb233_fizH1,            768
.equiv          nb233_fixH2,            784
.equiv          nb233_fiyH2,            800
.equiv          nb233_fizH2,            816
.equiv          nb233_fixM,             832
.equiv          nb233_fiyM,             848
.equiv          nb233_fizM,             864
.equiv          nb233_fjx,              880
.equiv          nb233_fjy,              896
.equiv          nb233_fjz,              912
.equiv          nb233_half,             928
.equiv          nb233_three,            944
.equiv          nb233_is3,              960
.equiv          nb233_ii3,              964
.equiv          nb233_ntia,             968
.equiv          nb233_innerjjnr,        972
.equiv          nb233_innerk,           976
.equiv          nb233_n,                980
.equiv          nb233_nn1,              984
.equiv          nb233_nri,              988
.equiv          nb233_nouter,           992
.equiv          nb233_ninner,           996
.equiv          nb233_salign,           1000
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx
	push esi
	push edi
	sub esp, 1004		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb233_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb233_p_nri]
	mov ecx, [ecx]
	mov [esp + nb233_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb233_nouter], eax
	mov [esp + nb233_ninner], eax

	mov eax, [ebp + nb233_p_tabscale]
	movss xmm3, [eax]
	shufps xmm3, xmm3, 0
	movaps [esp + nb233_tsc], xmm3

	mov esi, [ebp + nb233_argkrf]
	mov edi, [ebp + nb233_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	movaps [esp + nb233_krf], xmm5
	movaps [esp + nb233_crf], xmm6
	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb233_half], eax
	movss xmm1, [esp + nb233_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb233_half],  xmm1
	movaps [esp + nb233_two],  xmm2
	movaps [esp + nb233_three],  xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb233_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb233_charge]
	movss xmm4, [edx + ebx*4 + 4]	
	movss xmm3, [edx + ebx*4 + 12]	
	mov esi, [ebp + nb233_p_facel]
	movss xmm5, [esi]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [esp + nb233_iqM], xmm3
	movaps [esp + nb233_iqH], xmm4
	
	mov   edx, [ebp + nb233_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov edi, [ebp + nb233_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	mov   [esp + nb233_ntia], ecx		
.nb233_threadloop:
        mov   esi, [ebp + nb233_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb233_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb233_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb233_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb233_n], eax
        mov [esp + nb233_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb233_outerstart
        jmp .nb233_end
	
.nb233_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb233_nouter]
	mov [esp + nb233_nouter], ebx

.nb233_outer:
	mov   eax, [ebp + nb233_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]	;# ebx=3*is 
	mov   [esp + nb233_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb233_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb233_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]		;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2	
	movaps xmm6, xmm0
	movaps xmm7, xmm1
	
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb233_pos]	;# eax = base of pos[]  
	mov   [esp + nb233_ii3], ebx

	addss xmm3, [eax + ebx*4]  	;# ox
	addss xmm4, [eax + ebx*4 + 4]  ;# oy
	addss xmm5, [eax + ebx*4 + 8]  ;# oz
	addss xmm6, [eax + ebx*4 + 12] ;# h1x
	addss xmm7, [eax + ebx*4 + 16] ;# h1y
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	shufps xmm7, xmm7, 0
	movaps [esp + nb233_ixO], xmm3
	movaps [esp + nb233_iyO], xmm4
	movaps [esp + nb233_izO], xmm5
	movaps [esp + nb233_ixH1], xmm6
	movaps [esp + nb233_iyH1], xmm7

	movss xmm6, xmm2
	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm6, [eax + ebx*4 + 20] ;# h1z
	addss xmm0, [eax + ebx*4 + 24] ;# h2x
	addss xmm1, [eax + ebx*4 + 28] ;# h2y
	addss xmm2, [eax + ebx*4 + 32] ;# h2z
	addss xmm3, [eax + ebx*4 + 36] ;# mx
	addss xmm4, [eax + ebx*4 + 40] ;# my
	addss xmm5, [eax + ebx*4 + 44] ;# mz

	shufps xmm6, xmm6, 0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb233_izH1], xmm6
	movaps [esp + nb233_ixH2], xmm0
	movaps [esp + nb233_iyH2], xmm1
	movaps [esp + nb233_izH2], xmm2
	movaps [esp + nb233_ixM], xmm3
	movaps [esp + nb233_iyM], xmm4
	movaps [esp + nb233_izM], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb233_vctot], xmm4
	movaps [esp + nb233_Vvdwtot], xmm4
	movaps [esp + nb233_fixO], xmm4
	movaps [esp + nb233_fiyO], xmm4
	movaps [esp + nb233_fizO], xmm4
	movaps [esp + nb233_fixH1], xmm4
	movaps [esp + nb233_fiyH1], xmm4
	movaps [esp + nb233_fizH1], xmm4
	movaps [esp + nb233_fixH2], xmm4
	movaps [esp + nb233_fiyH2], xmm4
	movaps [esp + nb233_fizH2], xmm4
	movaps [esp + nb233_fixM], xmm4
	movaps [esp + nb233_fiyM], xmm4
	movaps [esp + nb233_fizM], xmm4
	
	mov   eax, [ebp + nb233_jindex]
	mov   ecx, [eax + esi*4]	 	;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   esi, [ebp + nb233_pos]
	mov   edi, [ebp + nb233_faction]	
	mov   eax, [ebp + nb233_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb233_innerjjnr], eax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb233_ninner]
	mov   [esp + nb233_ninner], ecx
	add   edx, 0
	mov   [esp + nb233_innerk], edx	;# number of innerloop atoms 
	jge   .nb233_unroll_loop
	jmp   .nb233_odd_inner
.nb233_unroll_loop:
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb233_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]     	;# eax-edx=jnr1-4 

	add dword ptr [esp + nb233_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb233_charge]	;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	 	;# and in xmm4 
	mulps  xmm3, [esp + nb233_iqM]
	mulps  xmm4, [esp + nb233_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [esp + nb233_qqM], xmm3
	movaps  [esp + nb233_qqH], xmm4
	
	mov esi, [ebp + nb233_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + nb233_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + nb233_ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm6, xmm7, 221  ;# constant 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + nb233_c6], xmm4
	movaps [esp + nb233_c12], xmm6

	mov esi, [ebp + nb233_pos]   	;# base of pos[] 

	lea   eax, [eax + eax*2] 	;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2] 	;# replace jnr with j3 
	lea   edx, [edx + edx*2]	

	;# move four coordinates to xmm0-xmm2 	
	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# constant 10001000
	
	shufps xmm0, xmm5, 136  ;# constant 10001000
	shufps xmm1, xmm5, 221  ;# constant 11011101		

	;# move ixO-izO to xmm4-xmm6 
	movaps xmm4, [esp + nb233_ixO]
	movaps xmm5, [esp + nb233_iyO]
	movaps xmm6, [esp + nb233_izO]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb233_dxO], xmm4
	movaps [esp + nb233_dyO], xmm5
	movaps [esp + nb233_dzO], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	;# rsqO in xmm7

	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [esp + nb233_ixH1]
	movaps xmm5, [esp + nb233_iyH1]
	movaps xmm6, [esp + nb233_izH1]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb233_dxH1], xmm4
	movaps [esp + nb233_dyH1], xmm5
	movaps [esp + nb233_dzH1], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	;# rsqH1 in xmm6 
		
	;# move ixH2-izH2 to xmm3-xmm5  
	movaps xmm3, [esp + nb233_ixH2]
	movaps xmm4, [esp + nb233_iyH2]
	movaps xmm5, [esp + nb233_izH2]

	;# calc dr 
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# store dr 
	movaps [esp + nb233_dxH2], xmm3
	movaps [esp + nb233_dyH2], xmm4
	movaps [esp + nb233_dzH2], xmm5
	;# square it 
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	
	;# move ixM-izM to xmm2-xmm4  
	movaps xmm3, [esp + nb233_iyM]
	movaps xmm4, [esp + nb233_izM]
	subps  xmm3, xmm1
	subps  xmm4, xmm2
	movaps xmm2, [esp + nb233_ixM]
	subps  xmm2, xmm0	

	;# store dr 
	movaps [esp + nb233_dxM], xmm2
	movaps [esp + nb233_dyM], xmm3
	movaps [esp + nb233_dzM], xmm4
	;# square it 
	mulps xmm2,xmm2
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	addps xmm4, xmm3
	addps xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 
	movaps xmm0, xmm4
	movaps xmm1, xmm5
	movaps xmm2, xmm6
	mulps  xmm0, [esp + nb233_krf]	
	mulps  xmm1, [esp + nb233_krf]	
	mulps  xmm2, [esp + nb233_krf]	
	movaps [esp + nb233_krsqM], xmm0
	movaps [esp + nb233_krsqH2], xmm1
	movaps [esp + nb233_krsqH1], xmm2

	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm0, [esp + nb233_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [esp + nb233_half]
	movaps  [esp + nb233_rinvH1], xmm0	;# rinvH1 

	;# rsqH2 - seed to xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm0, [esp + nb233_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [esp + nb233_half]
	movaps  [esp + nb233_rinvH2], xmm0	;# rinvH2 

	;# rsqM - seed to xmm2 
	rsqrtps xmm2, xmm4
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm0, [esp + nb233_three]
	mulps   xmm2, xmm4	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [esp + nb233_half]
	movaps  [esp + nb233_rinvM], xmm0
	
	;# Do the O LJ-only interaction directly.	
	;# rsqO is in xmm7
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb233_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb233_half]
	movaps  xmm0, xmm4
	;# xmm0=rinvO
	
	mulps xmm7, xmm0
	mulps xmm7, [esp + nb233_tsc] ;# rtab
	
	movhlps xmm5, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm5	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps  xmm7, xmm6	
	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3
	pslld mm7, 3

	movd mm0, eax	
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx	

	mov  esi, [ebp + nb233_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7
	
	;# dispersion 
	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half dispersion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of dispersion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm7, [esp + nb233_two]	;# two*Heps2 
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [esp + nb233_c6]
	mulps  xmm7, xmm4	 ;# fijD 
	mulps  xmm5, xmm4	 ;# Vvdw6 
	mulps  xmm7, [esp + nb233_tsc]
	;# put scalar force on stack Update Vvdwtot directly 
	addps  xmm5, [esp + nb233_Vvdwtot]
	movaps [esp + nb233_fstmp], xmm7
	movaps [esp + nb233_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [esi + eax*4 + 16]
	movlps xmm7, [esi + ecx*4 + 16]
	movhps xmm5, [esi + ebx*4 + 16]
	movhps xmm7, [esi + edx*4 + 16] ;# got half repulsion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 24]
	movlps xmm3, [esi + ecx*4 + 24]
	movhps xmm7, [esi + ebx*4 + 24]
	movhps xmm3, [esi + edx*4 + 24] ;# other half of repulsion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm7, [esp + nb233_two]	;# two*Heps2 
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 	
	movaps xmm4, [esp + nb233_c12]
	mulps  xmm7, xmm4 ;# fijR 
	mulps  xmm5, xmm4 ;# Vvdw12 
	mulps  xmm7, [esp + nb233_tsc]
	addps  xmm7, [esp + nb233_fstmp]

	addps  xmm5, [esp + nb233_Vvdwtot]
	movaps [esp + nb233_Vvdwtot], xmm5

	movd  eax, mm0
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	xorps  xmm1, xmm1
	mulps  xmm7, xmm0
	subps  xmm1, xmm7
	
	movaps xmm3, [esp + nb233_dxO]
	movaps xmm4, [esp + nb233_dyO]
	movaps xmm5, [esp + nb233_dzO]
	mulps  xmm3, xmm1
	mulps  xmm4, xmm1
	mulps  xmm5, xmm1	;# tx in xmm3-xmm5

	;# update O forces 
	movaps xmm0, [esp + nb233_fixO]
	movaps xmm1, [esp + nb233_fiyO]
	movaps xmm2, [esp + nb233_fizO]
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5
	movaps [esp + nb233_fixO], xmm0
	movaps [esp + nb233_fiyO], xmm1
	movaps [esp + nb233_fizO], xmm2
	;# update j forces with water O 
	movaps [esp + nb233_fjx], xmm3
	movaps [esp + nb233_fjy], xmm4
	movaps [esp + nb233_fjz], xmm5

	;# Do H1 interaction
	movaps  xmm7, [esp + nb233_rinvH1]
	movaps  xmm4, xmm7	
	mulps   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq
	movaps xmm0, xmm7
	movaps xmm1, [esp + nb233_krsqH1]
	addps  xmm0, xmm1
	subps  xmm0, [esp + nb233_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm1, [esp + nb233_two]
	subps  xmm7, xmm1
	mulps  xmm0, [esp + nb233_qqH]
	mulps  xmm7, [esp + nb233_qqH]

	mulps  xmm4, xmm7	;# total fs H1 in xmm4 

	addps  xmm0, [esp + nb233_vctot]	
	movaps [esp + nb233_vctot], xmm0

	movaps xmm0, [esp + nb233_dxH1]
	movaps xmm1, [esp + nb233_dyH1]
	movaps xmm2, [esp + nb233_dzH1]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	;# update H1 forces 
	movaps xmm3, [esp + nb233_fixH1]
	movaps xmm4, [esp + nb233_fiyH1]
	movaps xmm7, [esp + nb233_fizH1]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb233_fixH1], xmm3
	movaps [esp + nb233_fiyH1], xmm4
	movaps [esp + nb233_fizH1], xmm7
	;# update j forces with water H1 
	addps  xmm0, [esp + nb233_fjx]
	addps  xmm1, [esp + nb233_fjy]
	addps  xmm2, [esp + nb233_fjz]
	movaps [esp + nb233_fjx], xmm0
	movaps [esp + nb233_fjy], xmm1
	movaps [esp + nb233_fjz], xmm2

	;# Done with H1, do H2 interactions
	movaps  xmm7, [esp + nb233_rinvH2]
	movaps  xmm4, xmm7	
	mulps   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq
	movaps xmm0, xmm7
	movaps xmm1, [esp + nb233_krsqH2]
	addps  xmm0, xmm1
	subps  xmm0, [esp + nb233_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm1, [esp + nb233_two]
	subps  xmm7, xmm1
	mulps  xmm0, [esp + nb233_qqH]
	mulps  xmm7, [esp + nb233_qqH]

	mulps  xmm4, xmm7	;# total fs H2 in xmm4 

	addps  xmm0, [esp + nb233_vctot]	
	movaps [esp + nb233_vctot], xmm0

	movaps xmm0, [esp + nb233_dxH2]
	movaps xmm1, [esp + nb233_dyH2]
	movaps xmm2, [esp + nb233_dzH2]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	;# update H2 forces 
	movaps xmm3, [esp + nb233_fixH2]
	movaps xmm4, [esp + nb233_fiyH2]
	movaps xmm7, [esp + nb233_fizH2]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb233_fixH2], xmm3
	movaps [esp + nb233_fiyH2], xmm4
	movaps [esp + nb233_fizH2], xmm7
	addps xmm0, [esp + nb233_fjx]
        addps xmm1, [esp + nb233_fjy]
        addps xmm2, [esp + nb233_fjz]
	movaps [esp + nb233_fjx], xmm0
	movaps [esp + nb233_fjy], xmm1
	movaps [esp + nb233_fjz], xmm2

	;# Done with H2, do M interactions
	movaps  xmm7, [esp + nb233_rinvM]
	movaps  xmm4, xmm7	
	mulps   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq
	movaps xmm0, xmm7
	movaps xmm1, [esp + nb233_krsqM]
	addps  xmm0, xmm1
	subps  xmm0, [esp + nb233_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm1, [esp + nb233_two]
	subps  xmm7, xmm1
	mulps  xmm0, [esp + nb233_qqM]
	mulps  xmm7, [esp + nb233_qqM]

	mulps  xmm4, xmm7	;# total fs M in xmm4 

	addps  xmm0, [esp + nb233_vctot]	
	movaps [esp + nb233_vctot], xmm0

	movaps xmm0, [esp + nb233_dxM]
	movaps xmm1, [esp + nb233_dyM]
	movaps xmm2, [esp + nb233_dzM]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	
	;# update M forces 
	movaps xmm3, [esp + nb233_fixM]
	movaps xmm4, [esp + nb233_fiyM]
	movaps xmm7, [esp + nb233_fizM]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb233_fixM], xmm3
	movaps [esp + nb233_fiyM], xmm4
	movaps [esp + nb233_fizM], xmm7

	mov edi, [ebp + nb233_faction]
	;# update j forces from stored values
	addps xmm0, [esp + nb233_fjx]
	addps xmm1, [esp + nb233_fjy]
	addps xmm2, [esp + nb233_fjz]

	movlps xmm4, [edi + eax*4]
	movlps xmm7, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm7, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm7, 136  ;# constant 10001000
	shufps xmm4, xmm7, 221  ;# constant 11011101
	
	;# xmm3 has fjx, xmm4 has fjy 
	subps xmm3, xmm0
	subps xmm4, xmm1
	;# unpack them back for storing 
	movaps xmm7, xmm3
	unpcklps xmm7, xmm4
	unpckhps xmm3, xmm4	
	movlps [edi + eax*4], xmm7
	movlps [edi + ecx*4], xmm3
	movhps [edi + ebx*4], xmm7
	movhps [edi + edx*4], xmm3
	;# finally z forces 
	movss  xmm0, [edi + eax*4 + 8]
	movss  xmm1, [edi + ebx*4 + 8]
	movss  xmm3, [edi + ecx*4 + 8]
	movss  xmm4, [edi + edx*4 + 8]
	subss  xmm0, xmm2
	shufps xmm2, xmm2, 229  ;# constant 11100101
	subss  xmm1, xmm2
	shufps xmm2, xmm2, 234  ;# constant 11101010
	subss  xmm3, xmm2
	shufps xmm2, xmm2, 255  ;# constant 11111111
	subss  xmm4, xmm2
	movss  [edi + eax*4 + 8], xmm0
	movss  [edi + ebx*4 + 8], xmm1
	movss  [edi + ecx*4 + 8], xmm3
	movss  [edi + edx*4 + 8], xmm4
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb233_innerk],  4
	jl    .nb233_odd_inner
	jmp   .nb233_unroll_loop
.nb233_odd_inner:	
	add dword ptr [esp + nb233_innerk],  4
	jnz   .nb233_odd_loop
	jmp   .nb233_updateouterdata
.nb233_odd_loop:
	mov   edx, [esp + nb233_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb233_innerjjnr],  4	

 	xorps xmm4, xmm4  	;# clear reg.
	movss xmm4, [esp + nb233_iqM]
	mov esi, [ebp + nb233_charge] 
	movhps xmm4, [esp + nb233_iqH]  ;# [qM  0  qH  qH] 
	shufps xmm4, xmm4, 41	;# [0 qH qH qM]

	movss xmm3, [esi + eax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [esp + nb233_qqM], xmm3	;# use dummy qq for storage 
	
	xorps xmm6, xmm6
	mov esi, [ebp + nb233_type]
	mov ebx, [esi + eax*4]
	mov esi, [ebp + nb233_vdwparam]
	shl ebx, 1	
	add ebx, [esp + nb233_ntia]
	movlps xmm6, [esi + ebx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# constant 11111100
	shufps xmm7, xmm7, 253  ;# constant 11111101
	movaps [esp + nb233_c6], xmm6
	movaps [esp + nb233_c12], xmm7
	
	mov esi, [ebp + nb233_pos]
	lea eax, [eax + eax*2]  

	movss xmm3, [esp + nb233_ixO]
	movss xmm4, [esp + nb233_iyO]
	movss xmm5, [esp + nb233_izO]
	movss xmm0, [esp + nb233_ixH1]
	movss xmm1, [esp + nb233_iyH1]
	movss xmm2, [esp + nb233_izH1]
	unpcklps xmm3, [esp + nb233_ixH2] 	;# ixO ixH2 - -
	unpcklps xmm4, [esp + nb233_iyH2]  	;# iyO iyH2 - -
	unpcklps xmm5, [esp + nb233_izH2]	;# izO izH2 - -
	unpcklps xmm0, [esp + nb233_ixM] 	;# ixH1 ixM - -
	unpcklps xmm1, [esp + nb233_iyM]  	;# iyH1 iyM - -
	unpcklps xmm2, [esp + nb233_izM]	;# izH1 izM - -
	unpcklps xmm3, xmm0  	;# ixO ixH1 ixH2 ixM
	unpcklps xmm4, xmm1 	;# same for y
	unpcklps xmm5, xmm2 	;# same for z
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# use O distances for storage
	movaps [esp + nb233_dxO], xmm3
	movaps [esp + nb233_dyO], xmm4
	movaps [esp + nb233_dzO], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 
	movaps xmm0, xmm4
	mulps xmm0, [esp + nb233_krf]
	movaps [esp + nb233_krsqM], xmm0
		
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb233_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb233_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv, xmm4=rsq

	;# LJ table interaction
	mulps xmm4, xmm0
	mulps  xmm4, [esp + nb233_tsc] ;# rtab
	
	cvttps2pi mm6, xmm4
	cvtpi2ps xmm6, mm6
	subss  xmm4, xmm6	
	movss xmm1, xmm4	;# xmm1=eps 
	movss xmm2, xmm1	
	mulss  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3

	movd mm0, eax	

	mov  esi, [ebp + nb233_VFtab]
	movd eax, mm6
	
	;# dispersion 
	movlps xmm5, [esi + eax*4]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [esi + eax*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm7, [esp + nb233_two]	;# two*Heps2 
	addss  xmm7, xmm6
	addss  xmm7, xmm5 ;# xmm7=FF 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 

	movss xmm4, [esp + nb233_c6]
	mulss  xmm7, xmm4	 ;# fijD 
	mulss  xmm5, xmm4	 ;# Vvdw6 
	mulss  xmm7, [esp + nb233_tsc]
	;# put scalar force on stack Update Vvdwtot directly 
	addss  xmm5, [esp + nb233_Vvdwtot]
	movss [esp + nb233_fstmp], xmm7
	movss [esp + nb233_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [esi + eax*4 + 16]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 24]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm7, [esp + nb233_two]	;# two*Heps2 
	addss  xmm7, xmm6
	addss  xmm7, xmm5 ;# xmm7=FF 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
 	
	movss xmm4, [esp + nb233_c12]
	mulss  xmm7, xmm4 ;# fijR 
	mulss  xmm5, xmm4 ;# Vvdw12 
	mulss  xmm7, [esp + nb233_tsc]
	addss  xmm7, [esp + nb233_fstmp]
	movss [esp + nb233_fstmp], xmm7
	addss  xmm5, [esp + nb233_Vvdwtot]
	movss [esp + nb233_Vvdwtot], xmm5	

	movd eax, mm0	

	movaps xmm4, xmm0
	movaps xmm1, xmm0
	movaps xmm5, xmm0
	mulps  xmm4, xmm4	;# xmm1=rinv, xmm4=rinvsq
	movaps xmm3, [esp + nb233_krsqM]
	addps  xmm5, xmm3	;# xmm0=rinv+ krsq 
	subps  xmm5, [esp + nb233_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm3, [esp + nb233_two]
	subps  xmm1, xmm3	;# xmm1=rinv-2*krsq
	movaps xmm2, xmm5
	mulps  xmm2, [esp + nb233_qqM]	;# xmm2=vcoul 
	mulps  xmm1, [esp + nb233_qqM] 	;# xmm1=coul part of fs 
	mulps  xmm1, xmm0
	subss  xmm1, [esp + nb233_fstmp]
	mulps  xmm1, xmm0
	movaps xmm4, xmm1
	
	addps  xmm2, [esp + nb233_vctot]	
	movaps [esp + nb233_vctot], xmm2

		
	movaps xmm0, [esp + nb233_dxO]
	movaps xmm1, [esp + nb233_dyO]
	movaps xmm2, [esp + nb233_dzO]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4 ;# xmm0-xmm2 now contains tx-tz (partial force)
	
	movss  xmm3, [esp + nb233_fixO]	
	movss  xmm4, [esp + nb233_fiyO]	
	movss  xmm5, [esp + nb233_fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esp + nb233_fixO], xmm3	
	movss  [esp + nb233_fiyO], xmm4	
	movss  [esp + nb233_fizO], xmm5	;# updated the O force now do the H's
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2      
	shufps xmm3, xmm3, 0x39	;# shift right 
	shufps xmm4, xmm4, 0x39
	shufps xmm5, xmm5, 0x39
	addss  xmm3, [esp + nb233_fixH1]
	addss  xmm4, [esp + nb233_fiyH1]
	addss  xmm5, [esp + nb233_fizH1]
	movss  [esp + nb233_fixH1], xmm3	
	movss  [esp + nb233_fiyH1], xmm4	
	movss  [esp + nb233_fizH1], xmm5	;# updated the H1 force 

	shufps xmm3, xmm3, 0x39
	shufps xmm4, xmm4, 0x39
	shufps xmm5, xmm5, 0x39
	addss  xmm3, [esp + nb233_fixH2]
	addss  xmm4, [esp + nb233_fiyH2]
	addss  xmm5, [esp + nb233_fizH2]
	movss  [esp + nb233_fixH2], xmm3	
	movss  [esp + nb233_fiyH2], xmm4	
	movss  [esp + nb233_fizH2], xmm5	;# updated the H2 force 

	mov edi, [ebp + nb233_faction]
	shufps xmm3, xmm3, 0x39
	shufps xmm4, xmm4, 0x39
	shufps xmm5, xmm5, 0x39
	addss  xmm3, [esp + nb233_fixM]
	addss  xmm4, [esp + nb233_fiyM]
	addss  xmm5, [esp + nb233_fizM]
	movss  [esp + nb233_fixM], xmm3	
	movss  [esp + nb233_fiyM], xmm4	
	movss  [esp + nb233_fizM], xmm5	;# updated the M force 

	;# the fj's - move in from mem start by acc. tx/ty/tz in xmm0, xmm1
	movlps xmm6, [edi + eax*4]
	movss  xmm7, [edi + eax*4 + 8]
	
	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  xmm0, xmm3
	movaps  xmm1, xmm4
	movaps  xmm2, xmm5
		
	shufps xmm3, xmm3, 0x39	;# shift right 
	shufps xmm4, xmm4, 0x39
	shufps xmm5, xmm5, 0x39
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5
	unpcklps xmm0, xmm1 	;# x,y sum in xmm0, z sum in xmm2
	
	subps    xmm6, xmm0
	subss    xmm7, xmm2
	
	movlps [edi + eax*4],     xmm6
	movss  [edi + eax*4 + 8], xmm7

	dec dword ptr [esp + nb233_innerk]
	jz    .nb233_updateouterdata
	jmp   .nb233_odd_loop
.nb233_updateouterdata:
	mov   ecx, [esp + nb233_ii3]
	mov   edi, [ebp + nb233_faction]
	mov   esi, [ebp + nb233_fshift]
	mov   edx, [esp + nb233_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb233_fixO]
	movaps xmm1, [esp + nb233_fiyO]
	movaps xmm2, [esp + nb233_fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# constant 00001000	

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb233_fixH1]
	movaps xmm1, [esp + nb233_fiyH1]
	movaps xmm2, [esp + nb233_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb233_fixH2]
	movaps xmm1, [esp + nb233_fiyH2]
	movaps xmm2, [esp + nb233_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# accumulate Mi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb233_fixM]
	movaps xmm1, [esp + nb233_fiyM]
	movaps xmm2, [esp + nb233_fizM]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4 + 36]
	movss  xmm4, [edi + ecx*4 + 40]
	movss  xmm5, [edi + ecx*4 + 44]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 36], xmm3
	movss  [edi + ecx*4 + 40], xmm4
	movss  [edi + ecx*4 + 44], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	;# get n from stack
	mov esi, [esp + nb233_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb233_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb233_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb233_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + nb233_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb233_Vvdw]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb233_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb233_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb233_n], esi
        jmp .nb233_outer
.nb233_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb233_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb233_end
        ;# non-zero, do one more workunit
        jmp   .nb233_threadloop
.nb233_end:
	emms

	mov eax, [esp + nb233_nouter]
	mov ebx, [esp + nb233_ninner]
	mov ecx, [ebp + nb233_outeriter]
	mov edx, [ebp + nb233_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb233_salign]
	add esp, eax
	add esp, 1004
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret
	

.globl nb_kernel233nf_ia32_sse
.globl _nb_kernel233nf_ia32_sse
nb_kernel233nf_ia32_sse:	
_nb_kernel233nf_ia32_sse:	
.equiv          nb233nf_p_nri,            8
.equiv          nb233nf_iinr,             12
.equiv          nb233nf_jindex,           16
.equiv          nb233nf_jjnr,             20
.equiv          nb233nf_shift,            24
.equiv          nb233nf_shiftvec,         28
.equiv          nb233nf_fshift,           32
.equiv          nb233nf_gid,              36
.equiv          nb233nf_pos,              40
.equiv          nb233nf_faction,          44
.equiv          nb233nf_charge,           48
.equiv          nb233nf_p_facel,          52
.equiv          nb233nf_argkrf,           56
.equiv          nb233nf_argcrf,           60
.equiv          nb233nf_Vc,               64
.equiv          nb233nf_type,             68
.equiv          nb233nf_p_ntype,          72
.equiv          nb233nf_vdwparam,         76
.equiv          nb233nf_Vvdw,             80
.equiv          nb233nf_p_tabscale,       84
.equiv          nb233nf_VFtab,            88
.equiv          nb233nf_invsqrta,         92
.equiv          nb233nf_dvda,             96
.equiv          nb233nf_p_gbtabscale,     100
.equiv          nb233nf_GBtab,            104
.equiv          nb233nf_p_nthreads,       108
.equiv          nb233nf_count,            112
.equiv          nb233nf_mtx,              116
.equiv          nb233nf_outeriter,        120
.equiv          nb233nf_inneriter,        124
.equiv          nb233nf_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb233nf_ixO,              0
.equiv          nb233nf_iyO,              16
.equiv          nb233nf_izO,              32
.equiv          nb233nf_ixH1,             48
.equiv          nb233nf_iyH1,             64
.equiv          nb233nf_izH1,             80
.equiv          nb233nf_ixH2,             96
.equiv          nb233nf_iyH2,             112
.equiv          nb233nf_izH2,             128
.equiv          nb233nf_ixM,              144
.equiv          nb233nf_iyM,              160
.equiv          nb233nf_izM,              176
.equiv          nb233nf_iqM,              192
.equiv          nb233nf_iqH,              208
.equiv          nb233nf_qqH,              224
.equiv          nb233nf_rinvH1,           240
.equiv          nb233nf_rinvH2,           256
.equiv          nb233nf_rinvM,            272
.equiv          nb233nf_c6,               288
.equiv          nb233nf_c12,              304
.equiv          nb233nf_tsc,              320
.equiv          nb233nf_krf,              336
.equiv          nb233nf_crf,              352
.equiv          nb233nf_krsqH1,           368
.equiv          nb233nf_krsqH2,           384
.equiv          nb233nf_krsqM,            400
.equiv          nb233nf_vctot,            416
.equiv          nb233nf_Vvdwtot,          432
.equiv          nb233nf_half,             448
.equiv          nb233nf_three,            464
.equiv          nb233nf_qqM,              480
.equiv          nb233nf_is3,              496
.equiv          nb233nf_ii3,              500
.equiv          nb233nf_ntia,             504
.equiv          nb233nf_innerjjnr,        508
.equiv          nb233nf_innerk,           512
.equiv          nb233nf_n,                516
.equiv          nb233nf_nn1,              520
.equiv          nb233nf_nri,              524
.equiv          nb233nf_nouter,           528
.equiv          nb233nf_ninner,           532
.equiv          nb233nf_salign,           536
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx
	push esi
	push edi
	sub esp, 540		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb233nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb233nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb233nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb233nf_nouter], eax
	mov [esp + nb233nf_ninner], eax

	mov eax, [ebp + nb233nf_p_tabscale]
	movss xmm3, [eax]
	shufps xmm3, xmm3, 0
	movaps [esp + nb233nf_tsc], xmm3

	mov esi, [ebp + nb233nf_argkrf]
	mov edi, [ebp + nb233nf_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	movaps [esp + nb233nf_krf], xmm5
	movaps [esp + nb233nf_crf], xmm6
	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb233nf_half], eax
	movss xmm1, [esp + nb233nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb233nf_half],  xmm1
	movaps [esp + nb233nf_three],  xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb233nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb233nf_charge]
	movss xmm4, [edx + ebx*4 + 4]	
	movss xmm3, [edx + ebx*4 + 12]	
	mov esi, [ebp + nb233nf_p_facel]
	movss xmm5, [esi]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [esp + nb233nf_iqM], xmm3
	movaps [esp + nb233nf_iqH], xmm4
	
	mov   edx, [ebp + nb233nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov edi, [ebp + nb233nf_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	mov   [esp + nb233nf_ntia], ecx		
.nb233nf_threadloop:
        mov   esi, [ebp + nb233nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb233nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb233nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb233nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb233nf_n], eax
        mov [esp + nb233nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb233nf_outerstart
        jmp .nb233nf_end
	
.nb233nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb233nf_nouter]
	mov [esp + nb233nf_nouter], ebx

.nb233nf_outer:
	mov   eax, [ebp + nb233nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]	;# ebx=3*is 
	mov   [esp + nb233nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb233nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb233nf_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]		;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2	
	movaps xmm6, xmm0
	movaps xmm7, xmm1
	
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb233nf_pos]	;# eax = base of pos[]  
	mov   [esp + nb233nf_ii3], ebx

	addss xmm3, [eax + ebx*4]  	;# ox
	addss xmm4, [eax + ebx*4 + 4]  ;# oy
	addss xmm5, [eax + ebx*4 + 8]  ;# oz
	addss xmm6, [eax + ebx*4 + 12] ;# h1x
	addss xmm7, [eax + ebx*4 + 16] ;# h1y
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	shufps xmm7, xmm7, 0
	movaps [esp + nb233nf_ixO], xmm3
	movaps [esp + nb233nf_iyO], xmm4
	movaps [esp + nb233nf_izO], xmm5
	movaps [esp + nb233nf_ixH1], xmm6
	movaps [esp + nb233nf_iyH1], xmm7

	movss xmm6, xmm2
	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm6, [eax + ebx*4 + 20] ;# h1z
	addss xmm0, [eax + ebx*4 + 24] ;# h2x
	addss xmm1, [eax + ebx*4 + 28] ;# h2y
	addss xmm2, [eax + ebx*4 + 32] ;# h2z
	addss xmm3, [eax + ebx*4 + 36] ;# mx
	addss xmm4, [eax + ebx*4 + 40] ;# my
	addss xmm5, [eax + ebx*4 + 44] ;# mz

	shufps xmm6, xmm6, 0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb233nf_izH1], xmm6
	movaps [esp + nb233nf_ixH2], xmm0
	movaps [esp + nb233nf_iyH2], xmm1
	movaps [esp + nb233nf_izH2], xmm2
	movaps [esp + nb233nf_ixM], xmm3
	movaps [esp + nb233nf_iyM], xmm4
	movaps [esp + nb233nf_izM], xmm5
	
	;# clear vctot 
	xorps xmm4, xmm4
	movaps [esp + nb233nf_vctot], xmm4
	movaps [esp + nb233nf_Vvdwtot], xmm4

	mov   eax, [ebp + nb233nf_jindex]
	mov   ecx, [eax + esi*4]	 	;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   esi, [ebp + nb233nf_pos]
	mov   eax, [ebp + nb233nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb233nf_innerjjnr], eax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb233nf_ninner]
	mov   [esp + nb233nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb233nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb233nf_unroll_loop
	jmp   .nb233nf_odd_inner
.nb233nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb233nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]     	;# eax-edx=jnr1-4 

	add dword ptr [esp + nb233nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb233nf_charge]	;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	 	;# and in xmm4 
	mulps  xmm3, [esp + nb233nf_iqM]
	mulps  xmm4, [esp + nb233nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [esp + nb233nf_qqM], xmm3
	movaps  [esp + nb233nf_qqH], xmm4
	
	mov esi, [ebp + nb233nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + nb233nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + nb233nf_ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm6, xmm7, 221  ;# constant 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + nb233nf_c6], xmm4
	movaps [esp + nb233nf_c12], xmm6

	mov esi, [ebp + nb233nf_pos]   	;# base of pos[] 

	lea   eax, [eax + eax*2] 	;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2] 	;# replace jnr with j3 
	lea   edx, [edx + edx*2]	

	;# move four coordinates to xmm0-xmm2 	
	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# constant 10001000
	
	shufps xmm0, xmm5, 136  ;# constant 10001000
	shufps xmm1, xmm5, 221  ;# constant 11011101		

	;# move ixO-izO to xmm4-xmm6 
	movaps xmm4, [esp + nb233nf_ixO]
	movaps xmm5, [esp + nb233nf_iyO]
	movaps xmm6, [esp + nb233nf_izO]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	;# rsqO in xmm7

	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [esp + nb233nf_ixH1]
	movaps xmm5, [esp + nb233nf_iyH1]
	movaps xmm6, [esp + nb233nf_izH1]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	;# rsqH1 in xmm6 
		
	;# move ixH2-izH2 to xmm3-xmm5  
	movaps xmm3, [esp + nb233nf_ixH2]
	movaps xmm4, [esp + nb233nf_iyH2]
	movaps xmm5, [esp + nb233nf_izH2]

	;# calc dr 
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# square it 
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	
	;# move ixM-izM to xmm2-xmm4  
	movaps xmm3, [esp + nb233nf_iyM]
	movaps xmm4, [esp + nb233nf_izM]
	subps  xmm3, xmm1
	subps  xmm4, xmm2
	movaps xmm2, [esp + nb233nf_ixM]
	subps  xmm2, xmm0	

	;# square it 
	mulps xmm2,xmm2
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	addps xmm4, xmm3
	addps xmm4, xmm2	
	;# rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 
	movaps xmm0, xmm4
	movaps xmm1, xmm5
	movaps xmm2, xmm6
	mulps  xmm0, [esp + nb233nf_krf]	
	mulps  xmm1, [esp + nb233nf_krf]	
	mulps  xmm2, [esp + nb233nf_krf]	
	movaps [esp + nb233nf_krsqM], xmm0
	movaps [esp + nb233nf_krsqH2], xmm1
	movaps [esp + nb233nf_krsqH1], xmm2

	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm0, [esp + nb233nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [esp + nb233nf_half]
	movaps  [esp + nb233nf_rinvH1], xmm0	;# rinvH1 

	;# rsqH2 - seed to xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm0, [esp + nb233nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [esp + nb233nf_half]
	movaps  [esp + nb233nf_rinvH2], xmm0	;# rinvH2 

	;# rsqM - seed to xmm2 
	rsqrtps xmm2, xmm4
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm0, [esp + nb233nf_three]
	mulps   xmm2, xmm4	;# rsq*lu*lu 
	subps   xmm0, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm0, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm0, [esp + nb233nf_half]
	movaps  [esp + nb233nf_rinvM], xmm0
	
	;# Do the O LJ-only interaction directly.	
	;# rsqO is in xmm7
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb233nf_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb233nf_half]
	movaps  xmm0, xmm4
	;# xmm0=rinvO
	
	mulps xmm7, xmm0
	mulps xmm7, [esp + nb233nf_tsc] ;# rtab
	
	movhlps xmm5, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm5	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps  xmm7, xmm6	
	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3
	pslld mm7, 3

	mov  esi, [ebp + nb233nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7
	
	;# dispersion 
	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ;# got half dispersion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ;# other half of dispersion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [esp + nb233nf_c6]
	mulps  xmm5, xmm4	 ;# Vvdw6 

	addps  xmm5, [esp + nb233nf_Vvdwtot]
	movaps [esp + nb233nf_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [esi + eax*4 + 16]
	movlps xmm7, [esi + ecx*4 + 16]
	movhps xmm5, [esi + ebx*4 + 16]
	movhps xmm7, [esi + edx*4 + 16] ;# got half repulsion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 24]
	movlps xmm3, [esi + ecx*4 + 24]
	movhps xmm7, [esi + ebx*4 + 24]
	movhps xmm3, [esi + edx*4 + 24] ;# other half of repulsion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 	
	movaps xmm4, [esp + nb233nf_c12]
	mulps  xmm5, xmm4 ;# Vvdw12 

	addps  xmm5, [esp + nb233nf_Vvdwtot]
	movaps [esp + nb233nf_Vvdwtot], xmm5

	;# Do H1 interaction
	movaps  xmm7, [esp + nb233nf_rinvH1]
	movaps  xmm4, xmm7	
	mulps   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq
	movaps xmm0, xmm7
	movaps xmm1, [esp + nb233nf_krsqH1]
	addps  xmm0, xmm1
	subps  xmm0, [esp + nb233nf_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm0, [esp + nb233nf_qqH]

	addps  xmm0, [esp + nb233nf_vctot]	
	movaps [esp + nb233nf_vctot], xmm0

	;# Done with H1, do H2 interactions
	movaps  xmm7, [esp + nb233nf_rinvH2]
	movaps  xmm4, xmm7	
	mulps   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq
	movaps xmm0, xmm7
	movaps xmm1, [esp + nb233nf_krsqH2]
	addps  xmm0, xmm1
	subps  xmm0, [esp + nb233nf_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm0, [esp + nb233nf_qqH]

	addps  xmm0, [esp + nb233nf_vctot]	
	movaps [esp + nb233nf_vctot], xmm0

	;# Done with H2, do M interactions
	movaps  xmm7, [esp + nb233nf_rinvM]
	movaps  xmm4, xmm7	
	mulps   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq
	movaps xmm0, xmm7
	movaps xmm1, [esp + nb233nf_krsqM]
	addps  xmm0, xmm1
	subps  xmm0, [esp + nb233nf_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm0, [esp + nb233nf_qqM]

	addps  xmm0, [esp + nb233nf_vctot]	
	movaps [esp + nb233nf_vctot], xmm0

	;# should we do one more iteration? 
	sub dword ptr [esp + nb233nf_innerk],  4
	jl    .nb233nf_odd_inner
	jmp   .nb233nf_unroll_loop
.nb233nf_odd_inner:	
	add dword ptr [esp + nb233nf_innerk],  4
	jnz   .nb233nf_odd_loop
	jmp   .nb233nf_updateouterdata
.nb233nf_odd_loop:
	mov   edx, [esp + nb233nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb233nf_innerjjnr],  4	

 	xorps xmm4, xmm4  	;# clear reg.
	movss xmm4, [esp + nb233nf_iqM]
	mov esi, [ebp + nb233nf_charge] 
	movhps xmm4, [esp + nb233nf_iqH]  ;# [qM  0  qH  qH] 
	shufps xmm4, xmm4, 41	;# [0 qH qH qM]

	movss xmm3, [esi + eax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [esp + nb233nf_qqM], xmm3	;# use dummy qq for storage 
	
	xorps xmm6, xmm6
	mov esi, [ebp + nb233nf_type]
	mov ebx, [esi + eax*4]
	mov esi, [ebp + nb233nf_vdwparam]
	shl ebx, 1	
	add ebx, [esp + nb233nf_ntia]
	movlps xmm6, [esi + ebx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# constant 11111100
	shufps xmm7, xmm7, 253  ;# constant 11111101
	movaps [esp + nb233nf_c6], xmm6
	movaps [esp + nb233nf_c12], xmm7
	
	mov esi, [ebp + nb233nf_pos]
	lea eax, [eax + eax*2]  

	movss xmm3, [esp + nb233nf_ixO]
	movss xmm4, [esp + nb233nf_iyO]
	movss xmm5, [esp + nb233nf_izO]
	movss xmm0, [esp + nb233nf_ixH1]
	movss xmm1, [esp + nb233nf_iyH1]
	movss xmm2, [esp + nb233nf_izH1]
	unpcklps xmm3, [esp + nb233nf_ixH2] 	;# ixO ixH2 - -
	unpcklps xmm4, [esp + nb233nf_iyH2]  	;# iyO iyH2 - -
	unpcklps xmm5, [esp + nb233nf_izH2]	;# izO izH2 - -
	unpcklps xmm0, [esp + nb233nf_ixM] 	;# ixH1 ixM - -
	unpcklps xmm1, [esp + nb233nf_iyM]  	;# iyH1 iyM - -
	unpcklps xmm2, [esp + nb233nf_izM]	;# izH1 izM - -
	unpcklps xmm3, xmm0  	;# ixO ixH1 ixH2 ixM
	unpcklps xmm4, xmm1 	;# same for y
	unpcklps xmm5, xmm2 	;# same for z
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 
	movaps xmm0, xmm4
	mulps xmm0, [esp + nb233nf_krf]
	movaps [esp + nb233nf_krsqM], xmm0
		
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb233nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb233nf_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv, xmm4=rsq

	;# LJ table interaction
	mulps xmm4, xmm0
	mulps  xmm4, [esp + nb233nf_tsc] ;# rtab
	
	cvttps2pi mm6, xmm4
	cvtpi2ps xmm6, mm6
	subss  xmm4, xmm6	
	movss xmm1, xmm4	;# xmm1=eps 
	movss xmm2, xmm1	
	mulss  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3

	mov  esi, [ebp + nb233nf_VFtab]
	movd eax, mm6
	
	;# dispersion 
	movlps xmm5, [esi + eax*4]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [esi + eax*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 

	movss xmm4, [esp + nb233nf_c6]
	mulss  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addss  xmm5, [esp + nb233nf_Vvdwtot]
	movss [esp + nb233nf_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [esi + eax*4 + 16]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [esi + eax*4 + 24]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
 	
	movss xmm4, [esp + nb233nf_c12]
	mulss  xmm5, xmm4 ;# Vvdw12 

	addss  xmm5, [esp + nb233nf_Vvdwtot]
	movss [esp + nb233nf_Vvdwtot], xmm5	

	movaps xmm4, xmm0
	movaps xmm1, xmm0
	movaps xmm5, xmm0
	mulps  xmm4, xmm4	;# xmm1=rinv, xmm4=rinvsq
	movaps xmm3, [esp + nb233nf_krsqM]
	addps  xmm5, xmm3	;# xmm0=rinv+ krsq 
	subps  xmm5, [esp + nb233nf_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm5, [esp + nb233nf_qqM]	;# xmm2=vcoul 
	
	addps  xmm5, [esp + nb233nf_vctot]	
	movaps [esp + nb233nf_vctot], xmm5

	dec dword ptr [esp + nb233nf_innerk]
	jz    .nb233nf_updateouterdata
	jmp   .nb233nf_odd_loop
.nb233nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb233nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb233nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb233nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb233nf_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + nb233nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb233nf_Vvdw]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb233nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb233nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb233nf_n], esi
        jmp .nb233nf_outer
.nb233nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb233nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb233nf_end
        ;# non-zero, do one more workunit
        jmp   .nb233nf_threadloop
.nb233nf_end:
	emms

	mov eax, [esp + nb233nf_nouter]
	mov ebx, [esp + nb233nf_ninner]
	mov ecx, [ebp + nb233nf_outeriter]
	mov edx, [ebp + nb233nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb233nf_salign]
	add esp, eax
	add esp, 540
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret
	
