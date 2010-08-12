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



.globl nb_kernel203_ia32_sse
.globl _nb_kernel203_ia32_sse
nb_kernel203_ia32_sse:	
_nb_kernel203_ia32_sse:	
.equiv          nb203_p_nri,            8
.equiv          nb203_iinr,             12
.equiv          nb203_jindex,           16
.equiv          nb203_jjnr,             20
.equiv          nb203_shift,            24
.equiv          nb203_shiftvec,         28
.equiv          nb203_fshift,           32
.equiv          nb203_gid,              36
.equiv          nb203_pos,              40
.equiv          nb203_faction,          44
.equiv          nb203_charge,           48
.equiv          nb203_p_facel,          52
.equiv          nb203_argkrf,           56
.equiv          nb203_argcrf,           60
.equiv          nb203_Vc,               64
.equiv          nb203_type,             68
.equiv          nb203_p_ntype,          72
.equiv          nb203_vdwparam,         76
.equiv          nb203_Vvdw,             80
.equiv          nb203_p_tabscale,       84
.equiv          nb203_VFtab,            88
.equiv          nb203_invsqrta,         92
.equiv          nb203_dvda,             96
.equiv          nb203_p_gbtabscale,     100
.equiv          nb203_GBtab,            104
.equiv          nb203_p_nthreads,       108
.equiv          nb203_count,            112
.equiv          nb203_mtx,              116
.equiv          nb203_outeriter,        120
.equiv          nb203_inneriter,        124
.equiv          nb203_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb203_ixH1,             0
.equiv          nb203_iyH1,             16
.equiv          nb203_izH1,             32
.equiv          nb203_ixH2,             48
.equiv          nb203_iyH2,             64
.equiv          nb203_izH2,             80
.equiv          nb203_ixM,              96
.equiv          nb203_iyM,              112
.equiv          nb203_izM,              128
.equiv          nb203_iqH,              144
.equiv          nb203_iqM,              160
.equiv          nb203_dxH1,             176
.equiv          nb203_dyH1,             192
.equiv          nb203_dzH1,             208
.equiv          nb203_dxH2,             224
.equiv          nb203_dyH2,             240
.equiv          nb203_dzH2,             256
.equiv          nb203_dxM,              272
.equiv          nb203_dyM,              288
.equiv          nb203_dzM,              304
.equiv          nb203_qqH,              320
.equiv          nb203_qqM,              336
.equiv          nb203_vctot,            352
.equiv          nb203_fixH1,            384
.equiv          nb203_fiyH1,            400
.equiv          nb203_fizH1,            416
.equiv          nb203_fixH2,            432
.equiv          nb203_fiyH2,            448
.equiv          nb203_fizH2,            464
.equiv          nb203_fixM,             480
.equiv          nb203_fiyM,             496
.equiv          nb203_fizM,             512
.equiv          nb203_fjx,              528
.equiv          nb203_fjy,              544
.equiv          nb203_fjz,              560
.equiv          nb203_half,             576
.equiv          nb203_three,            592
.equiv          nb203_two,              608
.equiv          nb203_krf,              624
.equiv          nb203_crf,              640
.equiv          nb203_krsqH1,           656
.equiv          nb203_krsqH2,           672
.equiv          nb203_krsqM,            688
.equiv          nb203_is3,              704
.equiv          nb203_ii3,              708
.equiv          nb203_innerjjnr,        712
.equiv          nb203_innerk,           716
.equiv          nb203_n,                720
.equiv          nb203_nn1,              724
.equiv          nb203_nri,              728
.equiv          nb203_nouter,           732
.equiv          nb203_ninner,           736
.equiv          nb203_salign,           740
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx
	push esi
	push edi
	sub esp, 744		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb203_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb203_p_nri]
	mov ecx, [ecx]
	mov [esp + nb203_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb203_nouter], eax
	mov [esp + nb203_ninner], eax


	mov esi, [ebp + nb203_argkrf]
	mov edi, [ebp + nb203_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	movaps [esp + nb203_krf], xmm5
	movaps [esp + nb203_crf], xmm6
	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb203_half], eax
	movss xmm1, [esp + nb203_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb203_half],  xmm1
	movaps [esp + nb203_two],  xmm2
	movaps [esp + nb203_three],  xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb203_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]		;# ebx =ii 

	mov   edx, [ebp + nb203_charge]
	movss xmm3, [edx + ebx*4 + 4]	
	movss xmm4, [edx + ebx*4 + 12]	
	mov esi, [ebp + nb203_p_facel]
	movss xmm5, [esi]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [esp + nb203_iqH], xmm3
	movaps [esp + nb203_iqM], xmm4
			
.nb203_threadloop:
        mov   esi, [ebp + nb203_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb203_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb203_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb203_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb203_n], eax
        mov [esp + nb203_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb203_outerstart
        jmp .nb203_end
	
.nb203_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb203_nouter]
	mov [esp + nb203_nouter], ebx

.nb203_outer:
	mov   eax, [ebp + nb203_shift]  	;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]	;# ebx=3*is 
	mov   [esp + nb203_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb203_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb203_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]		;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb203_pos]	;# eax = base of pos[]  
	mov   [esp + nb203_ii3], ebx

	addss xmm3, [eax + ebx*4 + 12]
	addss xmm4, [eax + ebx*4 + 16]
	addss xmm5, [eax + ebx*4 + 20]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb203_ixH1], xmm3
	movaps [esp + nb203_iyH1], xmm4
	movaps [esp + nb203_izH1], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 24]
	addss xmm1, [eax + ebx*4 + 28]
	addss xmm2, [eax + ebx*4 + 32]		
	addss xmm3, [eax + ebx*4 + 36]
	addss xmm4, [eax + ebx*4 + 40]
	addss xmm5, [eax + ebx*4 + 44]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb203_ixH2], xmm0
	movaps [esp + nb203_iyH2], xmm1
	movaps [esp + nb203_izH2], xmm2
	movaps [esp + nb203_ixM], xmm3
	movaps [esp + nb203_iyM], xmm4
	movaps [esp + nb203_izM], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb203_vctot], xmm4
	movaps [esp + nb203_fixH1], xmm4
	movaps [esp + nb203_fiyH1], xmm4
	movaps [esp + nb203_fizH1], xmm4
	movaps [esp + nb203_fixH2], xmm4
	movaps [esp + nb203_fiyH2], xmm4
	movaps [esp + nb203_fizH2], xmm4
	movaps [esp + nb203_fixM], xmm4
	movaps [esp + nb203_fiyM], xmm4
	movaps [esp + nb203_fizM], xmm4
	
	mov   eax, [ebp + nb203_jindex]
	mov   ecx, [eax + esi*4]	 	;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   esi, [ebp + nb203_pos]
	mov   edi, [ebp + nb203_faction]	
	mov   eax, [ebp + nb203_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb203_innerjjnr], eax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb203_ninner]
	mov   [esp + nb203_ninner], ecx
	add   edx, 0
	mov   [esp + nb203_innerk], edx	;# number of innerloop atoms 
	jge   .nb203_unroll_loop
	jmp   .nb203_odd_inner
.nb203_unroll_loop:
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb203_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]     	;# eax-edx=jnr1-4 

	add dword ptr [esp + nb203_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb203_charge]	;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	 	;# and in xmm4 
	mulps  xmm3, [esp + nb203_iqH]
	mulps  xmm4, [esp + nb203_iqM]

	movaps  [esp + nb203_qqH], xmm3
	movaps  [esp + nb203_qqM], xmm4

	mov esi, [ebp + nb203_pos]   	;# base of pos[] 

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

	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [esp + nb203_ixH1]
	movaps xmm5, [esp + nb203_iyH1]
	movaps xmm6, [esp + nb203_izH1]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb203_dxH1], xmm4
	movaps [esp + nb203_dyH1], xmm5
	movaps [esp + nb203_dzH1], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	;# rsqH1 in xmm7 

	;# move ixH2-izH2 to xmm4-xmm6 
	movaps xmm4, [esp + nb203_ixH2]
	movaps xmm5, [esp + nb203_iyH2]
	movaps xmm6, [esp + nb203_izH2]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb203_dxH2], xmm4
	movaps [esp + nb203_dyH2], xmm5
	movaps [esp + nb203_dzH2], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	;# rsqH2 in xmm6 

	;# move ixM-izM to xmm3-xmm5  
	movaps xmm3, [esp + nb203_ixM]
	movaps xmm4, [esp + nb203_iyM]
	movaps xmm5, [esp + nb203_izM]

	;# calc dr 
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# store dr 
	movaps [esp + nb203_dxM], xmm3
	movaps [esp + nb203_dyM], xmm4
	movaps [esp + nb203_dzM], xmm5
	;# square it 
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	;# rsqM in xmm5, rsqH2 in xmm6, rsqH1 in xmm7 

	movaps xmm0, xmm5
	movaps xmm1, xmm6
	movaps xmm2, xmm7

	mulps  xmm0, [esp + nb203_krf]	
	mulps  xmm1, [esp + nb203_krf]	
	mulps  xmm2, [esp + nb203_krf]	

	movaps [esp + nb203_krsqM], xmm0
	movaps [esp + nb203_krsqH2], xmm1
	movaps [esp + nb203_krsqH1], xmm2
	
	;# start with rsqH1 - seed in xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb203_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb203_half]
	movaps  xmm7, xmm4	;# rinvH1 in xmm7 
	;# rsqH2 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb203_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb203_half]
	movaps  xmm6, xmm4	;# rinvH2 in xmm6 
	;# rsqM - seed in xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb203_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb203_half]
	movaps  xmm5, xmm4	;# rinvM in xmm5 

	;# do H1 interactions 
	movaps  xmm4, xmm7	
	mulps   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq 

	movaps xmm0, xmm7
	movaps xmm1, [esp + nb203_krsqH1]
	addps  xmm0, xmm1
	subps  xmm0, [esp + nb203_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm1, [esp + nb203_two]
	subps  xmm7, xmm1
	mulps  xmm0, [esp + nb203_qqH]
	mulps  xmm7, [esp + nb203_qqH]

	mulps  xmm4, xmm7	;# total fs H1 in xmm4 

	addps  xmm0, [esp + nb203_vctot]
	movaps [esp + nb203_vctot], xmm0

	movaps xmm0, [esp + nb203_dxH1]
	movaps xmm1, [esp + nb203_dyH1]
	movaps xmm2, [esp + nb203_dzH1]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	;# update H1 forces 
	movaps xmm3, [esp + nb203_fixH1]
	movaps xmm4, [esp + nb203_fiyH1]
	movaps xmm7, [esp + nb203_fizH1]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb203_fixH1], xmm3
	movaps [esp + nb203_fiyH1], xmm4
	movaps [esp + nb203_fizH1], xmm7
	;# update j forces with water O 
	movaps [esp + nb203_fjx], xmm0
	movaps [esp + nb203_fjy], xmm1
	movaps [esp + nb203_fjz], xmm2

	;# H2 interactions 
	movaps  xmm4, xmm6	
	mulps   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movaps  xmm7, xmm6
	movaps  xmm0, [esp + nb203_krsqH2]
	addps   xmm6, xmm0	;# xmm6=rinv+ krsq 
	subps   xmm6, [esp + nb203_crf] ;# xmm6=rinv+ krsq-crf 
	mulps   xmm0, [esp + nb203_two]
	subps   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulps   xmm6, [esp + nb203_qqH] ;# vcoul 
	mulps   xmm7, [esp + nb203_qqH]
	mulps  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addps  xmm6, [esp + nb203_vctot]

	movaps xmm0, [esp + nb203_dxH2]
	movaps xmm1, [esp + nb203_dyH2]
	movaps xmm2, [esp + nb203_dzH2]
	movaps [esp + nb203_vctot], xmm6
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	;# update H2 forces 
	movaps xmm3, [esp + nb203_fixH2]
	movaps xmm4, [esp + nb203_fiyH2]
	movaps xmm7, [esp + nb203_fizH2]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb203_fixH2], xmm3
	movaps [esp + nb203_fiyH2], xmm4
	movaps [esp + nb203_fizH2], xmm7
	;# update j forces with water H2
	addps  xmm0, [esp + nb203_fjx]
	addps  xmm1, [esp + nb203_fjy]
	addps  xmm2, [esp + nb203_fjz]
	movaps [esp + nb203_fjx], xmm0
	movaps [esp + nb203_fjy], xmm1
	movaps [esp + nb203_fjz], xmm2

	;# M interactions 
	movaps  xmm4, xmm5	
	mulps   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movaps  xmm7, xmm5
	movaps  xmm0, [esp + nb203_krsqM]
	addps   xmm5, xmm0	;# xmm6=rinv+ krsq 
	subps   xmm5, [esp + nb203_crf] ;# xmm5=rinv+ krsq-crf 
	mulps   xmm0, [esp + nb203_two]
	subps   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulps   xmm5, [esp + nb203_qqM] ;# vcoul 
	mulps   xmm7, [esp + nb203_qqM]
	mulps  xmm4, xmm7		;# total fsM in xmm4 
	
	addps  xmm5, [esp + nb203_vctot]

	movaps xmm0, [esp + nb203_dxM]
	movaps xmm1, [esp + nb203_dyM]
	movaps xmm2, [esp + nb203_dzM]
	movaps [esp + nb203_vctot], xmm5
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	;# update M forces 
	movaps xmm3, [esp + nb203_fixM]
	movaps xmm4, [esp + nb203_fiyM]
	movaps xmm7, [esp + nb203_fizM]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb203_fixM], xmm3
	movaps [esp + nb203_fiyM], xmm4
	movaps [esp + nb203_fizM], xmm7

	mov edi, [ebp + nb203_faction]
	
	;# update j forces 
	addps xmm0, [esp + nb203_fjx]
	addps xmm1, [esp + nb203_fjy]
	addps xmm2, [esp + nb203_fjz]

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
	sub dword ptr [esp + nb203_innerk],  4
	jl    .nb203_odd_inner
	jmp   .nb203_unroll_loop
.nb203_odd_inner:	
	add dword ptr [esp + nb203_innerk],  4
	jnz   .nb203_odd_loop
	jmp   .nb203_updateouterdata
.nb203_odd_loop:
	mov   edx, [esp + nb203_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb203_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + nb203_iqM]
	mov esi, [ebp + nb203_charge] 
	movhps xmm4, [esp + nb203_iqH]     
	movss xmm3, [esi + eax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [esp + nb203_qqM], xmm3	;# use dummy qq for storage 

	mov esi, [ebp + nb203_pos]
	lea   eax, [eax + eax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [esp + nb203_ixM]
	movss xmm4, [esp + nb203_iyM]
	movss xmm5, [esp + nb203_izM]
		
	movlps xmm6, [esp + nb203_ixH1]
	movlps xmm7, [esp + nb203_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + nb203_iyH1]
	movlps xmm7, [esp + nb203_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + nb203_izH1]
	movlps xmm7, [esp + nb203_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	;# use dummy dx for storage
	movaps [esp + nb203_dxM], xmm3
	movaps [esp + nb203_dyM], xmm4
	movaps [esp + nb203_dzM], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 

	movaps xmm0, xmm4
	mulps xmm0, [esp + nb203_krf]
	movaps [esp + nb203_krsqM], xmm0
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb203_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb203_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# constant 11100000

	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 

	movaps xmm1, xmm0	;# xmm1=rinv 
	movaps xmm3, [esp + nb203_krsqM]
	addps  xmm0, xmm3	;# xmm0=rinv+ krsq 
	subps  xmm0, [esp + nb203_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm3, [esp + nb203_two]
	subps  xmm1, xmm3	;# xmm1=rinv-2*krsq 
	mulps  xmm0, [esp + nb203_qqM]	;# xmm0=vcoul 
	mulps  xmm1, [esp + nb203_qqM] 	;# xmm1=coul part of fs 

	
	mulps  xmm4, xmm1	;# xmm4=total fscal 
	addps  xmm0, [esp + nb203_vctot]
	movaps [esp + nb203_vctot], xmm0
	
	movaps xmm0, [esp + nb203_dxM]
	movaps xmm1, [esp + nb203_dyM]
	movaps xmm2, [esp + nb203_dzM]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	movss  xmm3, [esp + nb203_fixM]	
	movss  xmm4, [esp + nb203_fiyM]	
	movss  xmm5, [esp + nb203_fizM]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esp + nb203_fixM], xmm3	
	movss  [esp + nb203_fiyM], xmm4	
	movss  [esp + nb203_fizM], xmm5	;# updated the O force now do the H's 
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 230 ;# constant 11100110	;# shift right 
	shufps xmm4, xmm4, 230 ;# constant 11100110
	shufps xmm5, xmm5, 230 ;# constant 11100110
	addss  xmm3, [esp + nb203_fixH1]
	addss  xmm4, [esp + nb203_fiyH1]
	addss  xmm5, [esp + nb203_fizH1]
	movss  [esp + nb203_fixH1], xmm3	
	movss  [esp + nb203_fiyH1], xmm4	
	movss  [esp + nb203_fizH1], xmm5	;# updated the H1 force 

	mov edi, [ebp + nb203_faction]
	shufps xmm3, xmm3, 231 ;# constant 11100111	;# shift right 
	shufps xmm4, xmm4, 231 ;# constant 11100111
	shufps xmm5, xmm5, 231 ;# constant 11100111
	addss  xmm3, [esp + nb203_fixH2]
	addss  xmm4, [esp + nb203_fiyH2]
	addss  xmm5, [esp + nb203_fizH2]
	movss  [esp + nb203_fixH2], xmm3	
	movss  [esp + nb203_fiyH2], xmm4	
	movss  [esp + nb203_fizH2], xmm5	;# updated the H2 force 

	;# the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1 
	xorps  xmm5, xmm5
	movaps xmm3, xmm0
	movlps xmm6, [edi + eax*4]
	movss  xmm7, [edi + eax*4 + 8]
	unpcklps xmm3, xmm1
	movlhps  xmm3, xmm5	
	unpckhps xmm0, xmm1		
	addps    xmm0, xmm3
	movhlps  xmm3, xmm0	
	addps    xmm0, xmm3	;# x,y sum in xmm0 

	movhlps  xmm1, xmm2
	addss    xmm2, xmm1
	shufps   xmm1, xmm1, 1 
	addss    xmm2, xmm1	;# z sum in xmm2 
	subps    xmm6, xmm0
	subss    xmm7, xmm2
	
	movlps [edi + eax*4],     xmm6
	movss  [edi + eax*4 + 8], xmm7

	dec dword ptr [esp + nb203_innerk]
	jz    .nb203_updateouterdata
	jmp   .nb203_odd_loop
.nb203_updateouterdata:
	mov   ecx, [esp + nb203_ii3]
	mov   edi, [ebp + nb203_faction]
	mov   esi, [ebp + nb203_fshift]
	mov   edx, [esp + nb203_is3]

	;# accumulate  H1 i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb203_fixH1]
	movaps xmm1, [esp + nb203_fiyH1]
	movaps xmm2, [esp + nb203_fizH1]

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
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# constant 00001000	

	;# accumulate H2 i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb203_fixH2]
	movaps xmm1, [esp + nb203_fiyH2]
	movaps xmm2, [esp + nb203_fizH2]

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

	;# accumulate m i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb203_fixM]
	movaps xmm1, [esp + nb203_fiyM]
	movaps xmm2, [esp + nb203_fizM]

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
	mov esi, [esp + nb203_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb203_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb203_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb203_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb203_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb203_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb203_n], esi
        jmp .nb203_outer
.nb203_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb203_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb203_end
        ;# non-zero, do one more workunit
        jmp   .nb203_threadloop
.nb203_end:
	emms

	mov eax, [esp + nb203_nouter]
	mov ebx, [esp + nb203_ninner]
	mov ecx, [ebp + nb203_outeriter]
	mov edx, [ebp + nb203_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb203_salign]
	add esp, eax
	add esp, 744
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret




.globl nb_kernel203nf_ia32_sse
.globl _nb_kernel203nf_ia32_sse
nb_kernel203nf_ia32_sse:	
_nb_kernel203nf_ia32_sse:	
.equiv          nb203nf_p_nri,          8
.equiv          nb203nf_iinr,           12
.equiv          nb203nf_jindex,         16
.equiv          nb203nf_jjnr,           20
.equiv          nb203nf_shift,          24
.equiv          nb203nf_shiftvec,       28
.equiv          nb203nf_fshift,         32
.equiv          nb203nf_gid,            36
.equiv          nb203nf_pos,            40
.equiv          nb203nf_faction,        44
.equiv          nb203nf_charge,         48
.equiv          nb203nf_p_facel,        52
.equiv          nb203nf_argkrf,         56
.equiv          nb203nf_argcrf,         60
.equiv          nb203nf_Vc,             64
.equiv          nb203nf_type,           68
.equiv          nb203nf_p_ntype,        72
.equiv          nb203nf_vdwparam,       76
.equiv          nb203nf_Vvdw,           80
.equiv          nb203nf_p_tabscale,     84
.equiv          nb203nf_VFtab,          88
.equiv          nb203nf_invsqrta,       92
.equiv          nb203nf_dvda,           96
.equiv          nb203nf_p_gbtabscale,   100
.equiv          nb203nf_GBtab,          104
.equiv          nb203nf_p_nthreads,     108
.equiv          nb203nf_count,          112
.equiv          nb203nf_mtx,            116
.equiv          nb203nf_outeriter,      120
.equiv          nb203nf_inneriter,      124
.equiv          nb203nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb203nf_ixH1,           0
.equiv          nb203nf_iyH1,           16
.equiv          nb203nf_izH1,           32
.equiv          nb203nf_ixH2,           48
.equiv          nb203nf_iyH2,           64
.equiv          nb203nf_izH2,           80
.equiv          nb203nf_ixM,            96
.equiv          nb203nf_iyM,            112
.equiv          nb203nf_izM,            128
.equiv          nb203nf_iqH,            144
.equiv          nb203nf_iqM,            160
.equiv          nb203nf_qqH,            176
.equiv          nb203nf_qqM,            192
.equiv          nb203nf_vctot,          208
.equiv          nb203nf_half,           224
.equiv          nb203nf_three,          240
.equiv          nb203nf_krf,            256
.equiv          nb203nf_crf,            272
.equiv          nb203nf_krsqH1,         288
.equiv          nb203nf_krsqH2,         304
.equiv          nb203nf_krsqM,          320
.equiv          nb203nf_is3,            336
.equiv          nb203nf_ii3,            340
.equiv          nb203nf_innerjjnr,      344
.equiv          nb203nf_innerk,         348
.equiv          nb203nf_n,              352
.equiv          nb203nf_nn1,            356
.equiv          nb203nf_nri,            360
.equiv          nb203nf_nouter,         364
.equiv          nb203nf_ninner,         368
.equiv          nb203nf_salign,         372
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx
	push esi
	push edi
	sub esp, 376		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb203nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb203nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb203nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb203nf_nouter], eax
	mov [esp + nb203nf_ninner], eax


	mov esi, [ebp + nb203nf_argkrf]
	mov edi, [ebp + nb203nf_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	movaps [esp + nb203nf_krf], xmm5
	movaps [esp + nb203nf_crf], xmm6
	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb203nf_half], eax
	movss xmm1, [esp + nb203nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb203nf_half],  xmm1
	movaps [esp + nb203nf_three],  xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb203nf_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]		;# ebx =ii 

	mov   edx, [ebp + nb203nf_charge]
	movss xmm3, [edx + ebx*4 + 4]	
	movss xmm4, [edx + ebx*4 + 12]	
	mov esi, [ebp + nb203nf_p_facel]
	movss xmm5, [esi]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [esp + nb203nf_iqH], xmm3
	movaps [esp + nb203nf_iqM], xmm4
			
.nb203nf_threadloop:
        mov   esi, [ebp + nb203nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb203nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb203nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb203nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb203nf_n], eax
        mov [esp + nb203nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb203nf_outerstart
        jmp .nb203nf_end
	
.nb203nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb203nf_nouter]
	mov [esp + nb203nf_nouter], ebx

.nb203nf_outer:
	mov   eax, [ebp + nb203nf_shift]  	;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]	;# ebx=3*is 
	mov   [esp + nb203nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb203nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb203nf_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]		;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb203nf_pos]	;# eax = base of pos[]  
	mov   [esp + nb203nf_ii3], ebx

	addss xmm3, [eax + ebx*4 + 12]
	addss xmm4, [eax + ebx*4 + 16]
	addss xmm5, [eax + ebx*4 + 20]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb203nf_ixH1], xmm3
	movaps [esp + nb203nf_iyH1], xmm4
	movaps [esp + nb203nf_izH1], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 24]
	addss xmm1, [eax + ebx*4 + 28]
	addss xmm2, [eax + ebx*4 + 32]		
	addss xmm3, [eax + ebx*4 + 36]
	addss xmm4, [eax + ebx*4 + 40]
	addss xmm5, [eax + ebx*4 + 44]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb203nf_ixH2], xmm0
	movaps [esp + nb203nf_iyH2], xmm1
	movaps [esp + nb203nf_izH2], xmm2
	movaps [esp + nb203nf_ixM], xmm3
	movaps [esp + nb203nf_iyM], xmm4
	movaps [esp + nb203nf_izM], xmm5
	
	;# clear vctot 
	xorps xmm4, xmm4
	movaps [esp + nb203nf_vctot], xmm4
	
	mov   eax, [ebp + nb203nf_jindex]
	mov   ecx, [eax + esi*4]	 	;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   esi, [ebp + nb203nf_pos]
	mov   edi, [ebp + nb203nf_faction]	
	mov   eax, [ebp + nb203nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb203nf_innerjjnr], eax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb203nf_ninner]
	mov   [esp + nb203nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb203nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb203nf_unroll_loop
	jmp   .nb203nf_odd_inner
.nb203nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb203nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]     	;# eax-edx=jnr1-4 

	add dword ptr [esp + nb203nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb203nf_charge]	;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	 	;# and in xmm4 
	mulps  xmm3, [esp + nb203nf_iqH]
	mulps  xmm4, [esp + nb203nf_iqM]

	movaps  [esp + nb203nf_qqH], xmm3
	movaps  [esp + nb203nf_qqM], xmm4

	mov esi, [ebp + nb203nf_pos]   	;# base of pos[] 

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

	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [esp + nb203nf_ixH1]
	movaps xmm5, [esp + nb203nf_iyH1]
	movaps xmm6, [esp + nb203nf_izH1]

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
	;# rsqH1 in xmm7 

	;# move ixH2-izH2 to xmm4-xmm6 
	movaps xmm4, [esp + nb203nf_ixH2]
	movaps xmm5, [esp + nb203nf_iyH2]
	movaps xmm6, [esp + nb203nf_izH2]

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
	;# rsqH2 in xmm6 

	;# move ixM-izM to xmm3-xmm5  
	movaps xmm3, [esp + nb203nf_ixM]
	movaps xmm4, [esp + nb203nf_iyM]
	movaps xmm5, [esp + nb203nf_izM]

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
	;# rsqM in xmm5, rsqH2 in xmm6, rsqH1 in xmm7 

	movaps xmm0, xmm5
	movaps xmm1, xmm6
	movaps xmm2, xmm7

	mulps  xmm0, [esp + nb203nf_krf]	
	mulps  xmm1, [esp + nb203nf_krf]	
	mulps  xmm2, [esp + nb203nf_krf]	

	movaps [esp + nb203nf_krsqM], xmm0
	movaps [esp + nb203nf_krsqH2], xmm1
	movaps [esp + nb203nf_krsqH1], xmm2
	
	;# start with rsqH1 - seed in xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb203nf_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb203nf_half]
	movaps  xmm7, xmm4	;# rinvH1 in xmm7 
	;# rsqH2 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb203nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb203nf_half]
	movaps  xmm6, xmm4	;# rinvH2 in xmm6 
	;# rsqM - seed in xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb203nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb203nf_half]
	movaps  xmm5, xmm4	;# rinvM in xmm5 

	;# do H1 interactions - xmm7=rinv
	addps xmm7, [esp + nb203nf_krsqH1]
	subps xmm7, [esp + nb203nf_crf] ;# xmm7=rinv+ krsq-crf 
	mulps xmm7, [esp + nb203nf_qqH]
	addps xmm7, [esp + nb203nf_vctot]

	;# H2 interactions - xmm6=rinv
	addps xmm6, [esp + nb203nf_krsqH2]
	subps xmm6, [esp + nb203nf_crf] ;# xmm6=rinv+ krsq-crf 
	mulps xmm6, [esp + nb203nf_qqH]
	addps xmm6, xmm7

	;# M interactions - xmm5=rinv
	addps xmm5, [esp + nb203nf_krsqM]
	subps xmm5, [esp + nb203nf_crf] ;# xmm5=rinv+ krsq-crf 
	mulps xmm5, [esp + nb203nf_qqM]
	addps xmm5, xmm6
	movaps [esp + nb203nf_vctot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb203nf_innerk],  4
	jl    .nb203nf_odd_inner
	jmp   .nb203nf_unroll_loop
.nb203nf_odd_inner:	
	add dword ptr [esp + nb203nf_innerk],  4
	jnz   .nb203nf_odd_loop
	jmp   .nb203nf_updateouterdata
.nb203nf_odd_loop:
	mov   edx, [esp + nb203nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb203nf_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + nb203nf_iqM]
	mov esi, [ebp + nb203nf_charge] 
	movhps xmm4, [esp + nb203nf_iqH]     
	movss xmm3, [esi + eax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [esp + nb203nf_qqM], xmm3	;# use dummy qq for storage 

	mov esi, [ebp + nb203nf_pos]
	lea   eax, [eax + eax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [esp + nb203nf_ixM]
	movss xmm4, [esp + nb203nf_iyM]
	movss xmm5, [esp + nb203nf_izM]
		
	movlps xmm6, [esp + nb203nf_ixH1]
	movlps xmm7, [esp + nb203nf_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + nb203nf_iyH1]
	movlps xmm7, [esp + nb203nf_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + nb203nf_izH1]
	movlps xmm7, [esp + nb203nf_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

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
	mulps xmm0, [esp + nb203nf_krf]
	movaps [esp + nb203nf_krsqM], xmm0
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb203nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb203nf_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# constant 11100000

	;# xmm0=rinv 
	addps  xmm0, [esp + nb203nf_krsqM]
	subps  xmm0, [esp + nb203nf_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm0, [esp + nb203nf_qqM]	;# xmm0=vcoul 
	addps  xmm0, [esp + nb203nf_vctot]
	movaps [esp + nb203nf_vctot], xmm0
	
	dec dword ptr [esp + nb203nf_innerk]
	jz    .nb203nf_updateouterdata
	jmp   .nb203nf_odd_loop
.nb203nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb203nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb203nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb203nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb203nf_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb203nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb203nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb203nf_n], esi
        jmp .nb203nf_outer
.nb203nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb203nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb203nf_end
        ;# non-zero, do one more workunit
        jmp   .nb203nf_threadloop
.nb203nf_end:
	emms

	mov eax, [esp + nb203nf_nouter]
	mov ebx, [esp + nb203nf_ninner]
	mov ecx, [ebp + nb203nf_outeriter]
	mov edx, [ebp + nb203nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb203nf_salign]
	add esp, eax
	add esp, 376
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret
