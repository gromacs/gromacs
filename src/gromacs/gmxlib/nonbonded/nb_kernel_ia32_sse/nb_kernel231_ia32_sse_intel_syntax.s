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



.globl nb_kernel231_ia32_sse
.globl _nb_kernel231_ia32_sse
nb_kernel231_ia32_sse:	
_nb_kernel231_ia32_sse:	
.equiv          nb231_p_nri,            8
.equiv          nb231_iinr,             12
.equiv          nb231_jindex,           16
.equiv          nb231_jjnr,             20
.equiv          nb231_shift,            24
.equiv          nb231_shiftvec,         28
.equiv          nb231_fshift,           32
.equiv          nb231_gid,              36
.equiv          nb231_pos,              40
.equiv          nb231_faction,          44
.equiv          nb231_charge,           48
.equiv          nb231_p_facel,          52
.equiv          nb231_argkrf,           56
.equiv          nb231_argcrf,           60
.equiv          nb231_Vc,               64
.equiv          nb231_type,             68
.equiv          nb231_p_ntype,          72
.equiv          nb231_vdwparam,         76
.equiv          nb231_Vvdw,             80
.equiv          nb231_p_tabscale,       84
.equiv          nb231_VFtab,            88
.equiv          nb231_invsqrta,         92
.equiv          nb231_dvda,             96
.equiv          nb231_p_gbtabscale,     100
.equiv          nb231_GBtab,            104
.equiv          nb231_p_nthreads,       108
.equiv          nb231_count,            112
.equiv          nb231_mtx,              116
.equiv          nb231_outeriter,        120
.equiv          nb231_inneriter,        124
.equiv          nb231_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb231_ixO,              0
.equiv          nb231_iyO,              16
.equiv          nb231_izO,              32
.equiv          nb231_ixH1,             48
.equiv          nb231_iyH1,             64
.equiv          nb231_izH1,             80
.equiv          nb231_ixH2,             96
.equiv          nb231_iyH2,             112
.equiv          nb231_izH2,             128
.equiv          nb231_iqO,              144
.equiv          nb231_iqH,              160
.equiv          nb231_dxO,              176
.equiv          nb231_dyO,              192
.equiv          nb231_dzO,              208
.equiv          nb231_dxH1,             224
.equiv          nb231_dyH1,             240
.equiv          nb231_dzH1,             256
.equiv          nb231_dxH2,             272
.equiv          nb231_dyH2,             288
.equiv          nb231_dzH2,             304
.equiv          nb231_qqO,              320
.equiv          nb231_qqH,              336
.equiv          nb231_c6,               352
.equiv          nb231_c12,              368
.equiv          nb231_tsc,              384
.equiv          nb231_fstmp,            400
.equiv          nb231_vctot,            416
.equiv          nb231_Vvdwtot,          432
.equiv          nb231_fixO,             448
.equiv          nb231_fiyO,             464
.equiv          nb231_fizO,             480
.equiv          nb231_fixH1,            496
.equiv          nb231_fiyH1,            512
.equiv          nb231_fizH1,            528
.equiv          nb231_fixH2,            544
.equiv          nb231_fiyH2,            560
.equiv          nb231_fizH2,            576
.equiv          nb231_fjx,              592
.equiv          nb231_fjy,              608
.equiv          nb231_fjz,              624
.equiv          nb231_half,             640
.equiv          nb231_three,            656
.equiv          nb231_two,              672
.equiv          nb231_krf,              688
.equiv          nb231_crf,              704
.equiv          nb231_krsqO,            720
.equiv          nb231_krsqH1,           736
.equiv          nb231_krsqH2,           752
.equiv          nb231_rinvO,            768
.equiv          nb231_rinvH1,           784
.equiv          nb231_rinvH2,           800
.equiv          nb231_is3,              816
.equiv          nb231_ii3,              820
.equiv          nb231_ntia,             824
.equiv          nb231_innerjjnr,        828
.equiv          nb231_innerk,           832
.equiv          nb231_n,                836
.equiv          nb231_nn1,              840
.equiv          nb231_nri,              844
.equiv          nb231_nouter,           848
.equiv          nb231_ninner,           852
.equiv          nb231_salign,           856
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 860		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb231_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb231_p_nri]
	mov ecx, [ecx]
	mov [esp + nb231_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb231_nouter], eax
	mov [esp + nb231_ninner], eax

	mov eax, [ebp + nb231_p_tabscale]
	movss xmm3, [eax]
	shufps xmm3, xmm3, 0
	movaps [esp + nb231_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb231_half], eax
	movss xmm1, [esp + nb231_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb231_half],  xmm1
	movaps [esp + nb231_two],  xmm2
	movaps [esp + nb231_three],  xmm3
	
	mov esi, [ebp + nb231_argkrf]
	mov edi, [ebp + nb231_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	movaps [esp + nb231_krf], xmm5
	movaps [esp + nb231_crf], xmm6
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb231_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb231_charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, [edx + ebx*4 + 4]	
	mov esi, [ebp + nb231_p_facel]
	movss xmm5, [esi]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [esp + nb231_iqO], xmm3
	movaps [esp + nb231_iqH], xmm4
	
	mov   edx, [ebp + nb231_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov edi, [ebp + nb231_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	mov   [esp + nb231_ntia], ecx		

.nb231_threadloop:
        mov   esi, [ebp + nb231_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb231_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb231_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb231_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb231_n], eax
        mov [esp + nb231_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb231_outerstart
        jmp .nb231_end

.nb231_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb231_nouter]
	mov [esp + nb231_nouter], ebx

.nb231_outer:
	mov   eax, [ebp + nb231_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax +esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb231_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb231_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb231_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb231_pos]    ;# eax = base of pos[]  
	mov   [esp + nb231_ii3], ebx

	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb231_ixO], xmm3
	movaps [esp + nb231_iyO], xmm4
	movaps [esp + nb231_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb231_ixH1], xmm0
	movaps [esp + nb231_iyH1], xmm1
	movaps [esp + nb231_izH1], xmm2
	movaps [esp + nb231_ixH2], xmm3
	movaps [esp + nb231_iyH2], xmm4
	movaps [esp + nb231_izH2], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb231_vctot], xmm4
	movaps [esp + nb231_Vvdwtot], xmm4
	movaps [esp + nb231_fixO], xmm4
	movaps [esp + nb231_fiyO], xmm4
	movaps [esp + nb231_fizO], xmm4
	movaps [esp + nb231_fixH1], xmm4
	movaps [esp + nb231_fiyH1], xmm4
	movaps [esp + nb231_fizH1], xmm4
	movaps [esp + nb231_fixH2], xmm4
	movaps [esp + nb231_fiyH2], xmm4
	movaps [esp + nb231_fizH2], xmm4
	
	mov   eax, [ebp + nb231_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb231_pos]
	mov   edi, [ebp + nb231_faction]	
	mov   eax, [ebp + nb231_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb231_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb231_ninner]
	mov   [esp + nb231_ninner], ecx
	add   edx, 0
	mov   [esp + nb231_innerk], edx    ;# number of innerloop atoms 
	jge   .nb231_unroll_loop
	jmp   .nb231_odd_inner
.nb231_unroll_loop:
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb231_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 

	add dword ptr [esp + nb231_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb231_charge]    ;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [esp + nb231_iqO]
	mulps  xmm4, [esp + nb231_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [esp + nb231_qqO], xmm3
	movaps  [esp + nb231_qqH], xmm4
	
	mov esi, [ebp + nb231_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + nb231_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + nb231_ntia]
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

	movaps [esp + nb231_c6], xmm4
	movaps [esp + nb231_c12], xmm6

	mov esi, [ebp + nb231_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
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
	movaps xmm4, [esp + nb231_ixO]
	movaps xmm5, [esp + nb231_iyO]
	movaps xmm6, [esp + nb231_izO]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb231_dxO], xmm4
	movaps [esp + nb231_dyO], xmm5
	movaps [esp + nb231_dzO], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	;# rsqO in xmm7 
	
	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [esp + nb231_ixH1]
	movaps xmm5, [esp + nb231_iyH1]
	movaps xmm6, [esp + nb231_izH1]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb231_dxH1], xmm4
	movaps [esp + nb231_dyH1], xmm5
	movaps [esp + nb231_dzH1], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movaps xmm3, [esp + nb231_ixH2]
	movaps xmm4, [esp + nb231_iyH2]
	movaps xmm5, [esp + nb231_izH2]

	;# calc dr 
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# store dr 
	movaps [esp + nb231_dxH2], xmm3
	movaps [esp + nb231_dyH2], xmm4
	movaps [esp + nb231_dzH2], xmm5
	;# square it 
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7

	
	movaps xmm0, xmm5
	movaps xmm1, xmm6
	movaps xmm2, xmm7
	
	mulps  xmm0, [esp + nb231_krf]	
	mulps  xmm1, [esp + nb231_krf]	
	mulps  xmm2, [esp + nb231_krf]	

	movaps [esp + nb231_krsqH2], xmm0
	movaps [esp + nb231_krsqH1], xmm1
	movaps [esp + nb231_krsqO], xmm2
	
	;# start with rsqO (still in xmm7) - seed to xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb231_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb231_half]
	movaps  [esp + nb231_rinvO], xmm4	;# rinvO in xmm0, rsqO in xmm7
	
	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb231_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb231_half]
	movaps  [esp + nb231_rinvH1], xmm4	;# rinvH1 in xmm6 

	;# rsqH2 - seed in xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb231_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb231_half]
	movaps  [esp + nb231_rinvH2], xmm4	;# rinvH2 in xmm5 
	
	
	;# do O table interactions - rsqO in xmm7.
	mulps xmm7, [esp + nb231_rinvO]
	mulps xmm7, [esp + nb231_tsc] ;# rtab
	
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

	mov  esi, [ebp + nb231_VFtab]
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
	mulps  xmm7, [esp + nb231_two]	;# two*Heps2 
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [esp + nb231_c6]
	mulps  xmm7, xmm4	 ;# fijD 
	mulps  xmm5, xmm4	 ;# Vvdw6 
	mulps  xmm7, [esp + nb231_tsc]
	;# put scalar force on stack Update Vvdwtot directly 
	addps  xmm5, [esp + nb231_Vvdwtot]
	movaps [esp + nb231_fstmp], xmm7
	movaps [esp + nb231_Vvdwtot], xmm5

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
	mulps  xmm7, [esp + nb231_two]	;# two*Heps2 
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 	
	movaps xmm4, [esp + nb231_c12]
	mulps  xmm7, xmm4 ;# fijR 
	mulps  xmm5, xmm4 ;# Vvdw12 
	mulps  xmm7, [esp + nb231_tsc]
	addps  xmm7, [esp + nb231_fstmp]
	movaps [esp + nb231_fstmp], xmm7

	addps  xmm5, [esp + nb231_Vvdwtot]
	movaps [esp + nb231_Vvdwtot], xmm5
	
	movaps xmm2, [esp + nb231_rinvO]
	movaps xmm1, [esp + nb231_krsqO]
	movaps xmm3, xmm2
	addps  xmm2, xmm1 ;# rinv+krsq
	subps  xmm2, [esp + nb231_crf] ;# rinv+krsq-crf
	movaps xmm0, xmm3
	subps  xmm3, xmm1
	subps  xmm3, xmm1 ;# rinv-2*krsq
	mulps  xmm2, [esp + nb231_qqO]
	mulps  xmm3, [esp + nb231_qqO]

	mulps  xmm3, xmm0
	subps  xmm3, [esp + nb231_fstmp]
	mulps  xmm3, xmm0

	addps  xmm2, [esp + nb231_vctot]
	movaps [esp + nb231_vctot], xmm2

	movaps xmm0, [esp + nb231_dxO]
	movaps xmm1, [esp + nb231_dyO]
	movaps xmm2, [esp + nb231_dzO]
	mulps  xmm0, xmm3
	mulps  xmm1, xmm3
	mulps  xmm2, xmm3

	movd eax, mm0	
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3	

	;# update O forces 
	movaps xmm3, [esp + nb231_fixO]
	movaps xmm4, [esp + nb231_fiyO]
	movaps xmm7, [esp + nb231_fizO]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb231_fixO], xmm3
	movaps [esp + nb231_fiyO], xmm4
	movaps [esp + nb231_fizO], xmm7
	;# update j forces with water O 
	movaps [esp + nb231_fjx], xmm0
	movaps [esp + nb231_fjy], xmm1
	movaps [esp + nb231_fjz], xmm2

	;# H1 interactions 
	movaps xmm2, [esp + nb231_rinvH1]
	movaps xmm1, [esp + nb231_krsqH1]
	movaps xmm3, xmm2
	addps  xmm2, xmm1 ;# rinv+krsq
	subps  xmm2, [esp + nb231_crf] ;# rinv+krsq-crf
	movaps xmm0, xmm3
	subps  xmm3, xmm1
	subps  xmm3, xmm1 ;# rinv-2*krsq
	mulps  xmm0, xmm0
	mulps  xmm2, [esp + nb231_qqH]
	mulps  xmm3, [esp + nb231_qqH]	
	
	mulps  xmm3, xmm0
	addps  xmm2, [esp + nb231_vctot]
	movaps [esp + nb231_vctot], xmm2

	movaps xmm0, [esp + nb231_dxH1]
	movaps xmm1, [esp + nb231_dyH1]
	movaps xmm2, [esp + nb231_dzH1]

	mulps  xmm0, xmm3
	mulps  xmm1, xmm3
	mulps  xmm2, xmm3

	;# update H1 forces 
	movaps xmm3, [esp + nb231_fixH1]
	movaps xmm4, [esp + nb231_fiyH1]
	movaps xmm7, [esp + nb231_fizH1]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb231_fixH1], xmm3
	movaps [esp + nb231_fiyH1], xmm4
	movaps [esp + nb231_fizH1], xmm7
	;# update j forces with water H1 
	addps  xmm0, [esp + nb231_fjx]
	addps  xmm1, [esp + nb231_fjy]
	addps  xmm2, [esp + nb231_fjz]
	movaps [esp + nb231_fjx], xmm0
	movaps [esp + nb231_fjy], xmm1
	movaps [esp + nb231_fjz], xmm2

	;# H2 interactions 
	movaps xmm2, [esp + nb231_rinvH2]
	movaps xmm1, [esp + nb231_krsqH2]
	movaps xmm3, xmm2
	addps  xmm2, xmm1 ;# rinv+krsq
	subps  xmm2, [esp + nb231_crf] ;# rinv+krsq-crf
	movaps xmm0, xmm3
	subps  xmm3, xmm1
	subps  xmm3, xmm1 ;# rinv-2*krsq
	mulps  xmm0, xmm0
	mulps  xmm2, [esp + nb231_qqH]
	mulps  xmm3, [esp + nb231_qqH]	
	
	mulps  xmm3, xmm0
	addps  xmm2, [esp + nb231_vctot]
	movaps [esp + nb231_vctot], xmm2

	movaps xmm0, [esp + nb231_dxH2]
	movaps xmm1, [esp + nb231_dyH2]
	movaps xmm2, [esp + nb231_dzH2]

	mulps  xmm0, xmm3
	mulps  xmm1, xmm3
	mulps  xmm2, xmm3
	;# update H2 forces 
	movaps xmm3, [esp + nb231_fixH2]
	movaps xmm4, [esp + nb231_fiyH2]
	movaps xmm7, [esp + nb231_fizH2]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb231_fixH2], xmm3
	movaps [esp + nb231_fiyH2], xmm4
	movaps [esp + nb231_fizH2], xmm7

	mov edi, [ebp + nb231_faction]
	;# update j forces 
	addps xmm0, [esp + nb231_fjx]
	addps xmm1, [esp + nb231_fjy]
	addps xmm2, [esp + nb231_fjz]

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
	sub dword ptr [esp + nb231_innerk],  4
	jl    .nb231_odd_inner
	jmp   .nb231_unroll_loop
.nb231_odd_inner:	
	add dword ptr [esp + nb231_innerk],  4
	jnz   .nb231_odd_loop
	jmp   .nb231_updateouterdata
.nb231_odd_loop:
	mov   edx, [esp + nb231_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb231_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + nb231_iqO]
	mov esi, [ebp + nb231_charge] 
	movhps xmm4, [esp + nb231_iqH]     
	movss xmm3, [esi + eax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [esp + nb231_qqO], xmm3	;# use oxygen qq for storage 

	xorps xmm6, xmm6
	mov esi, [ebp + nb231_type]
	mov ebx, [esi + eax*4]
	mov esi, [ebp + nb231_vdwparam]
	shl ebx, 1	
	add ebx, [esp + nb231_ntia]
	movlps xmm6, [esi + ebx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# constant 11111100
	shufps xmm7, xmm7, 253  ;# constant 11111101
	movaps [esp + nb231_c6], xmm6
	movaps [esp + nb231_c12], xmm7

	mov esi, [ebp + nb231_pos]
	lea   eax, [eax + eax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [esp + nb231_ixO]
	movss xmm4, [esp + nb231_iyO]
	movss xmm5, [esp + nb231_izO]
		
	movlps xmm6, [esp + nb231_ixH1]
	movlps xmm7, [esp + nb231_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + nb231_iyH1]
	movlps xmm7, [esp + nb231_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + nb231_izH1]
	movlps xmm7, [esp + nb231_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	movaps [esp + nb231_dxO], xmm3
	movaps [esp + nb231_dyO], xmm4
	movaps [esp + nb231_dzO], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 

	movaps xmm0, xmm4
	mulps xmm0, [esp + nb231_krf]
	movaps [esp + nb231_krsqO], xmm0
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb231_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb231_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# constant 11100000
		
	;# rsq still in xmm4, rinv in xmm0.
	mulps  xmm4, xmm0
	mulps  xmm4, [esp + nb231_tsc] ;# rtab
	
	cvttps2pi mm6, xmm4
	cvtpi2ps xmm6, mm6
	subss  xmm4, xmm6	
	movss xmm1, xmm4	;# xmm1=eps 
	movss xmm2, xmm1	
	mulss  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3

	movd mm0, eax	

	mov  esi, [ebp + nb231_VFtab]
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
	mulss  xmm7, [esp + nb231_two]	;# two*Heps2 
	addss  xmm7, xmm6
	addss  xmm7, xmm5 ;# xmm7=FF 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 

	movss xmm4, [esp + nb231_c6]
	mulss  xmm7, xmm4	 ;# fijD 
	mulss  xmm5, xmm4	 ;# Vvdw6 
	mulss  xmm7, [esp + nb231_tsc]
	;# put scalar force on stack Update Vvdwtot directly 
	addss  xmm5, [esp + nb231_Vvdwtot]
	movss [esp + nb231_fstmp], xmm7
	movss [esp + nb231_Vvdwtot], xmm5

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
	mulss  xmm7, [esp + nb231_two]	;# two*Heps2 
	addss  xmm7, xmm6
	addss  xmm7, xmm5 ;# xmm7=FF 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
 	
	movss xmm4, [esp + nb231_c12]
	mulss  xmm7, xmm4 ;# fijR 
	mulss  xmm5, xmm4 ;# Vvdw12 
	mulss  xmm7, [esp + nb231_tsc]
	addss  xmm7, [esp + nb231_fstmp]
	movss [esp + nb231_fstmp], xmm7

	addss  xmm5, [esp + nb231_Vvdwtot]
	movss [esp + nb231_Vvdwtot], xmm5

	movaps xmm1, xmm0	;# xmm1=rinv 
	movaps xmm4, xmm0
	movaps xmm3, [esp + nb231_krsqO]
	addps  xmm0, xmm3	;# xmm0=rinv+ krsq 
	mulps  xmm3, [esp + nb231_two]
	subps  xmm0, [esp + nb231_crf] ;# xmm0=rinv+ krsq-crf 
	subps  xmm1, xmm3	;# xmm1=rinv-2*krsq 
	mulps  xmm0, [esp + nb231_qqO]	;# xmm0=vcoul 
	mulps  xmm1, [esp + nb231_qqO] 	;# xmm1=coul part of fs 

	mulps  xmm1, xmm4
	subss  xmm1, [esp + nb231_fstmp]
	mulps  xmm4, xmm1
	
	addps  xmm0, [esp + nb231_vctot]
	movaps [esp + nb231_vctot], xmm0
	
	movaps xmm0, [esp + nb231_dxO]
	movaps xmm1, [esp + nb231_dyO]
	movaps xmm2, [esp + nb231_dzO]

	movd eax, mm0	

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	movss  xmm3, [esp + nb231_fixO]	
	movss  xmm4, [esp + nb231_fiyO]	
	movss  xmm5, [esp + nb231_fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esp + nb231_fixO], xmm3	
	movss  [esp + nb231_fiyO], xmm4	
	movss  [esp + nb231_fizO], xmm5	;# updated the O force now do the H's 
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 230 ;# constant 11100110	;# shift right 
	shufps xmm4, xmm4, 230 ;# constant 11100110
	shufps xmm5, xmm5, 230 ;# constant 11100110
	addss  xmm3, [esp + nb231_fixH1]
	addss  xmm4, [esp + nb231_fiyH1]
	addss  xmm5, [esp + nb231_fizH1]
	movss  [esp + nb231_fixH1], xmm3	
	movss  [esp + nb231_fiyH1], xmm4	
	movss  [esp + nb231_fizH1], xmm5	;# updated the H1 force 

	mov edi, [ebp + nb231_faction]
	shufps xmm3, xmm3, 231 ;# constant 11100111	;# shift right 
	shufps xmm4, xmm4, 231 ;# constant 11100111
	shufps xmm5, xmm5, 231 ;# constant 11100111
	addss  xmm3, [esp + nb231_fixH2]
	addss  xmm4, [esp + nb231_fiyH2]
	addss  xmm5, [esp + nb231_fizH2]
	movss  [esp + nb231_fixH2], xmm3	
	movss  [esp + nb231_fiyH2], xmm4	
	movss  [esp + nb231_fizH2], xmm5	;# updated the H2 force 

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
	addss    xmm2, xmm1    ;# z sum in xmm2 
	subps    xmm6, xmm0
	subss    xmm7, xmm2
	
	movlps [edi + eax*4],     xmm6
	movss  [edi + eax*4 + 8], xmm7

	dec dword ptr [esp + nb231_innerk]
	jz    .nb231_updateouterdata
	jmp   .nb231_odd_loop
.nb231_updateouterdata:
	mov   ecx, [esp + nb231_ii3]
	mov   edi, [ebp + nb231_faction]
	mov   esi, [ebp + nb231_fshift]
	mov   edx, [esp + nb231_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb231_fixO]
	movaps xmm1, [esp + nb231_fiyO]
	movaps xmm2, [esp + nb231_fizO]

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
	movaps xmm0, [esp + nb231_fixH1]
	movaps xmm1, [esp + nb231_fiyH1]
	movaps xmm2, [esp + nb231_fizH1]

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
	movaps xmm0, [esp + nb231_fixH2]
	movaps xmm1, [esp + nb231_fiyH2]
	movaps xmm2, [esp + nb231_fizH2]

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

	;# increment fshift force  
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	;# get n from stack
	mov esi, [esp + nb231_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb231_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb231_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb231_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + nb231_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb231_Vvdw]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb231_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb231_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb231_n], esi
        jmp .nb231_outer
.nb231_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb231_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb231_end
        ;# non-zero, do one more workunit
        jmp   .nb231_threadloop
.nb231_end:
	emms

	mov eax, [esp + nb231_nouter]
	mov ebx, [esp + nb231_ninner]
	mov ecx, [ebp + nb231_outeriter]
	mov edx, [ebp + nb231_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb231_salign]
	add esp, eax
	add esp, 860
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret



.globl nb_kernel231nf_ia32_sse
.globl _nb_kernel231nf_ia32_sse
nb_kernel231nf_ia32_sse:	
_nb_kernel231nf_ia32_sse:	
.equiv          nb231nf_p_nri,            8
.equiv          nb231nf_iinr,             12
.equiv          nb231nf_jindex,           16
.equiv          nb231nf_jjnr,             20
.equiv          nb231nf_shift,            24
.equiv          nb231nf_shiftvec,         28
.equiv          nb231nf_fshift,           32
.equiv          nb231nf_gid,              36
.equiv          nb231nf_pos,              40
.equiv          nb231nf_faction,          44
.equiv          nb231nf_charge,           48
.equiv          nb231nf_p_facel,          52
.equiv          nb231nf_argkrf,           56
.equiv          nb231nf_argcrf,           60
.equiv          nb231nf_Vc,               64
.equiv          nb231nf_type,             68
.equiv          nb231nf_p_ntype,          72
.equiv          nb231nf_vdwparam,         76
.equiv          nb231nf_Vvdw,             80
.equiv          nb231nf_p_tabscale,       84
.equiv          nb231nf_VFtab,            88
.equiv          nb231nf_invsqrta,         92
.equiv          nb231nf_dvda,             96
.equiv          nb231nf_p_gbtabscale,     100
.equiv          nb231nf_GBtab,            104
.equiv          nb231nf_p_nthreads,       108
.equiv          nb231nf_count,            112
.equiv          nb231nf_mtx,              116
.equiv          nb231nf_outeriter,        120
.equiv          nb231nf_inneriter,        124
.equiv          nb231nf_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb231nf_ixO,              0
.equiv          nb231nf_iyO,              16
.equiv          nb231nf_izO,              32
.equiv          nb231nf_ixH1,             48
.equiv          nb231nf_iyH1,             64
.equiv          nb231nf_izH1,             80
.equiv          nb231nf_ixH2,             96
.equiv          nb231nf_iyH2,             112
.equiv          nb231nf_izH2,             128
.equiv          nb231nf_iqO,              144
.equiv          nb231nf_iqH,              160
.equiv          nb231nf_qqO,              176
.equiv          nb231nf_qqH,              192
.equiv          nb231nf_c6,               208
.equiv          nb231nf_c12,              224
.equiv          nb231nf_tsc,              240
.equiv          nb231nf_vctot,            256
.equiv          nb231nf_Vvdwtot,          272
.equiv          nb231nf_half,             288
.equiv          nb231nf_three,            304
.equiv          nb231nf_krf,              320
.equiv          nb231nf_crf,              336
.equiv          nb231nf_krsqO,            352
.equiv          nb231nf_krsqH1,           368
.equiv          nb231nf_krsqH2,           384
.equiv          nb231nf_rinvO,            400
.equiv          nb231nf_rinvH1,           416
.equiv          nb231nf_rinvH2,           432
.equiv          nb231nf_is3,              448
.equiv          nb231nf_ii3,              452
.equiv          nb231nf_ntia,             456
.equiv          nb231nf_innerjjnr,        460
.equiv          nb231nf_innerk,           464
.equiv          nb231nf_n,                468
.equiv          nb231nf_nn1,              472
.equiv          nb231nf_nri,              476
.equiv          nb231nf_nouter,           480
.equiv          nb231nf_ninner,           484
.equiv          nb231nf_salign,           488
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 492		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb231nf_salign], eax


	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb231nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb231nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb231nf_nouter], eax
	mov [esp + nb231nf_ninner], eax

	mov eax, [ebp + nb231nf_p_tabscale]
	movss xmm3, [eax]
	shufps xmm3, xmm3, 0
	movaps [esp + nb231nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb231nf_half], eax
	movss xmm1, [esp + nb231nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb231nf_half],  xmm1
	movaps [esp + nb231nf_three],  xmm3
	
	mov esi, [ebp + nb231nf_argkrf]
	mov edi, [ebp + nb231nf_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	movaps [esp + nb231nf_krf], xmm5
	movaps [esp + nb231nf_crf], xmm6
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb231nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb231nf_charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, [edx + ebx*4 + 4]	
	mov esi, [ebp + nb231nf_p_facel]
	movss xmm5, [esi]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [esp + nb231nf_iqO], xmm3
	movaps [esp + nb231nf_iqH], xmm4
	
	mov   edx, [ebp + nb231nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov edi, [ebp + nb231nf_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	mov   [esp + nb231nf_ntia], ecx		

.nb231nf_threadloop:
        mov   esi, [ebp + nb231nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb231nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb231nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb231nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb231nf_n], eax
        mov [esp + nb231nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb231nf_outerstart
        jmp .nb231nf_end

.nb231nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb231nf_nouter]
	mov [esp + nb231nf_nouter], ebx

.nb231nf_outer:
	mov   eax, [ebp + nb231nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax +esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb231nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb231nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb231nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb231nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb231nf_ii3], ebx

	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb231nf_ixO], xmm3
	movaps [esp + nb231nf_iyO], xmm4
	movaps [esp + nb231nf_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb231nf_ixH1], xmm0
	movaps [esp + nb231nf_iyH1], xmm1
	movaps [esp + nb231nf_izH1], xmm2
	movaps [esp + nb231nf_ixH2], xmm3
	movaps [esp + nb231nf_iyH2], xmm4
	movaps [esp + nb231nf_izH2], xmm5
	
	;# clear vctot
	xorps xmm4, xmm4
	movaps [esp + nb231nf_vctot], xmm4
	movaps [esp + nb231nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb231nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb231nf_pos]
	mov   eax, [ebp + nb231nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb231nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb231nf_ninner]
	mov   [esp + nb231nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb231nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb231nf_unroll_loop
	jmp   .nb231nf_odd_inner
.nb231nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb231nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 

	add dword ptr [esp + nb231nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb231nf_charge]    ;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [esp + nb231nf_iqO]
	mulps  xmm4, [esp + nb231nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [esp + nb231nf_qqO], xmm3
	movaps  [esp + nb231nf_qqH], xmm4
	
	mov esi, [ebp + nb231nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + nb231nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + nb231nf_ntia]
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

	movaps [esp + nb231nf_c6], xmm4
	movaps [esp + nb231nf_c12], xmm6

	mov esi, [ebp + nb231nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
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
	movaps xmm4, [esp + nb231nf_ixO]
	movaps xmm5, [esp + nb231nf_iyO]
	movaps xmm6, [esp + nb231nf_izO]

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
	movaps xmm4, [esp + nb231nf_ixH1]
	movaps xmm5, [esp + nb231nf_iyH1]
	movaps xmm6, [esp + nb231nf_izH1]

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
	movaps xmm3, [esp + nb231nf_ixH2]
	movaps xmm4, [esp + nb231nf_iyH2]
	movaps xmm5, [esp + nb231nf_izH2]

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
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7

	
	movaps xmm0, xmm5
	movaps xmm1, xmm6
	movaps xmm2, xmm7
	
	mulps  xmm0, [esp + nb231nf_krf]	
	mulps  xmm1, [esp + nb231nf_krf]	
	mulps  xmm2, [esp + nb231nf_krf]	

	movaps [esp + nb231nf_krsqH2], xmm0
	movaps [esp + nb231nf_krsqH1], xmm1
	movaps [esp + nb231nf_krsqO], xmm2
	
	;# start with rsqO (still in xmm7) - seed to xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb231nf_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb231nf_half]
	movaps  [esp + nb231nf_rinvO], xmm4	;# rinvO in xmm0, rsqO in xmm7
	
	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb231nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb231nf_half]
	movaps  [esp + nb231nf_rinvH1], xmm4	;# rinvH1 in xmm6 

	;# rsqH2 - seed in xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb231nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb231nf_half]
	movaps  [esp + nb231nf_rinvH2], xmm4	;# rinvH2 in xmm5 
	
	
	;# do O table interactions - rsqO in xmm7.
	mulps xmm7, [esp + nb231nf_rinvO]
	mulps xmm7, [esp + nb231nf_tsc] ;# rtab
	
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

	mov  esi, [ebp + nb231nf_VFtab]
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

	movaps xmm4, [esp + nb231nf_c6]
	mulps  xmm5, xmm4	 ;# Vvdw6 

	addps  xmm5, [esp + nb231nf_Vvdwtot]
	movaps [esp + nb231nf_Vvdwtot], xmm5

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
 	
	movaps xmm4, [esp + nb231nf_c12]
	mulps  xmm5, xmm4 ;# Vvdw12 

	addps  xmm5, [esp + nb231nf_Vvdwtot]
	movaps [esp + nb231nf_Vvdwtot], xmm5
	
	movaps xmm2, [esp + nb231nf_rinvO]
	movaps xmm1, [esp + nb231nf_krsqO]
	addps  xmm2, xmm1 ;# rinv+krsq
	subps  xmm2, [esp + nb231nf_crf] ;# rinv+krsq-crf
	mulps  xmm2, [esp + nb231nf_qqO]

	addps  xmm2, [esp + nb231nf_vctot]
	movaps [esp + nb231nf_vctot], xmm2

	;# H1 interactions 
	movaps xmm2, [esp + nb231nf_rinvH1]
	movaps xmm1, [esp + nb231nf_krsqH1]
	addps  xmm2, xmm1 ;# rinv+krsq
	subps  xmm2, [esp + nb231nf_crf] ;# rinv+krsq-crf
	mulps  xmm2, [esp + nb231nf_qqH]
	
	addps  xmm2, [esp + nb231nf_vctot]
	movaps [esp + nb231nf_vctot], xmm2

	;# H2 interactions 
	movaps xmm2, [esp + nb231nf_rinvH2]
	movaps xmm1, [esp + nb231nf_krsqH2]
	addps  xmm2, xmm1 ;# rinv+krsq
	subps  xmm2, [esp + nb231nf_crf] ;# rinv+krsq-crf
	mulps  xmm2, [esp + nb231nf_qqH]
	addps  xmm2, [esp + nb231nf_vctot]
	movaps [esp + nb231nf_vctot], xmm2

	;# should we do one more iteration? 
	sub dword ptr [esp + nb231nf_innerk],  4
	jl    .nb231nf_odd_inner
	jmp   .nb231nf_unroll_loop
.nb231nf_odd_inner:	
	add dword ptr [esp + nb231nf_innerk],  4
	jnz   .nb231nf_odd_loop
	jmp   .nb231nf_updateouterdata
.nb231nf_odd_loop:
	mov   edx, [esp + nb231nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb231nf_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + nb231nf_iqO]
	mov esi, [ebp + nb231nf_charge] 
	movhps xmm4, [esp + nb231nf_iqH]     
	movss xmm3, [esi + eax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [esp + nb231nf_qqO], xmm3	;# use oxygen qq for storage 

	xorps xmm6, xmm6
	mov esi, [ebp + nb231nf_type]
	mov ebx, [esi + eax*4]
	mov esi, [ebp + nb231nf_vdwparam]
	shl ebx, 1	
	add ebx, [esp + nb231nf_ntia]
	movlps xmm6, [esi + ebx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# constant 11111100
	shufps xmm7, xmm7, 253  ;# constant 11111101
	movaps [esp + nb231nf_c6], xmm6
	movaps [esp + nb231nf_c12], xmm7

	mov esi, [ebp + nb231nf_pos]
	lea   eax, [eax + eax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [esp + nb231nf_ixO]
	movss xmm4, [esp + nb231nf_iyO]
	movss xmm5, [esp + nb231nf_izO]
		
	movlps xmm6, [esp + nb231nf_ixH1]
	movlps xmm7, [esp + nb231nf_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + nb231nf_iyH1]
	movlps xmm7, [esp + nb231nf_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + nb231nf_izH1]
	movlps xmm7, [esp + nb231nf_izH2]
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
	mulps xmm0, [esp + nb231nf_krf]
	movaps [esp + nb231nf_krsqO], xmm0
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb231nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb231nf_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# constant 11100000
		
	;# rsq still in xmm4, rinv in xmm0.
	mulps  xmm4, xmm0
	mulps  xmm4, [esp + nb231nf_tsc] ;# rtab
	
	cvttps2pi mm6, xmm4
	cvtpi2ps xmm6, mm6
	subss  xmm4, xmm6	
	movss xmm1, xmm4	;# xmm1=eps 
	movss xmm2, xmm1	
	mulss  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3

	mov  esi, [ebp + nb231nf_VFtab]
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

	movss xmm4, [esp + nb231nf_c6]
	mulss  xmm5, xmm4	 ;# Vvdw6 

	addss  xmm5, [esp + nb231nf_Vvdwtot]
	movss [esp + nb231nf_Vvdwtot], xmm5

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
 	
	movss xmm4, [esp + nb231nf_c12]
	mulss  xmm5, xmm4 ;# Vvdw12 
	addss  xmm5, [esp + nb231nf_Vvdwtot]
	movss [esp + nb231nf_Vvdwtot], xmm5

	movaps xmm1, xmm0	;# xmm1=rinv 
	movaps xmm4, xmm0
	movaps xmm3, [esp + nb231nf_krsqO]
	addps  xmm0, xmm3	;# xmm0=rinv+ krsq 

	subps  xmm0, [esp + nb231nf_crf] ;# xmm0=rinv+ krsq-crf 

	mulps  xmm0, [esp + nb231nf_qqO]	;# xmm0=vcoul 
	
	addps  xmm0, [esp + nb231nf_vctot]
	movaps [esp + nb231nf_vctot], xmm0
	
	dec dword ptr [esp + nb231nf_innerk]
	jz    .nb231nf_updateouterdata
	jmp   .nb231nf_odd_loop
.nb231nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb231nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb231nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb231nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb231nf_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + nb231nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb231nf_Vvdw]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb231nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb231nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb231nf_n], esi
        jmp .nb231nf_outer
.nb231nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb231nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb231nf_end
        ;# non-zero, do one more workunit
        jmp   .nb231nf_threadloop
.nb231nf_end:
	emms

	mov eax, [esp + nb231nf_nouter]
	mov ebx, [esp + nb231nf_ninner]
	mov ecx, [ebp + nb231nf_outeriter]
	mov edx, [ebp + nb231nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb231nf_salign]
	add esp, eax
	add esp, 492
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


