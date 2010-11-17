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

	
	

.globl nb_kernel103_ia32_sse2
.globl _nb_kernel103_ia32_sse2
nb_kernel103_ia32_sse2:	
_nb_kernel103_ia32_sse2:	
.equiv          nb103_p_nri,            8
.equiv          nb103_iinr,             12
.equiv          nb103_jindex,           16
.equiv          nb103_jjnr,             20
.equiv          nb103_shift,            24
.equiv          nb103_shiftvec,         28
.equiv          nb103_fshift,           32
.equiv          nb103_gid,              36
.equiv          nb103_pos,              40
.equiv          nb103_faction,          44
.equiv          nb103_charge,           48
.equiv          nb103_p_facel,          52
.equiv          nb103_argkrf,           56
.equiv          nb103_argcrf,           60
.equiv          nb103_Vc,               64
.equiv          nb103_type,             68
.equiv          nb103_p_ntype,          72
.equiv          nb103_vdwparam,         76
.equiv          nb103_Vvdw,             80
.equiv          nb103_p_tabscale,       84
.equiv          nb103_VFtab,            88
.equiv          nb103_invsqrta,         92
.equiv          nb103_dvda,             96
.equiv          nb103_p_gbtabscale,     100
.equiv          nb103_GBtab,            104
.equiv          nb103_p_nthreads,       108
.equiv          nb103_count,            112
.equiv          nb103_mtx,              116
.equiv          nb103_outeriter,        120
.equiv          nb103_inneriter,        124
.equiv          nb103_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb103_ixH1,             0
.equiv          nb103_iyH1,             16
.equiv          nb103_izH1,             32
.equiv          nb103_ixH2,             48
.equiv          nb103_iyH2,             64
.equiv          nb103_izH2,             80
.equiv          nb103_ixM,              96
.equiv          nb103_iyM,              112
.equiv          nb103_izM,              128
.equiv          nb103_iqM,              144
.equiv          nb103_iqH,              160
.equiv          nb103_dxH1,             176
.equiv          nb103_dyH1,             192
.equiv          nb103_dzH1,             208
.equiv          nb103_dxH2,             224
.equiv          nb103_dyH2,             240
.equiv          nb103_dzH2,             256
.equiv          nb103_dxM,              272
.equiv          nb103_dyM,              288
.equiv          nb103_dzM,              304
.equiv          nb103_qqM,              320
.equiv          nb103_qqH,              336
.equiv          nb103_vctot,            352
.equiv          nb103_fixM,             368
.equiv          nb103_fiyM,             384
.equiv          nb103_fizM,             400
.equiv          nb103_fixH1,            416
.equiv          nb103_fiyH1,            432
.equiv          nb103_fizH1,            448
.equiv          nb103_fixH2,            464
.equiv          nb103_fiyH2,            480
.equiv          nb103_fizH2,            496
.equiv          nb103_fjx,              512
.equiv          nb103_fjy,              528
.equiv          nb103_fjz,              544
.equiv          nb103_half,             560
.equiv          nb103_three,            576
.equiv          nb103_is3,              592
.equiv          nb103_ii3,              596
.equiv          nb103_innerjjnr,        600
.equiv          nb103_innerk,           604
.equiv          nb103_n,                608
.equiv          nb103_nn1,              612
.equiv          nb103_nri,              616
.equiv          nb103_nouter,           620
.equiv          nb103_ninner,           624
.equiv          nb103_salign,           628
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 632		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb103_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb103_p_nri]
	mov ecx, [ecx]
	mov [esp + nb103_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb103_nouter], eax
	mov [esp + nb103_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb103_half], eax
	mov [esp + nb103_half+4], ebx
	movsd xmm1, [esp + nb103_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb103_half], xmm1
	movapd [esp + nb103_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb103_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb103_charge]
	movsd xmm3, [edx + ebx*8 + 8]	
	movsd xmm4, [edx + ebx*8 + 24]	
	mov esi, [ebp + nb103_p_facel]
	movsd xmm5, [esi]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb103_iqH], xmm3
	movapd [esp + nb103_iqM], xmm4
	
.nb103_threadloop:
        mov   esi, [ebp + nb103_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb103_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb103_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb103_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb103_n], eax
        mov [esp + nb103_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb103_outerstart
        jmp .nb103_end

.nb103_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb103_nouter]
	mov [esp + nb103_nouter], ebx

.nb103_outer:
	mov   eax, [ebp + nb103_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb103_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb103_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb103_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb103_pos]    ;# eax = base of pos[]  
	mov   [esp + nb103_ii3], ebx

	addsd xmm3, [eax + ebx*8 + 24]
	addsd xmm4, [eax + ebx*8 + 32]
	addsd xmm5, [eax + ebx*8 + 40]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb103_ixH1], xmm3
	movapd [esp + nb103_iyH1], xmm4
	movapd [esp + nb103_izH1], xmm5

	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm0, [eax + ebx*8 + 48]
	addsd xmm1, [eax + ebx*8 + 56]
	addsd xmm2, [eax + ebx*8 + 64]		
	addsd xmm3, [eax + ebx*8 + 72]
	addsd xmm4, [eax + ebx*8 + 80]
	addsd xmm5, [eax + ebx*8 + 88]		

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb103_ixH2], xmm0
	movapd [esp + nb103_iyH2], xmm1
	movapd [esp + nb103_izH2], xmm2
	movapd [esp + nb103_ixM], xmm3
	movapd [esp + nb103_iyM], xmm4
	movapd [esp + nb103_izM], xmm5
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb103_vctot], xmm4
	movapd [esp + nb103_fixM], xmm4
	movapd [esp + nb103_fiyM], xmm4
	movapd [esp + nb103_fizM], xmm4
	movapd [esp + nb103_fixH1], xmm4
	movapd [esp + nb103_fiyH1], xmm4
	movapd [esp + nb103_fizH1], xmm4
	movapd [esp + nb103_fixH2], xmm4
	movapd [esp + nb103_fiyH2], xmm4
	movapd [esp + nb103_fizH2], xmm4
	
	mov   eax, [ebp + nb103_jindex]
	mov   ecx, [eax+esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb103_pos]
	mov   edi, [ebp + nb103_faction]	
	mov   eax, [ebp + nb103_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb103_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb103_ninner]
	mov   [esp + nb103_ninner], ecx
	add   edx, 0
	mov   [esp + nb103_innerk], edx    ;# number of innerloop atoms 
	jge   .nb103_unroll_loop
	jmp   .nb103_checksingle
.nb103_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb103_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              

	add dword ptr [esp + nb103_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb103_charge]    ;# base of charge[] 
	
	
	movlpd xmm6, [esi + eax*8]	;# jq A 
	movhpd xmm6, [esi + ebx*8]	;# jq B 
	movapd xmm3, [esp + nb103_iqM]
	movapd xmm4, [esp + nb103_iqH]
	mulpd xmm3, xmm6		;# qqM 
	mulpd xmm4, xmm6		;# qqH 
	
	movapd  [esp + nb103_qqM], xmm3
	movapd  [esp + nb103_qqH], xmm4	

	mov esi, [ebp + nb103_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		

	;# move ixM-izM to xmm4-xmm6 
	movapd xmm4, [esp + nb103_ixM]
	movapd xmm5, [esp + nb103_iyM]
	movapd xmm6, [esp + nb103_izM]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb103_dxM], xmm4
	movapd [esp + nb103_dyM], xmm5
	movapd [esp + nb103_dzM], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb103_ixH1]
	movapd xmm5, [esp + nb103_iyH1]
	movapd xmm6, [esp + nb103_izH1]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb103_dxH1], xmm4
	movapd [esp + nb103_dyH1], xmm5
	movapd [esp + nb103_dzH1], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm6, xmm5
	addpd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb103_ixH2]
	movapd xmm4, [esp + nb103_iyH2]
	movapd xmm5, [esp + nb103_izH2]

	;# calc dr 
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2

	;# store dr 
	movapd [esp + nb103_dxH2], xmm3
	movapd [esp + nb103_dyH2], xmm4
	movapd [esp + nb103_dzH2], xmm5
	;# square it 
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	addpd xmm5, xmm4
	addpd xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqM - put seed in xmm2 
	cvtpd2ps xmm2, xmm7	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb103_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb103_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb103_three]
	subpd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb103_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvM in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb103_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb103_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb103_three]
	subpd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb103_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb103_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb103_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb103_three]
	subpd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb103_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do M interactions 
	movapd  xmm4, xmm7	
	mulpd   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq 
	mulpd  xmm7, [esp + nb103_qqM]	;# xmm7=vcoul 
	
	mulpd  xmm4, xmm7	;# total fsO in xmm4 

	addpd  xmm7, [esp + nb103_vctot]
	
	movapd [esp + nb103_vctot], xmm7

	movapd xmm0, [esp + nb103_dxM]
	movapd xmm1, [esp + nb103_dyM]
	movapd xmm2, [esp + nb103_dzM]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update M forces 
	movapd xmm3, [esp + nb103_fixM]
	movapd xmm4, [esp + nb103_fiyM]
	movapd xmm7, [esp + nb103_fizM]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb103_fixM], xmm3
	movapd [esp + nb103_fiyM], xmm4
	movapd [esp + nb103_fizM], xmm7
	;# update j forces with water M
	movapd [esp + nb103_fjx], xmm0
	movapd [esp + nb103_fjy], xmm1
	movapd [esp + nb103_fjz], xmm2

	;# H1 interactions 
	movapd  xmm4, xmm6	
	mulpd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	mulpd  xmm6, [esp + nb103_qqH]	;# xmm6=vcoul 
	mulpd  xmm4, xmm6		;# total fsH1 in xmm4 
	
	addpd  xmm6, [esp + nb103_vctot]

	movapd xmm0, [esp + nb103_dxH1]
	movapd xmm1, [esp + nb103_dyH1]
	movapd xmm2, [esp + nb103_dzH1]
	movapd [esp + nb103_vctot], xmm6
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [esp + nb103_fixH1]
	movapd xmm4, [esp + nb103_fiyH1]
	movapd xmm7, [esp + nb103_fizH1]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb103_fixH1], xmm3
	movapd [esp + nb103_fiyH1], xmm4
	movapd [esp + nb103_fizH1], xmm7
	;# update j forces with water H1 
	addpd  xmm0, [esp + nb103_fjx]
	addpd  xmm1, [esp + nb103_fjy]
	addpd  xmm2, [esp + nb103_fjz]
	movapd [esp + nb103_fjx], xmm0
	movapd [esp + nb103_fjy], xmm1
	movapd [esp + nb103_fjz], xmm2

	;# H2 interactions 
	movapd  xmm4, xmm5	
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	mulpd  xmm5, [esp + nb103_qqH]	;# xmm5=vcoul 
	mulpd  xmm4, xmm5		;# total fsH1 in xmm4 
	
	addpd  xmm5, [esp + nb103_vctot]

	movapd xmm0, [esp + nb103_dxH2]
	movapd xmm1, [esp + nb103_dyH2]
	movapd xmm2, [esp + nb103_dzH2]
	movapd [esp + nb103_vctot], xmm5
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [esp + nb103_fixH2]
	movapd xmm4, [esp + nb103_fiyH2]
	movapd xmm7, [esp + nb103_fizH2]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb103_fixH2], xmm3
	movapd [esp + nb103_fiyH2], xmm4
	movapd [esp + nb103_fizH2], xmm7

	mov edi, [ebp + nb103_faction]
	;# update j forces 
	addpd  xmm0, [esp + nb103_fjx]
	addpd  xmm1, [esp + nb103_fjy]
	addpd  xmm2, [esp + nb103_fjz]

	movlpd xmm3, [edi + eax*8]
	movlpd xmm4, [edi + eax*8 + 8]
	movlpd xmm5, [edi + eax*8 + 16]
	movhpd xmm3, [edi + ebx*8]
	movhpd xmm4, [edi + ebx*8 + 8]
	movhpd xmm5, [edi + ebx*8 + 16]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	movlpd [edi + eax*8], xmm3
	movlpd [edi + eax*8 + 8], xmm4
	movlpd [edi + eax*8 + 16], xmm5
	movhpd [edi + ebx*8], xmm3
	movhpd [edi + ebx*8 + 8], xmm4
	movhpd [edi + ebx*8 + 16], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb103_innerk],  2
	jl    .nb103_checksingle
	jmp   .nb103_unroll_loop
.nb103_checksingle:				
	mov   edx, [esp + nb103_innerk]
	and   edx, 1
	jnz    .nb103_dosingle
	jmp    .nb103_updateouterdata
.nb103_dosingle:
	mov   edx, [esp + nb103_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	

	mov esi, [ebp + nb103_charge]    ;# base of charge[] 
	xorpd xmm6, xmm6
	movlpd xmm6, [esi + eax*8]	;# jq A 
	
	movapd xmm3, [esp + nb103_iqM]
	movapd xmm4, [esp + nb103_iqH]
	mulsd xmm3, xmm6		;# qqM
	mulsd xmm4, xmm6		;# qqH 
	
	movapd  [esp + nb103_qqM], xmm3
	movapd  [esp + nb103_qqH], xmm4	

	mov esi, [ebp + nb103_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move ixM-izM to xmm4-xmm6 
	movapd xmm4, [esp + nb103_ixM]
	movapd xmm5, [esp + nb103_iyM]
	movapd xmm6, [esp + nb103_izM]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb103_dxM], xmm4
	movapd [esp + nb103_dyM], xmm5
	movapd [esp + nb103_dzM], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqM in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb103_ixH1]
	movapd xmm5, [esp + nb103_iyH1]
	movapd xmm6, [esp + nb103_izH1]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb103_dxH1], xmm4
	movapd [esp + nb103_dyH1], xmm5
	movapd [esp + nb103_dzH1], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb103_ixH2]
	movapd xmm4, [esp + nb103_iyH2]
	movapd xmm5, [esp + nb103_izH2]

	;# calc dr 
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2

	;# store dr 
	movapd [esp + nb103_dxH2], xmm3
	movapd [esp + nb103_dyH2], xmm4
	movapd [esp + nb103_dzH2], xmm5
	;# square it 
	mulsd xmm3,xmm3
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	addsd xmm5, xmm4
	addsd xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqM - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb103_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb103_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb103_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb103_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvM in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb103_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb103_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb103_three]
	subsd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb103_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb103_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb103_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb103_three]
	subsd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb103_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do M interactions 
	movapd  xmm4, xmm7	
	mulsd   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq 
	mulsd  xmm7, [esp + nb103_qqM]	;# xmm7=vcoul 
	
	mulsd  xmm4, xmm7	;# total fsM in xmm4 

	addsd  xmm7, [esp + nb103_vctot]
	
	movlpd [esp + nb103_vctot], xmm7

	movapd xmm0, [esp + nb103_dxM]
	movapd xmm1, [esp + nb103_dyM]
	movapd xmm2, [esp + nb103_dzM]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update M forces 
	movapd xmm3, [esp + nb103_fixM]
	movapd xmm4, [esp + nb103_fiyM]
	movapd xmm7, [esp + nb103_fizM]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [esp + nb103_fixM], xmm3
	movlpd [esp + nb103_fiyM], xmm4
	movlpd [esp + nb103_fizM], xmm7
	;# update j forces with water M 
	movlpd [esp + nb103_fjx], xmm0
	movlpd [esp + nb103_fjy], xmm1
	movlpd [esp + nb103_fjz], xmm2

	;# H1 interactions 
	movapd  xmm4, xmm6	
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	mulsd  xmm6, [esp + nb103_qqH]	;# xmm6=vcoul 
	mulsd  xmm4, xmm6		;# total fsH1 in xmm4 
	
	addsd  xmm6, [esp + nb103_vctot]

	movapd xmm0, [esp + nb103_dxH1]
	movapd xmm1, [esp + nb103_dyH1]
	movapd xmm2, [esp + nb103_dzH1]
	movlpd [esp + nb103_vctot], xmm6
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [esp + nb103_fixH1]
	movapd xmm4, [esp + nb103_fiyH1]
	movapd xmm7, [esp + nb103_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [esp + nb103_fixH1], xmm3
	movlpd [esp + nb103_fiyH1], xmm4
	movlpd [esp + nb103_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [esp + nb103_fjx]
	addsd  xmm1, [esp + nb103_fjy]
	addsd  xmm2, [esp + nb103_fjz]
	movsd [esp + nb103_fjx], xmm0
	movsd [esp + nb103_fjy], xmm1
	movsd [esp + nb103_fjz], xmm2

	;# H2 interactions 
	movapd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	mulsd  xmm5, [esp + nb103_qqH]	;# xmm5=vcoul 
	mulsd  xmm4, xmm5		;# total fsH1 in xmm4 
	
	addsd  xmm5, [esp + nb103_vctot]

	movapd xmm0, [esp + nb103_dxH2]
	movapd xmm1, [esp + nb103_dyH2]
	movapd xmm2, [esp + nb103_dzH2]
	movlpd [esp + nb103_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [esp + nb103_fixH2]
	movapd xmm4, [esp + nb103_fiyH2]
	movapd xmm7, [esp + nb103_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [esp + nb103_fixH2], xmm3
	movlpd [esp + nb103_fiyH2], xmm4
	movlpd [esp + nb103_fizH2], xmm7

	mov edi, [ebp + nb103_faction]
	;# update j forces 
	addsd  xmm0, [esp + nb103_fjx]
	addsd  xmm1, [esp + nb103_fjy]
	addsd  xmm2, [esp + nb103_fjz]

	movlpd xmm3, [edi + eax*8]
	movlpd xmm4, [edi + eax*8 + 8]
	movlpd xmm5, [edi + eax*8 + 16]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	movlpd [edi + eax*8], xmm3
	movlpd [edi + eax*8 + 8], xmm4
	movlpd [edi + eax*8 + 16], xmm5

.nb103_updateouterdata:
	mov   ecx, [esp + nb103_ii3]
	mov   edi, [ebp + nb103_faction]
	mov   esi, [ebp + nb103_fshift]
	mov   edx, [esp + nb103_is3]

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb103_fixH1]
	movapd xmm1, [esp + nb103_fiyH1]
	movapd xmm2, [esp + nb103_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	;# increment i force 
	movsd  xmm3, [edi + ecx*8 + 24]
	movsd  xmm4, [edi + ecx*8 + 32]
	movsd  xmm5, [edi + ecx*8 + 40]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movsd  [edi + ecx*8 + 24],     xmm3
	movsd  [edi + ecx*8 + 32], xmm4
	movsd  [edi + ecx*8 + 40], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movapd xmm6, xmm0
	movsd xmm7, xmm2
	unpcklpd xmm6, xmm1

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb103_fixH2]
	movapd xmm1, [esp + nb103_fiyH2]
	movapd xmm2, [esp + nb103_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	;# increment i force 
	movsd  xmm3, [edi + ecx*8 + 48]
	movsd  xmm4, [edi + ecx*8 + 56]
	movsd  xmm5, [edi + ecx*8 + 64]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movsd  [edi + ecx*8 + 48], xmm3
	movsd  [edi + ecx*8 + 56], xmm4
	movsd  [edi + ecx*8 + 64], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb103_fixM]
	movapd xmm1, [esp + nb103_fiyM]
	movapd xmm2, [esp + nb103_fizM]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	movapd xmm3, xmm0	
	movapd xmm4, xmm1	
	movapd xmm5, xmm2	

	;# increment i force 
	movsd  xmm3, [edi + ecx*8 + 72]
	movsd  xmm4, [edi + ecx*8 + 80]
	movsd  xmm5, [edi + ecx*8 + 88]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movsd  [edi + ecx*8 + 72], xmm3
	movsd  [edi + ecx*8 + 80], xmm4
	movsd  [edi + ecx*8 + 88], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# increment fshift force 
	movlpd xmm3, [esi + edx*8]
	movhpd xmm3, [esi + edx*8 + 8]
	movsd  xmm4, [esi + edx*8 + 16]
	addpd  xmm3, xmm6
	addsd  xmm4, xmm7
	movlpd [esi + edx*8],      xmm3
	movhpd [esi + edx*8 + 8],  xmm3
	movsd  [esi + edx*8 + 16], xmm4

	;# get n from stack
	mov esi, [esp + nb103_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb103_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb103_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	        
	;# add earlier value from mem 
	mov   eax, [ebp + nb103_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 	
	
        ;# finish if last 
        mov ecx, [esp + nb103_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb103_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb103_n], esi
        jmp .nb103_outer
.nb103_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb103_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb103_end
        ;# non-zero, do one more workunit
        jmp   .nb103_threadloop
.nb103_end:
	emms

	mov eax, [esp + nb103_nouter]
	mov ebx, [esp + nb103_ninner]
	mov ecx, [ebp + nb103_outeriter]
	mov edx, [ebp + nb103_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb103_salign]
	add esp, eax
	add esp, 632
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret



.globl nb_kernel103nf_ia32_sse2
.globl _nb_kernel103nf_ia32_sse2
nb_kernel103nf_ia32_sse2:	
_nb_kernel103nf_ia32_sse2:	
.equiv          nb103nf_p_nri,          8
.equiv          nb103nf_iinr,           12
.equiv          nb103nf_jindex,         16
.equiv          nb103nf_jjnr,           20
.equiv          nb103nf_shift,          24
.equiv          nb103nf_shiftvec,       28
.equiv          nb103nf_fshift,         32
.equiv          nb103nf_gid,            36
.equiv          nb103nf_pos,            40
.equiv          nb103nf_faction,        44
.equiv          nb103nf_charge,         48
.equiv          nb103nf_p_facel,        52
.equiv          nb103nf_argkrf,         56
.equiv          nb103nf_argcrf,         60
.equiv          nb103nf_Vc,             64
.equiv          nb103nf_type,           68
.equiv          nb103nf_p_ntype,        72
.equiv          nb103nf_vdwparam,       76
.equiv          nb103nf_Vvdw,           80
.equiv          nb103nf_p_tabscale,     84
.equiv          nb103nf_VFtab,          88
.equiv          nb103nf_invsqrta,       92
.equiv          nb103nf_dvda,           96
.equiv          nb103nf_p_gbtabscale,   100
.equiv          nb103nf_GBtab,          104
.equiv          nb103nf_p_nthreads,     108
.equiv          nb103nf_count,          112
.equiv          nb103nf_mtx,            116
.equiv          nb103nf_outeriter,      120
.equiv          nb103nf_inneriter,      124
.equiv          nb103nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb103nf_ixM,            0
.equiv          nb103nf_iyM,            16
.equiv          nb103nf_izM,            32
.equiv          nb103nf_ixH1,           48
.equiv          nb103nf_iyH1,           64
.equiv          nb103nf_izH1,           80
.equiv          nb103nf_ixH2,           96
.equiv          nb103nf_iyH2,           112
.equiv          nb103nf_izH2,           128
.equiv          nb103nf_iqM,            144
.equiv          nb103nf_iqH,            160
.equiv          nb103nf_qqM,            176
.equiv          nb103nf_qqH,            192
.equiv          nb103nf_vctot,          208
.equiv          nb103nf_half,           224
.equiv          nb103nf_three,          240
.equiv          nb103nf_is3,            256
.equiv          nb103nf_ii3,            260
.equiv          nb103nf_innerjjnr,      264
.equiv          nb103nf_innerk,         268
.equiv          nb103nf_n,              272
.equiv          nb103nf_nn1,            276
.equiv          nb103nf_nri,            280
.equiv          nb103nf_nouter,         284
.equiv          nb103nf_ninner,         288
.equiv          nb103nf_salign,         292
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 296		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb103nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb103nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb103nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb103nf_nouter], eax
	mov [esp + nb103nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb103nf_half], eax
	mov [esp + nb103nf_half+4], ebx
	movsd xmm1, [esp + nb103nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb103nf_half], xmm1
	movapd [esp + nb103nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb103nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb103nf_charge]
	movsd xmm3, [edx + ebx*8 + 8]	
	movsd xmm4, [edx + ebx*8 + 24]	
	mov esi, [ebp + nb103nf_p_facel]
	movsd xmm5, [esi]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb103nf_iqH], xmm3
	movapd [esp + nb103nf_iqM], xmm4
	
.nb103nf_threadloop:
        mov   esi, [ebp + nb103nf_count]        ;# pointer to sync counter
        mov   eax, [esi]
.nb103nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb103nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb103nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb103nf_n], eax
        mov [esp + nb103nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb103nf_outerstart
        jmp .nb103nf_end

.nb103nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb103nf_nouter]
	mov [esp + nb103nf_nouter], ebx

.nb103nf_outer:
	mov   eax, [ebp + nb103nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 

	mov   eax, [ebp + nb103nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb103nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb103nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb103nf_ii3], ebx

	addsd xmm3, [eax + ebx*8 + 24]
	addsd xmm4, [eax + ebx*8 + 32]
	addsd xmm5, [eax + ebx*8 + 40]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb103nf_ixH1], xmm3
	movapd [esp + nb103nf_iyH1], xmm4
	movapd [esp + nb103nf_izH1], xmm5

	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm0, [eax + ebx*8 + 48]
	addsd xmm1, [eax + ebx*8 + 56]
	addsd xmm2, [eax + ebx*8 + 64]		
	addsd xmm3, [eax + ebx*8 + 72]
	addsd xmm4, [eax + ebx*8 + 80]
	addsd xmm5, [eax + ebx*8 + 88]		

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb103nf_ixH2], xmm0
	movapd [esp + nb103nf_iyH2], xmm1
	movapd [esp + nb103nf_izH2], xmm2
	movapd [esp + nb103nf_ixM], xmm3
	movapd [esp + nb103nf_iyM], xmm4
	movapd [esp + nb103nf_izM], xmm5
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb103nf_vctot], xmm4
	
	mov   eax, [ebp + nb103nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb103nf_pos]
	mov   eax, [ebp + nb103nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb103nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb103nf_ninner]
	mov   [esp + nb103nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb103nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb103nf_unroll_loop
	jmp   .nb103nf_checksingle
.nb103nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb103nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              

	add dword ptr [esp + nb103nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb103nf_charge]    ;# base of charge[] 
	
	
	movlpd xmm6, [esi + eax*8]	;# jq A 
	movhpd xmm6, [esi + ebx*8]	;# jq B 
	movapd xmm3, [esp + nb103nf_iqM]
	movapd xmm4, [esp + nb103nf_iqH]
	mulpd xmm3, xmm6		;# qqM 
	mulpd xmm4, xmm6		;# qqH 
	
	movapd  [esp + nb103nf_qqM], xmm3
	movapd  [esp + nb103nf_qqH], xmm4	

	mov esi, [ebp + nb103nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		

	;# move ixM-izM to xmm4-xmm6 
	movapd xmm4, [esp + nb103nf_ixM]
	movapd xmm5, [esp + nb103nf_iyM]
	movapd xmm6, [esp + nb103nf_izM]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqM in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb103nf_ixH1]
	movapd xmm5, [esp + nb103nf_iyH1]
	movapd xmm6, [esp + nb103nf_izH1]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm6, xmm5
	addpd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb103nf_ixH2]
	movapd xmm4, [esp + nb103nf_iyH2]
	movapd xmm5, [esp + nb103nf_izH2]

	;# calc dr 
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2

	;# square it 
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	addpd xmm5, xmm4
	addpd xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqM in xmm7 

	;# start with rsqM - put seed in xmm2 
	cvtpd2ps xmm2, xmm7	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb103nf_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb103nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb103nf_three]
	subpd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb103nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvM in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb103nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb103nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb103nf_three]
	subpd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb103nf_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb103nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb103nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb103nf_three]
	subpd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb103nf_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do M interactions 
	mulpd  xmm7, [esp + nb103nf_qqM]	;# xmm7=vcoul 
	addpd  xmm7, [esp + nb103nf_vctot]
	movapd [esp + nb103nf_vctot], xmm7

	;# H1 interactions 
	mulpd  xmm6, [esp + nb103nf_qqH]	;# xmm6=vcoul 
	addpd  xmm6, [esp + nb103nf_vctot]
	movapd [esp + nb103nf_vctot], xmm6

	;# H2 interactions 
	mulpd  xmm5, [esp + nb103nf_qqH]	;# xmm5=vcoul 
	addpd  xmm5, [esp + nb103nf_vctot]
	movapd [esp + nb103nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [esp + nb103nf_innerk],  2
	jl    .nb103nf_checksingle
	jmp   .nb103nf_unroll_loop
.nb103nf_checksingle:				
	mov   edx, [esp + nb103nf_innerk]
	and   edx, 1
	jnz   .nb103nf_dosingle
	jmp   .nb103nf_updateouterdata
.nb103nf_dosingle:
	mov   edx, [esp + nb103nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	

	mov esi, [ebp + nb103nf_charge]    ;# base of charge[] 
	xorpd xmm6, xmm6
	movlpd xmm6, [esi + eax*8]	;# jq A 
	
	movapd xmm3, [esp + nb103nf_iqM]
	movapd xmm4, [esp + nb103nf_iqH]
	mulsd xmm3, xmm6		;# qqM 
	mulsd xmm4, xmm6		;# qqH 
	
	movapd  [esp + nb103nf_qqM], xmm3
	movapd  [esp + nb103nf_qqH], xmm4	

	mov esi, [ebp + nb103nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move ixM-izM to xmm4-xmm6 
	movapd xmm4, [esp + nb103nf_ixM]
	movapd xmm5, [esp + nb103nf_iyM]
	movapd xmm6, [esp + nb103nf_izM]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqM in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb103nf_ixH1]
	movapd xmm5, [esp + nb103nf_iyH1]
	movapd xmm6, [esp + nb103nf_izH1]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb103nf_ixH2]
	movapd xmm4, [esp + nb103nf_iyH2]
	movapd xmm5, [esp + nb103nf_izH2]

	;# calc dr 
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2

	;# square it 
	mulsd xmm3,xmm3
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	addsd xmm5, xmm4
	addsd xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqM in xmm7 

	;# start with rsqM - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb103nf_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb103nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb103nf_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb103nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvM in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb103nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb103nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb103nf_three]
	subsd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb103nf_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb103nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb103nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb103nf_three]
	subsd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb103nf_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do M interactions 
	mulsd  xmm7, [esp + nb103nf_qqM]	;# xmm7=vcoul 
	addsd  xmm7, [esp + nb103nf_vctot]
	movlpd [esp + nb103nf_vctot], xmm7

	;# H1 interactions 
	mulsd  xmm6, [esp + nb103nf_qqH]	;# xmm6=vcoul 
	addsd  xmm6, [esp + nb103nf_vctot]
	movlpd [esp + nb103nf_vctot], xmm6

	;# H2 interactions 
	mulsd  xmm5, [esp + nb103nf_qqH]	;# xmm5=vcoul 
	addsd  xmm5, [esp + nb103nf_vctot]
	movlpd [esp + nb103nf_vctot], xmm5

.nb103nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb103nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb103nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	movapd xmm7, [esp + nb103nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	        
	;# add earlier value from mem 
	mov   eax, [ebp + nb103nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 	
	
        ;# finish if last 
        mov ecx, [esp + nb103nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb103nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb103nf_n], esi
        jmp .nb103nf_outer
.nb103nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb103nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb103nf_end
        ;# non-zero, do one more workunit
        jmp   .nb103nf_threadloop
.nb103nf_end:
	emms

	mov eax, [esp + nb103nf_nouter]
	mov ebx, [esp + nb103nf_ninner]
	mov ecx, [ebp + nb103nf_outeriter]
	mov edx, [ebp + nb103nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb103nf_salign]
	add esp, eax
	add esp, 296
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret

