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

	
.globl nb_kernel301_ia32_sse2
.globl _nb_kernel301_ia32_sse2
nb_kernel301_ia32_sse2:	
_nb_kernel301_ia32_sse2:	
.equiv          nb301_p_nri,            8
.equiv          nb301_iinr,             12
.equiv          nb301_jindex,           16
.equiv          nb301_jjnr,             20
.equiv          nb301_shift,            24
.equiv          nb301_shiftvec,         28
.equiv          nb301_fshift,           32
.equiv          nb301_gid,              36
.equiv          nb301_pos,              40
.equiv          nb301_faction,          44
.equiv          nb301_charge,           48
.equiv          nb301_p_facel,          52
.equiv          nb301_argkrf,           56
.equiv          nb301_argcrf,           60
.equiv          nb301_Vc,               64
.equiv          nb301_type,             68
.equiv          nb301_p_ntype,          72
.equiv          nb301_vdwparam,         76
.equiv          nb301_Vvdw,             80
.equiv          nb301_p_tabscale,       84
.equiv          nb301_VFtab,            88
.equiv          nb301_invsqrta,         92
.equiv          nb301_dvda,             96
.equiv          nb301_p_gbtabscale,     100
.equiv          nb301_GBtab,            104
.equiv          nb301_p_nthreads,       108
.equiv          nb301_count,            112
.equiv          nb301_mtx,              116
.equiv          nb301_outeriter,        120
.equiv          nb301_inneriter,        124
.equiv          nb301_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb301_ixO,              0
.equiv          nb301_iyO,              16
.equiv          nb301_izO,              32
.equiv          nb301_ixH1,             48
.equiv          nb301_iyH1,             64
.equiv          nb301_izH1,             80
.equiv          nb301_ixH2,             96
.equiv          nb301_iyH2,             112
.equiv          nb301_izH2,             128
.equiv          nb301_iqO,              144
.equiv          nb301_iqH,              160
.equiv          nb301_dxO,              176
.equiv          nb301_dyO,              192
.equiv          nb301_dzO,              208
.equiv          nb301_dxH1,             224
.equiv          nb301_dyH1,             240
.equiv          nb301_dzH1,             256
.equiv          nb301_dxH2,             272
.equiv          nb301_dyH2,             288
.equiv          nb301_dzH2,             304
.equiv          nb301_qqO,              320
.equiv          nb301_qqH,              336
.equiv          nb301_rinvO,            352
.equiv          nb301_rinvH1,           368
.equiv          nb301_rinvH2,           384
.equiv          nb301_rO,               400
.equiv          nb301_rH1,              416
.equiv          nb301_rH2,              432
.equiv          nb301_tsc,              448
.equiv          nb301_two,              464
.equiv          nb301_vctot,            480
.equiv          nb301_fixO,             496
.equiv          nb301_fiyO,             512
.equiv          nb301_fizO,             528
.equiv          nb301_fixH1,            544
.equiv          nb301_fiyH1,            560
.equiv          nb301_fizH1,            576
.equiv          nb301_fixH2,            592
.equiv          nb301_fiyH2,            608
.equiv          nb301_fizH2,            624
.equiv          nb301_fjx,              640
.equiv          nb301_fjy,              656
.equiv          nb301_fjz,              672
.equiv          nb301_half,             688
.equiv          nb301_three,            704
.equiv          nb301_is3,              720
.equiv          nb301_ii3,              724
.equiv          nb301_innerjjnr,        728
.equiv          nb301_innerk,           732
.equiv          nb301_n,                736
.equiv          nb301_nn1,              740
.equiv          nb301_nri,              744
.equiv          nb301_nouter,           748
.equiv          nb301_ninner,           752
.equiv          nb301_salign,           756
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 760		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb301_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb301_p_nri]
	mov ecx, [ecx]
	mov [esp + nb301_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb301_nouter], eax
	mov [esp + nb301_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb301_half], eax
	mov [esp + nb301_half+4], ebx
	movsd xmm1, [esp + nb301_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb301_half], xmm1
	movapd [esp + nb301_two], xmm2
	movapd [esp + nb301_three], xmm3

	mov eax, [ebp + nb301_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0 
	movapd [esp + nb301_tsc], xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb301_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb301_charge]
	movsd xmm3, [edx + ebx*8]	
	movsd xmm4, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb301_p_facel]
	movsd xmm5, [esi]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb301_iqO], xmm3
	movapd [esp + nb301_iqH], xmm4
	
.nb301_threadloop:
        mov   esi, [ebp + nb301_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb301_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
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
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb301_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb301_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb301_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb301_pos]    ;# eax = base of pos[]  
	mov   [esp + nb301_ii3], ebx

	addsd xmm3, [eax + ebx*8]
	addsd xmm4, [eax + ebx*8 + 8]
	addsd xmm5, [eax + ebx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb301_ixO], xmm3
	movapd [esp + nb301_iyO], xmm4
	movapd [esp + nb301_izO], xmm5

	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm0, [eax + ebx*8 + 24]
	addsd xmm1, [eax + ebx*8 + 32]
	addsd xmm2, [eax + ebx*8 + 40]		
	addsd xmm3, [eax + ebx*8 + 48]
	addsd xmm4, [eax + ebx*8 + 56]
	addsd xmm5, [eax + ebx*8 + 64]		

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb301_ixH1], xmm0
	movapd [esp + nb301_iyH1], xmm1
	movapd [esp + nb301_izH1], xmm2
	movapd [esp + nb301_ixH2], xmm3
	movapd [esp + nb301_iyH2], xmm4
	movapd [esp + nb301_izH2], xmm5
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb301_vctot], xmm4
	movapd [esp + nb301_fixO], xmm4
	movapd [esp + nb301_fiyO], xmm4
	movapd [esp + nb301_fizO], xmm4
	movapd [esp + nb301_fixH1], xmm4
	movapd [esp + nb301_fiyH1], xmm4
	movapd [esp + nb301_fizH1], xmm4
	movapd [esp + nb301_fixH2], xmm4
	movapd [esp + nb301_fiyH2], xmm4
	movapd [esp + nb301_fizH2], xmm4
	
	mov   eax, [ebp + nb301_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb301_pos]
	mov   edi, [ebp + nb301_faction]	
	mov   eax, [ebp + nb301_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb301_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb301_ninner]
	mov   [esp + nb301_ninner], ecx
	add   edx, 0
	mov   [esp + nb301_innerk], edx    ;# number of innerloop atoms 
	jge   .nb301_unroll_loop
	jmp   .nb301_checksingle
.nb301_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb301_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]

	add dword ptr [esp + nb301_innerjjnr],  8	;# advance pointer (unrolled 2) 
	mov esi, [ebp + nb301_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]
	movapd xmm4, xmm3	     
	mulpd  xmm3, [esp + nb301_iqO]
	mulpd  xmm4, [esp + nb301_iqH]

	movapd  [esp + nb301_qqO], xmm3
	movapd  [esp + nb301_qqH], xmm4	

	mov esi, [ebp + nb301_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [esp + nb301_ixO]
	movapd xmm5, [esp + nb301_iyO]
	movapd xmm6, [esp + nb301_izO]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb301_dxO], xmm4
	movapd [esp + nb301_dyO], xmm5
	movapd [esp + nb301_dzO], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb301_ixH1]
	movapd xmm5, [esp + nb301_iyH1]
	movapd xmm6, [esp + nb301_izH1]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb301_dxH1], xmm4
	movapd [esp + nb301_dyH1], xmm5
	movapd [esp + nb301_dzH1], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm6, xmm5
	addpd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb301_ixH2]
	movapd xmm4, [esp + nb301_iyH2]
	movapd xmm5, [esp + nb301_izH2]

	;# calc dr 
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2

	;# store dr 
	movapd [esp + nb301_dxH2], xmm3
	movapd [esp + nb301_dyH2], xmm4
	movapd [esp + nb301_dzH2], xmm5
	;# square it 
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	addpd xmm5, xmm4
	addpd xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqO - put seed in xmm2 
	cvtpd2ps xmm2, xmm7	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb301_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb301_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm7
	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb301_three]
	subpd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb301_half] ;# rinv 
	movapd  [esp + nb301_rinvO], xmm4	;# rinvO in xmm4 
	mulpd   xmm7, xmm4
	movapd  [esp + nb301_rO], xmm7	;# r in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb301_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb301_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm6
	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb301_three]
	subpd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb301_half] ;# rinv 
	movapd [esp + nb301_rinvH1], xmm4	;# rinvH1 
	mulpd  xmm6, xmm4
	movapd [esp + nb301_rH1], xmm6	;# rH1 
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb301_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb301_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm5
	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb301_three]
	subpd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb301_half] ;# rinv 
	movapd [esp + nb301_rinvH2], xmm4 ;# rinv 
	mulpd xmm5, xmm4
	movapd [esp + nb301_rH2], xmm5 ;# r 

	;# do O interactions 
	;# rO is still in xmm7 
	mulpd xmm7, [esp + nb301_tsc]
	cvttpd2pi mm6, xmm7	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	movd mm0, eax	
	movd mm1, ebx
	mov  esi, [ebp + nb301_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movlpd xmm4, [esi + eax*8]	;# Y1 F1 	
	movlpd xmm3, [esi + ebx*8]	;# Y2 F2 
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb301_two]	;# two*Heps2 
	movapd xmm3, [esp + nb301_qqO]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 
    ;# increment vcoul - then we can get rid of mm5 
    addpd  xmm5, [esp + nb301_vctot]
    movapd [esp + nb301_vctot], xmm5 
	xorpd  xmm4, xmm4

	mulpd  xmm3, [esp + nb301_tsc]
	mulpd  xmm3, [esp + nb301_rinvO]	
	subpd  xmm4, xmm3

	movapd xmm0, [esp + nb301_dxO]
	movapd xmm1, [esp + nb301_dyO]
	movapd xmm2, [esp + nb301_dzO]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4	;# tx in xmm0-xmm2 

	;# update O forces 
	movapd xmm3, [esp + nb301_fixO]
	movapd xmm4, [esp + nb301_fiyO]
	movapd xmm7, [esp + nb301_fizO]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb301_fixO], xmm3
	movapd [esp + nb301_fiyO], xmm4
	movapd [esp + nb301_fizO], xmm7
	;# update j forces with water O 
	movapd [esp + nb301_fjx], xmm0
	movapd [esp + nb301_fjy], xmm1
	movapd [esp + nb301_fjz], xmm2

	;# Done with O interactions - now H1! 
	movapd xmm7, [esp + nb301_rH1]
	mulpd xmm7, [esp + nb301_tsc]
	cvttpd2pi mm6, xmm7	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb301_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movlpd xmm4, [esi + eax*8]	;# Y1 F1 	
	movlpd xmm3, [esi + ebx*8]	;# Y2 F2 
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb301_two]	;# two*Heps2 
	movapd xmm3, [esp + nb301_qqH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 
    ;# increment vcoul 
	xorpd  xmm4, xmm4
    addpd  xmm5, [esp + nb301_vctot]
	mulpd  xmm3, [esp + nb301_rinvH1]
    movapd [esp + nb301_vctot], xmm5 
	mulpd  xmm3, [esp + nb301_tsc]
	subpd xmm4, xmm3

	movapd xmm0, [esp + nb301_dxH1]
	movapd xmm1, [esp + nb301_dyH1]
	movapd xmm2, [esp + nb301_dzH1]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [esp + nb301_fixH1]
	movapd xmm4, [esp + nb301_fiyH1]
	movapd xmm7, [esp + nb301_fizH1]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb301_fixH1], xmm3
	movapd [esp + nb301_fiyH1], xmm4
	movapd [esp + nb301_fizH1], xmm7
	;# update j forces with water H1 
	addpd  xmm0, [esp + nb301_fjx]
	addpd  xmm1, [esp + nb301_fjy]
	addpd  xmm2, [esp + nb301_fjz]
	movapd [esp + nb301_fjx], xmm0
	movapd [esp + nb301_fjy], xmm1
	movapd [esp + nb301_fjz], xmm2

	;# Done with H1, finally we do H2 interactions 
	movapd xmm7, [esp + nb301_rH2]
	mulpd   xmm7, [esp + nb301_tsc]
	cvttpd2pi mm6, xmm7	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb301_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movlpd xmm4, [esi + eax*8]	;# Y1 F1 	
	movlpd xmm3, [esi + ebx*8]	;# Y2 F2 
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb301_two]	;# two*Heps2 
	movapd xmm3, [esp + nb301_qqH]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 
    ;# increment vcoul 
	xorpd  xmm4, xmm4
    addpd  xmm5, [esp + nb301_vctot]
	mulpd  xmm3, [esp + nb301_rinvH2]
    movapd [esp + nb301_vctot], xmm5 
	mulpd  xmm3, [esp + nb301_tsc]
	subpd  xmm4, xmm3

	movapd xmm0, [esp + nb301_dxH2]
	movapd xmm1, [esp + nb301_dyH2]
	movapd xmm2, [esp + nb301_dzH2]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

    movd eax, mm0   
    movd ebx, mm1
	
	;# update H2 forces 
	movapd xmm3, [esp + nb301_fixH2]
	movapd xmm4, [esp + nb301_fiyH2]
	movapd xmm7, [esp + nb301_fizH2]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb301_fixH2], xmm3
	movapd [esp + nb301_fiyH2], xmm4
	movapd [esp + nb301_fizH2], xmm7

	mov edi, [ebp + nb301_faction]
	;# update j forces 
	;# update j forces with water H1 
	addpd  xmm0, [esp + nb301_fjx]
	addpd  xmm1, [esp + nb301_fjy]
	addpd  xmm2, [esp + nb301_fjz]

	;# the fj's - start by accumulating forces from memory 
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
	sub dword ptr [esp + nb301_innerk],  2
	jl    .nb301_checksingle
	jmp   .nb301_unroll_loop
.nb301_checksingle:	
	mov   edx, [esp + nb301_innerk]
	and   edx, 1
	jnz   .nb301_dosingle
	jmp   .nb301_updateouterdata
.nb301_dosingle:
	mov   edx, [esp + nb301_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	

	mov esi, [ebp + nb301_charge]    ;# base of charge[] 
	xorpd xmm3, xmm3
	movlpd xmm3, [esi + eax*8]
	movapd xmm4, xmm3	     
	mulpd  xmm3, [esp + nb301_iqO]
	mulpd  xmm4, [esp + nb301_iqH]

	movapd  [esp + nb301_qqO], xmm3
	movapd  [esp + nb301_qqH], xmm4	

	mov esi, [ebp + nb301_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	;# move coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [esp + nb301_ixO]
	movapd xmm5, [esp + nb301_iyO]
	movapd xmm6, [esp + nb301_izO]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb301_dxO], xmm4
	movapd [esp + nb301_dyO], xmm5
	movapd [esp + nb301_dzO], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb301_ixH1]
	movapd xmm5, [esp + nb301_iyH1]
	movapd xmm6, [esp + nb301_izH1]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb301_dxH1], xmm4
	movapd [esp + nb301_dyH1], xmm5
	movapd [esp + nb301_dzH1], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb301_ixH2]
	movapd xmm4, [esp + nb301_iyH2]
	movapd xmm5, [esp + nb301_izH2]

	;# calc dr 
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2

	;# store dr 
	movapd [esp + nb301_dxH2], xmm3
	movapd [esp + nb301_dyH2], xmm4
	movapd [esp + nb301_dzH2], xmm5
	;# square it 
	mulsd xmm3,xmm3
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	addsd xmm5, xmm4
	addsd xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb301_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb301_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm7
	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb301_three]
	subsd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb301_half] ;# rinv 
	movapd  [esp + nb301_rinvO], xmm4	;# rinvO in xmm4 
	mulsd   xmm7, xmm4
	movapd  [esp + nb301_rO], xmm7	;# r in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb301_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb301_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm6
	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb301_three]
	subsd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb301_half] ;# rinv 
	movapd [esp + nb301_rinvH1], xmm4	;# rinvH1 
	mulsd  xmm6, xmm4
	movapd [esp + nb301_rH1], xmm6	;# rH1 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb301_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb301_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm5
	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb301_three]
	subsd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb301_half] ;# rinv 
	movapd [esp + nb301_rinvH2], xmm4 ;# rinv 
	mulsd xmm5, xmm4
	movapd [esp + nb301_rH2], xmm5 ;# r 

	;# do O interactions 
	movd mm0, eax	
	;# rO is still in xmm7 
	mulsd   xmm7, [esp + nb301_tsc]
	cvttsd2si eax, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb301_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 F1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb301_two]	;# two*Heps2 
	movapd xmm3, [esp + nb301_qqO]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 
    ;# increment vcoul - then we can get rid of mm5 
    addsd  xmm5, [esp + nb301_vctot]
    movlpd [esp + nb301_vctot], xmm5 
	xorpd  xmm4, xmm4

	mulsd  xmm3, [esp + nb301_tsc]
	mulsd  xmm3, [esp + nb301_rinvO]	
	subsd  xmm4, xmm3

	movapd xmm0, [esp + nb301_dxO]
	movapd xmm1, [esp + nb301_dyO]
	movapd xmm2, [esp + nb301_dzO]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4	;# tx in xmm0-xmm2 

	;# update O forces 
	movapd xmm3, [esp + nb301_fixO]
	movapd xmm4, [esp + nb301_fiyO]
	movapd xmm7, [esp + nb301_fizO]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [esp + nb301_fixO], xmm3
	movlpd [esp + nb301_fiyO], xmm4
	movlpd [esp + nb301_fizO], xmm7
	;# update j forces with water O 
	movlpd [esp + nb301_fjx], xmm0
	movlpd [esp + nb301_fjy], xmm1
	movlpd [esp + nb301_fjz], xmm2

	;# Done with O interactions - now H1! 
	movapd xmm7, [esp + nb301_rH1]
	mulsd xmm7, [esp + nb301_tsc]
	cvttsd2si eax, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb301_VFtab]
	
	movlpd xmm4, [esi + eax*8]	;# Y1 F1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb301_two]	;# two*Heps2 
	movapd xmm3, [esp + nb301_qqH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 
    ;# increment vcoul 
	xorpd  xmm4, xmm4
    addsd  xmm5, [esp + nb301_vctot]
	mulsd  xmm3, [esp + nb301_rinvH1]
    movlpd [esp + nb301_vctot], xmm5 
	mulsd  xmm3, [esp + nb301_tsc]
	subsd xmm4, xmm3

	movapd xmm0, [esp + nb301_dxH1]
	movapd xmm1, [esp + nb301_dyH1]
	movapd xmm2, [esp + nb301_dzH1]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [esp + nb301_fixH1]
	movapd xmm4, [esp + nb301_fiyH1]
	movapd xmm7, [esp + nb301_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [esp + nb301_fixH1], xmm3
	movlpd [esp + nb301_fiyH1], xmm4
	movlpd [esp + nb301_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [esp + nb301_fjx]
	addsd  xmm1, [esp + nb301_fjy]
	addsd  xmm2, [esp + nb301_fjz]
	movlpd [esp + nb301_fjx], xmm0
	movlpd [esp + nb301_fjy], xmm1
	movlpd [esp + nb301_fjz], xmm2

	;# Done with H1, finally we do H2 interactions 
	movapd xmm7, [esp + nb301_rH2]
	mulsd   xmm7, [esp + nb301_tsc]
	cvttsd2si eax, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb301_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 F1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb301_two]	;# two*Heps2 
	movapd xmm3, [esp + nb301_qqH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 
    ;# increment vcoul 
	xorpd  xmm4, xmm4
    addsd  xmm5, [esp + nb301_vctot]
	mulsd  xmm3, [esp + nb301_rinvH2]
    movlpd [esp + nb301_vctot], xmm5 
	mulsd  xmm3, [esp + nb301_tsc]
	subsd  xmm4, xmm3

	movapd xmm0, [esp + nb301_dxH2]
	movapd xmm1, [esp + nb301_dyH2]
	movapd xmm2, [esp + nb301_dzH2]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

    movd eax, mm0   
	
	;# update H2 forces 
	movapd xmm3, [esp + nb301_fixH2]
	movapd xmm4, [esp + nb301_fiyH2]
	movapd xmm7, [esp + nb301_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [esp + nb301_fixH2], xmm3
	movlpd [esp + nb301_fiyH2], xmm4
	movlpd [esp + nb301_fizH2], xmm7

	mov edi, [ebp + nb301_faction]
	;# update j forces 
	;# update j forces with water H1 
	addsd  xmm0, [esp + nb301_fjx]
	addsd  xmm1, [esp + nb301_fjy]
	addsd  xmm2, [esp + nb301_fjz]

	;# the fj's - start by accumulating forces from memory 
	movlpd xmm3, [edi + eax*8]
	movlpd xmm4, [edi + eax*8 + 8]
	movlpd xmm5, [edi + eax*8 + 16]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	movlpd [edi + eax*8], xmm3
	movlpd [edi + eax*8 + 8], xmm4
	movlpd [edi + eax*8 + 16], xmm5

.nb301_updateouterdata:
	mov   ecx, [esp + nb301_ii3]
	mov   edi, [ebp + nb301_faction]
	mov   esi, [ebp + nb301_fshift]
	mov   edx, [esp + nb301_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb301_fixO]
	movapd xmm1, [esp + nb301_fiyO]
	movapd xmm2, [esp + nb301_fizO]

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
	movsd  xmm3, [edi + ecx*8]
	movsd  xmm4, [edi + ecx*8 + 8]
	movsd  xmm5, [edi + ecx*8 + 16]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movsd  [edi + ecx*8],     xmm3
	movsd  [edi + ecx*8 + 8], xmm4
	movsd  [edi + ecx*8 + 16], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movapd xmm6, xmm0
	movsd xmm7, xmm2
	unpcklpd xmm6, xmm1

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb301_fixH1]
	movapd xmm1, [esp + nb301_fiyH1]
	movapd xmm2, [esp + nb301_fizH1]

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
	movsd  [edi + ecx*8 + 24], xmm3
	movsd  [edi + ecx*8 + 32], xmm4
	movsd  [edi + ecx*8 + 40], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb301_fixH2]
	movapd xmm1, [esp + nb301_fiyH2]
	movapd xmm2, [esp + nb301_fizH2]

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
	mov esi, [esp + nb301_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb301_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb301_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb301_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb301_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb301_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb301_n], esi
        jmp .nb301_outer
.nb301_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb301_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb301_end
        ;# non-zero, do one more workunit
        jmp   .nb301_threadloop
.nb301_end:
	emms

	mov eax, [esp + nb301_nouter]
	mov ebx, [esp + nb301_ninner]
	mov ecx, [ebp + nb301_outeriter]
	mov edx, [ebp + nb301_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb301_salign]
	add esp, eax
	add esp, 760
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret
	


.globl nb_kernel301nf_ia32_sse2
.globl _nb_kernel301nf_ia32_sse2
nb_kernel301nf_ia32_sse2:	
_nb_kernel301nf_ia32_sse2:	
.equiv          nb301nf_p_nri,          8
.equiv          nb301nf_iinr,           12
.equiv          nb301nf_jindex,         16
.equiv          nb301nf_jjnr,           20
.equiv          nb301nf_shift,          24
.equiv          nb301nf_shiftvec,       28
.equiv          nb301nf_fshift,         32
.equiv          nb301nf_gid,            36
.equiv          nb301nf_pos,            40
.equiv          nb301nf_faction,        44
.equiv          nb301nf_charge,         48
.equiv          nb301nf_p_facel,        52
.equiv          nb301nf_argkrf,         56
.equiv          nb301nf_argcrf,         60
.equiv          nb301nf_Vc,             64
.equiv          nb301nf_type,           68
.equiv          nb301nf_p_ntype,        72
.equiv          nb301nf_vdwparam,       76
.equiv          nb301nf_Vvdw,           80
.equiv          nb301nf_p_tabscale,     84
.equiv          nb301nf_VFtab,          88
.equiv          nb301nf_invsqrta,       92
.equiv          nb301nf_dvda,           96
.equiv          nb301nf_p_gbtabscale,   100
.equiv          nb301nf_GBtab,          104
.equiv          nb301nf_p_nthreads,     108
.equiv          nb301nf_count,          112
.equiv          nb301nf_mtx,            116
.equiv          nb301nf_outeriter,      120
.equiv          nb301nf_inneriter,      124
.equiv          nb301nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb301nf_ixO,            0
.equiv          nb301nf_iyO,            16
.equiv          nb301nf_izO,            32
.equiv          nb301nf_ixH1,           48
.equiv          nb301nf_iyH1,           64
.equiv          nb301nf_izH1,           80
.equiv          nb301nf_ixH2,           96
.equiv          nb301nf_iyH2,           112
.equiv          nb301nf_izH2,           128
.equiv          nb301nf_iqO,            144
.equiv          nb301nf_iqH,            160
.equiv          nb301nf_qqO,            176
.equiv          nb301nf_qqH,            192
.equiv          nb301nf_rinvO,          208
.equiv          nb301nf_rinvH1,         224
.equiv          nb301nf_rinvH2,         240
.equiv          nb301nf_rO,             256
.equiv          nb301nf_rH1,            272
.equiv          nb301nf_rH2,            288
.equiv          nb301nf_tsc,            304
.equiv          nb301nf_vctot,          320
.equiv          nb301nf_half,           336
.equiv          nb301nf_three,          352
.equiv          nb301nf_is3,            368
.equiv          nb301nf_ii3,            372
.equiv          nb301nf_innerjjnr,      376
.equiv          nb301nf_innerk,         380
.equiv          nb301nf_n,              384
.equiv          nb301nf_nn1,            388
.equiv          nb301nf_nri,            392
.equiv          nb301nf_nouter,         396
.equiv          nb301nf_ninner,         400
.equiv          nb301nf_salign,         404
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 408		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb301nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb301nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb301nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb301nf_nouter], eax
	mov [esp + nb301nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb301nf_half], eax
	mov [esp + nb301nf_half+4], ebx
	movsd xmm1, [esp + nb301nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb301nf_half], xmm1
	movapd [esp + nb301nf_three], xmm3

	mov eax, [ebp + nb301nf_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0 
	movapd [esp + nb301nf_tsc], xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb301nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb301nf_charge]
	movsd xmm3, [edx + ebx*8]	
	movsd xmm4, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb301nf_p_facel]
	movsd xmm5, [esi]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb301nf_iqO], xmm3
	movapd [esp + nb301nf_iqH], xmm4
	
.nb301nf_threadloop:
        mov   esi, [ebp + nb301nf_count]        ;# pointer to sync counter
        mov   eax, [esi]
.nb301nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
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
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 

	mov   eax, [ebp + nb301nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb301nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb301nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb301nf_ii3], ebx

	addsd xmm3, [eax + ebx*8]
	addsd xmm4, [eax + ebx*8 + 8]
	addsd xmm5, [eax + ebx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb301nf_ixO], xmm3
	movapd [esp + nb301nf_iyO], xmm4
	movapd [esp + nb301nf_izO], xmm5

	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm0, [eax + ebx*8 + 24]
	addsd xmm1, [eax + ebx*8 + 32]
	addsd xmm2, [eax + ebx*8 + 40]		
	addsd xmm3, [eax + ebx*8 + 48]
	addsd xmm4, [eax + ebx*8 + 56]
	addsd xmm5, [eax + ebx*8 + 64]		

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb301nf_ixH1], xmm0
	movapd [esp + nb301nf_iyH1], xmm1
	movapd [esp + nb301nf_izH1], xmm2
	movapd [esp + nb301nf_ixH2], xmm3
	movapd [esp + nb301nf_iyH2], xmm4
	movapd [esp + nb301nf_izH2], xmm5
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb301nf_vctot], xmm4
	
	mov   eax, [ebp + nb301nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb301nf_pos]
	mov   eax, [ebp + nb301nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb301nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb301nf_ninner]
	mov   [esp + nb301nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb301nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb301nf_unroll_loop
	jmp   .nb301nf_checksingle
.nb301nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb301nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]

	add dword ptr [esp + nb301nf_innerjjnr],  8	;# advance pointer (unrolled 2) 
	mov esi, [ebp + nb301nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]
	movapd xmm4, xmm3	     
	mulpd  xmm3, [esp + nb301nf_iqO]
	mulpd  xmm4, [esp + nb301nf_iqH]

	movapd  [esp + nb301nf_qqO], xmm3
	movapd  [esp + nb301nf_qqH], xmm4	

	mov esi, [ebp + nb301nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [esp + nb301nf_ixO]
	movapd xmm5, [esp + nb301nf_iyO]
	movapd xmm6, [esp + nb301nf_izO]

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
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb301nf_ixH1]
	movapd xmm5, [esp + nb301nf_iyH1]
	movapd xmm6, [esp + nb301nf_izH1]

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
	movapd xmm3, [esp + nb301nf_ixH2]
	movapd xmm4, [esp + nb301nf_iyH2]
	movapd xmm5, [esp + nb301nf_izH2]

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
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqO - put seed in xmm2 
	cvtpd2ps xmm2, xmm7	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb301nf_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb301nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm7
	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb301nf_three]
	subpd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb301nf_half] ;# rinv 
	movapd  [esp + nb301nf_rinvO], xmm4	;# rinvO in xmm4 
	mulpd   xmm7, xmm4
	movapd  [esp + nb301nf_rO], xmm7	;# r in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb301nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb301nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm6
	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb301nf_three]
	subpd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb301nf_half] ;# rinv 
	movapd [esp + nb301nf_rinvH1], xmm4	;# rinvH1 
	mulpd  xmm6, xmm4
	movapd [esp + nb301nf_rH1], xmm6	;# rH1 
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb301nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb301nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm5
	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb301nf_three]
	subpd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb301nf_half] ;# rinv 
	movapd [esp + nb301nf_rinvH2], xmm4 ;# rinv 
	mulpd xmm5, xmm4
	movapd [esp + nb301nf_rH2], xmm5 ;# r 

	;# do O interactions 
	;# rO is still in xmm7 
	mulpd xmm7, [esp + nb301nf_tsc]
	cvttpd2pi mm6, xmm7	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb301nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movlpd xmm4, [esi + eax*8]	;# Y1 F1 	
	movlpd xmm3, [esi + ebx*8]	;# Y2 F2 
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb301nf_qqO]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
    ;# increment vcoul - then we can get rid of mm5 
    addpd  xmm5, [esp + nb301nf_vctot]
    movapd [esp + nb301nf_vctot], xmm5 

	;# Done with O interactions - now H1! 
	movapd xmm7, [esp + nb301nf_rH1]
	mulpd xmm7, [esp + nb301nf_tsc]
	cvttpd2pi mm6, xmm7	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb301nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movlpd xmm4, [esi + eax*8]	;# Y1 F1 	
	movlpd xmm3, [esi + ebx*8]	;# Y2 F2 
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb301nf_qqH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
    ;# increment vcoul 
    addpd  xmm5, [esp + nb301nf_vctot]
	movapd [esp + nb301nf_vctot], xmm5
	
	;# Done with H1, finally we do H2 interactions 
	movapd xmm7, [esp + nb301nf_rH2]
	mulpd   xmm7, [esp + nb301nf_tsc]
	cvttpd2pi mm6, xmm7	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  esi, [ebp + nb301nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movlpd xmm4, [esi + eax*8]	;# Y1 F1 	
	movlpd xmm3, [esi + ebx*8]	;# Y2 F2 
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb301nf_qqH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
    ;# increment vcoul 
    addpd  xmm5, [esp + nb301nf_vctot]
	movapd [esp + nb301nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [esp + nb301nf_innerk],  2
	jl    .nb301nf_checksingle
	jmp   .nb301nf_unroll_loop
.nb301nf_checksingle:	
	mov   edx, [esp + nb301nf_innerk]
	and   edx, 1
	jnz   .nb301nf_dosingle
	jmp   .nb301nf_updateouterdata
.nb301nf_dosingle:
	mov   edx, [esp + nb301nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	

	mov esi, [ebp + nb301nf_charge]    ;# base of charge[] 
	xorpd xmm3, xmm3
	movlpd xmm3, [esi + eax*8]
	movapd xmm4, xmm3	     
	mulpd  xmm3, [esp + nb301nf_iqO]
	mulpd  xmm4, [esp + nb301nf_iqH]

	movapd  [esp + nb301nf_qqO], xmm3
	movapd  [esp + nb301nf_qqH], xmm4	

	mov esi, [ebp + nb301nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	;# move coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [esp + nb301nf_ixO]
	movapd xmm5, [esp + nb301nf_iyO]
	movapd xmm6, [esp + nb301nf_izO]

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
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb301nf_ixH1]
	movapd xmm5, [esp + nb301nf_iyH1]
	movapd xmm6, [esp + nb301nf_izH1]

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
	movapd xmm3, [esp + nb301nf_ixH2]
	movapd xmm4, [esp + nb301nf_iyH2]
	movapd xmm5, [esp + nb301nf_izH2]

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
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb301nf_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb301nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm7
	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb301nf_three]
	subsd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb301nf_half] ;# rinv 
	movapd  [esp + nb301nf_rinvO], xmm4	;# rinvO in xmm4 
	mulsd   xmm7, xmm4
	movapd  [esp + nb301nf_rO], xmm7	;# r in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb301nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb301nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm6
	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb301nf_three]
	subsd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb301nf_half] ;# rinv 
	movapd [esp + nb301nf_rinvH1], xmm4	;# rinvH1 
	mulsd  xmm6, xmm4
	movapd [esp + nb301nf_rH1], xmm6	;# rH1 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb301nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb301nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm5
	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb301nf_three]
	subsd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb301nf_half] ;# rinv 
	movapd [esp + nb301nf_rinvH2], xmm4 ;# rinv 
	mulsd xmm5, xmm4
	movapd [esp + nb301nf_rH2], xmm5 ;# r 

	;# do O interactions 
	movd mm0, eax	
	;# rO is still in xmm7 
	mulsd   xmm7, [esp + nb301nf_tsc]
	cvttsd2si eax, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb301nf_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 F1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb301nf_qqO]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
    ;# increment vcoul - then we can get rid of mm5 
    addsd  xmm5, [esp + nb301nf_vctot]
    movlpd [esp + nb301nf_vctot], xmm5 

	;# Done with O interactions - now H1! 
	movapd xmm7, [esp + nb301nf_rH1]
	mulsd xmm7, [esp + nb301nf_tsc]
	cvttsd2si eax, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb301nf_VFtab]
	
	movlpd xmm4, [esi + eax*8]	;# Y1 F1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb301nf_qqH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
    ;# increment vcoul 
    addsd  xmm5, [esp + nb301nf_vctot]
    movlpd [esp + nb301nf_vctot], xmm5 


	;# Done with H1, finally we do H2 interactions 
	movapd xmm7, [esp + nb301nf_rH2]
	mulsd   xmm7, [esp + nb301nf_tsc]
	cvttsd2si eax, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  esi, [ebp + nb301nf_VFtab]

	movlpd xmm4, [esi + eax*8]	;# Y1 F1 	
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	

	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb301nf_qqH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
    ;# increment vcoul 
    addsd  xmm5, [esp + nb301nf_vctot]
    movlpd [esp + nb301nf_vctot], xmm5 

.nb301nf_updateouterdata:
	;# get group index for i particle 
	;# get n from stack
	mov esi, [esp + nb301nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb301nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb301nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb301nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb301nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb301nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb301nf_n], esi
        jmp .nb301nf_outer
.nb301nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb301nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb301nf_end
        ;# non-zero, do one more workunit
        jmp   .nb301nf_threadloop
.nb301nf_end:
	emms

	mov eax, [esp + nb301nf_nouter]
	mov ebx, [esp + nb301nf_ninner]
	mov ecx, [ebp + nb301nf_outeriter]
	mov edx, [ebp + nb301nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb301nf_salign]
	add esp, eax
	add esp, 408
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret
	
