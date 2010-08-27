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



.globl nb_kernel131_ia32_sse2
.globl _nb_kernel131_ia32_sse2
nb_kernel131_ia32_sse2:	
_nb_kernel131_ia32_sse2:	
.equiv          nb131_p_nri,            8
.equiv          nb131_iinr,             12
.equiv          nb131_jindex,           16
.equiv          nb131_jjnr,             20
.equiv          nb131_shift,            24
.equiv          nb131_shiftvec,         28
.equiv          nb131_fshift,           32
.equiv          nb131_gid,              36
.equiv          nb131_pos,              40
.equiv          nb131_faction,          44
.equiv          nb131_charge,           48
.equiv          nb131_p_facel,          52
.equiv          nb131_argkrf,           56
.equiv          nb131_argcrf,           60
.equiv          nb131_Vc,               64
.equiv          nb131_type,             68
.equiv          nb131_p_ntype,          72
.equiv          nb131_vdwparam,         76
.equiv          nb131_Vvdw,             80
.equiv          nb131_p_tabscale,       84
.equiv          nb131_VFtab,            88
.equiv          nb131_invsqrta,         92
.equiv          nb131_dvda,             96
.equiv          nb131_p_gbtabscale,     100
.equiv          nb131_GBtab,            104
.equiv          nb131_p_nthreads,       108
.equiv          nb131_count,            112
.equiv          nb131_mtx,              116
.equiv          nb131_outeriter,        120
.equiv          nb131_inneriter,        124
.equiv          nb131_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb131_ixO,              0
.equiv          nb131_iyO,              16
.equiv          nb131_izO,              32
.equiv          nb131_ixH1,             48
.equiv          nb131_iyH1,             64
.equiv          nb131_izH1,             80
.equiv          nb131_ixH2,             96
.equiv          nb131_iyH2,             112
.equiv          nb131_izH2,             128
.equiv          nb131_iqO,              144
.equiv          nb131_iqH,              160
.equiv          nb131_dxO,              176
.equiv          nb131_dyO,              192
.equiv          nb131_dzO,              208
.equiv          nb131_dxH1,             224
.equiv          nb131_dyH1,             240
.equiv          nb131_dzH1,             256
.equiv          nb131_dxH2,             272
.equiv          nb131_dyH2,             288
.equiv          nb131_dzH2,             304
.equiv          nb131_qqO,              320
.equiv          nb131_qqH,              336
.equiv          nb131_c6,               352
.equiv          nb131_c12,              368
.equiv          nb131_tsc,              384
.equiv          nb131_fstmp,            400
.equiv          nb131_vctot,            416
.equiv          nb131_Vvdwtot,          432
.equiv          nb131_fixO,             448
.equiv          nb131_fiyO,             464
.equiv          nb131_fizO,             480
.equiv          nb131_fixH1,            496
.equiv          nb131_fiyH1,            512
.equiv          nb131_fizH1,            528
.equiv          nb131_fixH2,            544
.equiv          nb131_fiyH2,            560
.equiv          nb131_fizH2,            576
.equiv          nb131_fjx,              592
.equiv          nb131_fjy,              608
.equiv          nb131_fjz,              624
.equiv          nb131_half,             640
.equiv          nb131_three,            656
.equiv          nb131_two,              672
.equiv          nb131_krsqO,            720
.equiv          nb131_krsqH1,           736
.equiv          nb131_krsqH2,           752
.equiv          nb131_rsqO,             768
.equiv          nb131_rinvH1,           784
.equiv          nb131_rinvH2,           800
.equiv          nb131_is3,              816
.equiv          nb131_ii3,              820
.equiv          nb131_ntia,             824
.equiv          nb131_innerjjnr,        828
.equiv          nb131_innerk,           832
.equiv          nb131_n,                836
.equiv          nb131_nn1,              840
.equiv          nb131_nri,              844
.equiv          nb131_nouter,           848
.equiv          nb131_ninner,           852
.equiv          nb131_salign,           856
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
	mov [esp + nb131_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb131_p_nri]
	mov ecx, [ecx]
	mov [esp + nb131_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb131_nouter], eax
	mov [esp + nb131_ninner], eax

	mov eax, [ebp + nb131_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb131_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb131_half], eax
	mov [esp + nb131_half+4], ebx
	movsd xmm1, [esp + nb131_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb131_half], xmm1
	movapd [esp + nb131_two], xmm2
	movapd [esp + nb131_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb131_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb131_charge]
	movsd xmm3, [edx + ebx*8]	
	movsd xmm4, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb131_p_facel]
	movsd xmm5, [esi]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb131_iqO], xmm3
	movapd [esp + nb131_iqH], xmm4
	
	mov   edx, [ebp + nb131_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov edi, [ebp + nb131_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	mov   [esp + nb131_ntia], ecx		
.nb131_threadloop:
        mov   esi, [ebp + nb131_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb131_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb131_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb131_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb131_n], eax
        mov [esp + nb131_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb131_outerstart
        jmp .nb131_end

.nb131_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb131_nouter]
	mov [esp + nb131_nouter], ebx

.nb131_outer:
	mov   eax, [ebp + nb131_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb131_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb131_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb131_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb131_pos]    ;# eax = base of pos[]  
	mov   [esp + nb131_ii3], ebx

	addsd xmm3, [eax + ebx*8]
	addsd xmm4, [eax + ebx*8 + 8]
	addsd xmm5, [eax + ebx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb131_ixO], xmm3
	movapd [esp + nb131_iyO], xmm4
	movapd [esp + nb131_izO], xmm5

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
	movapd [esp + nb131_ixH1], xmm0
	movapd [esp + nb131_iyH1], xmm1
	movapd [esp + nb131_izH1], xmm2
	movapd [esp + nb131_ixH2], xmm3
	movapd [esp + nb131_iyH2], xmm4
	movapd [esp + nb131_izH2], xmm5
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb131_vctot], xmm4
	movapd [esp + nb131_Vvdwtot], xmm4
	movapd [esp + nb131_fixO], xmm4
	movapd [esp + nb131_fiyO], xmm4
	movapd [esp + nb131_fizO], xmm4
	movapd [esp + nb131_fixH1], xmm4
	movapd [esp + nb131_fiyH1], xmm4
	movapd [esp + nb131_fizH1], xmm4
	movapd [esp + nb131_fixH2], xmm4
	movapd [esp + nb131_fiyH2], xmm4
	movapd [esp + nb131_fizH2], xmm4
	
	mov   eax, [ebp + nb131_jindex]
	mov   ecx, [eax+esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb131_pos]
	mov   edi, [ebp + nb131_faction]	
	mov   eax, [ebp + nb131_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb131_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb131_ninner]
	mov   [esp + nb131_ninner], ecx
	add   edx, 0
	mov   [esp + nb131_innerk], edx    ;# number of innerloop atoms 
	jge   .nb131_unroll_loop
	jmp   .nb131_checksingle
.nb131_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb131_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]

	add dword ptr [esp + nb131_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb131_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb131_iqO]
	mulpd  xmm4, [esp + nb131_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx

	movapd  [esp + nb131_qqO], xmm3
	movapd  [esp + nb131_qqH], xmm4
	
	mov esi, [ebp + nb131_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb131_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [esp + nb131_ntia]
	add eax, edi
	add ebx, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movlpd xmm7, [esi + ebx*8]	;# c6b
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 
	movhpd xmm7, [esi + ebx*8 + 8]	;# c6b c12b 
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [esp + nb131_c6], xmm4
	movapd [esp + nb131_c12], xmm6
	
	mov esi, [ebp + nb131_pos]       ;# base of pos[] 

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
	movapd xmm4, [esp + nb131_ixO]
	movapd xmm5, [esp + nb131_iyO]
	movapd xmm6, [esp + nb131_izO]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb131_dxO], xmm4
	movapd [esp + nb131_dyO], xmm5
	movapd [esp + nb131_dzO], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 
	movapd [esp + nb131_rsqO], xmm7
	
	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb131_ixH1]
	movapd xmm5, [esp + nb131_iyH1]
	movapd xmm6, [esp + nb131_izH1]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb131_dxH1], xmm4
	movapd [esp + nb131_dyH1], xmm5
	movapd [esp + nb131_dzH1], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm6, xmm5
	addpd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb131_ixH2]
	movapd xmm4, [esp + nb131_iyH2]
	movapd xmm5, [esp + nb131_izH2]

	;# calc dr 
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2

	;# store dr 
	movapd [esp + nb131_dxH2], xmm3
	movapd [esp + nb131_dyH2], xmm4
	movapd [esp + nb131_dzH2], xmm5
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
	movapd  xmm4, [esp + nb131_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb131_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb131_three]
	subpd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb131_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb131_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb131_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb131_three]
	subpd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb131_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	movapd  [esp + nb131_rinvH1], xmm6
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb131_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb131_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb131_three]
	subpd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb131_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 
	movapd  [esp + nb131_rinvH2], xmm5

	;# do O interactions 
	movapd xmm0, xmm7
	mulpd  xmm7, [esp + nb131_qqO] ;# vcoul
	movapd xmm6, xmm0
	mulpd  xmm6, xmm7 ;# vcoul*rinv

	movapd [esp + nb131_fstmp], xmm6 ;# save to temp. storage

	addpd  xmm7, [esp + nb131_vctot]
	movapd [esp + nb131_vctot], xmm7

	movapd xmm4, [esp + nb131_rsqO]	
	;# LJ table interaction. xmm0=rinv, xmm4=rsq
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb131_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  esi, [ebp + nb131_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6

	;# dispersion 
	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movlpd xmm3, [esi + ebx*8]	;# Y2 
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
	;# dispersion table ready, in xmm4-xmm7 	
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb131_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb131_c6]
	mulpd  xmm7, xmm4	 ;# fijD 
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addpd  xmm5, [esp + nb131_Vvdwtot]
	movapd xmm3, [esp + nb131_fstmp]
	mulpd  xmm7, [esp + nb131_tsc]
	subpd  xmm3, xmm7
	movapd [esp + nb131_fstmp], xmm3
	movapd [esp + nb131_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [esi + eax*8 + 32]	;# Y1 	
	movlpd xmm3, [esi + ebx*8 + 32]	;# Y2 
	movhpd xmm4, [esi + eax*8 + 40]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 40]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 48]	;# G1
	movlpd xmm3, [esi + ebx*8 + 48]	;# G2
	movhpd xmm6, [esi + eax*8 + 56]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 56]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	
	;# table ready, in xmm4-xmm7 	
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb131_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	
	movapd xmm4, [esp + nb131_c12]
	mulpd  xmm7, xmm4 
	mulpd  xmm5, xmm4  
	
	addpd  xmm5, [esp + nb131_Vvdwtot]
	movapd xmm3, [esp + nb131_fstmp]
	mulpd  xmm7, [esp + nb131_tsc]
	subpd  xmm3, xmm7
	movapd [esp + nb131_Vvdwtot], xmm5

	mulpd  xmm3, xmm0
		
	movapd xmm0, [esp + nb131_dxO]
	movapd xmm1, [esp + nb131_dyO]
	movapd xmm2, [esp + nb131_dzO]

	movd eax, mm0	
	movd ebx, mm1

	mov    edi, [ebp + nb131_faction]
	mulpd  xmm0, xmm3
	mulpd  xmm1, xmm3
	mulpd  xmm2, xmm3

	;# update O forces 
	movapd xmm3, [esp + nb131_fixO]
	movapd xmm4, [esp + nb131_fiyO]
	movapd xmm7, [esp + nb131_fizO]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb131_fixO], xmm3
	movapd [esp + nb131_fiyO], xmm4
	movapd [esp + nb131_fizO], xmm7
	;# update j forces with water O 
	movapd [esp + nb131_fjx], xmm0
	movapd [esp + nb131_fjy], xmm1
	movapd [esp + nb131_fjz], xmm2

	;# H1 interactions 
	movapd  xmm6, [esp + nb131_rinvH1]	
	movapd  xmm4, xmm6
	mulpd   xmm6, [esp + nb131_qqH] ;# vcoul 
	mulpd   xmm4, xmm4 ;# rinvsq
	mulpd   xmm4, xmm6 ;# vcoul*rinvsq
	
	addpd   xmm6, [esp + nb131_vctot]
	movapd [esp + nb131_vctot], xmm6

	movapd xmm0, [esp + nb131_dxH1]
	movapd xmm1, [esp + nb131_dyH1]
	movapd xmm2, [esp + nb131_dzH1]

	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [esp + nb131_fixH1]
	movapd xmm4, [esp + nb131_fiyH1]
	movapd xmm7, [esp + nb131_fizH1]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb131_fixH1], xmm3
	movapd [esp + nb131_fiyH1], xmm4
	movapd [esp + nb131_fizH1], xmm7
	;# update j forces with water H1 
	addpd  xmm0, [esp + nb131_fjx]
	addpd  xmm1, [esp + nb131_fjy]
	addpd  xmm2, [esp + nb131_fjz]
	movapd [esp + nb131_fjx], xmm0
	movapd [esp + nb131_fjy], xmm1
	movapd [esp + nb131_fjz], xmm2

	;# H2 interactions 
	movapd  xmm6, [esp + nb131_rinvH2]	
	movapd  xmm4, xmm6
	mulpd   xmm6, [esp + nb131_qqH] ;# vcoul 
	mulpd   xmm4, xmm4 ;# rinvsq
	mulpd   xmm4, xmm6 ;# vcoul*rinvsq

	addpd  xmm6, [esp + nb131_vctot]

	movapd xmm0, [esp + nb131_dxH2]
	movapd xmm1, [esp + nb131_dyH2]
	movapd xmm2, [esp + nb131_dzH2]
	movapd [esp + nb131_vctot], xmm6
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [esp + nb131_fixH2]
	movapd xmm4, [esp + nb131_fiyH2]
	movapd xmm7, [esp + nb131_fizH2]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb131_fixH2], xmm3
	movapd [esp + nb131_fiyH2], xmm4
	movapd [esp + nb131_fizH2], xmm7

	mov edi, [ebp + nb131_faction]
	;# update j forces 
	addpd  xmm0, [esp + nb131_fjx]
	addpd  xmm1, [esp + nb131_fjy]
	addpd  xmm2, [esp + nb131_fjz]
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
	sub dword ptr [esp + nb131_innerk],  2
	jl    .nb131_checksingle
	jmp   .nb131_unroll_loop
.nb131_checksingle:	
	mov   edx, [esp + nb131_innerk]
	and   edx, 1
	jnz   .nb131_dosingle
	jmp   .nb131_updateouterdata
.nb131_dosingle:
	mov   edx, [esp + nb131_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb131_innerjjnr],  4	

	mov esi, [ebp + nb131_charge]    ;# base of charge[] 
	xorpd xmm3, xmm3
	movlpd xmm3, [esi + eax*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb131_iqO]
	mulpd  xmm4, [esp + nb131_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 

	movapd  [esp + nb131_qqO], xmm3
	movapd  [esp + nb131_qqH], xmm4
	
	mov esi, [ebp + nb131_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb131_vdwparam]
	shl eax, 1	
	mov edi, [esp + nb131_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movapd [esp + nb131_c6], xmm4
	movapd [esp + nb131_c12], xmm6
	
	mov esi, [ebp + nb131_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [esp + nb131_ixO]
	movapd xmm5, [esp + nb131_iyO]
	movapd xmm6, [esp + nb131_izO]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb131_dxO], xmm4
	movapd [esp + nb131_dyO], xmm5
	movapd [esp + nb131_dzO], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 
	movapd [esp + nb131_rsqO], xmm7
	
	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb131_ixH1]
	movapd xmm5, [esp + nb131_iyH1]
	movapd xmm6, [esp + nb131_izH1]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb131_dxH1], xmm4
	movapd [esp + nb131_dyH1], xmm5
	movapd [esp + nb131_dzH1], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb131_ixH2]
	movapd xmm4, [esp + nb131_iyH2]
	movapd xmm5, [esp + nb131_izH2]

	;# calc dr 
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2

	;# store dr 
	movapd [esp + nb131_dxH2], xmm3
	movapd [esp + nb131_dyH2], xmm4
	movapd [esp + nb131_dzH2], xmm5
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
	movapd  xmm4, [esp + nb131_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb131_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb131_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb131_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb131_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb131_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb131_three]
	subsd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb131_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	movsd  [esp + nb131_rinvH1], xmm6
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb131_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb131_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb131_three]
	subsd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb131_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 
	movsd [esp + nb131_rinvH2], xmm5
	
	;# do O interactions 
	movsd xmm0, xmm7
	mulsd  xmm7, [esp + nb131_qqO] ;# vcoul
	movsd xmm6, xmm0
	mulsd  xmm6, xmm7 ;# vcoul*rinv

	movsd [esp + nb131_fstmp], xmm6 ;# save to temp. storage

	addsd  xmm7, [esp + nb131_vctot]
	movsd [esp + nb131_vctot], xmm7

	movsd xmm4, [esp + nb131_rsqO]

	;# LJ table interaction. xmm0=rinv, xmm4=rsq	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb131_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3
	
	mov  esi, [ebp + nb131_VFtab]

	;# dispersion 
	movlpd xmm4, [esi + ebx*8]	;# Y1 	
	movhpd xmm4, [esi + ebx*8 + 8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + ebx*8 + 16]	;# G1
	movhpd xmm6, [esi + ebx*8 + 24]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# dispersion table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb131_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movsd xmm4, [esp + nb131_c6]
	mulsd  xmm7, xmm4	 ;# fijD 
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addsd  xmm5, [esp + nb131_Vvdwtot]
	movsd xmm3, [esp + nb131_fstmp]
	mulsd  xmm7, [esp + nb131_tsc]
	subsd  xmm3, xmm7
	movsd [esp + nb131_fstmp], xmm3
	movsd [esp + nb131_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [esi + ebx*8 + 32]	;# Y1 	
	movhpd xmm4, [esi + ebx*8 + 40]	;# Y1 F1 	

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + ebx*8 + 48]	;# G1
	movhpd xmm6, [esi + ebx*8 + 56]	;# G1 H1 	

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	
	;# table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb131_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	
	movsd xmm4, [esp + nb131_c12]
	mulsd  xmm7, xmm4 
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [esp + nb131_Vvdwtot]
	movsd xmm3, [esp + nb131_fstmp]
	mulsd  xmm7, [esp + nb131_tsc]
	subsd  xmm3, xmm7
	movsd [esp + nb131_Vvdwtot], xmm5

	mulsd  xmm3, xmm0
		
	movsd xmm0, [esp + nb131_dxO]
	movsd xmm1, [esp + nb131_dyO]
	movsd xmm2, [esp + nb131_dzO]

	mov    edi, [ebp + nb131_faction]
	mulsd  xmm0, xmm3
	mulsd  xmm1, xmm3
	mulsd  xmm2, xmm3

	;# update O forces 
	movapd xmm3, [esp + nb131_fixO]
	movapd xmm4, [esp + nb131_fiyO]
	movapd xmm7, [esp + nb131_fizO]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [esp + nb131_fixO], xmm3
	movlpd [esp + nb131_fiyO], xmm4
	movlpd [esp + nb131_fizO], xmm7
	;# update j forces with water O 
	movlpd [esp + nb131_fjx], xmm0
	movlpd [esp + nb131_fjy], xmm1
	movlpd [esp + nb131_fjz], xmm2

	;# H1 interactions 
	movsd  xmm6, [esp + nb131_rinvH1]	
	movsd  xmm4, xmm6
	mulsd   xmm6, [esp + nb131_qqH] ;# vcoul 
	mulsd   xmm4, xmm4 ;# rinvsq
	mulsd   xmm4, xmm6 ;# vcoul*rinvsq

	addsd  xmm6, [esp + nb131_vctot]

	movapd xmm0, [esp + nb131_dxH1]
	movapd xmm1, [esp + nb131_dyH1]
	movapd xmm2, [esp + nb131_dzH1]
	movlpd [esp + nb131_vctot], xmm6
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [esp + nb131_fixH1]
	movapd xmm4, [esp + nb131_fiyH1]
	movapd xmm7, [esp + nb131_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [esp + nb131_fixH1], xmm3
	movlpd [esp + nb131_fiyH1], xmm4
	movlpd [esp + nb131_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [esp + nb131_fjx]
	addsd  xmm1, [esp + nb131_fjy]
	addsd  xmm2, [esp + nb131_fjz]
	movlpd [esp + nb131_fjx], xmm0
	movlpd [esp + nb131_fjy], xmm1
	movlpd [esp + nb131_fjz], xmm2

	;# H2 interactions 
	movsd  xmm6, [esp + nb131_rinvH2]	
	movsd  xmm4, xmm6
	mulsd   xmm6, [esp + nb131_qqH] ;# vcoul 
	mulsd   xmm4, xmm4 ;# rinvsq
	mulsd   xmm4, xmm6 ;# vcoul*rinvsq
	
	addsd  xmm6, [esp + nb131_vctot]

	movapd xmm0, [esp + nb131_dxH2]
	movapd xmm1, [esp + nb131_dyH2]
	movapd xmm2, [esp + nb131_dzH2]
	movlpd [esp + nb131_vctot], xmm6
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [esp + nb131_fixH2]
	movapd xmm4, [esp + nb131_fiyH2]
	movapd xmm7, [esp + nb131_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [esp + nb131_fixH2], xmm3
	movlpd [esp + nb131_fiyH2], xmm4
	movlpd [esp + nb131_fizH2], xmm7

	mov edi, [ebp + nb131_faction]
	;# update j forces 
	addsd  xmm0, [esp + nb131_fjx]
	addsd  xmm1, [esp + nb131_fjy]
	addsd  xmm2, [esp + nb131_fjz]
	movlpd xmm3, [edi + eax*8]
	movlpd xmm4, [edi + eax*8 + 8]
	movlpd xmm5, [edi + eax*8 + 16]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	movlpd [edi + eax*8], xmm3
	movlpd [edi + eax*8 + 8], xmm4
	movlpd [edi + eax*8 + 16], xmm5	

.nb131_updateouterdata:
	mov   ecx, [esp + nb131_ii3]
	mov   edi, [ebp + nb131_faction]
	mov   esi, [ebp + nb131_fshift]
	mov   edx, [esp + nb131_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb131_fixO]
	movapd xmm1, [esp + nb131_fiyO]
	movapd xmm2, [esp + nb131_fizO]

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
	movapd xmm0, [esp + nb131_fixH1]
	movapd xmm1, [esp + nb131_fiyH1]
	movapd xmm2, [esp + nb131_fizH1]

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
	movapd xmm0, [esp + nb131_fixH2]
	movapd xmm1, [esp + nb131_fiyH2]
	movapd xmm2, [esp + nb131_fizH2]

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
	mov esi, [esp + nb131_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb131_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb131_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb131_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb131_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb131_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [esp + nb131_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb131_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb131_n], esi
        jmp .nb131_outer
.nb131_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb131_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb131_end
        ;# non-zero, do one more workunit
        jmp   .nb131_threadloop
.nb131_end:
	emms

	mov eax, [esp + nb131_nouter]
	mov ebx, [esp + nb131_ninner]
	mov ecx, [ebp + nb131_outeriter]
	mov edx, [ebp + nb131_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb131_salign]
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





.globl nb_kernel131nf_ia32_sse2
.globl _nb_kernel131nf_ia32_sse2
nb_kernel131nf_ia32_sse2:	
_nb_kernel131nf_ia32_sse2:	
.equiv          nb131nf_p_nri,            8
.equiv          nb131nf_iinr,             12
.equiv          nb131nf_jindex,           16
.equiv          nb131nf_jjnr,             20
.equiv          nb131nf_shift,            24
.equiv          nb131nf_shiftvec,         28
.equiv          nb131nf_fshift,           32
.equiv          nb131nf_gid,              36
.equiv          nb131nf_pos,              40
.equiv          nb131nf_faction,          44
.equiv          nb131nf_charge,           48
.equiv          nb131nf_p_facel,          52
.equiv          nb131nf_argkrf,           56
.equiv          nb131nf_argcrf,           60
.equiv          nb131nf_Vc,               64
.equiv          nb131nf_type,             68
.equiv          nb131nf_p_ntype,          72
.equiv          nb131nf_vdwparam,         76
.equiv          nb131nf_Vvdw,             80
.equiv          nb131nf_p_tabscale,       84
.equiv          nb131nf_VFtab,            88
.equiv          nb131nf_invsqrta,         92
.equiv          nb131nf_dvda,             96
.equiv          nb131nf_p_gbtabscale,     100
.equiv          nb131nf_GBtab,            104
.equiv          nb131nf_p_nthreads,       108
.equiv          nb131nf_count,            112
.equiv          nb131nf_mtx,              116
.equiv          nb131nf_outeriter,        120
.equiv          nb131nf_inneriter,        124
.equiv          nb131nf_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb131nf_ixO,              0
.equiv          nb131nf_iyO,              16
.equiv          nb131nf_izO,              32
.equiv          nb131nf_ixH1,             48
.equiv          nb131nf_iyH1,             64
.equiv          nb131nf_izH1,             80
.equiv          nb131nf_ixH2,             96
.equiv          nb131nf_iyH2,             112
.equiv          nb131nf_izH2,             128
.equiv          nb131nf_iqO,              144
.equiv          nb131nf_iqH,              160
.equiv          nb131nf_qqO,              176
.equiv          nb131nf_qqH,              192
.equiv          nb131nf_c6,               208
.equiv          nb131nf_c12,              224
.equiv          nb131nf_tsc,              240
.equiv          nb131nf_vctot,            256
.equiv          nb131nf_Vvdwtot,          272
.equiv          nb131nf_half,             288
.equiv          nb131nf_three,            304
.equiv          nb131nf_two,              320
.equiv          nb131nf_krsqO,            368
.equiv          nb131nf_krsqH1,           384
.equiv          nb131nf_krsqH2,           400
.equiv          nb131nf_rsqO,             416
.equiv          nb131nf_rinvH1,           432
.equiv          nb131nf_rinvH2,           448
.equiv          nb131nf_is3,              464
.equiv          nb131nf_ii3,              468
.equiv          nb131nf_ntia,             472
.equiv          nb131nf_innerjjnr,        476
.equiv          nb131nf_innerk,           480
.equiv          nb131nf_n,                484
.equiv          nb131nf_nn1,              488
.equiv          nb131nf_nri,              492
.equiv          nb131nf_nouter,           496
.equiv          nb131nf_ninner,           500
.equiv          nb131nf_salign,           504
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 508		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb131nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb131nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb131nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb131nf_nouter], eax
	mov [esp + nb131nf_ninner], eax

	mov eax, [ebp + nb131nf_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb131nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb131nf_half], eax
	mov [esp + nb131nf_half+4], ebx
	movsd xmm1, [esp + nb131nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb131nf_half], xmm1
	movapd [esp + nb131nf_two], xmm2
	movapd [esp + nb131nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb131nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb131nf_charge]
	movsd xmm3, [edx + ebx*8]	
	movsd xmm4, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb131nf_p_facel]
	movsd xmm5, [esi]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb131nf_iqO], xmm3
	movapd [esp + nb131nf_iqH], xmm4
	
	mov   edx, [ebp + nb131nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov edi, [ebp + nb131nf_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	mov   [esp + nb131nf_ntia], ecx		
.nb131nf_threadloop:
        mov   esi, [ebp + nb131nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb131nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb131nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb131nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb131nf_n], eax
        mov [esp + nb131nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb131nf_outerstart
        jmp .nb131nf_end

.nb131nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb131nf_nouter]
	mov [esp + nb131nf_nouter], ebx

.nb131nf_outer:
	mov   eax, [ebp + nb131nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb131nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb131nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb131nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb131nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb131nf_ii3], ebx

	addsd xmm3, [eax + ebx*8]
	addsd xmm4, [eax + ebx*8 + 8]
	addsd xmm5, [eax + ebx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb131nf_ixO], xmm3
	movapd [esp + nb131nf_iyO], xmm4
	movapd [esp + nb131nf_izO], xmm5

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
	movapd [esp + nb131nf_ixH1], xmm0
	movapd [esp + nb131nf_iyH1], xmm1
	movapd [esp + nb131nf_izH1], xmm2
	movapd [esp + nb131nf_ixH2], xmm3
	movapd [esp + nb131nf_iyH2], xmm4
	movapd [esp + nb131nf_izH2], xmm5
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb131nf_vctot], xmm4
	movapd [esp + nb131nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb131nf_jindex]
	mov   ecx, [eax+esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb131nf_pos]
	mov   edi, [ebp + nb131nf_faction]	
	mov   eax, [ebp + nb131nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb131nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb131nf_ninner]
	mov   [esp + nb131nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb131nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb131nf_unroll_loop
	jmp   .nb131nf_checksingle
.nb131nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb131nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]

	add dword ptr [esp + nb131nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb131nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb131nf_iqO]
	mulpd  xmm4, [esp + nb131nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx

	movapd  [esp + nb131nf_qqO], xmm3
	movapd  [esp + nb131nf_qqH], xmm4
	
	mov esi, [ebp + nb131nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb131nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [esp + nb131nf_ntia]
	add eax, edi
	add ebx, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movlpd xmm7, [esi + ebx*8]	;# c6b
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 
	movhpd xmm7, [esi + ebx*8 + 8]	;# c6b c12b 
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [esp + nb131nf_c6], xmm4
	movapd [esp + nb131nf_c12], xmm6
	
	mov esi, [ebp + nb131nf_pos]       ;# base of pos[] 

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
	movapd xmm4, [esp + nb131nf_ixO]
	movapd xmm5, [esp + nb131nf_iyO]
	movapd xmm6, [esp + nb131nf_izO]

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
	movapd [esp + nb131nf_rsqO], xmm7
	
	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb131nf_ixH1]
	movapd xmm5, [esp + nb131nf_iyH1]
	movapd xmm6, [esp + nb131nf_izH1]

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
	movapd xmm3, [esp + nb131nf_ixH2]
	movapd xmm4, [esp + nb131nf_iyH2]
	movapd xmm5, [esp + nb131nf_izH2]

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
	movapd  xmm4, [esp + nb131nf_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb131nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb131nf_three]
	subpd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb131nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb131nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb131nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb131nf_three]
	subpd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb131nf_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	movapd  [esp + nb131nf_rinvH1], xmm6
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb131nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb131nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb131nf_three]
	subpd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb131nf_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 
	movapd  [esp + nb131nf_rinvH2], xmm5

	;# do O interactions 
	movapd xmm0, xmm7
	mulpd  xmm7, [esp + nb131nf_qqO]

	addpd  xmm7, [esp + nb131nf_vctot]
	movapd [esp + nb131nf_vctot], xmm7

	movapd xmm4, [esp + nb131nf_rsqO]
	;# LJ table interaction. xmm0=rinv, xmm4=rsq
	
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb131nf_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	mov  esi, [ebp + nb131nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6

	;# dispersion 
	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movlpd xmm3, [esi + ebx*8]	;# Y2 
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
	;# dispersion table ready, in xmm4-xmm7 	
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb131nf_c6]
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;#  Update Vvdwtot directly 
	addpd  xmm5, [esp + nb131nf_Vvdwtot]
	movapd [esp + nb131nf_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [esi + eax*8 + 32]	;# Y1 	
	movlpd xmm3, [esi + ebx*8 + 32]	;# Y2 
	movhpd xmm4, [esi + eax*8 + 40]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 40]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 48]	;# G1
	movlpd xmm3, [esi + ebx*8 + 48]	;# G2
	movhpd xmm6, [esi + eax*8 + 56]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 56]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	
	;# table ready, in xmm4-xmm7 	
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	
	movapd xmm4, [esp + nb131nf_c12]
	mulpd  xmm5, xmm4  
	
	addpd  xmm5, [esp + nb131nf_Vvdwtot]
	movapd [esp + nb131nf_Vvdwtot], xmm5

	;# H1 & H2 interactions 
	movapd  xmm6, [esp + nb131nf_rinvH1]	
	addpd   xmm6, [esp + nb131nf_rinvH2]
	mulpd   xmm6, [esp + nb131nf_qqH] ;# vcoul 
	
	addpd  xmm6, [esp + nb131nf_vctot]
	movapd [esp + nb131nf_vctot], xmm6

	;# should we do one more iteration? 
	sub dword ptr [esp + nb131nf_innerk],  2
	jl    .nb131nf_checksingle
	jmp   .nb131nf_unroll_loop
.nb131nf_checksingle:	
	mov   edx, [esp + nb131nf_innerk]
	and   edx, 1
	jnz   .nb131nf_dosingle
	jmp   .nb131nf_updateouterdata
.nb131nf_dosingle:
	mov   edx, [esp + nb131nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb131nf_innerjjnr],  4	

	mov esi, [ebp + nb131nf_charge]    ;# base of charge[] 
	xorpd xmm3, xmm3
	movlpd xmm3, [esi + eax*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb131nf_iqO]
	mulpd  xmm4, [esp + nb131nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 

	movapd  [esp + nb131nf_qqO], xmm3
	movapd  [esp + nb131nf_qqH], xmm4
	
	mov esi, [ebp + nb131nf_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb131nf_vdwparam]
	shl eax, 1	
	mov edi, [esp + nb131nf_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movapd [esp + nb131nf_c6], xmm4
	movapd [esp + nb131nf_c12], xmm6
	
	mov esi, [ebp + nb131nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [esp + nb131nf_ixO]
	movapd xmm5, [esp + nb131nf_iyO]
	movapd xmm6, [esp + nb131nf_izO]

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
	movapd [esp + nb131nf_rsqO], xmm7
	
	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb131nf_ixH1]
	movapd xmm5, [esp + nb131nf_iyH1]
	movapd xmm6, [esp + nb131nf_izH1]

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
	movapd xmm3, [esp + nb131nf_ixH2]
	movapd xmm4, [esp + nb131nf_iyH2]
	movapd xmm5, [esp + nb131nf_izH2]

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
	movapd  xmm4, [esp + nb131nf_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb131nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb131nf_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb131nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb131nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb131nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb131nf_three]
	subsd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb131nf_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	movsd  [esp + nb131nf_rinvH1], xmm6
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb131nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb131nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb131nf_three]
	subsd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb131nf_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 
	movsd [esp + nb131nf_rinvH2], xmm5
	
	;# do O interactions 
	movsd xmm0, xmm7
	mulsd  xmm7, [esp + nb131nf_qqO]

	addsd  xmm7, [esp + nb131nf_vctot]
	movsd [esp + nb131nf_vctot], xmm7

	movsd xmm4, [esp + nb131nf_rsqO]
	;# LJ table interaction. xmm0=rinv, xmm4=rsq
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb131nf_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3
	
	mov  esi, [ebp + nb131nf_VFtab]

	;# dispersion 
	movlpd xmm4, [esi + ebx*8]	;# Y1 	
	movhpd xmm4, [esi + ebx*8 + 8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + ebx*8 + 16]	;# G1
	movhpd xmm6, [esi + ebx*8 + 24]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# dispersion table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movsd xmm4, [esp + nb131nf_c6]
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# Update Vvdwtot directly 
	addsd  xmm5, [esp + nb131nf_Vvdwtot]
	movsd [esp + nb131nf_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [esi + ebx*8 + 32]	;# Y1 	
	movhpd xmm4, [esi + ebx*8 + 40]	;# Y1 F1 	

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + ebx*8 + 48]	;# G1
	movhpd xmm6, [esi + ebx*8 + 56]	;# G1 H1 	

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	
	;# table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	
	movsd xmm4, [esp + nb131nf_c12]
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [esp + nb131nf_Vvdwtot]
	movsd [esp + nb131nf_Vvdwtot], xmm5


	;# H1 & H2 interactions 
	movsd  xmm6, [esp + nb131nf_rinvH1]	
	addsd   xmm6, [esp + nb131nf_rinvH2]
	mulsd   xmm6, [esp + nb131nf_qqH] ;# vcoul 
	addsd   xmm6, [esp + nb131nf_vctot]
	movsd  [esp + nb131nf_vctot], xmm6
	
.nb131nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb131nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb131nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb131nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb131nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb131nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb131nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [esp + nb131nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb131nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb131nf_n], esi
        jmp .nb131nf_outer
.nb131nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb131nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb131nf_end
        ;# non-zero, do one more workunit
        jmp   .nb131nf_threadloop
.nb131nf_end:
	emms

	mov eax, [esp + nb131nf_nouter]
	mov ebx, [esp + nb131nf_ninner]
	mov ecx, [ebp + nb131nf_outeriter]
	mov edx, [ebp + nb131nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb131nf_salign]
	add esp, eax
	add esp, 508
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret


