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



.globl nb_kernel231_ia32_sse2
.globl _nb_kernel231_ia32_sse2
nb_kernel231_ia32_sse2:	
_nb_kernel231_ia32_sse2:	
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
	;# bottom of stack is cache-aligned for sse2 use 
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
.equiv          nb231_rsqO,             768
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
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb231_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb231_half], eax
	mov [esp + nb231_half+4], ebx
	movsd xmm1, [esp + nb231_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb231_half], xmm1
	movapd [esp + nb231_two], xmm2
	movapd [esp + nb231_three], xmm3

	mov esi, [ebp + nb231_argkrf]
	mov edi, [ebp + nb231_argcrf]
	movsd xmm5, [esi]
	movsd xmm6, [edi]
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	movapd [esp + nb231_krf], xmm5
	movapd [esp + nb231_crf], xmm6
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb231_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb231_charge]
	movsd xmm3, [edx + ebx*8]	
	movsd xmm4, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb231_p_facel]
	movsd xmm5, [esi]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb231_iqO], xmm3
	movapd [esp + nb231_iqH], xmm4
	
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
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb231_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb231_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb231_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb231_pos]    ;# eax = base of pos[]  
	mov   [esp + nb231_ii3], ebx

	addsd xmm3, [eax + ebx*8]
	addsd xmm4, [eax + ebx*8 + 8]
	addsd xmm5, [eax + ebx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb231_ixO], xmm3
	movapd [esp + nb231_iyO], xmm4
	movapd [esp + nb231_izO], xmm5

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
	movapd [esp + nb231_ixH1], xmm0
	movapd [esp + nb231_iyH1], xmm1
	movapd [esp + nb231_izH1], xmm2
	movapd [esp + nb231_ixH2], xmm3
	movapd [esp + nb231_iyH2], xmm4
	movapd [esp + nb231_izH2], xmm5
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb231_vctot], xmm4
	movapd [esp + nb231_Vvdwtot], xmm4
	movapd [esp + nb231_fixO], xmm4
	movapd [esp + nb231_fiyO], xmm4
	movapd [esp + nb231_fizO], xmm4
	movapd [esp + nb231_fixH1], xmm4
	movapd [esp + nb231_fiyH1], xmm4
	movapd [esp + nb231_fizH1], xmm4
	movapd [esp + nb231_fixH2], xmm4
	movapd [esp + nb231_fiyH2], xmm4
	movapd [esp + nb231_fizH2], xmm4
	
	mov   eax, [ebp + nb231_jindex]
	mov   ecx, [eax+esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb231_pos]
	mov   edi, [ebp + nb231_faction]	
	mov   eax, [ebp + nb231_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb231_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb231_ninner]
	mov   [esp + nb231_ninner], ecx
	add   edx, 0
	mov   [esp + nb231_innerk], edx    ;# number of innerloop atoms 
	jge   .nb231_unroll_loop
	jmp   .nb231_checksingle
.nb231_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb231_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]

	add dword ptr [esp + nb231_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb231_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb231_iqO]
	mulpd  xmm4, [esp + nb231_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx

	movapd  [esp + nb231_qqO], xmm3
	movapd  [esp + nb231_qqH], xmm4
	
	mov esi, [ebp + nb231_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb231_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [esp + nb231_ntia]
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
	movapd [esp + nb231_c6], xmm4
	movapd [esp + nb231_c12], xmm6
	
	mov esi, [ebp + nb231_pos]       ;# base of pos[] 

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
	movapd xmm4, [esp + nb231_ixO]
	movapd xmm5, [esp + nb231_iyO]
	movapd xmm6, [esp + nb231_izO]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb231_dxO], xmm4
	movapd [esp + nb231_dyO], xmm5
	movapd [esp + nb231_dzO], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 
	movapd [esp + nb231_rsqO], xmm7
	
	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb231_ixH1]
	movapd xmm5, [esp + nb231_iyH1]
	movapd xmm6, [esp + nb231_izH1]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb231_dxH1], xmm4
	movapd [esp + nb231_dyH1], xmm5
	movapd [esp + nb231_dzH1], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm6, xmm5
	addpd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb231_ixH2]
	movapd xmm4, [esp + nb231_iyH2]
	movapd xmm5, [esp + nb231_izH2]

	;# calc dr 
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2

	;# store dr 
	movapd [esp + nb231_dxH2], xmm3
	movapd [esp + nb231_dyH2], xmm4
	movapd [esp + nb231_dzH2], xmm5
	;# square it 
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	addpd xmm5, xmm4
	addpd xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	movapd xmm0, xmm5
	movapd xmm1, xmm6
	movapd xmm2, xmm7

	mulpd  xmm0, [esp + nb231_krf]	
	mulpd  xmm1, [esp + nb231_krf]	
	mulpd  xmm2, [esp + nb231_krf]	

	movapd [esp + nb231_krsqH2], xmm0
	movapd [esp + nb231_krsqH1], xmm1
	movapd [esp + nb231_krsqO], xmm2
		
	;# start with rsqO - put seed in xmm2 
	cvtpd2ps xmm2, xmm7	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb231_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb231_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb231_three]
	subpd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb231_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb231_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb231_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb231_three]
	subpd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb231_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	movapd  [esp + nb231_rinvH1], xmm6
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb231_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb231_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb231_three]
	subpd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb231_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 
	movapd  [esp + nb231_rinvH2], xmm5

	;# do O interactions 
	movapd xmm0, xmm7
	movapd xmm6, xmm7
	movapd xmm1, [esp + nb231_krsqO]
	addpd  xmm0, xmm1
	mulpd  xmm1, [esp + nb231_two]
	subpd  xmm0, [esp + nb231_crf] ;# xmm0=rinv+ krsq-crf 
	subpd  xmm6, xmm1
	mulpd  xmm0, [esp + nb231_qqO]
	mulpd  xmm6, [esp + nb231_qqO]

	mulpd  xmm6, xmm7
	movapd [esp + nb231_fstmp], xmm6 ;# save to temp. storage

	addpd  xmm0, [esp + nb231_vctot]
	movapd [esp + nb231_vctot], xmm0

	movapd xmm0, xmm7
	movapd xmm4, [esp + nb231_rsqO]
	;# LJ table interaction. xmm0=rinv, xmm4=rsq
	
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb231_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  esi, [ebp + nb231_VFtab]
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
	mulpd  xmm7, [esp + nb231_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb231_c6]
	mulpd  xmm7, xmm4	 ;# fijD 
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addpd  xmm5, [esp + nb231_Vvdwtot]
	movapd xmm3, [esp + nb231_fstmp]
	mulpd  xmm7, [esp + nb231_tsc]
	subpd  xmm3, xmm7
	movapd [esp + nb231_fstmp], xmm3
	movapd [esp + nb231_Vvdwtot], xmm5

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
	mulpd  xmm7, [esp + nb231_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	
	movapd xmm4, [esp + nb231_c12]
	mulpd  xmm7, xmm4 
	mulpd  xmm5, xmm4  
	
	addpd  xmm5, [esp + nb231_Vvdwtot]
	movapd xmm3, [esp + nb231_fstmp]
	mulpd  xmm7, [esp + nb231_tsc]
	subpd  xmm3, xmm7
	movapd [esp + nb231_Vvdwtot], xmm5

	mulpd  xmm3, xmm0
		
	movapd xmm0, [esp + nb231_dxO]
	movapd xmm1, [esp + nb231_dyO]
	movapd xmm2, [esp + nb231_dzO]

	movd eax, mm0	
	movd ebx, mm1

	mov    edi, [ebp + nb231_faction]
	mulpd  xmm0, xmm3
	mulpd  xmm1, xmm3
	mulpd  xmm2, xmm3

	;# update O forces 
	movapd xmm3, [esp + nb231_fixO]
	movapd xmm4, [esp + nb231_fiyO]
	movapd xmm7, [esp + nb231_fizO]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb231_fixO], xmm3
	movapd [esp + nb231_fiyO], xmm4
	movapd [esp + nb231_fizO], xmm7
	;# update j forces with water O 
	movapd [esp + nb231_fjx], xmm0
	movapd [esp + nb231_fjy], xmm1
	movapd [esp + nb231_fjz], xmm2

	;# H1 interactions 
	movapd  xmm6, [esp + nb231_rinvH1]	
	movapd  xmm4, xmm6
	mulpd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm6
	movapd  xmm0, [esp + nb231_krsqH1]
	addpd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	mulpd   xmm0, [esp + nb231_two]
	subpd   xmm6, [esp + nb231_crf]
	subpd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulpd   xmm6, [esp + nb231_qqH] ;# vcoul 
	mulpd   xmm7, [esp + nb231_qqH]
	mulpd  xmm4, xmm7		;# total fsH1 in xmm4 
	
	addpd  xmm6, [esp + nb231_vctot]

	movapd xmm0, [esp + nb231_dxH1]
	movapd xmm1, [esp + nb231_dyH1]
	movapd xmm2, [esp + nb231_dzH1]
	movapd [esp + nb231_vctot], xmm6
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [esp + nb231_fixH1]
	movapd xmm4, [esp + nb231_fiyH1]
	movapd xmm7, [esp + nb231_fizH1]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb231_fixH1], xmm3
	movapd [esp + nb231_fiyH1], xmm4
	movapd [esp + nb231_fizH1], xmm7
	;# update j forces with water H1 
	addpd  xmm0, [esp + nb231_fjx]
	addpd  xmm1, [esp + nb231_fjy]
	addpd  xmm2, [esp + nb231_fjz]
	movapd [esp + nb231_fjx], xmm0
	movapd [esp + nb231_fjy], xmm1
	movapd [esp + nb231_fjz], xmm2

	;# H2 interactions 
	movapd  xmm5, [esp + nb231_rinvH2]	
	movapd  xmm4, xmm5
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm5
	movapd  xmm0, [esp + nb231_krsqH2]
	addpd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulpd   xmm0, [esp + nb231_two]
	subpd   xmm5, [esp + nb231_crf]
	subpd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulpd   xmm5, [esp + nb231_qqH] ;# vcoul 
	mulpd   xmm7, [esp + nb231_qqH]
	mulpd  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addpd  xmm5, [esp + nb231_vctot]

	movapd xmm0, [esp + nb231_dxH2]
	movapd xmm1, [esp + nb231_dyH2]
	movapd xmm2, [esp + nb231_dzH2]
	movapd [esp + nb231_vctot], xmm5
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [esp + nb231_fixH2]
	movapd xmm4, [esp + nb231_fiyH2]
	movapd xmm7, [esp + nb231_fizH2]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb231_fixH2], xmm3
	movapd [esp + nb231_fiyH2], xmm4
	movapd [esp + nb231_fizH2], xmm7

	mov edi, [ebp + nb231_faction]
	;# update j forces 
	addpd  xmm0, [esp + nb231_fjx]
	addpd  xmm1, [esp + nb231_fjy]
	addpd  xmm2, [esp + nb231_fjz]
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
	sub dword ptr [esp + nb231_innerk],  2
	jl    .nb231_checksingle
	jmp   .nb231_unroll_loop
.nb231_checksingle:	
	mov   edx, [esp + nb231_innerk]
	and   edx, 1
	jnz   .nb231_dosingle
	jmp   .nb231_updateouterdata
.nb231_dosingle:
	mov   edx, [esp + nb231_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb231_innerjjnr],  4	

	mov esi, [ebp + nb231_charge]    ;# base of charge[] 
	xorpd xmm3, xmm3
	movlpd xmm3, [esi + eax*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb231_iqO]
	mulpd  xmm4, [esp + nb231_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 

	movapd  [esp + nb231_qqO], xmm3
	movapd  [esp + nb231_qqH], xmm4
	
	mov esi, [ebp + nb231_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb231_vdwparam]
	shl eax, 1	
	mov edi, [esp + nb231_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movapd [esp + nb231_c6], xmm4
	movapd [esp + nb231_c12], xmm6
	
	mov esi, [ebp + nb231_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [esp + nb231_ixO]
	movapd xmm5, [esp + nb231_iyO]
	movapd xmm6, [esp + nb231_izO]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb231_dxO], xmm4
	movapd [esp + nb231_dyO], xmm5
	movapd [esp + nb231_dzO], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 
	movapd [esp + nb231_rsqO], xmm7
	
	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb231_ixH1]
	movapd xmm5, [esp + nb231_iyH1]
	movapd xmm6, [esp + nb231_izH1]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb231_dxH1], xmm4
	movapd [esp + nb231_dyH1], xmm5
	movapd [esp + nb231_dzH1], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb231_ixH2]
	movapd xmm4, [esp + nb231_iyH2]
	movapd xmm5, [esp + nb231_izH2]

	;# calc dr 
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2

	;# store dr 
	movapd [esp + nb231_dxH2], xmm3
	movapd [esp + nb231_dyH2], xmm4
	movapd [esp + nb231_dzH2], xmm5
	;# square it 
	mulsd xmm3,xmm3
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	addsd xmm5, xmm4
	addsd xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	movapd xmm0, xmm5
	movapd xmm1, xmm6
	movapd xmm2, xmm7

	mulsd  xmm0, [esp + nb231_krf]	
	mulsd  xmm1, [esp + nb231_krf]	
	mulsd  xmm2, [esp + nb231_krf]	

	movapd [esp + nb231_krsqH2], xmm0
	movapd [esp + nb231_krsqH1], xmm1
	movapd [esp + nb231_krsqO], xmm2

	;# start with rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb231_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb231_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb231_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb231_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb231_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb231_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb231_three]
	subsd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb231_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	movsd  [esp + nb231_rinvH1], xmm6
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb231_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb231_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb231_three]
	subsd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb231_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 
	movsd [esp + nb231_rinvH2], xmm5
	
	;# do O interactions 
	movsd xmm0, xmm7
	movsd xmm6, xmm7
	movsd xmm1, [esp + nb231_krsqO]
	addsd  xmm0, xmm1
	mulsd  xmm1, [esp + nb231_two]
	subsd  xmm0, [esp + nb231_crf] ;# xmm0=rinv+ krsq-crf 
	subsd  xmm6, xmm1
	mulsd  xmm0, [esp + nb231_qqO]
	mulsd  xmm6, [esp + nb231_qqO]

	mulsd  xmm6, xmm7
	movsd [esp + nb231_fstmp], xmm6 ;# save to temp. storage

	addsd  xmm0, [esp + nb231_vctot]
	movsd [esp + nb231_vctot], xmm0

	movsd xmm0, xmm7
	movsd xmm4, [esp + nb231_rsqO]
	;# LJ table interaction. xmm0=rinv, xmm4=rsq
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb231_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3
	
	mov  esi, [ebp + nb231_VFtab]

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
	mulsd  xmm7, [esp + nb231_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movsd xmm4, [esp + nb231_c6]
	mulsd  xmm7, xmm4	 ;# fijD 
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addsd  xmm5, [esp + nb231_Vvdwtot]
	movsd xmm3, [esp + nb231_fstmp]
	mulsd  xmm7, [esp + nb231_tsc]
	subsd  xmm3, xmm7
	movsd [esp + nb231_fstmp], xmm3
	movsd [esp + nb231_Vvdwtot], xmm5

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
	mulsd  xmm7, [esp + nb231_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	
	movsd xmm4, [esp + nb231_c12]
	mulsd  xmm7, xmm4 
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [esp + nb231_Vvdwtot]
	movsd xmm3, [esp + nb231_fstmp]
	mulsd  xmm7, [esp + nb231_tsc]
	subsd  xmm3, xmm7
	movsd [esp + nb231_Vvdwtot], xmm5

	mulsd  xmm3, xmm0
		
	movsd xmm0, [esp + nb231_dxO]
	movsd xmm1, [esp + nb231_dyO]
	movsd xmm2, [esp + nb231_dzO]

	mov    edi, [ebp + nb231_faction]
	mulsd  xmm0, xmm3
	mulsd  xmm1, xmm3
	mulsd  xmm2, xmm3

	;# update O forces 
	movapd xmm3, [esp + nb231_fixO]
	movapd xmm4, [esp + nb231_fiyO]
	movapd xmm7, [esp + nb231_fizO]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [esp + nb231_fixO], xmm3
	movlpd [esp + nb231_fiyO], xmm4
	movlpd [esp + nb231_fizO], xmm7
	;# update j forces with water O 
	movlpd [esp + nb231_fjx], xmm0
	movlpd [esp + nb231_fjy], xmm1
	movlpd [esp + nb231_fjz], xmm2

	;# H1 interactions 
	movsd   xmm6, [esp + nb231_rinvH1]	
	movsd   xmm4, xmm6
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movsd   xmm7, xmm6
	movsd   xmm0, [esp + nb231_krsqH1]
	addsd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	mulsd   xmm0, [esp + nb231_two]
	subsd   xmm6, [esp + nb231_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm6, [esp + nb231_qqH] ;# vcoul 
	mulsd   xmm7, [esp + nb231_qqH]
	mulsd  xmm4, xmm7		;# total fsH1 in xmm4 
	
	addsd  xmm6, [esp + nb231_vctot]

	movapd xmm0, [esp + nb231_dxH1]
	movapd xmm1, [esp + nb231_dyH1]
	movapd xmm2, [esp + nb231_dzH1]
	movlpd [esp + nb231_vctot], xmm6
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [esp + nb231_fixH1]
	movapd xmm4, [esp + nb231_fiyH1]
	movapd xmm7, [esp + nb231_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [esp + nb231_fixH1], xmm3
	movlpd [esp + nb231_fiyH1], xmm4
	movlpd [esp + nb231_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [esp + nb231_fjx]
	addsd  xmm1, [esp + nb231_fjy]
	addsd  xmm2, [esp + nb231_fjz]
	movlpd [esp + nb231_fjx], xmm0
	movlpd [esp + nb231_fjy], xmm1
	movlpd [esp + nb231_fjz], xmm2

	;# H2 interactions 
	movapd  xmm5, [esp + nb231_rinvH2]	
	movapd  xmm4, xmm5
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm5
	movsd   xmm0, [esp + nb231_krsqH2]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulsd   xmm0, [esp + nb231_two]
	subsd   xmm5, [esp + nb231_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm5, [esp + nb231_qqH] ;# vcoul 
	mulsd   xmm7, [esp + nb231_qqH]
	mulsd  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addsd  xmm5, [esp + nb231_vctot]

	movapd xmm0, [esp + nb231_dxH2]
	movapd xmm1, [esp + nb231_dyH2]
	movapd xmm2, [esp + nb231_dzH2]
	movlpd [esp + nb231_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [esp + nb231_fixH2]
	movapd xmm4, [esp + nb231_fiyH2]
	movapd xmm7, [esp + nb231_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [esp + nb231_fixH2], xmm3
	movlpd [esp + nb231_fiyH2], xmm4
	movlpd [esp + nb231_fizH2], xmm7

	mov edi, [ebp + nb231_faction]
	;# update j forces 
	addsd  xmm0, [esp + nb231_fjx]
	addsd  xmm1, [esp + nb231_fjy]
	addsd  xmm2, [esp + nb231_fjz]
	movlpd xmm3, [edi + eax*8]
	movlpd xmm4, [edi + eax*8 + 8]
	movlpd xmm5, [edi + eax*8 + 16]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	movlpd [edi + eax*8], xmm3
	movlpd [edi + eax*8 + 8], xmm4
	movlpd [edi + eax*8 + 16], xmm5	

.nb231_updateouterdata:
	mov   ecx, [esp + nb231_ii3]
	mov   edi, [ebp + nb231_faction]
	mov   esi, [ebp + nb231_fshift]
	mov   edx, [esp + nb231_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb231_fixO]
	movapd xmm1, [esp + nb231_fiyO]
	movapd xmm2, [esp + nb231_fizO]

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
	movapd xmm0, [esp + nb231_fixH1]
	movapd xmm1, [esp + nb231_fiyH1]
	movapd xmm2, [esp + nb231_fizH1]

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
	movapd xmm0, [esp + nb231_fixH2]
	movapd xmm1, [esp + nb231_fiyH2]
	movapd xmm2, [esp + nb231_fizH2]

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
	mov esi, [esp + nb231_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb231_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb231_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb231_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb231_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb231_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
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





.globl nb_kernel231nf_ia32_sse2
.globl _nb_kernel231nf_ia32_sse2
nb_kernel231nf_ia32_sse2:	
_nb_kernel231nf_ia32_sse2:	
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
	;# bottom of stack is cache-aligned for sse2 use 
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
.equiv          nb231nf_two,              320
.equiv          nb231nf_krf,              336
.equiv          nb231nf_crf,              352
.equiv          nb231nf_krsqO,            368
.equiv          nb231nf_krsqH1,           384
.equiv          nb231nf_krsqH2,           400
.equiv          nb231nf_rsqO,             416
.equiv          nb231nf_rinvH1,           432
.equiv          nb231nf_rinvH2,           448
.equiv          nb231nf_is3,              464
.equiv          nb231nf_ii3,              468
.equiv          nb231nf_ntia,             472
.equiv          nb231nf_innerjjnr,        476
.equiv          nb231nf_innerk,           480
.equiv          nb231nf_n,                484
.equiv          nb231nf_nn1,              488
.equiv          nb231nf_nri,              492
.equiv          nb231nf_nouter,           496
.equiv          nb231nf_ninner,           500
.equiv          nb231nf_salign,           504
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
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb231nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb231nf_half], eax
	mov [esp + nb231nf_half+4], ebx
	movsd xmm1, [esp + nb231nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb231nf_half], xmm1
	movapd [esp + nb231nf_two], xmm2
	movapd [esp + nb231nf_three], xmm3

	mov esi, [ebp + nb231nf_argkrf]
	mov edi, [ebp + nb231nf_argcrf]
	movsd xmm5, [esi]
	movsd xmm6, [edi]
	shufpd xmm5, xmm5, 0
	shufpd xmm6, xmm6, 0
	movapd [esp + nb231nf_krf], xmm5
	movapd [esp + nb231nf_crf], xmm6
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb231nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb231nf_charge]
	movsd xmm3, [edx + ebx*8]	
	movsd xmm4, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb231nf_p_facel]
	movsd xmm5, [esi]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb231nf_iqO], xmm3
	movapd [esp + nb231nf_iqH], xmm4
	
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
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb231nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb231nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb231nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb231nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb231nf_ii3], ebx

	addsd xmm3, [eax + ebx*8]
	addsd xmm4, [eax + ebx*8 + 8]
	addsd xmm5, [eax + ebx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb231nf_ixO], xmm3
	movapd [esp + nb231nf_iyO], xmm4
	movapd [esp + nb231nf_izO], xmm5

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
	movapd [esp + nb231nf_ixH1], xmm0
	movapd [esp + nb231nf_iyH1], xmm1
	movapd [esp + nb231nf_izH1], xmm2
	movapd [esp + nb231nf_ixH2], xmm3
	movapd [esp + nb231nf_iyH2], xmm4
	movapd [esp + nb231nf_izH2], xmm5
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb231nf_vctot], xmm4
	movapd [esp + nb231nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb231nf_jindex]
	mov   ecx, [eax+esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb231nf_pos]
	mov   edi, [ebp + nb231nf_faction]	
	mov   eax, [ebp + nb231nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb231nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb231nf_ninner]
	mov   [esp + nb231nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb231nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb231nf_unroll_loop
	jmp   .nb231nf_checksingle
.nb231nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb231nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]

	add dword ptr [esp + nb231nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb231nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb231nf_iqO]
	mulpd  xmm4, [esp + nb231nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx

	movapd  [esp + nb231nf_qqO], xmm3
	movapd  [esp + nb231nf_qqH], xmm4
	
	mov esi, [ebp + nb231nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb231nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [esp + nb231nf_ntia]
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
	movapd [esp + nb231nf_c6], xmm4
	movapd [esp + nb231nf_c12], xmm6
	
	mov esi, [ebp + nb231nf_pos]       ;# base of pos[] 

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
	movapd xmm4, [esp + nb231nf_ixO]
	movapd xmm5, [esp + nb231nf_iyO]
	movapd xmm6, [esp + nb231nf_izO]

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
	movapd [esp + nb231nf_rsqO], xmm7
	
	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb231nf_ixH1]
	movapd xmm5, [esp + nb231nf_iyH1]
	movapd xmm6, [esp + nb231nf_izH1]

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
	movapd xmm3, [esp + nb231nf_ixH2]
	movapd xmm4, [esp + nb231nf_iyH2]
	movapd xmm5, [esp + nb231nf_izH2]

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

	movapd xmm0, xmm5
	movapd xmm1, xmm6
	movapd xmm2, xmm7

	mulpd  xmm0, [esp + nb231nf_krf]	
	mulpd  xmm1, [esp + nb231nf_krf]	
	mulpd  xmm2, [esp + nb231nf_krf]	

	movapd [esp + nb231nf_krsqH2], xmm0
	movapd [esp + nb231nf_krsqH1], xmm1
	movapd [esp + nb231nf_krsqO], xmm2
		
	;# start with rsqO - put seed in xmm2 
	cvtpd2ps xmm2, xmm7	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb231nf_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb231nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb231nf_three]
	subpd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb231nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb231nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb231nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb231nf_three]
	subpd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb231nf_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	movapd  [esp + nb231nf_rinvH1], xmm6
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb231nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb231nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb231nf_three]
	subpd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb231nf_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 
	movapd  [esp + nb231nf_rinvH2], xmm5

	;# do O interactions 
	movapd xmm0, xmm7
	movapd xmm1, [esp + nb231nf_krsqO]
	addpd  xmm0, xmm1
	subpd  xmm0, [esp + nb231nf_crf] ;# xmm0=rinv+ krsq-crf 
	mulpd  xmm0, [esp + nb231nf_qqO]

	mulpd  xmm6, xmm7

	addpd  xmm0, [esp + nb231nf_vctot]
	movapd [esp + nb231nf_vctot], xmm0

	movapd xmm0, xmm7
	movapd xmm4, [esp + nb231nf_rsqO]
	;# LJ table interaction. xmm0=rinv, xmm4=rsq
	
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb231nf_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	mov  esi, [ebp + nb231nf_VFtab]
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

	movapd xmm4, [esp + nb231nf_c6]
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;#  Update Vvdwtot directly 
	addpd  xmm5, [esp + nb231nf_Vvdwtot]
	movapd [esp + nb231nf_Vvdwtot], xmm5

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
	
	movapd xmm4, [esp + nb231nf_c12]
	mulpd  xmm5, xmm4  
	
	addpd  xmm5, [esp + nb231nf_Vvdwtot]
	movapd [esp + nb231nf_Vvdwtot], xmm5

	;# H1 interactions 
	movapd  xmm6, [esp + nb231nf_rinvH1]	
	movapd  xmm4, xmm6
	mulpd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm6
	movapd  xmm0, [esp + nb231nf_krsqH1]
	addpd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	subpd   xmm6, [esp + nb231nf_crf]
	mulpd   xmm6, [esp + nb231nf_qqH] ;# vcoul 
	
	addpd  xmm6, [esp + nb231nf_vctot]
	movapd [esp + nb231nf_vctot], xmm6

	;# H2 interactions 
	movapd  xmm5, [esp + nb231nf_rinvH2]	
	movapd  xmm4, xmm5
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm5
	movapd  xmm0, [esp + nb231nf_krsqH2]
	addpd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	subpd   xmm5, [esp + nb231nf_crf]
	mulpd   xmm5, [esp + nb231nf_qqH] ;# vcoul 
	
	addpd  xmm5, [esp + nb231nf_vctot]
	movapd [esp + nb231nf_vctot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb231nf_innerk],  2
	jl    .nb231nf_checksingle
	jmp   .nb231nf_unroll_loop
.nb231nf_checksingle:	
	mov   edx, [esp + nb231nf_innerk]
	and   edx, 1
	jnz   .nb231nf_dosingle
	jmp   .nb231nf_updateouterdata
.nb231nf_dosingle:
	mov   edx, [esp + nb231nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb231nf_innerjjnr],  4	

	mov esi, [ebp + nb231nf_charge]    ;# base of charge[] 
	xorpd xmm3, xmm3
	movlpd xmm3, [esi + eax*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb231nf_iqO]
	mulpd  xmm4, [esp + nb231nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 

	movapd  [esp + nb231nf_qqO], xmm3
	movapd  [esp + nb231nf_qqH], xmm4
	
	mov esi, [ebp + nb231nf_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb231nf_vdwparam]
	shl eax, 1	
	mov edi, [esp + nb231nf_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movapd [esp + nb231nf_c6], xmm4
	movapd [esp + nb231nf_c12], xmm6
	
	mov esi, [ebp + nb231nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [esp + nb231nf_ixO]
	movapd xmm5, [esp + nb231nf_iyO]
	movapd xmm6, [esp + nb231nf_izO]

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
	movapd [esp + nb231nf_rsqO], xmm7
	
	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb231nf_ixH1]
	movapd xmm5, [esp + nb231nf_iyH1]
	movapd xmm6, [esp + nb231nf_izH1]

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
	movapd xmm3, [esp + nb231nf_ixH2]
	movapd xmm4, [esp + nb231nf_iyH2]
	movapd xmm5, [esp + nb231nf_izH2]

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

	movapd xmm0, xmm5
	movapd xmm1, xmm6
	movapd xmm2, xmm7

	mulsd  xmm0, [esp + nb231nf_krf]	
	mulsd  xmm1, [esp + nb231nf_krf]	
	mulsd  xmm2, [esp + nb231nf_krf]	

	movapd [esp + nb231nf_krsqH2], xmm0
	movapd [esp + nb231nf_krsqH1], xmm1
	movapd [esp + nb231nf_krsqO], xmm2

	;# start with rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb231nf_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb231nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb231nf_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb231nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb231nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb231nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb231nf_three]
	subsd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb231nf_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	movsd  [esp + nb231nf_rinvH1], xmm6
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb231nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb231nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb231nf_three]
	subsd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb231nf_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 
	movsd [esp + nb231nf_rinvH2], xmm5
	
	;# do O interactions 
	movsd xmm0, xmm7
	movsd xmm1, [esp + nb231nf_krsqO]
	addsd  xmm0, xmm1
	subsd  xmm0, [esp + nb231nf_crf] ;# xmm0=rinv+ krsq-crf 
	mulsd  xmm0, [esp + nb231nf_qqO]

	addsd  xmm0, [esp + nb231nf_vctot]
	movsd [esp + nb231nf_vctot], xmm0

	movsd xmm0, xmm7
	movsd xmm4, [esp + nb231nf_rsqO]
	;# LJ table interaction. xmm0=rinv, xmm4=rsq
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb231nf_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3
	
	mov  esi, [ebp + nb231nf_VFtab]

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

	movsd xmm4, [esp + nb231nf_c6]
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# Update Vvdwtot directly 
	addsd  xmm5, [esp + nb231nf_Vvdwtot]
	movsd [esp + nb231nf_Vvdwtot], xmm5

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
	
	movsd xmm4, [esp + nb231nf_c12]
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [esp + nb231nf_Vvdwtot]
	movsd [esp + nb231nf_Vvdwtot], xmm5

	;# H1 interactions 
	movsd   xmm6, [esp + nb231nf_rinvH1]	
	movsd   xmm4, xmm6
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movsd   xmm0, [esp + nb231nf_krsqH1]
	addsd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	subsd   xmm6, [esp + nb231nf_crf]
	mulsd   xmm6, [esp + nb231nf_qqH] ;# vcoul 
	
	addsd  xmm6, [esp + nb231nf_vctot]
	movlpd [esp + nb231nf_vctot], xmm6

	;# H2 interactions 
	movapd  xmm5, [esp + nb231nf_rinvH2]	
	movapd  xmm4, xmm5
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm5
	movsd   xmm0, [esp + nb231nf_krsqH2]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	subsd   xmm5, [esp + nb231nf_crf]
	mulsd   xmm5, [esp + nb231nf_qqH] ;# vcoul 
	
	addsd  xmm5, [esp + nb231nf_vctot]
	movlpd [esp + nb231nf_vctot], xmm5

.nb231nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb231nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb231nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb231nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb231nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb231nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb231nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
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
	add esp, 508
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret


