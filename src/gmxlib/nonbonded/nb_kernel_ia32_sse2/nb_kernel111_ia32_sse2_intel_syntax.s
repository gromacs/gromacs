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


.globl nb_kernel111_ia32_sse2
.globl _nb_kernel111_ia32_sse2
nb_kernel111_ia32_sse2:	
_nb_kernel111_ia32_sse2:	
.equiv          nb111_p_nri,            8
.equiv          nb111_iinr,             12
.equiv          nb111_jindex,           16
.equiv          nb111_jjnr,             20
.equiv          nb111_shift,            24
.equiv          nb111_shiftvec,         28
.equiv          nb111_fshift,           32
.equiv          nb111_gid,              36
.equiv          nb111_pos,              40
.equiv          nb111_faction,          44
.equiv          nb111_charge,           48
.equiv          nb111_p_facel,          52
.equiv          nb111_argkrf,           56
.equiv          nb111_argcrf,           60
.equiv          nb111_Vc,               64
.equiv          nb111_type,             68
.equiv          nb111_p_ntype,          72
.equiv          nb111_vdwparam,         76
.equiv          nb111_Vvdw,             80
.equiv          nb111_p_tabscale,       84
.equiv          nb111_VFtab,            88
.equiv          nb111_invsqrta,         92
.equiv          nb111_dvda,             96
.equiv          nb111_p_gbtabscale,     100
.equiv          nb111_GBtab,            104
.equiv          nb111_p_nthreads,       108
.equiv          nb111_count,            112
.equiv          nb111_mtx,              116
.equiv          nb111_outeriter,        120
.equiv          nb111_inneriter,        124
.equiv          nb111_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb111_ixO,              0
.equiv          nb111_iyO,              16
.equiv          nb111_izO,              32
.equiv          nb111_ixH1,             48
.equiv          nb111_iyH1,             64
.equiv          nb111_izH1,             80
.equiv          nb111_ixH2,             96
.equiv          nb111_iyH2,             112
.equiv          nb111_izH2,             128
.equiv          nb111_iqO,              144
.equiv          nb111_iqH,              160
.equiv          nb111_dxO,              176
.equiv          nb111_dyO,              192
.equiv          nb111_dzO,              208
.equiv          nb111_dxH1,             224
.equiv          nb111_dyH1,             240
.equiv          nb111_dzH1,             256
.equiv          nb111_dxH2,             272
.equiv          nb111_dyH2,             288
.equiv          nb111_dzH2,             304
.equiv          nb111_qqO,              320
.equiv          nb111_qqH,              336
.equiv          nb111_c6,               352
.equiv          nb111_c12,              368
.equiv          nb111_six,              384
.equiv          nb111_twelve,           400
.equiv          nb111_vctot,            416
.equiv          nb111_Vvdwtot,          432
.equiv          nb111_fixO,             448
.equiv          nb111_fiyO,             464
.equiv          nb111_fizO,             480
.equiv          nb111_fixH1,            496
.equiv          nb111_fiyH1,            512
.equiv          nb111_fizH1,            528
.equiv          nb111_fixH2,            544
.equiv          nb111_fiyH2,            560
.equiv          nb111_fizH2,            576
.equiv          nb111_fjx,              592
.equiv          nb111_fjy,              608
.equiv          nb111_fjz,              624
.equiv          nb111_half,             640
.equiv          nb111_three,            656
.equiv          nb111_is3,              672
.equiv          nb111_ii3,              676
.equiv          nb111_ntia,             680
.equiv          nb111_innerjjnr,        684
.equiv          nb111_innerk,           688
.equiv          nb111_n,                692
.equiv          nb111_nn1,              696
.equiv          nb111_nri,              700
.equiv          nb111_nouter,           704
.equiv          nb111_ninner,           708
.equiv          nb111_salign,           712
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 716		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb111_salign], eax
	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb111_p_nri]
	mov ecx, [ecx]
	mov [esp + nb111_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb111_nouter], eax
	mov [esp + nb111_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb111_half], eax
	mov [esp + nb111_half+4], ebx
	movsd xmm1, [esp + nb111_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd xmm4, xmm3
	addpd  xmm4, xmm4       ;# 6.0
	movapd xmm5, xmm4
	addpd  xmm5, xmm5       ;# 12.0
	movapd [esp + nb111_half], xmm1
	movapd [esp + nb111_three], xmm3
	movapd [esp + nb111_six], xmm4
	movapd [esp + nb111_twelve], xmm5

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb111_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb111_charge]
	movsd xmm3, [edx + ebx*8]	
	movsd xmm4, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb111_p_facel]
	movsd xmm5, [esi]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb111_iqO], xmm3
	movapd [esp + nb111_iqH], xmm4
	
	mov   edx, [ebp + nb111_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov edi, [ebp + nb111_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	mov   [esp + nb111_ntia], ecx		
.nb111_threadloop:
        mov   esi, [ebp + nb111_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb111_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb111_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb111_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb111_n], eax
        mov [esp + nb111_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb111_outerstart
        jmp .nb111_end

.nb111_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb111_nouter]
	mov [esp + nb111_nouter], ebx

.nb111_outer:
	mov   eax, [ebp + nb111_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb111_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb111_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb111_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb111_pos]    ;# eax = base of pos[]  
	mov   [esp + nb111_ii3], ebx

	addsd xmm3, [eax + ebx*8]
	addsd xmm4, [eax + ebx*8 + 8]
	addsd xmm5, [eax + ebx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb111_ixO], xmm3
	movapd [esp + nb111_iyO], xmm4
	movapd [esp + nb111_izO], xmm5

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
	movapd [esp + nb111_ixH1], xmm0
	movapd [esp + nb111_iyH1], xmm1
	movapd [esp + nb111_izH1], xmm2
	movapd [esp + nb111_ixH2], xmm3
	movapd [esp + nb111_iyH2], xmm4
	movapd [esp + nb111_izH2], xmm5
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb111_vctot], xmm4
	movapd [esp + nb111_Vvdwtot], xmm4
	movapd [esp + nb111_fixO], xmm4
	movapd [esp + nb111_fiyO], xmm4
	movapd [esp + nb111_fizO], xmm4
	movapd [esp + nb111_fixH1], xmm4
	movapd [esp + nb111_fiyH1], xmm4
	movapd [esp + nb111_fizH1], xmm4
	movapd [esp + nb111_fixH2], xmm4
	movapd [esp + nb111_fiyH2], xmm4
	movapd [esp + nb111_fizH2], xmm4
	
	mov   eax, [ebp + nb111_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb111_pos]
	mov   edi, [ebp + nb111_faction]	
	mov   eax, [ebp + nb111_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb111_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb111_ninner]
	mov   [esp + nb111_ninner], ecx
	add   edx, 0
	mov   [esp + nb111_innerk], edx    ;# number of innerloop atoms 
	jge   .nb111_unroll_loop
	jmp   .nb111_checksingle
.nb111_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb111_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]

	add dword ptr [esp + nb111_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb111_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb111_iqO]
	mulpd  xmm4, [esp + nb111_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx

	movapd  [esp + nb111_qqO], xmm3
	movapd  [esp + nb111_qqH], xmm4
	
	mov esi, [ebp + nb111_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb111_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [esp + nb111_ntia]
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
	movapd [esp + nb111_c6], xmm4
	movapd [esp + nb111_c12], xmm6
	
	mov esi, [ebp + nb111_pos]       ;# base of pos[] 

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
	movapd xmm4, [esp + nb111_ixO]
	movapd xmm5, [esp + nb111_iyO]
	movapd xmm6, [esp + nb111_izO]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb111_dxO], xmm4
	movapd [esp + nb111_dyO], xmm5
	movapd [esp + nb111_dzO], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb111_ixH1]
	movapd xmm5, [esp + nb111_iyH1]
	movapd xmm6, [esp + nb111_izH1]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb111_dxH1], xmm4
	movapd [esp + nb111_dyH1], xmm5
	movapd [esp + nb111_dzH1], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm6, xmm5
	addpd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb111_ixH2]
	movapd xmm4, [esp + nb111_iyH2]
	movapd xmm5, [esp + nb111_izH2]

	;# calc dr 
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2

	;# store dr 
	movapd [esp + nb111_dxH2], xmm3
	movapd [esp + nb111_dyH2], xmm4
	movapd [esp + nb111_dzH2], xmm5
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
	movapd  xmm4, [esp + nb111_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb111_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb111_three]
	subpd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb111_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb111_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb111_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb111_three]
	subpd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb111_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb111_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb111_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb111_three]
	subpd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb111_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 
	movapd  xmm4, xmm7	
	mulpd   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq 
	movapd xmm1, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulpd  xmm7, [esp + nb111_qqO]	;# xmm7=vcoul 
	
	mulpd  xmm1, [esp + nb111_c6]
	mulpd  xmm2, [esp + nb111_c12]
	movapd xmm3, xmm2
	subpd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addpd  xmm3, [esp + nb111_Vvdwtot]
	mulpd  xmm1, [esp + nb111_six]
	mulpd  xmm2, [esp + nb111_twelve]
	subpd  xmm2, xmm1
	addpd  xmm2, xmm7	
	mulpd  xmm4, xmm2	;# total fsO in xmm4 

	addpd  xmm7, [esp + nb111_vctot]
	
	movapd [esp + nb111_Vvdwtot], xmm3
	movapd [esp + nb111_vctot], xmm7

	movapd xmm0, [esp + nb111_dxO]
	movapd xmm1, [esp + nb111_dyO]
	movapd xmm2, [esp + nb111_dzO]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update O forces 
	movapd xmm3, [esp + nb111_fixO]
	movapd xmm4, [esp + nb111_fiyO]
	movapd xmm7, [esp + nb111_fizO]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb111_fixO], xmm3
	movapd [esp + nb111_fiyO], xmm4
	movapd [esp + nb111_fizO], xmm7
	;# update j forces with water O 
	movapd [esp + nb111_fjx], xmm0
	movapd [esp + nb111_fjy], xmm1
	movapd [esp + nb111_fjz], xmm2

	;# H1 interactions 
	movapd  xmm4, xmm6	
	mulpd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	mulpd  xmm6, [esp + nb111_qqH]	;# xmm6=vcoul 
	mulpd  xmm4, xmm6		;# total fsH1 in xmm4 
	
	addpd  xmm6, [esp + nb111_vctot]

	movapd xmm0, [esp + nb111_dxH1]
	movapd xmm1, [esp + nb111_dyH1]
	movapd xmm2, [esp + nb111_dzH1]
	movapd [esp + nb111_vctot], xmm6
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [esp + nb111_fixH1]
	movapd xmm4, [esp + nb111_fiyH1]
	movapd xmm7, [esp + nb111_fizH1]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb111_fixH1], xmm3
	movapd [esp + nb111_fiyH1], xmm4
	movapd [esp + nb111_fizH1], xmm7
	;# update j forces with water H1 
	addpd  xmm0, [esp + nb111_fjx]
	addpd  xmm1, [esp + nb111_fjy]
	addpd  xmm2, [esp + nb111_fjz]
	movapd [esp + nb111_fjx], xmm0
	movapd [esp + nb111_fjy], xmm1
	movapd [esp + nb111_fjz], xmm2

	;# H2 interactions 
	movapd  xmm4, xmm5	
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	mulpd  xmm5, [esp + nb111_qqH]	;# xmm5=vcoul 
	mulpd  xmm4, xmm5		;# total fsH1 in xmm4 
	
	addpd  xmm5, [esp + nb111_vctot]

	movapd xmm0, [esp + nb111_dxH2]
	movapd xmm1, [esp + nb111_dyH2]
	movapd xmm2, [esp + nb111_dzH2]
	movapd [esp + nb111_vctot], xmm5
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [esp + nb111_fixH2]
	movapd xmm4, [esp + nb111_fiyH2]
	movapd xmm7, [esp + nb111_fizH2]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm7, xmm2
	movapd [esp + nb111_fixH2], xmm3
	movapd [esp + nb111_fiyH2], xmm4
	movapd [esp + nb111_fizH2], xmm7

	mov edi, [ebp + nb111_faction]
	;# update j forces 
	addpd  xmm0, [esp + nb111_fjx]
	addpd  xmm1, [esp + nb111_fjy]
	addpd  xmm2, [esp + nb111_fjz]
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
	sub dword ptr [esp + nb111_innerk],  2
	jl   .nb111_checksingle
	jmp  .nb111_unroll_loop
.nb111_checksingle:	
	mov   edx, [esp + nb111_innerk]
	and   edx, 1
	jnz  .nb111_dosingle
	jmp  .nb111_updateouterdata
.nb111_dosingle:
	mov   edx, [esp + nb111_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb111_innerjjnr],  4	

	mov esi, [ebp + nb111_charge]    ;# base of charge[] 

	xorpd xmm3, xmm3
	movlpd xmm3, [esi + eax*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb111_iqO]
	mulpd  xmm4, [esp + nb111_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 

	movapd  [esp + nb111_qqO], xmm3
	movapd  [esp + nb111_qqH], xmm4
	
	mov esi, [ebp + nb111_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb111_vdwparam]
	shl eax, 1	
	mov edi, [esp + nb111_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [esp + nb111_c6], xmm4
	movapd [esp + nb111_c12], xmm6
	
	mov esi, [ebp + nb111_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [esp + nb111_ixO]
	movapd xmm5, [esp + nb111_iyO]
	movapd xmm6, [esp + nb111_izO]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb111_dxO], xmm4
	movapd [esp + nb111_dyO], xmm5
	movapd [esp + nb111_dzO], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [esp + nb111_ixH1]
	movapd xmm5, [esp + nb111_iyH1]
	movapd xmm6, [esp + nb111_izH1]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb111_dxH1], xmm4
	movapd [esp + nb111_dyH1], xmm5
	movapd [esp + nb111_dzH1], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [esp + nb111_ixH2]
	movapd xmm4, [esp + nb111_iyH2]
	movapd xmm5, [esp + nb111_izH2]

	;# calc dr 
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2

	;# store dr 
	movapd [esp + nb111_dxH2], xmm3
	movapd [esp + nb111_dyH2], xmm4
	movapd [esp + nb111_dzH2], xmm5
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
	movapd  xmm4, [esp + nb111_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb111_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb111_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb111_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb111_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb111_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb111_three]
	subsd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb111_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb111_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb111_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb111_three]
	subsd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb111_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 
	movapd  xmm4, xmm7	
	mulsd   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq 
	movapd xmm1, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm7, [esp + nb111_qqO]	;# xmm7=vcoul 
	
	mulsd  xmm1, [esp + nb111_c6]
	mulsd  xmm2, [esp + nb111_c12]
	movapd xmm3, xmm2
	subsd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addsd  xmm3, [esp + nb111_Vvdwtot]
	mulsd  xmm1, [esp + nb111_six]
	mulsd  xmm2, [esp + nb111_twelve]
	subsd  xmm2, xmm1
	addsd  xmm2, xmm7	
	mulsd  xmm4, xmm2	;# total fsO in xmm4 

	addsd  xmm7, [esp + nb111_vctot]
	
	movsd [esp + nb111_Vvdwtot], xmm3
	movsd [esp + nb111_vctot], xmm7

	movapd xmm0, [esp + nb111_dxO]
	movapd xmm1, [esp + nb111_dyO]
	movapd xmm2, [esp + nb111_dzO]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update O forces 
	movapd xmm3, [esp + nb111_fixO]
	movapd xmm4, [esp + nb111_fiyO]
	movapd xmm7, [esp + nb111_fizO]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [esp + nb111_fixO], xmm3
	movsd [esp + nb111_fiyO], xmm4
	movsd [esp + nb111_fizO], xmm7
	;# update j forces with water O 
	movsd [esp + nb111_fjx], xmm0
	movsd [esp + nb111_fjy], xmm1
	movsd [esp + nb111_fjz], xmm2

	;# H1 interactions 
	movapd  xmm4, xmm6	
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	mulsd  xmm6, [esp + nb111_qqH]	;# xmm6=vcoul 
	mulsd  xmm4, xmm6		;# total fsH1 in xmm4 
	
	addsd  xmm6, [esp + nb111_vctot]

	movapd xmm0, [esp + nb111_dxH1]
	movapd xmm1, [esp + nb111_dyH1]
	movapd xmm2, [esp + nb111_dzH1]
	movsd [esp + nb111_vctot], xmm6
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [esp + nb111_fixH1]
	movapd xmm4, [esp + nb111_fiyH1]
	movapd xmm7, [esp + nb111_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [esp + nb111_fixH1], xmm3
	movsd [esp + nb111_fiyH1], xmm4
	movsd [esp + nb111_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [esp + nb111_fjx]
	addsd  xmm1, [esp + nb111_fjy]
	addsd  xmm2, [esp + nb111_fjz]
	movsd [esp + nb111_fjx], xmm0
	movsd [esp + nb111_fjy], xmm1
	movsd [esp + nb111_fjz], xmm2

	;# H2 interactions 
	movapd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	mulsd  xmm5, [esp + nb111_qqH]	;# xmm5=vcoul 
	mulsd  xmm4, xmm5		;# total fsH1 in xmm4 
	
	addsd  xmm5, [esp + nb111_vctot]

	movapd xmm0, [esp + nb111_dxH2]
	movapd xmm1, [esp + nb111_dyH2]
	movapd xmm2, [esp + nb111_dzH2]
	movsd [esp + nb111_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [esp + nb111_fixH2]
	movapd xmm4, [esp + nb111_fiyH2]
	movapd xmm7, [esp + nb111_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [esp + nb111_fixH2], xmm3
	movsd [esp + nb111_fiyH2], xmm4
	movsd [esp + nb111_fizH2], xmm7

	mov edi, [ebp + nb111_faction]
	;# update j forces 
	addsd  xmm0, [esp + nb111_fjx]
	addsd  xmm1, [esp + nb111_fjy]
	addsd  xmm2, [esp + nb111_fjz]
	movlpd xmm3, [edi + eax*8]
	movlpd xmm4, [edi + eax*8 + 8]
	movlpd xmm5, [edi + eax*8 + 16]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	movlpd [edi + eax*8], xmm3
	movlpd [edi + eax*8 + 8], xmm4
	movlpd [edi + eax*8 + 16], xmm5

.nb111_updateouterdata:
	mov   ecx, [esp + nb111_ii3]
	mov   edi, [ebp + nb111_faction]
	mov   esi, [ebp + nb111_fshift]
	mov   edx, [esp + nb111_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb111_fixO]
	movapd xmm1, [esp + nb111_fiyO]
	movapd xmm2, [esp + nb111_fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

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
	unpcklpd xmm6,xmm1 

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb111_fixH1]
	movapd xmm1, [esp + nb111_fiyH1]
	movapd xmm2, [esp + nb111_fizH1]

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
	movapd xmm0, [esp + nb111_fixH2]
	movapd xmm1, [esp + nb111_fiyH2]
	movapd xmm2, [esp + nb111_fizH2]

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
	mov esi, [esp + nb111_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb111_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb111_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb111_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb111_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb111_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [esp + nb111_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb111_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb111_n], esi
        jmp .nb111_outer
.nb111_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb111_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb111_end
        ;# non-zero, do one more workunit
        jmp   .nb111_threadloop
.nb111_end:
	emms

	mov eax, [esp + nb111_nouter]
	mov ebx, [esp + nb111_ninner]
	mov ecx, [ebp + nb111_outeriter]
	mov edx, [ebp + nb111_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb111_salign]
	add esp, eax
	add esp, 716
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


.globl nb_kernel111nf_ia32_sse2
.globl _nb_kernel111nf_ia32_sse2
nb_kernel111nf_ia32_sse2:	
_nb_kernel111nf_ia32_sse2:	
.equiv          nb111nf_p_nri,          8
.equiv          nb111nf_iinr,           12
.equiv          nb111nf_jindex,         16
.equiv          nb111nf_jjnr,           20
.equiv          nb111nf_shift,          24
.equiv          nb111nf_shiftvec,       28
.equiv          nb111nf_fshift,         32
.equiv          nb111nf_gid,            36
.equiv          nb111nf_pos,            40
.equiv          nb111nf_faction,        44
.equiv          nb111nf_charge,         48
.equiv          nb111nf_p_facel,        52
.equiv          nb111nf_argkrf,         56
.equiv          nb111nf_argcrf,         60
.equiv          nb111nf_Vc,             64
.equiv          nb111nf_type,           68
.equiv          nb111nf_p_ntype,        72
.equiv          nb111nf_vdwparam,       76
.equiv          nb111nf_Vvdw,           80
.equiv          nb111nf_p_tabscale,     84
.equiv          nb111nf_VFtab,          88
.equiv          nb111nf_invsqrta,       92
.equiv          nb111nf_dvda,           96
.equiv          nb111nf_p_gbtabscale,   100
.equiv          nb111nf_GBtab,          104
.equiv          nb111nf_p_nthreads,     108
.equiv          nb111nf_count,          112
.equiv          nb111nf_mtx,            116
.equiv          nb111nf_outeriter,      120
.equiv          nb111nf_inneriter,      124
.equiv          nb111nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb111nf_ixO,            0
.equiv          nb111nf_iyO,            16
.equiv          nb111nf_izO,            32
.equiv          nb111nf_ixH1,           48
.equiv          nb111nf_iyH1,           64
.equiv          nb111nf_izH1,           80
.equiv          nb111nf_ixH2,           96
.equiv          nb111nf_iyH2,           112
.equiv          nb111nf_izH2,           128
.equiv          nb111nf_iqO,            144
.equiv          nb111nf_iqH,            160
.equiv          nb111nf_qqO,            176
.equiv          nb111nf_qqH,            192
.equiv          nb111nf_c6,             208
.equiv          nb111nf_c12,            224
.equiv          nb111nf_vctot,          240
.equiv          nb111nf_Vvdwtot,        256
.equiv          nb111nf_half,           272
.equiv          nb111nf_three,          288
.equiv          nb111nf_is3,            304
.equiv          nb111nf_ii3,            308
.equiv          nb111nf_ntia,           312
.equiv          nb111nf_innerjjnr,      316
.equiv          nb111nf_innerk,         320
.equiv          nb111nf_n,              324
.equiv          nb111nf_nn1,            328
.equiv          nb111nf_nri,            332
.equiv          nb111nf_nouter,         336
.equiv          nb111nf_ninner,         340
.equiv          nb111nf_salign,         344
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 348		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb111nf_salign], eax
	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb111nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb111nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb111nf_nouter], eax
	mov [esp + nb111nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb111nf_half], eax
	mov [esp + nb111nf_half+4], ebx
	movsd xmm1, [esp + nb111nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb111nf_half], xmm1
	movapd [esp + nb111nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb111nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb111nf_charge]
	movsd xmm3, [edx + ebx*8]	
	movsd xmm4, [edx + ebx*8 + 8]	
	mov esi, [ebp + nb111nf_p_facel]
	movsd xmm5, [esi]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb111nf_iqO], xmm3
	movapd [esp + nb111nf_iqH], xmm4
	
	mov   edx, [ebp + nb111nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov edi, [ebp + nb111nf_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	mov   [esp + nb111nf_ntia], ecx		
.nb111nf_threadloop:
        mov   esi, [ebp + nb111nf_count]        ;# pointer to sync counter
        mov   eax, [esi]
.nb111nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb111nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb111nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb111nf_n], eax
        mov [esp + nb111nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb111nf_outerstart
        jmp .nb111nf_end

.nb111nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb111nf_nouter]
	mov [esp + nb111nf_nouter], ebx

.nb111nf_outer:
	mov   eax, [ebp + nb111nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 

	mov   eax, [ebp + nb111nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb111nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb111nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb111nf_ii3], ebx

	addsd xmm3, [eax + ebx*8]
	addsd xmm4, [eax + ebx*8 + 8]
	addsd xmm5, [eax + ebx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [esp + nb111nf_ixO], xmm3
	movapd [esp + nb111nf_iyO], xmm4
	movapd [esp + nb111nf_izO], xmm5

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
	movapd [esp + nb111nf_ixH1], xmm0
	movapd [esp + nb111nf_iyH1], xmm1
	movapd [esp + nb111nf_izH1], xmm2
	movapd [esp + nb111nf_ixH2], xmm3
	movapd [esp + nb111nf_iyH2], xmm4
	movapd [esp + nb111nf_izH2], xmm5
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb111nf_vctot], xmm4
	movapd [esp + nb111nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb111nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb111nf_pos]
	mov   eax, [ebp + nb111nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb111nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb111nf_ninner]
	mov   [esp + nb111nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb111nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb111nf_unroll_loop
	jmp   .nb111nf_checksingle
.nb111nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb111nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]

	add dword ptr [esp + nb111nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb111nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb111nf_iqO]
	mulpd  xmm4, [esp + nb111nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx

	movapd  [esp + nb111nf_qqO], xmm3
	movapd  [esp + nb111nf_qqH], xmm4
	
	mov esi, [ebp + nb111nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb111nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [esp + nb111nf_ntia]
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
	movapd [esp + nb111nf_c6], xmm4
	movapd [esp + nb111nf_c12], xmm6
	
	mov esi, [ebp + nb111nf_pos]       ;# base of pos[] 

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
	movapd xmm4, [esp + nb111nf_ixO]
	movapd xmm5, [esp + nb111nf_iyO]
	movapd xmm6, [esp + nb111nf_izO]

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
	movapd xmm4, [esp + nb111nf_ixH1]
	movapd xmm5, [esp + nb111nf_iyH1]
	movapd xmm6, [esp + nb111nf_izH1]

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
	movapd xmm3, [esp + nb111nf_ixH2]
	movapd xmm4, [esp + nb111nf_iyH2]
	movapd xmm5, [esp + nb111nf_izH2]

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
	movapd  xmm4, [esp + nb111nf_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb111nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb111nf_three]
	subpd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb111nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb111nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb111nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb111nf_three]
	subpd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb111nf_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [esp + nb111nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [esp + nb111nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb111nf_three]
	subpd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [esp + nb111nf_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 
	movapd  xmm4, xmm7	
	mulpd   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq 
	movapd xmm1, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulpd  xmm7, [esp + nb111nf_qqO]	;# xmm7=vcoul 
	
	mulpd  xmm1, [esp + nb111nf_c6]
	mulpd  xmm2, [esp + nb111nf_c12]
	movapd xmm3, xmm2
	subpd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addpd  xmm3, [esp + nb111nf_Vvdwtot]
	addpd  xmm7, [esp + nb111nf_vctot]	
	movapd [esp + nb111nf_Vvdwtot], xmm3
	movapd [esp + nb111nf_vctot], xmm7

	;# H1 interactions 
	mulpd  xmm6, [esp + nb111nf_qqH]	;# xmm6=vcoul 
	addpd  xmm6, [esp + nb111nf_vctot]
	movapd [esp + nb111nf_vctot], xmm6

	;# H2 interactions 
	mulpd  xmm5, [esp + nb111nf_qqH]	;# xmm5=vcoul 
	addpd  xmm5, [esp + nb111nf_vctot]
	movapd [esp + nb111nf_vctot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb111nf_innerk],  2
	jl    .nb111nf_checksingle
	jmp   .nb111nf_unroll_loop
.nb111nf_checksingle:	
	mov   edx, [esp + nb111nf_innerk]
	and   edx, 1
	jnz   .nb111nf_dosingle
	jmp   .nb111nf_updateouterdata
.nb111nf_dosingle:
	mov   edx, [esp + nb111nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb111nf_innerjjnr],  4	

	mov esi, [ebp + nb111nf_charge]    ;# base of charge[] 

	xorpd xmm3, xmm3
	movlpd xmm3, [esi + eax*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [esp + nb111nf_iqO]
	mulpd  xmm4, [esp + nb111nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 

	movapd  [esp + nb111nf_qqO], xmm3
	movapd  [esp + nb111nf_qqH], xmm4
	
	mov esi, [ebp + nb111nf_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb111nf_vdwparam]
	shl eax, 1	
	mov edi, [esp + nb111nf_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [esp + nb111nf_c6], xmm4
	movapd [esp + nb111nf_c12], xmm6
	
	mov esi, [ebp + nb111nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [esp + nb111nf_ixO]
	movapd xmm5, [esp + nb111nf_iyO]
	movapd xmm6, [esp + nb111nf_izO]

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
	movapd xmm4, [esp + nb111nf_ixH1]
	movapd xmm5, [esp + nb111nf_iyH1]
	movapd xmm6, [esp + nb111nf_izH1]

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
	movapd xmm3, [esp + nb111nf_ixH2]
	movapd xmm4, [esp + nb111nf_iyH2]
	movapd xmm5, [esp + nb111nf_izH2]

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
	movapd  xmm4, [esp + nb111nf_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb111nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb111nf_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb111nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb111nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb111nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb111nf_three]
	subsd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb111nf_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [esp + nb111nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [esp + nb111nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [esp + nb111nf_three]
	subsd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [esp + nb111nf_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 
	movapd  xmm4, xmm7	
	mulsd   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq 
	movapd xmm1, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm7, [esp + nb111nf_qqO]	;# xmm7=vcoul 
	
	mulsd  xmm1, [esp + nb111nf_c6]
	mulsd  xmm2, [esp + nb111nf_c12]
	movapd xmm3, xmm2
	subsd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addsd  xmm3, [esp + nb111nf_Vvdwtot]
	addsd  xmm7, [esp + nb111nf_vctot]
	movsd [esp + nb111nf_Vvdwtot], xmm3
	movsd [esp + nb111nf_vctot], xmm7

	;# H1 interactions 
	mulsd  xmm6, [esp + nb111nf_qqH]	;# xmm6=vcoul 
	addsd  xmm6, [esp + nb111nf_vctot]
	movsd [esp + nb111nf_vctot], xmm6

	;# H2 interactions 
	mulsd  xmm5, [esp + nb111nf_qqH]	;# xmm5=vcoul 
	addsd  xmm5, [esp + nb111nf_vctot]
	movsd [esp + nb111nf_vctot], xmm5
	
.nb111nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb111nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb111nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb111nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb111nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb111nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb111nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb111nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb111nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb111nf_n], esi
        jmp .nb111nf_outer
.nb111nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb111nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb111nf_end
        ;# non-zero, do one more workunit
        jmp   .nb111nf_threadloop
.nb111nf_end:
	emms

	mov eax, [esp + nb111nf_nouter]
	mov ebx, [esp + nb111nf_ninner]
	mov ecx, [ebp + nb111nf_outeriter]
	mov edx, [ebp + nb111nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb111nf_salign]
	add esp, eax
	add esp, 348
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret

