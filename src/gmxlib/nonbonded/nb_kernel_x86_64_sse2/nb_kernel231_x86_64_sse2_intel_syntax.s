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


.globl nb_kernel231_x86_64_sse2
.globl _nb_kernel231_x86_64_sse2
nb_kernel231_x86_64_sse2:	
_nb_kernel231_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb231_fshift,           16
.equiv          nb231_gid,              24
.equiv          nb231_pos,              32
.equiv          nb231_faction,          40
.equiv          nb231_charge,           48
.equiv          nb231_p_facel,          56
.equiv          nb231_argkrf,           64
.equiv          nb231_argcrf,           72
.equiv          nb231_Vc,               80
.equiv          nb231_type,             88
.equiv          nb231_p_ntype,          96
.equiv          nb231_vdwparam,         104
.equiv          nb231_Vvdw,             112
.equiv          nb231_p_tabscale,       120
.equiv          nb231_VFtab,            128
.equiv          nb231_invsqrta,         136
.equiv          nb231_dvda,             144
.equiv          nb231_p_gbtabscale,     152
.equiv          nb231_GBtab,            160
.equiv          nb231_p_nthreads,       168
.equiv          nb231_count,            176
.equiv          nb231_mtx,              184
.equiv          nb231_outeriter,        192
.equiv          nb231_inneriter,        200
.equiv          nb231_work,             208
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
.equiv          nb231_rsqO,             720
.equiv          nb231_rsqH1,            736
.equiv          nb231_rsqH2,            752
.equiv          nb231_rinvO,            768
.equiv          nb231_rinvH1,           784
.equiv          nb231_rinvH2,           800
.equiv          nb231_facel,            816
.equiv          nb231_iinr,             824
.equiv          nb231_jindex,           832
.equiv          nb231_jjnr,             840
.equiv          nb231_shift,            848
.equiv          nb231_shiftvec,         856
.equiv          nb231_innerjjnr,        864
.equiv          nb231_nri,              872
.equiv          nb231_is3,              876
.equiv          nb231_ii3,              880
.equiv          nb231_ntia,             884
.equiv          nb231_innerk,           888
.equiv          nb231_n,                892
.equiv          nb231_nn1,              896
.equiv          nb231_nouter,           900
.equiv          nb231_ninner,           904

	push rbp
	mov  rbp, rsp
    
    ;# Push integer registers on stack
	push rbx
    push rsi
    push rdi
    push r12
    push r13
    push r14
    push r15

    ;# Make room for registers xmm6-xmm15 (10 registers=160 bytes)
    sub rsp, 168
    
    ;# Save xmm registers to stack
    movaps [rsp      ], xmm6
    movaps [rsp + 16 ], xmm7
    movaps [rsp + 32 ], xmm8
    movaps [rsp + 48 ], xmm9
    movaps [rsp + 64 ], xmm10
    movaps [rsp + 80 ], xmm11
    movaps [rsp + 96 ], xmm12
    movaps [rsp + 112], xmm13
    movaps [rsp + 128], xmm14
    movaps [rsp + 144], xmm15
    
; .if 0    # block below only read by NASM - special calling convention on win64
%ifidn __OUTPUT_FORMAT__, win64
    ;# Adjust rbp to account for shadow space (32) & two extra args (2*8) on stack
    add rbp, 48
    ;# Adjust stack pointer for different alignment
    ;# Move around arguments to fit AMD64 convention below
    ;# AMD64 passes args in: rdi,rsi,rdx,rcx,r8,r9 + stack
    ;# win64 passes args in: rcx,rdx,r8,r9         + stack
    mov rdi, rcx
    mov rsi, rdx
    mov rdx, r8
    mov rcx, r9
    mov r8,  [rbp]
    mov r9,  [rbp + 8]
%endif
; .endif   # end NASM- and win64-specific block

	emms
    sub rsp, 912    ;  # local variable stack space (n*16+8)
	
	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb231_nouter], eax
	mov [rsp + nb231_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb231_nri], edi
	mov [rsp + nb231_iinr], rsi
	mov [rsp + nb231_jindex], rdx
	mov [rsp + nb231_jjnr], rcx
	mov [rsp + nb231_shift], r8
	mov [rsp + nb231_shiftvec], r9
	mov rsi, [rbp + nb231_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb231_facel], xmm0

	mov rax, [rbp + nb231_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb231_tsc], xmm3

	mov rsi, [rbp + nb231_argkrf]
	mov rdi, [rbp + nb231_argcrf]
	movsd xmm1, [rsi]
	movsd xmm2, [rdi]
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	movapd [rsp + nb231_krf], xmm1
	movapd [rsp + nb231_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb231_half], eax
	mov [rsp + nb231_half+4], ebx
	movsd xmm1, [rsp + nb231_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb231_half], xmm1
	movapd [rsp + nb231_two], xmm2
	movapd [rsp + nb231_three], xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb231_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb231_charge]
	movsd xmm3, [rdx + rbx*8]	
	movsd xmm4, [rdx + rbx*8 + 8]	

	movsd xmm5, [rsp + nb231_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb231_iqO], xmm3
	movapd [rsp + nb231_iqH], xmm4
	
	mov   rdx, [rbp + nb231_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb231_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb231_ntia], ecx		
.nb231_threadloop:
        mov   rsi, [rbp + nb231_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb231_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb231_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb231_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb231_n], eax
        mov [rsp + nb231_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb231_outerstart
        jmp .nb231_end

.nb231_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb231_nouter]
	mov [rsp + nb231_nouter], ebx

.nb231_outer:
	mov   rax, [rsp + nb231_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb231_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb231_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb231_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb231_pos]    ;# eax = base of pos[]  
	mov   [rsp + nb231_ii3], ebx

	addsd xmm3, [rax + rbx*8]
	addsd xmm4, [rax + rbx*8 + 8]
	addsd xmm5, [rax + rbx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb231_ixO], xmm3
	movapd [rsp + nb231_iyO], xmm4
	movapd [rsp + nb231_izO], xmm5

	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm0, [rax + rbx*8 + 24]
	addsd xmm1, [rax + rbx*8 + 32]
	addsd xmm2, [rax + rbx*8 + 40]		
	addsd xmm3, [rax + rbx*8 + 48]
	addsd xmm4, [rax + rbx*8 + 56]
	addsd xmm5, [rax + rbx*8 + 64]		

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb231_ixH1], xmm0
	movapd [rsp + nb231_iyH1], xmm1
	movapd [rsp + nb231_izH1], xmm2
	movapd [rsp + nb231_ixH2], xmm3
	movapd [rsp + nb231_iyH2], xmm4
	movapd [rsp + nb231_izH2], xmm5
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb231_vctot], xmm4
	movapd [rsp + nb231_Vvdwtot], xmm4
	movapd [rsp + nb231_fixO], xmm4
	movapd [rsp + nb231_fiyO], xmm4
	movapd [rsp + nb231_fizO], xmm4
	movapd [rsp + nb231_fixH1], xmm4
	movapd [rsp + nb231_fiyH1], xmm4
	movapd [rsp + nb231_fizH1], xmm4
	movapd [rsp + nb231_fixH2], xmm4
	movapd [rsp + nb231_fiyH2], xmm4
	movapd [rsp + nb231_fizH2], xmm4
	
	mov   rax, [rsp + nb231_jindex]
	mov   ecx, [rax+rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb231_pos]
	mov   rdi, [rbp + nb231_faction]	
	mov   rax, [rsp + nb231_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb231_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb231_ninner]
	mov   [rsp + nb231_ninner], ecx
	add   edx, 0
	mov   [rsp + nb231_innerk], edx    ;# number of innerloop atoms 
	jge   .nb231_unroll_loop
	jmp   .nb231_checksingle
.nb231_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb231_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]

	add qword ptr [rsp + nb231_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb231_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb231_iqO]
	mulpd  xmm4, [rsp + nb231_iqH]

	movapd  [rsp + nb231_qqO], xmm3
	movapd  [rsp + nb231_qqH], xmm4
	
	mov rsi, [rbp + nb231_type]
	mov r8d, [rsi + rax*4]
	mov r9d, [rsi + rbx*4]
	mov rsi, [rbp + nb231_vdwparam]
	shl r8d, 1	
	shl r9d, 1	
	mov edi, [rsp + nb231_ntia]
	add r8d, edi
	add r9d, edi

	movlpd xmm6, [rsi + r8*8]	;# c6a
	movlpd xmm7, [rsi + r9*8]	;# c6b
	movhpd xmm6, [rsi + r8*8 + 8]	;# c6a c12a 
	movhpd xmm7, [rsi + r9*8 + 8]	;# c6b c12b 
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7

	movapd [rsp + nb231_c6], xmm4
	movapd [rsp + nb231_c12], xmm6
	
	mov rsi, [rbp + nb231_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move j coordinates to local temp variables 
    movlpd xmm0, [rsi + rax*8] 
    movlpd xmm1, [rsi + rax*8 + 8] 
    movlpd xmm2, [rsi + rax*8 + 16] 
    movhpd xmm0, [rsi + rbx*8] 
    movhpd xmm1, [rsi + rbx*8 + 8] 
    movhpd xmm2, [rsi + rbx*8 + 16] 

    ;# xmm0 = jx
    ;# xmm1 = jy
    ;# xmm2 = jz
        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb231_ixO]
    subpd xmm1, [rsp + nb231_iyO]
    subpd xmm2, [rsp + nb231_izO]
    subpd xmm3, [rsp + nb231_ixH1]
    subpd xmm4, [rsp + nb231_iyH1]
    subpd xmm5, [rsp + nb231_izH1]
    subpd xmm6, [rsp + nb231_ixH2]
    subpd xmm7, [rsp + nb231_iyH2]
    subpd xmm8, [rsp + nb231_izH2]
    
	movapd [rsp + nb231_dxO], xmm0
	movapd [rsp + nb231_dyO], xmm1
	movapd [rsp + nb231_dzO], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb231_dxH1], xmm3
	movapd [rsp + nb231_dyH1], xmm4
	movapd [rsp + nb231_dzH1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb231_dxH2], xmm6
	movapd [rsp + nb231_dyH2], xmm7
	movapd [rsp + nb231_dzH2], xmm8
	mulpd  xmm6, xmm6
	mulpd  xmm7, xmm7
	mulpd  xmm8, xmm8
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    addpd  xmm6, xmm7
    addpd  xmm6, xmm8

	;# start doing invsqrt for j atoms
    cvtpd2ps xmm1, xmm0
    cvtpd2ps xmm4, xmm3
    cvtpd2ps xmm7, xmm6
	rsqrtps xmm1, xmm1
	rsqrtps xmm4, xmm4
    rsqrtps xmm7, xmm7
    cvtps2pd xmm1, xmm1
    cvtps2pd xmm4, xmm4
    cvtps2pd xmm7, xmm7
	
	movapd  xmm2, xmm1
	movapd  xmm5, xmm4
    movapd  xmm8, xmm7
    
	mulpd   xmm1, xmm1 ;# lu*lu
	mulpd   xmm4, xmm4 ;# lu*lu
    mulpd   xmm7, xmm7 ;# lu*lu
		
	movapd  xmm9, [rsp + nb231_three]
	movapd  xmm10, xmm9
    movapd  xmm11, xmm9

	mulpd   xmm1, xmm0 ;# rsq*lu*lu
	mulpd   xmm4, xmm3 ;# rsq*lu*lu 
    mulpd   xmm7, xmm6 ;# rsq*lu*lu
	
	subpd   xmm9, xmm1
	subpd   xmm10, xmm4
    subpd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulpd   xmm9, xmm2
	mulpd   xmm10, xmm5
    mulpd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb231_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvO
	mulpd   xmm10, xmm15 ;# first iteration for rinvH1
    mulpd   xmm11, xmm15 ;# first iteration for rinvH2

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb231_three]
	movapd  xmm4, xmm1
    movapd  xmm7, xmm1

	mulpd   xmm2, xmm0 ;# rsq*lu*lu
	mulpd   xmm5, xmm3 ;# rsq*lu*lu 
    mulpd   xmm8, xmm6 ;# rsq*lu*lu
	
	subpd   xmm1, xmm2
	subpd   xmm4, xmm5
    subpd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulpd   xmm9, xmm1
	mulpd   xmm10, xmm4
    mulpd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb231_half]
	mulpd   xmm9, xmm15  ;#  rinvO
	mulpd   xmm10, xmm15 ;#   rinvH1
    mulpd   xmm11, xmm15 ;#   rinvH2
	
	;# O interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11
    movapd [rsp + nb231_rsqO], xmm0
    movapd [rsp + nb231_rsqH1], xmm3
    movapd [rsp + nb231_rsqH2], xmm6
    movapd [rsp + nb231_rinvO], xmm9
    movapd [rsp + nb231_rinvH1], xmm10
    movapd [rsp + nb231_rinvH2], xmm11

    ;# table LJ interaction
    mulpd  xmm0, xmm9
    mulpd  xmm0, [rsp + nb231_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttpd2dq xmm1, xmm0

    ;# convert back to float
    cvtdq2pd  xmm2, xmm1
         
    ;# multiply by 8
    pslld   xmm1, 3

    ;# move to integer registers
    pshufd  xmm13, xmm1, 1
    movd    r8d, xmm1
    movd    r10d, xmm13
    
    ;# calculate eps
    subpd     xmm0, xmm2
    mov  rsi, [rbp + nb231_VFtab]
            
    movlpd xmm4, [rsi + r8*8]
   	movlpd xmm5, [rsi + r8*8 + 8]
	movlpd xmm6, [rsi + r8*8 + 16]
	movlpd xmm7, [rsi + r8*8 + 24]
    movlpd xmm8, [rsi + r8*8 + 32]
   	movlpd xmm9, [rsi + r8*8 + 40]
	movlpd xmm10, [rsi + r8*8 + 48]
	movlpd xmm11, [rsi + r8*8 + 56]
    
    movhpd xmm4, [rsi + r10*8]
   	movhpd xmm5, [rsi + r10*8 + 8]
	movhpd xmm6, [rsi + r10*8 + 16]
	movhpd xmm7, [rsi + r10*8 + 24]
    movhpd xmm8, [rsi + r10*8 + 32]
   	movhpd xmm9, [rsi + r10*8 + 40]
	movhpd xmm10, [rsi + r10*8 + 48]
	movhpd xmm11, [rsi + r10*8 + 56]
    ;# dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

    mulpd  xmm7, xmm0    ;# Heps
    mulpd  xmm11, xmm0 
    mulpd  xmm6, xmm0   ;# Geps
    mulpd  xmm10, xmm0 
    mulpd  xmm7, xmm0   ;# Heps2
    mulpd  xmm11, xmm0 
    addpd  xmm5, xmm6  ;# F+Geps
    addpd  xmm9, xmm10 
    addpd  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addpd  xmm9, xmm11 
    addpd  xmm7, xmm7    ;# 2*Heps2
    addpd  xmm11, xmm11
    addpd  xmm7, xmm6   ;# 2*Heps2+Geps
    addpd  xmm11, xmm10
    
    addpd  xmm7, xmm5  ;# FF = Fp + 2*Heps2 + Geps
    addpd  xmm11, xmm9
    mulpd  xmm5, xmm0  ;# eps*Fp
    mulpd  xmm9, xmm0
    movapd xmm12, [rsp + nb231_c6]
    movapd xmm13, [rsp + nb231_c12]
    addpd  xmm5, xmm4 ;# VV
    addpd  xmm9, xmm8

    mulpd  xmm5, xmm12  ;# VV*c6 = vnb6
    mulpd  xmm9, xmm13  ;# VV*c12 = vnb12
    addpd  xmm5, xmm9
    addpd  xmm5, [rsp + nb231_Vvdwtot]
    movapd [rsp + nb231_Vvdwtot], xmm5
        
    mulpd  xmm7, xmm12   ;# FF*c6 = fnb6
    mulpd  xmm11, xmm13   ;# FF*c12  = fnb12
    addpd  xmm7, xmm11
    mulpd  xmm7, [rsp + nb231_tsc]
    movapd [rsp + nb231_fstmp], xmm7 
           
    ;# Coulomb reaction-field interaction
    movapd xmm0, [rsp + nb231_rsqO]
    movapd xmm3, [rsp + nb231_rsqH1]
    movapd xmm6, [rsp + nb231_rsqH2]
    movapd xmm9, [rsp + nb231_rinvO]
    movapd xmm10, [rsp + nb231_rinvH1]
    movapd xmm11, [rsp + nb231_rinvH2]
    
    movapd xmm1, xmm9 ;# copy of rinv
    movapd xmm4, xmm10
    movapd xmm7, xmm11
    movapd xmm2, [rsp + nb231_krf]    
    mulpd  xmm10, xmm10 ;# rinvsq
    mulpd  xmm11, xmm11
    mulpd  xmm0, xmm2  ;# k*rsq
    mulpd  xmm3, xmm2
    mulpd  xmm6, xmm2
    movapd xmm2, xmm0 ;# copy of k*rsq
    movapd xmm5, xmm3
    movapd xmm8, xmm6
    addpd  xmm2, xmm1  ;# rinv+krsq
    addpd  xmm5, xmm4
    addpd  xmm8, xmm7
    subpd  xmm2, [rsp + nb231_crf]   ;# rinv+krsq-crf
    subpd  xmm5, [rsp + nb231_crf]
    subpd  xmm8, [rsp + nb231_crf]   
    mulpd  xmm2, [rsp + nb231_qqO] ;# voul=qq*(rinv+ krsq-crf)
    mulpd  xmm5, [rsp + nb231_qqH] ;# voul=qq*(rinv+ krsq-crf)
    mulpd  xmm8, [rsp + nb231_qqH] ;# voul=qq*(rinv+ krsq-crf)
    addpd  xmm0, xmm0 ;# 2*krsq
    addpd  xmm3, xmm3 
    addpd  xmm6, xmm6 
    subpd  xmm1, xmm0 ;# rinv-2*krsq
    subpd  xmm4, xmm3
    subpd  xmm7, xmm6
    mulpd  xmm1, [rsp + nb231_qqO]   ;# (rinv-2*krsq)*qq
    mulpd  xmm4, [rsp + nb231_qqH] 
    mulpd  xmm7, [rsp + nb231_qqH] 
    addpd  xmm2, [rsp + nb231_vctot]
    addpd  xmm5, xmm8
    addpd  xmm2, xmm5
    movapd [rsp + nb231_vctot], xmm2
    mulpd  xmm1, xmm9   ;# fijC
    mulpd  xmm10, xmm4
    mulpd  xmm11, xmm7    
    ;# xmm15 contains vctot
    subpd xmm1, [rsp + nb231_fstmp]
    mulpd xmm9, xmm1    

    ;# move j forces to xmm0-xmm2
    mov rdi, [rbp + nb231_faction]
	movlpd xmm0, [rdi + rax*8]
	movlpd xmm1, [rdi + rax*8 + 8]
	movlpd xmm2, [rdi + rax*8 + 16]
	movhpd xmm0, [rdi + rbx*8]
	movhpd xmm1, [rdi + rbx*8 + 8]
	movhpd xmm2, [rdi + rbx*8 + 16]

    movapd xmm7, xmm9
    movapd xmm8, xmm9
    movapd xmm13, xmm11
    movapd xmm14, xmm11
    movapd xmm15, xmm11
    movapd xmm11, xmm10
    movapd xmm12, xmm10

	mulpd xmm7, [rsp + nb231_dxO]
	mulpd xmm8, [rsp + nb231_dyO]
	mulpd xmm9, [rsp + nb231_dzO]
	mulpd xmm10, [rsp + nb231_dxH1]
	mulpd xmm11, [rsp + nb231_dyH1]
	mulpd xmm12, [rsp + nb231_dzH1]
	mulpd xmm13, [rsp + nb231_dxH2]
	mulpd xmm14, [rsp + nb231_dyH2]
	mulpd xmm15, [rsp + nb231_dzH2]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb231_fixO]
    addpd xmm8, [rsp + nb231_fiyO]
    addpd xmm9, [rsp + nb231_fizO]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb231_fixH1]
    addpd xmm11, [rsp + nb231_fiyH1]
    addpd xmm12, [rsp + nb231_fizH1]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb231_fixH2]
    addpd xmm14, [rsp + nb231_fiyH2]
    addpd xmm15, [rsp + nb231_fizH2]

    movapd [rsp + nb231_fixO], xmm7
    movapd [rsp + nb231_fiyO], xmm8
    movapd [rsp + nb231_fizO], xmm9
    movapd [rsp + nb231_fixH1], xmm10
    movapd [rsp + nb231_fiyH1], xmm11
    movapd [rsp + nb231_fizH1], xmm12
    movapd [rsp + nb231_fixH2], xmm13
    movapd [rsp + nb231_fiyH2], xmm14
    movapd [rsp + nb231_fizH2], xmm15
   
    ;# store back j forces from xmm0-xmm2
	movlpd [rdi + rax*8],      xmm0
	movlpd [rdi + rax*8 + 8],  xmm1
	movlpd [rdi + rax*8 + 16], xmm2
	movhpd [rdi + rbx*8],      xmm0
	movhpd [rdi + rbx*8 + 8],  xmm1
	movhpd [rdi + rbx*8 + 16], xmm2

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb231_innerk],  2
	jl    .nb231_checksingle
	jmp   .nb231_unroll_loop
.nb231_checksingle:	
	mov   edx, [rsp + nb231_innerk]
	and   edx, 1
	jnz   .nb231_dosingle
	jmp   .nb231_updateouterdata
.nb231_dosingle:
	mov   rdx, [rsp + nb231_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb231_innerjjnr],  4	

	mov rsi, [rbp + nb231_charge]    ;# base of charge[] 
	xorpd xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb231_iqO]
	mulpd  xmm4, [rsp + nb231_iqH]

	movapd  [rsp + nb231_qqO], xmm3
	movapd  [rsp + nb231_qqH], xmm4
	
	mov rsi, [rbp + nb231_type]
	mov r8d, [rsi + rax*4]
	mov rsi, [rbp + nb231_vdwparam]
	shl r8d, 1	
	mov edi, [rsp + nb231_ntia]
	add r8d, edi

	movlpd xmm6, [rsi + r8*8]	;# c6a
	movhpd xmm6, [rsi + r8*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movapd [rsp + nb231_c6], xmm4
	movapd [rsp + nb231_c12], xmm6
	
	mov rsi, [rbp + nb231_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	
	movlpd xmm4, [rsi + rax*8]
	movlpd xmm5, [rsi + rax*8 + 8]
	movlpd xmm6, [rsi + rax*8 + 16]
    movapd xmm0, xmm4
    movapd xmm1, xmm5
    movapd xmm2, xmm6

	;# calc dr 
	subsd xmm4, [rsp + nb231_ixO]
	subsd xmm5, [rsp + nb231_iyO]
	subsd xmm6, [rsp + nb231_izO]

	;# store dr 
	movapd [rsp + nb231_dxO], xmm4
	movapd [rsp + nb231_dyO], xmm5
	movapd [rsp + nb231_dzO], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 
	movapd [rsp + nb231_rsqO], xmm7
	
	;# move j coords to xmm4-xmm6 
	movapd xmm4, xmm0
	movapd xmm5, xmm1
	movapd xmm6, xmm2

	;# calc dr 
	subsd xmm4, [rsp + nb231_ixH1]
	subsd xmm5, [rsp + nb231_iyH1]
	subsd xmm6, [rsp + nb231_izH1]

	;# store dr 
	movapd [rsp + nb231_dxH1], xmm4
	movapd [rsp + nb231_dyH1], xmm5
	movapd [rsp + nb231_dzH1], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move j coords to xmm3-xmm5
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	;# calc dr 
	subsd xmm3, [rsp + nb231_ixH2]
	subsd xmm4, [rsp + nb231_iyH2]
	subsd xmm5, [rsp + nb231_izH2]

	;# store dr 
	movapd [rsp + nb231_dxH2], xmm3
	movapd [rsp + nb231_dyH2], xmm4
	movapd [rsp + nb231_dzH2], xmm5
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

	movapd [rsp + nb231_rsqH2], xmm0
	movapd [rsp + nb231_rsqH1], xmm1
	movapd [rsp + nb231_rsqO], xmm2

	;# start with rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb231_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb231_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb231_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb231_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb231_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb231_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb231_three]
	subsd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb231_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	movsd  [rsp + nb231_rinvH1], xmm6
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb231_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb231_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb231_three]
	subsd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb231_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 
	movsd [rsp + nb231_rinvH2], xmm5
	
	;# do O interactions 
	movsd xmm0, xmm7
	movsd xmm6, xmm7
	movsd xmm1, [rsp + nb231_rsqO]
	mulsd  xmm1, [rsp + nb231_krf]	

	addsd  xmm0, xmm1
	mulsd  xmm1, [rsp + nb231_two]
	subsd  xmm0, [rsp + nb231_crf] ;# xmm0=rinv+ krsq-crf 
	subsd  xmm6, xmm1
	mulsd  xmm0, [rsp + nb231_qqO]
	mulsd  xmm6, [rsp + nb231_qqO]

	mulsd  xmm6, xmm7
	movsd [rsp + nb231_fstmp], xmm6 ;# save to temp. storage

	addsd  xmm0, [rsp + nb231_vctot]
	movsd [rsp + nb231_vctot], xmm0

	movsd xmm0, xmm7
	movsd xmm4, [rsp + nb231_rsqO]
	;# LJ table interaction. xmm0=rinv, xmm4=rsq
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [rsp + nb231_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3
	
	mov  rsi, [rbp + nb231_VFtab]

	;# dispersion 
	movlpd xmm4, [rsi + rbx*8]	;# Y1 	
	movhpd xmm4, [rsi + rbx*8 + 8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rbx*8 + 16]	;# G1
	movhpd xmm6, [rsi + rbx*8 + 24]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# dispersion table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [rsp + nb231_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movsd xmm4, [rsp + nb231_c6]
	mulsd  xmm7, xmm4	 ;# fijD 
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addsd  xmm5, [rsp + nb231_Vvdwtot]
	movsd xmm3, [rsp + nb231_fstmp]
	mulsd  xmm7, [rsp + nb231_tsc]
	subsd  xmm3, xmm7
	movsd [rsp + nb231_fstmp], xmm3
	movsd [rsp + nb231_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [rsi + rbx*8 + 32]	;# Y1 	
	movhpd xmm4, [rsi + rbx*8 + 40]	;# Y1 F1 	

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rbx*8 + 48]	;# G1
	movhpd xmm6, [rsi + rbx*8 + 56]	;# G1 H1 	

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	
	;# table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [rsp + nb231_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	
	movsd xmm4, [rsp + nb231_c12]
	mulsd  xmm7, xmm4 
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [rsp + nb231_Vvdwtot]
	movsd xmm3, [rsp + nb231_fstmp]
	mulsd  xmm7, [rsp + nb231_tsc]
	subsd  xmm3, xmm7
	movsd [rsp + nb231_Vvdwtot], xmm5

	mulsd  xmm3, xmm0
		
	movsd xmm0, [rsp + nb231_dxO]
	movsd xmm1, [rsp + nb231_dyO]
	movsd xmm2, [rsp + nb231_dzO]

	mov    rdi, [rbp + nb231_faction]
	mulsd  xmm0, xmm3
	mulsd  xmm1, xmm3
	mulsd  xmm2, xmm3

	;# update O forces 
	movapd xmm3, [rsp + nb231_fixO]
	movapd xmm4, [rsp + nb231_fiyO]
	movapd xmm7, [rsp + nb231_fizO]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [rsp + nb231_fixO], xmm3
	movlpd [rsp + nb231_fiyO], xmm4
	movlpd [rsp + nb231_fizO], xmm7
	;# update j forces with water O 
	movlpd [rsp + nb231_fjx], xmm0
	movlpd [rsp + nb231_fjy], xmm1
	movlpd [rsp + nb231_fjz], xmm2

	;# H1 interactions 
	movsd   xmm6, [rsp + nb231_rinvH1]	
	movsd   xmm4, xmm6
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movsd   xmm7, xmm6
	movsd   xmm0, [rsp + nb231_rsqH1]
   	mulsd  xmm0, [rsp + nb231_krf]	

	addsd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	mulsd   xmm0, [rsp + nb231_two]
	subsd   xmm6, [rsp + nb231_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm6, [rsp + nb231_qqH] ;# vcoul 
	mulsd   xmm7, [rsp + nb231_qqH]
	mulsd  xmm4, xmm7		;# total fsH1 in xmm4 
	
	addsd  xmm6, [rsp + nb231_vctot]

	movapd xmm0, [rsp + nb231_dxH1]
	movapd xmm1, [rsp + nb231_dyH1]
	movapd xmm2, [rsp + nb231_dzH1]
	movlpd [rsp + nb231_vctot], xmm6
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [rsp + nb231_fixH1]
	movapd xmm4, [rsp + nb231_fiyH1]
	movapd xmm7, [rsp + nb231_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [rsp + nb231_fixH1], xmm3
	movlpd [rsp + nb231_fiyH1], xmm4
	movlpd [rsp + nb231_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [rsp + nb231_fjx]
	addsd  xmm1, [rsp + nb231_fjy]
	addsd  xmm2, [rsp + nb231_fjz]
	movlpd [rsp + nb231_fjx], xmm0
	movlpd [rsp + nb231_fjy], xmm1
	movlpd [rsp + nb231_fjz], xmm2

	;# H2 interactions 
	movapd  xmm5, [rsp + nb231_rinvH2]	
	movapd  xmm4, xmm5
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm5
	movsd   xmm0, [rsp + nb231_rsqH2]
	mulsd  xmm0, [rsp + nb231_krf]	

	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulsd   xmm0, [rsp + nb231_two]
	subsd   xmm5, [rsp + nb231_crf]
	subsd   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulsd   xmm5, [rsp + nb231_qqH] ;# vcoul 
	mulsd   xmm7, [rsp + nb231_qqH]
	mulsd  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addsd  xmm5, [rsp + nb231_vctot]

	movapd xmm0, [rsp + nb231_dxH2]
	movapd xmm1, [rsp + nb231_dyH2]
	movapd xmm2, [rsp + nb231_dzH2]
	movlpd [rsp + nb231_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [rsp + nb231_fixH2]
	movapd xmm4, [rsp + nb231_fiyH2]
	movapd xmm7, [rsp + nb231_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [rsp + nb231_fixH2], xmm3
	movlpd [rsp + nb231_fiyH2], xmm4
	movlpd [rsp + nb231_fizH2], xmm7

	mov rdi, [rbp + nb231_faction]
	;# update j forces 
	addsd  xmm0, [rsp + nb231_fjx]
	addsd  xmm1, [rsp + nb231_fjy]
	addsd  xmm2, [rsp + nb231_fjz]
	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]
	addsd xmm3, xmm0
	addsd xmm4, xmm1
	addsd xmm5, xmm2
	movlpd [rdi + rax*8], xmm3
	movlpd [rdi + rax*8 + 8], xmm4
	movlpd [rdi + rax*8 + 16], xmm5	

.nb231_updateouterdata:
	mov   ecx, [rsp + nb231_ii3]
	mov   rdi, [rbp + nb231_faction]
	mov   rsi, [rbp + nb231_fshift]
	mov   edx, [rsp + nb231_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb231_fixO]
	movapd xmm1, [rsp + nb231_fiyO]
	movapd xmm2, [rsp + nb231_fizO]

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
	movsd  xmm3, [rdi + rcx*8]
	movsd  xmm4, [rdi + rcx*8 + 8]
	movsd  xmm5, [rdi + rcx*8 + 16]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8],     xmm3
	movsd  [rdi + rcx*8 + 8], xmm4
	movsd  [rdi + rcx*8 + 16], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movapd xmm6, xmm0
	movsd xmm7, xmm2
	unpcklpd xmm6, xmm1

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb231_fixH1]
	movapd xmm1, [rsp + nb231_fiyH1]
	movapd xmm2, [rsp + nb231_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	;# increment i force 
	movsd  xmm3, [rdi + rcx*8 + 24]
	movsd  xmm4, [rdi + rcx*8 + 32]
	movsd  xmm5, [rdi + rcx*8 + 40]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8 + 24], xmm3
	movsd  [rdi + rcx*8 + 32], xmm4
	movsd  [rdi + rcx*8 + 40], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb231_fixH2]
	movapd xmm1, [rsp + nb231_fiyH2]
	movapd xmm2, [rsp + nb231_fizH2]

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
	movsd  xmm3, [rdi + rcx*8 + 48]
	movsd  xmm4, [rdi + rcx*8 + 56]
	movsd  xmm5, [rdi + rcx*8 + 64]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8 + 48], xmm3
	movsd  [rdi + rcx*8 + 56], xmm4
	movsd  [rdi + rcx*8 + 64], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# increment fshift force 
	movlpd xmm3, [rsi + rdx*8]
	movhpd xmm3, [rsi + rdx*8 + 8]
	movsd  xmm4, [rsi + rdx*8 + 16]
	subpd  xmm3, xmm6
	subsd  xmm4, xmm7
	movlpd [rsi + rdx*8],      xmm3
	movhpd [rsi + rdx*8 + 8],  xmm3
	movsd  [rsi + rdx*8 + 16], xmm4

	;# get n from stack
	mov esi, [rsp + nb231_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb231_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb231_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb231_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb231_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   rax, [rbp + nb231_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [rsp + nb231_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb231_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb231_n], esi
        jmp .nb231_outer
.nb231_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb231_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb231_end
        ;# non-zero, do one more workunit
        jmp   .nb231_threadloop
.nb231_end:
	mov eax, [rsp + nb231_nouter]
	mov ebx, [rsp + nb231_ninner]
	mov rcx, [rbp + nb231_outeriter]
	mov rdx, [rbp + nb231_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 912
	emms

    ;# Save xmm registers to stack
    movaps xmm6,  [rsp      ]
    movaps xmm7,  [rsp + 16 ]
    movaps xmm8,  [rsp + 32 ]
    movaps xmm9,  [rsp + 48 ]
    movaps xmm10, [rsp + 64 ]
    movaps xmm11, [rsp + 80 ]
    movaps xmm12, [rsp + 96 ]
    movaps xmm13, [rsp + 112]
    movaps xmm14, [rsp + 128]
    movaps xmm15, [rsp + 144]

    ;# Reset pointers after restoring xmm6-15
    add rsp, 168

    pop r15
    pop r14
    pop r13
    pop r12
    pop rdi
    pop rsi
    pop rbx
    
	pop	rbp
	ret




.globl nb_kernel231nf_x86_64_sse2
.globl _nb_kernel231nf_x86_64_sse2
nb_kernel231nf_x86_64_sse2:	
_nb_kernel231nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb231nf_fshift,         16
.equiv          nb231nf_gid,            24
.equiv          nb231nf_pos,            32
.equiv          nb231nf_faction,        40
.equiv          nb231nf_charge,         48
.equiv          nb231nf_p_facel,        56
.equiv          nb231nf_argkrf,         64
.equiv          nb231nf_argcrf,         72
.equiv          nb231nf_Vc,             80
.equiv          nb231nf_type,           88
.equiv          nb231nf_p_ntype,        96
.equiv          nb231nf_vdwparam,       104
.equiv          nb231nf_Vvdw,           112
.equiv          nb231nf_p_tabscale,     120
.equiv          nb231nf_VFtab,          128
.equiv          nb231nf_invsqrta,       136
.equiv          nb231nf_dvda,           144
.equiv          nb231nf_p_gbtabscale,   152
.equiv          nb231nf_GBtab,          160
.equiv          nb231nf_p_nthreads,     168
.equiv          nb231nf_count,          176
.equiv          nb231nf_mtx,            184
.equiv          nb231nf_outeriter,      192
.equiv          nb231nf_inneriter,      200
.equiv          nb231nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb231nf_ixO,            0
.equiv          nb231nf_iyO,            16
.equiv          nb231nf_izO,            32
.equiv          nb231nf_ixH1,           48
.equiv          nb231nf_iyH1,           64
.equiv          nb231nf_izH1,           80
.equiv          nb231nf_ixH2,           96
.equiv          nb231nf_iyH2,           112
.equiv          nb231nf_izH2,           128
.equiv          nb231nf_iqO,            144
.equiv          nb231nf_iqH,            160
.equiv          nb231nf_qqO,            176
.equiv          nb231nf_qqH,            192
.equiv          nb231nf_c6,             208
.equiv          nb231nf_c12,            224
.equiv          nb231nf_vctot,          240
.equiv          nb231nf_Vvdwtot,        256
.equiv          nb231nf_half,           272
.equiv          nb231nf_three,          288
.equiv          nb231nf_krf,            304
.equiv          nb231nf_crf,            320
.equiv          nb231nf_rsqO,           336
.equiv          nb231nf_rinvH1,         352
.equiv          nb231nf_rinvH2,         368
.equiv          nb231nf_krsqO,          384
.equiv          nb231nf_krsqH1,         400
.equiv          nb231nf_krsqH2,         416
.equiv          nb231nf_tsc,            432
.equiv          nb231nf_facel,          448
.equiv          nb231nf_iinr,           456
.equiv          nb231nf_jindex,         464
.equiv          nb231nf_jjnr,           472
.equiv          nb231nf_shift,          480
.equiv          nb231nf_shiftvec,       488
.equiv          nb231nf_innerjjnr,      496
.equiv          nb231nf_nri,            504
.equiv          nb231nf_ntia,           508
.equiv          nb231nf_is3,            512
.equiv          nb231nf_ii3,            516
.equiv          nb231nf_innerk,         520
.equiv          nb231nf_n,              524
.equiv          nb231nf_nn1,            528
.equiv          nb231nf_nouter,         532
.equiv          nb231nf_ninner,         536

	push rbp
	mov  rbp, rsp
    
    ;# Push integer registers on stack
	push rbx
    push rsi
    push rdi
    push r12
    push r13
    push r14
    push r15

    ;# Make room for registers xmm6-xmm15 (10 registers=160 bytes)
    sub rsp, 168
    
    ;# Save xmm registers to stack
    movaps [rsp      ], xmm6
    movaps [rsp + 16 ], xmm7
    movaps [rsp + 32 ], xmm8
    movaps [rsp + 48 ], xmm9
    movaps [rsp + 64 ], xmm10
    movaps [rsp + 80 ], xmm11
    movaps [rsp + 96 ], xmm12
    movaps [rsp + 112], xmm13
    movaps [rsp + 128], xmm14
    movaps [rsp + 144], xmm15
    
; .if 0    # block below only read by NASM - special calling convention on win64
%ifidn __OUTPUT_FORMAT__, win64
    ;# Adjust rbp to account for shadow space (32) & two extra args (2*8) on stack
    add rbp, 48
    ;# Adjust stack pointer for different alignment
    ;# Move around arguments to fit AMD64 convention below
    ;# AMD64 passes args in: rdi,rsi,rdx,rcx,r8,r9 + stack
    ;# win64 passes args in: rcx,rdx,r8,r9         + stack
    mov rdi, rcx
    mov rsi, rdx
    mov rdx, r8
    mov rcx, r9
    mov r8,  [rbp]
    mov r9,  [rbp + 8]
%endif
; .endif   # end NASM- and win64-specific block

	emms
	sub rsp, 544		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb231nf_nouter], eax
	mov [rsp + nb231nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb231nf_nri], edi
	mov [rsp + nb231nf_iinr], rsi
	mov [rsp + nb231nf_jindex], rdx
	mov [rsp + nb231nf_jjnr], rcx
	mov [rsp + nb231nf_shift], r8
	mov [rsp + nb231nf_shiftvec], r9
	mov rsi, [rbp + nb231nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb231nf_facel], xmm0

	mov rax, [rbp + nb231nf_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb231nf_tsc], xmm3

	mov rsi, [rbp + nb231nf_argkrf]
	mov rdi, [rbp + nb231nf_argcrf]
	movsd xmm1, [rsi]
	movsd xmm2, [rdi]
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	movapd [rsp + nb231nf_krf], xmm1
	movapd [rsp + nb231nf_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb231nf_half], eax
	mov [rsp + nb231nf_half+4], ebx
	movsd xmm1, [rsp + nb231nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb231nf_half], xmm1
	movapd [rsp + nb231nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb231nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb231nf_charge]
	movsd xmm3, [rdx + rbx*8]	
	movsd xmm4, [rdx + rbx*8 + 8]	
	movsd xmm5, [rsp + nb231nf_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb231nf_iqO], xmm3
	movapd [rsp + nb231nf_iqH], xmm4
	
	mov   rdx, [rbp + nb231nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb231nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb231nf_ntia], ecx		
.nb231nf_threadloop:
        mov   rsi, [rbp + nb231nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb231nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb231nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb231nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb231nf_n], eax
        mov [rsp + nb231nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb231nf_outerstart
        jmp .nb231nf_end

.nb231nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb231nf_nouter]
	mov [rsp + nb231nf_nouter], ebx

.nb231nf_outer:
	mov   rax, [rsp + nb231nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb231nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb231nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb231nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb231nf_pos]    ;# eax = base of pos[]  
	mov   [rsp + nb231nf_ii3], ebx

	addsd xmm3, [rax + rbx*8]
	addsd xmm4, [rax + rbx*8 + 8]
	addsd xmm5, [rax + rbx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb231nf_ixO], xmm3
	movapd [rsp + nb231nf_iyO], xmm4
	movapd [rsp + nb231nf_izO], xmm5

	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm0, [rax + rbx*8 + 24]
	addsd xmm1, [rax + rbx*8 + 32]
	addsd xmm2, [rax + rbx*8 + 40]		
	addsd xmm3, [rax + rbx*8 + 48]
	addsd xmm4, [rax + rbx*8 + 56]
	addsd xmm5, [rax + rbx*8 + 64]		

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb231nf_ixH1], xmm0
	movapd [rsp + nb231nf_iyH1], xmm1
	movapd [rsp + nb231nf_izH1], xmm2
	movapd [rsp + nb231nf_ixH2], xmm3
	movapd [rsp + nb231nf_iyH2], xmm4
	movapd [rsp + nb231nf_izH2], xmm5
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb231nf_vctot], xmm4
	movapd [rsp + nb231nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb231nf_jindex]
	mov   ecx, [rax+rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb231nf_pos]
	mov   rdi, [rbp + nb231nf_faction]	
	mov   rax, [rsp + nb231nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb231nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb231nf_ninner]
	mov   [rsp + nb231nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb231nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb231nf_unroll_loop
	jmp   .nb231nf_checksingle
.nb231nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb231nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]

	add qword ptr [rsp + nb231nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb231nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb231nf_iqO]
	mulpd  xmm4, [rsp + nb231nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx

	movapd  [rsp + nb231nf_qqO], xmm3
	movapd  [rsp + nb231nf_qqH], xmm4
	
	mov rsi, [rbp + nb231nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov rsi, [rbp + nb231nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [rsp + nb231nf_ntia]
	add eax, edi
	add ebx, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movlpd xmm7, [rsi + rbx*8]	;# c6b
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 
	movhpd xmm7, [rsi + rbx*8 + 8]	;# c6b c12b 
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [rsp + nb231nf_c6], xmm4
	movapd [rsp + nb231nf_c12], xmm6
	
	mov rsi, [rbp + nb231nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	movhpd xmm0, [rsi + rbx*8]
	movhpd xmm1, [rsi + rbx*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 16]		

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [rsp + nb231nf_ixO]
	movapd xmm5, [rsp + nb231nf_iyO]
	movapd xmm6, [rsp + nb231nf_izO]

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
	movapd [rsp + nb231nf_rsqO], xmm7
	
	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [rsp + nb231nf_ixH1]
	movapd xmm5, [rsp + nb231nf_iyH1]
	movapd xmm6, [rsp + nb231nf_izH1]

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
	movapd xmm3, [rsp + nb231nf_ixH2]
	movapd xmm4, [rsp + nb231nf_iyH2]
	movapd xmm5, [rsp + nb231nf_izH2]

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

	mulpd  xmm0, [rsp + nb231nf_krf]	
	mulpd  xmm1, [rsp + nb231nf_krf]	
	mulpd  xmm2, [rsp + nb231nf_krf]	

	movapd [rsp + nb231nf_krsqH2], xmm0
	movapd [rsp + nb231nf_krsqH1], xmm1
	movapd [rsp + nb231nf_krsqO], xmm2
		
	;# start with rsqO - put seed in xmm2 
	cvtpd2ps xmm2, xmm7	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [rsp + nb231nf_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [rsp + nb231nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb231nf_three]
	subpd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [rsp + nb231nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [rsp + nb231nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [rsp + nb231nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb231nf_three]
	subpd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [rsp + nb231nf_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	movapd  [rsp + nb231nf_rinvH1], xmm6
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [rsp + nb231nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [rsp + nb231nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb231nf_three]
	subpd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [rsp + nb231nf_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 
	movapd  [rsp + nb231nf_rinvH2], xmm5

	;# do O interactions 
	movapd xmm0, xmm7
	movapd xmm1, [rsp + nb231nf_krsqO]
	addpd  xmm0, xmm1
	subpd  xmm0, [rsp + nb231nf_crf] ;# xmm0=rinv+ krsq-crf 
	mulpd  xmm0, [rsp + nb231nf_qqO]

	mulpd  xmm6, xmm7

	addpd  xmm0, [rsp + nb231nf_vctot]
	movapd [rsp + nb231nf_vctot], xmm0

	movapd xmm0, xmm7
	movapd xmm4, [rsp + nb231nf_rsqO]
	;# LJ table interaction. xmm0=rinv, xmm4=rsq
	
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [rsp + nb231nf_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	mov  rsi, [rbp + nb231nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6

	;# dispersion 
	movlpd xmm4, [rsi + rax*8]	;# Y1 	
	movlpd xmm3, [rsi + rbx*8]	;# Y2 
	movhpd xmm4, [rsi + rax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [rsi + rbx*8 + 8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rax*8 + 16]	;# G1
	movlpd xmm3, [rsi + rbx*8 + 16]	;# G2
	movhpd xmm6, [rsi + rax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [rsi + rbx*8 + 24]	;# G2 H2 
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

	movapd xmm4, [rsp + nb231nf_c6]
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;#  Update Vvdwtot directly 
	addpd  xmm5, [rsp + nb231nf_Vvdwtot]
	movapd [rsp + nb231nf_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [rsi + rax*8 + 32]	;# Y1 	
	movlpd xmm3, [rsi + rbx*8 + 32]	;# Y2 
	movhpd xmm4, [rsi + rax*8 + 40]	;# Y1 F1 	
	movhpd xmm3, [rsi + rbx*8 + 40]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rax*8 + 48]	;# G1
	movlpd xmm3, [rsi + rbx*8 + 48]	;# G2
	movhpd xmm6, [rsi + rax*8 + 56]	;# G1 H1 	
	movhpd xmm3, [rsi + rbx*8 + 56]	;# G2 H2 

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
	
	movapd xmm4, [rsp + nb231nf_c12]
	mulpd  xmm5, xmm4  
	
	addpd  xmm5, [rsp + nb231nf_Vvdwtot]
	movapd [rsp + nb231nf_Vvdwtot], xmm5

	;# H1 interactions 
	movapd  xmm6, [rsp + nb231nf_rinvH1]	
	movapd  xmm4, xmm6
	mulpd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm6
	movapd  xmm0, [rsp + nb231nf_krsqH1]
	addpd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	subpd   xmm6, [rsp + nb231nf_crf]
	mulpd   xmm6, [rsp + nb231nf_qqH] ;# vcoul 
	
	addpd  xmm6, [rsp + nb231nf_vctot]
	movapd [rsp + nb231nf_vctot], xmm6

	;# H2 interactions 
	movapd  xmm5, [rsp + nb231nf_rinvH2]	
	movapd  xmm4, xmm5
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm5
	movapd  xmm0, [rsp + nb231nf_krsqH2]
	addpd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	subpd   xmm5, [rsp + nb231nf_crf]
	mulpd   xmm5, [rsp + nb231nf_qqH] ;# vcoul 
	
	addpd  xmm5, [rsp + nb231nf_vctot]
	movapd [rsp + nb231nf_vctot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb231nf_innerk],  2
	jl    .nb231nf_checksingle
	jmp   .nb231nf_unroll_loop
.nb231nf_checksingle:	
	mov   edx, [rsp + nb231nf_innerk]
	and   edx, 1
	jnz   .nb231nf_dosingle
	jmp   .nb231nf_updateouterdata
.nb231nf_dosingle:
	mov   rdx, [rsp + nb231nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb231nf_innerjjnr],  4

	mov rsi, [rbp + nb231nf_charge]    ;# base of charge[] 
	xorpd xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb231nf_iqO]
	mulpd  xmm4, [rsp + nb231nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 

	movapd  [rsp + nb231nf_qqO], xmm3
	movapd  [rsp + nb231nf_qqH], xmm4
	
	mov rsi, [rbp + nb231nf_type]
	mov eax, [rsi + rax*4]
	mov rsi, [rbp + nb231nf_vdwparam]
	shl eax, 1	
	mov edi, [rsp + nb231nf_ntia]
	add eax, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movapd [rsp + nb231nf_c6], xmm4
	movapd [rsp + nb231nf_c12], xmm6
	
	mov rsi, [rbp + nb231nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [rsp + nb231nf_ixO]
	movapd xmm5, [rsp + nb231nf_iyO]
	movapd xmm6, [rsp + nb231nf_izO]

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
	movapd [rsp + nb231nf_rsqO], xmm7
	
	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [rsp + nb231nf_ixH1]
	movapd xmm5, [rsp + nb231nf_iyH1]
	movapd xmm6, [rsp + nb231nf_izH1]

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
	movapd xmm3, [rsp + nb231nf_ixH2]
	movapd xmm4, [rsp + nb231nf_iyH2]
	movapd xmm5, [rsp + nb231nf_izH2]

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

	mulsd  xmm0, [rsp + nb231nf_krf]	
	mulsd  xmm1, [rsp + nb231nf_krf]	
	mulsd  xmm2, [rsp + nb231nf_krf]	

	movapd [rsp + nb231nf_krsqH2], xmm0
	movapd [rsp + nb231nf_krsqH1], xmm1
	movapd [rsp + nb231nf_krsqO], xmm2

	;# start with rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb231nf_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb231nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb231nf_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb231nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb231nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb231nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb231nf_three]
	subsd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb231nf_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	movsd  [rsp + nb231nf_rinvH1], xmm6
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb231nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb231nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb231nf_three]
	subsd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb231nf_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 
	movsd [rsp + nb231nf_rinvH2], xmm5
	
	;# do O interactions 
	movsd xmm0, xmm7
	movsd xmm1, [rsp + nb231nf_krsqO]
	addsd  xmm0, xmm1
	subsd  xmm0, [rsp + nb231nf_crf] ;# xmm0=rinv+ krsq-crf 
	mulsd  xmm0, [rsp + nb231nf_qqO]

	addsd  xmm0, [rsp + nb231nf_vctot]
	movsd [rsp + nb231nf_vctot], xmm0

	movsd xmm0, xmm7
	movsd xmm4, [rsp + nb231nf_rsqO]
	;# LJ table interaction. xmm0=rinv, xmm4=rsq
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [rsp + nb231nf_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3
	
	mov  rsi, [rbp + nb231nf_VFtab]

	;# dispersion 
	movlpd xmm4, [rsi + rbx*8]	;# Y1 	
	movhpd xmm4, [rsi + rbx*8 + 8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rbx*8 + 16]	;# G1
	movhpd xmm6, [rsi + rbx*8 + 24]	;# G1 H1 	
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

	movsd xmm4, [rsp + nb231nf_c6]
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# Update Vvdwtot directly 
	addsd  xmm5, [rsp + nb231nf_Vvdwtot]
	movsd [rsp + nb231nf_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [rsi + rbx*8 + 32]	;# Y1 	
	movhpd xmm4, [rsi + rbx*8 + 40]	;# Y1 F1 	

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rbx*8 + 48]	;# G1
	movhpd xmm6, [rsi + rbx*8 + 56]	;# G1 H1 	

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
	
	movsd xmm4, [rsp + nb231nf_c12]
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [rsp + nb231nf_Vvdwtot]
	movsd [rsp + nb231nf_Vvdwtot], xmm5

	;# H1 interactions 
	movsd   xmm6, [rsp + nb231nf_rinvH1]	
	movsd   xmm4, xmm6
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movsd   xmm0, [rsp + nb231nf_krsqH1]
	addsd   xmm6, xmm0	;# xmm6=rinv+ krsq 
	subsd   xmm6, [rsp + nb231nf_crf]
	mulsd   xmm6, [rsp + nb231nf_qqH] ;# vcoul 
	
	addsd  xmm6, [rsp + nb231nf_vctot]
	movlpd [rsp + nb231nf_vctot], xmm6

	;# H2 interactions 
	movapd  xmm5, [rsp + nb231nf_rinvH2]	
	movapd  xmm4, xmm5
	mulpd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movapd  xmm7, xmm5
	movsd   xmm0, [rsp + nb231nf_krsqH2]
	addsd   xmm5, xmm0	;# xmm5=rinv+ krsq 
	subsd   xmm5, [rsp + nb231nf_crf]
	mulsd   xmm5, [rsp + nb231nf_qqH] ;# vcoul 
	
	addsd  xmm5, [rsp + nb231nf_vctot]
	movlpd [rsp + nb231nf_vctot], xmm5

.nb231nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb231nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb231nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb231nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb231nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb231nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   rax, [rbp + nb231nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [rsp + nb231nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb231nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb231nf_n], esi
        jmp .nb231nf_outer
.nb231nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb231nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb231nf_end
        ;# non-zero, do one more workunit
        jmp   .nb231nf_threadloop
.nb231nf_end:
	mov eax, [rsp + nb231nf_nouter]
	mov ebx, [rsp + nb231nf_ninner]
	mov rcx, [rbp + nb231nf_outeriter]
	mov rdx, [rbp + nb231nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 544
	emms

    ;# Save xmm registers to stack
    movaps xmm6,  [rsp      ]
    movaps xmm7,  [rsp + 16 ]
    movaps xmm8,  [rsp + 32 ]
    movaps xmm9,  [rsp + 48 ]
    movaps xmm10, [rsp + 64 ]
    movaps xmm11, [rsp + 80 ]
    movaps xmm12, [rsp + 96 ]
    movaps xmm13, [rsp + 112]
    movaps xmm14, [rsp + 128]
    movaps xmm15, [rsp + 144]

    ;# Reset pointers after restoring xmm6-15
    add rsp, 168

    pop r15
    pop r14
    pop r13
    pop r12
    pop rdi
    pop rsi
    pop rbx
    
	pop	rbp
	ret
