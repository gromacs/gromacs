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



.globl nb_kernel104_x86_64_sse2
.globl _nb_kernel104_x86_64_sse2
nb_kernel104_x86_64_sse2:	
_nb_kernel104_x86_64_sse2:	
.equiv          nb104_fshift,           16
.equiv          nb104_gid,              24
.equiv          nb104_pos,              32
.equiv          nb104_faction,          40
.equiv          nb104_charge,           48
.equiv          nb104_p_facel,          56
.equiv          nb104_argkrf,           64
.equiv          nb104_argcrf,           72
.equiv          nb104_Vc,               80
.equiv          nb104_type,             88
.equiv          nb104_p_ntype,          96
.equiv          nb104_vdwparam,         104
.equiv          nb104_Vvdw,             112
.equiv          nb104_p_tabscale,       120
.equiv          nb104_VFtab,            128
.equiv          nb104_invsqrta,         136
.equiv          nb104_dvda,             144
.equiv          nb104_p_gbtabscale,     152
.equiv          nb104_GBtab,            160
.equiv          nb104_p_nthreads,       168
.equiv          nb104_count,            176
.equiv          nb104_mtx,              184
.equiv          nb104_outeriter,        192
.equiv          nb104_inneriter,        200
.equiv          nb104_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 	
.equiv          nb104_ixM,              0
.equiv          nb104_iyM,              16
.equiv          nb104_izM,              32
.equiv          nb104_ixH1,             48
.equiv          nb104_iyH1,             64
.equiv          nb104_izH1,             80
.equiv          nb104_ixH2,             96
.equiv          nb104_iyH2,             112
.equiv          nb104_izH2,             128
.equiv          nb104_jxM,              144
.equiv          nb104_jyM,              160
.equiv          nb104_jzM,              176
.equiv          nb104_jxH1,             192
.equiv          nb104_jyH1,             208
.equiv          nb104_jzH1,             224
.equiv          nb104_jxH2,             240
.equiv          nb104_jyH2,             256
.equiv          nb104_jzH2,             272
.equiv          nb104_dxMM,             288
.equiv          nb104_dyMM,             304
.equiv          nb104_dzMM,             320
.equiv          nb104_dxMH1,            336
.equiv          nb104_dyMH1,            352
.equiv          nb104_dzMH1,            368
.equiv          nb104_dxMH2,            384
.equiv          nb104_dyMH2,            400
.equiv          nb104_dzMH2,            416
.equiv          nb104_dxH1M,            432
.equiv          nb104_dyH1M,            448
.equiv          nb104_dzH1M,            464
.equiv          nb104_dxH1H1,           480
.equiv          nb104_dyH1H1,           496
.equiv          nb104_dzH1H1,           512
.equiv          nb104_dxH1H2,           528
.equiv          nb104_dyH1H2,           544
.equiv          nb104_dzH1H2,           560
.equiv          nb104_dxH2M,            576
.equiv          nb104_dyH2M,            592
.equiv          nb104_dzH2M,            608
.equiv          nb104_dxH2H1,           624
.equiv          nb104_dyH2H1,           640
.equiv          nb104_dzH2H1,           656
.equiv          nb104_dxH2H2,           672
.equiv          nb104_dyH2H2,           688
.equiv          nb104_dzH2H2,           704
.equiv          nb104_qqMM,             720
.equiv          nb104_qqMH,             736
.equiv          nb104_qqHH,             752
.equiv          nb104_vctot,            768
.equiv          nb104_fixM,             784
.equiv          nb104_fiyM,             800
.equiv          nb104_fizM,             816
.equiv          nb104_fixH1,            832
.equiv          nb104_fiyH1,            848
.equiv          nb104_fizH1,            864
.equiv          nb104_fixH2,            880
.equiv          nb104_fiyH2,            896
.equiv          nb104_fizH2,            912
.equiv          nb104_fjxM,             928
.equiv          nb104_fjyM,             944
.equiv          nb104_fjzM,             960
.equiv          nb104_fjxH1,            976
.equiv          nb104_fjyH1,            992
.equiv          nb104_fjzH1,            1008
.equiv          nb104_fjxH2,            1024
.equiv          nb104_fjyH2,            1040
.equiv          nb104_fjzH2,            1056
.equiv          nb104_half,             1072
.equiv          nb104_three,            1088
.equiv          nb104_rsqMM,            1104
.equiv          nb104_rsqMH1,           1120
.equiv          nb104_rsqMH2,           1136
.equiv          nb104_rsqH1M,           1152
.equiv          nb104_rsqH1H1,          1168
.equiv          nb104_rsqH1H2,          1184
.equiv          nb104_rsqH2M,           1200
.equiv          nb104_rsqH2H1,          1216
.equiv          nb104_rsqH2H2,          1232
.equiv          nb104_rinvMM,           1248
.equiv          nb104_rinvMH1,          1264
.equiv          nb104_rinvMH2,          1280
.equiv          nb104_rinvH1M,          1296
.equiv          nb104_rinvH1H1,         1312
.equiv          nb104_rinvH1H2,         1328
.equiv          nb104_rinvH2M,          1344
.equiv          nb104_rinvH2H1,         1360
.equiv          nb104_rinvH2H2,         1376
.equiv          nb104_nri,              1392
.equiv          nb104_innerjjnr,        1400
.equiv          nb104_iinr,             1408
.equiv          nb104_jindex,           1416
.equiv          nb104_jjnr,             1424
.equiv          nb104_shift,            1432
.equiv          nb104_shiftvec,         1440
.equiv          nb104_facel,            1448
.equiv          nb104_is3,              1456
.equiv          nb104_ii3,              1460
.equiv          nb104_innerk,           1464
.equiv          nb104_n,                1468
.equiv          nb104_nn1,              1472
.equiv          nb104_nouter,           1476
.equiv          nb104_ninner,           1480

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
	sub rsp, 1488		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb104_nouter], eax
	mov [rsp + nb104_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb104_nri], edi
	mov [rsp + nb104_iinr], rsi
	mov [rsp + nb104_jindex], rdx
	mov [rsp + nb104_jjnr], rcx
	mov [rsp + nb104_shift], r8
	mov [rsp + nb104_shiftvec], r9
	mov rsi, [rbp + nb104_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb104_facel], xmm0


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb104_half], eax
	mov [rsp + nb104_half+4], ebx
	movsd xmm1, [rsp + nb104_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb104_half], xmm1
	movapd [rsp + nb104_three], xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb104_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb104_charge]
	movsd xmm3, [rdx + rbx*8 + 24]	;# qM 
	movsd xmm4, xmm3		;# qM 
	movsd xmm5, [rdx + rbx*8 + 8]	;# qH 
	mov rsi, [rbp + nb104_p_facel]
	movsd xmm0, [rsi]
	movsd xmm6, [rsp + nb104_facel]	;# facel 
	mulsd  xmm3, xmm3		;# qM*qM 
	mulsd  xmm4, xmm5		;# qM*qH 
	mulsd  xmm5, xmm5		;# qH*qH 
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb104_qqMM], xmm3
	movapd [rsp + nb104_qqMH], xmm4
	movapd [rsp + nb104_qqHH], xmm5

.nb104_threadloop:
        mov   rsi, [rbp + nb104_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb104_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb104_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb104_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb104_n], eax
        mov [rsp + nb104_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb104_outerstart
        jmp .nb104_end

.nb104_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb104_nouter]
	mov [rsp + nb104_nouter], ebx

.nb104_outer:
	mov   rax, [rsp + nb104_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb104_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb104_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb104_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb104_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb104_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [rax + rbx*8 + 24]
	addsd xmm4, [rax + rbx*8 + 32]
	addsd xmm5, [rax + rbx*8 + 40]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb104_ixH1], xmm3
	movapd [rsp + nb104_iyH1], xmm4
	movapd [rsp + nb104_izH1], xmm5

	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm0, [rax + rbx*8 + 48]
	addsd xmm1, [rax + rbx*8 + 56]
	addsd xmm2, [rax + rbx*8 + 64]		
	addsd xmm3, [rax + rbx*8 + 72]
	addsd xmm4, [rax + rbx*8 + 80]
	addsd xmm5, [rax + rbx*8 + 88]		

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb104_ixH2], xmm0
	movapd [rsp + nb104_iyH2], xmm1
	movapd [rsp + nb104_izH2], xmm2
	movapd [rsp + nb104_ixM], xmm3
	movapd [rsp + nb104_iyM], xmm4
	movapd [rsp + nb104_izM], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb104_vctot], xmm4
	movapd [rsp + nb104_fixM], xmm4
	movapd [rsp + nb104_fiyM], xmm4
	movapd [rsp + nb104_fizM], xmm4
	movapd [rsp + nb104_fixH1], xmm4
	movapd [rsp + nb104_fiyH1], xmm4
	movapd [rsp + nb104_fizH1], xmm4
	movapd [rsp + nb104_fixH2], xmm4
	movapd [rsp + nb104_fiyH2], xmm4
	movapd [rsp + nb104_fizH2], xmm4
	
	mov   rax, [rsp + nb104_jindex]
	mov   ecx, [rax+rsi*4]	     ;# jindex[n] 
	mov   edx, [rax +rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb104_pos]
	mov   rdi, [rbp + nb104_faction]	
	mov   rax, [rsp + nb104_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb104_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb104_ninner]
	mov   [rsp + nb104_ninner], ecx
	add   edx, 0
	mov   [rsp + nb104_innerk], edx    ;# number of innerloop atoms 
	jge   .nb104_unroll_loop
	jmp   .nb104_checksingle
.nb104_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb104_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	
	add qword ptr [rsp + nb104_innerjjnr], 8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb104_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	
	;# move j H1 coordinates to local temp variables 
    movlpd xmm0, [rsi + rax*8 + 24] 
    movlpd xmm1, [rsi + rax*8 + 32] 
    movlpd xmm2, [rsi + rax*8 + 40] 
    movhpd xmm0, [rsi + rbx*8 + 24] 
    movhpd xmm1, [rsi + rbx*8 + 32] 
    movhpd xmm2, [rsi + rbx*8 + 40] 

    ;# xmm0 = H1x
    ;# xmm1 = H1y
    ;# xmm2 = H1z
        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb104_ixH1]
    subpd xmm1, [rsp + nb104_iyH1]
    subpd xmm2, [rsp + nb104_izH1]
    subpd xmm3, [rsp + nb104_ixH2]
    subpd xmm4, [rsp + nb104_iyH2]
    subpd xmm5, [rsp + nb104_izH2]
    subpd xmm6, [rsp + nb104_ixM]
    subpd xmm7, [rsp + nb104_iyM]
    subpd xmm8, [rsp + nb104_izM]
    
	movapd [rsp + nb104_dxH1H1], xmm0
	movapd [rsp + nb104_dyH1H1], xmm1
	movapd [rsp + nb104_dzH1H1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb104_dxH2H1], xmm3
	movapd [rsp + nb104_dyH2H1], xmm4
	movapd [rsp + nb104_dzH2H1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb104_dxMH1], xmm6
	movapd [rsp + nb104_dyMH1], xmm7
	movapd [rsp + nb104_dzMH1], xmm8
	mulpd  xmm6, xmm6
	mulpd  xmm7, xmm7
	mulpd  xmm8, xmm8
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    addpd  xmm6, xmm7
    addpd  xmm6, xmm8

	;# start doing invsqrt for jH1 atoms
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
		
	movapd  xmm9, [rsp + nb104_three]
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

	movapd  xmm15, [rsp + nb104_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvH1H1 
	mulpd   xmm10, xmm15 ;# first iteration for rinvH2H1
    mulpd   xmm11, xmm15 ;# first iteration for rinvMH1	

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb104_three]
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

	movapd  xmm15, [rsp + nb104_half]
	mulpd   xmm9, xmm15  ;#  rinvH1H1 
	mulpd   xmm10, xmm15 ;#   rinvH2H1
    mulpd   xmm11, xmm15 ;#   rinvMH1
	
	;# H1 interactions 
    movapd xmm0, xmm9
    movapd xmm1, xmm10
    movapd xmm2, xmm11
    mulpd  xmm9, xmm9
    mulpd  xmm10, xmm10
    mulpd  xmm11, xmm11
    mulpd  xmm0, [rsp + nb104_qqHH] 
    mulpd  xmm1, [rsp + nb104_qqHH] 
    mulpd  xmm2, [rsp + nb104_qqMH] 
    mulpd  xmm9, xmm0
    mulpd  xmm10, xmm1
    mulpd  xmm11, xmm2
    
    addpd xmm0, [rsp + nb104_vctot] 
    addpd xmm1, xmm2
    addpd xmm0, xmm1
    movapd [rsp + nb104_vctot], xmm0
    
    ;# move j H1 forces to xmm0-xmm2
	mov   rdi, [rbp + nb104_faction]	
	movlpd xmm0, [rdi + rax*8 + 24]
	movlpd xmm1, [rdi + rax*8 + 32]
	movlpd xmm2, [rdi + rax*8 + 40]
	movhpd xmm0, [rdi + rbx*8 + 24]
	movhpd xmm1, [rdi + rbx*8 + 32]
	movhpd xmm2, [rdi + rbx*8 + 40]

    movapd xmm7, xmm9
    movapd xmm8, xmm9
    movapd xmm13, xmm11
    movapd xmm14, xmm11
    movapd xmm15, xmm11
    movapd xmm11, xmm10
    movapd xmm12, xmm10

	mulpd xmm7, [rsp + nb104_dxH1H1]
	mulpd xmm8, [rsp + nb104_dyH1H1]
	mulpd xmm9, [rsp + nb104_dzH1H1]
	mulpd xmm10, [rsp + nb104_dxH2H1]
	mulpd xmm11, [rsp + nb104_dyH2H1]
	mulpd xmm12, [rsp + nb104_dzH2H1]
	mulpd xmm13, [rsp + nb104_dxMH1]
	mulpd xmm14, [rsp + nb104_dyMH1]
	mulpd xmm15, [rsp + nb104_dzMH1]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb104_fixH1]
    addpd xmm8, [rsp + nb104_fiyH1]
    addpd xmm9, [rsp + nb104_fizH1]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb104_fixH2]
    addpd xmm11, [rsp + nb104_fiyH2]
    addpd xmm12, [rsp + nb104_fizH2]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb104_fixM]
    addpd xmm14, [rsp + nb104_fiyM]
    addpd xmm15, [rsp + nb104_fizM]

    movapd [rsp + nb104_fixH1], xmm7
    movapd [rsp + nb104_fiyH1], xmm8
    movapd [rsp + nb104_fizH1], xmm9
    movapd [rsp + nb104_fixH2], xmm10
    movapd [rsp + nb104_fiyH2], xmm11
    movapd [rsp + nb104_fizH2], xmm12
    movapd [rsp + nb104_fixM], xmm13
    movapd [rsp + nb104_fiyM], xmm14
    movapd [rsp + nb104_fizM], xmm15
   
    ;# store back j H1 forces from xmm0-xmm2
	movlpd [rdi + rax*8 + 24], xmm0
	movlpd [rdi + rax*8 + 32], xmm1
	movlpd [rdi + rax*8 + 40], xmm2
	movhpd [rdi + rbx*8 + 24], xmm0
	movhpd [rdi + rbx*8 + 32], xmm1
	movhpd [rdi + rbx*8 + 40], xmm2

	;# move j H2 coordinates to local temp variables 
	mov   rsi, [rbp + nb104_pos]	
    movlpd xmm0, [rsi + rax*8 + 48] 
    movlpd xmm1, [rsi + rax*8 + 56] 
    movlpd xmm2, [rsi + rax*8 + 64] 
    movhpd xmm0, [rsi + rbx*8 + 48] 
    movhpd xmm1, [rsi + rbx*8 + 56] 
    movhpd xmm2, [rsi + rbx*8 + 64] 

    ;# xmm0 = H2x
    ;# xmm1 = H2y
    ;# xmm2 = H2z
        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb104_ixH1]
    subpd xmm1, [rsp + nb104_iyH1]
    subpd xmm2, [rsp + nb104_izH1]
    subpd xmm3, [rsp + nb104_ixH2]
    subpd xmm4, [rsp + nb104_iyH2]
    subpd xmm5, [rsp + nb104_izH2]
    subpd xmm6, [rsp + nb104_ixM]
    subpd xmm7, [rsp + nb104_iyM]
    subpd xmm8, [rsp + nb104_izM]
    
	movapd [rsp + nb104_dxH1H2], xmm0
	movapd [rsp + nb104_dyH1H2], xmm1
	movapd [rsp + nb104_dzH1H2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb104_dxH2H2], xmm3
	movapd [rsp + nb104_dyH2H2], xmm4
	movapd [rsp + nb104_dzH2H2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb104_dxMH2], xmm6
	movapd [rsp + nb104_dyMH2], xmm7
	movapd [rsp + nb104_dzMH2], xmm8
	mulpd  xmm6, xmm6
	mulpd  xmm7, xmm7
	mulpd  xmm8, xmm8
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    addpd  xmm6, xmm7
    addpd  xmm6, xmm8

	;# start doing invsqrt for jH2 atoms
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
		
	movapd  xmm9, [rsp + nb104_three]
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

	movapd  xmm15, [rsp + nb104_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvH1H2 
	mulpd   xmm10, xmm15 ;# first iteration for rinvH2H2
    mulpd   xmm11, xmm15 ;# first iteration for rinvMH1H2

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb104_three]
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

	movapd  xmm15, [rsp + nb104_half]
	mulpd   xmm9, xmm15  ;#  rinvH1H2
	mulpd   xmm10, xmm15 ;#   rinvH2H2
    mulpd   xmm11, xmm15 ;#   rinvMH2
	
	;# H2 interactions 
    movapd xmm0, xmm9
    movapd xmm1, xmm10
    movapd xmm2, xmm11
    mulpd  xmm9, xmm9
    mulpd  xmm10, xmm10
    mulpd  xmm11, xmm11
    mulpd  xmm0, [rsp + nb104_qqHH] 
    mulpd  xmm1, [rsp + nb104_qqHH] 
    mulpd  xmm2, [rsp + nb104_qqMH] 
    mulpd  xmm9, xmm0
    mulpd  xmm10, xmm1
    mulpd  xmm11, xmm2
    
    addpd xmm0, [rsp + nb104_vctot] 
    addpd xmm1, xmm2
    addpd xmm0, xmm1
    movapd [rsp + nb104_vctot], xmm0
    
    ;# move j H2 forces to xmm0-xmm2
	mov   rdi, [rbp + nb104_faction]	
	movlpd xmm0, [rdi + rax*8 + 48]
	movlpd xmm1, [rdi + rax*8 + 56]
	movlpd xmm2, [rdi + rax*8 + 64]
	movhpd xmm0, [rdi + rbx*8 + 48]
	movhpd xmm1, [rdi + rbx*8 + 56]
	movhpd xmm2, [rdi + rbx*8 + 64]

    movapd xmm7, xmm9
    movapd xmm8, xmm9
    movapd xmm13, xmm11
    movapd xmm14, xmm11
    movapd xmm15, xmm11
    movapd xmm11, xmm10
    movapd xmm12, xmm10

	mulpd xmm7, [rsp + nb104_dxH1H2]
	mulpd xmm8, [rsp + nb104_dyH1H2]
	mulpd xmm9, [rsp + nb104_dzH1H2]
	mulpd xmm10, [rsp + nb104_dxH2H2]
	mulpd xmm11, [rsp + nb104_dyH2H2]
	mulpd xmm12, [rsp + nb104_dzH2H2]
	mulpd xmm13, [rsp + nb104_dxMH2]
	mulpd xmm14, [rsp + nb104_dyMH2]
	mulpd xmm15, [rsp + nb104_dzMH2]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb104_fixH1]
    addpd xmm8, [rsp + nb104_fiyH1]
    addpd xmm9, [rsp + nb104_fizH1]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb104_fixH2]
    addpd xmm11, [rsp + nb104_fiyH2]
    addpd xmm12, [rsp + nb104_fizH2]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb104_fixM]
    addpd xmm14, [rsp + nb104_fiyM]
    addpd xmm15, [rsp + nb104_fizM]

    movapd [rsp + nb104_fixH1], xmm7
    movapd [rsp + nb104_fiyH1], xmm8
    movapd [rsp + nb104_fizH1], xmm9
    movapd [rsp + nb104_fixH2], xmm10
    movapd [rsp + nb104_fiyH2], xmm11
    movapd [rsp + nb104_fizH2], xmm12
    movapd [rsp + nb104_fixM], xmm13
    movapd [rsp + nb104_fiyM], xmm14
    movapd [rsp + nb104_fizM], xmm15
   
    ;# store back j H2 forces from xmm0-xmm2
	movlpd [rdi + rax*8 + 48], xmm0
	movlpd [rdi + rax*8 + 56], xmm1
	movlpd [rdi + rax*8 + 64], xmm2
	movhpd [rdi + rbx*8 + 48], xmm0
	movhpd [rdi + rbx*8 + 56], xmm1
	movhpd [rdi + rbx*8 + 64], xmm2
       
	;# move j M coordinates to local temp variables 
	mov   rsi, [rbp + nb104_pos]	
    movlpd xmm0, [rsi + rax*8 + 72] 
    movlpd xmm1, [rsi + rax*8 + 80] 
    movlpd xmm2, [rsi + rax*8 + 88] 
    movhpd xmm0, [rsi + rbx*8 + 72] 
    movhpd xmm1, [rsi + rbx*8 + 80] 
    movhpd xmm2, [rsi + rbx*8 + 88] 

    ;# xmm0 = Mx
    ;# xmm1 = My
    ;# xmm2 = Mz
        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb104_ixH1]
    subpd xmm1, [rsp + nb104_iyH1]
    subpd xmm2, [rsp + nb104_izH1]
    subpd xmm3, [rsp + nb104_ixH2]
    subpd xmm4, [rsp + nb104_iyH2]
    subpd xmm5, [rsp + nb104_izH2]
    subpd xmm6, [rsp + nb104_ixM]
    subpd xmm7, [rsp + nb104_iyM]
    subpd xmm8, [rsp + nb104_izM]
    
	movapd [rsp + nb104_dxH1M], xmm0
	movapd [rsp + nb104_dyH1M], xmm1
	movapd [rsp + nb104_dzH1M], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb104_dxH2M], xmm3
	movapd [rsp + nb104_dyH2M], xmm4
	movapd [rsp + nb104_dzH2M], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb104_dxMM], xmm6
	movapd [rsp + nb104_dyMM], xmm7
	movapd [rsp + nb104_dzMM], xmm8
	mulpd  xmm6, xmm6
	mulpd  xmm7, xmm7
	mulpd  xmm8, xmm8
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    addpd  xmm6, xmm7
    addpd  xmm6, xmm8

	;# start doing invsqrt for jM atoms
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
		
	movapd  xmm9, [rsp + nb104_three]
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

	movapd  xmm15, [rsp + nb104_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvH1M 
	mulpd   xmm10, xmm15 ;# first iteration for rinvH2M
    mulpd   xmm11, xmm15 ;# first iteration for rinvMM

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb104_three]
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

	movapd  xmm15, [rsp + nb104_half]
	mulpd   xmm9, xmm15  ;#  rinvH1M
	mulpd   xmm10, xmm15 ;#   rinvH2M
    mulpd   xmm11, xmm15 ;#   rinvMM
	
	;# M interactions 
    movapd xmm0, xmm9
    movapd xmm1, xmm10
    movapd xmm2, xmm11
    mulpd  xmm9, xmm9
    mulpd  xmm10, xmm10
    mulpd  xmm11, xmm11
    mulpd  xmm0, [rsp + nb104_qqMH] 
    mulpd  xmm1, [rsp + nb104_qqMH] 
    mulpd  xmm2, [rsp + nb104_qqMM] 
    mulpd  xmm9, xmm0
    mulpd  xmm10, xmm1
    mulpd  xmm11, xmm2
    
    addpd xmm0, [rsp + nb104_vctot] 
    addpd xmm1, xmm2
    addpd xmm0, xmm1
    movapd [rsp + nb104_vctot], xmm0
    
    ;# move j M forces to xmm0-xmm2
	mov   rdi, [rbp + nb104_faction]	
	movlpd xmm0, [rdi + rax*8 + 72]
	movlpd xmm1, [rdi + rax*8 + 80]
	movlpd xmm2, [rdi + rax*8 + 88]
	movhpd xmm0, [rdi + rbx*8 + 72]
	movhpd xmm1, [rdi + rbx*8 + 80]
	movhpd xmm2, [rdi + rbx*8 + 88]

    movapd xmm7, xmm9
    movapd xmm8, xmm9
    movapd xmm13, xmm11
    movapd xmm14, xmm11
    movapd xmm15, xmm11
    movapd xmm11, xmm10
    movapd xmm12, xmm10

	mulpd xmm7, [rsp + nb104_dxH1M]
	mulpd xmm8, [rsp + nb104_dyH1M]
	mulpd xmm9, [rsp + nb104_dzH1M]
	mulpd xmm10, [rsp + nb104_dxH2M]
	mulpd xmm11, [rsp + nb104_dyH2M]
	mulpd xmm12, [rsp + nb104_dzH2M]
	mulpd xmm13, [rsp + nb104_dxMM]
	mulpd xmm14, [rsp + nb104_dyMM]
	mulpd xmm15, [rsp + nb104_dzMM]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb104_fixH1]
    addpd xmm8, [rsp + nb104_fiyH1]
    addpd xmm9, [rsp + nb104_fizH1]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb104_fixH2]
    addpd xmm11, [rsp + nb104_fiyH2]
    addpd xmm12, [rsp + nb104_fizH2]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb104_fixM]
    addpd xmm14, [rsp + nb104_fiyM]
    addpd xmm15, [rsp + nb104_fizM]

    movapd [rsp + nb104_fixH1], xmm7
    movapd [rsp + nb104_fiyH1], xmm8
    movapd [rsp + nb104_fizH1], xmm9
    movapd [rsp + nb104_fixH2], xmm10
    movapd [rsp + nb104_fiyH2], xmm11
    movapd [rsp + nb104_fizH2], xmm12
    movapd [rsp + nb104_fixM], xmm13
    movapd [rsp + nb104_fiyM], xmm14
    movapd [rsp + nb104_fizM], xmm15
   
    ;# store back j M forces from xmm0-xmm2
	movlpd [rdi + rax*8 + 72], xmm0
	movlpd [rdi + rax*8 + 80], xmm1
	movlpd [rdi + rax*8 + 88], xmm2
	movhpd [rdi + rbx*8 + 72], xmm0
	movhpd [rdi + rbx*8 + 80], xmm1
	movhpd [rdi + rbx*8 + 88], xmm2
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb104_innerk],  2
	jl    .nb104_checksingle
	jmp   .nb104_unroll_loop
.nb104_checksingle:
	mov   edx, [rsp + nb104_innerk]
	and   edx, 1
	jnz   .nb104_dosingle
	jmp   .nb104_updateouterdata
.nb104_dosingle:
	mov   rdx, [rsp + nb104_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	

	mov rsi, [rbp + nb104_pos]
	lea   rax, [rax + rax*2]  

	
	;# move j H1 coordinates to local temp variables 
    movsd xmm0, [rsi + rax*8 + 24] 
    movsd xmm1, [rsi + rax*8 + 32] 
    movsd xmm2, [rsi + rax*8 + 40] 

    ;# xmm0 = H1x
    ;# xmm1 = H1y
    ;# xmm2 = H1z
        
    movsd xmm3, xmm0
    movsd xmm4, xmm1
    movsd xmm5, xmm2
    movsd xmm6, xmm0
    movsd xmm7, xmm1
    movsd xmm8, xmm2
    
    subsd xmm0, [rsp + nb104_ixH1]
    subsd xmm1, [rsp + nb104_iyH1]
    subsd xmm2, [rsp + nb104_izH1]
    subsd xmm3, [rsp + nb104_ixH2]
    subsd xmm4, [rsp + nb104_iyH2]
    subsd xmm5, [rsp + nb104_izH2]
    subsd xmm6, [rsp + nb104_ixM]
    subsd xmm7, [rsp + nb104_iyM]
    subsd xmm8, [rsp + nb104_izM]
    
	movsd [rsp + nb104_dxH1H1], xmm0
	movsd [rsp + nb104_dyH1H1], xmm1
	movsd [rsp + nb104_dzH1H1], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movsd [rsp + nb104_dxH2H1], xmm3
	movsd [rsp + nb104_dyH2H1], xmm4
	movsd [rsp + nb104_dzH2H1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movsd [rsp + nb104_dxMH1], xmm6
	movsd [rsp + nb104_dyMH1], xmm7
	movsd [rsp + nb104_dzMH1], xmm8
	mulsd  xmm6, xmm6
	mulsd  xmm7, xmm7
	mulsd  xmm8, xmm8
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
    addsd  xmm6, xmm7
    addsd  xmm6, xmm8

	;# start doing invsqrt for jH1 atoms
    cvtsd2ss xmm1, xmm0
    cvtsd2ss xmm4, xmm3
    cvtsd2ss xmm7, xmm6
	rsqrtss xmm1, xmm1
	rsqrtss xmm4, xmm4
    rsqrtss xmm7, xmm7
    cvtss2sd xmm1, xmm1
    cvtss2sd xmm4, xmm4
    cvtss2sd xmm7, xmm7
	
	movsd  xmm2, xmm1
	movsd  xmm5, xmm4
    movsd  xmm8, xmm7
    
	mulsd   xmm1, xmm1 ;# lu*lu
	mulsd   xmm4, xmm4 ;# lu*lu
    mulsd   xmm7, xmm7 ;# lu*lu
		
	movsd  xmm9, [rsp + nb104_three]
	movsd  xmm10, xmm9
    movsd  xmm11, xmm9

	mulsd   xmm1, xmm0 ;# rsq*lu*lu
	mulsd   xmm4, xmm3 ;# rsq*lu*lu 
    mulsd   xmm7, xmm6 ;# rsq*lu*lu
	
	subsd   xmm9, xmm1
	subsd   xmm10, xmm4
    subsd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm2
	mulsd   xmm10, xmm5
    mulsd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movsd  xmm15, [rsp + nb104_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvH1H1 
	mulsd   xmm10, xmm15 ;# first iteration for rinvH2H1
    mulsd   xmm11, xmm15 ;# first iteration for rinvMH1	

    ;# second iteration step    
	movsd  xmm2, xmm9
	movsd  xmm5, xmm10
    movsd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movsd  xmm1, [rsp + nb104_three]
	movsd  xmm4, xmm1
    movsd  xmm7, xmm1

	mulsd   xmm2, xmm0 ;# rsq*lu*lu
	mulsd   xmm5, xmm3 ;# rsq*lu*lu 
    mulsd   xmm8, xmm6 ;# rsq*lu*lu
	
	subsd   xmm1, xmm2
	subsd   xmm4, xmm5
    subsd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm1
	mulsd   xmm10, xmm4
    mulsd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movsd  xmm15, [rsp + nb104_half]
	mulsd   xmm9, xmm15  ;#  rinvH1H1 
	mulsd   xmm10, xmm15 ;#   rinvH2H1
    mulsd   xmm11, xmm15 ;#   rinvMH1
	
	;# H1 interactions 
    movsd xmm0, xmm9
    movsd xmm1, xmm10
    movsd xmm2, xmm11
    mulsd  xmm9, xmm9
    mulsd  xmm10, xmm10
    mulsd  xmm11, xmm11
    mulsd  xmm0, [rsp + nb104_qqHH] 
    mulsd  xmm1, [rsp + nb104_qqHH] 
    mulsd  xmm2, [rsp + nb104_qqMH] 
    mulsd  xmm9, xmm0
    mulsd  xmm10, xmm1
    mulsd  xmm11, xmm2
    
    addsd xmm0, [rsp + nb104_vctot] 
    addsd xmm1, xmm2
    addsd xmm0, xmm1
    movsd [rsp + nb104_vctot], xmm0
    
    ;# move j H1 forces to xmm0-xmm2
	movsd xmm0, [rdi + rax*8 + 24]
	movsd xmm1, [rdi + rax*8 + 32]
	movsd xmm2, [rdi + rax*8 + 40]

    movsd xmm7, xmm9
    movsd xmm8, xmm9
    movsd xmm13, xmm11
    movsd xmm14, xmm11
    movsd xmm15, xmm11
    movsd xmm11, xmm10
    movsd xmm12, xmm10

	mulsd xmm7, [rsp + nb104_dxH1H1]
	mulsd xmm8, [rsp + nb104_dyH1H1]
	mulsd xmm9, [rsp + nb104_dzH1H1]
	mulsd xmm10, [rsp + nb104_dxH2H1]
	mulsd xmm11, [rsp + nb104_dyH2H1]
	mulsd xmm12, [rsp + nb104_dzH2H1]
	mulsd xmm13, [rsp + nb104_dxMH1]
	mulsd xmm14, [rsp + nb104_dyMH1]
	mulsd xmm15, [rsp + nb104_dzMH1]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb104_fixH1]
    addsd xmm8, [rsp + nb104_fiyH1]
    addsd xmm9, [rsp + nb104_fizH1]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb104_fixH2]
    addsd xmm11, [rsp + nb104_fiyH2]
    addsd xmm12, [rsp + nb104_fizH2]

    addsd xmm0, xmm13
    addsd xmm1, xmm14
    addsd xmm2, xmm15
    addsd xmm13, [rsp + nb104_fixM]
    addsd xmm14, [rsp + nb104_fiyM]
    addsd xmm15, [rsp + nb104_fizM]

    movsd [rsp + nb104_fixH1], xmm7
    movsd [rsp + nb104_fiyH1], xmm8
    movsd [rsp + nb104_fizH1], xmm9
    movsd [rsp + nb104_fixH2], xmm10
    movsd [rsp + nb104_fiyH2], xmm11
    movsd [rsp + nb104_fizH2], xmm12
    movsd [rsp + nb104_fixM], xmm13
    movsd [rsp + nb104_fiyM], xmm14
    movsd [rsp + nb104_fizM], xmm15
   
    ;# store back j H1 forces from xmm0-xmm2
	movsd [rdi + rax*8 + 24], xmm0
	movsd [rdi + rax*8 + 32], xmm1
	movsd [rdi + rax*8 + 40], xmm2

	;# move j H2 coordinates to local temp variables 
    movsd xmm0, [rsi + rax*8 + 48] 
    movsd xmm1, [rsi + rax*8 + 56] 
    movsd xmm2, [rsi + rax*8 + 64] 

    ;# xmm0 = H2x
    ;# xmm1 = H2y
    ;# xmm2 = H2z
        
    movsd xmm3, xmm0
    movsd xmm4, xmm1
    movsd xmm5, xmm2
    movsd xmm6, xmm0
    movsd xmm7, xmm1
    movsd xmm8, xmm2
    
    subsd xmm0, [rsp + nb104_ixH1]
    subsd xmm1, [rsp + nb104_iyH1]
    subsd xmm2, [rsp + nb104_izH1]
    subsd xmm3, [rsp + nb104_ixH2]
    subsd xmm4, [rsp + nb104_iyH2]
    subsd xmm5, [rsp + nb104_izH2]
    subsd xmm6, [rsp + nb104_ixM]
    subsd xmm7, [rsp + nb104_iyM]
    subsd xmm8, [rsp + nb104_izM]
    
	movsd [rsp + nb104_dxH1H2], xmm0
	movsd [rsp + nb104_dyH1H2], xmm1
	movsd [rsp + nb104_dzH1H2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movsd [rsp + nb104_dxH2H2], xmm3
	movsd [rsp + nb104_dyH2H2], xmm4
	movsd [rsp + nb104_dzH2H2], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movsd [rsp + nb104_dxMH2], xmm6
	movsd [rsp + nb104_dyMH2], xmm7
	movsd [rsp + nb104_dzMH2], xmm8
	mulsd  xmm6, xmm6
	mulsd  xmm7, xmm7
	mulsd  xmm8, xmm8
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
    addsd  xmm6, xmm7
    addsd  xmm6, xmm8

	;# start doing invsqrt for jH2 atoms
    cvtsd2ss xmm1, xmm0
    cvtsd2ss xmm4, xmm3
    cvtsd2ss xmm7, xmm6
	rsqrtss xmm1, xmm1
	rsqrtss xmm4, xmm4
    rsqrtss xmm7, xmm7
    cvtss2sd xmm1, xmm1
    cvtss2sd xmm4, xmm4
    cvtss2sd xmm7, xmm7
	
	movsd  xmm2, xmm1
	movsd  xmm5, xmm4
    movsd  xmm8, xmm7
    
	mulsd   xmm1, xmm1 ;# lu*lu
	mulsd   xmm4, xmm4 ;# lu*lu
    mulsd   xmm7, xmm7 ;# lu*lu
		
	movsd  xmm9, [rsp + nb104_three]
	movsd  xmm10, xmm9
    movsd  xmm11, xmm9

	mulsd   xmm1, xmm0 ;# rsq*lu*lu
	mulsd   xmm4, xmm3 ;# rsq*lu*lu 
    mulsd   xmm7, xmm6 ;# rsq*lu*lu
	
	subsd   xmm9, xmm1
	subsd   xmm10, xmm4
    subsd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm2
	mulsd   xmm10, xmm5
    mulsd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movsd  xmm15, [rsp + nb104_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvH1H2 
	mulsd   xmm10, xmm15 ;# first iteration for rinvH2H2
    mulsd   xmm11, xmm15 ;# first iteration for rinvMH2

    ;# second iteration step    
	movsd  xmm2, xmm9
	movsd  xmm5, xmm10
    movsd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movsd  xmm1, [rsp + nb104_three]
	movsd  xmm4, xmm1
    movsd  xmm7, xmm1

	mulsd   xmm2, xmm0 ;# rsq*lu*lu
	mulsd   xmm5, xmm3 ;# rsq*lu*lu 
    mulsd   xmm8, xmm6 ;# rsq*lu*lu
	
	subsd   xmm1, xmm2
	subsd   xmm4, xmm5
    subsd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm1
	mulsd   xmm10, xmm4
    mulsd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movsd  xmm15, [rsp + nb104_half]
	mulsd   xmm9, xmm15  ;#  rinvH1H2
	mulsd   xmm10, xmm15 ;#   rinvH2H2
    mulsd   xmm11, xmm15 ;#   rinvMH2
	
	;# H2 interactions 
    movsd xmm0, xmm9
    movsd xmm1, xmm10
    movsd xmm2, xmm11
    mulsd  xmm9, xmm9
    mulsd  xmm10, xmm10
    mulsd  xmm11, xmm11
    mulsd  xmm0, [rsp + nb104_qqHH] 
    mulsd  xmm1, [rsp + nb104_qqHH] 
    mulsd  xmm2, [rsp + nb104_qqMH] 
    mulsd  xmm9, xmm0
    mulsd  xmm10, xmm1
    mulsd  xmm11, xmm2
    
    addsd xmm0, [rsp + nb104_vctot] 
    addsd xmm1, xmm2
    addsd xmm0, xmm1
    movsd [rsp + nb104_vctot], xmm0
    
    ;# move j H2 forces to xmm0-xmm2
	movsd xmm0, [rdi + rax*8 + 48]
	movsd xmm1, [rdi + rax*8 + 56]
	movsd xmm2, [rdi + rax*8 + 64]

    movsd xmm7, xmm9
    movsd xmm8, xmm9
    movsd xmm13, xmm11
    movsd xmm14, xmm11
    movsd xmm15, xmm11
    movsd xmm11, xmm10
    movsd xmm12, xmm10

	mulsd xmm7, [rsp + nb104_dxH1H2]
	mulsd xmm8, [rsp + nb104_dyH1H2]
	mulsd xmm9, [rsp + nb104_dzH1H2]
	mulsd xmm10, [rsp + nb104_dxH2H2]
	mulsd xmm11, [rsp + nb104_dyH2H2]
	mulsd xmm12, [rsp + nb104_dzH2H2]
	mulsd xmm13, [rsp + nb104_dxMH2]
	mulsd xmm14, [rsp + nb104_dyMH2]
	mulsd xmm15, [rsp + nb104_dzMH2]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb104_fixH1]
    addsd xmm8, [rsp + nb104_fiyH1]
    addsd xmm9, [rsp + nb104_fizH1]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb104_fixH2]
    addsd xmm11, [rsp + nb104_fiyH2]
    addsd xmm12, [rsp + nb104_fizH2]

    addsd xmm0, xmm13
    addsd xmm1, xmm14
    addsd xmm2, xmm15
    addsd xmm13, [rsp + nb104_fixM]
    addsd xmm14, [rsp + nb104_fiyM]
    addsd xmm15, [rsp + nb104_fizM]

    movsd [rsp + nb104_fixH1], xmm7
    movsd [rsp + nb104_fiyH1], xmm8
    movsd [rsp + nb104_fizH1], xmm9
    movsd [rsp + nb104_fixH2], xmm10
    movsd [rsp + nb104_fiyH2], xmm11
    movsd [rsp + nb104_fizH2], xmm12
    movsd [rsp + nb104_fixM], xmm13
    movsd [rsp + nb104_fiyM], xmm14
    movsd [rsp + nb104_fizM], xmm15
   
    ;# store back j H2 forces from xmm0-xmm2
	movsd [rdi + rax*8 + 48], xmm0
	movsd [rdi + rax*8 + 56], xmm1
	movsd [rdi + rax*8 + 64], xmm2
       
	;# move j M coordinates to local temp variables 
    movsd xmm0, [rsi + rax*8 + 72] 
    movsd xmm1, [rsi + rax*8 + 80] 
    movsd xmm2, [rsi + rax*8 + 88] 

    ;# xmm0 = Mx
    ;# xmm1 = My
    ;# xmm2 = Mz
        
    movsd xmm3, xmm0
    movsd xmm4, xmm1
    movsd xmm5, xmm2
    movsd xmm6, xmm0
    movsd xmm7, xmm1
    movsd xmm8, xmm2
    
    subsd xmm0, [rsp + nb104_ixH1]
    subsd xmm1, [rsp + nb104_iyH1]
    subsd xmm2, [rsp + nb104_izH1]
    subsd xmm3, [rsp + nb104_ixH2]
    subsd xmm4, [rsp + nb104_iyH2]
    subsd xmm5, [rsp + nb104_izH2]
    subsd xmm6, [rsp + nb104_ixM]
    subsd xmm7, [rsp + nb104_iyM]
    subsd xmm8, [rsp + nb104_izM]
    
	movsd [rsp + nb104_dxH1M], xmm0
	movsd [rsp + nb104_dyH1M], xmm1
	movsd [rsp + nb104_dzH1M], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movsd [rsp + nb104_dxH2M], xmm3
	movsd [rsp + nb104_dyH2M], xmm4
	movsd [rsp + nb104_dzH2M], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movsd [rsp + nb104_dxMM], xmm6
	movsd [rsp + nb104_dyMM], xmm7
	movsd [rsp + nb104_dzMM], xmm8
	mulsd  xmm6, xmm6
	mulsd  xmm7, xmm7
	mulsd  xmm8, xmm8
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
    addsd  xmm6, xmm7
    addsd  xmm6, xmm8

	;# start doing invsqrt for jM atoms
    cvtsd2ss xmm1, xmm0
    cvtsd2ss xmm4, xmm3
    cvtsd2ss xmm7, xmm6
	rsqrtss xmm1, xmm1
	rsqrtss xmm4, xmm4
    rsqrtss xmm7, xmm7
    cvtss2sd xmm1, xmm1
    cvtss2sd xmm4, xmm4
    cvtss2sd xmm7, xmm7
	
	movsd  xmm2, xmm1
	movsd  xmm5, xmm4
    movsd  xmm8, xmm7
    
	mulsd   xmm1, xmm1 ;# lu*lu
	mulsd   xmm4, xmm4 ;# lu*lu
    mulsd   xmm7, xmm7 ;# lu*lu
		
	movsd  xmm9, [rsp + nb104_three]
	movsd  xmm10, xmm9
    movsd  xmm11, xmm9

	mulsd   xmm1, xmm0 ;# rsq*lu*lu
	mulsd   xmm4, xmm3 ;# rsq*lu*lu 
    mulsd   xmm7, xmm6 ;# rsq*lu*lu
	
	subsd   xmm9, xmm1
	subsd   xmm10, xmm4
    subsd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm2
	mulsd   xmm10, xmm5
    mulsd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movsd  xmm15, [rsp + nb104_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvH1M 
	mulsd   xmm10, xmm15 ;# first iteration for rinvH2M
    mulsd   xmm11, xmm15 ;# first iteration for rinvMM

    ;# second iteration step    
	movsd  xmm2, xmm9
	movsd  xmm5, xmm10
    movsd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movsd  xmm1, [rsp + nb104_three]
	movsd  xmm4, xmm1
    movsd  xmm7, xmm1

	mulsd   xmm2, xmm0 ;# rsq*lu*lu
	mulsd   xmm5, xmm3 ;# rsq*lu*lu 
    mulsd   xmm8, xmm6 ;# rsq*lu*lu
	
	subsd   xmm1, xmm2
	subsd   xmm4, xmm5
    subsd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulsd   xmm9, xmm1
	mulsd   xmm10, xmm4
    mulsd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movsd  xmm15, [rsp + nb104_half]
	mulsd   xmm9, xmm15  ;#  rinvH1M
	mulsd   xmm10, xmm15 ;#   rinvH2M
    mulsd   xmm11, xmm15 ;#   rinvMM
	
	;# M interactions 
    movsd xmm0, xmm9
    movsd xmm1, xmm10
    movsd xmm2, xmm11
    mulsd  xmm9, xmm9
    mulsd  xmm10, xmm10
    mulsd  xmm11, xmm11
    mulsd  xmm0, [rsp + nb104_qqMH] 
    mulsd  xmm1, [rsp + nb104_qqMH] 
    mulsd  xmm2, [rsp + nb104_qqMM] 
    mulsd  xmm9, xmm0
    mulsd  xmm10, xmm1
    mulsd  xmm11, xmm2
    
    addsd xmm0, [rsp + nb104_vctot] 
    addsd xmm1, xmm2
    addsd xmm0, xmm1
    movsd [rsp + nb104_vctot], xmm0
    
    ;# move j M forces to xmm0-xmm2
	movsd xmm0, [rdi + rax*8 + 72]
	movsd xmm1, [rdi + rax*8 + 80]
	movsd xmm2, [rdi + rax*8 + 88]

    movsd xmm7, xmm9
    movsd xmm8, xmm9
    movsd xmm13, xmm11
    movsd xmm14, xmm11
    movsd xmm15, xmm11
    movsd xmm11, xmm10
    movsd xmm12, xmm10

	mulsd xmm7, [rsp + nb104_dxH1M]
	mulsd xmm8, [rsp + nb104_dyH1M]
	mulsd xmm9, [rsp + nb104_dzH1M]
	mulsd xmm10, [rsp + nb104_dxH2M]
	mulsd xmm11, [rsp + nb104_dyH2M]
	mulsd xmm12, [rsp + nb104_dzH2M]
	mulsd xmm13, [rsp + nb104_dxMM]
	mulsd xmm14, [rsp + nb104_dyMM]
	mulsd xmm15, [rsp + nb104_dzMM]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb104_fixH1]
    addsd xmm8, [rsp + nb104_fiyH1]
    addsd xmm9, [rsp + nb104_fizH1]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb104_fixH2]
    addsd xmm11, [rsp + nb104_fiyH2]
    addsd xmm12, [rsp + nb104_fizH2]

    addsd xmm0, xmm13
    addsd xmm1, xmm14
    addsd xmm2, xmm15
    addsd xmm13, [rsp + nb104_fixM]
    addsd xmm14, [rsp + nb104_fiyM]
    addsd xmm15, [rsp + nb104_fizM]

    movsd [rsp + nb104_fixH1], xmm7
    movsd [rsp + nb104_fiyH1], xmm8
    movsd [rsp + nb104_fizH1], xmm9
    movsd [rsp + nb104_fixH2], xmm10
    movsd [rsp + nb104_fiyH2], xmm11
    movsd [rsp + nb104_fizH2], xmm12
    movsd [rsp + nb104_fixM], xmm13
    movsd [rsp + nb104_fiyM], xmm14
    movsd [rsp + nb104_fizM], xmm15
   
    ;# store back j M forces from xmm0-xmm2
	movsd [rdi + rax*8 + 72], xmm0
	movsd [rdi + rax*8 + 80], xmm1
	movsd [rdi + rax*8 + 88], xmm2
	
.nb104_updateouterdata:
	mov   ecx, [rsp + nb104_ii3]
	mov   rdi, [rbp + nb104_faction]
	mov   rsi, [rbp + nb104_fshift]
	mov   edx, [rsp + nb104_is3]

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb104_fixH1]
	movapd xmm1, [rsp + nb104_fiyH1]
	movapd xmm2, [rsp + nb104_fizH1]

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
	movsd  [rdi + rcx*8 + 24],     xmm3
	movsd  [rdi + rcx*8 + 32], xmm4
	movsd  [rdi + rcx*8 + 40], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movapd xmm6, xmm0
	movsd xmm7, xmm2
	unpcklpd xmm6, xmm1

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb104_fixH2]
	movapd xmm1, [rsp + nb104_fiyH2]
	movapd xmm2, [rsp + nb104_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

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

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb104_fixM]
	movapd xmm1, [rsp + nb104_fiyM]
	movapd xmm2, [rsp + nb104_fizM]

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
	movsd  xmm3, [rdi + rcx*8 + 72]
	movsd  xmm4, [rdi + rcx*8 + 80]
	movsd  xmm5, [rdi + rcx*8 + 88]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8 + 72], xmm3
	movsd  [rdi + rcx*8 + 80], xmm4
	movsd  [rdi + rcx*8 + 88], xmm5

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
	mov esi, [rsp + nb104_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb104_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb104_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	        
	;# add earlier value from mem 
	mov   rax, [rbp + nb104_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 	
	
        ;# finish if last 
        mov ecx, [rsp + nb104_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb104_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb104_n], esi
        jmp .nb104_outer
.nb104_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb104_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb104_end
        ;# non-zero, do one more workunit
        jmp   .nb104_threadloop
.nb104_end:
	mov eax, [rsp + nb104_nouter]
	mov ebx, [rsp + nb104_ninner]
	mov rcx, [rbp + nb104_outeriter]
	mov rdx, [rbp + nb104_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1488
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



	
.globl nb_kernel104nf_x86_64_sse2
.globl _nb_kernel104nf_x86_64_sse2
nb_kernel104nf_x86_64_sse2:	
_nb_kernel104nf_x86_64_sse2:	
.equiv          nb104nf_fshift,         16
.equiv          nb104nf_gid,            24
.equiv          nb104nf_pos,            32
.equiv          nb104nf_faction,        40
.equiv          nb104nf_charge,         48
.equiv          nb104nf_p_facel,        56
.equiv          nb104nf_argkrf,         64
.equiv          nb104nf_argcrf,         72
.equiv          nb104nf_Vc,             80
.equiv          nb104nf_type,           88
.equiv          nb104nf_p_ntype,        96
.equiv          nb104nf_vdwparam,       104
.equiv          nb104nf_Vvdw,           112
.equiv          nb104nf_p_tabscale,     120
.equiv          nb104nf_VFtab,          128
.equiv          nb104nf_invsqrta,       136
.equiv          nb104nf_dvda,           144
.equiv          nb104nf_p_gbtabscale,   152
.equiv          nb104nf_GBtab,          160
.equiv          nb104nf_p_nthreads,     168
.equiv          nb104nf_count,          176
.equiv          nb104nf_mtx,            184
.equiv          nb104nf_outeriter,      192
.equiv          nb104nf_inneriter,      200
.equiv          nb104nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb104nf_ixM,            0
.equiv          nb104nf_iyM,            16
.equiv          nb104nf_izM,            32
.equiv          nb104nf_ixH1,           48
.equiv          nb104nf_iyH1,           64
.equiv          nb104nf_izH1,           80
.equiv          nb104nf_ixH2,           96
.equiv          nb104nf_iyH2,           112
.equiv          nb104nf_izH2,           128
.equiv          nb104nf_jxM,            144
.equiv          nb104nf_jyM,            160
.equiv          nb104nf_jzM,            176
.equiv          nb104nf_jxH1,           192
.equiv          nb104nf_jyH1,           208
.equiv          nb104nf_jzH1,           224
.equiv          nb104nf_jxH2,           240
.equiv          nb104nf_jyH2,           256
.equiv          nb104nf_jzH2,           272
.equiv          nb104nf_qqMM,           288
.equiv          nb104nf_qqMH,           304
.equiv          nb104nf_qqHH,           320
.equiv          nb104nf_vctot,          336
.equiv          nb104nf_half,           352
.equiv          nb104nf_three,          368
.equiv          nb104nf_rsqMM,          384
.equiv          nb104nf_rsqMH1,         400
.equiv          nb104nf_rsqMH2,         416
.equiv          nb104nf_rsqH1M,         432
.equiv          nb104nf_rsqH1H1,        448
.equiv          nb104nf_rsqH1H2,        464
.equiv          nb104nf_rsqH2M,         480
.equiv          nb104nf_rsqH2H1,        496
.equiv          nb104nf_rsqH2H2,        512
.equiv          nb104nf_rinvMM,         528
.equiv          nb104nf_rinvMH1,        544
.equiv          nb104nf_rinvMH2,        560
.equiv          nb104nf_rinvH1M,        576
.equiv          nb104nf_rinvH1H1,       592
.equiv          nb104nf_rinvH1H2,       608
.equiv          nb104nf_rinvH2M,        624
.equiv          nb104nf_rinvH2H1,       640
.equiv          nb104nf_rinvH2H2,       656
.equiv          nb104nf_is3,            672
.equiv          nb104nf_ii3,            676
.equiv          nb104nf_nri,            692
.equiv          nb104nf_iinr,           700
.equiv          nb104nf_jindex,         708
.equiv          nb104nf_jjnr,           716
.equiv          nb104nf_shift,          724
.equiv          nb104nf_shiftvec,       732
.equiv          nb104nf_facel,          740
.equiv          nb104nf_innerjjnr,      748
.equiv          nb104nf_innerk,         756
.equiv          nb104nf_n,              760
.equiv          nb104nf_nn1,            764
.equiv          nb104nf_nouter,         768
.equiv          nb104nf_ninner,         772

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
	sub rsp, 784		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb104nf_nouter], eax
	mov [rsp + nb104nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb104nf_nri], edi
	mov [rsp + nb104nf_iinr], rsi
	mov [rsp + nb104nf_jindex], rdx
	mov [rsp + nb104nf_jjnr], rcx
	mov [rsp + nb104nf_shift], r8
	mov [rsp + nb104nf_shiftvec], r9
	mov rsi, [rbp + nb104nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb104nf_facel], xmm0

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb104nf_half], eax
	mov [rsp + nb104nf_half+4], ebx
	movsd xmm1, [rsp + nb104nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb104nf_half], xmm1
	movapd [rsp + nb104nf_three], xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb104nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb104nf_charge]
	movsd xmm3, [rdx + rbx*8 + 24]	;# qM 
	movsd xmm4, xmm3		;# qM 
	movsd xmm5, [rdx + rbx*8 + 8]	;# qH 
	mov rsi, [rbp + nb104nf_p_facel]
	movsd xmm0, [rsi]
	movsd xmm6, [rsp + nb104nf_facel]	;# facel 
	mulsd  xmm3, xmm3		;# qM*qM 
	mulsd  xmm4, xmm5		;# qM*qH 
	mulsd  xmm5, xmm5		;# qH*qH 
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb104nf_qqMM], xmm3
	movapd [rsp + nb104nf_qqMH], xmm4
	movapd [rsp + nb104nf_qqHH], xmm5

.nb104nf_threadloop:
        mov   rsi, [rbp + nb104nf_count]        ;# pointer to sync counter
        mov   eax, [rsi]
.nb104nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb104nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb104nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb104nf_n], eax
        mov [rsp + nb104nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb104nf_outerstart
        jmp .nb104nf_end

.nb104nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb104nf_nouter]
	mov [rsp + nb104nf_nouter], ebx

.nb104nf_outer:
	mov   rax, [rsp + nb104nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 

	mov   rax, [rsp + nb104nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb104nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb104nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb104nf_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [rax + rbx*8 + 24]
	addsd xmm4, [rax + rbx*8 + 32]
	addsd xmm5, [rax + rbx*8 + 40]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb104nf_ixH1], xmm3
	movapd [rsp + nb104nf_iyH1], xmm4
	movapd [rsp + nb104nf_izH1], xmm5

	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm0, [rax + rbx*8 + 48]
	addsd xmm1, [rax + rbx*8 + 56]
	addsd xmm2, [rax + rbx*8 + 64]		
	addsd xmm3, [rax + rbx*8 + 72]
	addsd xmm4, [rax + rbx*8 + 80]
	addsd xmm5, [rax + rbx*8 + 88]		

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb104nf_ixH2], xmm0
	movapd [rsp + nb104nf_iyH2], xmm1
	movapd [rsp + nb104nf_izH2], xmm2
	movapd [rsp + nb104nf_ixM], xmm3
	movapd [rsp + nb104nf_iyM], xmm4
	movapd [rsp + nb104nf_izM], xmm5

	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb104nf_vctot], xmm4
	
	mov   rax, [rsp + nb104nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb104nf_pos]
	mov   rax, [rsp + nb104nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb104nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb104nf_ninner]
	mov   [rsp + nb104nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb104nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb104nf_unroll_loop
	jmp   .nb104nf_checksingle
.nb104nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb104nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	
	add qword ptr [rsp + nb104nf_innerjjnr], 8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb104nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	
	;# move j coordinates to local temp variables 
	movlpd xmm2, [rsi + rax*8 + 24]
	movlpd xmm3, [rsi + rax*8 + 32]
	movlpd xmm4, [rsi + rax*8 + 40]
	movlpd xmm5, [rsi + rax*8 + 48]
	movlpd xmm6, [rsi + rax*8 + 56]
	movlpd xmm7, [rsi + rax*8 + 64]
	movhpd xmm2, [rsi + rbx*8 + 24]
	movhpd xmm3, [rsi + rbx*8 + 32]
	movhpd xmm4, [rsi + rbx*8 + 40]
	movhpd xmm5, [rsi + rbx*8 + 48]
	movhpd xmm6, [rsi + rbx*8 + 56]
	movhpd xmm7, [rsi + rbx*8 + 64]
	movapd 	[rsp + nb104nf_jxH1], xmm2
	movapd 	[rsp + nb104nf_jyH1], xmm3
	movapd 	[rsp + nb104nf_jzH1], xmm4
	movapd 	[rsp + nb104nf_jxH2], xmm5
	movapd 	[rsp + nb104nf_jyH2], xmm6
	movapd 	[rsp + nb104nf_jzH2], xmm7
	movlpd xmm2, [rsi + rax*8 + 72]
	movlpd xmm3, [rsi + rax*8 + 80]
	movlpd xmm4, [rsi + rax*8 + 88]
	movhpd xmm2, [rsi + rbx*8 + 72]
	movhpd xmm3, [rsi + rbx*8 + 80]
	movhpd xmm4, [rsi + rbx*8 + 88]
	movapd 	[rsp + nb104nf_jxM], xmm2
	movapd 	[rsp + nb104nf_jyM], xmm3
	movapd 	[rsp + nb104nf_jzM], xmm4
	
	movapd xmm0, [rsp + nb104nf_ixM]
	movapd xmm1, [rsp + nb104nf_iyM]
	movapd xmm2, [rsp + nb104nf_izM]
	movapd xmm3, [rsp + nb104nf_ixM]
	movapd xmm4, [rsp + nb104nf_iyM]
	movapd xmm5, [rsp + nb104nf_izM]
	subpd  xmm0, [rsp + nb104nf_jxM]
	subpd  xmm1, [rsp + nb104nf_jyM]
	subpd  xmm2, [rsp + nb104nf_jzM]
	subpd  xmm3, [rsp + nb104nf_jxH1]
	subpd  xmm4, [rsp + nb104nf_jyH1]
	subpd  xmm5, [rsp + nb104nf_jzH1]
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [rsp + nb104nf_rsqMM], xmm0
	movapd [rsp + nb104nf_rsqMH1], xmm3

	movapd xmm0, [rsp + nb104nf_ixM]
	movapd xmm1, [rsp + nb104nf_iyM]
	movapd xmm2, [rsp + nb104nf_izM]
	movapd xmm3, [rsp + nb104nf_ixH1]
	movapd xmm4, [rsp + nb104nf_iyH1]
	movapd xmm5, [rsp + nb104nf_izH1]
	subpd  xmm0, [rsp + nb104nf_jxH2]
	subpd  xmm1, [rsp + nb104nf_jyH2]
	subpd  xmm2, [rsp + nb104nf_jzH2]
	subpd  xmm3, [rsp + nb104nf_jxM]
	subpd  xmm4, [rsp + nb104nf_jyM]
	subpd  xmm5, [rsp + nb104nf_jzM]
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [rsp + nb104nf_rsqMH2], xmm0
	movapd [rsp + nb104nf_rsqH1M], xmm3

	movapd xmm0, [rsp + nb104nf_ixH1]
	movapd xmm1, [rsp + nb104nf_iyH1]
	movapd xmm2, [rsp + nb104nf_izH1]
	movapd xmm3, [rsp + nb104nf_ixH1]
	movapd xmm4, [rsp + nb104nf_iyH1]
	movapd xmm5, [rsp + nb104nf_izH1]
	subpd  xmm0, [rsp + nb104nf_jxH1]
	subpd  xmm1, [rsp + nb104nf_jyH1]
	subpd  xmm2, [rsp + nb104nf_jzH1]
	subpd  xmm3, [rsp + nb104nf_jxH2]
	subpd  xmm4, [rsp + nb104nf_jyH2]
	subpd  xmm5, [rsp + nb104nf_jzH2]
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
	movapd [rsp + nb104nf_rsqH1H1], xmm0
	movapd [rsp + nb104nf_rsqH1H2], xmm3

	movapd xmm0, [rsp + nb104nf_ixH2]
	movapd xmm1, [rsp + nb104nf_iyH2]
	movapd xmm2, [rsp + nb104nf_izH2]
	movapd xmm3, [rsp + nb104nf_ixH2]
	movapd xmm4, [rsp + nb104nf_iyH2]
	movapd xmm5, [rsp + nb104nf_izH2]
	subpd  xmm0, [rsp + nb104nf_jxM]
	subpd  xmm1, [rsp + nb104nf_jyM]
	subpd  xmm2, [rsp + nb104nf_jzM]
	subpd  xmm3, [rsp + nb104nf_jxH1]
	subpd  xmm4, [rsp + nb104nf_jyH1]
	subpd  xmm5, [rsp + nb104nf_jzH1]
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm4, xmm3
	addpd  xmm4, xmm5
	movapd [rsp + nb104nf_rsqH2M], xmm0
	movapd [rsp + nb104nf_rsqH2H1], xmm4

	movapd xmm0, [rsp + nb104nf_ixH2]
	movapd xmm1, [rsp + nb104nf_iyH2]
	movapd xmm2, [rsp + nb104nf_izH2]
	subpd  xmm0, [rsp + nb104nf_jxH2]
	subpd  xmm1, [rsp + nb104nf_jyH2]
	subpd  xmm2, [rsp + nb104nf_jzH2]
	mulpd xmm0, xmm0
	mulpd xmm1, xmm1
	mulpd xmm2, xmm2
	addpd xmm0, xmm1
	addpd xmm0, xmm2
	movapd [rsp + nb104nf_rsqH2H2], xmm0
	
	;# start doing invsqrt use rsq values in xmm0 (h2h2) , xmm4 (h2h1) 
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [rsp + nb104nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb104nf_half] ;# iter1 
	mulpd   xmm7, [rsp + nb104nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb104nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb104nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb104nf_half] ;# rinv 
	movapd [rsp + nb104nf_rinvH2H2], xmm1
	movapd [rsp + nb104nf_rinvH2H1], xmm5

	movapd xmm0, [rsp + nb104nf_rsqMM]
	movapd xmm4, [rsp + nb104nf_rsqMH1]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [rsp + nb104nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb104nf_half] ;# iter1 of  
	mulpd   xmm7, [rsp + nb104nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb104nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb104nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb104nf_half] ;# rinv
	movapd [rsp + nb104nf_rinvMM], xmm1
	movapd [rsp + nb104nf_rinvMH1], xmm5

	movapd xmm0, [rsp + nb104nf_rsqMH2]
	movapd xmm4, [rsp + nb104nf_rsqH1M]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [rsp + nb104nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb104nf_half] ;# iter1 
	mulpd   xmm7, [rsp + nb104nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb104nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb104nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb104nf_half] ;# rinv 
	movapd [rsp + nb104nf_rinvMH2], xmm1
	movapd [rsp + nb104nf_rinvH1M], xmm5

	movapd xmm0, [rsp + nb104nf_rsqH1H1]
	movapd xmm4, [rsp + nb104nf_rsqH1H2]	
	cvtpd2ps xmm1, xmm0	
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm1, xmm1
	rsqrtps xmm5, xmm5
	cvtps2pd xmm1, xmm1
	cvtps2pd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulpd   xmm1, xmm1	;# luA*luA 
	mulpd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [rsp + nb104nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb104nf_half] ;# iter1a 
	mulpd   xmm7, [rsp + nb104nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb104nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb104nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb104nf_half] ;# rinv 
	movapd [rsp + nb104nf_rinvH1H1], xmm1
	movapd [rsp + nb104nf_rinvH1H2], xmm5

	movapd xmm0, [rsp + nb104nf_rsqH2M]
	cvtpd2ps xmm1, xmm0	
	rsqrtps xmm1, xmm1
	cvtps2pd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulpd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [rsp + nb104nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm3, [rsp + nb104nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulpd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [rsp + nb104nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm1, [rsp + nb104nf_half] ;# rinv 
	movapd [rsp + nb104nf_rinvH2M], xmm1

	;# start with MM interaction 
	movapd xmm0, [rsp + nb104nf_rinvMM]
	mulpd  xmm0, [rsp + nb104nf_qqMM]	
	addpd  xmm0, [rsp + nb104nf_vctot]
	
	;# other interactions 
	movapd xmm1, [rsp + nb104nf_rinvMH1]
	movapd xmm2, [rsp + nb104nf_rinvH1H1]
	
	addpd xmm1, [rsp + nb104nf_rinvMH2]
	addpd xmm2, [rsp + nb104nf_rinvH1H2]
	
	addpd xmm1, [rsp + nb104nf_rinvH1M]
	addpd xmm2, [rsp + nb104nf_rinvH2H1]

	addpd xmm1, [rsp + nb104nf_rinvH2M]
	addpd xmm2, [rsp + nb104nf_rinvH2H2]

	mulpd xmm1, [rsp + nb104nf_qqMH]
	mulpd xmm2, [rsp + nb104nf_qqHH]
	
	addpd xmm0, xmm1	
	addpd xmm0, xmm2

	movapd [rsp + nb104nf_vctot], xmm0
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb104nf_innerk],  2
	jl    .nb104nf_checksingle
	jmp   .nb104nf_unroll_loop
.nb104nf_checksingle:
	mov   edx, [rsp + nb104nf_innerk]
	and   edx, 1
	jnz   .nb104nf_dosingle
	jmp   .nb104nf_updateouterdata
.nb104nf_dosingle:
	mov   rdx, [rsp + nb104nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	

	mov rsi, [rbp + nb104nf_pos]
	lea   rax, [rax + rax*2]  

	;# move j coordinates to local temp variables 
	movlpd xmm2, [rsi + rax*8 + 24]
	movlpd xmm3, [rsi + rax*8 + 32]
	movlpd xmm4, [rsi + rax*8 + 40]
	movlpd xmm5, [rsi + rax*8 + 48]
	movlpd xmm6, [rsi + rax*8 + 56]
	movlpd xmm7, [rsi + rax*8 + 64]
	movapd 	[rsp + nb104nf_jxH1], xmm2
	movapd 	[rsp + nb104nf_jyH1], xmm3
	movapd 	[rsp + nb104nf_jzH1], xmm4
	movapd 	[rsp + nb104nf_jxH2], xmm5
	movapd 	[rsp + nb104nf_jyH2], xmm6
	movapd 	[rsp + nb104nf_jzH2], xmm7
	movlpd xmm2, [rsi + rax*8 + 72]
	movlpd xmm3, [rsi + rax*8 + 80]
	movlpd xmm4, [rsi + rax*8 + 88]
	movapd 	[rsp + nb104nf_jxM], xmm2
	movapd 	[rsp + nb104nf_jyM], xmm3
	movapd 	[rsp + nb104nf_jzM], xmm4
	
	movapd xmm0, [rsp + nb104nf_ixM]
	movapd xmm1, [rsp + nb104nf_iyM]
	movapd xmm2, [rsp + nb104nf_izM]
	movapd xmm3, [rsp + nb104nf_ixM]
	movapd xmm4, [rsp + nb104nf_iyM]
	movapd xmm5, [rsp + nb104nf_izM]
	subsd  xmm0, [rsp + nb104nf_jxM]
	subsd  xmm1, [rsp + nb104nf_jyM]
	subsd  xmm2, [rsp + nb104nf_jzM]
	subsd  xmm3, [rsp + nb104nf_jxH1]
	subsd  xmm4, [rsp + nb104nf_jyH1]
	subsd  xmm5, [rsp + nb104nf_jzH1]
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [rsp + nb104nf_rsqMM], xmm0
	movapd [rsp + nb104nf_rsqMH1], xmm3

	movapd xmm0, [rsp + nb104nf_ixM]
	movapd xmm1, [rsp + nb104nf_iyM]
	movapd xmm2, [rsp + nb104nf_izM]
	movapd xmm3, [rsp + nb104nf_ixH1]
	movapd xmm4, [rsp + nb104nf_iyH1]
	movapd xmm5, [rsp + nb104nf_izH1]
	subsd  xmm0, [rsp + nb104nf_jxH2]
	subsd  xmm1, [rsp + nb104nf_jyH2]
	subsd  xmm2, [rsp + nb104nf_jzH2]
	subsd  xmm3, [rsp + nb104nf_jxM]
	subsd  xmm4, [rsp + nb104nf_jyM]
	subsd  xmm5, [rsp + nb104nf_jzM]
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [rsp + nb104nf_rsqMH2], xmm0
	movapd [rsp + nb104nf_rsqH1M], xmm3

	movapd xmm0, [rsp + nb104nf_ixH1]
	movapd xmm1, [rsp + nb104nf_iyH1]
	movapd xmm2, [rsp + nb104nf_izH1]
	movapd xmm3, [rsp + nb104nf_ixH1]
	movapd xmm4, [rsp + nb104nf_iyH1]
	movapd xmm5, [rsp + nb104nf_izH1]
	subsd  xmm0, [rsp + nb104nf_jxH1]
	subsd  xmm1, [rsp + nb104nf_jyH1]
	subsd  xmm2, [rsp + nb104nf_jzH1]
	subsd  xmm3, [rsp + nb104nf_jxH2]
	subsd  xmm4, [rsp + nb104nf_jyH2]
	subsd  xmm5, [rsp + nb104nf_jzH2]
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
	movapd [rsp + nb104nf_rsqH1H1], xmm0
	movapd [rsp + nb104nf_rsqH1H2], xmm3

	movapd xmm0, [rsp + nb104nf_ixH2]
	movapd xmm1, [rsp + nb104nf_iyH2]
	movapd xmm2, [rsp + nb104nf_izH2]
	movapd xmm3, [rsp + nb104nf_ixH2]
	movapd xmm4, [rsp + nb104nf_iyH2]
	movapd xmm5, [rsp + nb104nf_izH2]
	subsd  xmm0, [rsp + nb104nf_jxM]
	subsd  xmm1, [rsp + nb104nf_jyM]
	subsd  xmm2, [rsp + nb104nf_jzM]
	subsd  xmm3, [rsp + nb104nf_jxH1]
	subsd  xmm4, [rsp + nb104nf_jyH1]
	subsd  xmm5, [rsp + nb104nf_jzH1]
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm4, xmm3
	addsd  xmm4, xmm5
	movapd [rsp + nb104nf_rsqH2M], xmm0
	movapd [rsp + nb104nf_rsqH2H1], xmm4

	movapd xmm0, [rsp + nb104nf_ixH2]
	movapd xmm1, [rsp + nb104nf_iyH2]
	movapd xmm2, [rsp + nb104nf_izH2]
	subsd  xmm0, [rsp + nb104nf_jxH2]
	subsd  xmm1, [rsp + nb104nf_jyH2]
	subsd  xmm2, [rsp + nb104nf_jzH2]
	mulsd xmm0, xmm0
	mulsd xmm1, xmm1
	mulsd xmm2, xmm2
	addsd xmm0, xmm1
	addsd xmm0, xmm2
	movapd [rsp + nb104nf_rsqH2H2], xmm0
	
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [rsp + nb104nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb104nf_half] ;# iter1 
	mulsd   xmm7, [rsp + nb104nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb104nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb104nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb104nf_half] ;# rinv 
	movapd [rsp + nb104nf_rinvH2H2], xmm1
	movapd [rsp + nb104nf_rinvH2H1], xmm5

	movapd xmm0, [rsp + nb104nf_rsqMM]
	movapd xmm4, [rsp + nb104nf_rsqMH1]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [rsp + nb104nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb104nf_half] ;# iter1 of  
	mulsd   xmm7, [rsp + nb104nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb104nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb104nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb104nf_half] ;# rinv
	movapd [rsp + nb104nf_rinvMM], xmm1
	movapd [rsp + nb104nf_rinvMH1], xmm5

	movapd xmm0, [rsp + nb104nf_rsqMH2]
	movapd xmm4, [rsp + nb104nf_rsqH1M]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [rsp + nb104nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb104nf_half] ;# iter1 
	mulsd   xmm7, [rsp + nb104nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb104nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb104nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb104nf_half] ;# rinv 
	movapd [rsp + nb104nf_rinvMH2], xmm1
	movapd [rsp + nb104nf_rinvH1M], xmm5

	movapd xmm0, [rsp + nb104nf_rsqH1H1]
	movapd xmm4, [rsp + nb104nf_rsqH1H2]	
	cvtsd2ss xmm1, xmm0	
	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm1, xmm1
	rsqrtss xmm5, xmm5
	cvtss2sd xmm1, xmm1
	cvtss2sd xmm5, xmm5
	
	movapd  xmm2, xmm1	;# copy of luA 
	movapd  xmm6, xmm5	;# copy of luB 
	mulsd   xmm1, xmm1	;# luA*luA 
	mulsd   xmm5, xmm5	;# luB*luB 
	movapd  xmm3, [rsp + nb104nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb104nf_half] ;# iter1a 
	mulsd   xmm7, [rsp + nb104nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb104nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb104nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb104nf_half] ;# rinv 
	movapd [rsp + nb104nf_rinvH1H1], xmm1
	movapd [rsp + nb104nf_rinvH1H2], xmm5

	movapd xmm0, [rsp + nb104nf_rsqH2M]
	cvtsd2ss xmm1, xmm0	
	rsqrtss xmm1, xmm1
	cvtss2sd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulsd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [rsp + nb104nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm3, [rsp + nb104nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulsd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [rsp + nb104nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm1, [rsp + nb104nf_half] ;# rinv 
	movapd [rsp + nb104nf_rinvH2M], xmm1

	;# start with MM interaction 
	movapd xmm0, [rsp + nb104nf_rinvMM]
	mulpd  xmm0, [rsp + nb104nf_qqMM]	
	addpd  xmm0, [rsp + nb104nf_vctot]
	
	;# other interactions 
	movapd xmm1, [rsp + nb104nf_rinvMH1]
	movapd xmm2, [rsp + nb104nf_rinvH1H1]
	
	addsd xmm1, [rsp + nb104nf_rinvMH2]
	addsd xmm2, [rsp + nb104nf_rinvH1H2]
	
	addsd xmm1, [rsp + nb104nf_rinvH1M]
	addsd xmm2, [rsp + nb104nf_rinvH2H1]

	addsd xmm1, [rsp + nb104nf_rinvH2M]
	addsd xmm2, [rsp + nb104nf_rinvH2H2]

	mulsd xmm1, [rsp + nb104nf_qqMH]
	mulsd xmm2, [rsp + nb104nf_qqHH]
	
	addsd xmm0, xmm1	
	addsd xmm0, xmm2

	movlpd [rsp + nb104nf_vctot], xmm0
	
.nb104nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb104nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb104nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	movapd xmm7, [rsp + nb104nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	        
	;# add earlier value from mem 
	mov   rax, [rbp + nb104nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 	
	
        ;# finish if last 
        mov ecx, [rsp + nb104nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb104nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb104nf_n], esi
        jmp .nb104nf_outer
.nb104nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb104nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb104nf_end
        ;# non-zero, do one more workunit
        jmp   .nb104nf_threadloop
.nb104nf_end:
	mov eax, [rsp + nb104nf_nouter]
	mov ebx, [rsp + nb104nf_ninner]
	mov rcx, [rbp + nb104nf_outeriter]
	mov rdx, [rbp + nb104nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 784
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
