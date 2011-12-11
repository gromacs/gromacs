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

	
.globl nb_kernel332_x86_64_sse2
.globl _nb_kernel332_x86_64_sse2
nb_kernel332_x86_64_sse2:	
_nb_kernel332_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb332_fshift,           16
.equiv          nb332_gid,              24
.equiv          nb332_pos,              32
.equiv          nb332_faction,          40
.equiv          nb332_charge,           48
.equiv          nb332_p_facel,          56
.equiv          nb332_argkrf,           64
.equiv          nb332_argcrf,           72
.equiv          nb332_Vc,               80
.equiv          nb332_type,             88
.equiv          nb332_p_ntype,          96
.equiv          nb332_vdwparam,         104
.equiv          nb332_Vvdw,             112
.equiv          nb332_p_tabscale,       120
.equiv          nb332_VFtab,            128
.equiv          nb332_invsqrta,         136
.equiv          nb332_dvda,             144
.equiv          nb332_p_gbtabscale,     152
.equiv          nb332_GBtab,            160
.equiv          nb332_p_nthreads,       168
.equiv          nb332_count,            176
.equiv          nb332_mtx,              184
.equiv          nb332_outeriter,        192
.equiv          nb332_inneriter,        200
.equiv          nb332_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb332_ixO,              0
.equiv          nb332_iyO,              16
.equiv          nb332_izO,              32
.equiv          nb332_ixH1,             48
.equiv          nb332_iyH1,             64
.equiv          nb332_izH1,             80
.equiv          nb332_ixH2,             96
.equiv          nb332_iyH2,             112
.equiv          nb332_izH2,             128
.equiv          nb332_jxO,              144
.equiv          nb332_jyO,              160
.equiv          nb332_jzO,              176
.equiv          nb332_jxH1,             192
.equiv          nb332_jyH1,             208
.equiv          nb332_jzH1,             224
.equiv          nb332_jxH2,             240
.equiv          nb332_jyH2,             256
.equiv          nb332_jzH2,             272
.equiv          nb332_dxOO,             288
.equiv          nb332_dyOO,             304
.equiv          nb332_dzOO,             320
.equiv          nb332_dxOH1,            336
.equiv          nb332_dyOH1,            352
.equiv          nb332_dzOH1,            368
.equiv          nb332_dxOH2,            384
.equiv          nb332_dyOH2,            400
.equiv          nb332_dzOH2,            416
.equiv          nb332_dxH1O,            432
.equiv          nb332_dyH1O,            448
.equiv          nb332_dzH1O,            464
.equiv          nb332_dxH1H1,           480
.equiv          nb332_dyH1H1,           496
.equiv          nb332_dzH1H1,           512
.equiv          nb332_dxH1H2,           528
.equiv          nb332_dyH1H2,           544
.equiv          nb332_dzH1H2,           560
.equiv          nb332_dxH2O,            576
.equiv          nb332_dyH2O,            592
.equiv          nb332_dzH2O,            608
.equiv          nb332_dxH2H1,           624
.equiv          nb332_dyH2H1,           640
.equiv          nb332_dzH2H1,           656
.equiv          nb332_dxH2H2,           672
.equiv          nb332_dyH2H2,           688
.equiv          nb332_dzH2H2,           704
.equiv          nb332_qqOO,             720
.equiv          nb332_qqOH,             736
.equiv          nb332_qqHH,             752
.equiv          nb332_two,              768
.equiv          nb332_tsc,              784
.equiv          nb332_c6,               800
.equiv          nb332_c12,              816
.equiv          nb332_vctot,            832
.equiv          nb332_Vvdwtot,          848
.equiv          nb332_fixO,             864
.equiv          nb332_fiyO,             880
.equiv          nb332_fizO,             896
.equiv          nb332_fixH1,            912
.equiv          nb332_fiyH1,            928
.equiv          nb332_fizH1,            944
.equiv          nb332_fixH2,            960
.equiv          nb332_fiyH2,            976
.equiv          nb332_fizH2,            992
.equiv          nb332_epsO,             1008
.equiv          nb332_epsH1,            1024
.equiv          nb332_epsH2,            1040
.equiv          nb332_fjxH1,            1056
.equiv          nb332_fjyH1,            1072
.equiv          nb332_fjzH1,            1088
.equiv          nb332_fjxH2,            1104
.equiv          nb332_fjyH2,            1120
.equiv          nb332_fjzH2,            1136
.equiv          nb332_half,             1152
.equiv          nb332_three,            1168
.equiv          nb332_rsqOO,            1184
.equiv          nb332_rsqOH1,           1200
.equiv          nb332_rsqOH2,           1216
.equiv          nb332_rsqH1O,           1232
.equiv          nb332_rsqH1H1,          1248
.equiv          nb332_rsqH1H2,          1264
.equiv          nb332_rsqH2O,           1280
.equiv          nb332_rsqH2H1,          1296
.equiv          nb332_rsqH2H2,          1312
.equiv          nb332_rinvOO,           1328
.equiv          nb332_rinvOH1,          1344
.equiv          nb332_rinvOH2,          1360
.equiv          nb332_rinvH1O,          1376
.equiv          nb332_rinvH1H1,         1392
.equiv          nb332_rinvH1H2,         1408
.equiv          nb332_rinvH2O,          1424
.equiv          nb332_rinvH2H1,         1440
.equiv          nb332_rinvH2H2,         1456
.equiv          nb332_fscal,            1472
.equiv          nb332_is3,              1488
.equiv          nb332_ii3,              1492
.equiv          nb332_nri,              1496
.equiv          nb332_iinr,             1504
.equiv          nb332_jindex,           1512
.equiv          nb332_jjnr,             1520
.equiv          nb332_shift,            1528
.equiv          nb332_shiftvec,         1536
.equiv          nb332_facel,            1544
.equiv          nb332_innerjjnr,        1552
.equiv          nb332_innerk,           1560
.equiv          nb332_n,                1564
.equiv          nb332_nn1,              1568
.equiv          nb332_nouter,           1572
.equiv          nb332_ninner,           1576

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
	sub rsp, 1584		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb332_nouter], eax
	mov [rsp + nb332_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb332_nri], edi
	mov [rsp + nb332_iinr], rsi
	mov [rsp + nb332_jindex], rdx
	mov [rsp + nb332_jjnr], rcx
	mov [rsp + nb332_shift], r8
	mov [rsp + nb332_shiftvec], r9
	mov rsi, [rbp + nb332_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb332_facel], xmm0

	mov rax, [rbp + nb332_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb332_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb332_half], eax
	mov [rsp + nb332_half+4], ebx
	movsd xmm1, [rsp + nb332_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb332_half], xmm1
	movapd [rsp + nb332_two], xmm2
	movapd [rsp + nb332_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb332_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb332_charge]
	movsd xmm3, [rdx + rbx*8]	
	movsd xmm4, xmm3	
	movsd xmm5, [rdx + rbx*8 + 8]	
	mov rsi, [rbp + nb332_p_facel]
	movsd xmm0, [rsi]
	movsd xmm6, [rsp + nb332_facel]
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm5
	mulsd  xmm5, xmm5
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb332_qqOO], xmm3
	movapd [rsp + nb332_qqOH], xmm4
	movapd [rsp + nb332_qqHH], xmm5
	
	xorpd xmm0, xmm0
	mov   rdx, [rbp + nb332_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb332_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   rax, [rbp + nb332_vdwparam]
	movsd  xmm0, [rax + rdx*8]
	movsd  xmm1, [rax + rdx*8 + 8]
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	movapd [rsp + nb332_c6], xmm0
	movapd [rsp + nb332_c12], xmm1

.nb332_threadloop:
        mov   rsi, [rbp + nb332_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb332_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb332_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb332_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb332_n], eax
        mov [rsp + nb332_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb332_outerstart
        jmp .nb332_end

.nb332_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb332_nouter]
	mov [rsp + nb332_nouter], ebx

.nb332_outer:
	mov   rax, [rsp + nb332_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb332_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb332_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb332_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb332_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb332_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [rax + rbx*8]
	addsd xmm4, [rax + rbx*8 + 8]
	addsd xmm5, [rax + rbx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb332_ixO], xmm3
	movapd [rsp + nb332_iyO], xmm4
	movapd [rsp + nb332_izO], xmm5

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
	movapd [rsp + nb332_ixH1], xmm0
	movapd [rsp + nb332_iyH1], xmm1
	movapd [rsp + nb332_izH1], xmm2
	movapd [rsp + nb332_ixH2], xmm3
	movapd [rsp + nb332_iyH2], xmm4
	movapd [rsp + nb332_izH2], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb332_vctot], xmm4
	movapd [rsp + nb332_Vvdwtot], xmm4
	movapd [rsp + nb332_fixO], xmm4
	movapd [rsp + nb332_fiyO], xmm4
	movapd [rsp + nb332_fizO], xmm4
	movapd [rsp + nb332_fixH1], xmm4
	movapd [rsp + nb332_fiyH1], xmm4
	movapd [rsp + nb332_fizH1], xmm4
	movapd [rsp + nb332_fixH2], xmm4
	movapd [rsp + nb332_fiyH2], xmm4
	movapd [rsp + nb332_fizH2], xmm4
	
	mov   rax, [rsp + nb332_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb332_pos] 
	mov   rdi, [rbp + nb332_faction]	
	mov   rax, [rsp + nb332_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb332_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb332_ninner]
	mov   [rsp + nb332_ninner], ecx
	add   edx, 0
	mov   [rsp + nb332_innerk], edx    ;# number of innerloop atoms 
	jge   .nb332_unroll_loop
	jmp   .nb332_checksingle
.nb332_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb332_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	
	add qword ptr [rsp + nb332_innerjjnr], 8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb332_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	
	;# move j O coordinates to local temp variables 
    movlpd xmm0, [rsi + rax*8] 
    movlpd xmm1, [rsi + rax*8 + 8] 
    movlpd xmm2, [rsi + rax*8 + 16] 
    movhpd xmm0, [rsi + rbx*8] 
    movhpd xmm1, [rsi + rbx*8 + 8] 
    movhpd xmm2, [rsi + rbx*8 + 16] 

    ;# xmm0 = Ox
    ;# xmm1 = Oy
    ;# xmm2 = Oz
        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb332_ixO]
    subpd xmm1, [rsp + nb332_iyO]
    subpd xmm2, [rsp + nb332_izO]
    subpd xmm3, [rsp + nb332_ixH1]
    subpd xmm4, [rsp + nb332_iyH1]
    subpd xmm5, [rsp + nb332_izH1]
    subpd xmm6, [rsp + nb332_ixH2]
    subpd xmm7, [rsp + nb332_iyH2]
    subpd xmm8, [rsp + nb332_izH2]
    
	movapd [rsp + nb332_dxOO], xmm0
	movapd [rsp + nb332_dyOO], xmm1
	movapd [rsp + nb332_dzOO], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb332_dxH1O], xmm3
	movapd [rsp + nb332_dyH1O], xmm4
	movapd [rsp + nb332_dzH1O], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb332_dxH2O], xmm6
	movapd [rsp + nb332_dyH2O], xmm7
	movapd [rsp + nb332_dzH2O], xmm8
	mulpd  xmm6, xmm6
	mulpd  xmm7, xmm7
	mulpd  xmm8, xmm8
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    addpd  xmm6, xmm7
    addpd  xmm6, xmm8

	;# start doing invsqrt for jO atoms
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
		
	movapd  xmm9, [rsp + nb332_three]
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

	movapd  xmm15, [rsp + nb332_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvOO 
	mulpd   xmm10, xmm15 ;# first iteration for rinvH1O
    mulpd   xmm11, xmm15 ;# first iteration for rinvH2O	

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb332_three]
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

	movapd  xmm15, [rsp + nb332_half]
	mulpd   xmm9, xmm15  ;#  rinvOO 
	mulpd   xmm10, xmm15 ;#   rinvH1O
    mulpd   xmm11, xmm15 ;#   rinvH2O
	
	movapd  [rsp + nb332_rinvOO], xmm9
	movapd  [rsp + nb332_rinvH1O], xmm10
	movapd  [rsp + nb332_rinvH2O], xmm11
	
	;# O interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11
    movapd xmm1, [rsp + nb332_tsc]
    
    mulpd  xmm0, xmm9  ;# r
    mulpd  xmm3, xmm10
    mulpd  xmm6, xmm11
    mulpd  xmm0, xmm1 ;# rtab
    mulpd  xmm3, xmm1
    mulpd  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttpd2dq xmm1, xmm0
    cvttpd2dq xmm4, xmm3
    cvttpd2dq xmm7, xmm6        

    ;# convert back to float
    cvtdq2pd  xmm2, xmm1
    cvtdq2pd  xmm5, xmm4
    cvtdq2pd  xmm8, xmm7
    
    ;# multiply by 4
    pslld   xmm1, 2
    pslld   xmm4, 2
    pslld   xmm7, 2
        
    ;# multiply by three (copy, mult. by two, add back)
    movapd  xmm10, xmm1
    movapd  xmm11, xmm4
    movapd  xmm12, xmm7
    pslld   xmm1, 1
    pslld   xmm4, 1
    pslld   xmm7, 1    
    paddd   xmm1, xmm10
    paddd   xmm4, xmm11
    paddd   xmm7, xmm12    
    
    ;# move to integer registers
    pshufd xmm13, xmm1, 1
    pshufd xmm14, xmm4, 1
    pshufd xmm15, xmm7, 1
    movd    r8d, xmm1
    movd    r10d, xmm4
    movd    r12d, xmm7
    movd    r9d, xmm13
    movd    r11d, xmm14
    movd    r13d, xmm15
        
    mov  rsi, [rbp + nb332_VFtab]

    ;# calculate eps
    subpd     xmm0, xmm2
    subpd     xmm3, xmm5
    subpd     xmm6, xmm8

    movapd    [rsp + nb332_epsO], xmm0
    movapd    [rsp + nb332_epsH1], xmm3
    movapd    [rsp + nb332_epsH2], xmm6

    ;# Load LOTS of table data
    movlpd xmm0,  [rsi + r8*8]
    movlpd xmm1,  [rsi + r8*8 + 8]
    movlpd xmm2,  [rsi + r8*8 + 16]
    movlpd xmm3,  [rsi + r8*8 + 24]
    movlpd xmm4,  [rsi + r10*8]
    movlpd xmm5,  [rsi + r10*8 + 8]
    movlpd xmm6,  [rsi + r10*8 + 16]
    movlpd xmm7,  [rsi + r10*8 + 24]
    movlpd xmm8,  [rsi + r12*8]
    movlpd xmm9,  [rsi + r12*8 + 8]
    movlpd xmm10, [rsi + r12*8 + 16]
    movlpd xmm11, [rsi + r12*8 + 24]
    movhpd xmm0,  [rsi + r9*8]
    movhpd xmm1,  [rsi + r9*8 + 8]
    movhpd xmm2,  [rsi + r9*8 + 16]
    movhpd xmm3,  [rsi + r9*8 + 24]
    movhpd xmm4,  [rsi + r11*8]
    movhpd xmm5,  [rsi + r11*8 + 8]
    movhpd xmm6,  [rsi + r11*8 + 16]
    movhpd xmm7,  [rsi + r11*8 + 24]
    movhpd xmm8,  [rsi + r13*8]
    movhpd xmm9,  [rsi + r13*8 + 8]
    movhpd xmm10, [rsi + r13*8 + 16]
    movhpd xmm11, [rsi + r13*8 + 24]
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    mulpd  xmm3, [rsp + nb332_epsO]   ;# Heps
    mulpd  xmm7, [rsp + nb332_epsH1]
    mulpd  xmm11, [rsp + nb332_epsH2]
    mulpd  xmm2, [rsp + nb332_epsO]   ;# Geps
    mulpd  xmm6, [rsp + nb332_epsH1]
    mulpd  xmm10, [rsp + nb332_epsH2]
    mulpd  xmm3, [rsp + nb332_epsO]   ;# Heps2
    mulpd  xmm7, [rsp + nb332_epsH1]
    mulpd  xmm11, [rsp + nb332_epsH2]

    addpd  xmm1, xmm2   ;# F+Geps
    addpd  xmm5, xmm6
    addpd  xmm9, xmm10 
    addpd  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addpd  xmm5, xmm7
    addpd  xmm9, xmm11 
    addpd  xmm3, xmm3    ;# 2*Heps2
    addpd  xmm7, xmm7
    addpd  xmm11, xmm11
    addpd  xmm3, xmm2    ;# 2*Heps2+Geps
    addpd  xmm7, xmm6  
    addpd  xmm11, xmm10
    addpd  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addpd  xmm7, xmm5
    addpd  xmm11, xmm9
    mulpd  xmm1, [rsp + nb332_epsO]   ;# eps*Fp
    mulpd  xmm5, [rsp + nb332_epsH1]
    mulpd  xmm9, [rsp + nb332_epsH2]
    addpd  xmm1, xmm0     ;# VV
    addpd  xmm5, xmm4
    addpd  xmm9, xmm8
    mulpd  xmm1, [rsp + nb332_qqOO]   ;# VV*qq = vcoul
    mulpd  xmm5, [rsp + nb332_qqOH]
    mulpd  xmm9, [rsp + nb332_qqOH]
    mulpd  xmm3, [rsp + nb332_qqOO]    ;# FF*qq = fij
    mulpd  xmm7, [rsp + nb332_qqOH]
    mulpd  xmm11, [rsp + nb332_qqOH] 

    ;# accumulate vctot
    addpd  xmm1, [rsp + nb332_vctot]
    addpd  xmm5, xmm9
    addpd  xmm1, xmm5
    movapd [rsp + nb332_vctot], xmm1

    movapd xmm2, xmm7
    movapd xmm1, xmm11

    ;# fij coul in xmm3, xmm2, xmm1    
    
    ;# calculate LJ table
    movlpd xmm4,  [rsi + r8*8 + 32]
    movlpd xmm5,  [rsi + r8*8 + 40]
    movlpd xmm6,  [rsi + r8*8 + 48]
    movlpd xmm7,  [rsi + r8*8 + 56]
    movlpd xmm8,  [rsi + r8*8 + 64]
    movlpd xmm9,  [rsi + r8*8 + 72]
    movlpd xmm10, [rsi + r8*8 + 80]
    movlpd xmm11, [rsi + r8*8 + 88]
    movhpd xmm4,  [rsi + r9*8 + 32]
    movhpd xmm5,  [rsi + r9*8 + 40]
    movhpd xmm6,  [rsi + r9*8 + 48]
    movhpd xmm7,  [rsi + r9*8 + 56]
    movhpd xmm8,  [rsi + r9*8 + 64]
    movhpd xmm9,  [rsi + r9*8 + 72]
    movhpd xmm10, [rsi + r9*8 + 80]
    movhpd xmm11, [rsi + r9*8 + 88]
    ;# dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11
    
    movapd xmm0, [rsp + nb332_epsO]
    
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
    movapd xmm12, [rsp + nb332_c6]
    movapd xmm13, [rsp + nb332_c12]
    addpd  xmm5, xmm4 ;# VV
    addpd  xmm9, xmm8

    mulpd  xmm5, xmm12  ;# VV*c6 = vnb6
    mulpd  xmm9, xmm13  ;# VV*c12 = vnb12
    addpd  xmm5, xmm9
    addpd  xmm5, [rsp + nb332_Vvdwtot]
    movapd [rsp + nb332_Vvdwtot], xmm5
        
    mulpd  xmm7, xmm12   ;# FF*c6 = fnb6
    mulpd  xmm11, xmm13   ;# FF*c12  = fnb12
    addpd  xmm7, xmm11
    
    addpd  xmm3, xmm7
    movapd xmm10, [rsp + nb332_tsc]

    mulpd  xmm3, xmm10  ;# fscal
    mulpd  xmm2, xmm10
    mulpd  xmm1, xmm10
        
    ;# move j O forces to xmm11-xmm13
    mov rdi, [rbp + nb332_faction]
	movlpd xmm11, [rdi + rax*8]
	movlpd xmm12, [rdi + rax*8 + 8]
	movlpd xmm13, [rdi + rax*8 + 16]
	movhpd xmm11, [rdi + rbx*8]
	movhpd xmm12, [rdi + rbx*8 + 8]
	movhpd xmm13, [rdi + rbx*8 + 16]

    xorpd  xmm0, xmm0
    xorpd  xmm4, xmm4
    xorpd  xmm8, xmm8
    
    subpd  xmm0, xmm3
    subpd  xmm4, xmm2
    subpd  xmm8, xmm1

    mulpd  xmm0, [rsp + nb332_rinvOO]
    mulpd  xmm4, [rsp + nb332_rinvH1O]
    mulpd  xmm8, [rsp + nb332_rinvH2O]
    
    movapd xmm1, xmm0
    movapd xmm2, xmm0
    movapd xmm3, xmm4
    movapd xmm5, xmm4
    movapd xmm6, xmm8
    movapd xmm7, xmm8

	mulpd xmm0, [rsp + nb332_dxOO]
	mulpd xmm1, [rsp + nb332_dyOO]
	mulpd xmm2, [rsp + nb332_dzOO]
	mulpd xmm3, [rsp + nb332_dxH1O]
	mulpd xmm4, [rsp + nb332_dyH1O]
	mulpd xmm5, [rsp + nb332_dzH1O]
	mulpd xmm6, [rsp + nb332_dxH2O]
	mulpd xmm7, [rsp + nb332_dyH2O]
	mulpd xmm8, [rsp + nb332_dzH2O]

    addpd xmm11,  xmm0
    addpd xmm12, xmm1
    addpd xmm13, xmm2
    addpd xmm0, [rsp + nb332_fixO]
    addpd xmm1, [rsp + nb332_fiyO]
    addpd xmm2, [rsp + nb332_fizO]

    addpd xmm11,  xmm3
    addpd xmm12, xmm4
    addpd xmm13, xmm5
    addpd xmm3, [rsp + nb332_fixH1]
    addpd xmm4, [rsp + nb332_fiyH1]
    addpd xmm5, [rsp + nb332_fizH1]

    addpd xmm11,  xmm6
    addpd xmm12, xmm7
    addpd xmm13, xmm8
    addpd xmm6, [rsp + nb332_fixH2]
    addpd xmm7, [rsp + nb332_fiyH2]
    addpd xmm8, [rsp + nb332_fizH2]

    movapd [rsp + nb332_fixO], xmm0
    movapd [rsp + nb332_fiyO], xmm1
    movapd [rsp + nb332_fizO], xmm2
    movapd [rsp + nb332_fixH1], xmm3
    movapd [rsp + nb332_fiyH1], xmm4
    movapd [rsp + nb332_fizH1], xmm5
    movapd [rsp + nb332_fixH2], xmm6
    movapd [rsp + nb332_fiyH2], xmm7
    movapd [rsp + nb332_fizH2], xmm8
       
    ;# store back j O forces from xmm11-xmm13
	movlpd [rdi + rax*8],      xmm11
	movlpd [rdi + rax*8 + 8],  xmm12
	movlpd [rdi + rax*8 + 16], xmm13
	movhpd [rdi + rbx*8],      xmm11
	movhpd [rdi + rbx*8 + 8],  xmm12
	movhpd [rdi + rbx*8 + 16], xmm13

	;# move j H1 coordinates to local temp variables 
    mov rsi, [rbp + nb332_pos]
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
    
    subpd xmm0, [rsp + nb332_ixO]
    subpd xmm1, [rsp + nb332_iyO]
    subpd xmm2, [rsp + nb332_izO]
    subpd xmm3, [rsp + nb332_ixH1]
    subpd xmm4, [rsp + nb332_iyH1]
    subpd xmm5, [rsp + nb332_izH1]
    subpd xmm6, [rsp + nb332_ixH2]
    subpd xmm7, [rsp + nb332_iyH2]
    subpd xmm8, [rsp + nb332_izH2]
    
	movapd [rsp + nb332_dxOH1], xmm0
	movapd [rsp + nb332_dyOH1], xmm1
	movapd [rsp + nb332_dzOH1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb332_dxH1H1], xmm3
	movapd [rsp + nb332_dyH1H1], xmm4
	movapd [rsp + nb332_dzH1H1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb332_dxH2H1], xmm6
	movapd [rsp + nb332_dyH2H1], xmm7
	movapd [rsp + nb332_dzH2H1], xmm8
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
		
	movapd  xmm9, [rsp + nb332_three]
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

	movapd  xmm15, [rsp + nb332_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvOH1 
	mulpd   xmm10, xmm15 ;# first iteration for rinvH1H1
    mulpd   xmm11, xmm15 ;# first iteration for rinvH2OH1

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb332_three]
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

	movapd  xmm15, [rsp + nb332_half]
	mulpd   xmm9, xmm15  ;#  rinvOH1
	mulpd   xmm10, xmm15 ;#   rinvH1H1
    mulpd   xmm11, xmm15 ;#   rinvH2H1
	
	movapd  [rsp + nb332_rinvOH1], xmm9
	movapd  [rsp + nb332_rinvH1H1], xmm10
	movapd  [rsp + nb332_rinvH2H1], xmm11
	
	;# H1 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movapd xmm1, [rsp + nb332_tsc]
    mulpd  xmm0, xmm9  ;# r
    mulpd  xmm3, xmm10
    mulpd  xmm6, xmm11
    mulpd  xmm0, xmm1 ;# rtab
    mulpd  xmm3, xmm1
    mulpd  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttpd2dq xmm1, xmm0
    cvttpd2dq xmm4, xmm3
    cvttpd2dq xmm7, xmm6        

    ;# convert back to float
    cvtdq2pd  xmm2, xmm1
    cvtdq2pd  xmm5, xmm4
    cvtdq2pd  xmm8, xmm7
    
    ;# multiply by 4
    pslld   xmm1, 2
    pslld   xmm4, 2
    pslld   xmm7, 2
    
    ;# multiply by three (copy, mult. by two, add back)
    movdqa  xmm10, xmm1
    movdqa  xmm11, xmm4
    movdqa  xmm12, xmm7
    pslld   xmm1, 1
    pslld   xmm4, 1
    pslld   xmm7, 1    
    paddd   xmm1, xmm10
    paddd   xmm4, xmm11
    paddd   xmm7, xmm12    
    
    ;# move to integer registers
    pshufd xmm13, xmm1, 1
    pshufd xmm14, xmm4, 1
    pshufd xmm15, xmm7, 1
    movd    r8d, xmm1
    movd    r10d, xmm4
    movd    r12d, xmm7
    movd    r9d, xmm13
    movd    r11d, xmm14
    movd    r13d, xmm15
        
    mov  rsi, [rbp + nb332_VFtab]

    ;# calculate eps
    subpd     xmm0, xmm2
    subpd     xmm3, xmm5
    subpd     xmm6, xmm8

    movapd    [rsp + nb332_epsO], xmm0
    movapd    [rsp + nb332_epsH1], xmm3
    movapd    [rsp + nb332_epsH2], xmm6

    ;# Load LOTS of table data
    movlpd xmm0,  [rsi + r8*8]
    movlpd xmm1,  [rsi + r8*8 + 8]
    movlpd xmm2,  [rsi + r8*8 + 16]
    movlpd xmm3,  [rsi + r8*8 + 24]
    movlpd xmm4,  [rsi + r10*8]
    movlpd xmm5,  [rsi + r10*8 + 8]
    movlpd xmm6,  [rsi + r10*8 + 16]
    movlpd xmm7,  [rsi + r10*8 + 24]
    movlpd xmm8,  [rsi + r12*8]
    movlpd xmm9,  [rsi + r12*8 + 8]
    movlpd xmm10, [rsi + r12*8 + 16]
    movlpd xmm11, [rsi + r12*8 + 24]
    movhpd xmm0,  [rsi + r9*8]
    movhpd xmm1,  [rsi + r9*8 + 8]
    movhpd xmm2,  [rsi + r9*8 + 16]
    movhpd xmm3,  [rsi + r9*8 + 24]
    movhpd xmm4,  [rsi + r11*8]
    movhpd xmm5,  [rsi + r11*8 + 8]
    movhpd xmm6,  [rsi + r11*8 + 16]
    movhpd xmm7,  [rsi + r11*8 + 24]
    movhpd xmm8,  [rsi + r13*8]
    movhpd xmm9,  [rsi + r13*8 + 8]
    movhpd xmm10, [rsi + r13*8 + 16]
    movhpd xmm11, [rsi + r13*8 + 24]
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movapd xmm12, [rsp + nb332_epsO]
    movapd xmm13, [rsp + nb332_epsH1]
    movapd xmm14, [rsp + nb332_epsH2]
    
    mulpd  xmm3, xmm12   ;# Heps
    mulpd  xmm7, xmm13
    mulpd  xmm11, xmm14 
    mulpd  xmm2, xmm12   ;# Geps
    mulpd  xmm6, xmm13
    mulpd  xmm10, xmm14 
    mulpd  xmm3, xmm12   ;# Heps2
    mulpd  xmm7, xmm13
    mulpd  xmm11, xmm14 

    addpd  xmm1, xmm2   ;# F+Geps
    addpd  xmm5, xmm6
    addpd  xmm9, xmm10 
    addpd  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addpd  xmm5, xmm7
    addpd  xmm9, xmm11 
    addpd  xmm3, xmm3    ;# 2*Heps2
    addpd  xmm7, xmm7
    addpd  xmm11, xmm11
    addpd  xmm3, xmm2    ;# 2*Heps2+Geps
    addpd  xmm7, xmm6  
    addpd  xmm11, xmm10
    addpd  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addpd  xmm7, xmm5
    addpd  xmm11, xmm9
    mulpd  xmm1, xmm12   ;# eps*Fp
    mulpd  xmm5, xmm13
    mulpd  xmm9, xmm14
    movapd xmm12, [rsp + nb332_qqOH]
    movapd xmm13, [rsp + nb332_qqHH]
    addpd  xmm1, xmm0     ;# VV
    addpd  xmm5, xmm4
    addpd  xmm9, xmm8
    mulpd  xmm1, xmm12   ;# VV*qq = vcoul
    mulpd  xmm5, xmm13
    mulpd  xmm9, xmm13
    mulpd  xmm3, xmm12    ;# FF*qq = fij
    mulpd  xmm7, xmm13
    mulpd  xmm11, xmm13
    
    ;# accumulate vctot
    addpd  xmm1, [rsp + nb332_vctot]
    addpd  xmm5, xmm9
    addpd  xmm1, xmm5
    movapd [rsp + nb332_vctot], xmm1
    
    movapd xmm10, [rsp + nb332_tsc]
    mulpd  xmm3, xmm10  ;# fscal
    mulpd  xmm7, xmm10
    mulpd  xmm10, xmm11
    
    xorpd  xmm4, xmm4
    xorpd  xmm8, xmm8
    xorpd  xmm11, xmm11
    
    subpd  xmm4, xmm3
    subpd  xmm8, xmm7
    subpd  xmm11, xmm10

    mulpd  xmm4, [rsp + nb332_rinvOH1]
    mulpd  xmm8, [rsp + nb332_rinvH1H1]
    mulpd  xmm11, [rsp + nb332_rinvH2H1]
    
    ;# move j H1 forces to xmm0-xmm2
    mov rdi, [rbp + nb332_faction]
	movlpd xmm0, [rdi + rax*8 + 24]
	movlpd xmm1, [rdi + rax*8 + 32]
	movlpd xmm2, [rdi + rax*8 + 40]
	movhpd xmm0, [rdi + rbx*8 + 24]
	movhpd xmm1, [rdi + rbx*8 + 32]
	movhpd xmm2, [rdi + rbx*8 + 40]

    movapd xmm3, xmm4
    movapd xmm5, xmm4
    movapd xmm7, xmm8
    movapd xmm9, xmm8
    movapd xmm10, xmm11
    movapd xmm12, xmm11

	mulpd xmm3, [rsp + nb332_dxOH1]
	mulpd xmm4, [rsp + nb332_dyOH1]
	mulpd xmm5, [rsp + nb332_dzOH1]
	mulpd xmm7, [rsp + nb332_dxH1H1]
	mulpd xmm8, [rsp + nb332_dyH1H1]
	mulpd xmm9, [rsp + nb332_dzH1H1]
	mulpd xmm10, [rsp + nb332_dxH2H1]
	mulpd xmm11, [rsp + nb332_dyH2H1]
	mulpd xmm12, [rsp + nb332_dzH2H1]

    addpd xmm0, xmm3
    addpd xmm1, xmm4
    addpd xmm2, xmm5
    addpd xmm3, [rsp + nb332_fixO]
    addpd xmm4, [rsp + nb332_fiyO]
    addpd xmm5, [rsp + nb332_fizO]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb332_fixH1]
    addpd xmm8, [rsp + nb332_fiyH1]
    addpd xmm9, [rsp + nb332_fizH1]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb332_fixH2]
    addpd xmm11, [rsp + nb332_fiyH2]
    addpd xmm12, [rsp + nb332_fizH2]

    movapd [rsp + nb332_fixO], xmm3
    movapd [rsp + nb332_fiyO], xmm4
    movapd [rsp + nb332_fizO], xmm5
    movapd [rsp + nb332_fixH1], xmm7
    movapd [rsp + nb332_fiyH1], xmm8
    movapd [rsp + nb332_fizH1], xmm9
    movapd [rsp + nb332_fixH2], xmm10
    movapd [rsp + nb332_fiyH2], xmm11
    movapd [rsp + nb332_fizH2], xmm12
   
    ;# store back j H1 forces from xmm0-xmm2
	movlpd [rdi + rax*8 + 24], xmm0
	movlpd [rdi + rax*8 + 32], xmm1
	movlpd [rdi + rax*8 + 40], xmm2
	movhpd [rdi + rbx*8 + 24], xmm0
	movhpd [rdi + rbx*8 + 32], xmm1
	movhpd [rdi + rbx*8 + 40], xmm2
       
	;# move j H2 coordinates to local temp variables 
    mov rsi, [rbp + nb332_pos]
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
    
    subpd xmm0, [rsp + nb332_ixO]
    subpd xmm1, [rsp + nb332_iyO]
    subpd xmm2, [rsp + nb332_izO]
    subpd xmm3, [rsp + nb332_ixH1]
    subpd xmm4, [rsp + nb332_iyH1]
    subpd xmm5, [rsp + nb332_izH1]
    subpd xmm6, [rsp + nb332_ixH2]
    subpd xmm7, [rsp + nb332_iyH2]
    subpd xmm8, [rsp + nb332_izH2]
    
	movapd [rsp + nb332_dxOH2], xmm0
	movapd [rsp + nb332_dyOH2], xmm1
	movapd [rsp + nb332_dzOH2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb332_dxH1H2], xmm3
	movapd [rsp + nb332_dyH1H2], xmm4
	movapd [rsp + nb332_dzH1H2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb332_dxH2H2], xmm6
	movapd [rsp + nb332_dyH2H2], xmm7
	movapd [rsp + nb332_dzH2H2], xmm8
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
		
	movapd  xmm9, [rsp + nb332_three]
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

	movapd  xmm15, [rsp + nb332_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvOH2 
	mulpd   xmm10, xmm15 ;# first iteration for rinvH1H2
    mulpd   xmm11, xmm15 ;# first iteration for rinvH2H2

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb332_three]
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

	movapd  xmm15, [rsp + nb332_half]
	mulpd   xmm9, xmm15  ;#  rinvOH2
	mulpd   xmm10, xmm15 ;#   rinvH1H2
    mulpd   xmm11, xmm15 ;#   rinvH2H2
    
	
	movapd  [rsp + nb332_rinvOH2], xmm9
	movapd  [rsp + nb332_rinvH1H2], xmm10
	movapd  [rsp + nb332_rinvH2H2], xmm11
	
	;# H2 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movapd xmm1, [rsp + nb332_tsc]
    mulpd  xmm0, xmm9  ;# r
    mulpd  xmm3, xmm10
    mulpd  xmm6, xmm11
    mulpd  xmm0, xmm1 ;# rtab
    mulpd  xmm3, xmm1
    mulpd  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttpd2dq xmm1, xmm0
    cvttpd2dq xmm4, xmm3
    cvttpd2dq xmm7, xmm6        

    ;# convert back to float
    cvtdq2pd  xmm2, xmm1
    cvtdq2pd  xmm5, xmm4
    cvtdq2pd  xmm8, xmm7
    
    ;# multiply by 4
    pslld   xmm1, 2
    pslld   xmm4, 2
    pslld   xmm7, 2
    
    ;# multiply by three (copy, mult. by two, add back)
    movdqa  xmm10, xmm1
    movdqa  xmm11, xmm4
    movdqa  xmm12, xmm7
    pslld   xmm1, 1
    pslld   xmm4, 1
    pslld   xmm7, 1    
    paddd   xmm1, xmm10
    paddd   xmm4, xmm11
    paddd   xmm7, xmm12    
    
    ;# move to integer registers
    pshufd xmm13, xmm1, 1
    pshufd xmm14, xmm4, 1
    pshufd xmm15, xmm7, 1
    movd    r8d, xmm1
    movd    r10d, xmm4
    movd    r12d, xmm7
    movd    r9d, xmm13
    movd    r11d, xmm14
    movd    r13d, xmm15
        
    mov  rsi, [rbp + nb332_VFtab]

    ;# calculate eps
    subpd     xmm0, xmm2
    subpd     xmm3, xmm5
    subpd     xmm6, xmm8

    movapd    [rsp + nb332_epsO], xmm0
    movapd    [rsp + nb332_epsH1], xmm3
    movapd    [rsp + nb332_epsH2], xmm6

    ;# Load LOTS of table data
    movlpd xmm0,  [rsi + r8*8]
    movlpd xmm1,  [rsi + r8*8 + 8]
    movlpd xmm2,  [rsi + r8*8 + 16]
    movlpd xmm3,  [rsi + r8*8 + 24]
    movlpd xmm4,  [rsi + r10*8]
    movlpd xmm5,  [rsi + r10*8 + 8]
    movlpd xmm6,  [rsi + r10*8 + 16]
    movlpd xmm7,  [rsi + r10*8 + 24]
    movlpd xmm8,  [rsi + r12*8]
    movlpd xmm9,  [rsi + r12*8 + 8]
    movlpd xmm10, [rsi + r12*8 + 16]
    movlpd xmm11, [rsi + r12*8 + 24]
    movhpd xmm0,  [rsi + r9*8]
    movhpd xmm1,  [rsi + r9*8 + 8]
    movhpd xmm2,  [rsi + r9*8 + 16]
    movhpd xmm3,  [rsi + r9*8 + 24]
    movhpd xmm4,  [rsi + r11*8]
    movhpd xmm5,  [rsi + r11*8 + 8]
    movhpd xmm6,  [rsi + r11*8 + 16]
    movhpd xmm7,  [rsi + r11*8 + 24]
    movhpd xmm8,  [rsi + r13*8]
    movhpd xmm9,  [rsi + r13*8 + 8]
    movhpd xmm10, [rsi + r13*8 + 16]
    movhpd xmm11, [rsi + r13*8 + 24]
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movapd xmm12, [rsp + nb332_epsO]
    movapd xmm13, [rsp + nb332_epsH1]
    movapd xmm14, [rsp + nb332_epsH2]
    
    mulpd  xmm3, xmm12   ;# Heps
    mulpd  xmm7, xmm13
    mulpd  xmm11, xmm14 
    mulpd  xmm2, xmm12   ;# Geps
    mulpd  xmm6, xmm13
    mulpd  xmm10, xmm14 
    mulpd  xmm3, xmm12   ;# Heps2
    mulpd  xmm7, xmm13
    mulpd  xmm11, xmm14 

    addpd  xmm1, xmm2   ;# F+Geps
    addpd  xmm5, xmm6
    addpd  xmm9, xmm10 
    addpd  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addpd  xmm5, xmm7
    addpd  xmm9, xmm11 
    addpd  xmm3, xmm3    ;# 2*Heps2
    addpd  xmm7, xmm7
    addpd  xmm11, xmm11
    addpd  xmm3, xmm2    ;# 2*Heps2+Geps
    addpd  xmm7, xmm6  
    addpd  xmm11, xmm10
    addpd  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addpd  xmm7, xmm5
    addpd  xmm11, xmm9
    mulpd  xmm1, xmm12   ;# eps*Fp
    mulpd  xmm5, xmm13
    mulpd  xmm9, xmm14
    movapd xmm12, [rsp + nb332_qqOH]
    movapd xmm13, [rsp + nb332_qqHH]
    addpd  xmm1, xmm0     ;# VV
    addpd  xmm5, xmm4
    addpd  xmm9, xmm8
    mulpd  xmm1, xmm12   ;# VV*qq = vcoul
    mulpd  xmm5, xmm13
    mulpd  xmm9, xmm13
    mulpd  xmm3, xmm12    ;# FF*qq = fij
    mulpd  xmm7, xmm13
    mulpd  xmm11, xmm13
    
    ;# accumulate vctot
    addpd  xmm1, [rsp + nb332_vctot]
    addpd  xmm5, xmm9
    addpd  xmm1, xmm5
    movapd [rsp + nb332_vctot], xmm1
    
    movapd xmm10, [rsp + nb332_tsc]
    mulpd  xmm3, xmm10  ;# fscal
    mulpd  xmm7, xmm10
    mulpd  xmm10, xmm11
        
    xorpd  xmm4, xmm4
    xorpd  xmm8, xmm8
    xorpd  xmm11, xmm11
    
    subpd  xmm4, xmm3
    subpd  xmm8, xmm7
    subpd  xmm11, xmm10

    mulpd  xmm4, [rsp + nb332_rinvOH2]
    mulpd  xmm8, [rsp + nb332_rinvH1H2]
    mulpd  xmm11, [rsp + nb332_rinvH2H2]
    
    ;# move j H2 forces to xmm0-xmm2
    mov rdi, [rbp + nb332_faction]
	movlpd xmm0, [rdi + rax*8 + 48]
	movlpd xmm1, [rdi + rax*8 + 56]
	movlpd xmm2, [rdi + rax*8 + 64]
	movhpd xmm0, [rdi + rbx*8 + 48]
	movhpd xmm1, [rdi + rbx*8 + 56]
	movhpd xmm2, [rdi + rbx*8 + 64]

    movapd xmm3, xmm4
    movapd xmm5, xmm4
    movapd xmm7, xmm8
    movapd xmm9, xmm8
    movapd xmm10, xmm11
    movapd xmm12, xmm11

	mulpd xmm3, [rsp + nb332_dxOH2]
	mulpd xmm4, [rsp + nb332_dyOH2]
	mulpd xmm5, [rsp + nb332_dzOH2]
	mulpd xmm7, [rsp + nb332_dxH1H2]
	mulpd xmm8, [rsp + nb332_dyH1H2]
	mulpd xmm9, [rsp + nb332_dzH1H2]
	mulpd xmm10, [rsp + nb332_dxH2H2]
	mulpd xmm11, [rsp + nb332_dyH2H2]
	mulpd xmm12, [rsp + nb332_dzH2H2]

    addpd xmm0, xmm3
    addpd xmm1, xmm4
    addpd xmm2, xmm5
    addpd xmm3, [rsp + nb332_fixO]
    addpd xmm4, [rsp + nb332_fiyO]
    addpd xmm5, [rsp + nb332_fizO]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb332_fixH1]
    addpd xmm8, [rsp + nb332_fiyH1]
    addpd xmm9, [rsp + nb332_fizH1]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb332_fixH2]
    addpd xmm11, [rsp + nb332_fiyH2]
    addpd xmm12, [rsp + nb332_fizH2]

    movapd [rsp + nb332_fixO], xmm3
    movapd [rsp + nb332_fiyO], xmm4
    movapd [rsp + nb332_fizO], xmm5
    movapd [rsp + nb332_fixH1], xmm7
    movapd [rsp + nb332_fiyH1], xmm8
    movapd [rsp + nb332_fizH1], xmm9
    movapd [rsp + nb332_fixH2], xmm10
    movapd [rsp + nb332_fiyH2], xmm11
    movapd [rsp + nb332_fizH2], xmm12
   
    ;# store back j H2 forces from xmm0-xmm2
	movlpd [rdi + rax*8 + 48], xmm0
	movlpd [rdi + rax*8 + 56], xmm1
	movlpd [rdi + rax*8 + 64], xmm2
	movhpd [rdi + rbx*8 + 48], xmm0
	movhpd [rdi + rbx*8 + 56], xmm1
	movhpd [rdi + rbx*8 + 64], xmm2
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb332_innerk],  2
	jl    .nb332_checksingle
	jmp   .nb332_unroll_loop
.nb332_checksingle:
	mov   edx, [rsp + nb332_innerk]
	and   edx, 1
	jnz   .nb332_dosingle
	jmp   .nb332_updateouterdata
.nb332_dosingle:
	mov   rdx, [rsp + nb332_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	

	mov rsi, [rbp + nb332_pos]
	lea   rax, [rax + rax*2]  

	
	;# move j O coordinates to local temp variables 
    movsd xmm0, [rsi + rax*8] 
    movsd xmm1, [rsi + rax*8 + 8] 
    movsd xmm2, [rsi + rax*8 + 16] 

    ;# xmm0 = Ox
    ;# xmm1 = Oy
    ;# xmm2 = Oz
        
    movsd xmm3, xmm0
    movsd xmm4, xmm1
    movsd xmm5, xmm2
    movsd xmm6, xmm0
    movsd xmm7, xmm1
    movsd xmm8, xmm2
    
    subsd xmm0, [rsp + nb332_ixO]
    subsd xmm1, [rsp + nb332_iyO]
    subsd xmm2, [rsp + nb332_izO]
    subsd xmm3, [rsp + nb332_ixH1]
    subsd xmm4, [rsp + nb332_iyH1]
    subsd xmm5, [rsp + nb332_izH1]
    subsd xmm6, [rsp + nb332_ixH2]
    subsd xmm7, [rsp + nb332_iyH2]
    subsd xmm8, [rsp + nb332_izH2]
    
	movsd [rsp + nb332_dxOO], xmm0
	movsd [rsp + nb332_dyOO], xmm1
	movsd [rsp + nb332_dzOO], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movsd [rsp + nb332_dxH1O], xmm3
	movsd [rsp + nb332_dyH1O], xmm4
	movsd [rsp + nb332_dzH1O], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movsd [rsp + nb332_dxH2O], xmm6
	movsd [rsp + nb332_dyH2O], xmm7
	movsd [rsp + nb332_dzH2O], xmm8
	mulsd  xmm6, xmm6
	mulsd  xmm7, xmm7
	mulsd  xmm8, xmm8
	addsd  xmm0, xmm1
	addsd  xmm0, xmm2
	addsd  xmm3, xmm4
	addsd  xmm3, xmm5
    addsd  xmm6, xmm7
    addsd  xmm6, xmm8

	;# start doing invsqrt for jO atoms
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
		
	movsd  xmm9, [rsp + nb332_three]
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

	movsd  xmm15, [rsp + nb332_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvOO 
	mulsd   xmm10, xmm15 ;# first iteration for rinvH1O
    mulsd   xmm11, xmm15 ;# first iteration for rinvH2O	

    ;# second iteration step    
	movsd  xmm2, xmm9
	movsd  xmm5, xmm10
    movsd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movsd  xmm1, [rsp + nb332_three]
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

	movsd  xmm15, [rsp + nb332_half]
	mulsd   xmm9, xmm15  ;#  rinvOO 
	mulsd   xmm10, xmm15 ;#   rinvH1O
    mulsd   xmm11, xmm15 ;#   rinvH2O
		
	movsd  [rsp + nb332_rinvOO], xmm9
	movsd  [rsp + nb332_rinvH1O], xmm10
	movsd  [rsp + nb332_rinvH2O], xmm11
	
	;# O interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movapd [rsp + nb332_rinvOO], xmm9
    movsd  xmm1, [rsp + nb332_tsc]

    mulsd  xmm0, xmm9  ;# r
    mulsd  xmm3, xmm10
    mulsd  xmm6, xmm11
    mulsd  xmm0, xmm1 ;# rtab
    mulsd  xmm3, xmm1
    mulsd  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttsd2si r8d, xmm0
    cvttsd2si r10d, xmm3
    cvttsd2si r12d, xmm6        

    ;# convert back to float
    cvtsi2sd  xmm2, r8d
    cvtsi2sd  xmm5, r10d
    cvtsi2sd  xmm8, r12d
    
    ;# multiply by 4
    shl   r8, 2
    shl   r10d, 2
    shl   r12d, 2
    
	;# multiply by 3
    lea   r8, [r8 + r8*2]        
    lea   r10, [r10 + r10*2]        
    lea   r12, [r12 + r12*2]        

    mov  rsi, [rbp + nb332_VFtab]

    ;# calculate eps
    subsd     xmm0, xmm2
    subsd     xmm3, xmm5
    subsd     xmm6, xmm8

    movsd    [rsp + nb332_epsO], xmm0
    movsd    [rsp + nb332_epsH1], xmm3
    movsd    [rsp + nb332_epsH2], xmm6

    ;# Load LOTS of table data
    movsd xmm0,  [rsi + r8*8]
    movsd xmm1,  [rsi + r8*8 + 8]
    movsd xmm2,  [rsi + r8*8 + 16]
    movsd xmm3,  [rsi + r8*8 + 24]
    movsd xmm4,  [rsi + r10*8]
    movsd xmm5,  [rsi + r10*8 + 8]
    movsd xmm6,  [rsi + r10*8 + 16]
    movsd xmm7,  [rsi + r10*8 + 24]
    movsd xmm8,  [rsi + r12*8]
    movsd xmm9,  [rsi + r12*8 + 8]
    movsd xmm10, [rsi + r12*8 + 16]
    movsd xmm11, [rsi + r12*8 + 24]
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movsd xmm12, [rsp + nb332_epsO]
    movsd xmm13, [rsp + nb332_epsH1]
    movsd xmm14, [rsp + nb332_epsH2]
    
    mulsd  xmm3, xmm12   ;# Heps
    mulsd  xmm7, xmm13
    mulsd  xmm11, xmm14 
    mulsd  xmm2, xmm12   ;# Geps
    mulsd  xmm6, xmm13
    mulsd  xmm10, xmm14 
    mulsd  xmm3, xmm12   ;# Heps2
    mulsd  xmm7, xmm13
    mulsd  xmm11, xmm14 

    addsd  xmm1, xmm2   ;# F+Geps
    addsd  xmm5, xmm6
    addsd  xmm9, xmm10 
    addsd  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addsd  xmm5, xmm7
    addsd  xmm9, xmm11 
    addsd  xmm3, xmm3    ;# 2*Heps2
    addsd  xmm7, xmm7
    addsd  xmm11, xmm11
    addsd  xmm3, xmm2    ;# 2*Heps2+Geps
    addsd  xmm7, xmm6  
    addsd  xmm11, xmm10
    addsd  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addsd  xmm7, xmm5
    addsd  xmm11, xmm9
    mulsd  xmm1, [rsp + nb332_epsO]   ;# eps*Fp
    mulsd  xmm5, [rsp + nb332_epsH1]
    mulsd  xmm9, [rsp + nb332_epsH2]
    movsd xmm12, [rsp + nb332_qqOO]
    movsd xmm13, [rsp + nb332_qqOH]
    addsd  xmm1, xmm0     ;# VV
    addsd  xmm5, xmm4
    addsd  xmm9, xmm8
    mulsd  xmm1, [rsp + nb332_qqOO]   ;# VV*qq = vcoul
    mulsd  xmm5, [rsp + nb332_qqOH]
    mulsd  xmm9, [rsp + nb332_qqOH]
    mulsd  xmm3, [rsp + nb332_qqOO]    ;# FF*qq = fij
    mulsd  xmm7, [rsp + nb332_qqOH]
    mulsd  xmm11, [rsp + nb332_qqOH] 
    
    ;# accumulate vctot
    addsd  xmm1, [rsp + nb332_vctot]
    addsd  xmm5, xmm9
    addsd  xmm1, xmm5
    movsd [rsp + nb332_vctot], xmm1
    
    movsd xmm2, xmm7
    movsd xmm1, xmm11

    ;# fij coul in xmm3, xmm2, xmm1    
    
    ;# calculate LJ table
    movsd xmm4,  [rsi + r8*8 + 32]
    movsd xmm5,  [rsi + r8*8 + 40]
    movsd xmm6,  [rsi + r8*8 + 48]
    movsd xmm7,  [rsi + r8*8 + 56]
    movsd xmm8,  [rsi + r8*8 + 64]
    movsd xmm9,  [rsi + r8*8 + 72]
    movsd xmm10, [rsi + r8*8 + 80]
    movsd xmm11, [rsi + r8*8 + 88]
    ;# dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11
    
    movsd xmm0, [rsp + nb332_epsO]
    
    mulsd  xmm7, xmm0    ;# Heps
    mulsd  xmm11, xmm0 
    mulsd  xmm6, xmm0   ;# Geps
    mulsd  xmm10, xmm0 
    mulsd  xmm7, xmm0   ;# Heps2
    mulsd  xmm11, xmm0 
    addsd  xmm5, xmm6  ;# F+Geps
    addsd  xmm9, xmm10 
    addsd  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addsd  xmm9, xmm11 
    addsd  xmm7, xmm7    ;# 2*Heps2
    addsd  xmm11, xmm11
    addsd  xmm7, xmm6   ;# 2*Heps2+Geps
    addsd  xmm11, xmm10
    
    addsd  xmm7, xmm5  ;# FF = Fp + 2*Heps2 + Geps
    addsd  xmm11, xmm9
    mulsd  xmm5, xmm0  ;# eps*Fp
    mulsd  xmm9, xmm0
    movsd xmm12, [rsp + nb332_c6]
    movsd xmm13, [rsp + nb332_c12]
    addsd  xmm5, xmm4 ;# VV
    addsd  xmm9, xmm8

    mulsd  xmm5, xmm12  ;# VV*c6 = vnb6
    mulsd  xmm9, xmm13  ;# VV*c12 = vnb12
    addsd  xmm5, xmm9
    addsd  xmm5, [rsp + nb332_Vvdwtot]
    movsd [rsp + nb332_Vvdwtot], xmm5
        
    mulsd  xmm7, xmm12   ;# FF*c6 = fnb6
    mulsd  xmm11, xmm13   ;# FF*c12  = fnb12
    addsd  xmm7, xmm11
    
    addsd  xmm3, xmm7
    movsd xmm10, [rsp + nb332_tsc]

    mulsd  xmm3, xmm10  ;# fscal
    mulsd  xmm2, xmm10
    mulsd  xmm1, xmm10
        
    ;# move j O forces to xmm11-xmm13
    mov rdi, [rbp + nb332_faction]
	movsd xmm11, [rdi + rax*8]
	movsd xmm12, [rdi + rax*8 + 8]
	movsd xmm13, [rdi + rax*8 + 16]

    xorpd  xmm0, xmm0
    xorpd  xmm4, xmm4
    xorpd  xmm8, xmm8
    
    subsd  xmm0, xmm3
    subsd  xmm4, xmm2
    subsd  xmm8, xmm1

    mulsd  xmm0, [rsp + nb332_rinvOO]
    mulsd  xmm4, [rsp + nb332_rinvH1O]
    mulsd  xmm8, [rsp + nb332_rinvH2O]
    
    movsd xmm1, xmm0
    movsd xmm2, xmm0
    movsd xmm3, xmm4
    movsd xmm5, xmm4
    movsd xmm6, xmm8
    movsd xmm7, xmm8

	mulsd xmm0, [rsp + nb332_dxOO]
	mulsd xmm1, [rsp + nb332_dyOO]
	mulsd xmm2, [rsp + nb332_dzOO]
	mulsd xmm3, [rsp + nb332_dxH1O]
	mulsd xmm4, [rsp + nb332_dyH1O]
	mulsd xmm5, [rsp + nb332_dzH1O]
	mulsd xmm6, [rsp + nb332_dxH2O]
	mulsd xmm7, [rsp + nb332_dyH2O]
	mulsd xmm8, [rsp + nb332_dzH2O]

    addsd xmm11,  xmm0
    addsd xmm12, xmm1
    addsd xmm13, xmm2
    addsd xmm0, [rsp + nb332_fixO]
    addsd xmm1, [rsp + nb332_fiyO]
    addsd xmm2, [rsp + nb332_fizO]

    addsd xmm11,  xmm3
    addsd xmm12, xmm4
    addsd xmm13, xmm5
    addsd xmm3, [rsp + nb332_fixH1]
    addsd xmm4, [rsp + nb332_fiyH1]
    addsd xmm5, [rsp + nb332_fizH1]

    addsd xmm11,  xmm6
    addsd xmm12, xmm7
    addsd xmm13, xmm8
    addsd xmm6, [rsp + nb332_fixH2]
    addsd xmm7, [rsp + nb332_fiyH2]
    addsd xmm8, [rsp + nb332_fizH2]

    movsd [rsp + nb332_fixO], xmm0
    movsd [rsp + nb332_fiyO], xmm1
    movsd [rsp + nb332_fizO], xmm2
    movsd [rsp + nb332_fixH1], xmm3
    movsd [rsp + nb332_fiyH1], xmm4
    movsd [rsp + nb332_fizH1], xmm5
    movsd [rsp + nb332_fixH2], xmm6
    movsd [rsp + nb332_fiyH2], xmm7
    movsd [rsp + nb332_fizH2], xmm8
          
    ;# store back j O forces from xmm11-xmm13
	movsd [rdi + rax*8],      xmm11
	movsd [rdi + rax*8 + 8],  xmm12
	movsd [rdi + rax*8 + 16], xmm13

	;# move j H1 coordinates to local temp variables 
    mov rsi, [rbp + nb332_pos]
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
    
    subsd xmm0, [rsp + nb332_ixO]
    subsd xmm1, [rsp + nb332_iyO]
    subsd xmm2, [rsp + nb332_izO]
    subsd xmm3, [rsp + nb332_ixH1]
    subsd xmm4, [rsp + nb332_iyH1]
    subsd xmm5, [rsp + nb332_izH1]
    subsd xmm6, [rsp + nb332_ixH2]
    subsd xmm7, [rsp + nb332_iyH2]
    subsd xmm8, [rsp + nb332_izH2]
    
	movsd [rsp + nb332_dxOH1], xmm0
	movsd [rsp + nb332_dyOH1], xmm1
	movsd [rsp + nb332_dzOH1], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movsd [rsp + nb332_dxH1H1], xmm3
	movsd [rsp + nb332_dyH1H1], xmm4
	movsd [rsp + nb332_dzH1H1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movsd [rsp + nb332_dxH2H1], xmm6
	movsd [rsp + nb332_dyH2H1], xmm7
	movsd [rsp + nb332_dzH2H1], xmm8
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
		
	movsd  xmm9, [rsp + nb332_three]
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

	movsd  xmm15, [rsp + nb332_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvOH1 
	mulsd   xmm10, xmm15 ;# first iteration for rinvH1H1
    mulsd   xmm11, xmm15 ;# first iteration for rinvH2OH1

    ;# second iteration step    
	movsd  xmm2, xmm9
	movsd  xmm5, xmm10
    movsd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movsd  xmm1, [rsp + nb332_three]
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

	movsd  xmm15, [rsp + nb332_half]
	mulsd   xmm9, xmm15  ;#  rinvOH1
	mulsd   xmm10, xmm15 ;#   rinvH1H1
    mulsd   xmm11, xmm15 ;#   rinvH2H1
	
	movsd  [rsp + nb332_rinvOH1], xmm9
	movsd  [rsp + nb332_rinvH1H1], xmm10
	movsd  [rsp + nb332_rinvH2H1], xmm11
	
	;# H1 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movsd  xmm1, [rsp + nb332_tsc]
    mulsd  xmm0, xmm9  ;# r
    mulsd  xmm3, xmm10
    mulsd  xmm6, xmm11
    mulsd  xmm0, xmm1 ;# rtab
    mulsd  xmm3, xmm1
    mulsd  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttsd2si r8d, xmm0
    cvttsd2si r10d, xmm3
    cvttsd2si r12d, xmm6        

    ;# convert back to float
    cvtsi2sd  xmm2, r8d
    cvtsi2sd  xmm5, r10d
    cvtsi2sd  xmm8, r12d
    
    ;# multiply by 4
    shl   r8, 2
    shl   r10d, 2
    shl   r12d, 2
    
	;# multiply by 3
    lea   r8, [r8 + r8*2]        
    lea   r10, [r10 + r10*2]        
    lea   r12, [r12 + r12*2]        

    mov  rsi, [rbp + nb332_VFtab]

    ;# calculate eps
    subsd     xmm0, xmm2
    subsd     xmm3, xmm5
    subsd     xmm6, xmm8

    movsd    [rsp + nb332_epsO], xmm0
    movsd    [rsp + nb332_epsH1], xmm3
    movsd    [rsp + nb332_epsH2], xmm6

    ;# Load LOTS of table data
    movsd xmm0,  [rsi + r8*8]
    movsd xmm1,  [rsi + r8*8 + 8]
    movsd xmm2,  [rsi + r8*8 + 16]
    movsd xmm3,  [rsi + r8*8 + 24]
    movsd xmm4,  [rsi + r10*8]
    movsd xmm5,  [rsi + r10*8 + 8]
    movsd xmm6,  [rsi + r10*8 + 16]
    movsd xmm7,  [rsi + r10*8 + 24]
    movsd xmm8,  [rsi + r12*8]
    movsd xmm9,  [rsi + r12*8 + 8]
    movsd xmm10, [rsi + r12*8 + 16]
    movsd xmm11, [rsi + r12*8 + 24]
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movsd xmm12, [rsp + nb332_epsO]
    movsd xmm13, [rsp + nb332_epsH1]
    movsd xmm14, [rsp + nb332_epsH2]
    
    mulsd  xmm3, xmm12   ;# Heps
    mulsd  xmm7, xmm13
    mulsd  xmm11, xmm14 
    mulsd  xmm2, xmm12   ;# Geps
    mulsd  xmm6, xmm13
    mulsd  xmm10, xmm14 
    mulsd  xmm3, xmm12   ;# Heps2
    mulsd  xmm7, xmm13
    mulsd  xmm11, xmm14 

    addsd  xmm1, xmm2   ;# F+Geps
    addsd  xmm5, xmm6
    addsd  xmm9, xmm10 
    addsd  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addsd  xmm5, xmm7
    addsd  xmm9, xmm11 
    addsd  xmm3, xmm3    ;# 2*Heps2
    addsd  xmm7, xmm7
    addsd  xmm11, xmm11
    addsd  xmm3, xmm2    ;# 2*Heps2+Geps
    addsd  xmm7, xmm6  
    addsd  xmm11, xmm10
    addsd  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addsd  xmm7, xmm5
    addsd  xmm11, xmm9
    mulsd  xmm1, xmm12   ;# eps*Fp
    mulsd  xmm5, xmm13
    mulsd  xmm9, xmm14
    movsd xmm12, [rsp + nb332_qqOH]
    movsd xmm13, [rsp + nb332_qqHH]
    addsd  xmm1, xmm0     ;# VV
    addsd  xmm5, xmm4
    addsd  xmm9, xmm8
    mulsd  xmm1, xmm12   ;# VV*qq = vcoul
    mulsd  xmm5, xmm13
    mulsd  xmm9, xmm13
    mulsd  xmm3, xmm12    ;# FF*qq = fij
    mulsd  xmm7, xmm13
    mulsd  xmm11, xmm13
    
    ;# accumulate vctot
    addsd  xmm1, [rsp + nb332_vctot]
    addsd  xmm5, xmm9
    addsd  xmm1, xmm5
    movsd [rsp + nb332_vctot], xmm1
    
    movsd xmm10, [rsp + nb332_tsc]
    mulsd  xmm3, xmm10  ;# fscal
    mulsd  xmm7, xmm10
    mulsd  xmm10, xmm11
    
    xorpd  xmm4, xmm4
    xorpd  xmm8, xmm8
    xorpd  xmm11, xmm11
    
    subsd  xmm4, xmm3
    subsd  xmm8, xmm7
    subsd  xmm11, xmm10

    mulsd  xmm4, [rsp + nb332_rinvOH1]
    mulsd  xmm8, [rsp + nb332_rinvH1H1]
    mulsd  xmm11, [rsp + nb332_rinvH2H1]
    
    ;# move j H1 forces to xmm0-xmm2
    mov rdi, [rbp + nb332_faction]
	movsd xmm0, [rdi + rax*8 + 24]
	movsd xmm1, [rdi + rax*8 + 32]
	movsd xmm2, [rdi + rax*8 + 40]

    movapd xmm3, xmm4
    movapd xmm5, xmm4
    movapd xmm7, xmm8
    movapd xmm9, xmm8
    movapd xmm10, xmm11
    movapd xmm12, xmm11

	mulsd xmm3, [rsp + nb332_dxOH1]
	mulsd xmm4, [rsp + nb332_dyOH1]
	mulsd xmm5, [rsp + nb332_dzOH1]
	mulsd xmm7, [rsp + nb332_dxH1H1]
	mulsd xmm8, [rsp + nb332_dyH1H1]
	mulsd xmm9, [rsp + nb332_dzH1H1]
	mulsd xmm10, [rsp + nb332_dxH2H1]
	mulsd xmm11, [rsp + nb332_dyH2H1]
	mulsd xmm12, [rsp + nb332_dzH2H1]

    addsd xmm0, xmm3
    addsd xmm1, xmm4
    addsd xmm2, xmm5
    addsd xmm3, [rsp + nb332_fixO]
    addsd xmm4, [rsp + nb332_fiyO]
    addsd xmm5, [rsp + nb332_fizO]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb332_fixH1]
    addsd xmm8, [rsp + nb332_fiyH1]
    addsd xmm9, [rsp + nb332_fizH1]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb332_fixH2]
    addsd xmm11, [rsp + nb332_fiyH2]
    addsd xmm12, [rsp + nb332_fizH2]

    movsd [rsp + nb332_fixO], xmm3
    movsd [rsp + nb332_fiyO], xmm4
    movsd [rsp + nb332_fizO], xmm5
    movsd [rsp + nb332_fixH1], xmm7
    movsd [rsp + nb332_fiyH1], xmm8
    movsd [rsp + nb332_fizH1], xmm9
    movsd [rsp + nb332_fixH2], xmm10
    movsd [rsp + nb332_fiyH2], xmm11
    movsd [rsp + nb332_fizH2], xmm12
   
    ;# store back j H1 forces from xmm0-xmm2
	movsd [rdi + rax*8 + 24], xmm0
	movsd [rdi + rax*8 + 32], xmm1
	movsd [rdi + rax*8 + 40], xmm2
       
	;# move j H2 coordinates to local temp variables 
    mov rsi, [rbp + nb332_pos]
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
    
    subsd xmm0, [rsp + nb332_ixO]
    subsd xmm1, [rsp + nb332_iyO]
    subsd xmm2, [rsp + nb332_izO]
    subsd xmm3, [rsp + nb332_ixH1]
    subsd xmm4, [rsp + nb332_iyH1]
    subsd xmm5, [rsp + nb332_izH1]
    subsd xmm6, [rsp + nb332_ixH2]
    subsd xmm7, [rsp + nb332_iyH2]
    subsd xmm8, [rsp + nb332_izH2]
    
	movsd [rsp + nb332_dxOH2], xmm0
	movsd [rsp + nb332_dyOH2], xmm1
	movsd [rsp + nb332_dzOH2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movsd [rsp + nb332_dxH1H2], xmm3
	movsd [rsp + nb332_dyH1H2], xmm4
	movsd [rsp + nb332_dzH1H2], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movsd [rsp + nb332_dxH2H2], xmm6
	movsd [rsp + nb332_dyH2H2], xmm7
	movsd [rsp + nb332_dzH2H2], xmm8
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
		
	movsd  xmm9, [rsp + nb332_three]
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

	movsd  xmm15, [rsp + nb332_half]
	mulsd   xmm9, xmm15  ;# first iteration for rinvOH2 
	mulsd   xmm10, xmm15 ;# first iteration for rinvH1H2
    mulsd   xmm11, xmm15 ;# first iteration for rinvH2H2

    ;# second iteration step    
	movsd  xmm2, xmm9
	movsd  xmm5, xmm10
    movsd  xmm8, xmm11
    
	mulsd   xmm2, xmm2 ;# lu*lu
	mulsd   xmm5, xmm5 ;# lu*lu
    mulsd   xmm8, xmm8 ;# lu*lu
		
	movsd  xmm1, [rsp + nb332_three]
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

	movsd  xmm15, [rsp + nb332_half]
	mulsd   xmm9, xmm15  ;#  rinvOH2
	mulsd   xmm10, xmm15 ;#   rinvH1H2
    mulsd   xmm11, xmm15 ;#   rinvH2H2
	
	movsd  [rsp + nb332_rinvOH2], xmm9
	movsd  [rsp + nb332_rinvH1H2], xmm10
	movsd  [rsp + nb332_rinvH2H2], xmm11
	
	;# H2 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movsd  xmm1, [rsp + nb332_tsc]
    mulsd  xmm0, xmm9  ;# r
    mulsd  xmm3, xmm10
    mulsd  xmm6, xmm11
    mulsd  xmm0, xmm1 ;# rtab
    mulsd  xmm3, xmm1
    mulsd  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttsd2si r8d, xmm0
    cvttsd2si r10d, xmm3
    cvttsd2si r12d, xmm6        

    ;# convert back to float
    cvtsi2sd  xmm2, r8d
    cvtsi2sd  xmm5, r10d
    cvtsi2sd  xmm8, r12d
    
    ;# multiply by 4
    shl   r8, 2
    shl   r10d, 2
    shl   r12d, 2

	;# multiply by 3
    lea   r8, [r8 + r8*2]        
    lea   r10, [r10 + r10*2]        
    lea   r12, [r12 + r12*2]        

    mov  rsi, [rbp + nb332_VFtab]

    ;# calculate eps
    subsd     xmm0, xmm2
    subsd     xmm3, xmm5
    subsd     xmm6, xmm8

    movsd    [rsp + nb332_epsO], xmm0
    movsd    [rsp + nb332_epsH1], xmm3
    movsd    [rsp + nb332_epsH2], xmm6

    ;# Load LOTS of table data
    movsd xmm0,  [rsi + r8*8]
    movsd xmm1,  [rsi + r8*8 + 8]
    movsd xmm2,  [rsi + r8*8 + 16]
    movsd xmm3,  [rsi + r8*8 + 24]
    movsd xmm4,  [rsi + r10*8]
    movsd xmm5,  [rsi + r10*8 + 8]
    movsd xmm6,  [rsi + r10*8 + 16]
    movsd xmm7,  [rsi + r10*8 + 24]
    movsd xmm8,  [rsi + r12*8]
    movsd xmm9,  [rsi + r12*8 + 8]
    movsd xmm10, [rsi + r12*8 + 16]
    movsd xmm11, [rsi + r12*8 + 24]
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movsd xmm12, [rsp + nb332_epsO]
    movsd xmm13, [rsp + nb332_epsH1]
    movsd xmm14, [rsp + nb332_epsH2]
    
    mulsd  xmm3, xmm12   ;# Heps
    mulsd  xmm7, xmm13
    mulsd  xmm11, xmm14 
    mulsd  xmm2, xmm12   ;# Geps
    mulsd  xmm6, xmm13
    mulsd  xmm10, xmm14 
    mulsd  xmm3, xmm12   ;# Heps2
    mulsd  xmm7, xmm13
    mulsd  xmm11, xmm14 

    addsd  xmm1, xmm2   ;# F+Geps
    addsd  xmm5, xmm6
    addsd  xmm9, xmm10 
    addsd  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addsd  xmm5, xmm7
    addsd  xmm9, xmm11 
    addsd  xmm3, xmm3    ;# 2*Heps2
    addsd  xmm7, xmm7
    addsd  xmm11, xmm11
    addsd  xmm3, xmm2    ;# 2*Heps2+Geps
    addsd  xmm7, xmm6  
    addsd  xmm11, xmm10
    addsd  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addsd  xmm7, xmm5
    addsd  xmm11, xmm9
    mulsd  xmm1, xmm12   ;# eps*Fp
    mulsd  xmm5, xmm13
    mulsd  xmm9, xmm14
    movsd xmm12, [rsp + nb332_qqOH]
    movsd xmm13, [rsp + nb332_qqHH]
    addsd  xmm1, xmm0     ;# VV
    addsd  xmm5, xmm4
    addsd  xmm9, xmm8
    mulsd  xmm1, xmm12   ;# VV*qq = vcoul
    mulsd  xmm5, xmm13
    mulsd  xmm9, xmm13
    mulsd  xmm3, xmm12    ;# FF*qq = fij
    mulsd  xmm7, xmm13
    mulsd  xmm11, xmm13
    
    ;# accumulate vctot
    addsd  xmm1, [rsp + nb332_vctot]
    addsd  xmm5, xmm9
    addsd  xmm1, xmm5
    movsd [rsp + nb332_vctot], xmm1
    
    movsd xmm10, [rsp + nb332_tsc]
    mulsd  xmm3, xmm10  ;# fscal
    mulsd  xmm7, xmm10
    mulsd  xmm10, xmm11
    
    xorpd  xmm4, xmm4
    xorpd  xmm8, xmm8
    xorpd  xmm11, xmm11
    
    subsd  xmm4, xmm3
    subsd  xmm8, xmm7
    subsd  xmm11, xmm10

    mulsd  xmm4, [rsp + nb332_rinvOH2]
    mulsd  xmm8, [rsp + nb332_rinvH1H2]
    mulsd  xmm11, [rsp + nb332_rinvH2H2]
    
    ;# move j H2 forces to xmm0-xmm2
    mov rdi, [rbp + nb332_faction]
	movsd xmm0, [rdi + rax*8 + 48]
	movsd xmm1, [rdi + rax*8 + 56]
	movsd xmm2, [rdi + rax*8 + 64]

    movapd xmm3, xmm4
    movapd xmm5, xmm4
    movapd xmm7, xmm8
    movapd xmm9, xmm8
    movapd xmm10, xmm11
    movapd xmm12, xmm11

	mulsd xmm3, [rsp + nb332_dxOH2]
	mulsd xmm4, [rsp + nb332_dyOH2]
	mulsd xmm5, [rsp + nb332_dzOH2]
	mulsd xmm7, [rsp + nb332_dxH1H2]
	mulsd xmm8, [rsp + nb332_dyH1H2]
	mulsd xmm9, [rsp + nb332_dzH1H2]
	mulsd xmm10, [rsp + nb332_dxH2H2]
	mulsd xmm11, [rsp + nb332_dyH2H2]
	mulsd xmm12, [rsp + nb332_dzH2H2]

    addsd xmm0, xmm3
    addsd xmm1, xmm4
    addsd xmm2, xmm5
    addsd xmm3, [rsp + nb332_fixO]
    addsd xmm4, [rsp + nb332_fiyO]
    addsd xmm5, [rsp + nb332_fizO]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb332_fixH1]
    addsd xmm8, [rsp + nb332_fiyH1]
    addsd xmm9, [rsp + nb332_fizH1]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb332_fixH2]
    addsd xmm11, [rsp + nb332_fiyH2]
    addsd xmm12, [rsp + nb332_fizH2]

    movsd [rsp + nb332_fixO], xmm3
    movsd [rsp + nb332_fiyO], xmm4
    movsd [rsp + nb332_fizO], xmm5
    movsd [rsp + nb332_fixH1], xmm7
    movsd [rsp + nb332_fiyH1], xmm8
    movsd [rsp + nb332_fizH1], xmm9
    movsd [rsp + nb332_fixH2], xmm10
    movsd [rsp + nb332_fiyH2], xmm11
    movsd [rsp + nb332_fizH2], xmm12
   
    ;# store back j H2 forces from xmm0-xmm2
	movsd [rdi + rax*8 + 48], xmm0
	movsd [rdi + rax*8 + 56], xmm1
	movsd [rdi + rax*8 + 64], xmm2
	
.nb332_updateouterdata:
	mov   ecx, [rsp + nb332_ii3]
	mov   rdi, [rbp + nb332_faction]
	mov   rsi, [rbp + nb332_fshift]
	mov   edx, [rsp + nb332_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb332_fixO]
	movapd xmm1, [rsp + nb332_fiyO]
	movapd xmm2, [rsp + nb332_fizO]

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
	movapd xmm0, [rsp + nb332_fixH1]
	movapd xmm1, [rsp + nb332_fiyH1]
	movapd xmm2, [rsp + nb332_fizH1]

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
	movapd xmm0, [rsp + nb332_fixH2]
	movapd xmm1, [rsp + nb332_fiyH2]
	movapd xmm2, [rsp + nb332_fizH2]

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
	mov esi, [rsp + nb332_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb332_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb332_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb332_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb332_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb332_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb332_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb332_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb332_n], esi
        jmp .nb332_outer
.nb332_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb332_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb332_end
        ;# non-zero, do one more workunit
        jmp   .nb332_threadloop
.nb332_end:
	mov eax, [rsp + nb332_nouter]
	mov ebx, [rsp + nb332_ninner]
	mov rcx, [rbp + nb332_outeriter]
	mov rdx, [rbp + nb332_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1584
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




.globl nb_kernel332nf_x86_64_sse2
.globl _nb_kernel332nf_x86_64_sse2
nb_kernel332nf_x86_64_sse2:	
_nb_kernel332nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb332nf_fshift,         16
.equiv          nb332nf_gid,            24
.equiv          nb332nf_pos,            32
.equiv          nb332nf_faction,        40
.equiv          nb332nf_charge,         48
.equiv          nb332nf_p_facel,        56
.equiv          nb332nf_argkrf,         64
.equiv          nb332nf_argcrf,         72
.equiv          nb332nf_Vc,             80
.equiv          nb332nf_type,           88
.equiv          nb332nf_p_ntype,        96
.equiv          nb332nf_vdwparam,       104
.equiv          nb332nf_Vvdw,           112
.equiv          nb332nf_p_tabscale,     120
.equiv          nb332nf_VFtab,          128
.equiv          nb332nf_invsqrta,       136
.equiv          nb332nf_dvda,           144
.equiv          nb332nf_p_gbtabscale,   152
.equiv          nb332nf_GBtab,          160
.equiv          nb332nf_p_nthreads,     168
.equiv          nb332nf_count,          176
.equiv          nb332nf_mtx,            184
.equiv          nb332nf_outeriter,      192
.equiv          nb332nf_inneriter,      200
.equiv          nb332nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb332nf_ixO,            0
.equiv          nb332nf_iyO,            16
.equiv          nb332nf_izO,            32
.equiv          nb332nf_ixH1,           48
.equiv          nb332nf_iyH1,           64
.equiv          nb332nf_izH1,           80
.equiv          nb332nf_ixH2,           96
.equiv          nb332nf_iyH2,           112
.equiv          nb332nf_izH2,           128
.equiv          nb332nf_jxO,            144
.equiv          nb332nf_jyO,            160
.equiv          nb332nf_jzO,            176
.equiv          nb332nf_jxH1,           192
.equiv          nb332nf_jyH1,           208
.equiv          nb332nf_jzH1,           224
.equiv          nb332nf_jxH2,           240
.equiv          nb332nf_jyH2,           256
.equiv          nb332nf_jzH2,           272
.equiv          nb332nf_qqOO,           288
.equiv          nb332nf_qqOH,           304
.equiv          nb332nf_qqHH,           320
.equiv          nb332nf_tsc,            336
.equiv          nb332nf_c6,             352
.equiv          nb332nf_c12,            368
.equiv          nb332nf_vctot,          384
.equiv          nb332nf_Vvdwtot,        400
.equiv          nb332nf_half,           416
.equiv          nb332nf_three,          432
.equiv          nb332nf_rsqOO,          448
.equiv          nb332nf_rsqOH1,         464
.equiv          nb332nf_rsqOH2,         480
.equiv          nb332nf_rsqH1O,         496
.equiv          nb332nf_rsqH1H1,        512
.equiv          nb332nf_rsqH1H2,        528
.equiv          nb332nf_rsqH2O,         544
.equiv          nb332nf_rsqH2H1,        560
.equiv          nb332nf_rsqH2H2,        576
.equiv          nb332nf_rinvOO,         592
.equiv          nb332nf_rinvOH1,        608
.equiv          nb332nf_rinvOH2,        624
.equiv          nb332nf_rinvH1O,        640
.equiv          nb332nf_rinvH1H1,       656
.equiv          nb332nf_rinvH1H2,       672
.equiv          nb332nf_rinvH2O,        688
.equiv          nb332nf_rinvH2H1,       704
.equiv          nb332nf_rinvH2H2,       720
.equiv          nb332nf_is3,            736
.equiv          nb332nf_ii3,            740
.equiv          nb332nf_nri,            744
.equiv          nb332nf_iinr,           752
.equiv          nb332nf_jindex,         760
.equiv          nb332nf_jjnr,           768
.equiv          nb332nf_shift,          776
.equiv          nb332nf_shiftvec,       784
.equiv          nb332nf_facel,          792
.equiv          nb332nf_innerjjnr,      800
.equiv          nb332nf_innerk,         808
.equiv          nb332nf_n,              812
.equiv          nb332nf_nn1,            816
.equiv          nb332nf_nouter,         820
.equiv          nb332nf_ninner,         824

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
	sub rsp, 832		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb332nf_nouter], eax
	mov [rsp + nb332nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb332nf_nri], edi
	mov [rsp + nb332nf_iinr], rsi
	mov [rsp + nb332nf_jindex], rdx
	mov [rsp + nb332nf_jjnr], rcx
	mov [rsp + nb332nf_shift], r8
	mov [rsp + nb332nf_shiftvec], r9
	mov rsi, [rbp + nb332nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb332nf_facel], xmm0

	mov rax, [rbp + nb332nf_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb332nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb332nf_half], eax
	mov [rsp + nb332nf_half+4], ebx
	movsd xmm1, [rsp + nb332nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb332nf_half], xmm1
	movapd [rsp + nb332nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb332nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb332nf_charge]
	movsd xmm3, [rdx + rbx*8]	
	movsd xmm4, xmm3	
	movsd xmm5, [rdx + rbx*8 + 8]	
	mov rsi, [rbp + nb332nf_p_facel]
	movsd xmm0, [rsi]
	movsd xmm6, [rsp + nb332nf_facel]
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm5
	mulsd  xmm5, xmm5
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb332nf_qqOO], xmm3
	movapd [rsp + nb332nf_qqOH], xmm4
	movapd [rsp + nb332nf_qqHH], xmm5
	
	xorpd xmm0, xmm0
	mov   rdx, [rbp + nb332nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb332nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   rax, [rbp + nb332nf_vdwparam]
	movlpd xmm0, [rax + rdx*8]
	movlpd xmm1, [rax + rdx*8 + 8]
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	movapd [rsp + nb332nf_c6], xmm0
	movapd [rsp + nb332nf_c12], xmm1

.nb332nf_threadloop:
        mov   rsi, [rbp + nb332nf_count]        ;# pointer to sync counter
        mov   eax, [rsi]
.nb332nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb332nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb332nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb332nf_n], eax
        mov [rsp + nb332nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb332nf_outerstart
        jmp .nb332nf_end

.nb332nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb332nf_nouter]
	mov [rsp + nb332nf_nouter], ebx

.nb332nf_outer:
	mov   rax, [rsp + nb332nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 

	mov   rax, [rsp + nb332nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb332nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb332nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb332nf_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [rax + rbx*8]
	addsd xmm4, [rax + rbx*8 + 8]
	addsd xmm5, [rax + rbx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb332nf_ixO], xmm3
	movapd [rsp + nb332nf_iyO], xmm4
	movapd [rsp + nb332nf_izO], xmm5

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
	movapd [rsp + nb332nf_ixH1], xmm0
	movapd [rsp + nb332nf_iyH1], xmm1
	movapd [rsp + nb332nf_izH1], xmm2
	movapd [rsp + nb332nf_ixH2], xmm3
	movapd [rsp + nb332nf_iyH2], xmm4
	movapd [rsp + nb332nf_izH2], xmm5

	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb332nf_vctot], xmm4
	movapd [rsp + nb332nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb332nf_jindex]
	mov   ecx, [rax+rsi*4]	     ;# jindex[n] 
	mov   edx, [rax+rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb332nf_pos] 
	mov   rax, [rsp + nb332nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb332nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb332nf_ninner]
	mov   [rsp + nb332nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb332nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb332nf_unroll_loop
	jmp   .nb332nf_checksingle
.nb332nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb332nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	
	add qword ptr [rsp + nb332nf_innerjjnr], 8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb332nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	
	;# move j coordinates to local temp variables 
	movlpd xmm2, [rsi + rax*8]
	movlpd xmm3, [rsi + rax*8 + 8]
	movlpd xmm4, [rsi + rax*8 + 16]
	movlpd xmm5, [rsi + rax*8 + 24]
	movlpd xmm6, [rsi + rax*8 + 32]
	movlpd xmm7, [rsi + rax*8 + 40]
	movhpd xmm2, [rsi + rbx*8]
	movhpd xmm3, [rsi + rbx*8 + 8]
	movhpd xmm4, [rsi + rbx*8 + 16]
	movhpd xmm5, [rsi + rbx*8 + 24]
	movhpd xmm6, [rsi + rbx*8 + 32]
	movhpd xmm7, [rsi + rbx*8 + 40]
	movapd 	[rsp + nb332nf_jxO], xmm2
	movapd 	[rsp + nb332nf_jyO], xmm3
	movapd 	[rsp + nb332nf_jzO], xmm4
	movapd 	[rsp + nb332nf_jxH1], xmm5
	movapd 	[rsp + nb332nf_jyH1], xmm6
	movapd 	[rsp + nb332nf_jzH1], xmm7
	movlpd xmm2, [rsi + rax*8 + 48]
	movlpd xmm3, [rsi + rax*8 + 56]
	movlpd xmm4, [rsi + rax*8 + 64]
	movhpd xmm2, [rsi + rbx*8 + 48]
	movhpd xmm3, [rsi + rbx*8 + 56]
	movhpd xmm4, [rsi + rbx*8 + 64]
	movapd 	[rsp + nb332nf_jxH2], xmm2
	movapd 	[rsp + nb332nf_jyH2], xmm3
	movapd 	[rsp + nb332nf_jzH2], xmm4
	
	movapd xmm0, [rsp + nb332nf_ixO]
	movapd xmm1, [rsp + nb332nf_iyO]
	movapd xmm2, [rsp + nb332nf_izO]
	movapd xmm3, [rsp + nb332nf_ixO]
	movapd xmm4, [rsp + nb332nf_iyO]
	movapd xmm5, [rsp + nb332nf_izO]
	subpd  xmm0, [rsp + nb332nf_jxO]
	subpd  xmm1, [rsp + nb332nf_jyO]
	subpd  xmm2, [rsp + nb332nf_jzO]
	subpd  xmm3, [rsp + nb332nf_jxH1]
	subpd  xmm4, [rsp + nb332nf_jyH1]
	subpd  xmm5, [rsp + nb332nf_jzH1]
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
	movapd [rsp + nb332nf_rsqOO], xmm0
	movapd [rsp + nb332nf_rsqOH1], xmm3

	movapd xmm0, [rsp + nb332nf_ixO]
	movapd xmm1, [rsp + nb332nf_iyO]
	movapd xmm2, [rsp + nb332nf_izO]
	movapd xmm3, [rsp + nb332nf_ixH1]
	movapd xmm4, [rsp + nb332nf_iyH1]
	movapd xmm5, [rsp + nb332nf_izH1]
	subpd  xmm0, [rsp + nb332nf_jxH2]
	subpd  xmm1, [rsp + nb332nf_jyH2]
	subpd  xmm2, [rsp + nb332nf_jzH2]
	subpd  xmm3, [rsp + nb332nf_jxO]
	subpd  xmm4, [rsp + nb332nf_jyO]
	subpd  xmm5, [rsp + nb332nf_jzO]
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
	movapd [rsp + nb332nf_rsqOH2], xmm0
	movapd [rsp + nb332nf_rsqH1O], xmm3

	movapd xmm0, [rsp + nb332nf_ixH1]
	movapd xmm1, [rsp + nb332nf_iyH1]
	movapd xmm2, [rsp + nb332nf_izH1]
	movapd xmm3, [rsp + nb332nf_ixH1]
	movapd xmm4, [rsp + nb332nf_iyH1]
	movapd xmm5, [rsp + nb332nf_izH1]
	subpd  xmm0, [rsp + nb332nf_jxH1]
	subpd  xmm1, [rsp + nb332nf_jyH1]
	subpd  xmm2, [rsp + nb332nf_jzH1]
	subpd  xmm3, [rsp + nb332nf_jxH2]
	subpd  xmm4, [rsp + nb332nf_jyH2]
	subpd  xmm5, [rsp + nb332nf_jzH2]
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
	movapd [rsp + nb332nf_rsqH1H1], xmm0
	movapd [rsp + nb332nf_rsqH1H2], xmm3

	movapd xmm0, [rsp + nb332nf_ixH2]
	movapd xmm1, [rsp + nb332nf_iyH2]
	movapd xmm2, [rsp + nb332nf_izH2]
	movapd xmm3, [rsp + nb332nf_ixH2]
	movapd xmm4, [rsp + nb332nf_iyH2]
	movapd xmm5, [rsp + nb332nf_izH2]
	subpd  xmm0, [rsp + nb332nf_jxO]
	subpd  xmm1, [rsp + nb332nf_jyO]
	subpd  xmm2, [rsp + nb332nf_jzO]
	subpd  xmm3, [rsp + nb332nf_jxH1]
	subpd  xmm4, [rsp + nb332nf_jyH1]
	subpd  xmm5, [rsp + nb332nf_jzH1]
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
	movapd [rsp + nb332nf_rsqH2O], xmm0
	movapd [rsp + nb332nf_rsqH2H1], xmm4

	movapd xmm0, [rsp + nb332nf_ixH2]
	movapd xmm1, [rsp + nb332nf_iyH2]
	movapd xmm2, [rsp + nb332nf_izH2]
	subpd  xmm0, [rsp + nb332nf_jxH2]
	subpd  xmm1, [rsp + nb332nf_jyH2]
	subpd  xmm2, [rsp + nb332nf_jzH2]
	mulpd xmm0, xmm0
	mulpd xmm1, xmm1
	mulpd xmm2, xmm2
	addpd xmm0, xmm1
	addpd xmm0, xmm2
	movapd [rsp + nb332nf_rsqH2H2], xmm0
	
	;# start doing invsqrt use rsq values in xmm0, xmm4 
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
	movapd  xmm3, [rsp + nb332nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb332nf_half] ;# iter1 
	mulpd   xmm7, [rsp + nb332nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb332nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb332nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb332nf_half] ;# rinv 
	movapd [rsp + nb332nf_rinvH2H2], xmm1
	movapd [rsp + nb332nf_rinvH2H1], xmm5

	movapd xmm0, [rsp + nb332nf_rsqOO]
	movapd xmm4, [rsp + nb332nf_rsqOH1]	
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
	movapd  xmm3, [rsp + nb332nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb332nf_half] ;# iter1 of  
	mulpd   xmm7, [rsp + nb332nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb332nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb332nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb332nf_half] ;# rinv
	movapd [rsp + nb332nf_rinvOO], xmm1
	movapd [rsp + nb332nf_rinvOH1], xmm5

	movapd xmm0, [rsp + nb332nf_rsqOH2]
	movapd xmm4, [rsp + nb332nf_rsqH1O]	
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
	movapd  xmm3, [rsp + nb332nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb332nf_half] ;# iter1 
	mulpd   xmm7, [rsp + nb332nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb332nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb332nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb332nf_half] ;# rinv 
	movapd [rsp + nb332nf_rinvOH2], xmm1
	movapd [rsp + nb332nf_rinvH1O], xmm5

	movapd xmm0, [rsp + nb332nf_rsqH1H1]
	movapd xmm4, [rsp + nb332nf_rsqH1H2]	
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
	movapd  xmm3, [rsp + nb332nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb332nf_half] ;# iter1a 
	mulpd   xmm7, [rsp + nb332nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb332nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb332nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb332nf_half] ;# rinv 
	movapd [rsp + nb332nf_rinvH1H1], xmm1
	movapd [rsp + nb332nf_rinvH1H2], xmm5

	movapd xmm0, [rsp + nb332nf_rsqH2O]
	cvtpd2ps xmm1, xmm0	
	rsqrtps xmm1, xmm1
	cvtps2pd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulpd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [rsp + nb332nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm3, [rsp + nb332nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulpd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [rsp + nb332nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm1, [rsp + nb332nf_half] ;# rinv 
	movapd [rsp + nb332nf_rinvH2O], xmm1
	
	;# start with OO interaction 
	movapd xmm0, [rsp + nb332nf_rinvOO]
	movapd xmm1, xmm0
	mulpd  xmm1, [rsp + nb332nf_rsqOO] ;# xmm1=r 
	mulpd  xmm1, [rsp + nb332nf_tsc]

	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb332nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 
	lea   rbx, [rbx + rbx*2]	;# idx*=3 (12 total now) 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb332nf_qqOO]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV 
    ;# at this point mm5 contains vcoul 
    ;# increment vcoul - then we can get rid of mm5 
    ;# update vctot 
    addpd  xmm5, [rsp + nb332nf_vctot]
    movapd [rsp + nb332nf_vctot], xmm5

	;# Dispersion 
	movapd xmm4, [rsi + rax*8 + 32]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8 + 32]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 48]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 48]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	mulpd  xmm5, [rsp + nb332nf_c6] ;# Vvdw6 

	addpd  xmm5, [rsp + nb332nf_Vvdwtot]
	movapd [rsp + nb332nf_Vvdwtot], xmm5

	;# Repulsion 
	movapd xmm4, [rsi + rax*8 + 64]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8 + 64]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 80]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 80]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	mulpd  xmm5, [rsp + nb332nf_c12] ;# Vvdw12 

	addpd  xmm5, [rsp + nb332nf_Vvdwtot]
	movapd [rsp + nb332nf_Vvdwtot], xmm5

	;# O-H1 interaction 
	movapd xmm0, [rsp + nb332nf_rinvOH1]
	movapd xmm1, xmm0
	mulpd  xmm1, [rsp + nb332nf_rsqOH1] ;# xmm1=r 
	mulpd  xmm1, [rsp + nb332nf_tsc]

	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb332nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 
	lea   rbx, [rbx + rbx*2]	;# idx*=3 (12 total now) 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb332nf_qqOH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [rsp + nb332nf_vctot]
    movapd [rsp + nb332nf_vctot], xmm5

	;# O-H2 interaction  
	movapd xmm0, [rsp + nb332nf_rinvOH2]
	movapd xmm1, xmm0
	mulpd  xmm1, [rsp + nb332nf_rsqOH2] ;# xmm1=r 
	mulpd  xmm1, [rsp + nb332nf_tsc]
	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb332nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 
	lea   rbx, [rbx + rbx*2]	;# idx*=3 (12 total now) 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb332nf_qqOH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [rsp + nb332nf_vctot]
    movapd [rsp + nb332nf_vctot], xmm5

	;# H1-O interaction 
	movapd xmm0, [rsp + nb332nf_rinvH1O]
	movapd xmm1, xmm0
	mulpd  xmm1, [rsp + nb332nf_rsqH1O] ;# xmm1=r 
	mulpd  xmm1, [rsp + nb332nf_tsc]
	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb332nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 
	lea   rbx, [rbx + rbx*2]	;# idx*=3 (12 total now) 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb332nf_qqOH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV 
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [rsp + nb332nf_vctot]
    movapd [rsp + nb332nf_vctot], xmm5

	;# H1-H1 interaction 
	movapd xmm0, [rsp + nb332nf_rinvH1H1]
	movapd xmm1, xmm0
	mulpd  xmm1, [rsp + nb332nf_rsqH1H1] ;# xmm1=r 
	mulpd  xmm1, [rsp + nb332nf_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb332nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 
	lea   rbx, [rbx + rbx*2]	;# idx*=3 (12 total now) 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb332nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [rsp + nb332nf_vctot]
    movapd [rsp + nb332nf_vctot], xmm5

	;# H1-H2 interaction 
	movapd xmm0, [rsp + nb332nf_rinvH1H2]
	movapd xmm1, xmm0
	mulpd  xmm1, [rsp + nb332nf_rsqH1H2] ;# xmm1=r 
	mulpd  xmm1, [rsp + nb332nf_tsc]
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb332nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 
	lea   rbx, [rbx + rbx*2]	;# idx*=3 (12 total now) 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb332nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV 
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [rsp + nb332nf_vctot]
    movapd [rsp + nb332nf_vctot], xmm5

	;# H2-O interaction 
	movapd xmm0, [rsp + nb332nf_rinvH2O]
	movapd xmm1, xmm0
	mulpd  xmm1, [rsp + nb332nf_rsqH2O] ;# xmm1=r 
	mulpd  xmm1, [rsp + nb332nf_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb332nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 
	lea   rbx, [rbx + rbx*2]	;# idx*=3 (12 total now) 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb332nf_qqOH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [rsp + nb332nf_vctot]
    movapd [rsp + nb332nf_vctot], xmm5
	
	;# H2-H1 interaction 
	movapd xmm0, [rsp + nb332nf_rinvH2H1]
	movapd xmm1, xmm0
	mulpd  xmm1, [rsp + nb332nf_rsqH2H1] ;# xmm1=r 
	mulpd  xmm1, [rsp + nb332nf_tsc]
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb332nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 
	lea   rbx, [rbx + rbx*2]	;# idx*=3 (12 total now) 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb332nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addpd  xmm5, [rsp + nb332nf_vctot]
    movapd [rsp + nb332nf_vctot], xmm5

	;# H2-H2 interaction 
	movapd xmm0, [rsp + nb332nf_rinvH2H2]
	movapd xmm1, xmm0
	mulpd  xmm1, [rsp + nb332nf_rsqH2H2] ;# xmm1=r 
	mulpd  xmm1, [rsp + nb332nf_tsc]	
	cvttpd2pi mm6, xmm1	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb332nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 
	lea   rbx, [rbx + rbx*2]	;# idx*=3 (12 total now) 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb332nf_qqHH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
	
    addpd  xmm5, [rsp + nb332nf_vctot]
    movapd [rsp + nb332nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb332nf_innerk],  2
	jl    .nb332nf_checksingle
	jmp   .nb332nf_unroll_loop
.nb332nf_checksingle:
	mov   edx, [rsp + nb332nf_innerk]
	and   edx, 1
	jnz   .nb332nf_dosingle
	jmp   .nb332nf_updateouterdata
.nb332nf_dosingle:
	mov   rdx, [rsp + nb332nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]

	mov rsi, [rbp + nb332nf_pos]
	lea   rax, [rax + rax*2]  

	;# fetch j coordinates 
	movlpd xmm2, [rsi + rax*8]
	movlpd xmm3, [rsi + rax*8 + 8]
	movlpd xmm4, [rsi + rax*8 + 16]
	movlpd xmm5, [rsi + rax*8 + 24]
	movlpd xmm6, [rsi + rax*8 + 32]
	movlpd xmm7, [rsi + rax*8 + 40]
	movapd 	[rsp + nb332nf_jxO], xmm2
	movapd 	[rsp + nb332nf_jyO], xmm3
	movapd 	[rsp + nb332nf_jzO], xmm4
	movapd 	[rsp + nb332nf_jxH1], xmm5
	movapd 	[rsp + nb332nf_jyH1], xmm6
	movapd 	[rsp + nb332nf_jzH1], xmm7
	movlpd xmm2, [rsi + rax*8 + 48]
	movlpd xmm3, [rsi + rax*8 + 56]
	movlpd xmm4, [rsi + rax*8 + 64]
	movapd 	[rsp + nb332nf_jxH2], xmm2
	movapd 	[rsp + nb332nf_jyH2], xmm3
	movapd 	[rsp + nb332nf_jzH2], xmm4
	
	movapd xmm0, [rsp + nb332nf_ixO]
	movapd xmm1, [rsp + nb332nf_iyO]
	movapd xmm2, [rsp + nb332nf_izO]
	movapd xmm3, [rsp + nb332nf_ixO]
	movapd xmm4, [rsp + nb332nf_iyO]
	movapd xmm5, [rsp + nb332nf_izO]
	subsd  xmm0, [rsp + nb332nf_jxO]
	subsd  xmm1, [rsp + nb332nf_jyO]
	subsd  xmm2, [rsp + nb332nf_jzO]
	subsd  xmm3, [rsp + nb332nf_jxH1]
	subsd  xmm4, [rsp + nb332nf_jyH1]
	subsd  xmm5, [rsp + nb332nf_jzH1]
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
	movapd [rsp + nb332nf_rsqOO], xmm0
	movapd [rsp + nb332nf_rsqOH1], xmm3

	movapd xmm0, [rsp + nb332nf_ixO]
	movapd xmm1, [rsp + nb332nf_iyO]
	movapd xmm2, [rsp + nb332nf_izO]
	movapd xmm3, [rsp + nb332nf_ixH1]
	movapd xmm4, [rsp + nb332nf_iyH1]
	movapd xmm5, [rsp + nb332nf_izH1]
	subsd  xmm0, [rsp + nb332nf_jxH2]
	subsd  xmm1, [rsp + nb332nf_jyH2]
	subsd  xmm2, [rsp + nb332nf_jzH2]
	subsd  xmm3, [rsp + nb332nf_jxO]
	subsd  xmm4, [rsp + nb332nf_jyO]
	subsd  xmm5, [rsp + nb332nf_jzO]
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
	movapd [rsp + nb332nf_rsqOH2], xmm0
	movapd [rsp + nb332nf_rsqH1O], xmm3

	movapd xmm0, [rsp + nb332nf_ixH1]
	movapd xmm1, [rsp + nb332nf_iyH1]
	movapd xmm2, [rsp + nb332nf_izH1]
	movapd xmm3, [rsp + nb332nf_ixH1]
	movapd xmm4, [rsp + nb332nf_iyH1]
	movapd xmm5, [rsp + nb332nf_izH1]
	subsd  xmm0, [rsp + nb332nf_jxH1]
	subsd  xmm1, [rsp + nb332nf_jyH1]
	subsd  xmm2, [rsp + nb332nf_jzH1]
	subsd  xmm3, [rsp + nb332nf_jxH2]
	subsd  xmm4, [rsp + nb332nf_jyH2]
	subsd  xmm5, [rsp + nb332nf_jzH2]
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
	movapd [rsp + nb332nf_rsqH1H1], xmm0
	movapd [rsp + nb332nf_rsqH1H2], xmm3

	movapd xmm0, [rsp + nb332nf_ixH2]
	movapd xmm1, [rsp + nb332nf_iyH2]
	movapd xmm2, [rsp + nb332nf_izH2]
	movapd xmm3, [rsp + nb332nf_ixH2]
	movapd xmm4, [rsp + nb332nf_iyH2]
	movapd xmm5, [rsp + nb332nf_izH2]
	subsd  xmm0, [rsp + nb332nf_jxO]
	subsd  xmm1, [rsp + nb332nf_jyO]
	subsd  xmm2, [rsp + nb332nf_jzO]
	subsd  xmm3, [rsp + nb332nf_jxH1]
	subsd  xmm4, [rsp + nb332nf_jyH1]
	subsd  xmm5, [rsp + nb332nf_jzH1]
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
	movapd [rsp + nb332nf_rsqH2O], xmm0
	movapd [rsp + nb332nf_rsqH2H1], xmm4

	movapd xmm0, [rsp + nb332nf_ixH2]
	movapd xmm1, [rsp + nb332nf_iyH2]
	movapd xmm2, [rsp + nb332nf_izH2]
	subsd  xmm0, [rsp + nb332nf_jxH2]
	subsd  xmm1, [rsp + nb332nf_jyH2]
	subsd  xmm2, [rsp + nb332nf_jzH2]
	mulsd xmm0, xmm0
	mulsd xmm1, xmm1
	mulsd xmm2, xmm2
	addsd xmm0, xmm1
	addsd xmm0, xmm2
	movapd [rsp + nb332nf_rsqH2H2], xmm0
	
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
	movapd  xmm3, [rsp + nb332nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb332nf_half] ;# iter1 
	mulsd   xmm7, [rsp + nb332nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb332nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb332nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb332nf_half] ;# rinv 
	movapd [rsp + nb332nf_rinvH2H2], xmm1
	movapd [rsp + nb332nf_rinvH2H1], xmm5

	movapd xmm0, [rsp + nb332nf_rsqOO]
	movapd xmm4, [rsp + nb332nf_rsqOH1]	
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
	movapd  xmm3, [rsp + nb332nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb332nf_half] ;# iter1 of  
	mulsd   xmm7, [rsp + nb332nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb332nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb332nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb332nf_half] ;# rinv
	movapd [rsp + nb332nf_rinvOO], xmm1
	movapd [rsp + nb332nf_rinvOH1], xmm5

	movapd xmm0, [rsp + nb332nf_rsqOH2]
	movapd xmm4, [rsp + nb332nf_rsqH1O]	
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
	movapd  xmm3, [rsp + nb332nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb332nf_half] ;# iter1 
	mulsd   xmm7, [rsp + nb332nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb332nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb332nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb332nf_half] ;# rinv 
	movapd [rsp + nb332nf_rinvOH2], xmm1
	movapd [rsp + nb332nf_rinvH1O], xmm5

	movapd xmm0, [rsp + nb332nf_rsqH1H1]
	movapd xmm4, [rsp + nb332nf_rsqH1H2]	
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
	movapd  xmm3, [rsp + nb332nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb332nf_half] ;# iter1a 
	mulsd   xmm7, [rsp + nb332nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb332nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb332nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb332nf_half] ;# rinv 
	movapd [rsp + nb332nf_rinvH1H1], xmm1
	movapd [rsp + nb332nf_rinvH1H2], xmm5

	movapd xmm0, [rsp + nb332nf_rsqH2O]
	cvtsd2ss xmm1, xmm0	
	rsqrtss xmm1, xmm1
	cvtss2sd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulsd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [rsp + nb332nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm3, [rsp + nb332nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulsd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [rsp + nb332nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm1, [rsp + nb332nf_half] ;# rinv 
	movapd [rsp + nb332nf_rinvH2O], xmm1
	
	;# start with OO interaction 
	movapd xmm0, [rsp + nb332nf_rinvOO]
	movapd xmm1, xmm0
	mulsd  xmm1, [rsp + nb332nf_rsqOO] ;# xmm1=r 
	mulsd  xmm1, [rsp + nb332nf_tsc]

	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb332nf_VFtab]
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1  

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1  
	unpckhpd xmm7, xmm3	;# H1  
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb332nf_qqOO]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    ;# increment vcoul - then we can get rid of mm5 
    ;# update vctot 
    addsd  xmm5, [rsp + nb332nf_vctot]
    movlpd [rsp + nb332nf_vctot], xmm5

	;# Dispersion 
	movapd xmm4, [rsi + rax*8 + 32]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 48]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	mulsd  xmm5, [rsp + nb332nf_c6]	 ;# Vvdw6 

	addsd  xmm5, [rsp + nb332nf_Vvdwtot]
	movlpd [rsp + nb332nf_Vvdwtot], xmm5

	;# Repulsion 
	movapd xmm4, [rsi + rax*8 + 64]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 80]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# Repulsion table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	mulsd  xmm5, [rsp + nb332nf_c12] ;# Vvdw12 

	addsd  xmm5, [rsp + nb332nf_Vvdwtot]
	movlpd [rsp + nb332nf_Vvdwtot], xmm5

	;# O-H1 interaction 
	movapd xmm0, [rsp + nb332nf_rinvOH1]
	movapd xmm1, xmm0
	mulsd  xmm1, [rsp + nb332nf_rsqOH1] ;# xmm1=r 
	mulsd  xmm1, [rsp + nb332nf_tsc]

	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb332nf_VFtab]
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb332nf_qqOH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [rsp + nb332nf_vctot]
    movlpd [rsp + nb332nf_vctot], xmm5

	;# O-H2 interaction  
	movapd xmm0, [rsp + nb332nf_rinvOH2]
	movapd xmm1, xmm0
	mulsd  xmm1, [rsp + nb332nf_rsqOH2] ;# xmm1=r 
	mulsd  xmm1, [rsp + nb332nf_tsc]
	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb332nf_VFtab]
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb332nf_qqOH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [rsp + nb332nf_vctot]
    movlpd [rsp + nb332nf_vctot], xmm5

	;# H1-O interaction 
	movapd xmm0, [rsp + nb332nf_rinvH1O]
	movapd xmm1, xmm0
	mulsd  xmm1, [rsp + nb332nf_rsqH1O] ;# xmm1=r 
	mulsd  xmm1, [rsp + nb332nf_tsc]
	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb332nf_VFtab]
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb332nf_qqOH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [rsp + nb332nf_vctot]
    movlpd [rsp + nb332nf_vctot], xmm5

	;# H1-H1 interaction 
	movapd xmm0, [rsp + nb332nf_rinvH1H1]
	movapd xmm1, xmm0
	mulsd  xmm1, [rsp + nb332nf_rsqH1H1] ;# xmm1=r 
	mulsd  xmm1, [rsp + nb332nf_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb332nf_VFtab]
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb332nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV 
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [rsp + nb332nf_vctot]
    movlpd [rsp + nb332nf_vctot], xmm5

	;# H1-H2 interaction 
	movapd xmm0, [rsp + nb332nf_rinvH1H2]
	movapd xmm1, xmm0
	mulsd  xmm1, [rsp + nb332nf_rsqH1H2] ;# xmm1=r 
	mulsd  xmm1, [rsp + nb332nf_tsc]
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb332nf_VFtab]
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb332nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [rsp + nb332nf_vctot]
    movlpd [rsp + nb332nf_vctot], xmm5
	
	;# H2-O interaction 
	movapd xmm0, [rsp + nb332nf_rinvH2O]
	movapd xmm1, xmm0
	mulsd  xmm1, [rsp + nb332nf_rsqH2O] ;# xmm1=r 
	mulsd  xmm1, [rsp + nb332nf_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb332nf_VFtab]
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb332nf_qqOH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [rsp + nb332nf_vctot]
    movlpd [rsp + nb332nf_vctot], xmm5

	;# H2-H1 interaction 
	movapd xmm0, [rsp + nb332nf_rinvH2H1]
	movapd xmm1, xmm0
	mulsd  xmm1, [rsp + nb332nf_rsqH2H1] ;# xmm1=r 
	mulsd  xmm1, [rsp + nb332nf_tsc]
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subpd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb332nf_VFtab]
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb332nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [rsp + nb332nf_vctot]
    movlpd [rsp + nb332nf_vctot], xmm5

	;# H2-H2 interaction 
	movapd xmm0, [rsp + nb332nf_rinvH2H2]
	movapd xmm1, xmm0
	mulsd  xmm1, [rsp + nb332nf_rsqH2H2] ;# xmm1=r 
	mulsd  xmm1, [rsp + nb332nf_tsc]	
	cvttsd2si eax, xmm1	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm1, xmm6	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb332nf_VFtab]
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 

	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb332nf_qqHH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addsd  xmm5, [rsp + nb332nf_vctot]
    movlpd [rsp + nb332nf_vctot], xmm5
	
.nb332nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb332nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb332nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb332nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb332nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb332nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb332nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [rsp + nb332nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb332nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb332nf_n], esi
        jmp .nb332nf_outer
.nb332nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb332nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb332nf_end
        ;# non-zero, do one more workunit
        jmp   .nb332nf_threadloop
.nb332nf_end:
	mov eax, [rsp + nb332nf_nouter]
	mov ebx, [rsp + nb332nf_ninner]
	mov rcx, [rbp + nb332nf_outeriter]
	mov rdx, [rbp + nb332nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 832
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
