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



	
.globl nb_kernel112_x86_64_sse2
.globl _nb_kernel112_x86_64_sse2
nb_kernel112_x86_64_sse2:	
_nb_kernel112_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb112_fshift,           16
.equiv          nb112_gid,              24
.equiv          nb112_pos,              32
.equiv          nb112_faction,          40
.equiv          nb112_charge,           48
.equiv          nb112_p_facel,          56
.equiv          nb112_argkrf,           64
.equiv          nb112_argcrf,           72
.equiv          nb112_Vc,               80
.equiv          nb112_type,             88
.equiv          nb112_p_ntype,          96
.equiv          nb112_vdwparam,         104
.equiv          nb112_Vvdw,             112
.equiv          nb112_p_tabscale,       120
.equiv          nb112_VFtab,            128
.equiv          nb112_invsqrta,         136
.equiv          nb112_dvda,             144
.equiv          nb112_p_gbtabscale,     152
.equiv          nb112_GBtab,            160
.equiv          nb112_p_nthreads,       168
.equiv          nb112_count,            176
.equiv          nb112_mtx,              184
.equiv          nb112_outeriter,        192
.equiv          nb112_inneriter,        200
.equiv          nb112_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb112_ixO,              0
.equiv          nb112_iyO,              16
.equiv          nb112_izO,              32
.equiv          nb112_ixH1,             48
.equiv          nb112_iyH1,             64
.equiv          nb112_izH1,             80
.equiv          nb112_ixH2,             96
.equiv          nb112_iyH2,             112
.equiv          nb112_izH2,             128
.equiv          nb112_jxO,              144
.equiv          nb112_jyO,              160
.equiv          nb112_jzO,              176
.equiv          nb112_jxH1,             192
.equiv          nb112_jyH1,             208
.equiv          nb112_jzH1,             224
.equiv          nb112_jxH2,             240
.equiv          nb112_jyH2,             256
.equiv          nb112_jzH2,             272
.equiv          nb112_dxOO,             288
.equiv          nb112_dyOO,             304
.equiv          nb112_dzOO,             320
.equiv          nb112_dxOH1,            336
.equiv          nb112_dyOH1,            352
.equiv          nb112_dzOH1,            368
.equiv          nb112_dxOH2,            384
.equiv          nb112_dyOH2,            400
.equiv          nb112_dzOH2,            416
.equiv          nb112_dxH1O,            432
.equiv          nb112_dyH1O,            448
.equiv          nb112_dzH1O,            464
.equiv          nb112_dxH1H1,           480
.equiv          nb112_dyH1H1,           496
.equiv          nb112_dzH1H1,           512
.equiv          nb112_dxH1H2,           528
.equiv          nb112_dyH1H2,           544
.equiv          nb112_dzH1H2,           560
.equiv          nb112_dxH2O,            576
.equiv          nb112_dyH2O,            592
.equiv          nb112_dzH2O,            608
.equiv          nb112_dxH2H1,           624
.equiv          nb112_dyH2H1,           640
.equiv          nb112_dzH2H1,           656
.equiv          nb112_dxH2H2,           672
.equiv          nb112_dyH2H2,           688
.equiv          nb112_dzH2H2,           704
.equiv          nb112_qqOO,             720
.equiv          nb112_qqOH,             736
.equiv          nb112_qqHH,             752
.equiv          nb112_c6,               768
.equiv          nb112_c12,              784
.equiv          nb112_six,              800
.equiv          nb112_twelve,           816
.equiv          nb112_vctot,            832
.equiv          nb112_Vvdwtot,          848
.equiv          nb112_fixO,             864
.equiv          nb112_fiyO,             880
.equiv          nb112_fizO,             896
.equiv          nb112_fixH1,            912
.equiv          nb112_fiyH1,            928
.equiv          nb112_fizH1,            944
.equiv          nb112_fixH2,            960
.equiv          nb112_fiyH2,            976
.equiv          nb112_fizH2,            992
.equiv          nb112_fjxO,             1008
.equiv          nb112_fjyO,             1024
.equiv          nb112_fjzO,             1040
.equiv          nb112_fjxH1,            1056
.equiv          nb112_fjyH1,            1072
.equiv          nb112_fjzH1,            1088
.equiv          nb112_fjxH2,            1104
.equiv          nb112_fjyH2,            1120
.equiv          nb112_fjzH2,            1136
.equiv          nb112_half,             1152
.equiv          nb112_three,            1168
.equiv          nb112_rsqOO,            1184
.equiv          nb112_rsqOH1,           1200
.equiv          nb112_rsqOH2,           1216
.equiv          nb112_rsqH1O,           1232
.equiv          nb112_rsqH1H1,          1248
.equiv          nb112_rsqH1H2,          1264
.equiv          nb112_rsqH2O,           1280
.equiv          nb112_rsqH2H1,          1296
.equiv          nb112_rsqH2H2,          1312
.equiv          nb112_rinvOO,           1328
.equiv          nb112_rinvOH1,          1344
.equiv          nb112_rinvOH2,          1360
.equiv          nb112_rinvH1O,          1376
.equiv          nb112_rinvH1H1,         1392
.equiv          nb112_rinvH1H2,         1408
.equiv          nb112_rinvH2O,          1424
.equiv          nb112_rinvH2H1,         1440
.equiv          nb112_rinvH2H2,         1456
.equiv          nb112_is3,              1472
.equiv          nb112_ii3,              1476
.equiv          nb112_nri,              1492
.equiv          nb112_iinr,             1500
.equiv          nb112_jindex,           1508
.equiv          nb112_jjnr,             1516
.equiv          nb112_shift,            1524
.equiv          nb112_shiftvec,         1532
.equiv          nb112_facel,            1540
.equiv          nb112_innerjjnr,        1548
.equiv          nb112_innerk,           1556
.equiv          nb112_n,                1560
.equiv          nb112_nn1,              1564
.equiv          nb112_nouter,           1568
.equiv          nb112_ninner,           1572

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
	mov [rsp + nb112_nouter], eax
	mov [rsp + nb112_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb112_nri], edi
	mov [rsp + nb112_iinr], rsi
	mov [rsp + nb112_jindex], rdx
	mov [rsp + nb112_jjnr], rcx
	mov [rsp + nb112_shift], r8
	mov [rsp + nb112_shiftvec], r9
	mov rsi, [rbp + nb112_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb112_facel], xmm0


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb112_half], eax
	mov [rsp + nb112_half+4], ebx
	movsd xmm1, [rsp + nb112_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd xmm4, xmm3
	addpd  xmm4, xmm4       ;# six
	movapd xmm5, xmm4
	addpd  xmm5, xmm5       ;# twelve
	movapd [rsp + nb112_half], xmm1
	movapd [rsp + nb112_three], xmm3
	movapd [rsp + nb112_six], xmm4
	movapd [rsp + nb112_twelve], xmm5

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb112_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb112_charge]
	movsd xmm3, [rdx + rbx*8]	
	movsd xmm4, xmm3	
	movsd xmm5, [rdx + rbx*8 + 8]	
	mov rsi, [rbp + nb112_p_facel]
	movsd xmm0, [rsi]
	movsd xmm6, [rsp + nb112_facel]
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm5
	mulsd  xmm5, xmm5
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb112_qqOO], xmm3
	movapd [rsp + nb112_qqOH], xmm4
	movapd [rsp + nb112_qqHH], xmm5
	
	xorpd xmm0, xmm0
	mov   rdx, [rbp + nb112_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb112_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   rax, [rbp + nb112_vdwparam]
	movlpd xmm0, [rax + rdx*8]
	movlpd xmm1, [rax + rdx*8 + 8]
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0	
	movapd [rsp + nb112_c6], xmm0
	movapd [rsp + nb112_c12], xmm1

.nb112_threadloop:
        mov   rsi, [rbp + nb112_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb112_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb112_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb112_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb112_n], eax
        mov [rsp + nb112_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb112_outerstart
        jmp .nb112_end

.nb112_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb112_nouter]
	mov [rsp + nb112_nouter], ebx

.nb112_outer:
	mov   rax, [rsp + nb112_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb112_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb112_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb112_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb112_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb112_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [rax + rbx*8]
	addsd xmm4, [rax + rbx*8 + 8]
	addsd xmm5, [rax + rbx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb112_ixO], xmm3
	movapd [rsp + nb112_iyO], xmm4
	movapd [rsp + nb112_izO], xmm5

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
	movapd [rsp + nb112_ixH1], xmm0
	movapd [rsp + nb112_iyH1], xmm1
	movapd [rsp + nb112_izH1], xmm2
	movapd [rsp + nb112_ixH2], xmm3
	movapd [rsp + nb112_iyH2], xmm4
	movapd [rsp + nb112_izH2], xmm5

	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb112_vctot], xmm4
	movapd [rsp + nb112_Vvdwtot], xmm4
	movapd [rsp + nb112_fixO], xmm4
	movapd [rsp + nb112_fiyO], xmm4
	movapd [rsp + nb112_fizO], xmm4
	movapd [rsp + nb112_fixH1], xmm4
	movapd [rsp + nb112_fiyH1], xmm4
	movapd [rsp + nb112_fizH1], xmm4
	movapd [rsp + nb112_fixH2], xmm4
	movapd [rsp + nb112_fiyH2], xmm4
	movapd [rsp + nb112_fizH2], xmm4
	
	mov   rax, [rsp + nb112_jindex]
	mov   ecx, [rax+rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb112_pos]
	mov   rdi, [rbp + nb112_faction]	
	mov   rax, [rsp + nb112_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb112_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb112_ninner]
	mov   [rsp + nb112_ninner], ecx
	add   edx, 0
	mov   [rsp + nb112_innerk], edx    ;# number of innerloop atoms 
	jge  .nb112_unroll_loop
	jmp  .nb112_checksingle
.nb112_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb112_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	
	add qword ptr [rsp + nb112_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb112_pos]       ;# base of pos[] 

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
    
    subpd xmm0, [rsp + nb112_ixO]
    subpd xmm1, [rsp + nb112_iyO]
    subpd xmm2, [rsp + nb112_izO]
    subpd xmm3, [rsp + nb112_ixH1]
    subpd xmm4, [rsp + nb112_iyH1]
    subpd xmm5, [rsp + nb112_izH1]
    subpd xmm6, [rsp + nb112_ixH2]
    subpd xmm7, [rsp + nb112_iyH2]
    subpd xmm8, [rsp + nb112_izH2]
    
	movapd [rsp + nb112_dxOO], xmm0
	movapd [rsp + nb112_dyOO], xmm1
	movapd [rsp + nb112_dzOO], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb112_dxH1O], xmm3
	movapd [rsp + nb112_dyH1O], xmm4
	movapd [rsp + nb112_dzH1O], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb112_dxH2O], xmm6
	movapd [rsp + nb112_dyH2O], xmm7
	movapd [rsp + nb112_dzH2O], xmm8
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
		
	movapd  xmm9, [rsp + nb112_three]
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

	movapd  xmm15, [rsp + nb112_half]
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
		
	movapd  xmm1, [rsp + nb112_three]
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

	movapd  xmm15, [rsp + nb112_half]
	mulpd   xmm9, xmm15  ;#  rinvOO 
	mulpd   xmm10, xmm15 ;#   rinvH1O
    mulpd   xmm11, xmm15 ;#   rinvH2O
	
	
	;# O interactions 
    movapd xmm0, xmm9
    movapd xmm1, xmm10
    movapd xmm2, xmm11
    mulpd  xmm9, xmm9    ;# rinvsq
    mulpd  xmm10, xmm10
    mulpd  xmm11, xmm11
    movapd xmm12, xmm9
    mulpd  xmm12, xmm12 ;# rinv4
    mulpd  xmm12, xmm9  ;# rinv6
    mulpd  xmm0, [rsp + nb112_qqOO] 
    mulpd  xmm1, [rsp + nb112_qqOH] 
    mulpd  xmm2, [rsp + nb112_qqOH] 
    movapd xmm13, xmm12 ;# rinv6
    mulpd  xmm12, xmm12 ;# rinv12
	mulpd  xmm13, [rsp + nb112_c6]
	mulpd  xmm12, [rsp + nb112_c12]
    movapd xmm14, xmm12
    subpd  xmm14, xmm13
    
	addpd  xmm14, [rsp + nb112_Vvdwtot]
	mulpd  xmm13, [rsp + nb112_six]
	mulpd  xmm12, [rsp + nb112_twelve]
	movapd [rsp + nb112_Vvdwtot], xmm14
    subpd  xmm12, xmm13 ;# LJ fscal        

    addpd  xmm12, xmm0
    
    mulpd  xmm9, xmm12
    mulpd  xmm10, xmm1
    mulpd  xmm11, xmm2
    
    addpd xmm0, [rsp + nb112_vctot] 
    addpd xmm1, xmm2
    addpd xmm0, xmm1
    movapd [rsp + nb112_vctot], xmm0
    
    ;# move j O forces to xmm0-xmm2
    mov rdi, [rbp + nb112_faction]
	movlpd xmm0, [rdi + rax*8]
	movlpd xmm1, [rdi + rax*8 + 8]
	movlpd xmm2, [rdi + rax*8 + 16]

    movapd xmm7, xmm9
    movapd xmm8, xmm9
    movapd xmm13, xmm11
    movapd xmm14, xmm11
    movapd xmm15, xmm11
    movapd xmm11, xmm10
    movapd xmm12, xmm10

	movhpd xmm0, [rdi + rbx*8]
	movhpd xmm1, [rdi + rbx*8 + 8]
	movhpd xmm2, [rdi + rbx*8 + 16]

	mulpd xmm7, [rsp + nb112_dxOO]
	mulpd xmm8, [rsp + nb112_dyOO]
	mulpd xmm9, [rsp + nb112_dzOO]
	mulpd xmm10, [rsp + nb112_dxH1O]
	mulpd xmm11, [rsp + nb112_dyH1O]
	mulpd xmm12, [rsp + nb112_dzH1O]
	mulpd xmm13, [rsp + nb112_dxH2O]
	mulpd xmm14, [rsp + nb112_dyH2O]
	mulpd xmm15, [rsp + nb112_dzH2O]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb112_fixO]
    addpd xmm8, [rsp + nb112_fiyO]
    addpd xmm9, [rsp + nb112_fizO]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb112_fixH1]
    addpd xmm11, [rsp + nb112_fiyH1]
    addpd xmm12, [rsp + nb112_fizH1]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb112_fixH2]
    addpd xmm14, [rsp + nb112_fiyH2]
    addpd xmm15, [rsp + nb112_fizH2]

    movapd [rsp + nb112_fixO], xmm7
    movapd [rsp + nb112_fiyO], xmm8
    movapd [rsp + nb112_fizO], xmm9
    movapd [rsp + nb112_fixH1], xmm10
    movapd [rsp + nb112_fiyH1], xmm11
    movapd [rsp + nb112_fizH1], xmm12
    movapd [rsp + nb112_fixH2], xmm13
    movapd [rsp + nb112_fiyH2], xmm14
    movapd [rsp + nb112_fizH2], xmm15
   
    ;# store back j O forces from xmm0-xmm2
	movlpd [rdi + rax*8],      xmm0
	movlpd [rdi + rax*8 + 8],  xmm1
	movlpd [rdi + rax*8 + 16], xmm2
	movhpd [rdi + rbx*8],      xmm0
	movhpd [rdi + rbx*8 + 8],  xmm1
	movhpd [rdi + rbx*8 + 16], xmm2

	;# move j H1 coordinates to local temp variables 
    mov rsi, [rbp + nb112_pos]
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
    
    subpd xmm0, [rsp + nb112_ixO]
    subpd xmm1, [rsp + nb112_iyO]
    subpd xmm2, [rsp + nb112_izO]
    subpd xmm3, [rsp + nb112_ixH1]
    subpd xmm4, [rsp + nb112_iyH1]
    subpd xmm5, [rsp + nb112_izH1]
    subpd xmm6, [rsp + nb112_ixH2]
    subpd xmm7, [rsp + nb112_iyH2]
    subpd xmm8, [rsp + nb112_izH2]
    
	movapd [rsp + nb112_dxOH1], xmm0
	movapd [rsp + nb112_dyOH1], xmm1
	movapd [rsp + nb112_dzOH1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb112_dxH1H1], xmm3
	movapd [rsp + nb112_dyH1H1], xmm4
	movapd [rsp + nb112_dzH1H1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb112_dxH2H1], xmm6
	movapd [rsp + nb112_dyH2H1], xmm7
	movapd [rsp + nb112_dzH2H1], xmm8
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
		
	movapd  xmm9, [rsp + nb112_three]
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

	movapd  xmm15, [rsp + nb112_half]
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
		
	movapd  xmm1, [rsp + nb112_three]
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

	movapd  xmm15, [rsp + nb112_half]
	mulpd   xmm9, xmm15  ;#  rinvOH1
	mulpd   xmm10, xmm15 ;#   rinvH1H1
    mulpd   xmm11, xmm15 ;#   rinvH2H1
	
	;# H1 interactions 
    movapd xmm0, xmm9
    movapd xmm1, xmm10
    movapd xmm2, xmm11
    mulpd  xmm9, xmm9
    mulpd  xmm10, xmm10
    mulpd  xmm11, xmm11
    mulpd  xmm0, [rsp + nb112_qqOH] 
    mulpd  xmm1, [rsp + nb112_qqHH] 
    mulpd  xmm2, [rsp + nb112_qqHH] 
    mulpd  xmm9, xmm0
    mulpd  xmm10, xmm1
    mulpd  xmm11, xmm2
    
    addpd xmm0, [rsp + nb112_vctot] 
    addpd xmm1, xmm2
    addpd xmm0, xmm1
    movapd [rsp + nb112_vctot], xmm0
    
    ;# move j H1 forces to xmm0-xmm2
    mov rdi, [rbp + nb112_faction]
	movlpd xmm0, [rdi + rax*8 + 24]
	movlpd xmm1, [rdi + rax*8 + 32]
	movlpd xmm2, [rdi + rax*8 + 40]

    movapd xmm7, xmm9
    movapd xmm8, xmm9
    movapd xmm13, xmm11
    movapd xmm14, xmm11
    movapd xmm15, xmm11
    movapd xmm11, xmm10
    movapd xmm12, xmm10

	movhpd xmm0, [rdi + rbx*8 + 24]
	movhpd xmm1, [rdi + rbx*8 + 32]
	movhpd xmm2, [rdi + rbx*8 + 40]

	mulpd xmm7, [rsp + nb112_dxOH1]
	mulpd xmm8, [rsp + nb112_dyOH1]
	mulpd xmm9, [rsp + nb112_dzOH1]
	mulpd xmm10, [rsp + nb112_dxH1H1]
	mulpd xmm11, [rsp + nb112_dyH1H1]
	mulpd xmm12, [rsp + nb112_dzH1H1]
	mulpd xmm13, [rsp + nb112_dxH2H1]
	mulpd xmm14, [rsp + nb112_dyH2H1]
	mulpd xmm15, [rsp + nb112_dzH2H1]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb112_fixO]
    addpd xmm8, [rsp + nb112_fiyO]
    addpd xmm9, [rsp + nb112_fizO]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb112_fixH1]
    addpd xmm11, [rsp + nb112_fiyH1]
    addpd xmm12, [rsp + nb112_fizH1]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb112_fixH2]
    addpd xmm14, [rsp + nb112_fiyH2]
    addpd xmm15, [rsp + nb112_fizH2]

    movapd [rsp + nb112_fixO], xmm7
    movapd [rsp + nb112_fiyO], xmm8
    movapd [rsp + nb112_fizO], xmm9
    movapd [rsp + nb112_fixH1], xmm10
    movapd [rsp + nb112_fiyH1], xmm11
    movapd [rsp + nb112_fizH1], xmm12
    movapd [rsp + nb112_fixH2], xmm13
    movapd [rsp + nb112_fiyH2], xmm14
    movapd [rsp + nb112_fizH2], xmm15
   
    ;# store back j H1 forces from xmm0-xmm2
	movlpd [rdi + rax*8 + 24], xmm0
	movlpd [rdi + rax*8 + 32], xmm1
	movlpd [rdi + rax*8 + 40], xmm2
	movhpd [rdi + rbx*8 + 24], xmm0
	movhpd [rdi + rbx*8 + 32], xmm1
	movhpd [rdi + rbx*8 + 40], xmm2
       
	;# move j H2 coordinates to local temp variables 
    mov rsi, [rbp + nb112_pos]
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
    
    subpd xmm0, [rsp + nb112_ixO]
    subpd xmm1, [rsp + nb112_iyO]
    subpd xmm2, [rsp + nb112_izO]
    subpd xmm3, [rsp + nb112_ixH1]
    subpd xmm4, [rsp + nb112_iyH1]
    subpd xmm5, [rsp + nb112_izH1]
    subpd xmm6, [rsp + nb112_ixH2]
    subpd xmm7, [rsp + nb112_iyH2]
    subpd xmm8, [rsp + nb112_izH2]
    
	movapd [rsp + nb112_dxOH2], xmm0
	movapd [rsp + nb112_dyOH2], xmm1
	movapd [rsp + nb112_dzOH2], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb112_dxH1H2], xmm3
	movapd [rsp + nb112_dyH1H2], xmm4
	movapd [rsp + nb112_dzH1H2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb112_dxH2H2], xmm6
	movapd [rsp + nb112_dyH2H2], xmm7
	movapd [rsp + nb112_dzH2H2], xmm8
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
		
	movapd  xmm9, [rsp + nb112_three]
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

	movapd  xmm15, [rsp + nb112_half]
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
		
	movapd  xmm1, [rsp + nb112_three]
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

	movapd  xmm15, [rsp + nb112_half]
	mulpd   xmm9, xmm15  ;#  rinvOH2
	mulpd   xmm10, xmm15 ;#   rinvH1H2
    mulpd   xmm11, xmm15 ;#   rinvH2H2
	
	;# H2 interactions 
    movapd xmm0, xmm9
    movapd xmm1, xmm10
    movapd xmm2, xmm11
    mulpd  xmm9, xmm9
    mulpd  xmm10, xmm10
    mulpd  xmm11, xmm11
    mulpd  xmm0, [rsp + nb112_qqOH] 
    mulpd  xmm1, [rsp + nb112_qqHH] 
    mulpd  xmm2, [rsp + nb112_qqHH] 
    mulpd  xmm9, xmm0
    mulpd  xmm10, xmm1
    mulpd  xmm11, xmm2
    
    addpd xmm0, [rsp + nb112_vctot] 
    addpd xmm1, xmm2
    addpd xmm0, xmm1
    movapd [rsp + nb112_vctot], xmm0
    
    ;# move j H2 forces to xmm0-xmm2
    mov rdi, [rbp + nb112_faction]
	movlpd xmm0, [rdi + rax*8 + 48]
	movlpd xmm1, [rdi + rax*8 + 56]
	movlpd xmm2, [rdi + rax*8 + 64]

    movapd xmm7, xmm9
    movapd xmm8, xmm9
    movapd xmm13, xmm11
    movapd xmm14, xmm11
    movapd xmm15, xmm11
    movapd xmm11, xmm10
    movapd xmm12, xmm10

	movhpd xmm0, [rdi + rbx*8 + 48]
	movhpd xmm1, [rdi + rbx*8 + 56]
	movhpd xmm2, [rdi + rbx*8 + 64]

	mulpd xmm7, [rsp + nb112_dxOH2]
	mulpd xmm8, [rsp + nb112_dyOH2]
	mulpd xmm9, [rsp + nb112_dzOH2]
	mulpd xmm10, [rsp + nb112_dxH1H2]
	mulpd xmm11, [rsp + nb112_dyH1H2]
	mulpd xmm12, [rsp + nb112_dzH1H2]
	mulpd xmm13, [rsp + nb112_dxH2H2]
	mulpd xmm14, [rsp + nb112_dyH2H2]
	mulpd xmm15, [rsp + nb112_dzH2H2]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb112_fixO]
    addpd xmm8, [rsp + nb112_fiyO]
    addpd xmm9, [rsp + nb112_fizO]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb112_fixH1]
    addpd xmm11, [rsp + nb112_fiyH1]
    addpd xmm12, [rsp + nb112_fizH1]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb112_fixH2]
    addpd xmm14, [rsp + nb112_fiyH2]
    addpd xmm15, [rsp + nb112_fizH2]

    movapd [rsp + nb112_fixO], xmm7
    movapd [rsp + nb112_fiyO], xmm8
    movapd [rsp + nb112_fizO], xmm9
    movapd [rsp + nb112_fixH1], xmm10
    movapd [rsp + nb112_fiyH1], xmm11
    movapd [rsp + nb112_fizH1], xmm12
    movapd [rsp + nb112_fixH2], xmm13
    movapd [rsp + nb112_fiyH2], xmm14
    movapd [rsp + nb112_fizH2], xmm15
   
    ;# store back j H2 forces from xmm0-xmm2
	movlpd [rdi + rax*8 + 48], xmm0
	movlpd [rdi + rax*8 + 56], xmm1
	movlpd [rdi + rax*8 + 64], xmm2
	movhpd [rdi + rbx*8 + 48], xmm0
	movhpd [rdi + rbx*8 + 56], xmm1
	movhpd [rdi + rbx*8 + 64], xmm2
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb112_innerk],  2
	jl    .nb112_checksingle
	jmp   .nb112_unroll_loop
.nb112_checksingle:
	mov   edx, [rsp + nb112_innerk]
	and   edx, 1
	jnz   .nb112_dosingle
	jmp   .nb112_updateouterdata
.nb112_dosingle:
	mov   rdx, [rsp + nb112_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	

	mov rsi, [rbp + nb112_pos]
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
    
    subsd xmm0, [rsp + nb112_ixO]
    subsd xmm1, [rsp + nb112_iyO]
    subsd xmm2, [rsp + nb112_izO]
    subsd xmm3, [rsp + nb112_ixH1]
    subsd xmm4, [rsp + nb112_iyH1]
    subsd xmm5, [rsp + nb112_izH1]
    subsd xmm6, [rsp + nb112_ixH2]
    subsd xmm7, [rsp + nb112_iyH2]
    subsd xmm8, [rsp + nb112_izH2]
    
	movsd [rsp + nb112_dxOO], xmm0
	movsd [rsp + nb112_dyOO], xmm1
	movsd [rsp + nb112_dzOO], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movsd [rsp + nb112_dxH1O], xmm3
	movsd [rsp + nb112_dyH1O], xmm4
	movsd [rsp + nb112_dzH1O], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movsd [rsp + nb112_dxH2O], xmm6
	movsd [rsp + nb112_dyH2O], xmm7
	movsd [rsp + nb112_dzH2O], xmm8
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
		
	movsd  xmm9, [rsp + nb112_three]
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

	movsd  xmm15, [rsp + nb112_half]
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
		
	movsd  xmm1, [rsp + nb112_three]
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

	movsd  xmm15, [rsp + nb112_half]
	mulsd   xmm9, xmm15  ;#  rinvOO 
	mulsd   xmm10, xmm15 ;#   rinvH1O
    mulsd   xmm11, xmm15 ;#   rinvH2O
	
	;# O interactions 
    movsd xmm0, xmm9
    movsd xmm1, xmm10
    movsd xmm2, xmm11
    mulsd  xmm9, xmm9    ;# rinvsq
    mulsd  xmm10, xmm10
    mulsd  xmm11, xmm11
    movsd xmm12, xmm9
    mulsd  xmm12, xmm12 ;# rinv4
    mulsd  xmm12, xmm9  ;# rinv6
    mulsd  xmm0, [rsp + nb112_qqOO] 
    mulsd  xmm1, [rsp + nb112_qqOH] 
    mulsd  xmm2, [rsp + nb112_qqOH] 
    movsd xmm13, xmm12 ;# rinv6
    mulsd  xmm12, xmm12 ;# rinv12
	mulsd  xmm13, [rsp + nb112_c6]
	mulsd  xmm12, [rsp + nb112_c12]
    movsd xmm14, xmm12
    subsd  xmm14, xmm13
    
	addsd  xmm14, [rsp + nb112_Vvdwtot]
	mulsd  xmm13, [rsp + nb112_six]
	mulsd  xmm12, [rsp + nb112_twelve]
	movsd [rsp + nb112_Vvdwtot], xmm14
    subsd  xmm12, xmm13 ;# LJ fscal        

    addsd  xmm12, xmm0
    
    mulsd  xmm9, xmm12
    mulsd  xmm10, xmm1
    mulsd  xmm11, xmm2
    
    addsd xmm0, [rsp + nb112_vctot] 
    addsd xmm1, xmm2
    addsd xmm0, xmm1
    movsd [rsp + nb112_vctot], xmm0
    
    ;# move j O forces to xmm0-xmm2
	movsd xmm0, [rdi + rax*8]
	movsd xmm1, [rdi + rax*8 + 8]
	movsd xmm2, [rdi + rax*8 + 16]

    movsd xmm7, xmm9
    movsd xmm8, xmm9
    movsd xmm13, xmm11
    movsd xmm14, xmm11
    movsd xmm15, xmm11
    movsd xmm11, xmm10
    movsd xmm12, xmm10

	mulsd xmm7, [rsp + nb112_dxOO]
	mulsd xmm8, [rsp + nb112_dyOO]
	mulsd xmm9, [rsp + nb112_dzOO]
	mulsd xmm10, [rsp + nb112_dxH1O]
	mulsd xmm11, [rsp + nb112_dyH1O]
	mulsd xmm12, [rsp + nb112_dzH1O]
	mulsd xmm13, [rsp + nb112_dxH2O]
	mulsd xmm14, [rsp + nb112_dyH2O]
	mulsd xmm15, [rsp + nb112_dzH2O]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb112_fixO]
    addsd xmm8, [rsp + nb112_fiyO]
    addsd xmm9, [rsp + nb112_fizO]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb112_fixH1]
    addsd xmm11, [rsp + nb112_fiyH1]
    addsd xmm12, [rsp + nb112_fizH1]

    addsd xmm0, xmm13
    addsd xmm1, xmm14
    addsd xmm2, xmm15
    addsd xmm13, [rsp + nb112_fixH2]
    addsd xmm14, [rsp + nb112_fiyH2]
    addsd xmm15, [rsp + nb112_fizH2]

    movsd [rsp + nb112_fixO], xmm7
    movsd [rsp + nb112_fiyO], xmm8
    movsd [rsp + nb112_fizO], xmm9
    movsd [rsp + nb112_fixH1], xmm10
    movsd [rsp + nb112_fiyH1], xmm11
    movsd [rsp + nb112_fizH1], xmm12
    movsd [rsp + nb112_fixH2], xmm13
    movsd [rsp + nb112_fiyH2], xmm14
    movsd [rsp + nb112_fizH2], xmm15
   
    ;# store back j O forces from xmm0-xmm2
	movsd [rdi + rax*8],      xmm0
	movsd [rdi + rax*8 + 8],  xmm1
	movsd [rdi + rax*8 + 16], xmm2

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
    
    subsd xmm0, [rsp + nb112_ixO]
    subsd xmm1, [rsp + nb112_iyO]
    subsd xmm2, [rsp + nb112_izO]
    subsd xmm3, [rsp + nb112_ixH1]
    subsd xmm4, [rsp + nb112_iyH1]
    subsd xmm5, [rsp + nb112_izH1]
    subsd xmm6, [rsp + nb112_ixH2]
    subsd xmm7, [rsp + nb112_iyH2]
    subsd xmm8, [rsp + nb112_izH2]
    
	movsd [rsp + nb112_dxOH1], xmm0
	movsd [rsp + nb112_dyOH1], xmm1
	movsd [rsp + nb112_dzOH1], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movsd [rsp + nb112_dxH1H1], xmm3
	movsd [rsp + nb112_dyH1H1], xmm4
	movsd [rsp + nb112_dzH1H1], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movsd [rsp + nb112_dxH2H1], xmm6
	movsd [rsp + nb112_dyH2H1], xmm7
	movsd [rsp + nb112_dzH2H1], xmm8
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
		
	movsd  xmm9, [rsp + nb112_three]
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

	movsd  xmm15, [rsp + nb112_half]
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
		
	movsd  xmm1, [rsp + nb112_three]
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

	movsd  xmm15, [rsp + nb112_half]
	mulsd   xmm9, xmm15  ;#  rinvOH1
	mulsd   xmm10, xmm15 ;#   rinvH1H1
    mulsd   xmm11, xmm15 ;#   rinvH2H1
	
	;# H1 interactions 
    movsd xmm0, xmm9
    movsd xmm1, xmm10
    movsd xmm2, xmm11
    mulsd  xmm9, xmm9
    mulsd  xmm10, xmm10
    mulsd  xmm11, xmm11
    mulsd  xmm0, [rsp + nb112_qqOH] 
    mulsd  xmm1, [rsp + nb112_qqHH] 
    mulsd  xmm2, [rsp + nb112_qqHH] 
    mulsd  xmm9, xmm0
    mulsd  xmm10, xmm1
    mulsd  xmm11, xmm2
    
    addsd xmm0, [rsp + nb112_vctot] 
    addsd xmm1, xmm2
    addsd xmm0, xmm1
    movsd [rsp + nb112_vctot], xmm0
    
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

	mulsd xmm7, [rsp + nb112_dxOH1]
	mulsd xmm8, [rsp + nb112_dyOH1]
	mulsd xmm9, [rsp + nb112_dzOH1]
	mulsd xmm10, [rsp + nb112_dxH1H1]
	mulsd xmm11, [rsp + nb112_dyH1H1]
	mulsd xmm12, [rsp + nb112_dzH1H1]
	mulsd xmm13, [rsp + nb112_dxH2H1]
	mulsd xmm14, [rsp + nb112_dyH2H1]
	mulsd xmm15, [rsp + nb112_dzH2H1]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb112_fixO]
    addsd xmm8, [rsp + nb112_fiyO]
    addsd xmm9, [rsp + nb112_fizO]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb112_fixH1]
    addsd xmm11, [rsp + nb112_fiyH1]
    addsd xmm12, [rsp + nb112_fizH1]

    addsd xmm0, xmm13
    addsd xmm1, xmm14
    addsd xmm2, xmm15
    addsd xmm13, [rsp + nb112_fixH2]
    addsd xmm14, [rsp + nb112_fiyH2]
    addsd xmm15, [rsp + nb112_fizH2]

    movsd [rsp + nb112_fixO], xmm7
    movsd [rsp + nb112_fiyO], xmm8
    movsd [rsp + nb112_fizO], xmm9
    movsd [rsp + nb112_fixH1], xmm10
    movsd [rsp + nb112_fiyH1], xmm11
    movsd [rsp + nb112_fizH1], xmm12
    movsd [rsp + nb112_fixH2], xmm13
    movsd [rsp + nb112_fiyH2], xmm14
    movsd [rsp + nb112_fizH2], xmm15
   
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
    
    subsd xmm0, [rsp + nb112_ixO]
    subsd xmm1, [rsp + nb112_iyO]
    subsd xmm2, [rsp + nb112_izO]
    subsd xmm3, [rsp + nb112_ixH1]
    subsd xmm4, [rsp + nb112_iyH1]
    subsd xmm5, [rsp + nb112_izH1]
    subsd xmm6, [rsp + nb112_ixH2]
    subsd xmm7, [rsp + nb112_iyH2]
    subsd xmm8, [rsp + nb112_izH2]
    
	movsd [rsp + nb112_dxOH2], xmm0
	movsd [rsp + nb112_dyOH2], xmm1
	movsd [rsp + nb112_dzOH2], xmm2
	mulsd  xmm0, xmm0
	mulsd  xmm1, xmm1
	mulsd  xmm2, xmm2
	movsd [rsp + nb112_dxH1H2], xmm3
	movsd [rsp + nb112_dyH1H2], xmm4
	movsd [rsp + nb112_dzH1H2], xmm5
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm4
	mulsd  xmm5, xmm5
	movsd [rsp + nb112_dxH2H2], xmm6
	movsd [rsp + nb112_dyH2H2], xmm7
	movsd [rsp + nb112_dzH2H2], xmm8
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
		
	movsd  xmm9, [rsp + nb112_three]
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

	movsd  xmm15, [rsp + nb112_half]
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
		
	movsd  xmm1, [rsp + nb112_three]
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

	movsd  xmm15, [rsp + nb112_half]
	mulsd   xmm9, xmm15  ;#  rinvOH2
	mulsd   xmm10, xmm15 ;#   rinvH1H2
    mulsd   xmm11, xmm15 ;#   rinvH2H2
	
	;# H2 interactions 
    movsd xmm0, xmm9
    movsd xmm1, xmm10
    movsd xmm2, xmm11
    mulsd  xmm9, xmm9
    mulsd  xmm10, xmm10
    mulsd  xmm11, xmm11
    mulsd  xmm0, [rsp + nb112_qqOH] 
    mulsd  xmm1, [rsp + nb112_qqHH] 
    mulsd  xmm2, [rsp + nb112_qqHH] 
    mulsd  xmm9, xmm0
    mulsd  xmm10, xmm1
    mulsd  xmm11, xmm2
    
    addsd xmm0, [rsp + nb112_vctot] 
    addsd xmm1, xmm2
    addsd xmm0, xmm1
    movsd [rsp + nb112_vctot], xmm0
    
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

	mulsd xmm7, [rsp + nb112_dxOH2]
	mulsd xmm8, [rsp + nb112_dyOH2]
	mulsd xmm9, [rsp + nb112_dzOH2]
	mulsd xmm10, [rsp + nb112_dxH1H2]
	mulsd xmm11, [rsp + nb112_dyH1H2]
	mulsd xmm12, [rsp + nb112_dzH1H2]
	mulsd xmm13, [rsp + nb112_dxH2H2]
	mulsd xmm14, [rsp + nb112_dyH2H2]
	mulsd xmm15, [rsp + nb112_dzH2H2]

    addsd xmm0, xmm7
    addsd xmm1, xmm8
    addsd xmm2, xmm9
    addsd xmm7, [rsp + nb112_fixO]
    addsd xmm8, [rsp + nb112_fiyO]
    addsd xmm9, [rsp + nb112_fizO]

    addsd xmm0, xmm10
    addsd xmm1, xmm11
    addsd xmm2, xmm12
    addsd xmm10, [rsp + nb112_fixH1]
    addsd xmm11, [rsp + nb112_fiyH1]
    addsd xmm12, [rsp + nb112_fizH1]

    addsd xmm0, xmm13
    addsd xmm1, xmm14
    addsd xmm2, xmm15
    addsd xmm13, [rsp + nb112_fixH2]
    addsd xmm14, [rsp + nb112_fiyH2]
    addsd xmm15, [rsp + nb112_fizH2]

    movsd [rsp + nb112_fixO], xmm7
    movsd [rsp + nb112_fiyO], xmm8
    movsd [rsp + nb112_fizO], xmm9
    movsd [rsp + nb112_fixH1], xmm10
    movsd [rsp + nb112_fiyH1], xmm11
    movsd [rsp + nb112_fizH1], xmm12
    movsd [rsp + nb112_fixH2], xmm13
    movsd [rsp + nb112_fiyH2], xmm14
    movsd [rsp + nb112_fizH2], xmm15
   
    ;# store back j H2 forces from xmm0-xmm2
	movsd [rdi + rax*8 + 48], xmm0
	movsd [rdi + rax*8 + 56], xmm1
	movsd [rdi + rax*8 + 64], xmm2
	
.nb112_updateouterdata:
	mov   ecx, [rsp + nb112_ii3]
	mov   rdi, [rbp + nb112_faction]
	mov   rsi, [rbp + nb112_fshift]
	mov   edx, [rsp + nb112_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb112_fixO]
	movapd xmm1, [rsp + nb112_fiyO]
	movapd xmm2, [rsp + nb112_fizO]

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
	movapd xmm0, [rsp + nb112_fixH1]
	movapd xmm1, [rsp + nb112_fiyH1]
	movapd xmm2, [rsp + nb112_fizH1]

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
	movapd xmm0, [rsp + nb112_fixH2]
	movapd xmm1, [rsp + nb112_fiyH2]
	movapd xmm2, [rsp + nb112_fizH2]

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
	mov esi, [rsp + nb112_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb112_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb112_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb112_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb112_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   rax, [rbp + nb112_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [rsp + nb112_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb112_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb112_n], esi
        jmp .nb112_outer
.nb112_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb112_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb112_end
        ;# non-zero, do one more workunit
        jmp   .nb112_threadloop
.nb112_end:
	mov eax, [rsp + nb112_nouter]
	mov ebx, [rsp + nb112_ninner]
	mov rcx, [rbp + nb112_outeriter]
	mov rdx, [rbp + nb112_inneriter]
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



	
.globl nb_kernel112nf_x86_64_sse2
.globl _nb_kernel112nf_x86_64_sse2
nb_kernel112nf_x86_64_sse2:	
_nb_kernel112nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb112nf_fshift,         16
.equiv          nb112nf_gid,            24
.equiv          nb112nf_pos,            32
.equiv          nb112nf_faction,        40
.equiv          nb112nf_charge,         48
.equiv          nb112nf_p_facel,        56
.equiv          nb112nf_argkrf,         64
.equiv          nb112nf_argcrf,         72
.equiv          nb112nf_Vc,             80
.equiv          nb112nf_type,           88
.equiv          nb112nf_p_ntype,        96
.equiv          nb112nf_vdwparam,       104
.equiv          nb112nf_Vvdw,           112
.equiv          nb112nf_p_tabscale,     120
.equiv          nb112nf_VFtab,          128
.equiv          nb112nf_invsqrta,       136
.equiv          nb112nf_dvda,           144
.equiv          nb112nf_p_gbtabscale,   152
.equiv          nb112nf_GBtab,          160
.equiv          nb112nf_p_nthreads,     168
.equiv          nb112nf_count,          176
.equiv          nb112nf_mtx,            184
.equiv          nb112nf_outeriter,      192
.equiv          nb112nf_inneriter,      200
.equiv          nb112nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb112nf_ixO,            0
.equiv          nb112nf_iyO,            16
.equiv          nb112nf_izO,            32
.equiv          nb112nf_ixH1,           48
.equiv          nb112nf_iyH1,           64
.equiv          nb112nf_izH1,           80
.equiv          nb112nf_ixH2,           96
.equiv          nb112nf_iyH2,           112
.equiv          nb112nf_izH2,           128
.equiv          nb112nf_jxO,            144
.equiv          nb112nf_jyO,            160
.equiv          nb112nf_jzO,            176
.equiv          nb112nf_jxH1,           192
.equiv          nb112nf_jyH1,           208
.equiv          nb112nf_jzH1,           224
.equiv          nb112nf_jxH2,           240
.equiv          nb112nf_jyH2,           256
.equiv          nb112nf_jzH2,           272
.equiv          nb112nf_qqOO,           288
.equiv          nb112nf_qqOH,           304
.equiv          nb112nf_qqHH,           320
.equiv          nb112nf_c6,             336
.equiv          nb112nf_c12,            352
.equiv          nb112nf_vctot,          368
.equiv          nb112nf_Vvdwtot,        384
.equiv          nb112nf_half,           400
.equiv          nb112nf_three,          416
.equiv          nb112nf_rsqOO,          432
.equiv          nb112nf_rsqOH1,         448
.equiv          nb112nf_rsqOH2,         464
.equiv          nb112nf_rsqH1O,         480
.equiv          nb112nf_rsqH1H1,        496
.equiv          nb112nf_rsqH1H2,        512
.equiv          nb112nf_rsqH2O,         528
.equiv          nb112nf_rsqH2H1,        544
.equiv          nb112nf_rsqH2H2,        560
.equiv          nb112nf_rinvOO,         576
.equiv          nb112nf_rinvOH1,        592
.equiv          nb112nf_rinvOH2,        608
.equiv          nb112nf_rinvH1O,        624
.equiv          nb112nf_rinvH1H1,       640
.equiv          nb112nf_rinvH1H2,       656
.equiv          nb112nf_rinvH2O,        672
.equiv          nb112nf_rinvH2H1,       688
.equiv          nb112nf_rinvH2H2,       704
.equiv          nb112nf_is3,            720
.equiv          nb112nf_ii3,            724
.equiv          nb112nf_nri,            728
.equiv          nb112nf_iinr,           736
.equiv          nb112nf_jindex,         744
.equiv          nb112nf_jjnr,           752
.equiv          nb112nf_shift,          760
.equiv          nb112nf_shiftvec,       768
.equiv          nb112nf_facel,          776
.equiv          nb112nf_innerjjnr,      784
.equiv          nb112nf_innerk,         792
.equiv          nb112nf_n,              796
.equiv          nb112nf_nn1,            800
.equiv          nb112nf_nouter,         804
.equiv          nb112nf_ninner,         808

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
	sub rsp, 816		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb112nf_nouter], eax
	mov [rsp + nb112nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb112nf_nri], edi
	mov [rsp + nb112nf_iinr], rsi
	mov [rsp + nb112nf_jindex], rdx
	mov [rsp + nb112nf_jjnr], rcx
	mov [rsp + nb112nf_shift], r8
	mov [rsp + nb112nf_shiftvec], r9
	mov rsi, [rbp + nb112nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb112nf_facel], xmm0


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb112nf_half], eax
	mov [rsp + nb112nf_half+4], ebx
	movsd xmm1, [rsp + nb112nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb112nf_half], xmm1
	movapd [rsp + nb112nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb112nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb112nf_charge]
	movsd xmm3, [rdx + rbx*8]	
	movsd xmm4, xmm3	
	movsd xmm5, [rdx + rbx*8 + 8]	
	mov rsi, [rbp + nb112nf_p_facel]
	movsd xmm0, [rsi]
	movsd xmm6, [rsp + nb112nf_facel]
	mulsd  xmm3, xmm3
	mulsd  xmm4, xmm5
	mulsd  xmm5, xmm5
	mulsd  xmm3, xmm6
	mulsd  xmm4, xmm6
	mulsd  xmm5, xmm6
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb112nf_qqOO], xmm3
	movapd [rsp + nb112nf_qqOH], xmm4
	movapd [rsp + nb112nf_qqHH], xmm5
	
	xorpd xmm0, xmm0
	mov   rdx, [rbp + nb112nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb112nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   rax, [rbp + nb112nf_vdwparam]
	movlpd xmm0, [rax + rdx*8]
	movlpd xmm1, [rax + rdx*8 + 8]
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0	
	movapd [rsp + nb112nf_c6], xmm0
	movapd [rsp + nb112nf_c12], xmm1

.nb112nf_threadloop:
        mov   rsi, [rbp + nb112nf_count]        ;# pointer to sync counter
        mov   eax, [rsi]
.nb112nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb112nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb112nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb112nf_n], eax
        mov [rsp + nb112nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb112nf_outerstart
        jmp .nb112nf_end

.nb112nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb112nf_nouter]
	mov [rsp + nb112nf_nouter], ebx

.nb112nf_outer:
	mov   rax, [rsp + nb112nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 

	mov   rax, [rsp + nb112nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb112nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb112nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb112nf_ii3], ebx	
	
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2
	addsd xmm3, [rax + rbx*8]
	addsd xmm4, [rax + rbx*8 + 8]
	addsd xmm5, [rax + rbx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb112nf_ixO], xmm3
	movapd [rsp + nb112nf_iyO], xmm4
	movapd [rsp + nb112nf_izO], xmm5

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
	movapd [rsp + nb112nf_ixH1], xmm0
	movapd [rsp + nb112nf_iyH1], xmm1
	movapd [rsp + nb112nf_izH1], xmm2
	movapd [rsp + nb112nf_ixH2], xmm3
	movapd [rsp + nb112nf_iyH2], xmm4
	movapd [rsp + nb112nf_izH2], xmm5

	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb112nf_vctot], xmm4
	movapd [rsp + nb112nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb112nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb112nf_pos]
	mov   rax, [rsp + nb112nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb112nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb112nf_ninner]
	mov   [rsp + nb112nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb112nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb112nf_unroll_loop
	jmp   .nb112nf_checksingle
.nb112nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb112nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	
	add qword ptr [rsp + nb112nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb112nf_pos]       ;# base of pos[] 

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
	movapd 	[rsp + nb112nf_jxO], xmm2
	movapd 	[rsp + nb112nf_jyO], xmm3
	movapd 	[rsp + nb112nf_jzO], xmm4
	movapd 	[rsp + nb112nf_jxH1], xmm5
	movapd 	[rsp + nb112nf_jyH1], xmm6
	movapd 	[rsp + nb112nf_jzH1], xmm7
	movlpd xmm2, [rsi + rax*8 + 48]
	movlpd xmm3, [rsi + rax*8 + 56]
	movlpd xmm4, [rsi + rax*8 + 64]
	movhpd xmm2, [rsi + rbx*8 + 48]
	movhpd xmm3, [rsi + rbx*8 + 56]
	movhpd xmm4, [rsi + rbx*8 + 64]
	movapd 	[rsp + nb112nf_jxH2], xmm2
	movapd 	[rsp + nb112nf_jyH2], xmm3
	movapd 	[rsp + nb112nf_jzH2], xmm4
	
	movapd xmm0, [rsp + nb112nf_ixO]
	movapd xmm1, [rsp + nb112nf_iyO]
	movapd xmm2, [rsp + nb112nf_izO]
	movapd xmm3, [rsp + nb112nf_ixO]
	movapd xmm4, [rsp + nb112nf_iyO]
	movapd xmm5, [rsp + nb112nf_izO]
	subpd  xmm0, [rsp + nb112nf_jxO]
	subpd  xmm1, [rsp + nb112nf_jyO]
	subpd  xmm2, [rsp + nb112nf_jzO]
	subpd  xmm3, [rsp + nb112nf_jxH1]
	subpd  xmm4, [rsp + nb112nf_jyH1]
	subpd  xmm5, [rsp + nb112nf_jzH1]
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
	movapd [rsp + nb112nf_rsqOO], xmm0
	movapd [rsp + nb112nf_rsqOH1], xmm3

	movapd xmm0, [rsp + nb112nf_ixO]
	movapd xmm1, [rsp + nb112nf_iyO]
	movapd xmm2, [rsp + nb112nf_izO]
	movapd xmm3, [rsp + nb112nf_ixH1]
	movapd xmm4, [rsp + nb112nf_iyH1]
	movapd xmm5, [rsp + nb112nf_izH1]
	subpd  xmm0, [rsp + nb112nf_jxH2]
	subpd  xmm1, [rsp + nb112nf_jyH2]
	subpd  xmm2, [rsp + nb112nf_jzH2]
	subpd  xmm3, [rsp + nb112nf_jxO]
	subpd  xmm4, [rsp + nb112nf_jyO]
	subpd  xmm5, [rsp + nb112nf_jzO]
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
	movapd [rsp + nb112nf_rsqOH2], xmm0
	movapd [rsp + nb112nf_rsqH1O], xmm3

	movapd xmm0, [rsp + nb112nf_ixH1]
	movapd xmm1, [rsp + nb112nf_iyH1]
	movapd xmm2, [rsp + nb112nf_izH1]
	movapd xmm3, [rsp + nb112nf_ixH1]
	movapd xmm4, [rsp + nb112nf_iyH1]
	movapd xmm5, [rsp + nb112nf_izH1]
	subpd  xmm0, [rsp + nb112nf_jxH1]
	subpd  xmm1, [rsp + nb112nf_jyH1]
	subpd  xmm2, [rsp + nb112nf_jzH1]
	subpd  xmm3, [rsp + nb112nf_jxH2]
	subpd  xmm4, [rsp + nb112nf_jyH2]
	subpd  xmm5, [rsp + nb112nf_jzH2]
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
	movapd [rsp + nb112nf_rsqH1H1], xmm0
	movapd [rsp + nb112nf_rsqH1H2], xmm3

	movapd xmm0, [rsp + nb112nf_ixH2]
	movapd xmm1, [rsp + nb112nf_iyH2]
	movapd xmm2, [rsp + nb112nf_izH2]
	movapd xmm3, [rsp + nb112nf_ixH2]
	movapd xmm4, [rsp + nb112nf_iyH2]
	movapd xmm5, [rsp + nb112nf_izH2]
	subpd  xmm0, [rsp + nb112nf_jxO]
	subpd  xmm1, [rsp + nb112nf_jyO]
	subpd  xmm2, [rsp + nb112nf_jzO]
	subpd  xmm3, [rsp + nb112nf_jxH1]
	subpd  xmm4, [rsp + nb112nf_jyH1]
	subpd  xmm5, [rsp + nb112nf_jzH1]
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
	movapd [rsp + nb112nf_rsqH2O], xmm0
	movapd [rsp + nb112nf_rsqH2H1], xmm4

	movapd xmm0, [rsp + nb112nf_ixH2]
	movapd xmm1, [rsp + nb112nf_iyH2]
	movapd xmm2, [rsp + nb112nf_izH2]
	subpd  xmm0, [rsp + nb112nf_jxH2]
	subpd  xmm1, [rsp + nb112nf_jyH2]
	subpd  xmm2, [rsp + nb112nf_jzH2]
	mulpd xmm0, xmm0
	mulpd xmm1, xmm1
	mulpd xmm2, xmm2
	addpd xmm0, xmm1
	addpd xmm0, xmm2
	movapd [rsp + nb112nf_rsqH2H2], xmm0
	
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
	movapd  xmm3, [rsp + nb112nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb112nf_half] ;# iter1 
	mulpd   xmm7, [rsp + nb112nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb112nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb112nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb112nf_half] ;# rinv 
	movapd [rsp + nb112nf_rinvH2H2], xmm1
	movapd [rsp + nb112nf_rinvH2H1], xmm5

	movapd xmm0, [rsp + nb112nf_rsqOO]
	movapd xmm4, [rsp + nb112nf_rsqOH1]	
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
	movapd  xmm3, [rsp + nb112nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb112nf_half] ;# iter1 of  
	mulpd   xmm7, [rsp + nb112nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb112nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb112nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb112nf_half] ;# rinv
	movapd [rsp + nb112nf_rinvOO], xmm1
	movapd [rsp + nb112nf_rinvOH1], xmm5

	movapd xmm0, [rsp + nb112nf_rsqOH2]
	movapd xmm4, [rsp + nb112nf_rsqH1O]	
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
	movapd  xmm3, [rsp + nb112nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb112nf_half] ;# iter1 
	mulpd   xmm7, [rsp + nb112nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb112nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb112nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb112nf_half] ;# rinv 
	movapd [rsp + nb112nf_rinvOH2], xmm1
	movapd [rsp + nb112nf_rinvH1O], xmm5

	movapd xmm0, [rsp + nb112nf_rsqH1H1]
	movapd xmm4, [rsp + nb112nf_rsqH1H2]	
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
	movapd  xmm3, [rsp + nb112nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	mulpd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subpd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm3, [rsp + nb112nf_half] ;# iter1a 
	mulpd   xmm7, [rsp + nb112nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulpd   xmm3, xmm3	;# luA*luA 
	mulpd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb112nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	mulpd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subpd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulpd   xmm1, [rsp + nb112nf_half] ;# rinv 
	mulpd   xmm5, [rsp + nb112nf_half] ;# rinv 
	movapd [rsp + nb112nf_rinvH1H1], xmm1
	movapd [rsp + nb112nf_rinvH1H2], xmm5

	movapd xmm0, [rsp + nb112nf_rsqH2O]
	cvtpd2ps xmm1, xmm0	
	rsqrtps xmm1, xmm1
	cvtps2pd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulpd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [rsp + nb112nf_three]
	mulpd   xmm1, xmm0	;# rsqA*luA*luA 
	subpd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulpd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm3, [rsp + nb112nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulpd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [rsp + nb112nf_three]
	mulpd   xmm3, xmm0	;# rsqA*luA*luA 
	subpd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulpd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulpd   xmm1, [rsp + nb112nf_half] ;# rinv 
	movapd [rsp + nb112nf_rinvH2O], xmm1

	;# start with OO interaction 
	movapd xmm0, [rsp + nb112nf_rinvOO]
	movapd xmm7, xmm0
	mulpd  xmm0, xmm0
	movapd xmm1, xmm0
	mulpd  xmm1, xmm0
	mulpd  xmm1, xmm0	;# xmm1=rinvsix 
	mulpd  xmm7, [rsp + nb112nf_qqOO]
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulpd  xmm1, [rsp + nb112nf_c6]	
	mulpd  xmm2, [rsp + nb112nf_c12]	
	movapd xmm3, xmm2
	subpd  xmm3, xmm1	;# xmm3=Vvdw12-Vvdw6 
	addpd  xmm3, [rsp + nb112nf_Vvdwtot]
	movapd [rsp + nb112nf_Vvdwtot], xmm3
	addpd  xmm7, [rsp + nb112nf_vctot]

	;# other interactions 
	movapd xmm1, [rsp + nb112nf_rinvOH1]
	movapd xmm2, [rsp + nb112nf_rinvH1H1]
	
	addpd xmm1, [rsp + nb112nf_rinvOH2]
	addpd xmm2, [rsp + nb112nf_rinvH1H2]
	
	addpd xmm1, [rsp + nb112nf_rinvH1O]
	addpd xmm2, [rsp + nb112nf_rinvH2H1]

	addpd xmm1, [rsp + nb112nf_rinvH2O]
	addpd xmm2, [rsp + nb112nf_rinvH2H2]

	mulpd xmm1, [rsp + nb112nf_qqOH]
	mulpd xmm2, [rsp + nb112nf_qqHH]
	
	addpd xmm7, xmm1	
	addpd xmm7, xmm2

	movapd [rsp + nb112nf_vctot], xmm7
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb112nf_innerk],  2
	jl    .nb112nf_checksingle
	jmp   .nb112nf_unroll_loop
.nb112nf_checksingle:
	mov   edx, [rsp + nb112nf_innerk]
	and   edx, 1
	jnz   .nb112nf_dosingle
	jmp   .nb112nf_updateouterdata
.nb112nf_dosingle:
	mov   rdx, [rsp + nb112nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb112nf_innerjjnr],  4	

	mov rsi, [rbp + nb112nf_pos]
	lea   rax, [rax + rax*2]  

	;# fetch j coordinates 
	movlpd xmm2, [rsi + rax*8]
	movlpd xmm3, [rsi + rax*8 + 8]
	movlpd xmm4, [rsi + rax*8 + 16]
	movlpd xmm5, [rsi + rax*8 + 24]
	movlpd xmm6, [rsi + rax*8 + 32]
	movlpd xmm7, [rsi + rax*8 + 40]
	movapd 	[rsp + nb112nf_jxO], xmm2
	movapd 	[rsp + nb112nf_jyO], xmm3
	movapd 	[rsp + nb112nf_jzO], xmm4
	movapd 	[rsp + nb112nf_jxH1], xmm5
	movapd 	[rsp + nb112nf_jyH1], xmm6
	movapd 	[rsp + nb112nf_jzH1], xmm7
	movlpd xmm2, [rsi + rax*8 + 48]
	movlpd xmm3, [rsi + rax*8 + 56]
	movlpd xmm4, [rsi + rax*8 + 64]
	movapd 	[rsp + nb112nf_jxH2], xmm2
	movapd 	[rsp + nb112nf_jyH2], xmm3
	movapd 	[rsp + nb112nf_jzH2], xmm4
	
	movapd xmm0, [rsp + nb112nf_ixO]
	movapd xmm1, [rsp + nb112nf_iyO]
	movapd xmm2, [rsp + nb112nf_izO]
	movapd xmm3, [rsp + nb112nf_ixO]
	movapd xmm4, [rsp + nb112nf_iyO]
	movapd xmm5, [rsp + nb112nf_izO]
	subsd  xmm0, [rsp + nb112nf_jxO]
	subsd  xmm1, [rsp + nb112nf_jyO]
	subsd  xmm2, [rsp + nb112nf_jzO]
	subsd  xmm3, [rsp + nb112nf_jxH1]
	subsd  xmm4, [rsp + nb112nf_jyH1]
	subsd  xmm5, [rsp + nb112nf_jzH1]
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
	movapd [rsp + nb112nf_rsqOO], xmm0
	movapd [rsp + nb112nf_rsqOH1], xmm3

	movapd xmm0, [rsp + nb112nf_ixO]
	movapd xmm1, [rsp + nb112nf_iyO]
	movapd xmm2, [rsp + nb112nf_izO]
	movapd xmm3, [rsp + nb112nf_ixH1]
	movapd xmm4, [rsp + nb112nf_iyH1]
	movapd xmm5, [rsp + nb112nf_izH1]
	subsd  xmm0, [rsp + nb112nf_jxH2]
	subsd  xmm1, [rsp + nb112nf_jyH2]
	subsd  xmm2, [rsp + nb112nf_jzH2]
	subsd  xmm3, [rsp + nb112nf_jxO]
	subsd  xmm4, [rsp + nb112nf_jyO]
	subsd  xmm5, [rsp + nb112nf_jzO]
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
	movapd [rsp + nb112nf_rsqOH2], xmm0
	movapd [rsp + nb112nf_rsqH1O], xmm3

	movapd xmm0, [rsp + nb112nf_ixH1]
	movapd xmm1, [rsp + nb112nf_iyH1]
	movapd xmm2, [rsp + nb112nf_izH1]
	movapd xmm3, [rsp + nb112nf_ixH1]
	movapd xmm4, [rsp + nb112nf_iyH1]
	movapd xmm5, [rsp + nb112nf_izH1]
	subsd  xmm0, [rsp + nb112nf_jxH1]
	subsd  xmm1, [rsp + nb112nf_jyH1]
	subsd  xmm2, [rsp + nb112nf_jzH1]
	subsd  xmm3, [rsp + nb112nf_jxH2]
	subsd  xmm4, [rsp + nb112nf_jyH2]
	subsd  xmm5, [rsp + nb112nf_jzH2]
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
	movapd [rsp + nb112nf_rsqH1H1], xmm0
	movapd [rsp + nb112nf_rsqH1H2], xmm3

	movapd xmm0, [rsp + nb112nf_ixH2]
	movapd xmm1, [rsp + nb112nf_iyH2]
	movapd xmm2, [rsp + nb112nf_izH2]
	movapd xmm3, [rsp + nb112nf_ixH2]
	movapd xmm4, [rsp + nb112nf_iyH2]
	movapd xmm5, [rsp + nb112nf_izH2]
	subsd  xmm0, [rsp + nb112nf_jxO]
	subsd  xmm1, [rsp + nb112nf_jyO]
	subsd  xmm2, [rsp + nb112nf_jzO]
	subsd  xmm3, [rsp + nb112nf_jxH1]
	subsd  xmm4, [rsp + nb112nf_jyH1]
	subsd  xmm5, [rsp + nb112nf_jzH1]
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
	movapd [rsp + nb112nf_rsqH2O], xmm0
	movapd [rsp + nb112nf_rsqH2H1], xmm4

	movapd xmm0, [rsp + nb112nf_ixH2]
	movapd xmm1, [rsp + nb112nf_iyH2]
	movapd xmm2, [rsp + nb112nf_izH2]
	subsd  xmm0, [rsp + nb112nf_jxH2]
	subsd  xmm1, [rsp + nb112nf_jyH2]
	subsd  xmm2, [rsp + nb112nf_jzH2]
	mulsd xmm0, xmm0
	mulsd xmm1, xmm1
	mulsd xmm2, xmm2
	addsd xmm0, xmm1
	addsd xmm0, xmm2
	movapd [rsp + nb112nf_rsqH2H2], xmm0
	
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
	movapd  xmm3, [rsp + nb112nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb112nf_half] ;# iter1 
	mulsd   xmm7, [rsp + nb112nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb112nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb112nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb112nf_half] ;# rinv 
	movapd [rsp + nb112nf_rinvH2H2], xmm1
	movapd [rsp + nb112nf_rinvH2H1], xmm5

	movapd xmm0, [rsp + nb112nf_rsqOO]
	movapd xmm4, [rsp + nb112nf_rsqOH1]	
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
	movapd  xmm3, [rsp + nb112nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb112nf_half] ;# iter1 of  
	mulsd   xmm7, [rsp + nb112nf_half] ;# iter1 of  

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb112nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb112nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb112nf_half] ;# rinv
	movapd [rsp + nb112nf_rinvOO], xmm1
	movapd [rsp + nb112nf_rinvOH1], xmm5

	movapd xmm0, [rsp + nb112nf_rsqOH2]
	movapd xmm4, [rsp + nb112nf_rsqH1O]	
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
	movapd  xmm3, [rsp + nb112nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb112nf_half] ;# iter1 
	mulsd   xmm7, [rsp + nb112nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb112nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb112nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb112nf_half] ;# rinv 
	movapd [rsp + nb112nf_rinvOH2], xmm1
	movapd [rsp + nb112nf_rinvH1O], xmm5

	movapd xmm0, [rsp + nb112nf_rsqH1H1]
	movapd xmm4, [rsp + nb112nf_rsqH1H2]	
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
	movapd  xmm3, [rsp + nb112nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	mulsd   xmm5, xmm4	;# rsqB*luB*luB 	
	movapd  xmm7, xmm3
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	subsd   xmm7, xmm5	;# 3-rsqB*luB*luB 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm7, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm3, [rsp + nb112nf_half] ;# iter1a 
	mulsd   xmm7, [rsp + nb112nf_half] ;# iter1b 

	movapd  xmm2, xmm3	;# copy of luA 
	movapd  xmm6, xmm7	;# copy of luB 
	mulsd   xmm3, xmm3	;# luA*luA 
	mulsd   xmm7, xmm7	;# luB*luB 
	movapd  xmm1, [rsp + nb112nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	mulsd   xmm7, xmm4	;# rsqB*luB*luB 	
	movapd  xmm5, xmm1
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	subsd   xmm5, xmm7	;# 3-rsqB*luB*luB 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm5, xmm6	;# luB*(3-rsqB*luB*luB) 
	mulsd   xmm1, [rsp + nb112nf_half] ;# rinv 
	mulsd   xmm5, [rsp + nb112nf_half] ;# rinv 
	movapd [rsp + nb112nf_rinvH1H1], xmm1
	movapd [rsp + nb112nf_rinvH1H2], xmm5

	movapd xmm0, [rsp + nb112nf_rsqH2O]
	cvtsd2ss xmm1, xmm0	
	rsqrtss xmm1, xmm1
	cvtss2sd xmm1, xmm1
	
	movapd  xmm2, xmm1	;# copy of luA 
	mulsd   xmm1, xmm1	;# luA*luA 
	movapd  xmm3, [rsp + nb112nf_three]
	mulsd   xmm1, xmm0	;# rsqA*luA*luA 
	subsd   xmm3, xmm1	;# 3-rsqA*luA*luA 
	mulsd   xmm3, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm3, [rsp + nb112nf_half] ;# iter1 

	movapd  xmm2, xmm3	;# copy of luA 
	mulsd   xmm3, xmm3	;# luA*luA 
	movapd  xmm1, [rsp + nb112nf_three]
	mulsd   xmm3, xmm0	;# rsqA*luA*luA 
	subsd   xmm1, xmm3	;# 3-rsqA*luA*luA 
	mulsd   xmm1, xmm2	;# luA*(3-rsqA*luA*luA) 
	mulsd   xmm1, [rsp + nb112nf_half] ;# rinv 
	movapd [rsp + nb112nf_rinvH2O], xmm1

	;# start with OO interaction 
	movapd xmm0, [rsp + nb112nf_rinvOO]
	movapd xmm7, xmm0
	mulsd  xmm0, xmm0
	movapd xmm1, xmm0
	mulsd  xmm1, xmm0
	mulsd  xmm1, xmm0	;# xmm1=rinvsix 
	mulsd  xmm7, [rsp + nb112nf_qqOO]
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm1, [rsp + nb112nf_c6]	
	mulsd  xmm2, [rsp + nb112nf_c12]	
	movapd xmm3, xmm2
	subsd  xmm3, xmm1	;# xmm3=Vvdw12-Vvdw6 
	addsd  xmm3, [rsp + nb112nf_Vvdwtot]
	movlpd [rsp + nb112nf_Vvdwtot], xmm3
	addsd  xmm7, [rsp + nb112nf_vctot]

	;# other interactions 
	movapd xmm1, [rsp + nb112nf_rinvOH1]
	movapd xmm2, [rsp + nb112nf_rinvH1H1]
	
	addsd xmm1, [rsp + nb112nf_rinvOH2]
	addsd xmm2, [rsp + nb112nf_rinvH1H2]
	
	addsd xmm1, [rsp + nb112nf_rinvH1O]
	addsd xmm2, [rsp + nb112nf_rinvH2H1]

	addsd xmm1, [rsp + nb112nf_rinvH2O]
	addsd xmm2, [rsp + nb112nf_rinvH2H2]

	mulsd xmm1, [rsp + nb112nf_qqOH]
	mulsd xmm2, [rsp + nb112nf_qqHH]
	
	addsd xmm7, xmm1	
	addsd xmm7, xmm2

	movlpd [rsp + nb112nf_vctot], xmm7

.nb112nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb112nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb112nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb112nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb112nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb112nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   rax, [rbp + nb112nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb112nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb112nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb112nf_n], esi
        jmp .nb112nf_outer
.nb112nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb112nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb112nf_end
        ;# non-zero, do one more workunit
        jmp   .nb112nf_threadloop
.nb112nf_end:
	mov eax, [rsp + nb112nf_nouter]
	mov ebx, [rsp + nb112nf_ninner]
	mov rcx, [rbp + nb112nf_outeriter]
	mov rdx, [rbp + nb112nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 816
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
