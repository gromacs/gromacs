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


.globl nb_kernel332_x86_64_sse
.globl _nb_kernel332_x86_64_sse
nb_kernel332_x86_64_sse:	
_nb_kernel332_x86_64_sse:	
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
	;# bottom of stack is cache-aligned for sse use 
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
.equiv          nb332_fjxO,             1008
.equiv          nb332_fjyO,             1024
.equiv          nb332_fjzO,             1040
.equiv          nb332_fjxH1,            1056
.equiv          nb332_fjyH1,            1072
.equiv          nb332_fjzH1,            1088
.equiv          nb332_fjxH2,            1104
.equiv          nb332_fjyH2,            1120
.equiv          nb332_fjzH2,            1136
.equiv          nb332_half,             1152
.equiv          nb332_three,            1168
.equiv          nb332_epsO,             1184
.equiv          nb332_epsH1,            1200
.equiv          nb332_epsH2,            1216
.equiv          nb332_rsqH1O,           1232
.equiv          nb332_rsqOO,            1248
.equiv          nb332_rsqH1H2,          1264
.equiv          nb332_rsqH2O,           1280
.equiv          nb332_fstmpH1,          1296
.equiv          nb332_fstmpH2,          1312
.equiv          nb332_rinvOO,           1328
.equiv          nb332_rinvOH1,          1344
.equiv          nb332_rinvOH2,          1360
.equiv          nb332_rinvH1O,          1376
.equiv          nb332_rinvH1H1,         1392
.equiv          nb332_rinvH1H2,         1408
.equiv          nb332_rinvH2O,          1424
.equiv          nb332_rinvH2H1,         1440
.equiv          nb332_rinvH2H2,         1456
.equiv          nb332_fstmp,            1472
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
    
	emms

	sub rsp, 1584		;# local variable stack space (n*16+8)
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
	movss xmm0, [rsi]
	movss [rsp + nb332_facel], xmm0

	mov rax, [rbp + nb332_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb332_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb332_half], eax
	movss xmm1, [rsp + nb332_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb332_half],  xmm1
	movaps [rsp + nb332_two],  xmm2
	movaps [rsp + nb332_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb332_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb332_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, xmm3	
	movss xmm5, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb332_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb332_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb332_qqOO], xmm3
	movaps [rsp + nb332_qqOH], xmm4
	movaps [rsp + nb332_qqHH], xmm5
	
	xorps xmm0, xmm0
	mov   rdx, [rbp + nb332_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb332_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   rax, [rbp + nb332_vdwparam]
	movlps xmm0, [rax + rdx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 85  ;# 01010101
	movaps [rsp + nb332_c6], xmm0
	movaps [rsp + nb332_c12], xmm1

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
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb332_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb332_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb332_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb332_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb332_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb332_ixO], xmm3
	movaps [rsp + nb332_iyO], xmm4
	movaps [rsp + nb332_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 12]
	addss xmm1, [rax + rbx*4 + 16]
	addss xmm2, [rax + rbx*4 + 20]		
	addss xmm3, [rax + rbx*4 + 24]
	addss xmm4, [rax + rbx*4 + 28]
	addss xmm5, [rax + rbx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb332_ixH1], xmm0
	movaps [rsp + nb332_iyH1], xmm1
	movaps [rsp + nb332_izH1], xmm2
	movaps [rsp + nb332_ixH2], xmm3
	movaps [rsp + nb332_iyH2], xmm4
	movaps [rsp + nb332_izH2], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb332_vctot], xmm4
	movaps [rsp + nb332_Vvdwtot], xmm4
	movaps [rsp + nb332_fixO], xmm4
	movaps [rsp + nb332_fiyO], xmm4
	movaps [rsp + nb332_fizO], xmm4
	movaps [rsp + nb332_fixH1], xmm4
	movaps [rsp + nb332_fiyH1], xmm4
	movaps [rsp + nb332_fizH1], xmm4
	movaps [rsp + nb332_fixH2], xmm4
	movaps [rsp + nb332_fiyH2], xmm4
	movaps [rsp + nb332_fizH2], xmm4
	
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
	sub   edx,  4
	add   ecx, [rsp + nb332_ninner]
	mov   [rsp + nb332_ninner], ecx
	add   edx, 0
	mov   [rsp + nb332_innerk], edx    ;# number of innerloop atoms 
	jge   .nb332_unroll_loop
	jmp   .nb332_single_check
.nb332_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb332_innerjjnr]     ;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb332_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb332_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2]     ;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	

	;# move j O coordinates to local temp variables 
    movlps xmm0, [rsi + rax*4] ;# jxOa jyOa  -   -
    movlps xmm1, [rsi + rcx*4] ;# jxOc jyOc  -   -
    movhps xmm0, [rsi + rbx*4] ;# jxOa jyOa jxOb jyOb 
    movhps xmm1, [rsi + rdx*4] ;# jxOc jyOc jxOd jyOd 

    movss  xmm2, [rsi + rax*4 + 8] ;# jzOa  -  -  -
    movss  xmm3, [rsi + rcx*4 + 8] ;# jzOc  -  -  -
    movhps xmm2, [rsi + rbx*4 + 8] ;# jzOa  -  jzOb  -
    movhps xmm3, [rsi + rdx*4 + 8] ;# jzOc  -  jzOd -
    
    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxOa jxOc jyOa jyOc        
    unpckhps xmm4, xmm1  ;# jxOb jxOd jyOb jyOd
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps xmm2, xmm3,  136  ;# 10001000 => jzOa jzOb jzOc jzOd

    ;# xmm0 = Ox
    ;# xmm1 = Oy
    ;# xmm2 = Oz
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb332_ixO]
    subps xmm1, [rsp + nb332_iyO]
    subps xmm2, [rsp + nb332_izO]
    subps xmm3, [rsp + nb332_ixH1]
    subps xmm4, [rsp + nb332_iyH1]
    subps xmm5, [rsp + nb332_izH1]
    subps xmm6, [rsp + nb332_ixH2]
    subps xmm7, [rsp + nb332_iyH2]
    subps xmm8, [rsp + nb332_izH2]
    
    movd mm0, eax ;# save j3 in mm0-mm3
    movd mm1, ebx
    movd mm2, ecx
    movd mm3, edx

	movaps [rsp + nb332_dxOO], xmm0
	movaps [rsp + nb332_dyOO], xmm1
	movaps [rsp + nb332_dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb332_dxH1O], xmm3
	movaps [rsp + nb332_dyH1O], xmm4
	movaps [rsp + nb332_dzH1O], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb332_dxH2O], xmm6
	movaps [rsp + nb332_dyH2O], xmm7
	movaps [rsp + nb332_dzH2O], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for jO atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb332_three]
	movaps  xmm10, xmm9
    movaps  xmm11, xmm9

	mulps   xmm1, xmm0 ;# rsq*lu*lu
	mulps   xmm4, xmm3 ;# rsq*lu*lu 
    mulps   xmm7, xmm6 ;# rsq*lu*lu
	
	subps   xmm9, xmm1
	subps   xmm10, xmm4
    subps   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulps   xmm9, xmm2
	mulps   xmm10, xmm5
    mulps   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movaps  xmm4, [rsp + nb332_half]
	mulps   xmm9, xmm4  ;# rinvOO 
	mulps   xmm10, xmm4 ;# rinvH1O
    mulps   xmm11, xmm4 ;# rinvH2O

	movaps  [rsp + nb332_rinvOO], xmm9
	movaps  [rsp + nb332_rinvH1O], xmm10
	movaps  [rsp + nb332_rinvH2O], xmm11
	
	;# O interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11
    movaps [rsp + nb332_rinvOO], xmm9
    movaps xmm1, [rsp + nb332_tsc]
    
    mulps  xmm0, xmm9  ;# r
    mulps  xmm3, xmm10
    mulps  xmm6, xmm11
    mulps  xmm0, xmm1 ;# rtab
    mulps  xmm3, xmm1
    mulps  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttps2dq xmm1, xmm0
    cvttps2dq xmm4, xmm3
    cvttps2dq xmm7, xmm6        

    ;# convert back to float
    cvtdq2ps  xmm2, xmm1
    cvtdq2ps  xmm5, xmm4
    cvtdq2ps  xmm8, xmm7
    
    ;# multiply by 4
    pslld   xmm1, 2
    pslld   xmm4, 2
    pslld   xmm7, 2
    
    ;# multiply by three (copy, mult. by two, add back)
    movaps  xmm10, xmm1
    movaps  xmm11, xmm4
    movaps  xmm12, xmm7
    pslld   xmm1, 1
    pslld   xmm4, 1
    pslld   xmm7, 1    
    paddd   xmm1, xmm10
    paddd   xmm4, xmm11
    paddd   xmm7, xmm12    
    
    ;# move to integer registers
    movhlps xmm13, xmm1
    movhlps xmm14, xmm4
    movhlps xmm15, xmm7
    movd    eax, xmm1
    movd    r8d, xmm4
    movd    r12d, xmm7
    movd    ecx, xmm13
    movd    r10d, xmm14
    movd    r14d, xmm15
    pshufd  xmm1, xmm1, 1
    pshufd  xmm4, xmm4, 1
    pshufd  xmm7, xmm7, 1
    pshufd  xmm13, xmm13, 1
    pshufd  xmm14, xmm14, 1
    pshufd  xmm15, xmm15, 1
    movd    ebx, xmm1
    movd    r9d, xmm4
    movd    r13d, xmm7    
    movd    edx, xmm13
    movd    r11d, xmm14
    movd    r15d, xmm15   
        
    mov  rsi, [rbp + nb332_VFtab]

    ;# calculate eps
    subps     xmm0, xmm2
    subps     xmm3, xmm5
    subps     xmm6, xmm8

    movaps    [rsp + nb332_epsO], xmm0
    movaps    [rsp + nb332_epsH1], xmm3
    movaps    [rsp + nb332_epsH2], xmm6

    ;# Load LOTS of table data
   	movlps xmm1, [rsi + rax*4]
   	movlps xmm5, [rsi + r8*4]
   	movlps xmm9, [rsi + r12*4]

	movlps xmm3, [rsi + rcx*4]
	movlps xmm7, [rsi + r10*4]
	movlps xmm11, [rsi + r14*4]

	movhps xmm1, [rsi + rbx*4]
	movhps xmm5, [rsi + r9*4]
	movhps xmm9, [rsi + r13*4]

	movhps xmm3, [rsi + rdx*4]
	movhps xmm7, [rsi + r11*4]
	movhps xmm11, [rsi + r15*4]

    movaps xmm0, xmm1
    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm0, xmm3, 136  ;# 10001000
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm8, xmm11, 136  ;# 10001000
	shufps xmm1, xmm3, 221  ;# 11011101
	shufps xmm5, xmm7, 221  ;# 11011101
	shufps xmm9, xmm11, 221  ;# 11011101
    
	movlps xmm3, [rsi + rax*4 + 8]
	movlps xmm7, [rsi + r8*4 + 8]
	movlps xmm11, [rsi + r12*4 + 8]
    
	movlps xmm12, [rsi + rcx*4 + 8]
	movlps xmm13, [rsi + r10*4 + 8]
	movlps xmm14, [rsi + r14*4 + 8]

	movhps xmm3, [rsi + rbx*4 + 8]
	movhps xmm7, [rsi + r9*4 + 8]
	movhps xmm11, [rsi + r13*4 + 8]
    
	movhps xmm12, [rsi + rdx*4 + 8]
	movhps xmm13, [rsi + r11*4 + 8]
	movhps xmm14, [rsi + r15*4 + 8]

    movaps xmm2, xmm3
    movaps xmm6, xmm7
    movaps xmm10, xmm11
    
	shufps xmm2, xmm12, 136  ;# 10001000
	shufps xmm6, xmm13, 136  ;# 10001000
	shufps xmm10, xmm14, 136  ;# 10001000
	shufps xmm3, xmm12, 221  ;# 11011101
	shufps xmm7, xmm13, 221  ;# 11011101
	shufps xmm11, xmm14, 221  ;# 11011101
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    mulps  xmm3, [rsp + nb332_epsO]   ;# Heps
    mulps  xmm7, [rsp + nb332_epsH1]
    mulps  xmm11, [rsp + nb332_epsH2]
    mulps  xmm2, [rsp + nb332_epsO]   ;# Geps
    mulps  xmm6, [rsp + nb332_epsH1]
    mulps  xmm10, [rsp + nb332_epsH2]
    mulps  xmm3, [rsp + nb332_epsO]   ;# Heps2
    mulps  xmm7, [rsp + nb332_epsH1]
    mulps  xmm11, [rsp + nb332_epsH2]

    addps  xmm1, xmm2   ;# F+Geps
    addps  xmm5, xmm6
    addps  xmm9, xmm10 
    addps  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addps  xmm5, xmm7
    addps  xmm9, xmm11 
    addps  xmm3, xmm3    ;# 2*Heps2
    addps  xmm7, xmm7
    addps  xmm11, xmm11
    addps  xmm3, xmm2    ;# 2*Heps2+Geps
    addps  xmm7, xmm6  
    addps  xmm11, xmm10
    addps  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addps  xmm7, xmm5
    addps  xmm11, xmm9
    mulps  xmm1, [rsp + nb332_epsO]   ;# eps*Fp
    mulps  xmm5, [rsp + nb332_epsH1]
    mulps  xmm9, [rsp + nb332_epsH2]
    addps  xmm1, xmm0     ;# VV
    addps  xmm5, xmm4
    addps  xmm9, xmm8
    mulps  xmm1, [rsp + nb332_qqOO]   ;# VV*qq = vcoul
    mulps  xmm5, [rsp + nb332_qqOH]
    mulps  xmm9, [rsp + nb332_qqOH]
    mulps  xmm3, [rsp + nb332_qqOO]    ;# FF*qq = fij
    mulps  xmm7, [rsp + nb332_qqOH]
    mulps  xmm11, [rsp + nb332_qqOH] 

    ;# accumulate vctot
    addps  xmm1, [rsp + nb332_vctot]
    addps  xmm5, xmm9
    addps  xmm1, xmm5
    movaps [rsp + nb332_vctot], xmm1

    movaps xmm2, xmm7
    movaps xmm1, xmm11

    ;# fij coul in xmm3, xmm2, xmm1    
    
    ;# calculate LJ table
    movlps xmm5, [rsi + rax*4 + 16]
   	movlps xmm9, [rsi + rax*4 + 32]

	movlps xmm7, [rsi + rcx*4 + 16]
	movlps xmm11, [rsi + rcx*4 + 32]

	movhps xmm5, [rsi + rbx*4 + 16]
	movhps xmm9, [rsi + rbx*4 + 32]

	movhps xmm7, [rsi + rdx*4 + 16]
	movhps xmm11, [rsi + rdx*4 + 32]

    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm8, xmm11, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101
	shufps xmm9, xmm11, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 24]
	movlps xmm11, [rsi + rax*4 + 40]
    
	movlps xmm13, [rsi + rcx*4 + 24]
	movlps xmm14, [rsi + rcx*4 + 40]

	movhps xmm7, [rsi + rbx*4 + 24]
	movhps xmm11, [rsi + rbx*4 + 40]
    
	movhps xmm13, [rsi + rdx*4 + 24]
	movhps xmm14, [rsi + rdx*4 + 40]

    movaps xmm6, xmm7
    movaps xmm10, xmm11
    
	shufps xmm6, xmm13, 136  ;# 10001000
	shufps xmm10, xmm14, 136  ;# 10001000
	shufps xmm7, xmm13, 221  ;# 11011101
	shufps xmm11, xmm14, 221  ;# 11011101
    ;# dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11
    movaps xmm0, [rsp + nb332_epsO]
    
    mulps  xmm7, xmm0    ;# Heps
    mulps  xmm11, xmm0 
    mulps  xmm6, xmm0   ;# Geps
    mulps  xmm10, xmm0 
    mulps  xmm7, xmm0   ;# Heps2
    mulps  xmm11, xmm0 
    addps  xmm5, xmm6  ;# F+Geps
    addps  xmm9, xmm10 
    addps  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addps  xmm9, xmm11 
    addps  xmm7, xmm7    ;# 2*Heps2
    addps  xmm11, xmm11
    addps  xmm7, xmm6   ;# 2*Heps2+Geps
    addps  xmm11, xmm10
    
    addps  xmm7, xmm5  ;# FF = Fp + 2*Heps2 + Geps
    addps  xmm11, xmm9
    mulps  xmm5, xmm0  ;# eps*Fp
    mulps  xmm9, xmm0
    movaps xmm12, [rsp + nb332_c6]
    movaps xmm13, [rsp + nb332_c12]
    addps  xmm5, xmm4 ;# VV
    addps  xmm9, xmm8

    movd eax, mm0 ;# restore j3 from mm0-mm3
    movd ebx, mm1
    movd ecx, mm2
    movd edx, mm3

    mulps  xmm5, xmm12  ;# VV*c6 = vnb6
    mulps  xmm9, xmm13  ;# VV*c12 = vnb12
    addps  xmm5, xmm9
    addps  xmm5, [rsp + nb332_Vvdwtot]
    movaps [rsp + nb332_Vvdwtot], xmm5
        
    mulps  xmm7, xmm12   ;# FF*c6 = fnb6
    mulps  xmm11, xmm13   ;# FF*c12  = fnb12
    addps  xmm7, xmm11
    
    addps  xmm3, xmm7
    movaps xmm10, [rsp + nb332_tsc]

    mulps  xmm3, xmm10  ;# fscal
    mulps  xmm2, xmm10
    mulps  xmm1, xmm10
        
	;# move j O forces to local temp variables 
    mov rdi, [rbp + nb332_faction]
    movlps xmm11, [rdi + rax*4] ;# jxOa jyOa  -   -
    movlps xmm12, [rdi + rcx*4] ;# jxOc jyOc  -   -

    xorps  xmm0, xmm0
    xorps  xmm4, xmm4
    xorps  xmm8, xmm8
    
    movhps xmm11, [rdi + rbx*4] ;# jxOa jyOa jxOb jyOb 
    movhps xmm12, [rdi + rdx*4] ;# jxOc jyOc jxOd jyOd 

    subps  xmm0, xmm3
    subps  xmm4, xmm2
    subps  xmm8, xmm1

    movss  xmm13, [rdi + rax*4 + 8] ;# jzOa  -  -  -
    movss  xmm14, [rdi + rcx*4 + 8] ;# jzOc  -  -  -
    mulps  xmm0, [rsp + nb332_rinvOO]
    mulps  xmm4, [rsp + nb332_rinvH1O]
    mulps  xmm8, [rsp + nb332_rinvH2O]
    movhps xmm13, [rdi + rbx*4 + 8] ;# jzOa  -  jzOb  -
    movhps xmm14, [rdi + rdx*4 + 8] ;# jzOc  -  jzOd -    
    shufps xmm13, xmm14,  136  ;# 10001000 => jzOa jzOb jzOc jzOd

    ;# xmm11: jxOa jyOa jxOb jyOb 
    ;# xmm12: jxOc jyOc jxOd jyOd
    ;# xmm13: jzOa jzOb jzOc jzOd

    movaps xmm1, xmm0
    movaps xmm2, xmm0
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm8
    movaps xmm7, xmm8

	mulps xmm0, [rsp + nb332_dxOO]
	mulps xmm1, [rsp + nb332_dyOO]
	mulps xmm2, [rsp + nb332_dzOO]
	mulps xmm3, [rsp + nb332_dxH1O]
	mulps xmm4, [rsp + nb332_dyH1O]
	mulps xmm5, [rsp + nb332_dzH1O]
	mulps xmm6, [rsp + nb332_dxH2O]
	mulps xmm7, [rsp + nb332_dyH2O]
	mulps xmm8, [rsp + nb332_dzH2O]

    movaps xmm14,  xmm0
    movaps xmm15, xmm1
    addps xmm13, xmm2
    addps xmm0, [rsp + nb332_fixO]
    addps xmm1, [rsp + nb332_fiyO]
    addps xmm2, [rsp + nb332_fizO]

    addps xmm14,  xmm3
    addps xmm15, xmm4
    addps xmm13, xmm5
    addps xmm3, [rsp + nb332_fixH1]
    addps xmm4, [rsp + nb332_fiyH1]
    addps xmm5, [rsp + nb332_fizH1]

    addps xmm14,  xmm6
    addps xmm15, xmm7
    addps xmm13, xmm8
    addps xmm6, [rsp + nb332_fixH2]
    addps xmm7, [rsp + nb332_fiyH2]
    addps xmm8, [rsp + nb332_fizH2]

    movaps [rsp + nb332_fixO], xmm0
    movaps [rsp + nb332_fiyO], xmm1
    movaps [rsp + nb332_fizO], xmm2
    movaps [rsp + nb332_fixH1], xmm3
    movaps [rsp + nb332_fiyH1], xmm4
    movaps [rsp + nb332_fizH1], xmm5
    movaps [rsp + nb332_fixH2], xmm6
    movaps [rsp + nb332_fiyH2], xmm7
    movaps [rsp + nb332_fizH2], xmm8
    
    ;# xmm11 = fOx
    ;# xmm12 = fOy
    ;# xmm13 = fOz
    movaps xmm0, xmm14
    unpcklps xmm14, xmm15
    unpckhps xmm0,  xmm15
    
    addps  xmm11, xmm14
    addps  xmm12, xmm0
    
    movhlps  xmm14, xmm13 ;# fOzc fOzd
    
    movlps [rdi + rax*4], xmm11
    movhps [rdi + rbx*4], xmm11
    movlps [rdi + rcx*4], xmm12
    movhps [rdi + rdx*4], xmm12
    movss  [rdi + rax*4 + 8], xmm13
    movss  [rdi + rcx*4 + 8], xmm14
    shufps xmm13, xmm13, 1
    shufps xmm14, xmm14, 1
    movss  [rdi + rbx*4 + 8], xmm13
    movss  [rdi + rdx*4 + 8], xmm14
    
	;# move j H1 coordinates to local temp variables 
	mov   rsi, [rbp + nb332_pos]
    movlps xmm0, [rsi + rax*4 + 12] ;# jxH1a jyH1a  -   -
    movlps xmm1, [rsi + rcx*4 + 12] ;# jxH1c jyH1c  -   -
    movhps xmm0, [rsi + rbx*4 + 12] ;# jxH1a jyH1a jxH1b jyH1b 
    movhps xmm1, [rsi + rdx*4 + 12] ;# jxH1c jyH1c jxH1d jyH1d 

    movss  xmm2, [rsi + rax*4 + 20] ;# jzH1a  -  -  -
    movss  xmm3, [rsi + rcx*4 + 20] ;# jzH1c  -  -  -
    movhps xmm2, [rsi + rbx*4 + 20] ;# jzH1a  -  jzH1b  -
    movhps xmm3, [rsi + rdx*4 + 20] ;# jzH1c  -  jzH1d -

    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxH1a jxH1c jyH1a jyH1c        
    unpckhps xmm4, xmm1  ;# jxH1b jxH1d jyH1b jyH1d
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps   xmm2, xmm3,  136  ;# 10001000 => jzH1a jzH1b jzH1c jzH1d

    ;# xmm0 = H1x
    ;# xmm1 = H1y
    ;# xmm2 = H1z
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb332_ixO]
    subps xmm1, [rsp + nb332_iyO]
    subps xmm2, [rsp + nb332_izO]
    subps xmm3, [rsp + nb332_ixH1]
    subps xmm4, [rsp + nb332_iyH1]
    subps xmm5, [rsp + nb332_izH1]
    subps xmm6, [rsp + nb332_ixH2]
    subps xmm7, [rsp + nb332_iyH2]
    subps xmm8, [rsp + nb332_izH2]
    
	movaps [rsp + nb332_dxOH1], xmm0
	movaps [rsp + nb332_dyOH1], xmm1
	movaps [rsp + nb332_dzOH1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb332_dxH1H1], xmm3
	movaps [rsp + nb332_dyH1H1], xmm4
	movaps [rsp + nb332_dzH1H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb332_dxH2H1], xmm6
	movaps [rsp + nb332_dyH2H1], xmm7
	movaps [rsp + nb332_dzH2H1], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for jH1 atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb332_three]
	movaps  xmm10, xmm9
    movaps  xmm11, xmm9

	mulps   xmm1, xmm0 ;# rsq*lu*lu
	mulps   xmm4, xmm3 ;# rsq*lu*lu 
    mulps   xmm7, xmm6 ;# rsq*lu*lu
	
	subps   xmm9, xmm1
	subps   xmm10, xmm4
    subps   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulps   xmm9, xmm2
	mulps   xmm10, xmm5
    mulps   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movaps  xmm4, [rsp + nb332_half]
	mulps   xmm9, xmm4  ;# rinvOH1
	mulps   xmm10, xmm4 ;# rinvH1H1
    mulps   xmm11, xmm4 ;# rinvH2H1

	movaps  [rsp + nb332_rinvOH1], xmm9
	movaps  [rsp + nb332_rinvH1H1], xmm10
	movaps  [rsp + nb332_rinvH2H1], xmm11
	
	;# H1 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, [rsp + nb332_tsc]
    mulps  xmm0, xmm9  ;# r
    mulps  xmm3, xmm10
    mulps  xmm6, xmm11
    mulps  xmm0, xmm1 ;# rtab
    mulps  xmm3, xmm1
    mulps  xmm6, xmm1

    mov  rsi, [rbp + nb332_VFtab]
    
    ;# truncate and convert to integers
    cvttps2dq xmm1, xmm0
    cvttps2dq xmm4, xmm3
    cvttps2dq xmm7, xmm6        
    
    ;# convert back to float
    cvtdq2ps  xmm2, xmm1
    cvtdq2ps  xmm5, xmm4
    cvtdq2ps  xmm8, xmm7
         
    ;# multiply by 4
    pslld   xmm1, 2
    pslld   xmm4, 2
    pslld   xmm7, 2
    
    ;# multiply by three (copy, mult. by two, add back)
    movaps  xmm10, xmm1
    movaps  xmm11, xmm4
    movaps  xmm12, xmm7
    pslld   xmm1, 1
    pslld   xmm4, 1
    pslld   xmm7, 1    
    paddd   xmm1, xmm10
    paddd   xmm4, xmm11
    paddd   xmm7, xmm12    

    ;# move to integer registers
    movhlps xmm13, xmm1
    movhlps xmm14, xmm4
    movhlps xmm15, xmm7
    movd    eax, xmm1
    movd    r8d, xmm4
    movd    r12d, xmm7
    movd    ecx, xmm13
    movd    r10d, xmm14
    movd    r14d, xmm15
    pshufd  xmm1, xmm1, 1
    pshufd  xmm4, xmm4, 1
    pshufd  xmm7, xmm7, 1
    pshufd  xmm13, xmm13, 1
    pshufd  xmm14, xmm14, 1
    pshufd  xmm15, xmm15, 1
    movd    ebx, xmm1
    movd    r9d, xmm4
    movd    r13d, xmm7    
    movd    edx, xmm13
    movd    r11d, xmm14
    movd    r15d, xmm15   

    ;# calculate eps
    subps     xmm0, xmm2
    subps     xmm3, xmm5
    subps     xmm6, xmm8

    movaps    [rsp + nb332_epsO], xmm0
    movaps    [rsp + nb332_epsH1], xmm3
    movaps    [rsp + nb332_epsH2], xmm6


    ;# Load LOTS of table data
   	movlps xmm1, [rsi + rax*4]
   	movlps xmm5, [rsi + r8*4]
   	movlps xmm9, [rsi + r12*4]

	movlps xmm3, [rsi + rcx*4]
	movlps xmm7, [rsi + r10*4]
	movlps xmm11, [rsi + r14*4]

	movhps xmm1, [rsi + rbx*4]
	movhps xmm5, [rsi + r9*4]
	movhps xmm9, [rsi + r13*4]

	movhps xmm3, [rsi + rdx*4]
	movhps xmm7, [rsi + r11*4]
	movhps xmm11, [rsi + r15*4]

    movaps xmm0, xmm1
    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm0, xmm3, 136  ;# 10001000
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm8, xmm11, 136  ;# 10001000
	shufps xmm1, xmm3, 221  ;# 11011101
	shufps xmm5, xmm7, 221  ;# 11011101
	shufps xmm9, xmm11, 221  ;# 11011101
    
	movlps xmm3, [rsi + rax*4 + 8]
	movlps xmm7, [rsi + r8*4 + 8]
	movlps xmm11, [rsi + r12*4 + 8]
    
	movlps xmm12, [rsi + rcx*4 + 8]
	movlps xmm13, [rsi + r10*4 + 8]
	movlps xmm14, [rsi + r14*4 + 8]

	movhps xmm3, [rsi + rbx*4 + 8]
	movhps xmm7, [rsi + r9*4 + 8]
	movhps xmm11, [rsi + r13*4 + 8]
    
	movhps xmm12, [rsi + rdx*4 + 8]
	movhps xmm13, [rsi + r11*4 + 8]
	movhps xmm14, [rsi + r15*4 + 8]

    movaps xmm2, xmm3
    movaps xmm6, xmm7
    movaps xmm10, xmm11
    
	shufps xmm2, xmm12, 136  ;# 10001000
	shufps xmm6, xmm13, 136  ;# 10001000
	shufps xmm10, xmm14, 136  ;# 10001000
	shufps xmm3, xmm12, 221  ;# 11011101
	shufps xmm7, xmm13, 221  ;# 11011101
	shufps xmm11, xmm14, 221  ;# 11011101
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movaps xmm12, [rsp + nb332_epsO]
    movaps xmm13, [rsp + nb332_epsH1]
    movaps xmm14, [rsp + nb332_epsH2]
    
    mulps  xmm3, xmm12   ;# Heps
    mulps  xmm7, xmm13
    mulps  xmm11, xmm14 
    mulps  xmm2, xmm12   ;# Geps
    mulps  xmm6, xmm13
    mulps  xmm10, xmm14 
    mulps  xmm3, xmm12   ;# Heps2
    mulps  xmm7, xmm13
    mulps  xmm11, xmm14 

    addps  xmm1, xmm2   ;# F+Geps
    addps  xmm5, xmm6
    addps  xmm9, xmm10 
    addps  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addps  xmm5, xmm7
    addps  xmm9, xmm11 
    addps  xmm3, xmm3    ;# 2*Heps2
    addps  xmm7, xmm7
    addps  xmm11, xmm11
    addps  xmm3, xmm2    ;# 2*Heps2+Geps
    addps  xmm7, xmm6  
    addps  xmm11, xmm10
    addps  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addps  xmm7, xmm5
    addps  xmm11, xmm9
    mulps  xmm1, xmm12   ;# eps*Fp
    mulps  xmm5, xmm13
    mulps  xmm9, xmm14
    movaps xmm12, [rsp + nb332_qqOH]
    movaps xmm13, [rsp + nb332_qqHH]
    addps  xmm1, xmm0     ;# VV
    addps  xmm5, xmm4
    addps  xmm9, xmm8
    mulps  xmm1, xmm12   ;# VV*qq = vcoul
    mulps  xmm5, xmm13
    mulps  xmm9, xmm13
    mulps  xmm3, xmm12    ;# FF*qq = fij
    mulps  xmm7, xmm13
    mulps  xmm11, xmm13
    
    ;# accumulate vctot
    addps  xmm1, [rsp + nb332_vctot]
    addps  xmm5, xmm9
    addps  xmm1, xmm5
    movaps [rsp + nb332_vctot], xmm1
    
    movaps xmm10, [rsp + nb332_tsc]
    mulps  xmm3, xmm10  ;# fscal
    mulps  xmm7, xmm10
    mulps  xmm10, xmm11
        
    movd eax, mm0 ;# restore j3 from mm0-mm3
    movd ebx, mm1
    movd ecx, mm2
    movd edx, mm3
    
	;# move j H1 forces to local temp variables 
    mov rdi, [rbp + nb332_faction]
    movlps xmm11, [rdi + rax*4 + 12] ;# jxH1a jyH1a  -   -
    movlps xmm12, [rdi + rcx*4 + 12] ;# jxH1c jyH1c  -   -

    ;# xmm11: jxH1a jyH1a jxH1b jyH1b 
    ;# xmm12: jxH1c jyH1c jxH1d jyH1d
    ;# xmm13: jzH1a jzH1b jzH1c jzH1d

    xorps  xmm0, xmm0
    xorps  xmm4, xmm4
    xorps  xmm8, xmm8

    movhps xmm11, [rdi + rbx*4 + 12] ;# jxH1a jyH1a jxH1b jyH1b 
    movhps xmm12, [rdi + rdx*4 + 12] ;# jxH1c jyH1c jxH1d jyH1d 

    mulps  xmm3, [rsp + nb332_rinvOH1]
    mulps  xmm7, [rsp + nb332_rinvH1H1]
    mulps  xmm10, [rsp + nb332_rinvH2H1]
    
    movss  xmm13, [rdi + rax*4 + 20] ;# jzH1a  -  -  -
    movss  xmm14, [rdi + rcx*4 + 20] ;# jzH1c  -  -  -
    subps  xmm0, xmm3
    subps  xmm4, xmm7
    subps  xmm8, xmm10
    
    movhps xmm13, [rdi + rbx*4 + 20] ;# jzH1a  -  jzH1b  -
    movhps xmm14, [rdi + rdx*4 + 20] ;# jzH1c  -  jzH1d -
    
    movaps xmm1, xmm0
    movaps xmm2, xmm0
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm8
    movaps xmm7, xmm8

    shufps xmm13, xmm14,  136  ;# 10001000 => jzH1a jzH1b jzH1c jzH1d

	mulps xmm0, [rsp + nb332_dxOH1]
	mulps xmm1, [rsp + nb332_dyOH1]
	mulps xmm2, [rsp + nb332_dzOH1]
	mulps xmm3, [rsp + nb332_dxH1H1]
	mulps xmm4, [rsp + nb332_dyH1H1]
	mulps xmm5, [rsp + nb332_dzH1H1]
	mulps xmm6, [rsp + nb332_dxH2H1]
	mulps xmm7, [rsp + nb332_dyH2H1]
	mulps xmm8, [rsp + nb332_dzH2H1]

    movaps xmm14,  xmm0
    movaps xmm15, xmm1
    addps xmm13, xmm2
    addps xmm0, [rsp + nb332_fixO]
    addps xmm1, [rsp + nb332_fiyO]
    addps xmm2, [rsp + nb332_fizO]

    addps xmm14,  xmm3
    addps xmm15, xmm4
    addps xmm13, xmm5
    addps xmm3, [rsp + nb332_fixH1]
    addps xmm4, [rsp + nb332_fiyH1]
    addps xmm5, [rsp + nb332_fizH1]

    addps xmm14,  xmm6
    addps xmm15, xmm7
    addps xmm13, xmm8
    addps xmm6, [rsp + nb332_fixH2]
    addps xmm7, [rsp + nb332_fiyH2]
    addps xmm8, [rsp + nb332_fizH2]

    movaps [rsp + nb332_fixO], xmm0
    movaps [rsp + nb332_fiyO], xmm1
    movaps [rsp + nb332_fizO], xmm2
    movaps [rsp + nb332_fixH1], xmm3
    movaps [rsp + nb332_fiyH1], xmm4
    movaps [rsp + nb332_fizH1], xmm5
    movaps [rsp + nb332_fixH2], xmm6
    movaps [rsp + nb332_fiyH2], xmm7
    movaps [rsp + nb332_fizH2], xmm8
    
    ;# xmm11 = fH1x
    ;# xmm12 = fH1y
    ;# xmm13 = fH1z
    movaps xmm0, xmm14
    unpcklps xmm14, xmm15
    unpckhps xmm0,  xmm15
    
    addps  xmm11, xmm14
    addps  xmm12, xmm0
    
    movhlps  xmm14, xmm13 ;# fH1zc fH1zd
    
    movlps [rdi + rax*4 + 12], xmm11
    movhps [rdi + rbx*4 + 12], xmm11
    movlps [rdi + rcx*4 + 12], xmm12
    movhps [rdi + rdx*4 + 12], xmm12
    movss  [rdi + rax*4 + 20], xmm13
    movss  [rdi + rcx*4 + 20], xmm14
    shufps xmm13, xmm13, 1
    shufps xmm14, xmm14, 1
    movss  [rdi + rbx*4 + 20], xmm13
    movss  [rdi + rdx*4 + 20], xmm14
    
   	mov   rsi, [rbp + nb332_pos]
	;# move j H2 coordinates to local temp variables 
    movlps xmm0, [rsi + rax*4 + 24] ;# jxH2a jyH2a  -   -
    movlps xmm1, [rsi + rcx*4 + 24] ;# jxH2c jyH2c  -   -
    movhps xmm0, [rsi + rbx*4 + 24] ;# jxH2a jyH2a jxH2b jyH2b 
    movhps xmm1, [rsi + rdx*4 + 24] ;# jxH2c jyH2c jxH2d jyH2d 

    movss  xmm2, [rsi + rax*4 + 32] ;# jzH2a  -  -  -
    movss  xmm3, [rsi + rcx*4 + 32] ;# jzH2c  -  -  -
    movss  xmm5, [rsi + rbx*4 + 32] ;# jzH2b  -  -  -
    movss  xmm6, [rsi + rdx*4 + 32] ;# jzH2d  -  -  -
    movlhps xmm2, xmm5 ;# jzH2a  -  jzH2b  -
    movlhps xmm3, xmm6 ;# jzH2c  -  jzH2d -

    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxH2a jxH2c jyH2a jyH2c        
    unpckhps xmm4, xmm1  ;# jxH2b jxH2d jyH2b jyH2d
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps   xmm2, xmm3,  136  ;# 10001000 => jzH2a jzH2b jzH2c jzH2d

    ;# xmm0 = H2x
    ;# xmm1 = H2y
    ;# xmm2 = H2z
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb332_ixO]
    subps xmm1, [rsp + nb332_iyO]
    subps xmm2, [rsp + nb332_izO]
    subps xmm3, [rsp + nb332_ixH1]
    subps xmm4, [rsp + nb332_iyH1]
    subps xmm5, [rsp + nb332_izH1]
    subps xmm6, [rsp + nb332_ixH2]
    subps xmm7, [rsp + nb332_iyH2]
    subps xmm8, [rsp + nb332_izH2]
    
	movaps [rsp + nb332_dxOH2], xmm0
	movaps [rsp + nb332_dyOH2], xmm1
	movaps [rsp + nb332_dzOH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb332_dxH1H2], xmm3
	movaps [rsp + nb332_dyH1H2], xmm4
	movaps [rsp + nb332_dzH1H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb332_dxH2H2], xmm6
	movaps [rsp + nb332_dyH2H2], xmm7
	movaps [rsp + nb332_dzH2H2], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for jH2 atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb332_three]
	movaps  xmm10, xmm9
    movaps  xmm11, xmm9

	mulps   xmm1, xmm0 ;# rsq*lu*lu
	mulps   xmm4, xmm3 ;# rsq*lu*lu 
    mulps   xmm7, xmm6 ;# rsq*lu*lu
	
	subps   xmm9, xmm1
	subps   xmm10, xmm4
    subps   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulps   xmm9, xmm2
	mulps   xmm10, xmm5
    mulps   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movaps  xmm4, [rsp + nb332_half]
	mulps   xmm9, xmm4  ;# rinvOH2
	mulps   xmm10, xmm4 ;# rinvH1H2
    mulps   xmm11, xmm4 ;# rinvH2H2

	movaps  [rsp + nb332_rinvOH2], xmm9
	movaps  [rsp + nb332_rinvH1H2], xmm10
	movaps  [rsp + nb332_rinvH2H2], xmm11
	
	;# H2 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, [rsp + nb332_tsc]
    mulps  xmm0, xmm9  ;# r
    mulps  xmm3, xmm10
    mulps  xmm6, xmm11
    mulps  xmm0, xmm1 ;# rtab
    mulps  xmm3, xmm1
    mulps  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttps2dq xmm1, xmm0
    cvttps2dq xmm4, xmm3
    cvttps2dq xmm7, xmm6        
    
    ;# convert back to float
    cvtdq2ps  xmm2, xmm1
    cvtdq2ps  xmm5, xmm4
    cvtdq2ps  xmm8, xmm7
     
    ;# multiply by 4
    pslld   xmm1, 2
    pslld   xmm4, 2
    pslld   xmm7, 2
    
    ;# multiply by three (copy, mult. by two, add back)
    movaps  xmm10, xmm1
    movaps  xmm11, xmm4
    movaps  xmm12, xmm7
    pslld   xmm1, 1
    pslld   xmm4, 1
    pslld   xmm7, 1    
    paddd   xmm1, xmm10
    paddd   xmm4, xmm11
    paddd   xmm7, xmm12    

    ;# move to integer registers
    movhlps xmm13, xmm1
    movhlps xmm14, xmm4
    movhlps xmm15, xmm7
    movd    eax, xmm1
    movd    r8d, xmm4
    movd    r12d, xmm7
    movd    ecx, xmm13
    movd    r10d, xmm14
    movd    r14d, xmm15
    pshufd  xmm1, xmm1, 1
    pshufd  xmm4, xmm4, 1
    pshufd  xmm7, xmm7, 1
    pshufd  xmm13, xmm13, 1
    pshufd  xmm14, xmm14, 1
    pshufd  xmm15, xmm15, 1
    movd    ebx, xmm1
    movd    r9d, xmm4
    movd    r13d, xmm7    
    movd    edx, xmm13
    movd    r11d, xmm14
    movd    r15d, xmm15   
    
    mov  rsi, [rbp + nb332_VFtab]

    ;# calculate eps
    subps     xmm0, xmm2
    subps     xmm3, xmm5
    subps     xmm6, xmm8

    movaps    [rsp + nb332_epsO], xmm0
    movaps    [rsp + nb332_epsH1], xmm3
    movaps    [rsp + nb332_epsH2], xmm6

    ;# Load LOTS of table data
   	movlps xmm1, [rsi + rax*4]
   	movlps xmm5, [rsi + r8*4]
   	movlps xmm9, [rsi + r12*4]

	movlps xmm3, [rsi + rcx*4]
	movlps xmm7, [rsi + r10*4]
	movlps xmm11, [rsi + r14*4]

	movhps xmm1, [rsi + rbx*4]
	movhps xmm5, [rsi + r9*4]
	movhps xmm9, [rsi + r13*4]

	movhps xmm3, [rsi + rdx*4]
	movhps xmm7, [rsi + r11*4]
	movhps xmm11, [rsi + r15*4]

    movaps xmm0, xmm1
    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm0, xmm3, 136  ;# 10001000
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm8, xmm11, 136  ;# 10001000
	shufps xmm1, xmm3, 221  ;# 11011101
	shufps xmm5, xmm7, 221  ;# 11011101
	shufps xmm9, xmm11, 221  ;# 11011101
    
	movlps xmm3, [rsi + rax*4 + 8]
	movlps xmm7, [rsi + r8*4 + 8]
	movlps xmm11, [rsi + r12*4 + 8]
    
	movlps xmm12, [rsi + rcx*4 + 8]
	movlps xmm13, [rsi + r10*4 + 8]
	movlps xmm14, [rsi + r14*4 + 8]

	movhps xmm3, [rsi + rbx*4 + 8]
	movhps xmm7, [rsi + r9*4 + 8]
	movhps xmm11, [rsi + r13*4 + 8]
    
	movhps xmm12, [rsi + rdx*4 + 8]
	movhps xmm13, [rsi + r11*4 + 8]
	movhps xmm14, [rsi + r15*4 + 8]

    movaps xmm2, xmm3
    movaps xmm6, xmm7
    movaps xmm10, xmm11
    
	shufps xmm2, xmm12, 136  ;# 10001000
	shufps xmm6, xmm13, 136  ;# 10001000
	shufps xmm10, xmm14, 136  ;# 10001000
	shufps xmm3, xmm12, 221  ;# 11011101
	shufps xmm7, xmm13, 221  ;# 11011101
	shufps xmm11, xmm14, 221  ;# 11011101
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movaps xmm12, [rsp + nb332_epsO]
    movaps xmm13, [rsp + nb332_epsH1]
    movaps xmm14, [rsp + nb332_epsH2]
    
    mulps  xmm3, xmm12   ;# Heps
    mulps  xmm7, xmm13
    mulps  xmm11, xmm14 
    mulps  xmm2, xmm12   ;# Geps
    mulps  xmm6, xmm13
    mulps  xmm10, xmm14 
    mulps  xmm3, xmm12   ;# Heps2
    mulps  xmm7, xmm13
    mulps  xmm11, xmm14 

    addps  xmm1, xmm2   ;# F+Geps
    addps  xmm5, xmm6
    addps  xmm9, xmm10 
    addps  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addps  xmm5, xmm7
    addps  xmm9, xmm11 
    addps  xmm3, xmm3    ;# 2*Heps2
    addps  xmm7, xmm7
    addps  xmm11, xmm11
    addps  xmm3, xmm2    ;# 2*Heps2+Geps
    addps  xmm7, xmm6  
    addps  xmm11, xmm10
    addps  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addps  xmm7, xmm5
    addps  xmm11, xmm9
    mulps  xmm1, xmm12   ;# eps*Fp
    mulps  xmm5, xmm13
    mulps  xmm9, xmm14
    movaps xmm12, [rsp + nb332_qqOH]
    movaps xmm13, [rsp + nb332_qqHH]
    addps  xmm1, xmm0     ;# VV
    addps  xmm5, xmm4
    addps  xmm9, xmm8
    mulps  xmm1, xmm12   ;# VV*qq = vcoul
    mulps  xmm5, xmm13
    mulps  xmm9, xmm13
    mulps  xmm3, xmm12    ;# FF*qq = fij
    mulps  xmm7, xmm13
    mulps  xmm11, xmm13
    
    ;# accumulate vctot
    addps  xmm1, [rsp + nb332_vctot]
    addps  xmm5, xmm9
    addps  xmm1, xmm5
    movaps [rsp + nb332_vctot], xmm1
    
    movaps xmm10, [rsp + nb332_tsc]
    mulps  xmm3, xmm10  ;# fscal
    mulps  xmm7, xmm10
    mulps  xmm10, xmm11
        
    movd eax, mm0 ;# restore j3 from mm0-mm3
    movd ebx, mm1
    movd ecx, mm2
    movd edx, mm3

    
	;# move j H2 forces to local temp variables 
    mov rdi, [rbp + nb332_faction]
    movlps xmm11, [rdi + rax*4 + 24] ;# jxH2a jyH2a  -   -
    movlps xmm12, [rdi + rcx*4 + 24] ;# jxH2c jyH2c  -   -
    movhps xmm11, [rdi + rbx*4 + 24] ;# jxH2a jyH2a jxH2b jyH2b 
    movhps xmm12, [rdi + rdx*4 + 24] ;# jxH2c jyH2c jxH2d jyH2d 

    movss  xmm13, [rdi + rax*4 + 32] ;# jzH2a  -  -  -
    movss  xmm14, [rdi + rcx*4 + 32] ;# jzH2c  -  -  -
    movss  xmm1, [rdi + rbx*4 + 32] ;# jzH2b  -  -  -
    movss  xmm2, [rdi + rdx*4 + 32] ;# jzH2d  -  -  -
    movlhps xmm13, xmm1 ;# jzH2a  -  jzH2b  -
    movlhps xmm14, xmm2 ;# jzH2c  -  jzH2d -

    shufps xmm13, xmm14,  136  ;# 10001000 => jzH2a jzH2b jzH2c jzH2d

    ;# xmm11: jxH2a jyH2a jxH2b jyH2b 
    ;# xmm12: jxH2c jyH2c jxH2d jyH2d
    ;# xmm13: jzH2a jzH2b jzH2c jzH2d

    xorps  xmm0, xmm0
    xorps  xmm4, xmm4
    xorps  xmm8, xmm8

    mulps  xmm3, [rsp + nb332_rinvOH2]
    mulps  xmm7, [rsp + nb332_rinvH1H2]
    mulps  xmm10, [rsp + nb332_rinvH2H2]
    
    subps  xmm0, xmm3
    subps  xmm4, xmm7
    subps  xmm8, xmm10
    
    movaps xmm1, xmm0
    movaps xmm2, xmm0
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm8
    movaps xmm7, xmm8

	mulps xmm0, [rsp + nb332_dxOH2]
	mulps xmm1, [rsp + nb332_dyOH2]
	mulps xmm2, [rsp + nb332_dzOH2]
	mulps xmm3, [rsp + nb332_dxH1H2]
	mulps xmm4, [rsp + nb332_dyH1H2]
	mulps xmm5, [rsp + nb332_dzH1H2]
	mulps xmm6, [rsp + nb332_dxH2H2]
	mulps xmm7, [rsp + nb332_dyH2H2]
	mulps xmm8, [rsp + nb332_dzH2H2]

    movaps xmm14,  xmm0
    movaps xmm15, xmm1
    addps xmm13, xmm2
    addps xmm0, [rsp + nb332_fixO]
    addps xmm1, [rsp + nb332_fiyO]
    addps xmm2, [rsp + nb332_fizO]

    addps xmm14,  xmm3
    addps xmm15, xmm4
    addps xmm13, xmm5
    addps xmm3, [rsp + nb332_fixH1]
    addps xmm4, [rsp + nb332_fiyH1]
    addps xmm5, [rsp + nb332_fizH1]

    addps xmm14,  xmm6
    addps xmm15, xmm7
    addps xmm13, xmm8
    addps xmm6, [rsp + nb332_fixH2]
    addps xmm7, [rsp + nb332_fiyH2]
    addps xmm8, [rsp + nb332_fizH2]

    movaps [rsp + nb332_fixO], xmm0
    movaps [rsp + nb332_fiyO], xmm1
    movaps [rsp + nb332_fizO], xmm2
    movaps [rsp + nb332_fixH1], xmm3
    movaps [rsp + nb332_fiyH1], xmm4
    movaps [rsp + nb332_fizH1], xmm5
    movaps [rsp + nb332_fixH2], xmm6
    movaps [rsp + nb332_fiyH2], xmm7
    movaps [rsp + nb332_fizH2], xmm8
    
    ;# xmm11 = fH2x
    ;# xmm12 = fH2y
    ;# xmm13 = fH2z
    movaps xmm0, xmm14
    unpcklps xmm14, xmm15
    unpckhps xmm0,  xmm15
    
    addps  xmm11, xmm14
    addps  xmm12, xmm0
    
    movhlps  xmm14, xmm13 ;# fH2zc fH2zd
    
    movlps [rdi + rax*4 + 24], xmm11
    movhps [rdi + rbx*4 + 24], xmm11
    movlps [rdi + rcx*4 + 24], xmm12
    movhps [rdi + rdx*4 + 24], xmm12
    movss  [rdi + rax*4 + 32], xmm13
    movss  [rdi + rcx*4 + 32], xmm14
    shufps xmm13, xmm13, 1
    shufps xmm14, xmm14, 1
    movss  [rdi + rbx*4 + 32], xmm13
    movss  [rdi + rdx*4 + 32], xmm14
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb332_innerk],  4
	jl    .nb332_single_check
	jmp   .nb332_unroll_loop
.nb332_single_check:
	add dword ptr [rsp + nb332_innerk],  4
	jnz   .nb332_single_loop
	jmp   .nb332_updateouterdata
.nb332_single_loop:
	mov   rdx, [rsp + nb332_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb332_innerjjnr],  4	

	mov rsi, [rbp + nb332_pos]
	lea   rax, [rax + rax*2]  

	;# fetch j coordinates 
	xorps xmm0, xmm0
	xorps xmm1, xmm1
	xorps xmm2, xmm2
	
	movss xmm0, [rsi + rax*4]		;# jxO  -  -  -
	movss xmm1, [rsi + rax*4 + 4]		;# jyO  -  -  -
	movss xmm2, [rsi + rax*4 + 8]		;# jzO  -  -  -  

	movlps xmm6, [rsi + rax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, [rsi + rax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, [rsi + rax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm5, [rsi + rax*4 + 32]		;# xmm5 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm5			;# xmm7 = jzH1 jzH2   -    -
	movlhps xmm0, xmm6			;# xmm0 = jxO   0   jxH1 jxH2 
	shufps  xmm1, xmm6, 228 ;# 11100100	;# xmm1 = jyO   0   jyH1 jyH2 
	shufps  xmm2, xmm7, 68  ;# 01000100	;# xmm2 = jzO   0   jzH1 jzH2

	;# store all j coordinates in jO  
	movaps [rsp + nb332_jxO], xmm0
	movaps [rsp + nb332_jyO], xmm1
	movaps [rsp + nb332_jzO], xmm2
	subps  xmm0, [rsp + nb332_ixO]
	subps  xmm1, [rsp + nb332_iyO]
	subps  xmm2, [rsp + nb332_izO]
	movaps [rsp + nb332_dxOO], xmm0
	movaps [rsp + nb332_dyOO], xmm1
	movaps [rsp + nb332_dzOO], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1	
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb332_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [rsp + nb332_half] ;# rinv iO - j water 

	movaps  xmm1, xmm3
	mulps   xmm1, xmm0	;# xmm1=r 
	movaps  xmm0, xmm3	;# xmm0=rinv 
	mulps  xmm1, [rsp + nb332_tsc]
	
	movhlps xmm2, xmm1	
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 2
    pslld mm7, 2
    movd ebx, mm6
    movd ecx, mm7
    psrlq mm7, 32
    movd edx, mm7		;# table indices in ebx,ecx,edx 

	mov rsi, [rbp + nb332_VFtab]
	
    lea   rbx, [rbx + rbx*2]
    lea   rcx, [rcx + rcx*2]
    lea   rdx, [rdx + rdx*2]
	
    movlps xmm5, [rsi + rbx*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rbx*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    mulps  xmm7, [rsp + nb332_two]       ;# two*Heps2 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb332_qqOO]
	movhps  xmm3, [rsp + nb332_qqOH]
	
    addps  xmm7, xmm6
    addps  xmm7, xmm5 ;# xmm7=FF 
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    mulps  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point xmm5 contains vcoul and xmm3 fijC 
	
    addps  xmm5, [rsp + nb332_vctot]
    movaps [rsp + nb332_vctot], xmm5
    ;# put scalar force on stack temporarily 
    movaps [rsp + nb332_fstmp], xmm3

    ;# dispersion 
	movss  xmm4, [rsi + rbx*4 + 16]	
	movss  xmm5, [rsi + rbx*4 + 20]	
	movss  xmm6, [rsi + rbx*4 + 24]	
	movss  xmm7, [rsi + rbx*4 + 28]
    ;# dispersion table ready, in xmm4-xmm7 
    mulss  xmm6, xmm1       ;# xmm6=Geps 
    mulss  xmm7, xmm2       ;# xmm7=Heps2 
    addss  xmm5, xmm6
    addss  xmm5, xmm7       ;# xmm5=Fp 
    mulss  xmm7, [rsp + nb332_two]       ;# two*Heps2 
    addss  xmm7, xmm6
    addss  xmm7, xmm5 ;# xmm7=FF 
    mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
    addss  xmm5, xmm4 ;# xmm5=VV 
	xorps  xmm4, xmm4
    movss  xmm4, [rsp + nb332_c6]
    mulps  xmm7, xmm4    ;# fijD 
    mulps  xmm5, xmm4    ;# Vvdw6 
    addps  xmm7, [rsp + nb332_fstmp] ;# add to fscal 

    ;# put scalar force on stack Update Vvdwtot directly 
    addps  xmm5, [rsp + nb332_Vvdwtot]
    movaps [rsp + nb332_fstmp], xmm7
    movaps [rsp + nb332_Vvdwtot], xmm5

    ;# repulsion 
	movss  xmm4, [rsi + rbx*4 + 32]	
	movss  xmm5, [rsi + rbx*4 + 36]	
	movss  xmm6, [rsi + rbx*4 + 40]	
	movss  xmm7, [rsi + rbx*4 + 44]
    ;# table ready, in xmm4-xmm7 
    mulss  xmm6, xmm1       ;# xmm6=Geps 
    mulss  xmm7, xmm2       ;# xmm7=Heps2 
    addss  xmm5, xmm6
    addss  xmm5, xmm7       ;# xmm5=Fp 
    mulss  xmm7, [rsp + nb332_two]       ;# two*Heps2 
    addss  xmm7, xmm6
    addss  xmm7, xmm5 ;# xmm7=FF 
    mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
    addss  xmm5, xmm4 ;# xmm5=VV 

	xorps  xmm4, xmm4
    movss  xmm4, [rsp + nb332_c12]
    mulps  xmm7, xmm4 ;# fijR 
    mulps  xmm5, xmm4 ;# Vvdw12 
    addps  xmm7, [rsp + nb332_fstmp] 

    addps  xmm5, [rsp + nb332_Vvdwtot]
    movaps [rsp + nb332_Vvdwtot], xmm5
    xorps  xmm1, xmm1

    mulps xmm7, [rsp + nb332_tsc]
    mulps xmm7, xmm0
    subps  xmm1, xmm7

	movaps xmm0, xmm1
	movaps xmm2, xmm1		

	mulps   xmm0, [rsp + nb332_dxOO]
	mulps   xmm1, [rsp + nb332_dyOO]
	mulps   xmm2, [rsp + nb332_dzOO]
	;# initial update for j forces 
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb332_fjxO], xmm3
	movaps  [rsp + nb332_fjyO], xmm4
	movaps  [rsp + nb332_fjzO], xmm5
	addps   xmm0, [rsp + nb332_fixO]
	addps   xmm1, [rsp + nb332_fiyO]
	addps   xmm2, [rsp + nb332_fizO]
	movaps  [rsp + nb332_fixO], xmm0
	movaps  [rsp + nb332_fiyO], xmm1
	movaps  [rsp + nb332_fizO], xmm2

	
	;# done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
    movaps  xmm0, [rsp + nb332_jxO]
    movaps  xmm1, [rsp + nb332_jyO]
    movaps  xmm2, [rsp + nb332_jzO]
    movaps  xmm3, xmm0
    movaps  xmm4, xmm1
    movaps  xmm5, xmm2
	subps  xmm0, [rsp + nb332_ixH1]
	subps  xmm1, [rsp + nb332_iyH1]
	subps  xmm2, [rsp + nb332_izH1]
	subps  xmm3, [rsp + nb332_ixH2]
	subps  xmm4, [rsp + nb332_iyH2]
	subps  xmm5, [rsp + nb332_izH2]
    movaps [rsp + nb332_dxH1O], xmm0
	movaps [rsp + nb332_dyH1O], xmm1
	movaps [rsp + nb332_dzH1O], xmm2
	movaps [rsp + nb332_dxH2O], xmm3
	movaps [rsp + nb332_dyH2O], xmm4
	movaps [rsp + nb332_dzH2O], xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;# have rsqH1 in xmm0 
	addps xmm4, xmm5	;# have rsqH2 in xmm4 

	;# start with H1, save H2 data 
	movaps [rsp + nb332_rsqH2O], xmm4
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb332_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb332_half] ;# rinv H1 - j water 
	mulps   xmm7, [rsp + nb332_half] ;# rinv H2 - j water  

	;# start with H1, save H2 data 
	movaps [rsp + nb332_rinvH2O], xmm7

	movaps xmm1, xmm3
	mulps  xmm1, xmm0	;# xmm1=r 
	movaps xmm0, xmm3	;# xmm0=rinv 
	mulps  xmm1, [rsp + nb332_tsc]
	
	movhlps xmm2, xmm1	
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 2
    pslld mm7, 2
    movd ebx, mm6
    movd ecx, mm7
    psrlq mm7, 32
    movd edx, mm7		;# table indices in ebx,ecx,edx 

    lea   rbx, [rbx + rbx*2]
    lea   rcx, [rcx + rcx*2]
    lea   rdx, [rdx + rdx*2]
	
    movlps xmm5, [rsi + rbx*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rbx*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    mulps  xmm7, [rsp + nb332_two]       ;# two*Heps2 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb332_qqOH]
	movhps  xmm3, [rsp + nb332_qqHH]
	
    addps  xmm7, xmm6
    addps  xmm7, xmm5 ;# xmm7=FF 
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    mulps  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point xmm5 contains vcoul and xmm3 fijC 
    addps  xmm5, [rsp + nb332_vctot]
    movaps [rsp + nb332_vctot], xmm5	

    xorps  xmm1, xmm1

    mulps xmm3, [rsp + nb332_tsc]
    mulps xmm3, xmm0
    subps  xmm1, xmm3
	
	movaps  xmm0, xmm1
	movaps  xmm2, xmm1
	mulps   xmm0, [rsp + nb332_dxH1O]
	mulps   xmm1, [rsp + nb332_dyH1O]
	mulps   xmm2, [rsp + nb332_dzH1O]
	;# update forces H1 - j water 
	movaps  xmm3, [rsp + nb332_fjxO]
	movaps  xmm4, [rsp + nb332_fjyO]
	movaps  xmm5, [rsp + nb332_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb332_fjxO], xmm3
	movaps  [rsp + nb332_fjyO], xmm4
	movaps  [rsp + nb332_fjzO], xmm5
	addps   xmm0, [rsp + nb332_fixH1]
	addps   xmm1, [rsp + nb332_fiyH1]
	addps   xmm2, [rsp + nb332_fizH1]
	movaps  [rsp + nb332_fixH1], xmm0
	movaps  [rsp + nb332_fiyH1], xmm1
	movaps  [rsp + nb332_fizH1], xmm2
	;# do table for H2 - j water interaction 
	movaps xmm0, [rsp + nb332_rinvH2O]
	movaps xmm1, [rsp + nb332_rsqH2O]
	mulps  xmm1, xmm0	;# xmm0=rinv, xmm1=r 
	mulps  xmm1, [rsp + nb332_tsc]
	
	movhlps xmm2, xmm1	
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 2
    pslld mm7, 2
    movd ebx, mm6
    movd ecx, mm7
    psrlq mm7, 32
    movd edx, mm7		;# table indices in ebx,ecx,edx 

    lea   rbx, [rbx + rbx*2]
    lea   rcx, [rcx + rcx*2]
    lea   rdx, [rdx + rdx*2]
	
    movlps xmm5, [rsi + rbx*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rbx*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    mulps  xmm7, [rsp + nb332_two]       ;# two*Heps2 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb332_qqOH]
	movhps  xmm3, [rsp + nb332_qqHH]
	
    addps  xmm7, xmm6
    addps  xmm7, xmm5 ;# xmm7=FF 
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    mulps  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point xmm5 contains vcoul and xmm3 fijC 
    addps  xmm5, [rsp + nb332_vctot]
    movaps [rsp + nb332_vctot], xmm5	

    xorps  xmm1, xmm1

    mulps xmm3, [rsp + nb332_tsc]
    mulps xmm3, xmm0
    subps  xmm1, xmm3
	
	movaps  xmm0, xmm1
	movaps  xmm2, xmm1
	
	mulps   xmm0, [rsp + nb332_dxH2O]
	mulps   xmm1, [rsp + nb332_dyH2O]
	mulps   xmm2, [rsp + nb332_dzH2O]
	movaps  xmm3, [rsp + nb332_fjxO]
	movaps  xmm4, [rsp + nb332_fjyO]
	movaps  xmm5, [rsp + nb332_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	mov     rsi, [rbp + nb332_faction]
	movaps  [rsp + nb332_fjxO], xmm3
	movaps  [rsp + nb332_fjyO], xmm4
	movaps  [rsp + nb332_fjzO], xmm5
	addps   xmm0, [rsp + nb332_fixH2]
	addps   xmm1, [rsp + nb332_fiyH2]
	addps   xmm2, [rsp + nb332_fizH2]
	movaps  [rsp + nb332_fixH2], xmm0
	movaps  [rsp + nb332_fiyH2], xmm1
	movaps  [rsp + nb332_fizH2], xmm2

	;# update j water forces from local variables 
	movlps  xmm0, [rsi + rax*4]
	movlps  xmm1, [rsi + rax*4 + 12]
	movhps  xmm1, [rsi + rax*4 + 24]
	movaps  xmm3, [rsp + nb332_fjxO]
	movaps  xmm4, [rsp + nb332_fjyO]
	movaps  xmm5, [rsp + nb332_fjzO]
	movaps  xmm6, xmm5
	movaps  xmm7, xmm5
	shufps  xmm6, xmm6, 2 ;# 00000010
	shufps  xmm7, xmm7, 3 ;# 00000011
	addss   xmm5, [rsi + rax*4 + 8]
	addss   xmm6, [rsi + rax*4 + 20]
	addss   xmm7, [rsi + rax*4 + 32]
	movss   [rsi + rax*4 + 8], xmm5
	movss   [rsi + rax*4 + 20], xmm6
	movss   [rsi + rax*4 + 32], xmm7
	movaps   xmm5, xmm3
	unpcklps xmm3, xmm4
	unpckhps xmm5, xmm4
	addps    xmm0, xmm3
	addps    xmm1, xmm5
	movlps  [rsi + rax*4], xmm0 
	movlps  [rsi + rax*4 + 12], xmm1 
	movhps  [rsi + rax*4 + 24], xmm1 
	
	dec dword ptr [rsp + nb332_innerk]
	jz    .nb332_updateouterdata
	jmp   .nb332_single_loop
.nb332_updateouterdata:
	mov   ecx, [rsp + nb332_ii3]
	mov   rdi, [rbp + nb332_faction]
	mov   rsi, [rbp + nb332_fshift]
	mov   edx, [rsp + nb332_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb332_fixO]
	movaps xmm1, [rsp + nb332_fiyO] 
	movaps xmm2, [rsp + nb332_fizO]

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
	movss  xmm3, [rdi + rcx*4]
	movss  xmm4, [rdi + rcx*4 + 4]
	movss  xmm5, [rdi + rcx*4 + 8]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4],     xmm3
	movss  [rdi + rcx*4 + 4], xmm4
	movss  [rdi + rcx*4 + 8], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# 00001000	

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb332_fixH1]
	movaps xmm1, [rsp + nb332_fiyH1]
	movaps xmm2, [rsp + nb332_fizH1]

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
	movss  xmm3, [rdi + rcx*4 + 12]
	movss  xmm4, [rdi + rcx*4 + 16]
	movss  xmm5, [rdi + rcx*4 + 20]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 12], xmm3
	movss  [rdi + rcx*4 + 16], xmm4
	movss  [rdi + rcx*4 + 20], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb332_fixH2]
	movaps xmm1, [rsp + nb332_fiyH2]
	movaps xmm2, [rsp + nb332_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 i xmm0-xmm2 

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
	movss  xmm3, [rdi + rcx*4 + 24]
	movss  xmm4, [rdi + rcx*4 + 28]
	movss  xmm5, [rdi + rcx*4 + 32]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 24], xmm3
	movss  [rdi + rcx*4 + 28], xmm4
	movss  [rdi + rcx*4 + 32], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, [rsi + rdx*4]
	movss  xmm4, [rsi + rdx*4 + 8]
	subps  xmm3, xmm6
	subss  xmm4, xmm7
	movlps  [rsi + rdx*4],    xmm3
	movss  [rsi + rdx*4 + 8], xmm4

	;# get n from stack
	mov esi, [rsp + nb332_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb332_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb332_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb332_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb332_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb332_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
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



	

	
.globl nb_kernel332nf_x86_64_sse
.globl _nb_kernel332nf_x86_64_sse
nb_kernel332nf_x86_64_sse:	
_nb_kernel332nf_x86_64_sse:	
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
    
	emms

	sub rsp, 832		;# local variable stack space (n*16+8)
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
	movss xmm0, [rsi]
	movss [rsp + nb332nf_facel], xmm0

	mov rax, [rbp + nb332nf_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb332nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb332nf_half], eax
	movss xmm1, [rsp + nb332nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb332nf_half],  xmm1
	movaps [rsp + nb332nf_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb332nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb332nf_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, xmm3	
	movss xmm5, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb332nf_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb332nf_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb332nf_qqOO], xmm3
	movaps [rsp + nb332nf_qqOH], xmm4
	movaps [rsp + nb332nf_qqHH], xmm5
	
	xorps xmm0, xmm0
	mov   rdx, [rbp + nb332nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb332nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   rax, [rbp + nb332nf_vdwparam]
	movlps xmm0, [rax + rdx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 85   ;# 01010101
	movaps [rsp + nb332nf_c6], xmm0
	movaps [rsp + nb332nf_c12], xmm1

.nb332nf_threadloop:
        mov   rsi, [rbp + nb332nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb332nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
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
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb332nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb332nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb332nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb332nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb332nf_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb332nf_ixO], xmm3
	movaps [rsp + nb332nf_iyO], xmm4
	movaps [rsp + nb332nf_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 12]
	addss xmm1, [rax + rbx*4 + 16]
	addss xmm2, [rax + rbx*4 + 20]		
	addss xmm3, [rax + rbx*4 + 24]
	addss xmm4, [rax + rbx*4 + 28]
	addss xmm5, [rax + rbx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb332nf_ixH1], xmm0
	movaps [rsp + nb332nf_iyH1], xmm1
	movaps [rsp + nb332nf_izH1], xmm2
	movaps [rsp + nb332nf_ixH2], xmm3
	movaps [rsp + nb332nf_iyH2], xmm4
	movaps [rsp + nb332nf_izH2], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb332nf_vctot], xmm4
	movaps [rsp + nb332nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb332nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb332nf_pos]
	mov   rax, [rsp + nb332nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb332nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb332nf_ninner]
	mov   [rsp + nb332nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb332nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb332nf_unroll_loop
	jmp   .nb332nf_single_check
.nb332nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb332nf_innerjjnr]     ;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb332nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb332nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2]     ;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	
	
	;# move j coordinates to local temp variables 
	movlps xmm2, [rsi + rax*4]
	movlps xmm3, [rsi + rax*4 + 12]
	movlps xmm4, [rsi + rax*4 + 24]

	movlps xmm5, [rsi + rbx*4]
	movlps xmm6, [rsi + rbx*4 + 12]
	movlps xmm7, [rsi + rbx*4 + 24]

	movhps xmm2, [rsi + rcx*4]
	movhps xmm3, [rsi + rcx*4 + 12]
	movhps xmm4, [rsi + rcx*4 + 24]

	movhps xmm5, [rsi + rdx*4]
	movhps xmm6, [rsi + rdx*4 + 12]
	movhps xmm7, [rsi + rdx*4 + 24]

	;# current state: 	
	;# xmm2= jxOa  jyOa  jxOc  jyOc 
	;# xmm3= jxH1a jyH1a jxH1c jyH1c 
	;# xmm4= jxH2a jyH2a jxH2c jyH2c 
	;# xmm5= jxOb  jyOb  jxOd  jyOd 
	;# xmm6= jxH1b jyH1b jxH1d jyH1d 
	;# xmm7= jxH2b jyH2b jxH2d jyH2d 
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5	;# xmm0= jxOa  jxOb  jyOa  jyOb 
	unpcklps xmm1, xmm6	;# xmm1= jxH1a jxH1b jyH1a jyH1b 
	unpckhps xmm2, xmm5	;# xmm2= jxOc  jxOd  jyOc  jyOd 
	unpckhps xmm3, xmm6	;# xmm3= jxH1c jxH1d jyH1c jyH1d  
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7	;# xmm4= jxH2a jxH2b jyH2a jyH2b 		
	unpckhps xmm5, xmm7	;# xmm5= jxH2c jxH2d jyH2c jyH2d 
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2	;# xmm0= jxOa  jxOb  jxOc  jxOd  
	movaps [rsp + nb332nf_jxO], xmm0
	movhlps  xmm2, xmm6	;# xmm2= jyOa  jyOb  jyOc  jyOd 
	movaps [rsp + nb332nf_jyO], xmm2
	movlhps  xmm1, xmm3
	movaps [rsp + nb332nf_jxH1], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [rsp + nb332nf_jyH1], xmm3
	movlhps  xmm4, xmm5
	movaps [rsp + nb332nf_jxH2], xmm4
	movhlps  xmm5, xmm6
	movaps [rsp + nb332nf_jyH2], xmm5

	movss  xmm0, [rsi + rax*4 + 8]
	movss  xmm1, [rsi + rax*4 + 20]
	movss  xmm2, [rsi + rax*4 + 32]

	movss  xmm3, [rsi + rcx*4 + 8]
	movss  xmm4, [rsi + rcx*4 + 20]
	movss  xmm5, [rsi + rcx*4 + 32]

	movhps xmm0, [rsi + rbx*4 + 4]
	movhps xmm1, [rsi + rbx*4 + 16]
	movhps xmm2, [rsi + rbx*4 + 28]
	
	movhps xmm3, [rsi + rdx*4 + 4]
	movhps xmm4, [rsi + rdx*4 + 16]
	movhps xmm5, [rsi + rdx*4 + 28]
	
	shufps xmm0, xmm3, 204  ;# 11001100
	shufps xmm1, xmm4, 204  ;# 11001100
	shufps xmm2, xmm5, 204  ;# 11001100
	movaps [rsp + nb332nf_jzO],  xmm0
	movaps [rsp + nb332nf_jzH1],  xmm1
	movaps [rsp + nb332nf_jzH2],  xmm2

	movaps xmm0, [rsp + nb332nf_ixO]
	movaps xmm1, [rsp + nb332nf_iyO]
	movaps xmm2, [rsp + nb332nf_izO]
	movaps xmm3, [rsp + nb332nf_ixO]
	movaps xmm4, [rsp + nb332nf_iyO]
	movaps xmm5, [rsp + nb332nf_izO]
	subps  xmm0, [rsp + nb332nf_jxO]
	subps  xmm1, [rsp + nb332nf_jyO]
	subps  xmm2, [rsp + nb332nf_jzO]
	subps  xmm3, [rsp + nb332nf_jxH1]
	subps  xmm4, [rsp + nb332nf_jyH1]
	subps  xmm5, [rsp + nb332nf_jzH1]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb332nf_rsqOO], xmm0
	movaps [rsp + nb332nf_rsqOH1], xmm3

	movaps xmm0, [rsp + nb332nf_ixO]
	movaps xmm1, [rsp + nb332nf_iyO]
	movaps xmm2, [rsp + nb332nf_izO]
	movaps xmm3, [rsp + nb332nf_ixH1]
	movaps xmm4, [rsp + nb332nf_iyH1]
	movaps xmm5, [rsp + nb332nf_izH1]
	subps  xmm0, [rsp + nb332nf_jxH2]
	subps  xmm1, [rsp + nb332nf_jyH2]
	subps  xmm2, [rsp + nb332nf_jzH2]
	subps  xmm3, [rsp + nb332nf_jxO]
	subps  xmm4, [rsp + nb332nf_jyO]
	subps  xmm5, [rsp + nb332nf_jzO]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb332nf_rsqOH2], xmm0
	movaps [rsp + nb332nf_rsqH1O], xmm3

	movaps xmm0, [rsp + nb332nf_ixH1]
	movaps xmm1, [rsp + nb332nf_iyH1]
	movaps xmm2, [rsp + nb332nf_izH1]
	movaps xmm3, [rsp + nb332nf_ixH1]
	movaps xmm4, [rsp + nb332nf_iyH1]
	movaps xmm5, [rsp + nb332nf_izH1]
	subps  xmm0, [rsp + nb332nf_jxH1]
	subps  xmm1, [rsp + nb332nf_jyH1]
	subps  xmm2, [rsp + nb332nf_jzH1]
	subps  xmm3, [rsp + nb332nf_jxH2]
	subps  xmm4, [rsp + nb332nf_jyH2]
	subps  xmm5, [rsp + nb332nf_jzH2]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb332nf_rsqH1H1], xmm0
	movaps [rsp + nb332nf_rsqH1H2], xmm3

	movaps xmm0, [rsp + nb332nf_ixH2]
	movaps xmm1, [rsp + nb332nf_iyH2]
	movaps xmm2, [rsp + nb332nf_izH2]
	movaps xmm3, [rsp + nb332nf_ixH2]
	movaps xmm4, [rsp + nb332nf_iyH2]
	movaps xmm5, [rsp + nb332nf_izH2]
	subps  xmm0, [rsp + nb332nf_jxO]
	subps  xmm1, [rsp + nb332nf_jyO]
	subps  xmm2, [rsp + nb332nf_jzO]
	subps  xmm3, [rsp + nb332nf_jxH1]
	subps  xmm4, [rsp + nb332nf_jyH1]
	subps  xmm5, [rsp + nb332nf_jzH1]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [rsp + nb332nf_rsqH2O], xmm0
	movaps [rsp + nb332nf_rsqH2H1], xmm4

	movaps xmm0, [rsp + nb332nf_ixH2]
	movaps xmm1, [rsp + nb332nf_iyH2]
	movaps xmm2, [rsp + nb332nf_izH2]
	subps  xmm0, [rsp + nb332nf_jxH2]
	subps  xmm1, [rsp + nb332nf_jyH2]
	subps  xmm2, [rsp + nb332nf_jzH2]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [rsp + nb332nf_rsqH2H2], xmm0
	
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb332nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb332nf_half] ;# rinvH2H2 
	mulps   xmm7, [rsp + nb332nf_half] ;# rinvH2H1 
	movaps  [rsp + nb332nf_rinvH2H2], xmm3
	movaps  [rsp + nb332nf_rinvH2H1], xmm7
	
	rsqrtps xmm1, [rsp + nb332nf_rsqOO]
	rsqrtps xmm5, [rsp + nb332nf_rsqOH1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb332nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb332nf_rsqOO]
	mulps   xmm5, [rsp + nb332nf_rsqOH1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb332nf_half] 
	mulps   xmm7, [rsp + nb332nf_half]
	movaps  [rsp + nb332nf_rinvOO], xmm3
	movaps  [rsp + nb332nf_rinvOH1], xmm7
	
	rsqrtps xmm1, [rsp + nb332nf_rsqOH2]
	rsqrtps xmm5, [rsp + nb332nf_rsqH1O]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb332nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb332nf_rsqOH2]
	mulps   xmm5, [rsp + nb332nf_rsqH1O]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb332nf_half] 
	mulps   xmm7, [rsp + nb332nf_half]
	movaps  [rsp + nb332nf_rinvOH2], xmm3
	movaps  [rsp + nb332nf_rinvH1O], xmm7
	
	rsqrtps xmm1, [rsp + nb332nf_rsqH1H1]
	rsqrtps xmm5, [rsp + nb332nf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb332nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb332nf_rsqH1H1]
	mulps   xmm5, [rsp + nb332nf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb332nf_half] 
	mulps   xmm7, [rsp + nb332nf_half]
	movaps  [rsp + nb332nf_rinvH1H1], xmm3
	movaps  [rsp + nb332nf_rinvH1H2], xmm7
	
	rsqrtps xmm1, [rsp + nb332nf_rsqH2O]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb332nf_three]
	mulps   xmm1, [rsp + nb332nf_rsqH2O]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb332nf_half] 
	movaps  [rsp + nb332nf_rinvH2O], xmm3

	;# start with OO interaction 
	movaps xmm0, [rsp + nb332nf_rinvOO]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb332nf_rsqOO] ;# xmm1=r 
	mulps  xmm1, [rsp + nb332nf_tsc]
	
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 2
    pslld mm7, 2
	
    movd mm0, eax
    movd mm1, ebx
    movd mm2, ecx
    movd mm3, edx

    mov  rsi, [rbp + nb332nf_VFtab]
    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    lea   rax, [rax + rax*2]
    lea   rbx, [rbx + rbx*2]
    lea   rcx, [rcx + rcx*2]
    lea   rdx, [rdx + rdx*2]

    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  

    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    movaps xmm3, [rsp + nb332nf_qqOO]
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
    ;# increment vcoul - then we can get rid of mm5 
    ;# update vctot 
    addps  xmm5, [rsp + nb332nf_vctot]
    movaps [rsp + nb332nf_vctot], xmm5 

    ;# dispersion 
    movlps xmm5, [rsi + rax*4 + 16]
    movlps xmm7, [rsi + rcx*4 + 16]
    movhps xmm5, [rsi + rbx*4 + 16]
    movhps xmm7, [rsi + rdx*4 + 16] ;# got half dispersion table 
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 24]
    movlps xmm3, [rsi + rcx*4 + 24]
    movhps xmm7, [rsi + rbx*4 + 24]
    movhps xmm3, [rsi + rdx*4 + 24] ;# other half of dispersion table 
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# dispersion table ready, in xmm4-xmm7 
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 

    movaps xmm4, [rsp + nb332nf_c6]
    mulps  xmm5, xmm4    ;# Vvdw6 

    ;# put scalar force on stack Update Vvdwtot directly 
    addps  xmm5, [rsp + nb332nf_Vvdwtot]
    movaps [rsp + nb332nf_Vvdwtot], xmm5

    ;# repulsion 
    movlps xmm5, [rsi + rax*4 + 32]
    movlps xmm7, [rsi + rcx*4 + 32]
    movhps xmm5, [rsi + rbx*4 + 32]
    movhps xmm7, [rsi + rdx*4 + 32] ;# got half repulsion table 
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 40]
    movlps xmm3, [rsi + rcx*4 + 40]
    movhps xmm7, [rsi + rbx*4 + 40]
    movhps xmm3, [rsi + rdx*4 + 40] ;# other half of repulsion table 
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# table ready, in xmm4-xmm7 
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
 
    movaps xmm4, [rsp + nb332nf_c12]
    mulps  xmm5, xmm4 ;# Vvdw12 
    addps  xmm5, [rsp + nb332nf_Vvdwtot]
    movaps [rsp + nb332nf_Vvdwtot], xmm5
	
	;# O-H1 interaction 
	movaps xmm0, [rsp + nb332nf_rinvOH1]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb332nf_rsqOH1] ;# xmm1=r 
	mulps  xmm1, [rsp + nb332nf_tsc]	
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 2
    pslld mm7, 2

    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    lea   rax, [rax + rax*2]
    lea   rbx, [rbx + rbx*2]
    lea   rcx, [rcx + rcx*2]
    lea   rdx, [rdx + rdx*2]

    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  

    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    movaps xmm3, [rsp + nb332nf_qqOH]
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV 
    ;# at this point mm5 contains vcoul 

    addps  xmm5, [rsp + nb332nf_vctot]
    movaps [rsp + nb332nf_vctot], xmm5
	
	;# O-H2 interaction  
	movaps xmm0, [rsp + nb332nf_rinvOH2]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb332nf_rsqOH2] ;# xmm1=r 
	mulps  xmm1, [rsp + nb332nf_tsc]	
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 2
    pslld mm7, 2

    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    lea   rax, [rax + rax*2]
    lea   rbx, [rbx + rbx*2]
    lea   rcx, [rcx + rcx*2]
    lea   rdx, [rdx + rdx*2]

    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  

    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    movaps xmm3, [rsp + nb332nf_qqOH]
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV 
    ;# at this point mm5 contains vcoul 

    addps  xmm5, [rsp + nb332nf_vctot]
    movaps [rsp + nb332nf_vctot], xmm5
	
	;# H1-O interaction 
	movaps xmm0, [rsp + nb332nf_rinvH1O]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb332nf_rsqH1O] ;# xmm1=r 
	mulps  xmm1, [rsp + nb332nf_tsc]	
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 2
    pslld mm7, 2

    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    lea   rax, [rax + rax*2]
    lea   rbx, [rbx + rbx*2]
    lea   rcx, [rcx + rcx*2]
    lea   rdx, [rdx + rdx*2]

    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  

    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    movaps xmm3, [rsp + nb332nf_qqOH]
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV 
    ;# at this point mm5 contains vcoul 

    addps  xmm5, [rsp + nb332nf_vctot]
    movaps [rsp + nb332nf_vctot], xmm5
	
	;# H1-H1 interaction 
	movaps xmm0, [rsp + nb332nf_rinvH1H1]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb332nf_rsqH1H1] ;# xmm1=r 
	mulps  xmm1, [rsp + nb332nf_tsc]	
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 2
    pslld mm7, 2

    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    lea   rax, [rax + rax*2]
    lea   rbx, [rbx + rbx*2]
    lea   rcx, [rcx + rcx*2]
    lea   rdx, [rdx + rdx*2]

    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  

    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    movaps xmm3, [rsp + nb332nf_qqHH]
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV 
    ;# at this point mm5 contains vcoul 

    addps  xmm5, [rsp + nb332nf_vctot]
    movaps [rsp + nb332nf_vctot], xmm5
	
	;# H1-H2 interaction 
	movaps xmm0, [rsp + nb332nf_rinvH1H2]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb332nf_rsqH1H2] ;# xmm1=r 
	mulps  xmm1, [rsp + nb332nf_tsc]
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 2
    pslld mm7, 2

    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    lea   rax, [rax + rax*2]
    lea   rbx, [rbx + rbx*2]
    lea   rcx, [rcx + rcx*2]
    lea   rdx, [rdx + rdx*2]

    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  

    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    movaps xmm3, [rsp + nb332nf_qqHH]
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV 
    ;# at this point mm5 contains vcoul 

    addps  xmm5, [rsp + nb332nf_vctot]
    movaps [rsp + nb332nf_vctot], xmm5
	
	;# H2-O interaction 
	movaps xmm0, [rsp + nb332nf_rinvH2O]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb332nf_rsqH2O] ;# xmm1=r 
	mulps  xmm1, [rsp + nb332nf_tsc]	
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 2
    pslld mm7, 2

    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    lea   rax, [rax + rax*2]
    lea   rbx, [rbx + rbx*2]
    lea   rcx, [rcx + rcx*2]
    lea   rdx, [rdx + rdx*2]

    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  

    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    movaps xmm3, [rsp + nb332nf_qqOH]
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addps  xmm5, [rsp + nb332nf_vctot]
    movaps [rsp + nb332nf_vctot], xmm5

	;# H2-H1 interaction 
	movaps xmm0, [rsp + nb332nf_rinvH2H1]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb332nf_rsqH2H1] ;# xmm1=r 
	mulps  xmm1, [rsp + nb332nf_tsc]
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 2
    pslld mm7, 2

    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    lea   rax, [rax + rax*2]
    lea   rbx, [rbx + rbx*2]
    lea   rcx, [rcx + rcx*2]
    lea   rdx, [rdx + rdx*2]

    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  

    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    movaps xmm3, [rsp + nb332nf_qqHH]
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addps  xmm5, [rsp + nb332nf_vctot]
    movaps [rsp + nb332nf_vctot], xmm5

	;# H2-H2 interaction 
	movaps xmm0, [rsp + nb332nf_rinvH2H2]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb332nf_rsqH2H2] ;# xmm1=r 
	mulps  xmm1, [rsp + nb332nf_tsc]	
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 2
    pslld mm7, 2

    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    lea   rax, [rax + rax*2]
    lea   rbx, [rbx + rbx*2]
    lea   rcx, [rcx + rcx*2]
    lea   rdx, [rdx + rdx*2]

    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  

    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    movaps xmm3, [rsp + nb332nf_qqHH]
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addps  xmm5, [rsp + nb332nf_vctot]
    movaps [rsp + nb332nf_vctot], xmm5	
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb332nf_innerk],  4
	jl    .nb332nf_single_check
	jmp   .nb332nf_unroll_loop
.nb332nf_single_check:
	add dword ptr [rsp + nb332nf_innerk],  4
	jnz   .nb332nf_single_loop
	jmp   .nb332nf_updateouterdata
.nb332nf_single_loop:
	mov   rdx, [rsp + nb332nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb332nf_innerjjnr],  4	

	mov rsi, [rbp + nb332nf_pos]
	lea   rax, [rax + rax*2]  

	;# fetch j coordinates 
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	
	movss xmm3, [rsi + rax*4]		;# jxO  -  -  -
	movss xmm4, [rsi + rax*4 + 4]		;# jyO  -  -  -
	movss xmm5, [rsi + rax*4 + 8]		;# jzO  -  -  -  

	movlps xmm6, [rsi + rax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, [rsi + rax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, [rsi + rax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm2, [rsi + rax*4 + 32]		;# xmm2 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm2			;# xmm7 = jzH1 jzH2   -    -
	movaps  xmm0, [rsp + nb332nf_ixO]     
	movaps  xmm1, [rsp + nb332nf_iyO]
	movaps  xmm2, [rsp + nb332nf_izO]	
	movlhps xmm3, xmm6			;# xmm3 = jxO   0   jxH1 jxH2 
	shufps  xmm4, xmm6, 228 ;# 11100100	;# xmm4 = jyO   0   jyH1 jyH2 
	shufps  xmm5, xmm7, 68  ;# 01000100	;# xmm5 = jzO   0   jzH1 jzH2

	;# store all j coordinates in jO  
	movaps [rsp + nb332nf_jxO], xmm3
	movaps [rsp + nb332nf_jyO], xmm4
	movaps [rsp + nb332nf_jzO], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1	
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb332nf_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [rsp + nb332nf_half] ;# rinv iO - j water 

	movaps  xmm1, xmm3
	mulps   xmm1, xmm0	;# xmm1=r 
	movaps  xmm0, xmm3	;# xmm0=rinv 
	mulps  xmm1, [rsp + nb332nf_tsc]
	
	movhlps xmm2, xmm1	
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 2
    pslld mm7, 2
    movd ebx, mm6
    movd ecx, mm7
    psrlq mm7, 32
    movd edx, mm7		;# table indices in ebx,ecx,edx 

	mov rsi, [rbp + nb332nf_VFtab]
	
    lea   rbx, [rbx + rbx*2]
    lea   rcx, [rcx + rcx*2]
    lea   rdx, [rdx + rdx*2]
	
    movlps xmm5, [rsi + rbx*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rbx*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb332nf_qqOO]
	movhps  xmm3, [rsp + nb332nf_qqOH]
	
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV 
    ;# at this point xmm5 contains vcoul 
	
    addps  xmm5, [rsp + nb332nf_vctot]
    movaps [rsp + nb332nf_vctot], xmm5

    ;# dispersion 
	movss  xmm4, [rsi + rbx*4 + 16]	
	movss  xmm5, [rsi + rbx*4 + 20]	
	movss  xmm6, [rsi + rbx*4 + 24]	
	movss  xmm7, [rsi + rbx*4 + 28]
    ;# dispersion table ready, in xmm4-xmm7 
    mulss  xmm6, xmm1       ;# xmm6=Geps 
    mulss  xmm7, xmm2       ;# xmm7=Heps2 
    addss  xmm5, xmm6
    addss  xmm5, xmm7       ;# xmm5=Fp 
    mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
    addss  xmm5, xmm4 ;# xmm5=VV 
	xorps  xmm4, xmm4
    movss  xmm4, [rsp + nb332nf_c6]
    mulps  xmm5, xmm4    ;# Vvdw6 
    ;# put scalar force on stack 
    addps  xmm5, [rsp + nb332nf_Vvdwtot]
    movaps [rsp + nb332nf_Vvdwtot], xmm5

    ;# repulsion 
	movss  xmm4, [rsi + rbx*4 + 32]	
	movss  xmm5, [rsi + rbx*4 + 36]	
	movss  xmm6, [rsi + rbx*4 + 40]	
	movss  xmm7, [rsi + rbx*4 + 44]
    ;# table ready, in xmm4-xmm7 
    mulss  xmm6, xmm1       ;# xmm6=Geps 
    mulss  xmm7, xmm2       ;# xmm7=Heps2 
    addss  xmm5, xmm6
    addss  xmm5, xmm7       ;# xmm5=Fp 
    mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
    addss  xmm5, xmm4 ;# xmm5=VV 

	xorps  xmm4, xmm4
    movss  xmm4, [rsp + nb332nf_c12]
    mulps  xmm5, xmm4 ;# Vvdw12 
    addps  xmm5, [rsp + nb332nf_Vvdwtot]
    movaps [rsp + nb332nf_Vvdwtot], xmm5

	
	;# done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [rsp + nb332nf_ixH1]
	movaps  xmm1, [rsp + nb332nf_iyH1]
	movaps  xmm2, [rsp + nb332nf_izH1]	
	movaps  xmm3, [rsp + nb332nf_ixH2] 
	movaps  xmm4, [rsp + nb332nf_iyH2] 
	movaps  xmm5, [rsp + nb332nf_izH2] 
	subps   xmm0, [rsp + nb332nf_jxO]
	subps   xmm1, [rsp + nb332nf_jyO]
	subps   xmm2, [rsp + nb332nf_jzO]
	subps   xmm3, [rsp + nb332nf_jxO]
	subps   xmm4, [rsp + nb332nf_jyO]
	subps   xmm5, [rsp + nb332nf_jzO]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;# have rsqH1 in xmm0 
	addps xmm4, xmm5	;# have rsqH2 in xmm4 

	;# start with H1, save H2 data 
	movaps [rsp + nb332nf_rsqH2O], xmm4
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb332nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb332nf_half] ;# rinv H1 - j water 
	mulps   xmm7, [rsp + nb332nf_half] ;# rinv H2 - j water  

	;# start with H1, save H2 data 
	movaps [rsp + nb332nf_rinvH2O], xmm7

	movaps xmm1, xmm3
	mulps  xmm1, xmm0	;# xmm1=r 
	movaps xmm0, xmm3	;# xmm0=rinv 
	mulps  xmm1, [rsp + nb332nf_tsc]
	
	movhlps xmm2, xmm1	
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 2
    pslld mm7, 2
    movd ebx, mm6
    movd ecx, mm7
    psrlq mm7, 32
    movd edx, mm7		;# table indices in ebx,ecx,edx 

    lea   rbx, [rbx + rbx*2]
    lea   rcx, [rcx + rcx*2]
    lea   rdx, [rdx + rdx*2]
	
    movlps xmm5, [rsi + rbx*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rbx*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb332nf_qqOH]
	movhps  xmm3, [rsp + nb332nf_qqHH]
	
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV 
    ;# at this point xmm5 contains vcoul 
    addps  xmm5, [rsp + nb332nf_vctot]
    movaps [rsp + nb332nf_vctot], xmm5	

	
	;# do table for H2 - j water interaction 
	movaps xmm0, [rsp + nb332nf_rinvH2O]
	movaps xmm1, [rsp + nb332nf_rsqH2O]
	mulps  xmm1, xmm0	;# xmm0=rinv, xmm1=r 
	mulps  xmm1, [rsp + nb332nf_tsc]
	
	movhlps xmm2, xmm1	
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
    pslld mm6, 2
    pslld mm7, 2
    movd ebx, mm6
    movd ecx, mm7
    psrlq mm7, 32
    movd edx, mm7		;# table indices in ebx,ecx,edx 

    lea   rbx, [rbx + rbx*2]
    lea   rcx, [rcx + rcx*2]
    lea   rdx, [rdx + rdx*2]
	
    movlps xmm5, [rsi + rbx*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rbx*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb332nf_qqOH]
	movhps  xmm3, [rsp + nb332nf_qqHH]
	
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point xmm5 contains vcoul 
    addps  xmm5, [rsp + nb332nf_vctot]
    movaps [rsp + nb332nf_vctot], xmm5	
	
	dec dword ptr [rsp + nb332nf_innerk]
	jz    .nb332nf_updateouterdata
	jmp   .nb332nf_single_loop
.nb332nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb332nf_n]
    	;# get group index for i particle 
    	mov   rdx, [rbp + nb332nf_gid]      	;# base of gid[]
    	mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb332nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb332nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb332nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb332nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
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
