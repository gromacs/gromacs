##
##
## Gromacs 4.0                         Copyright (c) 1991-2003 
## David van der Spoel, Erik Lindahl
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## To help us fund GROMACS development, we humbly ask that you cite
## the research papers on the package. Check out http://www.gromacs.org
## 
## And Hey:
## Gnomes, ROck Monsters And Chili Sauce
##





.globl nb_kernel101_x86_64_sse
.globl _nb_kernel101_x86_64_sse
nb_kernel101_x86_64_sse:        
_nb_kernel101_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb101_fshift, 16
.set nb101_gid, 24
.set nb101_pos, 32
.set nb101_faction, 40
.set nb101_charge, 48
.set nb101_p_facel, 56
.set nb101_argkrf, 64
.set nb101_argcrf, 72
.set nb101_Vc, 80
.set nb101_type, 88
.set nb101_p_ntype, 96
.set nb101_vdwparam, 104
.set nb101_Vvdw, 112
.set nb101_p_tabscale, 120
.set nb101_VFtab, 128
.set nb101_invsqrta, 136
.set nb101_dvda, 144
.set nb101_p_gbtabscale, 152
.set nb101_GBtab, 160
.set nb101_p_nthreads, 168
.set nb101_count, 176
.set nb101_mtx, 184
.set nb101_outeriter, 192
.set nb101_inneriter, 200
.set nb101_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb101_ixO, 0
.set nb101_iyO, 16
.set nb101_izO, 32
.set nb101_ixH1, 48
.set nb101_iyH1, 64
.set nb101_izH1, 80
.set nb101_ixH2, 96
.set nb101_iyH2, 112
.set nb101_izH2, 128
.set nb101_iqO, 144
.set nb101_iqH, 160
.set nb101_dxO, 176
.set nb101_dyO, 192
.set nb101_dzO, 208
.set nb101_dxH1, 224
.set nb101_dyH1, 240
.set nb101_dzH1, 256
.set nb101_dxH2, 272
.set nb101_dyH2, 288
.set nb101_dzH2, 304
.set nb101_qqO, 320
.set nb101_qqH, 336
.set nb101_vctot, 352
.set nb101_fixO, 368
.set nb101_fiyO, 384
.set nb101_fizO, 400
.set nb101_fixH1, 416
.set nb101_fiyH1, 432
.set nb101_fizH1, 448
.set nb101_fixH2, 464
.set nb101_fiyH2, 480
.set nb101_fizH2, 496
.set nb101_fjx, 512
.set nb101_fjy, 528
.set nb101_fjz, 544
.set nb101_half, 560
.set nb101_three, 576
.set nb101_nri, 592
.set nb101_innerjjnr, 600
.set nb101_iinr, 608
.set nb101_jindex, 616
.set nb101_jjnr, 624
.set nb101_shift, 632
.set nb101_shiftvec, 640
.set nb101_facel, 648
.set nb101_is3, 656
.set nb101_ii3, 660
.set nb101_innerk, 664
.set nb101_n, 668
.set nb101_nn1, 672
.set nb101_nouter, 676
.set nb101_ninner, 680

        push %rbp
        movq %rsp,%rbp
        push %rbx

        push %r12
        push %r13
        push %r14
        push %r15


        subq $696,%rsp
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb101_nouter(%rsp)
        movl %eax,nb101_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb101_nri(%rsp)
        movq %rsi,nb101_iinr(%rsp)
        movq %rdx,nb101_jindex(%rsp)
        movq %rcx,nb101_jjnr(%rsp)
        movq %r8,nb101_shift(%rsp)
        movq %r9,nb101_shiftvec(%rsp)
        movq nb101_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb101_facel(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb101_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb101_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss 4(%rdx,%rbx,4),%xmm4
        movq nb101_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb101_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb101_iqO(%rsp)
        movaps %xmm4,nb101_iqH(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb101_half(%rsp)
        movss nb101_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb101_half(%rsp)
        movaps %xmm3,nb101_three(%rsp)

_nb_kernel101_x86_64_sse.nb101_threadloop: 
        movq  nb101_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel101_x86_64_sse.nb101_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel101_x86_64_sse.nb101_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb101_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb101_n(%rsp)
        movl %ebx,nb101_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel101_x86_64_sse.nb101_outerstart
        jmp _nb_kernel101_x86_64_sse.nb101_end

_nb_kernel101_x86_64_sse.nb101_outerstart: 
        ## ebx contains number of outer iterations
        addl nb101_nouter(%rsp),%ebx
        movl %ebx,nb101_nouter(%rsp)

_nb_kernel101_x86_64_sse.nb101_outer: 
        movq  nb101_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb101_is3(%rsp)      ## store is3 

        movq  nb101_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb101_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb101_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb101_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb101_ixO(%rsp)
        movaps %xmm4,nb101_iyO(%rsp)
        movaps %xmm5,nb101_izO(%rsp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%rax,%rbx,4),%xmm0
        addss 16(%rax,%rbx,4),%xmm1
        addss 20(%rax,%rbx,4),%xmm2
        addss 24(%rax,%rbx,4),%xmm3
        addss 28(%rax,%rbx,4),%xmm4
        addss 32(%rax,%rbx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb101_ixH1(%rsp)
        movaps %xmm1,nb101_iyH1(%rsp)
        movaps %xmm2,nb101_izH1(%rsp)
        movaps %xmm3,nb101_ixH2(%rsp)
        movaps %xmm4,nb101_iyH2(%rsp)
        movaps %xmm5,nb101_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb101_vctot(%rsp)
        movaps %xmm4,nb101_fixO(%rsp)
        movaps %xmm4,nb101_fiyO(%rsp)
        movaps %xmm4,nb101_fizO(%rsp)
        movaps %xmm4,nb101_fixH1(%rsp)
        movaps %xmm4,nb101_fiyH1(%rsp)
        movaps %xmm4,nb101_fizH1(%rsp)
        movaps %xmm4,nb101_fixH2(%rsp)
        movaps %xmm4,nb101_fiyH2(%rsp)
        movaps %xmm4,nb101_fizH2(%rsp)

        movq  nb101_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb101_pos(%rbp),%rsi
        movq  nb101_faction(%rbp),%rdi
        movq  nb101_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb101_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb101_ninner(%rsp),%ecx
        movl  %ecx,nb101_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb101_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel101_x86_64_sse.nb101_unroll_loop
        jmp   _nb_kernel101_x86_64_sse.nb101_odd_inner
_nb_kernel101_x86_64_sse.nb101_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb101_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%r8d
        movl  4(%rdx),%r9d
        movl  8(%rdx),%r10d
        movl  12(%rdx),%r11d           ## eax-edx=jnr1-4 

        addq $16,nb101_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb101_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%r8,4),%xmm12
        movss (%rsi,%r10,4),%xmm13
        movss (%rsi,%r9,4),%xmm14
        movss (%rsi,%r11,4),%xmm15

        shufps $0,%xmm14,%xmm12
        shufps $0,%xmm15,%xmm13
        shufps $136,%xmm13,%xmm12 ## 10001000 ;# all charges in xmm3  
        movaps %xmm12,%xmm13         ## and in xmm4 
        mulps  nb101_iqO(%rsp),%xmm12
        mulps  nb101_iqH(%rsp),%xmm13

        movaps  %xmm12,nb101_qqO(%rsp)
        movaps  %xmm13,nb101_qqH(%rsp)

        movq nb101_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%r8,%r8,2),%rax     ## j3 
        lea  (%r9,%r9,2),%rbx
        lea  (%r10,%r10,2),%rcx
        lea  (%r11,%r11,2),%rdx

        ## move four j coordinates to xmm0-xmm2         
        movlps (%rsi,%rax,4),%xmm4
        movlps (%rsi,%rcx,4),%xmm5
        movss 8(%rsi,%rax,4),%xmm2
        movss 8(%rsi,%rcx,4),%xmm6

        movhps (%rsi,%rbx,4),%xmm4
        movhps (%rsi,%rdx,4),%xmm5

        movss 8(%rsi,%rbx,4),%xmm0
        movss 8(%rsi,%rdx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## 10001000

        shufps $136,%xmm5,%xmm0 ## 10001000
        shufps $221,%xmm5,%xmm1 ## 11011101             

    ## xmm0 = jx
    ## xmm1 = jy
    ## xmm2 = jz

    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8

    subps nb101_ixO(%rsp),%xmm0
    subps nb101_iyO(%rsp),%xmm1
    subps nb101_izO(%rsp),%xmm2
    subps nb101_ixH1(%rsp),%xmm3
    subps nb101_iyH1(%rsp),%xmm4
    subps nb101_izH1(%rsp),%xmm5
    subps nb101_ixH2(%rsp),%xmm6
    subps nb101_iyH2(%rsp),%xmm7
    subps nb101_izH2(%rsp),%xmm8

        movaps %xmm0,nb101_dxO(%rsp)
        movaps %xmm1,nb101_dyO(%rsp)
        movaps %xmm2,nb101_dzO(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb101_dxH1(%rsp)
        movaps %xmm4,nb101_dyH1(%rsp)
        movaps %xmm5,nb101_dzH1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb101_dxH2(%rsp)
        movaps %xmm7,nb101_dyH2(%rsp)
        movaps %xmm8,nb101_dzH2(%rsp)
        mulps  %xmm6,%xmm6
        mulps  %xmm7,%xmm7
        mulps  %xmm8,%xmm8
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
    addps  %xmm7,%xmm6
    addps  %xmm8,%xmm6

        ## start doing invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm3,%xmm4
    rsqrtps %xmm6,%xmm7

        movaps  %xmm1,%xmm2
        movaps  %xmm4,%xmm5
    movaps  %xmm7,%xmm8

        mulps   %xmm1,%xmm1 ## lu*lu
        mulps   %xmm4,%xmm4 ## lu*lu
    mulps   %xmm7,%xmm7 ## lu*lu

        movaps  nb101_three(%rsp),%xmm9
        movaps  %xmm9,%xmm10
    movaps  %xmm9,%xmm11

        mulps   %xmm0,%xmm1 ## rsq*lu*lu
        mulps   %xmm3,%xmm4 ## rsq*lu*lu 
    mulps   %xmm6,%xmm7 ## rsq*lu*lu

        subps   %xmm1,%xmm9
        subps   %xmm4,%xmm10
    subps   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulps   %xmm2,%xmm9
        mulps   %xmm5,%xmm10
    mulps   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movaps  nb101_half(%rsp),%xmm0
        mulps   %xmm0,%xmm9 ## rinvO 
        mulps   %xmm0,%xmm10 ## rinvH1
    mulps   %xmm0,%xmm11 ## rinvH2

        ## interactions 
    movaps %xmm9,%xmm13
    movaps %xmm10,%xmm14
    movaps %xmm11,%xmm15
    mulps  %xmm9,%xmm9
    mulps  %xmm10,%xmm10
    mulps  %xmm11,%xmm11
    mulps  nb101_qqO(%rsp),%xmm13
    mulps  nb101_qqH(%rsp),%xmm14
    mulps  nb101_qqH(%rsp),%xmm15
    mulps  %xmm13,%xmm9
    mulps  %xmm14,%xmm10
    mulps  %xmm15,%xmm11

    addps nb101_vctot(%rsp),%xmm13
    addps %xmm15,%xmm14
    addps %xmm14,%xmm13
    movaps %xmm13,nb101_vctot(%rsp)

        ## move j forces to local temp variables 
    movlps (%rdi,%rax,4),%xmm0 ## jxa jya  -   -
    movlps (%rdi,%rcx,4),%xmm1 ## jxc jyc  -   -
    movhps (%rdi,%rbx,4),%xmm0 ## jxa jya jxb jyb 
    movhps (%rdi,%rdx,4),%xmm1 ## jxc jyc jxd jyd 

    movss  8(%rdi,%rax,4),%xmm2    ## jza  -  -  -
    movss  8(%rdi,%rcx,4),%xmm3    ## jzc  -  -  -
    movss  8(%rdi,%rbx,4),%xmm4     ## jzb
    movss  8(%rdi,%rdx,4),%xmm5     ## jzd
    movlhps %xmm4,%xmm2 ## jza  -  jzb  -
    movlhps %xmm5,%xmm3 ## jzc  -  jzd -

    shufps $136,%xmm3,%xmm2 ## 10001000 => jza jzb jzc jzd

    ## xmm0: jxa jya jxb jyb 
    ## xmm1: jxc jyc jxd jyd
    ## xmm2: jza jzb jzc jzd

    movaps %xmm9,%xmm7
    movaps %xmm9,%xmm8
    movaps %xmm11,%xmm13
    movaps %xmm11,%xmm14
    movaps %xmm11,%xmm15
    movaps %xmm10,%xmm11
    movaps %xmm10,%xmm12

        mulps nb101_dxO(%rsp),%xmm7
        mulps nb101_dyO(%rsp),%xmm8
        mulps nb101_dzO(%rsp),%xmm9
        mulps nb101_dxH1(%rsp),%xmm10
        mulps nb101_dyH1(%rsp),%xmm11
        mulps nb101_dzH1(%rsp),%xmm12
        mulps nb101_dxH2(%rsp),%xmm13
        mulps nb101_dyH2(%rsp),%xmm14
        mulps nb101_dzH2(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb101_fixO(%rsp),%xmm7
    addps nb101_fiyO(%rsp),%xmm8
    addps nb101_fizO(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb101_fixH1(%rsp),%xmm10
    addps nb101_fiyH1(%rsp),%xmm11
    addps nb101_fizH1(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb101_fixH2(%rsp),%xmm13
    addps nb101_fiyH2(%rsp),%xmm14
    addps nb101_fizH2(%rsp),%xmm15

    movaps %xmm7,nb101_fixO(%rsp)
    movaps %xmm8,nb101_fiyO(%rsp)
    movaps %xmm9,nb101_fizO(%rsp)
    movaps %xmm10,nb101_fixH1(%rsp)
    movaps %xmm11,nb101_fiyH1(%rsp)
    movaps %xmm12,nb101_fizH1(%rsp)
    movaps %xmm13,nb101_fixH2(%rsp)
    movaps %xmm14,nb101_fiyH2(%rsp)
    movaps %xmm15,nb101_fizH2(%rsp)

    ## xmm3 = fjx , xmm4 = fjy
    movaps %xmm3,%xmm5
    unpcklps %xmm4,%xmm3
    unpckhps %xmm4,%xmm5

    addps %xmm3,%xmm0
    addps %xmm5,%xmm1

    movhlps  %xmm2,%xmm3 ## fjzc fjzd

    movlps %xmm0,(%rdi,%rax,4)
    movhps %xmm0,(%rdi,%rbx,4)
    movlps %xmm1,(%rdi,%rcx,4)
    movhps %xmm1,(%rdi,%rdx,4)
    movss  %xmm2,8(%rdi,%rax,4)
    movss  %xmm3,8(%rdi,%rcx,4)
    shufps $1,%xmm2,%xmm2
    shufps $1,%xmm3,%xmm3
    movss  %xmm2,8(%rdi,%rbx,4)
    movss  %xmm3,8(%rdi,%rdx,4)

        ## should we do one more iteration? 
        subl $4,nb101_innerk(%rsp)
        jl    _nb_kernel101_x86_64_sse.nb101_odd_inner
        jmp   _nb_kernel101_x86_64_sse.nb101_unroll_loop
_nb_kernel101_x86_64_sse.nb101_odd_inner: 
        addl $4,nb101_innerk(%rsp)
        jnz   _nb_kernel101_x86_64_sse.nb101_odd_loop
        jmp   _nb_kernel101_x86_64_sse.nb101_updateouterdata
_nb_kernel101_x86_64_sse.nb101_odd_loop: 
        movq  nb101_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb101_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb101_iqO(%rsp),%xmm4
        movq nb101_charge(%rbp),%rsi
        movhps nb101_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb101_qqO(%rsp)    ## use oxygen qq for storage 

        movq nb101_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm3
        movss 4(%rsi,%rax,4),%xmm4
        movss 8(%rsi,%rax,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5

        movss nb101_ixO(%rsp),%xmm0
        movss nb101_iyO(%rsp),%xmm1
        movss nb101_izO(%rsp),%xmm2

        movlps nb101_ixH1(%rsp),%xmm6
        movlps nb101_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm0
        movlps nb101_iyH1(%rsp),%xmm6
        movlps nb101_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm1
        movlps nb101_izH1(%rsp),%xmm6
        movlps nb101_izH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm2

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        movaps %xmm3,nb101_dxO(%rsp)
        movaps %xmm4,nb101_dyO(%rsp)
        movaps %xmm5,nb101_dzO(%rsp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb101_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb101_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## 11100000      
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps nb101_qqO(%rsp),%xmm3

        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        mulps  %xmm3,%xmm4      ## xmm4=total fscal 
        addps  nb101_vctot(%rsp),%xmm3

        movaps nb101_dxO(%rsp),%xmm0
        movaps nb101_dyO(%rsp),%xmm1
        movaps nb101_dzO(%rsp),%xmm2

        movaps %xmm3,nb101_vctot(%rsp)

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        movss  nb101_fixO(%rsp),%xmm3
        movss  nb101_fiyO(%rsp),%xmm4
        movss  nb101_fizO(%rsp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb101_fixO(%rsp)
        movss  %xmm4,nb101_fiyO(%rsp)
        movss  %xmm5,nb101_fizO(%rsp)   ## updated the O force now do the H's 
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $230,%xmm3,%xmm3 ## 11100110      ;# shift right 
        shufps $230,%xmm4,%xmm4 ## 11100110
        shufps $230,%xmm5,%xmm5 ## 11100110
        addss  nb101_fixH1(%rsp),%xmm3
        addss  nb101_fiyH1(%rsp),%xmm4
        addss  nb101_fizH1(%rsp),%xmm5
        movss  %xmm3,nb101_fixH1(%rsp)
        movss  %xmm4,nb101_fiyH1(%rsp)
        movss  %xmm5,nb101_fizH1(%rsp)          ## updated the H1 force 

        movq nb101_faction(%rbp),%rdi
        shufps $231,%xmm3,%xmm3 ## 11100111      ;# shift right 
        shufps $231,%xmm4,%xmm4 ## 11100111
        shufps $231,%xmm5,%xmm5 ## 11100111
        addss  nb101_fixH2(%rsp),%xmm3
        addss  nb101_fiyH2(%rsp),%xmm4
        addss  nb101_fizH2(%rsp),%xmm5
        movss  %xmm3,nb101_fixH2(%rsp)
        movss  %xmm4,nb101_fiyH2(%rsp)
        movss  %xmm5,nb101_fizH2(%rsp)          ## updated the H2 force 

        ## the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1 
        xorps  %xmm5,%xmm5
        movaps %xmm0,%xmm3
        movlps (%rdi,%rax,4),%xmm6
        movss  8(%rdi,%rax,4),%xmm7
        unpcklps %xmm1,%xmm3
        movlhps  %xmm5,%xmm3
        unpckhps %xmm1,%xmm0
        addps    %xmm3,%xmm0
        movhlps  %xmm0,%xmm3
        addps    %xmm3,%xmm0    ## x,y sum in xmm0 

        movhlps  %xmm2,%xmm1
        addss    %xmm1,%xmm2
        shufps  $1,%xmm1,%xmm1
        addss    %xmm1,%xmm2   ## z sum in xmm2 
        addps    %xmm0,%xmm6
        addss    %xmm2,%xmm7

        movlps %xmm6,(%rdi,%rax,4)
        movss  %xmm7,8(%rdi,%rax,4)

        decl  nb101_innerk(%rsp)
        jz    _nb_kernel101_x86_64_sse.nb101_updateouterdata
        jmp   _nb_kernel101_x86_64_sse.nb101_odd_loop
_nb_kernel101_x86_64_sse.nb101_updateouterdata: 
        movl  nb101_ii3(%rsp),%ecx
        movq  nb101_faction(%rbp),%rdi
        movq  nb101_fshift(%rbp),%rsi
        movl  nb101_is3(%rsp),%edx

        ## accumulate Oi forces in xmm0, xmm1, xmm2 
        movaps nb101_fixO(%rsp),%xmm0
        movaps nb101_fiyO(%rsp),%xmm1
        movaps nb101_fizO(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  (%rdi,%rcx,4),%xmm3
        movss  4(%rdi,%rcx,4),%xmm4
        movss  8(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,(%rdi,%rcx,4)
        movss  %xmm4,4(%rdi,%rcx,4)
        movss  %xmm5,8(%rdi,%rcx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        movaps %xmm0,%xmm6
        movss %xmm2,%xmm7
        movlhps %xmm1,%xmm6
        shufps $8,%xmm6,%xmm6 ## 00001000       

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movaps nb101_fixH1(%rsp),%xmm0
        movaps nb101_fiyH1(%rsp),%xmm1
        movaps nb101_fizH1(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  12(%rdi,%rcx,4),%xmm3
        movss  16(%rdi,%rcx,4),%xmm4
        movss  20(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,12(%rdi,%rcx,4)
        movss  %xmm4,16(%rdi,%rcx,4)
        movss  %xmm5,20(%rdi,%rcx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## 00001000       
        addps   %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movaps nb101_fixH2(%rsp),%xmm0
        movaps nb101_fiyH2(%rsp),%xmm1
        movaps nb101_fizH2(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  24(%rdi,%rcx,4),%xmm3
        movss  28(%rdi,%rcx,4),%xmm4
        movss  32(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,24(%rdi,%rcx,4)
        movss  %xmm4,28(%rdi,%rcx,4)
        movss  %xmm5,32(%rdi,%rcx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## 00001000       
        addps   %xmm0,%xmm6

        ## increment fshift force  
        movlps  (%rsi,%rdx,4),%xmm3
        movss  8(%rsi,%rdx,4),%xmm4
        subps  %xmm6,%xmm3
        subss  %xmm7,%xmm4
        movlps  %xmm3,(%rsi,%rdx,4)
        movss  %xmm4,8(%rsi,%rdx,4)

        ## get n from stack
        movl nb101_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb101_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb101_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb101_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb101_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel101_x86_64_sse.nb101_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb101_n(%rsp)
        jmp _nb_kernel101_x86_64_sse.nb101_outer
_nb_kernel101_x86_64_sse.nb101_outerend: 
        ## check if more outer neighborlists remain
        movl  nb101_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel101_x86_64_sse.nb101_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel101_x86_64_sse.nb101_threadloop
_nb_kernel101_x86_64_sse.nb101_end: 

        movl nb101_nouter(%rsp),%eax
        movl nb101_ninner(%rsp),%ebx
        movq nb101_outeriter(%rbp),%rcx
        movq nb101_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $696,%rsp
        emms

        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret


.globl nb_kernel101nf_x86_64_sse
.globl _nb_kernel101nf_x86_64_sse
nb_kernel101nf_x86_64_sse:      
_nb_kernel101nf_x86_64_sse:     
.set nb101nf_fshift, 16
.set nb101nf_gid, 24
.set nb101nf_pos, 32
.set nb101nf_faction, 40
.set nb101nf_charge, 48
.set nb101nf_p_facel, 56
.set nb101nf_argkrf, 64
.set nb101nf_argcrf, 72
.set nb101nf_Vc, 80
.set nb101nf_type, 88
.set nb101nf_p_ntype, 96
.set nb101nf_vdwparam, 104
.set nb101nf_Vvdw, 112
.set nb101nf_p_tabscale, 120
.set nb101nf_VFtab, 128
.set nb101nf_invsqrta, 136
.set nb101nf_dvda, 144
.set nb101nf_p_gbtabscale, 152
.set nb101nf_GBtab, 160
.set nb101nf_p_nthreads, 168
.set nb101nf_count, 176
.set nb101nf_mtx, 184
.set nb101nf_outeriter, 192
.set nb101nf_inneriter, 200
.set nb101nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb101nf_ixO, 0
.set nb101nf_iyO, 16
.set nb101nf_izO, 32
.set nb101nf_ixH1, 48
.set nb101nf_iyH1, 64
.set nb101nf_izH1, 80
.set nb101nf_ixH2, 96
.set nb101nf_iyH2, 112
.set nb101nf_izH2, 128
.set nb101nf_iqO, 144
.set nb101nf_iqH, 160
.set nb101nf_qqO, 176
.set nb101nf_qqH, 192
.set nb101nf_vctot, 208
.set nb101nf_half, 224
.set nb101nf_three, 240
.set nb101nf_is3, 256
.set nb101nf_ii3, 260
.set nb101nf_nri, 264
.set nb101nf_iinr, 272
.set nb101nf_jindex, 280
.set nb101nf_jjnr, 288
.set nb101nf_shift, 296
.set nb101nf_shiftvec, 304
.set nb101nf_facel, 312
.set nb101nf_innerjjnr, 320
.set nb101nf_innerk, 328
.set nb101nf_n, 332
.set nb101nf_nn1, 336
.set nb101nf_nouter, 340
.set nb101nf_ninner, 344

        push %rbp
        movq %rsp,%rbp
        push %rbx

        subq $360,%rsp
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb101nf_nouter(%rsp)
        movl %eax,nb101nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb101nf_nri(%rsp)
        movq %rsi,nb101nf_iinr(%rsp)
        movq %rdx,nb101nf_jindex(%rsp)
        movq %rcx,nb101nf_jjnr(%rsp)
        movq %r8,nb101nf_shift(%rsp)
        movq %r9,nb101nf_shiftvec(%rsp)
        movq nb101nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb101nf_facel(%rsp)


        ## assume we have at least one i particle - start directly 
        movq  nb101nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb101nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss 4(%rdx,%rbx,4),%xmm4
        movq nb101nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb101nf_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb101nf_iqO(%rsp)
        movaps %xmm4,nb101nf_iqH(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb101nf_half(%rsp)
        movss nb101nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb101nf_half(%rsp)
        movaps %xmm3,nb101nf_three(%rsp)

_nb_kernel101nf_x86_64_sse.nb101nf_threadloop: 
        movq  nb101nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel101nf_x86_64_sse.nb101nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel101nf_x86_64_sse.nb101nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb101nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb101nf_n(%rsp)
        movl %ebx,nb101nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel101nf_x86_64_sse.nb101nf_outerstart
        jmp _nb_kernel101nf_x86_64_sse.nb101nf_end

_nb_kernel101nf_x86_64_sse.nb101nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb101nf_nouter(%rsp),%ebx
        movl %ebx,nb101nf_nouter(%rsp)

_nb_kernel101nf_x86_64_sse.nb101nf_outer: 
        movq  nb101nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb101nf_is3(%rsp)            ## store is3 

        movq  nb101nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb101nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb101nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb101nf_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb101nf_ixO(%rsp)
        movaps %xmm4,nb101nf_iyO(%rsp)
        movaps %xmm5,nb101nf_izO(%rsp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%rax,%rbx,4),%xmm0
        addss 16(%rax,%rbx,4),%xmm1
        addss 20(%rax,%rbx,4),%xmm2
        addss 24(%rax,%rbx,4),%xmm3
        addss 28(%rax,%rbx,4),%xmm4
        addss 32(%rax,%rbx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb101nf_ixH1(%rsp)
        movaps %xmm1,nb101nf_iyH1(%rsp)
        movaps %xmm2,nb101nf_izH1(%rsp)
        movaps %xmm3,nb101nf_ixH2(%rsp)
        movaps %xmm4,nb101nf_iyH2(%rsp)
        movaps %xmm5,nb101nf_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb101nf_vctot(%rsp)

        movq  nb101nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb101nf_pos(%rbp),%rsi
        movq  nb101nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb101nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb101nf_ninner(%rsp),%ecx
        movl  %ecx,nb101nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb101nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel101nf_x86_64_sse.nb101nf_unroll_loop
        jmp   _nb_kernel101nf_x86_64_sse.nb101nf_odd_inner
_nb_kernel101nf_x86_64_sse.nb101nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb101nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb101nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb101nf_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb101nf_iqO(%rsp),%xmm3
        mulps  nb101nf_iqH(%rsp),%xmm4

        movaps  %xmm3,nb101nf_qqO(%rsp)
        movaps  %xmm4,nb101nf_qqH(%rsp)

        movq nb101nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx     ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

        ## move four coordinates to xmm0-xmm2   
        movlps (%rsi,%rax,4),%xmm4
        movlps (%rsi,%rcx,4),%xmm5
        movss 8(%rsi,%rax,4),%xmm2
        movss 8(%rsi,%rcx,4),%xmm6

        movhps (%rsi,%rbx,4),%xmm4
        movhps (%rsi,%rdx,4),%xmm5

        movss 8(%rsi,%rbx,4),%xmm0
        movss 8(%rsi,%rdx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## 10001000

        shufps $136,%xmm5,%xmm0 ## 10001000
        shufps $221,%xmm5,%xmm1 ## 11011101             

        ## move ixO-izO to xmm4-xmm6 
        movaps nb101nf_ixO(%rsp),%xmm4
        movaps nb101nf_iyO(%rsp),%xmm5
        movaps nb101nf_izO(%rsp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        movaps %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb101nf_ixH1(%rsp),%xmm4
        movaps nb101nf_iyH1(%rsp),%xmm5
        movaps nb101nf_izH1(%rsp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movaps nb101nf_ixH2(%rsp),%xmm3
        movaps nb101nf_iyH2(%rsp),%xmm4
        movaps nb101nf_izH2(%rsp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## square it 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm4,%xmm5
        addps %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqO - seed in xmm2       
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb101nf_three(%rsp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb101nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm7     ## rinvO in xmm7 
        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb101nf_three(%rsp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb101nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm6     ## rinvH1 in xmm6 
        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb101nf_three(%rsp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb101nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        mulps  nb101nf_qqO(%rsp),%xmm7          ## xmm7=vcoul 
        addps  nb101nf_vctot(%rsp),%xmm7
        movaps %xmm7,nb101nf_vctot(%rsp)

        ## H1 interactions 
        mulps  nb101nf_qqH(%rsp),%xmm6          ## xmm6=vcoul 
        addps  nb101nf_vctot(%rsp),%xmm6
        movaps %xmm6,nb101nf_vctot(%rsp)

        ## H2 interactions 
        mulps  nb101nf_qqH(%rsp),%xmm5          ## xmm5=vcoul 
        addps  nb101nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb101nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb101nf_innerk(%rsp)
        jl    _nb_kernel101nf_x86_64_sse.nb101nf_odd_inner
        jmp   _nb_kernel101nf_x86_64_sse.nb101nf_unroll_loop
_nb_kernel101nf_x86_64_sse.nb101nf_odd_inner: 
        addl $4,nb101nf_innerk(%rsp)
        jnz   _nb_kernel101nf_x86_64_sse.nb101nf_odd_loop
        jmp   _nb_kernel101nf_x86_64_sse.nb101nf_updateouterdata
_nb_kernel101nf_x86_64_sse.nb101nf_odd_loop: 
        movq  nb101nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb101nf_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb101nf_iqO(%rsp),%xmm4
        movq nb101nf_charge(%rbp),%rsi
        movhps nb101nf_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb101nf_qqO(%rsp)          ## use oxygen qq for storage 

        movq nb101nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm0
        movss 4(%rsi,%rax,4),%xmm1
        movss 8(%rsi,%rax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb101nf_ixO(%rsp),%xmm3
        movss nb101nf_iyO(%rsp),%xmm4
        movss nb101nf_izO(%rsp),%xmm5

        movlps nb101nf_ixH1(%rsp),%xmm6
        movlps nb101nf_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb101nf_iyH1(%rsp),%xmm6
        movlps nb101nf_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb101nf_izH1(%rsp),%xmm6
        movlps nb101nf_izH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm5

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb101nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb101nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## 11100000      

        movaps nb101nf_qqO(%rsp),%xmm3

        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        addps  nb101nf_vctot(%rsp),%xmm3
        movaps %xmm3,nb101nf_vctot(%rsp)

        decl  nb101nf_innerk(%rsp)
        jz    _nb_kernel101nf_x86_64_sse.nb101nf_updateouterdata
        jmp   _nb_kernel101nf_x86_64_sse.nb101nf_odd_loop
_nb_kernel101nf_x86_64_sse.nb101nf_updateouterdata: 
        ## accumulate total potential energy and update it 
        ## get n from stack
        movl nb101nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb101nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        movaps nb101nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb101nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb101nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel101nf_x86_64_sse.nb101nf_outerend

        ## not last, iterate outer loop once more!
        movl %esi,nb101nf_n(%rsp)
        jmp _nb_kernel101nf_x86_64_sse.nb101nf_outer
_nb_kernel101nf_x86_64_sse.nb101nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb101nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel101nf_x86_64_sse.nb101nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel101nf_x86_64_sse.nb101nf_threadloop
_nb_kernel101nf_x86_64_sse.nb101nf_end: 

        movl nb101nf_nouter(%rsp),%eax
        movl nb101nf_ninner(%rsp),%ebx
        movq nb101nf_outeriter(%rbp),%rcx
        movq nb101nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $360,%rsp
        emms

        pop %rbx
        pop    %rbp
        ret

