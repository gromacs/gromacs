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





.globl nb_kernel111_x86_64_sse
.globl _nb_kernel111_x86_64_sse
nb_kernel111_x86_64_sse:        
_nb_kernel111_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb111_fshift, 16
.set nb111_gid, 24
.set nb111_pos, 32
.set nb111_faction, 40
.set nb111_charge, 48
.set nb111_p_facel, 56
.set nb111_argkrf, 64
.set nb111_argcrf, 72
.set nb111_Vc, 80
.set nb111_type, 88
.set nb111_p_ntype, 96
.set nb111_vdwparam, 104
.set nb111_Vvdw, 112
.set nb111_p_tabscale, 120
.set nb111_VFtab, 128
.set nb111_invsqrta, 136
.set nb111_dvda, 144
.set nb111_p_gbtabscale, 152
.set nb111_GBtab, 160
.set nb111_p_nthreads, 168
.set nb111_count, 176
.set nb111_mtx, 184
.set nb111_outeriter, 192
.set nb111_inneriter, 200
.set nb111_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb111_ixO, 0
.set nb111_iyO, 16
.set nb111_izO, 32
.set nb111_ixH1, 48
.set nb111_iyH1, 64
.set nb111_izH1, 80
.set nb111_ixH2, 96
.set nb111_iyH2, 112
.set nb111_izH2, 128
.set nb111_iqO, 144
.set nb111_iqH, 160
.set nb111_dxO, 176
.set nb111_dyO, 192
.set nb111_dzO, 208
.set nb111_dxH1, 224
.set nb111_dyH1, 240
.set nb111_dzH1, 256
.set nb111_dxH2, 272
.set nb111_dyH2, 288
.set nb111_dzH2, 304
.set nb111_qqO, 320
.set nb111_qqH, 336
.set nb111_c6, 352
.set nb111_c12, 368
.set nb111_six, 384
.set nb111_twelve, 400
.set nb111_vctot, 416
.set nb111_Vvdwtot, 432
.set nb111_fixO, 448
.set nb111_fiyO, 464
.set nb111_fizO, 480
.set nb111_fixH1, 496
.set nb111_fiyH1, 512
.set nb111_fizH1, 528
.set nb111_fixH2, 544
.set nb111_fiyH2, 560
.set nb111_fizH2, 576
.set nb111_fjx, 592
.set nb111_fjy, 608
.set nb111_fjz, 624
.set nb111_half, 640
.set nb111_three, 656
.set nb111_is3, 672
.set nb111_ii3, 676
.set nb111_nri, 692
.set nb111_iinr, 700
.set nb111_jindex, 708
.set nb111_jjnr, 716
.set nb111_shift, 724
.set nb111_shiftvec, 732
.set nb111_facel, 740
.set nb111_innerjjnr, 748
.set nb111_ntia, 756
.set nb111_innerk, 760
.set nb111_n, 764
.set nb111_nn1, 768
.set nb111_nouter, 772
.set nb111_ninner, 776


        push %rbp
        movq %rsp,%rbp
        push %rbx

        push %r12
        push %r13
        push %r14
        push %r15

        subq $792,%rsp
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb111_nouter(%rsp)
        movl %eax,nb111_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb111_nri(%rsp)
        movq %rsi,nb111_iinr(%rsp)
        movq %rdx,nb111_jindex(%rsp)
        movq %rcx,nb111_jjnr(%rsp)
        movq %r8,nb111_shift(%rsp)
        movq %r9,nb111_shiftvec(%rsp)
        movq nb111_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb111_facel(%rsp)


        ## assume we have at least one i particle - start directly 
        movq  nb111_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb111_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss 4(%rdx,%rbx,4),%xmm4
        movq nb111_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb111_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb111_iqO(%rsp)
        movaps %xmm4,nb111_iqH(%rsp)

        movq  nb111_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb111_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb111_ntia(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb111_half(%rsp)
        movss nb111_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm3,%xmm4
        addps  %xmm4,%xmm4      ## six
        movaps %xmm4,%xmm5
        addps  %xmm5,%xmm5      ## twelve
        movaps %xmm1,nb111_half(%rsp)
        movaps %xmm3,nb111_three(%rsp)
        movaps %xmm4,nb111_six(%rsp)
        movaps %xmm5,nb111_twelve(%rsp)

_nb_kernel111_x86_64_sse.nb111_threadloop: 
        movq  nb111_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel111_x86_64_sse.nb111_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel111_x86_64_sse.nb111_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb111_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb111_n(%rsp)
        movl %ebx,nb111_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel111_x86_64_sse.nb111_outerstart
        jmp _nb_kernel111_x86_64_sse.nb111_end

_nb_kernel111_x86_64_sse.nb111_outerstart: 
        ## ebx contains number of outer iterations
        addl nb111_nouter(%rsp),%ebx
        movl %ebx,nb111_nouter(%rsp)

_nb_kernel111_x86_64_sse.nb111_outer: 
        movq  nb111_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb111_is3(%rsp)      ## store is3 

        movq  nb111_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb111_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb111_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb111_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb111_ixO(%rsp)
        movaps %xmm4,nb111_iyO(%rsp)
        movaps %xmm5,nb111_izO(%rsp)

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
        movaps %xmm0,nb111_ixH1(%rsp)
        movaps %xmm1,nb111_iyH1(%rsp)
        movaps %xmm2,nb111_izH1(%rsp)
        movaps %xmm3,nb111_ixH2(%rsp)
        movaps %xmm4,nb111_iyH2(%rsp)
        movaps %xmm5,nb111_izH2(%rsp)

        ## clear potentials and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb111_vctot(%rsp)
        movaps %xmm4,nb111_Vvdwtot(%rsp)
        movaps %xmm4,nb111_fixO(%rsp)
        movaps %xmm4,nb111_fiyO(%rsp)
        movaps %xmm4,nb111_fizO(%rsp)
        movaps %xmm4,nb111_fixH1(%rsp)
        movaps %xmm4,nb111_fiyH1(%rsp)
        movaps %xmm4,nb111_fizH1(%rsp)
        movaps %xmm4,nb111_fixH2(%rsp)
        movaps %xmm4,nb111_fiyH2(%rsp)
        movaps %xmm4,nb111_fizH2(%rsp)

        movq  nb111_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb111_pos(%rbp),%rsi
        movq  nb111_faction(%rbp),%rdi
        movq  nb111_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb111_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb111_ninner(%rsp),%ecx
        movl  %ecx,nb111_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb111_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel111_x86_64_sse.nb111_unroll_loop
        jmp   _nb_kernel111_x86_64_sse.nb111_odd_inner
_nb_kernel111_x86_64_sse.nb111_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb111_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%r8d
        movl  4(%rdx),%r9d
        movl  8(%rdx),%r10d
        movl  12(%rdx),%r11d           ## eax-edx=jnr1-4 

        addq $16,nb111_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb111_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%r8,%r8,2),%rax     ## replace jnr with j3 
        lea  (%r9,%r9,2),%rbx
        lea  (%r10,%r10,2),%rcx     ## replace jnr with j3 
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
        movq nb111_charge(%rbp),%rsi     ## base of charge[] 

    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8

    subps nb111_ixO(%rsp),%xmm0
    subps nb111_iyO(%rsp),%xmm1
    subps nb111_izO(%rsp),%xmm2
    subps nb111_ixH1(%rsp),%xmm3
    subps nb111_iyH1(%rsp),%xmm4
    subps nb111_izH1(%rsp),%xmm5
    subps nb111_ixH2(%rsp),%xmm6
    subps nb111_iyH2(%rsp),%xmm7
    subps nb111_izH2(%rsp),%xmm8

        movaps %xmm0,nb111_dxO(%rsp)
        movaps %xmm1,nb111_dyO(%rsp)
        movaps %xmm2,nb111_dzO(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb111_dxH1(%rsp)
        movaps %xmm4,nb111_dyH1(%rsp)
        movaps %xmm5,nb111_dzH1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb111_dxH2(%rsp)
        movaps %xmm7,nb111_dyH2(%rsp)
        movaps %xmm8,nb111_dzH2(%rsp)
        mulps  %xmm6,%xmm6
        mulps  %xmm7,%xmm7
        mulps  %xmm8,%xmm8
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
    addps  %xmm7,%xmm6
    addps  %xmm8,%xmm6

        movss (%rsi,%r8,4),%xmm9
        movss (%rsi,%r10,4),%xmm10
        movss (%rsi,%r9,4),%xmm11
        movss (%rsi,%r11,4),%xmm12

    unpcklps %xmm10,%xmm9
    unpcklps %xmm12,%xmm11
    unpcklps %xmm11,%xmm9
        movaps %xmm9,%xmm7
        mulps  nb111_iqO(%rsp),%xmm9
        mulps  nb111_iqH(%rsp),%xmm7

        movaps  %xmm9,nb111_qqO(%rsp)
        movaps  %xmm7,nb111_qqH(%rsp)

        movq nb111_type(%rbp),%rsi
        movl (%rsi,%r8,4),%r8d
        movl (%rsi,%r9,4),%r9d
        movl (%rsi,%r10,4),%r10d
        movl (%rsi,%r11,4),%r11d

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

        movq nb111_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        shll %r10d
        shll %r11d

        movaps  nb111_three(%rsp),%xmm9
        movaps  %xmm9,%xmm10
    movaps  %xmm9,%xmm11

        mulps   %xmm0,%xmm1 ## rsq*lu*lu
        mulps   %xmm3,%xmm4 ## rsq*lu*lu 
    mulps   %xmm6,%xmm7 ## rsq*lu*lu

        subps   %xmm1,%xmm9
        subps   %xmm4,%xmm10
    subps   %xmm7,%xmm11 ## 3-rsq*lu*lu

        movl nb111_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d
        addl %edi,%r10d
        addl %edi,%r11d

        mulps   %xmm2,%xmm9
        mulps   %xmm5,%xmm10
    mulps   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movaps  nb111_half(%rsp),%xmm0
        mulps   %xmm0,%xmm9 ## rinvO 
        mulps   %xmm0,%xmm10 ## rinvH1
    mulps   %xmm0,%xmm11 ## rinvH2

    movlps (%rsi,%r8,4),%xmm3
        movlps (%rsi,%r10,4),%xmm2
        movhps (%rsi,%r9,4),%xmm3
        movhps (%rsi,%r11,4),%xmm2

        movaps %xmm3,%xmm5
        shufps $136,%xmm2,%xmm3 ## 10001000
        shufps $221,%xmm2,%xmm5 ## 11011101

        ## interactions 
    movaps %xmm9,%xmm0
    movaps %xmm10,%xmm1
    movaps %xmm11,%xmm2
    mulps  %xmm9,%xmm9   ## rinvsq
    mulps  %xmm10,%xmm10
    mulps  %xmm11,%xmm11
    movaps %xmm9,%xmm12
    mulps  %xmm12,%xmm12 ## rinv4
    mulps  %xmm9,%xmm12 ## rinv6
    mulps  nb111_qqO(%rsp),%xmm0
    mulps  nb111_qqH(%rsp),%xmm1
    mulps  nb111_qqH(%rsp),%xmm2
    movaps %xmm12,%xmm13 ## rinv6
    mulps %xmm12,%xmm12 ## rinv12
        mulps  %xmm3,%xmm13
        mulps  %xmm5,%xmm12
    movaps %xmm12,%xmm14
    subps  %xmm13,%xmm14

        addps  nb111_Vvdwtot(%rsp),%xmm14
        mulps  nb111_six(%rsp),%xmm13
        mulps  nb111_twelve(%rsp),%xmm12
        movaps %xmm14,nb111_Vvdwtot(%rsp)
    subps  %xmm13,%xmm12 ## LJ fscal        

    addps  %xmm0,%xmm12

    mulps  %xmm12,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps nb111_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb111_vctot(%rsp)

        movq  nb111_faction(%rbp),%rdi
        ## move j  forces to local temp variables 
    movlps (%rdi,%rax,4),%xmm0 ## jxa jya  -   -
    movlps (%rdi,%rcx,4),%xmm1 ## jxc jyc  -   -
    movhps (%rdi,%rbx,4),%xmm0 ## jxa jya jxb jyb 
    movhps (%rdi,%rdx,4),%xmm1 ## jxc jyc jxd jyd 

    movss  8(%rdi,%rax,4),%xmm2    ## jza  -  -  -
    movss  8(%rdi,%rcx,4),%xmm3    ## jzc  -  -  -
    movss  8(%rdi,%rbx,4),%xmm4    ## jzb
    movss  8(%rdi,%rdx,4),%xmm5    ## jzd
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

        mulps nb111_dxO(%rsp),%xmm7
        mulps nb111_dyO(%rsp),%xmm8
        mulps nb111_dzO(%rsp),%xmm9
        mulps nb111_dxH1(%rsp),%xmm10
        mulps nb111_dyH1(%rsp),%xmm11
        mulps nb111_dzH1(%rsp),%xmm12
        mulps nb111_dxH2(%rsp),%xmm13
        mulps nb111_dyH2(%rsp),%xmm14
        mulps nb111_dzH2(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb111_fixO(%rsp),%xmm7
    addps nb111_fiyO(%rsp),%xmm8
    addps nb111_fizO(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb111_fixH1(%rsp),%xmm10
    addps nb111_fiyH1(%rsp),%xmm11
    addps nb111_fizH1(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb111_fixH2(%rsp),%xmm13
    addps nb111_fiyH2(%rsp),%xmm14
    addps nb111_fizH2(%rsp),%xmm15

    movaps %xmm7,nb111_fixO(%rsp)
    movaps %xmm8,nb111_fiyO(%rsp)
    movaps %xmm9,nb111_fizO(%rsp)
    movaps %xmm10,nb111_fixH1(%rsp)
    movaps %xmm11,nb111_fiyH1(%rsp)
    movaps %xmm12,nb111_fizH1(%rsp)
    movaps %xmm13,nb111_fixH2(%rsp)
    movaps %xmm14,nb111_fiyH2(%rsp)
    movaps %xmm15,nb111_fizH2(%rsp)

    ## xmm3 = fjx , xmm4 = fjy, xmm5=fjz
    movaps %xmm3,%xmm5
    unpcklps %xmm4,%xmm3
    unpckhps %xmm4,%xmm5

    addps %xmm3,%xmm0
    addps %xmm5,%xmm1

    movhlps  %xmm2,%xmm3

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
        subl $4,nb111_innerk(%rsp)
        jl    _nb_kernel111_x86_64_sse.nb111_odd_inner
        jmp   _nb_kernel111_x86_64_sse.nb111_unroll_loop
_nb_kernel111_x86_64_sse.nb111_odd_inner: 
        addl $4,nb111_innerk(%rsp)
        jnz   _nb_kernel111_x86_64_sse.nb111_odd_loop
        jmp   _nb_kernel111_x86_64_sse.nb111_updateouterdata
_nb_kernel111_x86_64_sse.nb111_odd_loop: 
        movq  nb111_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb111_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb111_iqO(%rsp),%xmm4
        movq nb111_charge(%rbp),%rsi
        movhps nb111_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb111_qqO(%rsp)    ## use oxygen qq for storage 

        xorps %xmm6,%xmm6
        movq nb111_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ebx
        movq nb111_vdwparam(%rbp),%rsi
        shll %ebx
        addl nb111_ntia(%rsp),%ebx
        movlps (%rsi,%rbx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## 11111100
        shufps $253,%xmm7,%xmm7 ## 11111101
        movaps %xmm6,nb111_c6(%rsp)
        movaps %xmm7,nb111_c12(%rsp)

        movq nb111_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm3
        movss 4(%rsi,%rax,4),%xmm4
        movss 8(%rsi,%rax,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5

        movss nb111_ixO(%rsp),%xmm0
        movss nb111_iyO(%rsp),%xmm1
        movss nb111_izO(%rsp),%xmm2

        movlps nb111_ixH1(%rsp),%xmm6
        movlps nb111_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm0
        movlps nb111_iyH1(%rsp),%xmm6
        movlps nb111_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm1
        movlps nb111_izH1(%rsp),%xmm6
        movlps nb111_izH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm2

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        movaps %xmm3,nb111_dxO(%rsp)
        movaps %xmm4,nb111_dyO(%rsp)
        movaps %xmm5,nb111_dzO(%rsp)

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
        movaps nb111_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb111_half(%rsp),%xmm0
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
        movaps %xmm4,%xmm1
        mulss  %xmm4,%xmm1
        movaps nb111_qqO(%rsp),%xmm3
        mulss  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        mulps  nb111_c6(%rsp),%xmm1
        mulps  nb111_c12(%rsp),%xmm2
        movaps %xmm2,%xmm5
        subss  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb111_Vvdwtot(%rsp),%xmm5
        mulss  nb111_six(%rsp),%xmm1
        mulss  nb111_twelve(%rsp),%xmm2
        subss  %xmm1,%xmm2
        addps  %xmm3,%xmm2
        mulps  %xmm2,%xmm4      ## xmm4=total fscal 
        addps  nb111_vctot(%rsp),%xmm3

        movaps nb111_dxO(%rsp),%xmm0
        movaps nb111_dyO(%rsp),%xmm1
        movaps nb111_dzO(%rsp),%xmm2

        movaps %xmm3,nb111_vctot(%rsp)
        movaps %xmm5,nb111_Vvdwtot(%rsp)

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        movss  nb111_fixO(%rsp),%xmm3
        movss  nb111_fiyO(%rsp),%xmm4
        movss  nb111_fizO(%rsp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb111_fixO(%rsp)
        movss  %xmm4,nb111_fiyO(%rsp)
        movss  %xmm5,nb111_fizO(%rsp)   ## updated the O force now do the H's 
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $230,%xmm3,%xmm3 ## 11100110      ;# shift right 
        shufps $230,%xmm4,%xmm4 ## 11100110
        shufps $230,%xmm5,%xmm5 ## 11100110
        addss  nb111_fixH1(%rsp),%xmm3
        addss  nb111_fiyH1(%rsp),%xmm4
        addss  nb111_fizH1(%rsp),%xmm5
        movss  %xmm3,nb111_fixH1(%rsp)
        movss  %xmm4,nb111_fiyH1(%rsp)
        movss  %xmm5,nb111_fizH1(%rsp)          ## updated the H1 force 

        movq nb111_faction(%rbp),%rdi
        shufps $231,%xmm3,%xmm3 ## 11100111      ;# shift right 
        shufps $231,%xmm4,%xmm4 ## 11100111
        shufps $231,%xmm5,%xmm5 ## 11100111
        addss  nb111_fixH2(%rsp),%xmm3
        addss  nb111_fiyH2(%rsp),%xmm4
        addss  nb111_fizH2(%rsp),%xmm5
        movss  %xmm3,nb111_fixH2(%rsp)
        movss  %xmm4,nb111_fiyH2(%rsp)
        movss  %xmm5,nb111_fizH2(%rsp)          ## updated the H2 force 

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

        decl nb111_innerk(%rsp)
        jz    _nb_kernel111_x86_64_sse.nb111_updateouterdata
        jmp   _nb_kernel111_x86_64_sse.nb111_odd_loop
_nb_kernel111_x86_64_sse.nb111_updateouterdata: 
        movl  nb111_ii3(%rsp),%ecx
        movq  nb111_faction(%rbp),%rdi
        movq  nb111_fshift(%rbp),%rsi
        movl  nb111_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb111_fixO(%rsp),%xmm0
        movaps nb111_fiyO(%rsp),%xmm1
        movaps nb111_fizO(%rsp),%xmm2

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
        movaps nb111_fixH1(%rsp),%xmm0
        movaps nb111_fiyH1(%rsp),%xmm1
        movaps nb111_fizH1(%rsp),%xmm2

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
        movaps nb111_fixH2(%rsp),%xmm0
        movaps nb111_fiyH2(%rsp),%xmm1
        movaps nb111_fizH2(%rsp),%xmm2

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
        movl nb111_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb111_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb111_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb111_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb111_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb111_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb111_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel111_x86_64_sse.nb111_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb111_n(%rsp)
        jmp _nb_kernel111_x86_64_sse.nb111_outer
_nb_kernel111_x86_64_sse.nb111_outerend: 
        ## check if more outer neighborlists remain
        movl  nb111_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel111_x86_64_sse.nb111_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel111_x86_64_sse.nb111_threadloop
_nb_kernel111_x86_64_sse.nb111_end: 


        movl nb111_nouter(%rsp),%eax
        movl nb111_ninner(%rsp),%ebx
        movq nb111_outeriter(%rbp),%rcx
        movq nb111_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $792,%rsp
        emms

        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret





.globl nb_kernel111nf_x86_64_sse
.globl _nb_kernel111nf_x86_64_sse
nb_kernel111nf_x86_64_sse:      
_nb_kernel111nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb111nf_fshift, 16
.set nb111nf_gid, 24
.set nb111nf_pos, 32
.set nb111nf_faction, 40
.set nb111nf_charge, 48
.set nb111nf_p_facel, 56
.set nb111nf_argkrf, 64
.set nb111nf_argcrf, 72
.set nb111nf_Vc, 80
.set nb111nf_type, 88
.set nb111nf_p_ntype, 96
.set nb111nf_vdwparam, 104
.set nb111nf_Vvdw, 112
.set nb111nf_p_tabscale, 120
.set nb111nf_VFtab, 128
.set nb111nf_invsqrta, 136
.set nb111nf_dvda, 144
.set nb111nf_p_gbtabscale, 152
.set nb111nf_GBtab, 160
.set nb111nf_p_nthreads, 168
.set nb111nf_count, 176
.set nb111nf_mtx, 184
.set nb111nf_outeriter, 192
.set nb111nf_inneriter, 200
.set nb111nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb111nf_ixO, 0
.set nb111nf_iyO, 16
.set nb111nf_izO, 32
.set nb111nf_ixH1, 48
.set nb111nf_iyH1, 64
.set nb111nf_izH1, 80
.set nb111nf_ixH2, 96
.set nb111nf_iyH2, 112
.set nb111nf_izH2, 128
.set nb111nf_iqO, 144
.set nb111nf_iqH, 160
.set nb111nf_qqO, 176
.set nb111nf_qqH, 192
.set nb111nf_c6, 208
.set nb111nf_c12, 224
.set nb111nf_vctot, 240
.set nb111nf_Vvdwtot, 256
.set nb111nf_half, 272
.set nb111nf_three, 288
.set nb111nf_is3, 304
.set nb111nf_ii3, 308
.set nb111nf_nri, 324
.set nb111nf_iinr, 332
.set nb111nf_jindex, 340
.set nb111nf_jjnr, 348
.set nb111nf_shift, 356
.set nb111nf_shiftvec, 364
.set nb111nf_facel, 372
.set nb111nf_innerjjnr, 380
.set nb111nf_ntia, 388
.set nb111nf_innerk, 392
.set nb111nf_n, 396
.set nb111nf_nn1, 400
.set nb111nf_nouter, 404
.set nb111nf_ninner, 408

        push %rbp
        movq %rsp,%rbp
        push %rbx

        subq $424,%rsp
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb111nf_nouter(%rsp)
        movl %eax,nb111nf_ninner(%rsp)


        movl (%rdi),%edi
        movl %edi,nb111nf_nri(%rsp)
        movq %rsi,nb111nf_iinr(%rsp)
        movq %rdx,nb111nf_jindex(%rsp)
        movq %rcx,nb111nf_jjnr(%rsp)
        movq %r8,nb111nf_shift(%rsp)
        movq %r9,nb111nf_shiftvec(%rsp)
        movq nb111nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb111nf_facel(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb111nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb111nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss 4(%rdx,%rbx,4),%xmm4
        movq nb111nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb111nf_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb111nf_iqO(%rsp)
        movaps %xmm4,nb111nf_iqH(%rsp)

        movq  nb111nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb111nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb111nf_ntia(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb111nf_half(%rsp)
        movss nb111nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb111nf_half(%rsp)
        movaps %xmm3,nb111nf_three(%rsp)

_nb_kernel111nf_x86_64_sse.nb111nf_threadloop: 
        movq  nb111nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel111nf_x86_64_sse.nb111nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel111nf_x86_64_sse.nb111nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb111nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb111nf_n(%rsp)
        movl %ebx,nb111nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel111nf_x86_64_sse.nb111nf_outerstart
        jmp _nb_kernel111nf_x86_64_sse.nb111nf_end

_nb_kernel111nf_x86_64_sse.nb111nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb111nf_nouter(%rsp),%ebx
        movl %ebx,nb111nf_nouter(%rsp)

_nb_kernel111nf_x86_64_sse.nb111nf_outer: 
        movq  nb111nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb111nf_is3(%rsp)            ## store is3 

        movq  nb111nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb111nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb111nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb111nf_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb111nf_ixO(%rsp)
        movaps %xmm4,nb111nf_iyO(%rsp)
        movaps %xmm5,nb111nf_izO(%rsp)

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
        movaps %xmm0,nb111nf_ixH1(%rsp)
        movaps %xmm1,nb111nf_iyH1(%rsp)
        movaps %xmm2,nb111nf_izH1(%rsp)
        movaps %xmm3,nb111nf_ixH2(%rsp)
        movaps %xmm4,nb111nf_iyH2(%rsp)
        movaps %xmm5,nb111nf_izH2(%rsp)

        ## clear vctot and Vvdwtot
        xorps %xmm4,%xmm4
        movaps %xmm4,nb111nf_vctot(%rsp)
        movaps %xmm4,nb111nf_Vvdwtot(%rsp)

        movq  nb111nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb111nf_pos(%rbp),%rsi
        movq  nb111nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb111nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb111nf_ninner(%rsp),%ecx
        movl  %ecx,nb111nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb111nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel111nf_x86_64_sse.nb111nf_unroll_loop
        jmp   _nb_kernel111nf_x86_64_sse.nb111nf_odd_inner
_nb_kernel111nf_x86_64_sse.nb111nf_unroll_loop: 

        ## quad-unroll innerloop here 
        movq  nb111nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb111nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb111nf_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb111nf_iqO(%rsp),%xmm3
        mulps  nb111nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb111nf_qqO(%rsp)
        movaps  %xmm4,nb111nf_qqH(%rsp)

        movq nb111nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb111nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb111nf_ntia(%rsp),%edi
        addl %edi,%eax
        addl %edi,%ebx
        addl %edi,%ecx
        addl %edi,%edx

        movlps (%rsi,%rax,4),%xmm6
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rbx,4),%xmm6
        movhps (%rsi,%rdx,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm6 ## 11011101

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        movaps %xmm4,nb111nf_c6(%rsp)
        movaps %xmm6,nb111nf_c12(%rsp)

        movq nb111nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movaps nb111nf_ixO(%rsp),%xmm4
        movaps nb111nf_iyO(%rsp),%xmm5
        movaps nb111nf_izO(%rsp),%xmm6

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
        movaps nb111nf_ixH1(%rsp),%xmm4
        movaps nb111nf_iyH1(%rsp),%xmm5
        movaps nb111nf_izH1(%rsp),%xmm6

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
        movaps nb111nf_ixH2(%rsp),%xmm3
        movaps nb111nf_iyH2(%rsp),%xmm4
        movaps nb111nf_izH2(%rsp),%xmm5

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
        movaps  nb111nf_three(%rsp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb111nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm7     ## rinvO in xmm7 
        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb111nf_three(%rsp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb111nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm6     ## rinvH1 in xmm6 
        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb111nf_three(%rsp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb111nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  nb111nf_qqO(%rsp),%xmm7          ## xmm7=vcoul 

        mulps  nb111nf_c6(%rsp),%xmm1
        mulps  nb111nf_c12(%rsp),%xmm2
        movaps %xmm2,%xmm3
        subps  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addps  nb111nf_Vvdwtot(%rsp),%xmm3
        addps  nb111nf_vctot(%rsp),%xmm7
        movaps %xmm3,nb111nf_Vvdwtot(%rsp)
        movaps %xmm7,nb111nf_vctot(%rsp)

        ## H1 & H2 interactions 
        addps  %xmm5,%xmm6          ## add H2 rinv 
        mulps  nb111nf_qqH(%rsp),%xmm6          ## xmm6=vcoul 
        addps  nb111nf_vctot(%rsp),%xmm6
        movaps %xmm6,nb111nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb111nf_innerk(%rsp)
        jl    _nb_kernel111nf_x86_64_sse.nb111nf_odd_inner
        jmp   _nb_kernel111nf_x86_64_sse.nb111nf_unroll_loop
_nb_kernel111nf_x86_64_sse.nb111nf_odd_inner: 
        addl $4,nb111nf_innerk(%rsp)
        jnz   _nb_kernel111nf_x86_64_sse.nb111nf_odd_loop
        jmp   _nb_kernel111nf_x86_64_sse.nb111nf_updateouterdata
_nb_kernel111nf_x86_64_sse.nb111nf_odd_loop: 
        movq  nb111nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb111nf_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb111nf_iqO(%rsp),%xmm4
        movq nb111nf_charge(%rbp),%rsi
        movhps nb111nf_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb111nf_qqO(%rsp)          ## use oxygen qq for storage 

        xorps %xmm6,%xmm6
        movq nb111nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ebx
        movq nb111nf_vdwparam(%rbp),%rsi
        shll %ebx
        addl nb111nf_ntia(%rsp),%ebx
        movlps (%rsi,%rbx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## 11111100
        shufps $253,%xmm7,%xmm7 ## 11111101
        movaps %xmm6,nb111nf_c6(%rsp)
        movaps %xmm7,nb111nf_c12(%rsp)

        movq nb111nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm0
        movss 4(%rsi,%rax,4),%xmm1
        movss 8(%rsi,%rax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb111nf_ixO(%rsp),%xmm3
        movss nb111nf_iyO(%rsp),%xmm4
        movss nb111nf_izO(%rsp),%xmm5

        movlps nb111nf_ixH1(%rsp),%xmm6
        movlps nb111nf_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb111nf_iyH1(%rsp),%xmm6
        movlps nb111nf_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb111nf_izH1(%rsp),%xmm6
        movlps nb111nf_izH2(%rsp),%xmm7
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
        movaps nb111nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb111nf_half(%rsp),%xmm0
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
        movaps %xmm4,%xmm1
        mulss  %xmm4,%xmm1
        movaps nb111nf_qqO(%rsp),%xmm3
        mulss  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        mulps  nb111nf_c6(%rsp),%xmm1
        mulps  nb111nf_c12(%rsp),%xmm2
        movaps %xmm2,%xmm5
        subss  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb111nf_Vvdwtot(%rsp),%xmm5
        addps  nb111nf_vctot(%rsp),%xmm3
        movaps %xmm3,nb111nf_vctot(%rsp)
        movaps %xmm5,nb111nf_Vvdwtot(%rsp)

        decl nb111nf_innerk(%rsp)
        jz    _nb_kernel111nf_x86_64_sse.nb111nf_updateouterdata
        jmp   _nb_kernel111nf_x86_64_sse.nb111nf_odd_loop
_nb_kernel111nf_x86_64_sse.nb111nf_updateouterdata: 
        ## get n from stack
        movl nb111nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb111nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        movaps nb111nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb111nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb111nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb111nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb111nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel111nf_x86_64_sse.nb111nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb111nf_n(%rsp)
        jmp _nb_kernel111nf_x86_64_sse.nb111nf_outer
_nb_kernel111nf_x86_64_sse.nb111nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb111nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel111nf_x86_64_sse.nb111nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel111nf_x86_64_sse.nb111nf_threadloop
_nb_kernel111nf_x86_64_sse.nb111nf_end: 


        movl nb111nf_nouter(%rsp),%eax
        movl nb111nf_ninner(%rsp),%ebx
        movq nb111nf_outeriter(%rbp),%rcx
        movq nb111nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $424,%rsp
        emms
        pop %rbx
        pop    %rbp
        ret



