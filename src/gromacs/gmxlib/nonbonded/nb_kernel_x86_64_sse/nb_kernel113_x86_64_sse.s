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





.globl nb_kernel113_x86_64_sse
.globl _nb_kernel113_x86_64_sse
nb_kernel113_x86_64_sse:        
_nb_kernel113_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb113_fshift, 16
.set nb113_gid, 24
.set nb113_pos, 32
.set nb113_faction, 40
.set nb113_charge, 48
.set nb113_p_facel, 56
.set nb113_argkrf, 64
.set nb113_argcrf, 72
.set nb113_Vc, 80
.set nb113_type, 88
.set nb113_p_ntype, 96
.set nb113_vdwparam, 104
.set nb113_Vvdw, 112
.set nb113_p_tabscale, 120
.set nb113_VFtab, 128
.set nb113_invsqrta, 136
.set nb113_dvda, 144
.set nb113_p_gbtabscale, 152
.set nb113_GBtab, 160
.set nb113_p_nthreads, 168
.set nb113_count, 176
.set nb113_mtx, 184
.set nb113_outeriter, 192
.set nb113_inneriter, 200
.set nb113_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb113_ixO, 0
.set nb113_iyO, 16
.set nb113_izO, 32
.set nb113_ixH1, 48
.set nb113_iyH1, 64
.set nb113_izH1, 80
.set nb113_ixH2, 96
.set nb113_iyH2, 112
.set nb113_izH2, 128
.set nb113_ixM, 144
.set nb113_iyM, 160
.set nb113_izM, 176
.set nb113_iqM, 192
.set nb113_iqH, 208
.set nb113_dxO, 224
.set nb113_dyO, 240
.set nb113_dzO, 256
.set nb113_dxH1, 272
.set nb113_dyH1, 288
.set nb113_dzH1, 304
.set nb113_dxH2, 320
.set nb113_dyH2, 336
.set nb113_dzH2, 352
.set nb113_dxM, 368
.set nb113_dyM, 384
.set nb113_dzM, 400
.set nb113_qqM, 416
.set nb113_qqH, 432
.set nb113_rinvH1, 448
.set nb113_rinvH2, 464
.set nb113_rinvM, 480
.set nb113_two, 496
.set nb113_c6, 512
.set nb113_c12, 528
.set nb113_six, 544
.set nb113_twelve, 560
.set nb113_vctot, 576
.set nb113_Vvdwtot, 592
.set nb113_fixO, 608
.set nb113_fiyO, 624
.set nb113_fizO, 640
.set nb113_fixH1, 656
.set nb113_fiyH1, 672
.set nb113_fizH1, 688
.set nb113_fixH2, 704
.set nb113_fiyH2, 720
.set nb113_fizH2, 736
.set nb113_fixM, 752
.set nb113_fiyM, 768
.set nb113_fizM, 784
.set nb113_fjx, 800
.set nb113_fjy, 816
.set nb113_fjz, 832
.set nb113_half, 848
.set nb113_three, 864
.set nb113_nri, 880
.set nb113_iinr, 888
.set nb113_jindex, 896
.set nb113_jjnr, 904
.set nb113_shift, 912
.set nb113_shiftvec, 920
.set nb113_facel, 928
.set nb113_innerjjnr, 936
.set nb113_ntia, 944
.set nb113_is3, 948
.set nb113_ii3, 952
.set nb113_innerk, 956
.set nb113_n, 960
.set nb113_nn1, 964
.set nb113_nouter, 968
.set nb113_ninner, 972

        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $984,%rsp

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb113_nouter(%rsp)
        movl %eax,nb113_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb113_nri(%rsp)
        movq %rsi,nb113_iinr(%rsp)
        movq %rdx,nb113_jindex(%rsp)
        movq %rcx,nb113_jjnr(%rsp)
        movq %r8,nb113_shift(%rsp)
        movq %r9,nb113_shiftvec(%rsp)
        movq nb113_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb113_facel(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb113_half(%rsp)
        movss nb113_half(%rsp),%xmm1
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
        movaps %xmm1,nb113_half(%rsp)
        movaps %xmm2,nb113_two(%rsp)
        movaps %xmm3,nb113_three(%rsp)
        movaps %xmm4,nb113_six(%rsp)
        movaps %xmm5,nb113_twelve(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb113_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb113_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm4
        movss 12(%rdx,%rbx,4),%xmm3
        movq nb113_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb113_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb113_iqM(%rsp)
        movaps %xmm4,nb113_iqH(%rsp)

        movq  nb113_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb113_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb113_ntia(%rsp)

_nb_kernel113_x86_64_sse.nb113_threadloop: 
        movq  nb113_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel113_x86_64_sse.nb113_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel113_x86_64_sse.nb113_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb113_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb113_n(%rsp)
        movl %ebx,nb113_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel113_x86_64_sse.nb113_outerstart
        jmp _nb_kernel113_x86_64_sse.nb113_end

_nb_kernel113_x86_64_sse.nb113_outerstart: 
        ## ebx contains number of outer iterations
        addl nb113_nouter(%rsp),%ebx
        movl %ebx,nb113_nouter(%rsp)

_nb_kernel113_x86_64_sse.nb113_outer: 
        movq  nb113_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb113_is3(%rsp)      ## store is3 

        movq  nb113_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb113_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb113_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb113_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3       ## ox
        addss 4(%rax,%rbx,4),%xmm4     ## oy
        addss 8(%rax,%rbx,4),%xmm5     ## oz
        addss 12(%rax,%rbx,4),%xmm6    ## h1x
        addss 16(%rax,%rbx,4),%xmm7    ## h1y
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm7,%xmm7
        movaps %xmm3,nb113_ixO(%rsp)
        movaps %xmm4,nb113_iyO(%rsp)
        movaps %xmm5,nb113_izO(%rsp)
        movaps %xmm6,nb113_ixH1(%rsp)
        movaps %xmm7,nb113_iyH1(%rsp)

        movss %xmm2,%xmm6
        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 20(%rax,%rbx,4),%xmm6    ## h1z
        addss 24(%rax,%rbx,4),%xmm0    ## h2x
        addss 28(%rax,%rbx,4),%xmm1    ## h2y
        addss 32(%rax,%rbx,4),%xmm2    ## h2z
        addss 36(%rax,%rbx,4),%xmm3    ## mx
        addss 40(%rax,%rbx,4),%xmm4    ## my
        addss 44(%rax,%rbx,4),%xmm5    ## mz

        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm6,nb113_izH1(%rsp)
        movaps %xmm0,nb113_ixH2(%rsp)
        movaps %xmm1,nb113_iyH2(%rsp)
        movaps %xmm2,nb113_izH2(%rsp)
        movaps %xmm3,nb113_ixM(%rsp)
        movaps %xmm4,nb113_iyM(%rsp)
        movaps %xmm5,nb113_izM(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb113_vctot(%rsp)
        movaps %xmm4,nb113_Vvdwtot(%rsp)
        movaps %xmm4,nb113_fixO(%rsp)
        movaps %xmm4,nb113_fiyO(%rsp)
        movaps %xmm4,nb113_fizO(%rsp)
        movaps %xmm4,nb113_fixH1(%rsp)
        movaps %xmm4,nb113_fiyH1(%rsp)
        movaps %xmm4,nb113_fizH1(%rsp)
        movaps %xmm4,nb113_fixH2(%rsp)
        movaps %xmm4,nb113_fiyH2(%rsp)
        movaps %xmm4,nb113_fizH2(%rsp)
        movaps %xmm4,nb113_fixM(%rsp)
        movaps %xmm4,nb113_fiyM(%rsp)
        movaps %xmm4,nb113_fizM(%rsp)

        movq  nb113_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb113_pos(%rbp),%rsi
        movq  nb113_faction(%rbp),%rdi
        movq  nb113_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb113_innerjjnr(%rsp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb113_ninner(%rsp),%ecx
        movl  %ecx,nb113_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb113_innerk(%rsp)   ## number of innerloop atoms 
        jge   _nb_kernel113_x86_64_sse.nb113_unroll_loop
        jmp   _nb_kernel113_x86_64_sse.nb113_odd_inner
_nb_kernel113_x86_64_sse.nb113_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb113_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb113_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb113_charge(%rbp),%rsi    ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb113_iqM(%rsp),%xmm3
        mulps  nb113_iqH(%rsp),%xmm4

        movaps  %xmm3,nb113_qqM(%rsp)
        movaps  %xmm4,nb113_qqH(%rsp)

        movq nb113_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movl (%rsi,%rcx,4),%r10d
        movl (%rsi,%rdx,4),%r11d
        movq nb113_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        shll %r10d
        shll %r11d
        movl nb113_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d
        addl %edi,%r10d
        addl %edi,%r11d

        movlps (%rsi,%r8,4),%xmm6
        movlps (%rsi,%r10,4),%xmm7
        movhps (%rsi,%r9,4),%xmm6
        movhps (%rsi,%r11,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm6 ## 11011101

        movaps %xmm4,nb113_c6(%rsp)
        movaps %xmm6,nb113_c12(%rsp)

        movq nb113_pos(%rbp),%rsi       ## base of pos[] 

        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx        ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

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

    ## interaction
    ## copy to xmm3-xmm5
    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5

    subps nb113_ixO(%rsp),%xmm3
    subps nb113_iyO(%rsp),%xmm4
    subps nb113_izO(%rsp),%xmm5

    movaps %xmm3,%xmm13
    movaps %xmm4,%xmm14
    movaps %xmm5,%xmm15

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3

    ## calc 1/rsq
    rcpps %xmm3,%xmm5
    movaps nb113_two(%rsp),%xmm4
    mulps %xmm5,%xmm3
    subps %xmm3,%xmm4
    mulps %xmm5,%xmm4       ## xmm4=rinvsq

    movaps %xmm4,%xmm3      ## rinvsq
    mulps  %xmm4,%xmm4      ## rinv4
    mulps  %xmm3,%xmm4      ## rinv6
    movaps %xmm4,%xmm5
    mulps  %xmm5,%xmm5      ## rinv12
    mulps  nb113_c6(%rsp),%xmm4
    mulps  nb113_c12(%rsp),%xmm5
    movaps %xmm5,%xmm6
    subps  %xmm4,%xmm6 ## Vvdw=vvdw12-vvdw6
    mulps  nb113_six(%rsp),%xmm4
    mulps  nb113_twelve(%rsp),%xmm5
    subps  %xmm4,%xmm5
    mulps  %xmm5,%xmm3  ## fscal

    addps  nb113_Vvdwtot(%rsp),%xmm6
    movaps %xmm6,nb113_Vvdwtot(%rsp)

    mulps  %xmm3,%xmm13 ## fx
    mulps  %xmm3,%xmm14 ## fy
    mulps  %xmm3,%xmm15 ## fz

    ## save j force temporarily
    movaps %xmm13,nb113_fjx(%rsp)
    movaps %xmm14,nb113_fjy(%rsp)
    movaps %xmm15,nb113_fjz(%rsp)

    ## increment i O force
    addps nb113_fixO(%rsp),%xmm13
    addps nb113_fiyO(%rsp),%xmm14
    addps nb113_fizO(%rsp),%xmm15
    movaps %xmm13,nb113_fixO(%rsp)
    movaps %xmm14,nb113_fiyO(%rsp)
    movaps %xmm15,nb113_fizO(%rsp)
    ## finished O LJ interaction.


    ## do H1, H2, and M interactions in parallel.
    ## xmm0-xmm2 still contain j coordinates.
    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8


    subps nb113_ixH1(%rsp),%xmm0
    subps nb113_iyH1(%rsp),%xmm1
    subps nb113_izH1(%rsp),%xmm2
    subps nb113_ixH2(%rsp),%xmm3
    subps nb113_iyH2(%rsp),%xmm4
    subps nb113_izH2(%rsp),%xmm5
    subps nb113_ixM(%rsp),%xmm6
    subps nb113_iyM(%rsp),%xmm7
    subps nb113_izM(%rsp),%xmm8

        movaps %xmm0,nb113_dxH1(%rsp)
        movaps %xmm1,nb113_dyH1(%rsp)
        movaps %xmm2,nb113_dzH1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb113_dxH2(%rsp)
        movaps %xmm4,nb113_dyH2(%rsp)
        movaps %xmm5,nb113_dzH2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb113_dxM(%rsp)
        movaps %xmm7,nb113_dyM(%rsp)
        movaps %xmm8,nb113_dzM(%rsp)
        mulps  %xmm6,%xmm6
        mulps  %xmm7,%xmm7
        mulps  %xmm8,%xmm8
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
    addps  %xmm7,%xmm6
    addps  %xmm8,%xmm6

        ## start doing invsqrt for j atoms
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm3,%xmm4
    rsqrtps %xmm6,%xmm7

        movaps  %xmm1,%xmm2
        movaps  %xmm4,%xmm5
    movaps  %xmm7,%xmm8

        mulps   %xmm1,%xmm1 ## lu*lu
        mulps   %xmm4,%xmm4 ## lu*lu
    mulps   %xmm7,%xmm7 ## lu*lu

        movaps  nb113_three(%rsp),%xmm9
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

        movaps  nb113_half(%rsp),%xmm0
        mulps   %xmm0,%xmm9 ## rinvH1
        mulps   %xmm0,%xmm10 ## rinvH2
    mulps   %xmm0,%xmm11 ## rinvM

        ## interactions 
    movaps %xmm9,%xmm0
    movaps %xmm10,%xmm1
    movaps %xmm11,%xmm2
    mulps  %xmm9,%xmm9
    mulps  %xmm10,%xmm10
    mulps  %xmm11,%xmm11
    mulps  nb113_qqH(%rsp),%xmm0
    mulps  nb113_qqH(%rsp),%xmm1
    mulps  nb113_qqM(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps nb113_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb113_vctot(%rsp)

        movq  nb113_faction(%rbp),%rdi
        ## move j forces to local temp variables 
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

        mulps nb113_dxH1(%rsp),%xmm7
        mulps nb113_dyH1(%rsp),%xmm8
        mulps nb113_dzH1(%rsp),%xmm9
        mulps nb113_dxH2(%rsp),%xmm10
        mulps nb113_dyH2(%rsp),%xmm11
        mulps nb113_dzH2(%rsp),%xmm12
        mulps nb113_dxM(%rsp),%xmm13
        mulps nb113_dyM(%rsp),%xmm14
        mulps nb113_dzM(%rsp),%xmm15

    ## fetch forces from O interaction
    movaps nb113_fjx(%rsp),%xmm3
    movaps nb113_fjy(%rsp),%xmm4
    addps  nb113_fjz(%rsp),%xmm2

    addps %xmm7,%xmm3
    addps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb113_fixH1(%rsp),%xmm7
    addps nb113_fiyH1(%rsp),%xmm8
    addps nb113_fizH1(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb113_fixH2(%rsp),%xmm10
    addps nb113_fiyH2(%rsp),%xmm11
    addps nb113_fizH2(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb113_fixM(%rsp),%xmm13
    addps nb113_fiyM(%rsp),%xmm14
    addps nb113_fizM(%rsp),%xmm15

    movaps %xmm7,nb113_fixH1(%rsp)
    movaps %xmm8,nb113_fiyH1(%rsp)
    movaps %xmm9,nb113_fizH1(%rsp)
    movaps %xmm10,nb113_fixH2(%rsp)
    movaps %xmm11,nb113_fiyH2(%rsp)
    movaps %xmm12,nb113_fizH2(%rsp)
    movaps %xmm13,nb113_fixM(%rsp)
    movaps %xmm14,nb113_fiyM(%rsp)
    movaps %xmm15,nb113_fizM(%rsp)

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
        subl $4,nb113_innerk(%rsp)
        jl    _nb_kernel113_x86_64_sse.nb113_odd_inner
        jmp   _nb_kernel113_x86_64_sse.nb113_unroll_loop
_nb_kernel113_x86_64_sse.nb113_odd_inner: 
        addl $4,nb113_innerk(%rsp)
        jnz   _nb_kernel113_x86_64_sse.nb113_odd_loop
        jmp   _nb_kernel113_x86_64_sse.nb113_updateouterdata
_nb_kernel113_x86_64_sse.nb113_odd_loop: 
        movq  nb113_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb113_innerjjnr(%rsp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb113_iqM(%rsp),%xmm4
        movq nb113_charge(%rbp),%rsi
        movhps nb113_iqH(%rsp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb113_qqM(%rsp)    ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movq nb113_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ebx
        movq nb113_vdwparam(%rbp),%rsi
        shll %ebx
        addl nb113_ntia(%rsp),%ebx
        movlps (%rsi,%rbx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## 11111100
        shufps $253,%xmm7,%xmm7 ## 11111101
        movaps %xmm6,nb113_c6(%rsp)
        movaps %xmm7,nb113_c12(%rsp)

        movq nb113_pos(%rbp),%rsi
        lea    (%rax,%rax,2),%rax

        movss nb113_ixO(%rsp),%xmm0
        movss nb113_iyO(%rsp),%xmm1
        movss nb113_izO(%rsp),%xmm2
        movss nb113_ixH1(%rsp),%xmm3
        movss nb113_iyH1(%rsp),%xmm4
        movss nb113_izH1(%rsp),%xmm5
        unpcklps nb113_ixH2(%rsp),%xmm0         ## ixO ixH2 - -
        unpcklps nb113_iyH2(%rsp),%xmm1         ## iyO iyH2 - -
        unpcklps nb113_izH2(%rsp),%xmm2         ## izO izH2 - -
        unpcklps nb113_ixM(%rsp),%xmm3          ## ixH1 ixM - -
        unpcklps nb113_iyM(%rsp),%xmm4          ## iyH1 iyM - -
        unpcklps nb113_izM(%rsp),%xmm5          ## izH1 izM - -
        unpcklps %xmm3,%xmm0    ## ixO ixH1 ixH2 ixM
        unpcklps %xmm4,%xmm1    ## same for y
        unpcklps %xmm5,%xmm2    ## same for z

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm3
        movss 4(%rsi,%rax,4),%xmm4
        movss 8(%rsi,%rax,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## use O distances for storage
        movaps %xmm3,nb113_dxO(%rsp)
        movaps %xmm4,nb113_dyO(%rsp)
        movaps %xmm5,nb113_dzO(%rsp)

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
        movaps nb113_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb113_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv

        movaps %xmm0,%xmm1
        mulps %xmm1,%xmm1       ## rinvsq 
        movaps %xmm0,%xmm7
        mulps  nb113_qqM(%rsp),%xmm7   ## vcoul
        movaps %xmm7,%xmm6

        addps  nb113_vctot(%rsp),%xmm7
        movaps %xmm7,nb113_vctot(%rsp)

        movaps %xmm1,%xmm2
        mulss  %xmm1,%xmm1
        mulss  %xmm2,%xmm1      ## xmm1=rinvsix
        xorps  %xmm4,%xmm4
        movss  %xmm1,%xmm4
        mulss  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulss  nb113_c6(%rsp),%xmm1
        mulss  nb113_c12(%rsp),%xmm4
        movaps %xmm4,%xmm3
        subss  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        mulss  nb113_six(%rsp),%xmm1
        mulss  nb113_twelve(%rsp),%xmm4
        subss  %xmm1,%xmm4
        addss  nb113_Vvdwtot(%rsp),%xmm3
        movss  %xmm3,nb113_Vvdwtot(%rsp)
        addps  %xmm6,%xmm4
        mulps  %xmm0,%xmm4
        mulps  %xmm0,%xmm4      ## fscal

        movaps nb113_dxO(%rsp),%xmm0
        movaps nb113_dyO(%rsp),%xmm1
        movaps nb113_dzO(%rsp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2 ## xmm0-xmm2 now contains tx-tz (partial force)

        movss  nb113_fixO(%rsp),%xmm3
        movss  nb113_fiyO(%rsp),%xmm4
        movss  nb113_fizO(%rsp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb113_fixO(%rsp)
        movss  %xmm4,nb113_fiyO(%rsp)
        movss  %xmm5,nb113_fizO(%rsp)   ## updated the O force now do the H's

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $0x39,%xmm3,%xmm3 ## shift right 
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb113_fixH1(%rsp),%xmm3
        addss  nb113_fiyH1(%rsp),%xmm4
        addss  nb113_fizH1(%rsp),%xmm5
        movss  %xmm3,nb113_fixH1(%rsp)
        movss  %xmm4,nb113_fiyH1(%rsp)
        movss  %xmm5,nb113_fizH1(%rsp)          ## updated the H1 force 

        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb113_fixH2(%rsp),%xmm3
        addss  nb113_fiyH2(%rsp),%xmm4
        addss  nb113_fizH2(%rsp),%xmm5
        movss  %xmm3,nb113_fixH2(%rsp)
        movss  %xmm4,nb113_fiyH2(%rsp)
        movss  %xmm5,nb113_fizH2(%rsp)          ## updated the H2 force 

        movq nb113_faction(%rbp),%rdi
        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb113_fixM(%rsp),%xmm3
        addss  nb113_fiyM(%rsp),%xmm4
        addss  nb113_fizM(%rsp),%xmm5
        movss  %xmm3,nb113_fixM(%rsp)
        movss  %xmm4,nb113_fiyM(%rsp)
        movss  %xmm5,nb113_fizM(%rsp)   ## updated the M force 

        ## the fj's - move in from mem start by acc. tx/ty/tz in xmm0, xmm1
        movlps (%rdi,%rax,4),%xmm6
        movss  8(%rdi,%rax,4),%xmm7

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,%xmm0
        movaps  %xmm4,%xmm1
        movaps  %xmm5,%xmm2

        shufps $0x39,%xmm3,%xmm3 ## shift right 
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2
        unpcklps %xmm1,%xmm0    ## x,y sum in xmm0, z sum in xmm2

        addps    %xmm0,%xmm6
        addss    %xmm2,%xmm7

        movlps %xmm6,(%rdi,%rax,4)
        movss  %xmm7,8(%rdi,%rax,4)

        decl nb113_innerk(%rsp)
        jz    _nb_kernel113_x86_64_sse.nb113_updateouterdata
        jmp   _nb_kernel113_x86_64_sse.nb113_odd_loop
_nb_kernel113_x86_64_sse.nb113_updateouterdata: 
        movl  nb113_ii3(%rsp),%ecx
        movq  nb113_faction(%rbp),%rdi
        movq  nb113_fshift(%rbp),%rsi
        movl  nb113_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb113_fixO(%rsp),%xmm0
        movaps nb113_fiyO(%rsp),%xmm1
        movaps nb113_fizO(%rsp),%xmm2

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
        movaps nb113_fixH1(%rsp),%xmm0
        movaps nb113_fiyH1(%rsp),%xmm1
        movaps nb113_fizH1(%rsp),%xmm2

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
        movaps nb113_fixH2(%rsp),%xmm0
        movaps nb113_fiyH2(%rsp),%xmm1
        movaps nb113_fizH2(%rsp),%xmm2

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

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movaps nb113_fixM(%rsp),%xmm0
        movaps nb113_fiyM(%rsp),%xmm1
        movaps nb113_fizM(%rsp),%xmm2

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
        movss  36(%rdi,%rcx,4),%xmm3
        movss  40(%rdi,%rcx,4),%xmm4
        movss  44(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,36(%rdi,%rcx,4)
        movss  %xmm4,40(%rdi,%rcx,4)
        movss  %xmm5,44(%rdi,%rcx,4)

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
        movl nb113_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb113_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb113_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb113_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb113_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb113_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb113_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel113_x86_64_sse.nb113_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb113_n(%rsp)
        jmp _nb_kernel113_x86_64_sse.nb113_outer
_nb_kernel113_x86_64_sse.nb113_outerend: 
        ## check if more outer neighborlists remain
        movl  nb113_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel113_x86_64_sse.nb113_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel113_x86_64_sse.nb113_threadloop
_nb_kernel113_x86_64_sse.nb113_end: 


        movl nb113_nouter(%rsp),%eax
        movl nb113_ninner(%rsp),%ebx
        movq nb113_outeriter(%rbp),%rcx
        movq nb113_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $984,%rsp
        emms

        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret



.globl nb_kernel113nf_x86_64_sse
.globl _nb_kernel113nf_x86_64_sse
nb_kernel113nf_x86_64_sse:      
_nb_kernel113nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb113nf_fshift, 16
.set nb113nf_gid, 24
.set nb113nf_pos, 32
.set nb113nf_faction, 40
.set nb113nf_charge, 48
.set nb113nf_p_facel, 56
.set nb113nf_argkrf, 64
.set nb113nf_argcrf, 72
.set nb113nf_Vc, 80
.set nb113nf_type, 88
.set nb113nf_p_ntype, 96
.set nb113nf_vdwparam, 104
.set nb113nf_Vvdw, 112
.set nb113nf_p_tabscale, 120
.set nb113nf_VFtab, 128
.set nb113nf_invsqrta, 136
.set nb113nf_dvda, 144
.set nb113nf_p_gbtabscale, 152
.set nb113nf_GBtab, 160
.set nb113nf_p_nthreads, 168
.set nb113nf_count, 176
.set nb113nf_mtx, 184
.set nb113nf_outeriter, 192
.set nb113nf_inneriter, 200
.set nb113nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb113nf_ixO, 0
.set nb113nf_iyO, 16
.set nb113nf_izO, 32
.set nb113nf_ixH1, 48
.set nb113nf_iyH1, 64
.set nb113nf_izH1, 80
.set nb113nf_ixH2, 96
.set nb113nf_iyH2, 112
.set nb113nf_izH2, 128
.set nb113nf_ixM, 144
.set nb113nf_iyM, 160
.set nb113nf_izM, 176
.set nb113nf_iqM, 192
.set nb113nf_iqH, 208
.set nb113nf_qqM, 224
.set nb113nf_qqH, 240
.set nb113nf_rinvH1, 256
.set nb113nf_rinvH2, 272
.set nb113nf_rinvM, 288
.set nb113nf_c6, 304
.set nb113nf_c12, 320
.set nb113nf_vctot, 336
.set nb113nf_Vvdwtot, 352
.set nb113nf_half, 368
.set nb113nf_three, 384
.set nb113nf_two, 400
.set nb113nf_facel, 416
.set nb113nf_iinr, 424
.set nb113nf_jindex, 432
.set nb113nf_jjnr, 440
.set nb113nf_shift, 448
.set nb113nf_shiftvec, 456
.set nb113nf_innerjjnr, 464
.set nb113nf_ntia, 472
.set nb113nf_nri, 476
.set nb113nf_is3, 480
.set nb113nf_ii3, 484
.set nb113nf_innerk, 488
.set nb113nf_n, 492
.set nb113nf_nn1, 496
.set nb113nf_nouter, 500
.set nb113nf_ninner, 504

        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $520,%rsp

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb113nf_nouter(%rsp)
        movl %eax,nb113nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb113nf_nri(%rsp)
        movq %rsi,nb113nf_iinr(%rsp)
        movq %rdx,nb113nf_jindex(%rsp)
        movq %rcx,nb113nf_jjnr(%rsp)
        movq %r8,nb113nf_shift(%rsp)
        movq %r9,nb113nf_shiftvec(%rsp)
        movq nb113nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb113nf_facel(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb113nf_half(%rsp)
        movss nb113nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb113nf_half(%rsp)
        movaps %xmm2,nb113nf_two(%rsp)
        movaps %xmm3,nb113nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb113nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb113nf_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm4
        movss 12(%rdx,%rbx,4),%xmm3
        movq nb113nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb113nf_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb113nf_iqM(%rsp)
        movaps %xmm4,nb113nf_iqH(%rsp)

        movq  nb113nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb113nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb113nf_ntia(%rsp)

_nb_kernel113nf_x86_64_sse.nb113nf_threadloop: 
        movq  nb113nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel113nf_x86_64_sse.nb113nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel113nf_x86_64_sse.nb113nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb113nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb113nf_n(%rsp)
        movl %ebx,nb113nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel113nf_x86_64_sse.nb113nf_outerstart
        jmp _nb_kernel113nf_x86_64_sse.nb113nf_end

_nb_kernel113nf_x86_64_sse.nb113nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb113nf_nouter(%rsp),%ebx
        movl %ebx,nb113nf_nouter(%rsp)

_nb_kernel113nf_x86_64_sse.nb113nf_outer: 
        movq  nb113nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb113nf_is3(%rsp)            ## store is3 

        movq  nb113nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb113nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb113nf_pos(%rbp),%rax    ## rax = base of pos[]  
        movl  %ebx,nb113nf_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3       ## ox
        addss 4(%rax,%rbx,4),%xmm4     ## oy
        addss 8(%rax,%rbx,4),%xmm5     ## oz
        addss 12(%rax,%rbx,4),%xmm6    ## h1x
        addss 16(%rax,%rbx,4),%xmm7    ## h1y
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm7,%xmm7
        movaps %xmm3,nb113nf_ixO(%rsp)
        movaps %xmm4,nb113nf_iyO(%rsp)
        movaps %xmm5,nb113nf_izO(%rsp)
        movaps %xmm6,nb113nf_ixH1(%rsp)
        movaps %xmm7,nb113nf_iyH1(%rsp)

        movss %xmm2,%xmm6
        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 20(%rax,%rbx,4),%xmm6    ## h1z
        addss 24(%rax,%rbx,4),%xmm0    ## h2x
        addss 28(%rax,%rbx,4),%xmm1    ## h2y
        addss 32(%rax,%rbx,4),%xmm2    ## h2z
        addss 36(%rax,%rbx,4),%xmm3    ## mx
        addss 40(%rax,%rbx,4),%xmm4    ## my
        addss 44(%rax,%rbx,4),%xmm5    ## mz

        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm6,nb113nf_izH1(%rsp)
        movaps %xmm0,nb113nf_ixH2(%rsp)
        movaps %xmm1,nb113nf_iyH2(%rsp)
        movaps %xmm2,nb113nf_izH2(%rsp)
        movaps %xmm3,nb113nf_ixM(%rsp)
        movaps %xmm4,nb113nf_iyM(%rsp)
        movaps %xmm5,nb113nf_izM(%rsp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb113nf_vctot(%rsp)
        movaps %xmm4,nb113nf_Vvdwtot(%rsp)

        movq  nb113nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb113nf_pos(%rbp),%rsi
        movq  nb113nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb113nf_innerjjnr(%rsp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb113nf_ninner(%rsp),%ecx
        movl  %ecx,nb113nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb113nf_innerk(%rsp)         ## number of innerloop atoms 
        jge   _nb_kernel113nf_x86_64_sse.nb113nf_unroll_loop
        jmp   _nb_kernel113nf_x86_64_sse.nb113nf_odd_inner
_nb_kernel113nf_x86_64_sse.nb113nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb113nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb113nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb113nf_charge(%rbp),%rsi  ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb113nf_iqM(%rsp),%xmm3
        mulps  nb113nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb113nf_qqM(%rsp)
        movaps  %xmm4,nb113nf_qqH(%rsp)

        movq nb113nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb113nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb113nf_ntia(%rsp),%edi
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

        movaps %xmm4,nb113nf_c6(%rsp)
        movaps %xmm6,nb113nf_c12(%rsp)

        movq nb113nf_pos(%rbp),%rsi     ## base of pos[] 

        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx        ## replace jnr with j3 
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
        movaps nb113nf_ixO(%rsp),%xmm4
        movaps nb113nf_iyO(%rsp),%xmm5
        movaps nb113nf_izO(%rsp),%xmm6

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
        movaps nb113nf_ixH1(%rsp),%xmm4
        movaps nb113nf_iyH1(%rsp),%xmm5
        movaps nb113nf_izH1(%rsp),%xmm6

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
        movaps nb113nf_ixH2(%rsp),%xmm3
        movaps nb113nf_iyH2(%rsp),%xmm4
        movaps nb113nf_izH2(%rsp),%xmm5

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

        ## move ixM-izM to xmm2-xmm4  
        movaps nb113nf_iyM(%rsp),%xmm3
        movaps nb113nf_izM(%rsp),%xmm4
        subps  %xmm1,%xmm3
        subps  %xmm2,%xmm4
        movaps nb113nf_ixM(%rsp),%xmm2
        subps  %xmm0,%xmm2

        ## square it 
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        addps %xmm3,%xmm4
        addps %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb113nf_three(%rsp),%xmm0
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb113nf_half(%rsp),%xmm0
        movaps  %xmm0,nb113nf_rinvH1(%rsp)      ## rinvH1 

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb113nf_three(%rsp),%xmm0
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb113nf_half(%rsp),%xmm0
        movaps  %xmm0,nb113nf_rinvH2(%rsp)      ## rinvH2 

        ## rsqM - seed to xmm2 
        rsqrtps %xmm4,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb113nf_three(%rsp),%xmm0
        mulps   %xmm4,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb113nf_half(%rsp),%xmm0
        movaps  %xmm0,nb113nf_rinvM(%rsp)

        ## Do the O LJ-only interaction directly.       
        rcpps   %xmm7,%xmm2
        movaps  nb113nf_two(%rsp),%xmm1
        mulps   %xmm2,%xmm7
        subps   %xmm7,%xmm1
        mulps   %xmm1,%xmm2 ## rinvsq 
        movaps  %xmm2,%xmm0
        mulps   %xmm2,%xmm0     ## r4
        mulps   %xmm2,%xmm0     ## r6
        movaps  %xmm0,%xmm1
        mulps   %xmm1,%xmm1     ## r12
        mulps   nb113nf_c6(%rsp),%xmm0
        mulps   nb113nf_c12(%rsp),%xmm1
        movaps  %xmm1,%xmm3
        subps   %xmm0,%xmm3     ## Vvdw12-Vvdw6
        addps   nb113nf_Vvdwtot(%rsp),%xmm3
        movaps  %xmm3,nb113nf_Vvdwtot(%rsp)

        ## Do H1,H2,M interactions
        movaps  nb113nf_rinvH1(%rsp),%xmm7
        movaps  nb113nf_rinvM(%rsp),%xmm6
        addps   nb113nf_rinvH2(%rsp),%xmm7
        mulps  nb113nf_qqH(%rsp),%xmm7
        mulps  nb113nf_qqM(%rsp),%xmm6
        addps  %xmm6,%xmm7

        addps  nb113nf_vctot(%rsp),%xmm7
        movaps %xmm7,nb113nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb113nf_innerk(%rsp)
        jl    _nb_kernel113nf_x86_64_sse.nb113nf_odd_inner
        jmp   _nb_kernel113nf_x86_64_sse.nb113nf_unroll_loop
_nb_kernel113nf_x86_64_sse.nb113nf_odd_inner: 
        addl $4,nb113nf_innerk(%rsp)
        jnz   _nb_kernel113nf_x86_64_sse.nb113nf_odd_loop
        jmp   _nb_kernel113nf_x86_64_sse.nb113nf_updateouterdata
_nb_kernel113nf_x86_64_sse.nb113nf_odd_loop: 
        movq  nb113nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb113nf_innerjjnr(%rsp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb113nf_iqM(%rsp),%xmm4
        movq nb113nf_charge(%rbp),%rsi
        movhps nb113nf_iqH(%rsp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb113nf_qqM(%rsp)          ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movq nb113nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ebx
        movq nb113nf_vdwparam(%rbp),%rsi
        shll %ebx
        addl nb113nf_ntia(%rsp),%ebx
        movlps (%rsi,%rbx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## 11111100
        shufps $253,%xmm7,%xmm7 ## 11111101
        movaps %xmm6,nb113nf_c6(%rsp)
        movaps %xmm7,nb113nf_c12(%rsp)

        movq nb113nf_pos(%rbp),%rsi
        lea (%rax,%rax,2),%rax

        movss nb113nf_ixO(%rsp),%xmm3
        movss nb113nf_iyO(%rsp),%xmm4
        movss nb113nf_izO(%rsp),%xmm5
        movss nb113nf_ixH1(%rsp),%xmm0
        movss nb113nf_iyH1(%rsp),%xmm1
        movss nb113nf_izH1(%rsp),%xmm2
        unpcklps nb113nf_ixH2(%rsp),%xmm3       ## ixO ixH2 - -
        unpcklps nb113nf_iyH2(%rsp),%xmm4       ## iyO iyH2 - -
        unpcklps nb113nf_izH2(%rsp),%xmm5       ## izO izH2 - -
        unpcklps nb113nf_ixM(%rsp),%xmm0        ## ixH1 ixM - -
        unpcklps nb113nf_iyM(%rsp),%xmm1        ## iyH1 iyM - -
        unpcklps nb113nf_izM(%rsp),%xmm2        ## izH1 izM - -
        unpcklps %xmm0,%xmm3    ## ixO ixH1 ixH2 ixM
        unpcklps %xmm1,%xmm4    ## same for y
        unpcklps %xmm2,%xmm5    ## same for z

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm0
        movss 4(%rsi,%rax,4),%xmm1
        movss 8(%rsi,%rax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

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
        movaps nb113nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb113nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm7
        mulps  nb113nf_qqM(%rsp),%xmm7   ## vcoul
        addps  nb113nf_vctot(%rsp),%xmm7
        movaps %xmm7,nb113nf_vctot(%rsp)

        mulss  %xmm1,%xmm1
        movaps %xmm1,%xmm2
        mulss  %xmm1,%xmm1
        mulss  %xmm2,%xmm1      ## xmm1=rinvsix
        xorps  %xmm4,%xmm4
        movss  %xmm1,%xmm4
        mulss  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulss  nb113nf_c6(%rsp),%xmm1
        mulss  nb113nf_c12(%rsp),%xmm4
        movaps %xmm4,%xmm3
        subss  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addss  nb113nf_Vvdwtot(%rsp),%xmm3
        movss  %xmm3,nb113nf_Vvdwtot(%rsp)

        decl nb113nf_innerk(%rsp)
        jz    _nb_kernel113nf_x86_64_sse.nb113nf_updateouterdata
        jmp   _nb_kernel113nf_x86_64_sse.nb113nf_odd_loop
_nb_kernel113nf_x86_64_sse.nb113nf_updateouterdata: 
        ## get n from stack
        movl nb113nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb113nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb113nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb113nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb113nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb113nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb113nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel113nf_x86_64_sse.nb113nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb113nf_n(%rsp)
        jmp _nb_kernel113nf_x86_64_sse.nb113nf_outer
_nb_kernel113nf_x86_64_sse.nb113nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb113nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel113nf_x86_64_sse.nb113nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel113nf_x86_64_sse.nb113nf_threadloop
_nb_kernel113nf_x86_64_sse.nb113nf_end: 


        movl nb113nf_nouter(%rsp),%eax
        movl nb113nf_ninner(%rsp),%ebx
        movq nb113nf_outeriter(%rbp),%rcx
        movq nb113nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $520,%rsp
        emms

        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret


