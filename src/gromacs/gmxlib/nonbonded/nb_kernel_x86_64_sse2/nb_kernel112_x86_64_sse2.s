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







.globl nb_kernel112_x86_64_sse2
.globl _nb_kernel112_x86_64_sse2
nb_kernel112_x86_64_sse2:       
_nb_kernel112_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb112_fshift, 16
.set nb112_gid, 24
.set nb112_pos, 32
.set nb112_faction, 40
.set nb112_charge, 48
.set nb112_p_facel, 56
.set nb112_argkrf, 64
.set nb112_argcrf, 72
.set nb112_Vc, 80
.set nb112_type, 88
.set nb112_p_ntype, 96
.set nb112_vdwparam, 104
.set nb112_Vvdw, 112
.set nb112_p_tabscale, 120
.set nb112_VFtab, 128
.set nb112_invsqrta, 136
.set nb112_dvda, 144
.set nb112_p_gbtabscale, 152
.set nb112_GBtab, 160
.set nb112_p_nthreads, 168
.set nb112_count, 176
.set nb112_mtx, 184
.set nb112_outeriter, 192
.set nb112_inneriter, 200
.set nb112_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb112_ixO, 0
.set nb112_iyO, 16
.set nb112_izO, 32
.set nb112_ixH1, 48
.set nb112_iyH1, 64
.set nb112_izH1, 80
.set nb112_ixH2, 96
.set nb112_iyH2, 112
.set nb112_izH2, 128
.set nb112_jxO, 144
.set nb112_jyO, 160
.set nb112_jzO, 176
.set nb112_jxH1, 192
.set nb112_jyH1, 208
.set nb112_jzH1, 224
.set nb112_jxH2, 240
.set nb112_jyH2, 256
.set nb112_jzH2, 272
.set nb112_dxOO, 288
.set nb112_dyOO, 304
.set nb112_dzOO, 320
.set nb112_dxOH1, 336
.set nb112_dyOH1, 352
.set nb112_dzOH1, 368
.set nb112_dxOH2, 384
.set nb112_dyOH2, 400
.set nb112_dzOH2, 416
.set nb112_dxH1O, 432
.set nb112_dyH1O, 448
.set nb112_dzH1O, 464
.set nb112_dxH1H1, 480
.set nb112_dyH1H1, 496
.set nb112_dzH1H1, 512
.set nb112_dxH1H2, 528
.set nb112_dyH1H2, 544
.set nb112_dzH1H2, 560
.set nb112_dxH2O, 576
.set nb112_dyH2O, 592
.set nb112_dzH2O, 608
.set nb112_dxH2H1, 624
.set nb112_dyH2H1, 640
.set nb112_dzH2H1, 656
.set nb112_dxH2H2, 672
.set nb112_dyH2H2, 688
.set nb112_dzH2H2, 704
.set nb112_qqOO, 720
.set nb112_qqOH, 736
.set nb112_qqHH, 752
.set nb112_c6, 768
.set nb112_c12, 784
.set nb112_six, 800
.set nb112_twelve, 816
.set nb112_vctot, 832
.set nb112_Vvdwtot, 848
.set nb112_fixO, 864
.set nb112_fiyO, 880
.set nb112_fizO, 896
.set nb112_fixH1, 912
.set nb112_fiyH1, 928
.set nb112_fizH1, 944
.set nb112_fixH2, 960
.set nb112_fiyH2, 976
.set nb112_fizH2, 992
.set nb112_fjxO, 1008
.set nb112_fjyO, 1024
.set nb112_fjzO, 1040
.set nb112_fjxH1, 1056
.set nb112_fjyH1, 1072
.set nb112_fjzH1, 1088
.set nb112_fjxH2, 1104
.set nb112_fjyH2, 1120
.set nb112_fjzH2, 1136
.set nb112_half, 1152
.set nb112_three, 1168
.set nb112_rsqOO, 1184
.set nb112_rsqOH1, 1200
.set nb112_rsqOH2, 1216
.set nb112_rsqH1O, 1232
.set nb112_rsqH1H1, 1248
.set nb112_rsqH1H2, 1264
.set nb112_rsqH2O, 1280
.set nb112_rsqH2H1, 1296
.set nb112_rsqH2H2, 1312
.set nb112_rinvOO, 1328
.set nb112_rinvOH1, 1344
.set nb112_rinvOH2, 1360
.set nb112_rinvH1O, 1376
.set nb112_rinvH1H1, 1392
.set nb112_rinvH1H2, 1408
.set nb112_rinvH2O, 1424
.set nb112_rinvH2H1, 1440
.set nb112_rinvH2H2, 1456
.set nb112_is3, 1472
.set nb112_ii3, 1476
.set nb112_nri, 1492
.set nb112_iinr, 1500
.set nb112_jindex, 1508
.set nb112_jjnr, 1516
.set nb112_shift, 1524
.set nb112_shiftvec, 1532
.set nb112_facel, 1540
.set nb112_innerjjnr, 1548
.set nb112_innerk, 1556
.set nb112_n, 1560
.set nb112_nn1, 1564
.set nb112_nouter, 1568
.set nb112_ninner, 1572
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1592,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb112_nouter(%rsp)
        movl %eax,nb112_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb112_nri(%rsp)
        movq %rsi,nb112_iinr(%rsp)
        movq %rdx,nb112_jindex(%rsp)
        movq %rcx,nb112_jjnr(%rsp)
        movq %r8,nb112_shift(%rsp)
        movq %r9,nb112_shiftvec(%rsp)
        movq nb112_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb112_facel(%rsp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb112_half(%rsp)
        movl %ebx,nb112_half+4(%rsp)
        movsd nb112_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm3,%xmm4
        addpd  %xmm4,%xmm4      ## six
        movapd %xmm4,%xmm5
        addpd  %xmm5,%xmm5      ## twelve
        movapd %xmm1,nb112_half(%rsp)
        movapd %xmm3,nb112_three(%rsp)
        movapd %xmm4,nb112_six(%rsp)
        movapd %xmm5,nb112_twelve(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb112_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb112_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%rdx,%rbx,8),%xmm5
        movq nb112_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb112_facel(%rsp),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb112_qqOO(%rsp)
        movapd %xmm4,nb112_qqOH(%rsp)
        movapd %xmm5,nb112_qqHH(%rsp)

        xorpd %xmm0,%xmm0
        movq  nb112_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb112_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  nb112_vdwparam(%rbp),%rax
        movlpd (%rax,%rdx,8),%xmm0
        movlpd 8(%rax,%rdx,8),%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb112_c6(%rsp)
        movapd %xmm1,nb112_c12(%rsp)

_nb_kernel112_x86_64_sse2.nb112_threadloop: 
        movq  nb112_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel112_x86_64_sse2.nb112_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel112_x86_64_sse2.nb112_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb112_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb112_n(%rsp)
        movl %ebx,nb112_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel112_x86_64_sse2.nb112_outerstart
        jmp _nb_kernel112_x86_64_sse2.nb112_end

_nb_kernel112_x86_64_sse2.nb112_outerstart: 
        ## ebx contains number of outer iterations
        addl nb112_nouter(%rsp),%ebx
        movl %ebx,nb112_nouter(%rsp)

_nb_kernel112_x86_64_sse2.nb112_outer: 
        movq  nb112_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb112_is3(%rsp)      ## store is3 

        movq  nb112_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb112_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb112_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb112_ii3(%rsp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%rax,%rbx,8),%xmm3
        addsd 8(%rax,%rbx,8),%xmm4
        addsd 16(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb112_ixO(%rsp)
        movapd %xmm4,nb112_iyO(%rsp)
        movapd %xmm5,nb112_izO(%rsp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%rax,%rbx,8),%xmm0
        addsd 32(%rax,%rbx,8),%xmm1
        addsd 40(%rax,%rbx,8),%xmm2
        addsd 48(%rax,%rbx,8),%xmm3
        addsd 56(%rax,%rbx,8),%xmm4
        addsd 64(%rax,%rbx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb112_ixH1(%rsp)
        movapd %xmm1,nb112_iyH1(%rsp)
        movapd %xmm2,nb112_izH1(%rsp)
        movapd %xmm3,nb112_ixH2(%rsp)
        movapd %xmm4,nb112_iyH2(%rsp)
        movapd %xmm5,nb112_izH2(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb112_vctot(%rsp)
        movapd %xmm4,nb112_Vvdwtot(%rsp)
        movapd %xmm4,nb112_fixO(%rsp)
        movapd %xmm4,nb112_fiyO(%rsp)
        movapd %xmm4,nb112_fizO(%rsp)
        movapd %xmm4,nb112_fixH1(%rsp)
        movapd %xmm4,nb112_fiyH1(%rsp)
        movapd %xmm4,nb112_fizH1(%rsp)
        movapd %xmm4,nb112_fixH2(%rsp)
        movapd %xmm4,nb112_fiyH2(%rsp)
        movapd %xmm4,nb112_fizH2(%rsp)

        movq  nb112_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx     ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb112_pos(%rbp),%rsi
        movq  nb112_faction(%rbp),%rdi
        movq  nb112_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb112_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb112_ninner(%rsp),%ecx
        movl  %ecx,nb112_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb112_innerk(%rsp)      ## number of innerloop atoms 
        jge  _nb_kernel112_x86_64_sse2.nb112_unroll_loop
        jmp  _nb_kernel112_x86_64_sse2.nb112_checksingle
_nb_kernel112_x86_64_sse2.nb112_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb112_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb112_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 

        movq nb112_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move j O coordinates to local temp variables 
    movlpd (%rsi,%rax,8),%xmm0
    movlpd 8(%rsi,%rax,8),%xmm1
    movlpd 16(%rsi,%rax,8),%xmm2
    movhpd (%rsi,%rbx,8),%xmm0
    movhpd 8(%rsi,%rbx,8),%xmm1
    movhpd 16(%rsi,%rbx,8),%xmm2

    ## xmm0 = Ox
    ## xmm1 = Oy
    ## xmm2 = Oz

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb112_ixO(%rsp),%xmm0
    subpd nb112_iyO(%rsp),%xmm1
    subpd nb112_izO(%rsp),%xmm2
    subpd nb112_ixH1(%rsp),%xmm3
    subpd nb112_iyH1(%rsp),%xmm4
    subpd nb112_izH1(%rsp),%xmm5
    subpd nb112_ixH2(%rsp),%xmm6
    subpd nb112_iyH2(%rsp),%xmm7
    subpd nb112_izH2(%rsp),%xmm8

        movapd %xmm0,nb112_dxOO(%rsp)
        movapd %xmm1,nb112_dyOO(%rsp)
        movapd %xmm2,nb112_dzOO(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb112_dxH1O(%rsp)
        movapd %xmm4,nb112_dyH1O(%rsp)
        movapd %xmm5,nb112_dzH1O(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb112_dxH2O(%rsp)
        movapd %xmm7,nb112_dyH2O(%rsp)
        movapd %xmm8,nb112_dzH2O(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for jO atoms
    cvtpd2ps %xmm0,%xmm1
    cvtpd2ps %xmm3,%xmm4
    cvtpd2ps %xmm6,%xmm7
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm4,%xmm4
    rsqrtps %xmm7,%xmm7
    cvtps2pd %xmm1,%xmm1
    cvtps2pd %xmm4,%xmm4
    cvtps2pd %xmm7,%xmm7

        movapd  %xmm1,%xmm2
        movapd  %xmm4,%xmm5
    movapd  %xmm7,%xmm8

        mulpd   %xmm1,%xmm1 ## lu*lu
        mulpd   %xmm4,%xmm4 ## lu*lu
    mulpd   %xmm7,%xmm7 ## lu*lu

        movapd  nb112_three(%rsp),%xmm9
        movapd  %xmm9,%xmm10
    movapd  %xmm9,%xmm11

        mulpd   %xmm0,%xmm1 ## rsq*lu*lu
        mulpd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm7 ## rsq*lu*lu

        subpd   %xmm1,%xmm9
        subpd   %xmm4,%xmm10
    subpd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulpd   %xmm2,%xmm9
        mulpd   %xmm5,%xmm10
    mulpd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb112_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvOO 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH1O
    mulpd   %xmm15,%xmm11 ## first iteration for rinvH2O 

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb112_three(%rsp),%xmm1
        movapd  %xmm1,%xmm4
    movapd  %xmm1,%xmm7

        mulpd   %xmm0,%xmm2 ## rsq*lu*lu
        mulpd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm8 ## rsq*lu*lu

        subpd   %xmm2,%xmm1
        subpd   %xmm5,%xmm4
    subpd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulpd   %xmm1,%xmm9
        mulpd   %xmm4,%xmm10
    mulpd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb112_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvOO 
        mulpd   %xmm15,%xmm10 ##   rinvH1O
    mulpd   %xmm15,%xmm11 ##   rinvH2O


        ## O interactions 
    movapd %xmm9,%xmm0
    movapd %xmm10,%xmm1
    movapd %xmm11,%xmm2
    mulpd  %xmm9,%xmm9   ## rinvsq
    mulpd  %xmm10,%xmm10
    mulpd  %xmm11,%xmm11
    movapd %xmm9,%xmm12
    mulpd  %xmm12,%xmm12 ## rinv4
    mulpd  %xmm9,%xmm12 ## rinv6
    mulpd  nb112_qqOO(%rsp),%xmm0
    mulpd  nb112_qqOH(%rsp),%xmm1
    mulpd  nb112_qqOH(%rsp),%xmm2
    movapd %xmm12,%xmm13 ## rinv6
    mulpd  %xmm12,%xmm12 ## rinv12
        mulpd  nb112_c6(%rsp),%xmm13
        mulpd  nb112_c12(%rsp),%xmm12
    movapd %xmm12,%xmm14
    subpd  %xmm13,%xmm14

        addpd  nb112_Vvdwtot(%rsp),%xmm14
        mulpd  nb112_six(%rsp),%xmm13
        mulpd  nb112_twelve(%rsp),%xmm12
        movapd %xmm14,nb112_Vvdwtot(%rsp)
    subpd  %xmm13,%xmm12 ## LJ fscal        

    addpd  %xmm0,%xmm12

    mulpd  %xmm12,%xmm9
    mulpd  %xmm1,%xmm10
    mulpd  %xmm2,%xmm11

    addpd nb112_vctot(%rsp),%xmm0
    addpd %xmm2,%xmm1
    addpd %xmm1,%xmm0
    movapd %xmm0,nb112_vctot(%rsp)

    ## move j O forces to xmm0-xmm2
    movq nb112_faction(%rbp),%rdi
        movlpd (%rdi,%rax,8),%xmm0
        movlpd 8(%rdi,%rax,8),%xmm1
        movlpd 16(%rdi,%rax,8),%xmm2

    movapd %xmm9,%xmm7
    movapd %xmm9,%xmm8
    movapd %xmm11,%xmm13
    movapd %xmm11,%xmm14
    movapd %xmm11,%xmm15
    movapd %xmm10,%xmm11
    movapd %xmm10,%xmm12

        movhpd (%rdi,%rbx,8),%xmm0
        movhpd 8(%rdi,%rbx,8),%xmm1
        movhpd 16(%rdi,%rbx,8),%xmm2

        mulpd nb112_dxOO(%rsp),%xmm7
        mulpd nb112_dyOO(%rsp),%xmm8
        mulpd nb112_dzOO(%rsp),%xmm9
        mulpd nb112_dxH1O(%rsp),%xmm10
        mulpd nb112_dyH1O(%rsp),%xmm11
        mulpd nb112_dzH1O(%rsp),%xmm12
        mulpd nb112_dxH2O(%rsp),%xmm13
        mulpd nb112_dyH2O(%rsp),%xmm14
        mulpd nb112_dzH2O(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb112_fixO(%rsp),%xmm7
    addpd nb112_fiyO(%rsp),%xmm8
    addpd nb112_fizO(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb112_fixH1(%rsp),%xmm10
    addpd nb112_fiyH1(%rsp),%xmm11
    addpd nb112_fizH1(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb112_fixH2(%rsp),%xmm13
    addpd nb112_fiyH2(%rsp),%xmm14
    addpd nb112_fizH2(%rsp),%xmm15

    movapd %xmm7,nb112_fixO(%rsp)
    movapd %xmm8,nb112_fiyO(%rsp)
    movapd %xmm9,nb112_fizO(%rsp)
    movapd %xmm10,nb112_fixH1(%rsp)
    movapd %xmm11,nb112_fiyH1(%rsp)
    movapd %xmm12,nb112_fizH1(%rsp)
    movapd %xmm13,nb112_fixH2(%rsp)
    movapd %xmm14,nb112_fiyH2(%rsp)
    movapd %xmm15,nb112_fizH2(%rsp)

    ## store back j O forces from xmm0-xmm2
        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)
        movhpd %xmm0,(%rdi,%rbx,8)
        movhpd %xmm1,8(%rdi,%rbx,8)
        movhpd %xmm2,16(%rdi,%rbx,8)

        ## move j H1 coordinates to local temp variables 
    movq nb112_pos(%rbp),%rsi
    movlpd 24(%rsi,%rax,8),%xmm0
    movlpd 32(%rsi,%rax,8),%xmm1
    movlpd 40(%rsi,%rax,8),%xmm2
    movhpd 24(%rsi,%rbx,8),%xmm0
    movhpd 32(%rsi,%rbx,8),%xmm1
    movhpd 40(%rsi,%rbx,8),%xmm2

    ## xmm0 = H1x
    ## xmm1 = H1y
    ## xmm2 = H1z

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb112_ixO(%rsp),%xmm0
    subpd nb112_iyO(%rsp),%xmm1
    subpd nb112_izO(%rsp),%xmm2
    subpd nb112_ixH1(%rsp),%xmm3
    subpd nb112_iyH1(%rsp),%xmm4
    subpd nb112_izH1(%rsp),%xmm5
    subpd nb112_ixH2(%rsp),%xmm6
    subpd nb112_iyH2(%rsp),%xmm7
    subpd nb112_izH2(%rsp),%xmm8

        movapd %xmm0,nb112_dxOH1(%rsp)
        movapd %xmm1,nb112_dyOH1(%rsp)
        movapd %xmm2,nb112_dzOH1(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb112_dxH1H1(%rsp)
        movapd %xmm4,nb112_dyH1H1(%rsp)
        movapd %xmm5,nb112_dzH1H1(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb112_dxH2H1(%rsp)
        movapd %xmm7,nb112_dyH2H1(%rsp)
        movapd %xmm8,nb112_dzH2H1(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for jH1 atoms
    cvtpd2ps %xmm0,%xmm1
    cvtpd2ps %xmm3,%xmm4
    cvtpd2ps %xmm6,%xmm7
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm4,%xmm4
    rsqrtps %xmm7,%xmm7
    cvtps2pd %xmm1,%xmm1
    cvtps2pd %xmm4,%xmm4
    cvtps2pd %xmm7,%xmm7

        movapd  %xmm1,%xmm2
        movapd  %xmm4,%xmm5
    movapd  %xmm7,%xmm8

        mulpd   %xmm1,%xmm1 ## lu*lu
        mulpd   %xmm4,%xmm4 ## lu*lu
    mulpd   %xmm7,%xmm7 ## lu*lu

        movapd  nb112_three(%rsp),%xmm9
        movapd  %xmm9,%xmm10
    movapd  %xmm9,%xmm11

        mulpd   %xmm0,%xmm1 ## rsq*lu*lu
        mulpd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm7 ## rsq*lu*lu

        subpd   %xmm1,%xmm9
        subpd   %xmm4,%xmm10
    subpd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulpd   %xmm2,%xmm9
        mulpd   %xmm5,%xmm10
    mulpd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb112_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvOH1 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH1H1
    mulpd   %xmm15,%xmm11 ## first iteration for rinvH2OH1

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb112_three(%rsp),%xmm1
        movapd  %xmm1,%xmm4
    movapd  %xmm1,%xmm7

        mulpd   %xmm0,%xmm2 ## rsq*lu*lu
        mulpd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm8 ## rsq*lu*lu

        subpd   %xmm2,%xmm1
        subpd   %xmm5,%xmm4
    subpd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulpd   %xmm1,%xmm9
        mulpd   %xmm4,%xmm10
    mulpd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb112_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvOH1
        mulpd   %xmm15,%xmm10 ##   rinvH1H1
    mulpd   %xmm15,%xmm11 ##   rinvH2H1

        ## H1 interactions 
    movapd %xmm9,%xmm0
    movapd %xmm10,%xmm1
    movapd %xmm11,%xmm2
    mulpd  %xmm9,%xmm9
    mulpd  %xmm10,%xmm10
    mulpd  %xmm11,%xmm11
    mulpd  nb112_qqOH(%rsp),%xmm0
    mulpd  nb112_qqHH(%rsp),%xmm1
    mulpd  nb112_qqHH(%rsp),%xmm2
    mulpd  %xmm0,%xmm9
    mulpd  %xmm1,%xmm10
    mulpd  %xmm2,%xmm11

    addpd nb112_vctot(%rsp),%xmm0
    addpd %xmm2,%xmm1
    addpd %xmm1,%xmm0
    movapd %xmm0,nb112_vctot(%rsp)

    ## move j H1 forces to xmm0-xmm2
    movq nb112_faction(%rbp),%rdi
        movlpd 24(%rdi,%rax,8),%xmm0
        movlpd 32(%rdi,%rax,8),%xmm1
        movlpd 40(%rdi,%rax,8),%xmm2

    movapd %xmm9,%xmm7
    movapd %xmm9,%xmm8
    movapd %xmm11,%xmm13
    movapd %xmm11,%xmm14
    movapd %xmm11,%xmm15
    movapd %xmm10,%xmm11
    movapd %xmm10,%xmm12

        movhpd 24(%rdi,%rbx,8),%xmm0
        movhpd 32(%rdi,%rbx,8),%xmm1
        movhpd 40(%rdi,%rbx,8),%xmm2

        mulpd nb112_dxOH1(%rsp),%xmm7
        mulpd nb112_dyOH1(%rsp),%xmm8
        mulpd nb112_dzOH1(%rsp),%xmm9
        mulpd nb112_dxH1H1(%rsp),%xmm10
        mulpd nb112_dyH1H1(%rsp),%xmm11
        mulpd nb112_dzH1H1(%rsp),%xmm12
        mulpd nb112_dxH2H1(%rsp),%xmm13
        mulpd nb112_dyH2H1(%rsp),%xmm14
        mulpd nb112_dzH2H1(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb112_fixO(%rsp),%xmm7
    addpd nb112_fiyO(%rsp),%xmm8
    addpd nb112_fizO(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb112_fixH1(%rsp),%xmm10
    addpd nb112_fiyH1(%rsp),%xmm11
    addpd nb112_fizH1(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb112_fixH2(%rsp),%xmm13
    addpd nb112_fiyH2(%rsp),%xmm14
    addpd nb112_fizH2(%rsp),%xmm15

    movapd %xmm7,nb112_fixO(%rsp)
    movapd %xmm8,nb112_fiyO(%rsp)
    movapd %xmm9,nb112_fizO(%rsp)
    movapd %xmm10,nb112_fixH1(%rsp)
    movapd %xmm11,nb112_fiyH1(%rsp)
    movapd %xmm12,nb112_fizH1(%rsp)
    movapd %xmm13,nb112_fixH2(%rsp)
    movapd %xmm14,nb112_fiyH2(%rsp)
    movapd %xmm15,nb112_fizH2(%rsp)

    ## store back j H1 forces from xmm0-xmm2
        movlpd %xmm0,24(%rdi,%rax,8)
        movlpd %xmm1,32(%rdi,%rax,8)
        movlpd %xmm2,40(%rdi,%rax,8)
        movhpd %xmm0,24(%rdi,%rbx,8)
        movhpd %xmm1,32(%rdi,%rbx,8)
        movhpd %xmm2,40(%rdi,%rbx,8)

        ## move j H2 coordinates to local temp variables 
    movq nb112_pos(%rbp),%rsi
    movlpd 48(%rsi,%rax,8),%xmm0
    movlpd 56(%rsi,%rax,8),%xmm1
    movlpd 64(%rsi,%rax,8),%xmm2
    movhpd 48(%rsi,%rbx,8),%xmm0
    movhpd 56(%rsi,%rbx,8),%xmm1
    movhpd 64(%rsi,%rbx,8),%xmm2

    ## xmm0 = H2x
    ## xmm1 = H2y
    ## xmm2 = H2z

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb112_ixO(%rsp),%xmm0
    subpd nb112_iyO(%rsp),%xmm1
    subpd nb112_izO(%rsp),%xmm2
    subpd nb112_ixH1(%rsp),%xmm3
    subpd nb112_iyH1(%rsp),%xmm4
    subpd nb112_izH1(%rsp),%xmm5
    subpd nb112_ixH2(%rsp),%xmm6
    subpd nb112_iyH2(%rsp),%xmm7
    subpd nb112_izH2(%rsp),%xmm8

        movapd %xmm0,nb112_dxOH2(%rsp)
        movapd %xmm1,nb112_dyOH2(%rsp)
        movapd %xmm2,nb112_dzOH2(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb112_dxH1H2(%rsp)
        movapd %xmm4,nb112_dyH1H2(%rsp)
        movapd %xmm5,nb112_dzH1H2(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb112_dxH2H2(%rsp)
        movapd %xmm7,nb112_dyH2H2(%rsp)
        movapd %xmm8,nb112_dzH2H2(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for jH2 atoms
    cvtpd2ps %xmm0,%xmm1
    cvtpd2ps %xmm3,%xmm4
    cvtpd2ps %xmm6,%xmm7
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm4,%xmm4
    rsqrtps %xmm7,%xmm7
    cvtps2pd %xmm1,%xmm1
    cvtps2pd %xmm4,%xmm4
    cvtps2pd %xmm7,%xmm7

        movapd  %xmm1,%xmm2
        movapd  %xmm4,%xmm5
    movapd  %xmm7,%xmm8

        mulpd   %xmm1,%xmm1 ## lu*lu
        mulpd   %xmm4,%xmm4 ## lu*lu
    mulpd   %xmm7,%xmm7 ## lu*lu

        movapd  nb112_three(%rsp),%xmm9
        movapd  %xmm9,%xmm10
    movapd  %xmm9,%xmm11

        mulpd   %xmm0,%xmm1 ## rsq*lu*lu
        mulpd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm7 ## rsq*lu*lu

        subpd   %xmm1,%xmm9
        subpd   %xmm4,%xmm10
    subpd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulpd   %xmm2,%xmm9
        mulpd   %xmm5,%xmm10
    mulpd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb112_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvOH2 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH1H2
    mulpd   %xmm15,%xmm11 ## first iteration for rinvH2H2

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb112_three(%rsp),%xmm1
        movapd  %xmm1,%xmm4
    movapd  %xmm1,%xmm7

        mulpd   %xmm0,%xmm2 ## rsq*lu*lu
        mulpd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm8 ## rsq*lu*lu

        subpd   %xmm2,%xmm1
        subpd   %xmm5,%xmm4
    subpd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulpd   %xmm1,%xmm9
        mulpd   %xmm4,%xmm10
    mulpd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb112_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvOH2
        mulpd   %xmm15,%xmm10 ##   rinvH1H2
    mulpd   %xmm15,%xmm11 ##   rinvH2H2

        ## H2 interactions 
    movapd %xmm9,%xmm0
    movapd %xmm10,%xmm1
    movapd %xmm11,%xmm2
    mulpd  %xmm9,%xmm9
    mulpd  %xmm10,%xmm10
    mulpd  %xmm11,%xmm11
    mulpd  nb112_qqOH(%rsp),%xmm0
    mulpd  nb112_qqHH(%rsp),%xmm1
    mulpd  nb112_qqHH(%rsp),%xmm2
    mulpd  %xmm0,%xmm9
    mulpd  %xmm1,%xmm10
    mulpd  %xmm2,%xmm11

    addpd nb112_vctot(%rsp),%xmm0
    addpd %xmm2,%xmm1
    addpd %xmm1,%xmm0
    movapd %xmm0,nb112_vctot(%rsp)

    ## move j H2 forces to xmm0-xmm2
    movq nb112_faction(%rbp),%rdi
        movlpd 48(%rdi,%rax,8),%xmm0
        movlpd 56(%rdi,%rax,8),%xmm1
        movlpd 64(%rdi,%rax,8),%xmm2

    movapd %xmm9,%xmm7
    movapd %xmm9,%xmm8
    movapd %xmm11,%xmm13
    movapd %xmm11,%xmm14
    movapd %xmm11,%xmm15
    movapd %xmm10,%xmm11
    movapd %xmm10,%xmm12

        movhpd 48(%rdi,%rbx,8),%xmm0
        movhpd 56(%rdi,%rbx,8),%xmm1
        movhpd 64(%rdi,%rbx,8),%xmm2

        mulpd nb112_dxOH2(%rsp),%xmm7
        mulpd nb112_dyOH2(%rsp),%xmm8
        mulpd nb112_dzOH2(%rsp),%xmm9
        mulpd nb112_dxH1H2(%rsp),%xmm10
        mulpd nb112_dyH1H2(%rsp),%xmm11
        mulpd nb112_dzH1H2(%rsp),%xmm12
        mulpd nb112_dxH2H2(%rsp),%xmm13
        mulpd nb112_dyH2H2(%rsp),%xmm14
        mulpd nb112_dzH2H2(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb112_fixO(%rsp),%xmm7
    addpd nb112_fiyO(%rsp),%xmm8
    addpd nb112_fizO(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb112_fixH1(%rsp),%xmm10
    addpd nb112_fiyH1(%rsp),%xmm11
    addpd nb112_fizH1(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb112_fixH2(%rsp),%xmm13
    addpd nb112_fiyH2(%rsp),%xmm14
    addpd nb112_fizH2(%rsp),%xmm15

    movapd %xmm7,nb112_fixO(%rsp)
    movapd %xmm8,nb112_fiyO(%rsp)
    movapd %xmm9,nb112_fizO(%rsp)
    movapd %xmm10,nb112_fixH1(%rsp)
    movapd %xmm11,nb112_fiyH1(%rsp)
    movapd %xmm12,nb112_fizH1(%rsp)
    movapd %xmm13,nb112_fixH2(%rsp)
    movapd %xmm14,nb112_fiyH2(%rsp)
    movapd %xmm15,nb112_fizH2(%rsp)

    ## store back j H2 forces from xmm0-xmm2
        movlpd %xmm0,48(%rdi,%rax,8)
        movlpd %xmm1,56(%rdi,%rax,8)
        movlpd %xmm2,64(%rdi,%rax,8)
        movhpd %xmm0,48(%rdi,%rbx,8)
        movhpd %xmm1,56(%rdi,%rbx,8)
        movhpd %xmm2,64(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb112_innerk(%rsp)
        jl    _nb_kernel112_x86_64_sse2.nb112_checksingle
        jmp   _nb_kernel112_x86_64_sse2.nb112_unroll_loop
_nb_kernel112_x86_64_sse2.nb112_checksingle: 
        movl  nb112_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel112_x86_64_sse2.nb112_dosingle
        jmp   _nb_kernel112_x86_64_sse2.nb112_updateouterdata
_nb_kernel112_x86_64_sse2.nb112_dosingle: 
        movq  nb112_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb112_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax


        ## move j O coordinates to local temp variables 
    movsd (%rsi,%rax,8),%xmm0
    movsd 8(%rsi,%rax,8),%xmm1
    movsd 16(%rsi,%rax,8),%xmm2

    ## xmm0 = Ox
    ## xmm1 = Oy
    ## xmm2 = Oz

    movsd %xmm0,%xmm3
    movsd %xmm1,%xmm4
    movsd %xmm2,%xmm5
    movsd %xmm0,%xmm6
    movsd %xmm1,%xmm7
    movsd %xmm2,%xmm8

    subsd nb112_ixO(%rsp),%xmm0
    subsd nb112_iyO(%rsp),%xmm1
    subsd nb112_izO(%rsp),%xmm2
    subsd nb112_ixH1(%rsp),%xmm3
    subsd nb112_iyH1(%rsp),%xmm4
    subsd nb112_izH1(%rsp),%xmm5
    subsd nb112_ixH2(%rsp),%xmm6
    subsd nb112_iyH2(%rsp),%xmm7
    subsd nb112_izH2(%rsp),%xmm8

        movsd %xmm0,nb112_dxOO(%rsp)
        movsd %xmm1,nb112_dyOO(%rsp)
        movsd %xmm2,nb112_dzOO(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb112_dxH1O(%rsp)
        movsd %xmm4,nb112_dyH1O(%rsp)
        movsd %xmm5,nb112_dzH1O(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb112_dxH2O(%rsp)
        movsd %xmm7,nb112_dyH2O(%rsp)
        movsd %xmm8,nb112_dzH2O(%rsp)
        mulsd  %xmm6,%xmm6
        mulsd  %xmm7,%xmm7
        mulsd  %xmm8,%xmm8
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    addsd  %xmm7,%xmm6
    addsd  %xmm8,%xmm6

        ## start doing invsqrt for jO atoms
    cvtsd2ss %xmm0,%xmm1
    cvtsd2ss %xmm3,%xmm4
    cvtsd2ss %xmm6,%xmm7
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm4,%xmm4
    rsqrtss %xmm7,%xmm7
    cvtss2sd %xmm1,%xmm1
    cvtss2sd %xmm4,%xmm4
    cvtss2sd %xmm7,%xmm7

        movsd  %xmm1,%xmm2
        movsd  %xmm4,%xmm5
    movsd  %xmm7,%xmm8

        mulsd   %xmm1,%xmm1 ## lu*lu
        mulsd   %xmm4,%xmm4 ## lu*lu
    mulsd   %xmm7,%xmm7 ## lu*lu

        movsd  nb112_three(%rsp),%xmm9
        movsd  %xmm9,%xmm10
    movsd  %xmm9,%xmm11

        mulsd   %xmm0,%xmm1 ## rsq*lu*lu
        mulsd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm7 ## rsq*lu*lu

        subsd   %xmm1,%xmm9
        subsd   %xmm4,%xmm10
    subsd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulsd   %xmm2,%xmm9
        mulsd   %xmm5,%xmm10
    mulsd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb112_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvOO 
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH1O
    mulsd   %xmm15,%xmm11 ## first iteration for rinvH2O 

    ## second iteration step    
        movsd  %xmm9,%xmm2
        movsd  %xmm10,%xmm5
    movsd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movsd  nb112_three(%rsp),%xmm1
        movsd  %xmm1,%xmm4
    movsd  %xmm1,%xmm7

        mulsd   %xmm0,%xmm2 ## rsq*lu*lu
        mulsd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm8 ## rsq*lu*lu

        subsd   %xmm2,%xmm1
        subsd   %xmm5,%xmm4
    subsd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulsd   %xmm1,%xmm9
        mulsd   %xmm4,%xmm10
    mulsd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb112_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvOO 
        mulsd   %xmm15,%xmm10 ##   rinvH1O
    mulsd   %xmm15,%xmm11 ##   rinvH2O

        ## O interactions 
    movsd %xmm9,%xmm0
    movsd %xmm10,%xmm1
    movsd %xmm11,%xmm2
    mulsd  %xmm9,%xmm9   ## rinvsq
    mulsd  %xmm10,%xmm10
    mulsd  %xmm11,%xmm11
    movsd %xmm9,%xmm12
    mulsd  %xmm12,%xmm12 ## rinv4
    mulsd  %xmm9,%xmm12 ## rinv6
    mulsd  nb112_qqOO(%rsp),%xmm0
    mulsd  nb112_qqOH(%rsp),%xmm1
    mulsd  nb112_qqOH(%rsp),%xmm2
    movsd %xmm12,%xmm13 ## rinv6
    mulsd  %xmm12,%xmm12 ## rinv12
        mulsd  nb112_c6(%rsp),%xmm13
        mulsd  nb112_c12(%rsp),%xmm12
    movsd %xmm12,%xmm14
    subsd  %xmm13,%xmm14

        addsd  nb112_Vvdwtot(%rsp),%xmm14
        mulsd  nb112_six(%rsp),%xmm13
        mulsd  nb112_twelve(%rsp),%xmm12
        movsd %xmm14,nb112_Vvdwtot(%rsp)
    subsd  %xmm13,%xmm12 ## LJ fscal        

    addsd  %xmm0,%xmm12

    mulsd  %xmm12,%xmm9
    mulsd  %xmm1,%xmm10
    mulsd  %xmm2,%xmm11

    addsd nb112_vctot(%rsp),%xmm0
    addsd %xmm2,%xmm1
    addsd %xmm1,%xmm0
    movsd %xmm0,nb112_vctot(%rsp)

    ## move j O forces to xmm0-xmm2
        movsd (%rdi,%rax,8),%xmm0
        movsd 8(%rdi,%rax,8),%xmm1
        movsd 16(%rdi,%rax,8),%xmm2

    movsd %xmm9,%xmm7
    movsd %xmm9,%xmm8
    movsd %xmm11,%xmm13
    movsd %xmm11,%xmm14
    movsd %xmm11,%xmm15
    movsd %xmm10,%xmm11
    movsd %xmm10,%xmm12

        mulsd nb112_dxOO(%rsp),%xmm7
        mulsd nb112_dyOO(%rsp),%xmm8
        mulsd nb112_dzOO(%rsp),%xmm9
        mulsd nb112_dxH1O(%rsp),%xmm10
        mulsd nb112_dyH1O(%rsp),%xmm11
        mulsd nb112_dzH1O(%rsp),%xmm12
        mulsd nb112_dxH2O(%rsp),%xmm13
        mulsd nb112_dyH2O(%rsp),%xmm14
        mulsd nb112_dzH2O(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb112_fixO(%rsp),%xmm7
    addsd nb112_fiyO(%rsp),%xmm8
    addsd nb112_fizO(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb112_fixH1(%rsp),%xmm10
    addsd nb112_fiyH1(%rsp),%xmm11
    addsd nb112_fizH1(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb112_fixH2(%rsp),%xmm13
    addsd nb112_fiyH2(%rsp),%xmm14
    addsd nb112_fizH2(%rsp),%xmm15

    movsd %xmm7,nb112_fixO(%rsp)
    movsd %xmm8,nb112_fiyO(%rsp)
    movsd %xmm9,nb112_fizO(%rsp)
    movsd %xmm10,nb112_fixH1(%rsp)
    movsd %xmm11,nb112_fiyH1(%rsp)
    movsd %xmm12,nb112_fizH1(%rsp)
    movsd %xmm13,nb112_fixH2(%rsp)
    movsd %xmm14,nb112_fiyH2(%rsp)
    movsd %xmm15,nb112_fizH2(%rsp)

    ## store back j O forces from xmm0-xmm2
        movsd %xmm0,(%rdi,%rax,8)
        movsd %xmm1,8(%rdi,%rax,8)
        movsd %xmm2,16(%rdi,%rax,8)

        ## move j H1 coordinates to local temp variables 
    movsd 24(%rsi,%rax,8),%xmm0
    movsd 32(%rsi,%rax,8),%xmm1
    movsd 40(%rsi,%rax,8),%xmm2

    ## xmm0 = H1x
    ## xmm1 = H1y
    ## xmm2 = H1z

    movsd %xmm0,%xmm3
    movsd %xmm1,%xmm4
    movsd %xmm2,%xmm5
    movsd %xmm0,%xmm6
    movsd %xmm1,%xmm7
    movsd %xmm2,%xmm8

    subsd nb112_ixO(%rsp),%xmm0
    subsd nb112_iyO(%rsp),%xmm1
    subsd nb112_izO(%rsp),%xmm2
    subsd nb112_ixH1(%rsp),%xmm3
    subsd nb112_iyH1(%rsp),%xmm4
    subsd nb112_izH1(%rsp),%xmm5
    subsd nb112_ixH2(%rsp),%xmm6
    subsd nb112_iyH2(%rsp),%xmm7
    subsd nb112_izH2(%rsp),%xmm8

        movsd %xmm0,nb112_dxOH1(%rsp)
        movsd %xmm1,nb112_dyOH1(%rsp)
        movsd %xmm2,nb112_dzOH1(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb112_dxH1H1(%rsp)
        movsd %xmm4,nb112_dyH1H1(%rsp)
        movsd %xmm5,nb112_dzH1H1(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb112_dxH2H1(%rsp)
        movsd %xmm7,nb112_dyH2H1(%rsp)
        movsd %xmm8,nb112_dzH2H1(%rsp)
        mulsd  %xmm6,%xmm6
        mulsd  %xmm7,%xmm7
        mulsd  %xmm8,%xmm8
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    addsd  %xmm7,%xmm6
    addsd  %xmm8,%xmm6

        ## start doing invsqrt for jH1 atoms
    cvtsd2ss %xmm0,%xmm1
    cvtsd2ss %xmm3,%xmm4
    cvtsd2ss %xmm6,%xmm7
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm4,%xmm4
    rsqrtss %xmm7,%xmm7
    cvtss2sd %xmm1,%xmm1
    cvtss2sd %xmm4,%xmm4
    cvtss2sd %xmm7,%xmm7

        movsd  %xmm1,%xmm2
        movsd  %xmm4,%xmm5
    movsd  %xmm7,%xmm8

        mulsd   %xmm1,%xmm1 ## lu*lu
        mulsd   %xmm4,%xmm4 ## lu*lu
    mulsd   %xmm7,%xmm7 ## lu*lu

        movsd  nb112_three(%rsp),%xmm9
        movsd  %xmm9,%xmm10
    movsd  %xmm9,%xmm11

        mulsd   %xmm0,%xmm1 ## rsq*lu*lu
        mulsd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm7 ## rsq*lu*lu

        subsd   %xmm1,%xmm9
        subsd   %xmm4,%xmm10
    subsd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulsd   %xmm2,%xmm9
        mulsd   %xmm5,%xmm10
    mulsd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb112_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvOH1 
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH1H1
    mulsd   %xmm15,%xmm11 ## first iteration for rinvH2OH1

    ## second iteration step    
        movsd  %xmm9,%xmm2
        movsd  %xmm10,%xmm5
    movsd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movsd  nb112_three(%rsp),%xmm1
        movsd  %xmm1,%xmm4
    movsd  %xmm1,%xmm7

        mulsd   %xmm0,%xmm2 ## rsq*lu*lu
        mulsd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm8 ## rsq*lu*lu

        subsd   %xmm2,%xmm1
        subsd   %xmm5,%xmm4
    subsd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulsd   %xmm1,%xmm9
        mulsd   %xmm4,%xmm10
    mulsd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb112_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvOH1
        mulsd   %xmm15,%xmm10 ##   rinvH1H1
    mulsd   %xmm15,%xmm11 ##   rinvH2H1

        ## H1 interactions 
    movsd %xmm9,%xmm0
    movsd %xmm10,%xmm1
    movsd %xmm11,%xmm2
    mulsd  %xmm9,%xmm9
    mulsd  %xmm10,%xmm10
    mulsd  %xmm11,%xmm11
    mulsd  nb112_qqOH(%rsp),%xmm0
    mulsd  nb112_qqHH(%rsp),%xmm1
    mulsd  nb112_qqHH(%rsp),%xmm2
    mulsd  %xmm0,%xmm9
    mulsd  %xmm1,%xmm10
    mulsd  %xmm2,%xmm11

    addsd nb112_vctot(%rsp),%xmm0
    addsd %xmm2,%xmm1
    addsd %xmm1,%xmm0
    movsd %xmm0,nb112_vctot(%rsp)

    ## move j H1 forces to xmm0-xmm2
        movsd 24(%rdi,%rax,8),%xmm0
        movsd 32(%rdi,%rax,8),%xmm1
        movsd 40(%rdi,%rax,8),%xmm2

    movsd %xmm9,%xmm7
    movsd %xmm9,%xmm8
    movsd %xmm11,%xmm13
    movsd %xmm11,%xmm14
    movsd %xmm11,%xmm15
    movsd %xmm10,%xmm11
    movsd %xmm10,%xmm12

        mulsd nb112_dxOH1(%rsp),%xmm7
        mulsd nb112_dyOH1(%rsp),%xmm8
        mulsd nb112_dzOH1(%rsp),%xmm9
        mulsd nb112_dxH1H1(%rsp),%xmm10
        mulsd nb112_dyH1H1(%rsp),%xmm11
        mulsd nb112_dzH1H1(%rsp),%xmm12
        mulsd nb112_dxH2H1(%rsp),%xmm13
        mulsd nb112_dyH2H1(%rsp),%xmm14
        mulsd nb112_dzH2H1(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb112_fixO(%rsp),%xmm7
    addsd nb112_fiyO(%rsp),%xmm8
    addsd nb112_fizO(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb112_fixH1(%rsp),%xmm10
    addsd nb112_fiyH1(%rsp),%xmm11
    addsd nb112_fizH1(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb112_fixH2(%rsp),%xmm13
    addsd nb112_fiyH2(%rsp),%xmm14
    addsd nb112_fizH2(%rsp),%xmm15

    movsd %xmm7,nb112_fixO(%rsp)
    movsd %xmm8,nb112_fiyO(%rsp)
    movsd %xmm9,nb112_fizO(%rsp)
    movsd %xmm10,nb112_fixH1(%rsp)
    movsd %xmm11,nb112_fiyH1(%rsp)
    movsd %xmm12,nb112_fizH1(%rsp)
    movsd %xmm13,nb112_fixH2(%rsp)
    movsd %xmm14,nb112_fiyH2(%rsp)
    movsd %xmm15,nb112_fizH2(%rsp)

    ## store back j H1 forces from xmm0-xmm2
        movsd %xmm0,24(%rdi,%rax,8)
        movsd %xmm1,32(%rdi,%rax,8)
        movsd %xmm2,40(%rdi,%rax,8)

        ## move j H2 coordinates to local temp variables 
    movsd 48(%rsi,%rax,8),%xmm0
    movsd 56(%rsi,%rax,8),%xmm1
    movsd 64(%rsi,%rax,8),%xmm2

    ## xmm0 = H2x
    ## xmm1 = H2y
    ## xmm2 = H2z

    movsd %xmm0,%xmm3
    movsd %xmm1,%xmm4
    movsd %xmm2,%xmm5
    movsd %xmm0,%xmm6
    movsd %xmm1,%xmm7
    movsd %xmm2,%xmm8

    subsd nb112_ixO(%rsp),%xmm0
    subsd nb112_iyO(%rsp),%xmm1
    subsd nb112_izO(%rsp),%xmm2
    subsd nb112_ixH1(%rsp),%xmm3
    subsd nb112_iyH1(%rsp),%xmm4
    subsd nb112_izH1(%rsp),%xmm5
    subsd nb112_ixH2(%rsp),%xmm6
    subsd nb112_iyH2(%rsp),%xmm7
    subsd nb112_izH2(%rsp),%xmm8

        movsd %xmm0,nb112_dxOH2(%rsp)
        movsd %xmm1,nb112_dyOH2(%rsp)
        movsd %xmm2,nb112_dzOH2(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb112_dxH1H2(%rsp)
        movsd %xmm4,nb112_dyH1H2(%rsp)
        movsd %xmm5,nb112_dzH1H2(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb112_dxH2H2(%rsp)
        movsd %xmm7,nb112_dyH2H2(%rsp)
        movsd %xmm8,nb112_dzH2H2(%rsp)
        mulsd  %xmm6,%xmm6
        mulsd  %xmm7,%xmm7
        mulsd  %xmm8,%xmm8
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    addsd  %xmm7,%xmm6
    addsd  %xmm8,%xmm6

        ## start doing invsqrt for jH2 atoms
    cvtsd2ss %xmm0,%xmm1
    cvtsd2ss %xmm3,%xmm4
    cvtsd2ss %xmm6,%xmm7
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm4,%xmm4
    rsqrtss %xmm7,%xmm7
    cvtss2sd %xmm1,%xmm1
    cvtss2sd %xmm4,%xmm4
    cvtss2sd %xmm7,%xmm7

        movsd  %xmm1,%xmm2
        movsd  %xmm4,%xmm5
    movsd  %xmm7,%xmm8

        mulsd   %xmm1,%xmm1 ## lu*lu
        mulsd   %xmm4,%xmm4 ## lu*lu
    mulsd   %xmm7,%xmm7 ## lu*lu

        movsd  nb112_three(%rsp),%xmm9
        movsd  %xmm9,%xmm10
    movsd  %xmm9,%xmm11

        mulsd   %xmm0,%xmm1 ## rsq*lu*lu
        mulsd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm7 ## rsq*lu*lu

        subsd   %xmm1,%xmm9
        subsd   %xmm4,%xmm10
    subsd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulsd   %xmm2,%xmm9
        mulsd   %xmm5,%xmm10
    mulsd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb112_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvOH2 
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH1H2
    mulsd   %xmm15,%xmm11 ## first iteration for rinvH2H2

    ## second iteration step    
        movsd  %xmm9,%xmm2
        movsd  %xmm10,%xmm5
    movsd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movsd  nb112_three(%rsp),%xmm1
        movsd  %xmm1,%xmm4
    movsd  %xmm1,%xmm7

        mulsd   %xmm0,%xmm2 ## rsq*lu*lu
        mulsd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm8 ## rsq*lu*lu

        subsd   %xmm2,%xmm1
        subsd   %xmm5,%xmm4
    subsd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulsd   %xmm1,%xmm9
        mulsd   %xmm4,%xmm10
    mulsd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb112_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvOH2
        mulsd   %xmm15,%xmm10 ##   rinvH1H2
    mulsd   %xmm15,%xmm11 ##   rinvH2H2

        ## H2 interactions 
    movsd %xmm9,%xmm0
    movsd %xmm10,%xmm1
    movsd %xmm11,%xmm2
    mulsd  %xmm9,%xmm9
    mulsd  %xmm10,%xmm10
    mulsd  %xmm11,%xmm11
    mulsd  nb112_qqOH(%rsp),%xmm0
    mulsd  nb112_qqHH(%rsp),%xmm1
    mulsd  nb112_qqHH(%rsp),%xmm2
    mulsd  %xmm0,%xmm9
    mulsd  %xmm1,%xmm10
    mulsd  %xmm2,%xmm11

    addsd nb112_vctot(%rsp),%xmm0
    addsd %xmm2,%xmm1
    addsd %xmm1,%xmm0
    movsd %xmm0,nb112_vctot(%rsp)

    ## move j H2 forces to xmm0-xmm2
        movsd 48(%rdi,%rax,8),%xmm0
        movsd 56(%rdi,%rax,8),%xmm1
        movsd 64(%rdi,%rax,8),%xmm2

    movsd %xmm9,%xmm7
    movsd %xmm9,%xmm8
    movsd %xmm11,%xmm13
    movsd %xmm11,%xmm14
    movsd %xmm11,%xmm15
    movsd %xmm10,%xmm11
    movsd %xmm10,%xmm12

        mulsd nb112_dxOH2(%rsp),%xmm7
        mulsd nb112_dyOH2(%rsp),%xmm8
        mulsd nb112_dzOH2(%rsp),%xmm9
        mulsd nb112_dxH1H2(%rsp),%xmm10
        mulsd nb112_dyH1H2(%rsp),%xmm11
        mulsd nb112_dzH1H2(%rsp),%xmm12
        mulsd nb112_dxH2H2(%rsp),%xmm13
        mulsd nb112_dyH2H2(%rsp),%xmm14
        mulsd nb112_dzH2H2(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb112_fixO(%rsp),%xmm7
    addsd nb112_fiyO(%rsp),%xmm8
    addsd nb112_fizO(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb112_fixH1(%rsp),%xmm10
    addsd nb112_fiyH1(%rsp),%xmm11
    addsd nb112_fizH1(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb112_fixH2(%rsp),%xmm13
    addsd nb112_fiyH2(%rsp),%xmm14
    addsd nb112_fizH2(%rsp),%xmm15

    movsd %xmm7,nb112_fixO(%rsp)
    movsd %xmm8,nb112_fiyO(%rsp)
    movsd %xmm9,nb112_fizO(%rsp)
    movsd %xmm10,nb112_fixH1(%rsp)
    movsd %xmm11,nb112_fiyH1(%rsp)
    movsd %xmm12,nb112_fizH1(%rsp)
    movsd %xmm13,nb112_fixH2(%rsp)
    movsd %xmm14,nb112_fiyH2(%rsp)
    movsd %xmm15,nb112_fizH2(%rsp)

    ## store back j H2 forces from xmm0-xmm2
        movsd %xmm0,48(%rdi,%rax,8)
        movsd %xmm1,56(%rdi,%rax,8)
        movsd %xmm2,64(%rdi,%rax,8)

_nb_kernel112_x86_64_sse2.nb112_updateouterdata: 
        movl  nb112_ii3(%rsp),%ecx
        movq  nb112_faction(%rbp),%rdi
        movq  nb112_fshift(%rbp),%rsi
        movl  nb112_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb112_fixO(%rsp),%xmm0
        movapd nb112_fiyO(%rsp),%xmm1
        movapd nb112_fizO(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        ## increment i force 
        movsd  (%rdi,%rcx,8),%xmm3
        movsd  8(%rdi,%rcx,8),%xmm4
        movsd  16(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,(%rdi,%rcx,8)
        movsd  %xmm4,8(%rdi,%rcx,8)
        movsd  %xmm5,16(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        movapd %xmm0,%xmm6
        movsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm6

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb112_fixH1(%rsp),%xmm0
        movapd nb112_fiyH1(%rsp),%xmm1
        movapd nb112_fizH1(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  24(%rdi,%rcx,8),%xmm3
        movsd  32(%rdi,%rcx,8),%xmm4
        movsd  40(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,24(%rdi,%rcx,8)
        movsd  %xmm4,32(%rdi,%rcx,8)
        movsd  %xmm5,40(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb112_fixH2(%rsp),%xmm0
        movapd nb112_fiyH2(%rsp),%xmm1
        movapd nb112_fizH2(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        ## increment i force 
        movsd  48(%rdi,%rcx,8),%xmm3
        movsd  56(%rdi,%rcx,8),%xmm4
        movsd  64(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,48(%rdi,%rcx,8)
        movsd  %xmm4,56(%rdi,%rcx,8)
        movsd  %xmm5,64(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## increment fshift force 
        movlpd (%rsi,%rdx,8),%xmm3
        movhpd 8(%rsi,%rdx,8),%xmm3
        movsd  16(%rsi,%rdx,8),%xmm4
        subpd  %xmm6,%xmm3
        subsd  %xmm7,%xmm4
        movlpd %xmm3,(%rsi,%rdx,8)
        movhpd %xmm3,8(%rsi,%rdx,8)
        movsd  %xmm4,16(%rsi,%rdx,8)

        ## get n from stack
        movl nb112_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb112_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb112_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb112_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb112_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb112_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

       ## finish if last 
        movl nb112_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel112_x86_64_sse2.nb112_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb112_n(%rsp)
        jmp _nb_kernel112_x86_64_sse2.nb112_outer
_nb_kernel112_x86_64_sse2.nb112_outerend: 
        ## check if more outer neighborlists remain
        movl  nb112_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel112_x86_64_sse2.nb112_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel112_x86_64_sse2.nb112_threadloop
_nb_kernel112_x86_64_sse2.nb112_end: 
        movl nb112_nouter(%rsp),%eax
        movl nb112_ninner(%rsp),%ebx
        movq nb112_outeriter(%rbp),%rcx
        movq nb112_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1592,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




.globl nb_kernel112nf_x86_64_sse2
.globl _nb_kernel112nf_x86_64_sse2
nb_kernel112nf_x86_64_sse2:     
_nb_kernel112nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb112nf_fshift, 16
.set nb112nf_gid, 24
.set nb112nf_pos, 32
.set nb112nf_faction, 40
.set nb112nf_charge, 48
.set nb112nf_p_facel, 56
.set nb112nf_argkrf, 64
.set nb112nf_argcrf, 72
.set nb112nf_Vc, 80
.set nb112nf_type, 88
.set nb112nf_p_ntype, 96
.set nb112nf_vdwparam, 104
.set nb112nf_Vvdw, 112
.set nb112nf_p_tabscale, 120
.set nb112nf_VFtab, 128
.set nb112nf_invsqrta, 136
.set nb112nf_dvda, 144
.set nb112nf_p_gbtabscale, 152
.set nb112nf_GBtab, 160
.set nb112nf_p_nthreads, 168
.set nb112nf_count, 176
.set nb112nf_mtx, 184
.set nb112nf_outeriter, 192
.set nb112nf_inneriter, 200
.set nb112nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb112nf_ixO, 0
.set nb112nf_iyO, 16
.set nb112nf_izO, 32
.set nb112nf_ixH1, 48
.set nb112nf_iyH1, 64
.set nb112nf_izH1, 80
.set nb112nf_ixH2, 96
.set nb112nf_iyH2, 112
.set nb112nf_izH2, 128
.set nb112nf_jxO, 144
.set nb112nf_jyO, 160
.set nb112nf_jzO, 176
.set nb112nf_jxH1, 192
.set nb112nf_jyH1, 208
.set nb112nf_jzH1, 224
.set nb112nf_jxH2, 240
.set nb112nf_jyH2, 256
.set nb112nf_jzH2, 272
.set nb112nf_qqOO, 288
.set nb112nf_qqOH, 304
.set nb112nf_qqHH, 320
.set nb112nf_c6, 336
.set nb112nf_c12, 352
.set nb112nf_vctot, 368
.set nb112nf_Vvdwtot, 384
.set nb112nf_half, 400
.set nb112nf_three, 416
.set nb112nf_rsqOO, 432
.set nb112nf_rsqOH1, 448
.set nb112nf_rsqOH2, 464
.set nb112nf_rsqH1O, 480
.set nb112nf_rsqH1H1, 496
.set nb112nf_rsqH1H2, 512
.set nb112nf_rsqH2O, 528
.set nb112nf_rsqH2H1, 544
.set nb112nf_rsqH2H2, 560
.set nb112nf_rinvOO, 576
.set nb112nf_rinvOH1, 592
.set nb112nf_rinvOH2, 608
.set nb112nf_rinvH1O, 624
.set nb112nf_rinvH1H1, 640
.set nb112nf_rinvH1H2, 656
.set nb112nf_rinvH2O, 672
.set nb112nf_rinvH2H1, 688
.set nb112nf_rinvH2H2, 704
.set nb112nf_is3, 720
.set nb112nf_ii3, 724
.set nb112nf_nri, 728
.set nb112nf_iinr, 736
.set nb112nf_jindex, 744
.set nb112nf_jjnr, 752
.set nb112nf_shift, 760
.set nb112nf_shiftvec, 768
.set nb112nf_facel, 776
.set nb112nf_innerjjnr, 784
.set nb112nf_innerk, 792
.set nb112nf_n, 796
.set nb112nf_nn1, 800
.set nb112nf_nouter, 804
.set nb112nf_ninner, 808
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $824,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb112nf_nouter(%rsp)
        movl %eax,nb112nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb112nf_nri(%rsp)
        movq %rsi,nb112nf_iinr(%rsp)
        movq %rdx,nb112nf_jindex(%rsp)
        movq %rcx,nb112nf_jjnr(%rsp)
        movq %r8,nb112nf_shift(%rsp)
        movq %r9,nb112nf_shiftvec(%rsp)
        movq nb112nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb112nf_facel(%rsp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb112nf_half(%rsp)
        movl %ebx,nb112nf_half+4(%rsp)
        movsd nb112nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb112nf_half(%rsp)
        movapd %xmm3,nb112nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb112nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb112nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%rdx,%rbx,8),%xmm5
        movq nb112nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb112nf_facel(%rsp),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb112nf_qqOO(%rsp)
        movapd %xmm4,nb112nf_qqOH(%rsp)
        movapd %xmm5,nb112nf_qqHH(%rsp)

        xorpd %xmm0,%xmm0
        movq  nb112nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb112nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  nb112nf_vdwparam(%rbp),%rax
        movlpd (%rax,%rdx,8),%xmm0
        movlpd 8(%rax,%rdx,8),%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb112nf_c6(%rsp)
        movapd %xmm1,nb112nf_c12(%rsp)

_nb_kernel112nf_x86_64_sse2.nb112nf_threadloop: 
        movq  nb112nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel112nf_x86_64_sse2.nb112nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel112nf_x86_64_sse2.nb112nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb112nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb112nf_n(%rsp)
        movl %ebx,nb112nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel112nf_x86_64_sse2.nb112nf_outerstart
        jmp _nb_kernel112nf_x86_64_sse2.nb112nf_end

_nb_kernel112nf_x86_64_sse2.nb112nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb112nf_nouter(%rsp),%ebx
        movl %ebx,nb112nf_nouter(%rsp)

_nb_kernel112nf_x86_64_sse2.nb112nf_outer: 
        movq  nb112nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb112nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb112nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb112nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb112nf_ii3(%rsp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%rax,%rbx,8),%xmm3
        addsd 8(%rax,%rbx,8),%xmm4
        addsd 16(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb112nf_ixO(%rsp)
        movapd %xmm4,nb112nf_iyO(%rsp)
        movapd %xmm5,nb112nf_izO(%rsp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%rax,%rbx,8),%xmm0
        addsd 32(%rax,%rbx,8),%xmm1
        addsd 40(%rax,%rbx,8),%xmm2
        addsd 48(%rax,%rbx,8),%xmm3
        addsd 56(%rax,%rbx,8),%xmm4
        addsd 64(%rax,%rbx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb112nf_ixH1(%rsp)
        movapd %xmm1,nb112nf_iyH1(%rsp)
        movapd %xmm2,nb112nf_izH1(%rsp)
        movapd %xmm3,nb112nf_ixH2(%rsp)
        movapd %xmm4,nb112nf_iyH2(%rsp)
        movapd %xmm5,nb112nf_izH2(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb112nf_vctot(%rsp)
        movapd %xmm4,nb112nf_Vvdwtot(%rsp)

        movq  nb112nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb112nf_pos(%rbp),%rsi
        movq  nb112nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb112nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb112nf_ninner(%rsp),%ecx
        movl  %ecx,nb112nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb112nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel112nf_x86_64_sse2.nb112nf_unroll_loop
        jmp   _nb_kernel112nf_x86_64_sse2.nb112nf_checksingle
_nb_kernel112nf_x86_64_sse2.nb112nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb112nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb112nf_innerjjnr(%rsp)                 ## advance pointer (unrolled 2) 

        movq nb112nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move j coordinates to local temp variables 
        movlpd (%rsi,%rax,8),%xmm2
        movlpd 8(%rsi,%rax,8),%xmm3
        movlpd 16(%rsi,%rax,8),%xmm4
        movlpd 24(%rsi,%rax,8),%xmm5
        movlpd 32(%rsi,%rax,8),%xmm6
        movlpd 40(%rsi,%rax,8),%xmm7
        movhpd (%rsi,%rbx,8),%xmm2
        movhpd 8(%rsi,%rbx,8),%xmm3
        movhpd 16(%rsi,%rbx,8),%xmm4
        movhpd 24(%rsi,%rbx,8),%xmm5
        movhpd 32(%rsi,%rbx,8),%xmm6
        movhpd 40(%rsi,%rbx,8),%xmm7
        movapd  %xmm2,nb112nf_jxO(%rsp)
        movapd  %xmm3,nb112nf_jyO(%rsp)
        movapd  %xmm4,nb112nf_jzO(%rsp)
        movapd  %xmm5,nb112nf_jxH1(%rsp)
        movapd  %xmm6,nb112nf_jyH1(%rsp)
        movapd  %xmm7,nb112nf_jzH1(%rsp)
        movlpd 48(%rsi,%rax,8),%xmm2
        movlpd 56(%rsi,%rax,8),%xmm3
        movlpd 64(%rsi,%rax,8),%xmm4
        movhpd 48(%rsi,%rbx,8),%xmm2
        movhpd 56(%rsi,%rbx,8),%xmm3
        movhpd 64(%rsi,%rbx,8),%xmm4
        movapd  %xmm2,nb112nf_jxH2(%rsp)
        movapd  %xmm3,nb112nf_jyH2(%rsp)
        movapd  %xmm4,nb112nf_jzH2(%rsp)

        movapd nb112nf_ixO(%rsp),%xmm0
        movapd nb112nf_iyO(%rsp),%xmm1
        movapd nb112nf_izO(%rsp),%xmm2
        movapd nb112nf_ixO(%rsp),%xmm3
        movapd nb112nf_iyO(%rsp),%xmm4
        movapd nb112nf_izO(%rsp),%xmm5
        subpd  nb112nf_jxO(%rsp),%xmm0
        subpd  nb112nf_jyO(%rsp),%xmm1
        subpd  nb112nf_jzO(%rsp),%xmm2
        subpd  nb112nf_jxH1(%rsp),%xmm3
        subpd  nb112nf_jyH1(%rsp),%xmm4
        subpd  nb112nf_jzH1(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb112nf_rsqOO(%rsp)
        movapd %xmm3,nb112nf_rsqOH1(%rsp)

        movapd nb112nf_ixO(%rsp),%xmm0
        movapd nb112nf_iyO(%rsp),%xmm1
        movapd nb112nf_izO(%rsp),%xmm2
        movapd nb112nf_ixH1(%rsp),%xmm3
        movapd nb112nf_iyH1(%rsp),%xmm4
        movapd nb112nf_izH1(%rsp),%xmm5
        subpd  nb112nf_jxH2(%rsp),%xmm0
        subpd  nb112nf_jyH2(%rsp),%xmm1
        subpd  nb112nf_jzH2(%rsp),%xmm2
        subpd  nb112nf_jxO(%rsp),%xmm3
        subpd  nb112nf_jyO(%rsp),%xmm4
        subpd  nb112nf_jzO(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb112nf_rsqOH2(%rsp)
        movapd %xmm3,nb112nf_rsqH1O(%rsp)

        movapd nb112nf_ixH1(%rsp),%xmm0
        movapd nb112nf_iyH1(%rsp),%xmm1
        movapd nb112nf_izH1(%rsp),%xmm2
        movapd nb112nf_ixH1(%rsp),%xmm3
        movapd nb112nf_iyH1(%rsp),%xmm4
        movapd nb112nf_izH1(%rsp),%xmm5
        subpd  nb112nf_jxH1(%rsp),%xmm0
        subpd  nb112nf_jyH1(%rsp),%xmm1
        subpd  nb112nf_jzH1(%rsp),%xmm2
        subpd  nb112nf_jxH2(%rsp),%xmm3
        subpd  nb112nf_jyH2(%rsp),%xmm4
        subpd  nb112nf_jzH2(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb112nf_rsqH1H1(%rsp)
        movapd %xmm3,nb112nf_rsqH1H2(%rsp)

        movapd nb112nf_ixH2(%rsp),%xmm0
        movapd nb112nf_iyH2(%rsp),%xmm1
        movapd nb112nf_izH2(%rsp),%xmm2
        movapd nb112nf_ixH2(%rsp),%xmm3
        movapd nb112nf_iyH2(%rsp),%xmm4
        movapd nb112nf_izH2(%rsp),%xmm5
        subpd  nb112nf_jxO(%rsp),%xmm0
        subpd  nb112nf_jyO(%rsp),%xmm1
        subpd  nb112nf_jzO(%rsp),%xmm2
        subpd  nb112nf_jxH1(%rsp),%xmm3
        subpd  nb112nf_jyH1(%rsp),%xmm4
        subpd  nb112nf_jzH1(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb112nf_rsqH2O(%rsp)
        movapd %xmm4,nb112nf_rsqH2H1(%rsp)

        movapd nb112nf_ixH2(%rsp),%xmm0
        movapd nb112nf_iyH2(%rsp),%xmm1
        movapd nb112nf_izH2(%rsp),%xmm2
        subpd  nb112nf_jxH2(%rsp),%xmm0
        subpd  nb112nf_jyH2(%rsp),%xmm1
        subpd  nb112nf_jzH2(%rsp),%xmm2
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb112nf_rsqH2H2(%rsp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb112nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112nf_half(%rsp),%xmm3   ## iter1 
        mulpd   nb112nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb112nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb112nf_rinvH2H2(%rsp)
        movapd %xmm5,nb112nf_rinvH2H1(%rsp)

        movapd nb112nf_rsqOO(%rsp),%xmm0
        movapd nb112nf_rsqOH1(%rsp),%xmm4
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb112nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112nf_half(%rsp),%xmm3   ## iter1 of  
        mulpd   nb112nf_half(%rsp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb112nf_half(%rsp),%xmm5   ## rinv
        movapd %xmm1,nb112nf_rinvOO(%rsp)
        movapd %xmm5,nb112nf_rinvOH1(%rsp)

        movapd nb112nf_rsqOH2(%rsp),%xmm0
        movapd nb112nf_rsqH1O(%rsp),%xmm4
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb112nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112nf_half(%rsp),%xmm3   ## iter1 
        mulpd   nb112nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb112nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb112nf_rinvOH2(%rsp)
        movapd %xmm5,nb112nf_rinvH1O(%rsp)

        movapd nb112nf_rsqH1H1(%rsp),%xmm0
        movapd nb112nf_rsqH1H2(%rsp),%xmm4
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb112nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112nf_half(%rsp),%xmm3   ## iter1a 
        mulpd   nb112nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb112nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb112nf_rinvH1H1(%rsp)
        movapd %xmm5,nb112nf_rinvH1H2(%rsp)

        movapd nb112nf_rsqH2O(%rsp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb112nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb112nf_half(%rsp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb112nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb112nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb112nf_rinvH2O(%rsp)

        ## start with OO interaction 
        movapd nb112nf_rinvOO(%rsp),%xmm0
        movapd %xmm0,%xmm7
        mulpd  %xmm0,%xmm0
        movapd %xmm0,%xmm1
        mulpd  %xmm0,%xmm1
        mulpd  %xmm0,%xmm1      ## xmm1=rinvsix 
        mulpd  nb112nf_qqOO(%rsp),%xmm7
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulpd  nb112nf_c6(%rsp),%xmm1
        mulpd  nb112nf_c12(%rsp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addpd  nb112nf_Vvdwtot(%rsp),%xmm3
        movapd %xmm3,nb112nf_Vvdwtot(%rsp)
        addpd  nb112nf_vctot(%rsp),%xmm7

        ## other interactions 
        movapd nb112nf_rinvOH1(%rsp),%xmm1
        movapd nb112nf_rinvH1H1(%rsp),%xmm2

        addpd nb112nf_rinvOH2(%rsp),%xmm1
        addpd nb112nf_rinvH1H2(%rsp),%xmm2

        addpd nb112nf_rinvH1O(%rsp),%xmm1
        addpd nb112nf_rinvH2H1(%rsp),%xmm2

        addpd nb112nf_rinvH2O(%rsp),%xmm1
        addpd nb112nf_rinvH2H2(%rsp),%xmm2

        mulpd nb112nf_qqOH(%rsp),%xmm1
        mulpd nb112nf_qqHH(%rsp),%xmm2

        addpd %xmm1,%xmm7
        addpd %xmm2,%xmm7

        movapd %xmm7,nb112nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb112nf_innerk(%rsp)
        jl    _nb_kernel112nf_x86_64_sse2.nb112nf_checksingle
        jmp   _nb_kernel112nf_x86_64_sse2.nb112nf_unroll_loop
_nb_kernel112nf_x86_64_sse2.nb112nf_checksingle: 
        movl  nb112nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel112nf_x86_64_sse2.nb112nf_dosingle
        jmp   _nb_kernel112nf_x86_64_sse2.nb112nf_updateouterdata
_nb_kernel112nf_x86_64_sse2.nb112nf_dosingle: 
        movq  nb112nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb112nf_innerjjnr(%rsp)

        movq nb112nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## fetch j coordinates 
        movlpd (%rsi,%rax,8),%xmm2
        movlpd 8(%rsi,%rax,8),%xmm3
        movlpd 16(%rsi,%rax,8),%xmm4
        movlpd 24(%rsi,%rax,8),%xmm5
        movlpd 32(%rsi,%rax,8),%xmm6
        movlpd 40(%rsi,%rax,8),%xmm7
        movapd  %xmm2,nb112nf_jxO(%rsp)
        movapd  %xmm3,nb112nf_jyO(%rsp)
        movapd  %xmm4,nb112nf_jzO(%rsp)
        movapd  %xmm5,nb112nf_jxH1(%rsp)
        movapd  %xmm6,nb112nf_jyH1(%rsp)
        movapd  %xmm7,nb112nf_jzH1(%rsp)
        movlpd 48(%rsi,%rax,8),%xmm2
        movlpd 56(%rsi,%rax,8),%xmm3
        movlpd 64(%rsi,%rax,8),%xmm4
        movapd  %xmm2,nb112nf_jxH2(%rsp)
        movapd  %xmm3,nb112nf_jyH2(%rsp)
        movapd  %xmm4,nb112nf_jzH2(%rsp)

        movapd nb112nf_ixO(%rsp),%xmm0
        movapd nb112nf_iyO(%rsp),%xmm1
        movapd nb112nf_izO(%rsp),%xmm2
        movapd nb112nf_ixO(%rsp),%xmm3
        movapd nb112nf_iyO(%rsp),%xmm4
        movapd nb112nf_izO(%rsp),%xmm5
        subsd  nb112nf_jxO(%rsp),%xmm0
        subsd  nb112nf_jyO(%rsp),%xmm1
        subsd  nb112nf_jzO(%rsp),%xmm2
        subsd  nb112nf_jxH1(%rsp),%xmm3
        subsd  nb112nf_jyH1(%rsp),%xmm4
        subsd  nb112nf_jzH1(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb112nf_rsqOO(%rsp)
        movapd %xmm3,nb112nf_rsqOH1(%rsp)

        movapd nb112nf_ixO(%rsp),%xmm0
        movapd nb112nf_iyO(%rsp),%xmm1
        movapd nb112nf_izO(%rsp),%xmm2
        movapd nb112nf_ixH1(%rsp),%xmm3
        movapd nb112nf_iyH1(%rsp),%xmm4
        movapd nb112nf_izH1(%rsp),%xmm5
        subsd  nb112nf_jxH2(%rsp),%xmm0
        subsd  nb112nf_jyH2(%rsp),%xmm1
        subsd  nb112nf_jzH2(%rsp),%xmm2
        subsd  nb112nf_jxO(%rsp),%xmm3
        subsd  nb112nf_jyO(%rsp),%xmm4
        subsd  nb112nf_jzO(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb112nf_rsqOH2(%rsp)
        movapd %xmm3,nb112nf_rsqH1O(%rsp)

        movapd nb112nf_ixH1(%rsp),%xmm0
        movapd nb112nf_iyH1(%rsp),%xmm1
        movapd nb112nf_izH1(%rsp),%xmm2
        movapd nb112nf_ixH1(%rsp),%xmm3
        movapd nb112nf_iyH1(%rsp),%xmm4
        movapd nb112nf_izH1(%rsp),%xmm5
        subsd  nb112nf_jxH1(%rsp),%xmm0
        subsd  nb112nf_jyH1(%rsp),%xmm1
        subsd  nb112nf_jzH1(%rsp),%xmm2
        subsd  nb112nf_jxH2(%rsp),%xmm3
        subsd  nb112nf_jyH2(%rsp),%xmm4
        subsd  nb112nf_jzH2(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb112nf_rsqH1H1(%rsp)
        movapd %xmm3,nb112nf_rsqH1H2(%rsp)

        movapd nb112nf_ixH2(%rsp),%xmm0
        movapd nb112nf_iyH2(%rsp),%xmm1
        movapd nb112nf_izH2(%rsp),%xmm2
        movapd nb112nf_ixH2(%rsp),%xmm3
        movapd nb112nf_iyH2(%rsp),%xmm4
        movapd nb112nf_izH2(%rsp),%xmm5
        subsd  nb112nf_jxO(%rsp),%xmm0
        subsd  nb112nf_jyO(%rsp),%xmm1
        subsd  nb112nf_jzO(%rsp),%xmm2
        subsd  nb112nf_jxH1(%rsp),%xmm3
        subsd  nb112nf_jyH1(%rsp),%xmm4
        subsd  nb112nf_jzH1(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb112nf_rsqH2O(%rsp)
        movapd %xmm4,nb112nf_rsqH2H1(%rsp)

        movapd nb112nf_ixH2(%rsp),%xmm0
        movapd nb112nf_iyH2(%rsp),%xmm1
        movapd nb112nf_izH2(%rsp),%xmm2
        subsd  nb112nf_jxH2(%rsp),%xmm0
        subsd  nb112nf_jyH2(%rsp),%xmm1
        subsd  nb112nf_jzH2(%rsp),%xmm2
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb112nf_rsqH2H2(%rsp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb112nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112nf_half(%rsp),%xmm3   ## iter1 
        mulsd   nb112nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb112nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb112nf_rinvH2H2(%rsp)
        movapd %xmm5,nb112nf_rinvH2H1(%rsp)

        movapd nb112nf_rsqOO(%rsp),%xmm0
        movapd nb112nf_rsqOH1(%rsp),%xmm4
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb112nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112nf_half(%rsp),%xmm3   ## iter1 of  
        mulsd   nb112nf_half(%rsp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb112nf_half(%rsp),%xmm5   ## rinv
        movapd %xmm1,nb112nf_rinvOO(%rsp)
        movapd %xmm5,nb112nf_rinvOH1(%rsp)

        movapd nb112nf_rsqOH2(%rsp),%xmm0
        movapd nb112nf_rsqH1O(%rsp),%xmm4
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb112nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112nf_half(%rsp),%xmm3   ## iter1 
        mulsd   nb112nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb112nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb112nf_rinvOH2(%rsp)
        movapd %xmm5,nb112nf_rinvH1O(%rsp)

        movapd nb112nf_rsqH1H1(%rsp),%xmm0
        movapd nb112nf_rsqH1H2(%rsp),%xmm4
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb112nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112nf_half(%rsp),%xmm3   ## iter1a 
        mulsd   nb112nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb112nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb112nf_rinvH1H1(%rsp)
        movapd %xmm5,nb112nf_rinvH1H2(%rsp)

        movapd nb112nf_rsqH2O(%rsp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb112nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb112nf_half(%rsp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb112nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb112nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb112nf_rinvH2O(%rsp)

        ## start with OO interaction 
        movapd nb112nf_rinvOO(%rsp),%xmm0
        movapd %xmm0,%xmm7
        mulsd  %xmm0,%xmm0
        movapd %xmm0,%xmm1
        mulsd  %xmm0,%xmm1
        mulsd  %xmm0,%xmm1      ## xmm1=rinvsix 
        mulsd  nb112nf_qqOO(%rsp),%xmm7
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  nb112nf_c6(%rsp),%xmm1
        mulsd  nb112nf_c12(%rsp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addsd  nb112nf_Vvdwtot(%rsp),%xmm3
        movlpd %xmm3,nb112nf_Vvdwtot(%rsp)
        addsd  nb112nf_vctot(%rsp),%xmm7

        ## other interactions 
        movapd nb112nf_rinvOH1(%rsp),%xmm1
        movapd nb112nf_rinvH1H1(%rsp),%xmm2

        addsd nb112nf_rinvOH2(%rsp),%xmm1
        addsd nb112nf_rinvH1H2(%rsp),%xmm2

        addsd nb112nf_rinvH1O(%rsp),%xmm1
        addsd nb112nf_rinvH2H1(%rsp),%xmm2

        addsd nb112nf_rinvH2O(%rsp),%xmm1
        addsd nb112nf_rinvH2H2(%rsp),%xmm2

        mulsd nb112nf_qqOH(%rsp),%xmm1
        mulsd nb112nf_qqHH(%rsp),%xmm2

        addsd %xmm1,%xmm7
        addsd %xmm2,%xmm7

        movlpd %xmm7,nb112nf_vctot(%rsp)

_nb_kernel112nf_x86_64_sse2.nb112nf_updateouterdata: 
        ## get n from stack
        movl nb112nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb112nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb112nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb112nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb112nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb112nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb112nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel112nf_x86_64_sse2.nb112nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb112nf_n(%rsp)
        jmp _nb_kernel112nf_x86_64_sse2.nb112nf_outer
_nb_kernel112nf_x86_64_sse2.nb112nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb112nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel112nf_x86_64_sse2.nb112nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel112nf_x86_64_sse2.nb112nf_threadloop
_nb_kernel112nf_x86_64_sse2.nb112nf_end: 
        movl nb112nf_nouter(%rsp),%eax
        movl nb112nf_ninner(%rsp),%ebx
        movq nb112nf_outeriter(%rbp),%rcx
        movq nb112nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $824,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret



