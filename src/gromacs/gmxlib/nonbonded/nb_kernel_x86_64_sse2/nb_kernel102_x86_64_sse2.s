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






.globl nb_kernel102_x86_64_sse2
.globl _nb_kernel102_x86_64_sse2
nb_kernel102_x86_64_sse2:       
_nb_kernel102_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb102_fshift, 16
.set nb102_gid, 24
.set nb102_pos, 32
.set nb102_faction, 40
.set nb102_charge, 48
.set nb102_p_facel, 56
.set nb102_argkrf, 64
.set nb102_argcrf, 72
.set nb102_Vc, 80
.set nb102_type, 88
.set nb102_p_ntype, 96
.set nb102_vdwparam, 104
.set nb102_Vvdw, 112
.set nb102_p_tabscale, 120
.set nb102_VFtab, 128
.set nb102_invsqrta, 136
.set nb102_dvda, 144
.set nb102_p_gbtabscale, 152
.set nb102_GBtab, 160
.set nb102_p_nthreads, 168
.set nb102_count, 176
.set nb102_mtx, 184
.set nb102_outeriter, 192
.set nb102_inneriter, 200
.set nb102_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use        
.set nb102_ixO, 0
.set nb102_iyO, 16
.set nb102_izO, 32
.set nb102_ixH1, 48
.set nb102_iyH1, 64
.set nb102_izH1, 80
.set nb102_ixH2, 96
.set nb102_iyH2, 112
.set nb102_izH2, 128
.set nb102_jxO, 144
.set nb102_jyO, 160
.set nb102_jzO, 176
.set nb102_jxH1, 192
.set nb102_jyH1, 208
.set nb102_jzH1, 224
.set nb102_jxH2, 240
.set nb102_jyH2, 256
.set nb102_jzH2, 272
.set nb102_dxOO, 288
.set nb102_dyOO, 304
.set nb102_dzOO, 320
.set nb102_dxOH1, 336
.set nb102_dyOH1, 352
.set nb102_dzOH1, 368
.set nb102_dxOH2, 384
.set nb102_dyOH2, 400
.set nb102_dzOH2, 416
.set nb102_dxH1O, 432
.set nb102_dyH1O, 448
.set nb102_dzH1O, 464
.set nb102_dxH1H1, 480
.set nb102_dyH1H1, 496
.set nb102_dzH1H1, 512
.set nb102_dxH1H2, 528
.set nb102_dyH1H2, 544
.set nb102_dzH1H2, 560
.set nb102_dxH2O, 576
.set nb102_dyH2O, 592
.set nb102_dzH2O, 608
.set nb102_dxH2H1, 624
.set nb102_dyH2H1, 640
.set nb102_dzH2H1, 656
.set nb102_dxH2H2, 672
.set nb102_dyH2H2, 688
.set nb102_dzH2H2, 704
.set nb102_qqOO, 720
.set nb102_qqOH, 736
.set nb102_qqHH, 752
.set nb102_vctot, 768
.set nb102_fixO, 784
.set nb102_fiyO, 800
.set nb102_fizO, 816
.set nb102_fixH1, 832
.set nb102_fiyH1, 848
.set nb102_fizH1, 864
.set nb102_fixH2, 880
.set nb102_fiyH2, 896
.set nb102_fizH2, 912
.set nb102_fjxO, 928
.set nb102_fjyO, 944
.set nb102_fjzO, 960
.set nb102_fjxH1, 976
.set nb102_fjyH1, 992
.set nb102_fjzH1, 1008
.set nb102_fjxH2, 1024
.set nb102_fjyH2, 1040
.set nb102_fjzH2, 1056
.set nb102_half, 1072
.set nb102_three, 1088
.set nb102_rsqOO, 1104
.set nb102_rsqOH1, 1120
.set nb102_rsqOH2, 1136
.set nb102_rsqH1O, 1152
.set nb102_rsqH1H1, 1168
.set nb102_rsqH1H2, 1184
.set nb102_rsqH2O, 1200
.set nb102_rsqH2H1, 1216
.set nb102_rsqH2H2, 1232
.set nb102_rinvOO, 1248
.set nb102_rinvOH1, 1264
.set nb102_rinvOH2, 1280
.set nb102_rinvH1O, 1296
.set nb102_rinvH1H1, 1312
.set nb102_rinvH1H2, 1328
.set nb102_rinvH2O, 1344
.set nb102_rinvH2H1, 1360
.set nb102_rinvH2H2, 1376
.set nb102_is3, 1392
.set nb102_ii3, 1396
.set nb102_innerjjnr, 1400
.set nb102_nri, 1408
.set nb102_iinr, 1416
.set nb102_jindex, 1424
.set nb102_jjnr, 1432
.set nb102_shift, 1440
.set nb102_shiftvec, 1448
.set nb102_facel, 1456
.set nb102_innerk, 1464
.set nb102_n, 1468
.set nb102_nn1, 1472
.set nb102_nouter, 1476
.set nb102_ninner, 1480
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1496,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb102_nouter(%rsp)
        movl %eax,nb102_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb102_nri(%rsp)
        movq %rsi,nb102_iinr(%rsp)
        movq %rdx,nb102_jindex(%rsp)
        movq %rcx,nb102_jjnr(%rsp)
        movq %r8,nb102_shift(%rsp)
        movq %r9,nb102_shiftvec(%rsp)
        movq nb102_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb102_facel(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb102_half(%rsp)
        movl %ebx,nb102_half+4(%rsp)
        movsd nb102_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb102_half(%rsp)
        movapd %xmm3,nb102_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb102_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb102_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3       ## qO 
        movsd %xmm3,%xmm4               ## qO 
        movsd 8(%rdx,%rbx,8),%xmm5      ## qH 
        movq nb102_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb102_facel(%rsp),%xmm6   ## facel 
        mulsd  %xmm3,%xmm3              ## qO*qO 
        mulsd  %xmm5,%xmm4              ## qO*qH 
        mulsd  %xmm5,%xmm5              ## qH*qH 
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb102_qqOO(%rsp)
        movapd %xmm4,nb102_qqOH(%rsp)
        movapd %xmm5,nb102_qqHH(%rsp)

_nb_kernel102_x86_64_sse2.nb102_threadloop: 
        movq  nb102_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel102_x86_64_sse2.nb102_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel102_x86_64_sse2.nb102_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb102_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb102_n(%rsp)
        movl %ebx,nb102_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel102_x86_64_sse2.nb102_outerstart
        jmp _nb_kernel102_x86_64_sse2.nb102_end

_nb_kernel102_x86_64_sse2.nb102_outerstart: 
        ## ebx contains number of outer iterations
        addl nb102_nouter(%rsp),%ebx
        movl %ebx,nb102_nouter(%rsp)

_nb_kernel102_x86_64_sse2.nb102_outer: 
        movq  nb102_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb102_is3(%rsp)      ## store is3 

        movq  nb102_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb102_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb102_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb102_ii3(%rsp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%rax,%rbx,8),%xmm3
        addsd 8(%rax,%rbx,8),%xmm4
        addsd 16(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb102_ixO(%rsp)
        movapd %xmm4,nb102_iyO(%rsp)
        movapd %xmm5,nb102_izO(%rsp)

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
        movapd %xmm0,nb102_ixH1(%rsp)
        movapd %xmm1,nb102_iyH1(%rsp)
        movapd %xmm2,nb102_izH1(%rsp)
        movapd %xmm3,nb102_ixH2(%rsp)
        movapd %xmm4,nb102_iyH2(%rsp)
        movapd %xmm5,nb102_izH2(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb102_vctot(%rsp)
        movapd %xmm4,nb102_fixO(%rsp)
        movapd %xmm4,nb102_fiyO(%rsp)
        movapd %xmm4,nb102_fizO(%rsp)
        movapd %xmm4,nb102_fixH1(%rsp)
        movapd %xmm4,nb102_fiyH1(%rsp)
        movapd %xmm4,nb102_fizH1(%rsp)
        movapd %xmm4,nb102_fixH2(%rsp)
        movapd %xmm4,nb102_fiyH2(%rsp)
        movapd %xmm4,nb102_fizH2(%rsp)

        movq  nb102_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx     ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb102_pos(%rbp),%rsi
        movq  nb102_faction(%rbp),%rdi
        movq  nb102_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb102_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb102_ninner(%rsp),%ecx
        movl  %ecx,nb102_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb102_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel102_x86_64_sse2.nb102_unroll_loop
        jmp   _nb_kernel102_x86_64_sse2.nb102_checksingle
_nb_kernel102_x86_64_sse2.nb102_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb102_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb102_innerjjnr(%rsp)            ## advance pointer (unrolled 2) 

        movq nb102_pos(%rbp),%rsi        ## base of pos[] 

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

    subpd nb102_ixO(%rsp),%xmm0
    subpd nb102_iyO(%rsp),%xmm1
    subpd nb102_izO(%rsp),%xmm2
    subpd nb102_ixH1(%rsp),%xmm3
    subpd nb102_iyH1(%rsp),%xmm4
    subpd nb102_izH1(%rsp),%xmm5
    subpd nb102_ixH2(%rsp),%xmm6
    subpd nb102_iyH2(%rsp),%xmm7
    subpd nb102_izH2(%rsp),%xmm8

        movapd %xmm0,nb102_dxOO(%rsp)
        movapd %xmm1,nb102_dyOO(%rsp)
        movapd %xmm2,nb102_dzOO(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb102_dxH1O(%rsp)
        movapd %xmm4,nb102_dyH1O(%rsp)
        movapd %xmm5,nb102_dzH1O(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb102_dxH2O(%rsp)
        movapd %xmm7,nb102_dyH2O(%rsp)
        movapd %xmm8,nb102_dzH2O(%rsp)
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

        movapd  nb102_three(%rsp),%xmm9
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

        movapd  nb102_half(%rsp),%xmm15
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

        movapd  nb102_three(%rsp),%xmm1
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

        movapd  nb102_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvOO 
        mulpd   %xmm15,%xmm10 ##   rinvH1O
    mulpd   %xmm15,%xmm11 ##   rinvH2O

        ## O interactions 
    movapd %xmm9,%xmm0
    movapd %xmm10,%xmm1
    movapd %xmm11,%xmm2
    mulpd  %xmm9,%xmm9
    mulpd  %xmm10,%xmm10
    mulpd  %xmm11,%xmm11
    mulpd  nb102_qqOO(%rsp),%xmm0
    mulpd  nb102_qqOH(%rsp),%xmm1
    mulpd  nb102_qqOH(%rsp),%xmm2
    mulpd  %xmm0,%xmm9
    mulpd  %xmm1,%xmm10
    mulpd  %xmm2,%xmm11

    addpd nb102_vctot(%rsp),%xmm0
    addpd %xmm2,%xmm1
    addpd %xmm1,%xmm0
    movapd %xmm0,nb102_vctot(%rsp)

    ## move j O forces to xmm0-xmm2
        movq  nb102_faction(%rbp),%rdi
        movlpd (%rdi,%rax,8),%xmm0
        movlpd 8(%rdi,%rax,8),%xmm1
        movlpd 16(%rdi,%rax,8),%xmm2
        movhpd (%rdi,%rbx,8),%xmm0
        movhpd 8(%rdi,%rbx,8),%xmm1
        movhpd 16(%rdi,%rbx,8),%xmm2

    movapd %xmm9,%xmm7
    movapd %xmm9,%xmm8
    movapd %xmm11,%xmm13
    movapd %xmm11,%xmm14
    movapd %xmm11,%xmm15
    movapd %xmm10,%xmm11
    movapd %xmm10,%xmm12

        mulpd nb102_dxOO(%rsp),%xmm7
        mulpd nb102_dyOO(%rsp),%xmm8
        mulpd nb102_dzOO(%rsp),%xmm9
        mulpd nb102_dxH1O(%rsp),%xmm10
        mulpd nb102_dyH1O(%rsp),%xmm11
        mulpd nb102_dzH1O(%rsp),%xmm12
        mulpd nb102_dxH2O(%rsp),%xmm13
        mulpd nb102_dyH2O(%rsp),%xmm14
        mulpd nb102_dzH2O(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb102_fixO(%rsp),%xmm7
    addpd nb102_fiyO(%rsp),%xmm8
    addpd nb102_fizO(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb102_fixH1(%rsp),%xmm10
    addpd nb102_fiyH1(%rsp),%xmm11
    addpd nb102_fizH1(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb102_fixH2(%rsp),%xmm13
    addpd nb102_fiyH2(%rsp),%xmm14
    addpd nb102_fizH2(%rsp),%xmm15

    movapd %xmm7,nb102_fixO(%rsp)
    movapd %xmm8,nb102_fiyO(%rsp)
    movapd %xmm9,nb102_fizO(%rsp)
    movapd %xmm10,nb102_fixH1(%rsp)
    movapd %xmm11,nb102_fiyH1(%rsp)
    movapd %xmm12,nb102_fizH1(%rsp)
    movapd %xmm13,nb102_fixH2(%rsp)
    movapd %xmm14,nb102_fiyH2(%rsp)
    movapd %xmm15,nb102_fizH2(%rsp)

    ## store back j O forces from xmm0-xmm2
        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)
        movhpd %xmm0,(%rdi,%rbx,8)
        movhpd %xmm1,8(%rdi,%rbx,8)
        movhpd %xmm2,16(%rdi,%rbx,8)

        ## move j H1 coordinates to local temp variables 
        movq  nb102_pos(%rbp),%rsi
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

    subpd nb102_ixO(%rsp),%xmm0
    subpd nb102_iyO(%rsp),%xmm1
    subpd nb102_izO(%rsp),%xmm2
    subpd nb102_ixH1(%rsp),%xmm3
    subpd nb102_iyH1(%rsp),%xmm4
    subpd nb102_izH1(%rsp),%xmm5
    subpd nb102_ixH2(%rsp),%xmm6
    subpd nb102_iyH2(%rsp),%xmm7
    subpd nb102_izH2(%rsp),%xmm8

        movapd %xmm0,nb102_dxOH1(%rsp)
        movapd %xmm1,nb102_dyOH1(%rsp)
        movapd %xmm2,nb102_dzOH1(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb102_dxH1H1(%rsp)
        movapd %xmm4,nb102_dyH1H1(%rsp)
        movapd %xmm5,nb102_dzH1H1(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb102_dxH2H1(%rsp)
        movapd %xmm7,nb102_dyH2H1(%rsp)
        movapd %xmm8,nb102_dzH2H1(%rsp)
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

        movapd  nb102_three(%rsp),%xmm9
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

        movapd  nb102_half(%rsp),%xmm15
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

        movapd  nb102_three(%rsp),%xmm1
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

        movapd  nb102_half(%rsp),%xmm15
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
    mulpd  nb102_qqOH(%rsp),%xmm0
    mulpd  nb102_qqHH(%rsp),%xmm1
    mulpd  nb102_qqHH(%rsp),%xmm2
    mulpd  %xmm0,%xmm9
    mulpd  %xmm1,%xmm10
    mulpd  %xmm2,%xmm11

    addpd nb102_vctot(%rsp),%xmm0
    addpd %xmm2,%xmm1
    addpd %xmm1,%xmm0
    movapd %xmm0,nb102_vctot(%rsp)

    ## move j H1 forces to xmm0-xmm2
        movq  nb102_faction(%rbp),%rdi
        movlpd 24(%rdi,%rax,8),%xmm0
        movlpd 32(%rdi,%rax,8),%xmm1
        movlpd 40(%rdi,%rax,8),%xmm2
        movhpd 24(%rdi,%rbx,8),%xmm0
        movhpd 32(%rdi,%rbx,8),%xmm1
        movhpd 40(%rdi,%rbx,8),%xmm2

    movapd %xmm9,%xmm7
    movapd %xmm9,%xmm8
    movapd %xmm11,%xmm13
    movapd %xmm11,%xmm14
    movapd %xmm11,%xmm15
    movapd %xmm10,%xmm11
    movapd %xmm10,%xmm12

        mulpd nb102_dxOH1(%rsp),%xmm7
        mulpd nb102_dyOH1(%rsp),%xmm8
        mulpd nb102_dzOH1(%rsp),%xmm9
        mulpd nb102_dxH1H1(%rsp),%xmm10
        mulpd nb102_dyH1H1(%rsp),%xmm11
        mulpd nb102_dzH1H1(%rsp),%xmm12
        mulpd nb102_dxH2H1(%rsp),%xmm13
        mulpd nb102_dyH2H1(%rsp),%xmm14
        mulpd nb102_dzH2H1(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb102_fixO(%rsp),%xmm7
    addpd nb102_fiyO(%rsp),%xmm8
    addpd nb102_fizO(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb102_fixH1(%rsp),%xmm10
    addpd nb102_fiyH1(%rsp),%xmm11
    addpd nb102_fizH1(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb102_fixH2(%rsp),%xmm13
    addpd nb102_fiyH2(%rsp),%xmm14
    addpd nb102_fizH2(%rsp),%xmm15

    movapd %xmm7,nb102_fixO(%rsp)
    movapd %xmm8,nb102_fiyO(%rsp)
    movapd %xmm9,nb102_fizO(%rsp)
    movapd %xmm10,nb102_fixH1(%rsp)
    movapd %xmm11,nb102_fiyH1(%rsp)
    movapd %xmm12,nb102_fizH1(%rsp)
    movapd %xmm13,nb102_fixH2(%rsp)
    movapd %xmm14,nb102_fiyH2(%rsp)
    movapd %xmm15,nb102_fizH2(%rsp)

    ## store back j H1 forces from xmm0-xmm2
        movlpd %xmm0,24(%rdi,%rax,8)
        movlpd %xmm1,32(%rdi,%rax,8)
        movlpd %xmm2,40(%rdi,%rax,8)
        movhpd %xmm0,24(%rdi,%rbx,8)
        movhpd %xmm1,32(%rdi,%rbx,8)
        movhpd %xmm2,40(%rdi,%rbx,8)

        ## move j H2 coordinates to local temp variables 
        movq  nb102_pos(%rbp),%rsi
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

    subpd nb102_ixO(%rsp),%xmm0
    subpd nb102_iyO(%rsp),%xmm1
    subpd nb102_izO(%rsp),%xmm2
    subpd nb102_ixH1(%rsp),%xmm3
    subpd nb102_iyH1(%rsp),%xmm4
    subpd nb102_izH1(%rsp),%xmm5
    subpd nb102_ixH2(%rsp),%xmm6
    subpd nb102_iyH2(%rsp),%xmm7
    subpd nb102_izH2(%rsp),%xmm8

        movapd %xmm0,nb102_dxOH2(%rsp)
        movapd %xmm1,nb102_dyOH2(%rsp)
        movapd %xmm2,nb102_dzOH2(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb102_dxH1H2(%rsp)
        movapd %xmm4,nb102_dyH1H2(%rsp)
        movapd %xmm5,nb102_dzH1H2(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb102_dxH2H2(%rsp)
        movapd %xmm7,nb102_dyH2H2(%rsp)
        movapd %xmm8,nb102_dzH2H2(%rsp)
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

        movapd  nb102_three(%rsp),%xmm9
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

        movapd  nb102_half(%rsp),%xmm15
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

        movapd  nb102_three(%rsp),%xmm1
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

        movapd  nb102_half(%rsp),%xmm15
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
    mulpd  nb102_qqOH(%rsp),%xmm0
    mulpd  nb102_qqHH(%rsp),%xmm1
    mulpd  nb102_qqHH(%rsp),%xmm2
    mulpd  %xmm0,%xmm9
    mulpd  %xmm1,%xmm10
    mulpd  %xmm2,%xmm11

    addpd nb102_vctot(%rsp),%xmm0
    addpd %xmm2,%xmm1
    addpd %xmm1,%xmm0
    movapd %xmm0,nb102_vctot(%rsp)

    ## move j H2 forces to xmm0-xmm2
        movq  nb102_faction(%rbp),%rdi
        movlpd 48(%rdi,%rax,8),%xmm0
        movlpd 56(%rdi,%rax,8),%xmm1
        movlpd 64(%rdi,%rax,8),%xmm2
        movhpd 48(%rdi,%rbx,8),%xmm0
        movhpd 56(%rdi,%rbx,8),%xmm1
        movhpd 64(%rdi,%rbx,8),%xmm2

    movapd %xmm9,%xmm7
    movapd %xmm9,%xmm8
    movapd %xmm11,%xmm13
    movapd %xmm11,%xmm14
    movapd %xmm11,%xmm15
    movapd %xmm10,%xmm11
    movapd %xmm10,%xmm12

        mulpd nb102_dxOH2(%rsp),%xmm7
        mulpd nb102_dyOH2(%rsp),%xmm8
        mulpd nb102_dzOH2(%rsp),%xmm9
        mulpd nb102_dxH1H2(%rsp),%xmm10
        mulpd nb102_dyH1H2(%rsp),%xmm11
        mulpd nb102_dzH1H2(%rsp),%xmm12
        mulpd nb102_dxH2H2(%rsp),%xmm13
        mulpd nb102_dyH2H2(%rsp),%xmm14
        mulpd nb102_dzH2H2(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb102_fixO(%rsp),%xmm7
    addpd nb102_fiyO(%rsp),%xmm8
    addpd nb102_fizO(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb102_fixH1(%rsp),%xmm10
    addpd nb102_fiyH1(%rsp),%xmm11
    addpd nb102_fizH1(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb102_fixH2(%rsp),%xmm13
    addpd nb102_fiyH2(%rsp),%xmm14
    addpd nb102_fizH2(%rsp),%xmm15

    movapd %xmm7,nb102_fixO(%rsp)
    movapd %xmm8,nb102_fiyO(%rsp)
    movapd %xmm9,nb102_fizO(%rsp)
    movapd %xmm10,nb102_fixH1(%rsp)
    movapd %xmm11,nb102_fiyH1(%rsp)
    movapd %xmm12,nb102_fizH1(%rsp)
    movapd %xmm13,nb102_fixH2(%rsp)
    movapd %xmm14,nb102_fiyH2(%rsp)
    movapd %xmm15,nb102_fizH2(%rsp)

    ## store back j H2 forces from xmm0-xmm2
        movlpd %xmm0,48(%rdi,%rax,8)
        movlpd %xmm1,56(%rdi,%rax,8)
        movlpd %xmm2,64(%rdi,%rax,8)
        movhpd %xmm0,48(%rdi,%rbx,8)
        movhpd %xmm1,56(%rdi,%rbx,8)
        movhpd %xmm2,64(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb102_innerk(%rsp)
        jl    _nb_kernel102_x86_64_sse2.nb102_checksingle
        jmp   _nb_kernel102_x86_64_sse2.nb102_unroll_loop
_nb_kernel102_x86_64_sse2.nb102_checksingle: 
        movl  nb102_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel102_x86_64_sse2.nb102_dosingle
        jmp   _nb_kernel102_x86_64_sse2.nb102_updateouterdata
_nb_kernel102_x86_64_sse2.nb102_dosingle: 
        movq  nb102_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb102_pos(%rbp),%rsi
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

    subsd nb102_ixO(%rsp),%xmm0
    subsd nb102_iyO(%rsp),%xmm1
    subsd nb102_izO(%rsp),%xmm2
    subsd nb102_ixH1(%rsp),%xmm3
    subsd nb102_iyH1(%rsp),%xmm4
    subsd nb102_izH1(%rsp),%xmm5
    subsd nb102_ixH2(%rsp),%xmm6
    subsd nb102_iyH2(%rsp),%xmm7
    subsd nb102_izH2(%rsp),%xmm8

        movsd %xmm0,nb102_dxOO(%rsp)
        movsd %xmm1,nb102_dyOO(%rsp)
        movsd %xmm2,nb102_dzOO(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb102_dxH1O(%rsp)
        movsd %xmm4,nb102_dyH1O(%rsp)
        movsd %xmm5,nb102_dzH1O(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb102_dxH2O(%rsp)
        movsd %xmm7,nb102_dyH2O(%rsp)
        movsd %xmm8,nb102_dzH2O(%rsp)
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

        movsd  nb102_three(%rsp),%xmm9
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

        movsd  nb102_half(%rsp),%xmm15
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

        movsd  nb102_three(%rsp),%xmm1
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

        movsd  nb102_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvOO 
        mulsd   %xmm15,%xmm10 ##   rinvH1O
    mulsd   %xmm15,%xmm11 ##   rinvH2O

        ## O interactions 
    movsd %xmm9,%xmm0
    movsd %xmm10,%xmm1
    movsd %xmm11,%xmm2
    mulsd  %xmm9,%xmm9
    mulsd  %xmm10,%xmm10
    mulsd  %xmm11,%xmm11
    mulsd  nb102_qqOO(%rsp),%xmm0
    mulsd  nb102_qqOH(%rsp),%xmm1
    mulsd  nb102_qqOH(%rsp),%xmm2
    mulsd  %xmm0,%xmm9
    mulsd  %xmm1,%xmm10
    mulsd  %xmm2,%xmm11

    addsd nb102_vctot(%rsp),%xmm0
    addsd %xmm2,%xmm1
    addsd %xmm1,%xmm0
    movsd %xmm0,nb102_vctot(%rsp)

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

        mulsd nb102_dxOO(%rsp),%xmm7
        mulsd nb102_dyOO(%rsp),%xmm8
        mulsd nb102_dzOO(%rsp),%xmm9
        mulsd nb102_dxH1O(%rsp),%xmm10
        mulsd nb102_dyH1O(%rsp),%xmm11
        mulsd nb102_dzH1O(%rsp),%xmm12
        mulsd nb102_dxH2O(%rsp),%xmm13
        mulsd nb102_dyH2O(%rsp),%xmm14
        mulsd nb102_dzH2O(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb102_fixO(%rsp),%xmm7
    addsd nb102_fiyO(%rsp),%xmm8
    addsd nb102_fizO(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb102_fixH1(%rsp),%xmm10
    addsd nb102_fiyH1(%rsp),%xmm11
    addsd nb102_fizH1(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb102_fixH2(%rsp),%xmm13
    addsd nb102_fiyH2(%rsp),%xmm14
    addsd nb102_fizH2(%rsp),%xmm15

    movsd %xmm7,nb102_fixO(%rsp)
    movsd %xmm8,nb102_fiyO(%rsp)
    movsd %xmm9,nb102_fizO(%rsp)
    movsd %xmm10,nb102_fixH1(%rsp)
    movsd %xmm11,nb102_fiyH1(%rsp)
    movsd %xmm12,nb102_fizH1(%rsp)
    movsd %xmm13,nb102_fixH2(%rsp)
    movsd %xmm14,nb102_fiyH2(%rsp)
    movsd %xmm15,nb102_fizH2(%rsp)

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

    subsd nb102_ixO(%rsp),%xmm0
    subsd nb102_iyO(%rsp),%xmm1
    subsd nb102_izO(%rsp),%xmm2
    subsd nb102_ixH1(%rsp),%xmm3
    subsd nb102_iyH1(%rsp),%xmm4
    subsd nb102_izH1(%rsp),%xmm5
    subsd nb102_ixH2(%rsp),%xmm6
    subsd nb102_iyH2(%rsp),%xmm7
    subsd nb102_izH2(%rsp),%xmm8

        movsd %xmm0,nb102_dxOH1(%rsp)
        movsd %xmm1,nb102_dyOH1(%rsp)
        movsd %xmm2,nb102_dzOH1(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb102_dxH1H1(%rsp)
        movsd %xmm4,nb102_dyH1H1(%rsp)
        movsd %xmm5,nb102_dzH1H1(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb102_dxH2H1(%rsp)
        movsd %xmm7,nb102_dyH2H1(%rsp)
        movsd %xmm8,nb102_dzH2H1(%rsp)
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

        movsd  nb102_three(%rsp),%xmm9
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

        movsd  nb102_half(%rsp),%xmm15
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

        movsd  nb102_three(%rsp),%xmm1
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

        movsd  nb102_half(%rsp),%xmm15
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
    mulsd  nb102_qqOH(%rsp),%xmm0
    mulsd  nb102_qqHH(%rsp),%xmm1
    mulsd  nb102_qqHH(%rsp),%xmm2
    mulsd  %xmm0,%xmm9
    mulsd  %xmm1,%xmm10
    mulsd  %xmm2,%xmm11

    addsd nb102_vctot(%rsp),%xmm0
    addsd %xmm2,%xmm1
    addsd %xmm1,%xmm0
    movsd %xmm0,nb102_vctot(%rsp)

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

        mulsd nb102_dxOH1(%rsp),%xmm7
        mulsd nb102_dyOH1(%rsp),%xmm8
        mulsd nb102_dzOH1(%rsp),%xmm9
        mulsd nb102_dxH1H1(%rsp),%xmm10
        mulsd nb102_dyH1H1(%rsp),%xmm11
        mulsd nb102_dzH1H1(%rsp),%xmm12
        mulsd nb102_dxH2H1(%rsp),%xmm13
        mulsd nb102_dyH2H1(%rsp),%xmm14
        mulsd nb102_dzH2H1(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb102_fixO(%rsp),%xmm7
    addsd nb102_fiyO(%rsp),%xmm8
    addsd nb102_fizO(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb102_fixH1(%rsp),%xmm10
    addsd nb102_fiyH1(%rsp),%xmm11
    addsd nb102_fizH1(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb102_fixH2(%rsp),%xmm13
    addsd nb102_fiyH2(%rsp),%xmm14
    addsd nb102_fizH2(%rsp),%xmm15

    movsd %xmm7,nb102_fixO(%rsp)
    movsd %xmm8,nb102_fiyO(%rsp)
    movsd %xmm9,nb102_fizO(%rsp)
    movsd %xmm10,nb102_fixH1(%rsp)
    movsd %xmm11,nb102_fiyH1(%rsp)
    movsd %xmm12,nb102_fizH1(%rsp)
    movsd %xmm13,nb102_fixH2(%rsp)
    movsd %xmm14,nb102_fiyH2(%rsp)
    movsd %xmm15,nb102_fizH2(%rsp)

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

    subsd nb102_ixO(%rsp),%xmm0
    subsd nb102_iyO(%rsp),%xmm1
    subsd nb102_izO(%rsp),%xmm2
    subsd nb102_ixH1(%rsp),%xmm3
    subsd nb102_iyH1(%rsp),%xmm4
    subsd nb102_izH1(%rsp),%xmm5
    subsd nb102_ixH2(%rsp),%xmm6
    subsd nb102_iyH2(%rsp),%xmm7
    subsd nb102_izH2(%rsp),%xmm8

        movsd %xmm0,nb102_dxOH2(%rsp)
        movsd %xmm1,nb102_dyOH2(%rsp)
        movsd %xmm2,nb102_dzOH2(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb102_dxH1H2(%rsp)
        movsd %xmm4,nb102_dyH1H2(%rsp)
        movsd %xmm5,nb102_dzH1H2(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb102_dxH2H2(%rsp)
        movsd %xmm7,nb102_dyH2H2(%rsp)
        movsd %xmm8,nb102_dzH2H2(%rsp)
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

        movsd  nb102_three(%rsp),%xmm9
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

        movsd  nb102_half(%rsp),%xmm15
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

        movsd  nb102_three(%rsp),%xmm1
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

        movsd  nb102_half(%rsp),%xmm15
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
    mulsd  nb102_qqOH(%rsp),%xmm0
    mulsd  nb102_qqHH(%rsp),%xmm1
    mulsd  nb102_qqHH(%rsp),%xmm2
    mulsd  %xmm0,%xmm9
    mulsd  %xmm1,%xmm10
    mulsd  %xmm2,%xmm11

    addsd nb102_vctot(%rsp),%xmm0
    addsd %xmm2,%xmm1
    addsd %xmm1,%xmm0
    movsd %xmm0,nb102_vctot(%rsp)

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

        mulsd nb102_dxOH2(%rsp),%xmm7
        mulsd nb102_dyOH2(%rsp),%xmm8
        mulsd nb102_dzOH2(%rsp),%xmm9
        mulsd nb102_dxH1H2(%rsp),%xmm10
        mulsd nb102_dyH1H2(%rsp),%xmm11
        mulsd nb102_dzH1H2(%rsp),%xmm12
        mulsd nb102_dxH2H2(%rsp),%xmm13
        mulsd nb102_dyH2H2(%rsp),%xmm14
        mulsd nb102_dzH2H2(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb102_fixO(%rsp),%xmm7
    addsd nb102_fiyO(%rsp),%xmm8
    addsd nb102_fizO(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb102_fixH1(%rsp),%xmm10
    addsd nb102_fiyH1(%rsp),%xmm11
    addsd nb102_fizH1(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb102_fixH2(%rsp),%xmm13
    addsd nb102_fiyH2(%rsp),%xmm14
    addsd nb102_fizH2(%rsp),%xmm15

    movsd %xmm7,nb102_fixO(%rsp)
    movsd %xmm8,nb102_fiyO(%rsp)
    movsd %xmm9,nb102_fizO(%rsp)
    movsd %xmm10,nb102_fixH1(%rsp)
    movsd %xmm11,nb102_fiyH1(%rsp)
    movsd %xmm12,nb102_fizH1(%rsp)
    movsd %xmm13,nb102_fixH2(%rsp)
    movsd %xmm14,nb102_fiyH2(%rsp)
    movsd %xmm15,nb102_fizH2(%rsp)

    ## store back j H2 forces from xmm0-xmm2
        movsd %xmm0,48(%rdi,%rax,8)
        movsd %xmm1,56(%rdi,%rax,8)
        movsd %xmm2,64(%rdi,%rax,8)

_nb_kernel102_x86_64_sse2.nb102_updateouterdata: 
        movl  nb102_ii3(%rsp),%ecx
        movq  nb102_faction(%rbp),%rdi
        movq  nb102_fshift(%rbp),%rsi
        movl  nb102_is3(%rsp),%edx

        ## accumulate Oi forces in xmm0, xmm1, xmm2 
        movapd nb102_fixO(%rsp),%xmm0
        movapd nb102_fiyO(%rsp),%xmm1
        movapd nb102_fizO(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

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
        movapd nb102_fixH1(%rsp),%xmm0
        movapd nb102_fiyH1(%rsp),%xmm1
        movapd nb102_fizH1(%rsp),%xmm2

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
        movapd nb102_fixH2(%rsp),%xmm0
        movapd nb102_fiyH2(%rsp),%xmm1
        movapd nb102_fizH2(%rsp),%xmm2

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
        movl nb102_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb102_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb102_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 

        ## add earlier value from mem 
        movq  nb102_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb102_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel102_x86_64_sse2.nb102_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb102_n(%rsp)
        jmp _nb_kernel102_x86_64_sse2.nb102_outer
_nb_kernel102_x86_64_sse2.nb102_outerend: 
        ## check if more outer neighborlists remain
        movl  nb102_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel102_x86_64_sse2.nb102_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel102_x86_64_sse2.nb102_threadloop
_nb_kernel102_x86_64_sse2.nb102_end: 
        movl nb102_nouter(%rsp),%eax
        movl nb102_ninner(%rsp),%ebx
        movq nb102_outeriter(%rbp),%rcx
        movq nb102_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1496,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret


.globl nb_kernel102nf_x86_64_sse2
.globl _nb_kernel102nf_x86_64_sse2
nb_kernel102nf_x86_64_sse2:     
_nb_kernel102nf_x86_64_sse2:    
.set nb102nf_fshift, 16
.set nb102nf_gid, 24
.set nb102nf_pos, 32
.set nb102nf_faction, 40
.set nb102nf_charge, 48
.set nb102nf_p_facel, 56
.set nb102nf_argkrf, 64
.set nb102nf_argcrf, 72
.set nb102nf_Vc, 80
.set nb102nf_type, 88
.set nb102nf_p_ntype, 96
.set nb102nf_vdwparam, 104
.set nb102nf_Vvdw, 112
.set nb102nf_p_tabscale, 120
.set nb102nf_VFtab, 128
.set nb102nf_invsqrta, 136
.set nb102nf_dvda, 144
.set nb102nf_p_gbtabscale, 152
.set nb102nf_GBtab, 160
.set nb102nf_p_nthreads, 168
.set nb102nf_count, 176
.set nb102nf_mtx, 184
.set nb102nf_outeriter, 192
.set nb102nf_inneriter, 200
.set nb102nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb102nf_ixO, 0
.set nb102nf_iyO, 16
.set nb102nf_izO, 32
.set nb102nf_ixH1, 48
.set nb102nf_iyH1, 64
.set nb102nf_izH1, 80
.set nb102nf_ixH2, 96
.set nb102nf_iyH2, 112
.set nb102nf_izH2, 128
.set nb102nf_jxO, 144
.set nb102nf_jyO, 160
.set nb102nf_jzO, 176
.set nb102nf_jxH1, 192
.set nb102nf_jyH1, 208
.set nb102nf_jzH1, 224
.set nb102nf_jxH2, 240
.set nb102nf_jyH2, 256
.set nb102nf_jzH2, 272
.set nb102nf_qqOO, 288
.set nb102nf_qqOH, 304
.set nb102nf_qqHH, 320
.set nb102nf_vctot, 336
.set nb102nf_half, 352
.set nb102nf_three, 368
.set nb102nf_rsqOO, 384
.set nb102nf_rsqOH1, 400
.set nb102nf_rsqOH2, 416
.set nb102nf_rsqH1O, 432
.set nb102nf_rsqH1H1, 448
.set nb102nf_rsqH1H2, 464
.set nb102nf_rsqH2O, 480
.set nb102nf_rsqH2H1, 496
.set nb102nf_rsqH2H2, 512
.set nb102nf_rinvOO, 528
.set nb102nf_rinvOH1, 544
.set nb102nf_rinvOH2, 560
.set nb102nf_rinvH1O, 576
.set nb102nf_rinvH1H1, 592
.set nb102nf_rinvH1H2, 608
.set nb102nf_rinvH2O, 624
.set nb102nf_rinvH2H1, 640
.set nb102nf_rinvH2H2, 656
.set nb102nf_is3, 672
.set nb102nf_ii3, 676
.set nb102nf_nri, 680
.set nb102nf_iinr, 688
.set nb102nf_jindex, 696
.set nb102nf_jjnr, 704
.set nb102nf_shift, 712
.set nb102nf_shiftvec, 720
.set nb102nf_facel, 728
.set nb102nf_innerjjnr, 736
.set nb102nf_innerk, 744
.set nb102nf_n, 748
.set nb102nf_nn1, 752
.set nb102nf_nouter, 756
.set nb102nf_ninner, 760
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $792,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb102nf_nouter(%rsp)
        movl %eax,nb102nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb102nf_nri(%rsp)
        movq %rsi,nb102nf_iinr(%rsp)
        movq %rdx,nb102nf_jindex(%rsp)
        movq %rcx,nb102nf_jjnr(%rsp)
        movq %r8,nb102nf_shift(%rsp)
        movq %r9,nb102nf_shiftvec(%rsp)
        movq nb102nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb102nf_facel(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb102nf_half(%rsp)
        movl %ebx,nb102nf_half+4(%rsp)
        movsd nb102nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb102nf_half(%rsp)
        movapd %xmm3,nb102nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb102nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb102nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3       ## qO 
        movsd %xmm3,%xmm4               ## qO 
        movsd 8(%rdx,%rbx,8),%xmm5      ## qH 
        movq nb102nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb102nf_facel(%rsp),%xmm6         ## facel 
        mulsd  %xmm3,%xmm3              ## qO*qO 
        mulsd  %xmm5,%xmm4              ## qO*qH 
        mulsd  %xmm5,%xmm5              ## qH*qH 
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb102nf_qqOO(%rsp)
        movapd %xmm4,nb102nf_qqOH(%rsp)
        movapd %xmm5,nb102nf_qqHH(%rsp)

_nb_kernel102nf_x86_64_sse2.nb102nf_threadloop: 
        movq  nb102nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel102nf_x86_64_sse2.nb102nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel102nf_x86_64_sse2.nb102nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb102nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb102nf_n(%rsp)
        movl %ebx,nb102nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel102nf_x86_64_sse2.nb102nf_outerstart
        jmp _nb_kernel102nf_x86_64_sse2.nb102nf_end

_nb_kernel102nf_x86_64_sse2.nb102nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb102nf_nouter(%rsp),%ebx
        movl %ebx,nb102nf_nouter(%rsp)

_nb_kernel102nf_x86_64_sse2.nb102nf_outer: 
        movq  nb102nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb102nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb102nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb102nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb102nf_ii3(%rsp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%rax,%rbx,8),%xmm3
        addsd 8(%rax,%rbx,8),%xmm4
        addsd 16(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb102nf_ixO(%rsp)
        movapd %xmm4,nb102nf_iyO(%rsp)
        movapd %xmm5,nb102nf_izO(%rsp)

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
        movapd %xmm0,nb102nf_ixH1(%rsp)
        movapd %xmm1,nb102nf_iyH1(%rsp)
        movapd %xmm2,nb102nf_izH1(%rsp)
        movapd %xmm3,nb102nf_ixH2(%rsp)
        movapd %xmm4,nb102nf_iyH2(%rsp)
        movapd %xmm5,nb102nf_izH2(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb102nf_vctot(%rsp)

        movq  nb102nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb102nf_pos(%rbp),%rsi
        movq  nb102nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb102nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb102nf_ninner(%rsp),%ecx
        movl  %ecx,nb102nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb102nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel102nf_x86_64_sse2.nb102nf_unroll_loop
        jmp   _nb_kernel102nf_x86_64_sse2.nb102nf_checksingle
_nb_kernel102nf_x86_64_sse2.nb102nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb102nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb102nf_innerjjnr(%rsp)            ## advance pointer (unrolled 2) 

        movq nb102nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movapd  %xmm2,nb102nf_jxO(%rsp)
        movapd  %xmm3,nb102nf_jyO(%rsp)
        movapd  %xmm4,nb102nf_jzO(%rsp)
        movapd  %xmm5,nb102nf_jxH1(%rsp)
        movapd  %xmm6,nb102nf_jyH1(%rsp)
        movapd  %xmm7,nb102nf_jzH1(%rsp)
        movlpd 48(%rsi,%rax,8),%xmm2
        movlpd 56(%rsi,%rax,8),%xmm3
        movlpd 64(%rsi,%rax,8),%xmm4
        movhpd 48(%rsi,%rbx,8),%xmm2
        movhpd 56(%rsi,%rbx,8),%xmm3
        movhpd 64(%rsi,%rbx,8),%xmm4
        movapd  %xmm2,nb102nf_jxH2(%rsp)
        movapd  %xmm3,nb102nf_jyH2(%rsp)
        movapd  %xmm4,nb102nf_jzH2(%rsp)

        movapd nb102nf_ixO(%rsp),%xmm0
        movapd nb102nf_iyO(%rsp),%xmm1
        movapd nb102nf_izO(%rsp),%xmm2
        movapd nb102nf_ixO(%rsp),%xmm3
        movapd nb102nf_iyO(%rsp),%xmm4
        movapd nb102nf_izO(%rsp),%xmm5
        subpd  nb102nf_jxO(%rsp),%xmm0
        subpd  nb102nf_jyO(%rsp),%xmm1
        subpd  nb102nf_jzO(%rsp),%xmm2
        subpd  nb102nf_jxH1(%rsp),%xmm3
        subpd  nb102nf_jyH1(%rsp),%xmm4
        subpd  nb102nf_jzH1(%rsp),%xmm5
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
        movapd %xmm0,nb102nf_rsqOO(%rsp)
        movapd %xmm3,nb102nf_rsqOH1(%rsp)

        movapd nb102nf_ixO(%rsp),%xmm0
        movapd nb102nf_iyO(%rsp),%xmm1
        movapd nb102nf_izO(%rsp),%xmm2
        movapd nb102nf_ixH1(%rsp),%xmm3
        movapd nb102nf_iyH1(%rsp),%xmm4
        movapd nb102nf_izH1(%rsp),%xmm5
        subpd  nb102nf_jxH2(%rsp),%xmm0
        subpd  nb102nf_jyH2(%rsp),%xmm1
        subpd  nb102nf_jzH2(%rsp),%xmm2
        subpd  nb102nf_jxO(%rsp),%xmm3
        subpd  nb102nf_jyO(%rsp),%xmm4
        subpd  nb102nf_jzO(%rsp),%xmm5
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
        movapd %xmm0,nb102nf_rsqOH2(%rsp)
        movapd %xmm3,nb102nf_rsqH1O(%rsp)

        movapd nb102nf_ixH1(%rsp),%xmm0
        movapd nb102nf_iyH1(%rsp),%xmm1
        movapd nb102nf_izH1(%rsp),%xmm2
        movapd nb102nf_ixH1(%rsp),%xmm3
        movapd nb102nf_iyH1(%rsp),%xmm4
        movapd nb102nf_izH1(%rsp),%xmm5
        subpd  nb102nf_jxH1(%rsp),%xmm0
        subpd  nb102nf_jyH1(%rsp),%xmm1
        subpd  nb102nf_jzH1(%rsp),%xmm2
        subpd  nb102nf_jxH2(%rsp),%xmm3
        subpd  nb102nf_jyH2(%rsp),%xmm4
        subpd  nb102nf_jzH2(%rsp),%xmm5
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
        movapd %xmm0,nb102nf_rsqH1H1(%rsp)
        movapd %xmm3,nb102nf_rsqH1H2(%rsp)

        movapd nb102nf_ixH2(%rsp),%xmm0
        movapd nb102nf_iyH2(%rsp),%xmm1
        movapd nb102nf_izH2(%rsp),%xmm2
        movapd nb102nf_ixH2(%rsp),%xmm3
        movapd nb102nf_iyH2(%rsp),%xmm4
        movapd nb102nf_izH2(%rsp),%xmm5
        subpd  nb102nf_jxO(%rsp),%xmm0
        subpd  nb102nf_jyO(%rsp),%xmm1
        subpd  nb102nf_jzO(%rsp),%xmm2
        subpd  nb102nf_jxH1(%rsp),%xmm3
        subpd  nb102nf_jyH1(%rsp),%xmm4
        subpd  nb102nf_jzH1(%rsp),%xmm5
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
        movapd %xmm0,nb102nf_rsqH2O(%rsp)
        movapd %xmm4,nb102nf_rsqH2H1(%rsp)

        movapd nb102nf_ixH2(%rsp),%xmm0
        movapd nb102nf_iyH2(%rsp),%xmm1
        movapd nb102nf_izH2(%rsp),%xmm2
        subpd  nb102nf_jxH2(%rsp),%xmm0
        subpd  nb102nf_jyH2(%rsp),%xmm1
        subpd  nb102nf_jzH2(%rsp),%xmm2
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb102nf_rsqH2H2(%rsp)

        ## start doing invsqrt use rsq values in xmm0 (h2h2) , xmm4 (h2h1) 
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
        movapd  nb102nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102nf_half(%rsp),%xmm3   ## iter1 
        mulpd   nb102nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb102nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb102nf_rinvH2H2(%rsp)
        movapd %xmm5,nb102nf_rinvH2H1(%rsp)

        movapd nb102nf_rsqOO(%rsp),%xmm0
        movapd nb102nf_rsqOH1(%rsp),%xmm4
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
        movapd  nb102nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102nf_half(%rsp),%xmm3   ## iter1 of  
        mulpd   nb102nf_half(%rsp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb102nf_half(%rsp),%xmm5   ## rinv
        movapd %xmm1,nb102nf_rinvOO(%rsp)
        movapd %xmm5,nb102nf_rinvOH1(%rsp)

        movapd nb102nf_rsqOH2(%rsp),%xmm0
        movapd nb102nf_rsqH1O(%rsp),%xmm4
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
        movapd  nb102nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102nf_half(%rsp),%xmm3   ## iter1 
        mulpd   nb102nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb102nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb102nf_rinvOH2(%rsp)
        movapd %xmm5,nb102nf_rinvH1O(%rsp)

        movapd nb102nf_rsqH1H1(%rsp),%xmm0
        movapd nb102nf_rsqH1H2(%rsp),%xmm4
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
        movapd  nb102nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102nf_half(%rsp),%xmm3   ## iter1a 
        mulpd   nb102nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb102nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb102nf_rinvH1H1(%rsp)
        movapd %xmm5,nb102nf_rinvH1H2(%rsp)

        movapd nb102nf_rsqH2O(%rsp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb102nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb102nf_half(%rsp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb102nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb102nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb102nf_rinvH2O(%rsp)

        ## start with OO interaction 
        movapd nb102nf_rinvOO(%rsp),%xmm0
        mulpd  nb102nf_qqOO(%rsp),%xmm0
        addpd  nb102nf_vctot(%rsp),%xmm0

        ## other interactions 
        movapd nb102nf_rinvOH1(%rsp),%xmm1
        movapd nb102nf_rinvH1H1(%rsp),%xmm2

        addpd nb102nf_rinvOH2(%rsp),%xmm1
        addpd nb102nf_rinvH1H2(%rsp),%xmm2

        addpd nb102nf_rinvH1O(%rsp),%xmm1
        addpd nb102nf_rinvH2H1(%rsp),%xmm2

        addpd nb102nf_rinvH2O(%rsp),%xmm1
        addpd nb102nf_rinvH2H2(%rsp),%xmm2

        mulpd nb102nf_qqOH(%rsp),%xmm1
        mulpd nb102nf_qqHH(%rsp),%xmm2

        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0

        movapd %xmm0,nb102nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb102nf_innerk(%rsp)
        jl    _nb_kernel102nf_x86_64_sse2.nb102nf_checksingle
        jmp   _nb_kernel102nf_x86_64_sse2.nb102nf_unroll_loop
_nb_kernel102nf_x86_64_sse2.nb102nf_checksingle: 
        movl  nb102nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel102nf_x86_64_sse2.nb102nf_dosingle
        jmp   _nb_kernel102nf_x86_64_sse2.nb102nf_updateouterdata
_nb_kernel102nf_x86_64_sse2.nb102nf_dosingle: 
        movq  nb102nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb102nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coordinates to local temp variables 
        movlpd (%rsi,%rax,8),%xmm2
        movlpd 8(%rsi,%rax,8),%xmm3
        movlpd 16(%rsi,%rax,8),%xmm4
        movlpd 24(%rsi,%rax,8),%xmm5
        movlpd 32(%rsi,%rax,8),%xmm6
        movlpd 40(%rsi,%rax,8),%xmm7
        movapd  %xmm2,nb102nf_jxO(%rsp)
        movapd  %xmm3,nb102nf_jyO(%rsp)
        movapd  %xmm4,nb102nf_jzO(%rsp)
        movapd  %xmm5,nb102nf_jxH1(%rsp)
        movapd  %xmm6,nb102nf_jyH1(%rsp)
        movapd  %xmm7,nb102nf_jzH1(%rsp)
        movlpd 48(%rsi,%rax,8),%xmm2
        movlpd 56(%rsi,%rax,8),%xmm3
        movlpd 64(%rsi,%rax,8),%xmm4
        movapd  %xmm2,nb102nf_jxH2(%rsp)
        movapd  %xmm3,nb102nf_jyH2(%rsp)
        movapd  %xmm4,nb102nf_jzH2(%rsp)

        movapd nb102nf_ixO(%rsp),%xmm0
        movapd nb102nf_iyO(%rsp),%xmm1
        movapd nb102nf_izO(%rsp),%xmm2
        movapd nb102nf_ixO(%rsp),%xmm3
        movapd nb102nf_iyO(%rsp),%xmm4
        movapd nb102nf_izO(%rsp),%xmm5
        subsd  nb102nf_jxO(%rsp),%xmm0
        subsd  nb102nf_jyO(%rsp),%xmm1
        subsd  nb102nf_jzO(%rsp),%xmm2
        subsd  nb102nf_jxH1(%rsp),%xmm3
        subsd  nb102nf_jyH1(%rsp),%xmm4
        subsd  nb102nf_jzH1(%rsp),%xmm5
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
        movapd %xmm0,nb102nf_rsqOO(%rsp)
        movapd %xmm3,nb102nf_rsqOH1(%rsp)

        movapd nb102nf_ixO(%rsp),%xmm0
        movapd nb102nf_iyO(%rsp),%xmm1
        movapd nb102nf_izO(%rsp),%xmm2
        movapd nb102nf_ixH1(%rsp),%xmm3
        movapd nb102nf_iyH1(%rsp),%xmm4
        movapd nb102nf_izH1(%rsp),%xmm5
        subsd  nb102nf_jxH2(%rsp),%xmm0
        subsd  nb102nf_jyH2(%rsp),%xmm1
        subsd  nb102nf_jzH2(%rsp),%xmm2
        subsd  nb102nf_jxO(%rsp),%xmm3
        subsd  nb102nf_jyO(%rsp),%xmm4
        subsd  nb102nf_jzO(%rsp),%xmm5
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
        movapd %xmm0,nb102nf_rsqOH2(%rsp)
        movapd %xmm3,nb102nf_rsqH1O(%rsp)

        movapd nb102nf_ixH1(%rsp),%xmm0
        movapd nb102nf_iyH1(%rsp),%xmm1
        movapd nb102nf_izH1(%rsp),%xmm2
        movapd nb102nf_ixH1(%rsp),%xmm3
        movapd nb102nf_iyH1(%rsp),%xmm4
        movapd nb102nf_izH1(%rsp),%xmm5
        subsd  nb102nf_jxH1(%rsp),%xmm0
        subsd  nb102nf_jyH1(%rsp),%xmm1
        subsd  nb102nf_jzH1(%rsp),%xmm2
        subsd  nb102nf_jxH2(%rsp),%xmm3
        subsd  nb102nf_jyH2(%rsp),%xmm4
        subsd  nb102nf_jzH2(%rsp),%xmm5
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
        movapd %xmm0,nb102nf_rsqH1H1(%rsp)
        movapd %xmm3,nb102nf_rsqH1H2(%rsp)

        movapd nb102nf_ixH2(%rsp),%xmm0
        movapd nb102nf_iyH2(%rsp),%xmm1
        movapd nb102nf_izH2(%rsp),%xmm2
        movapd nb102nf_ixH2(%rsp),%xmm3
        movapd nb102nf_iyH2(%rsp),%xmm4
        movapd nb102nf_izH2(%rsp),%xmm5
        subsd  nb102nf_jxO(%rsp),%xmm0
        subsd  nb102nf_jyO(%rsp),%xmm1
        subsd  nb102nf_jzO(%rsp),%xmm2
        subsd  nb102nf_jxH1(%rsp),%xmm3
        subsd  nb102nf_jyH1(%rsp),%xmm4
        subsd  nb102nf_jzH1(%rsp),%xmm5
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
        movapd %xmm0,nb102nf_rsqH2O(%rsp)
        movapd %xmm4,nb102nf_rsqH2H1(%rsp)

        movapd nb102nf_ixH2(%rsp),%xmm0
        movapd nb102nf_iyH2(%rsp),%xmm1
        movapd nb102nf_izH2(%rsp),%xmm2
        subsd  nb102nf_jxH2(%rsp),%xmm0
        subsd  nb102nf_jyH2(%rsp),%xmm1
        subsd  nb102nf_jzH2(%rsp),%xmm2
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb102nf_rsqH2H2(%rsp)

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
        movapd  nb102nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102nf_half(%rsp),%xmm3   ## iter1 
        mulsd   nb102nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb102nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb102nf_rinvH2H2(%rsp)
        movapd %xmm5,nb102nf_rinvH2H1(%rsp)

        movapd nb102nf_rsqOO(%rsp),%xmm0
        movapd nb102nf_rsqOH1(%rsp),%xmm4
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
        movapd  nb102nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102nf_half(%rsp),%xmm3   ## iter1 of  
        mulsd   nb102nf_half(%rsp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb102nf_half(%rsp),%xmm5   ## rinv
        movapd %xmm1,nb102nf_rinvOO(%rsp)
        movapd %xmm5,nb102nf_rinvOH1(%rsp)

        movapd nb102nf_rsqOH2(%rsp),%xmm0
        movapd nb102nf_rsqH1O(%rsp),%xmm4
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
        movapd  nb102nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102nf_half(%rsp),%xmm3   ## iter1 
        mulsd   nb102nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb102nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb102nf_rinvOH2(%rsp)
        movapd %xmm5,nb102nf_rinvH1O(%rsp)

        movapd nb102nf_rsqH1H1(%rsp),%xmm0
        movapd nb102nf_rsqH1H2(%rsp),%xmm4
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
        movapd  nb102nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102nf_half(%rsp),%xmm3   ## iter1a 
        mulsd   nb102nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb102nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb102nf_rinvH1H1(%rsp)
        movapd %xmm5,nb102nf_rinvH1H2(%rsp)

        movapd nb102nf_rsqH2O(%rsp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb102nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb102nf_half(%rsp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb102nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb102nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb102nf_rinvH2O(%rsp)

        ## start with OO interaction 
        movapd nb102nf_rinvOO(%rsp),%xmm0
        mulpd  nb102nf_qqOO(%rsp),%xmm0
        addpd  nb102nf_vctot(%rsp),%xmm0

        ## other interactions 
        movapd nb102nf_rinvOH1(%rsp),%xmm1
        movapd nb102nf_rinvH1H1(%rsp),%xmm2

        addsd nb102nf_rinvOH2(%rsp),%xmm1
        addsd nb102nf_rinvH1H2(%rsp),%xmm2

        addsd nb102nf_rinvH1O(%rsp),%xmm1
        addsd nb102nf_rinvH2H1(%rsp),%xmm2

        addsd nb102nf_rinvH2O(%rsp),%xmm1
        addsd nb102nf_rinvH2H2(%rsp),%xmm2

        mulsd nb102nf_qqOH(%rsp),%xmm1
        mulsd nb102nf_qqHH(%rsp),%xmm2

        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0

        movlpd %xmm0,nb102nf_vctot(%rsp)

_nb_kernel102nf_x86_64_sse2.nb102nf_updateouterdata: 
        ## get n from stack
        movl nb102nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb102nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        movapd nb102nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 

        ## add earlier value from mem 
        movq  nb102nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb102nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel102nf_x86_64_sse2.nb102nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb102nf_n(%rsp)
        jmp _nb_kernel102nf_x86_64_sse2.nb102nf_outer
_nb_kernel102nf_x86_64_sse2.nb102nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb102nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel102nf_x86_64_sse2.nb102nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel102nf_x86_64_sse2.nb102nf_threadloop
_nb_kernel102nf_x86_64_sse2.nb102nf_end: 
        movl nb102nf_nouter(%rsp),%eax
        movl nb102nf_ninner(%rsp),%ebx
        movq nb102nf_outeriter(%rbp),%rcx
        movq nb102nf_inneriter(%rbp),%rdx
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

