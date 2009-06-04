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


.globl nb_kernel312_ia32_sse2
.globl _nb_kernel312_ia32_sse2
nb_kernel312_ia32_sse2: 
_nb_kernel312_ia32_sse2:        
.set nb312_p_nri, 8
.set nb312_iinr, 12
.set nb312_jindex, 16
.set nb312_jjnr, 20
.set nb312_shift, 24
.set nb312_shiftvec, 28
.set nb312_fshift, 32
.set nb312_gid, 36
.set nb312_pos, 40
.set nb312_faction, 44
.set nb312_charge, 48
.set nb312_p_facel, 52
.set nb312_argkrf, 56
.set nb312_argcrf, 60
.set nb312_Vc, 64
.set nb312_type, 68
.set nb312_p_ntype, 72
.set nb312_vdwparam, 76
.set nb312_Vvdw, 80
.set nb312_p_tabscale, 84
.set nb312_VFtab, 88
.set nb312_invsqrta, 92
.set nb312_dvda, 96
.set nb312_p_gbtabscale, 100
.set nb312_GBtab, 104
.set nb312_p_nthreads, 108
.set nb312_count, 112
.set nb312_mtx, 116
.set nb312_outeriter, 120
.set nb312_inneriter, 124
.set nb312_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb312_ixO, 0
.set nb312_iyO, 16
.set nb312_izO, 32
.set nb312_ixH1, 48
.set nb312_iyH1, 64
.set nb312_izH1, 80
.set nb312_ixH2, 96
.set nb312_iyH2, 112
.set nb312_izH2, 128
.set nb312_jxO, 144
.set nb312_jyO, 160
.set nb312_jzO, 176
.set nb312_jxH1, 192
.set nb312_jyH1, 208
.set nb312_jzH1, 224
.set nb312_jxH2, 240
.set nb312_jyH2, 256
.set nb312_jzH2, 272
.set nb312_dxOO, 288
.set nb312_dyOO, 304
.set nb312_dzOO, 320
.set nb312_dxOH1, 336
.set nb312_dyOH1, 352
.set nb312_dzOH1, 368
.set nb312_dxOH2, 384
.set nb312_dyOH2, 400
.set nb312_dzOH2, 416
.set nb312_dxH1O, 432
.set nb312_dyH1O, 448
.set nb312_dzH1O, 464
.set nb312_dxH1H1, 480
.set nb312_dyH1H1, 496
.set nb312_dzH1H1, 512
.set nb312_dxH1H2, 528
.set nb312_dyH1H2, 544
.set nb312_dzH1H2, 560
.set nb312_dxH2O, 576
.set nb312_dyH2O, 592
.set nb312_dzH2O, 608
.set nb312_dxH2H1, 624
.set nb312_dyH2H1, 640
.set nb312_dzH2H1, 656
.set nb312_dxH2H2, 672
.set nb312_dyH2H2, 688
.set nb312_dzH2H2, 704
.set nb312_qqOO, 720
.set nb312_qqOH, 736
.set nb312_qqHH, 752
.set nb312_two, 768
.set nb312_tsc, 784
.set nb312_c6, 800
.set nb312_c12, 816
.set nb312_six, 832
.set nb312_twelve, 848
.set nb312_vctot, 864
.set nb312_Vvdwtot, 880
.set nb312_fixO, 896
.set nb312_fiyO, 912
.set nb312_fizO, 928
.set nb312_fixH1, 944
.set nb312_fiyH1, 960
.set nb312_fizH1, 976
.set nb312_fixH2, 992
.set nb312_fiyH2, 1008
.set nb312_fizH2, 1024
.set nb312_fjxO, 1040
.set nb312_fjyO, 1056
.set nb312_fjzO, 1072
.set nb312_fjxH1, 1088
.set nb312_fjyH1, 1104
.set nb312_fjzH1, 1120
.set nb312_fjxH2, 1136
.set nb312_fjyH2, 1152
.set nb312_fjzH2, 1168
.set nb312_half, 1184
.set nb312_three, 1200
.set nb312_rsqOO, 1216
.set nb312_rsqOH1, 1232
.set nb312_rsqOH2, 1248
.set nb312_rsqH1O, 1264
.set nb312_rsqH1H1, 1280
.set nb312_rsqH1H2, 1296
.set nb312_rsqH2O, 1312
.set nb312_rsqH2H1, 1328
.set nb312_rsqH2H2, 1344
.set nb312_rinvOO, 1360
.set nb312_rinvOH1, 1376
.set nb312_rinvOH2, 1392
.set nb312_rinvH1O, 1408
.set nb312_rinvH1H1, 1424
.set nb312_rinvH1H2, 1440
.set nb312_rinvH2O, 1456
.set nb312_rinvH2H1, 1472
.set nb312_rinvH2H2, 1488
.set nb312_fstmp, 1504
.set nb312_is3, 1520
.set nb312_ii3, 1524
.set nb312_innerjjnr, 1528
.set nb312_innerk, 1532
.set nb312_n, 1536
.set nb312_nn1, 1540
.set nb312_nri, 1544
.set nb312_nouter, 1548
.set nb312_ninner, 1552
.set nb312_salign, 1556
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $1560,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb312_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb312_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb312_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb312_nouter(%esp)
        movl %eax,nb312_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb312_half(%esp)
        movl %ebx,nb312_half+4(%esp)
        movsd nb312_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm3,%xmm4
        addpd  %xmm4,%xmm4      ## 6.0
        movapd %xmm4,%xmm5
        addpd  %xmm5,%xmm5      ## 12.0
        movapd %xmm1,nb312_half(%esp)
        movapd %xmm2,nb312_two(%esp)
        movapd %xmm3,nb312_three(%esp)
        movapd %xmm4,nb312_six(%esp)
        movapd %xmm5,nb312_twelve(%esp)

        movl nb312_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm5
        shufpd $0,%xmm5,%xmm5
        movapd %xmm5,nb312_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb312_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb312_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb312_p_facel(%ebp),%esi
        movsd (%esi),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb312_qqOO(%esp)
        movapd %xmm4,nb312_qqOH(%esp)
        movapd %xmm5,nb312_qqHH(%esp)

        xorpd %xmm0,%xmm0
        movl  nb312_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb312_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb312_vdwparam(%ebp),%eax
        movlpd (%eax,%edx,8),%xmm0
        movlpd 8(%eax,%edx,8),%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb312_c6(%esp)
        movapd %xmm1,nb312_c12(%esp)

_nb_kernel312_ia32_sse2.nb312_threadloop: 
        movl  nb312_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel312_ia32_sse2.nb312_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel312_ia32_sse2.nb312_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb312_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb312_n(%esp)
        movl %ebx,nb312_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel312_ia32_sse2.nb312_outerstart
        jmp _nb_kernel312_ia32_sse2.nb312_end

_nb_kernel312_ia32_sse2.nb312_outerstart: 
        ## ebx contains number of outer iterations
        addl nb312_nouter(%esp),%ebx
        movl %ebx,nb312_nouter(%esp)

_nb_kernel312_ia32_sse2.nb312_outer: 
        movl  nb312_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb312_is3(%esp)      ## store is3 

        movl  nb312_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb312_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb312_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb312_ii3(%esp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb312_ixO(%esp)
        movapd %xmm4,nb312_iyO(%esp)
        movapd %xmm5,nb312_izO(%esp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%eax,%ebx,8),%xmm0
        addsd 32(%eax,%ebx,8),%xmm1
        addsd 40(%eax,%ebx,8),%xmm2
        addsd 48(%eax,%ebx,8),%xmm3
        addsd 56(%eax,%ebx,8),%xmm4
        addsd 64(%eax,%ebx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb312_ixH1(%esp)
        movapd %xmm1,nb312_iyH1(%esp)
        movapd %xmm2,nb312_izH1(%esp)
        movapd %xmm3,nb312_ixH2(%esp)
        movapd %xmm4,nb312_iyH2(%esp)
        movapd %xmm5,nb312_izH2(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb312_vctot(%esp)
        movapd %xmm4,nb312_Vvdwtot(%esp)
        movapd %xmm4,nb312_fixO(%esp)
        movapd %xmm4,nb312_fiyO(%esp)
        movapd %xmm4,nb312_fizO(%esp)
        movapd %xmm4,nb312_fixH1(%esp)
        movapd %xmm4,nb312_fiyH1(%esp)
        movapd %xmm4,nb312_fizH1(%esp)
        movapd %xmm4,nb312_fixH2(%esp)
        movapd %xmm4,nb312_fiyH2(%esp)
        movapd %xmm4,nb312_fizH2(%esp)

        movl  nb312_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb312_pos(%ebp),%esi
        movl  nb312_faction(%ebp),%edi
        movl  nb312_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb312_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb312_ninner(%esp),%ecx
        movl  %ecx,nb312_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb312_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel312_ia32_sse2.nb312_unroll_loop
        jmp   _nb_kernel312_ia32_sse2.nb312_checksingle
_nb_kernel312_ia32_sse2.nb312_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb312_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb312_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb312_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move j coordinates to local temp variables 
        movlpd (%esi,%eax,8),%xmm2
        movlpd 8(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%eax,8),%xmm4
        movlpd 24(%esi,%eax,8),%xmm5
        movlpd 32(%esi,%eax,8),%xmm6
        movlpd 40(%esi,%eax,8),%xmm7
        movhpd (%esi,%ebx,8),%xmm2
        movhpd 8(%esi,%ebx,8),%xmm3
        movhpd 16(%esi,%ebx,8),%xmm4
        movhpd 24(%esi,%ebx,8),%xmm5
        movhpd 32(%esi,%ebx,8),%xmm6
        movhpd 40(%esi,%ebx,8),%xmm7
        movapd  %xmm2,nb312_jxO(%esp)
        movapd  %xmm3,nb312_jyO(%esp)
        movapd  %xmm4,nb312_jzO(%esp)
        movapd  %xmm5,nb312_jxH1(%esp)
        movapd  %xmm6,nb312_jyH1(%esp)
        movapd  %xmm7,nb312_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 48(%esi,%ebx,8),%xmm2
        movhpd 56(%esi,%ebx,8),%xmm3
        movhpd 64(%esi,%ebx,8),%xmm4
        movapd  %xmm2,nb312_jxH2(%esp)
        movapd  %xmm3,nb312_jyH2(%esp)
        movapd  %xmm4,nb312_jzH2(%esp)

        movapd nb312_ixO(%esp),%xmm0
        movapd nb312_iyO(%esp),%xmm1
        movapd nb312_izO(%esp),%xmm2
        movapd nb312_ixO(%esp),%xmm3
        movapd nb312_iyO(%esp),%xmm4
        movapd nb312_izO(%esp),%xmm5
        subpd  nb312_jxO(%esp),%xmm0
        subpd  nb312_jyO(%esp),%xmm1
        subpd  nb312_jzO(%esp),%xmm2
        subpd  nb312_jxH1(%esp),%xmm3
        subpd  nb312_jyH1(%esp),%xmm4
        subpd  nb312_jzH1(%esp),%xmm5
        movapd %xmm0,nb312_dxOO(%esp)
        movapd %xmm1,nb312_dyOO(%esp)
        movapd %xmm2,nb312_dzOO(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb312_dxOH1(%esp)
        movapd %xmm4,nb312_dyOH1(%esp)
        movapd %xmm5,nb312_dzOH1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb312_rsqOO(%esp)
        movapd %xmm3,nb312_rsqOH1(%esp)

        movapd nb312_ixO(%esp),%xmm0
        movapd nb312_iyO(%esp),%xmm1
        movapd nb312_izO(%esp),%xmm2
        movapd nb312_ixH1(%esp),%xmm3
        movapd nb312_iyH1(%esp),%xmm4
        movapd nb312_izH1(%esp),%xmm5
        subpd  nb312_jxH2(%esp),%xmm0
        subpd  nb312_jyH2(%esp),%xmm1
        subpd  nb312_jzH2(%esp),%xmm2
        subpd  nb312_jxO(%esp),%xmm3
        subpd  nb312_jyO(%esp),%xmm4
        subpd  nb312_jzO(%esp),%xmm5
        movapd %xmm0,nb312_dxOH2(%esp)
        movapd %xmm1,nb312_dyOH2(%esp)
        movapd %xmm2,nb312_dzOH2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb312_dxH1O(%esp)
        movapd %xmm4,nb312_dyH1O(%esp)
        movapd %xmm5,nb312_dzH1O(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb312_rsqOH2(%esp)
        movapd %xmm3,nb312_rsqH1O(%esp)

        movapd nb312_ixH1(%esp),%xmm0
        movapd nb312_iyH1(%esp),%xmm1
        movapd nb312_izH1(%esp),%xmm2
        movapd nb312_ixH1(%esp),%xmm3
        movapd nb312_iyH1(%esp),%xmm4
        movapd nb312_izH1(%esp),%xmm5
        subpd  nb312_jxH1(%esp),%xmm0
        subpd  nb312_jyH1(%esp),%xmm1
        subpd  nb312_jzH1(%esp),%xmm2
        subpd  nb312_jxH2(%esp),%xmm3
        subpd  nb312_jyH2(%esp),%xmm4
        subpd  nb312_jzH2(%esp),%xmm5
        movapd %xmm0,nb312_dxH1H1(%esp)
        movapd %xmm1,nb312_dyH1H1(%esp)
        movapd %xmm2,nb312_dzH1H1(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb312_dxH1H2(%esp)
        movapd %xmm4,nb312_dyH1H2(%esp)
        movapd %xmm5,nb312_dzH1H2(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb312_rsqH1H1(%esp)
        movapd %xmm3,nb312_rsqH1H2(%esp)

        movapd nb312_ixH2(%esp),%xmm0
        movapd nb312_iyH2(%esp),%xmm1
        movapd nb312_izH2(%esp),%xmm2
        movapd nb312_ixH2(%esp),%xmm3
        movapd nb312_iyH2(%esp),%xmm4
        movapd nb312_izH2(%esp),%xmm5
        subpd  nb312_jxO(%esp),%xmm0
        subpd  nb312_jyO(%esp),%xmm1
        subpd  nb312_jzO(%esp),%xmm2
        subpd  nb312_jxH1(%esp),%xmm3
        subpd  nb312_jyH1(%esp),%xmm4
        subpd  nb312_jzH1(%esp),%xmm5
        movapd %xmm0,nb312_dxH2O(%esp)
        movapd %xmm1,nb312_dyH2O(%esp)
        movapd %xmm2,nb312_dzH2O(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb312_dxH2H1(%esp)
        movapd %xmm4,nb312_dyH2H1(%esp)
        movapd %xmm5,nb312_dzH2H1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb312_rsqH2O(%esp)
        movapd %xmm4,nb312_rsqH2H1(%esp)

        movapd nb312_ixH2(%esp),%xmm0
        movapd nb312_iyH2(%esp),%xmm1
        movapd nb312_izH2(%esp),%xmm2
        subpd  nb312_jxH2(%esp),%xmm0
        subpd  nb312_jyH2(%esp),%xmm1
        subpd  nb312_jzH2(%esp),%xmm2
        movapd %xmm0,nb312_dxH2H2(%esp)
        movapd %xmm1,nb312_dyH2H2(%esp)
        movapd %xmm2,nb312_dzH2H2(%esp)
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb312_rsqH2H2(%esp)

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
        movapd  nb312_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb312_half(%esp),%xmm3   ## iter1 
        mulpd   nb312_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb312_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb312_half(%esp),%xmm1   ## rinv 
        mulpd   nb312_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb312_rinvH2H2(%esp)
        movapd %xmm5,nb312_rinvH2H1(%esp)

        movapd nb312_rsqOO(%esp),%xmm0
        movapd nb312_rsqOH1(%esp),%xmm4
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
        movapd  nb312_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb312_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb312_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb312_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb312_half(%esp),%xmm1   ## rinv 
        mulpd   nb312_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb312_rinvOO(%esp)
        movapd %xmm5,nb312_rinvOH1(%esp)

        movapd nb312_rsqOH2(%esp),%xmm0
        movapd nb312_rsqH1O(%esp),%xmm4
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
        movapd  nb312_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb312_half(%esp),%xmm3   ## iter1 
        mulpd   nb312_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb312_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb312_half(%esp),%xmm1   ## rinv 
        mulpd   nb312_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb312_rinvOH2(%esp)
        movapd %xmm5,nb312_rinvH1O(%esp)

        movapd nb312_rsqH1H1(%esp),%xmm0
        movapd nb312_rsqH1H2(%esp),%xmm4
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
        movapd  nb312_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb312_half(%esp),%xmm3   ## iter1a 
        mulpd   nb312_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb312_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb312_half(%esp),%xmm1   ## rinv 
        mulpd   nb312_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb312_rinvH1H1(%esp)
        movapd %xmm5,nb312_rinvH1H2(%esp)

        movapd nb312_rsqH2O(%esp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb312_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb312_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb312_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb312_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb312_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb312_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb312_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulpd  nb312_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movd %eax,%mm0
        movd %ebx,%mm1
        movl nb312_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb312_two(%esp),%xmm7    ## two*Heps2 
        movapd nb312_qqOO(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul - then we can get rid of mm5 
    ## update vctot 
    addpd  nb312_vctot(%esp),%xmm5
        xorpd  %xmm2,%xmm2
    movapd %xmm5,nb312_vctot(%esp)
        mulpd  nb312_tsc(%esp),%xmm3

        ## start doing lj 
        movapd %xmm0,%xmm2
        mulpd  %xmm2,%xmm2
        movapd %xmm2,%xmm1
        mulpd  %xmm2,%xmm1
        mulpd  %xmm2,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulpd  nb312_c6(%esp),%xmm1
        mulpd  nb312_c12(%esp),%xmm2
        movapd %xmm2,%xmm4
        subpd  %xmm1,%xmm4
        addpd  nb312_Vvdwtot(%esp),%xmm4
        mulpd  nb312_six(%esp),%xmm1
        mulpd  nb312_twelve(%esp),%xmm2
        movapd %xmm4,nb312_Vvdwtot(%esp)
        subpd  %xmm1,%xmm2
        mulpd  %xmm0,%xmm2

        subpd  %xmm3,%xmm2
        mulpd  %xmm2,%xmm0

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb312_dxOO(%esp),%xmm0
        mulpd nb312_dyOO(%esp),%xmm1
        mulpd nb312_dzOO(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb312_fixO(%esp),%xmm0
        addpd nb312_fiyO(%esp),%xmm1
        addpd nb312_fizO(%esp),%xmm2
        movapd %xmm3,nb312_fjxO(%esp)
        movapd %xmm4,nb312_fjyO(%esp)
        movapd %xmm5,nb312_fjzO(%esp)
        movapd %xmm0,nb312_fixO(%esp)
        movapd %xmm1,nb312_fiyO(%esp)
        movapd %xmm2,nb312_fizO(%esp)

        ## O-H1 interaction 
        movapd nb312_rinvOH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb312_rsqOH1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb312_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb312_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb312_two(%esp),%xmm7    ## two*Heps2 
        movapd nb312_qqOH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb312_vctot(%esp),%xmm5
    movapd %xmm5,nb312_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb312_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb312_dxOH1(%esp),%xmm0
        mulpd nb312_dyOH1(%esp),%xmm1
        mulpd nb312_dzOH1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb312_fixO(%esp),%xmm0
        addpd nb312_fiyO(%esp),%xmm1
        addpd nb312_fizO(%esp),%xmm2
        movapd %xmm3,nb312_fjxH1(%esp)
        movapd %xmm4,nb312_fjyH1(%esp)
        movapd %xmm5,nb312_fjzH1(%esp)
        movapd %xmm0,nb312_fixO(%esp)
        movapd %xmm1,nb312_fiyO(%esp)
        movapd %xmm2,nb312_fizO(%esp)

        ## O-H2 interaction  
        movapd nb312_rinvOH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb312_rsqOH2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb312_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb312_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb312_two(%esp),%xmm7    ## two*Heps2 
        movapd nb312_qqOH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb312_vctot(%esp),%xmm5
    movapd %xmm5,nb312_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb312_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb312_dxOH2(%esp),%xmm0
        mulpd nb312_dyOH2(%esp),%xmm1
        mulpd nb312_dzOH2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb312_fixO(%esp),%xmm0
        addpd nb312_fiyO(%esp),%xmm1
        addpd nb312_fizO(%esp),%xmm2
        movapd %xmm3,nb312_fjxH2(%esp)
        movapd %xmm4,nb312_fjyH2(%esp)
        movapd %xmm5,nb312_fjzH2(%esp)
        movapd %xmm0,nb312_fixO(%esp)
        movapd %xmm1,nb312_fiyO(%esp)
        movapd %xmm2,nb312_fizO(%esp)

        ## H1-O interaction 
        movapd nb312_rinvH1O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb312_rsqH1O(%esp),%xmm1   ## xmm1=r 
        mulpd  nb312_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb312_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb312_two(%esp),%xmm7    ## two*Heps2 
        movapd nb312_qqOH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb312_vctot(%esp),%xmm5
    movapd %xmm5,nb312_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb312_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb312_fjxO(%esp),%xmm3
        movapd nb312_fjyO(%esp),%xmm4
        movapd nb312_fjzO(%esp),%xmm5
        mulpd nb312_dxH1O(%esp),%xmm0
        mulpd nb312_dyH1O(%esp),%xmm1
        mulpd nb312_dzH1O(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb312_fixH1(%esp),%xmm0
        addpd nb312_fiyH1(%esp),%xmm1
        addpd nb312_fizH1(%esp),%xmm2
        movapd %xmm3,nb312_fjxO(%esp)
        movapd %xmm4,nb312_fjyO(%esp)
        movapd %xmm5,nb312_fjzO(%esp)
        movapd %xmm0,nb312_fixH1(%esp)
        movapd %xmm1,nb312_fiyH1(%esp)
        movapd %xmm2,nb312_fizH1(%esp)

        ## H1-H1 interaction 
        movapd nb312_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb312_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb312_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb312_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb312_two(%esp),%xmm7    ## two*Heps2 
        movapd nb312_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb312_vctot(%esp),%xmm5
    movapd %xmm5,nb312_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb312_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb312_fjxH1(%esp),%xmm3
        movapd nb312_fjyH1(%esp),%xmm4
        movapd nb312_fjzH1(%esp),%xmm5
        mulpd nb312_dxH1H1(%esp),%xmm0
        mulpd nb312_dyH1H1(%esp),%xmm1
        mulpd nb312_dzH1H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb312_fixH1(%esp),%xmm0
        addpd nb312_fiyH1(%esp),%xmm1
        addpd nb312_fizH1(%esp),%xmm2
        movapd %xmm3,nb312_fjxH1(%esp)
        movapd %xmm4,nb312_fjyH1(%esp)
        movapd %xmm5,nb312_fjzH1(%esp)
        movapd %xmm0,nb312_fixH1(%esp)
        movapd %xmm1,nb312_fiyH1(%esp)
        movapd %xmm2,nb312_fizH1(%esp)

        ## H1-H2 interaction 
        movapd nb312_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb312_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb312_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb312_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb312_two(%esp),%xmm7    ## two*Heps2 
        movapd nb312_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb312_vctot(%esp),%xmm5
    movapd %xmm5,nb312_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb312_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb312_fjxH2(%esp),%xmm3
        movapd nb312_fjyH2(%esp),%xmm4
        movapd nb312_fjzH2(%esp),%xmm5
        mulpd nb312_dxH1H2(%esp),%xmm0
        mulpd nb312_dyH1H2(%esp),%xmm1
        mulpd nb312_dzH1H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb312_fixH1(%esp),%xmm0
        addpd nb312_fiyH1(%esp),%xmm1
        addpd nb312_fizH1(%esp),%xmm2
        movapd %xmm3,nb312_fjxH2(%esp)
        movapd %xmm4,nb312_fjyH2(%esp)
        movapd %xmm5,nb312_fjzH2(%esp)
        movapd %xmm0,nb312_fixH1(%esp)
        movapd %xmm1,nb312_fiyH1(%esp)
        movapd %xmm2,nb312_fizH1(%esp)

        ## H2-O interaction 
        movapd nb312_rinvH2O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb312_rsqH2O(%esp),%xmm1   ## xmm1=r 
        mulpd  nb312_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb312_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb312_two(%esp),%xmm7    ## two*Heps2 
        movapd nb312_qqOH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb312_vctot(%esp),%xmm5
    movapd %xmm5,nb312_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb312_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb312_fjxO(%esp),%xmm3
        movapd nb312_fjyO(%esp),%xmm4
        movapd nb312_fjzO(%esp),%xmm5
        mulpd nb312_dxH2O(%esp),%xmm0
        mulpd nb312_dyH2O(%esp),%xmm1
        mulpd nb312_dzH2O(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb312_fixH2(%esp),%xmm0
        addpd nb312_fiyH2(%esp),%xmm1
        addpd nb312_fizH2(%esp),%xmm2
        movapd %xmm3,nb312_fjxO(%esp)
        movapd %xmm4,nb312_fjyO(%esp)
        movapd %xmm5,nb312_fjzO(%esp)
        movapd %xmm0,nb312_fixH2(%esp)
        movapd %xmm1,nb312_fiyH2(%esp)
        movapd %xmm2,nb312_fizH2(%esp)

        ## H2-H1 interaction 
        movapd nb312_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb312_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb312_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb312_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb312_two(%esp),%xmm7    ## two*Heps2 
        movapd nb312_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb312_vctot(%esp),%xmm5
    movapd %xmm5,nb312_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb312_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb312_fjxH1(%esp),%xmm3
        movapd nb312_fjyH1(%esp),%xmm4
        movapd nb312_fjzH1(%esp),%xmm5
        mulpd nb312_dxH2H1(%esp),%xmm0
        mulpd nb312_dyH2H1(%esp),%xmm1
        mulpd nb312_dzH2H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb312_fixH2(%esp),%xmm0
        addpd nb312_fiyH2(%esp),%xmm1
        addpd nb312_fizH2(%esp),%xmm2
        movapd %xmm3,nb312_fjxH1(%esp)
        movapd %xmm4,nb312_fjyH1(%esp)
        movapd %xmm5,nb312_fjzH1(%esp)
        movapd %xmm0,nb312_fixH2(%esp)
        movapd %xmm1,nb312_fiyH2(%esp)
        movapd %xmm2,nb312_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb312_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb312_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb312_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb312_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb312_two(%esp),%xmm7    ## two*Heps2 
        movapd nb312_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb312_vctot(%esp),%xmm5
    movapd %xmm5,nb312_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb312_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb312_fjxH2(%esp),%xmm3
        movapd nb312_fjyH2(%esp),%xmm4
        movapd nb312_fjzH2(%esp),%xmm5
        mulpd nb312_dxH2H2(%esp),%xmm0
        mulpd nb312_dyH2H2(%esp),%xmm1
        mulpd nb312_dzH2H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb312_fixH2(%esp),%xmm0
        addpd nb312_fiyH2(%esp),%xmm1
        addpd nb312_fizH2(%esp),%xmm2
        movapd %xmm3,nb312_fjxH2(%esp)
        movapd %xmm4,nb312_fjyH2(%esp)
        movapd %xmm5,nb312_fjzH2(%esp)
        movapd %xmm0,nb312_fixH2(%esp)
        movapd %xmm1,nb312_fiyH2(%esp)
        movapd %xmm2,nb312_fizH2(%esp)

        movl nb312_faction(%ebp),%edi

        movd %mm0,%eax
        movd %mm1,%ebx

        ## Did all interactions - now update j forces 
        movlpd (%edi,%eax,8),%xmm0
        movlpd 8(%edi,%eax,8),%xmm1
        movlpd 16(%edi,%eax,8),%xmm2
        movlpd 24(%edi,%eax,8),%xmm3
        movlpd 32(%edi,%eax,8),%xmm4
        movlpd 40(%edi,%eax,8),%xmm5
        movlpd 48(%edi,%eax,8),%xmm6
        movlpd 56(%edi,%eax,8),%xmm7
        movhpd (%edi,%ebx,8),%xmm0
        movhpd 8(%edi,%ebx,8),%xmm1
        movhpd 16(%edi,%ebx,8),%xmm2
        movhpd 24(%edi,%ebx,8),%xmm3
        movhpd 32(%edi,%ebx,8),%xmm4
        movhpd 40(%edi,%ebx,8),%xmm5
        movhpd 48(%edi,%ebx,8),%xmm6
        movhpd 56(%edi,%ebx,8),%xmm7
        addpd nb312_fjxO(%esp),%xmm0
        addpd nb312_fjyO(%esp),%xmm1
        addpd nb312_fjzO(%esp),%xmm2
        addpd nb312_fjxH1(%esp),%xmm3
        addpd nb312_fjyH1(%esp),%xmm4
        addpd nb312_fjzH1(%esp),%xmm5
        addpd nb312_fjxH2(%esp),%xmm6
        addpd nb312_fjyH2(%esp),%xmm7
        movlpd %xmm0,(%edi,%eax,8)
        movlpd %xmm1,8(%edi,%eax,8)
        movlpd %xmm2,16(%edi,%eax,8)
        movlpd %xmm3,24(%edi,%eax,8)
        movlpd %xmm4,32(%edi,%eax,8)
        movlpd %xmm5,40(%edi,%eax,8)
        movlpd %xmm6,48(%edi,%eax,8)
        movlpd %xmm7,56(%edi,%eax,8)
        movhpd %xmm0,(%edi,%ebx,8)
        movhpd %xmm1,8(%edi,%ebx,8)
        movhpd %xmm2,16(%edi,%ebx,8)
        movhpd %xmm3,24(%edi,%ebx,8)
        movhpd %xmm4,32(%edi,%ebx,8)
        movhpd %xmm5,40(%edi,%ebx,8)
        movhpd %xmm6,48(%edi,%ebx,8)
        movhpd %xmm7,56(%edi,%ebx,8)

        movlpd 64(%edi,%eax,8),%xmm0
        movhpd 64(%edi,%ebx,8),%xmm0
        addpd nb312_fjzH2(%esp),%xmm0
        movlpd %xmm0,64(%edi,%eax,8)
        movhpd %xmm0,64(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl $2,nb312_innerk(%esp)
        jl    _nb_kernel312_ia32_sse2.nb312_checksingle
        jmp   _nb_kernel312_ia32_sse2.nb312_unroll_loop
_nb_kernel312_ia32_sse2.nb312_checksingle: 
        movl  nb312_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel312_ia32_sse2.nb312_dosingle
        jmp   _nb_kernel312_ia32_sse2.nb312_updateouterdata
_nb_kernel312_ia32_sse2.nb312_dosingle: 
        movl  nb312_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb312_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        movlpd (%esi,%eax,8),%xmm2
        movlpd 8(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%eax,8),%xmm4
        movlpd 24(%esi,%eax,8),%xmm5
        movlpd 32(%esi,%eax,8),%xmm6
        movlpd 40(%esi,%eax,8),%xmm7
        movapd  %xmm2,nb312_jxO(%esp)
        movapd  %xmm3,nb312_jyO(%esp)
        movapd  %xmm4,nb312_jzO(%esp)
        movapd  %xmm5,nb312_jxH1(%esp)
        movapd  %xmm6,nb312_jyH1(%esp)
        movapd  %xmm7,nb312_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movapd  %xmm2,nb312_jxH2(%esp)
        movapd  %xmm3,nb312_jyH2(%esp)
        movapd  %xmm4,nb312_jzH2(%esp)

        movapd nb312_ixO(%esp),%xmm0
        movapd nb312_iyO(%esp),%xmm1
        movapd nb312_izO(%esp),%xmm2
        movapd nb312_ixO(%esp),%xmm3
        movapd nb312_iyO(%esp),%xmm4
        movapd nb312_izO(%esp),%xmm5
        subsd  nb312_jxO(%esp),%xmm0
        subsd  nb312_jyO(%esp),%xmm1
        subsd  nb312_jzO(%esp),%xmm2
        subsd  nb312_jxH1(%esp),%xmm3
        subsd  nb312_jyH1(%esp),%xmm4
        subsd  nb312_jzH1(%esp),%xmm5
        movapd %xmm0,nb312_dxOO(%esp)
        movapd %xmm1,nb312_dyOO(%esp)
        movapd %xmm2,nb312_dzOO(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb312_dxOH1(%esp)
        movapd %xmm4,nb312_dyOH1(%esp)
        movapd %xmm5,nb312_dzOH1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb312_rsqOO(%esp)
        movapd %xmm3,nb312_rsqOH1(%esp)

        movapd nb312_ixO(%esp),%xmm0
        movapd nb312_iyO(%esp),%xmm1
        movapd nb312_izO(%esp),%xmm2
        movapd nb312_ixH1(%esp),%xmm3
        movapd nb312_iyH1(%esp),%xmm4
        movapd nb312_izH1(%esp),%xmm5
        subsd  nb312_jxH2(%esp),%xmm0
        subsd  nb312_jyH2(%esp),%xmm1
        subsd  nb312_jzH2(%esp),%xmm2
        subsd  nb312_jxO(%esp),%xmm3
        subsd  nb312_jyO(%esp),%xmm4
        subsd  nb312_jzO(%esp),%xmm5
        movapd %xmm0,nb312_dxOH2(%esp)
        movapd %xmm1,nb312_dyOH2(%esp)
        movapd %xmm2,nb312_dzOH2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb312_dxH1O(%esp)
        movapd %xmm4,nb312_dyH1O(%esp)
        movapd %xmm5,nb312_dzH1O(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb312_rsqOH2(%esp)
        movapd %xmm3,nb312_rsqH1O(%esp)

        movapd nb312_ixH1(%esp),%xmm0
        movapd nb312_iyH1(%esp),%xmm1
        movapd nb312_izH1(%esp),%xmm2
        movapd nb312_ixH1(%esp),%xmm3
        movapd nb312_iyH1(%esp),%xmm4
        movapd nb312_izH1(%esp),%xmm5
        subsd  nb312_jxH1(%esp),%xmm0
        subsd  nb312_jyH1(%esp),%xmm1
        subsd  nb312_jzH1(%esp),%xmm2
        subsd  nb312_jxH2(%esp),%xmm3
        subsd  nb312_jyH2(%esp),%xmm4
        subsd  nb312_jzH2(%esp),%xmm5
        movapd %xmm0,nb312_dxH1H1(%esp)
        movapd %xmm1,nb312_dyH1H1(%esp)
        movapd %xmm2,nb312_dzH1H1(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb312_dxH1H2(%esp)
        movapd %xmm4,nb312_dyH1H2(%esp)
        movapd %xmm5,nb312_dzH1H2(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb312_rsqH1H1(%esp)
        movapd %xmm3,nb312_rsqH1H2(%esp)

        movapd nb312_ixH2(%esp),%xmm0
        movapd nb312_iyH2(%esp),%xmm1
        movapd nb312_izH2(%esp),%xmm2
        movapd nb312_ixH2(%esp),%xmm3
        movapd nb312_iyH2(%esp),%xmm4
        movapd nb312_izH2(%esp),%xmm5
        subsd  nb312_jxO(%esp),%xmm0
        subsd  nb312_jyO(%esp),%xmm1
        subsd  nb312_jzO(%esp),%xmm2
        subsd  nb312_jxH1(%esp),%xmm3
        subsd  nb312_jyH1(%esp),%xmm4
        subsd  nb312_jzH1(%esp),%xmm5
        movapd %xmm0,nb312_dxH2O(%esp)
        movapd %xmm1,nb312_dyH2O(%esp)
        movapd %xmm2,nb312_dzH2O(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb312_dxH2H1(%esp)
        movapd %xmm4,nb312_dyH2H1(%esp)
        movapd %xmm5,nb312_dzH2H1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb312_rsqH2O(%esp)
        movapd %xmm4,nb312_rsqH2H1(%esp)

        movapd nb312_ixH2(%esp),%xmm0
        movapd nb312_iyH2(%esp),%xmm1
        movapd nb312_izH2(%esp),%xmm2
        subsd  nb312_jxH2(%esp),%xmm0
        subsd  nb312_jyH2(%esp),%xmm1
        subsd  nb312_jzH2(%esp),%xmm2
        movapd %xmm0,nb312_dxH2H2(%esp)
        movapd %xmm1,nb312_dyH2H2(%esp)
        movapd %xmm2,nb312_dzH2H2(%esp)
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb312_rsqH2H2(%esp)

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
        movapd  nb312_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb312_half(%esp),%xmm3   ## iter1 
        mulsd   nb312_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb312_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb312_half(%esp),%xmm1   ## rinv 
        mulsd   nb312_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb312_rinvH2H2(%esp)
        movapd %xmm5,nb312_rinvH2H1(%esp)

        movapd nb312_rsqOO(%esp),%xmm0
        movapd nb312_rsqOH1(%esp),%xmm4
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
        movapd  nb312_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb312_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb312_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb312_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb312_half(%esp),%xmm1   ## rinv 
        mulsd   nb312_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb312_rinvOO(%esp)
        movapd %xmm5,nb312_rinvOH1(%esp)

        movapd nb312_rsqOH2(%esp),%xmm0
        movapd nb312_rsqH1O(%esp),%xmm4
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
        movapd  nb312_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb312_half(%esp),%xmm3   ## iter1 
        mulsd   nb312_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb312_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb312_half(%esp),%xmm1   ## rinv 
        mulsd   nb312_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb312_rinvOH2(%esp)
        movapd %xmm5,nb312_rinvH1O(%esp)

        movapd nb312_rsqH1H1(%esp),%xmm0
        movapd nb312_rsqH1H2(%esp),%xmm4
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
        movapd  nb312_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb312_half(%esp),%xmm3   ## iter1a 
        mulsd   nb312_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb312_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb312_half(%esp),%xmm1   ## rinv 
        mulsd   nb312_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb312_rinvH1H1(%esp)
        movapd %xmm5,nb312_rinvH1H2(%esp)

        movapd nb312_rsqH2O(%esp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb312_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb312_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb312_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb312_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb312_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb312_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb312_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulsd  nb312_tsc(%esp),%xmm1

        movd %eax,%mm0
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb312_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb312_two(%esp),%xmm7    ## two*Heps2 
        movapd nb312_qqOO(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul - then we can get rid of mm5 
    ## update vctot 
    addsd  nb312_vctot(%esp),%xmm5
        xorpd  %xmm2,%xmm2
    movlpd %xmm5,nb312_vctot(%esp)
        mulsd  nb312_tsc(%esp),%xmm3

        ## start doing lj 
        movapd %xmm0,%xmm2
        mulsd  %xmm2,%xmm2
        movapd %xmm2,%xmm1
        mulsd  %xmm2,%xmm1
        mulsd  %xmm2,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  nb312_c6(%esp),%xmm1
        mulsd  nb312_c12(%esp),%xmm2
        movapd %xmm2,%xmm4
        subsd  %xmm1,%xmm4
        addsd  nb312_Vvdwtot(%esp),%xmm4
        mulsd  nb312_six(%esp),%xmm1
        mulsd  nb312_twelve(%esp),%xmm2
        movlpd %xmm4,nb312_Vvdwtot(%esp)
        subsd  %xmm1,%xmm2
        mulsd  %xmm0,%xmm2

        subsd  %xmm3,%xmm2
        mulsd  %xmm2,%xmm0

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb312_dxOO(%esp),%xmm0
        mulpd nb312_dyOO(%esp),%xmm1
        mulpd nb312_dzOO(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb312_fixO(%esp),%xmm0
        addsd nb312_fiyO(%esp),%xmm1
        addsd nb312_fizO(%esp),%xmm2
        movlpd %xmm3,nb312_fjxO(%esp)
        movlpd %xmm4,nb312_fjyO(%esp)
        movlpd %xmm5,nb312_fjzO(%esp)
        movlpd %xmm0,nb312_fixO(%esp)
        movlpd %xmm1,nb312_fiyO(%esp)
        movlpd %xmm2,nb312_fizO(%esp)

        ## O-H1 interaction 
        movapd nb312_rinvOH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb312_rsqOH1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb312_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb312_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb312_two(%esp),%xmm7    ## two*Heps2 
        movapd nb312_qqOH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb312_vctot(%esp),%xmm5
    movlpd %xmm5,nb312_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb312_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb312_dxOH1(%esp),%xmm0
        mulsd nb312_dyOH1(%esp),%xmm1
        mulsd nb312_dzOH1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb312_fixO(%esp),%xmm0
        addsd nb312_fiyO(%esp),%xmm1
        addsd nb312_fizO(%esp),%xmm2
        movlpd %xmm3,nb312_fjxH1(%esp)
        movlpd %xmm4,nb312_fjyH1(%esp)
        movlpd %xmm5,nb312_fjzH1(%esp)
        movlpd %xmm0,nb312_fixO(%esp)
        movlpd %xmm1,nb312_fiyO(%esp)
        movlpd %xmm2,nb312_fizO(%esp)

        ## O-H2 interaction  
        movapd nb312_rinvOH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb312_rsqOH2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb312_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb312_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb312_two(%esp),%xmm7    ## two*Heps2 
        movapd nb312_qqOH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb312_vctot(%esp),%xmm5
    movlpd %xmm5,nb312_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb312_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb312_dxOH2(%esp),%xmm0
        mulsd nb312_dyOH2(%esp),%xmm1
        mulsd nb312_dzOH2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb312_fixO(%esp),%xmm0
        addsd nb312_fiyO(%esp),%xmm1
        addsd nb312_fizO(%esp),%xmm2
        movlpd %xmm3,nb312_fjxH2(%esp)
        movlpd %xmm4,nb312_fjyH2(%esp)
        movlpd %xmm5,nb312_fjzH2(%esp)
        movlpd %xmm0,nb312_fixO(%esp)
        movlpd %xmm1,nb312_fiyO(%esp)
        movlpd %xmm2,nb312_fizO(%esp)

        ## H1-O interaction 
        movapd nb312_rinvH1O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb312_rsqH1O(%esp),%xmm1   ## xmm1=r 
        mulsd  nb312_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb312_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb312_two(%esp),%xmm7    ## two*Heps2 
        movapd nb312_qqOH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb312_vctot(%esp),%xmm5
    movlpd %xmm5,nb312_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb312_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb312_fjxO(%esp),%xmm3
        movapd nb312_fjyO(%esp),%xmm4
        movapd nb312_fjzO(%esp),%xmm5
        mulsd nb312_dxH1O(%esp),%xmm0
        mulsd nb312_dyH1O(%esp),%xmm1
        mulsd nb312_dzH1O(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb312_fixH1(%esp),%xmm0
        addsd nb312_fiyH1(%esp),%xmm1
        addsd nb312_fizH1(%esp),%xmm2
        movlpd %xmm3,nb312_fjxO(%esp)
        movlpd %xmm4,nb312_fjyO(%esp)
        movlpd %xmm5,nb312_fjzO(%esp)
        movlpd %xmm0,nb312_fixH1(%esp)
        movlpd %xmm1,nb312_fiyH1(%esp)
        movlpd %xmm2,nb312_fizH1(%esp)

        ## H1-H1 interaction 
        movapd nb312_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb312_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb312_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb312_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb312_two(%esp),%xmm7    ## two*Heps2 
        movapd nb312_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb312_vctot(%esp),%xmm5
    movlpd %xmm5,nb312_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb312_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb312_fjxH1(%esp),%xmm3
        movapd nb312_fjyH1(%esp),%xmm4
        movapd nb312_fjzH1(%esp),%xmm5
        mulsd nb312_dxH1H1(%esp),%xmm0
        mulsd nb312_dyH1H1(%esp),%xmm1
        mulsd nb312_dzH1H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb312_fixH1(%esp),%xmm0
        addsd nb312_fiyH1(%esp),%xmm1
        addsd nb312_fizH1(%esp),%xmm2
        movlpd %xmm3,nb312_fjxH1(%esp)
        movlpd %xmm4,nb312_fjyH1(%esp)
        movlpd %xmm5,nb312_fjzH1(%esp)
        movlpd %xmm0,nb312_fixH1(%esp)
        movlpd %xmm1,nb312_fiyH1(%esp)
        movlpd %xmm2,nb312_fizH1(%esp)

        ## H1-H2 interaction 
        movapd nb312_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb312_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb312_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb312_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb312_two(%esp),%xmm7    ## two*Heps2 
        movapd nb312_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb312_vctot(%esp),%xmm5
    movlpd %xmm5,nb312_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb312_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb312_fjxH2(%esp),%xmm3
        movapd nb312_fjyH2(%esp),%xmm4
        movapd nb312_fjzH2(%esp),%xmm5
        mulsd nb312_dxH1H2(%esp),%xmm0
        mulsd nb312_dyH1H2(%esp),%xmm1
        mulsd nb312_dzH1H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb312_fixH1(%esp),%xmm0
        addsd nb312_fiyH1(%esp),%xmm1
        addsd nb312_fizH1(%esp),%xmm2
        movlpd %xmm3,nb312_fjxH2(%esp)
        movlpd %xmm4,nb312_fjyH2(%esp)
        movlpd %xmm5,nb312_fjzH2(%esp)
        movlpd %xmm0,nb312_fixH1(%esp)
        movlpd %xmm1,nb312_fiyH1(%esp)
        movlpd %xmm2,nb312_fizH1(%esp)

        ## H2-O interaction 
        movapd nb312_rinvH2O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb312_rsqH2O(%esp),%xmm1   ## xmm1=r 
        mulsd  nb312_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb312_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb312_two(%esp),%xmm7    ## two*Heps2 
        movapd nb312_qqOH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb312_vctot(%esp),%xmm5
    movlpd %xmm5,nb312_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb312_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb312_fjxO(%esp),%xmm3
        movapd nb312_fjyO(%esp),%xmm4
        movapd nb312_fjzO(%esp),%xmm5
        mulsd nb312_dxH2O(%esp),%xmm0
        mulsd nb312_dyH2O(%esp),%xmm1
        mulsd nb312_dzH2O(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb312_fixH2(%esp),%xmm0
        addsd nb312_fiyH2(%esp),%xmm1
        addsd nb312_fizH2(%esp),%xmm2
        movlpd %xmm3,nb312_fjxO(%esp)
        movlpd %xmm4,nb312_fjyO(%esp)
        movlpd %xmm5,nb312_fjzO(%esp)
        movlpd %xmm0,nb312_fixH2(%esp)
        movlpd %xmm1,nb312_fiyH2(%esp)
        movlpd %xmm2,nb312_fizH2(%esp)

        ## H2-H1 interaction 
        movapd nb312_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb312_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb312_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb312_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb312_two(%esp),%xmm7    ## two*Heps2 
        movapd nb312_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb312_vctot(%esp),%xmm5
    movlpd %xmm5,nb312_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb312_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb312_fjxH1(%esp),%xmm3
        movapd nb312_fjyH1(%esp),%xmm4
        movapd nb312_fjzH1(%esp),%xmm5
        mulsd nb312_dxH2H1(%esp),%xmm0
        mulsd nb312_dyH2H1(%esp),%xmm1
        mulsd nb312_dzH2H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb312_fixH2(%esp),%xmm0
        addsd nb312_fiyH2(%esp),%xmm1
        addsd nb312_fizH2(%esp),%xmm2
        movlpd %xmm3,nb312_fjxH1(%esp)
        movlpd %xmm4,nb312_fjyH1(%esp)
        movlpd %xmm5,nb312_fjzH1(%esp)
        movlpd %xmm0,nb312_fixH2(%esp)
        movlpd %xmm1,nb312_fiyH2(%esp)
        movlpd %xmm2,nb312_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb312_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb312_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb312_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb312_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb312_two(%esp),%xmm7    ## two*Heps2 
        movapd nb312_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb312_vctot(%esp),%xmm5
    movlpd %xmm5,nb312_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb312_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb312_fjxH2(%esp),%xmm3
        movapd nb312_fjyH2(%esp),%xmm4
        movapd nb312_fjzH2(%esp),%xmm5
        mulsd nb312_dxH2H2(%esp),%xmm0
        mulsd nb312_dyH2H2(%esp),%xmm1
        mulsd nb312_dzH2H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb312_fixH2(%esp),%xmm0
        addsd nb312_fiyH2(%esp),%xmm1
        addsd nb312_fizH2(%esp),%xmm2
        movlpd %xmm3,nb312_fjxH2(%esp)
        movlpd %xmm4,nb312_fjyH2(%esp)
        movlpd %xmm5,nb312_fjzH2(%esp)
        movlpd %xmm0,nb312_fixH2(%esp)
        movlpd %xmm1,nb312_fiyH2(%esp)
        movlpd %xmm2,nb312_fizH2(%esp)

        movl nb312_faction(%ebp),%edi

        movd %mm0,%eax

        ## Did all interactions - now update j forces 
        movlpd (%edi,%eax,8),%xmm0
        movlpd 8(%edi,%eax,8),%xmm1
        movlpd 16(%edi,%eax,8),%xmm2
        movlpd 24(%edi,%eax,8),%xmm3
        movlpd 32(%edi,%eax,8),%xmm4
        movlpd 40(%edi,%eax,8),%xmm5
        movlpd 48(%edi,%eax,8),%xmm6
        movlpd 56(%edi,%eax,8),%xmm7
        addsd nb312_fjxO(%esp),%xmm0
        addsd nb312_fjyO(%esp),%xmm1
        addsd nb312_fjzO(%esp),%xmm2
        addsd nb312_fjxH1(%esp),%xmm3
        addsd nb312_fjyH1(%esp),%xmm4
        addsd nb312_fjzH1(%esp),%xmm5
        addsd nb312_fjxH2(%esp),%xmm6
        addsd nb312_fjyH2(%esp),%xmm7
        movlpd %xmm0,(%edi,%eax,8)
        movlpd %xmm1,8(%edi,%eax,8)
        movlpd %xmm2,16(%edi,%eax,8)
        movlpd %xmm3,24(%edi,%eax,8)
        movlpd %xmm4,32(%edi,%eax,8)
        movlpd %xmm5,40(%edi,%eax,8)
        movlpd %xmm6,48(%edi,%eax,8)
        movlpd %xmm7,56(%edi,%eax,8)

        movlpd 64(%edi,%eax,8),%xmm0
        addsd nb312_fjzH2(%esp),%xmm0
        movlpd %xmm0,64(%edi,%eax,8)

_nb_kernel312_ia32_sse2.nb312_updateouterdata: 
        movl  nb312_ii3(%esp),%ecx
        movl  nb312_faction(%ebp),%edi
        movl  nb312_fshift(%ebp),%esi
        movl  nb312_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb312_fixO(%esp),%xmm0
        movapd nb312_fiyO(%esp),%xmm1
        movapd nb312_fizO(%esp),%xmm2

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
        movsd  (%edi,%ecx,8),%xmm3
        movsd  8(%edi,%ecx,8),%xmm4
        movsd  16(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,(%edi,%ecx,8)
        movsd  %xmm4,8(%edi,%ecx,8)
        movsd  %xmm5,16(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        movapd %xmm0,%xmm6
        movsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm6

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb312_fixH1(%esp),%xmm0
        movapd nb312_fiyH1(%esp),%xmm1
        movapd nb312_fizH1(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  24(%edi,%ecx,8),%xmm3
        movsd  32(%edi,%ecx,8),%xmm4
        movsd  40(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,24(%edi,%ecx,8)
        movsd  %xmm4,32(%edi,%ecx,8)
        movsd  %xmm5,40(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb312_fixH2(%esp),%xmm0
        movapd nb312_fiyH2(%esp),%xmm1
        movapd nb312_fizH2(%esp),%xmm2

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
        movsd  48(%edi,%ecx,8),%xmm3
        movsd  56(%edi,%ecx,8),%xmm4
        movsd  64(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,48(%edi,%ecx,8)
        movsd  %xmm4,56(%edi,%ecx,8)
        movsd  %xmm5,64(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## increment fshift force 
        movlpd (%esi,%edx,8),%xmm3
        movhpd 8(%esi,%edx,8),%xmm3
        movsd  16(%esi,%edx,8),%xmm4
        addpd  %xmm6,%xmm3
        addsd  %xmm7,%xmm4
        movlpd %xmm3,(%esi,%edx,8)
        movhpd %xmm3,8(%esi,%edx,8)
        movsd  %xmm4,16(%esi,%edx,8)

        ## get n from stack
        movl nb312_n(%esp),%esi
        ## get group index for i particle 
        movl  nb312_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb312_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb312_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb312_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb312_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb312_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel312_ia32_sse2.nb312_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb312_n(%esp)
        jmp _nb_kernel312_ia32_sse2.nb312_outer
_nb_kernel312_ia32_sse2.nb312_outerend: 
        ## check if more outer neighborlists remain
        movl  nb312_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel312_ia32_sse2.nb312_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel312_ia32_sse2.nb312_threadloop
_nb_kernel312_ia32_sse2.nb312_end: 
        emms

        movl nb312_nouter(%esp),%eax
        movl nb312_ninner(%esp),%ebx
        movl nb312_outeriter(%ebp),%ecx
        movl nb312_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb312_salign(%esp),%eax
        addl %eax,%esp
        addl $1560,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret






.globl nb_kernel312nf_ia32_sse2
.globl _nb_kernel312nf_ia32_sse2
nb_kernel312nf_ia32_sse2:       
_nb_kernel312nf_ia32_sse2:      
.set nb312nf_p_nri, 8
.set nb312nf_iinr, 12
.set nb312nf_jindex, 16
.set nb312nf_jjnr, 20
.set nb312nf_shift, 24
.set nb312nf_shiftvec, 28
.set nb312nf_fshift, 32
.set nb312nf_gid, 36
.set nb312nf_pos, 40
.set nb312nf_faction, 44
.set nb312nf_charge, 48
.set nb312nf_p_facel, 52
.set nb312nf_argkrf, 56
.set nb312nf_argcrf, 60
.set nb312nf_Vc, 64
.set nb312nf_type, 68
.set nb312nf_p_ntype, 72
.set nb312nf_vdwparam, 76
.set nb312nf_Vvdw, 80
.set nb312nf_p_tabscale, 84
.set nb312nf_VFtab, 88
.set nb312nf_invsqrta, 92
.set nb312nf_dvda, 96
.set nb312nf_p_gbtabscale, 100
.set nb312nf_GBtab, 104
.set nb312nf_p_nthreads, 108
.set nb312nf_count, 112
.set nb312nf_mtx, 116
.set nb312nf_outeriter, 120
.set nb312nf_inneriter, 124
.set nb312nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb312nf_ixO, 0
.set nb312nf_iyO, 16
.set nb312nf_izO, 32
.set nb312nf_ixH1, 48
.set nb312nf_iyH1, 64
.set nb312nf_izH1, 80
.set nb312nf_ixH2, 96
.set nb312nf_iyH2, 112
.set nb312nf_izH2, 128
.set nb312nf_jxO, 144
.set nb312nf_jyO, 160
.set nb312nf_jzO, 176
.set nb312nf_jxH1, 192
.set nb312nf_jyH1, 208
.set nb312nf_jzH1, 224
.set nb312nf_jxH2, 240
.set nb312nf_jyH2, 256
.set nb312nf_jzH2, 272
.set nb312nf_qqOO, 288
.set nb312nf_qqOH, 304
.set nb312nf_qqHH, 320
.set nb312nf_tsc, 336
.set nb312nf_c6, 352
.set nb312nf_c12, 368
.set nb312nf_vctot, 384
.set nb312nf_Vvdwtot, 400
.set nb312nf_half, 416
.set nb312nf_three, 432
.set nb312nf_rsqOO, 448
.set nb312nf_rsqOH1, 464
.set nb312nf_rsqOH2, 480
.set nb312nf_rsqH1O, 496
.set nb312nf_rsqH1H1, 512
.set nb312nf_rsqH1H2, 528
.set nb312nf_rsqH2O, 544
.set nb312nf_rsqH2H1, 560
.set nb312nf_rsqH2H2, 576
.set nb312nf_rinvOO, 592
.set nb312nf_rinvOH1, 608
.set nb312nf_rinvOH2, 624
.set nb312nf_rinvH1O, 640
.set nb312nf_rinvH1H1, 656
.set nb312nf_rinvH1H2, 672
.set nb312nf_rinvH2O, 688
.set nb312nf_rinvH2H1, 704
.set nb312nf_rinvH2H2, 720
.set nb312nf_is3, 736
.set nb312nf_ii3, 740
.set nb312nf_innerjjnr, 744
.set nb312nf_innerk, 748
.set nb312nf_n, 752
.set nb312nf_nn1, 756
.set nb312nf_nri, 760
.set nb312nf_nouter, 764
.set nb312nf_ninner, 768
.set nb312nf_salign, 772
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $776,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb312nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb312nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb312nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb312nf_nouter(%esp)
        movl %eax,nb312nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb312nf_half(%esp)
        movl %ebx,nb312nf_half+4(%esp)
        movsd nb312nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb312nf_half(%esp)
        movapd %xmm3,nb312nf_three(%esp)

        movl nb312nf_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm5
        shufpd $0,%xmm5,%xmm5
        movapd %xmm5,nb312nf_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb312nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb312nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb312nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb312nf_qqOO(%esp)
        movapd %xmm4,nb312nf_qqOH(%esp)
        movapd %xmm5,nb312nf_qqHH(%esp)

        xorpd %xmm0,%xmm0
        movl  nb312nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb312nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb312nf_vdwparam(%ebp),%eax
        movlpd (%eax,%edx,8),%xmm0
        movlpd 8(%eax,%edx,8),%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb312nf_c6(%esp)
        movapd %xmm1,nb312nf_c12(%esp)

_nb_kernel312nf_ia32_sse2.nb312nf_threadloop: 
        movl  nb312nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel312nf_ia32_sse2.nb312nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel312nf_ia32_sse2.nb312nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb312nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb312nf_n(%esp)
        movl %ebx,nb312nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel312nf_ia32_sse2.nb312nf_outerstart
        jmp _nb_kernel312nf_ia32_sse2.nb312nf_end

_nb_kernel312nf_ia32_sse2.nb312nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb312nf_nouter(%esp),%ebx
        movl %ebx,nb312nf_nouter(%esp)

_nb_kernel312nf_ia32_sse2.nb312nf_outer: 
        movl  nb312nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb312nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb312nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb312nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb312nf_ii3(%esp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb312nf_ixO(%esp)
        movapd %xmm4,nb312nf_iyO(%esp)
        movapd %xmm5,nb312nf_izO(%esp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%eax,%ebx,8),%xmm0
        addsd 32(%eax,%ebx,8),%xmm1
        addsd 40(%eax,%ebx,8),%xmm2
        addsd 48(%eax,%ebx,8),%xmm3
        addsd 56(%eax,%ebx,8),%xmm4
        addsd 64(%eax,%ebx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb312nf_ixH1(%esp)
        movapd %xmm1,nb312nf_iyH1(%esp)
        movapd %xmm2,nb312nf_izH1(%esp)
        movapd %xmm3,nb312nf_ixH2(%esp)
        movapd %xmm4,nb312nf_iyH2(%esp)
        movapd %xmm5,nb312nf_izH2(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb312nf_vctot(%esp)
        movapd %xmm4,nb312nf_Vvdwtot(%esp)

        movl  nb312nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb312nf_pos(%ebp),%esi
        movl  nb312nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb312nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb312nf_ninner(%esp),%ecx
        movl  %ecx,nb312nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb312nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel312nf_ia32_sse2.nb312nf_unroll_loop
        jmp   _nb_kernel312nf_ia32_sse2.nb312nf_checksingle
_nb_kernel312nf_ia32_sse2.nb312nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb312nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb312nf_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb312nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move j coordinates to local temp variables 
        movlpd (%esi,%eax,8),%xmm2
        movlpd 8(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%eax,8),%xmm4
        movlpd 24(%esi,%eax,8),%xmm5
        movlpd 32(%esi,%eax,8),%xmm6
        movlpd 40(%esi,%eax,8),%xmm7
        movhpd (%esi,%ebx,8),%xmm2
        movhpd 8(%esi,%ebx,8),%xmm3
        movhpd 16(%esi,%ebx,8),%xmm4
        movhpd 24(%esi,%ebx,8),%xmm5
        movhpd 32(%esi,%ebx,8),%xmm6
        movhpd 40(%esi,%ebx,8),%xmm7
        movapd  %xmm2,nb312nf_jxO(%esp)
        movapd  %xmm3,nb312nf_jyO(%esp)
        movapd  %xmm4,nb312nf_jzO(%esp)
        movapd  %xmm5,nb312nf_jxH1(%esp)
        movapd  %xmm6,nb312nf_jyH1(%esp)
        movapd  %xmm7,nb312nf_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 48(%esi,%ebx,8),%xmm2
        movhpd 56(%esi,%ebx,8),%xmm3
        movhpd 64(%esi,%ebx,8),%xmm4
        movapd  %xmm2,nb312nf_jxH2(%esp)
        movapd  %xmm3,nb312nf_jyH2(%esp)
        movapd  %xmm4,nb312nf_jzH2(%esp)

        movapd nb312nf_ixO(%esp),%xmm0
        movapd nb312nf_iyO(%esp),%xmm1
        movapd nb312nf_izO(%esp),%xmm2
        movapd nb312nf_ixO(%esp),%xmm3
        movapd nb312nf_iyO(%esp),%xmm4
        movapd nb312nf_izO(%esp),%xmm5
        subpd  nb312nf_jxO(%esp),%xmm0
        subpd  nb312nf_jyO(%esp),%xmm1
        subpd  nb312nf_jzO(%esp),%xmm2
        subpd  nb312nf_jxH1(%esp),%xmm3
        subpd  nb312nf_jyH1(%esp),%xmm4
        subpd  nb312nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb312nf_rsqOO(%esp)
        movapd %xmm3,nb312nf_rsqOH1(%esp)

        movapd nb312nf_ixO(%esp),%xmm0
        movapd nb312nf_iyO(%esp),%xmm1
        movapd nb312nf_izO(%esp),%xmm2
        movapd nb312nf_ixH1(%esp),%xmm3
        movapd nb312nf_iyH1(%esp),%xmm4
        movapd nb312nf_izH1(%esp),%xmm5
        subpd  nb312nf_jxH2(%esp),%xmm0
        subpd  nb312nf_jyH2(%esp),%xmm1
        subpd  nb312nf_jzH2(%esp),%xmm2
        subpd  nb312nf_jxO(%esp),%xmm3
        subpd  nb312nf_jyO(%esp),%xmm4
        subpd  nb312nf_jzO(%esp),%xmm5
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
        movapd %xmm0,nb312nf_rsqOH2(%esp)
        movapd %xmm3,nb312nf_rsqH1O(%esp)

        movapd nb312nf_ixH1(%esp),%xmm0
        movapd nb312nf_iyH1(%esp),%xmm1
        movapd nb312nf_izH1(%esp),%xmm2
        movapd nb312nf_ixH1(%esp),%xmm3
        movapd nb312nf_iyH1(%esp),%xmm4
        movapd nb312nf_izH1(%esp),%xmm5
        subpd  nb312nf_jxH1(%esp),%xmm0
        subpd  nb312nf_jyH1(%esp),%xmm1
        subpd  nb312nf_jzH1(%esp),%xmm2
        subpd  nb312nf_jxH2(%esp),%xmm3
        subpd  nb312nf_jyH2(%esp),%xmm4
        subpd  nb312nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb312nf_rsqH1H1(%esp)
        movapd %xmm3,nb312nf_rsqH1H2(%esp)

        movapd nb312nf_ixH2(%esp),%xmm0
        movapd nb312nf_iyH2(%esp),%xmm1
        movapd nb312nf_izH2(%esp),%xmm2
        movapd nb312nf_ixH2(%esp),%xmm3
        movapd nb312nf_iyH2(%esp),%xmm4
        movapd nb312nf_izH2(%esp),%xmm5
        subpd  nb312nf_jxO(%esp),%xmm0
        subpd  nb312nf_jyO(%esp),%xmm1
        subpd  nb312nf_jzO(%esp),%xmm2
        subpd  nb312nf_jxH1(%esp),%xmm3
        subpd  nb312nf_jyH1(%esp),%xmm4
        subpd  nb312nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb312nf_rsqH2O(%esp)
        movapd %xmm4,nb312nf_rsqH2H1(%esp)

        movapd nb312nf_ixH2(%esp),%xmm0
        movapd nb312nf_iyH2(%esp),%xmm1
        movapd nb312nf_izH2(%esp),%xmm2
        subpd  nb312nf_jxH2(%esp),%xmm0
        subpd  nb312nf_jyH2(%esp),%xmm1
        subpd  nb312nf_jzH2(%esp),%xmm2
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb312nf_rsqH2H2(%esp)

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
        movapd  nb312nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb312nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb312nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb312nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb312nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb312nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb312nf_rinvH2H2(%esp)
        movapd %xmm5,nb312nf_rinvH2H1(%esp)

        movapd nb312nf_rsqOO(%esp),%xmm0
        movapd nb312nf_rsqOH1(%esp),%xmm4
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
        movapd  nb312nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb312nf_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb312nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb312nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb312nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb312nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb312nf_rinvOO(%esp)
        movapd %xmm5,nb312nf_rinvOH1(%esp)

        movapd nb312nf_rsqOH2(%esp),%xmm0
        movapd nb312nf_rsqH1O(%esp),%xmm4
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
        movapd  nb312nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb312nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb312nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb312nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb312nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb312nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb312nf_rinvOH2(%esp)
        movapd %xmm5,nb312nf_rinvH1O(%esp)

        movapd nb312nf_rsqH1H1(%esp),%xmm0
        movapd nb312nf_rsqH1H2(%esp),%xmm4
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
        movapd  nb312nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb312nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb312nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb312nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb312nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb312nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb312nf_rinvH1H1(%esp)
        movapd %xmm5,nb312nf_rinvH1H2(%esp)

        movapd nb312nf_rsqH2O(%esp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb312nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb312nf_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb312nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb312nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb312nf_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb312nf_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb312nf_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulpd  nb312nf_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb312nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb312nf_qqOO(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul - then we can get rid of mm5 
    ## update vctot 
    addpd  nb312nf_vctot(%esp),%xmm5
    movapd %xmm5,nb312nf_vctot(%esp)

        ## start doing lj 
        movapd %xmm0,%xmm2
        mulpd  %xmm2,%xmm2
        movapd %xmm2,%xmm1
        mulpd  %xmm2,%xmm1
        mulpd  %xmm2,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulpd  nb312nf_c6(%esp),%xmm1
        mulpd  nb312nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm4
        subpd  %xmm1,%xmm4
        addpd  nb312nf_Vvdwtot(%esp),%xmm4
        movapd %xmm4,nb312nf_Vvdwtot(%esp)

        ## O-H1 interaction 
        movapd nb312nf_rinvOH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb312nf_rsqOH1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb312nf_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb312nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb312nf_qqOH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb312nf_vctot(%esp),%xmm5
    movapd %xmm5,nb312nf_vctot(%esp)

        ## O-H2 interaction  
        movapd nb312nf_rinvOH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb312nf_rsqOH2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb312nf_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb312nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb312nf_qqOH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb312nf_vctot(%esp),%xmm5
    movapd %xmm5,nb312nf_vctot(%esp)

        ## H1-O interaction 
        movapd nb312nf_rinvH1O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb312nf_rsqH1O(%esp),%xmm1   ## xmm1=r 
        mulpd  nb312nf_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb312nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb312nf_qqOH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb312nf_vctot(%esp),%xmm5
    movapd %xmm5,nb312nf_vctot(%esp)

        ## H1-H1 interaction 
        movapd nb312nf_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb312nf_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb312nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb312nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb312nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb312nf_vctot(%esp),%xmm5
    movapd %xmm5,nb312nf_vctot(%esp)

        ## H1-H2 interaction 
        movapd nb312nf_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb312nf_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb312nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb312nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb312nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb312nf_vctot(%esp),%xmm5
    movapd %xmm5,nb312nf_vctot(%esp)

        ## H2-O interaction 
        movapd nb312nf_rinvH2O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb312nf_rsqH2O(%esp),%xmm1   ## xmm1=r 
        mulpd  nb312nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb312nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb312nf_qqOH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb312nf_vctot(%esp),%xmm5
    movapd %xmm5,nb312nf_vctot(%esp)

        ## H2-H1 interaction 
        movapd nb312nf_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb312nf_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb312nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb312nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb312nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb312nf_vctot(%esp),%xmm5
    movapd %xmm5,nb312nf_vctot(%esp)

        ## H2-H2 interaction 
        movapd nb312nf_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb312nf_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb312nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb312nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb312nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb312nf_vctot(%esp),%xmm5
    movapd %xmm5,nb312nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb312nf_innerk(%esp)
        jl    _nb_kernel312nf_ia32_sse2.nb312nf_checksingle
        jmp   _nb_kernel312nf_ia32_sse2.nb312nf_unroll_loop
_nb_kernel312nf_ia32_sse2.nb312nf_checksingle: 
        movl  nb312nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel312nf_ia32_sse2.nb312nf_dosingle
        jmp   _nb_kernel312nf_ia32_sse2.nb312nf_updateouterdata
_nb_kernel312nf_ia32_sse2.nb312nf_dosingle: 
        movl  nb312nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb312nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        movlpd (%esi,%eax,8),%xmm2
        movlpd 8(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%eax,8),%xmm4
        movlpd 24(%esi,%eax,8),%xmm5
        movlpd 32(%esi,%eax,8),%xmm6
        movlpd 40(%esi,%eax,8),%xmm7
        movapd  %xmm2,nb312nf_jxO(%esp)
        movapd  %xmm3,nb312nf_jyO(%esp)
        movapd  %xmm4,nb312nf_jzO(%esp)
        movapd  %xmm5,nb312nf_jxH1(%esp)
        movapd  %xmm6,nb312nf_jyH1(%esp)
        movapd  %xmm7,nb312nf_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movapd  %xmm2,nb312nf_jxH2(%esp)
        movapd  %xmm3,nb312nf_jyH2(%esp)
        movapd  %xmm4,nb312nf_jzH2(%esp)

        movapd nb312nf_ixO(%esp),%xmm0
        movapd nb312nf_iyO(%esp),%xmm1
        movapd nb312nf_izO(%esp),%xmm2
        movapd nb312nf_ixO(%esp),%xmm3
        movapd nb312nf_iyO(%esp),%xmm4
        movapd nb312nf_izO(%esp),%xmm5
        subsd  nb312nf_jxO(%esp),%xmm0
        subsd  nb312nf_jyO(%esp),%xmm1
        subsd  nb312nf_jzO(%esp),%xmm2
        subsd  nb312nf_jxH1(%esp),%xmm3
        subsd  nb312nf_jyH1(%esp),%xmm4
        subsd  nb312nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb312nf_rsqOO(%esp)
        movapd %xmm3,nb312nf_rsqOH1(%esp)

        movapd nb312nf_ixO(%esp),%xmm0
        movapd nb312nf_iyO(%esp),%xmm1
        movapd nb312nf_izO(%esp),%xmm2
        movapd nb312nf_ixH1(%esp),%xmm3
        movapd nb312nf_iyH1(%esp),%xmm4
        movapd nb312nf_izH1(%esp),%xmm5
        subsd  nb312nf_jxH2(%esp),%xmm0
        subsd  nb312nf_jyH2(%esp),%xmm1
        subsd  nb312nf_jzH2(%esp),%xmm2
        subsd  nb312nf_jxO(%esp),%xmm3
        subsd  nb312nf_jyO(%esp),%xmm4
        subsd  nb312nf_jzO(%esp),%xmm5
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
        movapd %xmm0,nb312nf_rsqOH2(%esp)
        movapd %xmm3,nb312nf_rsqH1O(%esp)

        movapd nb312nf_ixH1(%esp),%xmm0
        movapd nb312nf_iyH1(%esp),%xmm1
        movapd nb312nf_izH1(%esp),%xmm2
        movapd nb312nf_ixH1(%esp),%xmm3
        movapd nb312nf_iyH1(%esp),%xmm4
        movapd nb312nf_izH1(%esp),%xmm5
        subsd  nb312nf_jxH1(%esp),%xmm0
        subsd  nb312nf_jyH1(%esp),%xmm1
        subsd  nb312nf_jzH1(%esp),%xmm2
        subsd  nb312nf_jxH2(%esp),%xmm3
        subsd  nb312nf_jyH2(%esp),%xmm4
        subsd  nb312nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb312nf_rsqH1H1(%esp)
        movapd %xmm3,nb312nf_rsqH1H2(%esp)

        movapd nb312nf_ixH2(%esp),%xmm0
        movapd nb312nf_iyH2(%esp),%xmm1
        movapd nb312nf_izH2(%esp),%xmm2
        movapd nb312nf_ixH2(%esp),%xmm3
        movapd nb312nf_iyH2(%esp),%xmm4
        movapd nb312nf_izH2(%esp),%xmm5
        subsd  nb312nf_jxO(%esp),%xmm0
        subsd  nb312nf_jyO(%esp),%xmm1
        subsd  nb312nf_jzO(%esp),%xmm2
        subsd  nb312nf_jxH1(%esp),%xmm3
        subsd  nb312nf_jyH1(%esp),%xmm4
        subsd  nb312nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb312nf_rsqH2O(%esp)
        movapd %xmm4,nb312nf_rsqH2H1(%esp)

        movapd nb312nf_ixH2(%esp),%xmm0
        movapd nb312nf_iyH2(%esp),%xmm1
        movapd nb312nf_izH2(%esp),%xmm2
        subsd  nb312nf_jxH2(%esp),%xmm0
        subsd  nb312nf_jyH2(%esp),%xmm1
        subsd  nb312nf_jzH2(%esp),%xmm2
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb312nf_rsqH2H2(%esp)

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
        movapd  nb312nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb312nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb312nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb312nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb312nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb312nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb312nf_rinvH2H2(%esp)
        movapd %xmm5,nb312nf_rinvH2H1(%esp)

        movapd nb312nf_rsqOO(%esp),%xmm0
        movapd nb312nf_rsqOH1(%esp),%xmm4
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
        movapd  nb312nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb312nf_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb312nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb312nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb312nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb312nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb312nf_rinvOO(%esp)
        movapd %xmm5,nb312nf_rinvOH1(%esp)

        movapd nb312nf_rsqOH2(%esp),%xmm0
        movapd nb312nf_rsqH1O(%esp),%xmm4
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
        movapd  nb312nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb312nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb312nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb312nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb312nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb312nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb312nf_rinvOH2(%esp)
        movapd %xmm5,nb312nf_rinvH1O(%esp)

        movapd nb312nf_rsqH1H1(%esp),%xmm0
        movapd nb312nf_rsqH1H2(%esp),%xmm4
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
        movapd  nb312nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb312nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb312nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb312nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb312nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb312nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb312nf_rinvH1H1(%esp)
        movapd %xmm5,nb312nf_rinvH1H2(%esp)

        movapd nb312nf_rsqH2O(%esp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb312nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb312nf_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb312nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb312nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb312nf_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb312nf_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb312nf_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulsd  nb312nf_tsc(%esp),%xmm1

        movd %eax,%mm0
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb312nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb312nf_qqOO(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul - then we can get rid of mm5 
    ## update vctot 
    addsd  nb312nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb312nf_vctot(%esp)

        ## start doing lj 
        movapd %xmm0,%xmm2
        mulsd  %xmm2,%xmm2
        movapd %xmm2,%xmm1
        mulsd  %xmm2,%xmm1
        mulsd  %xmm2,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  nb312nf_c6(%esp),%xmm1
        mulsd  nb312nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm4
        subsd  %xmm1,%xmm4
        addsd  nb312nf_Vvdwtot(%esp),%xmm4
        movlpd %xmm4,nb312nf_Vvdwtot(%esp)

        ## O-H1 interaction 
        movapd nb312nf_rinvOH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb312nf_rsqOH1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb312nf_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb312nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb312nf_qqOH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV 
    ## at this point mm5 contains vcoul 

    addsd  nb312nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb312nf_vctot(%esp)

        ## O-H2 interaction  
        movapd nb312nf_rinvOH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb312nf_rsqOH2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb312nf_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb312nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb312nf_qqOH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb312nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb312nf_vctot(%esp)

        ## H1-O interaction 
        movapd nb312nf_rinvH1O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb312nf_rsqH1O(%esp),%xmm1   ## xmm1=r 
        mulsd  nb312nf_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb312nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb312nf_qqOH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb312nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb312nf_vctot(%esp)

        ## H1-H1 interaction 
        movapd nb312nf_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb312nf_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb312nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb312nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb312nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb312nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb312nf_vctot(%esp)

        ## H1-H2 interaction 
        movapd nb312nf_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb312nf_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb312nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb312nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb312nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb312nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb312nf_vctot(%esp)

        ## H2-O interaction 
        movapd nb312nf_rinvH2O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb312nf_rsqH2O(%esp),%xmm1   ## xmm1=r 
        mulsd  nb312nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb312nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb312nf_qqOH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb312nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb312nf_vctot(%esp)

        ## H2-H1 interaction 
        movapd nb312nf_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb312nf_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb312nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb312nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb312nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb312nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb312nf_vctot(%esp)

        ## H2-H2 interaction 
        movapd nb312nf_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb312nf_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb312nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb312nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb312nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb312nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb312nf_vctot(%esp)

_nb_kernel312nf_ia32_sse2.nb312nf_updateouterdata: 
        ## get n from stack
        movl nb312nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb312nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb312nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb312nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb312nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb312nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb312nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel312nf_ia32_sse2.nb312nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb312nf_n(%esp)
        jmp _nb_kernel312nf_ia32_sse2.nb312nf_outer
_nb_kernel312nf_ia32_sse2.nb312nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb312nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel312nf_ia32_sse2.nb312nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel312nf_ia32_sse2.nb312nf_threadloop
_nb_kernel312nf_ia32_sse2.nb312nf_end: 
        emms

        movl nb312nf_nouter(%esp),%eax
        movl nb312nf_ninner(%esp),%ebx
        movl nb312nf_outeriter(%ebp),%ecx
        movl nb312nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb312nf_salign(%esp),%eax
        addl %eax,%esp
        addl $776,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret

