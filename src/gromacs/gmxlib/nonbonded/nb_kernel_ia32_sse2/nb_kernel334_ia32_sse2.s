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


.globl nb_kernel334_ia32_sse2
.globl _nb_kernel334_ia32_sse2
nb_kernel334_ia32_sse2: 
_nb_kernel334_ia32_sse2:        
.set nb334_p_nri, 8
.set nb334_iinr, 12
.set nb334_jindex, 16
.set nb334_jjnr, 20
.set nb334_shift, 24
.set nb334_shiftvec, 28
.set nb334_fshift, 32
.set nb334_gid, 36
.set nb334_pos, 40
.set nb334_faction, 44
.set nb334_charge, 48
.set nb334_p_facel, 52
.set nb334_argkrf, 56
.set nb334_argcrf, 60
.set nb334_Vc, 64
.set nb334_type, 68
.set nb334_p_ntype, 72
.set nb334_vdwparam, 76
.set nb334_Vvdw, 80
.set nb334_p_tabscale, 84
.set nb334_VFtab, 88
.set nb334_invsqrta, 92
.set nb334_dvda, 96
.set nb334_p_gbtabscale, 100
.set nb334_GBtab, 104
.set nb334_p_nthreads, 108
.set nb334_count, 112
.set nb334_mtx, 116
.set nb334_outeriter, 120
.set nb334_inneriter, 124
.set nb334_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb334_ixO, 0
.set nb334_iyO, 16
.set nb334_izO, 32
.set nb334_ixH1, 48
.set nb334_iyH1, 64
.set nb334_izH1, 80
.set nb334_ixH2, 96
.set nb334_iyH2, 112
.set nb334_izH2, 128
.set nb334_ixM, 144
.set nb334_iyM, 160
.set nb334_izM, 176
.set nb334_jxO, 192
.set nb334_jyO, 208
.set nb334_jzO, 224
.set nb334_jxH1, 240
.set nb334_jyH1, 256
.set nb334_jzH1, 272
.set nb334_jxH2, 288
.set nb334_jyH2, 304
.set nb334_jzH2, 320
.set nb334_jxM, 336
.set nb334_jyM, 352
.set nb334_jzM, 368
.set nb334_dxOO, 384
.set nb334_dyOO, 400
.set nb334_dzOO, 416
.set nb334_dxH1H1, 432
.set nb334_dyH1H1, 448
.set nb334_dzH1H1, 464
.set nb334_dxH1H2, 480
.set nb334_dyH1H2, 496
.set nb334_dzH1H2, 512
.set nb334_dxH1M, 528
.set nb334_dyH1M, 544
.set nb334_dzH1M, 560
.set nb334_dxH2H1, 576
.set nb334_dyH2H1, 592
.set nb334_dzH2H1, 608
.set nb334_dxH2H2, 624
.set nb334_dyH2H2, 640
.set nb334_dzH2H2, 656
.set nb334_dxH2M, 672
.set nb334_dyH2M, 688
.set nb334_dzH2M, 704
.set nb334_dxMH1, 720
.set nb334_dyMH1, 736
.set nb334_dzMH1, 752
.set nb334_dxMH2, 768
.set nb334_dyMH2, 784
.set nb334_dzMH2, 800
.set nb334_dxMM, 816
.set nb334_dyMM, 832
.set nb334_dzMM, 848
.set nb334_qqMM, 864
.set nb334_qqMH, 880
.set nb334_qqHH, 896
.set nb334_two, 912
.set nb334_tsc, 928
.set nb334_c6, 944
.set nb334_c12, 960
.set nb334_vctot, 976
.set nb334_Vvdwtot, 992
.set nb334_fixO, 1008
.set nb334_fiyO, 1024
.set nb334_fizO, 1040
.set nb334_fixH1, 1056
.set nb334_fiyH1, 1072
.set nb334_fizH1, 1088
.set nb334_fixH2, 1104
.set nb334_fiyH2, 1120
.set nb334_fizH2, 1136
.set nb334_fixM, 1152
.set nb334_fiyM, 1168
.set nb334_fizM, 1184
.set nb334_fjxO, 1200
.set nb334_fjyO, 1216
.set nb334_fjzO, 1232
.set nb334_fjxH1, 1248
.set nb334_fjyH1, 1264
.set nb334_fjzH1, 1280
.set nb334_fjxH2, 1296
.set nb334_fjyH2, 1312
.set nb334_fjzH2, 1328
.set nb334_fjxM, 1344
.set nb334_fjyM, 1360
.set nb334_fjzM, 1376
.set nb334_half, 1392
.set nb334_three, 1408
.set nb334_rsqOO, 1424
.set nb334_rsqH1H1, 1440
.set nb334_rsqH1H2, 1456
.set nb334_rsqH1M, 1472
.set nb334_rsqH2H1, 1488
.set nb334_rsqH2H2, 1504
.set nb334_rsqH2M, 1520
.set nb334_rsqMH1, 1536
.set nb334_rsqMH2, 1552
.set nb334_rsqMM, 1568
.set nb334_rinvOO, 1584
.set nb334_rinvH1H1, 1600
.set nb334_rinvH1H2, 1616
.set nb334_rinvH1M, 1632
.set nb334_rinvH2H1, 1648
.set nb334_rinvH2H2, 1664
.set nb334_rinvH2M, 1680
.set nb334_rinvMH1, 1696
.set nb334_rinvMH2, 1712
.set nb334_rinvMM, 1728
.set nb334_fscal, 1744
.set nb334_is3, 1760
.set nb334_ii3, 1764
.set nb334_innerjjnr, 1768
.set nb334_innerk, 1772
.set nb334_n, 1776
.set nb334_nn1, 1780
.set nb334_nri, 1784
.set nb334_nouter, 1788
.set nb334_ninner, 1792
.set nb334_salign, 1796
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $1800,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb334_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb334_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb334_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb334_nouter(%esp)
        movl %eax,nb334_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb334_half(%esp)
        movl %ebx,nb334_half+4(%esp)
        movsd nb334_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb334_half(%esp)
        movapd %xmm2,nb334_two(%esp)
        movapd %xmm3,nb334_three(%esp)
        movl nb334_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3

        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb334_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb334_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb334_charge(%ebp),%edx
        movsd 24(%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb334_p_facel(%ebp),%esi
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
        movapd %xmm3,nb334_qqMM(%esp)
        movapd %xmm4,nb334_qqMH(%esp)
        movapd %xmm5,nb334_qqHH(%esp)

        xorpd %xmm0,%xmm0
        movl  nb334_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb334_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb334_vdwparam(%ebp),%eax
        movlpd (%eax,%edx,8),%xmm0
        movhpd 8(%eax,%edx,8),%xmm0
        movhlps %xmm0,%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb334_c6(%esp)
        movapd %xmm1,nb334_c12(%esp)

_nb_kernel334_ia32_sse2.nb334_threadloop: 
        movl  nb334_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel334_ia32_sse2.nb334_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel334_ia32_sse2.nb334_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb334_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb334_n(%esp)
        movl %ebx,nb334_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel334_ia32_sse2.nb334_outerstart
        jmp _nb_kernel334_ia32_sse2.nb334_end

_nb_kernel334_ia32_sse2.nb334_outerstart: 
        ## ebx contains number of outer iterations
        addl nb334_nouter(%esp),%ebx
        movl %ebx,nb334_nouter(%esp)

_nb_kernel334_ia32_sse2.nb334_outer: 
        movl  nb334_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb334_is3(%esp)      ## store is3 

        movl  nb334_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb334_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb334_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb334_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3       ## ox
        addsd 8(%eax,%ebx,8),%xmm4      ## oy
        addsd 16(%eax,%ebx,8),%xmm5     ## oz   
        addsd 24(%eax,%ebx,8),%xmm6     ## h1x
        addsd 32(%eax,%ebx,8),%xmm7     ## h1y
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm7,%xmm7
        movapd %xmm3,nb334_ixO(%esp)
        movapd %xmm4,nb334_iyO(%esp)
        movapd %xmm5,nb334_izO(%esp)
        movapd %xmm6,nb334_ixH1(%esp)
        movapd %xmm7,nb334_iyH1(%esp)

        movsd %xmm2,%xmm6
        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 40(%eax,%ebx,8),%xmm6    ## h1z
        addsd 48(%eax,%ebx,8),%xmm0    ## h2x
        addsd 56(%eax,%ebx,8),%xmm1    ## h2y
        addsd 64(%eax,%ebx,8),%xmm2    ## h2z
        addsd 72(%eax,%ebx,8),%xmm3    ## mx
        addsd 80(%eax,%ebx,8),%xmm4    ## my
        addsd 88(%eax,%ebx,8),%xmm5    ## mz

        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm6,nb334_izH1(%esp)
        movapd %xmm0,nb334_ixH2(%esp)
        movapd %xmm1,nb334_iyH2(%esp)
        movapd %xmm2,nb334_izH2(%esp)
        movapd %xmm3,nb334_ixM(%esp)
        movapd %xmm4,nb334_iyM(%esp)
        movapd %xmm5,nb334_izM(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb334_vctot(%esp)
        movapd %xmm4,nb334_Vvdwtot(%esp)
        movapd %xmm4,nb334_fixO(%esp)
        movapd %xmm4,nb334_fiyO(%esp)
        movapd %xmm4,nb334_fizO(%esp)
        movapd %xmm4,nb334_fixH1(%esp)
        movapd %xmm4,nb334_fiyH1(%esp)
        movapd %xmm4,nb334_fizH1(%esp)
        movapd %xmm4,nb334_fixH2(%esp)
        movapd %xmm4,nb334_fiyH2(%esp)
        movapd %xmm4,nb334_fizH2(%esp)
        movapd %xmm4,nb334_fixM(%esp)
        movapd %xmm4,nb334_fiyM(%esp)
        movapd %xmm4,nb334_fizM(%esp)

        movl  nb334_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb334_pos(%ebp),%esi
        movl  nb334_faction(%ebp),%edi
        movl  nb334_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb334_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb334_ninner(%esp),%ecx
        movl  %ecx,nb334_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb334_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel334_ia32_sse2.nb334_unroll_loop
        jmp   _nb_kernel334_ia32_sse2.nb334_checksingle
_nb_kernel334_ia32_sse2.nb334_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb334_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb334_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb334_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%esi,%eax,8),%xmm0
        movlpd (%esi,%ebx,8),%xmm2
        movhpd 8(%esi,%eax,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm2
        movlpd 16(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%ebx,8),%xmm5
        movhpd 24(%esi,%eax,8),%xmm3
        movhpd 24(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## ox 
        unpckhpd %xmm2,%xmm1 ## oy
        unpcklpd %xmm5,%xmm3 ## ox 
        unpckhpd %xmm5,%xmm4 ## oy
        movapd  %xmm0,nb334_jxO(%esp)
        movapd  %xmm1,nb334_jyO(%esp)
        movapd  %xmm3,nb334_jzO(%esp)
        movapd  %xmm4,nb334_jxH1(%esp)

        ## load h1y, h1z, h2x, h2y 
        movlpd 32(%esi,%eax,8),%xmm0
        movlpd 32(%esi,%ebx,8),%xmm2
        movhpd 40(%esi,%eax,8),%xmm0
        movhpd 40(%esi,%ebx,8),%xmm2
        movlpd 48(%esi,%eax,8),%xmm3
        movlpd 48(%esi,%ebx,8),%xmm5
        movhpd 56(%esi,%eax,8),%xmm3
        movhpd 56(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h1y
        unpckhpd %xmm2,%xmm1 ## h1z
        unpcklpd %xmm5,%xmm3 ## h2x
        unpckhpd %xmm5,%xmm4 ## h2y
        movapd  %xmm0,nb334_jyH1(%esp)
        movapd  %xmm1,nb334_jzH1(%esp)
        movapd  %xmm3,nb334_jxH2(%esp)
        movapd  %xmm4,nb334_jyH2(%esp)

        ## load h2z, mx, my, mz
        movlpd 64(%esi,%eax,8),%xmm0
        movlpd 64(%esi,%ebx,8),%xmm2
        movhpd 72(%esi,%eax,8),%xmm0
        movhpd 72(%esi,%ebx,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 80(%esi,%ebx,8),%xmm5
        movhpd 88(%esi,%eax,8),%xmm3
        movhpd 88(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h2z
        unpckhpd %xmm2,%xmm1 ## mx
        unpcklpd %xmm5,%xmm3 ## my
        unpckhpd %xmm5,%xmm4 ## mz
        movapd  %xmm0,nb334_jzH2(%esp)
        movapd  %xmm1,nb334_jxM(%esp)
        movapd  %xmm3,nb334_jyM(%esp)
        movapd  %xmm4,nb334_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb334_ixO(%esp),%xmm0
        movapd nb334_iyO(%esp),%xmm1
        movapd nb334_izO(%esp),%xmm2
        movapd nb334_ixH1(%esp),%xmm3
        movapd nb334_iyH1(%esp),%xmm4
        movapd nb334_izH1(%esp),%xmm5
        subpd  nb334_jxO(%esp),%xmm0
        subpd  nb334_jyO(%esp),%xmm1
        subpd  nb334_jzO(%esp),%xmm2
        subpd  nb334_jxH1(%esp),%xmm3
        subpd  nb334_jyH1(%esp),%xmm4
        subpd  nb334_jzH1(%esp),%xmm5
        movapd %xmm0,nb334_dxOO(%esp)
        movapd %xmm1,nb334_dyOO(%esp)
        movapd %xmm2,nb334_dzOO(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb334_dxH1H1(%esp)
        movapd %xmm4,nb334_dyH1H1(%esp)
        movapd %xmm5,nb334_dzH1H1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb334_rsqOO(%esp)
        movapd %xmm3,nb334_rsqH1H1(%esp)

        movapd nb334_ixH1(%esp),%xmm0
        movapd nb334_iyH1(%esp),%xmm1
        movapd nb334_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb334_jxH2(%esp),%xmm0
        subpd  nb334_jyH2(%esp),%xmm1
        subpd  nb334_jzH2(%esp),%xmm2
        subpd  nb334_jxM(%esp),%xmm3
        subpd  nb334_jyM(%esp),%xmm4
        subpd  nb334_jzM(%esp),%xmm5
        movapd %xmm0,nb334_dxH1H2(%esp)
        movapd %xmm1,nb334_dyH1H2(%esp)
        movapd %xmm2,nb334_dzH1H2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb334_dxH1M(%esp)
        movapd %xmm4,nb334_dyH1M(%esp)
        movapd %xmm5,nb334_dzH1M(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb334_rsqH1H2(%esp)
        movapd %xmm3,nb334_rsqH1M(%esp)

        movapd nb334_ixH2(%esp),%xmm0
        movapd nb334_iyH2(%esp),%xmm1
        movapd nb334_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb334_jxH1(%esp),%xmm0
        subpd  nb334_jyH1(%esp),%xmm1
        subpd  nb334_jzH1(%esp),%xmm2
        subpd  nb334_jxH2(%esp),%xmm3
        subpd  nb334_jyH2(%esp),%xmm4
        subpd  nb334_jzH2(%esp),%xmm5
        movapd %xmm0,nb334_dxH2H1(%esp)
        movapd %xmm1,nb334_dyH2H1(%esp)
        movapd %xmm2,nb334_dzH2H1(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb334_dxH2H2(%esp)
        movapd %xmm4,nb334_dyH2H2(%esp)
        movapd %xmm5,nb334_dzH2H2(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb334_rsqH2H1(%esp)
        movapd %xmm3,nb334_rsqH2H2(%esp)

        movapd nb334_ixH2(%esp),%xmm0
        movapd nb334_iyH2(%esp),%xmm1
        movapd nb334_izH2(%esp),%xmm2
        movapd nb334_ixM(%esp),%xmm3
        movapd nb334_iyM(%esp),%xmm4
        movapd nb334_izM(%esp),%xmm5
        subpd  nb334_jxM(%esp),%xmm0
        subpd  nb334_jyM(%esp),%xmm1
        subpd  nb334_jzM(%esp),%xmm2
        subpd  nb334_jxH1(%esp),%xmm3
        subpd  nb334_jyH1(%esp),%xmm4
        subpd  nb334_jzH1(%esp),%xmm5
        movapd %xmm0,nb334_dxH2M(%esp)
        movapd %xmm1,nb334_dyH2M(%esp)
        movapd %xmm2,nb334_dzH2M(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb334_dxMH1(%esp)
        movapd %xmm4,nb334_dyMH1(%esp)
        movapd %xmm5,nb334_dzMH1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb334_rsqH2M(%esp)
        movapd %xmm4,nb334_rsqMH1(%esp)

        movapd nb334_ixM(%esp),%xmm0
        movapd nb334_iyM(%esp),%xmm1
        movapd nb334_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb334_jxH2(%esp),%xmm0
        subpd  nb334_jyH2(%esp),%xmm1
        subpd  nb334_jzH2(%esp),%xmm2
        subpd  nb334_jxM(%esp),%xmm3
        subpd  nb334_jyM(%esp),%xmm4
        subpd  nb334_jzM(%esp),%xmm5
        movapd %xmm0,nb334_dxMH2(%esp)
        movapd %xmm1,nb334_dyMH2(%esp)
        movapd %xmm2,nb334_dzMH2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb334_dxMM(%esp)
        movapd %xmm4,nb334_dyMM(%esp)
        movapd %xmm5,nb334_dzMM(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb334_rsqMH2(%esp)
        movapd %xmm4,nb334_rsqMM(%esp)

        ## Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1  ## luA
        cvtps2pd %xmm5,%xmm5  ## luB

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb334_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334_half(%esp),%xmm3   ## iter1 
        mulpd   nb334_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334_half(%esp),%xmm1   ## rinv 
        mulpd   nb334_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb334_rinvMH2(%esp)
        movapd %xmm5,nb334_rinvMM(%esp)

        movapd nb334_rsqOO(%esp),%xmm0
        movapd nb334_rsqH1H1(%esp),%xmm4
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
        movapd  nb334_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb334_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334_half(%esp),%xmm1   ## rinv 
        mulpd   nb334_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb334_rinvOO(%esp)
        movapd %xmm5,nb334_rinvH1H1(%esp)

        movapd nb334_rsqH1H2(%esp),%xmm0
        movapd nb334_rsqH1M(%esp),%xmm4
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
        movapd  nb334_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334_half(%esp),%xmm3   ## iter1 
        mulpd   nb334_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334_half(%esp),%xmm1   ## rinv 
        mulpd   nb334_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb334_rinvH1H2(%esp)
        movapd %xmm5,nb334_rinvH1M(%esp)

        movapd nb334_rsqH2H1(%esp),%xmm0
        movapd nb334_rsqH2H2(%esp),%xmm4
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
        movapd  nb334_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334_half(%esp),%xmm3   ## iter1a 
        mulpd   nb334_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334_half(%esp),%xmm1   ## rinv 
        mulpd   nb334_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb334_rinvH2H1(%esp)
        movapd %xmm5,nb334_rinvH2H2(%esp)

        movapd nb334_rsqMH1(%esp),%xmm0
        movapd nb334_rsqH2M(%esp),%xmm4
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
        movapd  nb334_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334_half(%esp),%xmm3   ## iter1a 
        mulpd   nb334_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334_half(%esp),%xmm1   ## rinv 
        mulpd   nb334_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb334_rinvMH1(%esp)
        movapd %xmm5,nb334_rinvH2M(%esp)

        ## start with OO interaction 
        movapd nb334_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movd %eax,%mm0
        movd %ebx,%mm1
        movl nb334_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

        ## Dispersion 
        movlpd 32(%esi,%eax,8),%xmm4    ## Y1
        movlpd 32(%esi,%ebx,8),%xmm3    ## Y2
        movhpd 40(%esi,%eax,8),%xmm4    ## Y1 F1        
        movhpd 40(%esi,%ebx,8),%xmm3    ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%esi,%eax,8),%xmm6    ## G1
        movlpd 48(%esi,%ebx,8),%xmm3    ## G2
        movhpd 56(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 56(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## Dispersion table ready, in xmm4-xmm7                 

                mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb334_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb334_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm7       ## fijD 
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack. Update Vvdwtot directly 
        addpd  nb334_Vvdwtot(%esp),%xmm5
        movapd %xmm7,nb334_fscal(%esp)
        movapd %xmm5,nb334_Vvdwtot(%esp)

        ## Repulsion 
        movlpd 64(%esi,%eax,8),%xmm4    ## Y1
        movlpd 64(%esi,%ebx,8),%xmm3    ## Y2
        movhpd 72(%esi,%eax,8),%xmm4    ## Y1 F1        
        movhpd 72(%esi,%ebx,8),%xmm3    ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 80(%esi,%eax,8),%xmm6    ## G1
        movlpd 80(%esi,%ebx,8),%xmm3    ## G2
        movhpd 88(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 88(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## Repulsion table ready, in xmm4-xmm7                  
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb334_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb334_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm7 ## fijR 
        mulpd  %xmm4,%xmm5 ## Vvdw12 
        addpd  nb334_fscal(%esp),%xmm7

        addpd  nb334_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb334_Vvdwtot(%esp)
        xorpd  %xmm4,%xmm4

        mulpd nb334_tsc(%esp),%xmm7
        mulpd nb334_rinvOO(%esp),%xmm7
        subpd %xmm7,%xmm4

        movapd %xmm4,%xmm0
        movapd %xmm4,%xmm1
        movapd %xmm4,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb334_dxOO(%esp),%xmm0
        mulpd nb334_dyOO(%esp),%xmm1
        mulpd nb334_dzOO(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb334_fixO(%esp),%xmm0
        addpd nb334_fiyO(%esp),%xmm1
        addpd nb334_fizO(%esp),%xmm2
        movapd %xmm3,nb334_fjxO(%esp)
        movapd %xmm4,nb334_fjyO(%esp)
        movapd %xmm5,nb334_fjzO(%esp)
        movapd %xmm0,nb334_fixO(%esp)
        movapd %xmm1,nb334_fiyO(%esp)
        movapd %xmm2,nb334_fizO(%esp)

        ## H1-H1 interaction 
        movapd nb334_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

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
        mulpd  nb334_two(%esp),%xmm7    ## two*Heps2 
        movapd nb334_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addpd  nb334_vctot(%esp),%xmm5
        movapd %xmm5,nb334_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb334_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb334_dxH1H1(%esp),%xmm0
        mulpd nb334_dyH1H1(%esp),%xmm1
        mulpd nb334_dzH1H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb334_fixH1(%esp),%xmm0
        addpd nb334_fiyH1(%esp),%xmm1
        addpd nb334_fizH1(%esp),%xmm2
        movapd %xmm3,nb334_fjxH1(%esp)
        movapd %xmm4,nb334_fjyH1(%esp)
        movapd %xmm5,nb334_fjzH1(%esp)
        movapd %xmm0,nb334_fixH1(%esp)
        movapd %xmm1,nb334_fiyH1(%esp)
        movapd %xmm2,nb334_fizH1(%esp)

        ## H1-H2 interaction  
        movapd nb334_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

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
        mulpd  nb334_two(%esp),%xmm7    ## two*Heps2 
        movapd nb334_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addpd  nb334_vctot(%esp),%xmm5
        movapd %xmm5,nb334_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb334_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb334_dxH1H2(%esp),%xmm0
        mulpd nb334_dyH1H2(%esp),%xmm1
        mulpd nb334_dzH1H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb334_fixH1(%esp),%xmm0
        addpd nb334_fiyH1(%esp),%xmm1
        addpd nb334_fizH1(%esp),%xmm2
        movapd %xmm3,nb334_fjxH2(%esp)
        movapd %xmm4,nb334_fjyH2(%esp)
        movapd %xmm5,nb334_fjzH2(%esp)
        movapd %xmm0,nb334_fixH1(%esp)
        movapd %xmm1,nb334_fiyH1(%esp)
        movapd %xmm2,nb334_fizH1(%esp)

        ## H1-M interaction 
        movapd nb334_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334_rsqH1M(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

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
        mulpd  nb334_two(%esp),%xmm7    ## two*Heps2 
        movapd nb334_qqMH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addpd  nb334_vctot(%esp),%xmm5
        movapd %xmm5,nb334_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb334_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb334_dxH1M(%esp),%xmm0
        mulpd nb334_dyH1M(%esp),%xmm1
        mulpd nb334_dzH1M(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb334_fixH1(%esp),%xmm0
        addpd nb334_fiyH1(%esp),%xmm1
        addpd nb334_fizH1(%esp),%xmm2
        movapd %xmm3,nb334_fjxM(%esp)
        movapd %xmm4,nb334_fjyM(%esp)
        movapd %xmm5,nb334_fjzM(%esp)
        movapd %xmm0,nb334_fixH1(%esp)
        movapd %xmm1,nb334_fiyH1(%esp)
        movapd %xmm2,nb334_fizH1(%esp)

        ## H2-H1 interaction 
        movapd nb334_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

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
        mulpd  nb334_two(%esp),%xmm7    ## two*Heps2 
        movapd nb334_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addpd  nb334_vctot(%esp),%xmm5
        movapd %xmm5,nb334_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb334_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb334_fjxH1(%esp),%xmm3
        movapd nb334_fjyH1(%esp),%xmm4
        movapd nb334_fjzH1(%esp),%xmm5
        mulpd nb334_dxH2H1(%esp),%xmm0
        mulpd nb334_dyH2H1(%esp),%xmm1
        mulpd nb334_dzH2H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb334_fixH2(%esp),%xmm0
        addpd nb334_fiyH2(%esp),%xmm1
        addpd nb334_fizH2(%esp),%xmm2
        movapd %xmm3,nb334_fjxH1(%esp)
        movapd %xmm4,nb334_fjyH1(%esp)
        movapd %xmm5,nb334_fjzH1(%esp)
        movapd %xmm0,nb334_fixH2(%esp)
        movapd %xmm1,nb334_fiyH2(%esp)
        movapd %xmm2,nb334_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb334_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

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
        mulpd  nb334_two(%esp),%xmm7    ## two*Heps2 
        movapd nb334_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addpd  nb334_vctot(%esp),%xmm5
        movapd %xmm5,nb334_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb334_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb334_fjxH2(%esp),%xmm3
        movapd nb334_fjyH2(%esp),%xmm4
        movapd nb334_fjzH2(%esp),%xmm5
        mulpd nb334_dxH2H2(%esp),%xmm0
        mulpd nb334_dyH2H2(%esp),%xmm1
        mulpd nb334_dzH2H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb334_fixH2(%esp),%xmm0
        addpd nb334_fiyH2(%esp),%xmm1
        addpd nb334_fizH2(%esp),%xmm2
        movapd %xmm3,nb334_fjxH2(%esp)
        movapd %xmm4,nb334_fjyH2(%esp)
        movapd %xmm5,nb334_fjzH2(%esp)
        movapd %xmm0,nb334_fixH2(%esp)
        movapd %xmm1,nb334_fiyH2(%esp)
        movapd %xmm2,nb334_fizH2(%esp)

        ## H2-M interaction 
        movapd nb334_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334_rsqH2M(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

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
        mulpd  nb334_two(%esp),%xmm7    ## two*Heps2 
        movapd nb334_qqMH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addpd  nb334_vctot(%esp),%xmm5
        movapd %xmm5,nb334_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb334_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb334_fjxM(%esp),%xmm3
        movapd nb334_fjyM(%esp),%xmm4
        movapd nb334_fjzM(%esp),%xmm5
        mulpd nb334_dxH2M(%esp),%xmm0
        mulpd nb334_dyH2M(%esp),%xmm1
        mulpd nb334_dzH2M(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb334_fixH2(%esp),%xmm0
        addpd nb334_fiyH2(%esp),%xmm1
        addpd nb334_fizH2(%esp),%xmm2
        movapd %xmm3,nb334_fjxM(%esp)
        movapd %xmm4,nb334_fjyM(%esp)
        movapd %xmm5,nb334_fjzM(%esp)
        movapd %xmm0,nb334_fixH2(%esp)
        movapd %xmm1,nb334_fiyH2(%esp)
        movapd %xmm2,nb334_fizH2(%esp)

        ## M-H1 interaction 
        movapd nb334_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334_rsqMH1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

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
        mulpd  nb334_two(%esp),%xmm7    ## two*Heps2 
        movapd nb334_qqMH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addpd  nb334_vctot(%esp),%xmm5
        movapd %xmm5,nb334_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb334_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb334_fjxH1(%esp),%xmm3
        movapd nb334_fjyH1(%esp),%xmm4
        movapd nb334_fjzH1(%esp),%xmm5
        mulpd nb334_dxMH1(%esp),%xmm0
        mulpd nb334_dyMH1(%esp),%xmm1
        mulpd nb334_dzMH1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb334_fixM(%esp),%xmm0
        addpd nb334_fiyM(%esp),%xmm1
        addpd nb334_fizM(%esp),%xmm2
        movapd %xmm3,nb334_fjxH1(%esp)
        movapd %xmm4,nb334_fjyH1(%esp)
        movapd %xmm5,nb334_fjzH1(%esp)
        movapd %xmm0,nb334_fixM(%esp)
        movapd %xmm1,nb334_fiyM(%esp)
        movapd %xmm2,nb334_fizM(%esp)

        ## M-H2 interaction 
        movapd nb334_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334_rsqMH2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

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
        mulpd  nb334_two(%esp),%xmm7    ## two*Heps2 
        movapd nb334_qqMH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addpd  nb334_vctot(%esp),%xmm5
        movapd %xmm5,nb334_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb334_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb334_fjxH2(%esp),%xmm3
        movapd nb334_fjyH2(%esp),%xmm4
        movapd nb334_fjzH2(%esp),%xmm5
        mulpd nb334_dxMH2(%esp),%xmm0
        mulpd nb334_dyMH2(%esp),%xmm1
        mulpd nb334_dzMH2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb334_fixM(%esp),%xmm0
        addpd nb334_fiyM(%esp),%xmm1
        addpd nb334_fizM(%esp),%xmm2
        movapd %xmm3,nb334_fjxH2(%esp)
        movapd %xmm4,nb334_fjyH2(%esp)
        movapd %xmm5,nb334_fjzH2(%esp)
        movapd %xmm0,nb334_fixM(%esp)
        movapd %xmm1,nb334_fiyM(%esp)
        movapd %xmm2,nb334_fizM(%esp)

        ## M-M interaction 
        movapd nb334_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334_rsqMM(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

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
        mulpd  nb334_two(%esp),%xmm7    ## two*Heps2 
        movapd nb334_qqMM(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addpd  nb334_vctot(%esp),%xmm5
        movapd %xmm5,nb334_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb334_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb334_fjxM(%esp),%xmm3
        movapd nb334_fjyM(%esp),%xmm4
        movapd nb334_fjzM(%esp),%xmm5
        mulpd nb334_dxMM(%esp),%xmm0
        mulpd nb334_dyMM(%esp),%xmm1
        mulpd nb334_dzMM(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb334_fixM(%esp),%xmm0
        addpd nb334_fiyM(%esp),%xmm1
        addpd nb334_fizM(%esp),%xmm2
        movapd %xmm3,nb334_fjxM(%esp)
        movapd %xmm4,nb334_fjyM(%esp)
        movapd %xmm5,nb334_fjzM(%esp)
        movapd %xmm0,nb334_fixM(%esp)
        movapd %xmm1,nb334_fiyM(%esp)
        movapd %xmm2,nb334_fizM(%esp)

        movl nb334_faction(%ebp),%edi

        movd %mm0,%eax
        movd %mm1,%ebx

        ## Did all interactions - now update j forces 
        ## Step1 - transpose fjxO, fjyO and fjzO, fjxH1
        movapd nb334_fjxO(%esp),%xmm0
        movapd nb334_fjyO(%esp),%xmm1
        movapd nb334_fjzO(%esp),%xmm2
        movapd nb334_fjxH1(%esp),%xmm3
        movapd %xmm0,%xmm4
        movapd %xmm2,%xmm5
        unpcklpd %xmm1,%xmm0  ## fjOxA fjOyA
        unpckhpd %xmm1,%xmm4  ## fjOxB fjOyB
        unpcklpd %xmm3,%xmm2  ## fjOzA fjH1xA
        unpckhpd %xmm3,%xmm5  ## fjOzB fjH1xB
        movlpd (%edi,%eax,8),%xmm1
        movlpd (%edi,%ebx,8),%xmm3
        movhpd 8(%edi,%eax,8),%xmm1
        movhpd 8(%edi,%ebx,8),%xmm3
        movlpd 16(%edi,%eax,8),%xmm6
        movlpd 16(%edi,%ebx,8),%xmm7
        movhpd 24(%edi,%eax,8),%xmm6
        movhpd 24(%edi,%ebx,8),%xmm7
        addpd  %xmm0,%xmm1
        addpd  %xmm4,%xmm3
        addpd  %xmm2,%xmm6
        addpd  %xmm5,%xmm7
        movlpd %xmm1,(%edi,%eax,8)
        movlpd %xmm3,(%edi,%ebx,8)
        movhpd %xmm1,8(%edi,%eax,8)
        movhpd %xmm3,8(%edi,%ebx,8)
        movlpd %xmm6,16(%edi,%eax,8)
        movlpd %xmm7,16(%edi,%ebx,8)
        movhpd %xmm6,24(%edi,%eax,8)
        movhpd %xmm7,24(%edi,%ebx,8)

        ## Step2 - transpose fjyH1, fjzH1 and fjxH2, fjyH2
        movapd nb334_fjyH1(%esp),%xmm0
        movapd nb334_fjzH1(%esp),%xmm1
        movapd nb334_fjxH2(%esp),%xmm2
        movapd nb334_fjyH2(%esp),%xmm3
        movapd %xmm0,%xmm4
        movapd %xmm2,%xmm5
        unpcklpd %xmm1,%xmm0  ## fjOxA fjOyA
        unpckhpd %xmm1,%xmm4  ## fjOxB fjOyB
        unpcklpd %xmm3,%xmm2  ## fjOzA fjH1xA
        unpckhpd %xmm3,%xmm5  ## fjOzB fjH1xB
        movlpd 32(%edi,%eax,8),%xmm1
        movlpd 32(%edi,%ebx,8),%xmm3
        movhpd 40(%edi,%eax,8),%xmm1
        movhpd 40(%edi,%ebx,8),%xmm3
        movlpd 48(%edi,%eax,8),%xmm6
        movlpd 48(%edi,%ebx,8),%xmm7
        movhpd 56(%edi,%eax,8),%xmm6
        movhpd 56(%edi,%ebx,8),%xmm7
        addpd  %xmm0,%xmm1
        addpd  %xmm4,%xmm3
        addpd  %xmm2,%xmm6
        addpd  %xmm5,%xmm7
        movlpd %xmm1,32(%edi,%eax,8)
        movlpd %xmm3,32(%edi,%ebx,8)
        movhpd %xmm1,40(%edi,%eax,8)
        movhpd %xmm3,40(%edi,%ebx,8)
        movlpd %xmm6,48(%edi,%eax,8)
        movlpd %xmm7,48(%edi,%ebx,8)
        movhpd %xmm6,56(%edi,%eax,8)
        movhpd %xmm7,56(%edi,%ebx,8)

        ## Step3 - transpose fjzH2, fjxM and fjyM, fjzM
        movapd nb334_fjzH2(%esp),%xmm0
        movapd nb334_fjxM(%esp),%xmm1
        movapd nb334_fjyM(%esp),%xmm2
        movapd nb334_fjzM(%esp),%xmm3
        movapd %xmm0,%xmm4
        movapd %xmm2,%xmm5
        unpcklpd %xmm1,%xmm0  ## fjOxA fjOyA
        unpckhpd %xmm1,%xmm4  ## fjOxB fjOyB
        unpcklpd %xmm3,%xmm2  ## fjOzA fjH1xA
        unpckhpd %xmm3,%xmm5  ## fjOzB fjH1xB
        movlpd 64(%edi,%eax,8),%xmm1
        movlpd 64(%edi,%ebx,8),%xmm3
        movhpd 72(%edi,%eax,8),%xmm1
        movhpd 72(%edi,%ebx,8),%xmm3
        movlpd 80(%edi,%eax,8),%xmm6
        movlpd 80(%edi,%ebx,8),%xmm7
        movhpd 88(%edi,%eax,8),%xmm6
        movhpd 88(%edi,%ebx,8),%xmm7
        addpd  %xmm0,%xmm1
        addpd  %xmm4,%xmm3
        addpd  %xmm2,%xmm6
        addpd  %xmm5,%xmm7
        movlpd %xmm1,64(%edi,%eax,8)
        movlpd %xmm3,64(%edi,%ebx,8)
        movhpd %xmm1,72(%edi,%eax,8)
        movhpd %xmm3,72(%edi,%ebx,8)
        movlpd %xmm6,80(%edi,%eax,8)
        movlpd %xmm7,80(%edi,%ebx,8)
        movhpd %xmm6,88(%edi,%eax,8)
        movhpd %xmm7,88(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl $2,nb334_innerk(%esp)
        jl    _nb_kernel334_ia32_sse2.nb334_checksingle
        jmp   _nb_kernel334_ia32_sse2.nb334_unroll_loop
_nb_kernel334_ia32_sse2.nb334_checksingle: 
        movl  nb334_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel334_ia32_sse2.nb334_dosingle
        jmp   _nb_kernel334_ia32_sse2.nb334_updateouterdata
_nb_kernel334_ia32_sse2.nb334_dosingle: 
        movl  nb334_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb334_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%esi,%eax,8),%xmm0
        movhpd 8(%esi,%eax,8),%xmm0
        movlpd 16(%esi,%eax,8),%xmm1
        movhpd 24(%esi,%eax,8),%xmm1
        movlpd 32(%esi,%eax,8),%xmm2
        movhpd 40(%esi,%eax,8),%xmm2
        movlpd 48(%esi,%eax,8),%xmm3
        movhpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 72(%esi,%eax,8),%xmm4
        movlpd 80(%esi,%eax,8),%xmm5
        movhpd 88(%esi,%eax,8),%xmm5
        movsd  %xmm0,nb334_jxO(%esp)
        movsd  %xmm1,nb334_jzO(%esp)
        movsd  %xmm2,nb334_jyH1(%esp)
        movsd  %xmm3,nb334_jxH2(%esp)
        movsd  %xmm4,nb334_jzH2(%esp)
        movsd  %xmm5,nb334_jyM(%esp)
        unpckhpd %xmm0,%xmm0
        unpckhpd %xmm1,%xmm1
        unpckhpd %xmm2,%xmm2
        unpckhpd %xmm3,%xmm3
        unpckhpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5
        movsd  %xmm0,nb334_jyO(%esp)
        movsd  %xmm1,nb334_jxH1(%esp)
        movsd  %xmm2,nb334_jzH1(%esp)
        movsd  %xmm3,nb334_jyH2(%esp)
        movsd  %xmm4,nb334_jxM(%esp)
        movsd  %xmm5,nb334_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb334_ixO(%esp),%xmm0
        movapd nb334_iyO(%esp),%xmm1
        movapd nb334_izO(%esp),%xmm2
        movapd nb334_ixH1(%esp),%xmm3
        movapd nb334_iyH1(%esp),%xmm4
        movapd nb334_izH1(%esp),%xmm5
        subsd  nb334_jxO(%esp),%xmm0
        subsd  nb334_jyO(%esp),%xmm1
        subsd  nb334_jzO(%esp),%xmm2
        subsd  nb334_jxH1(%esp),%xmm3
        subsd  nb334_jyH1(%esp),%xmm4
        subsd  nb334_jzH1(%esp),%xmm5
        movapd %xmm0,nb334_dxOO(%esp)
        movapd %xmm1,nb334_dyOO(%esp)
        movapd %xmm2,nb334_dzOO(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb334_dxH1H1(%esp)
        movapd %xmm4,nb334_dyH1H1(%esp)
        movapd %xmm5,nb334_dzH1H1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb334_rsqOO(%esp)
        movapd %xmm3,nb334_rsqH1H1(%esp)

        movapd nb334_ixH1(%esp),%xmm0
        movapd nb334_iyH1(%esp),%xmm1
        movapd nb334_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb334_jxH2(%esp),%xmm0
        subsd  nb334_jyH2(%esp),%xmm1
        subsd  nb334_jzH2(%esp),%xmm2
        subsd  nb334_jxM(%esp),%xmm3
        subsd  nb334_jyM(%esp),%xmm4
        subsd  nb334_jzM(%esp),%xmm5
        movapd %xmm0,nb334_dxH1H2(%esp)
        movapd %xmm1,nb334_dyH1H2(%esp)
        movapd %xmm2,nb334_dzH1H2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb334_dxH1M(%esp)
        movapd %xmm4,nb334_dyH1M(%esp)
        movapd %xmm5,nb334_dzH1M(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb334_rsqH1H2(%esp)
        movapd %xmm3,nb334_rsqH1M(%esp)

        movapd nb334_ixH2(%esp),%xmm0
        movapd nb334_iyH2(%esp),%xmm1
        movapd nb334_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb334_jxH1(%esp),%xmm0
        subsd  nb334_jyH1(%esp),%xmm1
        subsd  nb334_jzH1(%esp),%xmm2
        subsd  nb334_jxH2(%esp),%xmm3
        subsd  nb334_jyH2(%esp),%xmm4
        subsd  nb334_jzH2(%esp),%xmm5
        movapd %xmm0,nb334_dxH2H1(%esp)
        movapd %xmm1,nb334_dyH2H1(%esp)
        movapd %xmm2,nb334_dzH2H1(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb334_dxH2H2(%esp)
        movapd %xmm4,nb334_dyH2H2(%esp)
        movapd %xmm5,nb334_dzH2H2(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb334_rsqH2H1(%esp)
        movapd %xmm3,nb334_rsqH2H2(%esp)

        movapd nb334_ixH2(%esp),%xmm0
        movapd nb334_iyH2(%esp),%xmm1
        movapd nb334_izH2(%esp),%xmm2
        movapd nb334_ixM(%esp),%xmm3
        movapd nb334_iyM(%esp),%xmm4
        movapd nb334_izM(%esp),%xmm5
        subsd  nb334_jxM(%esp),%xmm0
        subsd  nb334_jyM(%esp),%xmm1
        subsd  nb334_jzM(%esp),%xmm2
        subsd  nb334_jxH1(%esp),%xmm3
        subsd  nb334_jyH1(%esp),%xmm4
        subsd  nb334_jzH1(%esp),%xmm5
        movapd %xmm0,nb334_dxH2M(%esp)
        movapd %xmm1,nb334_dyH2M(%esp)
        movapd %xmm2,nb334_dzH2M(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb334_dxMH1(%esp)
        movapd %xmm4,nb334_dyMH1(%esp)
        movapd %xmm5,nb334_dzMH1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb334_rsqH2M(%esp)
        movapd %xmm4,nb334_rsqMH1(%esp)

        movapd nb334_ixM(%esp),%xmm0
        movapd nb334_iyM(%esp),%xmm1
        movapd nb334_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb334_jxH2(%esp),%xmm0
        subsd  nb334_jyH2(%esp),%xmm1
        subsd  nb334_jzH2(%esp),%xmm2
        subsd  nb334_jxM(%esp),%xmm3
        subsd  nb334_jyM(%esp),%xmm4
        subsd  nb334_jzM(%esp),%xmm5
        movapd %xmm0,nb334_dxMH2(%esp)
        movapd %xmm1,nb334_dyMH2(%esp)
        movapd %xmm2,nb334_dzMH2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb334_dxMM(%esp)
        movapd %xmm4,nb334_dyMM(%esp)
        movapd %xmm5,nb334_dzMM(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb334_rsqMH2(%esp)
        movapd %xmm4,nb334_rsqMM(%esp)

        ## Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1  ## luA
        cvtss2sd %xmm5,%xmm5  ## luB

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb334_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334_half(%esp),%xmm3   ## iter1 
        mulsd   nb334_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334_half(%esp),%xmm1   ## rinv 
        mulsd   nb334_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb334_rinvMH2(%esp)
        movapd %xmm5,nb334_rinvMM(%esp)

        movapd nb334_rsqOO(%esp),%xmm0
        movapd nb334_rsqH1H1(%esp),%xmm4
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
        movapd  nb334_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb334_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334_half(%esp),%xmm1   ## rinv 
        mulsd   nb334_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb334_rinvOO(%esp)
        movapd %xmm5,nb334_rinvH1H1(%esp)

        movapd nb334_rsqH1H2(%esp),%xmm0
        movapd nb334_rsqH1M(%esp),%xmm4
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
        movapd  nb334_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334_half(%esp),%xmm3   ## iter1 
        mulsd   nb334_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334_half(%esp),%xmm1   ## rinv 
        mulsd   nb334_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb334_rinvH1H2(%esp)
        movapd %xmm5,nb334_rinvH1M(%esp)

        movapd nb334_rsqH2H1(%esp),%xmm0
        movapd nb334_rsqH2H2(%esp),%xmm4
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
        movapd  nb334_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334_half(%esp),%xmm3   ## iter1a 
        mulsd   nb334_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334_half(%esp),%xmm1   ## rinv 
        mulsd   nb334_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb334_rinvH2H1(%esp)
        movapd %xmm5,nb334_rinvH2H2(%esp)

        movapd nb334_rsqMH1(%esp),%xmm0
        movapd nb334_rsqH2M(%esp),%xmm4
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
        movapd  nb334_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334_half(%esp),%xmm3   ## iter1a 
        mulsd   nb334_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334_half(%esp),%xmm1   ## rinv 
        mulsd   nb334_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb334_rinvMH1(%esp)
        movapd %xmm5,nb334_rinvH2M(%esp)

        ## start with OO interaction 
        movapd nb334_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334_tsc(%esp),%xmm1

        movd %eax,%mm0


        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 


        shll $2,%eax            ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        ## Dispersion 
        movsd 32(%esi,%eax,8),%xmm4     ## Y1   
        movsd 40(%esi,%eax,8),%xmm5     ## F1   
        movsd 48(%esi,%eax,8),%xmm6     ## G1   
        movsd 56(%esi,%eax,8),%xmm7     ## H1   

        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb334_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb334_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm7       ## fijD 
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack. Update Vvdwtot directly 
        addsd  nb334_Vvdwtot(%esp),%xmm5
        movsd %xmm7,nb334_fscal(%esp)
        movsd %xmm5,nb334_Vvdwtot(%esp)

        ## Repulsion 
        movsd 64(%esi,%eax,8),%xmm4     ## Y1   
        movsd 72(%esi,%eax,8),%xmm5     ## F1   
        movsd 80(%esi,%eax,8),%xmm6     ## G1
        movsd 88(%esi,%eax,8),%xmm7     ## H1   
        ## Repulsion table ready, in xmm4-xmm7                  
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb334_two(%esp),%xmm7    ## two*Heps2 
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb334_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm7 ## fijR 
        mulsd  %xmm4,%xmm5 ## Vvdw12 
        addsd  nb334_fscal(%esp),%xmm7

        addsd  nb334_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb334_Vvdwtot(%esp)
        xorpd  %xmm4,%xmm4

        mulsd nb334_tsc(%esp),%xmm7
        mulsd nb334_rinvOO(%esp),%xmm7
        subsd %xmm7,%xmm4

        movapd %xmm4,%xmm0
        movapd %xmm4,%xmm1
        movapd %xmm4,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb334_dxOO(%esp),%xmm0
        mulsd nb334_dyOO(%esp),%xmm1
        mulsd nb334_dzOO(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb334_fixO(%esp),%xmm0
        addsd nb334_fiyO(%esp),%xmm1
        addsd nb334_fizO(%esp),%xmm2
        movsd %xmm3,nb334_fjxO(%esp)
        movsd %xmm4,nb334_fjyO(%esp)
        movsd %xmm5,nb334_fjzO(%esp)
        movsd %xmm0,nb334_fixO(%esp)
        movsd %xmm1,nb334_fiyO(%esp)
        movsd %xmm2,nb334_fizO(%esp)

        ## H1-H1 interaction 
        movapd nb334_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb334_two(%esp),%xmm7    ## two*Heps2 
        movapd nb334_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addsd  nb334_vctot(%esp),%xmm5
        movsd %xmm5,nb334_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb334_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb334_dxH1H1(%esp),%xmm0
        mulsd nb334_dyH1H1(%esp),%xmm1
        mulsd nb334_dzH1H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb334_fixH1(%esp),%xmm0
        addsd nb334_fiyH1(%esp),%xmm1
        addsd nb334_fizH1(%esp),%xmm2
        movsd %xmm3,nb334_fjxH1(%esp)
        movsd %xmm4,nb334_fjyH1(%esp)
        movsd %xmm5,nb334_fjzH1(%esp)
        movsd %xmm0,nb334_fixH1(%esp)
        movsd %xmm1,nb334_fiyH1(%esp)
        movsd %xmm2,nb334_fizH1(%esp)

        ## H1-H2 interaction  
        movapd nb334_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb334_two(%esp),%xmm7    ## two*Heps2 
        movapd nb334_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addsd  nb334_vctot(%esp),%xmm5
        movsd %xmm5,nb334_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb334_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb334_dxH1H2(%esp),%xmm0
        mulsd nb334_dyH1H2(%esp),%xmm1
        mulsd nb334_dzH1H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb334_fixH1(%esp),%xmm0
        addsd nb334_fiyH1(%esp),%xmm1
        addsd nb334_fizH1(%esp),%xmm2
        movsd %xmm3,nb334_fjxH2(%esp)
        movsd %xmm4,nb334_fjyH2(%esp)
        movsd %xmm5,nb334_fjzH2(%esp)
        movsd %xmm0,nb334_fixH1(%esp)
        movsd %xmm1,nb334_fiyH1(%esp)
        movsd %xmm2,nb334_fizH1(%esp)

        ## H1-M interaction 
        movapd nb334_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334_rsqH1M(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        leal (%eax,%eax,2),%eax ## idx*=3 (12 total now) 

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb334_two(%esp),%xmm7    ## two*Heps2 
        movapd nb334_qqMH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addsd  nb334_vctot(%esp),%xmm5
        movsd %xmm5,nb334_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb334_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2


        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb334_dxH1M(%esp),%xmm0
        mulsd nb334_dyH1M(%esp),%xmm1
        mulsd nb334_dzH1M(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb334_fixH1(%esp),%xmm0
        addsd nb334_fiyH1(%esp),%xmm1
        addsd nb334_fizH1(%esp),%xmm2
        movsd %xmm3,nb334_fjxM(%esp)
        movsd %xmm4,nb334_fjyM(%esp)
        movsd %xmm5,nb334_fjzM(%esp)
        movsd %xmm0,nb334_fixH1(%esp)
        movsd %xmm1,nb334_fiyH1(%esp)
        movsd %xmm2,nb334_fizH1(%esp)

        ## H2-H1 interaction 
        movapd nb334_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb334_two(%esp),%xmm7    ## two*Heps2 
        movapd nb334_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addsd  nb334_vctot(%esp),%xmm5
        movsd %xmm5,nb334_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb334_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb334_fjxH1(%esp),%xmm3
        movapd nb334_fjyH1(%esp),%xmm4
        movapd nb334_fjzH1(%esp),%xmm5
        mulsd nb334_dxH2H1(%esp),%xmm0
        mulsd nb334_dyH2H1(%esp),%xmm1
        mulsd nb334_dzH2H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb334_fixH2(%esp),%xmm0
        addsd nb334_fiyH2(%esp),%xmm1
        addsd nb334_fizH2(%esp),%xmm2
        movsd %xmm3,nb334_fjxH1(%esp)
        movsd %xmm4,nb334_fjyH1(%esp)
        movsd %xmm5,nb334_fjzH1(%esp)
        movsd %xmm0,nb334_fixH2(%esp)
        movsd %xmm1,nb334_fiyH2(%esp)
        movsd %xmm2,nb334_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb334_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb334_two(%esp),%xmm7    ## two*Heps2 
        movapd nb334_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addsd  nb334_vctot(%esp),%xmm5
        movsd %xmm5,nb334_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb334_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb334_fjxH2(%esp),%xmm3
        movapd nb334_fjyH2(%esp),%xmm4
        movapd nb334_fjzH2(%esp),%xmm5
        mulsd nb334_dxH2H2(%esp),%xmm0
        mulsd nb334_dyH2H2(%esp),%xmm1
        mulsd nb334_dzH2H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb334_fixH2(%esp),%xmm0
        addsd nb334_fiyH2(%esp),%xmm1
        addsd nb334_fizH2(%esp),%xmm2
        movsd %xmm3,nb334_fjxH2(%esp)
        movsd %xmm4,nb334_fjyH2(%esp)
        movsd %xmm5,nb334_fjzH2(%esp)
        movsd %xmm0,nb334_fixH2(%esp)
        movsd %xmm1,nb334_fiyH2(%esp)
        movsd %xmm2,nb334_fizH2(%esp)

        ## H2-M interaction 
        movapd nb334_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334_rsqH2M(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb334_two(%esp),%xmm7    ## two*Heps2 
        movapd nb334_qqMH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addsd  nb334_vctot(%esp),%xmm5
        movsd %xmm5,nb334_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb334_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb334_fjxM(%esp),%xmm3
        movapd nb334_fjyM(%esp),%xmm4
        movapd nb334_fjzM(%esp),%xmm5
        mulsd nb334_dxH2M(%esp),%xmm0
        mulsd nb334_dyH2M(%esp),%xmm1
        mulsd nb334_dzH2M(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb334_fixH2(%esp),%xmm0
        addsd nb334_fiyH2(%esp),%xmm1
        addsd nb334_fizH2(%esp),%xmm2
        movsd %xmm3,nb334_fjxM(%esp)
        movsd %xmm4,nb334_fjyM(%esp)
        movsd %xmm5,nb334_fjzM(%esp)
        movsd %xmm0,nb334_fixH2(%esp)
        movsd %xmm1,nb334_fiyH2(%esp)
        movsd %xmm2,nb334_fizH2(%esp)

        ## M-H1 interaction 
        movapd nb334_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334_rsqMH1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb334_two(%esp),%xmm7    ## two*Heps2 
        movapd nb334_qqMH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addsd  nb334_vctot(%esp),%xmm5
        movsd %xmm5,nb334_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb334_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb334_fjxH1(%esp),%xmm3
        movapd nb334_fjyH1(%esp),%xmm4
        movapd nb334_fjzH1(%esp),%xmm5
        mulsd nb334_dxMH1(%esp),%xmm0
        mulsd nb334_dyMH1(%esp),%xmm1
        mulsd nb334_dzMH1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb334_fixM(%esp),%xmm0
        addsd nb334_fiyM(%esp),%xmm1
        addsd nb334_fizM(%esp),%xmm2
        movsd %xmm3,nb334_fjxH1(%esp)
        movsd %xmm4,nb334_fjyH1(%esp)
        movsd %xmm5,nb334_fjzH1(%esp)
        movsd %xmm0,nb334_fixM(%esp)
        movsd %xmm1,nb334_fiyM(%esp)
        movsd %xmm2,nb334_fizM(%esp)

        ## M-H2 interaction 
        movapd nb334_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334_rsqMH2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb334_two(%esp),%xmm7    ## two*Heps2 
        movapd nb334_qqMH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addsd  nb334_vctot(%esp),%xmm5
        movsd %xmm5,nb334_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb334_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb334_fjxH2(%esp),%xmm3
        movapd nb334_fjyH2(%esp),%xmm4
        movapd nb334_fjzH2(%esp),%xmm5
        mulsd nb334_dxMH2(%esp),%xmm0
        mulsd nb334_dyMH2(%esp),%xmm1
        mulsd nb334_dzMH2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb334_fixM(%esp),%xmm0
        addsd nb334_fiyM(%esp),%xmm1
        addsd nb334_fizM(%esp),%xmm2
        movsd %xmm3,nb334_fjxH2(%esp)
        movsd %xmm4,nb334_fjyH2(%esp)
        movsd %xmm5,nb334_fjzH2(%esp)
        movsd %xmm0,nb334_fixM(%esp)
        movsd %xmm1,nb334_fiyM(%esp)
        movsd %xmm2,nb334_fizM(%esp)

        ## M-M interaction 
        movapd nb334_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334_rsqMM(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb334_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb334_two(%esp),%xmm7    ## two*Heps2 
        movapd nb334_qqMM(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addsd  nb334_vctot(%esp),%xmm5
        movsd %xmm5,nb334_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb334_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb334_fjxM(%esp),%xmm3
        movapd nb334_fjyM(%esp),%xmm4
        movapd nb334_fjzM(%esp),%xmm5
        mulsd nb334_dxMM(%esp),%xmm0
        mulsd nb334_dyMM(%esp),%xmm1
        mulsd nb334_dzMM(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb334_fixM(%esp),%xmm0
        addsd nb334_fiyM(%esp),%xmm1
        addsd nb334_fizM(%esp),%xmm2
        movsd %xmm3,nb334_fjxM(%esp)
        movsd %xmm4,nb334_fjyM(%esp)
        movsd %xmm5,nb334_fjzM(%esp)
        movsd %xmm0,nb334_fixM(%esp)
        movsd %xmm1,nb334_fiyM(%esp)
        movsd %xmm2,nb334_fizM(%esp)

        movl nb334_faction(%ebp),%edi

        movd %mm0,%eax

        ## Did all interactions - now update j forces 
        ## Step1 - merge forces
        movlpd nb334_fjxO(%esp),%xmm0
        movlpd nb334_fjzO(%esp),%xmm1
        movlpd nb334_fjyH1(%esp),%xmm2
        movlpd nb334_fjxH2(%esp),%xmm3
        movlpd nb334_fjzH2(%esp),%xmm4
        movlpd nb334_fjyM(%esp),%xmm5

        movhpd nb334_fjyO(%esp),%xmm0
        movhpd nb334_fjxH1(%esp),%xmm1
        movhpd nb334_fjzH1(%esp),%xmm2
        movhpd nb334_fjyH2(%esp),%xmm3
        movhpd nb334_fjxM(%esp),%xmm4
        movhpd nb334_fjzM(%esp),%xmm5

        movlps (%edi,%eax,8),%xmm6
        movhps 8(%edi,%eax,8),%xmm6
        movlps 16(%edi,%eax,8),%xmm7
        movhps 24(%edi,%eax,8),%xmm7
        addpd  %xmm6,%xmm0
        addpd  %xmm7,%xmm1
        movlps 32(%edi,%eax,8),%xmm6
        movhps 40(%edi,%eax,8),%xmm6
        movlps 48(%edi,%eax,8),%xmm7
        movhps 56(%edi,%eax,8),%xmm7
        addpd  %xmm6,%xmm2
        addpd  %xmm7,%xmm3
        movlps 64(%edi,%eax,8),%xmm6
        movhps 72(%edi,%eax,8),%xmm6
        movlps 80(%edi,%eax,8),%xmm7
        movhps 88(%edi,%eax,8),%xmm7
        addpd  %xmm6,%xmm4
        addpd  %xmm7,%xmm5

        movlpd %xmm0,(%edi,%eax,8)
        movhpd %xmm0,8(%edi,%eax,8)
        movlpd %xmm1,16(%edi,%eax,8)
        movhpd %xmm1,24(%edi,%eax,8)
        movlpd %xmm2,32(%edi,%eax,8)
        movhpd %xmm2,40(%edi,%eax,8)
        movlpd %xmm3,48(%edi,%eax,8)
        movhpd %xmm3,56(%edi,%eax,8)
        movlpd %xmm4,64(%edi,%eax,8)
        movhpd %xmm4,72(%edi,%eax,8)
        movlpd %xmm5,80(%edi,%eax,8)
        movhpd %xmm5,88(%edi,%eax,8)


_nb_kernel334_ia32_sse2.nb334_updateouterdata: 
        movl  nb334_ii3(%esp),%ecx
        movl  nb334_faction(%ebp),%edi
        movl  nb334_fshift(%ebp),%esi
        movl  nb334_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb334_fixO(%esp),%xmm0
        movapd nb334_fiyO(%esp),%xmm1
        movapd nb334_fizO(%esp),%xmm2

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
        movapd nb334_fixH1(%esp),%xmm0
        movapd nb334_fiyH1(%esp),%xmm1
        movapd nb334_fizH1(%esp),%xmm2

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
        movapd nb334_fixH2(%esp),%xmm0
        movapd nb334_fiyH2(%esp),%xmm1
        movapd nb334_fizH2(%esp),%xmm2

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

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movapd nb334_fixM(%esp),%xmm0
        movapd nb334_fiyM(%esp),%xmm1
        movapd nb334_fizM(%esp),%xmm2

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
        movsd  72(%edi,%ecx,8),%xmm3
        movsd  80(%edi,%ecx,8),%xmm4
        movsd  88(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,72(%edi,%ecx,8)
        movsd  %xmm4,80(%edi,%ecx,8)
        movsd  %xmm5,88(%edi,%ecx,8)

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
        movl nb334_n(%esp),%esi
        ## get group index for i particle 
        movl  nb334_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb334_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb334_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb334_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb334_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb334_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel334_ia32_sse2.nb334_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb334_n(%esp)
        jmp _nb_kernel334_ia32_sse2.nb334_outer
_nb_kernel334_ia32_sse2.nb334_outerend: 
        ## check if more outer neighborlists remain
        movl  nb334_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel334_ia32_sse2.nb334_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel334_ia32_sse2.nb334_threadloop
_nb_kernel334_ia32_sse2.nb334_end: 
        emms

        movl nb334_nouter(%esp),%eax
        movl nb334_ninner(%esp),%ebx
        movl nb334_outeriter(%ebp),%ecx
        movl nb334_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb334_salign(%esp),%eax
        addl %eax,%esp
        addl $1800,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



.globl nb_kernel334nf_ia32_sse2
.globl _nb_kernel334nf_ia32_sse2
nb_kernel334nf_ia32_sse2:       
_nb_kernel334nf_ia32_sse2:      
.set nb334nf_p_nri, 8
.set nb334nf_iinr, 12
.set nb334nf_jindex, 16
.set nb334nf_jjnr, 20
.set nb334nf_shift, 24
.set nb334nf_shiftvec, 28
.set nb334nf_fshift, 32
.set nb334nf_gid, 36
.set nb334nf_pos, 40
.set nb334nf_faction, 44
.set nb334nf_charge, 48
.set nb334nf_p_facel, 52
.set nb334nf_argkrf, 56
.set nb334nf_argcrf, 60
.set nb334nf_Vc, 64
.set nb334nf_type, 68
.set nb334nf_p_ntype, 72
.set nb334nf_vdwparam, 76
.set nb334nf_Vvdw, 80
.set nb334nf_p_tabscale, 84
.set nb334nf_VFtab, 88
.set nb334nf_invsqrta, 92
.set nb334nf_dvda, 96
.set nb334nf_p_gbtabscale, 100
.set nb334nf_GBtab, 104
.set nb334nf_p_nthreads, 108
.set nb334nf_count, 112
.set nb334nf_mtx, 116
.set nb334nf_outeriter, 120
.set nb334nf_inneriter, 124
.set nb334nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb334nf_ixO, 0
.set nb334nf_iyO, 16
.set nb334nf_izO, 32
.set nb334nf_ixH1, 48
.set nb334nf_iyH1, 64
.set nb334nf_izH1, 80
.set nb334nf_ixH2, 96
.set nb334nf_iyH2, 112
.set nb334nf_izH2, 128
.set nb334nf_ixM, 144
.set nb334nf_iyM, 160
.set nb334nf_izM, 176
.set nb334nf_jxO, 192
.set nb334nf_jyO, 208
.set nb334nf_jzO, 224
.set nb334nf_jxH1, 240
.set nb334nf_jyH1, 256
.set nb334nf_jzH1, 272
.set nb334nf_jxH2, 288
.set nb334nf_jyH2, 304
.set nb334nf_jzH2, 320
.set nb334nf_jxM, 336
.set nb334nf_jyM, 352
.set nb334nf_jzM, 368
.set nb334nf_qqMM, 384
.set nb334nf_qqMH, 400
.set nb334nf_qqHH, 416
.set nb334nf_tsc, 432
.set nb334nf_c6, 448
.set nb334nf_c12, 464
.set nb334nf_vctot, 480
.set nb334nf_Vvdwtot, 496
.set nb334nf_half, 512
.set nb334nf_three, 528
.set nb334nf_rsqOO, 544
.set nb334nf_rsqH1H1, 560
.set nb334nf_rsqH1H2, 576
.set nb334nf_rsqH1M, 592
.set nb334nf_rsqH2H1, 608
.set nb334nf_rsqH2H2, 624
.set nb334nf_rsqH2M, 640
.set nb334nf_rsqMH1, 656
.set nb334nf_rsqMH2, 672
.set nb334nf_rsqMM, 688
.set nb334nf_rinvOO, 704
.set nb334nf_rinvH1H1, 720
.set nb334nf_rinvH1H2, 736
.set nb334nf_rinvH1M, 752
.set nb334nf_rinvH2H1, 768
.set nb334nf_rinvH2H2, 784
.set nb334nf_rinvH2M, 800
.set nb334nf_rinvMH1, 816
.set nb334nf_rinvMH2, 832
.set nb334nf_rinvMM, 848
.set nb334nf_is3, 864
.set nb334nf_ii3, 868
.set nb334nf_innerjjnr, 872
.set nb334nf_innerk, 876
.set nb334nf_n, 880
.set nb334nf_nn1, 884
.set nb334nf_nri, 888
.set nb334nf_nouter, 892
.set nb334nf_ninner, 896
.set nb334nf_salign, 900
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $904,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb334nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb334nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb334nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb334nf_nouter(%esp)
        movl %eax,nb334nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb334nf_half(%esp)
        movl %ebx,nb334nf_half+4(%esp)
        movsd nb334nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb334nf_half(%esp)
        movapd %xmm3,nb334nf_three(%esp)
        movl nb334nf_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3

        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb334nf_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb334nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb334nf_charge(%ebp),%edx
        movsd 24(%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb334nf_p_facel(%ebp),%esi
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
        movapd %xmm3,nb334nf_qqMM(%esp)
        movapd %xmm4,nb334nf_qqMH(%esp)
        movapd %xmm5,nb334nf_qqHH(%esp)

        xorpd %xmm0,%xmm0
        movl  nb334nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb334nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb334nf_vdwparam(%ebp),%eax
        movlpd (%eax,%edx,8),%xmm0
        movlpd 8(%eax,%edx,8),%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb334nf_c6(%esp)
        movapd %xmm1,nb334nf_c12(%esp)

_nb_kernel334nf_ia32_sse2.nb334nf_threadloop: 
        movl  nb334nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel334nf_ia32_sse2.nb334nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel334nf_ia32_sse2.nb334nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb334nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb334nf_n(%esp)
        movl %ebx,nb334nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel334nf_ia32_sse2.nb334nf_outerstart
        jmp _nb_kernel334nf_ia32_sse2.nb334nf_end

_nb_kernel334nf_ia32_sse2.nb334nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb334nf_nouter(%esp),%ebx
        movl %ebx,nb334nf_nouter(%esp)

_nb_kernel334nf_ia32_sse2.nb334nf_outer: 
        movl  nb334nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb334nf_is3(%esp)            ## store is3 

        movl  nb334nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb334nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb334nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb334nf_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3       ## ox
        addsd 8(%eax,%ebx,8),%xmm4      ## oy
        addsd 16(%eax,%ebx,8),%xmm5     ## oz   
        addsd 24(%eax,%ebx,8),%xmm6     ## h1x
        addsd 32(%eax,%ebx,8),%xmm7     ## h1y
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm7,%xmm7
        movapd %xmm3,nb334nf_ixO(%esp)
        movapd %xmm4,nb334nf_iyO(%esp)
        movapd %xmm5,nb334nf_izO(%esp)
        movapd %xmm6,nb334nf_ixH1(%esp)
        movapd %xmm7,nb334nf_iyH1(%esp)

        movsd %xmm2,%xmm6
        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 40(%eax,%ebx,8),%xmm6    ## h1z
        addsd 48(%eax,%ebx,8),%xmm0    ## h2x
        addsd 56(%eax,%ebx,8),%xmm1    ## h2y
        addsd 64(%eax,%ebx,8),%xmm2    ## h2z
        addsd 72(%eax,%ebx,8),%xmm3    ## mx
        addsd 80(%eax,%ebx,8),%xmm4    ## my
        addsd 88(%eax,%ebx,8),%xmm5    ## mz

        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm6,nb334nf_izH1(%esp)
        movapd %xmm0,nb334nf_ixH2(%esp)
        movapd %xmm1,nb334nf_iyH2(%esp)
        movapd %xmm2,nb334nf_izH2(%esp)
        movapd %xmm3,nb334nf_ixM(%esp)
        movapd %xmm4,nb334nf_iyM(%esp)
        movapd %xmm5,nb334nf_izM(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb334nf_vctot(%esp)
        movapd %xmm4,nb334nf_Vvdwtot(%esp)

        movl  nb334nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb334nf_pos(%ebp),%esi
        movl  nb334nf_faction(%ebp),%edi
        movl  nb334nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb334nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb334nf_ninner(%esp),%ecx
        movl  %ecx,nb334nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb334nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel334nf_ia32_sse2.nb334nf_unroll_loop
        jmp   _nb_kernel334nf_ia32_sse2.nb334nf_checksingle
_nb_kernel334nf_ia32_sse2.nb334nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb334nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb334nf_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb334nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%esi,%eax,8),%xmm0
        movlpd (%esi,%ebx,8),%xmm2
        movhpd 8(%esi,%eax,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm2
        movlpd 16(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%ebx,8),%xmm5
        movhpd 24(%esi,%eax,8),%xmm3
        movhpd 24(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## ox 
        unpckhpd %xmm2,%xmm1 ## oy
        unpcklpd %xmm5,%xmm3 ## ox 
        unpckhpd %xmm5,%xmm4 ## oy
        movapd  %xmm0,nb334nf_jxO(%esp)
        movapd  %xmm1,nb334nf_jyO(%esp)
        movapd  %xmm3,nb334nf_jzO(%esp)
        movapd  %xmm4,nb334nf_jxH1(%esp)

        ## load h1y, h1z, h2x, h2y 
        movlpd 32(%esi,%eax,8),%xmm0
        movlpd 32(%esi,%ebx,8),%xmm2
        movhpd 40(%esi,%eax,8),%xmm0
        movhpd 40(%esi,%ebx,8),%xmm2
        movlpd 48(%esi,%eax,8),%xmm3
        movlpd 48(%esi,%ebx,8),%xmm5
        movhpd 56(%esi,%eax,8),%xmm3
        movhpd 56(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h1y
        unpckhpd %xmm2,%xmm1 ## h1z
        unpcklpd %xmm5,%xmm3 ## h2x
        unpckhpd %xmm5,%xmm4 ## h2y
        movapd  %xmm0,nb334nf_jyH1(%esp)
        movapd  %xmm1,nb334nf_jzH1(%esp)
        movapd  %xmm3,nb334nf_jxH2(%esp)
        movapd  %xmm4,nb334nf_jyH2(%esp)

        ## load h2z, mx, my, mz
        movlpd 64(%esi,%eax,8),%xmm0
        movlpd 64(%esi,%ebx,8),%xmm2
        movhpd 72(%esi,%eax,8),%xmm0
        movhpd 72(%esi,%ebx,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 80(%esi,%ebx,8),%xmm5
        movhpd 88(%esi,%eax,8),%xmm3
        movhpd 88(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h2z
        unpckhpd %xmm2,%xmm1 ## mx
        unpcklpd %xmm5,%xmm3 ## my
        unpckhpd %xmm5,%xmm4 ## mz
        movapd  %xmm0,nb334nf_jzH2(%esp)
        movapd  %xmm1,nb334nf_jxM(%esp)
        movapd  %xmm3,nb334nf_jyM(%esp)
        movapd  %xmm4,nb334nf_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb334nf_ixO(%esp),%xmm0
        movapd nb334nf_iyO(%esp),%xmm1
        movapd nb334nf_izO(%esp),%xmm2
        movapd nb334nf_ixH1(%esp),%xmm3
        movapd nb334nf_iyH1(%esp),%xmm4
        movapd nb334nf_izH1(%esp),%xmm5
        subpd  nb334nf_jxO(%esp),%xmm0
        subpd  nb334nf_jyO(%esp),%xmm1
        subpd  nb334nf_jzO(%esp),%xmm2
        subpd  nb334nf_jxH1(%esp),%xmm3
        subpd  nb334nf_jyH1(%esp),%xmm4
        subpd  nb334nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb334nf_rsqOO(%esp)
        movapd %xmm3,nb334nf_rsqH1H1(%esp)

        movapd nb334nf_ixH1(%esp),%xmm0
        movapd nb334nf_iyH1(%esp),%xmm1
        movapd nb334nf_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb334nf_jxH2(%esp),%xmm0
        subpd  nb334nf_jyH2(%esp),%xmm1
        subpd  nb334nf_jzH2(%esp),%xmm2
        subpd  nb334nf_jxM(%esp),%xmm3
        subpd  nb334nf_jyM(%esp),%xmm4
        subpd  nb334nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb334nf_rsqH1H2(%esp)
        movapd %xmm3,nb334nf_rsqH1M(%esp)

        movapd nb334nf_ixH2(%esp),%xmm0
        movapd nb334nf_iyH2(%esp),%xmm1
        movapd nb334nf_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb334nf_jxH1(%esp),%xmm0
        subpd  nb334nf_jyH1(%esp),%xmm1
        subpd  nb334nf_jzH1(%esp),%xmm2
        subpd  nb334nf_jxH2(%esp),%xmm3
        subpd  nb334nf_jyH2(%esp),%xmm4
        subpd  nb334nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb334nf_rsqH2H1(%esp)
        movapd %xmm3,nb334nf_rsqH2H2(%esp)

        movapd nb334nf_ixH2(%esp),%xmm0
        movapd nb334nf_iyH2(%esp),%xmm1
        movapd nb334nf_izH2(%esp),%xmm2
        movapd nb334nf_ixM(%esp),%xmm3
        movapd nb334nf_iyM(%esp),%xmm4
        movapd nb334nf_izM(%esp),%xmm5
        subpd  nb334nf_jxM(%esp),%xmm0
        subpd  nb334nf_jyM(%esp),%xmm1
        subpd  nb334nf_jzM(%esp),%xmm2
        subpd  nb334nf_jxH1(%esp),%xmm3
        subpd  nb334nf_jyH1(%esp),%xmm4
        subpd  nb334nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb334nf_rsqH2M(%esp)
        movapd %xmm4,nb334nf_rsqMH1(%esp)

        movapd nb334nf_ixM(%esp),%xmm0
        movapd nb334nf_iyM(%esp),%xmm1
        movapd nb334nf_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb334nf_jxH2(%esp),%xmm0
        subpd  nb334nf_jyH2(%esp),%xmm1
        subpd  nb334nf_jzH2(%esp),%xmm2
        subpd  nb334nf_jxM(%esp),%xmm3
        subpd  nb334nf_jyM(%esp),%xmm4
        subpd  nb334nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb334nf_rsqMH2(%esp)
        movapd %xmm4,nb334nf_rsqMM(%esp)

        ## Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1  ## luA
        cvtps2pd %xmm5,%xmm5  ## luB

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb334nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb334nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb334nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb334nf_rinvMH2(%esp)
        movapd %xmm5,nb334nf_rinvMM(%esp)

        movapd nb334nf_rsqOO(%esp),%xmm0
        movapd nb334nf_rsqH1H1(%esp),%xmm4
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
        movapd  nb334nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334nf_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb334nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb334nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb334nf_rinvOO(%esp)
        movapd %xmm5,nb334nf_rinvH1H1(%esp)

        movapd nb334nf_rsqH1H2(%esp),%xmm0
        movapd nb334nf_rsqH1M(%esp),%xmm4
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
        movapd  nb334nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb334nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb334nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb334nf_rinvH1H2(%esp)
        movapd %xmm5,nb334nf_rinvH1M(%esp)

        movapd nb334nf_rsqH2H1(%esp),%xmm0
        movapd nb334nf_rsqH2H2(%esp),%xmm4
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
        movapd  nb334nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb334nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb334nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb334nf_rinvH2H1(%esp)
        movapd %xmm5,nb334nf_rinvH2H2(%esp)

        movapd nb334nf_rsqMH1(%esp),%xmm0
        movapd nb334nf_rsqH2M(%esp),%xmm4
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
        movapd  nb334nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb334nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb334nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb334nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb334nf_rinvMH1(%esp)
        movapd %xmm5,nb334nf_rinvH2M(%esp)

        ## start with OO interaction 
        movapd nb334nf_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334nf_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334nf_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movd %eax,%mm0
        movd %ebx,%mm1
        movl nb334nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

        ## Dispersion 
        movlpd 32(%esi,%eax,8),%xmm4    ## Y1
        movlpd 32(%esi,%ebx,8),%xmm3    ## Y2
        movhpd 40(%esi,%eax,8),%xmm4    ## Y1 F1        
        movhpd 40(%esi,%ebx,8),%xmm3    ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%esi,%eax,8),%xmm6    ## G1
        movlpd 48(%esi,%ebx,8),%xmm3    ## G2
        movhpd 56(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 56(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb334nf_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm5       ## Vvdw6 
        addpd  nb334nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb334nf_Vvdwtot(%esp)

        ## Repulsion 
        movlpd 64(%esi,%eax,8),%xmm4    ## Y1
        movlpd 64(%esi,%ebx,8),%xmm3    ## Y2
        movhpd 72(%esi,%eax,8),%xmm4    ## Y1 F1        
        movhpd 72(%esi,%ebx,8),%xmm3    ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 80(%esi,%eax,8),%xmm6    ## G1
        movlpd 80(%esi,%ebx,8),%xmm3    ## G2
        movhpd 88(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 88(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb334nf_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm5 ## Vvdw12 
        addpd  nb334nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb334nf_Vvdwtot(%esp)

        ## H1-H1 interaction 
        movapd nb334nf_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334nf_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334nf_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

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
        movapd nb334nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        addpd  nb334nf_vctot(%esp),%xmm5
        movapd %xmm5,nb334nf_vctot(%esp)

        ## H1-H2 interaction  
        movapd nb334nf_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334nf_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334nf_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

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
        movapd nb334nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        addpd  nb334nf_vctot(%esp),%xmm5
        movapd %xmm5,nb334nf_vctot(%esp)

        ## H1-M interaction 
        movapd nb334nf_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334nf_rsqH1M(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334nf_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

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
        movapd nb334nf_qqMH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb334nf_vctot(%esp),%xmm5
        movapd %xmm5,nb334nf_vctot(%esp)

        ## H2-H1 interaction 
        movapd nb334nf_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334nf_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

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
        movapd nb334nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb334nf_vctot(%esp),%xmm5
        movapd %xmm5,nb334nf_vctot(%esp)

        ## H2-H2 interaction 
        movapd nb334nf_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334nf_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

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
        movapd nb334nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb334nf_vctot(%esp),%xmm5
        movapd %xmm5,nb334nf_vctot(%esp)

        ## H2-M interaction 
        movapd nb334nf_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334nf_rsqH2M(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

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
        movapd nb334nf_qqMH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb334nf_vctot(%esp),%xmm5
        movapd %xmm5,nb334nf_vctot(%esp)

        ## M-H1 interaction 
        movapd nb334nf_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334nf_rsqMH1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

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
        movapd nb334nf_qqMH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb334nf_vctot(%esp),%xmm5
        movapd %xmm5,nb334nf_vctot(%esp)

        ## M-H2 interaction 
        movapd nb334nf_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334nf_rsqMH2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

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
        movapd nb334nf_qqMH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb334nf_vctot(%esp),%xmm5
        movapd %xmm5,nb334nf_vctot(%esp)

        ## M-M interaction 
        movapd nb334nf_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb334nf_rsqMM(%esp),%xmm1   ## xmm1=r 
        mulpd  nb334nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 
        leal  (%ebx,%ebx,2),%ebx        ## idx*=3 (12 total now) 

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
        movapd nb334nf_qqMM(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb334nf_vctot(%esp),%xmm5
        movapd %xmm5,nb334nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb334nf_innerk(%esp)
        jl    _nb_kernel334nf_ia32_sse2.nb334nf_checksingle
        jmp   _nb_kernel334nf_ia32_sse2.nb334nf_unroll_loop
_nb_kernel334nf_ia32_sse2.nb334nf_checksingle: 
        movl  nb334nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel334nf_ia32_sse2.nb334nf_dosingle
        jmp   _nb_kernel334nf_ia32_sse2.nb334nf_updateouterdata
_nb_kernel334nf_ia32_sse2.nb334nf_dosingle: 
        movl  nb334nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb334nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%esi,%eax,8),%xmm0
        movhpd 8(%esi,%eax,8),%xmm0
        movlpd 16(%esi,%eax,8),%xmm1
        movhpd 24(%esi,%eax,8),%xmm1
        movlpd 32(%esi,%eax,8),%xmm2
        movhpd 40(%esi,%eax,8),%xmm2
        movlpd 48(%esi,%eax,8),%xmm3
        movhpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 72(%esi,%eax,8),%xmm4
        movlpd 80(%esi,%eax,8),%xmm5
        movhpd 88(%esi,%eax,8),%xmm5
        movsd  %xmm0,nb334nf_jxO(%esp)
        movsd  %xmm1,nb334nf_jzO(%esp)
        movsd  %xmm2,nb334nf_jyH1(%esp)
        movsd  %xmm3,nb334nf_jxH2(%esp)
        movsd  %xmm4,nb334nf_jzH2(%esp)
        movsd  %xmm5,nb334nf_jyM(%esp)
        unpckhpd %xmm0,%xmm0
        unpckhpd %xmm1,%xmm1
        unpckhpd %xmm2,%xmm2
        unpckhpd %xmm3,%xmm3
        unpckhpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5
        movsd  %xmm0,nb334nf_jyO(%esp)
        movsd  %xmm1,nb334nf_jxH1(%esp)
        movsd  %xmm2,nb334nf_jzH1(%esp)
        movsd  %xmm3,nb334nf_jyH2(%esp)
        movsd  %xmm4,nb334nf_jxM(%esp)
        movsd  %xmm5,nb334nf_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb334nf_ixO(%esp),%xmm0
        movapd nb334nf_iyO(%esp),%xmm1
        movapd nb334nf_izO(%esp),%xmm2
        movapd nb334nf_ixH1(%esp),%xmm3
        movapd nb334nf_iyH1(%esp),%xmm4
        movapd nb334nf_izH1(%esp),%xmm5
        subsd  nb334nf_jxO(%esp),%xmm0
        subsd  nb334nf_jyO(%esp),%xmm1
        subsd  nb334nf_jzO(%esp),%xmm2
        subsd  nb334nf_jxH1(%esp),%xmm3
        subsd  nb334nf_jyH1(%esp),%xmm4
        subsd  nb334nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb334nf_rsqOO(%esp)
        movapd %xmm3,nb334nf_rsqH1H1(%esp)

        movapd nb334nf_ixH1(%esp),%xmm0
        movapd nb334nf_iyH1(%esp),%xmm1
        movapd nb334nf_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb334nf_jxH2(%esp),%xmm0
        subsd  nb334nf_jyH2(%esp),%xmm1
        subsd  nb334nf_jzH2(%esp),%xmm2
        subsd  nb334nf_jxM(%esp),%xmm3
        subsd  nb334nf_jyM(%esp),%xmm4
        subsd  nb334nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb334nf_rsqH1H2(%esp)
        movapd %xmm3,nb334nf_rsqH1M(%esp)

        movapd nb334nf_ixH2(%esp),%xmm0
        movapd nb334nf_iyH2(%esp),%xmm1
        movapd nb334nf_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb334nf_jxH1(%esp),%xmm0
        subsd  nb334nf_jyH1(%esp),%xmm1
        subsd  nb334nf_jzH1(%esp),%xmm2
        subsd  nb334nf_jxH2(%esp),%xmm3
        subsd  nb334nf_jyH2(%esp),%xmm4
        subsd  nb334nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb334nf_rsqH2H1(%esp)
        movapd %xmm3,nb334nf_rsqH2H2(%esp)

        movapd nb334nf_ixH2(%esp),%xmm0
        movapd nb334nf_iyH2(%esp),%xmm1
        movapd nb334nf_izH2(%esp),%xmm2
        movapd nb334nf_ixM(%esp),%xmm3
        movapd nb334nf_iyM(%esp),%xmm4
        movapd nb334nf_izM(%esp),%xmm5
        subsd  nb334nf_jxM(%esp),%xmm0
        subsd  nb334nf_jyM(%esp),%xmm1
        subsd  nb334nf_jzM(%esp),%xmm2
        subsd  nb334nf_jxH1(%esp),%xmm3
        subsd  nb334nf_jyH1(%esp),%xmm4
        subsd  nb334nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb334nf_rsqH2M(%esp)
        movapd %xmm4,nb334nf_rsqMH1(%esp)

        movapd nb334nf_ixM(%esp),%xmm0
        movapd nb334nf_iyM(%esp),%xmm1
        movapd nb334nf_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb334nf_jxH2(%esp),%xmm0
        subsd  nb334nf_jyH2(%esp),%xmm1
        subsd  nb334nf_jzH2(%esp),%xmm2
        subsd  nb334nf_jxM(%esp),%xmm3
        subsd  nb334nf_jyM(%esp),%xmm4
        subsd  nb334nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb334nf_rsqMH2(%esp)
        movapd %xmm4,nb334nf_rsqMM(%esp)

        ## Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1  ## luA
        cvtss2sd %xmm5,%xmm5  ## luB

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb334nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb334nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb334nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb334nf_rinvMH2(%esp)
        movapd %xmm5,nb334nf_rinvMM(%esp)

        movapd nb334nf_rsqOO(%esp),%xmm0
        movapd nb334nf_rsqH1H1(%esp),%xmm4
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
        movapd  nb334nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334nf_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb334nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb334nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb334nf_rinvOO(%esp)
        movapd %xmm5,nb334nf_rinvH1H1(%esp)

        movapd nb334nf_rsqH1H2(%esp),%xmm0
        movapd nb334nf_rsqH1M(%esp),%xmm4
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
        movapd  nb334nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb334nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb334nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb334nf_rinvH1H2(%esp)
        movapd %xmm5,nb334nf_rinvH1M(%esp)

        movapd nb334nf_rsqH2H1(%esp),%xmm0
        movapd nb334nf_rsqH2H2(%esp),%xmm4
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
        movapd  nb334nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb334nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb334nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb334nf_rinvH2H1(%esp)
        movapd %xmm5,nb334nf_rinvH2H2(%esp)

        movapd nb334nf_rsqMH1(%esp),%xmm0
        movapd nb334nf_rsqH2M(%esp),%xmm4
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
        movapd  nb334nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb334nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb334nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb334nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb334nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb334nf_rinvMH1(%esp)
        movapd %xmm5,nb334nf_rinvH2M(%esp)

        ## start with OO interaction 
        movapd nb334nf_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334nf_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334nf_tsc(%esp),%xmm1

        movd %eax,%mm0


        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 


        shll $2,%eax            ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        ## Dispersion 
        movsd 32(%esi,%eax,8),%xmm4     ## Y1   
        movsd 40(%esi,%eax,8),%xmm5     ## F1   
        movsd 48(%esi,%eax,8),%xmm6     ## G1   
        movsd 56(%esi,%eax,8),%xmm7     ## H1   
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb334nf_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        addsd  nb334nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb334nf_Vvdwtot(%esp)

        ## Repulsion 
        movsd 64(%esi,%eax,8),%xmm4     ## Y1   
        movsd 72(%esi,%eax,8),%xmm5     ## F1   
        movsd 80(%esi,%eax,8),%xmm6     ## G1   
        movsd 88(%esi,%eax,8),%xmm7     ## H1   

        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb334nf_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm5 ## Vvdw12 


        addsd  nb334nf_Vvdwtot(%esp),%xmm5
        movsd %xmm5,nb334nf_Vvdwtot(%esp)

        ## H1-H1 interaction 
        movapd nb334nf_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334nf_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334nf_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb334nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  


        addsd  nb334nf_vctot(%esp),%xmm5
        movsd %xmm5,nb334nf_vctot(%esp)

        ## H1-H2 interaction  
        movapd nb334nf_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334nf_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334nf_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb334nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        addsd  nb334nf_vctot(%esp),%xmm5
        movsd %xmm5,nb334nf_vctot(%esp)

        ## H1-M interaction 
        movapd nb334nf_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334nf_rsqH1M(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334nf_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        leal (%eax,%eax,2),%eax ## idx*=3 (12 total now) 

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb334nf_qqMH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb334nf_vctot(%esp),%xmm5
        movsd %xmm5,nb334nf_vctot(%esp)

        ## H2-H1 interaction 
        movapd nb334nf_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334nf_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334nf_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb334nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb334nf_vctot(%esp),%xmm5
        movsd %xmm5,nb334nf_vctot(%esp)

        ## H2-H2 interaction 
        movapd nb334nf_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334nf_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb334nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb334nf_vctot(%esp),%xmm5
        movsd %xmm5,nb334nf_vctot(%esp)

        ## H2-M interaction 
        movapd nb334nf_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334nf_rsqH2M(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb334nf_qqMH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb334nf_vctot(%esp),%xmm5
        movsd %xmm5,nb334nf_vctot(%esp)

        ## M-H1 interaction 
        movapd nb334nf_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334nf_rsqMH1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb334nf_qqMH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb334nf_vctot(%esp),%xmm5
        movsd %xmm5,nb334nf_vctot(%esp)

        ## M-H2 interaction 
        movapd nb334nf_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334nf_rsqMH2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb334nf_qqMH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb334nf_vctot(%esp),%xmm5
        movsd %xmm5,nb334nf_vctot(%esp)

        ## M-M interaction 
        movapd nb334nf_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb334nf_rsqMM(%esp),%xmm1   ## xmm1=r 
        mulsd  nb334nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb334nf_VFtab(%ebp),%esi
        leal  (%eax,%eax,2),%eax        ## idx*=3 (12 total now) 

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb334nf_qqMM(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb334nf_vctot(%esp),%xmm5
        movsd %xmm5,nb334nf_vctot(%esp)

_nb_kernel334nf_ia32_sse2.nb334nf_updateouterdata: 
        ## get n from stack
        movl nb334nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb334nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb334nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb334nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb334nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb334nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb334nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel334nf_ia32_sse2.nb334nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb334nf_n(%esp)
        jmp _nb_kernel334nf_ia32_sse2.nb334nf_outer
_nb_kernel334nf_ia32_sse2.nb334nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb334nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel334nf_ia32_sse2.nb334nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel334nf_ia32_sse2.nb334nf_threadloop
_nb_kernel334nf_ia32_sse2.nb334nf_end: 
        emms

        movl nb334nf_nouter(%esp),%eax
        movl nb334nf_ninner(%esp),%ebx
        movl nb334nf_outeriter(%ebp),%ecx
        movl nb334nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb334nf_salign(%esp),%eax
        addl %eax,%esp
        addl $904,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret

