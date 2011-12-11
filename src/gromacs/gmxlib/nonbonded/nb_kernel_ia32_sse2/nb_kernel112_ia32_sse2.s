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



.globl nb_kernel112_ia32_sse2
.globl _nb_kernel112_ia32_sse2
nb_kernel112_ia32_sse2: 
_nb_kernel112_ia32_sse2:        
.set nb112_p_nri, 8
.set nb112_iinr, 12
.set nb112_jindex, 16
.set nb112_jjnr, 20
.set nb112_shift, 24
.set nb112_shiftvec, 28
.set nb112_fshift, 32
.set nb112_gid, 36
.set nb112_pos, 40
.set nb112_faction, 44
.set nb112_charge, 48
.set nb112_p_facel, 52
.set nb112_argkrf, 56
.set nb112_argcrf, 60
.set nb112_Vc, 64
.set nb112_type, 68
.set nb112_p_ntype, 72
.set nb112_vdwparam, 76
.set nb112_Vvdw, 80
.set nb112_p_tabscale, 84
.set nb112_VFtab, 88
.set nb112_invsqrta, 92
.set nb112_dvda, 96
.set nb112_p_gbtabscale, 100
.set nb112_GBtab, 104
.set nb112_p_nthreads, 108
.set nb112_count, 112
.set nb112_mtx, 116
.set nb112_outeriter, 120
.set nb112_inneriter, 124
.set nb112_work, 128
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
.set nb112_innerjjnr, 1480
.set nb112_innerk, 1484
.set nb112_n, 1488
.set nb112_nn1, 1492
.set nb112_nri, 1496
.set nb112_nouter, 1500
.set nb112_ninner, 1504
.set nb112_salign, 1508
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $1512,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb112_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb112_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb112_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb112_nouter(%esp)
        movl %eax,nb112_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb112_half(%esp)
        movl %ebx,nb112_half+4(%esp)
        movsd nb112_half(%esp),%xmm1
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
        movapd %xmm1,nb112_half(%esp)
        movapd %xmm3,nb112_three(%esp)
        movapd %xmm4,nb112_six(%esp)
        movapd %xmm5,nb112_twelve(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb112_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb112_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb112_p_facel(%ebp),%esi
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
        movapd %xmm3,nb112_qqOO(%esp)
        movapd %xmm4,nb112_qqOH(%esp)
        movapd %xmm5,nb112_qqHH(%esp)

        xorpd %xmm0,%xmm0
        movl  nb112_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb112_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb112_vdwparam(%ebp),%eax
        movlpd (%eax,%edx,8),%xmm0
        movhpd 8(%eax,%edx,8),%xmm0
        movhlps %xmm0,%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb112_c6(%esp)
        movapd %xmm1,nb112_c12(%esp)

_nb_kernel112_ia32_sse2.nb112_threadloop: 
        movl  nb112_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel112_ia32_sse2.nb112_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel112_ia32_sse2.nb112_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb112_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb112_n(%esp)
        movl %ebx,nb112_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel112_ia32_sse2.nb112_outerstart
        jmp _nb_kernel112_ia32_sse2.nb112_end

_nb_kernel112_ia32_sse2.nb112_outerstart: 
        ## ebx contains number of outer iterations
        addl nb112_nouter(%esp),%ebx
        movl %ebx,nb112_nouter(%esp)

_nb_kernel112_ia32_sse2.nb112_outer: 
        movl  nb112_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb112_is3(%esp)      ## store is3 

        movl  nb112_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb112_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb112_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb112_ii3(%esp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb112_ixO(%esp)
        movapd %xmm4,nb112_iyO(%esp)
        movapd %xmm5,nb112_izO(%esp)

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
        movapd %xmm0,nb112_ixH1(%esp)
        movapd %xmm1,nb112_iyH1(%esp)
        movapd %xmm2,nb112_izH1(%esp)
        movapd %xmm3,nb112_ixH2(%esp)
        movapd %xmm4,nb112_iyH2(%esp)
        movapd %xmm5,nb112_izH2(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb112_vctot(%esp)
        movapd %xmm4,nb112_Vvdwtot(%esp)
        movapd %xmm4,nb112_fixO(%esp)
        movapd %xmm4,nb112_fiyO(%esp)
        movapd %xmm4,nb112_fizO(%esp)
        movapd %xmm4,nb112_fixH1(%esp)
        movapd %xmm4,nb112_fiyH1(%esp)
        movapd %xmm4,nb112_fizH1(%esp)
        movapd %xmm4,nb112_fixH2(%esp)
        movapd %xmm4,nb112_fiyH2(%esp)
        movapd %xmm4,nb112_fizH2(%esp)

        movl  nb112_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb112_pos(%ebp),%esi
        movl  nb112_faction(%ebp),%edi
        movl  nb112_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb112_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb112_ninner(%esp),%ecx
        movl  %ecx,nb112_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb112_innerk(%esp)      ## number of innerloop atoms 
        jge  _nb_kernel112_ia32_sse2.nb112_unroll_loop
        jmp  _nb_kernel112_ia32_sse2.nb112_checksingle
_nb_kernel112_ia32_sse2.nb112_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb112_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb112_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb112_pos(%ebp),%esi        ## base of pos[] 

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
        movapd  %xmm2,nb112_jxO(%esp)
        movapd  %xmm3,nb112_jyO(%esp)
        movapd  %xmm4,nb112_jzO(%esp)
        movapd  %xmm5,nb112_jxH1(%esp)
        movapd  %xmm6,nb112_jyH1(%esp)
        movapd  %xmm7,nb112_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 48(%esi,%ebx,8),%xmm2
        movhpd 56(%esi,%ebx,8),%xmm3
        movhpd 64(%esi,%ebx,8),%xmm4
        movapd  %xmm2,nb112_jxH2(%esp)
        movapd  %xmm3,nb112_jyH2(%esp)
        movapd  %xmm4,nb112_jzH2(%esp)

        movapd nb112_ixO(%esp),%xmm0
        movapd nb112_iyO(%esp),%xmm1
        movapd nb112_izO(%esp),%xmm2
        movapd nb112_ixO(%esp),%xmm3
        movapd nb112_iyO(%esp),%xmm4
        movapd nb112_izO(%esp),%xmm5
        subpd  nb112_jxO(%esp),%xmm0
        subpd  nb112_jyO(%esp),%xmm1
        subpd  nb112_jzO(%esp),%xmm2
        subpd  nb112_jxH1(%esp),%xmm3
        subpd  nb112_jyH1(%esp),%xmm4
        subpd  nb112_jzH1(%esp),%xmm5
        movapd %xmm0,nb112_dxOO(%esp)
        movapd %xmm1,nb112_dyOO(%esp)
        movapd %xmm2,nb112_dzOO(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb112_dxOH1(%esp)
        movapd %xmm4,nb112_dyOH1(%esp)
        movapd %xmm5,nb112_dzOH1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb112_rsqOO(%esp)
        movapd %xmm3,nb112_rsqOH1(%esp)

        movapd nb112_ixO(%esp),%xmm0
        movapd nb112_iyO(%esp),%xmm1
        movapd nb112_izO(%esp),%xmm2
        movapd nb112_ixH1(%esp),%xmm3
        movapd nb112_iyH1(%esp),%xmm4
        movapd nb112_izH1(%esp),%xmm5
        subpd  nb112_jxH2(%esp),%xmm0
        subpd  nb112_jyH2(%esp),%xmm1
        subpd  nb112_jzH2(%esp),%xmm2
        subpd  nb112_jxO(%esp),%xmm3
        subpd  nb112_jyO(%esp),%xmm4
        subpd  nb112_jzO(%esp),%xmm5
        movapd %xmm0,nb112_dxOH2(%esp)
        movapd %xmm1,nb112_dyOH2(%esp)
        movapd %xmm2,nb112_dzOH2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb112_dxH1O(%esp)
        movapd %xmm4,nb112_dyH1O(%esp)
        movapd %xmm5,nb112_dzH1O(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb112_rsqOH2(%esp)
        movapd %xmm3,nb112_rsqH1O(%esp)

        movapd nb112_ixH1(%esp),%xmm0
        movapd nb112_iyH1(%esp),%xmm1
        movapd nb112_izH1(%esp),%xmm2
        movapd nb112_ixH1(%esp),%xmm3
        movapd nb112_iyH1(%esp),%xmm4
        movapd nb112_izH1(%esp),%xmm5
        subpd  nb112_jxH1(%esp),%xmm0
        subpd  nb112_jyH1(%esp),%xmm1
        subpd  nb112_jzH1(%esp),%xmm2
        subpd  nb112_jxH2(%esp),%xmm3
        subpd  nb112_jyH2(%esp),%xmm4
        subpd  nb112_jzH2(%esp),%xmm5
        movapd %xmm0,nb112_dxH1H1(%esp)
        movapd %xmm1,nb112_dyH1H1(%esp)
        movapd %xmm2,nb112_dzH1H1(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb112_dxH1H2(%esp)
        movapd %xmm4,nb112_dyH1H2(%esp)
        movapd %xmm5,nb112_dzH1H2(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb112_rsqH1H1(%esp)
        movapd %xmm3,nb112_rsqH1H2(%esp)

        movapd nb112_ixH2(%esp),%xmm0
        movapd nb112_iyH2(%esp),%xmm1
        movapd nb112_izH2(%esp),%xmm2
        movapd nb112_ixH2(%esp),%xmm3
        movapd nb112_iyH2(%esp),%xmm4
        movapd nb112_izH2(%esp),%xmm5
        subpd  nb112_jxO(%esp),%xmm0
        subpd  nb112_jyO(%esp),%xmm1
        subpd  nb112_jzO(%esp),%xmm2
        subpd  nb112_jxH1(%esp),%xmm3
        subpd  nb112_jyH1(%esp),%xmm4
        subpd  nb112_jzH1(%esp),%xmm5
        movapd %xmm0,nb112_dxH2O(%esp)
        movapd %xmm1,nb112_dyH2O(%esp)
        movapd %xmm2,nb112_dzH2O(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb112_dxH2H1(%esp)
        movapd %xmm4,nb112_dyH2H1(%esp)
        movapd %xmm5,nb112_dzH2H1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb112_rsqH2O(%esp)
        movapd %xmm4,nb112_rsqH2H1(%esp)

        movapd nb112_ixH2(%esp),%xmm0
        movapd nb112_iyH2(%esp),%xmm1
        movapd nb112_izH2(%esp),%xmm2
        subpd  nb112_jxH2(%esp),%xmm0
        subpd  nb112_jyH2(%esp),%xmm1
        subpd  nb112_jzH2(%esp),%xmm2
        movapd %xmm0,nb112_dxH2H2(%esp)
        movapd %xmm1,nb112_dyH2H2(%esp)
        movapd %xmm2,nb112_dzH2H2(%esp)
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb112_rsqH2H2(%esp)

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
        movapd  nb112_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112_half(%esp),%xmm3   ## iter1 
        mulpd   nb112_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112_half(%esp),%xmm1   ## rinv 
        mulpd   nb112_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb112_rinvH2H2(%esp)
        movapd %xmm5,nb112_rinvH2H1(%esp)

        movapd nb112_rsqOO(%esp),%xmm0
        movapd nb112_rsqOH1(%esp),%xmm4
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
        movapd  nb112_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb112_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112_half(%esp),%xmm1   ## rinv 
        mulpd   nb112_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb112_rinvOO(%esp)
        movapd %xmm5,nb112_rinvOH1(%esp)

        movapd nb112_rsqOH2(%esp),%xmm0
        movapd nb112_rsqH1O(%esp),%xmm4
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
        movapd  nb112_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112_half(%esp),%xmm3   ## iter1 
        mulpd   nb112_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112_half(%esp),%xmm1   ## rinv 
        mulpd   nb112_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb112_rinvOH2(%esp)
        movapd %xmm5,nb112_rinvH1O(%esp)

        movapd nb112_rsqH1H1(%esp),%xmm0
        movapd nb112_rsqH1H2(%esp),%xmm4
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
        movapd  nb112_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112_half(%esp),%xmm3   ## iter1a 
        mulpd   nb112_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112_half(%esp),%xmm1   ## rinv 
        mulpd   nb112_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb112_rinvH1H1(%esp)
        movapd %xmm5,nb112_rinvH1H2(%esp)

        movapd nb112_rsqH2O(%esp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb112_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb112_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb112_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb112_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb112_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb112_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm7
        mulpd  %xmm0,%xmm0
        movapd %xmm0,%xmm1
        mulpd  %xmm0,%xmm1
        mulpd  %xmm0,%xmm1      ## xmm1=rinvsix 
        mulpd  nb112_qqOO(%esp),%xmm7
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulpd  nb112_c6(%esp),%xmm1
        mulpd  nb112_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addpd  nb112_Vvdwtot(%esp),%xmm3
        mulpd  nb112_six(%esp),%xmm1
        mulpd  nb112_twelve(%esp),%xmm2
        movapd %xmm3,nb112_Vvdwtot(%esp)
        subpd  %xmm1,%xmm2
        addpd  %xmm7,%xmm2
        addpd  nb112_vctot(%esp),%xmm7
        mulpd  %xmm2,%xmm0

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb112_dxOO(%esp),%xmm0
        mulpd nb112_dyOO(%esp),%xmm1
        mulpd nb112_dzOO(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb112_fixO(%esp),%xmm0
        addpd nb112_fiyO(%esp),%xmm1
        addpd nb112_fizO(%esp),%xmm2
        movapd %xmm3,nb112_fjxO(%esp)
        movapd %xmm4,nb112_fjyO(%esp)
        movapd %xmm5,nb112_fjzO(%esp)
        movapd %xmm0,nb112_fixO(%esp)
        movapd %xmm1,nb112_fiyO(%esp)
        movapd %xmm2,nb112_fizO(%esp)

        ## O-H1 interaction 
        movapd nb112_rinvOH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb112_qqOH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsOH1  
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb112_dxOH1(%esp),%xmm0
        mulpd nb112_dyOH1(%esp),%xmm1
        mulpd nb112_dzOH1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb112_fixO(%esp),%xmm0
        addpd nb112_fiyO(%esp),%xmm1
        addpd nb112_fizO(%esp),%xmm2
        movapd %xmm3,nb112_fjxH1(%esp)
        movapd %xmm4,nb112_fjyH1(%esp)
        movapd %xmm5,nb112_fjzH1(%esp)
        movapd %xmm0,nb112_fixO(%esp)
        movapd %xmm1,nb112_fiyO(%esp)
        movapd %xmm2,nb112_fizO(%esp)

        ## O-H2 interaction  
        movapd nb112_rinvOH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb112_qqOH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsOH2  
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb112_dxOH2(%esp),%xmm0
        mulpd nb112_dyOH2(%esp),%xmm1
        mulpd nb112_dzOH2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb112_fixO(%esp),%xmm0
        addpd nb112_fiyO(%esp),%xmm1
        addpd nb112_fizO(%esp),%xmm2
        movapd %xmm3,nb112_fjxH2(%esp)
        movapd %xmm4,nb112_fjyH2(%esp)
        movapd %xmm5,nb112_fjzH2(%esp)
        movapd %xmm0,nb112_fixO(%esp)
        movapd %xmm1,nb112_fiyO(%esp)
        movapd %xmm2,nb112_fizO(%esp)

        ## H1-O interaction 
        movapd nb112_rinvH1O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb112_qqOH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsH1O 
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb112_fjxO(%esp),%xmm3
        movapd nb112_fjyO(%esp),%xmm4
        movapd nb112_fjzO(%esp),%xmm5
        mulpd nb112_dxH1O(%esp),%xmm0
        mulpd nb112_dyH1O(%esp),%xmm1
        mulpd nb112_dzH1O(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb112_fixH1(%esp),%xmm0
        addpd nb112_fiyH1(%esp),%xmm1
        addpd nb112_fizH1(%esp),%xmm2
        movapd %xmm3,nb112_fjxO(%esp)
        movapd %xmm4,nb112_fjyO(%esp)
        movapd %xmm5,nb112_fjzO(%esp)
        movapd %xmm0,nb112_fixH1(%esp)
        movapd %xmm1,nb112_fiyH1(%esp)
        movapd %xmm2,nb112_fizH1(%esp)

        ## H1-H1 interaction 
        movapd nb112_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb112_qqHH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsH1H1 
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb112_fjxH1(%esp),%xmm3
        movapd nb112_fjyH1(%esp),%xmm4
        movapd nb112_fjzH1(%esp),%xmm5
        mulpd nb112_dxH1H1(%esp),%xmm0
        mulpd nb112_dyH1H1(%esp),%xmm1
        mulpd nb112_dzH1H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb112_fixH1(%esp),%xmm0
        addpd nb112_fiyH1(%esp),%xmm1
        addpd nb112_fizH1(%esp),%xmm2
        movapd %xmm3,nb112_fjxH1(%esp)
        movapd %xmm4,nb112_fjyH1(%esp)
        movapd %xmm5,nb112_fjzH1(%esp)
        movapd %xmm0,nb112_fixH1(%esp)
        movapd %xmm1,nb112_fiyH1(%esp)
        movapd %xmm2,nb112_fizH1(%esp)

        ## H1-H2 interaction 
        movapd nb112_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb112_qqHH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsOH2  
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb112_fjxH2(%esp),%xmm3
        movapd nb112_fjyH2(%esp),%xmm4
        movapd nb112_fjzH2(%esp),%xmm5
        mulpd nb112_dxH1H2(%esp),%xmm0
        mulpd nb112_dyH1H2(%esp),%xmm1
        mulpd nb112_dzH1H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb112_fixH1(%esp),%xmm0
        addpd nb112_fiyH1(%esp),%xmm1
        addpd nb112_fizH1(%esp),%xmm2
        movapd %xmm3,nb112_fjxH2(%esp)
        movapd %xmm4,nb112_fjyH2(%esp)
        movapd %xmm5,nb112_fjzH2(%esp)
        movapd %xmm0,nb112_fixH1(%esp)
        movapd %xmm1,nb112_fiyH1(%esp)
        movapd %xmm2,nb112_fizH1(%esp)

        ## H2-O interaction 
        movapd nb112_rinvH2O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb112_qqOH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsH2O 
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb112_fjxO(%esp),%xmm3
        movapd nb112_fjyO(%esp),%xmm4
        movapd nb112_fjzO(%esp),%xmm5
        mulpd nb112_dxH2O(%esp),%xmm0
        mulpd nb112_dyH2O(%esp),%xmm1
        mulpd nb112_dzH2O(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb112_fixH2(%esp),%xmm0
        addpd nb112_fiyH2(%esp),%xmm1
        addpd nb112_fizH2(%esp),%xmm2
        movapd %xmm3,nb112_fjxO(%esp)
        movapd %xmm4,nb112_fjyO(%esp)
        movapd %xmm5,nb112_fjzO(%esp)
        movapd %xmm0,nb112_fixH2(%esp)
        movapd %xmm1,nb112_fiyH2(%esp)
        movapd %xmm2,nb112_fizH2(%esp)

        ## H2-H1 interaction 
        movapd nb112_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb112_qqHH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsH2H1 
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb112_fjxH1(%esp),%xmm3
        movapd nb112_fjyH1(%esp),%xmm4
        movapd nb112_fjzH1(%esp),%xmm5
        mulpd nb112_dxH2H1(%esp),%xmm0
        mulpd nb112_dyH2H1(%esp),%xmm1
        mulpd nb112_dzH2H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb112_fixH2(%esp),%xmm0
        addpd nb112_fiyH2(%esp),%xmm1
        addpd nb112_fizH2(%esp),%xmm2
        movapd %xmm3,nb112_fjxH1(%esp)
        movapd %xmm4,nb112_fjyH1(%esp)
        movapd %xmm5,nb112_fjzH1(%esp)
        movapd %xmm0,nb112_fixH2(%esp)
        movapd %xmm1,nb112_fiyH2(%esp)
        movapd %xmm2,nb112_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb112_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb112_qqHH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsH2H2 
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm7,nb112_vctot(%esp)
        movapd %xmm0,%xmm2
        movapd nb112_fjxH2(%esp),%xmm3
        movapd nb112_fjyH2(%esp),%xmm4
        movapd nb112_fjzH2(%esp),%xmm5
        mulpd nb112_dxH2H2(%esp),%xmm0
        mulpd nb112_dyH2H2(%esp),%xmm1
        mulpd nb112_dzH2H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb112_fixH2(%esp),%xmm0
        addpd nb112_fiyH2(%esp),%xmm1
        addpd nb112_fizH2(%esp),%xmm2
        movapd %xmm3,nb112_fjxH2(%esp)
        movapd %xmm4,nb112_fjyH2(%esp)
        movapd %xmm5,nb112_fjzH2(%esp)
        movapd %xmm0,nb112_fixH2(%esp)
        movapd %xmm1,nb112_fiyH2(%esp)
        movapd %xmm2,nb112_fizH2(%esp)

        movl nb112_faction(%ebp),%edi

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
        addpd nb112_fjxO(%esp),%xmm0
        addpd nb112_fjyO(%esp),%xmm1
        addpd nb112_fjzO(%esp),%xmm2
        addpd nb112_fjxH1(%esp),%xmm3
        addpd nb112_fjyH1(%esp),%xmm4
        addpd nb112_fjzH1(%esp),%xmm5
        addpd nb112_fjxH2(%esp),%xmm6
        addpd nb112_fjyH2(%esp),%xmm7
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
        addpd nb112_fjzH2(%esp),%xmm0
        movlpd %xmm0,64(%edi,%eax,8)
        movhpd %xmm0,64(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl $2,nb112_innerk(%esp)
        jl   _nb_kernel112_ia32_sse2.nb112_checksingle
        jmp  _nb_kernel112_ia32_sse2.nb112_unroll_loop
_nb_kernel112_ia32_sse2.nb112_checksingle: 
        movl nb112_innerk(%esp),%edx
        andl $1,%edx
        jnz  _nb_kernel112_ia32_sse2.nb112_dosingle
        jmp  _nb_kernel112_ia32_sse2.nb112_updateouterdata
_nb_kernel112_ia32_sse2.nb112_dosingle: 
        movl  nb112_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb112_innerjjnr(%esp)

        movl nb112_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        movlpd (%esi,%eax,8),%xmm2
        movlpd 8(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%eax,8),%xmm4
        movlpd 24(%esi,%eax,8),%xmm5
        movlpd 32(%esi,%eax,8),%xmm6
        movlpd 40(%esi,%eax,8),%xmm7
        movapd  %xmm2,nb112_jxO(%esp)
        movapd  %xmm3,nb112_jyO(%esp)
        movapd  %xmm4,nb112_jzO(%esp)
        movapd  %xmm5,nb112_jxH1(%esp)
        movapd  %xmm6,nb112_jyH1(%esp)
        movapd  %xmm7,nb112_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movapd  %xmm2,nb112_jxH2(%esp)
        movapd  %xmm3,nb112_jyH2(%esp)
        movapd  %xmm4,nb112_jzH2(%esp)

        movapd nb112_ixO(%esp),%xmm0
        movapd nb112_iyO(%esp),%xmm1
        movapd nb112_izO(%esp),%xmm2
        movapd nb112_ixO(%esp),%xmm3
        movapd nb112_iyO(%esp),%xmm4
        movapd nb112_izO(%esp),%xmm5
        subsd  nb112_jxO(%esp),%xmm0
        subsd  nb112_jyO(%esp),%xmm1
        subsd  nb112_jzO(%esp),%xmm2
        subsd  nb112_jxH1(%esp),%xmm3
        subsd  nb112_jyH1(%esp),%xmm4
        subsd  nb112_jzH1(%esp),%xmm5
        movapd %xmm0,nb112_dxOO(%esp)
        movapd %xmm1,nb112_dyOO(%esp)
        movapd %xmm2,nb112_dzOO(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb112_dxOH1(%esp)
        movapd %xmm4,nb112_dyOH1(%esp)
        movapd %xmm5,nb112_dzOH1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb112_rsqOO(%esp)
        movapd %xmm3,nb112_rsqOH1(%esp)

        movapd nb112_ixO(%esp),%xmm0
        movapd nb112_iyO(%esp),%xmm1
        movapd nb112_izO(%esp),%xmm2
        movapd nb112_ixH1(%esp),%xmm3
        movapd nb112_iyH1(%esp),%xmm4
        movapd nb112_izH1(%esp),%xmm5
        subsd  nb112_jxH2(%esp),%xmm0
        subsd  nb112_jyH2(%esp),%xmm1
        subsd  nb112_jzH2(%esp),%xmm2
        subsd  nb112_jxO(%esp),%xmm3
        subsd  nb112_jyO(%esp),%xmm4
        subsd  nb112_jzO(%esp),%xmm5
        movapd %xmm0,nb112_dxOH2(%esp)
        movapd %xmm1,nb112_dyOH2(%esp)
        movapd %xmm2,nb112_dzOH2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb112_dxH1O(%esp)
        movapd %xmm4,nb112_dyH1O(%esp)
        movapd %xmm5,nb112_dzH1O(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb112_rsqOH2(%esp)
        movapd %xmm3,nb112_rsqH1O(%esp)

        movapd nb112_ixH1(%esp),%xmm0
        movapd nb112_iyH1(%esp),%xmm1
        movapd nb112_izH1(%esp),%xmm2
        movapd nb112_ixH1(%esp),%xmm3
        movapd nb112_iyH1(%esp),%xmm4
        movapd nb112_izH1(%esp),%xmm5
        subsd  nb112_jxH1(%esp),%xmm0
        subsd  nb112_jyH1(%esp),%xmm1
        subsd  nb112_jzH1(%esp),%xmm2
        subsd  nb112_jxH2(%esp),%xmm3
        subsd  nb112_jyH2(%esp),%xmm4
        subsd  nb112_jzH2(%esp),%xmm5
        movapd %xmm0,nb112_dxH1H1(%esp)
        movapd %xmm1,nb112_dyH1H1(%esp)
        movapd %xmm2,nb112_dzH1H1(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb112_dxH1H2(%esp)
        movapd %xmm4,nb112_dyH1H2(%esp)
        movapd %xmm5,nb112_dzH1H2(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb112_rsqH1H1(%esp)
        movapd %xmm3,nb112_rsqH1H2(%esp)

        movapd nb112_ixH2(%esp),%xmm0
        movapd nb112_iyH2(%esp),%xmm1
        movapd nb112_izH2(%esp),%xmm2
        movapd nb112_ixH2(%esp),%xmm3
        movapd nb112_iyH2(%esp),%xmm4
        movapd nb112_izH2(%esp),%xmm5
        subsd  nb112_jxO(%esp),%xmm0
        subsd  nb112_jyO(%esp),%xmm1
        subsd  nb112_jzO(%esp),%xmm2
        subsd  nb112_jxH1(%esp),%xmm3
        subsd  nb112_jyH1(%esp),%xmm4
        subsd  nb112_jzH1(%esp),%xmm5
        movapd %xmm0,nb112_dxH2O(%esp)
        movapd %xmm1,nb112_dyH2O(%esp)
        movapd %xmm2,nb112_dzH2O(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb112_dxH2H1(%esp)
        movapd %xmm4,nb112_dyH2H1(%esp)
        movapd %xmm5,nb112_dzH2H1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb112_rsqH2O(%esp)
        movapd %xmm4,nb112_rsqH2H1(%esp)

        movapd nb112_ixH2(%esp),%xmm0
        movapd nb112_iyH2(%esp),%xmm1
        movapd nb112_izH2(%esp),%xmm2
        subsd  nb112_jxH2(%esp),%xmm0
        subsd  nb112_jyH2(%esp),%xmm1
        subsd  nb112_jzH2(%esp),%xmm2
        movapd %xmm0,nb112_dxH2H2(%esp)
        movapd %xmm1,nb112_dyH2H2(%esp)
        movapd %xmm2,nb112_dzH2H2(%esp)
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb112_rsqH2H2(%esp)

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
        movapd  nb112_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112_half(%esp),%xmm3   ## iter1 
        mulsd   nb112_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112_half(%esp),%xmm1   ## rinv 
        mulsd   nb112_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb112_rinvH2H2(%esp)
        movapd %xmm5,nb112_rinvH2H1(%esp)

        movapd nb112_rsqOO(%esp),%xmm0
        movapd nb112_rsqOH1(%esp),%xmm4
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
        movapd  nb112_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb112_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112_half(%esp),%xmm1   ## rinv 
        mulsd   nb112_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb112_rinvOO(%esp)
        movapd %xmm5,nb112_rinvOH1(%esp)

        movapd nb112_rsqOH2(%esp),%xmm0
        movapd nb112_rsqH1O(%esp),%xmm4
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
        movapd  nb112_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112_half(%esp),%xmm3   ## iter1 
        mulsd   nb112_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112_half(%esp),%xmm1   ## rinv 
        mulsd   nb112_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb112_rinvOH2(%esp)
        movapd %xmm5,nb112_rinvH1O(%esp)

        movapd nb112_rsqH1H1(%esp),%xmm0
        movapd nb112_rsqH1H2(%esp),%xmm4
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
        movapd  nb112_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112_half(%esp),%xmm3   ## iter1a 
        mulsd   nb112_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112_half(%esp),%xmm1   ## rinv 
        mulsd   nb112_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb112_rinvH1H1(%esp)
        movapd %xmm5,nb112_rinvH1H2(%esp)

        movapd nb112_rsqH2O(%esp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb112_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb112_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb112_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb112_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb112_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb112_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm7
        mulsd  %xmm0,%xmm0
        movapd %xmm0,%xmm1
        mulsd  %xmm0,%xmm1
        mulsd  %xmm0,%xmm1      ## xmm1=rinvsix 
        mulsd  nb112_qqOO(%esp),%xmm7
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  nb112_c6(%esp),%xmm1
        mulsd  nb112_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addsd  nb112_Vvdwtot(%esp),%xmm3
        mulsd  nb112_six(%esp),%xmm1
        mulsd  nb112_twelve(%esp),%xmm2
        movlpd %xmm3,nb112_Vvdwtot(%esp)
        subsd  %xmm1,%xmm2
        addsd  %xmm7,%xmm2
        addsd  nb112_vctot(%esp),%xmm7
        mulsd  %xmm2,%xmm0

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb112_dxOO(%esp),%xmm0
        mulsd nb112_dyOO(%esp),%xmm1
        mulsd nb112_dzOO(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb112_fixO(%esp),%xmm0
        addsd nb112_fiyO(%esp),%xmm1
        addsd nb112_fizO(%esp),%xmm2
        movlpd %xmm3,nb112_fjxO(%esp)
        movlpd %xmm4,nb112_fjyO(%esp)
        movlpd %xmm5,nb112_fjzO(%esp)
        movlpd %xmm0,nb112_fixO(%esp)
        movlpd %xmm1,nb112_fiyO(%esp)
        movlpd %xmm2,nb112_fizO(%esp)

        ## O-H1 interaction 
        movapd nb112_rinvOH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb112_qqOH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsOH1  
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb112_dxOH1(%esp),%xmm0
        mulsd nb112_dyOH1(%esp),%xmm1
        mulsd nb112_dzOH1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb112_fixO(%esp),%xmm0
        addsd nb112_fiyO(%esp),%xmm1
        addsd nb112_fizO(%esp),%xmm2
        movlpd %xmm3,nb112_fjxH1(%esp)
        movlpd %xmm4,nb112_fjyH1(%esp)
        movlpd %xmm5,nb112_fjzH1(%esp)
        movlpd %xmm0,nb112_fixO(%esp)
        movlpd %xmm1,nb112_fiyO(%esp)
        movlpd %xmm2,nb112_fizO(%esp)

        ## O-H2 interaction  
        movapd nb112_rinvOH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb112_qqOH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsOH2  
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb112_dxOH2(%esp),%xmm0
        mulsd nb112_dyOH2(%esp),%xmm1
        mulsd nb112_dzOH2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb112_fixO(%esp),%xmm0
        addsd nb112_fiyO(%esp),%xmm1
        addsd nb112_fizO(%esp),%xmm2
        movlpd %xmm3,nb112_fjxH2(%esp)
        movlpd %xmm4,nb112_fjyH2(%esp)
        movlpd %xmm5,nb112_fjzH2(%esp)
        movlpd %xmm0,nb112_fixO(%esp)
        movlpd %xmm1,nb112_fiyO(%esp)
        movlpd %xmm2,nb112_fizO(%esp)

        ## H1-O interaction 
        movapd nb112_rinvH1O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb112_qqOH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsH1O 
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb112_fjxO(%esp),%xmm3
        movapd nb112_fjyO(%esp),%xmm4
        movapd nb112_fjzO(%esp),%xmm5
        mulsd nb112_dxH1O(%esp),%xmm0
        mulsd nb112_dyH1O(%esp),%xmm1
        mulsd nb112_dzH1O(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb112_fixH1(%esp),%xmm0
        addsd nb112_fiyH1(%esp),%xmm1
        addsd nb112_fizH1(%esp),%xmm2
        movlpd %xmm3,nb112_fjxO(%esp)
        movlpd %xmm4,nb112_fjyO(%esp)
        movlpd %xmm5,nb112_fjzO(%esp)
        movlpd %xmm0,nb112_fixH1(%esp)
        movlpd %xmm1,nb112_fiyH1(%esp)
        movlpd %xmm2,nb112_fizH1(%esp)

        ## H1-H1 interaction 
        movapd nb112_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb112_qqHH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsH1H1 
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb112_fjxH1(%esp),%xmm3
        movapd nb112_fjyH1(%esp),%xmm4
        movapd nb112_fjzH1(%esp),%xmm5
        mulsd nb112_dxH1H1(%esp),%xmm0
        mulsd nb112_dyH1H1(%esp),%xmm1
        mulsd nb112_dzH1H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb112_fixH1(%esp),%xmm0
        addsd nb112_fiyH1(%esp),%xmm1
        addsd nb112_fizH1(%esp),%xmm2
        movlpd %xmm3,nb112_fjxH1(%esp)
        movlpd %xmm4,nb112_fjyH1(%esp)
        movlpd %xmm5,nb112_fjzH1(%esp)
        movlpd %xmm0,nb112_fixH1(%esp)
        movlpd %xmm1,nb112_fiyH1(%esp)
        movlpd %xmm2,nb112_fizH1(%esp)

        ## H1-H2 interaction 
        movapd nb112_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb112_qqHH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsOH2  
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb112_fjxH2(%esp),%xmm3
        movapd nb112_fjyH2(%esp),%xmm4
        movapd nb112_fjzH2(%esp),%xmm5
        mulsd nb112_dxH1H2(%esp),%xmm0
        mulsd nb112_dyH1H2(%esp),%xmm1
        mulsd nb112_dzH1H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb112_fixH1(%esp),%xmm0
        addsd nb112_fiyH1(%esp),%xmm1
        addsd nb112_fizH1(%esp),%xmm2
        movlpd %xmm3,nb112_fjxH2(%esp)
        movlpd %xmm4,nb112_fjyH2(%esp)
        movlpd %xmm5,nb112_fjzH2(%esp)
        movlpd %xmm0,nb112_fixH1(%esp)
        movlpd %xmm1,nb112_fiyH1(%esp)
        movlpd %xmm2,nb112_fizH1(%esp)

        ## H2-O interaction 
        movapd nb112_rinvH2O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb112_qqOH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsH2O 
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb112_fjxO(%esp),%xmm3
        movapd nb112_fjyO(%esp),%xmm4
        movapd nb112_fjzO(%esp),%xmm5
        mulsd nb112_dxH2O(%esp),%xmm0
        mulsd nb112_dyH2O(%esp),%xmm1
        mulsd nb112_dzH2O(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb112_fixH2(%esp),%xmm0
        addsd nb112_fiyH2(%esp),%xmm1
        addsd nb112_fizH2(%esp),%xmm2
        movlpd %xmm3,nb112_fjxO(%esp)
        movlpd %xmm4,nb112_fjyO(%esp)
        movlpd %xmm5,nb112_fjzO(%esp)
        movlpd %xmm0,nb112_fixH2(%esp)
        movlpd %xmm1,nb112_fiyH2(%esp)
        movlpd %xmm2,nb112_fizH2(%esp)

        ## H2-H1 interaction 
        movapd nb112_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb112_qqHH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsH2H1 
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb112_fjxH1(%esp),%xmm3
        movapd nb112_fjyH1(%esp),%xmm4
        movapd nb112_fjzH1(%esp),%xmm5
        mulsd nb112_dxH2H1(%esp),%xmm0
        mulsd nb112_dyH2H1(%esp),%xmm1
        mulsd nb112_dzH2H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb112_fixH2(%esp),%xmm0
        addsd nb112_fiyH2(%esp),%xmm1
        addsd nb112_fizH2(%esp),%xmm2
        movlpd %xmm3,nb112_fjxH1(%esp)
        movlpd %xmm4,nb112_fjyH1(%esp)
        movlpd %xmm5,nb112_fjzH1(%esp)
        movlpd %xmm0,nb112_fixH2(%esp)
        movlpd %xmm1,nb112_fiyH2(%esp)
        movlpd %xmm2,nb112_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb112_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb112_qqHH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsH2H2 
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movlpd %xmm7,nb112_vctot(%esp)
        movapd %xmm0,%xmm2
        movlpd nb112_fjxH2(%esp),%xmm3
        movlpd nb112_fjyH2(%esp),%xmm4
        movlpd nb112_fjzH2(%esp),%xmm5
        mulsd nb112_dxH2H2(%esp),%xmm0
        mulsd nb112_dyH2H2(%esp),%xmm1
        mulsd nb112_dzH2H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb112_fixH2(%esp),%xmm0
        addsd nb112_fiyH2(%esp),%xmm1
        addsd nb112_fizH2(%esp),%xmm2
        movlpd %xmm3,nb112_fjxH2(%esp)
        movlpd %xmm4,nb112_fjyH2(%esp)
        movlpd %xmm5,nb112_fjzH2(%esp)
        movlpd %xmm0,nb112_fixH2(%esp)
        movlpd %xmm1,nb112_fiyH2(%esp)
        movlpd %xmm2,nb112_fizH2(%esp)

        movl nb112_faction(%ebp),%edi

        ## Did all interactions - now update j forces 
        movlpd (%edi,%eax,8),%xmm0
        movlpd 8(%edi,%eax,8),%xmm1
        movlpd 16(%edi,%eax,8),%xmm2
        movlpd 24(%edi,%eax,8),%xmm3
        movlpd 32(%edi,%eax,8),%xmm4
        movlpd 40(%edi,%eax,8),%xmm5
        movlpd 48(%edi,%eax,8),%xmm6
        movlpd 56(%edi,%eax,8),%xmm7
        addsd nb112_fjxO(%esp),%xmm0
        addsd nb112_fjyO(%esp),%xmm1
        addsd nb112_fjzO(%esp),%xmm2
        addsd nb112_fjxH1(%esp),%xmm3
        addsd nb112_fjyH1(%esp),%xmm4
        addsd nb112_fjzH1(%esp),%xmm5
        addsd nb112_fjxH2(%esp),%xmm6
        addsd nb112_fjyH2(%esp),%xmm7
        movlpd %xmm0,(%edi,%eax,8)
        movlpd %xmm1,8(%edi,%eax,8)
        movlpd %xmm2,16(%edi,%eax,8)
        movlpd %xmm3,24(%edi,%eax,8)
        movlpd %xmm4,32(%edi,%eax,8)
        movlpd %xmm5,40(%edi,%eax,8)
        movlpd %xmm6,48(%edi,%eax,8)
        movlpd %xmm7,56(%edi,%eax,8)

        movlpd 64(%edi,%eax,8),%xmm0
        addsd nb112_fjzH2(%esp),%xmm0
        movlpd %xmm0,64(%edi,%eax,8)
_nb_kernel112_ia32_sse2.nb112_updateouterdata: 
        movl  nb112_ii3(%esp),%ecx
        movl  nb112_faction(%ebp),%edi
        movl  nb112_fshift(%ebp),%esi
        movl  nb112_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb112_fixO(%esp),%xmm0
        movapd nb112_fiyO(%esp),%xmm1
        movapd nb112_fizO(%esp),%xmm2

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
        movapd nb112_fixH1(%esp),%xmm0
        movapd nb112_fiyH1(%esp),%xmm1
        movapd nb112_fizH1(%esp),%xmm2

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
        movapd nb112_fixH2(%esp),%xmm0
        movapd nb112_fiyH2(%esp),%xmm1
        movapd nb112_fizH2(%esp),%xmm2

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
        movl nb112_n(%esp),%esi
        ## get group index for i particle 
        movl  nb112_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb112_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb112_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb112_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb112_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

       ## finish if last 
        movl nb112_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel112_ia32_sse2.nb112_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb112_n(%esp)
        jmp _nb_kernel112_ia32_sse2.nb112_outer
_nb_kernel112_ia32_sse2.nb112_outerend: 
        ## check if more outer neighborlists remain
        movl  nb112_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel112_ia32_sse2.nb112_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel112_ia32_sse2.nb112_threadloop
_nb_kernel112_ia32_sse2.nb112_end: 
        emms

        movl nb112_nouter(%esp),%eax
        movl nb112_ninner(%esp),%ebx
        movl nb112_outeriter(%ebp),%ecx
        movl nb112_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb112_salign(%esp),%eax
        addl %eax,%esp
        addl $1512,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret





.globl nb_kernel112nf_ia32_sse2
.globl _nb_kernel112nf_ia32_sse2
nb_kernel112nf_ia32_sse2:       
_nb_kernel112nf_ia32_sse2:      
.set nb112nf_p_nri, 8
.set nb112nf_iinr, 12
.set nb112nf_jindex, 16
.set nb112nf_jjnr, 20
.set nb112nf_shift, 24
.set nb112nf_shiftvec, 28
.set nb112nf_fshift, 32
.set nb112nf_gid, 36
.set nb112nf_pos, 40
.set nb112nf_faction, 44
.set nb112nf_charge, 48
.set nb112nf_p_facel, 52
.set nb112nf_argkrf, 56
.set nb112nf_argcrf, 60
.set nb112nf_Vc, 64
.set nb112nf_type, 68
.set nb112nf_p_ntype, 72
.set nb112nf_vdwparam, 76
.set nb112nf_Vvdw, 80
.set nb112nf_p_tabscale, 84
.set nb112nf_VFtab, 88
.set nb112nf_invsqrta, 92
.set nb112nf_dvda, 96
.set nb112nf_p_gbtabscale, 100
.set nb112nf_GBtab, 104
.set nb112nf_p_nthreads, 108
.set nb112nf_count, 112
.set nb112nf_mtx, 116
.set nb112nf_outeriter, 120
.set nb112nf_inneriter, 124
.set nb112nf_work, 128
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
.set nb112nf_innerjjnr, 728
.set nb112nf_innerk, 732
.set nb112nf_n, 736
.set nb112nf_nn1, 740
.set nb112nf_nri, 744
.set nb112nf_nouter, 748
.set nb112nf_ninner, 752
.set nb112nf_salign, 756
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $760,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb112nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb112nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb112nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb112nf_nouter(%esp)
        movl %eax,nb112nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb112nf_half(%esp)
        movl %ebx,nb112nf_half+4(%esp)
        movsd nb112nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb112nf_half(%esp)
        movapd %xmm3,nb112nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb112nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb112nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb112nf_p_facel(%ebp),%esi
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
        movapd %xmm3,nb112nf_qqOO(%esp)
        movapd %xmm4,nb112nf_qqOH(%esp)
        movapd %xmm5,nb112nf_qqHH(%esp)

        xorpd %xmm0,%xmm0
        movl  nb112nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb112nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb112nf_vdwparam(%ebp),%eax
        movlpd (%eax,%edx,8),%xmm0
        movhpd 8(%eax,%edx,8),%xmm0
        movhlps %xmm0,%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb112nf_c6(%esp)
        movapd %xmm1,nb112nf_c12(%esp)

_nb_kernel112nf_ia32_sse2.nb112nf_threadloop: 
        movl  nb112nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel112nf_ia32_sse2.nb112nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel112nf_ia32_sse2.nb112nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb112nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb112nf_n(%esp)
        movl %ebx,nb112nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel112nf_ia32_sse2.nb112nf_outerstart
        jmp _nb_kernel112nf_ia32_sse2.nb112nf_end

_nb_kernel112nf_ia32_sse2.nb112nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb112nf_nouter(%esp),%ebx
        movl %ebx,nb112nf_nouter(%esp)

_nb_kernel112nf_ia32_sse2.nb112nf_outer: 
        movl  nb112nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb112nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb112nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb112nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb112nf_ii3(%esp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb112nf_ixO(%esp)
        movapd %xmm4,nb112nf_iyO(%esp)
        movapd %xmm5,nb112nf_izO(%esp)

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
        movapd %xmm0,nb112nf_ixH1(%esp)
        movapd %xmm1,nb112nf_iyH1(%esp)
        movapd %xmm2,nb112nf_izH1(%esp)
        movapd %xmm3,nb112nf_ixH2(%esp)
        movapd %xmm4,nb112nf_iyH2(%esp)
        movapd %xmm5,nb112nf_izH2(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb112nf_vctot(%esp)
        movapd %xmm4,nb112nf_Vvdwtot(%esp)

        movl  nb112nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb112nf_pos(%ebp),%esi
        movl  nb112nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb112nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb112nf_ninner(%esp),%ecx
        movl  %ecx,nb112nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb112nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel112nf_ia32_sse2.nb112nf_unroll_loop
        jmp   _nb_kernel112nf_ia32_sse2.nb112nf_checksingle
_nb_kernel112nf_ia32_sse2.nb112nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb112nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb112nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb112nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd  %xmm2,nb112nf_jxO(%esp)
        movapd  %xmm3,nb112nf_jyO(%esp)
        movapd  %xmm4,nb112nf_jzO(%esp)
        movapd  %xmm5,nb112nf_jxH1(%esp)
        movapd  %xmm6,nb112nf_jyH1(%esp)
        movapd  %xmm7,nb112nf_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 48(%esi,%ebx,8),%xmm2
        movhpd 56(%esi,%ebx,8),%xmm3
        movhpd 64(%esi,%ebx,8),%xmm4
        movapd  %xmm2,nb112nf_jxH2(%esp)
        movapd  %xmm3,nb112nf_jyH2(%esp)
        movapd  %xmm4,nb112nf_jzH2(%esp)

        movapd nb112nf_ixO(%esp),%xmm0
        movapd nb112nf_iyO(%esp),%xmm1
        movapd nb112nf_izO(%esp),%xmm2
        movapd nb112nf_ixO(%esp),%xmm3
        movapd nb112nf_iyO(%esp),%xmm4
        movapd nb112nf_izO(%esp),%xmm5
        subpd  nb112nf_jxO(%esp),%xmm0
        subpd  nb112nf_jyO(%esp),%xmm1
        subpd  nb112nf_jzO(%esp),%xmm2
        subpd  nb112nf_jxH1(%esp),%xmm3
        subpd  nb112nf_jyH1(%esp),%xmm4
        subpd  nb112nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb112nf_rsqOO(%esp)
        movapd %xmm3,nb112nf_rsqOH1(%esp)

        movapd nb112nf_ixO(%esp),%xmm0
        movapd nb112nf_iyO(%esp),%xmm1
        movapd nb112nf_izO(%esp),%xmm2
        movapd nb112nf_ixH1(%esp),%xmm3
        movapd nb112nf_iyH1(%esp),%xmm4
        movapd nb112nf_izH1(%esp),%xmm5
        subpd  nb112nf_jxH2(%esp),%xmm0
        subpd  nb112nf_jyH2(%esp),%xmm1
        subpd  nb112nf_jzH2(%esp),%xmm2
        subpd  nb112nf_jxO(%esp),%xmm3
        subpd  nb112nf_jyO(%esp),%xmm4
        subpd  nb112nf_jzO(%esp),%xmm5
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
        movapd %xmm0,nb112nf_rsqOH2(%esp)
        movapd %xmm3,nb112nf_rsqH1O(%esp)

        movapd nb112nf_ixH1(%esp),%xmm0
        movapd nb112nf_iyH1(%esp),%xmm1
        movapd nb112nf_izH1(%esp),%xmm2
        movapd nb112nf_ixH1(%esp),%xmm3
        movapd nb112nf_iyH1(%esp),%xmm4
        movapd nb112nf_izH1(%esp),%xmm5
        subpd  nb112nf_jxH1(%esp),%xmm0
        subpd  nb112nf_jyH1(%esp),%xmm1
        subpd  nb112nf_jzH1(%esp),%xmm2
        subpd  nb112nf_jxH2(%esp),%xmm3
        subpd  nb112nf_jyH2(%esp),%xmm4
        subpd  nb112nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb112nf_rsqH1H1(%esp)
        movapd %xmm3,nb112nf_rsqH1H2(%esp)

        movapd nb112nf_ixH2(%esp),%xmm0
        movapd nb112nf_iyH2(%esp),%xmm1
        movapd nb112nf_izH2(%esp),%xmm2
        movapd nb112nf_ixH2(%esp),%xmm3
        movapd nb112nf_iyH2(%esp),%xmm4
        movapd nb112nf_izH2(%esp),%xmm5
        subpd  nb112nf_jxO(%esp),%xmm0
        subpd  nb112nf_jyO(%esp),%xmm1
        subpd  nb112nf_jzO(%esp),%xmm2
        subpd  nb112nf_jxH1(%esp),%xmm3
        subpd  nb112nf_jyH1(%esp),%xmm4
        subpd  nb112nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb112nf_rsqH2O(%esp)
        movapd %xmm4,nb112nf_rsqH2H1(%esp)

        movapd nb112nf_ixH2(%esp),%xmm0
        movapd nb112nf_iyH2(%esp),%xmm1
        movapd nb112nf_izH2(%esp),%xmm2
        subpd  nb112nf_jxH2(%esp),%xmm0
        subpd  nb112nf_jyH2(%esp),%xmm1
        subpd  nb112nf_jzH2(%esp),%xmm2
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb112nf_rsqH2H2(%esp)

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
        movapd  nb112nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb112nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb112nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb112nf_rinvH2H2(%esp)
        movapd %xmm5,nb112nf_rinvH2H1(%esp)

        movapd nb112nf_rsqOO(%esp),%xmm0
        movapd nb112nf_rsqOH1(%esp),%xmm4
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
        movapd  nb112nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112nf_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb112nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb112nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb112nf_rinvOO(%esp)
        movapd %xmm5,nb112nf_rinvOH1(%esp)

        movapd nb112nf_rsqOH2(%esp),%xmm0
        movapd nb112nf_rsqH1O(%esp),%xmm4
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
        movapd  nb112nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb112nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb112nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb112nf_rinvOH2(%esp)
        movapd %xmm5,nb112nf_rinvH1O(%esp)

        movapd nb112nf_rsqH1H1(%esp),%xmm0
        movapd nb112nf_rsqH1H2(%esp),%xmm4
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
        movapd  nb112nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb112nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb112nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb112nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb112nf_rinvH1H1(%esp)
        movapd %xmm5,nb112nf_rinvH1H2(%esp)

        movapd nb112nf_rsqH2O(%esp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb112nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb112nf_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb112nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb112nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb112nf_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb112nf_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm7
        mulpd  %xmm0,%xmm0
        movapd %xmm0,%xmm1
        mulpd  %xmm0,%xmm1
        mulpd  %xmm0,%xmm1      ## xmm1=rinvsix 
        mulpd  nb112nf_qqOO(%esp),%xmm7
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulpd  nb112nf_c6(%esp),%xmm1
        mulpd  nb112nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addpd  nb112nf_Vvdwtot(%esp),%xmm3
        movapd %xmm3,nb112nf_Vvdwtot(%esp)
        addpd  nb112nf_vctot(%esp),%xmm7

        ## other interactions 
        movapd nb112nf_rinvOH1(%esp),%xmm1
        movapd nb112nf_rinvH1H1(%esp),%xmm2

        addpd nb112nf_rinvOH2(%esp),%xmm1
        addpd nb112nf_rinvH1H2(%esp),%xmm2

        addpd nb112nf_rinvH1O(%esp),%xmm1
        addpd nb112nf_rinvH2H1(%esp),%xmm2

        addpd nb112nf_rinvH2O(%esp),%xmm1
        addpd nb112nf_rinvH2H2(%esp),%xmm2

        mulpd nb112nf_qqOH(%esp),%xmm1
        mulpd nb112nf_qqHH(%esp),%xmm2

        addpd %xmm1,%xmm7
        addpd %xmm2,%xmm7

        movapd %xmm7,nb112nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb112nf_innerk(%esp)
        jl    _nb_kernel112nf_ia32_sse2.nb112nf_checksingle
        jmp   _nb_kernel112nf_ia32_sse2.nb112nf_unroll_loop
_nb_kernel112nf_ia32_sse2.nb112nf_checksingle: 
        movl  nb112nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel112nf_ia32_sse2.nb112nf_dosingle
        jmp   _nb_kernel112nf_ia32_sse2.nb112nf_updateouterdata
_nb_kernel112nf_ia32_sse2.nb112nf_dosingle: 
        movl  nb112nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb112nf_innerjjnr(%esp)

        movl nb112nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        movlpd (%esi,%eax,8),%xmm2
        movlpd 8(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%eax,8),%xmm4
        movlpd 24(%esi,%eax,8),%xmm5
        movlpd 32(%esi,%eax,8),%xmm6
        movlpd 40(%esi,%eax,8),%xmm7
        movapd  %xmm2,nb112nf_jxO(%esp)
        movapd  %xmm3,nb112nf_jyO(%esp)
        movapd  %xmm4,nb112nf_jzO(%esp)
        movapd  %xmm5,nb112nf_jxH1(%esp)
        movapd  %xmm6,nb112nf_jyH1(%esp)
        movapd  %xmm7,nb112nf_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movapd  %xmm2,nb112nf_jxH2(%esp)
        movapd  %xmm3,nb112nf_jyH2(%esp)
        movapd  %xmm4,nb112nf_jzH2(%esp)

        movapd nb112nf_ixO(%esp),%xmm0
        movapd nb112nf_iyO(%esp),%xmm1
        movapd nb112nf_izO(%esp),%xmm2
        movapd nb112nf_ixO(%esp),%xmm3
        movapd nb112nf_iyO(%esp),%xmm4
        movapd nb112nf_izO(%esp),%xmm5
        subsd  nb112nf_jxO(%esp),%xmm0
        subsd  nb112nf_jyO(%esp),%xmm1
        subsd  nb112nf_jzO(%esp),%xmm2
        subsd  nb112nf_jxH1(%esp),%xmm3
        subsd  nb112nf_jyH1(%esp),%xmm4
        subsd  nb112nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb112nf_rsqOO(%esp)
        movapd %xmm3,nb112nf_rsqOH1(%esp)

        movapd nb112nf_ixO(%esp),%xmm0
        movapd nb112nf_iyO(%esp),%xmm1
        movapd nb112nf_izO(%esp),%xmm2
        movapd nb112nf_ixH1(%esp),%xmm3
        movapd nb112nf_iyH1(%esp),%xmm4
        movapd nb112nf_izH1(%esp),%xmm5
        subsd  nb112nf_jxH2(%esp),%xmm0
        subsd  nb112nf_jyH2(%esp),%xmm1
        subsd  nb112nf_jzH2(%esp),%xmm2
        subsd  nb112nf_jxO(%esp),%xmm3
        subsd  nb112nf_jyO(%esp),%xmm4
        subsd  nb112nf_jzO(%esp),%xmm5
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
        movapd %xmm0,nb112nf_rsqOH2(%esp)
        movapd %xmm3,nb112nf_rsqH1O(%esp)

        movapd nb112nf_ixH1(%esp),%xmm0
        movapd nb112nf_iyH1(%esp),%xmm1
        movapd nb112nf_izH1(%esp),%xmm2
        movapd nb112nf_ixH1(%esp),%xmm3
        movapd nb112nf_iyH1(%esp),%xmm4
        movapd nb112nf_izH1(%esp),%xmm5
        subsd  nb112nf_jxH1(%esp),%xmm0
        subsd  nb112nf_jyH1(%esp),%xmm1
        subsd  nb112nf_jzH1(%esp),%xmm2
        subsd  nb112nf_jxH2(%esp),%xmm3
        subsd  nb112nf_jyH2(%esp),%xmm4
        subsd  nb112nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb112nf_rsqH1H1(%esp)
        movapd %xmm3,nb112nf_rsqH1H2(%esp)

        movapd nb112nf_ixH2(%esp),%xmm0
        movapd nb112nf_iyH2(%esp),%xmm1
        movapd nb112nf_izH2(%esp),%xmm2
        movapd nb112nf_ixH2(%esp),%xmm3
        movapd nb112nf_iyH2(%esp),%xmm4
        movapd nb112nf_izH2(%esp),%xmm5
        subsd  nb112nf_jxO(%esp),%xmm0
        subsd  nb112nf_jyO(%esp),%xmm1
        subsd  nb112nf_jzO(%esp),%xmm2
        subsd  nb112nf_jxH1(%esp),%xmm3
        subsd  nb112nf_jyH1(%esp),%xmm4
        subsd  nb112nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb112nf_rsqH2O(%esp)
        movapd %xmm4,nb112nf_rsqH2H1(%esp)

        movapd nb112nf_ixH2(%esp),%xmm0
        movapd nb112nf_iyH2(%esp),%xmm1
        movapd nb112nf_izH2(%esp),%xmm2
        subsd  nb112nf_jxH2(%esp),%xmm0
        subsd  nb112nf_jyH2(%esp),%xmm1
        subsd  nb112nf_jzH2(%esp),%xmm2
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb112nf_rsqH2H2(%esp)

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
        movapd  nb112nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb112nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb112nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb112nf_rinvH2H2(%esp)
        movapd %xmm5,nb112nf_rinvH2H1(%esp)

        movapd nb112nf_rsqOO(%esp),%xmm0
        movapd nb112nf_rsqOH1(%esp),%xmm4
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
        movapd  nb112nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112nf_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb112nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb112nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb112nf_rinvOO(%esp)
        movapd %xmm5,nb112nf_rinvOH1(%esp)

        movapd nb112nf_rsqOH2(%esp),%xmm0
        movapd nb112nf_rsqH1O(%esp),%xmm4
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
        movapd  nb112nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb112nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb112nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb112nf_rinvOH2(%esp)
        movapd %xmm5,nb112nf_rinvH1O(%esp)

        movapd nb112nf_rsqH1H1(%esp),%xmm0
        movapd nb112nf_rsqH1H2(%esp),%xmm4
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
        movapd  nb112nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb112nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb112nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb112nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb112nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb112nf_rinvH1H1(%esp)
        movapd %xmm5,nb112nf_rinvH1H2(%esp)

        movapd nb112nf_rsqH2O(%esp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb112nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb112nf_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb112nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb112nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb112nf_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb112nf_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm7
        mulsd  %xmm0,%xmm0
        movapd %xmm0,%xmm1
        mulsd  %xmm0,%xmm1
        mulsd  %xmm0,%xmm1      ## xmm1=rinvsix 
        mulsd  nb112nf_qqOO(%esp),%xmm7
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  nb112nf_c6(%esp),%xmm1
        mulsd  nb112nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addsd  nb112nf_Vvdwtot(%esp),%xmm3
        movlpd %xmm3,nb112nf_Vvdwtot(%esp)
        addsd  nb112nf_vctot(%esp),%xmm7

        ## other interactions 
        movapd nb112nf_rinvOH1(%esp),%xmm1
        movapd nb112nf_rinvH1H1(%esp),%xmm2

        addsd nb112nf_rinvOH2(%esp),%xmm1
        addsd nb112nf_rinvH1H2(%esp),%xmm2

        addsd nb112nf_rinvH1O(%esp),%xmm1
        addsd nb112nf_rinvH2H1(%esp),%xmm2

        addsd nb112nf_rinvH2O(%esp),%xmm1
        addsd nb112nf_rinvH2H2(%esp),%xmm2

        mulsd nb112nf_qqOH(%esp),%xmm1
        mulsd nb112nf_qqHH(%esp),%xmm2

        addsd %xmm1,%xmm7
        addsd %xmm2,%xmm7

        movlpd %xmm7,nb112nf_vctot(%esp)

_nb_kernel112nf_ia32_sse2.nb112nf_updateouterdata: 
        ## get n from stack
        movl nb112nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb112nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb112nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb112nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb112nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb112nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb112nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel112nf_ia32_sse2.nb112nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb112nf_n(%esp)
        jmp _nb_kernel112nf_ia32_sse2.nb112nf_outer
_nb_kernel112nf_ia32_sse2.nb112nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb112nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel112nf_ia32_sse2.nb112nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel112nf_ia32_sse2.nb112nf_threadloop
_nb_kernel112nf_ia32_sse2.nb112nf_end: 
        emms

        movl nb112nf_nouter(%esp),%eax
        movl nb112nf_ninner(%esp),%ebx
        movl nb112nf_outeriter(%ebp),%ecx
        movl nb112nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb112nf_salign(%esp),%eax
        addl %eax,%esp
        addl $760,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret

