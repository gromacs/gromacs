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


.globl nb_kernel102_ia32_sse2
.globl _nb_kernel102_ia32_sse2
nb_kernel102_ia32_sse2: 
_nb_kernel102_ia32_sse2:        
.set nb102_p_nri, 8
.set nb102_iinr, 12
.set nb102_jindex, 16
.set nb102_jjnr, 20
.set nb102_shift, 24
.set nb102_shiftvec, 28
.set nb102_fshift, 32
.set nb102_gid, 36
.set nb102_pos, 40
.set nb102_faction, 44
.set nb102_charge, 48
.set nb102_p_facel, 52
.set nb102_argkrf, 56
.set nb102_argcrf, 60
.set nb102_Vc, 64
.set nb102_type, 68
.set nb102_p_ntype, 72
.set nb102_vdwparam, 76
.set nb102_Vvdw, 80
.set nb102_p_tabscale, 84
.set nb102_VFtab, 88
.set nb102_invsqrta, 92
.set nb102_dvda, 96
.set nb102_p_gbtabscale, 100
.set nb102_GBtab, 104
.set nb102_p_nthreads, 108
.set nb102_count, 112
.set nb102_mtx, 116
.set nb102_outeriter, 120
.set nb102_inneriter, 124
.set nb102_work, 128
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
.set nb102_innerk, 1404
.set nb102_n, 1408
.set nb102_nn1, 1412
.set nb102_nri, 1416
.set nb102_nouter, 1420
.set nb102_ninner, 1424
.set nb102_salign, 1428
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $1432,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb102_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb102_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb102_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb102_nouter(%esp)
        movl %eax,nb102_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb102_half(%esp)
        movl %ebx,nb102_half+4(%esp)
        movsd nb102_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb102_half(%esp)
        movapd %xmm3,nb102_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb102_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb102_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3       ## qO 
        movsd %xmm3,%xmm4               ## qO 
        movsd 8(%edx,%ebx,8),%xmm5      ## qH 
        movl nb102_p_facel(%ebp),%esi
        movsd (%esi),%xmm6      ## facel 
        mulsd  %xmm3,%xmm3              ## qO*qO 
        mulsd  %xmm5,%xmm4              ## qO*qH 
        mulsd  %xmm5,%xmm5              ## qH*qH 
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb102_qqOO(%esp)
        movapd %xmm4,nb102_qqOH(%esp)
        movapd %xmm5,nb102_qqHH(%esp)

_nb_kernel102_ia32_sse2.nb102_threadloop: 
        movl  nb102_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel102_ia32_sse2.nb102_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel102_ia32_sse2.nb102_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb102_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb102_n(%esp)
        movl %ebx,nb102_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel102_ia32_sse2.nb102_outerstart
        jmp _nb_kernel102_ia32_sse2.nb102_end

_nb_kernel102_ia32_sse2.nb102_outerstart: 
        ## ebx contains number of outer iterations
        addl nb102_nouter(%esp),%ebx
        movl %ebx,nb102_nouter(%esp)

_nb_kernel102_ia32_sse2.nb102_outer: 
        movl  nb102_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb102_is3(%esp)      ## store is3 

        movl  nb102_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb102_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb102_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb102_ii3(%esp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb102_ixO(%esp)
        movapd %xmm4,nb102_iyO(%esp)
        movapd %xmm5,nb102_izO(%esp)

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
        movapd %xmm0,nb102_ixH1(%esp)
        movapd %xmm1,nb102_iyH1(%esp)
        movapd %xmm2,nb102_izH1(%esp)
        movapd %xmm3,nb102_ixH2(%esp)
        movapd %xmm4,nb102_iyH2(%esp)
        movapd %xmm5,nb102_izH2(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb102_vctot(%esp)
        movapd %xmm4,nb102_fixO(%esp)
        movapd %xmm4,nb102_fiyO(%esp)
        movapd %xmm4,nb102_fizO(%esp)
        movapd %xmm4,nb102_fixH1(%esp)
        movapd %xmm4,nb102_fiyH1(%esp)
        movapd %xmm4,nb102_fizH1(%esp)
        movapd %xmm4,nb102_fixH2(%esp)
        movapd %xmm4,nb102_fiyH2(%esp)
        movapd %xmm4,nb102_fizH2(%esp)

        movl  nb102_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb102_pos(%ebp),%esi
        movl  nb102_faction(%ebp),%edi
        movl  nb102_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb102_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb102_ninner(%esp),%ecx
        movl  %ecx,nb102_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb102_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel102_ia32_sse2.nb102_unroll_loop
        jmp   _nb_kernel102_ia32_sse2.nb102_checksingle
_nb_kernel102_ia32_sse2.nb102_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb102_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb102_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb102_pos(%ebp),%esi        ## base of pos[] 

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
        movapd  %xmm2,nb102_jxO(%esp)
        movapd  %xmm3,nb102_jyO(%esp)
        movapd  %xmm4,nb102_jzO(%esp)
        movapd  %xmm5,nb102_jxH1(%esp)
        movapd  %xmm6,nb102_jyH1(%esp)
        movapd  %xmm7,nb102_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 48(%esi,%ebx,8),%xmm2
        movhpd 56(%esi,%ebx,8),%xmm3
        movhpd 64(%esi,%ebx,8),%xmm4
        movapd  %xmm2,nb102_jxH2(%esp)
        movapd  %xmm3,nb102_jyH2(%esp)
        movapd  %xmm4,nb102_jzH2(%esp)

        movapd nb102_ixO(%esp),%xmm0
        movapd nb102_iyO(%esp),%xmm1
        movapd nb102_izO(%esp),%xmm2
        movapd nb102_ixO(%esp),%xmm3
        movapd nb102_iyO(%esp),%xmm4
        movapd nb102_izO(%esp),%xmm5
        subpd  nb102_jxO(%esp),%xmm0
        subpd  nb102_jyO(%esp),%xmm1
        subpd  nb102_jzO(%esp),%xmm2
        subpd  nb102_jxH1(%esp),%xmm3
        subpd  nb102_jyH1(%esp),%xmm4
        subpd  nb102_jzH1(%esp),%xmm5
        movapd %xmm0,nb102_dxOO(%esp)
        movapd %xmm1,nb102_dyOO(%esp)
        movapd %xmm2,nb102_dzOO(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb102_dxOH1(%esp)
        movapd %xmm4,nb102_dyOH1(%esp)
        movapd %xmm5,nb102_dzOH1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb102_rsqOO(%esp)
        movapd %xmm3,nb102_rsqOH1(%esp)

        movapd nb102_ixO(%esp),%xmm0
        movapd nb102_iyO(%esp),%xmm1
        movapd nb102_izO(%esp),%xmm2
        movapd nb102_ixH1(%esp),%xmm3
        movapd nb102_iyH1(%esp),%xmm4
        movapd nb102_izH1(%esp),%xmm5
        subpd  nb102_jxH2(%esp),%xmm0
        subpd  nb102_jyH2(%esp),%xmm1
        subpd  nb102_jzH2(%esp),%xmm2
        subpd  nb102_jxO(%esp),%xmm3
        subpd  nb102_jyO(%esp),%xmm4
        subpd  nb102_jzO(%esp),%xmm5
        movapd %xmm0,nb102_dxOH2(%esp)
        movapd %xmm1,nb102_dyOH2(%esp)
        movapd %xmm2,nb102_dzOH2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb102_dxH1O(%esp)
        movapd %xmm4,nb102_dyH1O(%esp)
        movapd %xmm5,nb102_dzH1O(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb102_rsqOH2(%esp)
        movapd %xmm3,nb102_rsqH1O(%esp)

        movapd nb102_ixH1(%esp),%xmm0
        movapd nb102_iyH1(%esp),%xmm1
        movapd nb102_izH1(%esp),%xmm2
        movapd nb102_ixH1(%esp),%xmm3
        movapd nb102_iyH1(%esp),%xmm4
        movapd nb102_izH1(%esp),%xmm5
        subpd  nb102_jxH1(%esp),%xmm0
        subpd  nb102_jyH1(%esp),%xmm1
        subpd  nb102_jzH1(%esp),%xmm2
        subpd  nb102_jxH2(%esp),%xmm3
        subpd  nb102_jyH2(%esp),%xmm4
        subpd  nb102_jzH2(%esp),%xmm5
        movapd %xmm0,nb102_dxH1H1(%esp)
        movapd %xmm1,nb102_dyH1H1(%esp)
        movapd %xmm2,nb102_dzH1H1(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb102_dxH1H2(%esp)
        movapd %xmm4,nb102_dyH1H2(%esp)
        movapd %xmm5,nb102_dzH1H2(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb102_rsqH1H1(%esp)
        movapd %xmm3,nb102_rsqH1H2(%esp)

        movapd nb102_ixH2(%esp),%xmm0
        movapd nb102_iyH2(%esp),%xmm1
        movapd nb102_izH2(%esp),%xmm2
        movapd nb102_ixH2(%esp),%xmm3
        movapd nb102_iyH2(%esp),%xmm4
        movapd nb102_izH2(%esp),%xmm5
        subpd  nb102_jxO(%esp),%xmm0
        subpd  nb102_jyO(%esp),%xmm1
        subpd  nb102_jzO(%esp),%xmm2
        subpd  nb102_jxH1(%esp),%xmm3
        subpd  nb102_jyH1(%esp),%xmm4
        subpd  nb102_jzH1(%esp),%xmm5
        movapd %xmm0,nb102_dxH2O(%esp)
        movapd %xmm1,nb102_dyH2O(%esp)
        movapd %xmm2,nb102_dzH2O(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb102_dxH2H1(%esp)
        movapd %xmm4,nb102_dyH2H1(%esp)
        movapd %xmm5,nb102_dzH2H1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb102_rsqH2O(%esp)
        movapd %xmm4,nb102_rsqH2H1(%esp)

        movapd nb102_ixH2(%esp),%xmm0
        movapd nb102_iyH2(%esp),%xmm1
        movapd nb102_izH2(%esp),%xmm2
        subpd  nb102_jxH2(%esp),%xmm0
        subpd  nb102_jyH2(%esp),%xmm1
        subpd  nb102_jzH2(%esp),%xmm2
        movapd %xmm0,nb102_dxH2H2(%esp)
        movapd %xmm1,nb102_dyH2H2(%esp)
        movapd %xmm2,nb102_dzH2H2(%esp)
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb102_rsqH2H2(%esp)

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
        movapd  nb102_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102_half(%esp),%xmm3   ## iter1 
        mulpd   nb102_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102_half(%esp),%xmm1   ## rinv 
        mulpd   nb102_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb102_rinvH2H2(%esp)
        movapd %xmm5,nb102_rinvH2H1(%esp)

        movapd nb102_rsqOO(%esp),%xmm0
        movapd nb102_rsqOH1(%esp),%xmm4
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
        movapd  nb102_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb102_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102_half(%esp),%xmm1   ## rinv 
        mulpd   nb102_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb102_rinvOO(%esp)
        movapd %xmm5,nb102_rinvOH1(%esp)

        movapd nb102_rsqOH2(%esp),%xmm0
        movapd nb102_rsqH1O(%esp),%xmm4
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
        movapd  nb102_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102_half(%esp),%xmm3   ## iter1 
        mulpd   nb102_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102_half(%esp),%xmm1   ## rinv 
        mulpd   nb102_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb102_rinvOH2(%esp)
        movapd %xmm5,nb102_rinvH1O(%esp)

        movapd nb102_rsqH1H1(%esp),%xmm0
        movapd nb102_rsqH1H2(%esp),%xmm4
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
        movapd  nb102_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102_half(%esp),%xmm3   ## iter1a 
        mulpd   nb102_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102_half(%esp),%xmm1   ## rinv 
        mulpd   nb102_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb102_rinvH1H1(%esp)
        movapd %xmm5,nb102_rinvH1H2(%esp)

        movapd nb102_rsqH2O(%esp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb102_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb102_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb102_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb102_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb102_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb102_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm7
        mulpd  %xmm0,%xmm0              ## rinvsq 
        mulpd  nb102_qqOO(%esp),%xmm7
        mulpd  %xmm7,%xmm0
        addpd  nb102_vctot(%esp),%xmm7
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb102_dxOO(%esp),%xmm0
        mulpd nb102_dyOO(%esp),%xmm1
        mulpd nb102_dzOO(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb102_fixO(%esp),%xmm0
        addpd nb102_fiyO(%esp),%xmm1
        addpd nb102_fizO(%esp),%xmm2
        movapd %xmm3,nb102_fjxO(%esp)
        movapd %xmm4,nb102_fjyO(%esp)
        movapd %xmm5,nb102_fjzO(%esp)
        movapd %xmm0,nb102_fixO(%esp)
        movapd %xmm1,nb102_fiyO(%esp)
        movapd %xmm2,nb102_fizO(%esp)

        ## O-H1 interaction 
        movapd nb102_rinvOH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb102_qqOH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsOH1  
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb102_dxOH1(%esp),%xmm0
        mulpd nb102_dyOH1(%esp),%xmm1
        mulpd nb102_dzOH1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb102_fixO(%esp),%xmm0
        addpd nb102_fiyO(%esp),%xmm1
        addpd nb102_fizO(%esp),%xmm2
        movapd %xmm3,nb102_fjxH1(%esp)
        movapd %xmm4,nb102_fjyH1(%esp)
        movapd %xmm5,nb102_fjzH1(%esp)
        movapd %xmm0,nb102_fixO(%esp)
        movapd %xmm1,nb102_fiyO(%esp)
        movapd %xmm2,nb102_fizO(%esp)

        ## O-H2 interaction  
        movapd nb102_rinvOH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb102_qqOH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsOH2  
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb102_dxOH2(%esp),%xmm0
        mulpd nb102_dyOH2(%esp),%xmm1
        mulpd nb102_dzOH2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb102_fixO(%esp),%xmm0
        addpd nb102_fiyO(%esp),%xmm1
        addpd nb102_fizO(%esp),%xmm2
        movapd %xmm3,nb102_fjxH2(%esp)
        movapd %xmm4,nb102_fjyH2(%esp)
        movapd %xmm5,nb102_fjzH2(%esp)
        movapd %xmm0,nb102_fixO(%esp)
        movapd %xmm1,nb102_fiyO(%esp)
        movapd %xmm2,nb102_fizO(%esp)

        ## H1-O interaction 
        movapd nb102_rinvH1O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb102_qqOH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsH1O 
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb102_fjxO(%esp),%xmm3
        movapd nb102_fjyO(%esp),%xmm4
        movapd nb102_fjzO(%esp),%xmm5
        mulpd nb102_dxH1O(%esp),%xmm0
        mulpd nb102_dyH1O(%esp),%xmm1
        mulpd nb102_dzH1O(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb102_fixH1(%esp),%xmm0
        addpd nb102_fiyH1(%esp),%xmm1
        addpd nb102_fizH1(%esp),%xmm2
        movapd %xmm3,nb102_fjxO(%esp)
        movapd %xmm4,nb102_fjyO(%esp)
        movapd %xmm5,nb102_fjzO(%esp)
        movapd %xmm0,nb102_fixH1(%esp)
        movapd %xmm1,nb102_fiyH1(%esp)
        movapd %xmm2,nb102_fizH1(%esp)

        ## H1-H1 interaction 
        movapd nb102_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb102_qqHH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsH1H1 
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb102_fjxH1(%esp),%xmm3
        movapd nb102_fjyH1(%esp),%xmm4
        movapd nb102_fjzH1(%esp),%xmm5
        mulpd nb102_dxH1H1(%esp),%xmm0
        mulpd nb102_dyH1H1(%esp),%xmm1
        mulpd nb102_dzH1H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb102_fixH1(%esp),%xmm0
        addpd nb102_fiyH1(%esp),%xmm1
        addpd nb102_fizH1(%esp),%xmm2
        movapd %xmm3,nb102_fjxH1(%esp)
        movapd %xmm4,nb102_fjyH1(%esp)
        movapd %xmm5,nb102_fjzH1(%esp)
        movapd %xmm0,nb102_fixH1(%esp)
        movapd %xmm1,nb102_fiyH1(%esp)
        movapd %xmm2,nb102_fizH1(%esp)

        ## H1-H2 interaction 
        movapd nb102_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb102_qqHH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsOH2  
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb102_fjxH2(%esp),%xmm3
        movapd nb102_fjyH2(%esp),%xmm4
        movapd nb102_fjzH2(%esp),%xmm5
        mulpd nb102_dxH1H2(%esp),%xmm0
        mulpd nb102_dyH1H2(%esp),%xmm1
        mulpd nb102_dzH1H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb102_fixH1(%esp),%xmm0
        addpd nb102_fiyH1(%esp),%xmm1
        addpd nb102_fizH1(%esp),%xmm2
        movapd %xmm3,nb102_fjxH2(%esp)
        movapd %xmm4,nb102_fjyH2(%esp)
        movapd %xmm5,nb102_fjzH2(%esp)
        movapd %xmm0,nb102_fixH1(%esp)
        movapd %xmm1,nb102_fiyH1(%esp)
        movapd %xmm2,nb102_fizH1(%esp)

        ## H2-O interaction 
        movapd nb102_rinvH2O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb102_qqOH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsH2O 
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb102_fjxO(%esp),%xmm3
        movapd nb102_fjyO(%esp),%xmm4
        movapd nb102_fjzO(%esp),%xmm5
        mulpd nb102_dxH2O(%esp),%xmm0
        mulpd nb102_dyH2O(%esp),%xmm1
        mulpd nb102_dzH2O(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb102_fixH2(%esp),%xmm0
        addpd nb102_fiyH2(%esp),%xmm1
        addpd nb102_fizH2(%esp),%xmm2
        movapd %xmm3,nb102_fjxO(%esp)
        movapd %xmm4,nb102_fjyO(%esp)
        movapd %xmm5,nb102_fjzO(%esp)
        movapd %xmm0,nb102_fixH2(%esp)
        movapd %xmm1,nb102_fiyH2(%esp)
        movapd %xmm2,nb102_fizH2(%esp)

        ## H2-H1 interaction 
        movapd nb102_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb102_qqHH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsH2H1 
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb102_fjxH1(%esp),%xmm3
        movapd nb102_fjyH1(%esp),%xmm4
        movapd nb102_fjzH1(%esp),%xmm5
        mulpd nb102_dxH2H1(%esp),%xmm0
        mulpd nb102_dyH2H1(%esp),%xmm1
        mulpd nb102_dzH2H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb102_fixH2(%esp),%xmm0
        addpd nb102_fiyH2(%esp),%xmm1
        addpd nb102_fizH2(%esp),%xmm2
        movapd %xmm3,nb102_fjxH1(%esp)
        movapd %xmm4,nb102_fjyH1(%esp)
        movapd %xmm5,nb102_fjzH1(%esp)
        movapd %xmm0,nb102_fixH2(%esp)
        movapd %xmm1,nb102_fiyH2(%esp)
        movapd %xmm2,nb102_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb102_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb102_qqHH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsH2H2 
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm7,nb102_vctot(%esp)
        movapd %xmm0,%xmm2
        movapd nb102_fjxH2(%esp),%xmm3
        movapd nb102_fjyH2(%esp),%xmm4
        movapd nb102_fjzH2(%esp),%xmm5
        mulpd nb102_dxH2H2(%esp),%xmm0
        mulpd nb102_dyH2H2(%esp),%xmm1
        mulpd nb102_dzH2H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb102_fixH2(%esp),%xmm0
        addpd nb102_fiyH2(%esp),%xmm1
        addpd nb102_fizH2(%esp),%xmm2
        movapd %xmm3,nb102_fjxH2(%esp)
        movapd %xmm4,nb102_fjyH2(%esp)
        movapd %xmm5,nb102_fjzH2(%esp)
        movapd %xmm0,nb102_fixH2(%esp)
        movapd %xmm1,nb102_fiyH2(%esp)
        movapd %xmm2,nb102_fizH2(%esp)

        movl nb102_faction(%ebp),%edi

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
        addpd nb102_fjxO(%esp),%xmm0
        addpd nb102_fjyO(%esp),%xmm1
        addpd nb102_fjzO(%esp),%xmm2
        addpd nb102_fjxH1(%esp),%xmm3
        addpd nb102_fjyH1(%esp),%xmm4
        addpd nb102_fjzH1(%esp),%xmm5
        addpd nb102_fjxH2(%esp),%xmm6
        addpd nb102_fjyH2(%esp),%xmm7
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
        addpd nb102_fjzH2(%esp),%xmm0
        movlpd %xmm0,64(%edi,%eax,8)
        movhpd %xmm0,64(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl $2,nb102_innerk(%esp)
        jl    _nb_kernel102_ia32_sse2.nb102_checksingle
        jmp   _nb_kernel102_ia32_sse2.nb102_unroll_loop
_nb_kernel102_ia32_sse2.nb102_checksingle: 
        movl  nb102_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel102_ia32_sse2.nb102_dosingle
        jmp   _nb_kernel102_ia32_sse2.nb102_updateouterdata
_nb_kernel102_ia32_sse2.nb102_dosingle: 
        movl  nb102_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb102_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coordinates to local temp variables 
        movlpd (%esi,%eax,8),%xmm2
        movlpd 8(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%eax,8),%xmm4
        movlpd 24(%esi,%eax,8),%xmm5
        movlpd 32(%esi,%eax,8),%xmm6
        movlpd 40(%esi,%eax,8),%xmm7
        movapd  %xmm2,nb102_jxO(%esp)
        movapd  %xmm3,nb102_jyO(%esp)
        movapd  %xmm4,nb102_jzO(%esp)
        movapd  %xmm5,nb102_jxH1(%esp)
        movapd  %xmm6,nb102_jyH1(%esp)
        movapd  %xmm7,nb102_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movapd  %xmm2,nb102_jxH2(%esp)
        movapd  %xmm3,nb102_jyH2(%esp)
        movapd  %xmm4,nb102_jzH2(%esp)

        movapd nb102_ixO(%esp),%xmm0
        movapd nb102_iyO(%esp),%xmm1
        movapd nb102_izO(%esp),%xmm2
        movapd nb102_ixO(%esp),%xmm3
        movapd nb102_iyO(%esp),%xmm4
        movapd nb102_izO(%esp),%xmm5
        subsd  nb102_jxO(%esp),%xmm0
        subsd  nb102_jyO(%esp),%xmm1
        subsd  nb102_jzO(%esp),%xmm2
        subsd  nb102_jxH1(%esp),%xmm3
        subsd  nb102_jyH1(%esp),%xmm4
        subsd  nb102_jzH1(%esp),%xmm5
        movapd %xmm0,nb102_dxOO(%esp)
        movapd %xmm1,nb102_dyOO(%esp)
        movapd %xmm2,nb102_dzOO(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb102_dxOH1(%esp)
        movapd %xmm4,nb102_dyOH1(%esp)
        movapd %xmm5,nb102_dzOH1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb102_rsqOO(%esp)
        movapd %xmm3,nb102_rsqOH1(%esp)

        movapd nb102_ixO(%esp),%xmm0
        movapd nb102_iyO(%esp),%xmm1
        movapd nb102_izO(%esp),%xmm2
        movapd nb102_ixH1(%esp),%xmm3
        movapd nb102_iyH1(%esp),%xmm4
        movapd nb102_izH1(%esp),%xmm5
        subsd  nb102_jxH2(%esp),%xmm0
        subsd  nb102_jyH2(%esp),%xmm1
        subsd  nb102_jzH2(%esp),%xmm2
        subsd  nb102_jxO(%esp),%xmm3
        subsd  nb102_jyO(%esp),%xmm4
        subsd  nb102_jzO(%esp),%xmm5
        movapd %xmm0,nb102_dxOH2(%esp)
        movapd %xmm1,nb102_dyOH2(%esp)
        movapd %xmm2,nb102_dzOH2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb102_dxH1O(%esp)
        movapd %xmm4,nb102_dyH1O(%esp)
        movapd %xmm5,nb102_dzH1O(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb102_rsqOH2(%esp)
        movapd %xmm3,nb102_rsqH1O(%esp)

        movapd nb102_ixH1(%esp),%xmm0
        movapd nb102_iyH1(%esp),%xmm1
        movapd nb102_izH1(%esp),%xmm2
        movapd nb102_ixH1(%esp),%xmm3
        movapd nb102_iyH1(%esp),%xmm4
        movapd nb102_izH1(%esp),%xmm5
        subsd  nb102_jxH1(%esp),%xmm0
        subsd  nb102_jyH1(%esp),%xmm1
        subsd  nb102_jzH1(%esp),%xmm2
        subsd  nb102_jxH2(%esp),%xmm3
        subsd  nb102_jyH2(%esp),%xmm4
        subsd  nb102_jzH2(%esp),%xmm5
        movapd %xmm0,nb102_dxH1H1(%esp)
        movapd %xmm1,nb102_dyH1H1(%esp)
        movapd %xmm2,nb102_dzH1H1(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb102_dxH1H2(%esp)
        movapd %xmm4,nb102_dyH1H2(%esp)
        movapd %xmm5,nb102_dzH1H2(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb102_rsqH1H1(%esp)
        movapd %xmm3,nb102_rsqH1H2(%esp)

        movapd nb102_ixH2(%esp),%xmm0
        movapd nb102_iyH2(%esp),%xmm1
        movapd nb102_izH2(%esp),%xmm2
        movapd nb102_ixH2(%esp),%xmm3
        movapd nb102_iyH2(%esp),%xmm4
        movapd nb102_izH2(%esp),%xmm5
        subsd  nb102_jxO(%esp),%xmm0
        subsd  nb102_jyO(%esp),%xmm1
        subsd  nb102_jzO(%esp),%xmm2
        subsd  nb102_jxH1(%esp),%xmm3
        subsd  nb102_jyH1(%esp),%xmm4
        subsd  nb102_jzH1(%esp),%xmm5
        movapd %xmm0,nb102_dxH2O(%esp)
        movapd %xmm1,nb102_dyH2O(%esp)
        movapd %xmm2,nb102_dzH2O(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb102_dxH2H1(%esp)
        movapd %xmm4,nb102_dyH2H1(%esp)
        movapd %xmm5,nb102_dzH2H1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb102_rsqH2O(%esp)
        movapd %xmm4,nb102_rsqH2H1(%esp)

        movapd nb102_ixH2(%esp),%xmm0
        movapd nb102_iyH2(%esp),%xmm1
        movapd nb102_izH2(%esp),%xmm2
        subsd  nb102_jxH2(%esp),%xmm0
        subsd  nb102_jyH2(%esp),%xmm1
        subsd  nb102_jzH2(%esp),%xmm2
        movapd %xmm0,nb102_dxH2H2(%esp)
        movapd %xmm1,nb102_dyH2H2(%esp)
        movapd %xmm2,nb102_dzH2H2(%esp)
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb102_rsqH2H2(%esp)

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
        movapd  nb102_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102_half(%esp),%xmm3   ## iter1 
        mulsd   nb102_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102_half(%esp),%xmm1   ## rinv 
        mulsd   nb102_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb102_rinvH2H2(%esp)
        movapd %xmm5,nb102_rinvH2H1(%esp)

        movapd nb102_rsqOO(%esp),%xmm0
        movapd nb102_rsqOH1(%esp),%xmm4
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
        movapd  nb102_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb102_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102_half(%esp),%xmm1   ## rinv 
        mulsd   nb102_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb102_rinvOO(%esp)
        movapd %xmm5,nb102_rinvOH1(%esp)

        movapd nb102_rsqOH2(%esp),%xmm0
        movapd nb102_rsqH1O(%esp),%xmm4
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
        movapd  nb102_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102_half(%esp),%xmm3   ## iter1 
        mulsd   nb102_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102_half(%esp),%xmm1   ## rinv 
        mulsd   nb102_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb102_rinvOH2(%esp)
        movapd %xmm5,nb102_rinvH1O(%esp)

        movapd nb102_rsqH1H1(%esp),%xmm0
        movapd nb102_rsqH1H2(%esp),%xmm4
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
        movapd  nb102_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102_half(%esp),%xmm3   ## iter1a 
        mulsd   nb102_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102_half(%esp),%xmm1   ## rinv 
        mulsd   nb102_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb102_rinvH1H1(%esp)
        movapd %xmm5,nb102_rinvH1H2(%esp)

        movapd nb102_rsqH2O(%esp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb102_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb102_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb102_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb102_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb102_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb102_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm7
        mulsd  %xmm0,%xmm0
        mulsd  nb102_qqOO(%esp),%xmm7
        mulsd  %xmm7,%xmm0
        addsd  nb102_vctot(%esp),%xmm7
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb102_dxOO(%esp),%xmm0
        mulsd nb102_dyOO(%esp),%xmm1
        mulsd nb102_dzOO(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb102_fixO(%esp),%xmm0
        addsd nb102_fiyO(%esp),%xmm1
        addsd nb102_fizO(%esp),%xmm2
        movlpd %xmm3,nb102_fjxO(%esp)
        movlpd %xmm4,nb102_fjyO(%esp)
        movlpd %xmm5,nb102_fjzO(%esp)
        movlpd %xmm0,nb102_fixO(%esp)
        movlpd %xmm1,nb102_fiyO(%esp)
        movlpd %xmm2,nb102_fizO(%esp)

        ## O-H1 interaction 
        movapd nb102_rinvOH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb102_qqOH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsOH1  
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb102_dxOH1(%esp),%xmm0
        mulsd nb102_dyOH1(%esp),%xmm1
        mulsd nb102_dzOH1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb102_fixO(%esp),%xmm0
        addsd nb102_fiyO(%esp),%xmm1
        addsd nb102_fizO(%esp),%xmm2
        movlpd %xmm3,nb102_fjxH1(%esp)
        movlpd %xmm4,nb102_fjyH1(%esp)
        movlpd %xmm5,nb102_fjzH1(%esp)
        movlpd %xmm0,nb102_fixO(%esp)
        movlpd %xmm1,nb102_fiyO(%esp)
        movlpd %xmm2,nb102_fizO(%esp)

        ## O-H2 interaction  
        movapd nb102_rinvOH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb102_qqOH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsOH2  
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb102_dxOH2(%esp),%xmm0
        mulsd nb102_dyOH2(%esp),%xmm1
        mulsd nb102_dzOH2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb102_fixO(%esp),%xmm0
        addsd nb102_fiyO(%esp),%xmm1
        addsd nb102_fizO(%esp),%xmm2
        movlpd %xmm3,nb102_fjxH2(%esp)
        movlpd %xmm4,nb102_fjyH2(%esp)
        movlpd %xmm5,nb102_fjzH2(%esp)
        movlpd %xmm0,nb102_fixO(%esp)
        movlpd %xmm1,nb102_fiyO(%esp)
        movlpd %xmm2,nb102_fizO(%esp)

        ## H1-O interaction 
        movapd nb102_rinvH1O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb102_qqOH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsH1O 
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb102_fjxO(%esp),%xmm3
        movapd nb102_fjyO(%esp),%xmm4
        movapd nb102_fjzO(%esp),%xmm5
        mulsd nb102_dxH1O(%esp),%xmm0
        mulsd nb102_dyH1O(%esp),%xmm1
        mulsd nb102_dzH1O(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb102_fixH1(%esp),%xmm0
        addsd nb102_fiyH1(%esp),%xmm1
        addsd nb102_fizH1(%esp),%xmm2
        movlpd %xmm3,nb102_fjxO(%esp)
        movlpd %xmm4,nb102_fjyO(%esp)
        movlpd %xmm5,nb102_fjzO(%esp)
        movlpd %xmm0,nb102_fixH1(%esp)
        movlpd %xmm1,nb102_fiyH1(%esp)
        movlpd %xmm2,nb102_fizH1(%esp)

        ## H1-H1 interaction 
        movapd nb102_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb102_qqHH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsH1H1 
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb102_fjxH1(%esp),%xmm3
        movapd nb102_fjyH1(%esp),%xmm4
        movapd nb102_fjzH1(%esp),%xmm5
        mulsd nb102_dxH1H1(%esp),%xmm0
        mulsd nb102_dyH1H1(%esp),%xmm1
        mulsd nb102_dzH1H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb102_fixH1(%esp),%xmm0
        addsd nb102_fiyH1(%esp),%xmm1
        addsd nb102_fizH1(%esp),%xmm2
        movlpd %xmm3,nb102_fjxH1(%esp)
        movlpd %xmm4,nb102_fjyH1(%esp)
        movlpd %xmm5,nb102_fjzH1(%esp)
        movlpd %xmm0,nb102_fixH1(%esp)
        movlpd %xmm1,nb102_fiyH1(%esp)
        movlpd %xmm2,nb102_fizH1(%esp)

        ## H1-H2 interaction 
        movapd nb102_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb102_qqHH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsOH2  
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb102_fjxH2(%esp),%xmm3
        movapd nb102_fjyH2(%esp),%xmm4
        movapd nb102_fjzH2(%esp),%xmm5
        mulsd nb102_dxH1H2(%esp),%xmm0
        mulsd nb102_dyH1H2(%esp),%xmm1
        mulsd nb102_dzH1H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb102_fixH1(%esp),%xmm0
        addsd nb102_fiyH1(%esp),%xmm1
        addsd nb102_fizH1(%esp),%xmm2
        movlpd %xmm3,nb102_fjxH2(%esp)
        movlpd %xmm4,nb102_fjyH2(%esp)
        movlpd %xmm5,nb102_fjzH2(%esp)
        movlpd %xmm0,nb102_fixH1(%esp)
        movlpd %xmm1,nb102_fiyH1(%esp)
        movlpd %xmm2,nb102_fizH1(%esp)

        ## H2-O interaction 
        movapd nb102_rinvH2O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb102_qqOH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsH2O 
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb102_fjxO(%esp),%xmm3
        movapd nb102_fjyO(%esp),%xmm4
        movapd nb102_fjzO(%esp),%xmm5
        mulsd nb102_dxH2O(%esp),%xmm0
        mulsd nb102_dyH2O(%esp),%xmm1
        mulsd nb102_dzH2O(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb102_fixH2(%esp),%xmm0
        addsd nb102_fiyH2(%esp),%xmm1
        addsd nb102_fizH2(%esp),%xmm2
        movlpd %xmm3,nb102_fjxO(%esp)
        movlpd %xmm4,nb102_fjyO(%esp)
        movlpd %xmm5,nb102_fjzO(%esp)
        movlpd %xmm0,nb102_fixH2(%esp)
        movlpd %xmm1,nb102_fiyH2(%esp)
        movlpd %xmm2,nb102_fizH2(%esp)

        ## H2-H1 interaction 
        movapd nb102_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb102_qqHH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsH2H1 
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb102_fjxH1(%esp),%xmm3
        movapd nb102_fjyH1(%esp),%xmm4
        movapd nb102_fjzH1(%esp),%xmm5
        mulsd nb102_dxH2H1(%esp),%xmm0
        mulsd nb102_dyH2H1(%esp),%xmm1
        mulsd nb102_dzH2H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb102_fixH2(%esp),%xmm0
        addsd nb102_fiyH2(%esp),%xmm1
        addsd nb102_fizH2(%esp),%xmm2
        movlpd %xmm3,nb102_fjxH1(%esp)
        movlpd %xmm4,nb102_fjyH1(%esp)
        movlpd %xmm5,nb102_fjzH1(%esp)
        movlpd %xmm0,nb102_fixH2(%esp)
        movlpd %xmm1,nb102_fiyH2(%esp)
        movlpd %xmm2,nb102_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb102_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb102_qqHH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsH2H2 
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movsd %xmm7,nb102_vctot(%esp)
        movapd %xmm0,%xmm2
        movapd nb102_fjxH2(%esp),%xmm3
        movapd nb102_fjyH2(%esp),%xmm4
        movapd nb102_fjzH2(%esp),%xmm5
        mulsd nb102_dxH2H2(%esp),%xmm0
        mulsd nb102_dyH2H2(%esp),%xmm1
        mulsd nb102_dzH2H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb102_fixH2(%esp),%xmm0
        addsd nb102_fiyH2(%esp),%xmm1
        addsd nb102_fizH2(%esp),%xmm2
        movlpd %xmm3,nb102_fjxH2(%esp)
        movlpd %xmm4,nb102_fjyH2(%esp)
        movlpd %xmm5,nb102_fjzH2(%esp)
        movlpd %xmm0,nb102_fixH2(%esp)
        movlpd %xmm1,nb102_fiyH2(%esp)
        movlpd %xmm2,nb102_fizH2(%esp)

        movl nb102_faction(%ebp),%edi

        ## Did all interactions - now update j forces 
        movlpd (%edi,%eax,8),%xmm0
        movlpd 8(%edi,%eax,8),%xmm1
        movlpd 16(%edi,%eax,8),%xmm2
        movlpd 24(%edi,%eax,8),%xmm3
        movlpd 32(%edi,%eax,8),%xmm4
        movlpd 40(%edi,%eax,8),%xmm5
        movlpd 48(%edi,%eax,8),%xmm6
        movlpd 56(%edi,%eax,8),%xmm7
        addsd nb102_fjxO(%esp),%xmm0
        addsd nb102_fjyO(%esp),%xmm1
        addsd nb102_fjzO(%esp),%xmm2
        addsd nb102_fjxH1(%esp),%xmm3
        addsd nb102_fjyH1(%esp),%xmm4
        addsd nb102_fjzH1(%esp),%xmm5
        addsd nb102_fjxH2(%esp),%xmm6
        addsd nb102_fjyH2(%esp),%xmm7
        movlpd %xmm0,(%edi,%eax,8)
        movlpd %xmm1,8(%edi,%eax,8)
        movlpd %xmm2,16(%edi,%eax,8)
        movlpd %xmm3,24(%edi,%eax,8)
        movlpd %xmm4,32(%edi,%eax,8)
        movlpd %xmm5,40(%edi,%eax,8)
        movlpd %xmm6,48(%edi,%eax,8)
        movlpd %xmm7,56(%edi,%eax,8)

        movlpd 64(%edi,%eax,8),%xmm0
        addsd nb102_fjzH2(%esp),%xmm0
        movlpd %xmm0,64(%edi,%eax,8)

_nb_kernel102_ia32_sse2.nb102_updateouterdata: 
        movl  nb102_ii3(%esp),%ecx
        movl  nb102_faction(%ebp),%edi
        movl  nb102_fshift(%ebp),%esi
        movl  nb102_is3(%esp),%edx

        ## accumulate Oi forces in xmm0, xmm1, xmm2 
        movapd nb102_fixO(%esp),%xmm0
        movapd nb102_fiyO(%esp),%xmm1
        movapd nb102_fizO(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

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
        movapd nb102_fixH1(%esp),%xmm0
        movapd nb102_fiyH1(%esp),%xmm1
        movapd nb102_fizH1(%esp),%xmm2

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
        movapd nb102_fixH2(%esp),%xmm0
        movapd nb102_fiyH2(%esp),%xmm1
        movapd nb102_fizH2(%esp),%xmm2

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
        movl nb102_n(%esp),%esi
        ## get group index for i particle 
        movl  nb102_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb102_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 

        ## add earlier value from mem 
        movl  nb102_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb102_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel102_ia32_sse2.nb102_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb102_n(%esp)
        jmp _nb_kernel102_ia32_sse2.nb102_outer
_nb_kernel102_ia32_sse2.nb102_outerend: 
        ## check if more outer neighborlists remain
        movl  nb102_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel102_ia32_sse2.nb102_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel102_ia32_sse2.nb102_threadloop
_nb_kernel102_ia32_sse2.nb102_end: 
        emms

        movl nb102_nouter(%esp),%eax
        movl nb102_ninner(%esp),%ebx
        movl nb102_outeriter(%ebp),%ecx
        movl nb102_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb102_salign(%esp),%eax
        addl %eax,%esp
        addl $1432,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



.globl nb_kernel102nf_ia32_sse2
.globl _nb_kernel102nf_ia32_sse2
nb_kernel102nf_ia32_sse2:       
_nb_kernel102nf_ia32_sse2:      
.set nb102nf_p_nri, 8
.set nb102nf_iinr, 12
.set nb102nf_jindex, 16
.set nb102nf_jjnr, 20
.set nb102nf_shift, 24
.set nb102nf_shiftvec, 28
.set nb102nf_fshift, 32
.set nb102nf_gid, 36
.set nb102nf_pos, 40
.set nb102nf_faction, 44
.set nb102nf_charge, 48
.set nb102nf_p_facel, 52
.set nb102nf_argkrf, 56
.set nb102nf_argcrf, 60
.set nb102nf_Vc, 64
.set nb102nf_type, 68
.set nb102nf_p_ntype, 72
.set nb102nf_vdwparam, 76
.set nb102nf_Vvdw, 80
.set nb102nf_p_tabscale, 84
.set nb102nf_VFtab, 88
.set nb102nf_invsqrta, 92
.set nb102nf_dvda, 96
.set nb102nf_p_gbtabscale, 100
.set nb102nf_GBtab, 104
.set nb102nf_p_nthreads, 108
.set nb102nf_count, 112
.set nb102nf_mtx, 116
.set nb102nf_outeriter, 120
.set nb102nf_inneriter, 124
.set nb102nf_work, 128
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
.set nb102nf_innerjjnr, 680
.set nb102nf_innerk, 684
.set nb102nf_n, 688
.set nb102nf_nn1, 692
.set nb102nf_nri, 696
.set nb102nf_nouter, 700
.set nb102nf_ninner, 704
.set nb102nf_salign, 708
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $712,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb102nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb102nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb102nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb102nf_nouter(%esp)
        movl %eax,nb102nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb102nf_half(%esp)
        movl %ebx,nb102nf_half+4(%esp)
        movsd nb102nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb102nf_half(%esp)
        movapd %xmm3,nb102nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb102nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb102nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3       ## qO 
        movsd %xmm3,%xmm4               ## qO 
        movsd 8(%edx,%ebx,8),%xmm5      ## qH 
        movl nb102nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm6      ## facel 
        mulsd  %xmm3,%xmm3              ## qO*qO 
        mulsd  %xmm5,%xmm4              ## qO*qH 
        mulsd  %xmm5,%xmm5              ## qH*qH 
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb102nf_qqOO(%esp)
        movapd %xmm4,nb102nf_qqOH(%esp)
        movapd %xmm5,nb102nf_qqHH(%esp)

_nb_kernel102nf_ia32_sse2.nb102nf_threadloop: 
        movl  nb102nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel102nf_ia32_sse2.nb102nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel102nf_ia32_sse2.nb102nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb102nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb102nf_n(%esp)
        movl %ebx,nb102nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel102nf_ia32_sse2.nb102nf_outerstart
        jmp _nb_kernel102nf_ia32_sse2.nb102nf_end

_nb_kernel102nf_ia32_sse2.nb102nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb102nf_nouter(%esp),%ebx
        movl %ebx,nb102nf_nouter(%esp)

_nb_kernel102nf_ia32_sse2.nb102nf_outer: 
        movl  nb102nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb102nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb102nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb102nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb102nf_ii3(%esp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb102nf_ixO(%esp)
        movapd %xmm4,nb102nf_iyO(%esp)
        movapd %xmm5,nb102nf_izO(%esp)

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
        movapd %xmm0,nb102nf_ixH1(%esp)
        movapd %xmm1,nb102nf_iyH1(%esp)
        movapd %xmm2,nb102nf_izH1(%esp)
        movapd %xmm3,nb102nf_ixH2(%esp)
        movapd %xmm4,nb102nf_iyH2(%esp)
        movapd %xmm5,nb102nf_izH2(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb102nf_vctot(%esp)

        movl  nb102nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb102nf_pos(%ebp),%esi
        movl  nb102nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb102nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb102nf_ninner(%esp),%ecx
        movl  %ecx,nb102nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb102nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel102nf_ia32_sse2.nb102nf_unroll_loop
        jmp   _nb_kernel102nf_ia32_sse2.nb102nf_checksingle
_nb_kernel102nf_ia32_sse2.nb102nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb102nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb102nf_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb102nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd  %xmm2,nb102nf_jxO(%esp)
        movapd  %xmm3,nb102nf_jyO(%esp)
        movapd  %xmm4,nb102nf_jzO(%esp)
        movapd  %xmm5,nb102nf_jxH1(%esp)
        movapd  %xmm6,nb102nf_jyH1(%esp)
        movapd  %xmm7,nb102nf_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 48(%esi,%ebx,8),%xmm2
        movhpd 56(%esi,%ebx,8),%xmm3
        movhpd 64(%esi,%ebx,8),%xmm4
        movapd  %xmm2,nb102nf_jxH2(%esp)
        movapd  %xmm3,nb102nf_jyH2(%esp)
        movapd  %xmm4,nb102nf_jzH2(%esp)

        movapd nb102nf_ixO(%esp),%xmm0
        movapd nb102nf_iyO(%esp),%xmm1
        movapd nb102nf_izO(%esp),%xmm2
        movapd nb102nf_ixO(%esp),%xmm3
        movapd nb102nf_iyO(%esp),%xmm4
        movapd nb102nf_izO(%esp),%xmm5
        subpd  nb102nf_jxO(%esp),%xmm0
        subpd  nb102nf_jyO(%esp),%xmm1
        subpd  nb102nf_jzO(%esp),%xmm2
        subpd  nb102nf_jxH1(%esp),%xmm3
        subpd  nb102nf_jyH1(%esp),%xmm4
        subpd  nb102nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb102nf_rsqOO(%esp)
        movapd %xmm3,nb102nf_rsqOH1(%esp)

        movapd nb102nf_ixO(%esp),%xmm0
        movapd nb102nf_iyO(%esp),%xmm1
        movapd nb102nf_izO(%esp),%xmm2
        movapd nb102nf_ixH1(%esp),%xmm3
        movapd nb102nf_iyH1(%esp),%xmm4
        movapd nb102nf_izH1(%esp),%xmm5
        subpd  nb102nf_jxH2(%esp),%xmm0
        subpd  nb102nf_jyH2(%esp),%xmm1
        subpd  nb102nf_jzH2(%esp),%xmm2
        subpd  nb102nf_jxO(%esp),%xmm3
        subpd  nb102nf_jyO(%esp),%xmm4
        subpd  nb102nf_jzO(%esp),%xmm5
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
        movapd %xmm0,nb102nf_rsqOH2(%esp)
        movapd %xmm3,nb102nf_rsqH1O(%esp)

        movapd nb102nf_ixH1(%esp),%xmm0
        movapd nb102nf_iyH1(%esp),%xmm1
        movapd nb102nf_izH1(%esp),%xmm2
        movapd nb102nf_ixH1(%esp),%xmm3
        movapd nb102nf_iyH1(%esp),%xmm4
        movapd nb102nf_izH1(%esp),%xmm5
        subpd  nb102nf_jxH1(%esp),%xmm0
        subpd  nb102nf_jyH1(%esp),%xmm1
        subpd  nb102nf_jzH1(%esp),%xmm2
        subpd  nb102nf_jxH2(%esp),%xmm3
        subpd  nb102nf_jyH2(%esp),%xmm4
        subpd  nb102nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb102nf_rsqH1H1(%esp)
        movapd %xmm3,nb102nf_rsqH1H2(%esp)

        movapd nb102nf_ixH2(%esp),%xmm0
        movapd nb102nf_iyH2(%esp),%xmm1
        movapd nb102nf_izH2(%esp),%xmm2
        movapd nb102nf_ixH2(%esp),%xmm3
        movapd nb102nf_iyH2(%esp),%xmm4
        movapd nb102nf_izH2(%esp),%xmm5
        subpd  nb102nf_jxO(%esp),%xmm0
        subpd  nb102nf_jyO(%esp),%xmm1
        subpd  nb102nf_jzO(%esp),%xmm2
        subpd  nb102nf_jxH1(%esp),%xmm3
        subpd  nb102nf_jyH1(%esp),%xmm4
        subpd  nb102nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb102nf_rsqH2O(%esp)
        movapd %xmm4,nb102nf_rsqH2H1(%esp)

        movapd nb102nf_ixH2(%esp),%xmm0
        movapd nb102nf_iyH2(%esp),%xmm1
        movapd nb102nf_izH2(%esp),%xmm2
        subpd  nb102nf_jxH2(%esp),%xmm0
        subpd  nb102nf_jyH2(%esp),%xmm1
        subpd  nb102nf_jzH2(%esp),%xmm2
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb102nf_rsqH2H2(%esp)

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
        movapd  nb102nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb102nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb102nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb102nf_rinvH2H2(%esp)
        movapd %xmm5,nb102nf_rinvH2H1(%esp)

        movapd nb102nf_rsqOO(%esp),%xmm0
        movapd nb102nf_rsqOH1(%esp),%xmm4
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
        movapd  nb102nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102nf_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb102nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb102nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb102nf_rinvOO(%esp)
        movapd %xmm5,nb102nf_rinvOH1(%esp)

        movapd nb102nf_rsqOH2(%esp),%xmm0
        movapd nb102nf_rsqH1O(%esp),%xmm4
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
        movapd  nb102nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb102nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb102nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb102nf_rinvOH2(%esp)
        movapd %xmm5,nb102nf_rinvH1O(%esp)

        movapd nb102nf_rsqH1H1(%esp),%xmm0
        movapd nb102nf_rsqH1H2(%esp),%xmm4
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
        movapd  nb102nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb102nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb102nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb102nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb102nf_rinvH1H1(%esp)
        movapd %xmm5,nb102nf_rinvH1H2(%esp)

        movapd nb102nf_rsqH2O(%esp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb102nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb102nf_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb102nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb102nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb102nf_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb102nf_rinvOO(%esp),%xmm0
        mulpd  nb102nf_qqOO(%esp),%xmm0
        addpd  nb102nf_vctot(%esp),%xmm0

        ## other interactions 
        movapd nb102nf_rinvOH1(%esp),%xmm1
        movapd nb102nf_rinvH1H1(%esp),%xmm2

        addpd nb102nf_rinvOH2(%esp),%xmm1
        addpd nb102nf_rinvH1H2(%esp),%xmm2

        addpd nb102nf_rinvH1O(%esp),%xmm1
        addpd nb102nf_rinvH2H1(%esp),%xmm2

        addpd nb102nf_rinvH2O(%esp),%xmm1
        addpd nb102nf_rinvH2H2(%esp),%xmm2

        mulpd nb102nf_qqOH(%esp),%xmm1
        mulpd nb102nf_qqHH(%esp),%xmm2

        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0

        movapd %xmm0,nb102nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb102nf_innerk(%esp)
        jl    _nb_kernel102nf_ia32_sse2.nb102nf_checksingle
        jmp   _nb_kernel102nf_ia32_sse2.nb102nf_unroll_loop
_nb_kernel102nf_ia32_sse2.nb102nf_checksingle: 
        movl  nb102nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel102nf_ia32_sse2.nb102nf_dosingle
        jmp   _nb_kernel102nf_ia32_sse2.nb102nf_updateouterdata
_nb_kernel102nf_ia32_sse2.nb102nf_dosingle: 
        movl  nb102nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb102nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coordinates to local temp variables 
        movlpd (%esi,%eax,8),%xmm2
        movlpd 8(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%eax,8),%xmm4
        movlpd 24(%esi,%eax,8),%xmm5
        movlpd 32(%esi,%eax,8),%xmm6
        movlpd 40(%esi,%eax,8),%xmm7
        movapd  %xmm2,nb102nf_jxO(%esp)
        movapd  %xmm3,nb102nf_jyO(%esp)
        movapd  %xmm4,nb102nf_jzO(%esp)
        movapd  %xmm5,nb102nf_jxH1(%esp)
        movapd  %xmm6,nb102nf_jyH1(%esp)
        movapd  %xmm7,nb102nf_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movapd  %xmm2,nb102nf_jxH2(%esp)
        movapd  %xmm3,nb102nf_jyH2(%esp)
        movapd  %xmm4,nb102nf_jzH2(%esp)

        movapd nb102nf_ixO(%esp),%xmm0
        movapd nb102nf_iyO(%esp),%xmm1
        movapd nb102nf_izO(%esp),%xmm2
        movapd nb102nf_ixO(%esp),%xmm3
        movapd nb102nf_iyO(%esp),%xmm4
        movapd nb102nf_izO(%esp),%xmm5
        subsd  nb102nf_jxO(%esp),%xmm0
        subsd  nb102nf_jyO(%esp),%xmm1
        subsd  nb102nf_jzO(%esp),%xmm2
        subsd  nb102nf_jxH1(%esp),%xmm3
        subsd  nb102nf_jyH1(%esp),%xmm4
        subsd  nb102nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb102nf_rsqOO(%esp)
        movapd %xmm3,nb102nf_rsqOH1(%esp)

        movapd nb102nf_ixO(%esp),%xmm0
        movapd nb102nf_iyO(%esp),%xmm1
        movapd nb102nf_izO(%esp),%xmm2
        movapd nb102nf_ixH1(%esp),%xmm3
        movapd nb102nf_iyH1(%esp),%xmm4
        movapd nb102nf_izH1(%esp),%xmm5
        subsd  nb102nf_jxH2(%esp),%xmm0
        subsd  nb102nf_jyH2(%esp),%xmm1
        subsd  nb102nf_jzH2(%esp),%xmm2
        subsd  nb102nf_jxO(%esp),%xmm3
        subsd  nb102nf_jyO(%esp),%xmm4
        subsd  nb102nf_jzO(%esp),%xmm5
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
        movapd %xmm0,nb102nf_rsqOH2(%esp)
        movapd %xmm3,nb102nf_rsqH1O(%esp)

        movapd nb102nf_ixH1(%esp),%xmm0
        movapd nb102nf_iyH1(%esp),%xmm1
        movapd nb102nf_izH1(%esp),%xmm2
        movapd nb102nf_ixH1(%esp),%xmm3
        movapd nb102nf_iyH1(%esp),%xmm4
        movapd nb102nf_izH1(%esp),%xmm5
        subsd  nb102nf_jxH1(%esp),%xmm0
        subsd  nb102nf_jyH1(%esp),%xmm1
        subsd  nb102nf_jzH1(%esp),%xmm2
        subsd  nb102nf_jxH2(%esp),%xmm3
        subsd  nb102nf_jyH2(%esp),%xmm4
        subsd  nb102nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb102nf_rsqH1H1(%esp)
        movapd %xmm3,nb102nf_rsqH1H2(%esp)

        movapd nb102nf_ixH2(%esp),%xmm0
        movapd nb102nf_iyH2(%esp),%xmm1
        movapd nb102nf_izH2(%esp),%xmm2
        movapd nb102nf_ixH2(%esp),%xmm3
        movapd nb102nf_iyH2(%esp),%xmm4
        movapd nb102nf_izH2(%esp),%xmm5
        subsd  nb102nf_jxO(%esp),%xmm0
        subsd  nb102nf_jyO(%esp),%xmm1
        subsd  nb102nf_jzO(%esp),%xmm2
        subsd  nb102nf_jxH1(%esp),%xmm3
        subsd  nb102nf_jyH1(%esp),%xmm4
        subsd  nb102nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb102nf_rsqH2O(%esp)
        movapd %xmm4,nb102nf_rsqH2H1(%esp)

        movapd nb102nf_ixH2(%esp),%xmm0
        movapd nb102nf_iyH2(%esp),%xmm1
        movapd nb102nf_izH2(%esp),%xmm2
        subsd  nb102nf_jxH2(%esp),%xmm0
        subsd  nb102nf_jyH2(%esp),%xmm1
        subsd  nb102nf_jzH2(%esp),%xmm2
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb102nf_rsqH2H2(%esp)

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
        movapd  nb102nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb102nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb102nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb102nf_rinvH2H2(%esp)
        movapd %xmm5,nb102nf_rinvH2H1(%esp)

        movapd nb102nf_rsqOO(%esp),%xmm0
        movapd nb102nf_rsqOH1(%esp),%xmm4
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
        movapd  nb102nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102nf_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb102nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb102nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb102nf_rinvOO(%esp)
        movapd %xmm5,nb102nf_rinvOH1(%esp)

        movapd nb102nf_rsqOH2(%esp),%xmm0
        movapd nb102nf_rsqH1O(%esp),%xmm4
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
        movapd  nb102nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb102nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb102nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb102nf_rinvOH2(%esp)
        movapd %xmm5,nb102nf_rinvH1O(%esp)

        movapd nb102nf_rsqH1H1(%esp),%xmm0
        movapd nb102nf_rsqH1H2(%esp),%xmm4
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
        movapd  nb102nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb102nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb102nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb102nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb102nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb102nf_rinvH1H1(%esp)
        movapd %xmm5,nb102nf_rinvH1H2(%esp)

        movapd nb102nf_rsqH2O(%esp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb102nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb102nf_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb102nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb102nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb102nf_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb102nf_rinvOO(%esp),%xmm0
        mulpd  nb102nf_qqOO(%esp),%xmm0
        addpd  nb102nf_vctot(%esp),%xmm0

        ## other interactions 
        movapd nb102nf_rinvOH1(%esp),%xmm1
        movapd nb102nf_rinvH1H1(%esp),%xmm2

        addsd nb102nf_rinvOH2(%esp),%xmm1
        addsd nb102nf_rinvH1H2(%esp),%xmm2

        addsd nb102nf_rinvH1O(%esp),%xmm1
        addsd nb102nf_rinvH2H1(%esp),%xmm2

        addsd nb102nf_rinvH2O(%esp),%xmm1
        addsd nb102nf_rinvH2H2(%esp),%xmm2

        mulsd nb102nf_qqOH(%esp),%xmm1
        mulsd nb102nf_qqHH(%esp),%xmm2

        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0

        movlpd %xmm0,nb102nf_vctot(%esp)

_nb_kernel102nf_ia32_sse2.nb102nf_updateouterdata: 
        ## get n from stack
        movl nb102nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb102nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movapd nb102nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 

        ## add earlier value from mem 
        movl  nb102nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb102nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel102nf_ia32_sse2.nb102nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb102nf_n(%esp)
        jmp _nb_kernel102nf_ia32_sse2.nb102nf_outer
_nb_kernel102nf_ia32_sse2.nb102nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb102nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel102nf_ia32_sse2.nb102nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel102nf_ia32_sse2.nb102nf_threadloop
_nb_kernel102nf_ia32_sse2.nb102nf_end: 
        emms

        movl nb102nf_nouter(%esp),%eax
        movl nb102nf_ninner(%esp),%ebx
        movl nb102nf_outeriter(%ebp),%ecx
        movl nb102nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb102nf_salign(%esp),%eax
        addl %eax,%esp
        addl $712,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



