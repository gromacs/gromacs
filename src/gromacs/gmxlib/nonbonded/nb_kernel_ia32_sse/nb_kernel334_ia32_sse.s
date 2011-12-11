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

.globl nb_kernel334_ia32_sse
.globl _nb_kernel334_ia32_sse
nb_kernel334_ia32_sse:  
_nb_kernel334_ia32_sse: 
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
        ## bottom of stack is cache-aligned for sse use 
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
.set nb334_fstmp, 1744
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

        movl nb334_p_tabscale(%ebp),%eax
        movss (%eax),%xmm5
        shufps $0,%xmm5,%xmm5
        movaps %xmm5,nb334_tsc(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb334_half(%esp)
        movss nb334_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## 2.0
        addps  %xmm2,%xmm3      ## 3.0
        movaps %xmm1,nb334_half(%esp)
        movaps %xmm2,nb334_two(%esp)
        movaps %xmm3,nb334_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb334_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb334_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm5
        movss 12(%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movl nb334_p_facel(%ebp),%esi
        movss (%esi),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb334_qqMM(%esp)
        movaps %xmm4,nb334_qqMH(%esp)
        movaps %xmm5,nb334_qqHH(%esp)

        xorps %xmm0,%xmm0
        movl  nb334_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb334_p_ntype(%ebp),%edi
        imull (%edi),%ecx       ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb334_vdwparam(%ebp),%eax
        movlps (%eax,%edx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $0x55,%xmm1,%xmm1
        movaps %xmm0,nb334_c6(%esp)
        movaps %xmm1,nb334_c12(%esp)

_nb_kernel334_ia32_sse.nb334_threadloop: 
        movl  nb334_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel334_ia32_sse.nb334_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel334_ia32_sse.nb334_spinlock

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
        jg  _nb_kernel334_ia32_sse.nb334_outerstart
        jmp _nb_kernel334_ia32_sse.nb334_end

_nb_kernel334_ia32_sse.nb334_outerstart: 
        ## ebx contains number of outer iterations
        addl nb334_nouter(%esp),%ebx
        movl %ebx,nb334_nouter(%esp)

_nb_kernel334_ia32_sse.nb334_outer: 
        movl  nb334_shift(%ebp),%eax            ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb334_is3(%esp)      ## store is3 

        movl  nb334_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb334_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb334_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb334_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        addss (%eax,%ebx,4),%xmm3       ## ox
        addss 4(%eax,%ebx,4),%xmm4     ## oy
        addss 8(%eax,%ebx,4),%xmm5     ## oz
        addss 12(%eax,%ebx,4),%xmm6    ## h1x
        addss 16(%eax,%ebx,4),%xmm7    ## h1y
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm7,%xmm7
        movaps %xmm3,nb334_ixO(%esp)
        movaps %xmm4,nb334_iyO(%esp)
        movaps %xmm5,nb334_izO(%esp)
        movaps %xmm6,nb334_ixH1(%esp)
        movaps %xmm7,nb334_iyH1(%esp)

        movss %xmm2,%xmm6
        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 20(%eax,%ebx,4),%xmm6    ## h1z
        addss 24(%eax,%ebx,4),%xmm0    ## h2x
        addss 28(%eax,%ebx,4),%xmm1    ## h2y
        addss 32(%eax,%ebx,4),%xmm2    ## h2z
        addss 36(%eax,%ebx,4),%xmm3    ## mx
        addss 40(%eax,%ebx,4),%xmm4    ## my
        addss 44(%eax,%ebx,4),%xmm5    ## mz

        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm6,nb334_izH1(%esp)
        movaps %xmm0,nb334_ixH2(%esp)
        movaps %xmm1,nb334_iyH2(%esp)
        movaps %xmm2,nb334_izH2(%esp)
        movaps %xmm3,nb334_ixM(%esp)
        movaps %xmm4,nb334_iyM(%esp)
        movaps %xmm5,nb334_izM(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb334_vctot(%esp)
        movaps %xmm4,nb334_Vvdwtot(%esp)
        movaps %xmm4,nb334_fixO(%esp)
        movaps %xmm4,nb334_fiyO(%esp)
        movaps %xmm4,nb334_fizO(%esp)
        movaps %xmm4,nb334_fixH1(%esp)
        movaps %xmm4,nb334_fiyH1(%esp)
        movaps %xmm4,nb334_fizH1(%esp)
        movaps %xmm4,nb334_fixH2(%esp)
        movaps %xmm4,nb334_fiyH2(%esp)
        movaps %xmm4,nb334_fizH2(%esp)
        movaps %xmm4,nb334_fixM(%esp)
        movaps %xmm4,nb334_fiyM(%esp)
        movaps %xmm4,nb334_fizM(%esp)

        movl  nb334_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb334_pos(%ebp),%esi
        movl  nb334_faction(%ebp),%edi
        movl  nb334_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb334_innerjjnr(%esp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb334_ninner(%esp),%ecx
        movl  %ecx,nb334_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb334_innerk(%esp)   ## number of innerloop atoms 
        jge   _nb_kernel334_ia32_sse.nb334_unroll_loop
        jmp   _nb_kernel334_ia32_sse.nb334_single_check
_nb_kernel334_ia32_sse.nb334_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb334_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb334_innerjjnr(%esp)             ## advance pointer (unroll 4) 

        movl nb334_pos(%ebp),%esi       ## base of pos[] 

        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx        ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables
        ## Load Ox, Oy, Oz, H1x 
        movlps (%esi,%eax,4),%xmm1      ##  Oxa   Oya    -    -
        movlps (%esi,%ecx,4),%xmm4      ##  Oxc   Oyc    -    -
        movhps (%esi,%ebx,4),%xmm1      ##  Oxa   Oya   Oxb   Oyb 
        movhps (%esi,%edx,4),%xmm4      ##  Oxc   Oyc   Oxd   Oyd 
        movaps %xmm1,%xmm0              ##  Oxa   Oya   Oxb   Oyb 
        shufps $0x88,%xmm4,%xmm0       ##  Oxa   Oxb   Oxc   Oxd
        shufps $0xDD,%xmm4,%xmm1       ##  Oya   Oyb   Oyc   Oyd
        movlps 8(%esi,%eax,4),%xmm3     ##  Oza  H1xa    -    -
        movlps 8(%esi,%ecx,4),%xmm5     ##  Ozc  H1xc    -    -
        movhps 8(%esi,%ebx,4),%xmm3     ##  Oza  H1xa   Ozb  H1xb 
        movhps 8(%esi,%edx,4),%xmm5     ##  Ozc  H1xc   Ozd  H1xd 
        movaps %xmm3,%xmm2              ##  Oza  H1xa   Ozb  H1xb 
        shufps $0x88,%xmm5,%xmm2       ##  Oza   Ozb   Ozc   Ozd
        shufps $0xDD,%xmm5,%xmm3       ## H1xa  H1xb  H1xc  H1xd
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb334_jxO(%esp)
        movaps %xmm1,nb334_jyO(%esp)
        movaps %xmm2,nb334_jzO(%esp)
        movaps %xmm3,nb334_jxH1(%esp)

        ## Load H1y H1z H2x H2y 
        movlps 16(%esi,%eax,4),%xmm1
        movlps 16(%esi,%ecx,4),%xmm4
        movhps 16(%esi,%ebx,4),%xmm1
        movhps 16(%esi,%edx,4),%xmm4
        movaps %xmm1,%xmm0
        shufps $0x88,%xmm4,%xmm0
        shufps $0xDD,%xmm4,%xmm1
        movlps 24(%esi,%eax,4),%xmm3
        movlps 24(%esi,%ecx,4),%xmm5
        movhps 24(%esi,%ebx,4),%xmm3
        movhps 24(%esi,%edx,4),%xmm5
        movaps %xmm3,%xmm2
        shufps $0x88,%xmm5,%xmm2
        shufps $0xDD,%xmm5,%xmm3
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb334_jyH1(%esp)
        movaps %xmm1,nb334_jzH1(%esp)
        movaps %xmm2,nb334_jxH2(%esp)
        movaps %xmm3,nb334_jyH2(%esp)

        ## Load H2z Mx My Mz 
        movlps 32(%esi,%eax,4),%xmm1
        movlps 32(%esi,%ecx,4),%xmm4
        movhps 32(%esi,%ebx,4),%xmm1
        movhps 32(%esi,%edx,4),%xmm4
        movaps %xmm1,%xmm0
        shufps $0x88,%xmm4,%xmm0
        shufps $0xDD,%xmm4,%xmm1
        movlps 40(%esi,%eax,4),%xmm3
        movlps 40(%esi,%ecx,4),%xmm5
        movhps 40(%esi,%ebx,4),%xmm3
        movhps 40(%esi,%edx,4),%xmm5
        movaps %xmm3,%xmm2
        shufps $0x88,%xmm5,%xmm2
        shufps $0xDD,%xmm5,%xmm3
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb334_jzH2(%esp)
        movaps %xmm1,nb334_jxM(%esp)
        movaps %xmm2,nb334_jyM(%esp)
        movaps %xmm3,nb334_jzM(%esp)

        ## start calculating pairwise distances
        movaps nb334_ixO(%esp),%xmm0
        movaps nb334_iyO(%esp),%xmm1
        movaps nb334_izO(%esp),%xmm2
        movaps nb334_ixH1(%esp),%xmm3
        movaps nb334_iyH1(%esp),%xmm4
        movaps nb334_izH1(%esp),%xmm5
        subps  nb334_jxO(%esp),%xmm0
        subps  nb334_jyO(%esp),%xmm1
        subps  nb334_jzO(%esp),%xmm2
        subps  nb334_jxH1(%esp),%xmm3
        subps  nb334_jyH1(%esp),%xmm4
        subps  nb334_jzH1(%esp),%xmm5
        movaps %xmm0,nb334_dxOO(%esp)
        movaps %xmm1,nb334_dyOO(%esp)
        movaps %xmm2,nb334_dzOO(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb334_dxH1H1(%esp)
        movaps %xmm4,nb334_dyH1H1(%esp)
        movaps %xmm5,nb334_dzH1H1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb334_rsqOO(%esp)
        movaps %xmm3,nb334_rsqH1H1(%esp)

        movaps nb334_ixH1(%esp),%xmm0
        movaps nb334_iyH1(%esp),%xmm1
        movaps nb334_izH1(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb334_jxH2(%esp),%xmm0
        subps  nb334_jyH2(%esp),%xmm1
        subps  nb334_jzH2(%esp),%xmm2
        subps  nb334_jxM(%esp),%xmm3
        subps  nb334_jyM(%esp),%xmm4
        subps  nb334_jzM(%esp),%xmm5
        movaps %xmm0,nb334_dxH1H2(%esp)
        movaps %xmm1,nb334_dyH1H2(%esp)
        movaps %xmm2,nb334_dzH1H2(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb334_dxH1M(%esp)
        movaps %xmm4,nb334_dyH1M(%esp)
        movaps %xmm5,nb334_dzH1M(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb334_rsqH1H2(%esp)
        movaps %xmm3,nb334_rsqH1M(%esp)

        movaps nb334_ixH2(%esp),%xmm0
        movaps nb334_iyH2(%esp),%xmm1
        movaps nb334_izH2(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb334_jxH1(%esp),%xmm0
        subps  nb334_jyH1(%esp),%xmm1
        subps  nb334_jzH1(%esp),%xmm2
        subps  nb334_jxH2(%esp),%xmm3
        subps  nb334_jyH2(%esp),%xmm4
        subps  nb334_jzH2(%esp),%xmm5
        movaps %xmm0,nb334_dxH2H1(%esp)
        movaps %xmm1,nb334_dyH2H1(%esp)
        movaps %xmm2,nb334_dzH2H1(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb334_dxH2H2(%esp)
        movaps %xmm4,nb334_dyH2H2(%esp)
        movaps %xmm5,nb334_dzH2H2(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb334_rsqH2H1(%esp)
        movaps %xmm3,nb334_rsqH2H2(%esp)

        movaps nb334_ixH2(%esp),%xmm0
        movaps nb334_iyH2(%esp),%xmm1
        movaps nb334_izH2(%esp),%xmm2
        movaps nb334_ixM(%esp),%xmm3
        movaps nb334_iyM(%esp),%xmm4
        movaps nb334_izM(%esp),%xmm5
        subps  nb334_jxM(%esp),%xmm0
        subps  nb334_jyM(%esp),%xmm1
        subps  nb334_jzM(%esp),%xmm2
        subps  nb334_jxH1(%esp),%xmm3
        subps  nb334_jyH1(%esp),%xmm4
        subps  nb334_jzH1(%esp),%xmm5
        movaps %xmm0,nb334_dxH2M(%esp)
        movaps %xmm1,nb334_dyH2M(%esp)
        movaps %xmm2,nb334_dzH2M(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb334_dxMH1(%esp)
        movaps %xmm4,nb334_dyMH1(%esp)
        movaps %xmm5,nb334_dzMH1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb334_rsqH2M(%esp)
        movaps %xmm4,nb334_rsqMH1(%esp)

        movaps nb334_ixM(%esp),%xmm0
        movaps nb334_iyM(%esp),%xmm1
        movaps nb334_izM(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb334_jxH2(%esp),%xmm0
        subps  nb334_jyH2(%esp),%xmm1
        subps  nb334_jzH2(%esp),%xmm2
        subps  nb334_jxM(%esp),%xmm3
        subps  nb334_jyM(%esp),%xmm4
        subps  nb334_jzM(%esp),%xmm5
        movaps %xmm0,nb334_dxMH2(%esp)
        movaps %xmm1,nb334_dyMH2(%esp)
        movaps %xmm2,nb334_dzMH2(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb334_dxMM(%esp)
        movaps %xmm4,nb334_dyMM(%esp)
        movaps %xmm5,nb334_dzMM(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb334_rsqMH2(%esp)
        movaps %xmm4,nb334_rsqMM(%esp)

        ## Invsqrt for O-O
        rsqrtps  nb334_rsqOO(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb334_three(%esp),%xmm3
        mulps   nb334_rsqOO(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb334_half(%esp),%xmm3   ## rinvOO
        movaps %xmm3,nb334_rinvOO(%esp)

        ## Invsqrt for H1-H1 and H1-H2
        rsqrtps nb334_rsqH1H1(%esp),%xmm1
        rsqrtps nb334_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb334_rsqH1H1(%esp),%xmm1
        mulps   nb334_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334_half(%esp),%xmm3   ## rinvH1H1 
        mulps   nb334_half(%esp),%xmm7   ## rinvH1H2 
        movaps  %xmm3,nb334_rinvH1H1(%esp)
        movaps  %xmm7,nb334_rinvH1H2(%esp)

        rsqrtps nb334_rsqH1M(%esp),%xmm1
        rsqrtps nb334_rsqH2H1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb334_rsqH1M(%esp),%xmm1
        mulps   nb334_rsqH2H1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334_half(%esp),%xmm3
        mulps   nb334_half(%esp),%xmm7
        movaps  %xmm3,nb334_rinvH1M(%esp)
        movaps  %xmm7,nb334_rinvH2H1(%esp)

        rsqrtps nb334_rsqH2H2(%esp),%xmm1
        rsqrtps nb334_rsqH2M(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb334_rsqH2H2(%esp),%xmm1
        mulps   nb334_rsqH2M(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334_half(%esp),%xmm3
        mulps   nb334_half(%esp),%xmm7
        movaps  %xmm3,nb334_rinvH2H2(%esp)
        movaps  %xmm7,nb334_rinvH2M(%esp)

        rsqrtps nb334_rsqMH1(%esp),%xmm1
        rsqrtps nb334_rsqMH2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb334_rsqMH1(%esp),%xmm1
        mulps   nb334_rsqMH2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334_half(%esp),%xmm3
        mulps   nb334_half(%esp),%xmm7
        movaps  %xmm3,nb334_rinvMH1(%esp)
        movaps  %xmm7,nb334_rinvMH2(%esp)

        rsqrtps nb334_rsqMM(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb334_three(%esp),%xmm3
        mulps   nb334_rsqMM(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb334_half(%esp),%xmm3
        movaps  %xmm3,nb334_rinvMM(%esp)

        ## start with OO table interaction
        movaps nb334_rinvOO(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334_rsqOO(%esp),%xmm1   ## xmm1=r
        mulps  nb334_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movl nb334_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        ## load dispersion table data into xmm4-xmm7
        movlps 16(%esi,%eax,4),%xmm5
        movlps 16(%esi,%ecx,4),%xmm7
        movhps 16(%esi,%ebx,4),%xmm5
        movhps 16(%esi,%edx,4),%xmm7    ## got half table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 24(%esi,%eax,4),%xmm7
        movlps 24(%esi,%ecx,4),%xmm3
        movhps 24(%esi,%ebx,4),%xmm7
        movhps 24(%esi,%edx,4),%xmm3    ## other half of table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## dispersion table YFGH ready in xmm4-xmm7
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5      ## F+Geps 
        addps  %xmm7,%xmm5      ## xmm5=Fp=F+Geps+Heps2 
        mulps  nb334_two(%esp),%xmm7            ## 2*Heps2 
        addps  %xmm6,%xmm7      ## Geps+2*Heps2
        addps  %xmm5,%xmm7 ## xmm7=FF = Fp+Geps+2*Heps2
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV = Y+eps*Fp

        movaps nb334_c6(%esp),%xmm4
        mulps  %xmm4,%xmm7      ## fijD 
        mulps  %xmm4,%xmm5      ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addps  nb334_Vvdwtot(%esp),%xmm5
        movaps %xmm7,nb334_fstmp(%esp)   ## fscal 
        movaps %xmm5,nb334_Vvdwtot(%esp)

        ## load repulsion table data into xmm4-xmm7
        movlps 32(%esi,%eax,4),%xmm5
        movlps 32(%esi,%ecx,4),%xmm7
        movhps 32(%esi,%ebx,4),%xmm5
        movhps 32(%esi,%edx,4),%xmm7    ## got half table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 40(%esi,%eax,4),%xmm7
        movlps 40(%esi,%ecx,4),%xmm3
        movhps 40(%esi,%ebx,4),%xmm7
        movhps 40(%esi,%edx,4),%xmm3    ## other half of table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 #
        ## repulsion table YFGH ready in xmm4-xmm7# 11011101

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb334_two(%esp),%xmm7            ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb334_c12(%esp),%xmm4
        mulps  %xmm4,%xmm7 ## fijR 
        mulps  %xmm4,%xmm5 ## Vvdw12 
        addps  nb334_fstmp(%esp),%xmm7

        addps  nb334_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb334_Vvdwtot(%esp)

        xorps %xmm1,%xmm1
        mulps nb334_tsc(%esp),%xmm7
        mulps %xmm0,%xmm7
        subps  %xmm7,%xmm1


        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb334_dxOO(%esp),%xmm0
        mulps nb334_dyOO(%esp),%xmm1
        mulps nb334_dzOO(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb334_fixO(%esp),%xmm0
        addps nb334_fiyO(%esp),%xmm1
        addps nb334_fizO(%esp),%xmm2
        movaps %xmm3,nb334_fjxO(%esp)
        movaps %xmm4,nb334_fjyO(%esp)
        movaps %xmm5,nb334_fjzO(%esp)
        movaps %xmm0,nb334_fixO(%esp)
        movaps %xmm1,nb334_fiyO(%esp)
        movaps %xmm2,nb334_fizO(%esp)

        ## Coulomb interactions - first H1H1
        movaps nb334_rinvH1H1(%esp),%xmm0

        movaps %xmm0,%xmm1
        mulps  nb334_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulps  nb334_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        movl nb334_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb334_two(%esp),%xmm7            ## two*Heps2 
        movaps nb334_qqHH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 
        ## update vctot 
        addps  nb334_vctot(%esp),%xmm5
        movaps %xmm5,nb334_vctot(%esp)

        xorps  %xmm1,%xmm1
        mulps  nb334_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5

        mulps nb334_dxH1H1(%esp),%xmm0
        mulps nb334_dyH1H1(%esp),%xmm1
        mulps nb334_dzH1H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb334_fixH1(%esp),%xmm0
        addps nb334_fiyH1(%esp),%xmm1
        addps nb334_fizH1(%esp),%xmm2
        movaps %xmm3,nb334_fjxH1(%esp)
        movaps %xmm4,nb334_fjyH1(%esp)
        movaps %xmm5,nb334_fjzH1(%esp)
        movaps %xmm0,nb334_fixH1(%esp)
        movaps %xmm1,nb334_fiyH1(%esp)
        movaps %xmm2,nb334_fizH1(%esp)

        ## H1-H2 interaction 
        movaps nb334_rinvH1H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulps  nb334_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb334_two(%esp),%xmm7            ## two*Heps2 
        movaps nb334_qqHH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb334_vctot(%esp),%xmm5
        movaps %xmm5,nb334_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb334_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb334_dxH1H2(%esp),%xmm0
        mulps nb334_dyH1H2(%esp),%xmm1
        mulps nb334_dzH1H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb334_fixH1(%esp),%xmm0
        addps nb334_fiyH1(%esp),%xmm1
        addps nb334_fizH1(%esp),%xmm2
        movaps %xmm3,nb334_fjxH2(%esp)
        movaps %xmm4,nb334_fjyH2(%esp)
        movaps %xmm5,nb334_fjzH2(%esp)
        movaps %xmm0,nb334_fixH1(%esp)
        movaps %xmm1,nb334_fiyH1(%esp)
        movaps %xmm2,nb334_fizH1(%esp)

        ## H1-M interaction  
        movaps nb334_rinvH1M(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334_rsqH1M(%esp),%xmm1   ## xmm1=r 
        mulps  nb334_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx


        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb334_two(%esp),%xmm7            ## two*Heps2 
        movaps nb334_qqMH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb334_vctot(%esp),%xmm5
        movaps %xmm5,nb334_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb334_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb334_dxH1M(%esp),%xmm0
        mulps nb334_dyH1M(%esp),%xmm1
        mulps nb334_dzH1M(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb334_fixH1(%esp),%xmm0
        addps nb334_fiyH1(%esp),%xmm1
        addps nb334_fizH1(%esp),%xmm2
        movaps %xmm3,nb334_fjxM(%esp)
        movaps %xmm4,nb334_fjyM(%esp)
        movaps %xmm5,nb334_fjzM(%esp)
        movaps %xmm0,nb334_fixH1(%esp)
        movaps %xmm1,nb334_fiyH1(%esp)
        movaps %xmm2,nb334_fizH1(%esp)

        ## H2-H1 interaction 
        movaps nb334_rinvH2H1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulps  nb334_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb334_two(%esp),%xmm7            ## two*Heps2 
        movaps nb334_qqHH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb334_vctot(%esp),%xmm5
        movaps %xmm5,nb334_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb334_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        movaps nb334_fjxH1(%esp),%xmm3
        movaps nb334_fjyH1(%esp),%xmm4
        movaps nb334_fjzH1(%esp),%xmm5
        mulps nb334_dxH2H1(%esp),%xmm0
        mulps nb334_dyH2H1(%esp),%xmm1
        mulps nb334_dzH2H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb334_fixH2(%esp),%xmm0
        addps nb334_fiyH2(%esp),%xmm1
        addps nb334_fizH2(%esp),%xmm2
        movaps %xmm3,nb334_fjxH1(%esp)
        movaps %xmm4,nb334_fjyH1(%esp)
        movaps %xmm5,nb334_fjzH1(%esp)
        movaps %xmm0,nb334_fixH2(%esp)
        movaps %xmm1,nb334_fiyH2(%esp)
        movaps %xmm2,nb334_fizH2(%esp)

        ## H2-H2 interaction 
        movaps nb334_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulps  nb334_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb334_two(%esp),%xmm7            ## two*Heps2 
        movaps nb334_qqHH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb334_vctot(%esp),%xmm5
        movaps %xmm5,nb334_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb334_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        movaps nb334_fjxH2(%esp),%xmm3
        movaps nb334_fjyH2(%esp),%xmm4
        movaps nb334_fjzH2(%esp),%xmm5
        mulps nb334_dxH2H2(%esp),%xmm0
        mulps nb334_dyH2H2(%esp),%xmm1
        mulps nb334_dzH2H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb334_fixH2(%esp),%xmm0
        addps nb334_fiyH2(%esp),%xmm1
        addps nb334_fizH2(%esp),%xmm2
        movaps %xmm3,nb334_fjxH2(%esp)
        movaps %xmm4,nb334_fjyH2(%esp)
        movaps %xmm5,nb334_fjzH2(%esp)
        movaps %xmm0,nb334_fixH2(%esp)
        movaps %xmm1,nb334_fiyH2(%esp)
        movaps %xmm2,nb334_fizH2(%esp)

        ## H2-M interaction 
        movaps nb334_rinvH2M(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334_rsqH2M(%esp),%xmm1   ## xmm1=r 
        mulps  nb334_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb334_two(%esp),%xmm7            ## two*Heps2 
        movaps nb334_qqMH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb334_vctot(%esp),%xmm5
        movaps %xmm5,nb334_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb334_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        movaps nb334_fjxM(%esp),%xmm3
        movaps nb334_fjyM(%esp),%xmm4
        movaps nb334_fjzM(%esp),%xmm5
        mulps nb334_dxH2M(%esp),%xmm0
        mulps nb334_dyH2M(%esp),%xmm1
        mulps nb334_dzH2M(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb334_fixH2(%esp),%xmm0
        addps nb334_fiyH2(%esp),%xmm1
        addps nb334_fizH2(%esp),%xmm2
        movaps %xmm3,nb334_fjxM(%esp)
        movaps %xmm4,nb334_fjyM(%esp)
        movaps %xmm5,nb334_fjzM(%esp)
        movaps %xmm0,nb334_fixH2(%esp)
        movaps %xmm1,nb334_fiyH2(%esp)
        movaps %xmm2,nb334_fizH2(%esp)

        ## M-H1 interaction 
        movaps nb334_rinvMH1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334_rsqMH1(%esp),%xmm1   ## xmm1=r 
        mulps  nb334_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb334_two(%esp),%xmm7            ## two*Heps2 
        movaps nb334_qqMH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb334_vctot(%esp),%xmm5
        movaps %xmm5,nb334_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb334_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        movaps nb334_fjxH1(%esp),%xmm3
        movaps nb334_fjyH1(%esp),%xmm4
        movaps nb334_fjzH1(%esp),%xmm5
        mulps nb334_dxMH1(%esp),%xmm0
        mulps nb334_dyMH1(%esp),%xmm1
        mulps nb334_dzMH1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb334_fixM(%esp),%xmm0
        addps nb334_fiyM(%esp),%xmm1
        addps nb334_fizM(%esp),%xmm2
        movaps %xmm3,nb334_fjxH1(%esp)
        movaps %xmm4,nb334_fjyH1(%esp)
        movaps %xmm5,nb334_fjzH1(%esp)
        movaps %xmm0,nb334_fixM(%esp)
        movaps %xmm1,nb334_fiyM(%esp)
        movaps %xmm2,nb334_fizM(%esp)

        ## M-H2 interaction 
        movaps nb334_rinvMH2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334_rsqMH2(%esp),%xmm1   ## xmm1=r 
        mulps  nb334_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb334_two(%esp),%xmm7            ## two*Heps2 
        movaps nb334_qqMH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb334_vctot(%esp),%xmm5
        movaps %xmm5,nb334_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb334_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        movaps nb334_fjxH2(%esp),%xmm3
        movaps nb334_fjyH2(%esp),%xmm4
        movaps nb334_fjzH2(%esp),%xmm5
        mulps nb334_dxMH2(%esp),%xmm0
        mulps nb334_dyMH2(%esp),%xmm1
        mulps nb334_dzMH2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb334_fixM(%esp),%xmm0
        addps nb334_fiyM(%esp),%xmm1
        addps nb334_fizM(%esp),%xmm2
        movaps %xmm3,nb334_fjxH2(%esp)
        movaps %xmm4,nb334_fjyH2(%esp)
        movaps %xmm5,nb334_fjzH2(%esp)
        movaps %xmm0,nb334_fixM(%esp)
        movaps %xmm1,nb334_fiyM(%esp)
        movaps %xmm2,nb334_fizM(%esp)

        ## M-M interaction 
        movaps nb334_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334_rsqMM(%esp),%xmm1   ## xmm1=r 
        mulps  nb334_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb334_two(%esp),%xmm7            ## two*Heps2 
        movaps nb334_qqMM(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb334_vctot(%esp),%xmm5
        movaps %xmm5,nb334_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb334_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        movaps nb334_fjxM(%esp),%xmm3
        movaps nb334_fjyM(%esp),%xmm4
        movaps nb334_fjzM(%esp),%xmm5
        mulps nb334_dxMM(%esp),%xmm0
        mulps nb334_dyMM(%esp),%xmm1
        mulps nb334_dzMM(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb334_fixM(%esp),%xmm0
        addps nb334_fiyM(%esp),%xmm1
        addps nb334_fizM(%esp),%xmm2
        movaps %xmm3,nb334_fjxM(%esp)
        movaps %xmm4,nb334_fjyM(%esp)
        movaps %xmm5,nb334_fjzM(%esp)
        movaps %xmm0,nb334_fixM(%esp)
        movaps %xmm1,nb334_fiyM(%esp)
        movaps %xmm2,nb334_fizM(%esp)

        movl nb334_faction(%ebp),%edi

        movd %mm0,%eax
        movd %mm1,%ebx
        movd %mm2,%ecx
        movd %mm3,%edx

        ## Did all interactions - now update j forces 
        ## 4 j waters with four atoms each.
        ## step 1 : transpose fjxO, fjyO, fjzO, fjxH1
        movaps nb334_fjxO(%esp),%xmm0
        movaps nb334_fjyO(%esp),%xmm1
        movaps nb334_fjzO(%esp),%xmm2
        movaps nb334_fjxH1(%esp),%xmm3
        movaps %xmm0,%xmm4
        movaps %xmm1,%xmm5
        unpcklps %xmm2,%xmm4
        unpcklps %xmm3,%xmm5
        unpckhps %xmm2,%xmm0
        unpckhps %xmm3,%xmm1
        movaps %xmm4,%xmm2
        movaps %xmm0,%xmm3
        unpcklps %xmm5,%xmm4
        unpckhps %xmm5,%xmm2

        unpcklps %xmm1,%xmm0
        unpckhps %xmm1,%xmm3
        ## results are now in xmm4, xmm2, xmm0, xmm3
        ## load the corresponding j forces from memory
        movlps   (%edi,%eax,4),%xmm1
        movlps   (%edi,%ebx,4),%xmm5
        movlps   (%edi,%ecx,4),%xmm6
        movlps   (%edi,%edx,4),%xmm7
        movhps   8(%edi,%eax,4),%xmm1
        movhps   8(%edi,%ebx,4),%xmm5
        movhps   8(%edi,%ecx,4),%xmm6
        movhps   8(%edi,%edx,4),%xmm7
        ## add
        addps    %xmm4,%xmm1
        addps    %xmm2,%xmm5
        addps    %xmm0,%xmm6
        addps    %xmm3,%xmm7
        ## store back
        movlps   %xmm1,(%edi,%eax,4)
        movlps   %xmm5,(%edi,%ebx,4)
        movlps   %xmm6,(%edi,%ecx,4)
        movlps   %xmm7,(%edi,%edx,4)
        movhps   %xmm1,8(%edi,%eax,4)
        movhps   %xmm5,8(%edi,%ebx,4)
        movhps   %xmm6,8(%edi,%ecx,4)
        movhps   %xmm7,8(%edi,%edx,4)

        ## step 2 : transpose fjyH1, fjzH1, fjxH2, fjyH2
        movaps nb334_fjyH1(%esp),%xmm0
        movaps nb334_fjzH1(%esp),%xmm1
        movaps nb334_fjxH2(%esp),%xmm2
        movaps nb334_fjyH2(%esp),%xmm3
        movaps %xmm0,%xmm4
        movaps %xmm1,%xmm5
        unpcklps %xmm2,%xmm4
        unpcklps %xmm3,%xmm5
        unpckhps %xmm2,%xmm0
        unpckhps %xmm3,%xmm1
        movaps %xmm4,%xmm2
        movaps %xmm0,%xmm3
        unpcklps %xmm5,%xmm4
        unpckhps %xmm5,%xmm2

        unpcklps %xmm1,%xmm0
        unpckhps %xmm1,%xmm3
        ## results are now in xmm4, xmm2, xmm0, xmm3
        ## load the corresponding j forces from memory
        movlps   16(%edi,%eax,4),%xmm1
        movlps   16(%edi,%ebx,4),%xmm5
        movlps   16(%edi,%ecx,4),%xmm6
        movlps   16(%edi,%edx,4),%xmm7
        movhps   24(%edi,%eax,4),%xmm1
        movhps   24(%edi,%ebx,4),%xmm5
        movhps   24(%edi,%ecx,4),%xmm6
        movhps   24(%edi,%edx,4),%xmm7
        ## add
        addps    %xmm4,%xmm1
        addps    %xmm2,%xmm5
        addps    %xmm0,%xmm6
        addps    %xmm3,%xmm7
        ## store back
        movlps   %xmm1,16(%edi,%eax,4)
        movlps   %xmm5,16(%edi,%ebx,4)
        movlps   %xmm6,16(%edi,%ecx,4)
        movlps   %xmm7,16(%edi,%edx,4)
        movhps   %xmm1,24(%edi,%eax,4)
        movhps   %xmm5,24(%edi,%ebx,4)
        movhps   %xmm6,24(%edi,%ecx,4)
        movhps   %xmm7,24(%edi,%edx,4)

        ## step 3 : transpose fjzH2, fjxM, fjyM, fjzM
        movaps nb334_fjzH2(%esp),%xmm0
        movaps nb334_fjxM(%esp),%xmm1
        movaps nb334_fjyM(%esp),%xmm2
        movaps nb334_fjzM(%esp),%xmm3
        movaps %xmm0,%xmm4
        movaps %xmm1,%xmm5
        unpcklps %xmm2,%xmm4
        unpcklps %xmm3,%xmm5
        unpckhps %xmm2,%xmm0
        unpckhps %xmm3,%xmm1
        movaps %xmm4,%xmm2
        movaps %xmm0,%xmm3
        unpcklps %xmm5,%xmm4
        unpckhps %xmm5,%xmm2

        unpcklps %xmm1,%xmm0
        unpckhps %xmm1,%xmm3
        ## results are now in xmm4, xmm2, xmm0, xmm3
        ## load the corresponding j forces from memory
        movlps   32(%edi,%eax,4),%xmm1
        movlps   32(%edi,%ebx,4),%xmm5
        movlps   32(%edi,%ecx,4),%xmm6
        movlps   32(%edi,%edx,4),%xmm7
        movhps   40(%edi,%eax,4),%xmm1
        movhps   40(%edi,%ebx,4),%xmm5
        movhps   40(%edi,%ecx,4),%xmm6
        movhps   40(%edi,%edx,4),%xmm7
        ## add
        addps    %xmm4,%xmm1
        addps    %xmm2,%xmm5
        addps    %xmm0,%xmm6
        addps    %xmm3,%xmm7
        ## store back
        movlps   %xmm1,32(%edi,%eax,4)
        movlps   %xmm5,32(%edi,%ebx,4)
        movlps   %xmm6,32(%edi,%ecx,4)
        movlps   %xmm7,32(%edi,%edx,4)
        movhps   %xmm1,40(%edi,%eax,4)
        movhps   %xmm5,40(%edi,%ebx,4)
        movhps   %xmm6,40(%edi,%ecx,4)
        movhps   %xmm7,40(%edi,%edx,4)

        ## should we do one more iteration? 
        subl $4,nb334_innerk(%esp)
        jl    _nb_kernel334_ia32_sse.nb334_single_check
        jmp   _nb_kernel334_ia32_sse.nb334_unroll_loop
_nb_kernel334_ia32_sse.nb334_single_check: 
        addl $4,nb334_innerk(%esp)
        jnz   _nb_kernel334_ia32_sse.nb334_single_loop
        jmp   _nb_kernel334_ia32_sse.nb334_updateouterdata
_nb_kernel334_ia32_sse.nb334_single_loop: 
        movl  nb334_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb334_innerjjnr(%esp)

        movl nb334_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates
        movlps (%esi,%eax,4),%xmm3              ##  Ox  Oy  
        movlps 16(%esi,%eax,4),%xmm4            ## H1y H1z 
        movlps 32(%esi,%eax,4),%xmm5            ## H2z  Mx 
        movhps 8(%esi,%eax,4),%xmm3             ##  Ox  Oy  Oz H1x
        movhps 24(%esi,%eax,4),%xmm4            ## H1y H1z H2x H2y
        movhps 40(%esi,%eax,4),%xmm5            ## H2z  Mx  My  Mz
        ## transpose
        movaps %xmm4,%xmm0
        movaps %xmm3,%xmm1
        movaps %xmm4,%xmm2
        movaps %xmm3,%xmm6
        shufps $18,%xmm5,%xmm4 ## (00010010)  h2x - Mx  - 
        shufps $193,%xmm0,%xmm3 ## (11000001)  Oy  - H1y - 
        shufps $35,%xmm5,%xmm2 ## (00100011) H2y - My  - 
        shufps $18,%xmm0,%xmm1 ## (00010010)  Oz  - H1z - 
        ##  xmm6: Ox - - H1x   xmm5: H2z - - Mz 
        shufps $140,%xmm4,%xmm6 ## (10001100) Ox H1x H2x Mx 
        shufps $136,%xmm2,%xmm3 ## (10001000) Oy H1y H2y My 
        shufps $200,%xmm5,%xmm1 ## (11001000) Oz H1z H2z Mz

        ## store all j coordinates in jO  
        movaps %xmm6,nb334_jxO(%esp)
        movaps %xmm3,nb334_jyO(%esp)
        movaps %xmm1,nb334_jzO(%esp)

        ## do O and H1 in parallel
        movaps nb334_ixO(%esp),%xmm0
        movaps nb334_iyO(%esp),%xmm1
        movaps nb334_izO(%esp),%xmm2
        movaps nb334_ixH1(%esp),%xmm3
        movaps nb334_iyH1(%esp),%xmm4
        movaps nb334_izH1(%esp),%xmm5
        subps  nb334_jxO(%esp),%xmm0
        subps  nb334_jyO(%esp),%xmm1
        subps  nb334_jzO(%esp),%xmm2
        subps  nb334_jxO(%esp),%xmm3
        subps  nb334_jyO(%esp),%xmm4
        subps  nb334_jzO(%esp),%xmm5

        movaps %xmm0,nb334_dxOO(%esp)
        movaps %xmm1,nb334_dyOO(%esp)
        movaps %xmm2,nb334_dzOO(%esp)
        movaps %xmm3,nb334_dxH1H1(%esp)
        movaps %xmm4,nb334_dyH1H1(%esp)
        movaps %xmm5,nb334_dzH1H1(%esp)

        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm3,%xmm4
        addps %xmm5,%xmm4       ## have rsq in xmm4
        movaps %xmm0,nb334_rsqOO(%esp)
        movaps %xmm4,nb334_rsqH1H1(%esp)

        ## do 1/sqrt(x) for O and  H1
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334_half(%esp),%xmm3   ## rinv O - j water 
        mulps   nb334_half(%esp),%xmm7   ## rinv H1 - j water  

        movaps %xmm3,nb334_rinvOO(%esp)
        movaps %xmm7,nb334_rinvH1H1(%esp)

        movl nb334_VFtab(%ebp),%esi

        ## do O table LJ interaction
        movaps %xmm3,%xmm0
        movaps %xmm0,%xmm1
        mulss  nb334_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulss  nb334_tsc(%esp),%xmm1

        cvttps2pi %xmm1,%mm6
        cvtpi2ps %mm6,%xmm3
        subss    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6

        movd %mm6,%ebx
        leal  (%ebx,%ebx,2),%ebx

        ## load dispersion table data into xmm4
        movlps 16(%esi,%ebx,4),%xmm4
        movlps 24(%esi,%ebx,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $0x1,%xmm5,%xmm5
        shufps $0x1,%xmm7,%xmm7
        ## dispersion table YFGH ready in xmm4-xmm7
        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp 
        mulss  nb334_two(%esp),%xmm7            ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb334_c6(%esp),%xmm4
        mulss  %xmm4,%xmm7      ## fijD 
        mulss  %xmm4,%xmm5      ## Vvdw6 

        ## save scalar force in xmm3. Update Vvdwtot directly 
        addss  nb334_Vvdwtot(%esp),%xmm5
        movaps %xmm7,%xmm3 ## fscal 
        movss %xmm5,nb334_Vvdwtot(%esp)

        ## load repulsion table data into xmm4
        movlps 32(%esi,%ebx,4),%xmm4
        movlps 40(%esi,%ebx,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $0x1,%xmm5,%xmm5
        shufps $0x1,%xmm7,%xmm7
        ## repulsion table YFGH ready in xmm4-xmm7

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp 
        mulss  nb334_two(%esp),%xmm7            ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb334_c12(%esp),%xmm4
        mulss  %xmm4,%xmm7 ## fijR 
        mulss  %xmm4,%xmm5 ## Vvdw12 
        addss  %xmm3,%xmm7

        addss  nb334_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb334_Vvdwtot(%esp)

        xorps  %xmm1,%xmm1
        mulss nb334_tsc(%esp),%xmm7
        mulss %xmm0,%xmm7
        subss  %xmm7,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        mulss  nb334_dxOO(%esp),%xmm0
        mulss  nb334_dyOO(%esp),%xmm1
        mulss  nb334_dzOO(%esp),%xmm2
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        subss   %xmm0,%xmm3
        subss   %xmm1,%xmm4
        subss   %xmm2,%xmm5
        movaps  %xmm3,nb334_fjxO(%esp)
        movaps  %xmm4,nb334_fjyO(%esp)
        movaps  %xmm5,nb334_fjzO(%esp)
        addss   nb334_fixO(%esp),%xmm0
        addss   nb334_fiyO(%esp),%xmm1
        addss   nb334_fizO(%esp),%xmm2
        movss  %xmm0,nb334_fixO(%esp)
        movss  %xmm1,nb334_fiyO(%esp)
        movss  %xmm2,nb334_fizO(%esp)

        ## do  H1 coulomb interaction
        movaps nb334_rinvH1H1(%esp),%xmm0   ## rinv 
        movaps %xmm0,%xmm1
        mulps  nb334_rsqH1H1(%esp),%xmm1        ## r
        mulps nb334_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movl nb334_VFtab(%ebp),%esi

        movlps (%esi,%ebx,4),%xmm4
        movlps (%esi,%ecx,4),%xmm3
        movlps (%esi,%edx,4),%xmm7
        movhps 8(%esi,%ebx,4),%xmm4
        movhps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm7
        movaps %xmm3,%xmm6
        unpcklps %xmm7,%xmm6
        unpckhps %xmm7,%xmm3
        movaps %xmm4,%xmm5
        movaps %xmm4,%xmm7
        shufps $0x40,%xmm6,%xmm4
        shufps $0xE4,%xmm6,%xmm5
        movaps %xmm7,%xmm6
        shufps $0x48,%xmm3,%xmm6
        shufps $0xEC,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb334_two(%esp),%xmm7            ## two*Heps2 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb334_qqHH(%esp),%xmm3
        movhps  nb334_qqMH(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001 

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC 

        addps  nb334_vctot(%esp),%xmm5
        movaps %xmm5,nb334_vctot(%esp)

        mulps  nb334_tsc(%esp),%xmm3
        xorps  %xmm2,%xmm2
        subps  %xmm3,%xmm2
        mulps  %xmm2,%xmm0
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb334_dxH1H1(%esp),%xmm0
        mulps   nb334_dyH1H1(%esp),%xmm1
        mulps   nb334_dzH1H1(%esp),%xmm2
        ## update forces H1 - j water 
        movaps  nb334_fjxO(%esp),%xmm3
        movaps  nb334_fjyO(%esp),%xmm4
        movaps  nb334_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb334_fjxO(%esp)
        movaps  %xmm4,nb334_fjyO(%esp)
        movaps  %xmm5,nb334_fjzO(%esp)
        addps   nb334_fixH1(%esp),%xmm0
        addps   nb334_fiyH1(%esp),%xmm1
        addps   nb334_fizH1(%esp),%xmm2
        movaps  %xmm0,nb334_fixH1(%esp)
        movaps  %xmm1,nb334_fiyH1(%esp)
        movaps  %xmm2,nb334_fizH1(%esp)

        ## i H2 & M simultaneously first get i particle coords: 
        movaps  nb334_ixH2(%esp),%xmm0
        movaps  nb334_iyH2(%esp),%xmm1
        movaps  nb334_izH2(%esp),%xmm2
        movaps  nb334_ixM(%esp),%xmm3
        movaps  nb334_iyM(%esp),%xmm4
        movaps  nb334_izM(%esp),%xmm5
        subps   nb334_jxO(%esp),%xmm0
        subps   nb334_jyO(%esp),%xmm1
        subps   nb334_jzO(%esp),%xmm2
        subps   nb334_jxO(%esp),%xmm3
        subps   nb334_jyO(%esp),%xmm4
        subps   nb334_jzO(%esp),%xmm5
        movaps %xmm0,nb334_dxH2H2(%esp)
        movaps %xmm1,nb334_dyH2H2(%esp)
        movaps %xmm2,nb334_dzH2H2(%esp)
        movaps %xmm3,nb334_dxMM(%esp)
        movaps %xmm4,nb334_dyMM(%esp)
        movaps %xmm5,nb334_dzMM(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm1,%xmm0
        addps %xmm3,%xmm4
        addps %xmm2,%xmm0       ## have rsqH2 in xmm0 
        addps %xmm5,%xmm4       ## have rsqM in xmm4 

        ## start with H2, save data 
        movaps %xmm0,nb334_rsqH2H2(%esp)
        movaps %xmm4,nb334_rsqMM(%esp)
        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334_half(%esp),%xmm3   ## rinv H2 - j water 
        mulps   nb334_half(%esp),%xmm7   ## rinv M - j water  

        movaps %xmm3,nb334_rinvH2H2(%esp)
        movaps %xmm7,nb334_rinvMM(%esp)

        movaps %xmm3,%xmm1
        mulps  nb334_rsqH2H2(%esp),%xmm1        ## xmm1=r 
        movaps %xmm3,%xmm0      ## xmm0=rinv 
        mulps  nb334_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%ebx,4),%xmm4
        movlps (%esi,%ecx,4),%xmm3
        movlps (%esi,%edx,4),%xmm7
        movhps 8(%esi,%ebx,4),%xmm4
        movhps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm7
        movaps %xmm3,%xmm6
        unpcklps %xmm7,%xmm6
        unpckhps %xmm7,%xmm3
        movaps %xmm4,%xmm5
        movaps %xmm4,%xmm7
        shufps $0x40,%xmm6,%xmm4
        shufps $0xE4,%xmm6,%xmm5
        movaps %xmm7,%xmm6
        shufps $0x48,%xmm3,%xmm6
        shufps $0xEC,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb334_two(%esp),%xmm7            ## two*Heps2 

        xorps  %xmm3,%xmm3

        ## fetch charges to xmm3 (temporary) 
        movss   nb334_qqHH(%esp),%xmm3
        movhps  nb334_qqMH(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC 
        addps  nb334_vctot(%esp),%xmm5
        movaps %xmm5,nb334_vctot(%esp)

        xorps  %xmm1,%xmm1

        mulps nb334_tsc(%esp),%xmm3
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps  %xmm1,%xmm0
        movaps  %xmm1,%xmm2
        mulps   nb334_dxH2H2(%esp),%xmm0
        mulps   nb334_dyH2H2(%esp),%xmm1
        mulps   nb334_dzH2H2(%esp),%xmm2
        ## update forces H1 - j water 
        movaps  nb334_fjxO(%esp),%xmm3
        movaps  nb334_fjyO(%esp),%xmm4
        movaps  nb334_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb334_fjxO(%esp)
        movaps  %xmm4,nb334_fjyO(%esp)
        movaps  %xmm5,nb334_fjzO(%esp)
        addps   nb334_fixH2(%esp),%xmm0
        addps   nb334_fiyH2(%esp),%xmm1
        addps   nb334_fizH2(%esp),%xmm2
        movaps  %xmm0,nb334_fixH2(%esp)
        movaps  %xmm1,nb334_fiyH2(%esp)
        movaps  %xmm2,nb334_fizH2(%esp)

        ## do table for i M - j water interaction 
        movaps nb334_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334_rsqMM(%esp),%xmm1          ## xmm0=rinv, xmm1=r 
        mulps  nb334_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movl nb334_VFtab(%ebp),%esi

        movlps (%esi,%ebx,4),%xmm4
        movlps (%esi,%ecx,4),%xmm3
        movlps (%esi,%edx,4),%xmm7
        movhps 8(%esi,%ebx,4),%xmm4
        movhps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm7
        movaps %xmm3,%xmm6
        unpcklps %xmm7,%xmm6
        unpckhps %xmm7,%xmm3
        movaps %xmm4,%xmm5
        movaps %xmm4,%xmm7
        shufps $0x40,%xmm6,%xmm4
        shufps $0xE4,%xmm6,%xmm5
        movaps %xmm7,%xmm6
        shufps $0x48,%xmm3,%xmm6
        shufps $0xEC,%xmm3,%xmm7
        ## # coulomb table ready, in xmm4-xmm7


        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb334_two(%esp),%xmm7            ## two*Heps2 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb334_qqMH(%esp),%xmm3
        movhps  nb334_qqMM(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC 
        addps  nb334_vctot(%esp),%xmm5
        movaps %xmm5,nb334_vctot(%esp)

        xorps  %xmm1,%xmm1

        mulps nb334_tsc(%esp),%xmm3
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps  %xmm1,%xmm0
        movaps  %xmm1,%xmm2

        mulps   nb334_dxMM(%esp),%xmm0
        mulps   nb334_dyMM(%esp),%xmm1
        mulps   nb334_dzMM(%esp),%xmm2
        movaps  nb334_fjxO(%esp),%xmm3
        movaps  nb334_fjyO(%esp),%xmm4
        movaps  nb334_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movl    nb334_faction(%ebp),%esi
        movaps  %xmm3,nb334_fjxO(%esp)
        movaps  %xmm4,nb334_fjyO(%esp)
        movaps  %xmm5,nb334_fjzO(%esp)
        addps   nb334_fixM(%esp),%xmm0
        addps   nb334_fiyM(%esp),%xmm1
        addps   nb334_fizM(%esp),%xmm2
        movaps  %xmm0,nb334_fixM(%esp)
        movaps  %xmm1,nb334_fiyM(%esp)
        movaps  %xmm2,nb334_fizM(%esp)

        ## update j water forces from local variables.
        ## transpose back first
        movaps  nb334_fjxO(%esp),%xmm0   ## Ox H1x H2x Mx 
        movaps  nb334_fjyO(%esp),%xmm1   ## Oy H1y H2y My
        movaps  nb334_fjzO(%esp),%xmm2   ## Oz H1z H2z Mz

        movaps  %xmm0,%xmm3
        movaps  %xmm0,%xmm4
        unpcklps %xmm1,%xmm3            ## Ox Oy - -
        shufps $0x1,%xmm2,%xmm4        ## h1x - Oz -
        movaps  %xmm1,%xmm5
        movaps  %xmm0,%xmm6
        unpcklps %xmm2,%xmm5            ## - - H1y H1z
        unpckhps %xmm1,%xmm6            ## h2x h2y - - 
        unpckhps %xmm2,%xmm1            ## - - My Mz

        shufps  $0x32,%xmm0,%xmm2 ## (00110010) h2z - Mx -
        shufps  $36,%xmm4,%xmm3 ## constant 00100100 ;# Ox Oy Oz H1x 
        shufps  $78,%xmm6,%xmm5 ## constant 01001110 ;# h1y h1z h2x h2y
        shufps  $232,%xmm1,%xmm2 ## 11101000 ;# h2z mx my mz

        movlps  (%esi,%eax,4),%xmm0
        movlps  16(%esi,%eax,4),%xmm1
        movlps  32(%esi,%eax,4),%xmm4
        movhps  8(%esi,%eax,4),%xmm0
        movhps  24(%esi,%eax,4),%xmm1
        movhps  40(%esi,%eax,4),%xmm4
        addps   %xmm3,%xmm0
        addps   %xmm5,%xmm1
        addps   %xmm2,%xmm4
        movlps   %xmm0,(%esi,%eax,4)
        movlps   %xmm1,16(%esi,%eax,4)
        movlps   %xmm4,32(%esi,%eax,4)
        movhps   %xmm0,8(%esi,%eax,4)
        movhps   %xmm1,24(%esi,%eax,4)
        movhps   %xmm4,40(%esi,%eax,4)

        decl nb334_innerk(%esp)
        jz    _nb_kernel334_ia32_sse.nb334_updateouterdata
        jmp   _nb_kernel334_ia32_sse.nb334_single_loop
_nb_kernel334_ia32_sse.nb334_updateouterdata: 
        movl  nb334_ii3(%esp),%ecx
        movl  nb334_faction(%ebp),%edi
        movl  nb334_fshift(%ebp),%esi
        movl  nb334_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb334_fixO(%esp),%xmm0
        movaps nb334_fiyO(%esp),%xmm1
        movaps nb334_fizO(%esp),%xmm2

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
        movss  (%edi,%ecx,4),%xmm3
        movss  4(%edi,%ecx,4),%xmm4
        movss  8(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,(%edi,%ecx,4)
        movss  %xmm4,4(%edi,%ecx,4)
        movss  %xmm5,8(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        movaps %xmm0,%xmm6
        movss %xmm2,%xmm7
        movlhps %xmm1,%xmm6
        shufps $8,%xmm6,%xmm6 ## constant 00001000      

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movaps nb334_fixH1(%esp),%xmm0
        movaps nb334_fiyH1(%esp),%xmm1
        movaps nb334_fizH1(%esp),%xmm2

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
        movss  12(%edi,%ecx,4),%xmm3
        movss  16(%edi,%ecx,4),%xmm4
        movss  20(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,12(%edi,%ecx,4)
        movss  %xmm4,16(%edi,%ecx,4)
        movss  %xmm5,20(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movaps nb334_fixH2(%esp),%xmm0
        movaps nb334_fiyH2(%esp),%xmm1
        movaps nb334_fizH2(%esp),%xmm2

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
        movss  24(%edi,%ecx,4),%xmm3
        movss  28(%edi,%ecx,4),%xmm4
        movss  32(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,24(%edi,%ecx,4)
        movss  %xmm4,28(%edi,%ecx,4)
        movss  %xmm5,32(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movaps nb334_fixM(%esp),%xmm0
        movaps nb334_fiyM(%esp),%xmm1
        movaps nb334_fizM(%esp),%xmm2

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
        movss  36(%edi,%ecx,4),%xmm3
        movss  40(%edi,%ecx,4),%xmm4
        movss  44(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,36(%edi,%ecx,4)
        movss  %xmm4,40(%edi,%ecx,4)
        movss  %xmm5,44(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## increment fshift force  
        movlps  (%esi,%edx,4),%xmm3
        movss  8(%esi,%edx,4),%xmm4
        addps  %xmm6,%xmm3
        addss  %xmm7,%xmm4
        movlps  %xmm3,(%esi,%edx,4)
        movss  %xmm4,8(%esi,%edx,4)

        ## get n from stack
        movl nb334_n(%esp),%esi
        ## get group index for i particle 
        movl  nb334_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb334_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb334_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb334_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb334_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb334_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel334_ia32_sse.nb334_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb334_n(%esp)
        jmp _nb_kernel334_ia32_sse.nb334_outer
_nb_kernel334_ia32_sse.nb334_outerend: 
        ## check if more outer neighborlists remain
        movl  nb334_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel334_ia32_sse.nb334_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel334_ia32_sse.nb334_threadloop
_nb_kernel334_ia32_sse.nb334_end: 
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
        popl    %ecx
        popl %ebx
        popl %eax
        leave
        ret




.globl nb_kernel334nf_ia32_sse
.globl _nb_kernel334nf_ia32_sse
nb_kernel334nf_ia32_sse:        
_nb_kernel334nf_ia32_sse:       
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
        ## bottom of stack is cache-aligned for sse use 
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
.set nb334nf_rsqH2H2, 704
.set nb334nf_rsqH2M, 720
.set nb334nf_rsqMH1, 736
.set nb334nf_rsqMH2, 752
.set nb334nf_rsqMM, 768
.set nb334nf_rinvOO, 784
.set nb334nf_rinvH1H1, 800
.set nb334nf_rinvH1H2, 816
.set nb334nf_rinvH1M, 832
.set nb334nf_rinvH2H1, 848
.set nb334nf_rinvH2H2, 864
.set nb334nf_rinvH2M, 880
.set nb334nf_rinvMH1, 896
.set nb334nf_rinvMH2, 912
.set nb334nf_rinvMM, 928
.set nb334nf_is3, 944
.set nb334nf_ii3, 948
.set nb334nf_innerjjnr, 952
.set nb334nf_innerk, 956
.set nb334nf_n, 960
.set nb334nf_nn1, 964
.set nb334nf_nri, 968
.set nb334nf_nouter, 972
.set nb334nf_ninner, 976
.set nb334nf_salign, 980
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $984,%esp          ## local stack space 
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


        movl nb334nf_p_tabscale(%ebp),%eax
        movss (%eax),%xmm5
        shufps $0,%xmm5,%xmm5
        movaps %xmm5,nb334nf_tsc(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb334nf_half(%esp)
        movss nb334nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## 2.0
        addps  %xmm2,%xmm3      ## 3.0
        movaps %xmm1,nb334nf_half(%esp)
        movaps %xmm3,nb334nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb334nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb334nf_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm5
        movss 12(%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movl nb334nf_p_facel(%ebp),%esi
        movss (%esi),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb334nf_qqMM(%esp)
        movaps %xmm4,nb334nf_qqMH(%esp)
        movaps %xmm5,nb334nf_qqHH(%esp)

        xorps %xmm0,%xmm0
        movl  nb334nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb334nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx       ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb334nf_vdwparam(%ebp),%eax
        movlps (%eax,%edx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $0x55,%xmm1,%xmm1
        movaps %xmm0,nb334nf_c6(%esp)
        movaps %xmm1,nb334nf_c12(%esp)

_nb_kernel334nf_ia32_sse.nb334nf_threadloop: 
        movl  nb334nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel334nf_ia32_sse.nb334nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel334nf_ia32_sse.nb334nf_spinlock

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
        jg  _nb_kernel334nf_ia32_sse.nb334nf_outerstart
        jmp _nb_kernel334nf_ia32_sse.nb334nf_end

_nb_kernel334nf_ia32_sse.nb334nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb334nf_nouter(%esp),%ebx
        movl %ebx,nb334nf_nouter(%esp)

_nb_kernel334nf_ia32_sse.nb334nf_outer: 
        movl  nb334nf_shift(%ebp),%eax          ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb334nf_is3(%esp)            ## store is3 

        movl  nb334nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb334nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb334nf_pos(%ebp),%eax    ## eax = base of pos[]  
        movl  %ebx,nb334nf_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        addss (%eax,%ebx,4),%xmm3       ## ox
        addss 4(%eax,%ebx,4),%xmm4     ## oy
        addss 8(%eax,%ebx,4),%xmm5     ## oz
        addss 12(%eax,%ebx,4),%xmm6    ## h1x
        addss 16(%eax,%ebx,4),%xmm7    ## h1y
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm7,%xmm7
        movaps %xmm3,nb334nf_ixO(%esp)
        movaps %xmm4,nb334nf_iyO(%esp)
        movaps %xmm5,nb334nf_izO(%esp)
        movaps %xmm6,nb334nf_ixH1(%esp)
        movaps %xmm7,nb334nf_iyH1(%esp)

        movss %xmm2,%xmm6
        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 20(%eax,%ebx,4),%xmm6    ## h1z
        addss 24(%eax,%ebx,4),%xmm0    ## h2x
        addss 28(%eax,%ebx,4),%xmm1    ## h2y
        addss 32(%eax,%ebx,4),%xmm2    ## h2z
        addss 36(%eax,%ebx,4),%xmm3    ## mx
        addss 40(%eax,%ebx,4),%xmm4    ## my
        addss 44(%eax,%ebx,4),%xmm5    ## mz

        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm6,nb334nf_izH1(%esp)
        movaps %xmm0,nb334nf_ixH2(%esp)
        movaps %xmm1,nb334nf_iyH2(%esp)
        movaps %xmm2,nb334nf_izH2(%esp)
        movaps %xmm3,nb334nf_ixM(%esp)
        movaps %xmm4,nb334nf_iyM(%esp)
        movaps %xmm5,nb334nf_izM(%esp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb334nf_vctot(%esp)
        movaps %xmm4,nb334nf_Vvdwtot(%esp)

        movl  nb334nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb334nf_pos(%ebp),%esi
        movl  nb334nf_faction(%ebp),%edi
        movl  nb334nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb334nf_innerjjnr(%esp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb334nf_ninner(%esp),%ecx
        movl  %ecx,nb334nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb334nf_innerk(%esp)         ## number of innerloop atoms 
        jge   _nb_kernel334nf_ia32_sse.nb334nf_unroll_loop
        jmp   _nb_kernel334nf_ia32_sse.nb334nf_single_check
_nb_kernel334nf_ia32_sse.nb334nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb334nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb334nf_innerjjnr(%esp)             ## advance pointer (unroll 4) 

        movl nb334nf_pos(%ebp),%esi     ## base of pos[] 

        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx        ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables
        ## Load Ox, Oy, Oz, H1x 
        movlps (%esi,%eax,4),%xmm1      ##  Oxa   Oya    -    -
        movlps (%esi,%ecx,4),%xmm4      ##  Oxc   Oyc    -    -
        movhps (%esi,%ebx,4),%xmm1      ##  Oxa   Oya   Oxb   Oyb 
        movhps (%esi,%edx,4),%xmm4      ##  Oxc   Oyc   Oxd   Oyd 
        movaps %xmm1,%xmm0              ##  Oxa   Oya   Oxb   Oyb 
        shufps $0x88,%xmm4,%xmm0       ##  Oxa   Oxb   Oxc   Oxd
        shufps $0xDD,%xmm4,%xmm1       ##  Oya   Oyb   Oyc   Oyd
        movlps 8(%esi,%eax,4),%xmm3     ##  Oza  H1xa    -    -
        movlps 8(%esi,%ecx,4),%xmm5     ##  Ozc  H1xc    -    -
        movhps 8(%esi,%ebx,4),%xmm3     ##  Oza  H1xa   Ozb  H1xb 
        movhps 8(%esi,%edx,4),%xmm5     ##  Ozc  H1xc   Ozd  H1xd 
        movaps %xmm3,%xmm2              ##  Oza  H1xa   Ozb  H1xb 
        shufps $0x88,%xmm5,%xmm2       ##  Oza   Ozb   Ozc   Ozd
        shufps $0xDD,%xmm5,%xmm3       ## H1xa  H1xb  H1xc  H1xd
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb334nf_jxO(%esp)
        movaps %xmm1,nb334nf_jyO(%esp)
        movaps %xmm2,nb334nf_jzO(%esp)
        movaps %xmm3,nb334nf_jxH1(%esp)

        ## Load H1y H1z H2x H2y 
        movlps 16(%esi,%eax,4),%xmm1
        movlps 16(%esi,%ecx,4),%xmm4
        movhps 16(%esi,%ebx,4),%xmm1
        movhps 16(%esi,%edx,4),%xmm4
        movaps %xmm1,%xmm0
        shufps $0x88,%xmm4,%xmm0
        shufps $0xDD,%xmm4,%xmm1
        movlps 24(%esi,%eax,4),%xmm3
        movlps 24(%esi,%ecx,4),%xmm5
        movhps 24(%esi,%ebx,4),%xmm3
        movhps 24(%esi,%edx,4),%xmm5
        movaps %xmm3,%xmm2
        shufps $0x88,%xmm5,%xmm2
        shufps $0xDD,%xmm5,%xmm3
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb334nf_jyH1(%esp)
        movaps %xmm1,nb334nf_jzH1(%esp)
        movaps %xmm2,nb334nf_jxH2(%esp)
        movaps %xmm3,nb334nf_jyH2(%esp)

        ## Load H2z Mx My Mz 
        movlps 32(%esi,%eax,4),%xmm1
        movlps 32(%esi,%ecx,4),%xmm4
        movhps 32(%esi,%ebx,4),%xmm1
        movhps 32(%esi,%edx,4),%xmm4
        movaps %xmm1,%xmm0
        shufps $0x88,%xmm4,%xmm0
        shufps $0xDD,%xmm4,%xmm1
        movlps 40(%esi,%eax,4),%xmm3
        movlps 40(%esi,%ecx,4),%xmm5
        movhps 40(%esi,%ebx,4),%xmm3
        movhps 40(%esi,%edx,4),%xmm5
        movaps %xmm3,%xmm2
        shufps $0x88,%xmm5,%xmm2
        shufps $0xDD,%xmm5,%xmm3
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb334nf_jzH2(%esp)
        movaps %xmm1,nb334nf_jxM(%esp)
        movaps %xmm2,nb334nf_jyM(%esp)
        movaps %xmm3,nb334nf_jzM(%esp)

        ## start calculating pairwise distances
        movaps nb334nf_ixO(%esp),%xmm0
        movaps nb334nf_iyO(%esp),%xmm1
        movaps nb334nf_izO(%esp),%xmm2
        movaps nb334nf_ixH1(%esp),%xmm3
        movaps nb334nf_iyH1(%esp),%xmm4
        movaps nb334nf_izH1(%esp),%xmm5
        subps  nb334nf_jxO(%esp),%xmm0
        subps  nb334nf_jyO(%esp),%xmm1
        subps  nb334nf_jzO(%esp),%xmm2
        subps  nb334nf_jxH1(%esp),%xmm3
        subps  nb334nf_jyH1(%esp),%xmm4
        subps  nb334nf_jzH1(%esp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb334nf_rsqOO(%esp)
        movaps %xmm3,nb334nf_rsqH1H1(%esp)

        movaps nb334nf_ixH1(%esp),%xmm0
        movaps nb334nf_iyH1(%esp),%xmm1
        movaps nb334nf_izH1(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb334nf_jxH2(%esp),%xmm0
        subps  nb334nf_jyH2(%esp),%xmm1
        subps  nb334nf_jzH2(%esp),%xmm2
        subps  nb334nf_jxM(%esp),%xmm3
        subps  nb334nf_jyM(%esp),%xmm4
        subps  nb334nf_jzM(%esp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb334nf_rsqH1H2(%esp)
        movaps %xmm3,nb334nf_rsqH1M(%esp)

        movaps nb334nf_ixH2(%esp),%xmm0
        movaps nb334nf_iyH2(%esp),%xmm1
        movaps nb334nf_izH2(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb334nf_jxH1(%esp),%xmm0
        subps  nb334nf_jyH1(%esp),%xmm1
        subps  nb334nf_jzH1(%esp),%xmm2
        subps  nb334nf_jxH2(%esp),%xmm3
        subps  nb334nf_jyH2(%esp),%xmm4
        subps  nb334nf_jzH2(%esp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb334nf_rsqH2H1(%esp)
        movaps %xmm3,nb334nf_rsqH2H2(%esp)

        movaps nb334nf_ixH2(%esp),%xmm0
        movaps nb334nf_iyH2(%esp),%xmm1
        movaps nb334nf_izH2(%esp),%xmm2
        movaps nb334nf_ixM(%esp),%xmm3
        movaps nb334nf_iyM(%esp),%xmm4
        movaps nb334nf_izM(%esp),%xmm5
        subps  nb334nf_jxM(%esp),%xmm0
        subps  nb334nf_jyM(%esp),%xmm1
        subps  nb334nf_jzM(%esp),%xmm2
        subps  nb334nf_jxH1(%esp),%xmm3
        subps  nb334nf_jyH1(%esp),%xmm4
        subps  nb334nf_jzH1(%esp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb334nf_rsqH2M(%esp)
        movaps %xmm4,nb334nf_rsqMH1(%esp)

        movaps nb334nf_ixM(%esp),%xmm0
        movaps nb334nf_iyM(%esp),%xmm1
        movaps nb334nf_izM(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb334nf_jxH2(%esp),%xmm0
        subps  nb334nf_jyH2(%esp),%xmm1
        subps  nb334nf_jzH2(%esp),%xmm2
        subps  nb334nf_jxM(%esp),%xmm3
        subps  nb334nf_jyM(%esp),%xmm4
        subps  nb334nf_jzM(%esp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb334nf_rsqMH2(%esp)
        movaps %xmm4,nb334nf_rsqMM(%esp)

        ## Invsqrt for O-O
        rsqrtps  nb334nf_rsqOO(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb334nf_three(%esp),%xmm3
        mulps   nb334nf_rsqOO(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb334nf_half(%esp),%xmm3   ## rinvOO
        movaps %xmm3,nb334nf_rinvOO(%esp)

        ## Invsqrt for H1-H1 and H1-H2
        rsqrtps nb334nf_rsqH1H1(%esp),%xmm1
        rsqrtps nb334nf_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb334nf_rsqH1H1(%esp),%xmm1
        mulps   nb334nf_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334nf_half(%esp),%xmm3   ## rinvH1H1 
        mulps   nb334nf_half(%esp),%xmm7   ## rinvH1H2 
        movaps  %xmm3,nb334nf_rinvH1H1(%esp)
        movaps  %xmm7,nb334nf_rinvH1H2(%esp)

        rsqrtps nb334nf_rsqH1M(%esp),%xmm1
        rsqrtps nb334nf_rsqH2H1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb334nf_rsqH1M(%esp),%xmm1
        mulps   nb334nf_rsqH2H1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334nf_half(%esp),%xmm3
        mulps   nb334nf_half(%esp),%xmm7
        movaps  %xmm3,nb334nf_rinvH1M(%esp)
        movaps  %xmm7,nb334nf_rinvH2H1(%esp)

        rsqrtps nb334nf_rsqH2H2(%esp),%xmm1
        rsqrtps nb334nf_rsqH2M(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb334nf_rsqH2H2(%esp),%xmm1
        mulps   nb334nf_rsqH2M(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334nf_half(%esp),%xmm3
        mulps   nb334nf_half(%esp),%xmm7
        movaps  %xmm3,nb334nf_rinvH2H2(%esp)
        movaps  %xmm7,nb334nf_rinvH2M(%esp)

        rsqrtps nb334nf_rsqMH1(%esp),%xmm1
        rsqrtps nb334nf_rsqMH2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb334nf_rsqMH1(%esp),%xmm1
        mulps   nb334nf_rsqMH2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334nf_half(%esp),%xmm3
        mulps   nb334nf_half(%esp),%xmm7
        movaps  %xmm3,nb334nf_rinvMH1(%esp)
        movaps  %xmm7,nb334nf_rinvMH2(%esp)

        rsqrtps nb334nf_rsqMM(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb334nf_three(%esp),%xmm3
        mulps   nb334nf_rsqMM(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb334nf_half(%esp),%xmm3
        movaps  %xmm3,nb334nf_rinvMM(%esp)

        ## start with OO table interaction
        movaps nb334nf_rinvOO(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqOO(%esp),%xmm1   ## xmm1=r
        mulps  nb334nf_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movl nb334nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        ## load dispersion table data into xmm4-xmm7
        movlps 16(%esi,%eax,4),%xmm5
        movlps 16(%esi,%ecx,4),%xmm7
        movhps 16(%esi,%ebx,4),%xmm5
        movhps 16(%esi,%edx,4),%xmm7    ## got half table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 24(%esi,%eax,4),%xmm7
        movlps 24(%esi,%ecx,4),%xmm3
        movhps 24(%esi,%ebx,4),%xmm7
        movhps 24(%esi,%edx,4),%xmm3    ## other half of table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## dispersion table YFGH ready in xmm4-xmm7
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb334nf_c6(%esp),%xmm4
        mulps  %xmm4,%xmm5      ## Vvdw6 

        ## Update Vvdwtot directly 
        addps  nb334nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb334nf_Vvdwtot(%esp)

        ## load repulsion table data into xmm4-xmm7
        movlps 32(%esi,%eax,4),%xmm5
        movlps 32(%esi,%ecx,4),%xmm7
        movhps 32(%esi,%ebx,4),%xmm5
        movhps 32(%esi,%edx,4),%xmm7    ## got half table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 40(%esi,%eax,4),%xmm7
        movlps 40(%esi,%ecx,4),%xmm3
        movhps 40(%esi,%ebx,4),%xmm7
        movhps 40(%esi,%edx,4),%xmm3    ## other half of table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## repulsion table YFGH ready in xmm4-xmm7

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb334nf_c12(%esp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 

        addps  nb334nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb334nf_Vvdwtot(%esp)

        ## Coulomb interactions - first H1H1
        movaps nb334nf_rinvH1H1(%esp),%xmm0

        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulps  nb334nf_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        movl nb334nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb334nf_qqHH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## update vctot 
        addps  nb334nf_vctot(%esp),%xmm5
        movaps %xmm5,nb334nf_vctot(%esp)

        ## H1-H2 interaction 
        movaps nb334nf_rinvH1H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulps  nb334nf_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb334nf_qqHH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb334nf_vctot(%esp),%xmm5
        movaps %xmm5,nb334nf_vctot(%esp)

        ## H1-M interaction  
        movaps nb334nf_rinvH1M(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqH1M(%esp),%xmm1   ## xmm1=r 
        mulps  nb334nf_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx


        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb334nf_qqMH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb334nf_vctot(%esp),%xmm5
        movaps %xmm5,nb334nf_vctot(%esp)

        ## H2-H1 interaction 
        movaps nb334nf_rinvH2H1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulps  nb334nf_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb334nf_qqHH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb334nf_vctot(%esp),%xmm5
        movaps %xmm5,nb334nf_vctot(%esp)

        ## H2-H2 interaction 
        movaps nb334nf_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulps  nb334nf_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb334nf_qqHH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb334nf_vctot(%esp),%xmm5
        movaps %xmm5,nb334nf_vctot(%esp)

        ## H2-M interaction 
        movaps nb334nf_rinvH2M(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqH2M(%esp),%xmm1   ## xmm1=r 
        mulps  nb334nf_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb334nf_qqMH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb334nf_vctot(%esp),%xmm5
        movaps %xmm5,nb334nf_vctot(%esp)

        ## M-H1 interaction 
        movaps nb334nf_rinvMH1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqMH1(%esp),%xmm1   ## xmm1=r 
        mulps  nb334nf_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb334nf_qqMH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb334nf_vctot(%esp),%xmm5
        movaps %xmm5,nb334nf_vctot(%esp)

        ## M-H2 interaction 
        movaps nb334nf_rinvMH2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqMH2(%esp),%xmm1   ## xmm1=r 
        mulps  nb334nf_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb334nf_qqMH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb334nf_vctot(%esp),%xmm5
        movaps %xmm5,nb334nf_vctot(%esp)

        ## M-M interaction 
        movaps nb334nf_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqMM(%esp),%xmm1   ## xmm1=r 
        mulps  nb334nf_tsc(%esp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb334nf_qqMM(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb334nf_vctot(%esp),%xmm5
        movaps %xmm5,nb334nf_vctot(%esp)
        ## should we do one more iteration? 
        subl $4,nb334nf_innerk(%esp)
        jl    _nb_kernel334nf_ia32_sse.nb334nf_single_check
        jmp   _nb_kernel334nf_ia32_sse.nb334nf_unroll_loop
_nb_kernel334nf_ia32_sse.nb334nf_single_check: 
        addl $4,nb334nf_innerk(%esp)
        jnz   _nb_kernel334nf_ia32_sse.nb334nf_single_loop
        jmp   _nb_kernel334nf_ia32_sse.nb334nf_updateouterdata
_nb_kernel334nf_ia32_sse.nb334nf_single_loop: 
        movl  nb334nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb334nf_innerjjnr(%esp)

        movl nb334nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates
        movlps (%esi,%eax,4),%xmm3              ##  Ox  Oy  
        movlps 16(%esi,%eax,4),%xmm4            ## H1y H1z 
        movlps 32(%esi,%eax,4),%xmm5            ## H2z  Mx 
        movhps 8(%esi,%eax,4),%xmm3             ##  Ox  Oy  Oz H1x
        movhps 24(%esi,%eax,4),%xmm4            ## H1y H1z H2x H2y
        movhps 40(%esi,%eax,4),%xmm5            ## H2z  Mx  My  Mz
        ## transpose
        movaps %xmm4,%xmm0
        movaps %xmm3,%xmm1
        movaps %xmm4,%xmm2
        movaps %xmm3,%xmm6
        shufps $18,%xmm5,%xmm4 ## (00010010)  h2x - Mx  - 
        shufps $193,%xmm0,%xmm3 ## (11000001)  Oy  - H1y - 
        shufps $35,%xmm5,%xmm2 ## (00100011) H2y - My  - 
        shufps $18,%xmm0,%xmm1 ## (00010010)  Oz  - H1z - 
        ##  xmm6: Ox - - H1x   xmm5: H2z - - Mz 
        shufps $140,%xmm4,%xmm6 ## (10001100) Ox H1x H2x Mx 
        shufps $136,%xmm2,%xmm3 ## (10001000) Oy H1y H2y My 
        shufps $200,%xmm5,%xmm1 ## (11001000) Oz H1z H2z Mz

        ## store all j coordinates in jO  
        movaps %xmm6,nb334nf_jxO(%esp)
        movaps %xmm3,nb334nf_jyO(%esp)
        movaps %xmm1,nb334nf_jzO(%esp)

        ## do O and H1 in parallel
        movaps nb334nf_ixO(%esp),%xmm0
        movaps nb334nf_iyO(%esp),%xmm1
        movaps nb334nf_izO(%esp),%xmm2
        movaps nb334nf_ixH1(%esp),%xmm3
        movaps nb334nf_iyH1(%esp),%xmm4
        movaps nb334nf_izH1(%esp),%xmm5
        subps  nb334nf_jxO(%esp),%xmm0
        subps  nb334nf_jyO(%esp),%xmm1
        subps  nb334nf_jzO(%esp),%xmm2
        subps  nb334nf_jxO(%esp),%xmm3
        subps  nb334nf_jyO(%esp),%xmm4
        subps  nb334nf_jzO(%esp),%xmm5

        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm3,%xmm4
        addps %xmm5,%xmm4       ## have rsq in xmm4
        movaps %xmm0,nb334nf_rsqOO(%esp)
        movaps %xmm4,nb334nf_rsqH1H1(%esp)

        ## do 1/sqrt(x) for O and  H1
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334nf_half(%esp),%xmm3   ## rinv O - j water 
        mulps   nb334nf_half(%esp),%xmm7   ## rinv H1 - j water  

        movaps %xmm3,nb334nf_rinvOO(%esp)
        movaps %xmm7,nb334nf_rinvH1H1(%esp)

        movl nb334nf_VFtab(%ebp),%esi

        ## do O table LJ interaction
        movaps %xmm3,%xmm0
        movaps %xmm0,%xmm1
        mulss  nb334nf_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulss  nb334nf_tsc(%esp),%xmm1

        cvttps2pi %xmm1,%mm6
        cvtpi2ps %mm6,%xmm3
        subss    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6

        movd %mm6,%ebx
        leal  (%ebx,%ebx,2),%ebx

        ## load dispersion table data into xmm4
        movlps 16(%esi,%ebx,4),%xmm4
        movlps 24(%esi,%ebx,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $0x1,%xmm5,%xmm5
        shufps $0x1,%xmm7,%xmm7
        ## dispersion table YFGH ready in xmm4-xmm7
        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb334nf_c6(%esp),%xmm4
        mulss  %xmm4,%xmm5      ## Vvdw6 

        ## Update Vvdwtot directly 
        addss  nb334nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb334nf_Vvdwtot(%esp)

        ## load repulsion table data into xmm4
        movlps 32(%esi,%ebx,4),%xmm4
        movlps 40(%esi,%ebx,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $0x1,%xmm5,%xmm5
        shufps $0x1,%xmm7,%xmm7
        ## repulsion table YFGH ready in xmm4-xmm7

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb334nf_c12(%esp),%xmm4
        mulss  %xmm4,%xmm5 ## Vvdw12 

        addss  nb334nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb334nf_Vvdwtot(%esp)

        ## do  H1 coulomb interaction
        movaps nb334nf_rinvH1H1(%esp),%xmm0   ## rinv 
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqH1H1(%esp),%xmm1      ## r
        mulps nb334nf_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movl nb334_VFtab(%ebp),%esi

        movlps (%esi,%ebx,4),%xmm4
        movlps (%esi,%ecx,4),%xmm3
        movlps (%esi,%edx,4),%xmm7
        movhps 8(%esi,%ebx,4),%xmm4
        movhps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm7
        movaps %xmm3,%xmm6
        unpcklps %xmm7,%xmm6
        unpckhps %xmm7,%xmm3
        movaps %xmm4,%xmm5
        movaps %xmm4,%xmm7
        shufps $0x40,%xmm6,%xmm4
        shufps $0xE4,%xmm6,%xmm5
        movaps %xmm7,%xmm6
        shufps $0x48,%xmm3,%xmm6
        shufps $0xEC,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb334nf_qqHH(%esp),%xmm3
        movhps  nb334nf_qqMH(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001 

        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 

        addps  nb334nf_vctot(%esp),%xmm5
        movaps %xmm5,nb334nf_vctot(%esp)

        ## i H2 & M simultaneously first get i particle coords: 
        movaps  nb334nf_ixH2(%esp),%xmm0
        movaps  nb334nf_iyH2(%esp),%xmm1
        movaps  nb334nf_izH2(%esp),%xmm2
        movaps  nb334nf_ixM(%esp),%xmm3
        movaps  nb334nf_iyM(%esp),%xmm4
        movaps  nb334nf_izM(%esp),%xmm5
        subps   nb334nf_jxO(%esp),%xmm0
        subps   nb334nf_jyO(%esp),%xmm1
        subps   nb334nf_jzO(%esp),%xmm2
        subps   nb334nf_jxO(%esp),%xmm3
        subps   nb334nf_jyO(%esp),%xmm4
        subps   nb334nf_jzO(%esp),%xmm5
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm1,%xmm0
        addps %xmm3,%xmm4
        addps %xmm2,%xmm0       ## have rsqH2 in xmm0 
        addps %xmm5,%xmm4       ## have rsqM in xmm4 

        ## start with H2, save data 
        movaps %xmm0,nb334nf_rsqH2H2(%esp)
        movaps %xmm4,nb334nf_rsqMM(%esp)
        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334nf_half(%esp),%xmm3   ## rinv H2 - j water 
        mulps   nb334nf_half(%esp),%xmm7   ## rinv M - j water  

        movaps %xmm3,nb334nf_rinvH2H2(%esp)
        movaps %xmm7,nb334nf_rinvMM(%esp)

        movaps %xmm3,%xmm1
        mulps  nb334nf_rsqH2H2(%esp),%xmm1      ## xmm1=r 
        movaps %xmm3,%xmm0      ## xmm0=rinv 
        mulps  nb334nf_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%ebx,4),%xmm4
        movlps (%esi,%ecx,4),%xmm3
        movlps (%esi,%edx,4),%xmm7
        movhps 8(%esi,%ebx,4),%xmm4
        movhps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm7
        movaps %xmm3,%xmm6
        unpcklps %xmm7,%xmm6
        unpckhps %xmm7,%xmm3
        movaps %xmm4,%xmm5
        movaps %xmm4,%xmm7
        shufps $0x40,%xmm6,%xmm4
        shufps $0xE4,%xmm6,%xmm5
        movaps %xmm7,%xmm6
        shufps $0x48,%xmm3,%xmm6
        shufps $0xEC,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 

        xorps  %xmm3,%xmm3

        ## fetch charges to xmm3 (temporary) 
        movss   nb334nf_qqHH(%esp),%xmm3
        movhps  nb334nf_qqMH(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001

        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 
        addps  nb334nf_vctot(%esp),%xmm5
        movaps %xmm5,nb334nf_vctot(%esp)

        ## do table for i M - j water interaction 
        movaps nb334nf_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqMM(%esp),%xmm1        ## xmm0=rinv, xmm1=r 
        mulps  nb334nf_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movl nb334_VFtab(%ebp),%esi

        movlps (%esi,%ebx,4),%xmm4
        movlps (%esi,%ecx,4),%xmm3
        movlps (%esi,%edx,4),%xmm7
        movhps 8(%esi,%ebx,4),%xmm4
        movhps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm7
        movaps %xmm3,%xmm6
        unpcklps %xmm7,%xmm6
        unpckhps %xmm7,%xmm3
        movaps %xmm4,%xmm5
        movaps %xmm4,%xmm7
        shufps $0x40,%xmm6,%xmm4
        shufps $0xE4,%xmm6,%xmm5
        movaps %xmm7,%xmm6
        shufps $0x48,%xmm3,%xmm6
        shufps $0xEC,%xmm3,%xmm7
        ## # coulomb table ready, in xmm4-xmm7

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb334nf_qqMH(%esp),%xmm3
        movhps  nb334nf_qqMM(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001

        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul
        addps  nb334nf_vctot(%esp),%xmm5
        movaps %xmm5,nb334nf_vctot(%esp)

        decl nb334nf_innerk(%esp)
        jz    _nb_kernel334nf_ia32_sse.nb334nf_updateouterdata
        jmp   _nb_kernel334nf_ia32_sse.nb334nf_single_loop
_nb_kernel334nf_ia32_sse.nb334nf_updateouterdata: 
        ## get n from stack
        movl nb334nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb334nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb334nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb334nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb334nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb334nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb334nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel334nf_ia32_sse.nb334nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb334nf_n(%esp)
        jmp _nb_kernel334nf_ia32_sse.nb334nf_outer
_nb_kernel334nf_ia32_sse.nb334nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb334nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel334nf_ia32_sse.nb334nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel334nf_ia32_sse.nb334nf_threadloop
_nb_kernel334nf_ia32_sse.nb334nf_end: 
        emms

        movl nb334nf_nouter(%esp),%eax
        movl nb334nf_ninner(%esp),%ebx
        movl nb334nf_outeriter(%ebp),%ecx
        movl nb334nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb334nf_salign(%esp),%eax
        addl %eax,%esp
        addl $984,%esp
        popl %edi
        popl %esi
        popl %edx
        popl    %ecx
        popl %ebx
        popl %eax
        leave
        ret

