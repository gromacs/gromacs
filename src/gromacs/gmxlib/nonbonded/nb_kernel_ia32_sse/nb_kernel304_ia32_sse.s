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




.globl nb_kernel304_ia32_sse
.globl _nb_kernel304_ia32_sse
nb_kernel304_ia32_sse:  
_nb_kernel304_ia32_sse: 
.set nb304_p_nri, 8
.set nb304_iinr, 12
.set nb304_jindex, 16
.set nb304_jjnr, 20
.set nb304_shift, 24
.set nb304_shiftvec, 28
.set nb304_fshift, 32
.set nb304_gid, 36
.set nb304_pos, 40
.set nb304_faction, 44
.set nb304_charge, 48
.set nb304_p_facel, 52
.set nb304_argkrf, 56
.set nb304_argcrf, 60
.set nb304_Vc, 64
.set nb304_type, 68
.set nb304_p_ntype, 72
.set nb304_vdwparam, 76
.set nb304_Vvdw, 80
.set nb304_p_tabscale, 84
.set nb304_VFtab, 88
.set nb304_invsqrta, 92
.set nb304_dvda, 96
.set nb304_p_gbtabscale, 100
.set nb304_GBtab, 104
.set nb304_p_nthreads, 108
.set nb304_count, 112
.set nb304_mtx, 116
.set nb304_outeriter, 120
.set nb304_inneriter, 124
.set nb304_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb304_ixH1, 0
.set nb304_iyH1, 16
.set nb304_izH1, 32
.set nb304_ixH2, 48
.set nb304_iyH2, 64
.set nb304_izH2, 80
.set nb304_ixM, 96
.set nb304_iyM, 112
.set nb304_izM, 128
.set nb304_jxH1, 144
.set nb304_jyH1, 160
.set nb304_jzH1, 176
.set nb304_jxH2, 192
.set nb304_jyH2, 208
.set nb304_jzH2, 224
.set nb304_jxM, 240
.set nb304_jyM, 256
.set nb304_jzM, 272
.set nb304_dxH1H1, 288
.set nb304_dyH1H1, 304
.set nb304_dzH1H1, 320
.set nb304_dxH1H2, 336
.set nb304_dyH1H2, 352
.set nb304_dzH1H2, 368
.set nb304_dxH1M, 384
.set nb304_dyH1M, 400
.set nb304_dzH1M, 416
.set nb304_dxH2H1, 432
.set nb304_dyH2H1, 448
.set nb304_dzH2H1, 464
.set nb304_dxH2H2, 480
.set nb304_dyH2H2, 496
.set nb304_dzH2H2, 512
.set nb304_dxH2M, 528
.set nb304_dyH2M, 544
.set nb304_dzH2M, 560
.set nb304_dxMH1, 576
.set nb304_dyMH1, 592
.set nb304_dzMH1, 608
.set nb304_dxMH2, 624
.set nb304_dyMH2, 640
.set nb304_dzMH2, 656
.set nb304_dxMM, 672
.set nb304_dyMM, 688
.set nb304_dzMM, 704
.set nb304_qqHH, 720
.set nb304_qqMH, 736
.set nb304_qqMM, 752
.set nb304_two, 768
.set nb304_tsc, 784
.set nb304_vctot, 800
.set nb304_fixH1, 816
.set nb304_fiyH1, 832
.set nb304_fizH1, 848
.set nb304_fixH2, 864
.set nb304_fiyH2, 880
.set nb304_fizH2, 896
.set nb304_fixM, 912
.set nb304_fiyM, 928
.set nb304_fizM, 944
.set nb304_fjxH1, 960
.set nb304_fjyH1, 976
.set nb304_fjzH1, 992
.set nb304_fjxH2, 1008
.set nb304_fjyH2, 1024
.set nb304_fjzH2, 1040
.set nb304_fjxM, 1056
.set nb304_fjyM, 1072
.set nb304_fjzM, 1088
.set nb304_fjzMb, 1092
.set nb304_fjzMc, 1096
.set nb304_fjzMd, 1100
.set nb304_half, 1104
.set nb304_three, 1120
.set nb304_rsqH1H1, 1136
.set nb304_rsqH1H2, 1152
.set nb304_rsqH1M, 1168
.set nb304_rsqH2H1, 1184
.set nb304_rsqH2H2, 1200
.set nb304_rsqH2M, 1216
.set nb304_rsqMH1, 1232
.set nb304_rsqMH2, 1248
.set nb304_rsqMM, 1264
.set nb304_rinvH1H1, 1280
.set nb304_rinvH1H2, 1296
.set nb304_rinvH1M, 1312
.set nb304_rinvH2H1, 1328
.set nb304_rinvH2H2, 1344
.set nb304_rinvH2M, 1360
.set nb304_rinvMH1, 1376
.set nb304_rinvMH2, 1392
.set nb304_rinvMM, 1408
.set nb304_is3, 1424
.set nb304_ii3, 1428
.set nb304_innerjjnr, 1432
.set nb304_innerk, 1436
.set nb304_n, 1440
.set nb304_nn1, 1444
.set nb304_nri, 1448
.set nb304_nouter, 1452
.set nb304_ninner, 1456
.set nb304_salign, 1460
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $1464,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb304_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb304_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb304_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb304_nouter(%esp)
        movl %eax,nb304_ninner(%esp)


        movl nb304_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb304_tsc(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb304_half(%esp)
        movss nb304_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb304_half(%esp)
        movaps %xmm2,nb304_two(%esp)
        movaps %xmm3,nb304_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb304_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb304_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movss 12(%edx,%ebx,4),%xmm5
        movl nb304_p_facel(%ebp),%esi
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
        movaps %xmm3,nb304_qqHH(%esp)
        movaps %xmm4,nb304_qqMH(%esp)
        movaps %xmm5,nb304_qqMM(%esp)

_nb_kernel304_ia32_sse.nb304_threadloop: 
        movl  nb304_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel304_ia32_sse.nb304_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel304_ia32_sse.nb304_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb304_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb304_n(%esp)
        movl %ebx,nb304_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel304_ia32_sse.nb304_outerstart
        jmp _nb_kernel304_ia32_sse.nb304_end

_nb_kernel304_ia32_sse.nb304_outerstart: 
        ## ebx contains number of outer iterations
        addl nb304_nouter(%esp),%ebx
        movl %ebx,nb304_nouter(%esp)

_nb_kernel304_ia32_sse.nb304_outer: 
        movl  nb304_shift(%ebp),%eax            ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb304_is3(%esp)      ## store is3 

        movl  nb304_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb304_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb304_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb304_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss 12(%eax,%ebx,4),%xmm3
        addss 16(%eax,%ebx,4),%xmm4
        addss 20(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb304_ixH1(%esp)
        movaps %xmm4,nb304_iyH1(%esp)
        movaps %xmm5,nb304_izH1(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 24(%eax,%ebx,4),%xmm0
        addss 28(%eax,%ebx,4),%xmm1
        addss 32(%eax,%ebx,4),%xmm2
        addss 36(%eax,%ebx,4),%xmm3
        addss 40(%eax,%ebx,4),%xmm4
        addss 44(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb304_ixH2(%esp)
        movaps %xmm1,nb304_iyH2(%esp)
        movaps %xmm2,nb304_izH2(%esp)
        movaps %xmm3,nb304_ixM(%esp)
        movaps %xmm4,nb304_iyM(%esp)
        movaps %xmm5,nb304_izM(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb304_vctot(%esp)
        movaps %xmm4,nb304_fixH1(%esp)
        movaps %xmm4,nb304_fiyH1(%esp)
        movaps %xmm4,nb304_fizH1(%esp)
        movaps %xmm4,nb304_fixH2(%esp)
        movaps %xmm4,nb304_fiyH2(%esp)
        movaps %xmm4,nb304_fizH2(%esp)
        movaps %xmm4,nb304_fixM(%esp)
        movaps %xmm4,nb304_fiyM(%esp)
        movaps %xmm4,nb304_fizM(%esp)

        movl  nb304_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb304_pos(%ebp),%esi
        movl  nb304_faction(%ebp),%edi
        movl  nb304_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb304_innerjjnr(%esp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb304_ninner(%esp),%ecx
        movl  %ecx,nb304_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb304_innerk(%esp)   ## number of innerloop atoms 
        jge   _nb_kernel304_ia32_sse.nb304_unroll_loop
        jmp   _nb_kernel304_ia32_sse.nb304_single_check
_nb_kernel304_ia32_sse.nb304_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb304_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb304_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb304_pos(%ebp),%esi       ## base of pos[] 

        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx        ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables 
        movlps 12(%esi,%eax,4),%xmm2
        movlps 24(%esi,%eax,4),%xmm3
        movlps 36(%esi,%eax,4),%xmm4

        movlps 12(%esi,%ebx,4),%xmm5
        movlps 24(%esi,%ebx,4),%xmm6
        movlps 36(%esi,%ebx,4),%xmm7

        movhps 12(%esi,%ecx,4),%xmm2
        movhps 24(%esi,%ecx,4),%xmm3
        movhps 36(%esi,%ecx,4),%xmm4

        movhps 12(%esi,%edx,4),%xmm5
        movhps 24(%esi,%edx,4),%xmm6
        movhps 36(%esi,%edx,4),%xmm7

        movaps %xmm2,%xmm0
        movaps %xmm3,%xmm1
        unpcklps %xmm5,%xmm0
        unpcklps %xmm6,%xmm1
        unpckhps %xmm5,%xmm2
        unpckhps %xmm6,%xmm3
        movaps %xmm4,%xmm5
        movaps   %xmm0,%xmm6
        unpcklps %xmm7,%xmm4
        unpckhps %xmm7,%xmm5
        movaps   %xmm1,%xmm7
        movlhps  %xmm2,%xmm0
        movaps %xmm0,nb304_jxH1(%esp)
        movhlps  %xmm6,%xmm2
        movaps %xmm2,nb304_jyH1(%esp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb304_jxH2(%esp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb304_jyH2(%esp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb304_jxM(%esp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb304_jyM(%esp)

        movss  20(%esi,%eax,4),%xmm0
        movss  32(%esi,%eax,4),%xmm1
        movss  44(%esi,%eax,4),%xmm2

        movss  20(%esi,%ecx,4),%xmm3
        movss  32(%esi,%ecx,4),%xmm4
        movss  44(%esi,%ecx,4),%xmm5

        movhps 16(%esi,%ebx,4),%xmm0
        movhps 28(%esi,%ebx,4),%xmm1
        movhps 40(%esi,%ebx,4),%xmm2

        movhps 16(%esi,%edx,4),%xmm3
        movhps 28(%esi,%edx,4),%xmm4
        movhps 40(%esi,%edx,4),%xmm5

        shufps $204,%xmm3,%xmm0 ## constant 11001100
        shufps $204,%xmm4,%xmm1 ## constant 11001100
        shufps $204,%xmm5,%xmm2 ## constant 11001100
        movaps %xmm0,nb304_jzH1(%esp)
        movaps %xmm1,nb304_jzH2(%esp)
        movaps %xmm2,nb304_jzM(%esp)

        movaps nb304_ixH1(%esp),%xmm0
        movaps nb304_iyH1(%esp),%xmm1
        movaps nb304_izH1(%esp),%xmm2
        movaps nb304_ixH1(%esp),%xmm3
        movaps nb304_iyH1(%esp),%xmm4
        movaps nb304_izH1(%esp),%xmm5
        subps  nb304_jxH1(%esp),%xmm0
        subps  nb304_jyH1(%esp),%xmm1
        subps  nb304_jzH1(%esp),%xmm2
        subps  nb304_jxH2(%esp),%xmm3
        subps  nb304_jyH2(%esp),%xmm4
        subps  nb304_jzH2(%esp),%xmm5
        movaps %xmm0,nb304_dxH1H1(%esp)
        movaps %xmm1,nb304_dyH1H1(%esp)
        movaps %xmm2,nb304_dzH1H1(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb304_dxH1H2(%esp)
        movaps %xmm4,nb304_dyH1H2(%esp)
        movaps %xmm5,nb304_dzH1H2(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb304_rsqH1H1(%esp)
        movaps %xmm3,nb304_rsqH1H2(%esp)

        movaps nb304_ixH1(%esp),%xmm0
        movaps nb304_iyH1(%esp),%xmm1
        movaps nb304_izH1(%esp),%xmm2
        movaps nb304_ixH2(%esp),%xmm3
        movaps nb304_iyH2(%esp),%xmm4
        movaps nb304_izH2(%esp),%xmm5
        subps  nb304_jxM(%esp),%xmm0
        subps  nb304_jyM(%esp),%xmm1
        subps  nb304_jzM(%esp),%xmm2
        subps  nb304_jxH1(%esp),%xmm3
        subps  nb304_jyH1(%esp),%xmm4
        subps  nb304_jzH1(%esp),%xmm5
        movaps %xmm0,nb304_dxH1M(%esp)
        movaps %xmm1,nb304_dyH1M(%esp)
        movaps %xmm2,nb304_dzH1M(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb304_dxH2H1(%esp)
        movaps %xmm4,nb304_dyH2H1(%esp)
        movaps %xmm5,nb304_dzH2H1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb304_rsqH1M(%esp)
        movaps %xmm3,nb304_rsqH2H1(%esp)

        movaps nb304_ixH2(%esp),%xmm0
        movaps nb304_iyH2(%esp),%xmm1
        movaps nb304_izH2(%esp),%xmm2
        movaps nb304_ixH2(%esp),%xmm3
        movaps nb304_iyH2(%esp),%xmm4
        movaps nb304_izH2(%esp),%xmm5
        subps  nb304_jxH2(%esp),%xmm0
        subps  nb304_jyH2(%esp),%xmm1
        subps  nb304_jzH2(%esp),%xmm2
        subps  nb304_jxM(%esp),%xmm3
        subps  nb304_jyM(%esp),%xmm4
        subps  nb304_jzM(%esp),%xmm5
        movaps %xmm0,nb304_dxH2H2(%esp)
        movaps %xmm1,nb304_dyH2H2(%esp)
        movaps %xmm2,nb304_dzH2H2(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb304_dxH2M(%esp)
        movaps %xmm4,nb304_dyH2M(%esp)
        movaps %xmm5,nb304_dzH2M(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb304_rsqH2H2(%esp)
        movaps %xmm3,nb304_rsqH2M(%esp)

        movaps nb304_ixM(%esp),%xmm0
        movaps nb304_iyM(%esp),%xmm1
        movaps nb304_izM(%esp),%xmm2
        movaps nb304_ixM(%esp),%xmm3
        movaps nb304_iyM(%esp),%xmm4
        movaps nb304_izM(%esp),%xmm5
        subps  nb304_jxH1(%esp),%xmm0
        subps  nb304_jyH1(%esp),%xmm1
        subps  nb304_jzH1(%esp),%xmm2
        subps  nb304_jxH2(%esp),%xmm3
        subps  nb304_jyH2(%esp),%xmm4
        subps  nb304_jzH2(%esp),%xmm5
        movaps %xmm0,nb304_dxMH1(%esp)
        movaps %xmm1,nb304_dyMH1(%esp)
        movaps %xmm2,nb304_dzMH1(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb304_dxMH2(%esp)
        movaps %xmm4,nb304_dyMH2(%esp)
        movaps %xmm5,nb304_dzMH2(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb304_rsqMH1(%esp)
        movaps %xmm4,nb304_rsqMH2(%esp)

        movaps nb304_ixM(%esp),%xmm0
        movaps nb304_iyM(%esp),%xmm1
        movaps nb304_izM(%esp),%xmm2
        subps  nb304_jxM(%esp),%xmm0
        subps  nb304_jyM(%esp),%xmm1
        subps  nb304_jzM(%esp),%xmm2
        movaps %xmm0,nb304_dxMM(%esp)
        movaps %xmm1,nb304_dyMM(%esp)
        movaps %xmm2,nb304_dzMM(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb304_rsqMM(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb304_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb304_half(%esp),%xmm3   ## rinvMM
        mulps   nb304_half(%esp),%xmm7   ## rinvMH2 
        movaps  %xmm3,nb304_rinvMM(%esp)
        movaps  %xmm7,nb304_rinvMH2(%esp)

        rsqrtps nb304_rsqH1H1(%esp),%xmm1
        rsqrtps nb304_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb304_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb304_rsqH1H1(%esp),%xmm1
        mulps   nb304_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb304_half(%esp),%xmm3
        mulps   nb304_half(%esp),%xmm7
        movaps  %xmm3,nb304_rinvH1H1(%esp)
        movaps  %xmm7,nb304_rinvH1H2(%esp)

        rsqrtps nb304_rsqH1M(%esp),%xmm1
        rsqrtps nb304_rsqH2H1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb304_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb304_rsqH1M(%esp),%xmm1
        mulps   nb304_rsqH2H1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb304_half(%esp),%xmm3
        mulps   nb304_half(%esp),%xmm7
        movaps  %xmm3,nb304_rinvH1M(%esp)
        movaps  %xmm7,nb304_rinvH2H1(%esp)

        rsqrtps nb304_rsqH2H2(%esp),%xmm1
        rsqrtps nb304_rsqH2M(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb304_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb304_rsqH2H2(%esp),%xmm1
        mulps   nb304_rsqH2M(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb304_half(%esp),%xmm3
        mulps   nb304_half(%esp),%xmm7
        movaps  %xmm3,nb304_rinvH2H2(%esp)
        movaps  %xmm7,nb304_rinvH2M(%esp)

        rsqrtps nb304_rsqMH1(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb304_three(%esp),%xmm3
        mulps   nb304_rsqMH1(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb304_half(%esp),%xmm3
        movaps  %xmm3,nb304_rinvMH1(%esp)

        ## start with H1-H1 interaction 
        movaps nb304_rinvH1H1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulps  nb304_tsc(%esp),%xmm1

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

        movl nb304_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb304_two(%esp),%xmm7            ## two*Heps2 
        movaps nb304_qqHH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb304_vctot(%esp),%xmm5
        xorps  %xmm2,%xmm2
        movaps %xmm5,nb304_vctot(%esp)
        mulps  nb304_tsc(%esp),%xmm3

        subps  %xmm3,%xmm2
        mulps  %xmm2,%xmm0

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb304_dxH1H1(%esp),%xmm0
        mulps nb304_dyH1H1(%esp),%xmm1
        mulps nb304_dzH1H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb304_fixH1(%esp),%xmm0
        addps nb304_fiyH1(%esp),%xmm1
        addps nb304_fizH1(%esp),%xmm2
        movaps %xmm3,nb304_fjxH1(%esp)
        movaps %xmm4,nb304_fjyH1(%esp)
        movaps %xmm5,nb304_fjzH1(%esp)
        movaps %xmm0,nb304_fixH1(%esp)
        movaps %xmm1,nb304_fiyH1(%esp)
        movaps %xmm2,nb304_fizH1(%esp)

        ## H1-H2 interaction 
        movaps nb304_rinvH1H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulps  nb304_tsc(%esp),%xmm1
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

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb304_two(%esp),%xmm7            ## two*Heps2 
        movaps nb304_qqHH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb304_vctot(%esp),%xmm5
        movaps %xmm5,nb304_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb304_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb304_dxH1H2(%esp),%xmm0
        mulps nb304_dyH1H2(%esp),%xmm1
        mulps nb304_dzH1H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb304_fixH1(%esp),%xmm0
        addps nb304_fiyH1(%esp),%xmm1
        addps nb304_fizH1(%esp),%xmm2
        movaps %xmm3,nb304_fjxH2(%esp)
        movaps %xmm4,nb304_fjyH2(%esp)
        movaps %xmm5,nb304_fjzH2(%esp)
        movaps %xmm0,nb304_fixH1(%esp)
        movaps %xmm1,nb304_fiyH1(%esp)
        movaps %xmm2,nb304_fizH1(%esp)

        ## H1-M interaction  
        movaps nb304_rinvH1M(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304_rsqH1M(%esp),%xmm1   ## xmm1=r 
        mulps  nb304_tsc(%esp),%xmm1
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

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb304_two(%esp),%xmm7            ## two*Heps2 
        movaps nb304_qqMH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb304_vctot(%esp),%xmm5
        movaps %xmm5,nb304_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb304_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb304_dxH1M(%esp),%xmm0
        mulps nb304_dyH1M(%esp),%xmm1
        mulps nb304_dzH1M(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb304_fixH1(%esp),%xmm0
        addps nb304_fiyH1(%esp),%xmm1
        addps nb304_fizH1(%esp),%xmm2
        movaps %xmm3,nb304_fjxM(%esp)
        movaps %xmm4,nb304_fjyM(%esp)
        movaps %xmm5,nb304_fjzM(%esp)
        movaps %xmm0,nb304_fixH1(%esp)
        movaps %xmm1,nb304_fiyH1(%esp)
        movaps %xmm2,nb304_fizH1(%esp)

        ## H2-H1 interaction 
        movaps nb304_rinvH2H1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulps  nb304_tsc(%esp),%xmm1
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

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb304_two(%esp),%xmm7            ## two*Heps2 
        movaps nb304_qqHH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb304_vctot(%esp),%xmm5
        movaps %xmm5,nb304_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb304_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        movaps nb304_fjxH1(%esp),%xmm3
        movaps nb304_fjyH1(%esp),%xmm4
        movaps nb304_fjzH1(%esp),%xmm5
        mulps nb304_dxH2H1(%esp),%xmm0
        mulps nb304_dyH2H1(%esp),%xmm1
        mulps nb304_dzH2H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb304_fixH2(%esp),%xmm0
        addps nb304_fiyH2(%esp),%xmm1
        addps nb304_fizH2(%esp),%xmm2
        movaps %xmm3,nb304_fjxH1(%esp)
        movaps %xmm4,nb304_fjyH1(%esp)
        movaps %xmm5,nb304_fjzH1(%esp)
        movaps %xmm0,nb304_fixH2(%esp)
        movaps %xmm1,nb304_fiyH2(%esp)
        movaps %xmm2,nb304_fizH2(%esp)

        ## H2-H2 interaction 
        movaps nb304_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulps  nb304_tsc(%esp),%xmm1
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

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb304_two(%esp),%xmm7            ## two*Heps2 
        movaps nb304_qqHH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb304_vctot(%esp),%xmm5
        movaps %xmm5,nb304_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb304_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        movaps nb304_fjxH2(%esp),%xmm3
        movaps nb304_fjyH2(%esp),%xmm4
        movaps nb304_fjzH2(%esp),%xmm5
        mulps nb304_dxH2H2(%esp),%xmm0
        mulps nb304_dyH2H2(%esp),%xmm1
        mulps nb304_dzH2H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb304_fixH2(%esp),%xmm0
        addps nb304_fiyH2(%esp),%xmm1
        addps nb304_fizH2(%esp),%xmm2
        movaps %xmm3,nb304_fjxH2(%esp)
        movaps %xmm4,nb304_fjyH2(%esp)
        movaps %xmm5,nb304_fjzH2(%esp)
        movaps %xmm0,nb304_fixH2(%esp)
        movaps %xmm1,nb304_fiyH2(%esp)
        movaps %xmm2,nb304_fizH2(%esp)

        ## H2-M interaction 
        movaps nb304_rinvH2M(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304_rsqH2M(%esp),%xmm1   ## xmm1=r 
        mulps  nb304_tsc(%esp),%xmm1
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

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb304_two(%esp),%xmm7            ## two*Heps2 
        movaps nb304_qqMH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb304_vctot(%esp),%xmm5
        movaps %xmm5,nb304_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb304_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        movaps nb304_fjxM(%esp),%xmm3
        movaps nb304_fjyM(%esp),%xmm4
        movaps nb304_fjzM(%esp),%xmm5
        mulps nb304_dxH2M(%esp),%xmm0
        mulps nb304_dyH2M(%esp),%xmm1
        mulps nb304_dzH2M(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb304_fixH2(%esp),%xmm0
        addps nb304_fiyH2(%esp),%xmm1
        addps nb304_fizH2(%esp),%xmm2
        movaps %xmm3,nb304_fjxM(%esp)
        movaps %xmm4,nb304_fjyM(%esp)
        movaps %xmm5,nb304_fjzM(%esp)
        movaps %xmm0,nb304_fixH2(%esp)
        movaps %xmm1,nb304_fiyH2(%esp)
        movaps %xmm2,nb304_fizH2(%esp)

        ## M-H1 interaction 
        movaps nb304_rinvMH1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304_rsqMH1(%esp),%xmm1   ## xmm1=r 
        mulps  nb304_tsc(%esp),%xmm1
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

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb304_two(%esp),%xmm7            ## two*Heps2 
        movaps nb304_qqMH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb304_vctot(%esp),%xmm5
        movaps %xmm5,nb304_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb304_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        movaps nb304_fjxH1(%esp),%xmm3
        movaps nb304_fjyH1(%esp),%xmm4
        movaps nb304_fjzH1(%esp),%xmm5
        mulps nb304_dxMH1(%esp),%xmm0
        mulps nb304_dyMH1(%esp),%xmm1
        mulps nb304_dzMH1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb304_fixM(%esp),%xmm0
        addps nb304_fiyM(%esp),%xmm1
        addps nb304_fizM(%esp),%xmm2
        movaps %xmm3,nb304_fjxH1(%esp)
        movaps %xmm4,nb304_fjyH1(%esp)
        movaps %xmm5,nb304_fjzH1(%esp)
        movaps %xmm0,nb304_fixM(%esp)
        movaps %xmm1,nb304_fiyM(%esp)
        movaps %xmm2,nb304_fizM(%esp)

        ## M-H2 interaction 
        movaps nb304_rinvMH2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304_rsqMH2(%esp),%xmm1   ## xmm1=r 
        mulps  nb304_tsc(%esp),%xmm1
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

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb304_two(%esp),%xmm7            ## two*Heps2 
        movaps nb304_qqMH(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb304_vctot(%esp),%xmm5
        movaps %xmm5,nb304_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb304_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        movaps nb304_fjxH2(%esp),%xmm3
        movaps nb304_fjyH2(%esp),%xmm4
        movaps nb304_fjzH2(%esp),%xmm5
        mulps nb304_dxMH2(%esp),%xmm0
        mulps nb304_dyMH2(%esp),%xmm1
        mulps nb304_dzMH2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb304_fixM(%esp),%xmm0
        addps nb304_fiyM(%esp),%xmm1
        addps nb304_fizM(%esp),%xmm2
        movaps %xmm3,nb304_fjxH2(%esp)
        movaps %xmm4,nb304_fjyH2(%esp)
        movaps %xmm5,nb304_fjzH2(%esp)
        movaps %xmm0,nb304_fixM(%esp)
        movaps %xmm1,nb304_fiyM(%esp)
        movaps %xmm2,nb304_fizM(%esp)

        ## M-M interaction 
        movaps nb304_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304_rsqMM(%esp),%xmm1   ## xmm1=r 
        mulps  nb304_tsc(%esp),%xmm1
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

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb304_two(%esp),%xmm7            ## two*Heps2 
        movaps nb304_qqMM(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 

        addps  nb304_vctot(%esp),%xmm5
        movaps %xmm5,nb304_vctot(%esp)
        xorps  %xmm1,%xmm1
        mulps  nb304_tsc(%esp),%xmm3
        mulps  %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        movaps nb304_fjxM(%esp),%xmm3
        movaps nb304_fjyM(%esp),%xmm4
        movaps nb304_fjzM(%esp),%xmm5
        mulps nb304_dxMM(%esp),%xmm0
        mulps nb304_dyMM(%esp),%xmm1
        mulps nb304_dzMM(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb304_fixM(%esp),%xmm0
        addps nb304_fiyM(%esp),%xmm1
        addps nb304_fizM(%esp),%xmm2
        movaps %xmm3,nb304_fjxM(%esp)
        movaps %xmm4,nb304_fjyM(%esp)
        movaps %xmm5,nb304_fjzM(%esp)
        movaps %xmm0,nb304_fixM(%esp)
        movaps %xmm1,nb304_fiyM(%esp)
        movaps %xmm2,nb304_fizM(%esp)

        movl nb304_faction(%ebp),%edi

        movd %mm0,%eax
        movd %mm1,%ebx
        movd %mm2,%ecx
        movd %mm3,%edx

        ## Did all interactions - now update j forces 
        ## At this stage forces are still on the stack, in positions:
        ## fjxH1, fjyH1, fjzH1, ... , fjzM.
        ## Each position is a quadruplet of forces for the four 
        ## corresponding j waters, so we need to transpose them before
        ## adding to the memory positions.
        ## 
        ## This _used_ to be a simple transpose, but the resulting high number
        ## of unaligned 128-bit load/stores might trigger a possible hardware 
        ## bug on Athlon and Opteron chips, so I have worked around it
        ## to use 64-bit load/stores instead. The performance hit should be
        ## very modest, since the 128-bit unaligned memory instructions were
        ## slow anyway. 

        ## 4 j waters with three atoms each - first do 1st Hydrogen X & Y forces for 4 j particles 
        movaps nb304_fjxH1(%esp),%xmm0   ## xmm0= fjxH1a  fjxH1b  fjxH1c  fjxH1d 
        movaps nb304_fjyH1(%esp),%xmm2   ## xmm1= fjyH1a  fjyH1b  fjyH1c  fjyH1d
        movlps 12(%edi,%eax,4),%xmm3
        movlps 12(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0       ## xmm0= fjxH1a  fjyH1a  fjxH1b  fjyH1b
        unpckhps %xmm2,%xmm1       ## xmm1= fjxH1c  fjyH1c  fjxH1d  fjyH1d
        movhps 12(%edi,%ebx,4),%xmm3
        movhps 12(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,12(%edi,%eax,4)
        movlps %xmm4,12(%edi,%ecx,4)
        movhps %xmm3,12(%edi,%ebx,4)
        movhps %xmm4,12(%edi,%edx,4)

        ## constant 1st Hydrogen Z & 2nd hydrogen X forces for 4 j particles 
        movaps nb304_fjzH1(%esp),%xmm0    ## xmm0= fjzH1a   fjzH1b   fjzH1c   fjzH1d 
        movaps nb304_fjxH2(%esp),%xmm2   ## xmm1= fjxH2a  fjxH2b  fjxH2c  fjxH2d
        movlps 20(%edi,%eax,4),%xmm3
        movlps 20(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0       ## xmm0= fjzH1a  fjxH2a  fjzH1b  fjxH2b
        unpckhps %xmm2,%xmm1       ## xmm1= fjzH1c  fjxH2c  fjzH1d  fjxH2d
        movhps 20(%edi,%ebx,4),%xmm3
        movhps 20(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,20(%edi,%eax,4)
        movlps %xmm4,20(%edi,%ecx,4)
        movhps %xmm3,20(%edi,%ebx,4)
        movhps %xmm4,20(%edi,%edx,4)

        ## constant 2nd hydrogen Y & Z forces for 4 j particles 
        movaps nb304_fjyH2(%esp),%xmm0   ## xmm0= fjyH2a  fjyH2b  fjyH2c  fjyH2d 
        movaps nb304_fjzH2(%esp),%xmm2   ## xmm1= fjzH2a  fjzH2b  fjzH2c  fjzH2d
        movlps 28(%edi,%eax,4),%xmm3
        movlps 28(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0            ## xmm0= fjyH2a  fjzH2a  fjyH2b  fjzH2b
        unpckhps %xmm2,%xmm1            ## xmm1= fjyH2c  fjzH2c  fjyH2d  fjzH2d
        movhps 28(%edi,%ebx,4),%xmm3
        movhps 28(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,28(%edi,%eax,4)
        movlps %xmm4,28(%edi,%ecx,4)
        movhps %xmm3,28(%edi,%ebx,4)
        movhps %xmm4,28(%edi,%edx,4)

        ## Dummy (M) X & Y forces for 4 j particles 
        movaps nb304_fjxM(%esp),%xmm0   ## xmm0= fjxMa  fjxMb  fjxMc  fjxMd 
        movaps nb304_fjyM(%esp),%xmm2   ## xmm1= fjyMa  fjyMb  fjyMc  fjyMd
        movlps 36(%edi,%eax,4),%xmm3
        movlps 36(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0            ## xmm0= fjxMa  fjyMa  fjxMb  fjyMb
        unpckhps %xmm2,%xmm1            ## xmm1= fjxMc  fjyMc  fjxMd  fjyMd
        movhps 36(%edi,%ebx,4),%xmm3
        movhps 36(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,36(%edi,%eax,4)
        movlps %xmm4,36(%edi,%ecx,4)
        movhps %xmm3,36(%edi,%ebx,4)
        movhps %xmm4,36(%edi,%edx,4)


        ## Dummy (M) Z forces for 4 j particles 
        ## Just load the four Z coords into one reg. each
        movss 44(%edi,%eax,4),%xmm4
        movss 44(%edi,%ebx,4),%xmm5
        movss 44(%edi,%ecx,4),%xmm6
        movss 44(%edi,%edx,4),%xmm7
        ## add what we have on the stack
        addss nb304_fjzM(%esp),%xmm4
        addss nb304_fjzMb(%esp),%xmm5
        addss nb304_fjzMc(%esp),%xmm6
        addss nb304_fjzMd(%esp),%xmm7
        ## store back
        movss %xmm4,44(%edi,%eax,4)
        movss %xmm5,44(%edi,%ebx,4)
        movss %xmm6,44(%edi,%ecx,4)
        movss %xmm7,44(%edi,%edx,4)

        ## should we do one more iteration? 
        subl $4,nb304_innerk(%esp)
        jl    _nb_kernel304_ia32_sse.nb304_single_check
        jmp   _nb_kernel304_ia32_sse.nb304_unroll_loop
_nb_kernel304_ia32_sse.nb304_single_check: 
        addl $4,nb304_innerk(%esp)
        jnz   _nb_kernel304_ia32_sse.nb304_single_loop
        jmp   _nb_kernel304_ia32_sse.nb304_updateouterdata
_nb_kernel304_ia32_sse.nb304_single_loop: 
        movl  nb304_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb304_innerjjnr(%esp)

        movl nb304_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        xorps %xmm3,%xmm3
        xorps %xmm4,%xmm4
        xorps %xmm5,%xmm5
        movss 36(%esi,%eax,4),%xmm3             ## jxM  -  -  -
        movss 40(%esi,%eax,4),%xmm4             ## jyM  -  -  -
        movss 44(%esi,%eax,4),%xmm5             ## jzM  -  -  -  

        movlps 12(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%esi,%eax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%esi,%eax,4),%xmm2            ## xmm2 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## constant 11011000     ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm2,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movaps  nb304_ixM(%esp),%xmm0
        movaps  nb304_iyM(%esp),%xmm1
        movaps  nb304_izM(%esp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxM   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## constant 11100100    ;# xmm4 = jyM   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## constant 01000100    ;# xmm5 = jzM   0   jzH1 jzH2

        ## store all j coordinates in jM
        movaps %xmm3,nb304_jxM(%esp)
        movaps %xmm4,nb304_jyM(%esp)
        movaps %xmm5,nb304_jzM(%esp)
        subps  %xmm3,%xmm0
        subps  %xmm4,%xmm1
        subps  %xmm5,%xmm2
        movaps %xmm0,nb304_dxMM(%esp)
        movaps %xmm1,nb304_dyMM(%esp)
        movaps %xmm2,nb304_dzMM(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb304_three(%esp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb304_half(%esp),%xmm3   ## rinv iO - j water 

        movaps  %xmm3,%xmm1
        mulps   %xmm0,%xmm1     ## xmm1=r 
        movaps  %xmm3,%xmm0     ## xmm0=rinv 
        mulps  nb304_tsc(%esp),%xmm1

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

        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        movl nb304_VFtab(%ebp),%esi

        movlps (%esi,%ebx,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%ebx,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb304_two(%esp),%xmm7            ## two*Heps2 

        xorps  %xmm3,%xmm3

        ## fetch charges to xmm3 (temporary) 
        movss   nb304_qqMM(%esp),%xmm3
        movhps  nb304_qqMH(%esp),%xmm3

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC 

        addps  nb304_vctot(%esp),%xmm5
        movaps %xmm5,nb304_vctot(%esp)
        xorps  %xmm2,%xmm2
        mulps  nb304_tsc(%esp),%xmm3

        subps  %xmm3,%xmm2
        mulps  %xmm2,%xmm0

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb304_dxMM(%esp),%xmm0
        mulps   nb304_dyMM(%esp),%xmm1
        mulps   nb304_dzMM(%esp),%xmm2
        ## initial update for j forces 
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb304_fjxM(%esp)
        movaps  %xmm4,nb304_fjyM(%esp)
        movaps  %xmm5,nb304_fjzM(%esp)
        addps   nb304_fixM(%esp),%xmm0
        addps   nb304_fiyM(%esp),%xmm1
        addps   nb304_fizM(%esp),%xmm2
        movaps  %xmm0,nb304_fixM(%esp)
        movaps  %xmm1,nb304_fiyM(%esp)
        movaps  %xmm2,nb304_fizM(%esp)


        ## done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb304_ixH1(%esp),%xmm0
        movaps  nb304_iyH1(%esp),%xmm1
        movaps  nb304_izH1(%esp),%xmm2
        movaps  nb304_ixH2(%esp),%xmm3
        movaps  nb304_iyH2(%esp),%xmm4
        movaps  nb304_izH2(%esp),%xmm5
        subps   nb304_jxM(%esp),%xmm0
        subps   nb304_jyM(%esp),%xmm1
        subps   nb304_jzM(%esp),%xmm2
        subps   nb304_jxM(%esp),%xmm3
        subps   nb304_jyM(%esp),%xmm4
        subps   nb304_jzM(%esp),%xmm5
        movaps %xmm0,nb304_dxH1M(%esp)
        movaps %xmm1,nb304_dyH1M(%esp)
        movaps %xmm2,nb304_dzH1M(%esp)
        movaps %xmm3,nb304_dxH2M(%esp)
        movaps %xmm4,nb304_dyH2M(%esp)
        movaps %xmm5,nb304_dzH2M(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm1,%xmm0
        addps %xmm3,%xmm4
        addps %xmm2,%xmm0       ## have rsqH1 in xmm0 
        addps %xmm5,%xmm4       ## have rsqH2 in xmm4 

        ## start with H1, save H2 data 
        movaps %xmm4,nb304_rsqH2M(%esp)

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb304_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb304_half(%esp),%xmm3   ## rinv H1 - j water 
        mulps   nb304_half(%esp),%xmm7   ## rinv H2 - j water  

        ## start with H1, save H2 data 
        movaps %xmm7,nb304_rinvH2M(%esp)

        movaps %xmm3,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=r 
        movaps %xmm3,%xmm0      ## xmm0=rinv 
        mulps  nb304_tsc(%esp),%xmm1

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

        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        movlps (%esi,%ebx,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%ebx,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb304_two(%esp),%xmm7            ## two*Heps2 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb304_qqMH(%esp),%xmm3
        movhps  nb304_qqHH(%esp),%xmm3

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC 
        addps  nb304_vctot(%esp),%xmm5
        movaps %xmm5,nb304_vctot(%esp)

        xorps  %xmm1,%xmm1

        mulps nb304_tsc(%esp),%xmm3
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps  %xmm1,%xmm0
        movaps  %xmm1,%xmm2
        mulps   nb304_dxH1M(%esp),%xmm0
        mulps   nb304_dyH1M(%esp),%xmm1
        mulps   nb304_dzH1M(%esp),%xmm2
        ## update forces H1 - j water 
        movaps  nb304_fjxM(%esp),%xmm3
        movaps  nb304_fjyM(%esp),%xmm4
        movaps  nb304_fjzM(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb304_fjxM(%esp)
        movaps  %xmm4,nb304_fjyM(%esp)
        movaps  %xmm5,nb304_fjzM(%esp)
        addps   nb304_fixH1(%esp),%xmm0
        addps   nb304_fiyH1(%esp),%xmm1
        addps   nb304_fizH1(%esp),%xmm2
        movaps  %xmm0,nb304_fixH1(%esp)
        movaps  %xmm1,nb304_fiyH1(%esp)
        movaps  %xmm2,nb304_fizH1(%esp)
        ## do table for H2 - j water interaction 
        movaps nb304_rinvH2M(%esp),%xmm0
        movaps nb304_rsqH2M(%esp),%xmm1
        mulps  %xmm0,%xmm1      ## xmm0=rinv, xmm1=r 
        mulps  nb304_tsc(%esp),%xmm1

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

        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        movlps (%esi,%ebx,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%ebx,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb304_two(%esp),%xmm7            ## two*Heps2 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb304_qqMH(%esp),%xmm3
        movhps  nb304_qqHH(%esp),%xmm3

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC 
        addps  nb304_vctot(%esp),%xmm5
        movaps %xmm5,nb304_vctot(%esp)

        xorps  %xmm1,%xmm1

        mulps nb304_tsc(%esp),%xmm3
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps  %xmm1,%xmm0
        movaps  %xmm1,%xmm2

        mulps   nb304_dxH2M(%esp),%xmm0
        mulps   nb304_dyH2M(%esp),%xmm1
        mulps   nb304_dzH2M(%esp),%xmm2
        movaps  nb304_fjxM(%esp),%xmm3
        movaps  nb304_fjyM(%esp),%xmm4
        movaps  nb304_fjzM(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movl    nb304_faction(%ebp),%esi
        movaps  %xmm3,nb304_fjxM(%esp)
        movaps  %xmm4,nb304_fjyM(%esp)
        movaps  %xmm5,nb304_fjzM(%esp)
        addps   nb304_fixH2(%esp),%xmm0
        addps   nb304_fiyH2(%esp),%xmm1
        addps   nb304_fizH2(%esp),%xmm2
        movaps  %xmm0,nb304_fixH2(%esp)
        movaps  %xmm1,nb304_fiyH2(%esp)
        movaps  %xmm2,nb304_fizH2(%esp)

        ## update j water forces from local variables 
        movlps  36(%esi,%eax,4),%xmm0
        movlps  12(%esi,%eax,4),%xmm1
        movhps  24(%esi,%eax,4),%xmm1
        movaps  nb304_fjxM(%esp),%xmm3
        movaps  nb304_fjyM(%esp),%xmm4
        movaps  nb304_fjzM(%esp),%xmm5
        movaps  %xmm5,%xmm6
        movaps  %xmm5,%xmm7
        shufps $2,%xmm6,%xmm6 ## constant 00000010
        shufps $3,%xmm7,%xmm7 ## constant 00000011
        addss   44(%esi,%eax,4),%xmm5
        addss   20(%esi,%eax,4),%xmm6
        addss   32(%esi,%eax,4),%xmm7
        movss   %xmm5,44(%esi,%eax,4)
        movss   %xmm6,20(%esi,%eax,4)
        movss   %xmm7,32(%esi,%eax,4)
        movaps   %xmm3,%xmm5
        unpcklps %xmm4,%xmm3
        unpckhps %xmm4,%xmm5
        addps    %xmm3,%xmm0
        addps    %xmm5,%xmm1
        movlps  %xmm0,36(%esi,%eax,4)
        movlps  %xmm1,12(%esi,%eax,4)
        movhps  %xmm1,24(%esi,%eax,4)

        decl nb304_innerk(%esp)
        jz    _nb_kernel304_ia32_sse.nb304_updateouterdata
        jmp   _nb_kernel304_ia32_sse.nb304_single_loop
_nb_kernel304_ia32_sse.nb304_updateouterdata: 
        movl  nb304_ii3(%esp),%ecx
        movl  nb304_faction(%ebp),%edi
        movl  nb304_fshift(%ebp),%esi
        movl  nb304_is3(%esp),%edx

        ## accumulate  H1 i forces in xmm0, xmm1, xmm2 
        movaps nb304_fixH1(%esp),%xmm0
        movaps nb304_fiyH1(%esp),%xmm1
        movaps nb304_fizH1(%esp),%xmm2

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
        movaps %xmm0,%xmm6
        movss %xmm2,%xmm7
        movlhps %xmm1,%xmm6
        shufps $8,%xmm6,%xmm6 ## constant 00001000      

        ## accumulate H2 i forces in xmm0, xmm1, xmm2 
        movaps nb304_fixH2(%esp),%xmm0
        movaps nb304_fiyH2(%esp),%xmm1
        movaps nb304_fizH2(%esp),%xmm2

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

        ## accumulate M i forces in xmm0, xmm1, xmm2 
        movaps nb304_fixM(%esp),%xmm0
        movaps nb304_fiyM(%esp),%xmm1
        movaps nb304_fizM(%esp),%xmm2

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
        movl nb304_n(%esp),%esi
        ## get group index for i particle 
        movl  nb304_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb304_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb304_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb304_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel304_ia32_sse.nb304_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb304_n(%esp)
        jmp _nb_kernel304_ia32_sse.nb304_outer
_nb_kernel304_ia32_sse.nb304_outerend: 
        ## check if more outer neighborlists remain
        movl  nb304_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel304_ia32_sse.nb304_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel304_ia32_sse.nb304_threadloop
_nb_kernel304_ia32_sse.nb304_end: 
        emms

        movl nb304_nouter(%esp),%eax
        movl nb304_ninner(%esp),%ebx
        movl nb304_outeriter(%ebp),%ecx
        movl nb304_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb304_salign(%esp),%eax
        addl %eax,%esp
        addl $1464,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret





.globl nb_kernel304nf_ia32_sse
.globl _nb_kernel304nf_ia32_sse
nb_kernel304nf_ia32_sse:        
_nb_kernel304nf_ia32_sse:       
.set nb304nf_p_nri, 8
.set nb304nf_iinr, 12
.set nb304nf_jindex, 16
.set nb304nf_jjnr, 20
.set nb304nf_shift, 24
.set nb304nf_shiftvec, 28
.set nb304nf_fshift, 32
.set nb304nf_gid, 36
.set nb304nf_pos, 40
.set nb304nf_faction, 44
.set nb304nf_charge, 48
.set nb304nf_p_facel, 52
.set nb304nf_argkrf, 56
.set nb304nf_argcrf, 60
.set nb304nf_Vc, 64
.set nb304nf_type, 68
.set nb304nf_p_ntype, 72
.set nb304nf_vdwparam, 76
.set nb304nf_Vvdw, 80
.set nb304nf_p_tabscale, 84
.set nb304nf_VFtab, 88
.set nb304nf_invsqrta, 92
.set nb304nf_dvda, 96
.set nb304nf_p_gbtabscale, 100
.set nb304nf_GBtab, 104
.set nb304nf_p_nthreads, 108
.set nb304nf_count, 112
.set nb304nf_mtx, 116
.set nb304nf_outeriter, 120
.set nb304nf_inneriter, 124
.set nb304nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb304nf_ixH1, 0
.set nb304nf_iyH1, 16
.set nb304nf_izH1, 32
.set nb304nf_ixH2, 48
.set nb304nf_iyH2, 64
.set nb304nf_izH2, 80
.set nb304nf_ixM, 96
.set nb304nf_iyM, 112
.set nb304nf_izM, 128
.set nb304nf_jxH1, 144
.set nb304nf_jyH1, 160
.set nb304nf_jzH1, 176
.set nb304nf_jxH2, 192
.set nb304nf_jyH2, 208
.set nb304nf_jzH2, 224
.set nb304nf_jxM, 240
.set nb304nf_jyM, 256
.set nb304nf_jzM, 272
.set nb304nf_qqHH, 288
.set nb304nf_qqMH, 304
.set nb304nf_qqMM, 320
.set nb304nf_tsc, 336
.set nb304nf_vctot, 352
.set nb304nf_half, 368
.set nb304nf_three, 384
.set nb304nf_rsqH1H1, 400
.set nb304nf_rsqH1H2, 416
.set nb304nf_rsqH1M, 432
.set nb304nf_rsqH2H1, 448
.set nb304nf_rsqH2H2, 464
.set nb304nf_rsqH2M, 480
.set nb304nf_rsqMH1, 496
.set nb304nf_rsqMH2, 512
.set nb304nf_rsqMM, 528
.set nb304nf_rinvH1H1, 544
.set nb304nf_rinvH1H2, 560
.set nb304nf_rinvH1M, 576
.set nb304nf_rinvH2H1, 592
.set nb304nf_rinvH2H2, 608
.set nb304nf_rinvH2M, 624
.set nb304nf_rinvMH1, 640
.set nb304nf_rinvMH2, 656
.set nb304nf_rinvMM, 672
.set nb304nf_is3, 688
.set nb304nf_ii3, 692
.set nb304nf_innerjjnr, 696
.set nb304nf_innerk, 700
.set nb304nf_n, 704
.set nb304nf_nn1, 708
.set nb304nf_nri, 712
.set nb304nf_nouter, 716
.set nb304nf_ninner, 720
.set nb304nf_salign, 724
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $728,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb304nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb304nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb304nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb304nf_nouter(%esp)
        movl %eax,nb304nf_ninner(%esp)


        movl nb304nf_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb304nf_tsc(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb304nf_half(%esp)
        movss nb304nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb304nf_half(%esp)
        movaps %xmm3,nb304nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb304nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb304nf_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movss 12(%edx,%ebx,4),%xmm5
        movl nb304nf_p_facel(%ebp),%esi
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
        movaps %xmm3,nb304nf_qqHH(%esp)
        movaps %xmm4,nb304nf_qqMH(%esp)
        movaps %xmm5,nb304nf_qqMM(%esp)

_nb_kernel304nf_ia32_sse.nb304nf_threadloop: 
        movl  nb304nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel304nf_ia32_sse.nb304nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel304nf_ia32_sse.nb304nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb304nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb304nf_n(%esp)
        movl %ebx,nb304nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel304nf_ia32_sse.nb304nf_outerstart
        jmp _nb_kernel304nf_ia32_sse.nb304nf_end

_nb_kernel304nf_ia32_sse.nb304nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb304nf_nouter(%esp),%ebx
        movl %ebx,nb304nf_nouter(%esp)

_nb_kernel304nf_ia32_sse.nb304nf_outer: 
        movl  nb304nf_shift(%ebp),%eax          ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb304nf_is3(%esp)            ## store is3 

        movl  nb304nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb304nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb304nf_pos(%ebp),%eax    ## eax = base of pos[]  
        movl  %ebx,nb304nf_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss 12(%eax,%ebx,4),%xmm3
        addss 16(%eax,%ebx,4),%xmm4
        addss 20(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb304nf_ixH1(%esp)
        movaps %xmm4,nb304nf_iyH1(%esp)
        movaps %xmm5,nb304nf_izH1(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 24(%eax,%ebx,4),%xmm0
        addss 28(%eax,%ebx,4),%xmm1
        addss 32(%eax,%ebx,4),%xmm2
        addss 36(%eax,%ebx,4),%xmm3
        addss 40(%eax,%ebx,4),%xmm4
        addss 44(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb304nf_ixH2(%esp)
        movaps %xmm1,nb304nf_iyH2(%esp)
        movaps %xmm2,nb304nf_izH2(%esp)
        movaps %xmm3,nb304nf_ixM(%esp)
        movaps %xmm4,nb304nf_iyM(%esp)
        movaps %xmm5,nb304nf_izM(%esp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb304nf_vctot(%esp)

        movl  nb304nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb304nf_pos(%ebp),%esi
        movl  nb304nf_faction(%ebp),%edi
        movl  nb304nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb304nf_innerjjnr(%esp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb304nf_ninner(%esp),%ecx
        movl  %ecx,nb304nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb304nf_innerk(%esp)         ## number of innerloop atoms 
        jge   _nb_kernel304nf_ia32_sse.nb304nf_unroll_loop
        jmp   _nb_kernel304nf_ia32_sse.nb304nf_single_check
_nb_kernel304nf_ia32_sse.nb304nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb304nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb304nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb304nf_pos(%ebp),%esi     ## base of pos[] 

        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx        ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables 
        movlps 12(%esi,%eax,4),%xmm2
        movlps 24(%esi,%eax,4),%xmm3
        movlps 36(%esi,%eax,4),%xmm4

        movlps 12(%esi,%ebx,4),%xmm5
        movlps 24(%esi,%ebx,4),%xmm6
        movlps 36(%esi,%ebx,4),%xmm7

        movhps 12(%esi,%ecx,4),%xmm2
        movhps 24(%esi,%ecx,4),%xmm3
        movhps 36(%esi,%ecx,4),%xmm4

        movhps 12(%esi,%edx,4),%xmm5
        movhps 24(%esi,%edx,4),%xmm6
        movhps 36(%esi,%edx,4),%xmm7

        movaps %xmm2,%xmm0
        movaps %xmm3,%xmm1
        unpcklps %xmm5,%xmm0
        unpcklps %xmm6,%xmm1
        unpckhps %xmm5,%xmm2
        unpckhps %xmm6,%xmm3
        movaps %xmm4,%xmm5
        movaps   %xmm0,%xmm6
        unpcklps %xmm7,%xmm4
        unpckhps %xmm7,%xmm5
        movaps   %xmm1,%xmm7
        movlhps  %xmm2,%xmm0
        movaps %xmm0,nb304nf_jxH1(%esp)
        movhlps  %xmm6,%xmm2
        movaps %xmm2,nb304nf_jyH1(%esp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb304nf_jxH2(%esp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb304nf_jyH2(%esp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb304nf_jxM(%esp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb304nf_jyM(%esp)

        movss  20(%esi,%eax,4),%xmm0
        movss  32(%esi,%eax,4),%xmm1
        movss  44(%esi,%eax,4),%xmm2

        movss  20(%esi,%ecx,4),%xmm3
        movss  32(%esi,%ecx,4),%xmm4
        movss  44(%esi,%ecx,4),%xmm5

        movhps 16(%esi,%ebx,4),%xmm0
        movhps 28(%esi,%ebx,4),%xmm1
        movhps 40(%esi,%ebx,4),%xmm2

        movhps 16(%esi,%edx,4),%xmm3
        movhps 28(%esi,%edx,4),%xmm4
        movhps 40(%esi,%edx,4),%xmm5

        shufps $204,%xmm3,%xmm0 ## constant 11001100
        shufps $204,%xmm4,%xmm1 ## constant 11001100
        shufps $204,%xmm5,%xmm2 ## constant 11001100
        movaps %xmm0,nb304nf_jzH1(%esp)
        movaps %xmm1,nb304nf_jzH2(%esp)
        movaps %xmm2,nb304nf_jzM(%esp)

        movaps nb304nf_ixH1(%esp),%xmm0
        movaps nb304nf_iyH1(%esp),%xmm1
        movaps nb304nf_izH1(%esp),%xmm2
        movaps nb304nf_ixH1(%esp),%xmm3
        movaps nb304nf_iyH1(%esp),%xmm4
        movaps nb304nf_izH1(%esp),%xmm5
        subps  nb304nf_jxH1(%esp),%xmm0
        subps  nb304nf_jyH1(%esp),%xmm1
        subps  nb304nf_jzH1(%esp),%xmm2
        subps  nb304nf_jxH2(%esp),%xmm3
        subps  nb304nf_jyH2(%esp),%xmm4
        subps  nb304nf_jzH2(%esp),%xmm5
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
        movaps %xmm0,nb304nf_rsqH1H1(%esp)
        movaps %xmm3,nb304nf_rsqH1H2(%esp)

        movaps nb304nf_ixH1(%esp),%xmm0
        movaps nb304nf_iyH1(%esp),%xmm1
        movaps nb304nf_izH1(%esp),%xmm2
        movaps nb304nf_ixH2(%esp),%xmm3
        movaps nb304nf_iyH2(%esp),%xmm4
        movaps nb304nf_izH2(%esp),%xmm5
        subps  nb304nf_jxM(%esp),%xmm0
        subps  nb304nf_jyM(%esp),%xmm1
        subps  nb304nf_jzM(%esp),%xmm2
        subps  nb304nf_jxH1(%esp),%xmm3
        subps  nb304nf_jyH1(%esp),%xmm4
        subps  nb304nf_jzH1(%esp),%xmm5
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
        movaps %xmm0,nb304nf_rsqH1M(%esp)
        movaps %xmm3,nb304nf_rsqH2H1(%esp)

        movaps nb304nf_ixH2(%esp),%xmm0
        movaps nb304nf_iyH2(%esp),%xmm1
        movaps nb304nf_izH2(%esp),%xmm2
        movaps nb304nf_ixH2(%esp),%xmm3
        movaps nb304nf_iyH2(%esp),%xmm4
        movaps nb304nf_izH2(%esp),%xmm5
        subps  nb304nf_jxH2(%esp),%xmm0
        subps  nb304nf_jyH2(%esp),%xmm1
        subps  nb304nf_jzH2(%esp),%xmm2
        subps  nb304nf_jxM(%esp),%xmm3
        subps  nb304nf_jyM(%esp),%xmm4
        subps  nb304nf_jzM(%esp),%xmm5
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
        movaps %xmm0,nb304nf_rsqH2H2(%esp)
        movaps %xmm3,nb304nf_rsqH2M(%esp)

        movaps nb304nf_ixM(%esp),%xmm0
        movaps nb304nf_iyM(%esp),%xmm1
        movaps nb304nf_izM(%esp),%xmm2
        movaps nb304nf_ixM(%esp),%xmm3
        movaps nb304nf_iyM(%esp),%xmm4
        movaps nb304nf_izM(%esp),%xmm5
        subps  nb304nf_jxH1(%esp),%xmm0
        subps  nb304nf_jyH1(%esp),%xmm1
        subps  nb304nf_jzH1(%esp),%xmm2
        subps  nb304nf_jxH2(%esp),%xmm3
        subps  nb304nf_jyH2(%esp),%xmm4
        subps  nb304nf_jzH2(%esp),%xmm5
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
        movaps %xmm0,nb304nf_rsqMH1(%esp)
        movaps %xmm4,nb304nf_rsqMH2(%esp)

        movaps nb304nf_ixM(%esp),%xmm0
        movaps nb304nf_iyM(%esp),%xmm1
        movaps nb304nf_izM(%esp),%xmm2
        subps  nb304nf_jxM(%esp),%xmm0
        subps  nb304nf_jyM(%esp),%xmm1
        subps  nb304nf_jzM(%esp),%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb304nf_rsqMM(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb304nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb304nf_half(%esp),%xmm3   ## rinvMM
        mulps   nb304nf_half(%esp),%xmm7   ## rinvMH2 
        movaps  %xmm3,nb304nf_rinvMM(%esp)
        movaps  %xmm7,nb304nf_rinvMH2(%esp)

        rsqrtps nb304nf_rsqH1H1(%esp),%xmm1
        rsqrtps nb304nf_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb304nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb304nf_rsqH1H1(%esp),%xmm1
        mulps   nb304nf_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb304nf_half(%esp),%xmm3
        mulps   nb304nf_half(%esp),%xmm7
        movaps  %xmm3,nb304nf_rinvH1H1(%esp)
        movaps  %xmm7,nb304nf_rinvH1H2(%esp)

        rsqrtps nb304nf_rsqH1M(%esp),%xmm1
        rsqrtps nb304nf_rsqH2H1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb304nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb304nf_rsqH1M(%esp),%xmm1
        mulps   nb304nf_rsqH2H1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb304nf_half(%esp),%xmm3
        mulps   nb304nf_half(%esp),%xmm7
        movaps  %xmm3,nb304nf_rinvH1M(%esp)
        movaps  %xmm7,nb304nf_rinvH2H1(%esp)

        rsqrtps nb304nf_rsqH2H2(%esp),%xmm1
        rsqrtps nb304nf_rsqH2M(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb304nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb304nf_rsqH2H2(%esp),%xmm1
        mulps   nb304nf_rsqH2M(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb304nf_half(%esp),%xmm3
        mulps   nb304nf_half(%esp),%xmm7
        movaps  %xmm3,nb304nf_rinvH2H2(%esp)
        movaps  %xmm7,nb304nf_rinvH2M(%esp)

        rsqrtps nb304nf_rsqMH1(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb304nf_three(%esp),%xmm3
        mulps   nb304nf_rsqMH1(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb304nf_half(%esp),%xmm3
        movaps  %xmm3,nb304nf_rinvMH1(%esp)

        ## start with H1-H1 interaction 
        movaps nb304nf_rinvH1H1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304nf_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulps  nb304nf_tsc(%esp),%xmm1

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

        movl nb304nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb304nf_qqHH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## update vctot 
        addps  nb304nf_vctot(%esp),%xmm5
        movaps %xmm5,nb304nf_vctot(%esp)

        ## H1-H2 interaction 
        movaps nb304nf_rinvH1H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304nf_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulps  nb304nf_tsc(%esp),%xmm1
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

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb304nf_qqHH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb304nf_vctot(%esp),%xmm5
        movaps %xmm5,nb304nf_vctot(%esp)

        ## H1-M interaction  
        movaps nb304nf_rinvH1M(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304nf_rsqH1M(%esp),%xmm1   ## xmm1=r 
        mulps  nb304nf_tsc(%esp),%xmm1
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

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb304nf_qqMH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb304nf_vctot(%esp),%xmm5
        movaps %xmm5,nb304nf_vctot(%esp)

        ## H2-H1 interaction 
        movaps nb304nf_rinvH2H1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304nf_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulps  nb304nf_tsc(%esp),%xmm1
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

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb304nf_qqHH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb304nf_vctot(%esp),%xmm5
        movaps %xmm5,nb304nf_vctot(%esp)

        ## H2-H2 interaction 
        movaps nb304nf_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304nf_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulps  nb304nf_tsc(%esp),%xmm1
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

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb304nf_qqHH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb304nf_vctot(%esp),%xmm5
        movaps %xmm5,nb304nf_vctot(%esp)

        ## H2-M interaction 
        movaps nb304nf_rinvH2M(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304nf_rsqH2M(%esp),%xmm1   ## xmm1=r 
        mulps  nb304nf_tsc(%esp),%xmm1
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

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb304nf_qqMH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        addps  nb304nf_vctot(%esp),%xmm5
        movaps %xmm5,nb304nf_vctot(%esp)

        ## M-H1 interaction 
        movaps nb304nf_rinvMH1(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304nf_rsqMH1(%esp),%xmm1   ## xmm1=r 
        mulps  nb304nf_tsc(%esp),%xmm1
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

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb304nf_qqMH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb304nf_vctot(%esp),%xmm5
        movaps %xmm5,nb304nf_vctot(%esp)

        ## M-H2 interaction 
        movaps nb304nf_rinvMH2(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304nf_rsqMH2(%esp),%xmm1   ## xmm1=r 
        mulps  nb304nf_tsc(%esp),%xmm1
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

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb304nf_qqMH(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb304nf_vctot(%esp),%xmm5
        movaps %xmm5,nb304nf_vctot(%esp)

        ## M-M interaction 
        movaps nb304nf_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304nf_rsqMM(%esp),%xmm1   ## xmm1=r 
        mulps  nb304nf_tsc(%esp),%xmm1
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

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb304nf_qqMM(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb304nf_vctot(%esp),%xmm5
        movaps %xmm5,nb304nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb304nf_innerk(%esp)
        jl    _nb_kernel304nf_ia32_sse.nb304nf_single_check
        jmp   _nb_kernel304nf_ia32_sse.nb304nf_unroll_loop
_nb_kernel304nf_ia32_sse.nb304nf_single_check: 
        addl $4,nb304nf_innerk(%esp)
        jnz   _nb_kernel304nf_ia32_sse.nb304nf_single_loop
        jmp   _nb_kernel304nf_ia32_sse.nb304nf_updateouterdata
_nb_kernel304nf_ia32_sse.nb304nf_single_loop: 
        movl  nb304nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb304nf_innerjjnr(%esp)

        movl nb304nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        xorps %xmm3,%xmm3
        xorps %xmm4,%xmm4
        xorps %xmm5,%xmm5
        movss 36(%esi,%eax,4),%xmm3             ## jxM  -  -  -
        movss 40(%esi,%eax,4),%xmm4             ## jyM  -  -  -
        movss 44(%esi,%eax,4),%xmm5             ## jzM  -  -  -  

        movlps 12(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%esi,%eax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%esi,%eax,4),%xmm2            ## xmm2 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## constant 11011000     ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm2,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movaps  nb304nf_ixM(%esp),%xmm0
        movaps  nb304nf_iyM(%esp),%xmm1
        movaps  nb304nf_izM(%esp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxM   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## constant 11100100    ;# xmm4 = jyM   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## constant 01000100    ;# xmm5 = jzM   0   jzH1 jzH2

        ## store all j coordinates in jM
        movaps %xmm3,nb304nf_jxM(%esp)
        movaps %xmm4,nb304nf_jyM(%esp)
        movaps %xmm5,nb304nf_jzM(%esp)
        subps  %xmm3,%xmm0
        subps  %xmm4,%xmm1
        subps  %xmm5,%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb304nf_three(%esp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb304nf_half(%esp),%xmm3   ## rinv iO - j water 

        movaps  %xmm3,%xmm1
        mulps   %xmm0,%xmm1     ## xmm1=r 
        movaps  %xmm3,%xmm0     ## xmm0=rinv 
        mulps  nb304nf_tsc(%esp),%xmm1

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

        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        movl nb304nf_VFtab(%ebp),%esi

        movlps (%esi,%ebx,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%ebx,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 

        xorps  %xmm3,%xmm3

        ## fetch charges to xmm3 (temporary) 
        movss   nb304nf_qqMM(%esp),%xmm3
        movhps  nb304nf_qqMH(%esp),%xmm3

        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 

        addps  nb304nf_vctot(%esp),%xmm5
        movaps %xmm5,nb304nf_vctot(%esp)

        ## done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb304nf_ixH1(%esp),%xmm0
        movaps  nb304nf_iyH1(%esp),%xmm1
        movaps  nb304nf_izH1(%esp),%xmm2
        movaps  nb304nf_ixH2(%esp),%xmm3
        movaps  nb304nf_iyH2(%esp),%xmm4
        movaps  nb304nf_izH2(%esp),%xmm5
        subps   nb304nf_jxM(%esp),%xmm0
        subps   nb304nf_jyM(%esp),%xmm1
        subps   nb304nf_jzM(%esp),%xmm2
        subps   nb304nf_jxM(%esp),%xmm3
        subps   nb304nf_jyM(%esp),%xmm4
        subps   nb304nf_jzM(%esp),%xmm5
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm1,%xmm0
        addps %xmm3,%xmm4
        addps %xmm2,%xmm0       ## have rsqH1 in xmm0 
        addps %xmm5,%xmm4       ## have rsqH2 in xmm4 

        ## start with H1, save H2 data 
        movaps %xmm4,nb304nf_rsqH2M(%esp)

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb304nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb304nf_half(%esp),%xmm3   ## rinv H1 - j water 
        mulps   nb304nf_half(%esp),%xmm7   ## rinv H2 - j water  

        ## start with H1, save H2 data 
        movaps %xmm7,nb304nf_rinvH2M(%esp)

        movaps %xmm3,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=r 
        movaps %xmm3,%xmm0      ## xmm0=rinv 
        mulps  nb304nf_tsc(%esp),%xmm1

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

        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        movlps (%esi,%ebx,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%ebx,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb304nf_qqMH(%esp),%xmm3
        movhps  nb304nf_qqHH(%esp),%xmm3

        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 
        addps  nb304nf_vctot(%esp),%xmm5
        movaps %xmm5,nb304nf_vctot(%esp)

        ## do table for H2 - j water interaction 
        movaps nb304nf_rinvH2M(%esp),%xmm0
        movaps nb304nf_rsqH2M(%esp),%xmm1
        mulps  %xmm0,%xmm1      ## xmm0=rinv, xmm1=r 
        mulps  nb304nf_tsc(%esp),%xmm1

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

        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        movlps (%esi,%ebx,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%ebx,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7  
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb304nf_qqMH(%esp),%xmm3
        movhps  nb304nf_qqHH(%esp),%xmm3

        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 
        addps  nb304nf_vctot(%esp),%xmm5
        movaps %xmm5,nb304nf_vctot(%esp)

        decl nb304nf_innerk(%esp)
        jz    _nb_kernel304nf_ia32_sse.nb304nf_updateouterdata
        jmp   _nb_kernel304nf_ia32_sse.nb304nf_single_loop
_nb_kernel304nf_ia32_sse.nb304nf_updateouterdata: 
        ## get n from stack
        movl nb304nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb304nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb304nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb304nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb304nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel304nf_ia32_sse.nb304nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb304nf_n(%esp)
        jmp _nb_kernel304nf_ia32_sse.nb304nf_outer
_nb_kernel304nf_ia32_sse.nb304nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb304nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel304nf_ia32_sse.nb304nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel304nf_ia32_sse.nb304nf_threadloop
_nb_kernel304nf_ia32_sse.nb304nf_end: 
        emms

        movl nb304nf_nouter(%esp),%eax
        movl nb304nf_ninner(%esp),%ebx
        movl nb304nf_outeriter(%ebp),%ecx
        movl nb304nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb304nf_salign(%esp),%eax
        addl %eax,%esp
        addl $728,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


