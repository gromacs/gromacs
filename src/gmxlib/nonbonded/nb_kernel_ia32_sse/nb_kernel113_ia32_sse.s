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



.globl nb_kernel113_ia32_sse
.globl _nb_kernel113_ia32_sse
nb_kernel113_ia32_sse:  
_nb_kernel113_ia32_sse: 
.set nb113_p_nri, 8
.set nb113_iinr, 12
.set nb113_jindex, 16
.set nb113_jjnr, 20
.set nb113_shift, 24
.set nb113_shiftvec, 28
.set nb113_fshift, 32
.set nb113_gid, 36
.set nb113_pos, 40
.set nb113_faction, 44
.set nb113_charge, 48
.set nb113_p_facel, 52
.set nb113_p_krf, 56
.set nb113_p_crf, 60
.set nb113_Vc, 64
.set nb113_type, 68
.set nb113_p_ntype, 72
.set nb113_vdwparam, 76
.set nb113_Vvdw, 80
.set nb113_p_tabscale, 84
.set nb113_VFtab, 88
.set nb113_invsqrta, 92
.set nb113_dvda, 96
.set nb113_p_gbtabscale, 100
.set nb113_GBtab, 104
.set nb113_p_nthreads, 108
.set nb113_count, 112
.set nb113_mtx, 116
.set nb113_outeriter, 120
.set nb113_inneriter, 124
.set nb113_work, 128
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
.set nb113_is3, 880
.set nb113_ii3, 884
.set nb113_ntia, 888
.set nb113_innerjjnr, 892
.set nb113_innerk, 896
.set nb113_n, 900
.set nb113_nn1, 904
.set nb113_nri, 908
.set nb113_nouter, 912
.set nb113_ninner, 916
.set nb113_salign, 920
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $924,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb113_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb113_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb113_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb113_nouter(%esp)
        movl %eax,nb113_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb113_half(%esp)
        movss nb113_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm3,%xmm4
        addps  %xmm4,%xmm4      ## 6.0
        movaps %xmm4,%xmm5
        addps  %xmm5,%xmm5      ## constant 12.0
        movaps %xmm1,nb113_half(%esp)
        movaps %xmm2,nb113_two(%esp)
        movaps %xmm3,nb113_three(%esp)
        movaps %xmm4,nb113_six(%esp)
        movaps %xmm5,nb113_twelve(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb113_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb113_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm4
        movss 12(%edx,%ebx,4),%xmm3
        movl nb113_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb113_iqM(%esp)
        movaps %xmm4,nb113_iqH(%esp)

        movl  nb113_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb113_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb113_ntia(%esp)

_nb_kernel113_ia32_sse.nb113_threadloop: 
        movl  nb113_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel113_ia32_sse.nb113_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel113_ia32_sse.nb113_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb113_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb113_n(%esp)
        movl %ebx,nb113_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel113_ia32_sse.nb113_outerstart
        jmp _nb_kernel113_ia32_sse.nb113_end

_nb_kernel113_ia32_sse.nb113_outerstart: 
        ## ebx contains number of outer iterations
        addl nb113_nouter(%esp),%ebx
        movl %ebx,nb113_nouter(%esp)

_nb_kernel113_ia32_sse.nb113_outer: 
        movl  nb113_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb113_is3(%esp)      ## store is3 

        movl  nb113_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb113_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb113_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb113_ii3(%esp)

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
        movaps %xmm3,nb113_ixO(%esp)
        movaps %xmm4,nb113_iyO(%esp)
        movaps %xmm5,nb113_izO(%esp)
        movaps %xmm6,nb113_ixH1(%esp)
        movaps %xmm7,nb113_iyH1(%esp)

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
        movaps %xmm6,nb113_izH1(%esp)
        movaps %xmm0,nb113_ixH2(%esp)
        movaps %xmm1,nb113_iyH2(%esp)
        movaps %xmm2,nb113_izH2(%esp)
        movaps %xmm3,nb113_ixM(%esp)
        movaps %xmm4,nb113_iyM(%esp)
        movaps %xmm5,nb113_izM(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb113_vctot(%esp)
        movaps %xmm4,nb113_Vvdwtot(%esp)
        movaps %xmm4,nb113_fixO(%esp)
        movaps %xmm4,nb113_fiyO(%esp)
        movaps %xmm4,nb113_fizO(%esp)
        movaps %xmm4,nb113_fixH1(%esp)
        movaps %xmm4,nb113_fiyH1(%esp)
        movaps %xmm4,nb113_fizH1(%esp)
        movaps %xmm4,nb113_fixH2(%esp)
        movaps %xmm4,nb113_fiyH2(%esp)
        movaps %xmm4,nb113_fizH2(%esp)
        movaps %xmm4,nb113_fixM(%esp)
        movaps %xmm4,nb113_fiyM(%esp)
        movaps %xmm4,nb113_fizM(%esp)

        movl  nb113_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb113_pos(%ebp),%esi
        movl  nb113_faction(%ebp),%edi
        movl  nb113_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb113_innerjjnr(%esp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb113_ninner(%esp),%ecx
        movl  %ecx,nb113_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb113_innerk(%esp)   ## number of innerloop atoms 
        jge   _nb_kernel113_ia32_sse.nb113_unroll_loop
        jmp   _nb_kernel113_ia32_sse.nb113_odd_inner
_nb_kernel113_ia32_sse.nb113_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb113_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb113_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb113_charge(%ebp),%esi    ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb113_iqM(%esp),%xmm3
        mulps  nb113_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb113_qqM(%esp)
        movaps  %xmm4,nb113_qqH(%esp)

        movl nb113_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb113_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb113_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx
        addl %edi,%ecx
        addl %edi,%edx

        movlps (%esi,%eax,4),%xmm6
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm6
        movhps (%esi,%edx,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm6 ## constant 11011101

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        movaps %xmm4,nb113_c6(%esp)
        movaps %xmm6,nb113_c12(%esp)

        movl nb113_pos(%ebp),%esi       ## base of pos[] 

        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx        ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move four coordinates to xmm0-xmm2   
        movlps (%esi,%eax,4),%xmm4
        movlps (%esi,%ecx,4),%xmm5
        movss 8(%esi,%eax,4),%xmm2
        movss 8(%esi,%ecx,4),%xmm6

        movhps (%esi,%ebx,4),%xmm4
        movhps (%esi,%edx,4),%xmm5

        movss 8(%esi,%ebx,4),%xmm0
        movss 8(%esi,%edx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## constant 10001000

        shufps $136,%xmm5,%xmm0 ## constant 10001000
        shufps $221,%xmm5,%xmm1 ## constant 11011101            

        ## move ixO-izO to xmm4-xmm6 
        movaps nb113_ixO(%esp),%xmm4
        movaps nb113_iyO(%esp),%xmm5
        movaps nb113_izO(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb113_dxO(%esp)
        movaps %xmm5,nb113_dyO(%esp)
        movaps %xmm6,nb113_dzO(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        movaps %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb113_ixH1(%esp),%xmm4
        movaps nb113_iyH1(%esp),%xmm5
        movaps nb113_izH1(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb113_dxH1(%esp)
        movaps %xmm5,nb113_dyH1(%esp)
        movaps %xmm6,nb113_dzH1(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movaps nb113_ixH2(%esp),%xmm3
        movaps nb113_iyH2(%esp),%xmm4
        movaps nb113_izH2(%esp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## store dr 
        movaps %xmm3,nb113_dxH2(%esp)
        movaps %xmm4,nb113_dyH2(%esp)
        movaps %xmm5,nb113_dzH2(%esp)
        ## square it 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm4,%xmm5
        addps %xmm3,%xmm5

        ## move ixM-izM to xmm2-xmm4  
        movaps nb113_iyM(%esp),%xmm3
        movaps nb113_izM(%esp),%xmm4
        subps  %xmm1,%xmm3
        subps  %xmm2,%xmm4
        movaps nb113_ixM(%esp),%xmm2
        subps  %xmm0,%xmm2

        ## store dr 
        movaps %xmm2,nb113_dxM(%esp)
        movaps %xmm3,nb113_dyM(%esp)
        movaps %xmm4,nb113_dzM(%esp)
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
        movaps  nb113_three(%esp),%xmm0
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb113_half(%esp),%xmm0
        movaps  %xmm0,nb113_rinvH1(%esp)        ## rinvH1 

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb113_three(%esp),%xmm0
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb113_half(%esp),%xmm0
        movaps  %xmm0,nb113_rinvH2(%esp)        ## rinvH2 

        ## rsqM - seed to xmm2 
        rsqrtps %xmm4,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb113_three(%esp),%xmm0
        mulps   %xmm4,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb113_half(%esp),%xmm0
        movaps  %xmm0,nb113_rinvM(%esp)

        ## Do the O LJ-only interaction directly.       
        rcpps   %xmm7,%xmm2
        movaps  nb113_two(%esp),%xmm1
        mulps   %xmm2,%xmm7
        subps   %xmm7,%xmm1
        mulps   %xmm1,%xmm2 ## rinvsq 
        movaps  %xmm2,%xmm0
        mulps   %xmm2,%xmm0     ## r4
        mulps   %xmm2,%xmm0     ## r6
        movaps  %xmm0,%xmm1
        mulps   %xmm1,%xmm1     ## r12
        mulps   nb113_c6(%esp),%xmm0
        mulps   nb113_c12(%esp),%xmm1
        movaps  %xmm1,%xmm3
        subps   %xmm0,%xmm3     ## Vvdw12-Vvdw6
        addps   nb113_Vvdwtot(%esp),%xmm3
        movaps  %xmm3,nb113_Vvdwtot(%esp)
        mulps   nb113_six(%esp),%xmm0
        mulps   nb113_twelve(%esp),%xmm1
        subps   %xmm0,%xmm1
        mulps   %xmm2,%xmm1     ## fscal
        movaps nb113_dxO(%esp),%xmm3
        movaps nb113_dyO(%esp),%xmm4
        movaps nb113_dzO(%esp),%xmm5
        mulps  %xmm1,%xmm3
        mulps  %xmm1,%xmm4
        mulps  %xmm1,%xmm5      ## tx in xmm3-xmm5

        ## update O forces 
        movaps nb113_fixO(%esp),%xmm0
        movaps nb113_fiyO(%esp),%xmm1
        movaps nb113_fizO(%esp),%xmm2
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2
        movaps %xmm0,nb113_fixO(%esp)
        movaps %xmm1,nb113_fiyO(%esp)
        movaps %xmm2,nb113_fizO(%esp)
        ## update j forces with water O 
        movaps %xmm3,nb113_fjx(%esp)
        movaps %xmm4,nb113_fjy(%esp)
        movaps %xmm5,nb113_fjz(%esp)

        ## Do H1 interaction
        movaps  nb113_rinvH1(%esp),%xmm7
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        mulps  nb113_qqH(%esp),%xmm7    ## xmm7=vcoul 
        mulps  %xmm7,%xmm4      ## total fsH1 in xmm4 

        addps  nb113_vctot(%esp),%xmm7
        movaps %xmm7,nb113_vctot(%esp)

        movaps nb113_dxH1(%esp),%xmm0
        movaps nb113_dyH1(%esp),%xmm1
        movaps nb113_dzH1(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H1 forces 
        movaps nb113_fixH1(%esp),%xmm3
        movaps nb113_fiyH1(%esp),%xmm4
        movaps nb113_fizH1(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb113_fixH1(%esp)
        movaps %xmm4,nb113_fiyH1(%esp)
        movaps %xmm7,nb113_fizH1(%esp)
        ## update j forces with water H1 
        addps  nb113_fjx(%esp),%xmm0
        addps  nb113_fjy(%esp),%xmm1
        addps  nb113_fjz(%esp),%xmm2
        movaps %xmm0,nb113_fjx(%esp)
        movaps %xmm1,nb113_fjy(%esp)
        movaps %xmm2,nb113_fjz(%esp)

        ## Done with H1, do H2 interactions
        movaps  nb113_rinvH2(%esp),%xmm7
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        mulps  nb113_qqH(%esp),%xmm7    ## xmm7=vcoul 
        mulps  %xmm7,%xmm4      ## total fsH2 in xmm4 

        addps  nb113_vctot(%esp),%xmm7
        movaps %xmm7,nb113_vctot(%esp)

        movaps nb113_dxH2(%esp),%xmm0
        movaps nb113_dyH2(%esp),%xmm1
        movaps nb113_dzH2(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update H2 forces 
        movaps nb113_fixH2(%esp),%xmm3
        movaps nb113_fiyH2(%esp),%xmm4
        movaps nb113_fizH2(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb113_fixH2(%esp)
        movaps %xmm4,nb113_fiyH2(%esp)
        movaps %xmm7,nb113_fizH2(%esp)
        addps nb113_fjx(%esp),%xmm0
        addps nb113_fjy(%esp),%xmm1
        addps nb113_fjz(%esp),%xmm2
        movaps %xmm0,nb113_fjx(%esp)
        movaps %xmm1,nb113_fjy(%esp)
        movaps %xmm2,nb113_fjz(%esp)

        ## Done with H2, do M interactions
        movaps  nb113_rinvM(%esp),%xmm7
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        mulps  nb113_qqM(%esp),%xmm7    ## xmm7=vcoul 
        mulps  %xmm7,%xmm4      ## total fsM in xmm4 

        addps  nb113_vctot(%esp),%xmm7
        movaps %xmm7,nb113_vctot(%esp)

        movaps nb113_dxM(%esp),%xmm0
        movaps nb113_dyM(%esp),%xmm1
        movaps nb113_dzM(%esp),%xmm2
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2

        ## update M forces 
        movaps nb113_fixM(%esp),%xmm3
        movaps nb113_fiyM(%esp),%xmm4
        movaps nb113_fizM(%esp),%xmm7
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm7
        movaps %xmm3,nb113_fixM(%esp)
        movaps %xmm4,nb113_fiyM(%esp)
        movaps %xmm7,nb113_fizM(%esp)

        movl nb113_faction(%ebp),%edi
        ## update j forces from stored values
        addps nb113_fjx(%esp),%xmm0
        addps nb113_fjy(%esp),%xmm1
        addps nb113_fjz(%esp),%xmm2

        movlps (%edi,%eax,4),%xmm4
        movlps (%edi,%ecx,4),%xmm7
        movhps (%edi,%ebx,4),%xmm4
        movhps (%edi,%edx,4),%xmm7

        movaps %xmm4,%xmm3
        shufps $136,%xmm7,%xmm3 ## constant 10001000
        shufps $221,%xmm7,%xmm4 ## constant 11011101

        ## xmm3 has fjx, xmm4 has fjy 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        ## unpack them back for storing 
        movaps %xmm3,%xmm7
        unpcklps %xmm4,%xmm7
        unpckhps %xmm4,%xmm3
        movlps %xmm7,(%edi,%eax,4)
        movlps %xmm3,(%edi,%ecx,4)
        movhps %xmm7,(%edi,%ebx,4)
        movhps %xmm3,(%edi,%edx,4)
        ## finally z forces 
        movss  8(%edi,%eax,4),%xmm0
        movss  8(%edi,%ebx,4),%xmm1
        movss  8(%edi,%ecx,4),%xmm3
        movss  8(%edi,%edx,4),%xmm4
        subss  %xmm2,%xmm0
        shufps $229,%xmm2,%xmm2 ## constant 11100101
        subss  %xmm2,%xmm1
        shufps $234,%xmm2,%xmm2 ## constant 11101010
        subss  %xmm2,%xmm3
        shufps $255,%xmm2,%xmm2 ## constant 11111111
        subss  %xmm2,%xmm4
        movss  %xmm0,8(%edi,%eax,4)
        movss  %xmm1,8(%edi,%ebx,4)
        movss  %xmm3,8(%edi,%ecx,4)
        movss  %xmm4,8(%edi,%edx,4)

        ## should we do one more iteration? 
        subl $4,nb113_innerk(%esp)
        jl    _nb_kernel113_ia32_sse.nb113_odd_inner
        jmp   _nb_kernel113_ia32_sse.nb113_unroll_loop
_nb_kernel113_ia32_sse.nb113_odd_inner: 
        addl $4,nb113_innerk(%esp)
        jnz   _nb_kernel113_ia32_sse.nb113_odd_loop
        jmp   _nb_kernel113_ia32_sse.nb113_updateouterdata
_nb_kernel113_ia32_sse.nb113_odd_loop: 
        movl  nb113_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb113_innerjjnr(%esp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb113_iqM(%esp),%xmm4
        movl nb113_charge(%ebp),%esi
        movhps nb113_iqH(%esp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb113_qqM(%esp)    ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movl nb113_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb113_vdwparam(%ebp),%esi
        shll %ebx
        addl nb113_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb113_c6(%esp)
        movaps %xmm7,nb113_c12(%esp)

        movl nb113_pos(%ebp),%esi
        leal (%eax,%eax,2),%eax

        movss nb113_ixO(%esp),%xmm3
        movss nb113_iyO(%esp),%xmm4
        movss nb113_izO(%esp),%xmm5
        movss nb113_ixH1(%esp),%xmm0
        movss nb113_iyH1(%esp),%xmm1
        movss nb113_izH1(%esp),%xmm2
        unpcklps nb113_ixH2(%esp),%xmm3         ## ixO ixH2 - -
        unpcklps nb113_iyH2(%esp),%xmm4         ## iyO iyH2 - -
        unpcklps nb113_izH2(%esp),%xmm5         ## izO izH2 - -
        unpcklps nb113_ixM(%esp),%xmm0          ## ixH1 ixM - -
        unpcklps nb113_iyM(%esp),%xmm1          ## iyH1 iyM - -
        unpcklps nb113_izM(%esp),%xmm2          ## izH1 izM - -
        unpcklps %xmm0,%xmm3    ## ixO ixH1 ixH2 ixM
        unpcklps %xmm1,%xmm4    ## same for y
        unpcklps %xmm2,%xmm5    ## same for z

        ## move j coords to xmm0-xmm2 
        movss (%esi,%eax,4),%xmm0
        movss 4(%esi,%eax,4),%xmm1
        movss 8(%esi,%eax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## use O distances for storage
        movaps %xmm3,nb113_dxO(%esp)
        movaps %xmm4,nb113_dyO(%esp)
        movaps %xmm5,nb113_dzO(%esp)

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
        movaps nb113_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb113_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv

        movaps %xmm0,%xmm1
        mulps %xmm1,%xmm1       ## rinvsq 
        movaps %xmm0,%xmm7
        mulps  nb113_qqM(%esp),%xmm7   ## vcoul
        movaps %xmm7,%xmm6

        addps  nb113_vctot(%esp),%xmm7
        movaps %xmm7,nb113_vctot(%esp)

        movaps %xmm1,%xmm2
        mulss  %xmm1,%xmm1
        mulss  %xmm2,%xmm1      ## xmm1=rinvsix
        xorps  %xmm4,%xmm4
        movss  %xmm1,%xmm4
        mulss  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulss  nb113_c6(%esp),%xmm1
        mulss  nb113_c12(%esp),%xmm4
        movaps %xmm4,%xmm3
        subss  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        mulss  nb113_six(%esp),%xmm1
        mulss  nb113_twelve(%esp),%xmm4
        subss  %xmm1,%xmm4
        addss  nb113_Vvdwtot(%esp),%xmm3
        movss  %xmm3,nb113_Vvdwtot(%esp)
        addps  %xmm6,%xmm4
        mulps  %xmm0,%xmm4
        mulps  %xmm0,%xmm4      ## fscal

        movaps nb113_dxO(%esp),%xmm0
        movaps nb113_dyO(%esp),%xmm1
        movaps nb113_dzO(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2 ## xmm0-xmm2 now contains tx-tz (partial force)

        movss  nb113_fixO(%esp),%xmm3
        movss  nb113_fiyO(%esp),%xmm4
        movss  nb113_fizO(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb113_fixO(%esp)
        movss  %xmm4,nb113_fiyO(%esp)
        movss  %xmm5,nb113_fizO(%esp)   ## updated the O force now do the H's

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $0x39,%xmm3,%xmm3 ## shift right 
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb113_fixH1(%esp),%xmm3
        addss  nb113_fiyH1(%esp),%xmm4
        addss  nb113_fizH1(%esp),%xmm5
        movss  %xmm3,nb113_fixH1(%esp)
        movss  %xmm4,nb113_fiyH1(%esp)
        movss  %xmm5,nb113_fizH1(%esp)          ## updated the H1 force 

        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb113_fixH2(%esp),%xmm3
        addss  nb113_fiyH2(%esp),%xmm4
        addss  nb113_fizH2(%esp),%xmm5
        movss  %xmm3,nb113_fixH2(%esp)
        movss  %xmm4,nb113_fiyH2(%esp)
        movss  %xmm5,nb113_fizH2(%esp)          ## updated the H2 force 

        movl nb113_faction(%ebp),%edi
        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb113_fixM(%esp),%xmm3
        addss  nb113_fiyM(%esp),%xmm4
        addss  nb113_fizM(%esp),%xmm5
        movss  %xmm3,nb113_fixM(%esp)
        movss  %xmm4,nb113_fiyM(%esp)
        movss  %xmm5,nb113_fizM(%esp)   ## updated the M force 

        ## the fj's - move in from mem start by acc. tx/ty/tz in xmm0, xmm1
        movlps (%edi,%eax,4),%xmm6
        movss  8(%edi,%eax,4),%xmm7

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

        subps    %xmm0,%xmm6
        subss    %xmm2,%xmm7

        movlps %xmm6,(%edi,%eax,4)
        movss  %xmm7,8(%edi,%eax,4)

        decl nb113_innerk(%esp)
        jz    _nb_kernel113_ia32_sse.nb113_updateouterdata
        jmp   _nb_kernel113_ia32_sse.nb113_odd_loop
_nb_kernel113_ia32_sse.nb113_updateouterdata: 
        movl  nb113_ii3(%esp),%ecx
        movl  nb113_faction(%ebp),%edi
        movl  nb113_fshift(%ebp),%esi
        movl  nb113_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb113_fixO(%esp),%xmm0
        movaps nb113_fiyO(%esp),%xmm1
        movaps nb113_fizO(%esp),%xmm2

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
        movaps nb113_fixH1(%esp),%xmm0
        movaps nb113_fiyH1(%esp),%xmm1
        movaps nb113_fizH1(%esp),%xmm2

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
        movaps nb113_fixH2(%esp),%xmm0
        movaps nb113_fiyH2(%esp),%xmm1
        movaps nb113_fizH2(%esp),%xmm2

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
        movaps nb113_fixM(%esp),%xmm0
        movaps nb113_fiyM(%esp),%xmm1
        movaps nb113_fizM(%esp),%xmm2

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
        movl nb113_n(%esp),%esi
        ## get group index for i particle 
        movl  nb113_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb113_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb113_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb113_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb113_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb113_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel113_ia32_sse.nb113_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb113_n(%esp)
        jmp _nb_kernel113_ia32_sse.nb113_outer
_nb_kernel113_ia32_sse.nb113_outerend: 
        ## check if more outer neighborlists remain
        movl  nb113_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel113_ia32_sse.nb113_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel113_ia32_sse.nb113_threadloop
_nb_kernel113_ia32_sse.nb113_end: 
        emms

        movl nb113_nouter(%esp),%eax
        movl nb113_ninner(%esp),%ebx
        movl nb113_outeriter(%ebp),%ecx
        movl nb113_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb113_salign(%esp),%eax
        addl %eax,%esp
        addl $924,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



.globl nb_kernel113nf_ia32_sse
.globl _nb_kernel113nf_ia32_sse
nb_kernel113nf_ia32_sse:        
_nb_kernel113nf_ia32_sse:       
.set nb113nf_p_nri, 8
.set nb113nf_iinr, 12
.set nb113nf_jindex, 16
.set nb113nf_jjnr, 20
.set nb113nf_shift, 24
.set nb113nf_shiftvec, 28
.set nb113nf_fshift, 32
.set nb113nf_gid, 36
.set nb113nf_pos, 40
.set nb113nf_faction, 44
.set nb113nf_charge, 48
.set nb113nf_p_facel, 52
.set nb113nf_p_krf, 56
.set nb113nf_p_crf, 60
.set nb113nf_Vc, 64
.set nb113nf_type, 68
.set nb113nf_p_ntype, 72
.set nb113nf_vdwparam, 76
.set nb113nf_Vvdw, 80
.set nb113nf_p_tabscale, 84
.set nb113nf_VFtab, 88
.set nb113nf_invsqrta, 92
.set nb113nf_dvda, 96
.set nb113nf_p_gbtabscale, 100
.set nb113nf_GBtab, 104
.set nb113nf_p_nthreads, 108
.set nb113nf_count, 112
.set nb113nf_mtx, 116
.set nb113nf_outeriter, 120
.set nb113nf_inneriter, 124
.set nb113nf_work, 128
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
.set nb113nf_two, 368
.set nb113nf_half, 384
.set nb113nf_three, 400
.set nb113nf_is3, 416
.set nb113nf_ii3, 420
.set nb113nf_ntia, 424
.set nb113nf_innerjjnr, 428
.set nb113nf_innerk, 432
.set nb113nf_n, 436
.set nb113nf_nn1, 440
.set nb113nf_nri, 444
.set nb113nf_nouter, 448
.set nb113nf_ninner, 452
.set nb113nf_salign, 456
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $460,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb113nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb113nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb113nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb113nf_nouter(%esp)
        movl %eax,nb113nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb113nf_half(%esp)
        movss nb113nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb113nf_half(%esp)
        movaps %xmm2,nb113nf_two(%esp)
        movaps %xmm3,nb113nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb113nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb113nf_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm4
        movss 12(%edx,%ebx,4),%xmm3
        movl nb113nf_p_facel(%ebp),%esi
        movss (%esi),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb113nf_iqM(%esp)
        movaps %xmm4,nb113nf_iqH(%esp)

        movl  nb113nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb113nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb113nf_ntia(%esp)

_nb_kernel113nf_ia32_sse.nb113nf_threadloop: 
        movl  nb113nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel113nf_ia32_sse.nb113nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel113nf_ia32_sse.nb113nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb113nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb113nf_n(%esp)
        movl %ebx,nb113nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel113nf_ia32_sse.nb113nf_outerstart
        jmp _nb_kernel113nf_ia32_sse.nb113nf_end

_nb_kernel113nf_ia32_sse.nb113nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb113nf_nouter(%esp),%ebx
        movl %ebx,nb113nf_nouter(%esp)

_nb_kernel113nf_ia32_sse.nb113nf_outer: 
        movl  nb113nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb113nf_is3(%esp)            ## store is3 

        movl  nb113nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb113nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb113nf_pos(%ebp),%eax    ## eax = base of pos[]  
        movl  %ebx,nb113nf_ii3(%esp)

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
        movaps %xmm3,nb113nf_ixO(%esp)
        movaps %xmm4,nb113nf_iyO(%esp)
        movaps %xmm5,nb113nf_izO(%esp)
        movaps %xmm6,nb113nf_ixH1(%esp)
        movaps %xmm7,nb113nf_iyH1(%esp)

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
        movaps %xmm6,nb113nf_izH1(%esp)
        movaps %xmm0,nb113nf_ixH2(%esp)
        movaps %xmm1,nb113nf_iyH2(%esp)
        movaps %xmm2,nb113nf_izH2(%esp)
        movaps %xmm3,nb113nf_ixM(%esp)
        movaps %xmm4,nb113nf_iyM(%esp)
        movaps %xmm5,nb113nf_izM(%esp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb113nf_vctot(%esp)
        movaps %xmm4,nb113nf_Vvdwtot(%esp)

        movl  nb113nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb113nf_pos(%ebp),%esi
        movl  nb113nf_faction(%ebp),%edi
        movl  nb113nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb113nf_innerjjnr(%esp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb113nf_ninner(%esp),%ecx
        movl  %ecx,nb113nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb113nf_innerk(%esp)         ## number of innerloop atoms 
        jge   _nb_kernel113nf_ia32_sse.nb113nf_unroll_loop
        jmp   _nb_kernel113nf_ia32_sse.nb113nf_odd_inner
_nb_kernel113nf_ia32_sse.nb113nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb113nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb113nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb113nf_charge(%ebp),%esi  ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb113nf_iqM(%esp),%xmm3
        mulps  nb113nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb113nf_qqM(%esp)
        movaps  %xmm4,nb113nf_qqH(%esp)

        movl nb113nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb113nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb113nf_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx
        addl %edi,%ecx
        addl %edi,%edx

        movlps (%esi,%eax,4),%xmm6
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm6
        movhps (%esi,%edx,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm6 ## constant 11011101

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        movaps %xmm4,nb113nf_c6(%esp)
        movaps %xmm6,nb113nf_c12(%esp)

        movl nb113nf_pos(%ebp),%esi     ## base of pos[] 

        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx        ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move four coordinates to xmm0-xmm2   
        movlps (%esi,%eax,4),%xmm4
        movlps (%esi,%ecx,4),%xmm5
        movss 8(%esi,%eax,4),%xmm2
        movss 8(%esi,%ecx,4),%xmm6

        movhps (%esi,%ebx,4),%xmm4
        movhps (%esi,%edx,4),%xmm5

        movss 8(%esi,%ebx,4),%xmm0
        movss 8(%esi,%edx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## constant 10001000

        shufps $136,%xmm5,%xmm0 ## constant 10001000
        shufps $221,%xmm5,%xmm1 ## constant 11011101            

        ## move ixO-izO to xmm4-xmm6 
        movaps nb113nf_ixO(%esp),%xmm4
        movaps nb113nf_iyO(%esp),%xmm5
        movaps nb113nf_izO(%esp),%xmm6

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
        movaps nb113nf_ixH1(%esp),%xmm4
        movaps nb113nf_iyH1(%esp),%xmm5
        movaps nb113nf_izH1(%esp),%xmm6

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
        movaps nb113nf_ixH2(%esp),%xmm3
        movaps nb113nf_iyH2(%esp),%xmm4
        movaps nb113nf_izH2(%esp),%xmm5

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
        movaps nb113nf_iyM(%esp),%xmm3
        movaps nb113nf_izM(%esp),%xmm4
        subps  %xmm1,%xmm3
        subps  %xmm2,%xmm4
        movaps nb113nf_ixM(%esp),%xmm2
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
        movaps  nb113nf_three(%esp),%xmm0
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb113nf_half(%esp),%xmm0
        movaps  %xmm0,nb113nf_rinvH1(%esp)      ## rinvH1 

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb113nf_three(%esp),%xmm0
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb113nf_half(%esp),%xmm0
        movaps  %xmm0,nb113nf_rinvH2(%esp)      ## rinvH2 

        ## rsqM - seed to xmm2 
        rsqrtps %xmm4,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb113nf_three(%esp),%xmm0
        mulps   %xmm4,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb113nf_half(%esp),%xmm0
        movaps  %xmm0,nb113nf_rinvM(%esp)

        ## Do the O LJ-only interaction directly.       
        rcpps   %xmm7,%xmm2
        movaps  nb113nf_two(%esp),%xmm1
        mulps   %xmm2,%xmm7
        subps   %xmm7,%xmm1
        mulps   %xmm1,%xmm2 ## rinvsq 
        movaps  %xmm2,%xmm0
        mulps   %xmm2,%xmm0     ## r4
        mulps   %xmm2,%xmm0     ## r6
        movaps  %xmm0,%xmm1
        mulps   %xmm1,%xmm1     ## r12
        mulps   nb113nf_c6(%esp),%xmm0
        mulps   nb113nf_c12(%esp),%xmm1
        movaps  %xmm1,%xmm3
        subps   %xmm0,%xmm3     ## Vvdw12-Vvdw6
        addps   nb113nf_Vvdwtot(%esp),%xmm3
        movaps  %xmm3,nb113nf_Vvdwtot(%esp)

        ## Do H1 interaction
        movaps  nb113nf_rinvH1(%esp),%xmm7
        mulps  nb113nf_qqH(%esp),%xmm7          ## xmm7=vcoul 
        addps  nb113nf_vctot(%esp),%xmm7

        ## Done with H1, do H2 interactions
        movaps  nb113nf_rinvH2(%esp),%xmm6
        mulps  nb113nf_qqH(%esp),%xmm6          ## xmm6=vcoul 
        addps  %xmm7,%xmm6

        ## Done with H2, do M interactions
        movaps  nb113nf_rinvM(%esp),%xmm5
        mulps  nb113nf_qqM(%esp),%xmm5          ## xmm5=vcoul 
        addps  %xmm6,%xmm5
        movaps %xmm5,nb113nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb113nf_innerk(%esp)
        jl    _nb_kernel113nf_ia32_sse.nb113nf_odd_inner
        jmp   _nb_kernel113nf_ia32_sse.nb113nf_unroll_loop
_nb_kernel113nf_ia32_sse.nb113nf_odd_inner: 
        addl $4,nb113nf_innerk(%esp)
        jnz   _nb_kernel113nf_ia32_sse.nb113nf_odd_loop
        jmp   _nb_kernel113nf_ia32_sse.nb113nf_updateouterdata
_nb_kernel113nf_ia32_sse.nb113nf_odd_loop: 
        movl  nb113nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb113nf_innerjjnr(%esp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb113nf_iqM(%esp),%xmm4
        movl nb113nf_charge(%ebp),%esi
        movhps nb113nf_iqH(%esp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%esi,%eax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb113nf_qqM(%esp)          ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movl nb113nf_type(%ebp),%esi
        movl (%esi,%eax,4),%ebx
        movl nb113nf_vdwparam(%ebp),%esi
        shll %ebx
        addl nb113nf_ntia(%esp),%ebx
        movlps (%esi,%ebx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb113nf_c6(%esp)
        movaps %xmm7,nb113nf_c12(%esp)

        movl nb113nf_pos(%ebp),%esi
        leal (%eax,%eax,2),%eax

        movss nb113nf_ixO(%esp),%xmm3
        movss nb113nf_iyO(%esp),%xmm4
        movss nb113nf_izO(%esp),%xmm5
        movss nb113nf_ixH1(%esp),%xmm0
        movss nb113nf_iyH1(%esp),%xmm1
        movss nb113nf_izH1(%esp),%xmm2
        unpcklps nb113nf_ixH2(%esp),%xmm3       ## ixO ixH2 - -
        unpcklps nb113nf_iyH2(%esp),%xmm4       ## iyO iyH2 - -
        unpcklps nb113nf_izH2(%esp),%xmm5       ## izO izH2 - -
        unpcklps nb113nf_ixM(%esp),%xmm0        ## ixH1 ixM - -
        unpcklps nb113nf_iyM(%esp),%xmm1        ## iyH1 iyM - -
        unpcklps nb113nf_izM(%esp),%xmm2        ## izH1 izM - -
        unpcklps %xmm0,%xmm3    ## ixO ixH1 ixH2 ixM
        unpcklps %xmm1,%xmm4    ## same for y
        unpcklps %xmm2,%xmm5    ## same for z

        ## move j coords to xmm0-xmm2 
        movss (%esi,%eax,4),%xmm0
        movss 4(%esi,%eax,4),%xmm1
        movss 8(%esi,%eax,4),%xmm2
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
        movaps nb113nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb113nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv

        movaps %xmm0,%xmm1
        mulps %xmm1,%xmm1       ## rinvsq 
        movaps %xmm0,%xmm7
        mulps  nb113nf_qqM(%esp),%xmm7   ## vcoul

        addps  nb113nf_vctot(%esp),%xmm7
        movaps %xmm7,nb113nf_vctot(%esp)

        movaps %xmm1,%xmm2
        mulss  %xmm1,%xmm1
        mulss  %xmm2,%xmm1      ## xmm1=rinvsix
        xorps  %xmm4,%xmm4
        movss  %xmm1,%xmm4
        mulss  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulss  nb113nf_c6(%esp),%xmm1
        mulss  nb113nf_c12(%esp),%xmm4
        movaps %xmm4,%xmm3
        subss  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addss  nb113nf_Vvdwtot(%esp),%xmm3
        movss  %xmm3,nb113nf_Vvdwtot(%esp)

        decl nb113nf_innerk(%esp)
        jz    _nb_kernel113nf_ia32_sse.nb113nf_updateouterdata
        jmp   _nb_kernel113nf_ia32_sse.nb113nf_odd_loop
_nb_kernel113nf_ia32_sse.nb113nf_updateouterdata: 
        ## get n from stack
        movl nb113nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb113nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb113nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb113nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb113nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb113nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb113nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel113nf_ia32_sse.nb113nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb113nf_n(%esp)
        jmp _nb_kernel113nf_ia32_sse.nb113nf_outer
_nb_kernel113nf_ia32_sse.nb113nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb113nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel113nf_ia32_sse.nb113nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel113nf_ia32_sse.nb113nf_threadloop
_nb_kernel113nf_ia32_sse.nb113nf_end: 
        emms

        movl nb113nf_nouter(%esp),%eax
        movl nb113nf_ninner(%esp),%ebx
        movl nb113nf_outeriter(%ebp),%ecx
        movl nb113nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb113nf_salign(%esp),%eax
        addl %eax,%esp
        addl $460,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


